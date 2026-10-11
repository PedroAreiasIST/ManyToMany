#!/usr/bin/env bash
# 2026-10-11: follow-up of run-fin4-timed-amd2.sh (waits for it): interleaved queries and edits with the collector pauses recorded,
# concurrent readers with and without the writer (STEP2_NOWRITER=1 control of the per-reader rate), the fuzz gates with the destination check.
set -uo pipefail
cd "$(dirname "$0")/../src"
L=../benchmarks/logs_fin5; mkdir -p $L
D=$(cd ../benchmarks && pwd)/data_fin5; mkdir -p $D
export STEP2_DATA=$D
FOREIGN="SimPlas|xfortcom|ifx|gfortran|VBCSCompiler|simplas-amd/jobs"
foreign() { pgrep -f "$FOREIGN" > /dev/null; }
wait_idle() {
  local ok=0
  while [ $ok -lt 2 ]; do
    if foreign || awk -v l="$(cut -d' ' -f1 /proc/loadavg)" 'BEGIN{exit !(l > 1.0)}'; then ok=0; else ok=$((ok+1)); fi
    sleep 30
  done
}
R() {
  local try rc mon
  for try in 1 2 3 4 5 6; do
    wait_idle
    rm -f /tmp/step2_foreign.$$
    ( while true; do foreign && touch /tmp/step2_foreign.$$; sleep 5; done ) &
    mon=$!
    ./run-maxperf.sh "$@"; rc=$?
    kill $mon 2>/dev/null; wait $mon 2>/dev/null
    if [ ! -f /tmp/step2_foreign.$$ ]; then return $rc; fi
    echo "$(date) foreign load during: $*; repeating" >> $L/interference.log
  done
  return $rc
}
until [ -f ../benchmarks/logs_fin4/finished.txt ] || [ -f ../benchmarks/logs_fin4/failed.txt ]; do sleep 60; done
date > $L/started.txt
~/.dotnet/dotnet build Step2.Bench -c Release > $L/build.log 2>&1 || { echo "build failed" > $L/failed.txt; exit 1; }
~/.dotnet/dotnet build-server shutdown > /dev/null 2>&1
R interleaved 6,10,14,18,22,26    > $L/interleaved.log 2>&1
R concurrent 55 3 1,2,4,8,15      > $L/concurrent_55.log 2>&1
mkdir -p $D/nowriter; STEP2_NOWRITER=1 STEP2_DATA=$D/nowriter R concurrent 55 3 1,2,4,8,15 > $L/concurrent_55_nowriter.log 2>&1
R verify-fuzz 3 3000 5 12345 >  $L/fuzz.log 2>&1
R verify-fuzz 2 5000 3 1     >> $L/fuzz.log 2>&1
R verify-fuzz 4 2000 20 7    >> $L/fuzz.log 2>&1
R verify-fuzz 5 1000 50 99   >> $L/fuzz.log 2>&1
R verify-fuzz 3 20000 2 2024 >> $L/fuzz.log 2>&1
date > $L/finished.txt
