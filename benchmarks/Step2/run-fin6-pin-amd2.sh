#!/usr/bin/env bash
# 2026-10-11: control of the per-reader rate of the concurrent readers: writer on CPU 0 and reader i on CPU i+1 (STEP2_PINTHREADS=1), against the unpinned run.
set -uo pipefail
cd "$(dirname "$0")/../src"
L=../benchmarks/logs_fin6; mkdir -p $L
D=$(cd ../benchmarks && pwd)/data_fin6; mkdir -p $D
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
until [ -f ../benchmarks/logs_fin5/finished.txt ] || [ -f ../benchmarks/logs_fin5/failed.txt ]; do sleep 20; done
date > $L/started.txt
~/.dotnet/dotnet build Step2.Bench -c Release > $L/build.log 2>&1 || { echo "build failed" > $L/failed.txt; exit 1; }
~/.dotnet/dotnet build-server shutdown > /dev/null 2>&1
mkdir -p $D/pinned $D/unpinned
STEP2_PINTHREADS=1 STEP2_DATA=$D/pinned   R concurrent 55 3 1,2,4,8,15 > $L/concurrent_55_pinned.log 2>&1
                   STEP2_DATA=$D/unpinned R concurrent 55 3 1,2,4,8,15 > $L/concurrent_55_unpinned.log 2>&1
STEP2_PINTHREADS=1 STEP2_NOWRITER=1 STEP2_DATA=$D/pinned_nowriter R concurrent 55 3 1,15 > $L/concurrent_55_pinned_nowriter.log 2>&1
date > $L/finished.txt
