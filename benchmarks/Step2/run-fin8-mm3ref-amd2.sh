#!/usr/bin/env bash
# 2026-10-11: control of the per-reader rate of the concurrent readers: the MM3 cache-rebuild protocol with the collector pauses recorded (mm3ref, three repeats).
set -uo pipefail
cd "$(dirname "$0")/../src"
L=../benchmarks/logs_fin8; mkdir -p $L
D=$(cd ../benchmarks && pwd)/data_fin8; mkdir -p $D
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
until [ -f ../benchmarks/logs_fin7/finished.txt ] || [ -f ../benchmarks/logs_fin7/failed.txt ]; do sleep 20; done
date > $L/started.txt
~/.dotnet/dotnet build Step2.Bench -c Release > $L/build.log 2>&1 || { echo "build failed" > $L/failed.txt; exit 1; }
~/.dotnet/dotnet build-server shutdown > /dev/null 2>&1
for k in 1 2 3; do mkdir -p $D/mm3ref_$k; STEP2_DATA=$D/mm3ref_$k R mm3ref > $L/mm3ref_$k.log 2>&1; done
date > $L/finished.txt
