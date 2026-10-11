#!/usr/bin/env bash
# Measurements of the cutting extension (concurrent cutting) on amd2; output in data_cut/ and logs_cut/
# so that the measured data of the 2026-10-01 run stay untouched.
set -uo pipefail
cd "$(dirname "$0")/../src"
load=$(cut -d' ' -f1 /proc/loadavg)
busy=$(ps -eo pcpu= | awk '$1 > 50 {n++} END {print n+0}')
if awk -v l="$load" 'BEGIN {exit !(l > 2.0)}' || [ "$busy" -gt 0 ]; then
  echo "machine busy (load $load, $busy processes above 50% CPU); not measuring" >&2
  ps -eo pid,pcpu,comm --sort=-pcpu | head -8 >&2
  exit 3
fi
L=../benchmarks/logs_cut
export STEP2_DATA=$(cd ../benchmarks && pwd)/data_cut
mkdir -p $L $STEP2_DATA
R=./run-maxperf.sh
date > $L/started.txt
# correctness gate on this machine before measuring
$R verify-cutpar 16 12 1 16          > $L/verify_cutpar_16.log 2>&1
$R verify-fuzz 3 3000 5 12345        > $L/verify_fuzz.log 2>&1
# concurrent cutting: g^3 groups of sphere-crossing edges
for n in 55 119; do
  for g in 4 8; do
    $R cutpar $n $g 1,2,4,8,16 5     > $L/cutpar_${n}_${g}.log 2>&1
  done
done
# regression of the sequential hot path after the concurrent-mode changes (compare with data/)
$R work crack 55 0 7 1               > $L/work_crack_55_p.log 2>&1
$R work crack 55 0 1000000 0         > $L/work_crack_55_t.log 2>&1
$R work crack 119 0 15 1             > $L/work_crack_119_p.log 2>&1
date > $L/finished.txt
