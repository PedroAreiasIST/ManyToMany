#!/usr/bin/env bash
# concurrent-cut timings with the final code (exact slot bound), after run-final2-amd2.sh
set -uo pipefail
cd "$(dirname "$0")/../src"
until [ -f ../benchmarks/logs_final/finished.txt ]; do sleep 60; done
L=../benchmarks/logs_final; export STEP2_DATA=$(cd ../benchmarks && pwd)/data_final/cutpar; mkdir -p $STEP2_DATA
R=./run-maxperf.sh
for n in 55 119; do for g in 4 8; do $R cutpar $n $g 1,2,4,8,16 5 > $L/cutpar_${n}_${g}.log 2>&1; done; done
date > $L/cutpar_finished.txt
