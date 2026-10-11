#!/usr/bin/env bash
# penny cracks with the global (V = q N) interaction integral
set -uo pipefail
cd "$(dirname "$0")/../src"
L=../benchmarks/logs_cut/v2; mkdir -p $L
export STEP2_DATA=$(cd ../benchmarks && pwd)/data_cut/v2; mkdir -p $STEP2_DATA
R=./run-maxperf.sh
for d in 8 16; do for b in 45 0; do $R penny 20 $d $b 0.2 0.075 > $L/penny_20_${d}_${b}.log 2>&1; done; done
$R penny 20 16 45 0.2 0.05 > $L/penny_20_16_45_rq0.05.log 2>&1
date > $L/finished.txt
