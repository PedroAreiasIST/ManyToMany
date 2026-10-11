#!/usr/bin/env bash
# K verification (penny cracks) on amd2: Sneddon (beta = 0) and Kassir-Sih (beta = 45) under front refinement
# (div 4 and 8 at beta = 0 and the patch test were run before the parallel interaction integral; same results)
set -uo pipefail
cd "$(dirname "$0")/../src"
L=../benchmarks/logs_cut
export STEP2_DATA=$(cd ../benchmarks && pwd)/data_cut
mkdir -p $L $STEP2_DATA
R=./run-maxperf.sh
$R penny 20 16 0 0.2 0.075  > $L/penny_20_16_0.log 2>&1
for d in 4 8 16; do $R penny 20 $d 45 0.2 0.075 > $L/penny_20_${d}_45.log 2>&1; done
for rq in 0.05 0.11; do $R penny 20 16 0 0.2 $rq > $L/penny_20_16_0_rq${rq}.log 2>&1; done
date > $L/penny_finished.txt
