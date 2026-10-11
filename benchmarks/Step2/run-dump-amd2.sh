#!/usr/bin/env bash
# after the timing runs: dumps of the penny meshes (beta = 30, h0/8) for the crack-opening figure
set -uo pipefail
cd "$(dirname "$0")/../src"
until [ -f ../benchmarks/logs_final/cutpar_finished.txt ]; do sleep 60; done
L=../benchmarks/logs_final; export STEP2_DATA=$(cd ../benchmarks && pwd)/data_final/dump; mkdir -p $STEP2_DATA
./run-maxperf.sh penny 20 8 30 0.2 0.075 dump=$STEP2_DATA/penny_20_8_30 > $L/dump_20_8_30.log 2>&1
date > $L/dump_finished.txt
