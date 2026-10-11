#!/usr/bin/env bash
# finest penny fronts (h0/32) on amd2 after the other runs (Deucalion 1993350/1993351 may finish first)
set -uo pipefail
cd "$(dirname "$0")/../src"
until [ -f ../benchmarks/logs_final/dump_finished.txt ]; do sleep 60; done
L=../benchmarks/logs_final; export STEP2_DATA=$(cd ../benchmarks && pwd)/data_final; mkdir -p $STEP2_DATA
for b in 0 30; do /usr/bin/time -v ./run-maxperf.sh penny 20 32 $b 0.2 0.075 > $L/penny_20_32_${b}.log 2>&1; done
date > $L/penny32_finished.txt
