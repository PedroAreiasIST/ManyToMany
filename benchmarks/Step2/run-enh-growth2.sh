#!/usr/bin/env bash
# growth again with the deterministic K sum in blocks of 64 elements (first pass used 2048: K 1.5x slower)
cd ~/Desktop/NovosPapers/remesh/step2/src
L=../benchmarks/logs_enh
until [ -f $L/repeat/finished2.txt ]; do sleep 20; done
D=$(cd ../benchmarks && pwd)/data_enh/growth2; mkdir -p $D
STEP2_FRONTDBG=1 STEP2_DATA=$D ./run-maxperf.sh growth 20 8 0.25 0.15 14 0.025 4 0.04 0.06 > $L/growth2.log 2>&1
date > $L/growth2_finished.txt
