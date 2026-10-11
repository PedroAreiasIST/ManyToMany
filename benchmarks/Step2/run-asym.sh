#!/usr/bin/env bash
# exact near-tip fields of modes I, II, III imposed on the cracked mesh: interaction integral alone
set -uo pipefail
cd "$(dirname "$0")/../src"
L=../benchmarks/logs_cut; export STEP2_DATA=$(cd ../benchmarks && pwd)/data_cut
R=./run-maxperf.sh
for m in asym asym2 asym3; do for rq in 0.075 0.04; do $R penny 20 8 0 0.2 $rq $m > $L/${m}_rq${rq}.log 2>&1; done; done
date > $L/asym_finished.txt
