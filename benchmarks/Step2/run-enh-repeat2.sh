#!/usr/bin/env bash
# repeat of the 1e6 concurrent-cut cases (sequential reference 16 % above the batchcut run of the same marks)
cd ~/Desktop/NovosPapers/remesh/step2/src
L=../benchmarks/logs_enh/repeat; mkdir -p $L
until [ -f $L/finished.txt ]; do sleep 20; done
export STEP2_DATA=$(cd ../benchmarks && pwd)/data_enh/repeat_cutpar; mkdir -p $STEP2_DATA
for g in 4 8; do ./run-maxperf.sh cutpar 55 $g 1,2,4,8,16 5 > $L/cutpar_55_$g.log 2>&1; done
date > $L/finished2.txt
