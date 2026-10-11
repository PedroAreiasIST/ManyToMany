#!/usr/bin/env bash
# repeats of the 1e7 crack workload (the first enh pass showed a p90 outlier), after run-enh-amd2.sh
cd ~/Desktop/NovosPapers/remesh/step2/src
L=../benchmarks/logs_enh/repeat; mkdir -p $L
until [ -f ../benchmarks/logs_enh/finished.txt ]; do sleep 20; done
for r in 1 2; do
  mkdir -p ../benchmarks/data_enh/repeat$r
  STEP2_DATA=$(cd ../benchmarks && pwd)/data_enh/repeat$r ./run-maxperf.sh work crack 119 0 15 1 > $L/work_crack_119_p_$r.log 2>&1
  STEP2_DATA=$(cd ../benchmarks && pwd)/data_enh/repeat$r ./run-maxperf.sh work crack 119 0 1000000 0 > $L/work_crack_119_t_$r.log 2>&1
done
date > $L/finished.txt
