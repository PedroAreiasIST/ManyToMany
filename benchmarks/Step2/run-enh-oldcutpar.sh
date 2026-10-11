#!/usr/bin/env bash
# cutpar 55 with the pre-enhancement binary (step2_old) against the current one, same conditions
L=~/Desktop/NovosPapers/remesh/step2/benchmarks/logs_enh/oldcutpar; mkdir -p $L
until [ -f ~/Desktop/NovosPapers/remesh/step2/benchmarks/logs_enh/growth2_finished.txt ]; do sleep 20; done
export STEP2_DATA=$L
for impl in old new; do
  if [ $impl = old ]; then cd ~/step2_old/src; cp ~/Desktop/NovosPapers/remesh/step2/src/run-maxperf.sh .; else cd ~/Desktop/NovosPapers/remesh/step2/src; fi
  ./run-maxperf.sh cutpar 55 4 1,16 5 > $L/cutpar_55_4_$impl.log 2>&1
done
date > $L/finished.txt
