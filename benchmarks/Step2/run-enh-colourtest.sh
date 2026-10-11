#!/usr/bin/env bash
# why is the 1e6 sequential reference of cutpar 15 % slower than with the old binary? old / new / new with the
# sequential initial build (greedy colouring), two rounds
L=~/Desktop/NovosPapers/remesh/step2/benchmarks/logs_enh/colourtest; mkdir -p $L
until [ -f ~/Desktop/NovosPapers/remesh/step2/benchmarks/logs_enh/cut2/finished.txt ]; do sleep 20; done
export STEP2_DATA=$L
NEW=~/Desktop/NovosPapers/remesh/step2/src; BIN=$NEW/Step2.Bench/bin/Release/net9.0
cp $BIN/Step2.Bench.dll $BIN/Step2.Bench.dll.cut2 && cp ~/step2_seqbuild.dll $BIN/Step2.Bench.dll
for r in 1 2; do
  (cd ~/step2_old/src && ./run-maxperf.sh cutpar 55 4 1 3 > $L/old_$r.log 2>&1)
  (cd $NEW && ./run-maxperf.sh cutpar 55 4 1 3 > $L/new_$r.log 2>&1)
  (cd $NEW && STEP2_SEQBUILD=1 ./run-maxperf.sh cutpar 55 4 1 3 > $L/newseq_$r.log 2>&1)
done
date > $L/finished.txt
