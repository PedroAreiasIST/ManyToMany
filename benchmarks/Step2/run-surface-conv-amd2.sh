#!/usr/bin/env bash
# Surface crack (a = 0.15, c = 0.30, a/c = 0.5, a/t = 0.15, c/b = 0.3) against Newman-Raju under front refinement:
# one front, no growth. Accuracy runs (not timings): niced, serial, heap capped, exit codes logged.
cd ~/Desktop/NovosPapers/remesh/step2/src
L=../benchmarks/logs_surf; mkdir -p $L
export STEP2_DATA=$(cd ../benchmarks && pwd)/data_surf; mkdir -p $STEP2_DATA
export DOTNET_GCHeapHardLimit=0x700000000
R="nice -n 10 ./run-maxperf.sh"
run() { # name args...
  local name=$1; shift
  $R surface "$@" > $L/$name.log 2>&1; echo "$name rc=$?" >> $L/rc.txt
}
for d in 8 12 16; do run conv_div$d 20 $d 0.15 0.30 0 0.02 3 0.04 0.05 48; mv -f $STEP2_DATA/surface_K_20_$d.csv $STEP2_DATA/conv_K_20_${d}_rq0.050.csv; done
for rq in 0.035 0.07; do run rq_div8_$rq 20 8 0.15 0.30 0 0.02 3 0.04 $rq 48; mv -f $STEP2_DATA/surface_K_20_8.csv $STEP2_DATA/conv_K_20_8_rq$rq.csv; done
# the domain form alone (no free-surface term), same mesh as div 8
run nofs_div8 20 8 0.15 0.30 0 0.02 3 0.04 0.05 48 nofs; mv -f $STEP2_DATA/surface_K_20_8.csv $STEP2_DATA/conv_K_20_8_nofs.csv
date > $L/conv_finished.txt
