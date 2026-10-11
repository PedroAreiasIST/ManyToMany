#!/usr/bin/env bash
# Surface-crack K against Raju and Newman (NASA TP-1578 Table I), coarser meshes (the finest run on Deucalion).
# Waits for run-surface-conv-amd2.sh; accuracy runs (niced, heap capped, exit codes logged).
cd ~/Desktop/NovosPapers/remesh/step2/src
L=../benchmarks/logs_surf; mkdir -p $L
until [ -f $L/conv_finished.txt ]; do sleep 30; done
export STEP2_DATA=$(cd ../benchmarks && pwd)/data_surf_verify; mkdir -p $STEP2_DATA
export DOTNET_GCHeapHardLimit=0x700000000
export STEP2_QPLATEAU=0.5   # tube weight q = 1 up to 0.5 rq
R="nice -n 10 ./run-maxperf.sh"
run() { local name=$1; shift; $R surface "$@" > $L/$name.log 2>&1; echo "$name rc=$?" >> $L/rc_verify.txt; }
for d in 6 8 10; do run A1_div$d 20 $d 0.12 0.20 0 0.02 3 0.04 0.04 48 t=0.6 h=1; done
for d in 4 6 8; do run A2_div$d 20 $d 0.20 0.20 0 0.02 3 0.04 0.065 48 t=0.5 h=1; done
date > $L/verify_finished.txt
