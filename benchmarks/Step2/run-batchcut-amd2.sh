#!/usr/bin/env bash
# batch cutter vs local cut on the Table 6 meshes and marks; crossover rerun with the current code
cd ~/Desktop/NovosPapers/remesh/step2/src
L=../benchmarks/logs_final; export STEP2_DATA=$(cd ../benchmarks && pwd)/data_final/batchcut; mkdir -p $STEP2_DATA
~/.dotnet/dotnet build Step2.Bench -c Release > $L/batchcut_build.log 2>&1 || { echo "build failed" > $L/batchcut_finished.txt; exit 1; }
for n in 55 119; do for g in 4 8; do ./run-maxperf.sh batchcut $n $g 3 > $L/batchcut_${n}_${g}.log 2>&1; done; done
export STEP2_DATA=$(cd ../benchmarks && pwd)/data_final/crossover; mkdir -p $STEP2_DATA
./run-maxperf.sh crossover 55 > $L/crossover_55.log 2>&1
date > $L/batchcut_finished.txt
