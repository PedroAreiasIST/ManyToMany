#!/usr/bin/env bash
# 2026-10-11: collector-budget study for the cutting timings (amd2, idle machine). The sequential and the concurrent cut allocate
# several GB; the stop-the-world pauses of the server collector depend on the budget of generation 0 and on what else the process
# holds. Records time, collections and pauses of the same operations for three budgets (the protocol default is 0x20000000 = 512 MB).
set -uo pipefail
cd "$(dirname "$0")/../src"
L=../benchmarks/logs_gcstudy; mkdir -p $L
export STEP2_DATA=$(cd ../benchmarks && pwd)/data_gcstudy; mkdir -p $STEP2_DATA
until [ -f ~/step2diag/logs/finished.txt ]; do sleep 20; done
date > $L/started.txt
~/.dotnet/dotnet build Step2.Bench -c Release > $L/build.log 2>&1 || { echo "build failed" > $L/failed.txt; exit 1; }
~/.dotnet/dotnet build-server shutdown > /dev/null 2>&1
for b in 20000000 100000000 400000000; do
  export DOTNET_GCgen0size_OVERRIDE=$b
  # run-maxperf.sh sets DOTNET_GCgen0size itself; a modified copy takes the budget from STEP2_GEN0
  sed "s/DOTNET_GCgen0size=20000000/DOTNET_GCgen0size=\${STEP2_GEN0:-20000000}/" run-maxperf.sh > run-maxperf-gen0.sh; chmod +x run-maxperf-gen0.sh
  STEP2_GEN0=$b STEP2_DATA=$STEP2_DATA/gen0_$b ./run-maxperf-gen0.sh cutgc 119 4 3 > $L/cutgc_119_4_gen0_$b.log 2>&1
  STEP2_GEN0=$b STEP2_DATA=$STEP2_DATA/gen0_$b ./run-maxperf-gen0.sh cutpar 119 4 1,2,4,8,16 5 > $L/cutpar_119_4_gen0_$b.log 2>&1
done
date > $L/finished.txt
