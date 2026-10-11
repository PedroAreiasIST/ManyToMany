#!/usr/bin/env bash
# Energy balance on the geometry of case A1 (a/c = 0.6, a/t = 0.2): two fronts (a = 0.12 and the advance of 0.03), dPi = (f.u1 - f.u0)/2
# against the integral of G over the swept area (figs/src/surf/analyze_energy.py). Accuracy run: niced, heap capped.
cd ~/Desktop/NovosPapers/remesh/step2/src
L=../benchmarks/logs_surf; mkdir -p $L
export STEP2_DATA=$(cd ../benchmarks && pwd)/data_surf_energy; mkdir -p $STEP2_DATA
export DOTNET_GCHeapHardLimit=0x700000000 STEP2_QPLATEAU=0.5
for d in 8 12; do
  nice -n 10 ./run-maxperf.sh surface 20 $d 0.12 0.20 1 0.03 3 0.5 0.04 48 t=0.6 h=1 > $L/energy_A1_div$d.log 2>&1; echo "energy_A1_div$d rc=$?" >> $L/rc_energy.txt
done
date > $L/energy_finished.txt
