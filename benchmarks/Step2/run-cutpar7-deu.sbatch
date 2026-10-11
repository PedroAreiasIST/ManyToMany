#!/bin/bash
# concurrent cutting on one EPYC 7742 node, final library and harness with collector accounting (bin-gc): times with and without the pauses.
# 16 server-GC heaps (default) and 64 heaps; the sequential reference of the same process is the baseline of every speed-up.
#SBATCH --job-name=step2-cutpar7
#SBATCH --partition=dev-x86
#SBATCH --account=f2026002cncaist_ulx
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=128
#SBATCH --exclusive
#SBATCH --time=03:30:00
#SBATCH --output=/projects/F2026002CNCAIST_UL/step2/logs/cutpar7-%j.out
source /projects/F2026002CNCAIST_UL/step2/env.sh
B=/projects/F2026002CNCAIST_UL/step2/bin-gc/Step2.Bench
export STEP2_DATA=/projects/F2026002CNCAIST_UL/step2/data/cutpar7; mkdir -p $STEP2_DATA
echo "node $(hostname)"
$B cutpar 119 4 1,2,4,8,16,32,64,128 5
$B cutpar 119 8 1,2,4,8,16,32,64,128 5
$B cutgc 119 4 3
$B cutpar 160 8 1,16,32,64,128 3
mkdir -p $STEP2_DATA/gc64; STEP2_DATA=$STEP2_DATA/gc64 DOTNET_GCHeapCount=40 $B cutpar 119 8 1,16,32,64,128 3
mkdir -p $STEP2_DATA/gc64b; STEP2_DATA=$STEP2_DATA/gc64b DOTNET_GCHeapCount=40 $B cutpar 160 8 1,16,32,64,128 3
