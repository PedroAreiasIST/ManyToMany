#!/bin/bash
# controls requested by the audits (2026-10-11): (1) finite-body: the same penny crack at half the radius in the same cube (a/W from 0.2 to 0.1) at the same h_min/a;
# (2) tube factor: linear (STEP2_QPLATEAU=0) against the plateau 0.5 for three tube radii, penny beta = 30 deg and beta = 0.
#SBATCH --job-name=step2-pennyctl
#SBATCH --partition=dev-x86
#SBATCH --account=f2026002cncaist_ulx
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=128
#SBATCH --exclusive
#SBATCH --time=03:00:00
#SBATCH --output=/projects/F2026002CNCAIST_UL/step2/logs/pennyctl-%j.out
source /projects/F2026002CNCAIST_UL/step2/env.sh
B=/projects/F2026002CNCAIST_UL/step2/bin-gc/Step2.Bench
export STEP2_DATA=/projects/F2026002CNCAIST_UL/step2/data/pennyctl; mkdir -p $STEP2_DATA
echo "node $(hostname)"
echo "=== finite body: a = 0.1 (h_min/a = 1/div * 1/... see div), cube [-1,1]^3"
for beta in 0 30; do
  for div in 8 16; do echo "--- a=0.1 div=$div beta=$beta"; $B penny 20 $div $beta 0.1 0.0375; done
done
echo "=== tube factor: linear against plateau 0.5, a = 0.2, div 8"
for beta in 30 0; do
  for rq in 0.05 0.075 0.11; do
    for g in 0 0.5; do echo "--- a=0.2 div=8 beta=$beta rq=$rq plateau=$g"; STEP2_QPLATEAU=$g $B penny 20 8 $beta 0.2 $rq; done
  done
done
echo "=== tube factor at div 16 (rq 0.075)"
for beta in 30 0; do echo "--- a=0.2 div=16 beta=$beta rq=0.075 plateau=0.5"; STEP2_QPLATEAU=0.5 $B penny 20 16 $beta 0.2 0.075; done
