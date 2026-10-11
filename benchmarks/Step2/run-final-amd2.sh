#!/usr/bin/env bash
# 2026-10-08 final amd2 run: growth smoke, interfaces, penny convergence (beta 0, 30), all correctness gates
set -uo pipefail
cd "$(dirname "$0")/../src"
L=../benchmarks/logs_final; mkdir -p $L
export STEP2_DATA=$(cd ../benchmarks && pwd)/data_final; mkdir -p $STEP2_DATA
R=./run-maxperf.sh
mkdir -p $STEP2_DATA/growth_smoke; STEP2_DATA=$STEP2_DATA/growth_smoke $R growth 10 4 0.25 0.15 4 0.04 4 0.03 0.12 > $L/growth_smoke.log 2>&1
for g in 2 3; do $R interfaces 32 $g 16 $STEP2_DATA/interfaces_32_$g > $L/interfaces_32_$g.log 2>&1; done
$R penny 20 8 0 0.2 0.075 patch > $L/penny_patch.log 2>&1
for d in 4 8 16; do for b in 0 30; do $R penny 20 $d $b 0.2 0.075 > $L/penny_20_${d}_${b}.log 2>&1; done; done
for rq in 0.05 0.11; do for b in 0 30; do $R penny 20 16 $b 0.2 $rq > $L/penny_20_16_${b}_rq${rq}.log 2>&1; done; done
G=$L/gates.log; : > $G
for a in "3 3000 5 12345" "2 5000 3 1" "4 2000 20 7" "5 1000 50 99" "3 20000 2 2024"; do $R verify-fuzz $a 2>&1 | grep -E "^verify|Exception|Mismatch" >> $G; done
for n in 4 8 12; do $R verify-work $n 2>&1 | grep -E "^verify|Exception|Mismatch" >> $G; done
$R verify-cut 4 9 1 2>&1 | grep -E "^verify|Exception" >> $G; for s in 1 2 3 4 5; do $R verify-cut 5 12 $s 2>&1 | grep -E "^verify|Exception" >> $G; done
for s in 1 2 3; do $R verify-uncut 4 10 $s 2>&1 | grep -E "^verify|Exception" >> $G; done
for c in 0 1 2; do for n in 8 13 21; do $R verify-crack $n 10 $c 2>&1 | grep -E "^verify|Exception" >> $G; done; done
for a in "16 12 1 8" "16 12 2 16" "20 10 3 4" "24 10 4 16" "10 15 5 2"; do $R verify-cutpar $a 2>&1 | grep -E "^verify|Exception" >> $G; done
date > $L/finished.txt
