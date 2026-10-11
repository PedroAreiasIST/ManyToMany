#!/usr/bin/env bash
# 2026-10-08 (2): penny rq = 0.11, Paris growth (n = 20, h0/8), all correctness gates on the final code
set -uo pipefail
cd "$(dirname "$0")/../src"
L=../benchmarks/logs_final; mkdir -p $L
export STEP2_DATA=$(cd ../benchmarks && pwd)/data_final; mkdir -p $STEP2_DATA
R=./run-maxperf.sh
for b in 0 30; do $R penny 20 16 $b 0.2 0.11 > $L/penny_20_16_${b}_rq0.11.log 2>&1; done
mkdir -p $STEP2_DATA/growth; STEP2_FRONTDBG=1 STEP2_DATA=$STEP2_DATA/growth $R growth 20 8 0.25 0.15 14 0.025 4 0.04 0.06 > $L/growth.log 2>&1
G=$L/gates.log; : > $G
for a in "3 3000 5 12345" "2 5000 3 1" "4 2000 20 7" "5 1000 50 99" "3 20000 2 2024"; do $R verify-fuzz $a 2>&1 | grep -E "^verify|Exception|Mismatch" >> $G; done
for n in 4 8 12; do $R verify-work $n 2>&1 | grep -E "^verify|Exception|Mismatch" >> $G; done
$R verify-cut 4 9 1 2>&1 | grep -E "^verify|Exception" >> $G; for s in 1 2 3 4 5; do $R verify-cut 5 12 $s 2>&1 | grep -E "^verify|Exception" >> $G; done
for s in 1 2 3; do $R verify-uncut 4 10 $s 2>&1 | grep -E "^verify|Exception" >> $G; done
for c in 0 1 2; do for n in 8 13 21; do $R verify-crack $n 10 $c 2>&1 | grep -E "^verify|Exception" >> $G; done; done
for a in "16 12 1 8" "16 12 2 16" "20 10 3 4" "24 10 4 16" "10 15 5 2"; do $R verify-cutpar $a 2>&1 | grep -E "^verify|Exception" >> $G; done
date > $L/finished.txt
