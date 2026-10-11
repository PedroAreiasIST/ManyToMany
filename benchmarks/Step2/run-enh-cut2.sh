#!/usr/bin/env bash
# cut-dependent measurements again after the O(1) row-limit check (the per-node bound cost 17 % of the
# sequential cut at 1e6), then every gate
set -uo pipefail
cd ~/Desktop/NovosPapers/remesh/step2/src
L=../benchmarks/logs_enh/cut2; mkdir -p $L
D=$(cd ../benchmarks && pwd)/data_enh/cut2; mkdir -p $D
export STEP2_DATA=$D
R=./run-maxperf.sh
for n in 55 119; do for g in 4 8; do $R cutpar $n $g 1,2,4,8,16 5 > $L/cutpar_${n}_${g}.log 2>&1; done; done
for n in 55 119; do for g in 4 8; do $R batchcut $n $g 3 > $L/batchcut_${n}_${g}.log 2>&1; done; done
for g in 2 3; do $R interfaces 32 $g 16 $D/interfaces_32_$g > $L/interfaces_32_$g.log 2>&1; done
mkdir -p $D/growth; STEP2_FRONTDBG=1 STEP2_DATA=$D/growth $R growth 20 8 0.25 0.15 14 0.025 4 0.04 0.06 > $L/growth.log 2>&1
date > $L/measured.txt
G=$L/gates.log; : > $G
for a in "3 3000 5 12345" "2 5000 3 1" "4 2000 20 7" "5 1000 50 99" "3 20000 2 2024"; do $R verify-fuzz $a 2>&1 | grep -E "^verify|Exception|Mismatch" >> $G; done
for n in 4 8 12; do $R verify-work $n 2>&1 | grep -E "^verify|Exception|Mismatch" >> $G; done
$R verify-cut 4 9 1 2>&1 | grep -E "^verify|Exception" >> $G; for s in 1 2 3 4 5; do $R verify-cut 5 12 $s 2>&1 | grep -E "^verify|Exception" >> $G; done
for s in 1 2 3; do $R verify-uncut 4 10 $s 2>&1 | grep -E "^verify|Exception" >> $G; done
for c in 0 1 2; do for n in 8 13 21; do $R verify-crack $n 10 $c 2>&1 | grep -E "^verify|Exception" >> $G; done; done
for a in "16 12 1 8" "16 12 2 16" "20 10 3 4" "24 10 4 16" "10 15 5 2"; do $R verify-cutpar $a 2>&1 | grep -E "^verify|Exception" >> $G; done
for n in 30 60; do $R verify-flatbuild $n 2>&1 | grep -E "^verify|Exception" >> $G; STEP2_CPUS=0-1 $R verify-flatbuild $n 2>&1 | grep -E "^verify|Exception" >> $G; done
for n in 4 6 9; do $R verify-rowlimit $n 30 2>&1 | grep -E "^verify|Exception" >> $G; done
$R penny 20 8 0 0.2 0.075 patch > $L/penny_patch.log 2>&1
date > $L/finished.txt
