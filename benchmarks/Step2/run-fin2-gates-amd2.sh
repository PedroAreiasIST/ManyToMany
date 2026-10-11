#!/usr/bin/env bash
# 2026-10-09: final run, everything with the final code (initial states with the greedy colouring; (flat pattern copy, threaded
# rebuild with chunked transpose and deterministic colouring, gathered assembly, row-limit check).
# Gates first; every measurement after them, one at a time, on the otherwise idle machine.
set -uo pipefail
cd "$(dirname "$0")/../src"
L=../benchmarks/logs_fin; mkdir -p $L
D=$(cd ../benchmarks && pwd)/data_fin; mkdir -p $D
export STEP2_DATA=$D
GN() { nice -n 19 ./run-maxperf.sh "$@"; }
date > $L/gates_started.txt
G=$L/gates.log; : > $G
for a in "3 3000 5 12345" "2 5000 3 1" "4 2000 20 7" "5 1000 50 99" "3 20000 2 2024"; do GN verify-fuzz $a 2>&1 | grep -E "^verify|Exception|Mismatch" >> $G; done
for n in 4 8 12; do GN verify-work $n 2>&1 | grep -E "^verify|Exception|Mismatch" >> $G; done
GN verify-cut 4 9 1 2>&1 | grep -E "^verify|Exception" >> $G; for s in 1 2 3 4 5; do GN verify-cut 5 12 $s 2>&1 | grep -E "^verify|Exception" >> $G; done
for s in 1 2 3; do GN verify-uncut 4 10 $s 2>&1 | grep -E "^verify|Exception" >> $G; done
for c in 0 1 2; do for n in 8 13 21; do GN verify-crack $n 10 $c 2>&1 | grep -E "^verify|Exception" >> $G; done; done
for a in "16 12 1 8" "16 12 2 16" "20 10 3 4" "24 10 4 16" "10 15 5 2"; do GN verify-cutpar $a 2>&1 | grep -E "^verify|Exception" >> $G; done
for n in 30 60; do GN verify-flatbuild $n 2>&1 | grep -E "^verify|Exception" >> $G; STEP2_CPUS=0-1 GN verify-flatbuild $n 2>&1 | grep -E "^verify|Exception" >> $G; done
for n in 4 6 9; do GN verify-rowlimit $n 30 2>&1 | grep -E "^verify|Exception" >> $G; done
GN penny 20 8 0 0.2 0.075 patch > $L/penny_patch.log 2>&1
date > $L/gates_finished.txt

