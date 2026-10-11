#!/usr/bin/env bash
# 2026-10-08 enhancements: smoke gates (flat build, row limit)
cd ~/Desktop/NovosPapers/remesh/step2/src
L=../benchmarks/logs_enh
export DOTNET_gcServer=1
D=~/.dotnet/dotnet; B=Step2.Bench/bin/Release/net9.0/Step2.Bench.dll
$D $B verify-flatbuild 30 > $L/flatbuild_30.log 2>&1
taskset -c 0-1 $D $B verify-flatbuild 30 > $L/flatbuild_30_2cpu.log 2>&1
$D $B verify-rowlimit 6 20 > $L/rowlimit_6.log 2>&1
date > $L/smoke_finished.txt
