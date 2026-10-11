#!/usr/bin/env bash
# 2026-10-11: the final timed run (run-fin2-timed-amd2.sh) with the final library (SplitNode, corner crack) on the idle machine.
# Builds first, keeps the previous logs and data, then runs the timed set one measurement at a time.
set -uo pipefail
cd "$(dirname "$0")/../src"
L=../benchmarks/logs_fin; mkdir -p $L
[ -d ../benchmarks/logs_fin_before_20261011 ] || cp -a ../benchmarks/logs_fin ../benchmarks/logs_fin_before_20261011
[ -d ../benchmarks/data_fin_before_20261011 ] || cp -a ../benchmarks/data_fin ../benchmarks/data_fin_before_20261011
date > $L/fin3_started.txt
~/.dotnet/dotnet build Step2.Bench -c Release > $L/build_fin3.log 2>&1 || { echo "build failed" > $L/fin3_failed.txt; exit 1; }
~/.dotnet/dotnet build-server shutdown > /dev/null 2>&1
exec ../benchmarks/run-fin2-timed-amd2.sh
