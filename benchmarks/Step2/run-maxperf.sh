#!/usr/bin/env bash
# Runs the step2 measurements under the conditions of the released MM3 harness
# (benchmarks/Mm3Paper/run-maxperf.sh): pinned to the performance cores (0-15 on
# the i9-13900HX), tiered compilation and ReadyToRun off, server GC with
# concurrent GC off, ASLR off. Usage: ./run-maxperf.sh <mode> [args...]
set -euo pipefail
CPUS="${STEP2_CPUS:-0-15}"
cd "$(dirname "$0")"
export DOTNET_ROOT="$HOME/.dotnet" PATH="$HOME/.dotnet:$PATH"
export DOTNET_TieredCompilation=0
export DOTNET_TieredPGO=0
export DOTNET_ReadyToRun=0
export DOTNET_TC_QuickJitForLoops=0
export DOTNET_gcServer=1
export DOTNET_gcConcurrent=0
export DOTNET_GCgen0size="${STEP2_GEN0:-20000000}"
export DOTNET_GCRetainVM=1
echo "pinned to CPUs $CPUS; governor: $(cat /sys/devices/system/cpu/cpu0/cpufreq/scaling_governor 2>/dev/null || echo unknown)" >&2
exec taskset -c "$CPUS" setarch -R dotnet Step2.Bench/bin/Release/net9.0/Step2.Bench.dll "$@"
