#!/usr/bin/env bash
# Runs the Tables 10-14 re-measurement under the conditions stated in RESULTS.md:
# pinned to the performance cores (0-15 on the i9-13900HX), ASLR off, tiered
# compilation / TieredPGO / ReadyToRun off, server GC with concurrent GC off.
#
#   ./run-maxperf.sh [table10|table11|tables1213|table14|all]
set -euo pipefail
CPUS="${MM3_CPUS:-0-15}"
cd "$(dirname "$0")"
export DOTNET_TieredCompilation=0
export DOTNET_TieredPGO=0
export DOTNET_ReadyToRun=0
export DOTNET_TC_QuickJitForLoops=0
export DOTNET_gcServer=1
export DOTNET_gcConcurrent=0
export DOTNET_GCRetainVM=1
echo "pinned to CPUs $CPUS; governor: $(cat /sys/devices/system/cpu/cpu0/cpufreq/scaling_governor 2>/dev/null || echo unknown)" >&2
exec taskset -c "$CPUS" setarch -R dotnet run -c Release --no-build -- "$@"
