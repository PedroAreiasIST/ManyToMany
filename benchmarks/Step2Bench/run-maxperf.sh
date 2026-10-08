#!/usr/bin/env bash
# Runs the step2 measurements under the conditions of the released MM3 harness
# (benchmarks/Mm3Paper/run-maxperf.sh): pinned to the performance cores (0-15 on
# the i9-13900HX), tiered compilation and ReadyToRun off, server GC with
# concurrent GC off, ASLR off. Output files go to $STEP2_DATA (default ./data).
#
#   ./run-maxperf.sh <mode> [args...]      e.g.  ./run-maxperf.sh verify-work 6
#
# The bench must have been built first:  dotnet build -c Release
set -euo pipefail
CPUS="${STEP2_CPUS:-0-15}"
cd "$(dirname "$0")"
if [ -x "$HOME/.dotnet/dotnet" ]; then export DOTNET_ROOT="$HOME/.dotnet" PATH="$HOME/.dotnet:$PATH"; fi
export DOTNET_TieredCompilation=0
export DOTNET_TieredPGO=0
export DOTNET_ReadyToRun=0
export DOTNET_TC_QuickJitForLoops=0
export DOTNET_gcServer=1
export DOTNET_gcConcurrent=0
export DOTNET_GCgen0size=20000000
export DOTNET_GCRetainVM=1
# the output folder depends on how the bench was built (solution: bin/x64/Release, project: bin/Release)
DLL=$(ls -t bin/*/Release/net9.0/Step2.Bench.dll bin/Release/net9.0/Step2.Bench.dll 2>/dev/null | head -1 || true)
[ -n "$DLL" ] || { echo "Step2.Bench.dll not found: run 'dotnet build -c Release' first" >&2; exit 1; }
echo "pinned to CPUs $CPUS; governor: $(cat /sys/devices/system/cpu/cpu0/cpufreq/scaling_governor 2>/dev/null || echo unknown)" >&2
exec taskset -c "$CPUS" setarch -R dotnet "$DLL" "$@"
