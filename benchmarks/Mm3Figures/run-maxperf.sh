#!/usr/bin/env bash
# Runs the paper measurements under maximum-throughput conditions.
#
#   * pinned to the performance cores only (0-15 on a 13900HX): the hybrid
#     scheduler otherwise migrates work onto 3.9 GHz E-cores mid-measurement,
#     which shows up as 30-50% run-to-run swings;
#   * tiered compilation and ReadyToRun disabled so every method is fully
#     optimised from the first call instead of ramping through tier 0;
#   * server GC, concurrent GC off, 512 MB gen0 budget so measured regions are
#     not interrupted by collections;
#   * ASLR disabled (setarch -R) for run-to-run determinism.
#
# For the last few percent, also set the CPU governor to performance
# (needs root, and is NOT done here):
#     sudo cpupower frequency-set -g performance
set -euo pipefail

CPUS="${MM3_CPUS:-0-15}"
cd "$(dirname "$0")"

export DOTNET_TieredCompilation=0
export DOTNET_TieredPGO=0
export DOTNET_ReadyToRun=0
export DOTNET_TC_QuickJitForLoops=0
export DOTNET_gcServer=1
export DOTNET_gcConcurrent=0
export DOTNET_GCgen0size=20000000            # 512 MB gen0 (hex-free decimal is accepted)
export DOTNET_GCRetainVM=1
export DOTNET_TieredCompilation_BackgroundWorkerTimeoutMs=0

echo "pinned to CPUs $CPUS; governor: $(cat /sys/devices/system/cpu/cpu0/cpufreq/scaling_governor 2>/dev/null || echo unknown)"
exec taskset -c "$CPUS" setarch -R dotnet run -c Release --no-build -- "${1:-figdata}"
