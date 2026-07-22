#!/usr/bin/env bash
# BenchmarkDotNet runs under the same maximum-throughput conditions as Mm3Paper:
# pinned to the performance cores, server GC, concurrent GC off. BDN manages JIT
# warm-up itself, so tiering is left to its own harness configuration.
set -euo pipefail
CPUS="${MM3_CPUS:-0-15}"
cd "$(dirname "$0")"
export DOTNET_gcServer=1
export DOTNET_gcConcurrent=0
export DOTNET_GCRetainVM=1
exec taskset -c "$CPUS" dotnet run -c Release "$@"
