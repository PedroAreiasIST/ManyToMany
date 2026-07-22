#!/usr/bin/env bash
# Ablation: separate the contribution of (a) the runtime/compiler configuration
# from (b) pinning to the performance cores. Four conditions, same binary.
set -uo pipefail
OUT="${MM3_ABLATION_OUT:-$(pwd)/ablation-results}"
mkdir -p "$OUT"
cd "$(dirname "$0")"

flags_on() {
  export DOTNET_TieredCompilation=0 DOTNET_TieredPGO=0 DOTNET_ReadyToRun=0 \
         DOTNET_TC_QuickJitForLoops=0 DOTNET_gcServer=1 DOTNET_gcConcurrent=0 \
         DOTNET_GCgen0size=20000000 DOTNET_GCRetainVM=1
}
flags_off() {
  unset DOTNET_TieredCompilation DOTNET_TieredPGO DOTNET_ReadyToRun \
        DOTNET_TC_QuickJitForLoops DOTNET_gcServer DOTNET_gcConcurrent \
        DOTNET_GCgen0size DOTNET_GCRetainVM
  # csproj still requests server GC; force workstation+concurrent to isolate the flags
  export DOTNET_gcServer=0 DOTNET_gcConcurrent=1
}

run() { # $1=label $2=pin(yes/no)
  echo "===== $1 ====="
  if [ "$2" = yes ]; then
    taskset -c 0-15 setarch -R dotnet run -c Release --no-build 2>/dev/null
  else
    dotnet run -c Release --no-build 2>/dev/null
  fi
}

flags_off; run "A: baseline (no flags, no pinning)"      no  > "$OUT/A-baseline.log"
flags_off; run "B: pinning only"                          yes > "$OUT/B-pin-only.log"
flags_on;  run "C: flags only"                            no  > "$OUT/C-flags-only.log"
flags_on;  run "D: flags + pinning (MAX-PERF)"            yes > "$OUT/D-flags-and-pin.log"
echo "ABLATION_DONE -> $OUT"
