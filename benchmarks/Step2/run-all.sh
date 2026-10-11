#!/usr/bin/env bash
# Production measurements for the step2 paper (sequential; pinned to the P-cores).
set -uo pipefail
cd "$(dirname "$0")/../src"
# refuse to measure on a loaded machine: timings under foreign load are not reportable
load=$(cut -d' ' -f1 /proc/loadavg)
busy=$(ps -eo pcpu= | awk '$1 > 50 {n++} END {print n+0}')
if awk -v l="$load" 'BEGIN {exit !(l > 2.0)}' || [ "$busy" -gt 0 ]; then
  echo "machine busy (load $load, $busy processes above 50% CPU); not measuring" >&2
  ps -eo pid,pcpu,comm --sort=-pcpu | head -8 >&2
  exit 3
fi
L=../benchmarks/logs
R=./run-maxperf.sh
date > $L/started.txt
$R mm3ref                          > $L/mm3ref.log 2>&1
for n in 12 26 55 119; do
  every=$(( (n + 1) / 8 )); [ $every -lt 1 ] && every=1
  $R work crack  $n 0  $every 1    > $L/work_crack_${n}_p.log 2>&1
  $R work crack  $n 0  1000000 0   > $L/work_crack_${n}_t.log 2>&1
  $R work bisect $n 40 5 1         > $L/work_bisect_${n}_p.log 2>&1
  $R work erode  $n 40 5 1         > $L/work_erode_${n}_p.log 2>&1
done
$R interleaved 6,10,14,18,22,26    > $L/interleaved.log 2>&1
$R crossover 55                    > $L/crossover_55.log 2>&1
$R reads 55                        > $L/reads_55.log 2>&1
$R reads 119                       > $L/reads_119.log 2>&1
$R concurrent 55 3 1,2,4,8,15      > $L/concurrent_55.log 2>&1
$R persist 55                      > $L/persist_55.log 2>&1
$R persist 119                     > $L/persist_119.log 2>&1
$R memory 55                       > $L/memory_55.log 2>&1
$R memory 119                      > $L/memory_119.log 2>&1
$R endtoend 55 100                 > $L/endtoend_55.log 2>&1
STEP2_PROFILE=1 $R profile crack 55 0 > $L/profile_crack_55_t.log 2>&1
STEP2_PROFILE=1 $R profile crack 55 1 > $L/profile_crack_55_p.log 2>&1
STEP2_PROFILE=1 $R profile bisect 55 1 > $L/profile_bisect_55_p.log 2>&1
(uname -a; lscpu; free -g; cat /sys/devices/system/cpu/cpu0/cpufreq/scaling_driver \
  /sys/devices/system/cpu/cpu0/cpufreq/scaling_governor \
  /sys/devices/system/cpu/cpu0/cpufreq/energy_performance_preference; \
  $HOME/.dotnet/dotnet --info | head -30) > $L/machine.txt 2>&1
date > $L/finished.txt
