#!/usr/bin/env bash
# 2026-10-11: THE final timed run, one protocol for every timing of Sections 7-8: final library (SplitNode, corner crack), harness with
# collector accounting (every timed region records its collections and pauses; the cutting tables use the time without the pauses),
# one measurement at a time on the otherwise idle machine, previous logs and data kept in *_before_20261011.
set -uo pipefail
cd "$(dirname "$0")/../src"
L=../benchmarks/logs_fin4; mkdir -p $L
D=$(cd ../benchmarks && pwd)/data_fin4; mkdir -p $D
export STEP2_DATA=$D
FOREIGN="SimPlas|xfortcom|ifx|gfortran|VBCSCompiler|simplas-amd/jobs"
foreign() { pgrep -f "$FOREIGN" > /dev/null; }
wait_idle() {
  local ok=0
  while [ $ok -lt 2 ]; do
    if foreign || awk -v l="$(cut -d' ' -f1 /proc/loadavg)" 'BEGIN{exit !(l > 1.0)}'; then ok=0; else ok=$((ok+1)); fi
    sleep 30
  done
}
R() {
  local try rc mon
  for try in 1 2 3 4 5 6; do
    wait_idle
    rm -f /tmp/step2_foreign.$$
    ( while true; do foreign && touch /tmp/step2_foreign.$$; sleep 5; done ) &
    mon=$!
    ./run-maxperf.sh "$@"; rc=$?
    kill $mon 2>/dev/null; wait $mon 2>/dev/null
    if [ ! -f /tmp/step2_foreign.$$ ]; then return $rc; fi
    echo "$(date) foreign load during: $*; repeating" >> $L/interference.log
  done
  return $rc
}
date > $L/started.txt
~/.dotnet/dotnet build Step2.Bench -c Release > $L/build.log 2>&1 || { echo "build failed" > $L/failed.txt; exit 1; }
~/.dotnet/dotnet build-server shutdown > /dev/null 2>&1
(uname -a; lscpu; free -g; cat /sys/devices/system/cpu/cpu0/cpufreq/scaling_driver /sys/devices/system/cpu/cpu0/cpufreq/scaling_governor \
  /sys/devices/system/cpu/cpu0/cpufreq/energy_performance_preference; $HOME/.dotnet/dotnet --info | head -30) > $L/machine.txt 2>&1
# ---- cutting first (the tables most sensitive to the collector)
for n in 55 119; do for g in 4 8; do R batchcut $n $g 3 > $L/batchcut_${n}_${g}.log 2>&1; done; done
for n in 55 119; do for g in 4 8; do R cutpar $n $g 1,2,4,8,16 5 > $L/cutpar_${n}_${g}.log 2>&1; done; done
for n in 55 119; do R cutgc $n 4 5 > $L/cutgc_${n}_4.log 2>&1; done
for g in 2 3; do R interfaces 32 $g 16 $D/interfaces_32_$g > $L/interfaces_32_$g.log 2>&1; done
date > $L/cutting_finished.txt
# ---- Section 8 (as run-fin-amd2.sh)
R mm3ref                          > $L/mm3ref.log 2>&1
for n in 12 26 55 119; do
  every=$(( (n + 1) / 8 )); [ $every -lt 1 ] && every=1
  R work crack  $n 0  $every 1    > $L/work_crack_${n}_p.log 2>&1
  R work crack  $n 0  1000000 0   > $L/work_crack_${n}_t.log 2>&1
  R work bisect $n 40 5 1         > $L/work_bisect_${n}_p.log 2>&1
  R work erode  $n 40 5 1         > $L/work_erode_${n}_p.log 2>&1
done
R interleaved 6,10,14,18,22,26    > $L/interleaved.log 2>&1
R crossover 55                    > $L/crossover_55.log 2>&1
R reads 55                        > $L/reads_55.log 2>&1
R reads 119                       > $L/reads_119.log 2>&1
R concurrent 55 3 1,2,4,8,15      > $L/concurrent_55.log 2>&1
R persist 55                      > $L/persist_55.log 2>&1
R persist 119                     > $L/persist_119.log 2>&1
R memory 55                       > $L/memory_55.log 2>&1
R memory 119                      > $L/memory_119.log 2>&1
R endtoend 55 100                 > $L/endtoend_55.log 2>&1
STEP2_PROFILE=1 R profile crack 55 0 > $L/profile_crack_55_t.log 2>&1
STEP2_PROFILE=1 R profile crack 55 1 > $L/profile_crack_55_p.log 2>&1
STEP2_PROFILE=1 R profile bisect 55 1 > $L/profile_bisect_55_p.log 2>&1
date > $L/section8_finished.txt
# ---- growth (n = 20, h0/8)
mkdir -p $D/growth; STEP2_FRONTDBG=1 STEP2_DATA=$D/growth R growth 20 8 0.25 0.15 14 0.025 4 0.04 0.06 > $L/growth.log 2>&1
date > $L/finished.txt
