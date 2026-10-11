#!/usr/bin/env bash
# 2026-10-09: final run, everything with the final code (initial states with the greedy colouring; (flat pattern copy, threaded
# rebuild with chunked transpose and deterministic colouring, gathered assembly, row-limit check).
# Gates first; every measurement after them, one at a time, on the otherwise idle machine.
set -uo pipefail
cd "$(dirname "$0")/../src"
L=../benchmarks/logs_fin; mkdir -p $L
D=$(cd ../benchmarks && pwd)/data_fin; mkdir -p $D
export STEP2_DATA=$D
# idle gate: amd2 is shared with the SimPlas job queue; a measurement starts only when no foreign heavy
# process runs and the load stays below 1 for three minutes, and is repeated if one appears meanwhile
FOREIGN="SimPlas|xfortcom|ifx|gfortran|VBCSCompiler|simplas-amd/jobs"
foreign() { pgrep -f "$FOREIGN" > /dev/null; }
wait_idle() {
  local ok=0
  while [ $ok -lt 3 ]; do
    if foreign || awk -v l="$(cut -d' ' -f1 /proc/loadavg)" 'BEGIN{exit !(l > 1.0)}'; then ok=0; else ok=$((ok+1)); fi
    sleep 60
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
    echo "$(date) foreign load during: $*; repeating" >> ../benchmarks/logs_fin/interference.log
  done
  return $rc
}
date > $L/timed_resumed.txt
# resumed 2026-10-09 08:55 after the idle gate starved 5 h behind the SimPlas queue; erode 119 repeats (interfered 03:38)
R work erode 119 40 5 1         > $L/work_erode_119_p.log 2>&1
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
date > $L/section7_finished.txt
# ---- cutting
for n in 55 119; do for g in 4 8; do R batchcut $n $g 3 > $L/batchcut_${n}_${g}.log 2>&1; done; done
for n in 55 119; do for g in 4 8; do R cutpar $n $g 1,2,4,8,16 5 > $L/cutpar_${n}_${g}.log 2>&1; done; done
for g in 2 3; do R interfaces 32 $g 16 $D/interfaces_32_$g > $L/interfaces_32_$g.log 2>&1; done
# ---- growth (n = 20, h0/8)
mkdir -p $D/growth; STEP2_FRONTDBG=1 STEP2_DATA=$D/growth R growth 20 8 0.25 0.15 14 0.025 4 0.04 0.06 > $L/growth.log 2>&1

date > $L/timed_finished.txt
