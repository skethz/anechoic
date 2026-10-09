#!/bin/bash
# Questa RTL regression: unmodified v6.4b vs ASIC variant on all vector sets; compares the per-trial cycle lines.
export TMPDIR=/scratch/USER/snowball_asic_v64_20261007/tmp
cd /scratch/USER/snowball_asic_v64_20261007/gls
mkdir -p reg && cd reg
questa-2025.3 vlib work_v64b > vlog.log 2>&1; questa-2025.3 vlib work_asic >> vlog.log 2>&1
questa-2025.3 vlog -work work_v64b -sv ../../rtl/tb_core.v ../../rtl/sca_core_v64b.v >> vlog.log 2>&1
questa-2025.3 vlog -work work_asic -sv +define+SCA_SIM_CHECKS ../../rtl/tb_core.v ../../rtl/sca_core_asic.v ../../rtl/sca_mem_behav.v >> vlog.log 2>&1
grep -E "Errors:" vlog.log
for v in default O1 O2 O3 O4 P1 T1 X5; do
  for d in v64b asic; do
    ( mkdir -p ${d}_$v; cd ${d}_$v; cp ../../vecs/vec_$v/{stim.hex,expect.hex,counts.txt} .
      questa-2025.3 vsim -c -lib ../work_$d tb -do "run -all; quit -f" > sim.log 2>&1
      grep -E "^# (trial (start|out|end)|launch start|DONE)" sim.log > lines.txt ) &
  done
  wait
  if cmp -s v64b_$v/lines.txt asic_$v/lines.txt; then s=IDENTICAL; else s=DIFFER; fi
  echo "$v: v64b [$(grep DONE v64b_$v/lines.txt)] asic [$(grep DONE asic_$v/lines.txt)] | trial lines: $s | checks: $(grep -c CHECK: asic_$v/sim.log)"
done
