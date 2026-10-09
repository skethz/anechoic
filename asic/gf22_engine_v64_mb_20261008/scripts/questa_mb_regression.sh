#!/bin/bash
# Questa RTL regression of the multi-bit ASIC core (rtl/sca_core_mb_asic.v, behavioural SRAM models, SCA_SIM_CHECKS)
# against the reference vectors (sim/vec/*: expected words from sca_ref.hpp run_trial(..., dense) or sca_ref_bias.hpp
# run_trial_bias). Usage: questa_mb_regression.sh <set> [<set> ...]   (set = vec directory name, e.g. k2_k2000_short;
# the K of the build is the set's prefix). Runs the sets in parallel (one thread each).
R=/scratch/USER/snowball_asic_v64_mb_20261008
export TMPDIR=$R/tmp
cd $R/sim && mkdir -p reg && cd reg
if [ ! -d work ]; then
  questa-2025.3 vlib work > vlog.log 2>&1
  questa-2025.3 vlog -work work -sv +define+SCA_SIM_CHECKS ../tb_core_mbb.v ../../rtl/sca_core_mb_asic.v ../../rtl/sca_mem_behav.v >> vlog.log 2>&1
  grep -E "Errors:" vlog.log
fi
for v in "$@"; do
  k=${v:1:1}
  ( rm -rf $v; mkdir -p $v; cd $v; cp ../../vec/$v/{stim.hex,expect.hex,counts.txt} .
    questa-2025.3 vsim -c -lib ../work -gK=$k tb -do "run -all; quit -f" > sim.log 2>&1
    grep -E "^# (trial (start|out|end)|launch start|DONE)" sim.log > lines.txt
    echo "$v: K=$k $(grep -E '^# DONE' sim.log | sed 's/^# //') | $(grep -c 'TESTBENCH PASSED' sim.log) passed | checks $(grep -c 'CHECK:' sim.log)" > summary.txt ) &
done
wait
for v in "$@"; do cat $v/summary.txt; done
