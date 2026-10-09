#!/usr/bin/env bash
# r3 simulation driver (copy of sim_mb.sh): core src/v6/sca_core_mb_r3.v, testbench tb_core_mb3.v, vectors from
# build/gen_vectors_mb3 (golden model sca_ref_bias.hpp run_trial_bias) with a cases file; xelab generics K, BIAS, HBX, DRAIN.
# Usage: [BIAS=0|1] [HBX=<n>] [DRAIN=<n>] bash scripts/sim_mb3.sh <run_dir> <K> <matrix> <cases file> [bias 0|1|auto] [nfield]
# (Original sim_mb.sh notes follow.)
# RTL regression of the multi-bit core (src/v6/sca_core_mb.v) against the golden model, on fpga-host with the Vivado 2025.1
# simulator (xsim): generate vectors (build/gen_vectors_mb, sca_ref.hpp run_trial(..., dense)), compile, run tb_core_mb.v.
# Usage: [DRAIN=<n>] [SIM_CORE=<v> SIM_TB=<v>] bash scripts/sim_mb.sh <run_dir (fresh, under build/)> <K> <matrix spec> [suite [case]]
# (SIM_CORE/SIM_TB: simulate another core on the same vectors, e.g. v6.4b with K = 1 vectors for a cycle comparison.)
# Writes gen.log, cases.txt, xsim.log and status (PASS/FAIL) into the run directory. Never touches the board.
set -uo pipefail
source /scratch/USER/sca_v80_20261003/scripts/environment.sh
source "$XILINX_VIVADO/settings64.sh" > /dev/null
run=$1; K=$2; mat=$3; cases=$4; bm=${5:-auto}; nf=${6:-}; drain=${DRAIN:-2}; bias=${BIAS:-0}; hbx=${HBX:-0}
core=${SIM_CORE:-$TASK_ROOT/src/v6/sca_core_mb_r3.v}; tbf=${SIM_TB:-$TASK_ROOT/src/v6/tb_core_mb3.v}
case "$run" in "$TASK_ROOT"/build/*) ;; *) echo "run dir must be under $TASK_ROOT/build"; exit 2;; esac
[ -e "$run" ] && { echo "fresh run directory required: $run"; exit 2; }
mkdir -p "$run" && cd "$run" || exit 2
X="$SNOWBALL_ROOT/build/scratch_exec /scratch/USER"
echo running > status
{ date -u +%FT%TZ; echo "K=$K matrix=$mat cases=$cases bias_mode=$bm nfield=$nf BIAS=$bias HBX=$hbx DRAIN=$drain"; sha256sum "$core" "$tbf" "$TASK_ROOT/src/sca_ref.hpp" "$TASK_ROOT/src/sca_ref_bias.hpp" "$TASK_ROOT/src/mb_common3.hpp" "$TASK_ROOT/src/v6/gen_vectors_mb3.cpp" "$cases"; } > run.info
cp "$cases" cases.in
$X "$TASK_ROOT/build/gen_vectors_mb3" "$K" "$mat" . cases.in "$bm" $nf > gen.log 2>&1 || { echo FAIL > status; echo "vector generation failed"; exit 1; }
$X xvlog -nolog "$core" "$tbf" > xvlog.log 2>&1 || { echo FAIL > status; exit 1; }
$X xelab -nolog -O3 --mt 4 -generic_top "K=$K" -generic_top "DRAIN=$drain" -generic_top "BIAS=$bias" -generic_top "HBX=$hbx" tb -s tb_mb > xelab.log 2>&1 || { echo FAIL > status; exit 1; }
start=$(date +%s)
$X xsim -nolog tb_mb -R > xsim.log 2>&1
echo "sim_seconds $(( $(date +%s) - start ))" >> run.info
if grep -q "TESTBENCH PASSED" xsim.log; then echo PASS > status; else echo FAIL > status; fi
cat status
