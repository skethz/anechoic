#!/usr/bin/env bash
# After the multi-bit board build: per-engine utilization of the routed image (hier_util.tcl, as util_after_build.sh did for
# v6.2/v6.4), then the K = 8 out-of-context synthesis estimate. Sequential and low priority; read-only on the checkpoints.
# (The bring-up runs concurrently; its board runs use little CPU and memory.)
set -uo pipefail
cd /scratch/USER/sca_v80_20261003 && source scripts/environment.sh > /dev/null
T=${BOARD_TAG:-board_mb2_e12_250mhz}
until [ "$(cat results/pipeline_${T}_exit_status 2>/dev/null)" != running ]; do sleep 60; done
X="$SNOWBALL_ROOT/build/scratch_exec /scratch/USER"
out=results/util_mb_${T}; mkdir -p "$out"
d=$TASK_ROOT/build/$T/prj.runs/impl_1/top_wrapper_postroute_physopt.dcp
[ -f "$d" ] || d=$TASK_ROOT/build/$T/prj.runs/impl_1/top_wrapper_routed.dcp
if [ -f "$d" ]; then
  (cd "$out" && nice -n 10 $X "$XILINX_VIVADO/bin/vivado" -mode batch -nojournal -nolog -source "$TASK_ROOT/scripts/hier_util.tcl" \
     -tclargs "$d" mb2 > mb2_vivado.log 2>&1; echo "mb2 util exit $?" >> util.log)
fi
OOC_TAG=k8_synth OOC_PERIOD=4.0 OOC_PBLOCK=CLOCKREGION_X2Y5:CLOCKREGION_X4Y7 OOC_SYNTH_ONLY=1 MB_K=8 nice -n 10 bash scripts/ooc_mb_job.sh
echo "post-build done" >> "$out/util.log"
