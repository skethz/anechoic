#!/usr/bin/env bash
# Revision-2 multi-bit image (cut-shift fix, src/v6/sca_core_mb_r2.v): HLS + co-simulation of the unchanged front end under
# variant mb2r2, then the same board pipeline and settings as board_mb2_e12_250mhz. Never programs.
# Usage: [KERNEL_MHZ=250] bash scripts/build_mb2r2.sh   (WAIT_FOR=<status file>: start only when it no longer says running)
set -uo pipefail
cd /scratch/USER/sca_v80_20261003 && source scripts/environment.sh > /dev/null
[ -n "${WAIT_FOR:-}" ] && until [ -f "$WAIT_FOR" ] && [ "$(cat "$WAIT_FOR")" != running ]; do sleep 60; done
export SCA_VARIANT=mb2r2 SCA_KERNEL_SRC=v6/sca_mover_mb.cpp SCA_TB_SRC=v6/tb_mover_mb.cpp SCA_CFLAGS="-DSCA_LANES=256 -DMB_K=2"
[ "$(cat results/cosim_mb2r2_exit_status 2>/dev/null)" = 0 ] || { bash scripts/run_hls_job.sh && bash scripts/run_cosim_job.sh || exit 3; }
export SCA_HOST_SRC=host_sca_multi_mb.cpp NUM_ENGINES=12 KERNEL_MHZ=${KERNEL_MHZ:-250} \
  INTEGRATE_TCL=$TASK_ROOT/scripts/integrate_kernel_mb_r2.tcl PBLOCK_XDC=$TASK_ROOT/scripts/pblocks_v6_e12.xdc REGSLICE=1 POST_ROUTE_OPT=1
T=board_mb2r2_e12_${KERNEL_MHZ}mhz
[ -e "build/$T" ] && { echo "fresh build directory required: build/$T"; exit 2; }
env | grep -E "^(SCA_|NUM_ENGINES|KERNEL_MHZ|INTEGRATE_TCL|PBLOCK_XDC|REGSLICE|POST_ROUTE_OPT)" | sort > logs/pipeline_${T}.env
bash scripts/run_pipeline_job.sh > logs/pipeline_${T}.console 2>&1
