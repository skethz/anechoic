#!/usr/bin/env bash
# PROTOCOL_HW_MB fallback ("Fallback, fixed now"): rebuild the same multi-bit design at a lower kernel clock (default 225 MHz)
# with identical settings, then run the same unattended bring-up (protocol hashes re-checked there).
# Usage: KERNEL_MHZ=225 bash scripts/rebuild_mb2_fallback.sh
set -uo pipefail
cd /scratch/USER/sca_v80_20261003
source scripts/environment.sh > /dev/null
export SCA_VARIANT=mb2 SCA_CFLAGS="-DSCA_LANES=256 -DMB_K=2" SCA_HOST_SRC=host_sca_multi_mb.cpp NUM_ENGINES=12 KERNEL_MHZ=${KERNEL_MHZ:-225} \
  INTEGRATE_TCL=$TASK_ROOT/scripts/integrate_kernel_mb.tcl PBLOCK_XDC=$TASK_ROOT/scripts/pblocks_v6_e12.xdc REGSLICE=1 POST_ROUTE_OPT=1
T=board_mb2_e12_${KERNEL_MHZ}mhz
[ -e "build/$T" ] && { echo "fresh build directory required: build/$T"; exit 2; }
env | grep -E "^(SCA_|NUM_ENGINES|KERNEL_MHZ|INTEGRATE_TCL|PBLOCK_XDC|REGSLICE|POST_ROUTE_OPT)" | sort > logs/pipeline_${T}.env
bash scripts/run_pipeline_job.sh > logs/pipeline_${T}.console 2>&1
BOARD_TAG=$T GSET_CFG=$TASK_ROOT/PROTOCOL_HW_MB_gset.txt POWER_CFG=$TASK_ROOT/PROTOCOL_HW_MB_power.txt RESTORE_V64=1 \
  bash scripts/bringup_mb2.sh > logs/bringup_mb2_${T}.console 2>&1
