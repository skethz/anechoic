#!/usr/bin/env bash
# Revision-3/4 multi-bit images (src/v6/sca_core_mb_r3.v, or r4 for k4b: runtime active spins n, optional per-spin bias; front end
# src/v6/sca_mover_mb3.cpp, host src/host_sca_multi_mb3.cpp). Same board pipeline and settings as board_mb2_e12_250mhz
# (REGSLICE, post-route phys_opt), new variant names so that every image has its own HLS/IP/AMC/build directories.
# Never programs.
# Usage: bash scripts/build_mb3.sh k2   -> variant mb2n: K = 2, no bias, 12 engines (pblocks_v6_e12.xdc)
#        bash scripts/build_mb3.sh k4b  -> variant mb4b: K = 4, bias, HBX = 17, 6 engines (pblocks_v6_e6.xdc)
# [KERNEL_MHZ=250] [WAIT_FOR=<status file>: start only when it no longer says running]
set -uo pipefail
cd /scratch/USER/sca_v80_20261003 && source scripts/environment.sh > /dev/null
[ -n "${WAIT_FOR:-}" ] && until [ -f "$WAIT_FOR" ] && [ "$(cat "$WAIT_FOR")" != running ]; do sleep 60; done
case "$1" in
  k2)  export SCA_VARIANT=mb2n MB_K=2 MB_BIAS=0 MB_HBX=0 NUM_ENGINES=${NUM_ENGINES:-12} PBLOCK_XDC=$TASK_ROOT/scripts/pblocks_v6_e${NUM_ENGINES:-12}.xdc ;;
  k4b) export SCA_VARIANT=mb4b MB_K=4 MB_BIAS=1 MB_HBX=17 NUM_ENGINES=${NUM_ENGINES:-6} PBLOCK_XDC=$TASK_ROOT/scripts/pblocks_v6_e${NUM_ENGINES:-6}.xdc
       CORE_REV=4 ;;   # 8 Oct 03:10 UTC: K = 4 + bias uses revision 4 (S <= 8192 for the exported TSP schedules)
  *) echo "usage: build_mb3.sh k2|k4b"; exit 2 ;;
esac
export SCA_KERNEL_SRC=v6/sca_mover_mb3.cpp SCA_TB_SRC=v6/tb_mover_mb3.cpp SCA_CFLAGS="-DSCA_LANES=256 -DMB_K=$MB_K -DMB_BIAS=$MB_BIAS -DMB_HBX=$MB_HBX"
[ "$(cat results/cosim_${SCA_VARIANT}_exit_status 2>/dev/null)" = 0 ] || { bash scripts/run_hls_job.sh && bash scripts/run_cosim_job.sh || exit 3; }
if [ "${CORE_REV:-3}" = 4 ]; then HOSTSRC=host_sca_multi_mb4.cpp; ITCL=integrate_kernel_mb4.tcl; CORE=src/v6/sca_core_mb_r4.v
else HOSTSRC=host_sca_multi_mb3.cpp; ITCL=integrate_kernel_mb3.tcl; CORE=src/v6/sca_core_mb_r3.v; fi
export SCA_HOST_SRC=$HOSTSRC KERNEL_MHZ=${KERNEL_MHZ:-250} \
  INTEGRATE_TCL=$TASK_ROOT/scripts/$ITCL REGSLICE=1 POST_ROUTE_OPT=1
T=board_${SCA_VARIANT}_e${NUM_ENGINES}_${KERNEL_MHZ}mhz
[ -e "build/$T" ] && { echo "fresh build directory required: build/$T"; exit 2; }
env | grep -E "^(SCA_|MB_|NUM_ENGINES|KERNEL_MHZ|INTEGRATE_TCL|PBLOCK_XDC|REGSLICE|POST_ROUTE_OPT)" | sort > logs/pipeline_${T}.env
sha256sum $CORE src/v6/sca_mover_mb3.cpp src/$HOSTSRC scripts/$ITCL "$PBLOCK_XDC" >> logs/pipeline_${T}.env
bash scripts/run_pipeline_job.sh > logs/pipeline_${T}.console 2>&1
