#!/usr/bin/env bash
# Per-engine utilization of the v6.2 and v6.4 routed images, run only after the v6.4 pipeline has exited
# (two concurrent Vivado sessions exceed fpga-host's 62 GB RAM). Sequential, low priority, read-only on the checkpoints.
cd /scratch/USER/sca_v80_20261003 && source scripts/environment.sh
until [ "$(cat results/pipeline_board_v64_e12_250mhz_exit_status 2>/dev/null)" != running ]; do sleep 60; done
mkdir -p results/util_20261007 && cd results/util_20261007
for t in v62:board_v62_e12_275mhz v64:board_v64_e12_250mhz; do
  tag=${t%%:*}; b=${t#*:}; d=../../build/$b/prj.runs/impl_1/top_wrapper_postroute_physopt.dcp
  [ -f "$d" ] || d=../../build/$b/prj.runs/impl_1/top_wrapper_routed.dcp
  nice -n 10 vivado -mode batch -nojournal -nolog -source ../../scripts/hier_util.tcl -tclargs "$d" "$tag" > ${tag}_vivado.log 2>&1
  echo "$tag exit $?" >> util.log
done
