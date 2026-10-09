#!/usr/bin/env bash
# Sequential group-level OOC syntheses for the bias / width cost table (see ooc_group_mb3.tcl).
set -uo pipefail
source /scratch/USER/sca_v80_20261003/scripts/environment.sh > /dev/null
source "$XILINX_VIVADO/settings64.sh" > /dev/null
for v in "k2b0 2 12 0" "k2b1 2 12 1" "k4b0 4 15 0" "k4b0h17 4 17 0" "k4b1h17 4 17 1"; do
  set -- $v
  OOC_TAG=$1 G_K=$2 G_HB=$3 G_BIAS=$4 "$SNOWBALL_ROOT/build/scratch_exec" /scratch/USER "$XILINX_VIVADO/bin/vivado" -mode batch \
    -source "$TASK_ROOT/scripts/ooc_group_mb3.tcl" -log "$TASK_ROOT/logs/ooc_group_mb3_$1.log" -journal "$TASK_ROOT/logs/ooc_group_mb3_$1.jou" \
    > "$TASK_ROOT/logs/ooc_group_mb3_$1.console" 2>&1
  echo "$1 exit $?" >> "$TASK_ROOT/logs/ooc_group_mb3.status"
done
