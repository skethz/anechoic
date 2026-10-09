#!/usr/bin/env bash
set -euo pipefail
source /scratch/USER/sca_v80_20261003/scripts/environment.sh
source "$XILINX_VIVADO/settings64.sh"
export KERNEL_MHZ=${KERNEL_MHZ:-300}
export BOARD_TAG=${BOARD_TAG:-board_${KERNEL_MHZ}mhz}
[[ "$BOARD_TAG" =~ ^[a-zA-Z0-9_]+$ ]]
python3 scripts/prepare_ip.py
exec "$SNOWBALL_ROOT/build/scratch_exec" /scratch/USER "$XILINX_VIVADO/bin/vivado" \
  -mode batch -source "$TASK_ROOT/scripts/build_board.tcl" \
  -log "$TASK_ROOT/logs/${BOARD_TAG}_vivado.log" \
  -journal "$TASK_ROOT/logs/${BOARD_TAG}.jou"
