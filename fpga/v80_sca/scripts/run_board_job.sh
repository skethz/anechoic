#!/usr/bin/env bash
set -uo pipefail
source /scratch/USER/sca_v80_20261003/scripts/environment.sh
export KERNEL_MHZ=${KERNEL_MHZ:-300}
export BOARD_TAG=${BOARD_TAG:-board_${KERNEL_MHZ}mhz}
printf 'running\n' > "results/${BOARD_TAG}_exit_status"
python3 scripts/prepare_ip.py > "logs/prepare_ip_${BOARD_TAG}.log" 2>&1 && \
  bash scripts/run_board.sh > "logs/${BOARD_TAG}_console.log" 2>&1
result=$?
printf '%s\n' "$result" > "results/${BOARD_TAG}_exit_status"
exit "$result"
