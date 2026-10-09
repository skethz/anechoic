#!/usr/bin/env bash
set -uo pipefail
source /scratch/USER/sca_v80_20261003/scripts/environment.sh
export KERNEL_MHZ=${KERNEL_MHZ:-300}
printf 'running\n' > results/firmware_exit_status
bash scripts/build_firmware.sh > logs/firmware.log 2>&1
result=$?
printf '%s\n' "$result" > results/firmware_exit_status
exit "$result"
