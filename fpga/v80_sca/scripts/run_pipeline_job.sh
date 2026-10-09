#!/usr/bin/env bash
# Persistent end-to-end build for one kernel clock: board (prepare_ip + Vivado) -> AMC firmware -> host -> package + manifest.
# Each clock gets its own IP repository, AMC copy and status files, so several clocks can build in parallel. Never programs.
set -uo pipefail
source /scratch/USER/sca_v80_20261003/scripts/environment.sh
export KERNEL_MHZ=${KERNEL_MHZ:-250}
export BOARD_TAG=board${SCA_VARIANT:+_$SCA_VARIANT}${NUM_ENGINES:+_e$NUM_ENGINES}_${KERNEL_MHZ}mhz
export SNOWBALL_IP_REPO="$TASK_ROOT/build/iprepo_${BOARD_TAG}"
export AMC_BUILD_ROOT="$TASK_ROOT/build/AMC_${BOARD_TAG}"
export FIRMWARE_STATUS_FILE="results/firmware_${BOARD_TAG}_exit_status"
status="results/pipeline_${BOARD_TAG}_exit_status"
printf 'running\n' > "$status"
step() { printf '%s %s\n' "$(date -u +%FT%TZ)" "$1" >> "logs/pipeline_${BOARD_TAG}.log"; }
step start
bash scripts/run_board_job.sh && step board_ok && \
{ printf 'running\n' > "$FIRMWARE_STATUS_FILE"; bash scripts/build_firmware.sh > "logs/firmware_${BOARD_TAG}.log" 2>&1; r=$?; \
  printf '%s\n' "$r" > "$FIRMWARE_STATUS_FILE"; test "$r" = 0; } && step firmware_ok && \
bash scripts/build_host.sh > "logs/host_build_${BOARD_TAG}.log" 2>&1 && step host_ok && \
bash scripts/package_image.sh && step package_ok
result=$?
step "exit $result"
printf '%s\n' "$result" > "$status"
exit "$result"
