#!/usr/bin/env bash
# Resume a routed board build that missed timing: post-route phys_opt + finalize, then firmware, host, package.
# Usage: SCA_VARIANT=v53 NUM_ENGINES=3 BOARD_TAG=board_v53_e3_225mhz KERNEL_MHZ=225 bash scripts/resume_post_route_job.sh
set -uo pipefail
source /scratch/USER/sca_v80_20261003/scripts/environment.sh
source "$XILINX_VIVADO/settings64.sh"
[[ "$BOARD_TAG" =~ ^board_([a-z0-9]+_)*[0-9]+mhz$ ]] || exit 2
export SNOWBALL_IP_REPO="$TASK_ROOT/build/iprepo_${BOARD_TAG}"
export AMC_BUILD_ROOT="$TASK_ROOT/build/AMC_${BOARD_TAG}"
export FIRMWARE_STATUS_FILE="results/firmware_${BOARD_TAG}_exit_status"
status="results/pipeline_${BOARD_TAG}_exit_status"; printf 'running\n' > "$status"
step() { printf '%s %s\n' "$(date -u +%FT%TZ)" "$1" >> "logs/pipeline_${BOARD_TAG}.log"; }
step "resume: finish image (finalize only)"
printf 'running\n' > "results/${BOARD_TAG}_exit_status"
"$SNOWBALL_ROOT/build/scratch_exec" /scratch/USER "$XILINX_VIVADO/bin/vivado" -mode batch \
  -source "$TASK_ROOT/scripts/finish_image.tcl" -log "$TASK_ROOT/logs/${BOARD_TAG}_finish_vivado.log" \
  -journal "$TASK_ROOT/logs/${BOARD_TAG}_finish.jou" > "$TASK_ROOT/logs/${BOARD_TAG}_finish.log" 2>&1
r=$?; printf '%s\n' "$r" > "results/${BOARD_TAG}_exit_status"
[ "$r" = 0 ] && step board_ok && \
{ printf 'running\n' > "$FIRMWARE_STATUS_FILE"; bash scripts/build_firmware.sh > "logs/firmware_${BOARD_TAG}.log" 2>&1; f=$?; \
  printf '%s\n' "$f" > "$FIRMWARE_STATUS_FILE"; test "$f" = 0; } && step firmware_ok && \
bash scripts/build_host.sh > "logs/host_build_${BOARD_TAG}.log" 2>&1 && step host_ok && \
bash scripts/package_image.sh && step package_ok
result=$?; step "exit $result"; printf '%s\n' "$result" > "$status"; exit "$result"
