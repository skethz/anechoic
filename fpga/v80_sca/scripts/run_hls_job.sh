#!/usr/bin/env bash
# Persistent HLS job: survives SSH drops, records exit status. Writes confined to /scratch/USER by scratch_exec.
set -uo pipefail
source /scratch/USER/sca_v80_20261003/scripts/environment.sh
sfx=${SCA_VARIANT:+_$SCA_VARIANT}
printf 'running\n' > results/hls${sfx}_exit_status
(
  source /opt/tools/Xilinx/Vivado/2024.1/settings64.sh
  source /opt/tools/Xilinx/Vitis_HLS/2024.1/settings64.sh
  export LD_LIBRARY_PATH="$BASELINE_ROOT/lib/compat24${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"
  mkdir -p "$TASK_ROOT/build/hls_run${sfx}" && cd "$TASK_ROOT/build/hls_run${sfx}" && ulimit -s unlimited
  exec "$SNOWBALL_ROOT/build/scratch_exec" /scratch/USER /opt/tools/Xilinx/Vitis_HLS/2024.1/bin/vitis_hls -f "$TASK_ROOT/scripts/hls.tcl"
) > logs/hls${sfx}.log 2>&1
result=$?
printf '%s\n' "$result" > results/hls${sfx}_exit_status
exit "$result"
