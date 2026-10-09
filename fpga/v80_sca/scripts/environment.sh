#!/usr/bin/env bash
# Same validated toolchain as fpga/v80_snowball; all writes under this task root.
export TASK_ROOT=/scratch/USER/sca_v80_20261003
export BASELINE_ROOT=/scratch/USER/snowball_u250_v80_20260928
export SNOWBALL_ROOT=/scratch/USER/snowball_v80_20260928
export TMPDIR="$TASK_ROOT/tmp"
export TMP="$TMPDIR" TEMP="$TMPDIR"
export XDG_CACHE_HOME="$TASK_ROOT/cache" XDG_CONFIG_HOME="$TASK_ROOT/config" XDG_STATE_HOME="$TASK_ROOT/state"
export PYTHONPYCACHEPREFIX="$TASK_ROOT/cache/python"
export XILINX_LOCAL_USER_DATA=no
export XILINXD_LICENSE_FILE="/scratch/USER/licenses/smbus_tempo.lic${XILINXD_LICENSE_FILE:+:$XILINXD_LICENSE_FILE}"
export XILINX_VIVADO=/opt/xilinx/2025.1/Vivado
export XILINX_HLS=/opt/tools/Xilinx/Vitis_HLS/2024.1
export PATH="$XILINX_VIVADO/bin:$XILINX_HLS/bin:$PATH"
mkdir -p "$TMPDIR" "$XDG_CACHE_HOME" "$XDG_CONFIG_HOME" "$XDG_STATE_HOME" "$TASK_ROOT/logs" "$TASK_ROOT/results" "$TASK_ROOT/build"
cd "$TASK_ROOT" || return 1
export AVED_ROOT="$BASELINE_ROOT/third_party/AVED"
