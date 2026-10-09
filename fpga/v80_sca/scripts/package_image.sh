#!/usr/bin/env bash
set -euo pipefail
source /scratch/USER/sca_v80_20261003/scripts/environment.sh
source /opt/xilinx/2025.1/Vitis/settings64.sh
export KERNEL_MHZ=${KERNEL_MHZ:-300}
export BOARD_TAG=${BOARD_TAG:-board_${KERNEL_MHZ}mhz}
[[ "$BOARD_TAG" =~ ^[a-zA-Z0-9_]+$ ]]
test "$(cat "results/${BOARD_TAG}_exit_status")" = 0
firmware_status=${FIRMWARE_STATUS_FILE:-results/firmware_exit_status}
test "$(cat "$firmware_status")" = 0
board="$TASK_ROOT/build/$BOARD_TAG"
test -s "$board/prj.runs/impl_1/top_wrapper.pdi"
test -s "$board/amc.elf"
test -s "$board/kernel_clock.json"
python3 "$TASK_ROOT/scripts/check_timing_report.py" "$board/routed_timing.rpt"
cd "$board"
cat > snowball_combine.bif <<'BIF'
all: {
  image {
    { type=bootimage, file=./prj.runs/impl_1/top_wrapper.pdi }
  }
  image {
    id=0x1c000000, name=rpu_subsystem, delay_handoff
    { core=r5-0, file=./amc.elf }
  }
}
BIF
"$SNOWBALL_ROOT/build/scratch_exec" /scratch/USER bootgen -arch versal \
  -image snowball_combine.bif -w -o snowball_v80.pdi \
  > "$TASK_ROOT/logs/package_${BOARD_TAG}.log" 2>&1
"$SNOWBALL_ROOT/build/scratch_exec" /scratch/USER python3 "$TASK_ROOT/scripts/image_manifest.py"
