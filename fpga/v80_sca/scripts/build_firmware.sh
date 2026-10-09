#!/usr/bin/env bash
set -euo pipefail
source /scratch/USER/sca_v80_20261003/scripts/environment.sh
source /opt/xilinx/2025.1/Vitis/settings64.sh
export PATH="/opt/xilinx/2025.1/gnu/armr5/lin/gcc-arm-none-eabi/bin:$PATH"
command -v armr5-none-eabi-gcc >/dev/null
export KERNEL_MHZ=${KERNEL_MHZ:-300}
board_tag=${BOARD_TAG:-board_${KERNEL_MHZ}mhz}
[[ "$board_tag" =~ ^[a-zA-Z0-9_]+$ ]]
amc_build_root=${AMC_BUILD_ROOT:-$TASK_ROOT/build/AMC}
xsa="$TASK_ROOT/build/$board_tag/amd_v80_gen5x8_25.1.xsa"
test -s "$xsa"
python3 - "$xsa" <<'PY'
import sys, zipfile
with zipfile.ZipFile(sys.argv[1]) as archive:
    assert archive.testzip() is None, 'Damaged XSA archive'
    assert any(name.endswith('.hwh') for name in archive.namelist()), 'XSA has no hardware handoff'
PY
if [ ! -d "$amc_build_root" ]; then
  cp -a "$AVED_ROOT/fw/AMC" "$amc_build_root"
fi
cd "$amc_build_root"
"$SNOWBALL_ROOT/build/scratch_exec" /scratch/USER bash scripts/build.sh \
  -os freertos10_xilinx -profile v80 -xsa "$xsa"
# The upstream script can return zero after a piped build failure; check output.
test -s build/amc.elf
"$SNOWBALL_ROOT/build/scratch_exec" /scratch/USER python3 scripts/gen_fpt.py -f scripts/fpt.json
test -s fpt.bin
cp build/amc.elf "$TASK_ROOT/build/$board_tag/amc.elf"
cp fpt.bin "$TASK_ROOT/build/$board_tag/fpt.bin"
