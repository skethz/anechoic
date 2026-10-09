#!/usr/bin/env bash
# Program a manifest-verified SCA image with the scoped sudo command, then check UUID and BARs.
# Usage: BOARD_TAG=board_300mhz bash scripts/program_board.sh
set -euo pipefail
source /scratch/USER/sca_v80_20261003/scripts/environment.sh
[[ "$BOARD_TAG" =~ ^board_([a-z0-9]+_)*[0-9]+mhz$ ]]
manifest="$TASK_ROOT/build/$BOARD_TAG/image_manifest.json"
log="$TASK_ROOT/results/programming_${BOARD_TAG}_$(date -u +%Y%m%dT%H%M%SZ).log"
want=$(python3 -c "import json,sys; print(json.load(open(sys.argv[1]))['logic_uuid'])" "$manifest")
probe() { "$SNOWBALL_ROOT/build/scratch_exec" /scratch/USER --v80-ami "$TASK_ROOT/build/host_sca" --probe; }
bar2() { awk 'NR==3{print $1, $2}' /sys/bus/pci/devices/0000:01:00.0/resource; }
{
  echo "expected logic_uuid $want"
  echo "before: $(probe || echo probe-failed)"
  if pgrep -u "$(id -u)" -f 'host_sca --(verify|cohort)|run_physical_trials' >/dev/null; then echo "board busy"; exit 3; fi
  python3 -I "$TASK_ROOT/scripts/stage_program_image.py" --manifest "$manifest"
  sudo -n /usr/local/bin/ami_tool cfgmem_program -d 01:00.0 -t primary -p 1 \
    -i /scratch/USER/snowball_v80_20260928/programming/snowball_v80.pdi -y
  echo "after program: $(probe || echo probe-failed)"
  echo "BAR2 resource: $(bar2)"
  if [ "$(bar2 | awk '{print $1}')" = "0x0000000000000000" ]; then
    echo "BAR2 missing: running scoped rescan helper"
    sudo -n /usr/local/sbin/v80-rescan
    echo "BAR2 resource after rescan: $(bar2)"
  fi
  got=$(probe | python3 -c "import json,sys; print(json.load(sys.stdin)['logic_uuid'])")
  echo "after: $got"
  test "$got" = "$want" && echo "UUID MATCH" || { echo "UUID MISMATCH"; exit 4; }
  test "$(bar2 | awk '{print $1}')" != "0x0000000000000000" && echo "BAR2 OK"
} 2>&1 | tee "$log"
