#!/usr/bin/env bash
# PROTOCOL_HW_MB Amendment 1: the software G-set study's 12-engine selections on the multi-bit image (after the main protocol).
# Phase A1-1: exact gate on every N = 2000 instance (12 trials each, verified against sca_ref.hpp run_trial(..., dense)).
# Phase A1-2: one cohort of 2052 trials per frozen configuration (PROTOCOL_HW_MB_A1_gset.txt). Fresh result directory.
# Usage: SCA_VARIANT=mb2 NUM_ENGINES=12 BOARD_TAG=board_mb2_e12_250mhz bash scripts/run_hw_mb_a1.sh
set -uo pipefail
source /scratch/USER/sca_v80_20261003/scripts/environment.sh
[[ "$BOARD_TAG" =~ ^board_([a-z0-9]+_)*[0-9]+mhz$ ]] || exit 2
host="$TASK_ROOT/build/host_sca${SCA_VARIANT:+_$SCA_VARIANT}"; test -x "$host" || exit 2
E=${NUM_ENGINES:-12}
manifest="$TASK_ROOT/build/$BOARD_TAG/image_manifest.json"
uuid=$(python3 -c "import json,sys; print(json.load(open(sys.argv[1]))['logic_uuid'])" "$manifest")
mhz=$(python3 -c "import json,sys; print(json.load(open(sys.argv[1]))['clock']['actual_mhz'])" "$manifest")
out="$TASK_ROOT/results/hwmb_${BOARD_TAG}_a1"
if [ -e "$out" ]; then echo "fresh result directory required: $out"; exit 2; fi
mkdir -p "$out/verify" "$out/gset"
cp "$manifest" "$out/image_manifest.json"
cp "$TASK_ROOT/PROTOCOL_HW_MB_A1_gset.txt" "$TASK_ROOT/PROTOCOL_HW_MB_A1_verify.txt" "$out/"
sha256sum "$host" "$TASK_ROOT"/data/gset/G* > "$out/inputs.sha256"
log="$out/run.log"; echo running > "$out/status"
note() { printf '%s %s\n' "$(date -u +%FT%TZ)" "$*" >> "$log"; }
run() { "$SNOWBALL_ROOT/build/scratch_exec" /scratch/USER --v80-ami "$host" \
          --uuid "$uuid" --clock-mhz "$mhz" --seed 20261004 --engines $E --keep-tables "$@" >> "$log" 2>&1 < /dev/null; }
finish() { note "exit $1"; echo "$1" > "$out/status"; exit "$1"; }
note "uuid $uuid clock_mhz $mhz engines $E"
gate=0
while read -r key inst target trials off args; do
  case "$key" in ''|\#*) continue;; esac
  # shellcheck disable=SC2086
  run --verify --matrix gset:$TASK_ROOT/data/gset/$inst --target $target --output "$out/verify/$key" --trials $trials --per-launch 1 --offset $off $args
  r=$?; note "verify $key exit $r"; [ "$r" = 0 ] || gate=1
done < "$TASK_ROOT/PROTOCOL_HW_MB_A1_verify.txt"
[ "$gate" = 0 ] || finish 11
note "phase A1-1 PASS"
hard=0
while read -r key inst target trials args; do
  case "$key" in ''|\#*) continue;; esac
  # shellcheck disable=SC2086
  run --cohort --matrix gset:$TASK_ROOT/data/gset/$inst --target $target --output "$out/gset/$key" --trials $trials --per-launch 1 --offset 0 $args
  r=$?; note "gset $key exit $r"; [ "$r" = 1 ] && hard=1
done < "$TASK_ROOT/PROTOCOL_HW_MB_A1_gset.txt"
finish $((hard ? 21 : 0))
