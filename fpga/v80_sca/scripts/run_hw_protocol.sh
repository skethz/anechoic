#!/usr/bin/env bash
# PROTOCOL_HW.md: phase 1 (exact bring-up gate), then phase 2 (cohorts). Requires the image of BOARD_TAG to be programmed.
# Usage: BOARD_TAG=board_300mhz bash scripts/run_hw_protocol.sh
set -uo pipefail
source /scratch/USER/sca_v80_20261003/scripts/environment.sh
[[ "$BOARD_TAG" =~ ^board_([a-z0-9]+_)*[0-9]+mhz$ ]] || exit 2
host="$TASK_ROOT/build/host_sca${SCA_VARIANT:+_$SCA_VARIANT}"; test -x "$host" || exit 2
manifest="$TASK_ROOT/build/$BOARD_TAG/image_manifest.json"
uuid=$(python3 -c "import json,sys; print(json.load(open(sys.argv[1]))['logic_uuid'])" "$manifest")
mhz=$(python3 -c "import json,sys; print(json.load(open(sys.argv[1]))['clock']['actual_mhz'])" "$manifest")
out="$TASK_ROOT/results/hw_${BOARD_TAG}"
if [ -e "$out" ]; then echo "fresh result directory required: $out"; exit 2; fi
mkdir -p "$out/verify" "$out/cohort" "$out/single"
cp "$manifest" "$out/image_manifest.json"
log="$out/run.log"; echo running > "$out/status"
note() { printf '%s %s\n' "$(date -u +%FT%TZ)" "$*" >> "$log"; }
run() { "$SNOWBALL_ROOT/build/scratch_exec" /scratch/USER --v80-ami "$host" \
          --uuid "$uuid" --clock-mhz "$mhz" --seed 20261004 --t1 5 "$@" >> "$log" 2>&1; }
finish() { note "exit $1"; echo "$1" > "$out/status"; exit "$1"; }
note "uuid $uuid clock_mhz $mhz"

# Phase 1: every trial must match sca_ref exactly; any nonzero exit stops the protocol.
run --verify --output "$out/verify/V1" --q 4 --lambda 0 --t0 30 --steps 40 --trials 8 --per-launch 4 --offset 900000 || finish 11
run --verify --output "$out/verify/V2" --q 4 --lambda 0.7 --ramp --t0 20 --steps 560 --trials 8 --per-launch 8 --offset 900100 || finish 12
run --verify --output "$out/verify/V3" --q 8 --tec-jv -8 --t0 20 --steps 40 --trials 4 --per-launch 1 --offset 900200 || finish 13
run --verify --output "$out/verify/V4" --q 8 --lambda 1.05 --ramp --t0 12 --steps 960 --trials 4 --per-launch 2 --offset 900300 || finish 14
note "phase 1 PASS"

# Phase 2: 1,024 paired trials per configuration (offsets 0..1023), 64 per launch.
declare -A CFG=(
  [P1]="--q 4 --lambda 0 --t0 30 --steps 1560"   [P2]="--q 6 --lambda 0 --t0 30 --steps 1560"
  [P3]="--q 8 --lambda 0 --t0 30 --steps 1560"   [P4]="--q 8 --lambda 0 --t0 30 --steps 960"
  [T1]="--q 8 --tec-jv -4 --t0 30 --steps 1560"  [T2]="--q 8 --tec-jv -4 --t0 30 --steps 960"
  [T3]="--q 8 --tec-jv 4 --t0 30 --steps 1560"
  [O1]="--q 8 --lambda 1.05 --ramp --t0 12 --steps 960" [O2]="--q 6 --lambda 0.9 --ramp --t0 15 --steps 560"
  [O3]="--q 6 --lambda 0.9 --ramp --t0 12 --steps 560"  [O4]="--q 6 --lambda 0.9 --t0 12 --steps 360"
)
hard=0
for key in P1 P2 P3 P4 T1 T2 T3 O1 O2 O3 O4; do
  # shellcheck disable=SC2086
  run --cohort --output "$out/cohort/$key" ${CFG[$key]} --trials 1024 --per-launch 64 --offset 0; r=$?
  note "cohort $key exit $r"; [ "$r" = 1 ] && hard=1   # exit 2 = some trials failed integrity checks (recorded, counted as failures)
done
for key in P1 T1 O1; do
  # shellcheck disable=SC2086
  run --cohort --output "$out/single/$key" ${CFG[$key]} --trials 64 --per-launch 1 --offset 100000; r=$?
  note "single $key exit $r"; [ "$r" = 1 ] && hard=1
done
finish $((hard ? 21 : 0))
