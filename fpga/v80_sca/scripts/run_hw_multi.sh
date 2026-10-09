#!/usr/bin/env bash
# PROTOCOL_HW_MULTI.md: phase 1 (exact gate on all engines), 2a (concurrency), 2b (4-engine rounds). Requires the image programmed.
# Usage: SCA_VARIANT=v52 BOARD_TAG=board_v52_e4_300mhz bash scripts/run_hw_multi.sh
set -uo pipefail
source /scratch/USER/sca_v80_20261003/scripts/environment.sh
[[ "$BOARD_TAG" =~ ^board_([a-z0-9]+_)*[0-9]+mhz$ ]] || exit 2
host="$TASK_ROOT/build/host_sca${SCA_VARIANT:+_$SCA_VARIANT}"; test -x "$host" || exit 2
E=${NUM_ENGINES:-4}
manifest="$TASK_ROOT/build/$BOARD_TAG/image_manifest.json"
uuid=$(python3 -c "import json,sys; print(json.load(open(sys.argv[1]))['logic_uuid'])" "$manifest")
mhz=$(python3 -c "import json,sys; print(json.load(open(sys.argv[1]))['clock']['actual_mhz'])" "$manifest")
out="$TASK_ROOT/results/hwm_${BOARD_TAG}"
if [ -e "$out" ]; then echo "fresh result directory required: $out"; exit 2; fi
mkdir -p "$out/verify" "$out/concurrency" "$out/rounds"
cp "$manifest" "$out/image_manifest.json"
log="$out/run.log"; echo running > "$out/status"
note() { printf '%s %s\n' "$(date -u +%FT%TZ)" "$*" >> "$log"; }
run() { "$SNOWBALL_ROOT/build/scratch_exec" /scratch/USER --v80-ami "$host" \
          --uuid "$uuid" --clock-mhz "$mhz" --seed 20261004 --t1 5 --engines $E "$@" >> "$log" 2>&1; }
finish() { note "exit $1"; echo "$1" > "$out/status"; exit "$1"; }
note "uuid $uuid clock_mhz $mhz engines $E"
if [ "$E" = 3 ]; then W1=9; W2=12; W3=3; W4=9; NC=1152; NR=1026; else W1=8; W2=8; W3=4; W4=8; NC=1024; NR=1024; fi
run --verify --output "$out/verify/W1" --q 4 --lambda 0 --t0 30 --steps 40 --trials $W1 --per-launch 1 --offset 900000 || finish 11
run --verify --output "$out/verify/W2" --q 4 --lambda 0.7 --ramp --t0 20 --steps 560 --trials $W2 --per-launch 2 --offset 900100 || finish 12
run --verify --output "$out/verify/W3" --q 8 --tec-jv -8 --t0 20 --steps 40 --trials $W3 --per-launch 1 --offset 900200 || finish 13
run --verify --output "$out/verify/W4" --q 8 --lambda 1.05 --ramp --t0 12 --steps 960 --trials $W4 --per-launch 1 --offset 900300 || finish 14
note "phase 1 PASS"
run --cohort --output "$out/concurrency/O1" --q 8 --lambda 1.05 --ramp --t0 12 --steps 960 --trials $NC --per-launch 64 --offset 200000
note "phase 2a exit $?"
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
  run --cohort --output "$out/rounds/$key" ${CFG[$key]} --trials $NR --per-launch 1 --offset 0; r=$?
  note "rounds $key exit $r"; [ "$r" = 1 ] && hard=1
done
finish $((hard ? 21 : 0))
