#!/usr/bin/env bash
# PROTOCOL_HW_V6.md: phase 1 (exact gate on all engines), 2a (concurrency), 2b (E-engine rounds). Requires the image programmed.
# Usage: SCA_VARIANT=v6 NUM_ENGINES=9 BOARD_TAG=board_v6_e9_300mhz bash scripts/run_hw_v6.sh
# (Copy of run_hw_multi.sh: trial counts are derived from E, and results go to results/hwv6_<tag>.)
set -uo pipefail
source /scratch/USER/sca_v80_20261003/scripts/environment.sh
[[ "$BOARD_TAG" =~ ^board_([a-z0-9]+_)*[0-9]+mhz$ ]] || exit 2
host="$TASK_ROOT/build/host_sca${SCA_VARIANT:+_$SCA_VARIANT}"; test -x "$host" || exit 2
E=${NUM_ENGINES:-4}
manifest="$TASK_ROOT/build/$BOARD_TAG/image_manifest.json"
uuid=$(python3 -c "import json,sys; print(json.load(open(sys.argv[1]))['logic_uuid'])" "$manifest")
mhz=$(python3 -c "import json,sys; print(json.load(open(sys.argv[1]))['clock']['actual_mhz'])" "$manifest")
out="$TASK_ROOT/results/hwv6_${BOARD_TAG}"
if [ -e "$out" ]; then echo "fresh result directory required: $out"; exit 2; fi
mkdir -p "$out/verify" "$out/concurrency" "$out/rounds"
cp "$manifest" "$out/image_manifest.json"
log="$out/run.log"; echo running > "$out/status"
note() { printf '%s %s\n' "$(date -u +%FT%TZ)" "$*" >> "$log"; }
# KEEP_TABLES=1 (v6.2 images): tables stay resident in the cores after each invocation's first round
KT=${KEEP_TABLES:+--keep-tables}
run() { "$SNOWBALL_ROOT/build/scratch_exec" /scratch/USER --v80-ami "$host" \
          --uuid "$uuid" --clock-mhz "$mhz" --seed 20261004 --t1 5 --engines $E $KT "$@" >> "$log" 2>&1; }
finish() { note "exit $1"; echo "$1" > "$out/status"; exit "$1"; }
note "uuid $uuid clock_mhz $mhz engines $E keep_tables ${KEEP_TABLES:-0}"
# smallest multiples of E * per-launch covering the PROTOCOL_HW_MULTI sizes; rounds: at least 2048 trials (>= 128 rounds)
up() { echo $(( ( ($1 + $2 - 1) / $2 ) * $2 )); }
W1=$(up 8 $E); W2=$(up 8 $((2*E))); W3=$(up 4 $E); W4=$(up 8 $E); NC=$(up 1024 $((64*E))); NR=$(up 2048 $E)
note "W1 $W1 W2 $W2 W3 $W3 W4 $W4 NC $NC NR $NR"
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
  # Amendment 1 (PROTOCOL_HW_V6): schedules re-selected for v6 x 9 engines (research/v6_schedule_20261006), run after O4
  [X1]="--q 8 --tecT 1.75 --ramp --t0 15 --steps 320"    [X2]="--q 6 --lambda 0.9 --ramp --t0 15 --steps 320"
  [X3]="--q 8 --lambda 1.05 --ramp --t0 12 --steps 1260" [X4]="--q 8 --tecT 2.0 --ramp --t0 15 --steps 760"
  [B1]="--q 8 --lambda 0 --t0 40 --steps 960"            [B2]="--q 6 --lambda 0 --t0 40 --steps 1560"
  [B3]="--q 8 --tec-jv -4 --t0 30 --steps 960"           [B4]="--q 8 --tec-jv -4 --t0 30 --steps 1560"
  # Amendment 2: selected for 12 engines (primary objective); exploratory on the 9-engine image
  [X5]="--q 8 --tecT 1.75 --ramp --t0 15 --steps 280"
)
hard=0
for key in P1 P2 P3 P4 T1 T2 T3 O1 O2 O3 O4 X1 X2 X3 X4 B1 B2 B3 B4 X5; do
  # shellcheck disable=SC2086
  run --cohort --output "$out/rounds/$key" ${CFG[$key]} --trials $NR --per-launch 1 --offset 0; r=$?
  note "rounds $key exit $r"; [ "$r" = 1 ] && hard=1
done
finish $((hard ? 21 : 0))
