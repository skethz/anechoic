#!/usr/bin/env bash
# PROTOCOL_HW_MB Amendment 4 (revision-4 K = 4 + bias image, 6 engines; functional validation of the bias path): phase 1 exact gate on all engines (K2000 W1-W4 with
# bias mode off, K2000 with bias mode on and b = 0, random K = 4 matrices with random biases at n = 777 and 2048, and one set
# per problem of PROTOCOL_HW_MB_A4_verify.txt), phase 2 K2000 regression (X5, O4: trial ids and seed of the v6.4 run),
# phase 3 problem cohorts (PROTOCOL_HW_MB_A4_cohorts.txt). Requires the image programmed. Fresh result directory.
# Usage: SCA_VARIANT=mb4b NUM_ENGINES=6 BOARD_TAG=board_mb4b_e6_250mhz bash scripts/run_hw_mb4b.sh
# List formats: verify  key problem target_energy trials offset host-arguments...
#               cohorts key problem target_energy trials host-arguments...     (problem: matrix spec for --matrix)
set -uo pipefail
source /scratch/USER/sca_v80_20261003/scripts/environment.sh
[[ "$BOARD_TAG" =~ ^board_([a-z0-9]+_)*[0-9]+mhz$ ]] || exit 2
host="$TASK_ROOT/build/host_sca${SCA_VARIANT:+_$SCA_VARIANT}"; test -x "$host" || exit 2
E=${NUM_ENGINES:-6}
manifest="$TASK_ROOT/build/$BOARD_TAG/image_manifest.json"
uuid=$(python3 -c "import json,sys; print(json.load(open(sys.argv[1]))['logic_uuid'])" "$manifest")
mhz=$(python3 -c "import json,sys; print(json.load(open(sys.argv[1]))['clock']['actual_mhz'])" "$manifest")
out="$TASK_ROOT/results/hwmb_${BOARD_TAG}_a4"
if [ -e "$out" ]; then echo "fresh result directory required: $out"; exit 2; fi
mkdir -p "$out/verify" "$out/rounds" "$out/cohorts"
cp "$manifest" "$out/image_manifest.json"
cp "$TASK_ROOT/PROTOCOL_HW_MB_A4_verify.txt" "$TASK_ROOT/PROTOCOL_HW_MB_A4_cohorts.txt" "$out/"
sha256sum "$host" "$TASK_ROOT/src/host_sca_multi_mb4.cpp" "$TASK_ROOT/src/mb_common4.hpp" "$TASK_ROOT/src/sca_ref_bias.hpp" \
  "$TASK_ROOT/data/K2000.bin" "$TASK_ROOT"/data/a4_problems/* > "$out/inputs.sha256"
log="$out/run.log"; echo running > "$out/status"
note() { printf '%s %s\n' "$(date -u +%FT%TZ)" "$*" >> "$log"; }
run() { "$SNOWBALL_ROOT/build/scratch_exec" /scratch/USER --v80-ami "$host" \
          --uuid "$uuid" --clock-mhz "$mhz" --seed 20261004 --engines $E --keep-tables "$@" >> "$log" 2>&1 < /dev/null; }
finish() { note "exit $1"; echo "$1" > "$out/status"; exit "$1"; }
note "uuid $uuid clock_mhz $mhz engines $E host $(sha256sum "$host" | cut -c1-16)"
up() { echo $(( ( ($1 + $2 - 1) / $2 ) * $2 )); }
W1=$(up 8 $E); W2=$(up 8 $((2*E))); W3=$(up 4 $E); W4=$(up 8 $E); NR=$(up 2048 $E)
R777=rand:20261010:0.5:7:777:300; R2048=rand:20261011:0.5:7:2048:300
# ---- phase 1: exact gate (every trial bit-exact with sca_ref_bias.hpp run_trial_bias: spins, score word, flips, n_lin trace)
run --verify --t1 5 --output "$out/verify/W1" --q 4 --lambda 0 --t0 30 --steps 40 --trials $W1 --per-launch 1 --offset 900000 || finish 11
run --verify --t1 5 --output "$out/verify/W2" --q 4 --lambda 0.7 --ramp --t0 20 --steps 560 --trials $W2 --per-launch 2 --offset 900100 || finish 12
run --verify --t1 5 --output "$out/verify/W3" --q 8 --tec-jv -8 --t0 20 --steps 40 --trials $W3 --per-launch 1 --offset 900200 || finish 13
run --verify --t1 5 --output "$out/verify/W4" --q 8 --lambda 1.05 --ramp --t0 12 --steps 960 --trials $W4 --per-launch 1 --offset 900300 || finish 14
run --verify --t1 5 --bias-mode 1 --target-energy -1000000 --output "$out/verify/W1b0" --q 4 --lambda 0 --t0 30 --steps 40 --trials $W1 --per-launch 1 --offset 900000 || finish 15
run --verify --t1 5 --bias-mode 1 --target-energy -1000000 --output "$out/verify/W4b0" --q 8 --lambda 1.05 --ramp --t0 12 --steps 960 --trials $W4 --per-launch 1 --offset 900300 || finish 16
run --verify --matrix $R777 --target-energy -100000000 --output "$out/verify/R777_ons" --q 20 --lambda 0.7 --ramp --corr-scale auto --t0 100 --t1 20 --steps 40 --trials $((2*E)) --per-launch 2 --offset 960000 || finish 17
run --verify --matrix $R777 --target-energy -100000000 --output "$out/verify/R777_tec" --q 20 --tec-jv -20 --t0 100 --t1 20 --steps 400 --trials $E --per-launch 1 --offset 960100 || finish 18
run --verify --matrix $R2048 --target-energy -100000000 --output "$out/verify/R2048_ons" --q 10 --lambda 0.9 --ramp --corr-scale auto --t0 150 --t1 10 --steps 400 --trials $((2*E)) --per-launch 2 --offset 960200 || finish 19
run --verify --matrix $R2048 --target-energy -100000000 --output "$out/verify/R2048_tecT" --q 10 --tecT 1.0 --ramp --corr-scale auto --t0 150 --t1 10 --steps 300 --trials $E --per-launch 1 --offset 960300 || finish 20
gate=0
while read -r key prob target trials off args; do
  case "$key" in ''|\#*) continue;; esac
  # shellcheck disable=SC2086
  run --verify --matrix $prob --target-energy $target --output "$out/verify/$key" --trials $trials --per-launch 1 --offset $off $args
  r=$?; note "verify $key exit $r"; [ "$r" = 0 ] || gate=1
done < "$TASK_ROOT/PROTOCOL_HW_MB_A4_verify.txt"
[ "$gate" = 0 ] || finish 21
note "phase 1 PASS"
# ---- phase 2: K2000 regression (bias mode off), trial ids 0..NR-1 of the v6.4 run
run --cohort --t1 5 --output "$out/rounds/X5" --q 8 --tecT 1.75 --ramp --t0 15 --steps 280 --trials $NR --per-launch 1 --offset 0; note "rounds X5 exit $?"
run --cohort --t1 5 --output "$out/rounds/O4" --q 6 --lambda 0.9 --t0 12 --steps 360 --trials $NR --per-launch 1 --offset 0; note "rounds O4 exit $?"
note "phase 2 done"
# ---- phase 3: problem cohorts (frozen list)
hard=0
while read -r key prob target trials args; do
  case "$key" in ''|\#*) continue;; esac
  # shellcheck disable=SC2086
  run --cohort --matrix $prob --target-energy $target --output "$out/cohorts/$key" --trials $trials --per-launch 1 --offset 0 $args
  r=$?; note "cohort $key exit $r"; [ "$r" = 1 ] && hard=1
done < "$TASK_ROOT/PROTOCOL_HW_MB_A4_cohorts.txt"
finish $((hard ? 31 : 0))
