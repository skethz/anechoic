#!/usr/bin/env bash
# PROTOCOL_HW_MB Amendment 3 (revision-3 K = 2 image, runtime n): phase 1 exact gate on all engines (K2000 W1-W4 and the
# G22/G32 sets of the main protocol, the Amendment-2 negative-score sets, one set per G-set instance of every size), phase 2
# K2000 regression (X5, O4: trial ids, seed and launch shape of the v6.4 12-engine run), phase 3 G-set cohorts
# (PROTOCOL_HW_MB_A3_gset.txt, all 51 instances). Requires the image programmed. Fresh result directory.
# Usage: SCA_VARIANT=mb2n NUM_ENGINES=12 BOARD_TAG=board_mb2n_e12_250mhz bash scripts/run_hw_mb3.sh
set -uo pipefail
source /scratch/USER/sca_v80_20261003/scripts/environment.sh
[[ "$BOARD_TAG" =~ ^board_([a-z0-9]+_)*[0-9]+mhz$ ]] || exit 2
host="$TASK_ROOT/build/host_sca${SCA_VARIANT:+_$SCA_VARIANT}"; test -x "$host" || exit 2
E=${NUM_ENGINES:-12}
manifest="$TASK_ROOT/build/$BOARD_TAG/image_manifest.json"
uuid=$(python3 -c "import json,sys; print(json.load(open(sys.argv[1]))['logic_uuid'])" "$manifest")
mhz=$(python3 -c "import json,sys; print(json.load(open(sys.argv[1]))['clock']['actual_mhz'])" "$manifest")
out="$TASK_ROOT/results/hwmb_${BOARD_TAG}_a3"
if [ -e "$out" ]; then echo "fresh result directory required: $out"; exit 2; fi
mkdir -p "$out/verify" "$out/rounds" "$out/gset"
cp "$manifest" "$out/image_manifest.json"
cp "$TASK_ROOT/PROTOCOL_HW_MB_A3_gset.txt" "$TASK_ROOT/PROTOCOL_HW_MB_A3_verify.txt" "$out/"
sha256sum "$host" "$TASK_ROOT/src/host_sca_multi_mb3.cpp" "$TASK_ROOT/src/mb_common3.hpp" "$TASK_ROOT/src/sca_ref.hpp" "$TASK_ROOT/src/sca_ref_bias.hpp" \
  "$TASK_ROOT/data/K2000.bin" "$TASK_ROOT"/data/gset/G* > "$out/inputs.sha256"
log="$out/run.log"; echo running > "$out/status"
note() { printf '%s %s\n' "$(date -u +%FT%TZ)" "$*" >> "$log"; }
run() { "$SNOWBALL_ROOT/build/scratch_exec" /scratch/USER --v80-ami "$host" \
          --uuid "$uuid" --clock-mhz "$mhz" --seed 20261004 --engines $E --keep-tables "$@" >> "$log" 2>&1 < /dev/null; }
finish() { note "exit $1"; echo "$1" > "$out/status"; exit "$1"; }
note "uuid $uuid clock_mhz $mhz engines $E host $(sha256sum "$host" | cut -c1-16)"
up() { echo $(( ( ($1 + $2 - 1) / $2 ) * $2 )); }
W1=$(up 8 $E); W2=$(up 8 $((2*E))); W3=$(up 4 $E); W4=$(up 8 $E); NR=$(up 2048 $E)
G22=gset:$TASK_ROOT/data/gset/G22; G32=gset:$TASK_ROOT/data/gset/G32
# ---- phase 1: exact gate (every trial bit-exact with sca_ref_bias.hpp run_trial_bias: spins, score word, flips, n_lin trace)
run --verify --t1 5 --output "$out/verify/W1" --q 4 --lambda 0 --t0 30 --steps 40 --trials $W1 --per-launch 1 --offset 900000 || finish 11
run --verify --t1 5 --output "$out/verify/W2" --q 4 --lambda 0.7 --ramp --t0 20 --steps 560 --trials $W2 --per-launch 2 --offset 900100 || finish 12
run --verify --t1 5 --output "$out/verify/W3" --q 8 --tec-jv -8 --t0 20 --steps 40 --trials $W3 --per-launch 1 --offset 900200 || finish 13
run --verify --t1 5 --output "$out/verify/W4" --q 8 --lambda 1.05 --ramp --t0 12 --steps 960 --trials $W4 --per-launch 1 --offset 900300 || finish 14
run --verify --matrix $G22 --target 13226 --output "$out/verify/G22a" --q 7 --lambda 0.9 --ramp --t0 4 --t1 0.5 --steps 600 --corr-scale auto --trials $((2*E)) --per-launch 2 --offset 910000 || finish 15
run --verify --matrix $G22 --target 13226 --output "$out/verify/G22b" --q 7 --tecT 1.5 --ramp --t0 4 --t1 0.5 --steps 600 --corr-scale auto --trials $E --per-launch 1 --offset 910100 || finish 16
run --verify --matrix $G32 --target 1396 --output "$out/verify/G32a" --q 2 --lambda 0.9 --ramp --t0 4 --t1 0.5 --steps 600 --corr-scale auto --trials $((2*E)) --per-launch 2 --offset 920000 || finish 17
run --verify --matrix $G32 --target 1396 --output "$out/verify/G32b" --q 2 --tec-jv -1 --t0 4 --t1 0.5 --steps 40 --trials $E --per-launch 1 --offset 920100 || finish 18
run --verify --matrix $G22 --target 13226 --output "$out/verify/G22neg_plain" --q 1 --lambda 0 --t0 3 --t1 0.3 --steps 40 --trials $E --per-launch 1 --offset 940000 || finish 19
run --verify --matrix $G22 --target 13226 --output "$out/verify/G22neg_ons" --q 1 --lambda 0.9 --ramp --t0 2 --t1 0.3 --steps 40 --corr-scale auto --trials $E --per-launch 1 --offset 940100 || finish 20
gate=0
while read -r key inst target trials off args; do
  case "$key" in ''|\#*) continue;; esac
  # shellcheck disable=SC2086
  run --verify --matrix gset:$TASK_ROOT/data/gset/$inst --target $target --output "$out/verify/$key" --trials $trials --per-launch 1 --offset $off $args
  r=$?; note "verify $key exit $r"; [ "$r" = 0 ] || gate=1
done < "$TASK_ROOT/PROTOCOL_HW_MB_A3_verify.txt"
[ "$gate" = 0 ] || finish 21
note "phase 1 PASS"
# ---- phase 2: K2000 regression, same configuration strings, trial ids 0..NR-1 and launch shape as the v6.4 run
run --cohort --t1 5 --output "$out/rounds/X5" --q 8 --tecT 1.75 --ramp --t0 15 --steps 280 --trials $NR --per-launch 1 --offset 0; note "rounds X5 exit $?"
run --cohort --t1 5 --output "$out/rounds/O4" --q 6 --lambda 0.9 --t0 12 --steps 360 --trials $NR --per-launch 1 --offset 0; note "rounds O4 exit $?"
note "phase 2 done"
# ---- phase 3: G-set cohorts (frozen list, all instances)
hard=0
while read -r key inst target trials args; do
  case "$key" in ''|\#*) continue;; esac
  # shellcheck disable=SC2086
  run --cohort --matrix gset:$TASK_ROOT/data/gset/$inst --target $target --output "$out/gset/$key" --trials $trials --per-launch 1 --offset 0 $args
  r=$?; note "gset $key exit $r"; [ "$r" = 1 ] && hard=1
done < "$TASK_ROOT/PROTOCOL_HW_MB_A3_gset.txt"
finish $((hard ? 31 : 0))
