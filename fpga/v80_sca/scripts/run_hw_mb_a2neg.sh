#!/usr/bin/env bash
# PROTOCOL_HW_MB Amendment 2: negative-score exact gate. G22 with self-coupling q = 1 drives trials into the all-equal
# period-2 states (cut 0, sum_i s_i h_i < 0), where the first multi-bit image's device cut word was wrong. 12 trials per set,
# every one verified bit for bit (spins, cut, flips, n_lin trace) against sca_ref.hpp run_trial(..., dense).
# Usage: SCA_VARIANT=mb2r2 NUM_ENGINES=12 BOARD_TAG=board_mb2r2_e12_250mhz bash scripts/run_hw_mb_a2neg.sh
set -uo pipefail
source /scratch/USER/sca_v80_20261003/scripts/environment.sh
host="$TASK_ROOT/build/host_sca${SCA_VARIANT:+_$SCA_VARIANT}"; test -x "$host" || exit 2
E=${NUM_ENGINES:-12}; manifest="$TASK_ROOT/build/$BOARD_TAG/image_manifest.json"
uuid=$(python3 -c "import json,sys; print(json.load(open(sys.argv[1]))['logic_uuid'])" "$manifest")
mhz=$(python3 -c "import json,sys; print(json.load(open(sys.argv[1]))['clock']['actual_mhz'])" "$manifest")
out="$TASK_ROOT/results/hwmb_${BOARD_TAG}_a2"; [ -e "$out" ] && { echo "fresh result directory required"; exit 2; }
mkdir -p "$out/verify"; cp "$manifest" "$out/"; log="$out/run.log"; echo running > "$out/status"
run() { "$SNOWBALL_ROOT/build/scratch_exec" /scratch/USER --v80-ami "$host" --uuid "$uuid" --clock-mhz "$mhz" --seed 20261004 \
          --engines $E --keep-tables --matrix gset:$TASK_ROOT/data/gset/G22 --target 13226 "$@" >> "$log" 2>&1 < /dev/null; }
run --verify --output "$out/verify/G22neg_plain" --q 1 --lambda 0 --t0 3 --t1 0.3 --steps 40 --trials $E --per-launch 1 --offset 940000; r1=$?
run --verify --output "$out/verify/G22neg_ons" --q 1 --lambda 0.9 --ramp --t0 2 --t1 0.3 --steps 40 --corr-scale auto --trials $E --per-launch 1 --offset 940100; r2=$?
r=$(( r1 || r2 )); echo $r > "$out/status"; exit $r
