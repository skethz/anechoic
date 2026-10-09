#!/usr/bin/env bash
# PROTOCOL_HW_V6 Amendment 5: STATICA's published K2000 points (SP_long, SP_short) and a same-session X5 control (paper O1) on
# the programmed 12 x v6.4 image, with the frozen v6.4 host (build/host_sca_v64) and run_hw_v6.sh's launch shape: seed 20261004,
# --t1 5, tables resident, one trial per launch, trial ids 0..2051 (171 rounds of 12). No programming. Fresh result directory.
# Usage: bash scripts/run_hw_v6_a5.sh
set -uo pipefail
source /scratch/USER/sca_v80_20261003/scripts/environment.sh > /dev/null
cd "$TASK_ROOT"
sha256sum -c PROTOCOL_HW_V6_A5.sha256 > /dev/null 2>&1 || { echo "A5 files changed; stopping"; exit 7; }
T=board_v64_e12_250mhz; host=$TASK_ROOT/build/host_sca_v64; E=12; NR=2052
manifest=build/$T/image_manifest.json
uuid=$(python3 -c "import json,sys; print(json.load(open(sys.argv[1]))['logic_uuid'])" "$manifest")
mhz=$(python3 -c "import json,sys; print(json.load(open(sys.argv[1]))['clock']['actual_mhz'])" "$manifest")
out=$TASK_ROOT/results/hwv6_${T}_a5
[ -e "$out" ] && { echo "fresh result directory required: $out"; exit 2; }
got=$("$SNOWBALL_ROOT/build/scratch_exec" /scratch/USER --v80-ami "$host" --probe | python3 -c "import json,sys; print(json.load(sys.stdin)['logic_uuid'])")
[ "$got" = "$uuid" ] || { echo "board does not hold the v6.4 image ($got); stopping, no programming"; exit 3; }
pgrep -u "$(id -u)" -f "host_sca_[a-z0-9]* .*--(verify|cohort|power-run)" > /dev/null && { echo "board busy; stopping"; exit 8; }
mkdir -p "$out/rounds"; cp "$manifest" "$out/image_manifest.json"
sha256sum "$host" src/sca_ref.hpp data/K2000.bin PROTOCOL_HW_V6_A5.md scripts/run_hw_v6_a5.sh scripts/analyze_hw_v6_a5.py > "$out/inputs.sha256"
log=$out/run.log; echo running > "$out/status"
note() { printf '%s %s\n' "$(date -u +%FT%TZ)" "$*" >> "$log"; }
run() { "$SNOWBALL_ROOT/build/scratch_exec" /scratch/USER --v80-ami "$host" \
          --uuid "$uuid" --clock-mhz "$mhz" --seed 20261004 --t1 5 --engines $E --keep-tables "$@" >> "$log" 2>&1 < /dev/null; }
note "uuid $uuid clock_mhz $mhz engines $E host $(sha256sum "$host" | cut -c1-16)"
declare -A CFG=(
  [X5]="--q 8 --tecT 1.75 --ramp --t0 15 --steps 280"
  [SP_long]="--q 4 --lambda 0 --t0 40 --steps 1560"
  [SP_short]="--q 4 --lambda 0 --t0 30 --steps 560"
)
hard=0
for key in X5 SP_long SP_short; do
  # shellcheck disable=SC2086
  run --cohort --output "$out/rounds/$key" ${CFG[$key]} --trials $NR --per-launch 1 --offset 0; r=$?
  note "rounds $key exit $r"; [ "$r" = 1 ] && hard=1
done
note "exit $((hard ? 21 : 0))"; echo $((hard ? 21 : 0)) > "$out/status"
python3 scripts/analyze_hw_v6_a5.py "$out" results/hwv6_${T} > "$out/analysis.txt" 2>&1; note "analysis exit $?"
