#!/usr/bin/env bash
# Board power of a programmed v6 image (no sudo; ami_tool sensors is read-only). Phases: idle, then for each configuration a
# sustained back-to-back load on all engines (host --power-run, no readback) with sensors sampled throughout, then idle again.
# hwmon.txt columns: phase, unix time, Total_Power [uW], VCCINT current [mA], 12V_PEX current [mA]
# Usage: SCA_VARIANT=v62 NUM_ENGINES=12 BOARD_TAG=board_v62_e12_275mhz LOAD_S=150 [PER_LAUNCH=512] bash scripts/measure_power_v6.sh
# (First run 2026-10-06T23:37Z used 64 trials per launch: device busy only 71%, kept as a record.)
set -uo pipefail
source /scratch/USER/sca_v80_20261003/scripts/environment.sh
host="$TASK_ROOT/build/host_sca_${SCA_VARIANT}_power"; test -x "$host" || exit 2
E=${NUM_ENGINES:-12}; LOAD_S=${LOAD_S:-150}; IDLE_S=${IDLE_S:-60}; PL=${PER_LAUNCH:-512}   # 512: device busy >~95%
manifest="$TASK_ROOT/build/$BOARD_TAG/image_manifest.json"
uuid=$(python3 -c "import json,sys; print(json.load(open(sys.argv[1]))['logic_uuid'])" "$manifest")
mhz=$(python3 -c "import json,sys; print(json.load(open(sys.argv[1]))['clock']['actual_mhz'])" "$manifest")
out="$TASK_ROOT/results/power_${BOARD_TAG}_$(date -u +%Y%m%dT%H%M%SZ)"
mkdir "$out" || exit 2
cp "$manifest" "$out/image_manifest.json"; sha256sum "$host" "$TASK_ROOT/src/host_sca_multi_v6_power.cpp" > "$out/host.sha256"
echo running > "$out/status"
sample() {  # $1 = phase label, $2 = seconds. AMI hwmon (driver-cached, microsecond reads): ami_tool sensors stalls BAR access
  local end=$(( $(date +%s) + $2 )) h=/sys/class/hwmon/hwmon4
  [ "$(cat $h/name)" = Alveo ] || { echo "hwmon4 is not the Alveo sensor node" >> "$out/run.log"; exit 3; }
  while [ "$(date +%s)" -lt "$end" ]; do
    printf '%s %s %s %s %s\n' "$1" "$(date -u +%s.%N)" "$(cat $h/power1_input)" "$(cat $h/curr1_input)" "$(cat $h/curr7_input)" >> "$out/hwmon.txt"
    sleep 0.25
  done
}
load() {  # $1 = key, rest = schedule arguments
  local key=$1; shift
  "$SNOWBALL_ROOT/build/scratch_exec" /scratch/USER --v80-ami "$host" --uuid "$uuid" --clock-mhz "$mhz" --seed 20261004 \
      --t1 5 --engines "$E" ${KEEP_TABLES:+--keep-tables} --trials $((PL * E)) --per-launch "$PL" --offset 5000000 \
      --output "$out/unused_$key" --power-run "$LOAD_S" "$@" > "$out/load_$key.json" 2> "$out/load_$key.err" &
  local pid=$!
  sample "load_$key" "$LOAD_S"
  wait "$pid"; echo "$key exit $?" >> "$out/run.log"
}
sample idle0 "$IDLE_S"
load X5 --q 8 --tecT 1.75 --ramp --t0 15 --steps 280
sample idle1 "$IDLE_S"
load X4 --q 8 --tecT 2.0 --ramp --t0 15 --steps 760
sample idle2 "$IDLE_S"
load O1 --q 8 --lambda 1.05 --ramp --t0 12 --steps 960
sample idle3 "$IDLE_S"
echo 0 > "$out/status"; echo "$out"
