#!/usr/bin/env bash
# Board power of the programmed revision-3 multi-bit image (copy of measure_power_mb.sh recording the revision-3 host source;
# no sudo). (measure_power_mb.sh text follows) Board power of the programmed multi-bit image (copy of measure_power_v6.sh; no sudo). Phases: idle, then per configuration
# a sustained back-to-back load on all engines (host --power-run, no readback) with the AMI hwmon node sampled every 0.25 s
# (never `ami_tool sensors`: it stalls BAR access), idle again between loads.
# hwmon.txt columns: phase, unix time, Total_Power [uW], VCCINT current [mA], 12V_PEX current [mA]
# Usage: SCA_VARIANT=mb2 NUM_ENGINES=12 BOARD_TAG=board_mb2_e12_250mhz LOAD_S=150 [PER_LAUNCH=512] [POWER_CFG=<file>] bash scripts/measure_power_mb.sh
# POWER_CFG lines: key arguments...  (default: X5 on K2000)
set -uo pipefail
source /scratch/USER/sca_v80_20261003/scripts/environment.sh
host="$TASK_ROOT/build/host_sca_${SCA_VARIANT}"; test -x "$host" || exit 2
E=${NUM_ENGINES:-12}; LOAD_S=${LOAD_S:-150}; IDLE_S=${IDLE_S:-60}; PL=${PER_LAUNCH:-512}
manifest="$TASK_ROOT/build/$BOARD_TAG/image_manifest.json"
uuid=$(python3 -c "import json,sys; print(json.load(open(sys.argv[1]))['logic_uuid'])" "$manifest")
mhz=$(python3 -c "import json,sys; print(json.load(open(sys.argv[1]))['clock']['actual_mhz'])" "$manifest")
out="$TASK_ROOT/results/power_${BOARD_TAG}_$(date -u +%Y%m%dT%H%M%SZ)"
mkdir "$out" || exit 2
cp "$manifest" "$out/image_manifest.json"; sha256sum "$host" "$TASK_ROOT/src/host_sca_multi_mb3.cpp" "$TASK_ROOT/src/mb_common3.hpp" > "$out/host.sha256"
cfg=${POWER_CFG:-}; [ -n "$cfg" ] && cp "$cfg" "$out/power_configs.txt"
echo running > "$out/status"
sample() {  # $1 = phase label, $2 = seconds
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
      --t1 5 --engines "$E" --keep-tables --trials $((PL * E)) --per-launch "$PL" --offset 5000000 \
      --output "$out/unused_$key" --power-run "$LOAD_S" "$@" > "$out/load_$key.json" 2> "$out/load_$key.err" < /dev/null &
  local pid=$!
  sample "load_$key" "$LOAD_S"
  wait "$pid"; echo "$key exit $?" >> "$out/run.log"
}
sample idle0 "$IDLE_S"
if [ -z "$cfg" ]; then
  load X5 --q 8 --tecT 1.75 --ramp --t0 15 --steps 280
  sample idle1 "$IDLE_S"
else
  i=1
  while read -r key args; do
    case "$key" in ''|\#*) continue;; esac
    # shellcheck disable=SC2086
    load "$key" $args
    sample "idle$i" "$IDLE_S"; i=$((i + 1))
  done < "$cfg"
fi
echo 0 > "$out/status"; echo "$out"
