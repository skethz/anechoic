#!/usr/bin/env bash
# Unattended v6.4 bring-up (PROTOCOL_HW_V6 Amendment 4): wait for the build, require timing PASS, program (scoped sudo via
# program_board.sh), run the frozen protocol with resident tables, then measure board power. Stops at the first failure.
set -uo pipefail
cd /scratch/USER/sca_v80_20261003
source scripts/environment.sh
T=board_v64_e12_250mhz; st=results/pipeline_${T}_exit_status; log=logs/bringup_v64.log
note() { printf '%s %s\n' "$(date -u +%FT%TZ)" "$*" >> "$log"; }
note start
until [ "$(cat $st 2>/dev/null)" != running ]; do sleep 60; done
note "pipeline status $(cat $st)"
[ "$(cat $st)" = 0 ] || { note "pipeline failed; v6.2 result stands"; exit 1; }
python3 -c "import json,sys; d=json.load(open('build/$T/timing_signoff.json')); print(d['status'], d['setup_slack_ns'], d['hold_slack_ns']); sys.exit(0 if d['status'].startswith('PASS') else 1)" >> "$log" 2>&1 \
  || { note "timing not PASS; v6.2 result stands"; exit 2; }
test -x build/host_sca_v64 || { note "host_sca_v64 missing"; exit 3; }
note "programming"
SCA_VARIANT=v64 BOARD_TAG=$T bash scripts/program_board.sh >> "$log" 2>&1 || { note "programming failed"; exit 4; }
note "protocol run"
KEEP_TABLES=1 SCA_VARIANT=v64 NUM_ENGINES=12 BOARD_TAG=$T bash scripts/run_hw_v6.sh >> "$log" 2>&1; note "protocol exit $?"
sfx=_v64; DRV=$(dirname "$(find "$TASK_ROOT/build/hls_sca${sfx}/solution/impl/misc/drivers" -name xsca_v80_hw.h | head -1)")
"$SNOWBALL_ROOT/build/scratch_exec" /scratch/USER g++ -O3 -std=c++17 -Wall -Wextra -DSCA_LANES=256 \
  -include "$TASK_ROOT/build/sca_timer_offsets${sfx}.h" -I"$DRV" -I"$TASK_ROOT/src" -I"$AVED_ROOT/sw/AMI/api/include" \
  src/host_sca_multi_v6_power.cpp /usr/local/lib/libami.a -pthread -o build/host_sca_v64_power >> "$log" 2>&1 || { note "power host build failed"; exit 5; }
note "power measurement"
KEEP_TABLES=1 SCA_VARIANT=v64 NUM_ENGINES=12 BOARD_TAG=$T LOAD_S=150 IDLE_S=60 PER_LAUNCH=512 bash scripts/measure_power_v6.sh >> "$log" 2>&1
note "done $?"
