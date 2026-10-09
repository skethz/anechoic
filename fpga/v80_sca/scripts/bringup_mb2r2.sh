#!/usr/bin/env bash
# PROTOCOL_HW_MB Amendment 2: unattended validation of the revision-2 image (cut-shift fix). Waits for its build, requires
# timing PASS and all frozen hashes (main protocol, Amendments 1 and 2), backs up the staged image, programs board_mb2r2,
# then runs, unchanged: the main protocol phases 1-3 (run_hw_mb.sh), Amendment 1 (run_hw_mb_a1.sh), plus a negative-score
# exact-gate set; then restores the v6.4 headline image. Stops at the first failure (the r2 image then stays programmed).
# Usage: bash scripts/bringup_mb2r2.sh
set -uo pipefail
cd /scratch/USER/sca_v80_20261003
source scripts/environment.sh > /dev/null
T=board_mb2r2_e12_250mhz; st=results/pipeline_${T}_exit_status; log=logs/bringup_mb2r2_${T}.log
note() { printf '%s %s\n' "$(date -u +%FT%TZ)" "$*" >> "$log"; }
note start
for h in PROTOCOL_HW_MB.sha256 PROTOCOL_HW_MB_A1.sha256 PROTOCOL_HW_MB_A2.sha256; do
  sha256sum -c $h > /dev/null 2>&1 || { note "frozen files changed ($h); stopping"; exit 4; }; done
until [ -f "$st" ] && [ "$(cat $st)" != running ]; do sleep 60; done
# never overlap the first image's chain (it ends by restoring v6.4)
while pgrep -u "$(id -u)" -f "scripts/bringup_mb2_a1.sh|scripts/bringup_mb2.sh|scripts/run_hw_mb" > /dev/null; do sleep 60; done
note "pipeline status $(cat $st)"
[ "$(cat $st)" = 0 ] || { note "pipeline failed; board untouched"; exit 1; }
python3 -c "import json,sys; d=json.load(open('build/$T/timing_signoff.json')); print(d['status'], d['setup_slack_ns'], d['hold_slack_ns']); sys.exit(0 if d['status'].startswith('PASS') else 1)" >> "$log" 2>&1 \
  || { note "timing not PASS; board untouched"; exit 2; }
test -x build/host_sca_mb2r2 || { note "host_sca_mb2r2 missing"; exit 3; }
for h in PROTOCOL_HW_MB.sha256 PROTOCOL_HW_MB_A1.sha256 PROTOCOL_HW_MB_A2.sha256; do
  sha256sum -c $h >> "$log" 2>&1 || { note "frozen files changed ($h); stopping"; exit 4; }; done
bk=results/staged_image_backup_$(date -u +%Y%m%dT%H%M%SZ); mkdir -p "$bk"
cp -p /scratch/USER/snowball_v80_20260928/programming/snowball_v80.pdi /scratch/USER/snowball_v80_20260928/programming/staged_image.json "$bk/"
sha256sum "$bk/snowball_v80.pdi" > "$bk/SHA256SUMS"; note "staged image backed up to $bk: $(cut -c1-64 "$bk/SHA256SUMS")"
note "programming"
BOARD_TAG=$T bash scripts/program_board.sh >> "$log" 2>&1 || { note "programming failed"; exit 5; }
note "main protocol phases 1-3 on r2"
SCA_VARIANT=mb2r2 NUM_ENGINES=12 BOARD_TAG=$T GSET_CFG=$TASK_ROOT/PROTOCOL_HW_MB_gset.txt bash scripts/run_hw_mb.sh >> "$log" 2>&1; r=$?; note "main exit $r"
[ "$r" = 0 ] || [ "$r" = 21 ] || { note "main protocol phases failed; stopping (r2 stays programmed)"; exit 6; }
note "amendment 1 on r2"
SCA_VARIANT=mb2r2 NUM_ENGINES=12 BOARD_TAG=$T bash scripts/run_hw_mb_a1.sh >> "$log" 2>&1; note "A1 exit $?"
note "negative-score exact gate"
SCA_VARIANT=mb2r2 NUM_ENGINES=12 BOARD_TAG=$T bash scripts/run_hw_mb_a2neg.sh >> "$log" 2>&1; note "A2 neg exit $?"
note "restoring v6.4 headline image"
BOARD_TAG=board_v64_e12_250mhz bash scripts/program_board.sh >> "$log" 2>&1; note "restore exit $?"
note "done"
