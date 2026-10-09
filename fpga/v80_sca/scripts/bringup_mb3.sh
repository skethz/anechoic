#!/usr/bin/env bash
# Unattended chain for PROTOCOL_HW_MB Amendment 3 (revision-3 K = 2 image board_mb2n_e12_250mhz): wait for the build, require
# timing PASS and the frozen amendment hashes, back up the staged image (the v6.4 headline image, SHA-256 recorded), program
# the revision-3 image (scoped sudo via program_board.sh), run the amendment (exact gate, K2000 regression, G-set), measure
# board power, restore the v6.4 headline image, then run the frozen analysis. Stops at the first failure before programming;
# after programming, the v6.4 image is restored whatever the protocol outcome.
# Usage: BOARD_TAG=board_mb2n_e12_250mhz bash scripts/bringup_mb3.sh
set -uo pipefail
cd /scratch/USER/sca_v80_20261003
source scripts/environment.sh > /dev/null
T=${BOARD_TAG:-board_mb2n_e12_250mhz}; V=mb2n; st=results/pipeline_${T}_exit_status; log=logs/bringup_mb3_${T}.log
note() { printf '%s %s\n' "$(date -u +%FT%TZ)" "$*" >> "$log"; }
note start
sha256sum -c PROTOCOL_HW_MB_A3.sha256 > /dev/null 2>&1 || { note "A3 files changed before start; stopping"; exit 7; }
until [ -f "$st" ] && [ "$(cat $st)" != running ]; do sleep 60; done
note "pipeline status $(cat $st)"
[ "$(cat $st)" = 0 ] || { note "pipeline failed; board untouched (v6.4 image stays programmed)"; exit 1; }
python3 -c "import json,sys; d=json.load(open('build/$T/timing_signoff.json')); print(d['status'], d['setup_slack_ns'], d['hold_slack_ns']); sys.exit(0 if d['status'].startswith('PASS') else 1)" >> "$log" 2>&1 \
  || { note "timing not PASS; board untouched"; exit 2; }
test -x build/host_sca_$V || { note "host_sca_$V missing"; exit 3; }
sha256sum -c PROTOCOL_HW_MB_A3.sha256 >> "$log" 2>&1 || { note "A3 files changed; stopping"; exit 4; }
[ "$(date -u +%s)" -lt "$(date -u -d '2026-10-09T03:00:00Z' +%s)" ] || { note "too close to the 04:00 UTC board stop; board untouched"; exit 9; }
if pgrep -u "$(id -u)" -f "host_sca_[a-z0-9]* .*--(verify|cohort|power-run)" > /dev/null; then note "board busy; stopping"; exit 8; fi
bk=results/staged_image_backup_$(date -u +%Y%m%dT%H%M%SZ); mkdir -p "$bk"
cp -p /scratch/USER/snowball_v80_20260928/programming/snowball_v80.pdi /scratch/USER/snowball_v80_20260928/programming/staged_image.json "$bk/"
sha256sum "$bk/snowball_v80.pdi" > "$bk/SHA256SUMS"
note "staged image backed up to $bk: $(cut -c1-64 "$bk/SHA256SUMS")"
grep -q 67e591f351985531c3b8340f9c25bb80bc6add80c98e026ef26ca75f21f770b4 "$bk/SHA256SUMS" || note "WARNING: staged image is not the v6.4 headline image"
note "programming"
SCA_VARIANT=$V BOARD_TAG=$T bash scripts/program_board.sh >> "$log" 2>&1 || { note "programming failed; restoring v6.4"; BOARD_TAG=board_v64_e12_250mhz bash scripts/program_board.sh >> "$log" 2>&1; note "restore exit $?"; exit 5; }
note "A3 run"
SCA_VARIANT=$V NUM_ENGINES=12 BOARD_TAG=$T bash scripts/run_hw_mb3.sh >> "$log" 2>&1; r=$?; note "A3 exit $r"
if [ "$r" = 0 ] || [ "$r" = 31 ]; then
  note "power measurement"
  SCA_VARIANT=$V NUM_ENGINES=12 BOARD_TAG=$T LOAD_S=150 IDLE_S=60 PER_LAUNCH=512 POWER_CFG=$TASK_ROOT/PROTOCOL_HW_MB_A3_power.txt \
    bash scripts/measure_power_mb3.sh >> "$log" 2>&1
  note "power exit $?"
else
  note "A3 did not complete; skipping power"
fi
note "restoring v6.4 headline image"
BOARD_TAG=board_v64_e12_250mhz bash scripts/program_board.sh >> "$log" 2>&1; note "restore exit $?"
R=results/hwmb_${T}_a3
python3 scripts/analyze_hw_mb_a3.py $R results/hwv6_board_v64_e12_250mhz results/hwmb_board_mb2_e12_250mhz_a1 PROTOCOL_HW_MB_A3_predictions.json > $R/analysis.txt 2>&1
note "analysis exit $? ($(head -1 $R/analysis.txt))"
P=$(ls -td results/power_${T}_* 2>/dev/null | head -1)
[ -n "$P" ] && { python3 scripts/analyze_power_mb.py "$P" $R/a3_summary.json > "$P/analysis.txt" 2>&1; note "power analysis exit $?"; }
note "done"
