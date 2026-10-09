#!/usr/bin/env bash
# Unattended multi-bit bring-up (PROTOCOL_HW_MB.md), modelled on bringup_v64.sh: wait for the build, require timing PASS and
# the frozen protocol hashes, back up the currently staged image (the v6.4 headline image, SHA-256 recorded) so it can be
# restored, program the multi-bit image (scoped sudo via program_board.sh), run the protocol, then measure board power.
# Stops at the first failure. Restoring v6.4 afterwards: BOARD_TAG=board_v64_e12_250mhz bash scripts/program_board.sh
# (its manifest-verified PDI stays in build/board_v64_e12_250mhz; the backup copy is a second, independent record).
# Usage: BOARD_TAG=board_mb2_e12_250mhz GSET_CFG=<frozen list> POWER_CFG=<frozen list> bash scripts/bringup_mb2.sh
set -uo pipefail
cd /scratch/USER/sca_v80_20261003
source scripts/environment.sh
T=${BOARD_TAG:-board_mb2_e12_250mhz}; st=results/pipeline_${T}_exit_status; log=logs/bringup_mb2_${T}.log
note() { printf '%s %s\n' "$(date -u +%FT%TZ)" "$*" >> "$log"; }
note start
until [ "$(cat $st 2>/dev/null)" != running ]; do sleep 60; done
note "pipeline status $(cat $st)"
[ "$(cat $st)" = 0 ] || { note "pipeline failed; board untouched (v6.4 image stays programmed)"; exit 1; }
python3 -c "import json,sys; d=json.load(open('build/$T/timing_signoff.json')); print(d['status'], d['setup_slack_ns'], d['hold_slack_ns']); sys.exit(0 if d['status'].startswith('PASS') else 1)" >> "$log" 2>&1 \
  || { note "timing not PASS; board untouched"; exit 2; }
test -x build/host_sca_mb2 || { note "host_sca_mb2 missing"; exit 3; }
sha256sum -c PROTOCOL_HW_MB.sha256 >> "$log" 2>&1 || { note "frozen protocol files changed; stopping"; exit 4; }
# back up the staged image (fixed path of the sudo rule) before it is replaced
bk=results/staged_image_backup_$(date -u +%Y%m%dT%H%M%SZ); mkdir -p "$bk"
cp -p /scratch/USER/snowball_v80_20260928/programming/snowball_v80.pdi /scratch/USER/snowball_v80_20260928/programming/staged_image.json "$bk/"
sha256sum "$bk/snowball_v80.pdi" > "$bk/SHA256SUMS"
note "staged image backed up to $bk: $(cut -c1-64 "$bk/SHA256SUMS") (v6.4 headline 67e591f3...)"
grep -q 67e591f351985531c3b8340f9c25bb80bc6add80c98e026ef26ca75f21f770b4 "$bk/SHA256SUMS" || note "WARNING: staged image is not the v6.4 headline image"
note "programming"
SCA_VARIANT=mb2 BOARD_TAG=$T bash scripts/program_board.sh >> "$log" 2>&1 || { note "programming failed"; exit 5; }
note "protocol run"
SCA_VARIANT=mb2 NUM_ENGINES=12 BOARD_TAG=$T GSET_CFG=${GSET_CFG:-} bash scripts/run_hw_mb.sh >> "$log" 2>&1; r=$?; note "protocol exit $r"
[ "$r" = 0 ] || [ "$r" = 21 ] || { note "protocol did not complete; skipping power"; exit 6; }
note "power measurement"
SCA_VARIANT=mb2 NUM_ENGINES=12 BOARD_TAG=$T LOAD_S=150 IDLE_S=60 PER_LAUNCH=512 POWER_CFG=${POWER_CFG:-} bash scripts/measure_power_mb.sh >> "$log" 2>&1
note "power exit $?"
# RESTORE_V64=1: after a completed protocol, program the v6.4 headline image again (same scoped command, manifest-verified)
if [ "${RESTORE_V64:-0}" = 1 ]; then
  note "restoring v6.4 headline image"
  BOARD_TAG=board_v64_e12_250mhz bash scripts/program_board.sh >> "$log" 2>&1; note "restore exit $?"
fi
note "done"
