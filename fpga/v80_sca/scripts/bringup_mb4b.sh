#!/usr/bin/env bash
# Unattended chain for PROTOCOL_HW_MB Amendment 4 (revision-3 K = 4 + bias image, 6 engines): wait until the build queue
# (scripts/chain_builds_mb3.sh) is done and the Amendment-3 chain has finished with the board (it ends by restoring v6.4),
# take the K = 4 + bias image that passed (250 MHz preferred, else 225 MHz), check the frozen Amendment-4 hashes, back up the
# staged image, program, run the amendment, measure power, restore the v6.4 headline image, analyse. No programming after
# 02:00 UTC on 9 October (hard stop for board work 04:00 UTC).
# Usage: bash scripts/bringup_mb4b.sh
set -uo pipefail
cd /scratch/USER/sca_v80_20261003
source scripts/environment.sh > /dev/null
log=logs/bringup_mb4b.log; V=mb4b
note() { printf '%s %s\n' "$(date -u +%FT%TZ)" "$*" >> "$log"; }
ok() {
  [ "$(cat results/pipeline_$1_exit_status 2>/dev/null)" = 0 ] && \
  python3 -c "import json,sys; d=json.load(open('build/$1/timing_signoff.json')); sys.exit(0 if d['status'].startswith('PASS') else 1)" 2>/dev/null
}
note start
sha256sum -c PROTOCOL_HW_MB_A4.sha256 > /dev/null 2>&1 || { note "A4 files changed before start; stopping"; exit 7; }
until ok board_mb4b_e6_250mhz || grep -q " done$" logs/chain_builds_mb3.log 2>/dev/null; do sleep 60; done
if ok board_mb4b_e6_250mhz; then T=board_mb4b_e6_250mhz; elif ok board_mb4b_e6_225mhz; then T=board_mb4b_e6_225mhz
else note "no K=4+bias image passed (pipeline or timing); board untouched"; exit 1; fi
note "image $T"
# the board is free once no K = 2 (Amendment 3) chain is running; every such chain ends by restoring v6.4
while pgrep -u "$(id -u)" -f "bringup_mb3.sh" > /dev/null; do sleep 60; done
note "no Amendment-3 chain running; last A3 log line: $(tail -1 logs/bringup_mb3_board_mb2n_e12_250mhz.log 2>/dev/null)"
test -x build/host_sca_$V || { note "host_sca_$V missing"; exit 3; }
sha256sum -c PROTOCOL_HW_MB_A4.sha256 >> "$log" 2>&1 || { note "A4 files changed; stopping"; exit 4; }
[ "$(date -u +%s)" -lt "$(date -u -d '2026-10-09T02:00:00Z' +%s)" ] || { note "too close to the 04:00 UTC board stop (session plus restore can take 1.5 h); board untouched"; exit 9; }
if pgrep -u "$(id -u)" -f "host_sca_[a-z0-9]* .*--(verify|cohort|power-run)" > /dev/null; then note "board busy; stopping"; exit 8; fi
bk=results/staged_image_backup_$(date -u +%Y%m%dT%H%M%SZ); mkdir -p "$bk"
cp -p /scratch/USER/snowball_v80_20260928/programming/snowball_v80.pdi /scratch/USER/snowball_v80_20260928/programming/staged_image.json "$bk/"
sha256sum "$bk/snowball_v80.pdi" > "$bk/SHA256SUMS"
note "staged image backed up to $bk: $(cut -c1-64 "$bk/SHA256SUMS")"
grep -q 67e591f351985531c3b8340f9c25bb80bc6add80c98e026ef26ca75f21f770b4 "$bk/SHA256SUMS" || note "WARNING: staged image is not the v6.4 headline image"
note "programming"
SCA_VARIANT=$V BOARD_TAG=$T bash scripts/program_board.sh >> "$log" 2>&1 || { note "programming failed; restoring v6.4"; BOARD_TAG=board_v64_e12_250mhz bash scripts/program_board.sh >> "$log" 2>&1; note "restore exit $?"; exit 5; }
note "A4 run"
SCA_VARIANT=$V NUM_ENGINES=6 BOARD_TAG=$T bash scripts/run_hw_mb4b.sh >> "$log" 2>&1; r=$?; note "A4 exit $r"
if [ "$r" = 0 ] || [ "$r" = 31 ]; then
  note "power measurement"
  SCA_VARIANT=$V NUM_ENGINES=6 BOARD_TAG=$T LOAD_S=150 IDLE_S=60 PER_LAUNCH=512 POWER_CFG=$TASK_ROOT/PROTOCOL_HW_MB_A4_power.txt \
    bash scripts/measure_power_mb3.sh >> "$log" 2>&1
  note "power exit $?"
else
  note "A4 did not complete; skipping power"
fi
note "restoring v6.4 headline image"
BOARD_TAG=board_v64_e12_250mhz bash scripts/program_board.sh >> "$log" 2>&1; note "restore exit $?"
R=results/hwmb_${T}_a4
python3 scripts/analyze_hw_mb_a4.py $R results/hwv6_board_v64_e12_250mhz PROTOCOL_HW_MB_A4_predictions.json > $R/analysis.txt 2>&1
note "analysis exit $? ($(head -1 $R/analysis.txt))"
P=$(ls -td results/power_${T}_* 2>/dev/null | head -1)
[ -n "$P" ] && { python3 scripts/analyze_power_mb.py "$P" $R/a4_summary.json > "$P/analysis.txt" 2>&1; note "power analysis exit $?"; }
note "done"
