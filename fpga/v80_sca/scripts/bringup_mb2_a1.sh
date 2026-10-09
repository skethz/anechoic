#!/usr/bin/env bash
# Unattended chain for PROTOCOL_HW_MB with Amendment 1: the frozen bring-up (scripts/bringup_mb2.sh, unchanged, run with
# RESTORE_V64=0: wait for the build, timing and hash gates, back up the staged v6.4 image, program, main protocol, power),
# then Amendment 1 (hash-checked), then restore the v6.4 headline image. Stops at the first failure; if anything before the
# restore fails, the multi-bit image stays programmed (working, documented in the logs).
# Usage: BOARD_TAG=board_mb2_e12_250mhz bash scripts/bringup_mb2_a1.sh
set -uo pipefail
cd /scratch/USER/sca_v80_20261003
source scripts/environment.sh > /dev/null
T=${BOARD_TAG:-board_mb2_e12_250mhz}; log=logs/bringup_mb2_${T}.log
note() { printf '%s %s\n' "$(date -u +%FT%TZ)" "$*" >> "$log"; }
sha256sum -c PROTOCOL_HW_MB_A1.sha256 > /dev/null 2>&1 || { note "A1: amendment files changed before start; stopping"; exit 7; }
BOARD_TAG=$T GSET_CFG=$TASK_ROOT/PROTOCOL_HW_MB_gset.txt POWER_CFG=$TASK_ROOT/PROTOCOL_HW_MB_power.txt RESTORE_V64=0 \
  bash scripts/bringup_mb2.sh
tail -1 "$log" | grep -q " done$" || { note "A1: main bring-up did not finish; amendment not run, no restore"; exit 8; }
grep -qE "protocol exit (0|21)$" "$log" || { note "A1: main protocol not completed; amendment not run, no restore"; exit 9; }
sha256sum -c PROTOCOL_HW_MB_A1.sha256 >> "$log" 2>&1 || { note "A1: amendment files changed; stopping"; exit 10; }
note "A1: amendment run"
SCA_VARIANT=mb2 NUM_ENGINES=12 BOARD_TAG=$T bash scripts/run_hw_mb_a1.sh >> "$log" 2>&1; note "A1: amendment exit $?"
note "restoring v6.4 headline image"
BOARD_TAG=board_v64_e12_250mhz bash scripts/program_board.sh >> "$log" 2>&1; note "restore exit $?"
note "A1: chain done"
