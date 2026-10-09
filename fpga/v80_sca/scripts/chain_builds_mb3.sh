#!/usr/bin/env bash
# Build order for the revision-3 images on fpga-host (one Vivado build at a time: a 12-engine build peaks near 60 GB of 62 GB).
#  1. wait for board_mb2n_e12_250mhz (K = 2, 12 engines, already running);
#  2. if it failed (pipeline or timing): K = 2 again at 225 MHz first (the K = 2 board result has priority);
#  3. K = 4 + bias, 6 engines at 250 MHz; if that fails (pipeline or timing): K = 4 + bias at 225 MHz.
# Never programs. Usage: bash scripts/chain_builds_mb3.sh
set -uo pipefail
cd /scratch/USER/sca_v80_20261003 && source scripts/environment.sh > /dev/null
log=logs/chain_builds_mb3.log
note() { printf '%s %s\n' "$(date -u +%FT%TZ)" "$*" >> "$log"; }
ok() {  # $1 = board tag: pipeline status 0 and timing PASS
  [ "$(cat results/pipeline_$1_exit_status 2>/dev/null)" = 0 ] && \
  python3 -c "import json,sys; d=json.load(open('build/$1/timing_signoff.json')); sys.exit(0 if d['status'].startswith('PASS') else 1)" 2>/dev/null
}
wait_done() { until [ -f "results/pipeline_$1_exit_status" ] && [ "$(cat results/pipeline_$1_exit_status)" != running ]; do sleep 60; done; }
note start
wait_done board_mb2n_e12_250mhz
if ok board_mb2n_e12_250mhz; then note "K=2 250 MHz OK"
else
  note "K=2 250 MHz failed (status $(cat results/pipeline_board_mb2n_e12_250mhz_exit_status)); building K=2 at 225 MHz"
  KERNEL_MHZ=225 bash scripts/build_mb3.sh k2 >> "$log" 2>&1; note "K=2 225 MHz build exit $? ok=$(ok board_mb2n_e12_225mhz && echo 1 || echo 0)"
fi
note "building K=4+bias at 250 MHz"
KERNEL_MHZ=250 bash scripts/build_mb3.sh k4b >> "$log" 2>&1; note "K=4b 250 MHz build exit $?"
if ok board_mb4b_e6_250mhz; then note "K=4b 250 MHz OK"
else
  note "K=4b 250 MHz failed; building K=4+bias at 225 MHz"
  KERNEL_MHZ=225 bash scripts/build_mb3.sh k4b >> "$log" 2>&1; note "K=4b 225 MHz build exit $? ok=$(ok board_mb4b_e6_225mhz && echo 1 || echo 0)"
fi
note done
