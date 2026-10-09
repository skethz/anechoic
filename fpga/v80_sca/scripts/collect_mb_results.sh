#!/usr/bin/env bash
# After the multi-bit board chain: run the frozen analyses on fpga-host and pack the records (no PDI images) for the local
# results folder. Read-only on the result directories; writes tmp/mbpack_<tag>.tgz.
# Usage: BOARD_TAG=board_mb2_e12_250mhz bash scripts/collect_mb_results.sh
set -uo pipefail
cd /scratch/USER/sca_v80_20261003 && source scripts/environment.sh > /dev/null
T=${BOARD_TAG:-board_mb2_e12_250mhz}; R=results/hwmb_$T; A=results/hwmb_${T}_a1
P=$(ls -d results/power_${T}_* 2>/dev/null | tail -1)
[ -d "$R" ] && python3 scripts/analyze_hw_mb.py "$R" results/hwv6_board_v64_e12_250mhz --predictions PROTOCOL_HW_MB_predictions.json > "$R/analysis.txt" 2>&1
[ -n "$P" ] && [ -f "$R/mb_summary.json" ] && python3 scripts/analyze_power_mb.py "$P" "$R/mb_summary.json" > "$P/analysis.txt" 2>&1
[ -d "$A" ] && python3 scripts/analyze_hw_mb_a1.py "$A" PROTOCOL_HW_MB_A1_predictions.json > "$A/analysis.txt" 2>&1
pk=tmp/mbpack_$T; rm -rf "$pk"; mkdir -p "$pk/build"
for d in "$R" "$A" "$P" results/util_mb_$T; do [ -d "$d" ] && cp -r "$d" "$pk/"; done
cp logs/bringup_mb2_$T.log logs/pipeline_$T.log logs/pipeline_$T.env "$pk/" 2>/dev/null
cp results/programming_*_2026100[78]T*.log "$pk/" 2>/dev/null
for b in results/staged_image_backup_*; do [ -d "$b" ] && mkdir -p "$pk/$(basename $b)" && cp "$b/SHA256SUMS" "$b/staged_image.json" "$pk/$(basename $b)/"; done
for f in image_manifest.json timing_signoff.json kernel_clock.json timer_structure.json utilization.rpt route_status.rpt clocks.rpt; do
  cp "build/$T/$f" "$pk/build/" 2>/dev/null; done
cp build/$T/prj.runs/impl_1/place_report_utilization_0.rpt "$pk/build/" 2>/dev/null
grep -E "Timing Summary|SNOWBALL_|VERIFIED" build/$T/prj.runs/impl_1/runme.log logs/${T}_console.log 2>/dev/null | cut -c1-200 > "$pk/build/timing_lines.txt"
(cd tmp && tar czf "mbpack_$T.tgz" "mbpack_$T") && ls -la "tmp/mbpack_$T.tgz"
