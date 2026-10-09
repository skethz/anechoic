#!/bin/bash
# Copy results, logs and calibration outputs back from gpu-host into this folder (never deletes local files).
L="$(cd "$(dirname "$0")/." && pwd)"
R=USER@gpu-host:/scratch/USER/anechoic_cpu_20261008/research/reaim_benchmarks_20261008
for d in calib logs results hw_export precision; do
  rsync -a $R/$d/ $L/$d/ 2>/dev/null
done
rsync -a $R/verify_hope.json $R/verify_hw.json $R/repro_hope.json $L/ 2>/dev/null
rsync -a "$R/data/" "$L/data/" 2>/dev/null
true
