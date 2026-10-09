#!/bin/bash
# Copy results, logs and outputs back from gpu-host (never overwrites the frozen, hashed local files).
L="$(cd "$(dirname "$0")/." && pwd)"
R=USER@gpu-host:/scratch/USER/anechoic_cpu_20261008/research/reaim_benchmarks_20261008
for d in calib logs results hw_export precision; do
  rsync -a $R/$d/ $L/$d/ 2>/dev/null
done
rsync -a $R/verify_hw_frozen.json $L/ 2>/dev/null
true
