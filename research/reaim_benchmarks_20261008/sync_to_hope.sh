#!/bin/bash
# Mirror the code and data needed on gpu-host (scratch only). Hashed modules are copied unchanged.
L="$(cd "$(dirname "$0")/../.." && pwd)"
R=USER@gpu-host:/scratch/USER/anechoic_cpu_20261008
EXC="--exclude=__pycache__ --exclude=.numba_cache --exclude=*.pyc"
rsync -a $EXC $L/research/reaim_benchmarks_20261008/ $R/research/reaim_benchmarks_20261008/ --exclude=results --exclude=logs --exclude=remote_results
rsync -a $EXC --include='*.py' --include='bkv.json' --include='manifest.json' --include='data/***' --include='PROTOCOL*' --exclude='*' $L/research/gset_20261007/ $R/research/gset_20261007/
for d in algorithm_compare_20261007 ablation_20261005 reaim_reproduction_20261003 theory_ideas_20261003 statica_reproduction_20261003 iamp_evaluation_20261003; do
  rsync -a $EXC --include='*.py' --include='*.json' --exclude='*/' --exclude='*.npz' --exclude='*.log' $L/research/$d/ $R/research/$d/ 2>/dev/null
done
rsync -a $EXC $L/fpga/v80_sca/src/sca_ref.hpp $L/fpga/v80_sca/src/sca_ref_bias.hpp $L/fpga/v80_sca/src/test_ref_bias.cpp $L/fpga/v80_sca/src/mb_common.hpp $R/fpga/v80_sca/src/
rsync -a $EXC $L/gpu/multibit_20261007/src/jmat.hpp $R/gpu/multibit_20261007/src/
