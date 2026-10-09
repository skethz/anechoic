#!/usr/bin/env bash
# Adapted from /scratch/USER/anechoic_gpu_multibit_20261007/environment.sh (paths only).
R=/scratch/USER/anechoic_gpu_multibit_bias_20261008
export CUDA_VISIBLE_DEVICES=GPU-xxxxxxxx-xxxx-xxxx-xxxx-xxxxxxxxxxxx CUDA_DEVICE_ORDER=PCI_BUS_ID CUDA_CACHE_DISABLE=1
export CUDA_CACHE_PATH=$R/cache/cuda XDG_CACHE_HOME=$R/cache/xdg XDG_CONFIG_HOME=$R/cache/config XDG_STATE_HOME=$R/cache/state
export PYTHONPYCACHEPREFIX=$R/cache/python PYTHONDONTWRITEBYTECODE=1 TMPDIR=$R/tmp TMP=$R/tmp TEMP=$R/tmp HISTFILE=/dev/null
export PATH=/usr/local/cuda/bin:$PATH
mkdir -p $R/cache/cuda $R/cache/xdg $R/cache/config $R/cache/state $R/cache/python $R/tmp
cd $R
