#!/usr/bin/env bash
# Build on gpu-host (CUDA 13.3, GH200 sm_90). Run from the project root after sourcing environment.sh.
set -euo pipefail
mkdir -p build
nvcc -O3 -std=c++17 -arch=sm_90 -DSCA_LANES=256 -Xptxas -v -o build/sca_gpu src/sca_gpu.cu 2> build/ptxas_sca_gpu.log
grep -E "registers|spill" build/ptxas_sca_gpu.log | sed 's/^ptxas info    : //' | head -40
g++ -O3 -march=native -std=c++17 -DSCA_LANES=256 -pthread -o build/check_ref src/check_ref.cpp
sha256sum src/* build/sca_gpu build/check_ref > build/BUILD.sha256
