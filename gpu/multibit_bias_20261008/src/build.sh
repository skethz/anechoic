#!/usr/bin/env bash
# Build on gpu-host (CUDA 13.3, GH200 sm_90). Run from the project root after sourcing environment.sh.
set -euo pipefail
mkdir -p build
nvcc -O3 -std=c++17 -arch=sm_90 -DSCA_LANES=256 -Xptxas -v -o build/sca_gpu_bias src/sca_gpu_bias.cu 2> build/ptxas_sca_gpu_bias.log
grep -c "0 bytes spill stores" build/ptxas_sca_gpu_bias.log | sed 's/^/kernels without spills: /'
g++ -O3 -march=native -std=c++17 -DSCA_LANES=256 -pthread -o build/check_ref_bias src/check_ref_bias.cpp
g++ -O3 -march=native -std=c++17 -DSCA_LANES=256 -o build/gen_ising src/gen_ising.cpp
c++ -O2 -std=c++17 -DSCA_LANES=256 src/test_ref_bias.cpp -o build/test_ref_bias
sha256sum src/* build/sca_gpu_bias build/check_ref_bias build/gen_ising build/test_ref_bias > build/BUILD.sha256
