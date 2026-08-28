#!/bin/bash
# Builds the CUDA-C flex kernels (src/main/flex/cuda/) by setting USE_FLEX_CUDA=ON.
#
# This is the ONLY script that builds them. USE_FLEX_CUDA defaults to OFF and no
# other compile_*.sh sets it, so every other build leaves the device code out
# entirely -- no nvcc, no .cu in the source list, no CUDA symbol in the binary.
# In particular compile_gpu.sh is OpenMP nvptx offload, a different path that does
# NOT build these kernels.
#
# Toolchain overrides, for when nvcc/g++ are not on PATH or nvcc rejects the
# default host compiler:
#   NVCC=/usr/local/cuda/bin/nvcc CUDA_HOST=/usr/bin/g++ ./compile_flex_gpu.sh
# CUDA_ARCH defaults to 'native' (whatever card is in this machine); set it to
# cross-build for another card, e.g. CUDA_ARCH=75.
rm -rf build
mkdir build
cd build
cmake .. -DUSE_FLEX_CUDA=ON \
    ${NVCC:+-DCMAKE_CUDA_COMPILER=$NVCC} \
    ${CUDA_HOST:+-DCMAKE_CUDA_HOST_COMPILER=$CUDA_HOST} \
    ${CUDA_ARCH:+-DCMAKE_CUDA_ARCHITECTURES=$CUDA_ARCH}
make -j install
