#include "cuda.h"
/* based on NVIDIA reduction talk
 * https://developer.download.nvidia.com/compute/cuda/1.1-Beta/x86_website/projects/reduction/doc/reduction.pdf
 *
 */
/*
template <unsigned int blockSize>
__global__ void reduce_int(int *g_idata, int *g_odata, unsigned int n)
{
    extern __shared__ int sdata[];
    unsigned int tid = threadIdx.x;
    unsigned int i = blockIdx.x*(blockSize*2) + tid;
    unsigned int gridSize = blockSize*2*gridDim.x;
    sdata[tid] = 0;
    while (i < n) { sdata[tid] += g_idata[i] + g_idata[i+blockSize]; i += gridSize; }
    __syncthreads();
    if (blockSize >= 512) { if (tid < 256) { sdata[tid] += sdata[tid + 256]; } __syncthreads(); }
    if (blockSize >= 256) { if (tid < 128) { sdata[tid] += sdata[tid + 128]; } __syncthreads(); }
    if (blockSize >= 128) { if (tid < 64) { sdata[tid] += sdata[tid + 64]; } __syncthreads(); }
    if (tid < 32) {
       if (blockSize >= 64) sdata[tid] += sdata[tid + 32];
       if (blockSize >= 32) sdata[tid] += sdata[tid + 16];
       if (blockSize >= 16) sdata[tid] += sdata[tid + 8];
       if (blockSize >= 8) sdata[tid] += sdata[tid + 4];
       if (blockSize >= 4) sdata[tid] += sdata[tid + 2];
a       if (blockSize >= 2) sdata[tid] += sdata[tid + 1];
    }
    if (tid == 0) g_odata[blockIdx.x] = sdata[0];
}
*/
template <unsigned int blockSize>
__global__ void reduce_float(float *g_idata, float *g_odata, unsigned int n)
{
    extern __shared__ float sdata[];
    unsigned int tid = threadIdx.x;
    unsigned int i = blockIdx.x * (blockSize * 2) + tid;
    unsigned int gridSize = blockSize * 2 * gridDim.x;
    sdata[tid] = 0.;
    while(i < n) {
        sdata[tid] += g_idata[i] + g_idata[i + blockSize];
        i += gridSize;
    }
    __syncthreads();
    if(blockSize >= 512) {
        if(tid < 256) {
            sdata[tid] += sdata[tid + 256];
        } __syncthreads();
    }
    if(blockSize >= 256) {
        if(tid < 128) {
            sdata[tid] += sdata[tid + 128];
        } __syncthreads();
    }
    if(blockSize >= 128) {
        if(tid < 64) {
            sdata[tid] += sdata[tid + 64];
        } __syncthreads();
    }
    if(tid < 32) {
        if(blockSize >= 64) sdata[tid] += sdata[tid + 32];
        if(blockSize >= 32) sdata[tid] += sdata[tid + 16];
        if(blockSize >= 16) sdata[tid] += sdata[tid + 8];
        if(blockSize >= 8) sdata[tid] += sdata[tid + 4];
        if(blockSize >= 4) sdata[tid] += sdata[tid + 2];
        if(blockSize >= 2) sdata[tid] += sdata[tid + 1];
    }
    if(tid == 0) g_odata[blockIdx.x] = sdata[0];
}
