
#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <math.h>
#include "cuComplex.h"
#include "cuda.h"
#include "simple_cuda.h"

// Math_functions in cuda.h
__global__ void  normvecangle(float *in1, float *in2, float *out, int Nsize)
{

    unsigned int n = blockIdx.x * blockDim.x + threadIdx.x;

    if(n <= Nsize) {

        out[n] = acos(in1[n] * in2[n] + in1[n + 1] * in2[n + 1] + in1[n + 2] * in2[n + 2]);
    }
}


/*
   Fortran subroutine arguments are passed by references.
   call fun( array_a, array_b, N) will be mapped to
   function (*a, *b, *N);
*/
// C=acos (A.B)

extern "C" void kerneleulerdist_(float *a, float *b, float *c, int *Np, int *Ns, int *Bsize)
{
    int block_size = *Bsize;
    float *a_d, *b_d, *c_d;
    int N = *Np; int M = 3;
    cudaSetDevice(0);

    /* Allocate complex array on device */
    cudaMalloc((void **) &a_d , sizeof(float)*N * M);
    cudaMalloc((void **) &b_d , sizeof(float)*N * M);
    cudaMalloc((void **) &c_d , sizeof(float)*N * M);

    cudaCheckErrors(" kerneleulerdist Malloc failed ");

    /* Copy array from host memory to device memory */
    cudaMemcpy(a_d, a,  sizeof(float)*N * M   , cudaMemcpyHostToDevice);
    cudaMemcpy(b_d, b,  sizeof(float)*N * M   , cudaMemcpyHostToDevice);
    cudaCheckErrors(" kerneleulerdist Mencpy failed ");

    /* Compute execution configuration */
    dim3 dimBlock(block_size, 8);
    dim3 dimGrid ; //(N/dimBlock.x);
    dimGrid.x = (N + dimBlock.x - 1) / dimBlock.x;
    dimGrid.y = (M + dimBlock.y - 1) / dimBlock.y;


    if(N * M % block_size != 0) dimGrid.x += 1;

    /* Execute the kernel */
    normvecangle <<< dimGrid, dimBlock>>>(a_d, b_d, c_d, N);
    cudaCheckErrors(" acos_dot failed ");


    /* Copy the result back */
    cudaMemcpy(c, c_d, sizeof(float)*N * M, cudaMemcpyDeviceToHost);
    cudaCheckErrors(" kerneleulerdist Memcpy return failed ");

    /* Free memory on the device */
    cudaFree(a_d);
    cudaFree(b_d);
    cudaFree(c_d);
    cudaCheckErrors(" kerneleulerdist Free failed ");

    return;
}
