
#include <stdio.h>
#include <stdlib.h>
#include <errno.h>

#include "cuComplex.h"
#include "cuda.h"
#include "simple_cuda.h"

// Use cuCmulf from cuComplex.h
/* Define CUDA kernel that squares the input complex array */
__global__ void  mul1DComplex(cuComplex *in1, cuComplex *in2, cuComplex *out, int N)
{
    unsigned int index   = blockIdx.x * blockDim.x + threadIdx.x;
    if(index < N) {
        out[index] = cuCmulf(in1[index], in2[index]);
    }

}
/* Define CUDA kernel that squares the input complex array */
__global__ void  mul2DComplex(cuComplex *in1, cuComplex *in2, cuComplex *out, int Nsize, int Msize)
{

    unsigned int n = blockIdx.x * blockDim.x + threadIdx.x;
    unsigned int m = blockIdx.y * blockDim.y + threadIdx.y;

    if(n > Nsize || m > Msize) return;
    int grid_width = gridDim.x * blockDim.x;
    unsigned long index = m + (n * grid_width);

    out[index] = cuCmulf(in1[index] , in2[index]);

}
/* Define CUDA kernel that squares the input complex array */
__global__ void  muladd2DComplex(cuComplex *in1, cuComplex *in2, cuComplex *out, int Nsize, int Msize)
{

    unsigned int n = blockIdx.x * blockDim.x + threadIdx.x;
    unsigned int m = blockIdx.y * blockDim.y + threadIdx.y;

    if(n > Nsize || m > Msize) return;
    int grid_width = gridDim.x * blockDim.x;
    unsigned long index = m + (n * grid_width);

    out[index] = cuCfmaf(in1[index] , in2[index],  out[index]);

}

// math_functions in cuda.h
__global__ void  mul2DFloat(float *in1, float *in2, float *out, int Nsize, int Msize)
{

    unsigned int n = blockIdx.x * blockDim.x + threadIdx.x;
    unsigned int m = blockIdx.y * blockDim.y + threadIdx.y;

    if(n > Nsize || m > Msize) return;
    int grid_width = gridDim.x * blockDim.x;
    unsigned long index = m + (n * grid_width);

    out[index] = (in1[index] * in2[index]);

}
// math_functions in cuda.h
__global__ void  muladd2DFloat(float *in1, float *in2, float *out, int Nsize, int Msize)
{

    unsigned int n = blockIdx.x * blockDim.x + threadIdx.x;
    unsigned int m = blockIdx.y * blockDim.y + threadIdx.y;

    if(n > Nsize || m > Msize) return;
    int grid_width = gridDim.x * blockDim.x;
    unsigned long index = m + (n * grid_width);

    out[index] = (in1[index] * in2[index]) + out[index];

}


/*
   Fortran subroutine arguments are passed by references.
   call fun( array_a, array_b, N) will be mapped to
   function (*a, *b, *N);
*/


// C=A*B
//#define KERNELMUL2DCOMPLEX kernelmul2dcomplex_
extern "C" void kernelmul2dcomplex_(cuComplex *a, cuComplex *b, cuComplex *c, int *Ncol, int *Nrow, int *Bsize)
{
    int block_size = *Bsize;
    cuComplex *a_d, *b_d, *c_d;
    int N = *Ncol; int M = *Nrow;
//  printf("In kernelmul2dcomplex matsize %d, %d, matsize (bytes) %lu, %lu, %lu\n", N,M, sizeof(a),sizeof(b),sizeof(c));

    cudaSetDevice(0);

    /* Allocate complex array on device */
    cudaMalloc((void **) &a_d , sizeof(cuComplex)*N * M);
    cudaMalloc((void **) &b_d , sizeof(cuComplex)*N * M);
    cudaMalloc((void **) &c_d , sizeof(cuComplex)*N * M);
    cudaCheckErrors("mul2d Malloc failed.");


    /* Copy array from host memory to device memory */
    cudaMemcpy(a_d, a,  sizeof(cuComplex)*N * M   , cudaMemcpyHostToDevice);
    cudaMemcpy(b_d, b,  sizeof(cuComplex)*N * M   , cudaMemcpyHostToDevice);
    cudaCheckErrors("mul2d Memcpy failed.");

    /* Compute execution configuration */
    dim3 dimBlock(block_size, 8);
    dim3 dimGrid ;//(N/dimBlock.x);
    dimGrid.x = (N + dimBlock.x - 1) / dimBlock.x;
    dimGrid.y = (M + dimBlock.y - 1) / dimBlock.y;

    if(N * M % block_size != 0) dimGrid.x += 1;
//printf("dimGrid %d,%d,%d\n",dimGrid.x,dimGrid.y,dimGrid.z);
//printf("dimBlock %d,%d,%d\n",dimBlock.x,dimBlock.y,dimBlock.z);
    /* Execute the kernel */
    mul2DComplex <<< dimGrid, dimBlock>>>(a_d, b_d, c_d, N, M);
    cudaCheckErrors("mul2d mul2DComplex failed.");
    /* Copy the result back */
    cudaMemcpy(c, c_d, sizeof(cuComplex)*N * M, cudaMemcpyDeviceToHost);
    cudaCheckErrors("mul2d Memcpy failed.");
    /* Free memory on the device */
    cudaFree(a_d);
    cudaFree(b_d);
    cudaFree(c_d);

    return;
}

extern "C" void kernelmul1dcomplex_(cuComplex *a, cuComplex *b, cuComplex *c, int *Np)
{
    int block_size = 8;
    cuComplex *a_d, *b_d, *c_d;
    int N = *Np;
    cudaSetDevice(0);

    /* Allocate complex array on device */
    cudaMalloc((void **) &a_d , sizeof(cuComplex)*N);
    cudaMalloc((void **) &b_d , sizeof(cuComplex)*N);
    cudaMalloc((void **) &c_d , sizeof(cuComplex)*N);

    // if(cudaGetLastError() != cudaSuccess){
    // printf("%s\n",cudaGetErrorString(cudaGetLastError()));
    // exit(1);
    // }

    /* Copy array from host memory to device memory */
    cudaMemcpy(a_d, a,  sizeof(cuComplex)*N   , cudaMemcpyHostToDevice);
    cudaMemcpy(b_d, b,  sizeof(cuComplex)*N   , cudaMemcpyHostToDevice);
    cudaCheckErrors("mul1d Memcpy failed.");

    /* Compute execution configuration */
    dim3 dimBlock(block_size, 8);
    dim3 dimGrid(N / dimBlock.x);
    if(N % block_size != 0) dimGrid.x += 1;

    /* Execute the kernel */
    mul1DComplex <<< dimGrid, dimBlock>>>(a_d, b_d, c_d, N);
    cudaCheckErrors("mul1d mul1DComplex<> failed.");

    /* Copy the result back */
    cudaMemcpy(c, c_d, sizeof(cuComplex)*N, cudaMemcpyDeviceToHost);
    cudaCheckErrors("mul1d Memcpy failed.");

    /* Free memory on the device */
    cudaFree(a_d);
    cudaFree(b_d);
    cudaFree(c_d);

    return;
}




// C=A*B
//#define KERNELMUL2DFLOAT kernelmul2dfloat_
extern "C" void kernelmul2dfloat_(float *a, float *b, float *c, int *Np, int *Ns, int *Bsize)
{
    int block_size = *Bsize;
    float *a_d, *b_d, *c_d;
    int N = *Np; int M = *Ns;
    cudaSetDevice(0);

    /* Allocate complex array on device */
    cudaMalloc((void **) &a_d , sizeof(float)*N * M);
    cudaMalloc((void **) &b_d , sizeof(float)*N * M);
    cudaMalloc((void **) &c_d , sizeof(float)*N * M);
    cudaCheckErrors("mul2df Malloc failed.");


    /* Copy array from host memory to device memory */
    cudaMemcpy(a_d, a,  sizeof(float)*N * M   , cudaMemcpyHostToDevice);
    cudaMemcpy(b_d, b,  sizeof(float)*N * M   , cudaMemcpyHostToDevice);
    cudaCheckErrors("mul2df Memcpy failed.");


    /* Compute execution configuration */
    dim3 dimBlock(block_size, 8);
    dim3 dimGrid ; //(N/dimBlock.x);
    dimGrid.x = (N + dimBlock.x - 1) / dimBlock.x;
    dimGrid.y = (M + dimBlock.y - 1) / dimBlock.y;


    if(N * N % block_size != 0) dimGrid.x += 1;

    /* Execute the kernel */
    mul2DFloat <<< dimGrid, dimBlock>>>(a_d, b_d, c_d, N, M);
    cudaCheckErrors("mul2df kernel failed.");

    /* Copy the result back */
    cudaMemcpy(c, c_d, sizeof(float)*N * M, cudaMemcpyDeviceToHost);
    cudaCheckErrors("mul2df Memcpy failed.");

    /* Free memory on the device */
    cudaFree(a_d);
    cudaFree(b_d);
    cudaFree(c_d);

    return;
}

// D=A*B + C
//#define KERNELMULADD2DCOMPLEX kernelmul2dcomplex_
extern "C" void kernelmuladd2dcomplex_(cuComplex *a, cuComplex *b,cuComplex *c, cuComplex *d, int *Ncol, int *Nrow, int *Bsize)
{
    int block_size = *Bsize;
    cuComplex *a_d, *b_d, *c_d;
    int N = *Ncol; int M = *Nrow;
//  printf("In kernelmul2dcomplex matsize %d, %d, matsize (bytes) %lu, %lu, %lu\n", N,M, sizeof(a),sizeof(b),sizeof(c));

    cudaSetDevice(0);

    /* Allocate complex array on device */
    cudaMalloc((void **) &a_d , sizeof(cuComplex)*N * M);
    cudaMalloc((void **) &b_d , sizeof(cuComplex)*N * M);
    cudaMalloc((void **) &c_d , sizeof(cuComplex)*N * M);
    cudaCheckErrors("mul2d Malloc failed.");


    /* Copy array from host memory to device memory */
    cudaMemcpy(a_d, a,  sizeof(cuComplex)*N * M   , cudaMemcpyHostToDevice);
    cudaMemcpy(b_d, b,  sizeof(cuComplex)*N * M   , cudaMemcpyHostToDevice);
    cudaMemcpy(c_d, c,  sizeof(cuComplex)*N * M   , cudaMemcpyHostToDevice);

    cudaCheckErrors("mul2d Memcpy failed.");

    /* Compute execution configuration */
    dim3 dimBlock(block_size, 8);
    dim3 dimGrid ;//(N/dimBlock.x);
    dimGrid.x = (N + dimBlock.x - 1) / dimBlock.x;
    dimGrid.y = (M + dimBlock.y - 1) / dimBlock.y;

    if(N * M % block_size != 0) dimGrid.x += 1;
//printf("dimGrid %d,%d,%d\n",dimGrid.x,dimGrid.y,dimGrid.z);
//printf("dimBlock %d,%d,%d\n",dimBlock.x,dimBlock.y,dimBlock.z);
    /* Execute the kernel */
    muladd2DComplex <<< dimGrid, dimBlock>>>(a_d, b_d, c_d, N, M);
    cudaCheckErrors(" muladd2DComplex failed.");
    /* Copy the result back to D (not C)*/
    cudaMemcpy(d, c_d, sizeof(cuComplex)*N * M, cudaMemcpyDeviceToHost);
    cudaCheckErrors("mul2d Memcpy failed.");
    /* Free memory on the device */
    cudaFree(a_d);
    cudaFree(b_d);
    cudaFree(c_d);

    return;
}
