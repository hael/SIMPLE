
#include <stdio.h>
#include <stdlib.h>
#include <errno.h>

#include <cuComplex.h>
#include "cuda.h"


// Use cuCmulf from cuComplex.h

/* Define CUDA kernel that squares the input complex array */
__global__ void  mul2DComplex(cuFloatComplex *in1,cuFloatComplex *in2, cuFloatComplex *out, int Nsize, int Msize)
{

  unsigned int n = blockIdx.x * blockDim.x + threadIdx.x;
  unsigned int m = blockIdx.y * blockDim.y + threadIdx.y;

  if(n > Nsize || m > Msize) return;
  int grid_width = gridDim.x * blockDim.x;
  int index = m + (n * grid_width);

  out[index] = cuCmulf (in1[index] , in2[index]);

}
// math_functions in cuda.h
__global__ void  mul2DFloat(float *in1,float *in2, float *out, int Nsize, int Msize)
{

  unsigned int n = blockIdx.x * blockDim.x + threadIdx.x;
  unsigned int m = blockIdx.y * blockDim.y + threadIdx.y;

  if(n > Nsize || m > Msize) return;
  int grid_width = gridDim.x * blockDim.x;
  int index = m + (n * grid_width);

  out[index] = (in1[index] * in2[index]);

}


/*
   Fortran subroutine arguments are passed by references.
   call fun( array_a, array_b, N) will be mapped to
   function (*a, *b, *N);
*/
extern "C"
{
#define KERNELMUL2DFLOAT kernelMul2DFloat_
void kernelMul2DFloat(float *a, float *b, float *c, int *Np, int *Ns, int *Bsize)
{
  int block_size=*Bsize;
  cuFloat *a_d,*b_d,*c_d;
  int N=*Np;int M=*Ns;
  cudaSetDevice(0);

  /* Allocate complex array on device */
   cudaMalloc ((void **) &a_d , sizeof(cuFloat)*N*M );
   cudaMalloc ((void **) &b_d , sizeof(cuFloat)*N*M );
   cudaMalloc ((void **) &c_d , sizeof(cuFloat)*N*M );
  // if(cudaGetLastError() != cudaSuccess){
  // printf("%s\n",cudaGetErrorString(cudaGetLastError()));
  // exit(1);
  // }

  /* Copy array from host memory to device memory */
  cudaMemcpy( a_d, a,  sizeof(cuFloat)*N*M   ,cudaMemcpyHostToDevice);
  cudaMemcpy( b_d, b,  sizeof(cuFloat)*N*M   ,cudaMemcpyHostToDevice);
  // if( cudaMemcpy( a_d, a,  sizeof(cuFloatComplex)*N*M   ,cudaMemcpyHostToDevice) != cudaSuccess){
  // printf("%s\n",cudaGetErrorString(cudaGetLastError()));
  // exit(1);
  // }

  /* Compute execution configuration */
   dim3 dimBlock(block_size, 8);
   dim3 dimGrid ; //(N/dimBlock.x);
   dimGrid.x = (N + dimBlock.x - 1) / dimBlock.x;
   dimGrid.y = (M + dimBlock.y - 1) / dimBlock.y;


   if( N*N % block_size != 0 ) dimGrid.x+=1;

  /* Execute the kernel */
   mul2DFloat<<<dimGrid,dimBlock>>>(a_d,b_d,c_d,N,M);
  if( cudaGetLastError()!= cudaSuccess){
    printf("%s\n",cudaGetErrorString(cudaGetLastError()));
    exit(1);
  }

  /* Copy the result back */
  cudaMemcpy( c, c_d, sizeof(cuFloatComplex)*N*M,cudaMemcpyDeviceToHost);
  // if( cudaGetLastError()!= cudaSuccess){
  // printf("%s\n",cudaGetErrorString(cudaGetLastError()));
  // exit(1);
  // }

  /* Free memory on the device */
  cudaFree(a_d, b_d, c_d);

  return;
}
#define KERNELMUL2DCOMPLEX kernelMul2DComplex_
void kernelmul2dcomplex_(cuFloatComplex *a, cuFloatComplex *b, cuFloatComplex *c, int *Np, int *Ns, int *Bsize)
{
  int block_size=*Bsize;
  cuFloatComplex *a_d,*b_d;
  int N=*Np;int M=*Ns;
  cudaSetDevice(0);

  /* Allocate complex array on device */
   cudaMalloc ((void **) &a_d , sizeof(cuFloatComplex)*N*M*2 );
   cudaMalloc ((void **) &b_d , sizeof(cuFloatComplex)*N*M*2 );
  // if(cudaGetLastError() != cudaSuccess){
  // printf("%s\n",cudaGetErrorString(cudaGetLastError()));
  // exit(1);
  // }

  /* Copy array from host memory to device memory */
  cudaMemcpy( a_d, a,  sizeof(cuFloatComplex)*N*M   ,cudaMemcpyHostToDevice);
  cudaMemcpy( b_d, b,  sizeof(cuFloatComplex)*N*M   ,cudaMemcpyHostToDevice);
 // if( cudaMemcpy( a_d, a,  sizeof(cuFloatComplex)*N*M   ,cudaMemcpyHostToDevice) != cudaSuccess){
  // printf("%s\n",cudaGetErrorString(cudaGetLastError()));
  // exit(1);
  // }

  /* Compute execution configuration */
   dim3 dimBlock(block_size, 8);
   dim3 dimGrid ; //(N/dimBlock.x);
   dimGrid.x = (N + dimBlock.x - 1) / dimBlock.x;
   dimGrid.y = (M + dimBlock.y - 1) / dimBlock.y;


   if( N*N % block_size != 0 ) dimGrid.x+=1;

  /* Execute the kernel */
  mul2DComplex<<<dimGrid,dimBlock>>>(a_d,b_d,b_d,N,M);
  // if( cudaGetLastError()!= cudaSuccess){
  // printf("%s\n",cudaGetErrorString(cudaGetLastError()));
  // exit(1);
  // }

  /* Copy the result back */
  cudaMemcpy( c, b_d, sizeof(cuFloatComplex)*N*M,cudaMemcpyDeviceToHost);
  // if( cudaGetLastError()!= cudaSuccess){
  // printf("%s\n",cudaGetErrorString(cudaGetLastError()));
  // exit(1);
  // }

  /* Free memory on the device */
  cudaFree(a_d);

  return;
}

}
