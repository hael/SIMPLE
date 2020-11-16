#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <math.h>
#include "cuComplex.h"
#include "cuda.h"
#include "simple_cuda.h"

// Math_functions in cuda.h
// return real component of x * conjg(y)
__host__ __device__ static __inline__ float  rcorrelation(cuFloatComplex x, cuFloatComplex y)
{
     return (cuCrealf(x) * cuCrealf(y)) + (cuCimagf(x) * cuCimagf(y));

}

// return
__host__ __device__ static __inline__ float  complexsquared(cuFloatComplex a)
{
        float frac, sq;
        float x = fabsf(cuCrealf(a));
        float y = fabsf(cuCimagf(a));

        if( x < 1.0e-10f ) {
            sq = y*y;
        }else if( y < 1.0e-10f ) {
           sq = x*x;
        }else if( x > y ) {
            frac = y/x;
            sq = x*x*(1.+frac*frac);
        }else{
            frac = x/y;
            sq = y*y*(1.+frac*frac);
        }
    return sq;
}


// Reduction algorithm with csq in addition
template <unsigned int blockSize>
__global__ void sumComplexsquared(cuFloatComplex *g_idata, float *g_odata, unsigned int n)
{
    extern __shared__ float sdata[];
    unsigned int tid = threadIdx.x;

    unsigned int i = blockIdx.x*(blockSize*2) + tid;
    unsigned int gridSize = blockSize*2*gridDim.x;



    sdata[tid] = 0.;
    while (i < n) {
          sdata[tid] += complexsquared(g_idata[i]) + complexsquared(g_idata[i+blockSize]);
          i += gridSize;
    }
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
       if (blockSize >= 2) sdata[tid] += sdata[tid + 1];
    }
    if (tid == 0) g_odata[blockIdx.x] = sdata[0];
}


/*  Reduction algorithm with rcorrelation between two complex arrays in the addition */
template <unsigned int blockSize>
__global__ void sumrcorr(cuFloatComplex *g_idata1, cuFloatComplex *g_idata2, float *g_odata, unsigned int n)
{
    extern __shared__ float sdata[];
    unsigned int tid = threadIdx.x;
    unsigned int i = blockIdx.x*(blockSize*2) + tid;
    unsigned int gridSize = blockSize*2*gridDim.x;
    sdata[tid] = 0.;
    while (i < n) {
      sdata[tid] += rcorrelation(g_idata1[i], g_idata2[i]) + rcorrelation(g_idata1[i+blockSize], g_idata2[i+blockSize]);
          i += gridSize;
    }
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
       if (blockSize >= 2) sdata[tid] += sdata[tid + 1];
    }
    if (tid == 0) g_odata[blockIdx.x] = sdata[0];
}


/* Define CUDA kernel that squares the input complex array */
__global__ void  bpmask2DComplex(cuFloatComplex *in1,unsigned int Nsize, unsigned int Msize, unsigned int lp, unsigned int hp)
{

  unsigned int n = blockIdx.x * blockDim.x + threadIdx.x;
  unsigned int m = blockIdx.y * blockDim.y + threadIdx.y;

  if(n < Nsize-1 && m < Msize-1){
  int grid_width = gridDim.x * blockDim.x;
  unsigned long index = m + (n * grid_width);

  unsigned int r = n*n + m*m;
  if(r < lp || r > hp) in1[index] = make_cuFloatComplex(0.,0.);
  }
}


/*
   Fortran subroutine arguments are passed by references.
   call fun( array_a, array_b, N) will be mapped to
   function (*a, *b, *N);
*/
// C=corr (A,B)

extern "C" void kernelcrosscorr_(cuFloatComplex *a, cuFloatComplex *b, float *c, int *Ns, int*Np,int* sqlplim, int*sqhplim, int *Bsize, int threads)
{
  unsigned int block_size=(unsigned int)*Bsize;
  cuFloatComplex *a_d,*b_d;
  float *c_d;
  unsigned int N=(unsigned int)*Ns;int M=(unsigned int)*Np;

  int sqlp = (*sqlplim);
  int sqhp = (*sqhplim);

  float sumasq[1] = {0.};
  float sumbsq[1] = {0.};

  cudaSetDevice(0);

  /* Allocate complex array on device */
   cudaMalloc ((void **) &a_d , sizeof(cuFloatComplex)*N*M );
   cudaMalloc ((void **) &b_d , sizeof(cuFloatComplex)*N*M );
   cudaMalloc ((void **) &c_d , sizeof(float)*N*M );
   cudaCheckErrors(" kernelcrosscorr Malloc failed ");

  /* Copy array from host memory to device memory */
  cudaMemcpy( a_d, a,  sizeof(cuFloatComplex)*N*M   ,cudaMemcpyHostToDevice);
  cudaMemcpy( b_d, b,  sizeof(cuFloatComplex)*N*M   ,cudaMemcpyHostToDevice);
  cudaCheckErrors(" kernelcrosscorr Memcpy failed ");

  /* Compute execution configuration */
   dim3 dimBlock(block_size, 8);
   dim3 dimGrid ; //(N/dimBlock.x);
   dimGrid.x = (N + dimBlock.x - 1) / dimBlock.x;
   dimGrid.y = (M + dimBlock.y - 1) / dimBlock.y;
   int smemSize = (threads <= 32) ? 2 * threads * sizeof(float) : threads * sizeof(float);

   if( N*M % block_size != 0 ) dimGrid.x+=1;

   /* Mask A and B */
   bpmask2DComplex<<<dimGrid,dimBlock>>>(a_d, N, M, sqlp, sqhp); cudaCheckErrors("kernelcrosscorr masking A failed ");
   bpmask2DComplex<<<dimGrid,dimBlock>>>(b_d, N, M, sqlp, sqhp); cudaCheckErrors("kernelcrosscorr masking B failed ");

  /* Execute the kernels */

   sumComplexsquared<64><<<dimGrid,dimBlock,smemSize>>>(a_d,c_d,N*M); cudaCheckErrors(" A sumsq failed ");
   // copy device result back to host copy of c
   cudaMemcpy( sumasq, c_d, sizeof( float ) , cudaMemcpyDeviceToHost );

   sumComplexsquared<64><<<dimGrid,dimBlock,smemSize>>>(b_d,c_d,N*M);  cudaCheckErrors(" B sumsq failed ");
   // copy device result back to host copy of c
   cudaMemcpy( sumbsq, c_d, sizeof( float ) , cudaMemcpyDeviceToHost );


   if( sumasq[0] < (float)1e-8 || sumbsq[0] < (float)1e-8 ){
        *c = 0.;

   } else{
     sumrcorr<64><<<dimGrid,dimBlock,smemSize>>>(a_d,b_d,c_d,N*M);  cudaCheckErrors("kernelcrosscorr r sum failed ");
     cudaMemcpy( c, c_d, sizeof( float ) , cudaMemcpyDeviceToHost );

      *c = *c / sqrt(sumasq[0] * sumbsq[0]);
    }

  /* Free memory on the device */
  cudaFree(a_d);
  cudaFree(b_d);
  cudaFree(c_d);
  cudaCheckErrors(" kernelcrosscorr Free failed ");

  return;
}


extern "C" void kernelsumcsq_(cuFloatComplex *a,  float *c, int *Nsize, int *Bsize, int threads)
{
  unsigned int block_size=(unsigned int) *Bsize;
  cuFloatComplex *a_d;
  float *c_d;
  unsigned int N=(unsigned int)*Nsize;

  float sumasq[1] = {0.};
  float sumbsq[1] = {0.};

  cudaSetDevice(0);

  /* Allocate complex array on device */
   cudaMalloc ((void **) &a_d , sizeof(cuFloatComplex)*N );
   cudaMalloc ((void **) &c_d , sizeof(float)*N );
   cudaCheckErrors(" kernelcsq Malloc failed ");

  /* Copy array from host memory to device memory */
  cudaMemcpy( a_d, a,  sizeof(cuFloatComplex)*N   ,cudaMemcpyHostToDevice);
  cudaCheckErrors(" kernelcsq Memcpy failed ");

  /* Compute execution configuration */
   dim3 dimBlock(block_size, 8);
   dim3 dimGrid (N/dimBlock.x);
   //dimGrid.x = (N + dimBlock.x - 1) / dimBlock.x;
   //dimGrid.y = (M + dimBlock.y - 1) / dimBlock.y;
   int smemSize = (threads <= 32) ? 2 * threads * sizeof(float) : threads * sizeof(float);

   if( N % block_size != 0 ) dimGrid.x+=1;

   sumComplexsquared<64><<<dimGrid,dimBlock,smemSize>>>(a_d,c_d,N); cudaCheckErrors(" A sumsq failed ");
   // copy device result back to host copy of c
   cudaMemcpy( c, c_d, sizeof( float ) , cudaMemcpyDeviceToHost ); cudaCheckErrors(" Memcpy return failed ");


  /* Free memory on the device */
  cudaFree(a_d);
  cudaFree(c_d);
  cudaCheckErrors(" kernelcsq Free failed ");

  return;
}



extern "C" void kernelrcorr_(cuFloatComplex *a, cuFloatComplex *b, float *c, int *Nsize, int *Bsize, int threads)
{
  unsigned int block_size=(unsigned int)*Bsize;
  cuFloatComplex *a_d,*b_d;
  float *c_d;
  unsigned int N= (unsigned int)*Nsize;

  cudaSetDevice(0);

  /* Allocate complex array on device */
   cudaMalloc ((void **) &a_d , sizeof(cuFloatComplex)*N );
   cudaMalloc ((void **) &b_d , sizeof(cuFloatComplex)*N );
   cudaMalloc ((void **) &c_d , sizeof(float)*N );
   cudaCheckErrors(" kernelrcorr Malloc failed ");

  /* Copy array from host memory to device memory */
  cudaMemcpy( a_d, a,  sizeof(cuFloatComplex)*N   ,cudaMemcpyHostToDevice);
  cudaMemcpy( b_d, b,  sizeof(cuFloatComplex)*N   ,cudaMemcpyHostToDevice);
  cudaCheckErrors(" kernelrcorr Memcpy failed ");

  /* Compute execution configuration */
   dim3 dimBlock(block_size, 8);
   dim3 dimGrid (N/dimBlock.x);
   //dimGrid.x = (N + dimBlock.x - 1) / dimBlock.x;
   //dimGrid.y = (M + dimBlock.y - 1) / dimBlock.y;
   int smemSize = (threads <= 32) ? 2 * threads * sizeof(float) : threads * sizeof(float);

   if( N % block_size != 0 ) dimGrid.x+=1;


  /* Execute the kernel */
   sumrcorr<64><<<dimGrid,dimBlock,smemSize>>>(a_d,b_d,c_d,N);  cudaCheckErrors("kernelrcorr r sum failed ");
   cudaMemcpy( c, c_d, sizeof( float ) , cudaMemcpyDeviceToHost );

  /* Free memory on the device */
  cudaFree(a_d);
  cudaFree(b_d);
  cudaFree(c_d);
  cudaCheckErrors(" kernelrcorr Free failed ");

  return;
}
