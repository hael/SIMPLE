/*
 *   -- SIMPLE addon
 *      Author: Frederic Bonnet, Date: 27th of January 2016
 *      Monash University
 *      January 2016
 *
 *      Routine which calculates the product (element wise) of 
 *      cuDoubleComplex matrix A and cuDoubleComplex matrix B. and takes the
 *      real part of the product and puts it into a matrix of type double
 *      C = cuCreal(x) * cuCreal(y) + cuCimag(x) * cuCimag(y)
 *      C = Re( A * conjg(B) ) element wise
 *      The product is then summed over the 2nd and 3rd dimensions
 *      to produce a vector of correlators.
 *
 *      Non Special case
 * @precisions normal z -> s d c
*/

#include "commonblas_zz2d.h"
#include "simple_cuDoubleComplex.h"
#include "polarft_gpu.h"
#include "simple.h"

/* Conflicts with #include "common_magma.h" with the min argument*/
/* the cub callers  needs a not -DMAGMA */
#include <util_allocator.cuh>
#include <device/device_reduce.cuh>
using namespace cub;
//---------------------------------------------------------------------
// Globals, constants and typedefs
//---------------------------------------------------------------------
CachingDeviceAllocator gO_allocator(true); //Caching allocator for device memory
//---------------------------------------------------------------------
// Macros for this function kernel
//---------------------------------------------------------------------
#define imin(a,b) (a<b?a:b)
/* variables not used anymore replaced by the data structure debug_gpu_t */
//#define debug true
//#define debug_cpu false
//#define debug_high true
//#define debug_write false
#define debug_write_C false
/* */
#if defined (CUDA) /*preprossing for the CUDA environment */

/* doing the r product Re(A*conjg(B)) */
extern "C" __global__ void
zMXprodOpt_3D_mat( float *C,
               const cuFloatComplex *A,
               const cuFloatComplex *B,
               int d1lim0, int d1lim1, int d2lim0, int d2lim1,
               int npart, int nrot, int nk)
{

  int ipart = blockIdx.x * blockDim.x + threadIdx.x;
  int  irot = blockIdx.y * blockDim.y + threadIdx.y;

  int nhlfrot = nrot / 2;
  int n2rot = 2 * nrot;
  int n2part = 2 * npart;

  int jpart = ipart + d1lim0;
  int  jrot =  irot + d2lim0;
  int    jk = blockIdx.z * blockDim.z + threadIdx.z;

  if ( ipart < npart) {
    if (irot < nhlfrot ){
      if (jk < nk ){
        
        C[(irot    + nhlfrot * jk  ) *  npart + ipart    ] = cuReCCstarmulf(
        A[(irot    + nhlfrot * jk  ) *  npart + ipart    ],
        B[(jrot +   n2rot * jk ) * n2part + jpart] );

      }
    }
  }

  __syncthreads();
}
/* Normalising the cormmat3D with r/sqrt(sqsums_A*sqsums_B) */
extern "C" __global__ void
myONNormOpt_3D(float *cormat3d,
           const float *sqsums_A, const float *sqsums_B,
           float *r, int blocksPerGrid,
           int npart, int nrot, int ipart, int irot,  int d1lim0) {

  int jpart = blockIdx.x * blockDim.x + threadIdx.x;
  int kpart = jpart + d1lim0;
  
  if ( jpart < npart) {
    cormat3d[(ipart+npart*irot) * npart+jpart] = r[jpart] /
      cuDDsqrtf(sqsums_A[jpart],sqsums_B[kpart]);
  }
}
/* summing the 3D matrix treating it a 1D vec */
template <unsigned int blockSize>
__global__ void
myONsumOpt_23D(float *r, float *sumA, float *A, int chunk, int npart, int jpart,
           int blocksPerGrid){

  int i = threadIdx.x + blockIdx.x * blockDim.x;

  int sharedIdx = threadIdx.x;
  __shared__ float shared_A[blockSize];
  int offset = -1;
  
  float temp = 0.0;
  while ( i < chunk) {
    offset++;
    temp += A[i*npart + offset];//A[i];
    i += blockDim.x * gridDim.x;
  }
  shared_A[sharedIdx] = temp;
  __syncthreads();
  int j = blockDim.x / 2;
  while ( j != 0 ) {
    if ( sharedIdx < j ) shared_A[sharedIdx] += shared_A[sharedIdx + j];
    __syncthreads();
    j /= 2;
  }
  if ( sharedIdx == 0 ) sumA[blockIdx.x] = shared_A[0];

  //TODO: need to fix this for the over counting of the sum check junk.cu
  /*
    int ibloc = (i - threadIdx.x ) / blockDim.x;
    r[jpart] = 0.0;

    typedef cub::BlockReduce<float, BLOCKSIZE> BlockReduceT; 

    if ( ibloc < blocksPerGrid) {
    if ( blockIdx.x < blocksPerGrid) r[jpart] += sumA[blockIdx.x];
    }
    for (int s = 0 ; s < blocksPerGrid ; s++) {
    r[jpart] = s; //blocksPerGrid ;  //sumA[blockIdx.x];
    __syncthreads();
    }
  */
  //End TODO

}
/* summing the 3D matrix treating it a 1D vec */
/* using a slow modulo to select entries */
template <unsigned int blockSize>
__global__ void
myONsumOpt_23D_0(float *r, float *sumA, float *A, int chunk, int npart, int jpart,
           int blocksPerGrid){

  __shared__ float shared_A[blockSize];
  unsigned int tid = threadIdx.x;
  unsigned int i = threadIdx.x + blockIdx.x * blockDim.x;

  //shared_A[tid] = A[i*npart+(jpart)];
  //__syncthreads();
  int offset = -1;
  float temp = 0.0 ;
  while (i<chunk) {
    offset++;
    temp += A[i*npart + offset];
    i+= blockDim.x * gridDim.x;
  }
  shared_A[tid] = temp;
  __syncthreads();

  /* do the reduction in shared memory */
  for (unsigned int s=1 ; s< blockDim.x; s*= 2) {
    if (tid%(2*s)==0){shared_A[tid]+=shared_A[tid+s];}  
    __syncthreads();
  }
  if (tid==0) sumA[blockIdx.x]=shared_A[0];

  //r[jpart] = 0.0;
  //while ( i < blocksPerGrid ) {
    //atomicAdd(&r[jpart],shared_A[0]);
  //if (tid==0) atomicAdd(&(r[jpart]),sumA[blockIdx.x]);//shared_A[0]);
  //i++; //=blockDim.x * gridDim.x;;
  //}

}
/* summing the 3D matrix treating it a 1D vec, interleaved addressing */
/* no divergent branch in the for loop replacing the slow modulo */
template <unsigned int blockSize>
__global__ void
myONsumOpt_23D_1(float *r, float *sumA, float *A, int chunk, int npart, int jpart,
           int blocksPerGrid){

  __shared__ float shared_A[blockSize];
  unsigned int tid = threadIdx.x;
  unsigned int i = threadIdx.x + blockIdx.x * blockDim.x;

  int offset = -1;
  float temp = 0.0 ;
  while (i<chunk) {
    offset++;
    temp += A[i*npart + offset];
    i+= blockDim.x * gridDim.x;
  }
  shared_A[tid] = temp;
  __syncthreads();

  /* do the reduction in shared memory */
  for (unsigned int s=1 ; s< blockDim.x; s*= 2) {
    int index = 2*s*tid;
    if ( index < blockDim.x){shared_A[index]+=shared_A[index+s];}
    __syncthreads();
  }
  if (tid==0) sumA[blockIdx.x]=shared_A[0];

}
/* summing the 3D matrix treating it a 1D vec, sequential addressing */
/* which is similar to myONsumOpt_3D kernel above but some idle threads*/
template <unsigned int blockSize>
__global__ void
myONsumOpt_23D_2(float *r, float *sumA, float *A, int chunk, int npart, int jpart,
           int blocksPerGrid){

  __shared__ float shared_A[blockSize];
  unsigned int tid = threadIdx.x;
  unsigned int i = threadIdx.x + blockIdx.x * blockDim.x;

  int offset = -1;
  float temp = 0.0 ;
  while (i<chunk) {
    offset++;
    temp += A[i*npart + offset];
    i+= blockDim.x * gridDim.x;
  }
  shared_A[tid] = temp;
  __syncthreads();

  /* do the reduction in shared memory */
  for (unsigned int s=blockDim.x/2 ; s>0; s>>= 1) { 
    if ( tid < s){shared_A[tid]+=shared_A[tid+s];}
    __syncthreads();
  }
  if (tid==0) sumA[blockIdx.x]=shared_A[0];

}
/* unrolling the wrap as wrappe of 32 */
extern "C"
__device__ void
wrapMyONsum_23D_5(volatile float *shared_A, int tid) {
  shared_A[tid] += shared_A[tid + 32];
  shared_A[tid] += shared_A[tid + 16];
  shared_A[tid] += shared_A[tid +  8];
  shared_A[tid] += shared_A[tid +  4];
  shared_A[tid] += shared_A[tid +  2];
  shared_A[tid] += shared_A[tid +  1];  
}
/* summing the 3D matrix treating it a 1D vec, sequential addressing */
/* which is similar to myONsumOpt_3D kernel above but some idle threads*/
template <unsigned int blockSize>
__global__ void
myONsumOpt_23D_5(float *r, float *sumA, float *A, int chunk, int npart, int jpart,
           int blocksPerGrid){

 __shared__ float shared_A[blockSize];
  unsigned int tid = threadIdx.x;
  unsigned int i = threadIdx.x + blockIdx.x * blockDim.x;

  int offset = -1;
  float temp = 0.0 ;
  while (i<chunk) {
    offset++;
    temp += A[i*npart + offset];
    i+= blockDim.x * gridDim.x;
  }
  shared_A[tid] = temp;
  __syncthreads();

  /* do the reduction in shared memory */
  for (unsigned int s=blockDim.x/2 ; s>32; s>>= 1) { 
    if ( tid < s){shared_A[tid]+=shared_A[tid+s];}
    __syncthreads();
  }

  if (tid < 32 ) wrapMyONsum_23D_5(shared_A,tid);
  if (tid==0) sumA[blockIdx.x]=shared_A[0];

}
/* summing the 3D matrix treating it a 1D vec, competely unrolled */
/* which is similar to myONsumOpt_3D kernel but unrolled              */
template <unsigned int blockSize>
__device__ void
wrapMyONsum_23D_6(volatile float *shared_A, int tid) {
  if (blockSize >= 64 ) shared_A[tid] += shared_A[tid + 32];
  if (blockSize >= 32 ) shared_A[tid] += shared_A[tid + 16];
  if (blockSize >= 16 ) shared_A[tid] += shared_A[tid +  8];
  if (blockSize >=  8 ) shared_A[tid] += shared_A[tid +  4];
  if (blockSize >=  4 ) shared_A[tid] += shared_A[tid +  2];
  if (blockSize >=  2 ) shared_A[tid] += shared_A[tid +  1];  
}
/* summing the 3D matrix treating it a 1D vec, sequential addressing */
/* which is similar to myONsumOpt_3D kernel above but some idle threads*/
template <unsigned int blockSize>
__global__ void
myONsumOpt_23D_6(float *r, float *sumA, float *A, int chunk, int npart, int jpart,
        int blocksPerGrid){

  __shared__ float shared_A[blockSize];
  unsigned int tid = threadIdx.x;
  unsigned int i = threadIdx.x + blockIdx.x * blockDim.x;

  int offset = -1;
  float temp = 0.0 ;
  while (i<chunk) {
    offset++;
    temp += A[i*npart + offset];
    i+= blockDim.x * gridDim.x;
  }
  shared_A[tid] = temp;
  __syncthreads();

  /* do the reduction in shared memory */
  if(blockSize>=512){
    if(tid<256){shared_A[tid]+=shared_A[tid+256];}__syncthreads();}
  if(blockSize>=256){
    if(tid<128){shared_A[tid]+=shared_A[tid+128];}__syncthreads();}
  if(blockSize>=128){
    if(tid<64) {shared_A[tid]+=shared_A[tid+ 64];}__syncthreads();}
  
  if (tid < 32 ) wrapMyONsum_23D_6<blockSize>(shared_A,tid);
  if (tid==0) sumA[blockIdx.x]=shared_A[0];

}
/* getter function to get the kernel call for sum 2D for different kernels */
template <typename S>
int get_myONsumOpt_23D(S iKernel,
                   float *d_r, float *d_partial_1D, float *startptr,
                   int chunk, int npart, int jpart,
                   int blocksPerGrid, int threadsPerBlock) {
  int rc = RC_SUCCESS;
  
  switch (iKernel) {

  case 0:
    switch(threadsPerBlock){
    case 512: myONsumOpt_23D_0<512><<<blocksPerGrid,threadsPerBlock>>>
        (d_r,d_partial_1D,startptr, chunk,npart,jpart,blocksPerGrid); break;
    case 256:
      myONsumOpt_23D_0<256><<<blocksPerGrid,threadsPerBlock>>>
        (d_r,d_partial_1D,startptr, chunk,npart,jpart,blocksPerGrid); break;
    case 128:
      myONsumOpt_23D_0<128><<<blocksPerGrid,threadsPerBlock>>>
        (d_r,d_partial_1D,startptr, chunk,npart,jpart,blocksPerGrid); break;
    case 64:
      myONsumOpt_23D_0<64><<<blocksPerGrid,threadsPerBlock>>>
        (d_r,d_partial_1D,startptr, chunk,npart,jpart,blocksPerGrid); break;
    case 32:
      myONsumOpt_23D_0<32><<<blocksPerGrid,threadsPerBlock>>>
        (d_r,d_partial_1D,startptr, chunk,npart,jpart,blocksPerGrid); break;
    case 16:
      myONsumOpt_23D_0<16><<<blocksPerGrid,threadsPerBlock>>>
        (d_r,d_partial_1D,startptr, chunk,npart,jpart,blocksPerGrid); break;
    case 8:
      myONsumOpt_23D_0<8><<<blocksPerGrid,threadsPerBlock>>>
        (d_r,d_partial_1D,startptr, chunk,npart,jpart,blocksPerGrid); break;
    case 4:
      myONsumOpt_23D_0<4><<<blocksPerGrid,threadsPerBlock>>>
        (d_r,d_partial_1D,startptr, chunk,npart,jpart,blocksPerGrid); break;
    case 2:
      myONsumOpt_23D_0<2><<<blocksPerGrid,threadsPerBlock>>>
        (d_r,d_partial_1D,startptr, chunk,npart,jpart,blocksPerGrid); break;
    case 1:
      myONsumOpt_23D_0<1><<<blocksPerGrid,threadsPerBlock>>>
        (d_r,d_partial_1D,startptr, chunk,npart,jpart,blocksPerGrid); break;
    }
    break;
  case 1:
    switch(threadsPerBlock){
    case 512:
      myONsumOpt_23D_1<512><<<blocksPerGrid,threadsPerBlock>>>
        (d_r,d_partial_1D,startptr, chunk,npart,jpart,blocksPerGrid); break;
    case 256:
      myONsumOpt_23D_1<256><<<blocksPerGrid,threadsPerBlock>>>
        (d_r,d_partial_1D,startptr, chunk,npart,jpart,blocksPerGrid); break;
    case 128:
      myONsumOpt_23D_1<128><<<blocksPerGrid,threadsPerBlock>>>
        (d_r,d_partial_1D,startptr, chunk,npart,jpart,blocksPerGrid); break;
    case 64:
      myONsumOpt_23D_1<64><<<blocksPerGrid,threadsPerBlock>>>
        (d_r,d_partial_1D,startptr, chunk,npart,jpart,blocksPerGrid); break;
    case 32:
      myONsumOpt_23D_1<32><<<blocksPerGrid,threadsPerBlock>>>
        (d_r,d_partial_1D,startptr, chunk,npart,jpart,blocksPerGrid); break;
    case 16:
      myONsumOpt_23D_1<16><<<blocksPerGrid,threadsPerBlock>>>
        (d_r,d_partial_1D,startptr, chunk,npart,jpart,blocksPerGrid); break;
    case 8:
      myONsumOpt_23D_1<8><<<blocksPerGrid,threadsPerBlock>>>
        (d_r,d_partial_1D,startptr, chunk,npart,jpart,blocksPerGrid); break;
    case 4:
      myONsumOpt_23D_1<4><<<blocksPerGrid,threadsPerBlock>>>
        (d_r,d_partial_1D,startptr, chunk,npart,jpart,blocksPerGrid); break;
    case 2:
      myONsumOpt_23D_1<2><<<blocksPerGrid,threadsPerBlock>>>
        (d_r,d_partial_1D,startptr, chunk,npart,jpart,blocksPerGrid); break;
    case 1:
      myONsumOpt_23D_1<1><<<blocksPerGrid,threadsPerBlock>>>
        (d_r,d_partial_1D,startptr, chunk,npart,jpart,blocksPerGrid); break;
    }
    break;
  case 2:
    switch(threadsPerBlock){
    case 512:
      myONsumOpt_23D_2<512><<<blocksPerGrid,threadsPerBlock>>>
        (d_r,d_partial_1D,startptr, chunk,npart,jpart,blocksPerGrid); break;
    case 256:
      myONsumOpt_23D_2<256><<<blocksPerGrid,threadsPerBlock>>>
        (d_r,d_partial_1D,startptr, chunk,npart,jpart,blocksPerGrid); break;
    case 128:
      myONsumOpt_23D_2<128><<<blocksPerGrid,threadsPerBlock>>>
        (d_r,d_partial_1D,startptr, chunk,npart,jpart,blocksPerGrid); break;
    case 64:
      myONsumOpt_23D_2<64><<<blocksPerGrid,threadsPerBlock>>>
        (d_r,d_partial_1D,startptr, chunk,npart,jpart,blocksPerGrid); break;
    case 32:
      myONsumOpt_23D_2<32><<<blocksPerGrid,threadsPerBlock>>>
        (d_r,d_partial_1D,startptr, chunk,npart,jpart,blocksPerGrid); break;
    case 16:
      myONsumOpt_23D_2<16><<<blocksPerGrid,threadsPerBlock>>>
        (d_r,d_partial_1D,startptr, chunk,npart,jpart,blocksPerGrid); break;
    case 8:
      myONsumOpt_23D_2<8><<<blocksPerGrid,threadsPerBlock>>>
        (d_r,d_partial_1D,startptr, chunk,npart,jpart,blocksPerGrid); break;
    case 4:
      myONsumOpt_23D_2<4><<<blocksPerGrid,threadsPerBlock>>>
        (d_r,d_partial_1D,startptr, chunk,npart,jpart,blocksPerGrid); break;
    case 2:
      myONsumOpt_23D_2<2><<<blocksPerGrid,threadsPerBlock>>>
        (d_r,d_partial_1D,startptr, chunk,npart,jpart,blocksPerGrid); break;
    case 1:
      myONsumOpt_23D_2<1><<<blocksPerGrid,threadsPerBlock>>>
        (d_r,d_partial_1D,startptr, chunk,npart,jpart,blocksPerGrid); break;
    }
    break;
  case 5:
    switch(threadsPerBlock){
    case 512:
      myONsumOpt_23D_5<512><<<blocksPerGrid,threadsPerBlock>>>
        (d_r,d_partial_1D,startptr, chunk,npart,jpart,blocksPerGrid); break;
    case 256:
      myONsumOpt_23D_5<256><<<blocksPerGrid,threadsPerBlock>>>
        (d_r,d_partial_1D,startptr, chunk,npart,jpart,blocksPerGrid); break;
    case 128:
      myONsumOpt_23D_5<128><<<blocksPerGrid,threadsPerBlock>>>
        (d_r,d_partial_1D,startptr, chunk,npart,jpart,blocksPerGrid); break;
    case 64:
      myONsumOpt_23D_5<64><<<blocksPerGrid,threadsPerBlock>>>
        (d_r,d_partial_1D,startptr, chunk,npart,jpart,blocksPerGrid); break;
    case 32:
      myONsumOpt_23D_5<32><<<blocksPerGrid,threadsPerBlock>>>
        (d_r,d_partial_1D,startptr, chunk,npart,jpart,blocksPerGrid); break;
    case 16:
      myONsumOpt_23D_5<16><<<blocksPerGrid,threadsPerBlock>>>
        (d_r,d_partial_1D,startptr, chunk,npart,jpart,blocksPerGrid); break;
    case 8:
      myONsumOpt_23D_5<8><<<blocksPerGrid,threadsPerBlock>>>
        (d_r,d_partial_1D,startptr, chunk,npart,jpart,blocksPerGrid); break;
    case 4:
      myONsumOpt_23D_5<4><<<blocksPerGrid,threadsPerBlock>>>
        (d_r,d_partial_1D,startptr, chunk,npart,jpart,blocksPerGrid); break;
    case 2:
      myONsumOpt_23D_5<2><<<blocksPerGrid,threadsPerBlock>>>
        (d_r,d_partial_1D,startptr, chunk,npart,jpart,blocksPerGrid); break;
    case 1:
      myONsumOpt_23D_5<1><<<blocksPerGrid,threadsPerBlock>>>
        (d_r,d_partial_1D,startptr, chunk,npart,jpart,blocksPerGrid); break;
    }
    break;
  case 6:
    switch(threadsPerBlock){
    case 512:
      myONsumOpt_23D_6<512><<<blocksPerGrid,threadsPerBlock>>>
        (d_r,d_partial_1D,startptr,chunk,npart,jpart,blocksPerGrid); break;
    case 256:
      myONsumOpt_23D_6<256><<<blocksPerGrid,threadsPerBlock>>>
        (d_r,d_partial_1D,startptr,chunk,npart,jpart,blocksPerGrid); break;
    case 128:
      myONsumOpt_23D_6<128><<<blocksPerGrid,threadsPerBlock>>>
        (d_r,d_partial_1D,startptr,chunk,npart,jpart,blocksPerGrid); break;
    case 64:
      myONsumOpt_23D_6<64><<<blocksPerGrid,threadsPerBlock>>>
        (d_r,d_partial_1D,startptr,chunk,npart,jpart,blocksPerGrid); break;
    case 32:
      myONsumOpt_23D_6<32><<<blocksPerGrid,threadsPerBlock>>>
        (d_r,d_partial_1D,startptr,chunk,npart,jpart,blocksPerGrid); break;
    case 16:
      myONsumOpt_23D_6<16><<<blocksPerGrid,threadsPerBlock>>>
        (d_r,d_partial_1D,startptr,chunk,npart,jpart,blocksPerGrid); break;
    case 8:
      myONsumOpt_23D_6<8><<<blocksPerGrid,threadsPerBlock>>>
        (d_r,d_partial_1D,startptr,chunk,npart,jpart,blocksPerGrid); break;
    case 4:
      myONsumOpt_23D_6<4><<<blocksPerGrid,threadsPerBlock>>>
        (d_r,d_partial_1D,startptr,chunk,npart,jpart,blocksPerGrid); break;
    case 2:
      myONsumOpt_23D_6<2><<<blocksPerGrid,threadsPerBlock>>>
        (d_r,d_partial_1D,startptr,chunk,npart,jpart,blocksPerGrid); break;
    case 1:
      myONsumOpt_23D_6<1><<<blocksPerGrid,threadsPerBlock>>>
        (d_r,d_partial_1D,startptr,chunk,npart,jpart,blocksPerGrid); break;
    }
    break;
  default:
    switch(threadsPerBlock){
    case 512:
      myONsumOpt_23D<512><<<blocksPerGrid,threadsPerBlock>>>
        (d_r,d_partial_1D,startptr,chunk,npart,jpart,blocksPerGrid); break;
    case 256:
      myONsumOpt_23D<256><<<blocksPerGrid,threadsPerBlock>>>
        (d_r,d_partial_1D,startptr,chunk,npart,jpart,blocksPerGrid); break;
    case 128:
      myONsumOpt_23D<128><<<blocksPerGrid,threadsPerBlock>>>
        (d_r,d_partial_1D,startptr,chunk,npart,jpart,blocksPerGrid); break;
    case 64:
      myONsumOpt_23D<64><<<blocksPerGrid,threadsPerBlock>>>
        (d_r,d_partial_1D,startptr,chunk,npart,jpart,blocksPerGrid); break;
    case 32:
      myONsumOpt_23D<32><<<blocksPerGrid,threadsPerBlock>>>
        (d_r,d_partial_1D,startptr,chunk,npart,jpart,blocksPerGrid); break;
    case 16:
      myONsumOpt_23D<16><<<blocksPerGrid,threadsPerBlock>>>
        (d_r,d_partial_1D,startptr,chunk,npart,jpart,blocksPerGrid); break;
    case 8:
      myONsumOpt_23D<8><<<blocksPerGrid,threadsPerBlock>>>
        (d_r,d_partial_1D,startptr,chunk,npart,jpart,blocksPerGrid); break;
    case 4:
      myONsumOpt_23D<4><<<blocksPerGrid,threadsPerBlock>>>
        (d_r,d_partial_1D,startptr,chunk,npart,jpart,blocksPerGrid); break;
    case 2:
      myONsumOpt_23D<2><<<blocksPerGrid,threadsPerBlock>>>
        (d_r,d_partial_1D,startptr,chunk,npart,jpart,blocksPerGrid); break;
    case 1:
      myONsumOpt_23D<1><<<blocksPerGrid,threadsPerBlock>>>
        (d_r,d_partial_1D,startptr,chunk,npart,jpart,blocksPerGrid); break;
    }
    break;
  }

  return rc;
}

/* main kernel entry */
extern "C" int
polarft_krnl_Opti_O_N(deviceDetails_t * s_devD,
                      polar_corr_calc_t *s_polar,
                      float *r, float *cormat3d,
                      const cuFloatComplex *A,
                      const cuFloatComplex *B,
                      const float *sqsums_A,
                      const float *sqsums_B,
                      int npart, int nrot, int nk,
                      float alpha,
                      bench_t *s_bench, debug_gpu_t *s_debug_gpu)
{
  /* return code */
  int rc = RC_SUCCESS; 
  /* array dimension */
  int d1lim[2], d2lim[2];
  /* setting the values into the object */
  pfts_Sizes *p_pfts_Sizes = new pfts_Sizes(npart,nrot,nk);
  mesh_3D *p_mesh_3D = new mesh_3D(s_polar, p_pfts_Sizes->get_npart(),
                                   p_pfts_Sizes->get_nhlfrot(),
                                   p_pfts_Sizes->get_nk());
  mesh_1D *p_mesh_1D = new mesh_1D(s_polar, p_pfts_Sizes->get_npart(),
                                   p_pfts_Sizes->get_nhlfrot(),
                                   p_pfts_Sizes->get_nk());
  /* grid and threads block definition */
  int nx, ny, nz;               //Dimension of the threads
  int gridx, gridy, gridz;      //Mesh of the 3D Grid
  /* device allocation */
  cuFloatComplex *d_A = NULL;   //device pointer for A (const in)
  cuFloatComplex *d_B = NULL;   //device pointer for B (const in)
  float          *d_C = NULL;   //device pointer for C (re[A*Bstar])
  float          *d_r = NULL;   //device pointer for r (float in)
  float   *d_sqsums_A = NULL;
  float   *d_sqsums_B = NULL;
  float   *d_cormat3d = NULL;   //device pointer corrmat3D(ipart,ipart,nrot) 
  /* creating the working space on the GPU */
  float *d_partial_1D = NULL;
  /* size of the element in consideration */
  int    size_dA =    p_pfts_Sizes->get_npts_A() * sizeof(cuFloatComplex);
  int    size_dB =    p_pfts_Sizes->get_npts_B() * sizeof(cuFloatComplex);
  int    size_dC =    p_pfts_Sizes->get_npts_A() * sizeof(float);
  int    size_dr =     p_pfts_Sizes->get_npart() * sizeof(float);
  int  size_dsqA =     p_pfts_Sizes->get_npart() * sizeof(float);
  int  size_dsqB =    p_pfts_Sizes->get_n2part() * sizeof(float);
  int size_dcr3D = p_pfts_Sizes->get_npts_cr3D() * sizeof(float);

  /* indexer */
  int i, ipart, jpart, irot;
  /* the error handlers from the cuda library */
  cudaError_t err;
  /* printer index depth */
  float depth = 2;
  
  /*start of the execution commands */

  /* setting up the mesh using the s_polar structure */
  nx=s_polar->nx; ny=s_polar->ny; nz=s_polar->nz;
  dim3 threads(nx,ny,nz); /*mesh dimensions for the hadamard product*/
  gridx = p_mesh_3D->get_mesh3D_gridx(); gridy = p_mesh_3D->get_mesh3D_gridy();
  gridz = p_mesh_3D->get_mesh3D_gridz(); /*checking the mesh*/
  rc = check_grid3D(gridx,gridy,gridz,s_polar,p_pfts_Sizes);
  //if (rc != 0 ) {rc = 
  dim3 grid(gridx, gridy, gridz);
  if (s_debug_gpu->debug_i == true ) {
    rc = print_3D_mesh(1,p_mesh_3D,s_polar,p_pfts_Sizes,nx,ny,nz);//gridx, gridy, gridz);
    /* printing input variables on screen */
    if ( s_debug_gpu->debug_high_i == true ) {
      rc = print_s_debug_struct(s_debug_gpu);
      rc = print_s_bench_struct(s_bench);
      rc = print_function_header_Z_N_info(s_polar,
                                          r, cormat3d,
                                          A, B, 
                                          sqsums_A, sqsums_B,
                                          p_pfts_Sizes,
                                          npart, nrot, nk,
                                          alpha, depth);
    }
  }
  /* first on cpu */
  if (s_debug_gpu->debug_cpu_i == true ) {
    rc = polarft_gencorrAll_Z_N_cpu(s_polar, r, cormat3d,
                                    A, B, sqsums_A, sqsums_B,
                                    p_pfts_Sizes, npart, nrot, nk,
                                    alpha);

  }
  
  /* allocating the memory on GPU for d_{A,B,C}, d_sq{A,B} and d_cormat3d */
  err = cudaMalloc((void**)&d_A, size_dA);
  err = cudaMalloc((void**)&d_B, size_dB);
  err = cudaMalloc((void**)&d_C, size_dA);
  err = cudaMalloc((void**)&d_r, size_dr);
  err = cudaMalloc((void**)&d_sqsums_A, size_dsqA);
  err = cudaMalloc((void**)&d_sqsums_B, size_dsqB);
  err = cudaMalloc((void**)&d_cormat3d, size_dcr3D);
  if ( (int)err != CUDA_SUCCESS ) {rc = get_error_id_corr_Hadmr_gpu(err); }
  //uploading the matrices A, B to GPU
  err = cudaMemcpy(d_A, A, size_dA, cudaMemcpyHostToDevice);
  err = cudaMemcpy(d_B, B, size_dB, cudaMemcpyHostToDevice);
  //uploading the vectors sqsums A, B to GPU
  err = cudaMemcpy(d_sqsums_A, sqsums_A, size_dsqA, cudaMemcpyHostToDevice);
  err = cudaMemcpy(d_sqsums_B, sqsums_B, size_dsqB, cudaMemcpyHostToDevice);

  //Preparing mesh for summing the over the chunk and mapping it to vector r  
  //imin(32,(chunk+threadsPerBlock-1)/threadsPerBlock);
  int chunk = p_mesh_1D->get_mesh1D_chunk();//p_pfts_Sizes->get_nhlfrot() * p_pfts_Sizes->get_nk();
  int threadsPerBlock = p_mesh_1D->get_mesh1D_threadsPerBlock();//s_polar->threadsPerBlock; //32;
  int blocksPerGrid = p_mesh_1D->get_mesh1D_blocksPerGrid();//(chunk/threadsPerBlock)+(chunk%threadsPerBlock!=0);
  int size_pB = blocksPerGrid * sizeof(float);
  /* allocating the size of the partial sum on device */
  err = cudaMalloc(&d_partial_1D,size_pB);

  if ( s_debug_gpu->debug_i == true ) {
    rc = print_summing_method(1);
    rc = print_1D_mesh(1, p_mesh_1D, chunk, threadsPerBlock, blocksPerGrid);
    rc = print_iKernel_threadsPerBlock(s_polar, chunk,
                                       threadsPerBlock, blocksPerGrid);
  }

  // Request and allocate temporary storage for optimised size blocks
  void  *d_temp_storage = NULL;
  size_t temp_storage_bytes = 0;
  err = DeviceReduce::Sum(d_temp_storage, temp_storage_bytes,
                          d_partial_1D, d_r, blocksPerGrid);
  gO_allocator.DeviceAllocate(&d_temp_storage, temp_storage_bytes);

  /* the pft cross correlator calculation */
  /* initilizing the dim arrays */
  for (i = 0 ; i < 2 ; i++) {d1lim[i] = 0;} 
  for (i = 0 ; i < 2 ; i++) {d2lim[i] = 0;}
  /* start of the loop ipart irot */
  for (ipart = 0 ; ipart < p_pfts_Sizes->get_npart() ; ipart++) {
    d1lim[0] = ipart;
    d1lim[1] = ipart + p_pfts_Sizes->get_npart();
    for (irot = 0 ; irot < p_pfts_Sizes->get_nrot() ; irot++) {
      d2lim[0] = irot;
      d2lim[1] = irot + p_pfts_Sizes->get_nhlfrot();

      zMXprodOpt_3D_mat<<<grid,threads>>>(d_C, d_A, d_B,
                                          d1lim[0],d1lim[1],
                                          d2lim[0],d2lim[1],
                                          p_pfts_Sizes->get_npart(),
                                          p_pfts_Sizes->get_nrot(),
                                          p_pfts_Sizes->get_nk()     );

      //now summing the over the chunk and mapping it to vector r

      //Start loop here
      for (jpart = 0 ; jpart < p_pfts_Sizes->get_npart() ; jpart++ ) {
        float *startptr = d_C + jpart;

        //computing on sum of the last (jpart,:,:) on GPU //<threadsPerBlock>

        rc = get_myONsumOpt_23D<int>(s_polar->ikrnl, d_r,
                                     d_partial_1D, startptr,
                                     chunk, npart, jpart,
                                     blocksPerGrid,threadsPerBlock);

        err = DeviceReduce::Sum(d_temp_storage, temp_storage_bytes,
                                d_partial_1D, &d_r[jpart], blocksPerGrid);
      }//finish the loop here.

      myONNormOpt_3D<<<gridx,nx>>>(d_cormat3d, d_sqsums_A, d_sqsums_B, d_r,
                                   blocksPerGrid,
                                   p_pfts_Sizes->get_npart(),
                                   p_pfts_Sizes->get_nrot(),
                                   ipart, irot, d1lim[0]);
    }
  } /* end of the loop ipart irot */

  /* retrieving the corrmat 3D from device */
  err = cudaMemcpy(cormat3d,d_cormat3d,size_dcr3D,cudaMemcpyDeviceToHost);

  /* debugging for the hadamard product from device */
  if (s_debug_gpu->debug_write_C_i == true) {
    rc = print_dC_HadrProd(d_C,p_pfts_Sizes);}
  /* sanity check on the d{1,2}lim arrays */
  if ( s_debug_gpu->debug_i == true ) {rc = print_d12lim(d1lim, d2lim);}

  /* freeing the resources on the CPU */
  /* freeing the resources on the GPU */
  err = cudaFree(d_A);
  err = cudaFree(d_B);
  err = cudaFree(d_C);
  err = cudaFree(d_r);
  err = cudaFree(d_sqsums_A);
  err = cudaFree(d_sqsums_B);
  err = cudaFree(d_cormat3d);
  if ( (int)err != CUDA_SUCCESS ) {rc = get_error_id_corr_Hadmr_gpu(err); }
  /* calling the destructor for the pfts_Sizes object */
  p_pfts_Sizes->~pfts_Sizes();
  p_mesh_3D->~mesh_3D();
  p_mesh_1D->~mesh_1D();

  return rc;
 
} /* polarft_multi_GPUs_O_N */

#endif /* CUDA */
