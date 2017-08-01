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
CachingDeviceAllocator gC_allocator(true); //Caching allocator for device memory
//---------------------------------------------------------------------
// Macros for this function kernel
//---------------------------------------------------------------------
#define imin(a,b) (a<b?a:b)
/* variables not used anymore replaced by the data structure debug_gpu_t */
//#define debug true
//#define debug_cpu false
//#define debug_high true
//#define debug_write false
//#define debug_write_C false
/* */
#if defined (CUDA) /*preprossing for the CUDA environment */
//---------------------------------------------------------------------
// Hadamard products for the standard product Re(A*B)
//---------------------------------------------------------------------
/* doing the product Re(A*B) */
extern "C" __global__ void
NXprod_3D_mat(cuFloatComplex *C,
              const cuFloatComplex *A,
              const cuFloatComplex *B,
              int npart, int nrot, int nk)
{

  int ipart = blockIdx.x * blockDim.x + threadIdx.x;
  int  irot = blockIdx.y * blockDim.y + threadIdx.y;
  int    jk = blockIdx.z * blockDim.z + threadIdx.z;

  if ( ipart < npart) {
    if (irot < nrot ){
      if (jk < nk ){
        /* TODO include in the simple_cuDoubleComplex.h atomic function*/
        C[(irot + nrot * jk ) * npart + ipart  ] = cuCCmulf(
        A[(irot + nrot * jk ) * npart + ipart  ] ,
        B[(irot + nrot * jk ) * npart + ipart]   );
      }
    }
  }
  __syncthreads();
}
//---------------------------------------------------------------------
// Hadamard products for the standard product Re(A*conjg(B))
//---------------------------------------------------------------------
/* doing the r product Re(A*conjg(B)) */
extern "C" __global__ void
cNXprod_3D_mat( float *C,
               const cuFloatComplex *A,
               const cuFloatComplex *B,
               int npart, int nrot, int nk)
{

  int ipart = blockIdx.x * blockDim.x + threadIdx.x;
  int  irot = blockIdx.y * blockDim.y + threadIdx.y;
  int    jk = blockIdx.z * blockDim.z + threadIdx.z;

  if ( ipart < npart) {
    if (irot < nrot ){
      if (jk < nk ){
        
        C[(irot + nrot * jk ) * npart + ipart  ] = cuReCCstarmulf(
        A[(irot + nrot * jk ) * npart + ipart  ],
        B[(irot + nrot * jk ) * npart + ipart] );

      }
    }
  }
  __syncthreads();
}
//---------------------------------------------------------------------
// The common device only volatile kernels to the summing methods
//---------------------------------------------------------------------
/* unrolling the wrap as wrappe of 32 */
extern "C"
__device__ void
wrapMyCNsum_r_5(volatile float *shared_A, int tid) {
  shared_A[tid] += shared_A[tid + 32];
  shared_A[tid] += shared_A[tid + 16];
  shared_A[tid] += shared_A[tid +  8];
  shared_A[tid] += shared_A[tid +  4];
  shared_A[tid] += shared_A[tid +  2];
  shared_A[tid] += shared_A[tid +  1];  
}
/* summing the 3D matrix treating it a 1D vec, competely unrolled */
/* which is similar to myCNsum_3D kernel but unrolled              */
template <unsigned int blockSize>
__device__ void
wrapMyCNsum_r_6(volatile float *shared_A, int tid) {
  if (blockSize >= 64 ) shared_A[tid] += shared_A[tid + 32];
  if (blockSize >= 32 ) shared_A[tid] += shared_A[tid + 16];
  if (blockSize >= 16 ) shared_A[tid] += shared_A[tid +  8];
  if (blockSize >=  8 ) shared_A[tid] += shared_A[tid +  4];
  if (blockSize >=  4 ) shared_A[tid] += shared_A[tid +  2];
  if (blockSize >=  2 ) shared_A[tid] += shared_A[tid +  1];  
}
//---------------------------------------------------------------------
// The r sum including the atomic product Re(A*conjg(B))
//---------------------------------------------------------------------
/* summing the 3D matrix treating it a 1D vec */
template <unsigned int blockSize>
__global__ void
myCNsum_r(const cuFloatComplex *A, float *sumA, const cuFloatComplex *B, int chunk, int npart, int jpart,
          int blocksPerGrid){

  int i = threadIdx.x + blockIdx.x * blockDim.x;

  int sharedIdx = threadIdx.x;
  __shared__ float shared_A[blockSize];
  int offset = -1;
  
  float temp = 0.0;
  while ( i < chunk) {
    offset++;
    temp += cuReCCstarmulf(A[i],B[i]);//A[i];//A[i*npart + offset];
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

}
/* summing the 3D matrix treating it a 1D vec */
/* using a slow modulo to select entries */
template <unsigned int blockSize>
__global__ void
myCNsum_r_0(const cuFloatComplex *A, float *sumA, const cuFloatComplex *B, int chunk, int npart, int jpart,
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
    temp += cuReCCstarmulf(A[i],B[i]);//A[i];//A[i*npart + offset];
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

}
/* summing the 3D matrix treating it a 1D vec, interleaved addressing */
/* no divergent branch in the for loop replacing the slow modulo */
template <unsigned int blockSize>
__global__ void
myCNsum_r_1(const cuFloatComplex *A, float *sumA, const cuFloatComplex *B, int chunk, int npart, int jpart,
            int blocksPerGrid){

  __shared__ float shared_A[blockSize];
  unsigned int tid = threadIdx.x;
  unsigned int i = threadIdx.x + blockIdx.x * blockDim.x;

  int offset = -1;
  float temp = 0.0 ;
  while (i<chunk) {
    offset++;
    temp += cuReCCstarmulf(A[i],B[i]);//A[i];//A[i*npart + offset];
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
/* which is similar to myCNsum_3D kernel above but some idle threads*/
template <unsigned int blockSize>
__global__ void
myCNsum_r_2(const cuFloatComplex *A, float *sumA, const cuFloatComplex *B, int chunk, int npart, int jpart,
            int blocksPerGrid){

  __shared__ float shared_A[blockSize];
  unsigned int tid = threadIdx.x;
  unsigned int i = threadIdx.x + blockIdx.x * blockDim.x;

  int offset = -1;
  float temp = 0.0 ;
  while (i<chunk) {
    offset++;
    temp += cuReCCstarmulf(A[i],B[i]);//A[i];//A[i*npart + offset];
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
/* summing the 3D matrix treating it a 1D vec, sequential addressing */
/* which is similar to myCNsum_3D kernel above but some idle threads*/
template <unsigned int blockSize>
__global__ void
myCNsum_r_5(const cuFloatComplex *A, float *sumA, const cuFloatComplex *B, int chunk, int npart, int jpart,
            int blocksPerGrid){

 __shared__ float shared_A[blockSize];
  unsigned int tid = threadIdx.x;
  unsigned int i = threadIdx.x + blockIdx.x * blockDim.x;

  int offset = -1;
  float temp = 0.0 ;
  while (i<chunk) {
    offset++;
    temp += cuReCCstarmulf(A[i],B[i]);//A[i];//A[i*npart + offset];
    i+= blockDim.x * gridDim.x;
  }
  shared_A[tid] = temp;
  __syncthreads();

  /* do the reduction in shared memory */
  for (unsigned int s=blockDim.x/2 ; s>32; s>>= 1) { 
    if ( tid < s){shared_A[tid]+=shared_A[tid+s];}
    __syncthreads();
  }

  if (tid < 32 ) wrapMyCNsum_r_5(shared_A,tid);
  if (tid==0) sumA[blockIdx.x]=shared_A[0];

}
/* summing the 3D matrix treating it a 1D vec, sequential addressing */
/* which is similar to myCNsum_3D kernel above but some idle threads*/
template <unsigned int blockSize>
__global__ void
myCNsum_r_6(const cuFloatComplex *A, float *sumA, const cuFloatComplex *B, int chunk, int npart, int jpart,
            int blocksPerGrid){

  __shared__ float shared_A[blockSize];
  unsigned int tid = threadIdx.x;
  unsigned int i = threadIdx.x + blockIdx.x * blockDim.x;

  int offset = -1;
  float temp = 0.0 ;
  while (i<chunk) {
    offset++;
    temp += cuReCCstarmulf(A[i],B[i]);//A[i];//A[i*npart + offset];
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
  
  if (tid < 32 ) wrapMyCNsum_r_6<blockSize>(shared_A,tid);
  if (tid==0) sumA[blockIdx.x]=shared_A[0];

}
/* getter function to get the kernel call for sum 2D for different kernels */
/*float *d_r, float *d_partial_1D, float *startptr,*/
template <typename S>
int get_myCNsum_r(S iKernel,
                  cuFloatComplex *A, float *sumA, cuFloatComplex *B,
                  int chunk, int npart, int jpart,
                  int blocksPerGrid, int threadsPerBlock) {
  int rc = RC_SUCCESS;

  switch (iKernel) {

  case 0:
    switch(threadsPerBlock){
    case 512: myCNsum_r_0<512><<<blocksPerGrid,threadsPerBlock>>>
        (A,sumA,B,chunk,npart,jpart,blocksPerGrid); break;
    case 256:
      myCNsum_r_0<256><<<blocksPerGrid,threadsPerBlock>>>
        (A,sumA,B,chunk,npart,jpart,blocksPerGrid); break;
    case 128:
      myCNsum_r_0<128><<<blocksPerGrid,threadsPerBlock>>>
        (A,sumA,B,chunk,npart,jpart,blocksPerGrid); break;
    case 64:
      myCNsum_r_0<64><<<blocksPerGrid,threadsPerBlock>>>
        (A,sumA,B,chunk,npart,jpart,blocksPerGrid); break;
    case 32:
      myCNsum_r_0<32><<<blocksPerGrid,threadsPerBlock>>>
        (A,sumA,B,chunk,npart,jpart,blocksPerGrid); break;
    case 16:
      myCNsum_r_0<16><<<blocksPerGrid,threadsPerBlock>>>
        (A,sumA,B,chunk,npart,jpart,blocksPerGrid); break;
    case 8:
      myCNsum_r_0<8><<<blocksPerGrid,threadsPerBlock>>>
        (A,sumA,B,chunk,npart,jpart,blocksPerGrid); break;
    case 4:
      myCNsum_r_0<4><<<blocksPerGrid,threadsPerBlock>>>
        (A,sumA,B,chunk,npart,jpart,blocksPerGrid); break;
    case 2:
      myCNsum_r_0<2><<<blocksPerGrid,threadsPerBlock>>>
        (A,sumA,B,chunk,npart,jpart,blocksPerGrid); break;
    case 1:
      myCNsum_r_0<1><<<blocksPerGrid,threadsPerBlock>>>
        (A,sumA,B,chunk,npart,jpart,blocksPerGrid); break;
    }
    break;
  case 1:
    switch(threadsPerBlock){
    case 512:
      myCNsum_r_1<512><<<blocksPerGrid,threadsPerBlock>>>
        (A,sumA,B,chunk,npart,jpart,blocksPerGrid); break;
    case 256:
      myCNsum_r_1<256><<<blocksPerGrid,threadsPerBlock>>>
        (A,sumA,B,chunk,npart,jpart,blocksPerGrid); break;
    case 128:
      myCNsum_r_1<128><<<blocksPerGrid,threadsPerBlock>>>
        (A,sumA,B,chunk,npart,jpart,blocksPerGrid); break;
    case 64:
      myCNsum_r_1<64><<<blocksPerGrid,threadsPerBlock>>>
        (A,sumA,B,chunk,npart,jpart,blocksPerGrid); break;
    case 32:
      myCNsum_r_1<32><<<blocksPerGrid,threadsPerBlock>>>
        (A,sumA,B,chunk,npart,jpart,blocksPerGrid); break;
    case 16:
      myCNsum_r_1<16><<<blocksPerGrid,threadsPerBlock>>>
        (A,sumA,B,chunk,npart,jpart,blocksPerGrid); break;
    case 8:
      myCNsum_r_1<8><<<blocksPerGrid,threadsPerBlock>>>
        (A,sumA,B,chunk,npart,jpart,blocksPerGrid); break;
    case 4:
      myCNsum_r_1<4><<<blocksPerGrid,threadsPerBlock>>>
        (A,sumA,B,chunk,npart,jpart,blocksPerGrid); break;
    case 2:
      myCNsum_r_1<2><<<blocksPerGrid,threadsPerBlock>>>
        (A,sumA,B,chunk,npart,jpart,blocksPerGrid); break;
    case 1:
      myCNsum_r_1<1><<<blocksPerGrid,threadsPerBlock>>>
        (A,sumA,B,chunk,npart,jpart,blocksPerGrid); break;
    }
    break;
  case 2:
    switch(threadsPerBlock){
    case 512:
      myCNsum_r_2<512><<<blocksPerGrid,threadsPerBlock>>>
        (A,sumA,B,chunk,npart,jpart,blocksPerGrid); break;
    case 256:
      myCNsum_r_2<256><<<blocksPerGrid,threadsPerBlock>>>
        (A,sumA,B,chunk,npart,jpart,blocksPerGrid); break;
    case 128:
      myCNsum_r_2<128><<<blocksPerGrid,threadsPerBlock>>>
        (A,sumA,B,chunk,npart,jpart,blocksPerGrid); break;
    case 64:
      myCNsum_r_2<64><<<blocksPerGrid,threadsPerBlock>>>
        (A,sumA,B,chunk,npart,jpart,blocksPerGrid); break;
    case 32:
      myCNsum_r_2<32><<<blocksPerGrid,threadsPerBlock>>>
        (A,sumA,B,chunk,npart,jpart,blocksPerGrid); break;
    case 16:
      myCNsum_r_2<16><<<blocksPerGrid,threadsPerBlock>>>
        (A,sumA,B,chunk,npart,jpart,blocksPerGrid); break;
    case 8:
      myCNsum_r_2<8><<<blocksPerGrid,threadsPerBlock>>>
        (A,sumA,B,chunk,npart,jpart,blocksPerGrid); break;
    case 4:
      myCNsum_r_2<4><<<blocksPerGrid,threadsPerBlock>>>
        (A,sumA,B,chunk,npart,jpart,blocksPerGrid); break;
    case 2:
      myCNsum_r_2<2><<<blocksPerGrid,threadsPerBlock>>>
        (A,sumA,B,chunk,npart,jpart,blocksPerGrid); break;
    case 1:
      myCNsum_r_2<1><<<blocksPerGrid,threadsPerBlock>>>
        (A,sumA,B,chunk,npart,jpart,blocksPerGrid); break;
    }
    break;
  case 5:
    switch(threadsPerBlock){
    case 512:
      myCNsum_r_5<512><<<blocksPerGrid,threadsPerBlock>>>
        (A,sumA,B,chunk,npart,jpart,blocksPerGrid); break;
    case 256:
      myCNsum_r_5<256><<<blocksPerGrid,threadsPerBlock>>>
        (A,sumA,B,chunk,npart,jpart,blocksPerGrid); break;
    case 128:
      myCNsum_r_5<128><<<blocksPerGrid,threadsPerBlock>>>
        (A,sumA,B,chunk,npart,jpart,blocksPerGrid); break;
    case 64:
      myCNsum_r_5<64><<<blocksPerGrid,threadsPerBlock>>>
        (A,sumA,B,chunk,npart,jpart,blocksPerGrid); break;
    case 32:
      myCNsum_r_5<32><<<blocksPerGrid,threadsPerBlock>>>
        (A,sumA,B,chunk,npart,jpart,blocksPerGrid); break;
    case 16:
      myCNsum_r_5<16><<<blocksPerGrid,threadsPerBlock>>>
        (A,sumA,B,chunk,npart,jpart,blocksPerGrid); break;
    case 8:
      myCNsum_r_5<8><<<blocksPerGrid,threadsPerBlock>>>
        (A,sumA,B,chunk,npart,jpart,blocksPerGrid); break;
    case 4:
      myCNsum_r_5<4><<<blocksPerGrid,threadsPerBlock>>>
        (A,sumA,B,chunk,npart,jpart,blocksPerGrid); break;
    case 2:
      myCNsum_r_5<2><<<blocksPerGrid,threadsPerBlock>>>
        (A,sumA,B,chunk,npart,jpart,blocksPerGrid); break;
    case 1:
      myCNsum_r_5<1><<<blocksPerGrid,threadsPerBlock>>>
        (A,sumA,B,chunk,npart,jpart,blocksPerGrid); break;
    }
    break;
  case 6:
    switch(threadsPerBlock){
    case 512:
      myCNsum_r_6<512><<<blocksPerGrid,threadsPerBlock>>>
        (A,sumA,B,chunk,npart,jpart,blocksPerGrid); break;
    case 256:
      myCNsum_r_6<256><<<blocksPerGrid,threadsPerBlock>>>
        (A,sumA,B,chunk,npart,jpart,blocksPerGrid); break;
    case 128:
      myCNsum_r_6<128><<<blocksPerGrid,threadsPerBlock>>>
        (A,sumA,B,chunk,npart,jpart,blocksPerGrid); break;
    case 64:
      myCNsum_r_6<64><<<blocksPerGrid,threadsPerBlock>>>
        (A,sumA,B,chunk,npart,jpart,blocksPerGrid); break;
    case 32:
      myCNsum_r_6<32><<<blocksPerGrid,threadsPerBlock>>>
        (A,sumA,B,chunk,npart,jpart,blocksPerGrid); break;
    case 16:
      myCNsum_r_6<16><<<blocksPerGrid,threadsPerBlock>>>
        (A,sumA,B,chunk,npart,jpart,blocksPerGrid); break;
    case 8:
      myCNsum_r_6<8><<<blocksPerGrid,threadsPerBlock>>>
        (A,sumA,B,chunk,npart,jpart,blocksPerGrid); break;
    case 4:
      myCNsum_r_6<4><<<blocksPerGrid,threadsPerBlock>>>
        (A,sumA,B,chunk,npart,jpart,blocksPerGrid); break;
    case 2:
      myCNsum_r_6<2><<<blocksPerGrid,threadsPerBlock>>>
        (A,sumA,B,chunk,npart,jpart,blocksPerGrid); break;
    case 1:
      myCNsum_r_6<1><<<blocksPerGrid,threadsPerBlock>>>
        (A,sumA,B,chunk,npart,jpart,blocksPerGrid); break;
    }
    break;
  default:
    switch(threadsPerBlock){
    case 512:
      myCNsum_r<512><<<blocksPerGrid,threadsPerBlock>>>
        (A,sumA,B,chunk,npart,jpart,blocksPerGrid); break;
    case 256:
      myCNsum_r<256><<<blocksPerGrid,threadsPerBlock>>>
        (A,sumA,B,chunk,npart,jpart,blocksPerGrid); break;
    case 128:
      myCNsum_r<128><<<blocksPerGrid,threadsPerBlock>>>
        (A,sumA,B,chunk,npart,jpart,blocksPerGrid); break;
    case 64:
      myCNsum_r<64><<<blocksPerGrid,threadsPerBlock>>>
        (A,sumA,B,chunk,npart,jpart,blocksPerGrid); break;
    case 32:
      myCNsum_r<32><<<blocksPerGrid,threadsPerBlock>>>
        (A,sumA,B,chunk,npart,jpart,blocksPerGrid); break;
    case 16:
      myCNsum_r<16><<<blocksPerGrid,threadsPerBlock>>>
        (A,sumA,B,chunk,npart,jpart,blocksPerGrid); break;
    case 8:
      myCNsum_r<8><<<blocksPerGrid,threadsPerBlock>>>
        (A,sumA,B,chunk,npart,jpart,blocksPerGrid); break;
    case 4:
      myCNsum_r<4><<<blocksPerGrid,threadsPerBlock>>>
        (A,sumA,B,chunk,npart,jpart,blocksPerGrid); break;
    case 2:
      myCNsum_r<2><<<blocksPerGrid,threadsPerBlock>>>
        (A,sumA,B,chunk,npart,jpart,blocksPerGrid); break;
    case 1:
      myCNsum_r<1><<<blocksPerGrid,threadsPerBlock>>>
        (A,sumA,B,chunk,npart,jpart,blocksPerGrid); break;
    }
    break;
  }

  return rc;
}
//---------------------------------------------------------------------
// The sqsum including the atomic product Re(A*conjg(A))
//---------------------------------------------------------------------
/* summing the 3D matrix treating it a 1D vec */
template <unsigned int blockSize>
__global__ void
myCNsum_absq(const cuFloatComplex *A, float *sumA, int chunk, int npart, int jpart,
          int blocksPerGrid){

  int i = threadIdx.x + blockIdx.x * blockDim.x;

  int sharedIdx = threadIdx.x;
  __shared__ float shared_A[blockSize];
  int offset = -1;
  
  float temp = 0.0;
  while ( i < chunk) {
    offset++;
    temp += cuReCCstarf(A[i]);//A[i];//A[i*npart + offset];
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

}
/* summing the 3D matrix treating it a 1D vec */
/* using a slow modulo to select entries */
template <unsigned int blockSize>
__global__ void
myCNsum_absq_0(const cuFloatComplex *A, float *sumA, int chunk, int npart, int jpart,
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
    temp += cuReCCstarf(A[i]);//A[i];//A[i*npart + offset];
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

}
/* summing the 3D matrix treating it a 1D vec, interleaved addressing */
/* no divergent branch in the for loop replacing the slow modulo */
template <unsigned int blockSize>
__global__ void
myCNsum_absq_1(const cuFloatComplex *A, float *sumA, int chunk, int npart, int jpart,
            int blocksPerGrid){

  __shared__ float shared_A[blockSize];
  unsigned int tid = threadIdx.x;
  unsigned int i = threadIdx.x + blockIdx.x * blockDim.x;

  int offset = -1;
  float temp = 0.0 ;
  while (i<chunk) {
    offset++;
    temp += cuReCCstarf(A[i]);//A[i];//A[i*npart + offset];
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
/* which is similar to myCNsum_3D kernel above but some idle threads*/
template <unsigned int blockSize>
__global__ void
myCNsum_absq_2(const cuFloatComplex *A, float *sumA, int chunk, int npart, int jpart,
            int blocksPerGrid){

  __shared__ float shared_A[blockSize];
  unsigned int tid = threadIdx.x;
  unsigned int i = threadIdx.x + blockIdx.x * blockDim.x;

  int offset = -1;
  float temp = 0.0 ;
  while (i<chunk) {
    offset++;
    temp += cuReCCstarf(A[i]);//A[i];//A[i*npart + offset];
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
/* summing the 3D matrix treating it a 1D vec, sequential addressing */
/* which is similar to myCNsum_3D kernel above but some idle threads*/
template <unsigned int blockSize>
__global__ void
myCNsum_absq_5(const cuFloatComplex *A, float *sumA, int chunk, int npart, int jpart,
            int blocksPerGrid){

 __shared__ float shared_A[blockSize];
  unsigned int tid = threadIdx.x;
  unsigned int i = threadIdx.x + blockIdx.x * blockDim.x;

  int offset = -1;
  float temp = 0.0 ;
  while (i<chunk) {
    offset++;
    temp += cuReCCstarf(A[i]);//A[i];//A[i*npart + offset];
    i+= blockDim.x * gridDim.x;
  }
  shared_A[tid] = temp;
  __syncthreads();

  /* do the reduction in shared memory */
  for (unsigned int s=blockDim.x/2 ; s>32; s>>= 1) { 
    if ( tid < s){shared_A[tid]+=shared_A[tid+s];}
    __syncthreads();
  }

  if (tid < 32 ) wrapMyCNsum_r_5(shared_A,tid);
  if (tid==0) sumA[blockIdx.x]=shared_A[0];

}
/* summing the 3D matrix treating it a 1D vec, sequential addressing */
/* which is similar to myCNsum_3D kernel above but some idle threads*/
template <unsigned int blockSize>
__global__ void
myCNsum_absq_6(const cuFloatComplex *A, float *sumA, int chunk, int npart, int jpart,
            int blocksPerGrid){

  __shared__ float shared_A[blockSize];
  unsigned int tid = threadIdx.x;
  unsigned int i = threadIdx.x + blockIdx.x * blockDim.x;

  int offset = -1;
  float temp = 0.0 ;
  while (i<chunk) {
    offset++;
    temp += cuReCCstarf(A[i]);//A[i];//A[i*npart + offset];
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
  
  if (tid < 32 ) wrapMyCNsum_r_6<blockSize>(shared_A,tid);
  if (tid==0) sumA[blockIdx.x]=shared_A[0];

}
/* getter function to get the kernel call for sum 2D for different kernels */
/*float *d_r, float *d_partial_1D, float *startptr,*/
template <typename S>
int get_myCNsum_absq(S iKernel,
                  cuFloatComplex *A, float *sumA,
                  int chunk, int npart, int jpart,
                  int blocksPerGrid, int threadsPerBlock) {
  int rc = RC_SUCCESS;

  switch (iKernel) {

  case 0:
    switch(threadsPerBlock){
    case 512: myCNsum_absq_0<512><<<blocksPerGrid,threadsPerBlock>>>
        (A,sumA,chunk,npart,jpart,blocksPerGrid); break;
    case 256:
      myCNsum_absq_0<256><<<blocksPerGrid,threadsPerBlock>>>
        (A,sumA,chunk,npart,jpart,blocksPerGrid); break;
    case 128:
      myCNsum_absq_0<128><<<blocksPerGrid,threadsPerBlock>>>
        (A,sumA,chunk,npart,jpart,blocksPerGrid); break;
    case 64:
      myCNsum_absq_0<64><<<blocksPerGrid,threadsPerBlock>>>
        (A,sumA,chunk,npart,jpart,blocksPerGrid); break;
    case 32:
      myCNsum_absq_0<32><<<blocksPerGrid,threadsPerBlock>>>
        (A,sumA,chunk,npart,jpart,blocksPerGrid); break;
    case 16:
      myCNsum_absq_0<16><<<blocksPerGrid,threadsPerBlock>>>
        (A,sumA,chunk,npart,jpart,blocksPerGrid); break;
    case 8:
      myCNsum_absq_0<8><<<blocksPerGrid,threadsPerBlock>>>
        (A,sumA,chunk,npart,jpart,blocksPerGrid); break;
    case 4:
      myCNsum_absq_0<4><<<blocksPerGrid,threadsPerBlock>>>
        (A,sumA,chunk,npart,jpart,blocksPerGrid); break;
    case 2:
      myCNsum_absq_0<2><<<blocksPerGrid,threadsPerBlock>>>
        (A,sumA,chunk,npart,jpart,blocksPerGrid); break;
    case 1:
      myCNsum_absq_0<1><<<blocksPerGrid,threadsPerBlock>>>
        (A,sumA,chunk,npart,jpart,blocksPerGrid); break;
    }
    break;
  case 1:
    switch(threadsPerBlock){
    case 512:
      myCNsum_absq_1<512><<<blocksPerGrid,threadsPerBlock>>>
        (A,sumA,chunk,npart,jpart,blocksPerGrid); break;
    case 256:
      myCNsum_absq_1<256><<<blocksPerGrid,threadsPerBlock>>>
        (A,sumA,chunk,npart,jpart,blocksPerGrid); break;
    case 128:
      myCNsum_absq_1<128><<<blocksPerGrid,threadsPerBlock>>>
        (A,sumA,chunk,npart,jpart,blocksPerGrid); break;
    case 64:
      myCNsum_absq_1<64><<<blocksPerGrid,threadsPerBlock>>>
        (A,sumA,chunk,npart,jpart,blocksPerGrid); break;
    case 32:
      myCNsum_absq_1<32><<<blocksPerGrid,threadsPerBlock>>>
        (A,sumA,chunk,npart,jpart,blocksPerGrid); break;
    case 16:
      myCNsum_absq_1<16><<<blocksPerGrid,threadsPerBlock>>>
        (A,sumA,chunk,npart,jpart,blocksPerGrid); break;
    case 8:
      myCNsum_absq_1<8><<<blocksPerGrid,threadsPerBlock>>>
        (A,sumA,chunk,npart,jpart,blocksPerGrid); break;
    case 4:
      myCNsum_absq_1<4><<<blocksPerGrid,threadsPerBlock>>>
        (A,sumA,chunk,npart,jpart,blocksPerGrid); break;
    case 2:
      myCNsum_absq_1<2><<<blocksPerGrid,threadsPerBlock>>>
        (A,sumA,chunk,npart,jpart,blocksPerGrid); break;
    case 1:
      myCNsum_absq_1<1><<<blocksPerGrid,threadsPerBlock>>>
        (A,sumA,chunk,npart,jpart,blocksPerGrid); break;
    }
    break;
  case 2:
    switch(threadsPerBlock){
    case 512:
      myCNsum_absq_2<512><<<blocksPerGrid,threadsPerBlock>>>
        (A,sumA,chunk,npart,jpart,blocksPerGrid); break;
    case 256:
      myCNsum_absq_2<256><<<blocksPerGrid,threadsPerBlock>>>
        (A,sumA,chunk,npart,jpart,blocksPerGrid); break;
    case 128:
      myCNsum_absq_2<128><<<blocksPerGrid,threadsPerBlock>>>
        (A,sumA,chunk,npart,jpart,blocksPerGrid); break;
    case 64:
      myCNsum_absq_2<64><<<blocksPerGrid,threadsPerBlock>>>
        (A,sumA,chunk,npart,jpart,blocksPerGrid); break;
    case 32:
      myCNsum_absq_2<32><<<blocksPerGrid,threadsPerBlock>>>
        (A,sumA,chunk,npart,jpart,blocksPerGrid); break;
    case 16:
      myCNsum_absq_2<16><<<blocksPerGrid,threadsPerBlock>>>
        (A,sumA,chunk,npart,jpart,blocksPerGrid); break;
    case 8:
      myCNsum_absq_2<8><<<blocksPerGrid,threadsPerBlock>>>
        (A,sumA,chunk,npart,jpart,blocksPerGrid); break;
    case 4:
      myCNsum_absq_2<4><<<blocksPerGrid,threadsPerBlock>>>
        (A,sumA,chunk,npart,jpart,blocksPerGrid); break;
    case 2:
      myCNsum_absq_2<2><<<blocksPerGrid,threadsPerBlock>>>
        (A,sumA,chunk,npart,jpart,blocksPerGrid); break;
    case 1:
      myCNsum_absq_2<1><<<blocksPerGrid,threadsPerBlock>>>
        (A,sumA,chunk,npart,jpart,blocksPerGrid); break;
    }
    break;
  case 5:
    switch(threadsPerBlock){
    case 512:
      myCNsum_absq_5<512><<<blocksPerGrid,threadsPerBlock>>>
        (A,sumA,chunk,npart,jpart,blocksPerGrid); break;
    case 256:
      myCNsum_absq_5<256><<<blocksPerGrid,threadsPerBlock>>>
        (A,sumA,chunk,npart,jpart,blocksPerGrid); break;
    case 128:
      myCNsum_absq_5<128><<<blocksPerGrid,threadsPerBlock>>>
        (A,sumA,chunk,npart,jpart,blocksPerGrid); break;
    case 64:
      myCNsum_absq_5<64><<<blocksPerGrid,threadsPerBlock>>>
        (A,sumA,chunk,npart,jpart,blocksPerGrid); break;
    case 32:
      myCNsum_absq_5<32><<<blocksPerGrid,threadsPerBlock>>>
        (A,sumA,chunk,npart,jpart,blocksPerGrid); break;
    case 16:
      myCNsum_absq_5<16><<<blocksPerGrid,threadsPerBlock>>>
        (A,sumA,chunk,npart,jpart,blocksPerGrid); break;
    case 8:
      myCNsum_absq_5<8><<<blocksPerGrid,threadsPerBlock>>>
        (A,sumA,chunk,npart,jpart,blocksPerGrid); break;
    case 4:
      myCNsum_absq_5<4><<<blocksPerGrid,threadsPerBlock>>>
        (A,sumA,chunk,npart,jpart,blocksPerGrid); break;
    case 2:
      myCNsum_absq_5<2><<<blocksPerGrid,threadsPerBlock>>>
        (A,sumA,chunk,npart,jpart,blocksPerGrid); break;
    case 1:
      myCNsum_absq_5<1><<<blocksPerGrid,threadsPerBlock>>>
        (A,sumA,chunk,npart,jpart,blocksPerGrid); break;
    }
    break;
  case 6:
    switch(threadsPerBlock){
    case 512:
      myCNsum_absq_6<512><<<blocksPerGrid,threadsPerBlock>>>
        (A,sumA,chunk,npart,jpart,blocksPerGrid); break;
    case 256:
      myCNsum_absq_6<256><<<blocksPerGrid,threadsPerBlock>>>
        (A,sumA,chunk,npart,jpart,blocksPerGrid); break;
    case 128:
      myCNsum_absq_6<128><<<blocksPerGrid,threadsPerBlock>>>
        (A,sumA,chunk,npart,jpart,blocksPerGrid); break;
    case 64:
      myCNsum_absq_6<64><<<blocksPerGrid,threadsPerBlock>>>
        (A,sumA,chunk,npart,jpart,blocksPerGrid); break;
    case 32:
      myCNsum_absq_6<32><<<blocksPerGrid,threadsPerBlock>>>
        (A,sumA,chunk,npart,jpart,blocksPerGrid); break;
    case 16:
      myCNsum_absq_6<16><<<blocksPerGrid,threadsPerBlock>>>
        (A,sumA,chunk,npart,jpart,blocksPerGrid); break;
    case 8:
      myCNsum_absq_6<8><<<blocksPerGrid,threadsPerBlock>>>
        (A,sumA,chunk,npart,jpart,blocksPerGrid); break;
    case 4:
      myCNsum_absq_6<4><<<blocksPerGrid,threadsPerBlock>>>
        (A,sumA,chunk,npart,jpart,blocksPerGrid); break;
    case 2:
      myCNsum_absq_6<2><<<blocksPerGrid,threadsPerBlock>>>
        (A,sumA,chunk,npart,jpart,blocksPerGrid); break;
    case 1:
      myCNsum_absq_6<1><<<blocksPerGrid,threadsPerBlock>>>
        (A,sumA,chunk,npart,jpart,blocksPerGrid); break;
    }
    break;
  default:
    switch(threadsPerBlock){
    case 512:
      myCNsum_absq<512><<<blocksPerGrid,threadsPerBlock>>>
        (A,sumA,chunk,npart,jpart,blocksPerGrid); break;
    case 256:
      myCNsum_absq<256><<<blocksPerGrid,threadsPerBlock>>>
        (A,sumA,chunk,npart,jpart,blocksPerGrid); break;
    case 128:
      myCNsum_absq<128><<<blocksPerGrid,threadsPerBlock>>>
        (A,sumA,chunk,npart,jpart,blocksPerGrid); break;
    case 64:
      myCNsum_absq<64><<<blocksPerGrid,threadsPerBlock>>>
        (A,sumA,chunk,npart,jpart,blocksPerGrid); break;
    case 32:
      myCNsum_absq<32><<<blocksPerGrid,threadsPerBlock>>>
        (A,sumA,chunk,npart,jpart,blocksPerGrid); break;
    case 16:
      myCNsum_absq<16><<<blocksPerGrid,threadsPerBlock>>>
        (A,sumA,chunk,npart,jpart,blocksPerGrid); break;
    case 8:
      myCNsum_absq<8><<<blocksPerGrid,threadsPerBlock>>>
        (A,sumA,chunk,npart,jpart,blocksPerGrid); break;
    case 4:
      myCNsum_absq<4><<<blocksPerGrid,threadsPerBlock>>>
        (A,sumA,chunk,npart,jpart,blocksPerGrid); break;
    case 2:
      myCNsum_absq<2><<<blocksPerGrid,threadsPerBlock>>>
        (A,sumA,chunk,npart,jpart,blocksPerGrid); break;
    case 1:
      myCNsum_absq<1><<<blocksPerGrid,threadsPerBlock>>>
        (A,sumA,chunk,npart,jpart,blocksPerGrid); break;
    }
    break;
  }

  return rc;
}
//---------------------------------------------------------------------
// The main kernel driver  atomic product Re(A*conjg(B))
//---------------------------------------------------------------------
/* main kernel entry */
extern "C" int
carte2d_ftExt_corr_C_N(deviceDetails_t * s_devD,
                       polar_corr_calc_t *s_carte,
                       float *r,
                       cuFloatComplex *shmat,
                       const cuFloatComplex *A,
                       const cuFloatComplex *B,
                       int vx, int vy, int vz,
                       float alpha,
                       bench_t *s_bench, debug_gpu_t *s_debug_gpu)
{
  /* return code */
  int rc = RC_SUCCESS; 
  /* setting the values into the object */
  img_2D_cart_Sizes *p_img_2D_cart_Sizes = new img_2D_cart_Sizes(vx,vy,vz);
  mesh_3D           *p_mesh_3D           = new mesh_3D(s_carte,vx,vy,vz);
  mesh_1DV          *p_mesh_1DV          = new mesh_1DV(s_carte,vx,vy,vz);
  /* the error handlers from the cuda library */
  cudaError_t err;
  /* grid and threads block definition */
  int nx, ny, nz;               //Dimension of the threads
  int gridx, gridy, gridz;      //Mesh of the 3D Grid
  /* device allocation */
  cuFloatComplex *d_cmat2sh = NULL;   //device pointer for cmat2sh (const in)
  cuFloatComplex   *d_shmat = NULL;   //device pointer for shmat (const in)
  cuFloatComplex       *d_A = NULL;   //device pointer for A (const in)
  cuFloatComplex       *d_B = NULL;   //device pointer for B (const in)
  /* creating the working space on the GPU */
  float *d_partial_1D_r = NULL;
  float *d_partial_1D_a = NULL;
  float *d_partial_1D_b = NULL;
  /*The summed value mapped into the s_carte data structure*/
  float c_r;
  float c_a;
  float c_b;
  /* size of the element in consideration */
  int size_cmat2sh = p_img_2D_cart_Sizes->get_npts()  *sizeof(cuFloatComplex);
  int  size_dshmat = p_img_2D_cart_Sizes->get_npts()  *sizeof(cuFloatComplex);
  int      size_dA = p_img_2D_cart_Sizes->get_npts_A()*sizeof(cuFloatComplex);
  int      size_dB = p_img_2D_cart_Sizes->get_npts_B()*sizeof(cuFloatComplex);
  /* indexer */
  //int i,j,k;//, ipart, jpart, irot;
  /* printer index depth */
  float depth = 2;

  /*start of the execution commands */
  /* setting up the mesh using the s_polar structure */
  nx=s_carte->nx; ny=s_carte->ny; nz=s_carte->nz;
  dim3 threads(nx,ny,nz); /*mesh dimensions for the hadamard product*/
  gridx = p_mesh_3D->get_mesh3D_gridx();
  gridy = p_mesh_3D->get_mesh3D_gridy();
  gridz = p_mesh_3D->get_mesh3D_gridz(); /*checking the mesh*/
  rc = check_grid3D_V(gridx,gridy,gridz,s_carte,p_img_2D_cart_Sizes);
  //if (rc != 0 ) {rc = 
  dim3 grid(gridx, gridy, gridz);
  if (s_debug_gpu->debug_i == true ) {
    rc = print_3D_V_mesh(0,p_mesh_3D,s_carte,p_img_2D_cart_Sizes,nx,ny,nz);//gridx, gridy, gridz);
    /* printing input variables on screen */
    if ( s_debug_gpu->debug_high_i == true ) {
      rc = print_s_debug_struct(s_debug_gpu);
      rc = print_s_bench_struct(s_bench);
      rc = print_function_header_C_N_info(s_devD,
                                          s_carte,
                                          r,shmat,
                                          A, B, 
                                          p_img_2D_cart_Sizes,
                                          vx, vy, vz,
                                          alpha, depth);
    }
  }

  /* first on cpu */
  if (s_debug_gpu->debug_cpu_i == true ) {
    //TODO: need a cpu version to check result
    /* rc = polarft_gencorrAll_Z_N_cpu(....)*/}

  /* allocating the memory on GPU for d_{A,B,C}, d_{shmat,cmat2sh} */
  err = cudaMalloc((void**)&d_cmat2sh, size_cmat2sh);
  err = cudaMalloc((void**)&d_shmat, size_dshmat);
  err = cudaMalloc((void**)&d_A, size_dA);
  err = cudaMalloc((void**)&d_B, size_dB);
  if ( (int)err != CUDA_SUCCESS ) {rc = get_error_id_corr_Hadmr_gpu(err); }
  
  //uploading the matrices shmat, A, B to GPU
  err = cudaMemcpy(d_A, A, size_dA, cudaMemcpyHostToDevice);
  err = cudaMemcpy(d_B, B, size_dB, cudaMemcpyHostToDevice);
  err = cudaMemcpy(d_shmat, shmat, size_dshmat, cudaMemcpyHostToDevice);
  if ( (int)err != CUDA_SUCCESS ) {rc = get_error_id_corr_Hadmr_gpu(err); }

  /*First getting the product cmat2sh with shmat*/
  NXprod_3D_mat<<<grid,threads>>>(d_cmat2sh, d_B, d_shmat,
                                  p_img_2D_cart_Sizes->get_2D_vx(),
                                  p_img_2D_cart_Sizes->get_2D_vy(),
                                  p_img_2D_cart_Sizes->get_2D_vz() );
  cuCtxSynchronize();

  if (s_debug_gpu->debug_i == true ) {
    rc = print_cmat2sh_V_depth(d_cmat2sh, p_img_2D_cart_Sizes, vx, vy,vz,
                               depth);
  }

  //Preparing mesh for summing the over the chunk and mapping it to vector r  
  int N = p_mesh_1DV->get_mesh1DV_chunk();
  int threadsPerBlock = p_mesh_1DV->get_mesh1DV_threadsPerBlock();
  int blocksPerGrid = p_mesh_1DV->get_mesh1DV_blocksPerGrid();
  int size_p_c = blocksPerGrid*sizeof(float);
  if (s_debug_gpu->debug_i == true ) {
    rc = print_summing_method(0);
    rc = print_1D_mesh(0, NULL, N, threadsPerBlock, blocksPerGrid);
    rc = print_iKernel_threadsPerBlock(s_carte, N,
                                       threadsPerBlock, blocksPerGrid);
  }

  //computing the r product now summing the r 
  float *h_partial_1D_r = (float*)malloc(size_p_c);
  err = cudaMalloc((void**)&d_partial_1D_r, size_p_c);
  
  rc = get_myCNsum_r<int>(s_carte->ikrnl, d_A,d_partial_1D_r, d_cmat2sh,
                          N, vx, 0, blocksPerGrid,threadsPerBlock);
  cuCtxSynchronize();
  err = cudaMemcpy(h_partial_1D_r, d_partial_1D_r, size_p_c, cudaMemcpyDeviceToHost);

  c_r = 0.0;
  for ( int igrid = 0 ; igrid < blocksPerGrid ; igrid++ ) {  
    c_r += h_partial_1D_r[igrid];
  }
  s_carte->r_polar = c_r;
  free(h_partial_1D_r);
  cudaFree(d_partial_1D_r);

  //computing the sum_asq=sum(A*Astar)
  float *h_partial_1D_a = (float*)malloc(size_p_c);
  err = cudaMalloc((void**)&d_partial_1D_a, size_p_c);
  
  rc = get_myCNsum_absq<int>(s_carte->ikrnl, d_A, d_partial_1D_a,
                             N, vx, 0, blocksPerGrid,threadsPerBlock);
  cuCtxSynchronize();
  err = cudaMemcpy(h_partial_1D_a, d_partial_1D_a, size_p_c, cudaMemcpyDeviceToHost);

  c_a = 0.0;
  for ( int igrid = 0 ; igrid < blocksPerGrid ; igrid++ ) {  
    c_a += h_partial_1D_a[igrid];
  }
  s_carte->sumasq_polar = c_a;
  free(h_partial_1D_a);
  cudaFree(d_partial_1D_a);

  //computing the sum_asq=sum(cmat2sh*cmatshtar)
  
  //computing the sum_asq=sum(A*Astar)
  float *h_partial_1D_b = (float*)malloc(size_p_c);
  err = cudaMalloc((void**)&d_partial_1D_b, size_p_c);
  
  rc = get_myCNsum_absq<int>(s_carte->ikrnl, d_cmat2sh, d_partial_1D_b,
                             N, vx, 0, blocksPerGrid,threadsPerBlock);
  cuCtxSynchronize();
  err = cudaMemcpy(h_partial_1D_b, d_partial_1D_b, size_p_c, cudaMemcpyDeviceToHost);

  c_b = 0.0;
  for ( int igrid = 0 ; igrid < blocksPerGrid ; igrid++ ) {  
    c_b += h_partial_1D_b[igrid];
  }
  s_carte->sumbsq_polar = c_b;
  free(h_partial_1D_b);
  cudaFree(d_partial_1D_b);
  
  /* freeing the resources on the GPU */
  err = cudaFree(d_cmat2sh);
  err = cudaFree(d_shmat);
  err = cudaFree(d_A);
  err = cudaFree(d_B);
  //err = cudaFree(d_r);
  if ( (int)err != CUDA_SUCCESS ) {rc = get_error_id_corr_Hadmr_gpu(err); }
 
  /* calling the destructor for the pfts_Sizes object */
  p_img_2D_cart_Sizes->~img_2D_cart_Sizes();  
  p_mesh_3D->~mesh_3D();
  p_mesh_1DV->~mesh_1DV();

  return rc;
 
} /* End of carte2d_ftExt_corr_C_N */

#endif /* CUDA */
