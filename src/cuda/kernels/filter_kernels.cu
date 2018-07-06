

#include <stdio.h>
#include "cuda.h"
#include "cuda_runtime_api.h"


__global__ void gauss3DFTconvolve(float *dest, const float *src, const int *matSize, const float *s, const float factor )
{
    unsigned int idx = blockIdx.x*blockDim.x + threadIdx.x;
    unsigned int idy = blockIdx.y*blockDim.y + threadIdx.y;
    unsigned int idz = blockIdx.z*blockDim.z + threadIdx.z;
  float sigma[3];
  sigma[0]=s[0];sigma[1]=s[1];sigma[2]=s[2];
  const double TWOPISQ = 19.739208802178716;


  float i = (float)((idx / matSize[2]) / matSize[1]);
      i = (i - (float)floor((float)(matSize[0])/2.0))/(float)(matSize[0]);
  float j = (float)(idx / matSize[2]);
      if(j > matSize[1]) {j=(float)fmod(j, (float)matSize[1]);};
      j = (j - (float)floor((float)(matSize[1])/2.0f))/(float)(matSize[1]);
  //Account for large global index (stored as ulong) before performing modulus
  double pre_k=fmod((double)(idx) , (double) matSize[2]);
  float k = (float) pre_k;
      k = (k - (float)floor((float)(matSize[2])/2.0f))/(float)(matSize[2]);
  float weight = exp(-TWOPISQ*((i*i)*sigma[0]*sigma[0] + (j*j)*sigma[1]*sigma[1] + (k*k)*sigma[2]*sigma[2]));
  dest[idx] = src[idx] * weight;

}

__global__ void gaussElementwiseKernel(
     const float *u, const float *v, const float *w, float *z,  const float s, int N )
{
  const double TWOPISQ = 19.739208802178716;
  int i = threadIdx.x + blockIdx.x * blockDim.x;
  if ( i < N ){
    double x =  (u[i]*u[i] + v[i]* v[i] + w[i]*w[i]) * (s*s);
    z[i] = (float) exp( -( TWOPISQ * x) );
   }
}



extern "C"
{
#define GAUSSKERNELELEMENTWISE gausskernelelementwise_
void  gausskernelelementwise(float *A, float *B, float *C, float *Z, float *s, dim3 *dimGrid, dim3 *dimBlk,
                  int N, cudaStream_t *stream)
{
  gaussElementwiseKernel<<<*dimGrid, *dimBlk, 0, *stream>>>(A, B, C, Z,*s, N);
}
#define GAUSSCONVOLUTION3D gaussconvolution3d_
void  gaussconvolution3d(float *A, float *B,  float *sigma, dim3 *dimGrid, dim3 *dimBlk,
                         int* N, cudaStream_t *stream)
{
  gauss3DFTconvolve<<<*dimGrid, *dimBlk, 0, *stream>>>(A, B, N, sigma, 1.0);
}


}
