
#include "cuda.h"
#include <stdio.h>

__global__ void gauss_kernel(
const float *dest,
const float *src,
const int *dim,
const float *s, const float factor
)
{
  const ulong x = get_global_id(0);
  const SIZE_T y = get_local_id(0);
  const SIZE_T z = get_group_id(0);
  const SIZE_T dim1= dim[0];
  const SIZE_T dim2= dim[1];
  const SIZE_T dim3= dim[2];
  float sigma[3];
  sigma[0]=s[0];sigma[1]=s[1];sigma[2]=s[2];

  const double TWOPISQ = 19.739208802178716;
  const double SQRT2PI = 2.5066282746;
  const double CUBEDSQRT2PI = 15.749609945722419;
  const ulong idx = x; //z * dim2 * dim1  + y * dim1  + x;
                      //(SIZE_T)floor((float)(z * dim2 * dim1  + y * dim1  + x )/(float)dim1*dim2*dim3);
  //float i = ((float)(x) - floor((float)(dim1)/2.0))/(float)(dim1);
  //float j = ((float)(y) - floor((float)(dim2)/2.0))/(float)(dim2);
  //float k = ((float)(z) - floor((float)(dim3)/2.0))/(float)(dim3);
  float i = (float)((x / dim3) / dim2);
      i = (i - (float)floor((float)(dim1)/2.0))/(float)(dim1);
  float j = (float)(x / dim3);
      if((SIZE_T)j > dim2) {j=(float)fmod(j, (float)dim2);};
      j = (j - (float)floor((float)(dim2)/2.0f))/(float)(dim2);
  //Account for large global index (stored as ulong) before performing modulus
  double pre_k=fmod((double)(x) , (double) dim3);
  float k = (float) pre_k;
      k = (k - (float)floor((float)(dim3)/2.0f))/(float)(dim3);
  float weight = exp(-TWOPISQ*((i*i)*sigma[0]*sigma[0] + (j*j)*sigma[1]*sigma[1] + (k*k)*sigma[2]*sigma[2]));
  //float weight = expm1(-TWOPISQ*((i*i)*sigma[0]*sigma[0] + (j*j)*sigma[1]*sigma[1] + (k*k)*sigma[2]*sigma[2]))+1;
  //float weight= ${exp}(-TWOPISQ*((i*i)*sigma[0]*sigma[0] + (j*j)*sigma[1]*sigma[1] + (k*k)*sigma[2]*sigma[2]));

  dest[idx].x = src[idx].x * weight;
  dest[idx].y = src[idx].y * weight;

}

__global__ void gaussElementwiseKernel(
      const float s, const float *u, const float *v, const float *w, float *z, int N )
{
  int i = threadIdx.x + blockIdx.x * blockDim.x;
   if ( i < N ) {
 z[i] = exp(-2 * (PI ^ 2) * (u[i] ^ 2 + v[i] ^ 2 + w[i] ^ 2) * (s^2));
 }
}



extern "C"
{
  void  gaussKernel_(float *A, float *B, float *C, float *Z, float *s, dim3 *dimGrid, dim3 *dimBlk,
                   int N, cudaStream_t *stream)
  {
    gaussElementwiseKernel<<<*dimGrid, *dimBlk, 0, *stream>>>(*s,A, B, C, Z, N);
  }

}