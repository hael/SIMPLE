/*
 *   -- SIMPLE addon
 *      Author: Frederic Bonnet, Date: 29th of April 2015
 *      Monash University
 *      April 2015
 *
 *      Routine which calculates the product (element wise) of 
 *      cuDoubleComplex matrix A and cuDoubleComplex matrix B. and takes the
 *      real part of the product and puts it into a matrix of type double
 *      C = cuCreal(x) * cuCreal(y) + cuCimag(x) * cuCimag(y)
 *      C = Re( A * conjg(B) ) element wise
 *
 *      Non Special case
 * @precisions normal z -> s d c
*/
#include "common_magma.h"
#include "commonblas_zz2d.h"
#include "simple_cuDoubleComplex.h"
#include "simple.h"

#if defined (CUDA) /*preprossing for the OPENCL environment */

/*
//Slow kernel not used in computation only for error checking
extern "C" __global__  void
sumddaxpy(const cuDoubleComplex *A,
	  const cuDoubleComplex *B,
	  int lda, int ldb,
	  double *sumasq_gpu, double *sumbsq_gpu)
{
  int i;
  int j;

  sumasq_gpu[0] = 0.0;
  sumbsq_gpu[0] = 0.0;

  for( i = 0 ; i < lda ; i++ ) {
    for( j = 0 ; j < lda ; j++ ) {
      sumasq_gpu[0] += cuReCCstar( A[ i * lda + j] );
      sumbsq_gpu[0] += cuReCCstar( B[ i * ldb + j] );
    }
  }
  __syncthreads();

}
*/
//Kernel that takes element wise 
// C = alpha*A*B + beta*C
extern "C" __global__ void
main_zz2dgemm_kernel_T_T( double *C, 
			  const cuDoubleComplex *A, 
			  const cuDoubleComplex *B,
			  int m, int n,
			  int lda, int ldb, int ldc,
			  double alpha, double beta)
{

  int i = threadIdx.x + blockIdx.x * blockDim.x;
  int j = threadIdx.y + blockIdx.y * blockDim.y;

  C[i * ldc + j ] = cuDmul(alpha,cuReCCstarmul( A[ j * lda + i], B[ j * ldb + i] ));

  __syncthreads();

}

extern "C" int
zz2dgemm_kernel_T_T( double *C, 
		     const cuDoubleComplex *A, 
		     const cuDoubleComplex *B, 
		     int m, int n, int k, 
		     int lda, int ldb, int ldc, 
		     double alpha, double beta)
{
  int rc = 0; //return code
  dim3 threads(16,16);
  dim3 grid(m/16+(m%16!=0),n/16+(n%16!=0));

  main_zz2dgemm_kernel_T_T<<< grid, threads >>>(C, A, B, 
						m, n,
						lda, ldb, ldc,
						alpha, beta);
  
  //making sure that the GPU is synchronized with CPU before proceeding
  cuCtxSynchronize();
  return rc;
} /* End of zz2dgemm_kernel_T_T */

#endif /* CUDA */
