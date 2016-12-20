/*
 *   -- SIMPLE addon
 *      Author: Frederic Bonnet, Date: 17th Mar 2016
 *      Monash University
 *
 *      March 2016
 *
 *      wrappers for CUDA/MAGMA matrix multiplication with transfers.
 *
 *      c code for the MACOSX code this wrapps around the the cpp for LINUX
 *      files methods for the cpp files to extract the information from the GPU
 *      devices.
 *
 * @precisions normal z -> s d c
 */

#include "simple.h"

#if defined (CUDA) /*preprossing for the CUDA environment */

/**********************************************************************
 * SMATVEC.
 **********************************************************************/
int smatvec_gpu_(const char *trans, const int *m,
                 const int *n, const float *alpha,
                 const devptr_t *devPtrA, const int *lda, 
                 const devptr_t *devPtrX, const int *incx, const float *beta,
                 devptr_t *devPtrY, const int *incy) {

  int rc = RC_SUCCESS;
  
  float *A = (float*)*devPtrA;
  float *X = (float*)*devPtrX;
  float *Y = (float*)*devPtrY;

  cublasSgemv(*trans,*m,*n,*alpha,A,*lda,X,*incx,*beta,Y,*incy);

  return rc;
}
void SMATVEC_GPU_() __attribute__((weak,alias("smatvec_gpu_")));
void smatvec_gpu__() __attribute__((weak,alias("smatvec_gpu_")));
void SMATVEC_GPU__() __attribute__((weak,alias("smatvec_gpu_")));
/**********************************************************************
 * DMATVEC.
 **********************************************************************/
int dmatvec_gpu_(const char *trans, const int *m,
                 const int *n, const double *alpha,
                 const devptr_t *devPtrA, const int *lda, 
                 const devptr_t *devPtrX, const int *incx, const double *beta,
                 devptr_t *devPtrY, const int *incy) {

  int rc = RC_SUCCESS;
  
  double *A = (double*)*devPtrA;
  double *X = (double*)*devPtrX;
  double *Y = (double*)*devPtrY;

  cublasDgemv(*trans,*m,*n,*alpha,A,*lda,X,*incx,*beta,Y,*incy);

  return rc;
}
void DMATVEC_GPU_() __attribute__((weak,alias("dmatvec_gpu_")));
void dmatvec_gpu__() __attribute__((weak,alias("dmatvec_gpu_")));
void DMATVEC_GPU__() __attribute__((weak,alias("dmatvec_gpu_")));
/**********************************************************************
 * CMATVEC.
 **********************************************************************/
int cmatvec_gpu_(const char *trans, const int *m,
                 const int *n, const cuComplex *alpha,
                 const devptr_t *devPtrA, const int *lda, 
                 const devptr_t *devPtrX, const int *incx,
                 const cuComplex *beta,
                 devptr_t *devPtrY, const int *incy) {

  int rc = RC_SUCCESS;
  
  cuComplex *A = (cuComplex*)*devPtrA;
  cuComplex *X = (cuComplex*)*devPtrX;
  cuComplex *Y = (cuComplex*)*devPtrY;

  cublasCgemv(*trans,*m,*n,*alpha,A,*lda,X,*incx,*beta,Y,*incy);

  return rc;
}
void CMATVEC_GPU_() __attribute__((weak,alias("cmatvec_gpu_")));
void cmatvec_gpu__() __attribute__((weak,alias("cmatvec_gpu_")));
void CMATVEC_GPU__() __attribute__((weak,alias("cmatvec_gpu_")));
/**********************************************************************
 * ZMATVEC.
 **********************************************************************/
int zmatvec_gpu_(const char *trans, const int *m,
                 const int *n, const cuDoubleComplex *alpha,
                 const devptr_t *devPtrA, const int *lda, 
                 const devptr_t *devPtrX, const int *incx,
                 const cuDoubleComplex *beta,
                 devptr_t *devPtrY, const int *incy) {

  int rc = RC_SUCCESS;
  
  cuDoubleComplex *A = (cuDoubleComplex*)*devPtrA;
  cuDoubleComplex *X = (cuDoubleComplex*)*devPtrX;
  cuDoubleComplex *Y = (cuDoubleComplex*)*devPtrY;

  cublasZgemv(*trans,*m,*n,*alpha,A,*lda,X,*incx,*beta,Y,*incy);

  return rc;
}
void ZMATVEC_GPU_() __attribute__((weak,alias("zmatvec_gpu_")));
void zmatvec_gpu__() __attribute__((weak,alias("zmatvec_gpu_")));
void ZMATVEC_GPU__() __attribute__((weak,alias("zmatvec_gpu_")));

#endif /* CUDA */



