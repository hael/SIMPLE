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
 * SMATMUL.
 **********************************************************************/
int smatmul_gpu_(const char *transa, const char *transb, const int *m,
                 const int *n, const int *k, const float *alpha,
                 const devptr_t *devPtrA, const int *lda, 
                 const devptr_t *devPtrB, const int *ldb, const float *beta,
                 devptr_t *devPtrC, const int *ldc) {

  int rc = RC_SUCCESS;
  
  float *A = (float*)*devPtrA;
  float *B = (float*)*devPtrB;
  float *C = (float*)*devPtrC;

  cublasSgemm(*transa,*transb,*m,*n,*k,*alpha,A,*lda,B,*ldb,*beta,C,*ldc);

  return rc;
}
void SMATMUL_GPU_() __attribute__((weak,alias("smatmul_gpu_")));
void smatmul_gpu__() __attribute__((weak,alias("smatmul_gpu_")));
void SMATMUL_GPU__() __attribute__((weak,alias("smatmul_gpu_")));
/**********************************************************************
 * DMATMUL.
 **********************************************************************/
int dmatmul_gpu_(const char *transa, const char *transb, const int *m,
                 const int *n, const int *k, const double *alpha,
                 const devptr_t *devPtrA, const int *lda, 
                 const devptr_t *devPtrB, const int *ldb, const double *beta,
                 devptr_t *devPtrC, const int *ldc) {

  int rc = RC_SUCCESS;

  double *A = (double*)*devPtrA;
  double *B = (double*)*devPtrB;
  double *C = (double*)*devPtrC;

  cublasDgemm(*transa,*transb,*m,*n,*k,*alpha,A,*lda,B,*ldb,*beta,C,*ldc);

  return rc;
}
void DMATMUL_GPU_() __attribute__((weak,alias("dmatmul_gpu_")));
void dmatmul_gpu__() __attribute__((weak,alias("dmatmul_gpu_")));
void DMATMUL_GPU__() __attribute__((weak,alias("dmatmul_gpu_")));
/**********************************************************************
 * CMATMUL.
 **********************************************************************/
int cmatmul_gpu_(const char *transa, const char *transb, const int *m,
                 const int *n, const int *k, const cuComplex *alpha,
                 const devptr_t *devPtrA, const int *lda, 
                 const devptr_t *devPtrB, const int *ldb, const cuComplex *beta,
                 devptr_t *devPtrC, const int *ldc) {

  int rc = RC_SUCCESS;

  cuComplex *A = (cuComplex*)*devPtrA;
  cuComplex *B = (cuComplex*)*devPtrB;
  cuComplex *C = (cuComplex*)*devPtrC;

  cublasCgemm(*transa,*transb,*m,*n,*k,*alpha,A,*lda,B,*ldb,*beta,C,*ldc);

  return rc;
}
void CMATMUL_GPU_() __attribute__((weak,alias("cmatmul_gpu_")));
void cmatmul_gpu__() __attribute__((weak,alias("cmatmul_gpu_")));
void CMATMUL_GPU__() __attribute__((weak,alias("cmatmul_gpu_")));
/**********************************************************************
 * ZMATMUL.
 **********************************************************************/
int zmatmul_gpu_(const char *transa, const char *transb, const int *m,
                 const int *n, const int *k, const cuDoubleComplex *alpha,
                 const devptr_t *devPtrA, const int *lda, 
                 const devptr_t *devPtrB, const int *ldb,
                 const cuDoubleComplex *beta,
                 devptr_t *devPtrC, const int *ldc) {

  int rc = RC_SUCCESS;

  cuDoubleComplex *A = (cuDoubleComplex*)*devPtrA;
  cuDoubleComplex *B = (cuDoubleComplex*)*devPtrB;
  cuDoubleComplex *C = (cuDoubleComplex*)*devPtrC;

  cublasZgemm(*transa,*transb,*m,*n,*k,*alpha,A,*lda,B,*ldb,*beta,C,*ldc);

  return rc;
}
void ZMATMUL_GPU_() __attribute__((weak,alias("zmatmul_gpu_")));
void zmatmul_gpu__() __attribute__((weak,alias("zmatmul_gpu_")));
void ZMATMUL_GPU__() __attribute__((weak,alias("zmatmul_gpu_")));

#endif /* CUDA */
