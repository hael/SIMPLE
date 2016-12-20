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
int smatmul_tgpu_(const char *transa, const char *transb, const int *m,
                 const int *n, const int *k, const float *alpha,
                 float *A, const int *lda, 
                 float *B, const int *ldb, const float *beta,
                 float *C, const int *ldc){

  int rc = RC_SUCCESS;
  
  float *devPtrA = NULL;
  float *devPtrB = NULL;
  float *devPtrC = NULL;

  //input matrices
  rc = cublasAlloc(*m*(*k), size_of_float, (void**)&devPtrA);
  if (rc != RC_SUCCESS ) {rc = simple_cudblas_stat_return_c(rc);}
  rc = cublasAlloc(*k*(*n), size_of_float, (void**)&devPtrB);
  if (rc != RC_SUCCESS ) {rc = simple_cudblas_stat_return_c(rc);}

  rc = cublasSetMatrix (*m,*k,size_of_float,A,*lda,devPtrA,*lda);
  if (rc != RC_SUCCESS ) {rc = simple_cudblas_stat_return_c(rc);}
  rc = cublasSetMatrix (*k,*n,size_of_float,B,*ldb,devPtrB,*ldb);
  if (rc != RC_SUCCESS ) {rc = simple_cudblas_stat_return_c(rc);}
    
  /* !return matrix */
  rc = cublasAlloc(*m*(*n), size_of_float, (void**)&devPtrC);
  if (rc != RC_SUCCESS ) {rc = simple_cudblas_stat_return_c(rc);}

  cublasSgemm(*transa,*transb,*m,*n,*k,
              *alpha,
              devPtrA,*lda,
              devPtrB,*ldb,
              *beta,
              devPtrC,*ldc);

  /* retrieving back from GPU */
  rc = cublasGetMatrix( *m, *n, size_of_float, devPtrC, *lda, C, *ldc);
  if (rc != RC_SUCCESS ) {rc = simple_cudblas_stat_return_c(rc);}

  /* freeing ressources on the device */
  rc = cublasFree(devPtrA);
  rc = cublasFree(devPtrB);
  rc = cublasFree(devPtrC);

  return rc;
}
void SMATMUL_TGPU_() __attribute__((weak,alias("smatmul_tgpu_")));
void smatmul_tgpu__() __attribute__((weak,alias("smatmul_tgpu_")));
void SMATMUL_TGPU__() __attribute__((weak,alias("smatmul_tgpu_")));
/**********************************************************************
 * DMATMUL.
 **********************************************************************/
int dmatmul_tgpu_(const char *transa, const char *transb, const int *m,
                 const int *n, const int *k, const double *alpha,
                 double *A, const int *lda, 
                 double *B, const int *ldb, const double *beta,
                 double *C, const int *ldc){

  int rc = RC_SUCCESS;
  
  double *devPtrA = NULL;
  double *devPtrB = NULL;
  double *devPtrC = NULL;

  //input matrices
  rc = cublasAlloc(*m*(*k), size_of_double, (void**)&devPtrA);
  if (rc != RC_SUCCESS ) {rc = simple_cudblas_stat_return_c(rc);}
  rc = cublasAlloc(*k*(*n), size_of_double, (void**)&devPtrB);
  if (rc != RC_SUCCESS ) {rc = simple_cudblas_stat_return_c(rc);}

  rc = cublasSetMatrix (*m,*k,size_of_double,A,*lda,devPtrA,*lda);
  if (rc != RC_SUCCESS ) {rc = simple_cudblas_stat_return_c(rc);}
  rc = cublasSetMatrix (*k,*n,size_of_double,B,*ldb,devPtrB,*ldb);
  if (rc != RC_SUCCESS ) {rc = simple_cudblas_stat_return_c(rc);}
    
  /* !return matrix */
  rc = cublasAlloc(*m*(*n), size_of_double, (void**)&devPtrC);
  if (rc != RC_SUCCESS ) {rc = simple_cudblas_stat_return_c(rc);}

  cublasDgemm(*transa,*transb,*m,*n,*k,
              *alpha,
              devPtrA,*lda,
              devPtrB,*ldb,
              *beta,
              devPtrC,*ldc);

  /* retrieving back from GPU */
  rc = cublasGetMatrix( *m, *n, size_of_double, devPtrC, *lda, C, *ldc);
  if (rc != RC_SUCCESS ) {rc = simple_cudblas_stat_return_c(rc);}

  /* freeing ressources on the device */
  rc = cublasFree(devPtrA);
  rc = cublasFree(devPtrB);
  rc = cublasFree(devPtrC);

  return rc;
}
void DMATMUL_TGPU_() __attribute__((weak,alias("dmatmul_tgpu_")));
void dmatmul_tgpu__() __attribute__((weak,alias("dmatmul_tgpu_")));
void DMATMUL_TGPU__() __attribute__((weak,alias("dmatmul_tgpu_")));
/**********************************************************************
 * CMATMUL.
 **********************************************************************/
int cmatmul_tgpu_(const char *transa, const char *transb, const int *m,
                 const int *n, const int *k, const cuComplex *alpha,
                 cuComplex *A, const int *lda, 
                 cuComplex *B, const int *ldb, const cuComplex *beta,
                 cuComplex *C, const int *ldc){

  int rc = RC_SUCCESS;
  
  cuComplex *devPtrA = NULL;
  cuComplex *devPtrB = NULL;
  cuComplex *devPtrC = NULL;

  //input matrices
  rc = cublasAlloc(*m*(*k), size_of_complex, (void**)&devPtrA);
  if (rc != RC_SUCCESS ) {rc = simple_cudblas_stat_return_c(rc);}
  rc = cublasAlloc(*k*(*n), size_of_complex, (void**)&devPtrB);
  if (rc != RC_SUCCESS ) {rc = simple_cudblas_stat_return_c(rc);}

  rc = cublasSetMatrix (*m,*k,size_of_complex,A,*lda,devPtrA,*lda);
  if (rc != RC_SUCCESS ) {rc = simple_cudblas_stat_return_c(rc);}
  rc = cublasSetMatrix (*k,*n,size_of_complex,B,*ldb,devPtrB,*ldb);
  if (rc != RC_SUCCESS ) {rc = simple_cudblas_stat_return_c(rc);}
    
  /* !return matrix */
  rc = cublasAlloc(*m*(*n), size_of_complex, (void**)&devPtrC);
  if (rc != RC_SUCCESS ) {rc = simple_cudblas_stat_return_c(rc);}

  cublasCgemm(*transa,*transb,*m,*n,*k,
              *alpha,
              devPtrA,*lda,
              devPtrB,*ldb,
              *beta,
              devPtrC,*ldc);

  /* retrieving back from GPU */
  rc = cublasGetMatrix( *m, *n, size_of_complex, devPtrC, *lda, C, *ldc);
  if (rc != RC_SUCCESS ) {rc = simple_cudblas_stat_return_c(rc);}

  /* freeing ressources on the device */
  rc = cublasFree(devPtrA);
  rc = cublasFree(devPtrB);
  rc = cublasFree(devPtrC);

  return rc;
}
void CMATMUL_TGPU_() __attribute__((weak,alias("cmatmul_tgpu_")));
void cmatmul_tgpu__() __attribute__((weak,alias("cmatmul_tgpu_")));
void CMATMUL_TGPU__() __attribute__((weak,alias("cmatmul_tgpu_")));
/**********************************************************************
 * ZMATMUL.
 **********************************************************************/
int zmatmul_tgpu_(const char *transa, const char *transb, const int *m,
                 const int *n, const int *k, const cuDoubleComplex *alpha,
                 cuDoubleComplex *A, const int *lda, 
                 cuDoubleComplex *B, const int *ldb,
                  const cuDoubleComplex *beta,
                 cuDoubleComplex *C, const int *ldc){

  int rc = RC_SUCCESS;
  
  cuDoubleComplex *devPtrA = NULL;
  cuDoubleComplex *devPtrB = NULL;
  cuDoubleComplex *devPtrC = NULL;

  //input matrices
  rc = cublasAlloc(*m*(*k), size_of_double_complex, (void**)&devPtrA);
  if (rc != RC_SUCCESS ) {rc = simple_cudblas_stat_return_c(rc);}
  rc = cublasAlloc(*k*(*n), size_of_double_complex, (void**)&devPtrB);
  if (rc != RC_SUCCESS ) {rc = simple_cudblas_stat_return_c(rc);}

  rc = cublasSetMatrix (*m,*k,size_of_double_complex,A,*lda,devPtrA,*lda);
  if (rc != RC_SUCCESS ) {rc = simple_cudblas_stat_return_c(rc);}
  rc = cublasSetMatrix (*k,*n,size_of_double_complex,B,*ldb,devPtrB,*ldb);
  if (rc != RC_SUCCESS ) {rc = simple_cudblas_stat_return_c(rc);}
    
  /* !return matrix */
  rc = cublasAlloc(*m*(*n), size_of_double_complex, (void**)&devPtrC);
  if (rc != RC_SUCCESS ) {rc = simple_cudblas_stat_return_c(rc);}

  cublasZgemm(*transa,*transb,*m,*n,*k,
              *alpha,
              devPtrA,*lda,
              devPtrB,*ldb,
              *beta,
              devPtrC,*ldc);

  /* retrieving back from GPU */
  rc = cublasGetMatrix( *m, *n, size_of_double_complex, devPtrC, *lda, C, *ldc);
  if (rc != RC_SUCCESS ) {rc = simple_cudblas_stat_return_c(rc);}

  /* freeing ressources on the device */
  rc = cublasFree(devPtrA);
  rc = cublasFree(devPtrB);
  rc = cublasFree(devPtrC);

  return rc;
}
void ZMATMUL_TGPU_() __attribute__((weak,alias("zmatmul_tgpu_")));
void zmatmul_tgpu__() __attribute__((weak,alias("zmatmul_tgpu_")));
void ZMATMUL_TGPU__() __attribute__((weak,alias("zmatmul_tgpu_")));

#endif /* CUDA */
