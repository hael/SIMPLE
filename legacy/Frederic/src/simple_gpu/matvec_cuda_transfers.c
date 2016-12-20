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
int smatvec_tgpu_(const char *trans, const int *m,
                  const int *n, const float *alpha,
                  float *A, const int *lda, 
                  float *X, const int *incx, const float *beta,
                  float *Y, const int *incy){

  int rc = RC_SUCCESS;
  
  float *devPtrA = NULL;
  float *devPtrX = NULL;
  float *devPtrY = NULL;
  
  //input matrices
  rc = cublasAlloc(*m*(*n), size_of_float, (void**)&devPtrA);
  if (rc != RC_SUCCESS ) {rc = simple_cudblas_stat_return_c(rc);}
  rc = cublasAlloc((*n), size_of_float, (void**)&devPtrX);
  if (rc != RC_SUCCESS ) {rc = simple_cudblas_stat_return_c(rc);}

  rc = cublasSetMatrix (*m,*n,size_of_float,A,*lda,devPtrA,*lda);
  if (rc != RC_SUCCESS ) {rc = simple_cudblas_stat_return_c(rc);}
  rc = cublasSetVector (*n,size_of_float,X,*incx,devPtrX,*incx);
  if (rc != RC_SUCCESS ) {rc = simple_cudblas_stat_return_c(rc);}
    
  /* !return matrix */
  rc = cublasAlloc((*n), size_of_float, (void**)&devPtrY);
  if (rc != RC_SUCCESS ) {rc = simple_cudblas_stat_return_c(rc);}

  cublasSgemv(*trans,*m,*n,
              *alpha,
              devPtrA,*lda,
              devPtrX,*incx,
              *beta,
              devPtrY,*incy);

  /* retrieving back from GPU */
  rc = cublasGetVector(*n, size_of_float, devPtrY, *incy, Y, *incy);
  if (rc != RC_SUCCESS ) {rc = simple_cudblas_stat_return_c(rc);}

  /* freeing ressources on the device */
  rc = cublasFree(devPtrA);
  rc = cublasFree(devPtrX);
  rc = cublasFree(devPtrY);

  return rc;
}
void SMATVEC_TGPU_() __attribute__((weak,alias("smatvec_tgpu_")));
void smatvec_tgpu__() __attribute__((weak,alias("smatvec_tgpu_")));
void SMATVEC_TGPU__() __attribute__((weak,alias("smatvec_tgpu_")));
/**********************************************************************
 * DMATVEC.
 **********************************************************************/
int dmatvec_tgpu_(const char *trans, const int *m,
                  const int *n, const double *alpha,
                  double *A, const int *lda, 
                  double *X, const int *incx, const double *beta,
                  double *Y, const int *incy){

  int rc = RC_SUCCESS;
  
  double *devPtrA = NULL;
  double *devPtrX = NULL;
  double *devPtrY = NULL;

  //input matrices
  rc = cublasAlloc(*m*(*n), size_of_double, (void**)&devPtrA);
  if (rc != RC_SUCCESS ) {rc = simple_cudblas_stat_return_c(rc);}
  rc = cublasAlloc((*n), size_of_double, (void**)&devPtrX);
  if (rc != RC_SUCCESS ) {rc = simple_cudblas_stat_return_c(rc);}

  rc = cublasSetMatrix (*m,*n,size_of_double,A,*lda,devPtrA,*lda);
  if (rc != RC_SUCCESS ) {rc = simple_cudblas_stat_return_c(rc);}
  rc = cublasSetVector (*n,size_of_double,X,*incx,devPtrX,*incx);
  if (rc != RC_SUCCESS ) {rc = simple_cudblas_stat_return_c(rc);}
    
  /* !return matrix */
  rc = cublasAlloc((*n), size_of_double, (void**)&devPtrY);
  if (rc != RC_SUCCESS ) {rc = simple_cudblas_stat_return_c(rc);}

  cublasDgemv(*trans,*m,*n,
              *alpha,
              devPtrA,*lda,
              devPtrX,*incx,
              *beta,
              devPtrY,*incy);

  /* retrieving back from GPU */
  rc = cublasGetVector(*n, size_of_double, devPtrY, *incy, Y, *incy);
  if (rc != RC_SUCCESS ) {rc = simple_cudblas_stat_return_c(rc);}

  /* freeing ressources on the device */
  rc = cublasFree(devPtrA);
  rc = cublasFree(devPtrX);
  rc = cublasFree(devPtrY);

  return rc;
}
void DMATVEC_TGPU_() __attribute__((weak,alias("dmatvec_tgpu_")));
void dmatvec_tgpu__() __attribute__((weak,alias("dmatvec_tgpu_")));
void DMATVEC_TGPU__() __attribute__((weak,alias("dmatvec_tgpu_")));
/**********************************************************************
 * CMATVEC.
 **********************************************************************/
int cmatvec_tgpu_(const char *trans, const int *m,
                  const int *n, const cuComplex *alpha,
                  cuComplex *A, const int *lda, 
                  cuComplex *X, const int *incx, const cuComplex *beta,
                  cuComplex *Y, const int *incy){

  int rc = RC_SUCCESS;
  
  cuComplex *devPtrA = NULL;
  cuComplex *devPtrX = NULL;
  cuComplex *devPtrY = NULL;

  //input matrices
  rc = cublasAlloc(*m*(*n), size_of_complex, (void**)&devPtrA);
  if (rc != RC_SUCCESS ) {rc = simple_cudblas_stat_return_c(rc);}
  rc = cublasAlloc((*n), size_of_complex, (void**)&devPtrX);
  if (rc != RC_SUCCESS ) {rc = simple_cudblas_stat_return_c(rc);}

  rc = cublasSetMatrix (*m,*n,size_of_complex,A,*lda,devPtrA,*lda);
  if (rc != RC_SUCCESS ) {rc = simple_cudblas_stat_return_c(rc);}
  rc = cublasSetVector (*n,size_of_complex,X,*incx,devPtrX,*incx);
  if (rc != RC_SUCCESS ) {rc = simple_cudblas_stat_return_c(rc);}
    
  /* !return matrix */
  rc = cublasAlloc((*n), size_of_complex, (void**)&devPtrY);
  if (rc != RC_SUCCESS ) {rc = simple_cudblas_stat_return_c(rc);}

  cublasCgemv(*trans,*m,*n,
              *alpha,
              devPtrA,*lda,
              devPtrX,*incx,
              *beta,
              devPtrY,*incy);

  /* retrieving back from GPU */
  rc = cublasGetVector(*n, size_of_complex, devPtrY, *incy, Y, *incy);
  if (rc != RC_SUCCESS ) {rc = simple_cudblas_stat_return_c(rc);}

  /* freeing ressources on the device */
  rc = cublasFree(devPtrA);
  rc = cublasFree(devPtrX);
  rc = cublasFree(devPtrY);

  return rc;
}
void CMATVEC_TGPU_() __attribute__((weak,alias("cmatvec_tgpu_")));
void cmatvec_tgpu__() __attribute__((weak,alias("cmatvec_tgpu_")));
void CMATVEC_TGPU__() __attribute__((weak,alias("cmatvec_tgpu_")));
/**********************************************************************
 * ZMATVEC.
 **********************************************************************/
int zmatvec_tgpu_(const char *trans, const int *m,
                  const int *n, const cuDoubleComplex *alpha,
                  cuDoubleComplex *A, const int *lda, 
                  cuDoubleComplex *X, const int *incx,
                  const cuDoubleComplex *beta,
                  cuDoubleComplex *Y, const int *incy){

  int rc = RC_SUCCESS;
  
  cuDoubleComplex *devPtrA = NULL;
  cuDoubleComplex *devPtrX = NULL;
  cuDoubleComplex *devPtrY = NULL;

  //input matrices
  rc = cublasAlloc(*m*(*n), size_of_double_complex, (void**)&devPtrA);
  if (rc != RC_SUCCESS ) {rc = simple_cudblas_stat_return_c(rc);}
  rc = cublasAlloc((*n), size_of_double_complex, (void**)&devPtrX);
  if (rc != RC_SUCCESS ) {rc = simple_cudblas_stat_return_c(rc);}

  rc = cublasSetMatrix (*m,*n,size_of_double_complex,A,*lda,devPtrA,*lda);
  if (rc != RC_SUCCESS ) {rc = simple_cudblas_stat_return_c(rc);}
  rc = cublasSetVector (*n,size_of_double_complex,X,*incx,devPtrX,*incx);
  if (rc != RC_SUCCESS ) {rc = simple_cudblas_stat_return_c(rc);}
    
  /* !return matrix */
  rc = cublasAlloc((*n), size_of_double_complex, (void**)&devPtrY);
  if (rc != RC_SUCCESS ) {rc = simple_cudblas_stat_return_c(rc);}

  cublasZgemv(*trans,*m,*n,
              *alpha,
              devPtrA,*lda,
              devPtrX,*incx,
              *beta,
              devPtrY,*incy);

  /* retrieving back from GPU */
  rc = cublasGetVector(*n, size_of_double_complex, devPtrY, *incy, Y, *incy);
  if (rc != RC_SUCCESS ) {rc = simple_cudblas_stat_return_c(rc);}

  /* freeing ressources on the device */
  rc = cublasFree(devPtrA);
  rc = cublasFree(devPtrX);
  rc = cublasFree(devPtrY);

  return rc;
}
void ZMATVEC_TGPU_() __attribute__((weak,alias("zmatvec_tgpu_")));
void zmatvec_tgpu__() __attribute__((weak,alias("zmatvec_tgpu_")));
void ZMATVEC_TGPU__() __attribute__((weak,alias("zmatvec_tgpu_")));

#endif /* CUDA */
