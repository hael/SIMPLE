
/*
 * This file contains example Fortran bindings for the CUBLAS library, These
 * bindings have been tested with Intel Fortran 9.0 on 32-bit and 64-bit 
 * Windows, and with g77 3.4.5 on 32-bit and 64-bit Linux. They will likely
 * have to be adjusted for other Fortran compilers and platforms.
 */

#include <ctype.h>
#include <stdio.h>
#include <string.h>
#include <stddef.h>
#include <stdlib.h>
#if defined(__GNUC__)
#include <stdint.h>
#endif /* __GNUC__ */

#if defined (CUDA) /* prepropressing for the cuda environment */

#include "cublas.h"   /* CUBLAS public header file  */

#include "simple_fortran_common.h"
#include "simple_fortran.h"

//Customizable allocation wrappeprs
int SIMPLE_CUBLAS_ALLOC(const int *n, const int *elemSize, devptr_t *devicePtr)
{    
  void *tPtr;
  int rc;
  rc = (int)cublasAlloc (*n, *elemSize, &tPtr);
  *devicePtr = (devptr_t)tPtr;
  return rc;
}

/* memory synchronisation */
int SIMPLE_CUBLAS_FREE(const devptr_t *devicePtr) {
  cudaEvent_t start_event, stop_event;
  void *tPtr;

  cudaEventCreate(&start_event);
  cudaEventCreate(&stop_event);

  tPtr = (void *)(*devicePtr);
  cudaEventRecord(start_event,0);
  cublasFree (tPtr);
  cudaEventRecord(stop_event,0);

  cudaEventSynchronize(stop_event);

  return (int)cuCtxSynchronize();
}

#endif /* CUDA */
