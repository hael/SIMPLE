
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
#include "cudnn.h"    /* header file for the cudnn */

#include "simple.h"
#include "simple_cudnn_fortran_common.h"
#include "simple_cudnn_fortran.h"

//Customizable allocation wrappeprs

/*initialisation of the cuDNN environment */
int SIMPLE_CUDNN_CREATE(cudnnHandle_t *handle){
  int rc = RC_SUCCESS;
  cudnnStatus_t rc_cudnn;
  rc_cudnn = cudnnCreate(handle);
  if (rc_cudnn != CUDNN_STATUS_SUCCESS) {
    printf("error string: %s\n",cudnnGetErrorString(rc_cudnn));
  }

  return rc_cudnn;
}
/* Termination of the cuDNN environment */
int SIMPLE_CUDNN_DESTROY(cudnnHandle_t *handle){
  int rc = RC_SUCCESS;
  cudnnStatus_t rc_cudnn;
  rc_cudnn = cudnnDestroy(*handle);
  if (rc_cudnn != CUDNN_STATUS_SUCCESS) {
    printf("error string: %s\n",cudnnGetErrorString(rc_cudnn));
  }
  return rc_cudnn;
}

#endif /* CUDA */
