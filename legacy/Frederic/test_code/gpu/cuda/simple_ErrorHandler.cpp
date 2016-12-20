/*******************************************************************************
 * SIMPLE error Handler for the library                                        *
 *******************************************************************************
 *
 *   -- SIMPLE addon
 *      Author: Frederic Bonnet, Date: 17th Mar 2016
 *      Monash University
 *      March 2016
 *
 *      error hadler code for the library.
 *
 * @precisions normal z -> s d c
 */
#include "simple.h"

#if defined (CUDA) /*preprossing for the CUDA environment */

extern "C" int
simple_cudblas_stat_return_c(int err) {

  int rc = RC_SUCCESS;
    
#if defined (CUDA)

  switch(err) {

    case CUBLAS_STATUS_SUCCESS:
      printf("Error=%i : CUBLAS_STATUS_SUCCESS\n",err); rc = RC_FAIL;
      break;
    case  CUBLAS_STATUS_NOT_INITIALIZED:
      printf("Error=%i : CUBLAS_STATUS_NOT_INITIALIZED\n",err);
      printf("Error=%i : CUBLAS_STATUS_LICENSE_ERROR\n",err); rc = RC_FAIL;
      break;
    case  CUBLAS_STATUS_ALLOC_FAILED:
      printf("Error: device memory for devPtrA allocation failed\n");
      printf("Error=%i : CUBLAS_STATUS_ALLOC_FAILED\n",err); rc = RC_FAIL;
      break;
    case  CUBLAS_STATUS_INVALID_VALUE:
      printf("Error=%i : CUBLAS_STATUS_INVALID_VALUE\n",err); rc = RC_FAIL;
      break;
    case  CUBLAS_STATUS_ARCH_MISMATCH:
      printf("Error=%i : CUBLAS_STATUS_ARCH_MISMATCH\n",err); rc = RC_FAIL;
      break;
    case  CUBLAS_STATUS_MAPPING_ERROR:
      printf("Error: Data upload on the GPU failed\n"); rc = RC_FAIL;
      break;
      printf("Error=%i : CUBLAS_STATUS_MAPPING_ERROR\n",err); rc = RC_FAIL;
      break;
    case  CUBLAS_STATUS_EXECUTION_FAILED:
      printf("Error=%i : CUBLAS_STATUS_EXECUTION_FAILED\n",err); rc = RC_FAIL;
      break;
    case  CUBLAS_STATUS_INTERNAL_ERROR:
      printf("Error=%i : CUBLAS_STATUS_INTERNAL_ERROR\n",err); rc = RC_FAIL;
      break;
    case  CUBLAS_STATUS_NOT_SUPPORTED:
      printf("Error=%i : CUBLAS_STATUS_NOT_SUPPORTED\n",err); rc = RC_FAIL;
      break;
    default:
      printf("Error=%i : UNKOWN ERROR, check the argument list in code\n",err);
      rc = RC_FAIL;
      break;
    }

#else
    rc = -1;
    printf("**************************WARNING******************************\n");
    printf("You need to compile with -DCUDA to acces the CUDA environment  \n");
    printf("computation using GPU                                          \n");
    printf("***************************************************************\n");
#endif 

    return rc;
  }



#endif /* CUDA */
