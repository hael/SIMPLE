/*
 *   -- SIMPLE addon
 *      Author: Frederic Bonnet, Date: 06th April 2015
 *
 *      April 2015
 *
 *      c code for the MACOSX code this wrapps around the the cu for LINUX
 *      files methods for the cpp files to extract the information from the GPU
 *      devices.
 *
 * @precisions normal z -> s d c
 */

/* The Simple header */
#include "simple.h"

/* determine if peer to peer is allowed */
int print_s_devd_struct_c_(deviceDetails_t *devD) {
  int rc = RC_SUCCESS;
#if defined (CUDA) /*preprossing for the CUDA environment */
  rc = print_s_devD_struct(devD);
  if (rc == RC_FAIL) { rc = get_error_c(); }
#else
  rc = get_warning_message_cuda_c_();
  return rc;
#endif /* CUDA */
  return rc;
}
void PRINT_S_DEVD_STRUCT_C_() __attribute__((weak,alias("print_s_devd_struct_c_")));
void print_s_devd_struct_c__() __attribute__((weak,alias("print_s_devd_struct_c_")));
void PRINT_S_DEVD_STRUCT_C__() __attribute__((weak,alias("print_s_devd_struct_c_")));
/* getting the warning message from the cuda */
int get_warning_message_cuda_c_() {
  int rc = RC_SUCCESS;

  printf("***************************WARNING*****************************\n");
  printf("You need to compile with -DCUDA to acces the CUDA environment  \n");
  printf("computation using GPU ---> modify simple_user_input.pm         \n");
  printf("***************************************************************\n");
  printf("\n");
  printf("Line %i:%s ---> %s\n",__LINE__,__FILE__,__FUNCTION__);
  printf("\n");
  rc = RC_FAIL;

  return rc;
}
void GET_WARNING_MESSAGE_CUDA_C_() __attribute__((weak,alias("get_warning_message_cuda_c_")));
void get_warning_message_cuda_c__() __attribute__((weak,alias("get_warning_message_cuda_c_")));
void GET_WARNING_MESSAGE_CUDA_C__() __attribute__((weak,alias("get_warning_message_cuda_c_")));
