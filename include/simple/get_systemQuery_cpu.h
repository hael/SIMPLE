/*
 *   -- SIMPLE addon
 *      Author: Frederic Bonnet, Date: 21st May 2015
 *
 *      May 2015
 *
 *   Header files for cpp code to obtain information from the GPU cards
 *
 * @precisions normal z -> s d c
 */

/* The Simple header */
#include "simple.h"

#ifndef _GET_SYSTEMQUERY_CPU_H_
#define _GET_SYSTEMQUERY_CPU_H_
#define PRECISION_z

#ifdef __cplusplus
extern "C" {
#endif

/* ////////////////////////////////////////////////////////////////////////////
   -- routines used to interface back to fortran
   *
   *  FORTRAN API -
   * 
   */

  /*error and warning handlers methods */
  int get_error_cpu();
  int get_warning_message_linux_cpu();
  /* quering handlers methods */
  int get_CPU_cores(int *nCPU_cores, systemDetails_t *hstD);
  int get_memorySize_host(long long int *memSize, systemDetails_t *hstD);
  int get_available_memory_host(long long int *availMem, systemDetails_t *hstD);
  float convert_memory(long long int mem);

#ifdef __cplusplus
}
#endif

#undef PRECISION_z
#endif /* _GET_SYSTEMQUERY_CPU_H_ */
