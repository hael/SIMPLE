/*
 *   -- MAGMA addon
 *      Author: Frederic Bonnet, Date: 1st Apr 2015
 *      Monash University
 *      Apr 2015
 *
 * @precisions normal z -> s d c
 */

#include "simple.h"

#ifndef _SIMPLE_MAGMA_COMMON_
#define _SIMPLE_MAGMA_COMMON_

#ifdef __cplusplus
extern "C" {
#endif

/* //////////////////////////////////////////////////////////////////////////// 
 -- MAGMA function definitions / Some getters
*/

// ==== Definition of blocking sizes for Tesla ===============================
#if (GPUSHMEM < 200)
/* ////////////////////////////////////////////////////////////////////////////
   -- Return nb for zgetri based on n
*/
  int magma_get_getri_nb_gpu(int);
// ==== End Definition of blocking sizes for Tesla ===========================
#else
// ====     Definition of blocking sizes for Fermi ===========================
  int magma_get_getri_nb_gpu(int);
// ==== End Definition of blocking sizes for Fermi ===========================
#endif

#ifdef __cplusplus
}
#endif

#endif /* _SIMPLE_MAGMA_COMMON_ */
