/*
 *   -- MAGMA addon
 *      Author: Frederic Bonnet, Date: 22nd Apr 2015
 *      Monash University
 *      Apr 2015
 *
 *   code to get the identity in CUDA.
 *
 * @precisions normal z -> s d c
 */

/*preprossing for the CUDA and MAGMA environment */
#if defined (CUDA) && defined (MAGMA) 

//MAGMA include files from {path}/magma-1.6.1/include
#include "magma.h"

#include "gen_polar_coords_gpu.h"

//Magama definition of device pointer for a complex pointer
#define magma_zdevptr(ptr_) ((magmaDoubleComplex*)(uintptr_t)(*(ptr_)))
#define magma_ddevptr(ptr_) ((double*)(uintptr_t)(*(ptr_)))

/* 
 * typedef comming from fortran.h file provided in $CUDADIR/src directory
 * it will probably change with future release of cublas when they will use 64bits address
 */
typedef size_t devptr_t;

#ifdef __cplusplus
extern "C" {
#endif

/*****************************************************************************
 *  FORTRAN API - math functions (simple interface)
 *
*/

  void GEN_POLAR_COORDS_CUDA_GPU( devptr_t *dkfromto, int *ring2,
				  devptr_t *dcoords, devptr_t *dangtab ) {
    int rc; /* return code */ 
    float *kfromto = (float *)(*dkfromto);
    float *coords  = (float *)(*dcoords);
    float *angtab  = (float *)(*dangtab);
    rc = gen_polar_coords_cuda_gpu(kfromto, *ring2, coords, angtab);
  }

/* //////////////////////////////////////////////////////////////////////////// 
 -- MAGMA function definitions / Data on GPU
*/

  extern "C" void gen_polar_coords_cuda_gpu_( devptr_t, int, devptr_t, devptr_t) __attribute__((weak,alias("GEN_POLAR_COORDS_CUDA_GPU")));

#ifdef __cplusplus
}
#endif

#endif /* (CUDA) && (MAGMA) */
