/*
 *   -- SIMPLE addon
 *      Author: Frederic Bonnet, Date: 24th April 2015
 *      Monash University
 *      September 2015
 *
 *   code to get multiply a double matrix with a double complex matrix on GPU.
 *
 * @precisions normal z -> s d c
 */
//#include <ctype.h>
//#include <stdio.h>
//#include <string.h>
//#include <stddef.h>
//#include <stdlib.h>
//#include <stdlib.h>
//#if defined(__GNUC__)
//#include <stdint.h>
//#endif /* __GNUC__ */

#include "cublas.h"   /* CUBLAS public header file  */
#include <cuda.h>
#include <cuda_runtime.h>

#include "simple.h"

#include "gen_polar_coords_gpu.h"

#include "get_deviceQuery_gpu.h"
#include "get_systemQuery_cpu.h"

#if defined (CUDA) /*preprossing for the OPENCL environment */

/* ////////////////////////////////////////////////////////////////////////////
 *   -- The Kenels for the GPU
 *
 *      Author: Frederic Bonnet, Date: 1st September 2011
 *      CEA Grenoble
 *      September 2011
 *
 *      from ./cuda/include/cuComplex.h
 *      making use of:
 *      cuCadd(cuDoubleComplex x,cuDoubleComplex y)
 *      cuCmul(cuDoubleComplex x,cuDoubleComplex y)
*/
typedef struct {
  int ring2;
  float *A1;
  float *B1;
  float *C1;
} cublasDZgemm_gpu_params;

#define BLOCK_SIZE 16
/*///////////////////////////////////////////////////////////////////////////////////////////////////
 *   -- The Kenels for the GPU
 *      if ( toupper(transa[0]) == 'N' && toupper(transb[0]) == 'N' )
 *
 *      Author: Frederic Bonnet, Date: 1st September 2011
 *      CEA Grenoble
 *      September 2011
 *
 */
__global__ void matMull_DZgemm_NN_gpu(cublasDZgemm_gpu_params IDparams)
{
  /*  
  int i = threadIdx.x + blockIdx.x * blockDim.x;
  int j = threadIdx.y + blockIdx.y * blockDim.y;
  int k;

  cuDoubleComplex zero = make_cuDoubleComplex ( 0.0, 0.0 );
  
  int              r2 = IDparams.ring2;
  float      *IDmatA1 = IDparams.A1;
  float      *IDmatB1 = IDparams.B1;
  float      *IDmatC1 = IDparams.C1;
  */
}

/*////////////////////////////////////////////////////////////////////////////////////////////////////
 *   -- The wrappers for the Kenels for the GPU
 *      the gen_polar_coords_gpu
 *
 *      Author: Frederic Bonnet, Date: 24th April 2015
 *      Monash University
 *      April 2015
 */
extern "C" int
gen_polar_coords_cuda_gpu(float *kfromto, int ring2, float *coords, float *angtab) {

  int rc = 0; //return code

  //cuDoubleComplex *alpha_gpu; //uploaded value to GPU
  //cuDoubleComplex temp_alpha;
  cudaError_t err;

  //start of the execution commands.

  //dim3 grid( n / BLOCK_SIZE, n / BLOCK_SIZE, 1);
  //dim3 blocks(BLOCK_SIZE, BLOCK_SIZE, 1);

  size_t size_const = sizeof(cuDoubleComplex);

  //allocation GPU
  //cudaMalloc((void **) &alpha_gpu,size_const);
  //Copy alpha onto GPU
  //cudaMemcpy(alpha_gpu,&alpha,size_const,cudaMemcpyHostToDevice);

  //printf("sizeof(cuDoubleComplex: %i, at Line %i in file %s %s\n",size_const,__LINE__,__FILE__,__FUNCTION__);
  //cublasDZgemm_gpu_params params = {ring2, kfromto, coords, angtab};
  //matMull_DZgemm_NN_gpu<<<grid,blocks>>>(params);
  err=cudaThreadSynchronize();
  if(err != cudaSuccess)
    {
      printf("Error at execcution time:\n\t\t%s\n",cudaGetErrorString(err));
    }

  //getting alpha back to sync
  //cudaMemcpy(&temp_alpha,alpha_gpu,size_const,cudaMemcpyDeviceToHost);
  /*
  printf("ring2 = %i, at Line %i in file %s %s\n", ring2,__LINE__,__FILE__,__FUNCTION__);
  printf("kfromto = %d, at Line %i in file %s %s\n", kfromto,__LINE__,__FILE__,__FUNCTION__);
  printf("coords = %d, at Line %i in file %s %s\n", coords,__LINE__,__FILE__,__FUNCTION__);
  printf("angtab = %d, at Line %i in file %s %s\n", angtab,__LINE__,__FILE__,__FUNCTION__);
  */
  //freeing the memory on the GPU.
  //cudaFree(alpha_gpu);
  
  return rc;

} /* End of cublasDZgemm_gpu */

#endif /* CUDA */

