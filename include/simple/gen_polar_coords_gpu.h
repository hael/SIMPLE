/*
 *   -- DZgemm kernel addon for 
 *      Author: Frederic Bonnet, Date: 1st September 2011
 *      CEA Grenoble
 *      September 2011
 *
 * @precisions normal z -> s d c
 */

#ifndef _GEN_POLAR_COORDS_GPU_H_
#define _GEN_POLAR_COORDS_GPU_H_

#ifdef __cplusplus
extern "C" {
#endif

/* //////////////////////////////////////////////////////////////////////////// 
 -- gen_polar_coords_gpu kernel function definitions / Data on GPU
*/
#if defined (CUDA) /*preprossing for the OPENCL environment */

  int gen_polar_coords_cuda_gpu(float *kfromto, int ring2, float *coords, float *angtab);

#endif /* CUDA */

#ifdef __cplusplus
}
#endif

#endif /* _GEN_POLAR_COORDS_GPU_H_ */
