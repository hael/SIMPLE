/*******************************************************************************
 *     Author: Frederic D.R. Bonnet date: 12th of April 2011. 12:19am
 *
 * Name:
 * kernel_saX.cc - basic definitions used in all modules.
 *
 * Description:
 * OpenCL kernel saX to multiply a vector by a constant alpha
 *******************************************************************************
 */

#if defined (OPENCL) /*preprossing for the OPENCL environment */

__kernel void saX(
		  unsigned int N,  
		  float ALPHA,  
		  __global float* X,
		  unsigned int INCX
)
{
    // get index into global data array
  unsigned int tid = get_global_id(0);

  if (tid < N)
    X[tid*INCX] = ALPHA*X[tid*INCX];
}

#endif /* OPENCL */
