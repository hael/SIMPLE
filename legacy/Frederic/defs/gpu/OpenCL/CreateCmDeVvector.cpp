/*******************************************************************************
 *     Author: Frederic D.R. Bonnet date: 12th of April 2011.13:01am
 *
 * Name:
 * CreateCmDeVvector.cpp - basic definitions used in all modules.
 *
 * Description:
 * routine that creates, vectors and sets them to an object to be
 * passed in the routine that send data to the GPU. Used to 
 * create and launche the kernel on the GPU in OpenCL
 *******************************************************************************
 */

#if defined (OPENCL) /*preprossing for the OPENCL environment */

#include <global.h>
#include <openCL_GPU_funcDec.h>

//(cxGPUContext_2, CL_MEM_READ_WRITE|CL_MEM_COPY_HOST_PTR, sizeof(cl_float)*N, vector1, &ciErrNum);

//extern "C" void CreateCmDeVvector(float);

void CreateCmDeVvector(float *vector )
{
  cl_mem cmDeVvector;

  cl_context cxGPUContext_2 = context -> get_cxGPUContextGPU(); //getting the context
  unsigned int m = worksizeGPU -> get_szVectorGPU(); //getting the size of the vector

  cmDeVvector  = clCreateBuffer(cxGPUContext_2, CL_MEM_READ_WRITE|CL_MEM_COPY_HOST_PTR, sizeof(cl_float)*m, vector, &ciErrNum);

  //now setting this to the object to be retieved later in routin_...()

  DeVvectorGPU -> set_szVectorGPU(cmDeVvector);//return cmDeVvector;

}
extern "C" void createcmdevvector_gpu_(float *) __attribute__((weak,alias("CreateCmDeVvector")));

#endif /* OPENCL */
