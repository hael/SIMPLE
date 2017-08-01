/*******************************************************************************
 *     Author: Frederic D.R. Bonnet date: 12th of April 2011. 12:19am
 *
 * Name:
 * cxGPUContextGPU.cpp - basic definitions used in all modules.
 *
 * Description:
 * class to set and get the context on the GPU in OpenCL
 *******************************************************************************
 */
#include <iostream>

#if defined (OPENCL) /*preprossing for the OPENCL environment */

#include <openCL_GPU_funcDec.h>

cxGPUContextGPU::cxGPUContextGPU(cl_context GPUContext)
{
  cxGPUContext = GPUContext;
}
// setter
void cxGPUContextGPU::set_cxGPUContextGPU(cl_context GPUContext)
{
  cxGPUContext = GPUContext;
}
// getter
cl_context cxGPUContextGPU::get_cxGPUContextGPU()
{
  return cxGPUContext;
}
// the destructor
cxGPUContextGPU::~cxGPUContextGPU()
{
  std::cout << "Object cpPlatformGPU has been destroyed" << std::endl;
}

#endif /* OPENCL */
