/*******************************************************************************
 *     Author: Frederic D.R. Bonnet date: 08th of April 2011.11:19am
 *
 * Name:
 * cpPlatformGPU.cpp - basic definitions used in all modules.
 *
 * Description:
 * classes to set and get the cpPlatform on the GPU in OpenCL
 * routine to get the NVIDIA platform id and set to the class
 * set_cpPlatformGPU on the GPU in OpenCL 
 *******************************************************************************
 */
#include <iostream>

#if defined (OPENCL) /*preprossing for the OPENCL environment */

#include <openCL_GPU_funcDec.h>

cpPlatformGPU::cpPlatformGPU(cl_platform_id Platform)
{
  cpPlatform = Platform;
}
// setter
void cpPlatformGPU::set_cpPlatformGPU(cl_platform_id Platform)
{
  cpPlatform = Platform;
}
//getter
cl_platform_id cpPlatformGPU::get_cpPlatformGPU()
{
  return cpPlatform;
}
// the destructor
cpPlatformGPU::~cpPlatformGPU()
{
  std::cout << "Object cpPlatformGPU has been destroyed" << std::endl;
}

#endif /* OPENCL */
