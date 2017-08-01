/*******************************************************************************
 *     Author: Frederic D.R. Bonnet date: 08th of April 2011.11:19am
 *
 * Name:
 * ciDeviceCountGPU.cpp - basic definitions used in all modules.
 *
 * Description:
 * classes to set and get the ciDeviceCount on the GPU in OpenCL
 *******************************************************************************
 */
#include <iostream>

#if defined (OPENCL) /*preprossing for the OPENCL environment */

#include <openCL_GPU_funcDec.h>

ciDeviceCountGPU::ciDeviceCountGPU(cl_uint DeviceCount )
{
  ciDeviceCount = DeviceCount;
}
// setter
void ciDeviceCountGPU::set_ciDeviceCountGPU(cl_uint DeviceCount)
{
  ciDeviceCount = DeviceCount;
}
// getter
cl_uint ciDeviceCountGPU::get_ciDeviceCountGPU()
{
  return ciDeviceCount;
}
// the destructor
ciDeviceCountGPU::~ciDeviceCountGPU()
{
  std::cout << "Object cpPlatformGPU has been destroyed" << std::endl;
}

#endif /* OPENCL */
