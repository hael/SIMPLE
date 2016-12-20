/*******************************************************************************
 *     Author: Frederic D.R. Bonnet date: 12th of April 2011.11:48am
 *
 * Name:
 * cdDevicesGPU.cpp - basic definitions used in all modules.
 *
 * Description:
 * pclasses to set and get the cdDevices* on the GPU in OpenCL
 *******************************************************************************
 */
#include <iostream>

#if defined (OPENCL) /*preprossing for the OPENCL environment */

#include <openCL_GPU_funcDec.h>

cdDevicesGPU::cdDevicesGPU(cl_device_id* Devices )
{
  cdDevices = Devices;
}
// setter
void cdDevicesGPU::set_cdDevicesGPU(cl_device_id* Devices )
{
  cdDevices = Devices;
}
// getter
cl_device_id* cdDevicesGPU::get_cdDevicesGPU()
{
  return cdDevices;
}
// the destructor
cdDevicesGPU::~cdDevicesGPU()
{
  std::cout << "Object cdDevicesGPU has been destroyed" << std::endl;
}

#endif /* OPENCL */
