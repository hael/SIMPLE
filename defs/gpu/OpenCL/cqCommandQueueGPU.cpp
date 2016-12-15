/*******************************************************************************
 *     Author: Frederic D.R. Bonnet date: 12th of April 2011.15:01am
 *
 * Name:
 * cqCommandQueueGPU.cpp - basic definitions used in all modules.
 *
 * Description:
 * classes to set and get the cqCommandQueueGPU on the GPU in OpenCL
 *******************************************************************************
 */
#include <iostream>

#if defined (OPENCL) /*preprossing for the OPENCL environment */

#include <openCL_GPU_funcDec.h>

cqCommandQueueGPU::cqCommandQueueGPU(cl_command_queue* CommandQueue)
{
  cqCommandQueue = CommandQueue;
}
// setter
void cqCommandQueueGPU::set_cqCommandQueueGPU(cl_command_queue* CommandQueue)
{
  cqCommandQueue = CommandQueue;
}
// getter
cl_command_queue* cqCommandQueueGPU::get_cqCommandQueueGPU()
{
  return cqCommandQueue;
}
// the destructor
cqCommandQueueGPU::~cqCommandQueueGPU()
{
  std::cout << "Object cpPlatformGPU has been destroyed" << std::endl;
}

#endif /* OPENCL */
