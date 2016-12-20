/*******************************************************************************
 *     Author: Frederic D.R. Bonnet date: 12th of April 2011.15:01am
 *
 * Name:
 * cpProgramGPU.cpp - basic definitions used in all modules.
 *
 * Description:
 * classes to set and get the cpProgramGPU on the GPU in OpenCL
 *******************************************************************************
 */
#include <iostream>

#if defined (OPENCL) /*preprossing for the OPENCL environment */

#include <openCL_GPU_funcDec.h>

cpProgramGPU::cpProgramGPU(cl_program Prog)
{
  cpProgram = Prog;
}
// setter
void cpProgramGPU::set_cpProgramGPU(cl_program Prog)
{
  cpProgram = Prog;
}
// getter
cl_program cpProgramGPU::get_cpProgramGPU()
{
  return cpProgram;
}
// the destructor
cpProgramGPU::~cpProgramGPU()
{
  std::cout << "Object cpPlatformGPU has been destroyed" << std::endl;
}

#endif /* OPENCL */
