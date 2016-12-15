/*******************************************************************************
 *     Author: Frederic D.R. Bonnet date: 08th of April 2011. 15:01am
 *
 * Name:
 * SourcePathGPU.cpp - basic definitions used in all modules.
 *
 * Description:
 * classes to set and get the name of the source pathj
 *******************************************************************************
 */
#include <iostream>

#if defined (OPENCL) /*preprossing for the OPENCL environment */

#include <openCL_GPU_funcDec.h>

SourcePathGPU::SourcePathGPU(char* Path)
{
  cSourcePath = Path;
}

void SourcePathGPU::set_SourcePathGPU(char* Path)
{
  cSourcePath = Path;
}
char* SourcePathGPU::get_SourcePathGPU()
{
  return cSourcePath;
}

SourcePathGPU::~SourcePathGPU()
{
  std::cout << "Object SourcePathGPU has been destroyed" << std::endl;
}

#endif /* OPENCL */
