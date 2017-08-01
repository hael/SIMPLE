/*******************************************************************************
 *     Author: Frederic D.R. Bonnet date: 08th of April 2011. 15:01am
 *
 * Name:
 * SourceFileGPU.cpp - basic definitions used in all modules.
 *
 * Description:
 * classes to set and get the name of the source file
 *******************************************************************************
 */
#include <iostream>

#if defined (OPENCL) /*preprossing for the OPENCL environment */

#include <openCL_GPU_funcDec.h>

SourceFileGPU::SourceFileGPU(char* Filename)
{
  cSourceFile = Filename;
}

void SourceFileGPU::set_SourceFileGPU(char* Filename)
{
  cSourceFile = Filename;
}
char* SourceFileGPU::get_SourceFileGPU()
{
  return cSourceFile;
}

SourceFileGPU::~SourceFileGPU()
{
  std::cout << "Object SourceFileGPU has been destroyed" << std::endl;
}

#endif /* OPENCL */
