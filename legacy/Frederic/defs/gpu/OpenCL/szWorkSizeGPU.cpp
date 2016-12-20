/*******************************************************************************
 *     Author: Frederic D.R. Bonnet date: 08th of April 2011. 15:01am
 *
 * Name:
 * szWorkSizeGPU.cpp - basic definitions used in all modules.
 *
 * Description:
 * classes to set and get the szWorkSizeGPU
 *******************************************************************************
 */
#include <iostream>

#if defined (OPENCL) /*preprossing for the OPENCL environment */

#include <openCL_GPU_funcDec.h>

szWorkSizeGPU::szWorkSizeGPU(size_t LocalWorkSize, unsigned int nn )
{
  N = nn;
  szLocalWorkSize = LocalWorkSize;
}
//setter of the N
void szWorkSizeGPU::set_szVectorGPU(unsigned int nn)
{
  N = nn;
}
//getter of the szLocalWorkSize
unsigned int szWorkSizeGPU::get_szVectorGPU()
{
  return N;
}
//setter of the szLocalWorkSize
void szWorkSizeGPU::set_szLocalWorkSizeGPU(size_t LocalWorkSize)
{
  szLocalWorkSize = LocalWorkSize;
}
//getter of the szLocalWorkSize
size_t szWorkSizeGPU::get_szLocalWorkSizeGPU()
{
  return szLocalWorkSize;
}
//setter of the szGlobalWorkSizeGPU
void szWorkSizeGPU::set_szGlobalWorkSizeGPU()
{
  szGlobalWorkSize = shrRoundUp((int)szLocalWorkSize, N);
}
//getter of the szGlobalWorkSizeGPU
size_t szWorkSizeGPU::get_szGlobalWorkSizeGPU()
{
  return szGlobalWorkSize;
}
// the destructor
szWorkSizeGPU::~szWorkSizeGPU()
{
  std::cout << "Object szWorkSizeGPU has been destroyed" << std::endl;
}

#endif /* OPENCL */
