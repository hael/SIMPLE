/*******************************************************************************
 *     Author: Frederic D.R. Bonnet date: 12th of April 2011. 12:19am
 *
 * Name:
 * LaunchKernel.cpp - basic definitions used in all modules.
 *
 * Description:
 * routine that launches the kernel in question on the card
 *******************************************************************************
 */

#if defined (OPENCL) /*preprossing for the OPENCL environment */

#include <oclUtils.h>
#include <CL/cl.h>
#include <shrUtils.h>

#include <global.h>
#include <openCL_GPU_funcDec.h>

//extern "C" void LaunchKernel( cl_command_queue , cl_kernel , size_t, size_t );
//extern "C" void Cleanup(int);

void LaunchKernel(cl_command_queue commandQueue, cl_kernel ckKernel , size_t szGlobalWorkSize, size_t szLocalWorkSize)
{
  //  size_t szGlobalWorkSize;       // 1D var for Total # of work items
  //size_t szLocalWorkSize;	 // 1D var for # of work items in the work group
  cl_int ciErrNum = CL_SUCCESS;

  ciErrNum = clEnqueueNDRangeKernel(commandQueue, ckKernel, 1, NULL, &szGlobalWorkSize, &szLocalWorkSize, 0, NULL, NULL);
  if ( ciErrNum  != CL_SUCCESS)
    {
      shrLog("Error %i: error in clEnqueueNDRangeKernel, near Line %i in file %s %s\n",ciErrNum,__LINE__,__FILE__,__FUNCTION__);
      Cleanup(EXIT_FAILURE);
    }
  else
    {
      shrLog("Success %i: success in clEnqueueNDRangeKernel, near Line %i in file %s %s\n",ciErrNum,__LINE__,__FILE__,__FUNCTION__);
    }
} 

#endif /* OPENCL */
