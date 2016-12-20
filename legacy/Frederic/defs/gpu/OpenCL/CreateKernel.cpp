/*******************************************************************************
 *     Author: Frederic D.R. Bonnet date: 12th of April 2011. 11:19am
 *
 * Name:
 * Createkernel.cpp - basic definitions used in all modules.
 *
 * Description:
 * routine to create the kernel on the card and set the context to the class
 * set_cpPlatformGPU on the GPU in OpenCL 
 *******************************************************************************
 */

#if defined (OPENCL) /*preprossing for the OPENCL environment */

#include <oclUtils.h>
#include <CL/cl.h>
#include <shrUtils.h>

extern "C" cl_kernel Createkernel(cl_program , char *);
extern "C" void Cleanup(int);

cl_kernel Createkernel(cl_program cpProgram, char *kernelname)
{
  cl_kernel ckKernel;             // OpenCL kernel
  cl_int ciErrNum = CL_SUCCESS;

  //start of the execution commands

  shrLog("Creating the Kernel %s...\n",kernelname);

  //creatign the Kernel
  ckKernel = clCreateKernel(cpProgram, kernelname, &ciErrNum);
  if ( ciErrNum != CL_SUCCESS)
    {
      shrLog("Error %i: error in clCreateKernel, near Line %i in file %s %s\n",ciErrNum,__LINE__,__FILE__,__FUNCTION__);
      clReleaseKernel(ckKernel);
      Cleanup(EXIT_FAILURE);
    }
  else
    {
      shrLog("Success %i: success in clCreateKernel, near Line %i in file %s %s\n",ciErrNum,__LINE__,__FILE__,__FUNCTION__);
    }
  return ckKernel;               //returning the created kernel.
} 

#endif /* OPENCL */
