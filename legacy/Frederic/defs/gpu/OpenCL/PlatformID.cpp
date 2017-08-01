/*******************************************************************************
 *     Author: Frederic D.R. Bonnet date: 12th of April 2011. 12:19am
 *
 * Name:
 * PlatformID.cpp - basic definitions used in all modules.
 *
 * Description:
 * routine to get the NVIDIA platform id and set to the class
 * set_cpPlatformGPU on the GPU in OpenCL 
 *******************************************************************************
 */

#if defined (OPENCL) /*preprossing for the OPENCL environment */

#include <global.h>
#include <openCL_GPU_funcDec.h>

//extern "C" void PlatformID();

void PlatformID()
{

  //Get the NVIDIA platform
  ciErrNum = oclGetPlatformID(&cpPlatform);
  if (ciErrNum != CL_SUCCESS)
    {
      shrLog("Error %i: error in oclGetPlatformID, Failed to getPlatformID, near Line %i in file %s %s\n",ciErrNum,__LINE__,__FILE__,__FUNCTION__);
      //return ciErrNum;
    }
  else
    {
      shrLog("Success %i: success in oclGetPlatformID, -- The NVIDIA platform loaded successfully... near Line %i in file %s %s\n",ciErrNum,__LINE__,__FILE__,__FUNCTION__);
    }

  platform -> set_cpPlatformGPU(cpPlatform); //return cpPlatform;
}
extern "C" void platformid_gpu_() __attribute__((weak,alias("PlatformID")));

#endif /* OPENCL */
