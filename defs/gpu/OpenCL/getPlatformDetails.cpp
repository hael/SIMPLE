/*******************************************************************************
 *     Author: Frederic D.R. Bonnet date: 12th of April 2011. 12:19am
 *
 * Name:
 * getPlatformDetails.cpp - basic definitions used in all modules.
 *
 * Description:
 * routine to get the getPlatformDetails on the NVIDIA card
 *******************************************************************************
 */

#if defined (OPENCL) /*preprossing for the OPENCL environment */

#include <global.h>
#include <openCL_GPU_funcDec.h>

//extern "C" void getPlatformDetails();

void getPlatformDetails()
{
  //the object variables
  cl_platform_id cpPlatform_2 = NULL;
  cl_device_id *devices;

  // retrieving the cpPlatform details from the getter fo the files from the object

  shrLog(" in getplatformdetails(): near Line %i in file %s %s\n",__LINE__,__FILE__,__FUNCTION__);
  cpPlatform_2 = platform -> get_cpPlatformGPU();
  shrLog(" in getplatformdetails(): near Line %i in file %s %s\n",__LINE__,__FILE__,__FUNCTION__);
  ciErrNum = clGetPlatformInfo(cpPlatform_2, CL_PLATFORM_NAME, sizeof(cBuffer), cBuffer, NULL);
  if ( ciErrNum != CL_SUCCESS )
    {
      shrLog("Error %i: error in clGetPlatformInfo, Failed to getPlatformID, near Line %i in file %s %s\n",ciErrNum,__LINE__,__FILE__,__FUNCTION__);
      shrLog(" Error %i in getPlatformDetails Call !!!\n\n", ciErrNum);
      bPassed = false;
    }
  else 
    {
      shrLog("Success %i: success in clGetPlatformInfo, -- The NVIDIA platform loaded successfully... near Line %i in file %s %s\n",ciErrNum,__LINE__,__FILE__,__FUNCTION__);
      shrLog(" CL_PLATFORM_NAME: \t%s\n", cBuffer);
    }

  shrLog(" in getplatformdetails(): near Line %i in file %s %s\n",__LINE__,__FILE__,__FUNCTION__);

  ciErrNum = clGetPlatformInfo(cpPlatform, CL_PLATFORM_NAME, sizeof(cBuffer), cBuffer, NULL);
  if ( ciErrNum != CL_SUCCESS )
    {
      shrLog("Error %i: error in clGetPlatformInfo, Failed to getPlatformID, near Line %i in file %s %s\n",ciErrNum,__LINE__,__FILE__,__FUNCTION__);
      shrLog(" Error %i in getPlatformDetails Call !!!\n\n", ciErrNum);
      bPassed = false;
    }
  else 
    {
      shrLog("Success %i: success in clGetPlatformInfo, -- The NVIDIA platform loaded successfully... near Line %i in file %s %s\n",ciErrNum,__LINE__,__FILE__,__FUNCTION__);
      shrLog(" CL_PLATFORM_NAME: \t%s\n", cBuffer);
    }

  //Get the platform details, name and version...

  ciErrNum = clGetPlatformInfo(cpPlatform, CL_PLATFORM_NAME, sizeof(cBuffer), cBuffer, NULL);
  if ( ciErrNum == CL_SUCCESS )
    {
      shrLog(" CL_PLATFORM_NAME: \t%s\n", cBuffer);
    }
  else 
    {
      shrLog(" Error %i in getPlatformDetails Call !!!\n\n", ciErrNum);
      bPassed = false;
    }
  ciErrNum = clGetPlatformInfo(cpPlatform, CL_PLATFORM_VERSION, sizeof(cBuffer), cBuffer , NULL);
  if ( ciErrNum == CL_SUCCESS )
    {
      shrLog(" CL_PLATFORM_VERSION: \t%s\n", cBuffer);
    }
  else
    {
      shrLog(" Error %i in getPlatformDetails Call !!!\n\n", ciErrNum);
      bPassed = false;
    }
  
  ciErrNum = clGetPlatformInfo(cpPlatform, CL_PLATFORM_PROFILE, sizeof(cBuffer), cBuffer, NULL);
  if ( ciErrNum == CL_SUCCESS )
    {
      shrLog(" CL_PLATFORM_PROFILE: \t%s\n", cBuffer);
    }
  else 
    {
      shrLog(" Error %i in getPlatformDetails Call !!!\n\n", ciErrNum);
      bPassed = false;
    }
  shrLog(" OpenCL SDK Revision: \t%s\n\n\n", OCL_SDKREVISION);

  shrLog("OpenCL Device Info: \n");
  ciErrNum = clGetDeviceIDs(cpPlatform, CL_DEVICE_TYPE_ALL, 0, NULL, &ciDeviceCount);
  if ( ciDeviceCount == 0 )
    {
      shrLog("No device found supporting OpenCL (return code %i)\n\n",ciErrNum);
    }
  else if ( ciErrNum != CL_SUCCESS )
    {
      shrLog(" Error %i in getPlatformDetails Call !!!\n\n", ciErrNum);
      bPassed = false;      
    }
  else 
    {
      shrLog("\n");
      shrLog("%u unit(s) found to support the OpenCL\n",ciDeviceCount);
      shrLog("\n");
      
      if ((devices = (cl_device_id*)malloc(sizeof(cl_device_id) * ciDeviceCount)) == NULL)
	{
	  shrLog(" Failed to allocate memory for devices !!!\n\n");
	  bPassed = false;
	}
      
      ciErrNum = clGetDeviceIDs(cpPlatform, CL_DEVICE_TYPE_ALL, ciDeviceCount, devices, &ciDeviceCount);
      if ( ciErrNum == CL_SUCCESS )
	{
	  for ( unsigned int i = 0; i < ciDeviceCount; i++ )
	    {
	      shrLog(" ---------------------------------\n");
	      clGetDeviceInfo(devices[i], CL_DEVICE_NAME, sizeof(cBuffer), &cBuffer, NULL);
	      shrLog(" Device name: %s\n",cBuffer);
	      shrLog(" ---------------------------------\n");
	      oclPrintDevInfo(LOGBOTH, devices[i]);
	    }
	}
      else
	{
	  shrLog(" Error %i in getPlatformDetails Call !!!\n\n", ciErrNum);
	  bPassed = false;
	}
    }
  devicecount -> set_ciDeviceCountGPU(ciDeviceCount);//return ciDeviceCount;
}
extern "C" void getplatformdetails_gpu_() __attribute__((weak,alias("getPlatformDetails")));

#endif /* OPENCL */
