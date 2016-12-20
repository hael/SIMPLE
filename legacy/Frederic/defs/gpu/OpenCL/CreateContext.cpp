/*******************************************************************************
 *     Author: Frederic D.R. Bonnet date: 12th of April 2011. 11:19am
 *
 * Name:
 * CreateContext.cpp - basic definitions used in all modules.
 *
 * Description:
 * routine to create the context and set the context to the class
 * set_cpPlatformGPU on the GPU in OpenCL 
 *******************************************************************************
 */

#if defined (OPENCL) /*preprossing for the OPENCL environment */

#include <global.h>
#include <openCL_GPU_funcDec.h>

//extern "C" void CreateContext();
//extern "C" void Cleanup(int);

void CreateContext()
{
  //the object variables
  cl_uint ciDeviceCount_2;
  cl_platform_id cpPlatform_2 = NULL;
  cl_device_id* cdDevices_2 = NULL;

  //global variables: ./include/OpenCL/global.h

  //bringing the objects in for the ciDeviceCount and cpPlatform

  cpPlatform_2 = platform -> get_cpPlatformGPU();
  ciDeviceCount_2 = devicecount -> get_ciDeviceCountGPU();

  //Get the devicesId and compare with the passed object
  ciErrNum = clGetDeviceIDs(cpPlatform_2, CL_DEVICE_TYPE_GPU, 0, NULL, &ciDeviceCount);
  if( ciDeviceCount != ciDeviceCount_2 )
    {
      shrLog("Error %i: error in clGetDeviceIDs, DeviceCount do not match!!! , near Line %i in file %s %s\n",ciErrNum,__LINE__,__FILE__,__FUNCTION__);
    }

  cdDevices_2 = (cl_device_id *)malloc(ciDeviceCount_2 * sizeof(cl_device_id) );
  cdDevices_2 = devicesGPU -> get_cdDevicesGPU();

  //Create the context
  cxGPUContext = clCreateContext(0, ciDeviceCount_2, cdDevices_2, NULL, NULL, &ciErrNum);
  if (ciErrNum != CL_SUCCESS)
    {
      shrLog("Error %i: error in clCreateContext, near Line %i in file %s %s\n",ciErrNum,__LINE__,__FILE__,__FUNCTION__);
      clReleaseContext(cxGPUContext);
      Cleanup(EXIT_FAILURE);
    }
  else
    {
      shrLog("Success %i: success in clCreateContext, near Line %i in file %s %s\n",ciErrNum,__LINE__,__FILE__,__FUNCTION__);
    }

  context -> set_cxGPUContextGPU(cxGPUContext);//return cxGPUContext;

}
extern "C" void createcontext_gpu_() __attribute__((weak,alias("CreateContext")));

#endif /* OPENCL */
