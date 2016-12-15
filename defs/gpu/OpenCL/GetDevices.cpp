/*******************************************************************************
 *     Author: Frederic D.R. Bonnet date: 12th of April 2011. 12:19am
 *
 * Name:
 * GetDevices.cpp - basic definitions used in all modules.
 *
 * Description:
 * routine to get the devices cdDevices on the NVIDIA platform id
 * and set to the class set_cdDevicesGPU
 *******************************************************************************
 */

#if defined (OPENCL) /*preprossing for the OPENCL environment */

#include <global.h>
#include <openCL_GPU_funcDec.h>

//extern "C" void GetDevices();
//extern "C" void Cleanup(int);

void GetDevices()
{
  cl_uint ciDeviceCount_2;
  cl_platform_id cpPlatform_2 = NULL;

  //bringing the objects in for the ciDeviceCount and cpPlatform

  cpPlatform_2 = platform -> get_cpPlatformGPU();
  ciDeviceCount_2 = devicecount -> get_ciDeviceCountGPU();

  //Get the devices
  ciErrNum = clGetDeviceIDs(cpPlatform_2, CL_DEVICE_TYPE_GPU, 0, NULL, &ciDeviceCount);
  if( ciDeviceCount != ciDeviceCount_2 )
    {
      shrLog("Error %i: error in clGetDeviceIDs, DeviceCount do not match!!! , near Line %i in file %s %s\n",ciErrNum,__LINE__,__FILE__,__FUNCTION__);
    }

  cdDevices = (cl_device_id* )malloc(ciDeviceCount * sizeof(cl_device_id) );
  ciErrNum = clGetDeviceIDs(cpPlatform_2, CL_DEVICE_TYPE_GPU, ciDeviceCount, cdDevices, NULL);
  if (ciErrNum != CL_SUCCESS)
    {
      shrLog("Error %i: error in clGetDeviceIDs, Failed to GetDevices, near Line %i in file %s %s\n",ciErrNum,__LINE__,__FILE__,__FUNCTION__);
      free(cdDevices);
      //Cleanup(EXIT_FAILURE);
    }
  else
    {
      shrLog("Success %i: success in clGetDeviceIDs, -- The NVIDIA GetDevices loaded successfully... near Line %i in file %s %s\n",ciErrNum,__LINE__,__FILE__,__FUNCTION__);
    }

  devicesGPU -> set_cdDevicesGPU(cdDevices);//return cdDevices;

}
extern "C" void getdevices_gpu_() __attribute__((weak,alias("GetDevices")));

#endif /* OPENCL */
