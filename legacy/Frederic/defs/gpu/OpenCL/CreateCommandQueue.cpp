/*******************************************************************************
 *     Author: Frederic D.R. Bonnet date: 12th of April 2011.13:01am
 *
 * Name:
 * CreateCommandQueue.cpp - basic definitions used in all modules.
 *
 * Description:
 * routine to get the devices cdDevices on the NVIDIA platform id
 * and set to the class
 * set_cqCommandQueueGPU on the GPU in OpenCL.
 * included the OPENCL directives
 *******************************************************************************
 */

#if defined (OPENCL) /*preprossing for the OPENCL environment */

#include <global.h>
#include <openCL_GPU_funcDec.h>

#include <simul.h>

//extern "C" void CreateCommandQueue();
//extern "C" void Cleanup(int);

void CreateCommandQueue()
{
  cl_command_queue* commandQueue = new cl_command_queue[MAX_GPU_COUNT];   //the command queue array

  cl_device_id devices = NULL;
  //  cl_int ciErrNum = CL_SUCCESS;

  //local variabels

  cl_uint ciDeviceCount_2;
  cl_context cxGPUContext_2;

  //bringing the objects in for the ciDeviceCount

  ciDeviceCount_2 = devicecount -> get_ciDeviceCountGPU();
  cxGPUContext_2 = context -> get_cxGPUContextGPU();    

  // start of the execution commands

  shrLog("\n");
  shrLog("The number of devices is %i: \n", ciDeviceCount_2);
  shrLog("\n");

  for (unsigned int idev = 0 ; idev < ciDeviceCount_2 ; idev++)
    {
      shrLog("Devices number %d: ",idev);

      devices = oclGetDev(cxGPUContext_2,idev);
      
      oclPrintDevName(LOGBOTH, devices);            

      shrLog("\n");

      // create command queue
      commandQueue[idev] = clCreateCommandQueue(cxGPUContext_2, devices, 0, &ciErrNum);
      if (ciErrNum != CL_SUCCESS)
	{
	  shrLog("Error %i: error in clCreateCommandQueue, near Line %i in file %s %s\n",ciErrNum,__LINE__,__FILE__,__FUNCTION__);
	  clReleaseCommandQueue(commandQueue[idev]);
	  Cleanup(EXIT_FAILURE);
	}
      else
	{
	  shrLog("Success %i: success in clCreateCommandQueue, near Line %i in file %s %s\n",ciErrNum,__LINE__,__FILE__,__FUNCTION__);
	}
#ifdef GPU_PROFILING
      clSetCommandQueueProperty(commandQueue[idev], CL_QUEUE_PROFILING_ENABLE, CL_TRUE, NULL);
      if (ciErrNum != CL_SUCCESS)
	{
	  shrLog("Error %i: error in clSetCommandQueueProperty, near Line %i in file %s %s\n",ciErrNum,__LINE__,__FILE__,__FUNCTION__);
	  shrLog("Error %i in clSetCommandQueueProperty call !!!\n\n", ciErrNum);
	  clReleaseCommandQueue(commandQueue[idev]);
	  Cleanup(EXIT_FAILURE);
	}
      else
	{
	  shrLog("Success %i: success in clSetCommandQueueProperty, near Line %i in file %s %s\n",ciErrNum,__LINE__,__FILE__,__FUNCTION__);
	}
#endif

    }

  CommandQueueGPU -> set_cqCommandQueueGPU(commandQueue);//return commandQueue;

}
extern "C" void createcommandqueue_gpu_() __attribute__((weak,alias("CreateCommandQueue")));

#endif /* OPENCL */
