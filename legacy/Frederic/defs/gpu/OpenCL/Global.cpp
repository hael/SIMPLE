/*******************************************************************************
 *     Author: Frederic D.R. Bonnet date: 12th of April 2011. 12:19am
 *
 * Name:
 * Global.cpp - basic definitions used in all modules.
 *
 * Description:
 * file to englodbe all the global variables and
 * class declarion for the header file global.h
 *******************************************************************************
 */

#if defined (OPENCL) /*preprossing for the OPENCL environment */

#include <oclUtils.h>
#include <CL/cl.h>
#include <shrUtils.h>

//include file for the global objects 
#include <openCL_GPU_funcDec.h>
#include <simul.h>

//global variables used through out the code
char *cSourcePath = "../../../src/gpu/kernels/OpenCL/";
char *cSourceFile = "kernel_saX.cc";

// OpenCL global variables
char cBuffer[1024];
bool bPassed = true;

unsigned int nn;
size_t LocalWorkSize;

cl_int ciErrNum = CL_SUCCESS;
cl_platform_id cpPlatform = NULL;
3cl_uint ciDeviceCount;
cl_device_id* cdDevices = NULL;
cl_context cxGPUContext;
cl_command_queue* commandQueue;  //the command queue array
cl_program cpProgram;           // OpenCl program
cl_mem Vvector;

//global objects declared for OpenCl
cpPlatformGPU *platform = new cpPlatformGPU(cpPlatform);
ciDeviceCountGPU *devicecount = new ciDeviceCountGPU(ciDeviceCount);
cdDevicesGPU *devicesGPU = new cdDevicesGPU(cdDevices);
cxGPUContextGPU *context = new cxGPUContextGPU(cxGPUContext);
cqCommandQueueGPU *CommandQueueGPU = new cqCommandQueueGPU(commandQueue);
cpProgramGPU *programGPU = new cpProgramGPU(cpProgram);
szWorkSizeGPU *worksizeGPU = new szWorkSizeGPU(LocalWorkSize,nn);
cmDeVvectorGPU *DeVvectorGPU = new cmDeVvectorGPU(Vvector);

#endif /* OPENCL */
