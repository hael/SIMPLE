/*******************************************************************************
 *     Author: Frederic D.R. Bonnet date: 12th of April 2011. 12:19am
 *
 * Name:
 * global.h - basic definitions used in all modules.
 *
 * Description:
 * header file for the Global.cpp file to set the 
 * global variable swithin the cpp files and routines.
 * class declarion for the header file global.h
 *******************************************************************************
 */
#if defined (OPENCL) /*preprossing for the OPENCL environment */

#include <oclUtils.h>
#include <CL/cl.h>
#include <shrUtils.h>

//include file for the global objects 
#include <openCL_GPU_funcDec.h>

//global variables used through out the code
extern char* cSourcePath;
extern char* cSourceFile;

// OpenCL global variables
extern char cBuffer[1024];
extern bool bPassed;

extern cl_int ciErrNum;
extern cl_platform_id cpPlatform;
extern cl_uint ciDeviceCount;
extern cl_device_id* cdDevices;
extern cl_context cxGPUContext;

//global objects declared for OpenCl
extern cpPlatformGPU *platform;
extern ciDeviceCountGPU *devicecount;
extern cdDevicesGPU *devicesGPU;
extern cxGPUContextGPU *context;
extern cqCommandQueueGPU *CommandQueueGPU;
extern cpProgramGPU *programGPU;
extern szWorkSizeGPU *worksizeGPU;
extern cmDeVvectorGPU *DeVvectorGPU;

#endif /* OPENCL */
