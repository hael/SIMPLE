/*******************************************************************************
 *     Author: Frederic D.R. Bonnet date: 08th of April 2011. 12:19am
 *
 * Name:
 * routine_sAX.cpp - basic definitions used in all modules.
 *
 * Description:
 * routine that creates, alloacte the memeory on the buffer and
 * launches the kernel on the GPU in OpenCL routine that contraols the 
 * launch of the kernel on the card. Ity gets the context commandqueu etc and
 * strats the kernel on the card. This is the controller routine and the ones
 * that needs to be adapted to the problem.
 *******************************************************************************
 */

#if defined (OPENCL) /*preprossing for the OPENCL environment */

#include <global.h>
#include <openCL_GPU_funcDec.h>

#include <simul.h>
/*
#include <stdio.h>
#include <stdlib.h>
#include <oclUtils.h>
#include <CL/cl.h>
#include <shrUtils.h>

#include <global.h>
#include <simul.h>
*/
//-----------------The extern routines----------------------------//
//extern "C" void routine_sAX(int, float, cl_mem , int, size_t, size_t)
//extern "C" cl_kernel Createkernel(cl_program , char * );
//extern "C" void LaunchKernel( cl_command_queue , cl_kernel , size_t, size_t );
//extern "C" void Cleanup(int);
//----------------------------------------------------------------//

//void routine_sAX(int n, float alpha, cl_mem x, int incx )
// , size_t szLocalWorkSize, size_t szGlobalWorkSize)
void routine_saX(int n, float alpha, int incx )
{
  cl_mem   cmDevRET;                // OpenCL device destination buffer 
  cl_kernel *ckKernel;              // the created kernel
  //  float ret;
  float *ret;
  ret = (float *)malloc(sizeof(cl_float) * worksizeGPU -> get_szVectorGPU());

  //  cl_int ciErrNum = CL_SUCCESS;

  //retrieving the cxGPUContext, commandQueue, ciDeviceCount, cpProgram
  //from the object getters

  cl_context cxGPUContext_2        = context         -> get_cxGPUContextGPU();    //getting the context
  cl_command_queue* commandQueue_2 = CommandQueueGPU -> get_cqCommandQueueGPU();  //getting the command queue
  cl_uint ciDeviceCount_2          = devicecount     -> get_ciDeviceCountGPU();   //getting the device count
  cl_program cpProgram_2           = programGPU      -> get_cpProgramGPU();       //getting the OpenCl program
  size_t szLocalWorkSize           = worksizeGPU     -> get_szLocalWorkSizeGPU(); //getting the szLocalWorkSizeGPU 
  size_t szGlobalWorkSize          = worksizeGPU     -> get_szGlobalWorkSizeGPU();//getting the szGlobalWorkSizeGPU

  cl_mem x                         = DeVvectorGPU    -> get_szVectorGPU();        //getting the vector cmDeVvector1

  //allocating the memory for the kernel

  ckKernel =(cl_kernel *)malloc(sizeof(cl_kernel)*ciDeviceCount_2);

  for ( unsigned int idev = 0 ; idev < ciDeviceCount_2 ; idev++ )
    {
      //cmDevRET  = clCreateBuffer(cxGPUContext_2, CL_MEM_READ_WRITE, sizeof(cl_float), NULL, &ciErrNum);
      
      //Creating and getting the kernel.
      ckKernel[idev] = Createkernel(cpProgram_2, "saX");

      ciErrNum  = clSetKernelArg(ckKernel[idev], 0, sizeof(cl_int  ),  (void*)&n    );
      ciErrNum |= clSetKernelArg(ckKernel[idev], 1, sizeof(cl_float),  (void*)&alpha);
      ciErrNum |= clSetKernelArg(ckKernel[idev], 2, sizeof(cl_mem  ),  (void*)&x    );
      ciErrNum |= clSetKernelArg(ckKernel[idev], 3, sizeof(cl_int  ),  (void*)&incx );
      if (ciErrNum != CL_SUCCESS)
	{
	  shrLog("Error %i in routine_saX, clSetKernelArg, near Line %i in file %s %d\n",ciErrNum,__LINE__, __FILE__,__FUNCTION__);
	  clReleaseMemObject(x);
	  Cleanup(EXIT_FAILURE);
	}
      else
	{
	  shrLog("Success %i: success in clSetKernelArg, near Line %i in file %s %s\n",ciErrNum,__LINE__,__FILE__,__FUNCTION__);
	}

      //launch the kernel on the card

      LaunchKernel(commandQueue_2[idev], ckKernel[idev] , szGlobalWorkSize, szLocalWorkSize);

      //retrieving the data from the GPU

      //ciErrNum = clEnqueueReadBuffer(commandQueue_2[idev], cmDevRET, CL_TRUE, 0, sizeof(cl_float), &ret, 0, NULL, NULL);
      ciErrNum = clEnqueueReadBuffer(commandQueue_2[idev], x, CL_TRUE, 0, sizeof(cl_float)*10, ret, 0, NULL, NULL);
      if ( ciErrNum != CL_SUCCESS )
	{
	  shrLog("Error %i: error in clEnqueueReadBuffer, near Line %i in file %s %s\n",ciErrNum,__LINE__,__FILE__,__FUNCTION__);
	  clReleaseMemObject(x);
	  clReleaseKernel(ckKernel[idev]);
	  Cleanup(EXIT_FAILURE);
	}
      else
	{
	  shrLog("Success %i: success in clEnqueueReadBuffer, near Line %i in file %s %s\n",ciErrNum,__LINE__,__FILE__,__FUNCTION__);
	}

      printf("the vector output in routine_saX =[");
      for (int i = 0 ; i < 10 ; i++)
	{
	  printf(" %4.2f",ret[i]);
	}
      printf(" ]\n");

    }

  //cleaning up the ocl resources
  for (unsigned int idev = 0 ; idev < ciDeviceCount_2 ; idev++ )
    {
      clReleaseKernel( ckKernel[idev] );
      clReleaseCommandQueue( commandQueue_2[idev] );
    }
}
extern "C" void routine_sax_gpu_(int, float, int ) __attribute__((weak,alias("routine_saX")));

#endif /* OPENCL */
