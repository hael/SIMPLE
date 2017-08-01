/*
 *   -- SIMPLE addon
 *      Author: Frederic Bonnet, Date: 21st May 2015
 *
 *      May 2015
 *
 *   code to obtain information from the GPU cards
 *
 * @precisions normal z -> s d c
 */

#include "get_deviceQuery_gpu.h"

#if defined (CUDA)
cudaDeviceProp deviceProp;
#endif

#ifdef __cplusplus
extern "C" {
#endif
/***************************************************************************/
/**
 *  FORTRAN API - math functions (simple interface)
 **/

#if defined (CUDA) /*preprossing for the CUDA environment */

  /* sets the device for computation */
  int setDevice(int *idev) {
    int rc = RC_SUCCESS;
    int dev = *idev;
  
    cudaError_t error_id;
    cudaDeviceProp deviceProp_m;

    //cudaSetDevice(dev); here sets the device for computation
    error_id = cudaGetDeviceProperties(&deviceProp_m, dev);
    if (error_id != cudaSuccess) { rc = get_error_id (error_id); }

    error_id = cudaSetDevice(dev); //here sets the device for computation
    if (error_id != cudaSuccess) { rc = get_error_id (error_id); }

    return rc;
  }
  /* wrapper function to find the most suitable device on system */
  int get_findCudaDevice(int *idev) {
    int rc = RC_SUCCESS;

    int dev = 0;

    cudaError_t error_id;
    cudaDeviceProp deviceProp_m;

    //cudaSetDevice(dev); here sets the device for computation
    error_id = cudaGetDeviceProperties(&deviceProp_m, dev);
    if (error_id != cudaSuccess) { rc = get_error_id (error_id); }

    //method to find the device with the most capabilities
    *idev = findCudaDevice(deviceProp_m);

    return rc;
  }
  
  /* initialising the data structure */
  int deviceDetails_init(int *idev, deviceDetails_t *devD) {
    int rc = RC_SUCCESS; //return code
    int dev = *idev;
    cudaError_t error_id;

    //cudaSetDevice(dev); here sets the device for computation
    error_id = cudaGetDeviceProperties(&deviceProp, dev);
    if (error_id != cudaSuccess) { rc = get_error_id (error_id); }
    
    //method to find the device with the most capabilities
    //findCudaDevice(deviceProp);
    //TODO: add logic in case of additional multiple GPU selection process

    return rc;
  }
  /* freeing the data structure from memory */
  int deviceDetails_free(deviceDetails_t *devD) {
    int rc = RC_SUCCESS; //return code

    /*TODO: insert freeing statement in case of allocations */

    return rc;
  }
  /* determine if peer to peer is allowed */
  int get_peer_to_peer_capabilities(deviceDetails_t *devD) {
    int rc = RC_SUCCESS;

    //devD->is_p2p = false; /*initialises the p2p to false */

    cudaDeviceProp a_deviceProp[MAX_N_GPU];
    int gpuid[MAX_N_GPU];
    int ip2p = 0;
    for ( int i = 0 ; i < devD->ndev ; i++ ) {
      cudaGetDeviceProperties(&a_deviceProp[i], i);
     gpuid[ip2p++] = i;
    }
    /* getting the combinations for p2p support */
    int can_access_peer_0_1, can_access_peer_1_0;
    if ( devD->ndev >= 2) {
      if ( ip2p > 0 ) {

	for ( int i = 0 ; i < (ip2p-1) ; i++) {
	  for ( int j = 1 ; j < ip2p ; j++) {
	    cudaDeviceCanAccessPeer(&can_access_peer_0_1, gpuid[i], gpuid[j]);
	    printf(ANSI_COLOR_BRIGHT_BLUE"> Peer access from "
                   ANSI_COLOR_BRIGHT_YELLOW"%s (GPU%d)"
                   ANSI_COLOR_BRIGHT_BLUE" -> "
                   ANSI_COLOR_BRIGHT_YELLOW"%s (GPU%d) : %s\n"
                   ANSI_COLOR_RESET,
		   a_deviceProp[gpuid[i]].name, gpuid[i],
		   a_deviceProp[gpuid[j]].name, gpuid[j] ,
		   can_access_peer_0_1 ?
                   ANSI_COLOR_BRIGHT_GREEN"Yes" :
                   ANSI_COLOR_BRIGHT_RED"No"
                   ANSI_COLOR_RESET);
	    devD->is_p2p[i][j] = can_access_peer_0_1;
	  }
	}

	for ( int i = 1 ; i < ip2p ; i++) {
	  for ( int j = 0 ; j < (ip2p-1) ; j++) {
	    cudaDeviceCanAccessPeer(&can_access_peer_1_0, gpuid[i], gpuid[j]);
	    printf(ANSI_COLOR_BRIGHT_BLUE"> Peer access from "
                   ANSI_COLOR_BRIGHT_YELLOW"%s (GPU%d)"
                   ANSI_COLOR_BRIGHT_BLUE" -> "
                   ANSI_COLOR_BRIGHT_YELLOW"%s (GPU%d) : %s\n"
                   ANSI_COLOR_RESET,
		   a_deviceProp[gpuid[i]].name, gpuid[i],
		   a_deviceProp[gpuid[j]].name, gpuid[j] ,
		   can_access_peer_1_0 ?
                   ANSI_COLOR_BRIGHT_GREEN "Yes" :
                   ANSI_COLOR_BRIGHT_RED"No"
                   ANSI_COLOR_RESET);
	    devD->is_p2p[i][j] = can_access_peer_1_0;
	  }
	}

      }
    }

    return rc;
  }
  /* getting the error correcting code memory device true or flase*/
  int get_eec_support(int *ecc, int *idev, deviceDetails_t *devD) {
    int rc = RC_SUCCESS;
    
    *ecc = deviceProp.ECCEnabled ? 1 : get_ecc_warning_message(devD);
    devD->is_ecc = *ecc;
    return rc;
  }
  /* getting the threads details */
  int get_thread_details(int *warpSize, int *max_threads_per_mp,
                         int *max_threads_per_blk,
			 deviceDetails_t *devD) {
    int rc = RC_SUCCESS;

    *warpSize = deviceProp.warpSize;
    devD->warpSze = *warpSize;
    *max_threads_per_mp = deviceProp.maxThreadsPerMultiProcessor;
    devD->maxthreads_per_mp = *max_threads_per_mp;
    *max_threads_per_blk = deviceProp.maxThreadsPerBlock;
    devD->maxthreads_per_blk = *max_threads_per_blk;

    return rc;
  }
  /* getting the number of registers */
  int get_nregisters(int *nregisters_per_blk, deviceDetails_t *devD) {
    int rc = RC_SUCCESS;

    *nregisters_per_blk = deviceProp.regsPerBlock;
    devD->nregisters_per_blk = *nregisters_per_blk;

    return rc;
  }
  /* getting # of Multiprocessors, CUDA cores/MP and total # of CUDA cores */ 
  int get_CUDA_cores(int *idev, int *gpu_major, int *gpu_minor, int *nmp,
		     int *cuda_cores_per_mp, int *ncuda_cores,
		     deviceDetails_t *devD) {
    int rc = RC_SUCCESS;
    int dev = *idev;
    int major = *gpu_major; int minor = *gpu_minor;

    *nmp = deviceProp.multiProcessorCount;
    devD->nMultiProc = *nmp;
    *cuda_cores_per_mp = _ConvertSMVer2Cores(deviceProp.major, deviceProp.minor);
    devD->ncudacores_per_MultiProc = *cuda_cores_per_mp;
    *ncuda_cores = _ConvertSMVer2Cores(deviceProp.major, deviceProp.minor) *
                   deviceProp.multiProcessorCount;
    devD->ncudacores = *ncuda_cores;
    if (devD->ncudacores == 0 ) { rc = get_cuda_cores_error(devD->ncudacores); }

    /*requiring at leat a gpu_major.gpu_minor architecture Tesla card */ 
    bool bVal = checkCudaCapabilities(major,minor);
    /*    printf(ANSI_COLOR_BRIGHT_BLUE"Device %d: < "
           ANSI_COLOR_BRIGHT_YELLOW"%16s"ANSI_COLOR_BRIGHT_BLUE
           " >, Compute SM %d.%d detected, **suitable: %s**\n",
	   dev, deviceProp.name, deviceProp.major, deviceProp.minor,
           bVal?ANSI_COLOR_BRIGHT_GREEN"yes"ANSI_COLOR_RESET :
                  ANSI_COLOR_BRIGHT_RED"no" ANSI_COLOR_RESET);
    */
    devD->is_SMsuitable = bVal;

    return rc;
  }
  // General check for CUDA GPU SM Capabilities
  inline bool checkCudaCapabilities(int major_version, int minor_version) {
    if ((deviceProp.major > major_version) ||
	(deviceProp.major == major_version && deviceProp.minor >= minor_version)) {
        return true;
      } else { return false; }
  }
  /* getting the memory from the devices */
  int get_tot_global_mem_MB(float *totalGlobalMem_MB, deviceDetails_t *devD) {
    int rc = RC_SUCCESS;

    *totalGlobalMem_MB = (float)deviceProp.totalGlobalMem/1048576.0f;
    devD->tot_global_mem_MB = *totalGlobalMem_MB;

    return rc;
  }  
  /* getting the memory from the devices */
  int get_tot_global_mem_bytes(unsigned long long *totalGlobalMem_bytes,
                               deviceDetails_t *devD) {
    int rc = RC_SUCCESS;

    *totalGlobalMem_bytes = (unsigned long long) deviceProp.totalGlobalMem;
    devD->tot_global_mem_bytes = *totalGlobalMem_bytes;

    return rc;
  }  
  /* getting the driver version */
  int get_cuda_driverVersion(int *driverVersion , deviceDetails_t *devD) {
    int rc = RC_SUCCESS;
    cudaError_t error_id;

    error_id = cudaDriverGetVersion(driverVersion);
    if (error_id != cudaSuccess) { rc = get_error_id (error_id); }

    devD->d_ver = (float)*driverVersion/1000;

    return rc;
  }

  /* getting the runtime version */
  int get_cuda_runtimeVersion(int *runtimeVersion , deviceDetails_t *devD) {
    int rc = RC_SUCCESS;
    cudaError_t error_id;

    error_id = cudaRuntimeGetVersion(runtimeVersion);
    if (error_id != cudaSuccess) { rc = get_error_id (error_id); }

    devD->d_runver = (float)*runtimeVersion / 1000;

    return rc;
  }
  /* getting the name of the device */
  int get_dev_Name(int *idev, deviceDetails_t *devD) {
    int rc = RC_SUCCESS;

    devD->dev_name = deviceProp.name;

    return rc;
  }
  /* get the number of devices avaible on the system */
  int get_dev_count(deviceDetails_t *devD) {
    int rc = RC_SUCCESS; //return code
    int devCnt = 0;

#if defined (CUDA)
    cudaError_t error_id = cudaGetDeviceCount(&devCnt);
    devD->ndev = devCnt;
    //printf(ANSI_COLOR_BRIGHT_GREEN"devD->ndev = %i Devices\n"ANSI_COLOR_RESET
    //       ,devD->ndev);
    if (error_id != cudaSuccess) { rc = get_error_id (error_id); }
#else
    rc = get_warning_message();
#endif
    return rc;
  }
  /* getting the warning message from the cuda */
  int get_ecc_warning_message(deviceDetails_t *devD) {
    int rc = RC_SUCCESS;

    printf(ANSI_COLOR_BRIGHT_CYAN);
    printf("***************************WARNING*****************************\n");
    printf("Device: "ANSI_COLOR_BRIGHT_RED"%s"
           ANSI_COLOR_BRIGHT_CYAN
           " does not have error correcting code \n",devD->dev_name);
    printf("memory (ECC) enabled. This is probably because the device does \n");
    printf("not have CUDA computing capabilities. Check that the correct   \n");
    printf("device has been selected.                                      \n");
    printf("***************************************************************\n");
    printf(ANSI_COLOR_RESET);

    return rc;
  }
  /* getting the warning message from the cuda */
  int get_warning_message() {
    int rc = RC_SUCCESS;

    printf("***************************WARNING*****************************\n");
    printf("You need to compile with -DCUDA to acces the CUDA environment  \n");
    printf("computation using GPU                                          \n");
    printf("***************************************************************\n");
    printf("\n");
    printf("Exit at Line %i in file %s %s\n",__LINE__,__FILE__,__FUNCTION__);
    printf("\n");
    rc = RC_FAIL;

    return rc;
  }
  /* getting the error id from the cuda */
  int get_error_id (cudaError_t error_id) {
    int rc = RC_SUCCESS;
    printf("cudaDriverGetVersion returned %d\n-> %s\n", 
	     (int)error_id, cudaGetErrorString(error_id));
    printf("Result = FAIL\n");
    printf("Exit at Line %i in file %s %s\n",__LINE__,__FILE__,__FUNCTION__);
    rc = (int)error_id;
    exit(EXIT_FAILURE);
    return rc;
  }
  /* getting the cuda cores error */
  int get_cuda_cores_error(int ncores) {
    int rc = RC_SUCCESS;
    printf("There are no CUDA cores available on system %d\n",ncores);
    printf("Result = FAIL\n");
    printf("Exit at Line %i in file %s %s\n",__LINE__,__FILE__,__FUNCTION__);
    rc = RC_FAIL;
    exit(EXIT_FAILURE);
    return rc;
  }/* --end of get_cuda_cores_error--*/
  // Initialization code to find the best CUDA Device
  int findCudaDevice(cudaDeviceProp deviceProp) {
    
    int devID = 0;
    cudaError_t error_id;

    // pick the device with highest Gflops/s
    devID = gpuGetMaxGflopsDeviceId();
    error_id = cudaSetDevice(devID);
    error_id = cudaGetDeviceProperties(&deviceProp, devID);
    /*
    printf("Most suitable GPU Device %d: \"%s\" with compute capability %d.%d\n\n",
           devID, deviceProp.name, deviceProp.major, deviceProp.minor);
    */
    return devID;
  } /* --end of findCudaDevice */
  // This function returns the best GPU (with maximum GFLOPS)
  inline int gpuGetMaxGflopsDeviceId() {
    int current_device     = 0, sm_per_multiproc  = 0;
    int max_perf_device    = 0;
    int device_count       = 0, best_SM_arch      = 0;
    int devices_prohibited = 0;
    
    unsigned long long max_compute_perf = 0;
    cudaDeviceProp deviceProp;
    cudaError_t error_id;

    error_id = cudaGetDeviceCount(&device_count);

    if (device_count == 0) {
      fprintf(stderr, "gpuGetMaxGflopsDeviceId() CUDA error: no devices supporting CUDA.\n");
      exit(EXIT_FAILURE);
    }

    // Find the best major SM Architecture GPU device
    while (current_device < device_count) {
      cudaGetDeviceProperties(&deviceProp, current_device);

      // If this GPU is not running on Compute Mode prohibited,
      // then we can add it to the list
      if (deviceProp.computeMode != cudaComputeModeProhibited) {
        if (deviceProp.major > 0 && deviceProp.major < 9999) {
          best_SM_arch = MAX(best_SM_arch, deviceProp.major);
        }
      } else {devices_prohibited++;}
      current_device++;
    }

    if (devices_prohibited == device_count) {
      fprintf(stderr, "gpuGetMaxGflopsDeviceId() CUDA error: all devices have compute mode prohibited.\n");
      exit(EXIT_FAILURE);
    }

    // Find the best CUDA capable GPU device
    current_device = 0;

    while (current_device < device_count) {
      cudaGetDeviceProperties(&deviceProp, current_device);

      // If this GPU is not running on Compute Mode prohibited, then we can add it to the list
      if (deviceProp.computeMode != cudaComputeModeProhibited) {
        if (deviceProp.major == 9999 && deviceProp.minor == 9999) {
          sm_per_multiproc = 1;
        } else {
          sm_per_multiproc = _ConvertSMVer2Cores(deviceProp.major,
                                                 deviceProp.minor);
        }

        unsigned long long compute_perf  =
          (unsigned long long) deviceProp.multiProcessorCount *
                               sm_per_multiproc               *
                               deviceProp.clockRate;

        if (compute_perf  > max_compute_perf) {
          // If we find GPU with SM major > 2, search only these
          if (best_SM_arch > 2) {
            // If our device==dest_SM_arch, choose this, or else pass
            if (deviceProp.major == best_SM_arch) {
              max_compute_perf  = compute_perf;
              max_perf_device   = current_device;
            }
          } else {
            max_compute_perf  = compute_perf;
            max_perf_device   = current_device;
          }
        }
      }
      ++current_device;
    }

    return max_perf_device;
  } /* --end of gpuGetMaxGflopsDeviceId */
  /* Beginning of GPU Architecture definitions */
  inline int _ConvertSMVer2Cores(int major, int minor) {
    /* Defines for GPU Architecture types (using the SM version to determine
     * the # of cores per SM) 
     * 0xMm (hexidecimal notation),
     * M = SM Major version, and m = SM minor version */
    typedef struct {
      int SM; 
      int Cores;
    } sSMtoCores;

    sSMtoCores nGpuArchCoresPerSM[] = {
      { 0x20, 32 }, /* Fermi Generation (SM 2.0) GF100 class */
      { 0x21, 48 }, /* Fermi Generation (SM 2.1) GF10x class */
      { 0x30, 192}, /* Kepler Generation (SM 3.0) GK10x class */
      { 0x32, 192}, /* Kepler Generation (SM 3.2) GK10x class */
      { 0x35, 192}, /* Kepler Generation (SM 3.5) GK11x class */
      { 0x37, 192}, /* Kepler Generation (SM 3.7) GK21x class */
      { 0x50, 128}, /* Maxwell Generation (SM 5.0) GM10x class */
      { 0x52, 128}, /* Maxwell Generation (SM 5.2) GM20x class */
      {   -1, -1 }
    };
    
    int index = 0;

    while (nGpuArchCoresPerSM[index].SM != -1) {
        if (nGpuArchCoresPerSM[index].SM == ((major << 4) + minor)) {
            return nGpuArchCoresPerSM[index].Cores;
	  }	
        index++;
      }

    /* If we don't find the values, we default use the previous one to run properly */
    printf("MapSMtoCores for SM %d.%d is undefined. Default to use %d Cores/SM\n",
	   major, minor, nGpuArchCoresPerSM[index-1].Cores);
    return nGpuArchCoresPerSM[index-1].Cores;
  }

  /* the aliases for external access */

  extern "C" int setdevice_() __attribute__((weak,alias("setDevice")));
  extern "C" int devicedetails_init_() __attribute__((weak,alias("deviceDetails_init")));
  extern "C" int devicedetails_free_() __attribute__((weak,alias("deviceDetails_free")));

  extern "C" int findcudadevice_() __attribute__((weak,alias("findCudaDevice")));
  extern "C" int get_findcudadevice_() __attribute__((weak,alias("get_findCudaDevice")));

  extern "C" int get_dev_count_() __attribute__((weak,alias("get_dev_count")));
  extern "C" int get_dev_name_() __attribute__((weak,alias("get_dev_Name")));
  extern "C" int get_cuda_runtimeversion_() __attribute__((weak,alias("get_cuda_runtimeVersion")));
  extern "C" int get_cuda_driverversion_() __attribute__((weak,alias("get_cuda_driverVersion")));
  extern "C" int get_tot_global_mem_mb_() __attribute__((weak,alias("get_tot_global_mem_MB")));
  extern "C" int get_tot_global_mem_bytes_() __attribute__((weak,alias("get_tot_global_mem_bytes")));
  extern "C" int get_cuda_cores_() __attribute__((weak,alias("get_CUDA_cores")));
  extern "C" int get_nregisters_() __attribute__((weak,alias("get_nregisters")));
  extern "C" int get_thread_details_() __attribute__((weak,alias("get_thread_details")));
  extern "C" int get_eec_support_() __attribute__((weak,alias("get_eec_support")));
  extern "C" int get_peer_to_peer_capabilities_() __attribute__((weak,alias("get_peer_to_peer_capabilities")));

#endif /* CUDA */

#ifdef __cplusplus
}
#endif
