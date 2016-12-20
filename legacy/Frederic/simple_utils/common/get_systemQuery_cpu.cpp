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

#include "get_systemQuery_cpu.h"

#ifdef __cplusplus
extern "C" {
#endif
/***************************************************************************/
/**
 *  FORTRAN API - math functions (simple interface)
 **/
  float convert_memory(long long int mem) {
    int rc = RC_SUCCESS;
    float mem_conv = 0.0f;
    mem_conv = mem / (1024 * 1024);

    if (mem_conv == 0 ) { rc = RC_FAIL;}
    if (rc == RC_FAIL) { rc = get_error_cpu(); }

    return mem_conv;
  }

  /* determine the physical usable memory size on Linux */
  int get_memorySize_host(long long int *memSize, systemDetails_t *hstD) {
    int rc = RC_SUCCESS;
#if defined (LINUX)
    *memSize = 0;
    *memSize = sysconf(_SC_PHYS_PAGES) * sysconf(_SC_PAGESIZE);
    hstD->mem_Size = *memSize;
    printf(ANSI_COLOR_BRIGHT_BLUE"Sys total memory: "
           ANSI_COLOR_BRIGHT_YELLOW"%6.2f MBytes,"
           ANSI_COLOR_BRIGHT_GREEN" %lli bytes. \n"
           ANSI_COLOR_RESET,
      convert_memory(hstD->mem_Size), hstD->mem_Size);
#else
    rc = get_warning_message_linux_cpu();
#endif

    if (rc == RC_FAIL) { rc = get_error_cpu(); }
    return rc;
  }
  /* determine the available memory size on Linux */
  int get_available_memory_host(long long int *availMem,
				systemDetails_t *hstD) {
    int rc = RC_SUCCESS;
#if defined (LINUX)
    *availMem = 0;
    *availMem = sysconf(_SC_AVPHYS_PAGES) * sysconf(_SC_PAGESIZE);
    hstD->avail_Mem = *availMem;
    printf(ANSI_COLOR_BRIGHT_BLUE "Available memory: "
           ANSI_COLOR_BRIGHT_YELLOW"%6.2f MBytes,"
           ANSI_COLOR_BRIGHT_GREEN" %lli bytes\n"
           ANSI_COLOR_RESET,
      convert_memory(hstD->avail_Mem), hstD->avail_Mem);
#else
    rc = get_warning_message_linux_cpu();
#endif

    if (rc == RC_FAIL) { rc = get_error_cpu(); }
    return rc;
  }
  /* determine the number of cores on the host on Linux */
  int get_CPU_cores(int *nCPU_cores, systemDetails_t *hstD) {
    int rc = RC_SUCCESS;
#if defined (LINUX)
    *nCPU_cores = 0;
    *nCPU_cores = sysconf(_SC_NPROCESSORS_ONLN);
    hstD->nCPUcores = *nCPU_cores;
    if (hstD->nCPUcores < 1 ) {rc = RC_FAIL;}
    printf(ANSI_COLOR_BRIGHT_GREEN"CPU_cores: %i\n"ANSI_COLOR_RESET,
           hstD->nCPUcores);
#else
    rc = get_warning_message_linux_cpu();
#endif
    if (rc == RC_FAIL) { rc = get_error_cpu(); }
    return rc;
  }
  /* getting the warning message from the cuda */
  int get_warning_message_linux_cpu() {
    int rc = RC_SUCCESS;

    printf("***************************WARNING*****************************\n");
    printf("You need to compile with -DLINUX to acces the linux environment\n");
    printf("operating system                                               \n");
    printf("***************************************************************\n");
    printf("\n");
    printf("Exit at Line %i in file %s %s\n",__LINE__,__FILE__,__FUNCTION__);
    printf("\n");
    rc = RC_FAIL;

    return rc;
  }
  /*getting the error message and return code */
  int get_error_cpu() {
    int rc = RC_SUCCESS;
    printf("rc  = %i, at Line %i in file %s %s\n",
	   rc,__LINE__,__FILE__,__FUNCTION__);
    printf("Result = FAIL\n");
    rc = RC_FAIL;
    exit(EXIT_FAILURE);
    return rc;
  }

  /* the aliases for external access */
  extern "C" int get_cpu_cores_() __attribute__((weak,alias("get_CPU_cores")));
  extern "C" int get_memorysize_host_() __attribute__((weak,alias("get_memorySize_host")));
  extern "C" int get_available_memory_host_() __attribute__((weak,alias("get_available_memory_host")));

#ifdef __cplusplus
}
#endif
