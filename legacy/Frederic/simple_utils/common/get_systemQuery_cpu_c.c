/*
 *   -- SIMPLE addon
 *      Author: Frederic Bonnet, Date: 27st May 2015
 *
 *      May 2015
 *
 *      c code for the MACOSX code and LINUX to extract system information 
 *
 * @precisions normal z -> s d c
 */

/* The Simple header */
#include "simple.h"

/* determine the physical usable memory size on MacOSXC */
#if defined (MACOSX)
int get_memorysize_host_c_(uint64_t *memSize, systemDetails_t *hstD) {
  int rc = RC_SUCCESS;
  int nm[2];
  size_t len;
  uint64_t count;
  nm[0] = CTL_HW;      /* generic cpu/io */
  nm[1] = HW_MEMSIZE;   /* uint64_t: physical ram size */
  len = sizeof(uint64_t);
  printf(ANSI_COLOR_BRIGHT_BLUE"len: "
         ANSI_COLOR_BRIGHT_YELLOW"%zd\n"ANSI_COLOR_RESET,sizeof(uint64_t) );
  sysctl(nm,2,&count,&len,NULL,0);
  *memSize = count;
  hstD->mem_Size = *memSize;
  printf(ANSI_COLOR_BRIGHT_BLUE"Memory Size "
         ANSI_COLOR_BRIGHT_GREEN"%lli\n"ANSI_COLOR_RESET,hstD->mem_Size);

  if (rc == RC_FAIL) { rc = get_error_cpu_c(); }
  return rc;
}
void GET_MEMORYSIZE_HOST_C_() __attribute__((weak,alias("get_memorysize_host_c_")));
void get_memorysize_host_c__() __attribute__((weak,alias("get_memorysize_host_c_")));
void GET_MEMORYSIZE_HOST_C__() __attribute__((weak,alias("get_memorysize_host_c_")));
#endif
/* getting the user memeory */
int get_usermem_host_c_(int *userMem, systemDetails_t *hstD) {
  int rc = RC_SUCCESS;
  int nm[2];
  size_t len;
  int count;
  len = sizeof(int);
  *userMem = 0;
#if defined (MACOSX)
  nm[0] = CTL_HW;      /* generic cpu/io */
  nm[1] = HW_USERMEM;   /* int: non-kernel memory */
  sysctl(nm,2,&count,&len,NULL,0);
  *userMem = count;
  hstD->mem_User = *userMem;
  printf("User mem %i\n",hstD->mem_User);
#else
  rc = get_warning_message_macosx_cpu();
#endif
    
  if (rc == RC_FAIL) { rc = get_error_cpu_c(); }
  return rc;
}
void GET_USERMEM_HOST_C_() __attribute__((weak,alias("get_usermem_host_c_")));
void get_usermem_host_c__() __attribute__((weak,alias("get_usermem_host_c_")));
void GET_USERMEM_HOST_C__() __attribute__((weak,alias("get_usermem_host_c_")));
/* determine the number of cores on the host on MacOSX */
int get_cpu_cores_c_(int *nCPU_cores, systemDetails_t *hstD) {
  int rc = RC_SUCCESS;
  int nm[2];
  size_t len = 4;
  uint32_t count;
#if defined (MACOSX)
  nm[0] = CTL_HW;      /* generic cpu/io */
  nm[1] = HW_AVAILCPU; /* int: number of available CPUs */
  //nm[2] = HW_NCPU;     /* int: number of cpus */
  sysctl(nm,2,&count,&len,NULL,0);
  *nCPU_cores = (int)count;
  hstD->nCPUcores = *nCPU_cores;
  printf("Memory Size %lli\n",hstD->nCPUcores);
#else
    rc = get_warning_message_macosx_cpu();
#endif

  if (rc == RC_FAIL) { rc = get_error_cpu_c(); }
  return rc;
}
void GET_CPU_CORES_C_() __attribute__((weak,alias("get_cpu_cores_c_")));
void get_cpu_cores_c__() __attribute__((weak,alias("get_cpu_cores_c_")));
void GET_CPU_CORES_C__() __attribute__((weak,alias("get_cpu_cores_c_")));
/* getting the warning message from the cuda */
int get_warning_message_macosx_cpu() {
  int rc = RC_SUCCESS;
  printf("***************************WARNING*****************************\n");
  printf("You need to compile with -DMACOSX to acces the MacOSX          \n");
  printf("environment operating system                                   \n");
  printf("***************************************************************\n");
  printf("\n");
  printf("Exit at Line %i in file %s %s\n",__LINE__,__FILE__,__FUNCTION__);
  printf("\n");
  rc = RC_FAIL;
  return rc;
}
/*getting the error message and return code */
int get_error_cpu_c() {
  int rc = RC_SUCCESS;
  printf("rc  = %i, at Line %i in file %s %s\n",
	 rc,__LINE__,__FILE__,__FUNCTION__);
  printf("Result = FAIL\n");
  rc = RC_FAIL;
  exit(EXIT_FAILURE);
  return rc;
}
/* Other return varaibles that may be considered later on */
//#define HW_MACHINE       1     /* string: machine class */
//#define HW_MODEL         2     /* string: specific machine model */
//#define HW_PHYSMEM       5     /* int: total memory */
//#define HW_DISKNAMES     8     /* strings: disk drive names */
//#define HW_DISKSTATS     9     /* struct: diskstats[] */
//#define HW_CPU_FREQ     15     /* int: CPU Frequency */
//#define HW_MEMSIZE      24     /* uint64_t: physical ram size */
//#define HW_AVAILCPU     25     /* int: number of available CPUs */
