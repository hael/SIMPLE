/*
 *   -- SIMPLE addon
 *      Author: Frederic Bonnet, Date: 22nd Apr 2015
 *      Monash University
 *      Apr 2015
 *
 * @precisions normal z -> s d c
 */
/* The System headers */
#include <stdlib.h>
#include <stdio.h>
#include <ctype.h>
#include <string.h>
#include <stddef.h>
#include <math.h>
#include <stdbool.h> //TODO: check inm case it does not compile under MacOsX
#if defined(__GNUC__)
#include <stdint.h>
#endif /* __GNUC__ */

/* mathematical stuff */
#include <complex.h>

/*preprossing for the CUDA environment */
#if defined (CUDA)
#include <cuda.h>
#include <cuda_runtime_api.h>
#include <cublas.h>
#endif

#include <time.h>
#include <stdio.h>
#include <sys/types.h>
#include <sys/time.h>
#include <sys/resource.h>

#if defined (MACOSX)
#include <stdint.h>
#include <sys/param.h>
#include <sys/sysctl.h>
#include <sys/syscall.h>
#elif defined (LINUX)
#include <unistd.h>
#include <linux/param.h>
#include <linux/sysctl.h>
#endif

#ifndef _SIMPLE_H_
#define _SIMPLE_H_

/*some color definiitons for color output*/
#define ANSI_COLOR_BLACK   "\x1b[30m"
#define ANSI_COLOR_RED     "\x1b[31m"
#define ANSI_COLOR_GREEN   "\x1b[32m"
#define ANSI_COLOR_YELLOW  "\x1b[33m"
#define ANSI_COLOR_BLUE    "\x1b[34m"
#define ANSI_COLOR_MAGENTA "\x1b[35m"
#define ANSI_COLOR_CYAN    "\x1b[36m"
#define ANSI_COLOR_WHITE   "\x1b[37m"

#define ANSI_COLOR_BOLD_BLACK   "\x1b[1;30m"
#define ANSI_COLOR_BOLD_RED     "\x1b[1;31m"
#define ANSI_COLOR_BOLD_GREEN   "\x1b[1;32m"
#define ANSI_COLOR_BOLD_YELLOW  "\x1b[1;33m"
#define ANSI_COLOR_BOLD_BLUE    "\x1b[1;34m"
#define ANSI_COLOR_BOLD_MAGENTA "\x1b[1;35m"
#define ANSI_COLOR_BOLD_CYAN    "\x1b[1;36m"
#define ANSI_COLOR_BOLD_WHITE   "\x1b[1;37m"

#define ANSI_COLOR_BRIGHT_BLACK   "\x1b[0;90m"
#define ANSI_COLOR_BRIGHT_RED     "\x1b[0;91m"
#define ANSI_COLOR_BRIGHT_GREEN   "\x1b[0;92m"
#define ANSI_COLOR_BRIGHT_YELLOW  "\x1b[0;93m"
#define ANSI_COLOR_BRIGHT_BLUE    "\x1b[0;94m"
#define ANSI_COLOR_BRIGHT_MAGENTA "\x1b[0;95m"
#define ANSI_COLOR_BRIGHT_CYAN    "\x1b[0;96m"
#define ANSI_COLOR_BRIGHT_WHITE   "\x1b[0;97m"

#define ANSI_COLOR_BOLD_BRIGHT_BLACK   "\x1b[1;90m"
#define ANSI_COLOR_BOLD_BRIGHT_RED     "\x1b[1;91m"
#define ANSI_COLOR_BOLD_BRIGHT_GREEN   "\x1b[1;92m"
#define ANSI_COLOR_BOLD_BRIGHT_YELLOW  "\x1b[1;93m"
#define ANSI_COLOR_BOLD_BRIGHT_BLUE    "\x1b[1;94m"
#define ANSI_COLOR_BOLD_BRIGHT_MAGENTA "\x1b[1;95m"
#define ANSI_COLOR_BOLD_BRIGHT_CYAN    "\x1b[1;96m"
#define ANSI_COLOR_BOLD_BRIGHT_WHITE   "\x1b[1;97m"

#define ANSI_COLOR_RESET   "\x1b[0m"

#define SimpleLeft          'L'
#define SimpleRight         'R'
#define SimpleUpper         'U'
#define SimpleNoTrans       'N'
#define SimpleNonUnit       'N'
#define SimpleLower         'L'
#define SimpleUnit          'U'

/* Size of floating point types (in bytes).*/
#define size_of_float          4
#define size_of_double         8
#define size_of_complex        8
#define size_of_double_complex 16

/* Macro functions */
#define MAX(a,b) (a > b ? a : b)
#define MIN(a,b) (a < b ? a : b)
/* x:0 false 1: true */
#define get_bool(x)(x?"true":"false")

/* maximum number of GPU */
#define MAX_N_GPU            8

/* the return variables */
#define RC_SUCCESS           0
#define RC_FAIL             -1

/* CUBLAS status type returns */
#define  CUBLAS_STATUS_SUCCESS           0
#define  CUBLAS_STATUS_NOT_INITIALIZED   1
#define  CUBLAS_STATUS_ALLOC_FAILED      3
#define  CUBLAS_STATUS_INVALID_VALUE     7
#define  CUBLAS_STATUS_ARCH_MISMATCH     8
#define  CUBLAS_STATUS_MAPPING_ERROR     11
#define  CUBLAS_STATUS_EXECUTION_FAILED  13
#define  CUBLAS_STATUS_INTERNAL_ERROR    14
#define  CUBLAS_STATUS_NOT_SUPPORTED     15
#define  CUBLAS_STATUS_LICENSE_ERROR     1

/* The typedef struct definitions */
typedef size_t devptr_t;
/* system structure for details */
typedef struct systemDetails {
  int n_phys_proc;
  int nCPUcores;
#if defined (MACOSX)
  uint64_t mem_Size;
  int mem_User;
#elif defined (LINUX)
  long long int mem_Size;
  long long int avail_Mem;
#else
  int mem_size;
#endif 
} systemDetails_t;
/* device structure for details */
typedef struct deviceDetails {
  int ndev;
  char *dev_name;
  float d_ver;
  float d_runver;
  float tot_global_mem_MB;
  unsigned long long tot_global_mem_bytes;
  int nMultiProc;
  int ncudacores_per_MultiProc;
  int ncudacores;
  int is_SMsuitable;
  int nregisters_per_blk;
  int warpSze;
  int maxthreads_per_mp;
  int maxthreads_per_blk;
  int is_ecc;
  int is_p2p[MAX_N_GPU][MAX_N_GPU];
} deviceDetails_t;
/* polarFT structure, takers care of the mesh dimensions */
typedef struct polar_corr_calc {
  float r_polar;
  float sumasq_polar;
  float sumbsq_polar;
  int ikrnl;
  int threadsPerBlock;
  int nx,ny,nz;
} polar_corr_calc_t;
/* debugging structure used to substitute the #defined in the .cu files*/
typedef struct debug_gpu {
  int debug_i;         //true
  int debug_cpu_i;     //false
  int debug_high_i;    //true
  int debug_write_i;   //false
  int debug_write_C_i; //false
} debug_gpu_t;
/* benchmarking structure used in additions to DBENCH for more precisions */
typedef struct bench {
  int bench_i;
  int bench_write_i;
} bench_t;
/*data structure for the file_utils and handlers*/
typedef struct resources_avail {
  int nnodes;           //!N of nodes wanted to be used
  int size_socket_logi; //!N of logical cores on each sockets 
  int size_socket_phys; //!N physical cores on each socket
}resources_avail_t;

#ifdef __cplusplus
extern "C" {
#endif

/* ////////////////////////////////////////////////////////////////////////////
   -- Lapack routines used for the block code
*/

/**/
#if defined (CUDA) /* CUDA environment */

/* ////////////////////////////////////////////////////////////////////////////
   -- Error handlers for the cublas on GPU         
*/

  int simple_cudblas_stat_return_c(int err);

#endif /* CUDA */
/**/
#if defined (CUDA) && defined (MAGMA) /* CUDA and MAGMA environment */

/* ////////////////////////////////////////////////////////////////////////////
   -- MAGMA function definitions / Data on CPU
*/

/* //////////////////////////////////////////////////////////////////////////// 
   -- MAGMA function definitions / Data on GPU
*/

/* //////////////////////////////////////////////////////////////////////////// 
 -- MAGMA function definitions / Some getters
*/

#endif /* CUDA and MAGMA */

#ifdef __cplusplus
}
#endif

#undef PRECISION_z
#endif /* _SIMPLE_ */

