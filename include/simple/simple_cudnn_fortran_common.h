
/* Define CUBLAS_FORTRAN_COMPILER for Windows needed because 
   tbe compilation is done from the perl script fortran_nightly.pl which
   does not include fortran_tools.mk
   On Linux and Apple, CFLAGS is setup with fortran_tools.mk
*/
#if defined(_WIN32)
#ifndef CUBLAS_INTEL_FORTRAN
#define CUBLAS_INTEL_FORTRAN
#endif
#endif

#if defined(CUBLAS_GFORTRAN)
/* using option -ff2c make the ABI compatible with F77
   No need to define RETURN_COMPLEX, which cause problem 
   on Gfortran 4.x on 32 bit
*/
/* #define RETURN_COMPLEX 1 */
#endif

#if defined(CUBLAS_GFORTRAN) || defined (CUBLAS_G95)
/* NOTE: Must use -fno-second-underscore when building Fortran source with g77
 *       g77 invocation may not use -fno-f2c, which forces different return 
 *       type conventions than the one used below
 */
#define SIMPLE_CUDNN_CREATE             simple_cudnn_create_
#define SIMPLE_CUDNN_DESTROY            simple_cudnn_destroy_

#elif defined(CUBLAS_INTEL_FORTRAN)

#define SIMPLE_CUDNN_CREATE             simple_cudnn_create
#define SIMPLE_CUDNN_DESTROY            simple_cudnn_destroy

#else
#error unsupported Fortran compiler
#endif
