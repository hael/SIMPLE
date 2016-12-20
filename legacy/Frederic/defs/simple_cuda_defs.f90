!*******************************************************************************
!     Author: Frederic D.R. Bonnet date: 4th of March 2015.
!
! Name:
! simple_cuda_defs - basic definitions for cuda used in all modules.
!
! Description:
! simple_cuda_defs provides basic definitions for the types and declarations
! used in gpu calculations in modules using cuda calls. Using cuda-5.0
!*******************************************************************************
!
module simple_cuda_defs
  use, intrinsic :: iso_c_binding
  
  use simple_defs

  implicit none

!******
!DESCRIPTION
!defining the types of the for the GPU environment
#define devptr_t integer*8
#define size_t integer*8
#if defined (OPENCL)
#define cl_mem integer*8
#endif
  ! the return variables */
  integer, parameter :: RC_SUCCESS =  0
  integer, parameter :: RC_FAIL    = -1
  !polarFT data structure for the GPU environment
  !ikrnl          : selects the kernel in consideratrion
  !threadsPerBlock: number of threads in shared memory on GPU 1D mesh
  !nx,ny,nz       : Mesh diemensions in the 3D mesh  
  type, bind(c) :: polar_corr_calc
     real(c_float)  :: r_polar
     real(c_float)  :: sumasq_polar
     real(c_float)  :: sumbsq_polar
     integer(c_int) :: ikrnl
     integer(c_int) :: threadsPerBlock
     integer(c_int) :: nx,ny,nz
  end type polar_corr_calc
  !debugging derived type for the GPU environment to be passed in the kernels
  ! initialised to 0 and set in the code in question 0: false 1: true
  type, bind(c) :: t_debug_gpu
     integer(c_int) ::         debug_i = 1 !true/false
     integer(c_int) ::     debug_cpu_i = 0 !false/true
     integer(c_int) ::    debug_high_i = 1 !true/false
     integer(c_int) ::   debug_write_i = 0 !false/true
     integer(c_int) :: debug_write_C_i = 0 !false/true
  end type t_debug_gpu
  !Benchmarking type definition
  type, bind(c) :: t_bench
     integer(c_int) :: bench_i = 0
     integer(c_int) :: bench_write_i = 0
  end type t_bench

  !GPU parameters and settings
  integer, parameter :: R_POLAR_D = 1.321                 !kernel options
  integer, parameter :: SUMASQ_POLAR_D = 2.132            !kernel options
  integer, parameter :: SUMBSQ_POLAR_D = 3.123            !kernel options
  integer, parameter :: KERNEL_D = 6                      !kernel options
  integer, parameter :: THREADSPERBLOCK_D = 256           !kernel options
  integer, parameter :: NX_3D_D=16, NY_3D_D=16, NZ_3D_D=4 !3D mesh
  !logical variable (.true.) to determine if GPU(s) are available on the system.
  !.false. otherwise
  logical :: has_gpu = .false.
  logical :: has_multi_gpu = .false.
  logical :: use_gpu = .false.
  !logical variables (.true.) if BenchMarking is done. (.false.) otherwise
  logical :: ibench=.false.      !< benchmark result or not
  logical :: ibench_write=.false.!< write benchmark result or not
  ! The different kinds of precision
  !S: single D: double C: single complex Z: double complex
  integer, parameter :: kind_S = 1
  integer, parameter :: kind_D = 2
  integer, parameter :: kind_C = 3
  integer, parameter :: kind_Z = 4
  ! Size of floating point types (in bytes).
  integer, parameter :: size_of_real = 4
  integer, parameter :: size_of_double = 8
  integer, parameter :: size_of_complex = 8
  integer, parameter :: size_of_double_complex = 16
  !maximum number of GPU per nodes set
  integer, parameter :: MAX_N_GPU = 8
  integer, parameter :: DEVNAME_STRING_LENGTH = 80
  !Supported architecture set to 3.0 and above
  integer, parameter :: GPU_CARD_MAJOR = 3
  integer, parameter :: GPU_CARD_MINOR = 0

  !The c binding derived type for the device structure
  type, bind(c) :: deviceDetails
     integer(c_int)       :: ndev
     type(c_ptr)          :: dev_name
     real(c_float)        :: d_ver
     real(c_float)        :: d_runver
     real(c_float)        :: tot_global_mem_MB
     integer(c_long_long) :: tot_global_mem_bytes
     integer(c_int)       :: nMultiProc
     integer(c_int)       :: ncc_per_mp
     integer(c_int)       :: ncc
     integer(c_int)       :: is_SMsuitable
     integer(c_int)       :: nregisters_per_blk
     integer(c_int)       :: warpSze
     integer(c_int)       :: maxthreads_per_mp
     integer(c_int)       :: maxthreads_per_blk
     integer(c_int)       :: is_ecc
     integer(c_int)       :: is_p2p(0:MAX_N_GPU-1,0:MAX_N_GPU-1) = 999
  end type deviceDetails

  !fast fourier transform direction
  integer, parameter :: FFTW_FORWARD = -1
  integer, parameter :: FFTW_BACKWARD = +1

  integer, parameter :: CUFFT_FORWARD = -1
  integer, parameter :: CUFFT_INVERSE = +1
  !Type of transform single and double precision
  integer, parameter :: FFT_C2C = 11
  integer, parameter :: FFT_S2C = 12
  integer, parameter :: FFT_C2S = 13

  integer, parameter :: FFT_Z2Z = 14
  integer, parameter :: FFT_D2Z = 15
  integer, parameter :: FFT_Z2D = 16

#if defined (CUDA)
  integer, external :: cublas_init
  integer, external :: cublas_shutdown
  integer, external :: cublas_alloc         !memory allocation on GPU
  integer, external :: cublas_free          !memory deallocation on GPU

  integer,external :: cublas_set_vector     !setting the vector on the GPU
  integer,external :: cublas_get_vector     !getting the vector from the GPU

  integer,external :: cublas_set_matrix     !setting the matrix on the GPU
  integer,external :: cublas_get_matrix     !getting the matrix from the GPU

  integer,external :: simple_cublas_free    !simple version of cublas free synchronised
  integer,external :: simple_cublas_alloc   !simple version of cublas alloc synchronised

  integer,external :: cublas_get_error      !getting the error code from the cublas calls
  integer,external :: cublas_xerbla         !getting the error code from the xerbla message

  !* CUBLAS status type returns *!
  integer, parameter  ::    CUBLAS_STATUS_SUCCESS          = 0
  integer, parameter  ::    CUBLAS_STATUS_NOT_INITIALIZED  = 1
  integer, parameter  ::    CUBLAS_STATUS_ALLOC_FAILED     = 3
  integer, parameter  ::    CUBLAS_STATUS_INVALID_VALUE    = 7
  integer, parameter  ::    CUBLAS_STATUS_ARCH_MISMATCH    = 8
  integer, parameter  ::    CUBLAS_STATUS_MAPPING_ERROR    = 11
  integer, parameter  ::    CUBLAS_STATUS_EXECUTION_FAILED = 13
  integer, parameter  ::    CUBLAS_STATUS_INTERNAL_ERROR   = 14
  integer, parameter  ::    CUBLAS_STATUS_NOT_SUPPORTED    = 15
  integer, parameter  ::    CUBLAS_STATUS_LICENSE_ERROR    = 1

#endif

contains
!
! Error message from the cublas calls methods
!
  subroutine simple_cudblas_stat_return(err)
    implicit none

    !global varaiables
    integer :: err

    !start of the executions commands

#if defined (CUDA)
    
    select case(err)

    case( CUBLAS_STATUS_SUCCESS          )
       write(*,*) 'Error=',err,': CUBLAS_STATUS_SUCCESS'
    case( CUBLAS_STATUS_NOT_INITIALIZED  )
       write(*,*) 'Error=',err,': CUBLAS_STATUS_NOT_INITIALIZED'
       write(*,*) 'Error=',err,': CUBLAS_STATUS_LICENSE_ERROR'
    case( CUBLAS_STATUS_ALLOC_FAILED     )
       write(*,*) 'Error: device memory for devPtrA allocation failed'
       write(*,*) 'Error=',err,': CUBLAS_STATUS_ALLOC_FAILED'
    case( CUBLAS_STATUS_INVALID_VALUE    )
       write(*,*) 'Error=',err,': CUBLAS_STATUS_INVALID_VALUE'
    case( CUBLAS_STATUS_ARCH_MISMATCH    )
       write(*,*) 'Error=',err,': CUBLAS_STATUS_ARCH_MISMATCH'
    case( CUBLAS_STATUS_MAPPING_ERROR    )
       write(*,*) 'Error: Data upload on the GPU failed'
       write(*,*) 'Error=',err,': CUBLAS_STATUS_MAPPING_ERROR'
    case( CUBLAS_STATUS_EXECUTION_FAILED ) 
       write(*,*) 'Error=',err,': CUBLAS_STATUS_EXECUTION_FAILED'
    case( CUBLAS_STATUS_INTERNAL_ERROR   )
       write(*,*) 'Error=',err,': CUBLAS_STATUS_INTERNAL_ERROR'
    case( CUBLAS_STATUS_NOT_SUPPORTED    )
       write(*,*) 'Error=',err,': CUBLAS_STATUS_NOT_SUPPORTED'

    end select
#else
    err = -1
    write(*,*)"***************************WARNING******************************"
    write(*,*)"You need to compile with -DCUDA to acces the CUDA environment   "
    write(*,*)"computation using GPU                                           "
    write(*,*)"****************************************************************"
#endif

    return
  end subroutine simple_cudblas_stat_return

end module simple_cuda_defs
