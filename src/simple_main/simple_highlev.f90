! 
!*******************************************************************************
!     Author: Frederic D.R. Bonnet date: 29th of July 2016.
!
! Name:
! simple_highlev - High level routine to drive the library.
!
!*******************************************************************************
!
module simple_highlev
  use simple_defs                ! singleton
  use simple_cmdline             ! singleton
  use simple_jiffys,             only: simple_end
  use simple_hadamard3D_matcher, only: prime3D_exec, prime3D_find_resrange
  use simple_params,             only: params
  use simple_build,              only: build
  use simple_timing
  use simple_cuda_defs
  use simple_cuda
  use simple_file_utils
  use simple_file_defs
  !temporary
  !use simple_err_defs
  use simple_file_highlev
  !use simple_eglossary
  !use simple_error_handling
  !use simple_dynamic_memory
  !use simple_systemQuery_cpu
  !use simple_deviceQuery_gpu
  implicit none

contains
  !TODO: insert the higher level methods
end module simple_highlev
