!*******************************************************************************
!     Author: Frederic D.R. Bonnet date: 4th of March 2015.
!
! Name:
! simple_OpenCL_defs - basic definitions for OpenCL used in all modules.
!
! Description:
! simple_OpenCL_defs provides basic definitions for the types and declarations
! used in gpu calculations in modules using OpenCL calls. Using OpenCL
! Interface for OpenCl parallelization.
! Defines a few global variables (see below) and provides subroutines :
!  - Start the log files
!  - Defines the platform ID
!  - Gets the platform details
!  - get the devices 
!  - Defines and creates the context 
!  - Defines and creates the command queues
!  - Defines and creates the build of the opencl program
!
!*******************************************************************************
!
module simple_OpenCL_defs
  
  use simple_defs

  implicit none

  interface
#if defined (OPENCL)
     subroutine startlogs_gpu(filenamelog)
       implicit none
       character(len=80) :: filenamelog
     end subroutine startlogs_gpu

     subroutine platformid_gpu()
       implicit none
     end subroutine platformid_gpu

     subroutine getplatformdetails_gpu()
       implicit none
     end subroutine getplatformdetails_gpu

     subroutine getdevices_gpu()
       implicit none
     end subroutine getdevices_gpu

     subroutine createcontext_gpu()
       implicit none
     end subroutine createcontext_gpu

     subroutine createcommandqueue_gpu()
       implicit none
     end subroutine createcommandqueue_gpu

     subroutine createbuild_gpu(path, file)
       implicit none
       character(len=80) :: path, file
     end subroutine createbuild_gpu

     subroutine createcmdevvector_gpu(vector)
       use simple_defs
       implicit none
       real(sp) :: vector(*)
     end subroutine createcmdevvector_gpu

     subroutine routine_sax_gpu( n, alpha, incx )
       implicit none
       integer :: n
       real :: alpha
       integer :: incx
     end subroutine routine_sax_gpu

     subroutine getGPU_interface_gpu( nplat )
       implicit none
       integer :: nplat
     end subroutine getGPU_interface_gpu
#endif
  end interface

contains
!
! starts gpu initializing using OpenCL
!
  subroutine hello__gpu_OpenCL(string)
    implicit none

    !global variables
    
    character(len=80)        :: string

    !start of the execution commands

    write(*,*)
    write(*,*) "hello GPU OpenCL world"
    write(*,*)

  end subroutine hello__gpu_OpenCL

end module simple_OpenCL_defs
