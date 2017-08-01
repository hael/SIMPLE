!****h* TB_Sim/bonnet/modules/defs
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! NAME
! magma_defs.f90 - Various test to test the CUBLAS routines
!
! AUTHOR
! Frederic D.R. Bonnet
!
! COPYRIGHT
! Copyright (C) 2005-2010 Commissariat a l'Energie Atomique (CEA) de Grenoble,
! France. All rights reserved.
!
! DESCRIPTION
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
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! SOURCE

program testing_openCL_defs

  use simple_OpenCL_defs

  implicit none

  integer :: nplat = 1
  character(len=80) :: temp

  ! Main code.

  temp = "string"
  call hello(temp)

#if defined (OPENCL) /*preprossing for the OPENCL environment */

  call getGPU_interface_gpu( nplat ) !getPlatformDetails_interface( nplat )

#endif /* OPENCL */

end program testing_openCL_defs
