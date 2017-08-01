! ============================================================================
! Name        : testing_bench_analyser.f90
! Author      : Frederic Bonnet
! Version     :
! Date        : 19th of July 2016
! Description : Analyser code for the gpu benchmarking output
!             :
! ============================================================================
!
program testing_bench_analyser
  use, intrinsic :: iso_c_binding
  use simple_defs
  use simple_err_defs
  use simple_eglossary
  use simple_eglossary_lowlev
  use simple_file_highlev
  use simple_yaml_output
  use simple_yaml_strings
  use simple_file_utils
  use simple_file_defs
  use simple_timing
  use simple_systemQuery_cpu
  use simple_dynamic_memory
  use simple_error_handling
  implicit none

#define devptr_t integer*8
#define verbose .true.
  type(systemQuery_cpu)         :: sysQ
  type(systemDetails)           :: hstD
  type(polar_corr_calc)         :: s_polar

  !local variables
  !error eglossary
  type(eglossary), pointer      :: egloss1
  !invoice test for the eglossary
  real(kind=8) :: price
  type(eglossary), pointer :: egloss, egloss_tmp

  !for the eglossary
  integer :: ival
  !for the Yaml

!*******************************************************************************
!     start of the execution commands
!*******************************************************************************

  !start of the greeting message
  call timestamp()
  call start_Alltimers_cpu()

  call simple_file_lib_initialise

  write(*,*)'                                                                  '
  write(*,*)'******************************************************************'
  write(*,*)'   System fills in the object(sysQ) and the data structure(hstD)  '
  write(*,*)'******************************************************************'
  call sysQ%new_systemQuery_cpu(hstD)
  !call Sanity_check_cpu(sysQ, hstD)
  write(*,*)'******************************************************************'

  write(*,*)'                                                                  '
  write(*,*)'********************* YAML output ********************************'
  write(*,*)'                                                                  '
  
  write(*,*)'                                                                  '
  write(*,*)'********************* eglossary error egloss test*****************'
  write(*,*)'                                                                  '

  write(*,*)'                                                                  '
  write(*,*)'********************* eglossary invoice test**********************'
  write(*,*)'                                                                  '

 
  write(*,*)'******************************************************************'

  !shutting down the timers
  call stop_Alltimers_cpu()

end program testing_bench_analyser
!*******************************************************************************
!    Subroutine to run sanity checks on the data structure passed CPU
!
!*******************************************************************************
!
subroutine Sanity_check_cpu(sysQ, hstD)
  use simple_defs
  use greeting_version
  use simple_systemQuery_cpu
  implicit none
  type(systemQuery_cpu)           :: sysQ
  type(systemDetails)             :: hstD
  !local variables
  !cpu gear
  integer                         :: nCPU_cores
  integer*8                       :: h_TotalMemSize
  integer*8                       :: h_AvailMemSize

  !start of the execution commands

  write(*,*)'                                                                  '
  write(*,*)'******************************************************************'
  write(*,*)'  Sanity checks on the object(sysQ) and the data structure(hstD)  '
  write(*,*)'******************************************************************'
  nCPU_cores = sysQ%get_ncpu_cores()
  h_TotalMemSize = sysQ%get_SysTotalMem_size()
  h_AvailMemSize = sysQ%get_SysAvailMem_size()
  write(*,*)'Number of cores on system     (sysQ):',nCPU_cores
  write(*,*)'Number of cores               (hstD):',hstD%nCPUcores
  if ( nCPU_cores /= hstD%nCPUcores) call sysQ%get_warning_dataStruct()
  write(*,*)'Total Mem on system           (sysQ):',h_TotalMemSize
  write(*,*)'Total Mem on system           (hstD):',hstD%mem_Size
  if(h_TotalMemSize /= hstD%mem_Size)call sysQ%get_warning_dataStruct()
  write(*,*)'Total Available Mem on system (sysQ):',h_AvailMemSize
#if defined (MACOSX)
  write(*,*)'Total Mem on system           (hstD):',hstD%mem_User
  if(h_AvailMemSize /= hstD%mem_User)call sysQ%get_warning_dataStruct()
#elif defined (LINUX)
  write(*,*)'Total Available Mem on system (hstD):',hstD%avail_Mem
  if(h_AvailMemSize /= hstD%avail_Mem)call sysQ%get_warning_dataStruct()
#endif
  
  return
end subroutine Sanity_check_cpu
