! ============================================================================
! Name        : simple_systemQuery_cpu.f90
! Author      : Frederic Bonnet
! Version     :
! Date        : 27th of May 2015
! Description : To obtain infromation from the system on CPU mother boards
!             :
! ============================================================================
!
program testing_systemQuery_cpu

  use simple_defs
  use greeting_version
  use simple_timing
  use simple_systemQuery_cpu

  implicit none

  type(systemQuery_cpu)           :: sysQ
  type(systemDetails)             :: hstD
  integer                         :: err
  !local variables
  integer                         :: nCPU_cores
  integer*8                       :: h_TotalMemSize
  integer*8                       :: h_AvailMemSize
  !start of the execution commands
  !start of the greeting message
  call hello_systemQuery_cpu(err)
  call timestamp()
  call start_Alltimers_cpu()

  call sysQ%new_systemQuery_cpu(hstD)

  write(*,*)'                                                                  '
  write(*,*)'******************************************************************'
  write(*,*)'   Sanity checks on the object(sysQ) and the data structure(hstD) '
  write(*,*)'******************************************************************'
  nCPU_cores = sysQ%get_ncpu_cores()
  h_TotalMemSize = sysQ%get_SysTotalMem_size()
  h_AvailMemSize = sysQ%get_SysAvailMem_size()
  write(*,*)'Number of cores on the system                    : ',nCPU_cores
  !writting the number of cores avaible on the system
  write(*,*)'Number of cores returned from data structure hstD: ',hstD%nCPUcores
  if ( nCPU_cores /= hstD%nCPUcores) call sysQ%get_warning_dataStruct()

  write(*,*)'Total Mem on the system                          : ',h_TotalMemSize
  !writting the number of cores avaible on the system
  write(*,*)'Total Mem on the system returned from struct hstD: ',hstD%mem_Size
  if(h_TotalMemSize /= hstD%mem_Size)call sysQ%get_warning_dataStruct()

  write(*,*)'Total Available Mem on the system                : ',h_AvailMemSize
  !writting the number of cores avaible on the system
#if defined (MACOSX)
  write(*,*)'Total Mem on the system returned from struct hstD: ',hstD%mem_User
  if(h_AvailMemSize /= hstD%mem_User)call sysQ%get_warning_dataStruct()
#elif defined (LINUX)
  write(*,*)'Total Mem on the system returned from struct hstD: ',hstD%avail_Mem
  if(h_AvailMemSize /= hstD%avail_Mem)call sysQ%get_warning_dataStruct()
#endif
  write(*,*)'******************************************************************'

  !shutting down the timers
  call stop_Alltimers_cpu()

  !end of greeting message
  call bye_systemQuery_cpu()

end program testing_systemQuery_cpu
