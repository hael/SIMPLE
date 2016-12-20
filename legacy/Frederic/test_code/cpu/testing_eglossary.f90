! ============================================================================
! Name        : testing_eglossary.f90
! Author      : Frederic Bonnet
! Version     :
! Date        : 01st of July 2016
! Description : tests the glossary module and the yaml module
!             :
! ============================================================================
!
program testing_eglossary
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
  call Sanity_check_cpu(sysQ, hstD)
  write(*,*)'******************************************************************'

  write(*,*)'                                                                  '
  write(*,*)'********************* YAML output ********************************'
  write(*,*)'                                                                  '

  !memory_limit=
  !output_level=
  !logfile_name=
  call file_malloc_set_status(memory_limit=real(hstD%mem_Size),&
       output_level=2,logfile_name='memstatus.yaml')

  call yaml_comment('Now we test the eglossary inside yaml')
  call yaml_warning('warning test',8)
!  call yaml_map('mapname',hstD%nCPUcores,fmt='(es24.17)')
  call yaml_map('mapname',hstD%nCPUcores)
  
  write(*,*)'                                                                  '
  write(*,*)'********************* eglossary error egloss test*****************'
  write(*,*)'                                                                  '
  !Creating the new eglossary
  call yaml_comment('Now we test for just eglossary')
  egloss1=> egloss_new()
  call file_err_open_try()
  ival = egloss1//'Fredo'
  write(*,*) "ival: ",ival

  call yaml_map('ival not existing, fake value',ival)

  write(*,*)'                                                                  '
  write(*,*)'********************* eglossary invoice test**********************'
  write(*,*)'                                                                  '

  call yaml_comment('Yaml Invoice Example, using eglossaries',hfill='-')
  !setting the data in the fortran eglossary
  call egloss_init(egloss)
  call set(egloss//'invoice',34843)
  call set(egloss//'date',trim(yaml_date_toa()))

  call set(egloss//'bill-to'//'given','Chris')
  call set(egloss//'bill-to'//'family','Dumars')
  
  call egloss_init(egloss_tmp)
    call set(egloss_tmp//'lines','458 Walkman Dr. Suite #292')
    call set(egloss_tmp//'city','Royal Oak')
    call set(egloss_tmp//'state','MI')
    call set(egloss_tmp//'postal',48046)

  call set(egloss//'bill-to'//'address',egloss_tmp)
  !no need to free the eglossary after association
  
  !tagging of eglossary not yet implemented

  !products
  call egloss_init(egloss_tmp)
  call set(egloss_tmp//'sku','BL34D')
  call set(egloss_tmp//'quantity',4)
  call set(egloss_tmp//'description','Basketball')
  call set(egloss_tmp//'price',450.00)
  !adding to the item
  call add(egloss//'Product',egloss_tmp)

  call egloss_init(egloss_tmp)
  call set(egloss_tmp//'sku','BL4438')
  call set(egloss_tmp//'quantity',1)
  call set(egloss_tmp//'description','Super Hoop')
  call set(egloss_tmp//'price',2392.00)
  !TODO: need to fix the mapping value to real(8)
  price=egloss_tmp//'price'
  call yaml_map('Retrieve the price value',price)
  call add(egloss//'Product',egloss_tmp)

  call set(egloss//'Tax',251.42)
  call set(egloss//'Total',4443.52)
  call set(egloss//'Comments','Late afternoon is best. Backup contact is Nancy Billsmer @ 338-4338')

  !print invoice
  call yaml_egloss_write(egloss)

  call egloss_free(egloss)
  
  write(*,*)'******************************************************************'

  !shutting down the timers
  call stop_Alltimers_cpu()

end program testing_eglossary
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
