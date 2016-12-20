! ============================================================================
! Name        : testing_yaml.f90
! Author      : Frederic Bonnet
! Version     :
! Date        : 11th of February 2016
! Description : tests the yaml output and so on
!             :
! ============================================================================
!
program testing_yaml
  use, intrinsic :: iso_c_binding
  use simple_defs
  use simple_cuda_defs
  use matrixGetter
  use simple_testfunction
  use greeting_version
  use simple_timing
  use simple_systemQuery_cpu
  use simple_yaml_output
  use simple_yaml_strings
  use simple_file_utils
  use simple_file_defs
  use simple_err_defs
  use simple_eglossary
  use simple_eglossary_lowlev
  use simple_dynamic_memory

  implicit none

#define devptr_t integer*8
#define verbose .true.
  type(systemQuery_cpu)         :: sysQ
  type(systemDetails)           :: hstD
  type(polar_corr_calc)         :: s_polar

  integer, parameter            :: maxfls = 3

  integer, parameter            :: nptcls = 100000

  integer, parameter            :: npart=10000
  integer, parameter            :: start_npart=64, istep_npart=1

  integer, parameter            :: nrot =314
  integer, parameter            :: start_nrot=314, istep_nrot=1

  integer, parameter            :: nk = 59
  integer, parameter            :: start_nk=59, istep_nk=1
  !GPU optimastion problem
  integer, parameter            :: nthreads = 256 !number of threads
  !ressources available
  integer, parameter            :: nnodes = 14
  integer, parameter            :: ncores = 16

  integer                       :: err ! error code

  !local variables
  integer                       :: nCPU_cores
  integer                       :: ndev
  integer                       :: N !total number of elemts npart*nrot*nk 
  integer                       :: tcores = ncores*nnodes

  !Yaml imput parameters for the fucntion code
  integer                       :: unit = 1
  !filename string
  integer                       :: startfl
  character(len=3) :: char_out
  character(len=80) :: tmr_name
  !Files handlers
  integer      :: unt
  integer      :: length
  !function calls
  integer      :: get_length_of_string_c
  integer      :: convert_int2char_pos_c
  integer      :: convert_int2char_indexed_c
  integer      :: strlen
  !indexer
  integer      :: i, indx
  !  character(len=*)              :: mapname = 'mapname'
  !  character(len=*)              :: mapvalue
  !  character(len=*)              :: label
  !  character(len=*)              :: tag
  !  character(len=*)              :: advance

!*******************************************************************************
!     start of the execution commands
!*******************************************************************************

  !start of the greeting message
  call timestamp()
  call start_Alltimers_cpu()

  write(*,*)'                                                                  '
  write(*,*)'******************************************************************'
  write(*,*)'   System fills in the object(sysQ) and the data structure(hstD)  '
  write(*,*)'******************************************************************'
  call sysQ%new_systemQuery_cpu(hstD)
  call Sanity_check_cpu(sysQ, hstD)
  write(*,*)'******************************************************************'
  
!*******************************************************************************
!     now testing for the YAML output
!
!*******************************************************************************

  write(*,*)'                                                                  '
  write(*,*)'******************************************************************'
  write(*,*)'     now scanning for size of optimal size and factor of 256      '
  write(*,*)'******************************************************************'
  write(*,*)
  write(*,*)"Sumstack n particlesSize : ",nptcls
  write(*,*)"Number of nodes          : ",nnodes
  write(*,*)"Number of cores          : ",ncores
  write(*,*)"Total number of cores    : ",tcores
  write(*,*)"Size of blocks in threads: ",nthreads
  write(*,*)"Range for npart          : [1,",npart,"]"
  write(*,*)"Range for nrot           : [1,",nrot,"]"
  write(*,*)"Range for nk             : [1,",nk,"]"
  write(*,*)"In steps of istep_npart  :  ",istep_npart
  write(*,*)"In steps of istep_nrot   :  ",istep_nrot
  write(*,*)"In steps of istep_nk     :  ",istep_nk

  write(*,*)'                                                                  '
  write(*,*)'********************* file output ********************************'
  write(*,*)'                                                                  '

  !open(100,file='test.log',status='unknown',action='readwrite',iostat=err)
  
  do i=1,maxfls
     startfl = i
     unit = i
     err = convert_int2char_indexed_c(char_out,i,startfl,maxfls)
     tmr_name = 'yaml_map'
     tmr_name = tmr_name(1:strlen(tmr_name))//char_out(1:strlen(char_out))
     tmr_name = tmr_name(1:strlen(tmr_name))//".yaml"
     call file_open(tmr_name,unit,'unknown','asis','readwrite')
  end do

  write(*,*)'                                                                  '
  write(*,*)'********************* YAML output ********************************'
  write(*,*)'                                                                  '

  call yaml_comment('Now we test for just yaml output')
  call yaml_map('mapname',hstD%nCPUcores)!,fmt='(es24.17)')

  !call yaml_set_stream(tabbing=0)
  call yaml_comment('Yaml Invoice Example',hfill='-')
  call yaml_map('invoice',34843)
  call yaml_map('date',trim(yaml_date_toa()))
  call yaml_map_open('bill-to',label='id001')
   call yaml_map('given','Chris')
   call yaml_map_open('address')
      call yaml_map_open('lines')
      call yaml_scalar('458 Walkman Dr.')
      call yaml_scalar('Suite #292')
      call yaml_map_close()
   call yaml_map_close()
  call yaml_map_close()
  call yaml_map('ship_to','*id001')

  !next step: sequence elements
  call yaml_sequence_open('product')
  !call yaml_sequence_open()
    call yaml_sequence(advance='no')
!    call yaml_mapping_open()
      call yaml_map('sku','BL394D')
      call yaml_map('quantity',4)
      call yaml_map('description','Basketball')
      call yaml_map('price',450.,fmt='(f6.2)')
!    call yaml_mapping_close()
    !call yaml_newline() !new line in a flow 
     call yaml_sequence(advance='no')
!     call yaml_mapping_open()
     call yaml_map('sku','BL4438H')
     call yaml_map('quantity',1)
     call yaml_map('description','Super Hoop')
     call yaml_map('price',2392.,fmt='(f8.2)')
!     call yaml_mapping_close()
    call yaml_sequence_close()
    !final part
    call yaml_map('tax',251.42,fmt='(f6.2)')
    call yaml_map('total',4443.52d0,fmt='(f6.2)') !wrong format on purpose
    call yaml_map('comments','Late afternoon is best. Backup contact is Nancy Billsmer @ 338-4338.')

  write(*,*)'                                                                  '
  write(*,*)'******************************************************************'
  
  !shutting down the timers
  call stop_Alltimers_cpu()

end program testing_yaml
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
  write(*,*)'Total Mem on system           (hstD):',hstD%avail_Mem
  if(h_AvailMemSize /= hstD%avail_Mem)call sysQ%get_warning_dataStruct()
#endif
  
  return
end subroutine Sanity_check_cpu
