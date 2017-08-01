! ============================================================================
! Name        : testing_readWriteImg.f90
! Author      : Frederic Bonnet
! Version     :
! Date        : 01st of July 2016
! Description : tests the read write Imgage from stacks
!             :
! ============================================================================
!
program testing_readWriteImg
  use, intrinsic :: iso_c_binding
  use simple_defs
  use simple_cuda_defs
  use simple_cuda
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
  !image stuff
  use simple_jiffys
  use simple_image
  use simple_imgfile
  use simple_imghead
  use simple_params
  use simple_build
  use simple_cmdline
  !system stuff
  use simple_deviceQuery_gpu
  use simple_systemQuery_cpu
  use simple_commander_base, only: commander_base

  !$ use omp_lib
  !$ use omp_lib_kinds

  implicit none

#define devptr_t integer*8
#define verbose .true.
  type(polar_corr_calc)           :: s_polar
  type(systemQuery_cpu)           :: sysQ
  type(systemDetails)             :: hstD
  type(deviceQuery_gpu)           :: devQ
  type(deviceDetails)             :: devD
  type(params)                    :: p
  type(build)                     :: b
  type(imgfile)                   :: o_imgfile
  type(eglossary), pointer        :: egloss1
  type(cmdline)                   :: cline

  integer                         :: omp_nthr=8  !openMP variables default

  !local variables
  type(deviceDetails),allocatable :: a_devD(:)
  character(len=1),allocatable    :: devname(:,:)
  !image variables
  type(image)                     :: img
  type(imgfile)                   :: ioimg
  type(imghead)                   :: ioimghead
  integer                         :: nimg
  character(len=10)               :: imgkind
  real, allocatable               :: vec(:)
  !output files
  integer                         :: unit_vec = 2
  !filename string
  character(len=3)                :: char_out
  character(len=80)               :: tmr_name
  !gpu gear
  integer                         :: ndev
  !return code and err catchers
  integer                         :: err
  integer                         :: rc
  !indexers
  integer                         :: idev,iptcls,ivec_sz

  !functions calls
  integer :: get_dev_count_c
  integer :: findCudaDevice_c
  integer :: setDevice_c
  integer :: strlen
  integer :: convert_int2char_indexed_c

!*******************************************************************************
!     start of the execution commands, first the command line
!*******************************************************************************
  if( command_argument_count() < 2 )then
     write(*,'(a)',advance='yes') 'Usage: '
     write(*,'(a)',advance='no') 'testing_readWriteImg stk=<data1> crf_name1=<core_file_name1> [stk2=<data2>]'
     write(*,'(a)') ' [use_gpu=<{yes/no|no}>] [nthr=<nr of OpenMP threads{1}>]'
     stop
  endif
  call cline%parse
  call cline%checkvar('stk',     1)
  call cline%checkvar('crf_name1',  2)
  call cline%check
  !getting the parameters
  p=params(cline)
!*******************************************************************************
!     Initialising the environments
!*******************************************************************************
  !start of the greeting message
  call timestamp()
  call start_Alltimers_cpu()
  !initialising all of the initialisers high level for the file environment
  !call simple_file_lib_initialise
  !starting the cuda environment
  call simple_cuda_init(err,devQ,devD)
  if (err .ne. RC_SUCCESS ) write(*,*) 'cublas init failed'
!*******************************************************************************
!    Building the environments for image processing under SIMPLE
!
!*******************************************************************************
  ! general objects built
  call b%build_general_tbox(p,cline)
  call img%new([p%box,p%box,1],p%smpd,imgkind='em')
  imgkind = p%imgkind
  if ( omp_nthr /= p%nthr) then
     write(*,*) "The omp_nthr /= p%nthr from the command line!"
     write(*,*) "Default value for omp_nthr: ",omp_nthr
     write(*,*) "omp_nthr: ",omp_nthr," p%nthr: ",p%nthr
     write(*,*) "aligning the two for coherence..."
     omp_nthr = p%nthr
  end if
  !starting the openMP environemt and setting the number of threads to compute
  !$ call omp_set_num_threads(omp_nthr)

!*******************************************************************************
!    Getting Hardware information and setting up the data Structures using
!    the query modules 
!*******************************************************************************
  write(*,*)'                                                                  '
  write(*,*)'******************************************************************'
  write(*,*)'   System fills in the object(sysQ) and the data structure(hstD)  '
  write(*,*)'******************************************************************'
  call sysQ%new_systemQuery_cpu(hstD)
  !call Sanity_check_cpu(sysQ, hstD)
!*******************************************************************************
!     Printing some info on screen
!
!*******************************************************************************
  call find_ldim_nptcls(p%stk, p%ldim, p%nptcls)
  write(*,*)'                                                                  '
  write(*,*)'******************************************************************'
  write(*,*)'     Printing some info on screen                                 '
  write(*,*)'******************************************************************'
  write(*,*)
  write(*,*) "Stack of particles      : ",p%stk
  write(*,*) "Dimension of stack      : ",p%ldim
  write(*,*) "Sampling distance       : ",p%smpd
  write(*,*) "Image kind              : ",imgkind!(strlen(1:imgkind))
  write(*,*) "N images in stack (nimg): ",p%nptcls
  write(*,*) "OpenMP thrd             : ",omp_nthr
  write(*,*)'                                                                  '
  write(*,*)'******************************************************************'
  write(*,*)'                                                                  '
!*******************************************************************************
!    Loading up the data from stack and writting out to vector
!
!*******************************************************************************
  err = convert_int2char_indexed_c(char_out,1,p%nptcls,p%nptcls)
  tmr_name = p%crf_name1  !'stack_vec'
  !tmr_name = tmr_name(1:strlen(tmr_name))//char_out(1:strlen(char_out))
  tmr_name = tmr_name(1:strlen(tmr_name))//".asc"
  call file_open(tmr_name,unit_vec,'unknown','asis','readwrite')
  nimg = p%nptcls
!  nimg = 3
  do iptcls = 1,nimg
     call progress(iptcls,p%nptcls)
     call img%read(p%stk,iptcls)
     call img%serialize(pcavec=vec)
     do ivec_sz = 1, size(vec)-1
        write(unit_vec,'(f15.8)',advance='no')vec(ivec_sz)
     end do
     write(unit_vec,'(f15.8)',advance='yes')vec(size(vec))
     deallocate(vec)
  end do

  write(*,*) size(vec)

!*******************************************************************************
!    Starting the tests of the modules readwrite and file_utils           
!
!*******************************************************************************  
  write(*,*)'                                                                  '
  write(*,*)'********************* YAML output ********************************'
  write(*,*)'                                                                  '

  !call yaml_comment('Now we test the eglossary inside yaml')

  !call yaml_map('mapname',hstD%nCPUcores,fmt='(es24.17)')
  write(*,*)'                                                                  '

  write(*,*)'                                                                  '
  write(*,*)'********************* eglossary test******************************'
  write(*,*)'                                                                  '

  !egloss1=> egloss_new(egloss1)

  write(*,*)'******************************************************************'

!*******************************************************************************
!    Freeing resopurces                             
!
!*******************************************************************************

!*******************************************************************************
!    Shuting down the environments and finalisers
!
!*******************************************************************************

  !shutting down the environment
  call simple_cuda_shutdown()
  !shutting down the timers
  call stop_Alltimers_cpu()
  !end of greeting message
  call bye_deviceQuery_gpu()

end program testing_readWriteImg
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
