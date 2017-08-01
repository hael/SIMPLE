! ============================================================================
! Name        : simple_BenchMrkDtAnlzr.f90
! Author      : Frederic Bonnet
! Version     :
! Date        : 27th of July 2016
! Description : program to perform the data analysis from the BenchMarking runs
!             :
! ============================================================================
!
program simple_BenchMrkDtAnlzr
  use, intrinsic :: iso_c_binding
  use simple_defs
  use simple_cuda_defs
  use simple_cuda
  use matrixGetter
  use simple_testfunction
  use greeting_version
  use simple_timing
  use simple_systemQuery_cpu
  use simple_deviceQuery_gpu
  use simple_math, only: calc_corr, csq
  use simple_yaml_output
  use simple_yaml_strings
  use simple_file_utils
  use simple_file_defs
  use simple_err_defs
  use simple_eglossary
  use simple_eglossary_lowlev
  use simple_dynamic_memory
  use simple_params
  use simple_cmdline
  use simple_commander_base, only: commander_base
  implicit none
  include 'BenchMrkDtAnlzr-incFile.f90'
#define devptr_t integer*8
#define verbose .true.
  type(systemQuery_cpu)         :: sysQ
  type(systemDetails)           :: hstD
  type(polar_corr_calc)         :: s_polar
  type(params)                  :: p
  type(cmdline)                 :: cline

!  integer                       :: maxits = 8
!  integer   !in incFile         :: nchunk=1
!  integer                       :: nlp=10
!  integer                       :: nmsk=2
!  integer, parameter            :: nset = 2      #
  
  !statistical data
  integer,allocatable           ::     d_msk(:,:)
  integer,allocatable           ::   d_chunk(:,:)
  real(sp),allocatable          ::      d_lp(:,:)
  integer,allocatable           ::    d_nrot(:,:,:,:,:)
  integer,allocatable           ::      d_nk(:,:,:,:,:)
  real(sp),allocatable          ::    d_time(:,:,:,:,:)
  real(sp),allocatable          :: d_speedup(:,:,:,:)
  !local varaibales
  integer                       :: err ! error code
  !Yaml imput parameters for the fucntion code
  integer                       :: unit_yaml = 1
  integer                       :: unit_1 = 2
  integer                       :: unit_2 = 3
  integer                       :: unit_3 = 4
  integer                       :: unit_4 = 5
  !filename string
  character(len=3)  :: char_out
  character(len=80) :: wrk_name
  character(len=80) :: tmr_name_yaml
  character(len=80) :: tmr_name_1,tmr_name_2
  character(len=80) :: out_name_1,out_name_2
  !indexers
  integer :: imsk,ilp,ichunk
  !function calls
  integer      :: get_length_of_string_c
  integer      :: convert_int2char_pos_c
  integer      :: convert_int2char_indexed_c
  integer      :: strlen

!  TODO: further integration in module will require an interface
!  interface get_d_ave
!     module procedure get_d_ave_int , get_d_ave_real
!  end interface get_d_ave

  !start of the execution commands

!*******************************************************************************
!    Environment initilisers
!
!*******************************************************************************
  !start of the greeting message
  call hello_gpu_magma()
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
!    Using the command line object from simple
!
!*******************************************************************************
  if( command_argument_count() < 1 )then
     write(*,'(a)',advance='yes') 'Usage: '
     write(*,'(a)',advance='no') 'mcBenchMrkDtAnlzr crf_name1=<core_file_name1>'
     write(*,'(a)') ' crf_name2=<core_file_name2>'
     stop
  endif
  call cline%parse
  call cline%checkvar('crf_name1', 1)
  call cline%checkvar('crf_name2', 2)
! call cmdcheckvar('use_gpu', 2)
! call cmdcheckvar('nthr',    3)
  call cline%check
  p = params(cline)
  write(*,*) "input core filename1: ",p%crf_name1
  write(*,*) "input core filename2: ",p%crf_name2
!*******************************************************************************
!    Starting the YAML output environment     
!
!*******************************************************************************
  err = convert_int2char_indexed_c(char_out,0,1,1)
  tmr_name_yaml = 'benchMarking_yaml_map'
  tmr_name_yaml = tmr_name_yaml(1:strlen(tmr_name_yaml))//&
       char_out(1:strlen(char_out))
  tmr_name_yaml = tmr_name_yaml(1:strlen(tmr_name_yaml))//".yaml"
  call file_open(tmr_name_yaml,unit_yaml,'unknown','asis','readwrite')
  call yaml_comment('Now we output the benchmarking results in yaml file')
  call yaml_map('mapname',hstD%nCPUcores)!,fmt='(es24.17)')
!*******************************************************************************
!    Data allocation and initialisation
!
!*******************************************************************************
  !Data
  !making all arrays congruent, the entry 0 and maxits+1 in d_time
  !are used for the average and error based on simple std for now
  allocate( d_nrot(1:nset,1:nchunk,1:nlp,1:nmsk,0:maxits+1))
  allocate(  d_nk(1:nset,1:nchunk,1:nlp,1:nmsk,0:maxits+1))
  allocate(d_time(1:nset,1:nchunk,1:nlp,1:nmsk,0:maxits+1))
  d_nrot  = 0
  d_nk   = 0
  d_time = 0.0
  allocate( d_msk(1:nset,1:nmsk))
  allocate(  d_chunk(1:nset,1:nchunk))
  allocate(d_lp(1:nset,1:nlp))
  d_msk = 0
  d_chunk = 0
  d_lp = 0.0
  allocate( d_speedup(1:nchunk,1:nlp,1:nmsk,0:1) )
!*******************************************************************************
!    Data collection
!
!*******************************************************************************
 !Opening the input files
  !tmr_name_1 = 'gpu_no_bench_yes_fix_no_set_0_z_prime3D'
  tmr_name_1 = p%crf_name1(1:strlen(p%crf_name1))
  wrk_name = tmr_name_1(1:strlen(tmr_name_1))//".asc"
  call file_open(wrk_name,unit_1,'unknown','asis','read')
  tmr_name_2 = p%crf_name2(1:strlen(p%crf_name2)) !'gpu_yes_bench_no_fix_no_set_0_z_prime3D'
  wrk_name = tmr_name_2(1:strlen(tmr_name_2))//".asc"
  call file_open(wrk_name,unit_2,'unknown','asis','read')

  !reading in the data
  call readInData(tmr_name_1,unit_1, &
                   d_msk(1,1:nmsk),d_chunk(1,1:nchunk),d_lp(1,1:nlp),&
                  d_nrot(1,1:nchunk,1:nlp,1:nmsk,0:maxits+1),&
                    d_nk(1,1:nchunk,1:nlp,1:nmsk,0:maxits+1),&
                  d_time(1,1:nchunk,1:nlp,1:nmsk,0:maxits+1),&
                  maxits,nset,nchunk,nlp,nmsk)
  !reading in the data
  call readInData(tmr_name_2,unit_2,&
                   d_msk(2,1:nmsk),d_chunk(2,1:nchunk),d_lp(2,1:nlp),&
                  d_nrot(2,1:nchunk,1:nlp,1:nmsk,0:maxits+1),&
                    d_nk(2,1:nchunk,1:nlp,1:nmsk,0:maxits+1),&
                  d_time(2,1:nchunk,1:nlp,1:nmsk,0:maxits+1),&
                  maxits,nset,nchunk,nlp,nmsk)
!*******************************************************************************
!    Data processing
!
!*******************************************************************************
  !getting the sample mean
  call get_d_ave_real(d_time,maxits,nset,nchunk,nlp,nmsk)
  call get_d_ave_int(d_nrot,maxits,nset,nchunk,nlp,nmsk)
  call get_d_ave_int(d_nk,maxits,nset,nchunk,nlp,nmsk)
  !getting the std
  call get_d_std_real(d_time,maxits,nset,nchunk,nlp,nmsk)
  !getting the delat in times measurments
  call get_d_deltaBar_real(d_time,maxits,nset,nchunk,nlp,nmsk)
  !testing code
  do ichunk=1,nchunk
     do ilp=1,nlp
        do imsk=1,nmsk
           !internal loop of 8 iterations
           write(2500+unit_1,*) &
                   d_lp(1,ilp),&!lp
                  d_msk(1,imsk), &!,ring2,
                d_chunk(1,ichunk),&!nspace,&
                 d_nrot(1,ichunk,ilp,imsk,0),&!nrot,&
                   d_nk(1,ichunk,ilp,imsk,0),&!nk,&
                 d_time(1,ichunk,ilp,imsk,0),&!elps_r,&
                 d_time(1,ichunk,ilp,imsk,maxits+1)!elps_r,&
        end do
     end do
  end do
  write(*,*) d_nk(1,1,:,1,0)
  write(*,*) d_lp
  write(*,*) d_lp(1,1)
  !getting the speedup factors between the data sets
  call get_speedup(1,2,d_time,d_speedup)
!*******************************************************************************
!    writing data to disk
!
!*******************************************************************************
  !nrot and nk vs timeBar(err) for each data sets
  call write_nrotnk_time_err(tmr_name_1,unit_3,&
                             d_nrot(1,:,:,:,:),&
                               d_nk(1,:,:,:,:),&
                             d_time(1,:,:,:,:),&
                             maxits,nset,nchunk,nlp,nmsk)
  call write_nrotnk_time_err(tmr_name_2,unit_4,&
                             d_nrot(2,:,:,:,:),&
                               d_nk(2,:,:,:,:),&
                             d_time(2,:,:,:,:),&
                             maxits,nset,nchunk,nlp,nmsk)
  !nk vs timeBar(err) for each data sets
  call write_nk_time_err(tmr_name_1,unit_3,&
                           d_nk(1,:,:,:,:),&
                         d_time(1,:,:,:,:),&
                         maxits,nset,nchunk,nlp,nmsk)
  call write_nk_time_err(tmr_name_2,unit_4,&
                           d_nk(2,:,:,:,:),&
                         d_time(2,:,:,:,:),&
                         maxits,nset,nchunk,nlp,nmsk)
!*******************************************************************************
!    Data deallocation
!
!*******************************************************************************
  !making all arrays congruent, the entry 0 and maxits+1 in d_time
  !are used for the average and error based on simple std for now
  deallocate(   d_nrot)
  deallocate(     d_nk)
  deallocate(   d_time)
  deallocate(    d_msk)
  deallocate(  d_chunk)
  deallocate(     d_lp)
  deallocate(d_speedup)
!*******************************************************************************
!    Environment shutdown
!
!*******************************************************************************
  !shutting down the timers
  call stop_Alltimers_cpu()

  !end of greeting message
  call bye_gpu_magma()

contains
  
end program simple_BenchMrkDtAnlzr
!//////////////////// Helper subroutine and functions //////////////////////////
!//

!*******************************************************************************
!    Subroutine to get the speed up factors between the data            
!
!*******************************************************************************
!integer arrays
subroutine get_ratio(d_data)
  use simple_defs
  implicit none
  include 'BenchMrkDtAnlzr-incFile.f90'
  real(sp) :: d_data(1:nset,1:nchunk,1:nlp,1:nmsk,0:maxits+1)

  return
end subroutine get_ratio
!*******************************************************************************
!    Subroutine to get the speed up factors between the data            
!
!*******************************************************************************
!integer arrays
subroutine get_speedup(idata_1,idata_2,d_data,d_speedup)
  use simple_defs
  implicit none
  include 'BenchMrkDtAnlzr-incFile.f90'
  integer  :: idata_1, idata_2 
  real(sp) :: d_data(1:nset,1:nchunk,1:nlp,1:nmsk,0:maxits+1)
  real(sp) :: d_speedup(1:nchunk,1:nlp,1:nmsk,0:1)
  !local variables

  !start of the execution commands

  d_speedup(:,:,:,0) = d_data(idata_1,:,:,:,0) / d_data(idata_2,:,:,:,0)

  !getting the error bar based in the fractional error annalysis

  d_speedup(:,:,:,1) = d_speedup(:,:,:,0) * &
       sqrt( (d_data(idata_1,:,:,:,maxits+1)/d_data(idata_1,:,:,:,0))**2 + &
             (d_data(idata_2,:,:,:,maxits+1)/d_data(idata_2,:,:,:,0))**2   )

  return
end subroutine get_speedup
!*******************************************************************************
!    Subroutine to print to file the data in question for the 2 and 3D data
!    for a given data set
!*******************************************************************************
!just nrot and nk vs time(err)
subroutine write_nrotnk_time_err(tmr_name,unit,d_nrot,d_nk,d_time,&
                                 maxits,nset,nchunk,nlp,nmsk)
  use simple_defs
  use simple_file_utils
  implicit none
  integer           :: nset,nchunk,nlp,nmsk,maxits
  character(len=80) :: tmr_name
  integer           :: unit
  integer           ::  d_nrot(1:nchunk,1:nlp,1:nmsk,0:maxits+1)
  integer           ::    d_nk(1:nchunk,1:nlp,1:nmsk,0:maxits+1)
  real(sp)          ::  d_time(1:nchunk,1:nlp,1:nmsk,0:maxits+1)
  !local variables
  integer           :: err ! error code
  character(len=3)  :: char_out
  character(len=80) :: wrk_name
  character(len=80) :: out_name
  !indexers
  integer           :: imsk,ilp,ichunk
  !function calls  
  integer           :: get_length_of_string_c
  integer           :: convert_int2char_pos_c
  integer           :: convert_int2char_indexed_c
  integer           :: strlen

  out_name = tmr_name(1:strlen(tmr_name))//"_nrotnk"
  err = convert_int2char_indexed_c(char_out,0,1,1)
  wrk_name = out_name(1:strlen(out_name))//char_out(1:strlen(char_out))
  wrk_name = wrk_name(1:strlen(wrk_name))//".asc"
  call file_open(wrk_name,unit,'unknown','asis','readwrite')
  do ichunk=1,nchunk
     do ilp=1,nlp
        do imsk=1,nmsk
           !internal loop of 8 iterations
           write(unit,'(i3,a,i3,a,f12.4,a,f12.8)') &
                d_nrot(ichunk,ilp,imsk,0),",",&!nrot,&
                d_nk(ichunk,ilp,imsk,0),",",&!nk,&
                d_time(ichunk,ilp,imsk,0),",",&!elps_r,&
                d_time(ichunk,ilp,imsk,maxits+1)!delta elps_r,&
        end do
     end do
  end do

  return
end subroutine write_nrotnk_time_err
!just nk vs time(err)
subroutine write_nk_time_err(tmr_name,unit,d_nk,d_time,&
                                 maxits,nset,nchunk,nlp,nmsk)
  use simple_defs
  use simple_file_utils
  implicit none
  integer           :: nset,nchunk,nlp,nmsk,maxits
  character(len=80) :: tmr_name
  integer           :: unit
  integer           ::    d_nk(1:nchunk,1:nlp,1:nmsk,0:maxits+1)
  real(sp)          ::  d_time(1:nchunk,1:nlp,1:nmsk,0:maxits+1)
  !local variables
  integer           :: err ! error code
  character(len=3)  :: char_out
  character(len=80) :: wrk_name
  character(len=80) :: out_name
  !indexers
  integer           :: imsk,ilp,ichunk
  !function calls  
  integer           :: get_length_of_string_c
  integer           :: convert_int2char_pos_c
  integer           :: convert_int2char_indexed_c
  integer           :: strlen

  out_name = tmr_name(1:strlen(tmr_name))//"_nk"
  err = convert_int2char_indexed_c(char_out,0,1,1)
  wrk_name = out_name(1:strlen(out_name))//char_out(1:strlen(char_out))
  wrk_name = wrk_name(1:strlen(wrk_name))//".asc"
  call file_open(wrk_name,unit,'unknown','asis','readwrite')
  do ichunk=1,nchunk
     do ilp=1,nlp
        do imsk=1,nmsk
           !internal loop of 8 iterations
           write(unit,'(i3,a,f12.4,a,f12.8)') &
                d_nk(ichunk,ilp,imsk,0)," ",&!nk,&
                d_time(ichunk,ilp,imsk,0)," ",&!elps_r,&
                d_time(ichunk,ilp,imsk,maxits+1)!delat elps_r,&
        end do
     end do
  end do

  return
end subroutine write_nk_time_err
!*******************************************************************************
!    Subroutine to get the averages and error bars from the read in data
!
!*******************************************************************************
!integer arrays
subroutine get_d_ave_int(d_data,maxits,nset,nchunk,nlp,nmsk)
  use simple_defs
  implicit none
  integer           :: maxits
  integer           :: nset,nchunk,nlp,nmsk
  integer           :: d_data(1:nset,1:nchunk,1:nlp,1:nmsk,0:maxits+1)

  d_data(:,:,:,:,0) = sum(d_data,dim=5) / maxits

  return
end subroutine get_d_ave_int
!real arrays
subroutine get_d_ave_real(d_data,maxits,nset,nchunk,nlp,nmsk)
  use simple_defs
  implicit none
  integer           :: maxits
  integer           :: nset,nchunk,nlp,nmsk
  real(sp)          :: d_data(1:nset,1:nchunk,1:nlp,1:nmsk,0:maxits+1)

  d_data(:,:,:,:,0) = sum(d_data,dim=5) / maxits

  return
end subroutine get_d_ave_real
!real arrays
subroutine get_d_std_real(d_data,maxits,nset,nchunk,nlp,nmsk)
  use simple_defs
  implicit none
  integer           :: maxits
  integer           :: nset,nchunk,nlp,nmsk
  real(sp)          :: d_data(1:nset,1:nchunk,1:nlp,1:nmsk,0:maxits+1)
  !local variables
  !indexers
  integer           :: imaxits
  if (maxits > 1 )then
     do imaxits=1,maxits
        d_data(:,:,:,:,maxits+1) = d_data(:,:,:,:,maxits+1) + &
                                  (d_data(:,:,:,:,imaxits)-d_data(:,:,:,:,0))**2
     end do
  end if

  d_data(:,:,:,:,maxits+1) = sqrt((1/(maxits-1.0))*d_data(:,:,:,:,maxits+1))

  return
end subroutine get_d_std_real
!real arrays
subroutine get_d_deltaBar_real(d_data,maxits,nset,nchunk,nlp,nmsk)
  use simple_defs
  implicit none
  integer           :: maxits
  integer           :: nset,nchunk,nlp,nmsk
  real(sp)          :: d_data(1:nset,1:nchunk,1:nlp,1:nmsk,0:maxits+1)
  !local variables
  !indexers
  integer           :: imaxits

  d_data(:,:,:,:,maxits+1) = d_data(:,:,:,:,maxits+1)*sqrt(1/(maxits*1.0))

  return
end subroutine get_d_deltaBar_real

!*******************************************************************************
!    Subroutine to readin the data
!
!*******************************************************************************
subroutine readInData(tmr_name,unit,&
                      d_msk,d_chunk,d_lp,&
                      d_nrot,d_nk,d_time,&
                      maxits,nset,nchunk,nlp,nmsk)
  use simple_defs
  implicit none
  character(len=80) :: tmr_name
  integer           :: unit,maxits
  integer           :: nset,nchunk,nlp,nmsk
  integer           ::   d_msk(1:nmsk)
  integer           :: d_chunk(1:nchunk)
  real(sp)          ::    d_lp(1:nlp)
  integer           ::  d_nrot(1:nchunk,1:nlp,1:nmsk,*)
  integer           ::    d_nk(1:nchunk,1:nlp,1:nmsk,*)
  real(sp)          ::  d_time(1:nchunk,1:nlp,1:nmsk,*)
  !local varaiables
  character(len=3)  :: gpu
  character(len=5)  :: bench
  character(len=3)  :: fix
  character(len=3)  :: set
  character(len=3)  :: gpu_yn, bench_yn, fix_yn
  integer           :: set_id
  real(sp)          :: smpd,lp
  integer           :: ring2,nspace,nrot,nk
  !reading input in from generated file
  character(len=7) :: msg1
  character(len=4) :: msg2
  character(len=3) :: msg3
  character(len=12) :: msg4
  real(sp)          :: elps_r
  character(len=4)  :: t_unit
  character(len=20) :: molecule
  !indexers
  integer :: imaxits,imsk,ilp,ichunk
  !function call
  integer :: strlen

  do ichunk=1,nchunk
     do ilp=1,nlp
        do imsk=1,nmsk
           !internal loop of 8 iterations
           do imaxits=1,maxits
              !write(*,*) imaxits
              read(unit,*) &
                   gpu,   gpu_yn,&
                   bench, bench_yn,&
                   fix,   fix_yn, &
                   set,   set_id, &
                   smpd,lp,ring2,nspace,nrot,nk,&
                   msg1,msg2,msg3,msg4,elps_r,t_unit,&
                   molecule
              !filling the arrays
              d_lp(ilp) = lp
              d_msk(imsk) = ring2
              d_chunk(ichunk) = nspace
              d_nrot(ichunk,ilp,imsk,imaxits)  = nrot 
              d_nk(ichunk,ilp,imsk,imaxits)   = nk
              d_time(ichunk,ilp,imsk,imaxits) = elps_r
              write(1500+unit,'(7(a,x),i1,f7.4,x,f8.3,x,i5,x,i5,x,i5,x,i3,x,4(a,x),f12.4,2(x,a))') &
                   gpu, gpu_yn(1:strlen(gpu_yn)), &
                   bench, bench_yn(1:strlen(bench_yn)),&
                   fix,fix_yn(1:strlen(fix_yn)), &
                   set, set_id,&
                   smpd,&
                   d_lp(ilp),&!lp
                   d_msk(imsk), &!,ring2,
                   d_chunk(ichunk),&!nspace,&
                   d_nrot(ichunk,ilp,imsk,imaxits),&!nrot,&
                   d_nk(ichunk,ilp,imsk,imaxits),&!nk,&
                   msg1,msg2,msg3,msg4,&
                   d_time(ichunk,ilp,imsk,imaxits),&!elps_r,&
                   t_unit, molecule(1:strlen(molecule))
           end do

        end do
     end do
  end do
  
  return
end subroutine readInData
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
