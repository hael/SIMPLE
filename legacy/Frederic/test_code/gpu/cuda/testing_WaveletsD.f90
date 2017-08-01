! ============================================================================
! Name        : testing_WaveletsD.f90
! Author      : Frederic Bonnet
! Version     :
! Date        : 10th of August 2016
! Description : tests the wavelets implementation
!             :
! ============================================================================
!
program testing_WaveletsD
  use, intrinsic :: iso_c_binding
  use simple_defs
  use simple_cuda_defs
  use simple_cuda
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
  use simple_deviceQuery_gpu
  use fft123D_gpu

  implicit none

#define devptr_t integer*8
#define verbose .true.
  type(systemQuery_cpu)         :: sysQ
  type(systemDetails)           :: hstD
  type(deviceQuery_gpu)         :: devQ
  type(deviceDetails)           :: devD
  type(polar_corr_calc)         :: s_polar
  type(cuFFT_gpu)               :: t_cuFFT_gpu

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

  !functional details
  real(dp) :: xi,xf
  complex(dp) :: ui,uf
  integer  :: nx
  real(dp) :: a,b,c,na
  real(dp) :: delta
  complex(dp),allocatable :: u(:),fu(:)
  complex(dp),allocatable :: fh(:)
  !wavelet parameters
  integer      :: J
  !Yaml imput parameters for the fucntion code
  integer      :: unit = 1
  integer      :: unit_2 = 2
  integer      :: unit_3 = 3
  integer      :: unit_4 = 4
  !filename string
  integer      :: startfl
  character(len=3) :: char_out
  character(len=80) :: tmr_name
  !Files handlers
  integer      :: unt
  integer      :: length
  !function calls
  integer      :: get_wavelet_coef_c
  integer      :: test_transform_sig_1d_sym_c
  integer      :: test_transform_sig_2d_sym_c
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
  call hello_deviceQuery_gpu(err)
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
!     start of the CUDA environment  
!*******************************************************************************

  !starting the cuda environment
  call simple_cuda_init(err)
  if (err .ne. 0 ) write(*,*) 'cublas init failed'

  call devQ%new_deviceQuery_gpu(devD)

  ndev = devQ%get_devCnt()
  write(*,*)'cudaGetDeviceCount returned: ',ndev
  !writting the number of devices from data structure
  write(*,*)'cudaGetDeviceCount returned from data structure: ',devD%ndev
  if ( ndev /= devD%ndev) write(*,*)"data structure and getter do not match"

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

  unit = 1
  err = convert_int2char_indexed_c(char_out,unit,1,1)
  tmr_name = 'wavelets_yaml_map'
  tmr_name = tmr_name(1:strlen(tmr_name))//char_out(1:strlen(char_out))
  tmr_name = tmr_name(1:strlen(tmr_name))//".yaml"
  call file_open(tmr_name,unit,'unknown','asis','readwrite')

  write(*,*)'                                                                  '
  write(*,*)'********************* YAML output ********************************'
  write(*,*)'                                                                  '

  call yaml_comment('Now we test for just wavelet yaml output')
  call yaml_map('mapname',hstD%nCPUcores)!,fmt='(es24.17)')

  write(*,*)'                                                                  '
  write(*,*)'********* Constructing the functions and getting FT **************'
  write(*,*)'                                                                  '

  !real function 
  !test function details
  xi = -10
  xf = 10
  a = 2
  b = 2
  c = 1
  na = 1
  nx = 10000
  call yaml_map('x_i',xi)!,fmt='(es24.17)')
  call yaml_map('x_f',xf)!,fmt='(es24.17)')
  call yaml_map('a',a)!,fmt='(es24.17)')
  call yaml_map('b',b)!,fmt='(es24.17)')
  call yaml_map('c',c)!,fmt='(es24.17)')
  call yaml_map('power of n',na)!,fmt='(es24.17)')
  call yaml_map('Number of points in the function nx',nx)!,fmt='(es24.17)')
  !call t_cuFFT_gpu%new_cuFFT_1D_ZZ_gpu(nx)

  allocate(u(nx))
  allocate(fu(nx))
  allocate(fh(nx))

  !complex function
  ui = cmplx(-10.0, (4.0d0 * atan(1.0d0)/2.0d0)**(1/b) )
  uf = cmplx( 10.0, (4.0d0 * atan(1.0d0)/2.0d0)**(1/b) )
  fu = atom_like_complex(ui,uf,nx,a,b,c,u,na)
  !TODO: need to implement the list yaml output
  !call yaml_map('First 10 elemts of fu',real(fu(1:10)))!,fmt='(es24.17)')

  !now moving to the fast 1D fourier transform on GPU
  call hello_cuFFT_gpu(err)
  fh = 0.0d0
  call t_cuFFT_gpu%gather_fft_1D_Z2Z_gpu(nx,fu,fh,CUFFT_FORWARD)

  if ( nx >= 5 ) then
     write(*,*) "the first 5 entries of:"
     write(*,*) "u         data:f(u)     Fourier: f(h)"
     do i=1,nx-(nx-5)
        write(*,*)u(i), fu(i), fh(i)
     end do
  end if

  err = convert_int2char_indexed_c(char_out,unit,1,1)
  tmr_name = 'fu_real'
  tmr_name = tmr_name(1:strlen(tmr_name))//char_out(1:strlen(char_out))
  tmr_name = tmr_name(1:strlen(tmr_name))//".asc"
  call file_open(tmr_name,unit_2,'unknown','asis','readwrite')

  err = convert_int2char_indexed_c(char_out,unit,1,1)
  tmr_name = 'fh_real'
  tmr_name = tmr_name(1:strlen(tmr_name))//char_out(1:strlen(char_out))
  tmr_name = tmr_name(1:strlen(tmr_name))//".asc"
  call file_open(tmr_name,unit_3,'unknown','asis','readwrite')

  err = convert_int2char_indexed_c(char_out,unit,1,1)
  tmr_name = 'fh_imag'
  tmr_name = tmr_name(1:strlen(tmr_name))//char_out(1:strlen(char_out))
  tmr_name = tmr_name(1:strlen(tmr_name))//".asc"
  call file_open(tmr_name,unit_4,'unknown','asis','readwrite')
  
  !printing function to file for plotting
  if ( nx >= 5 ) then
     write(*,*) "the first 3x(i,j) = 9 entries of:"
     write(*,'(3(11x,a),4(15x,a))')"i","j","u","v","data:f(u,v)","Fourier: f(h,p)","Inverse Fourier"
     do i=1,nx
           write(unit_2,*)real(fu(i)) !real(u(i)),
           write(unit_3,*)real(u(i)),real(fh(i))
           write(unit_4,*)real(u(i)),imag(fh(i))
     end do
  end if

  write(*,*)'                                                                  '
  write(*,*)'********* Starting the wavelet analysis **************************'
  write(*,*)'                                                                  '

  err = get_wavelet_coef_c("db4")
  J = 6 !decompistion factor
  err = test_transform_sig_1d_sym_c("db4", J)
  J = 3
  err = test_transform_sig_2d_sym_c("db4", J)
  write(*,*)'                                                                  '
  write(*,*)'******************************************************************'

!*******************************************************************************
!    Freeing the resources                                                 
!
!*******************************************************************************
!
  deallocate(u)
  deallocate(fu)
  deallocate(fh)
!*******************************************************************************
!    Shuting down the environment                                     
!
!*******************************************************************************
!
  !shutting down the timers
  call stop_Alltimers_cpu()
  call bye_cuFFT_gpu()
  !shutting down the environment
  call simple_cuda_shutdown()

end program testing_WaveletsD
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
