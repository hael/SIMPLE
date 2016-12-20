! ============================================================================
! Name        : testing_cuFFT_1D_gpu.f90
! Author      : Frederic Bonnet
! Version     :
! Date        : 23rd of June 2015
! Description : tests the cuFFT library for 1D FFT on cuda
!             :
! ============================================================================
!
program testing_cuFFT_1D_Z2Z_gpu

  use simple_defs
  use simple_cuda_defs
  use simple_cuda
  use greeting_version
  use simple_timing
  use simple_deviceQuery_gpu
  use fft123D_gpu
  use simple_testfunction

  use simple_file_utils
  use simple_file_defs

  implicit none
#define devptr_t integer*8

  type(deviceQuery_gpu)           :: devQ
  type(deviceDetails)             :: devD
  type(cuFFT_gpu)                 :: t_cuFFT_gpu
  integer                         :: err
  !local variables
  integer                         :: ndev

  real(dp) :: xi,xf
  complex(dp) :: ui,uf
  integer  :: nx
  real(dp) :: a,b,c,n
  real(dp) :: delta
  complex(dp),allocatable :: u(:),fu(:)
  complex(dp),allocatable :: fh(:)
  complex(dp),allocatable :: fu_gpu(:)
  !function calls
  integer      :: get_length_of_string_c
  integer      :: convert_int2char_pos_c
  integer      :: convert_int2char_indexed_c
  integer      :: strlen
  !Yaml imput parameters for the fucntion code
  integer                       :: unit = 1
  !filename string
  character(len=3) :: char_out
  character(len=80) :: tmr_name
  !counters
  integer  :: i,j

  !start of the execution commands
  !start of the greeting message
  call hello_deviceQuery_gpu(err)
  call timestamp()
  call start_Alltimers_cpu()

#if defined (CUDA)

  !starting the cuda environment
  call simple_cuda_init(err)
  if (err .ne. 0 ) write(*,*) 'cublas init failed'

  call devQ%new_deviceQuery_gpu(devD)

  ndev = devQ%get_devCnt()
  write(*,*)'cudaGetDeviceCount returned: ',ndev
  !writting the number of devices from data structure
  write(*,*)'cudaGetDeviceCount returned from data structure: ',devD%ndev
  if ( ndev /= devD%ndev) write(*,*)"data structure and getter do not match"

  !now moving to the fast 1D fourier transform on GPU
  call hello_cuFFT_gpu(err)

  !real function 
  !test function details
  xi = -10
  xf = 10
  a = 2
  b = 2
  c = 1
  n = 1
  nx = 10000

  !call t_cuFFT_gpu%new_cuFFT_1D_ZZ_gpu(nx)

  allocate(u(nx))
  allocate(fu(nx))
  allocate(fh(nx))
  allocate(fu_gpu(nx))

  !complex function
  ui = cmplx(-10.0, (4.0d0 * atan(1.0d0)/2.0d0)**(1/b) )
  uf = cmplx( 10.0, (4.0d0 * atan(1.0d0)/2.0d0)**(1/b) )
  fu = atom_like_complex(ui,uf,nx,a,b,c,u,n)
  !call t_cuFFT_gpu%set_Zdata_1D(nx,u,fu)

  write(*,*)'                                                               '
  write(*,*)'************** GPU 1D cuFFT Z2Z -> Z2Z ************************'
  write(*,*)'                                                               '

  fh = 0.0d0
  call t_cuFFT_gpu%gather_fft_1D_Z2Z_gpu(nx,fu,fh,CUFFT_FORWARD)
  fu_gpu = 0.0d0
  call t_cuFFT_gpu%gather_fft_1D_Z2Z_gpu(nx,fh,fu_gpu,CUFFT_INVERSE)

  if ( nx >= 5 ) then
     write(*,*) "the first 5 entries of:"
     write(*,*) "u         data:f(u)     Fourier: f(h)      Inverse Fourier"
     do i=1,nx-(nx-5)
        write(*,*)u(i), fu(i), fh(i), fu_gpu(i)/nx
     end do
  end if

  unit = 1
  err = convert_int2char_indexed_c(char_out,unit,1,1)
  tmr_name = 'fu_real'
  tmr_name = tmr_name(1:strlen(tmr_name))//char_out(1:strlen(char_out))
  tmr_name = tmr_name(1:strlen(tmr_name))//".asc"
  call file_open(tmr_name,unit,'unknown','asis','readwrite')

  !printing function to file for plotting
  if ( nx >= 5 ) then
     write(*,*) "the first 3x(i,j) = 9 entries of:"
     write(*,'(3(11x,a),4(15x,a))')"i","j","u","v","data:f(u,v)","Fourier: f(h,p)","Inverse Fourier"
     do i=1,nx
           write(unit,*)real(u(i)),real(fu(i))
     end do
  end if

  !freeing the ressources
  deallocate(u)
  deallocate(fu)
  deallocate(fh)
  deallocate(fu_gpu)

  call bye_cuFFT_gpu()
  !shutting down the environment
  call simple_cuda_shutdown()

#else
  write(*,*)"**************************WARNING******************************"
  write(*,*)"You need to compile with -DCUDA                                "
  write(*,*)"to acces the CUDA environment computation using GPU            "
  write(*,*)"switching back to the CPU version of corr function             "
  write(*,*)"***************************************************************"
#endif

  !shutting down the timers
  call stop_Alltimers_cpu()

  !end of greeting message
  call bye_deviceQuery_gpu()
  
end program testing_cuFFT_1D_Z2Z_gpu
