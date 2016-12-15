! ============================================================================
! Name        : testing_cuFFT_2D_gpu.f90
! Author      : Frederic Bonnet
! Version     :
! Date        : 23rd of June 2015
! Description : tests the cuFFT library for 2D FFT on cuda
!             :
! ============================================================================
!
program testing_cuFFT_2D_S2C_gpu

  use simple_defs
  use simple_cuda_defs
  use simple_cuda
  use greeting_version
  use simple_timing
  use simple_deviceQuery_gpu
  use fft123D_gpu
  use simple_testfunction

  implicit none
#define devptr_t integer*8

  type(deviceQuery_gpu)           :: devQ
  type(deviceDetails)             :: devD
  type(cuFFT_gpu)                 :: t_cuFFT_gpu
  integer                         :: err
  !local variables
  integer                         :: ndev

  real(sp) :: a,b,c,n
  real(sp) :: delta
  !2D stuff
  integer :: nx,ny
  !real(dp)
  real(sp),allocatable :: x(:),y(:)
  real(sp),allocatable :: fxy(:,:)
  real(sp),allocatable :: fxy_gpu(:,:)
  !complex(dp)
  complex(sp) :: ui,uf
  complex(sp),allocatable :: fhp(:,:)
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
  ui = cmplx(-10.0,-10.0 )
  uf = cmplx( 10.0, 10.0 )
  a = 2
  b = 2
  c = 1
  n = 1
  nx = 8500
  ny = 8500

!  call t_cuFFT_gpu%new_cuFFT_2D_SC_gpu(nx,ny)

  write(*,*)'                                                               '
  write(*,*)'************** GPU 2D cuFFT S2C -> C2S ************************'
  write(*,*)'                                                               '

  allocate(fxy(nx,ny))
  allocate(x(nx))
  allocate(y(ny))

  x = 0.0
  y = 0.0
  fxy = 0.0

  call atom_like_2D_single(ui,uf,nx,ny,a,b,c,x,y,fxy,n)

  allocate(fhp(nx,(ny/2+1)))
  allocate(fxy_gpu(nx,ny))

  fhp = 0.0
  call start_timer_cpu("S2C_gpu")
  call t_cuFFT_gpu%gather_fft_2D_S2C_gpu(nx,ny,fxy,fhp) ! forward
  call stop_timer_cpu("S2C_gpu")

  fxy_gpu = 0.0
  call start_timer_cpu("C2S_gpu")
  call t_cuFFT_gpu%gather_fft_2D_C2S_gpu(nx,ny,fhp,fxy_gpu) ! backward
  call stop_timer_cpu("C2S_gpu")

  if ( nx >= 5 ) then
     write(*,*) "the first 3x(i,j) = 9 entries of:"
     write(*,'(2x,a,2x,a,5(15x,a))')"i","j","x","y","data:f(x,y)(D)","Fourier: f(h,p)(Z)","Inverse Fourier(D)"
     do i=1,nx-(nx-3)
        do j=1,ny-(ny-3)
           write(*,*)i,j,x(i),y(j),fxy(i,j),fhp(i,j),fxy_gpu(i,j)/(nx*ny)
        end do
     end do
  end if

  !freeing the ressources
  deallocate(fhp)
  deallocate(fxy_gpu)

  deallocate(x)
  deallocate(y)
  deallocate(fxy)

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
  
end program testing_cuFFT_2D_S2C_gpu
