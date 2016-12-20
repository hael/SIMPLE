! ============================================================================
! Name        : testing_cuFFT_2D_Z2D_gpu.f90
! Author      : Frederic Bonnet
! Version     :
! Date        : 23rd of June 2015
! Description : tests the cuFFT library for 2D FFT on cuda
!             :
! ============================================================================
!
program testing_cuFFT_2D_Z2D_gpu

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

  !real(dp) :: xi,xf
  real(dp) :: a,b,c,n
  real(dp) :: delta
  !2D stuff
  integer :: nx,ny
  !real(dp)
  real(dp),allocatable :: fhp(:,:)
  !complex(dp)
  complex(dp) :: ui,uf,vi,vf
  complex(dp),allocatable :: u(:),v(:),fuv(:,:)
  complex(dp),allocatable :: fuv_gpu(:,:)
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
  ui = cmplx(-10.0d0,-10.0d0 )
  uf = cmplx( 10.0d0, 10.0d0 )
  a = 2
  b = 2
  c = 1
  n = 1
  nx = 8500
  ny = 8500

!  call t_cuFFT_gpu%new_cuFFT_2D_DZ_gpu(nx,ny)

  write(*,*)'                                                               '
  write(*,*)'************** GPU 2D cuFFT Z2D -> D2Z ************************'
  write(*,*)'                                                               '

  !2D complex function
  ui = cmplx(-10.0, (4.0d0 * atan(1.0d0)/2.0d0)**(1/b) )
  uf = cmplx( 10.0, (4.0d0 * atan(1.0d0)/2.0d0)**(1/b) )
  vi = cmplx(-10.0, (4.0d0 * atan(1.0d0))**(1/b) )
  vf = cmplx( 10.0, (4.0d0 * atan(1.0d0))**(1/b) )

  allocate(u(nx))
  allocate(v(ny))
  allocate(fuv(nx,ny))
  allocate(fhp(nx,(ny/2+1)))
  allocate(fuv_gpu(nx,ny))
  
  call atom_like_2D_complex(ui,uf,vi,vf,nx,ny,a,b,c,u,v,fuv,n) 

  if ( nx >= 5 ) then
     write(*,*) "the first 3x(i,j) = 9 entries of:"
     write(*,'(2x,a,2x,a,3(15x,a))')"i","j","u","v","data:f(u,v)(Z)"
     do i=1,nx-(nx-3)
        do j=1,ny-(ny-3)
           write(*,*)i,j,u(i),v(j),fuv(i,j)
        end do
     end do
  end if

  fhp = 0.0d0
  call t_cuFFT_gpu%gather_fft_2D_Z2D_gpu(nx,ny,fuv,fhp)

  write(*,*)" !! Caution f(u,v)(Z) has been modified !! "

  fuv_gpu = 0.0d0
  call t_cuFFT_gpu%gather_fft_2D_D2Z_gpu(nx,ny,fhp,fuv_gpu)
  
  if ( nx >= 5 ) then
     write(*,*) "the first 3x(i,j) = 9 entries of:"
     write(*,'(2x,a,2x,a,5(15x,a))')"i","j","u","v","data:f(u,v)(Z)","Fourier: f(h,p)(D)","Inverse Fourier(Z)"
     do i=1,nx-(nx-3)
        do j=1,ny-(ny-3)
           write(*,*)i,j,u(i),v(j),fuv(i,j),fhp(i,j),fuv_gpu(i,j)/(nx*ny)
        end do
     end do
  end if

  deallocate(u)
  deallocate(v)
!  deallocate(fuv)
  deallocate(fhp)
  deallocate(fuv_gpu)

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
  
end program testing_cuFFT_2D_Z2D_gpu
