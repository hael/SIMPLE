! ============================================================================
! Name        : testing_cuFFT_2D_gpu.f90
! Author      : Frederic Bonnet
! Version     :
! Date        : 23rd of June 2015
! Description : tests the cuFFT library for 2D FFT on cuda
!             :
! ============================================================================
!
program testing_Functional3D_Z2Z

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

!  real(dp) :: xi,xf
  real(dp) :: a,b,c,n
  real(dp) :: delta
  !2D stuff
  integer :: nx,ny,nz
  !complex(dp)
  complex(dp) :: ui,uf,vi,vf,wi,wf
  complex(dp),allocatable :: u(:),v(:),w(:),fuvw(:,:,:)
  complex(dp),allocatable :: fuvw_gpu(:,:,:)
  complex(dp),allocatable :: fhpq(:,:,:)
  !counters
  integer  :: i,j,k

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
  a = 2
  b = 2
  c = 1
  n = 1
  nx = 130
  ny = 130
  nz = 130

!  call t_cuFFT_gpu%new_cuFFT_3D_ZZ_gpu(nx,ny,nz)

  write(*,*)'                                                               '
  write(*,*)'************** GPU 3D cuFFT Z2Z -> Z2Z ************************'
  write(*,*)'                                                               '

  !2D complex function
  ui = cmplx(-10.0, (4.0d0 * atan(1.0d0)/2.0d0)**(1/b) )
  uf = cmplx( 10.0, (4.0d0 * atan(1.0d0)/2.0d0)**(1/b) )
  vi = cmplx(-10.0, (4.0d0 * atan(1.0d0))**(1/b) )
  vf = cmplx( 10.0, (4.0d0 * atan(1.0d0))**(1/b) )
  wi = cmplx(-10.0, (4.0d0 * atan(1.0d0))**(1/b) )
  wf = cmplx( 10.0, (4.0d0 * atan(1.0d0))**(1/b) )

  allocate(u(nx))
  allocate(v(ny))
  allocate(w(nz))
  allocate(fuvw(nx,ny,nz))
  allocate(fhpq(nx,ny,nz))
  allocate(fuvw_gpu(nx,ny,nz))

  call atom_like_3D_complex(ui,uf,vi,vf,wi,wf,nx,ny,nz,a,b,c,u,v,w,fuvw,n) 
  !call getZRandomGaussianDistr_3D(nx,ny,nz,fuvw)

  fhpq = 0.0d0
  call t_cuFFT_gpu%gather_fft_3D_Z2Z_gpu(nx,ny,nz,fuvw,fhpq,CUFFT_FORWARD)

  fuvw_gpu = 0.0d0
  call t_cuFFT_gpu%gather_fft_3D_Z2Z_gpu(nx,ny,nz,fhpq,fuvw_gpu,CUFFT_INVERSE)
  
  if ( nx >= 5 ) then
     write(*,*) "the first 2x(i,j,k) = 8 entries of:"
     write(*,'(2x,a,2(10x,a))')"i","j","k","data:f(u,v,w)","Fourier: f(h,p,q)","Inverse Fourier"
     do i=1,nx-(nx-2)
        do j=1,ny-(ny-2)
           do k=1,nz-(nz-2)
              write(*,*)i,j,k,fuvw(i,j,k),fhpq(i,j,k),fuvw_gpu(i,j,k)/(nx*ny*nz)
           end do
        end do
     end do
  end if

  deallocate(u)
  deallocate(v)
  deallocate(w)
  deallocate(fuvw)
  deallocate(fhpq)
  deallocate(fuvw_gpu)

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
  
end program testing_Functional3D_Z2Z


