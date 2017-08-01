! ============================================================================
! Name        : testing_cuFFT_2D_gpu.f90
! Author      : Frederic Bonnet
! Version     :
! Date        : 23rd of June 2015
! Description : tests the cuFFT library for 2D FFT on cuda
!             :
! ============================================================================
!
program testing_cuFFT_3D_Z2D_gpu

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

  !2D stuff
  integer :: nx,ny,nz
  !real(dp)
  real(dp),allocatable :: fhpq(:,:,:)
  !complex(dp)
  complex(dp),allocatable :: fuvw(:,:,:)
  complex(dp),allocatable :: fuvw_gpu(:,:,:)
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
  nx = 330
  ny = 330
  nz = 330

!  call t_cuFFT_gpu%new_cuFFT_3D_DZ_gpu(nx,ny,nz)

  write(*,*)'                                                               '
  write(*,*)'************** GPU 3D cuFFT Z2D -> D2Z ************************'
  write(*,*)'                                                               '

  allocate(fuvw(nx,ny,((nz/2)+1)))
  allocate(fhpq(nx,ny,nz))
  allocate(fuvw_gpu(nx,ny,((nz/2)+1)))

  !call atom_like_3D_complex(ui,uf,vi,vf,wi,wf,nx,ny,nz,a,b,c,u,v,w,fuvw,n) 
  call getZRandomGaussianDistr_3D(nx,ny,nz/2+1,fuvw)

  fhpq = 0.0d0
  call t_cuFFT_gpu%gather_fft_3D_Z2D_gpu(nx,ny,nz,fuvw,fhpq)

  fuvw_gpu = 0.0d0
  call t_cuFFT_gpu%gather_fft_3D_D2Z_gpu(nx,ny,nz,fhpq,fuvw_gpu)
  
  if ( nx >= 5 ) then
     write(*,*) "the first 8x(i,j,k) = 8 entries of:"
     write(*,'(11x,a,11x,a,10x,a,2x,12x,a,2(20x,a))')"i","j","k","data:f(u,v,w)(Z)","Fourier: f(h,p,q)(D)","Inverse Fourier(Z)"
     do i=1,nx-(nx-2)
        do j=1,ny-(ny-2)
           do k=1,nz-(nz-2)
              write(*,*)i,j,k,fuvw(i,j,k),fhpq(i,j,k),fuvw_gpu(i,j,k)/(nx*ny*nz)
           end do
        end do
     end do
  end if

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
  
end program testing_cuFFT_3D_Z2D_gpu


