! ============================================================================
! Name        : testing_fft_3D_Z2Z_cpu.f90
! Author      : Frederic Bonnet
! Version     :
! Date        : 3rd of July 2015
! Description : tests the functionality of the 2D FFT library and calls
!             :
! ============================================================================
!
program testing_fft_3D_Z2Z_cpu

  use simple_defs
  use simple_cuda_defs
  use greeting_version
  use simple_timing
  use fft123D_cpu
  use simple_testfunction
  use simple_systemQuery_cpu

  implicit none

  type(systemQuery_cpu)           :: sysQ
  type(systemDetails)             :: hstD
  type(fft_cpu)                   :: t_fft_cpu
  integer                         :: err
  !local variables
  integer                         :: nCPU_cores
  real(dp) :: a,b,c,n
  real(dp) :: delta
  !2D stuff
  integer  :: nx,ny,nz
  !complex(dp)
  complex(dp) :: ui,uf,vi,vf,wi,wf
  complex(dp),allocatable :: u(:),v(:),w(:),fuvw(:,:,:)
  complex(dp),allocatable :: fuvw_cpu(:,:,:)
  complex(dp),allocatable :: fhpq(:,:,:)
  !counters
  integer  :: i,j,k

  !start of the execution commands
  !start of the greeting message
  call timestamp()
  call start_Alltimers_cpu()

  call hello_fft_cpu(err)

  call sysQ%new_systemQuery_cpu(hstD)

  write(*,*)'                                                                  '
  write(*,*)'******************************************************************'
  write(*,*)'   Sanity checks on the object(sysQ) and the data structure(hstD) '
  write(*,*)'******************************************************************'
  nCPU_cores = sysQ%get_ncpu_cores()
  write(*,*)'Number of cores on the system                    : ',nCPU_cores
  !writting the number of cores avaible on the system
  write(*,*)'Number of cores returned from data structure hstD: ',hstD%nCPUcores
  if ( nCPU_cores /= hstD%nCPUcores) call sysQ%get_warning_dataStruct()
  write(*,*)'******************************************************************'

  !test function details
  a = 2
  b = 2
  c = 1
  n = 1
  nx = 230
  ny = 230
  nz = 230

!  call t_fft_cpu%new_fft_3D_ZZ_cpu(nx,ny,nz)

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
  allocate(fuvw_cpu(nx,ny,nz))

  !call atom_like_3D_complex(ui,uf,vi,vf,wi,wf,nx,ny,nz,a,b,c,u,v,w,fuvw,n) 
  call getZRandomGaussianDistr_3D(nx,ny,nz,fuvw)

  fhpq = 0.0d0
  call t_fft_cpu%gather_fft_3D_Z2Z_cpu(nx,ny,nz,fuvw,fhpq,FFTW_FORWARD)

  fuvw_cpu = 0.0d0
  call t_fft_cpu%gather_fft_3D_Z2Z_cpu(nx,ny,nz,fhpq,fuvw_cpu,FFTW_BACKWARD)
  
  if ( nx >= 5 ) then
     write(*,*) "the first 2x(i,j,k) = 8 entries of:"
     write(*,'(2x,a,2(10x,a))')"i","j","k","data:f(u,v,w)","Fourier: f(h,p,q)","Inverse Fourier"
     do i=1,nx-(nx-2)
        do j=1,ny-(ny-2)
           do k=1,nz-(nz-2)
              write(*,*)i,j,k,fuvw(i,j,k),fhpq(i,j,k),fuvw_cpu(i,j,k)/(nx*ny*nz)
           end do
        end do
     end do
  end if

  deallocate(u)
  deallocate(v)
  deallocate(w)
  deallocate(fuvw)
  deallocate(fhpq)
  deallocate(fuvw_cpu)

  call bye_fft_cpu()
  !shutting down the timers
  call stop_Alltimers_cpu()

end program testing_fft_3D_Z2Z_cpu




