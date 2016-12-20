! ============================================================================
! Name        : testing_FFT_3D_D2Z_cpu.f90
! Author      : Frederic Bonnet
! Version     :
! Date        : 3rd of July 2015
! Description : tests the functionality of the 3D FFT library and calls
!             :
! ============================================================================
!
program testing_fft_3D_D2Z_cpu

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
  !3D stuff
  integer  :: nx,ny,nz
  !real(sp)
  real(sp),allocatable :: fxyz(:,:,:)
  real(sp),allocatable :: fxyz_cpu(:,:,:)
  !complex(sp)
  complex(sp),allocatable :: fhpq(:,:,:)
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
  nx = 230
  ny = 230
  nz = 230

  !now moving to the fast 3D fourier transform on CPU
!  call t_fft_cpu%new_fft_3D_SC_cpu(nx,ny,nz)

  write(*,*)'                                                               '
  write(*,*)'************** CPU 3D FFTW S2C -> C2S *************************'
  write(*,*)'                                                               '

  !3D complex random matrix

  allocate(fxyz(nx,ny,nz))
  allocate(fhpq(nx,ny,((nz/2)+1)))
  allocate(fxyz_cpu(nx,ny,nz))

  !call atom_like_3D_complex(ui,uf,vi,vf,wi,wf,nx,ny,nz,a,b,c,u,v,w,fuvw,n) 
  call getSRandomGaussianDistr_3D(nx,ny,nz,fxyz)

  fhpq = 0.0
  call start_timer_cpu("S2C_cpu")
  call t_fft_cpu%gather_fft_3D_S2C_cpu(nx,ny,nz,fxyz,fhpq)
  call stop_timer_cpu("S2C_cpu")

  fxyz_cpu = 0.0
  call start_timer_cpu("C2S_cpu")
  call t_fft_cpu%gather_fft_3D_C2S_cpu(nx,ny,nz,fhpq,fxyz_cpu)
  call stop_timer_cpu("C2S_cpu")
  
  if ( nx >= 5 ) then
     write(*,*) "the first 2x(i,j,k) = 8 entries of:"
     write(*,'(2x,a,2(10x,a))')"i","j","k","data:f(x,y,z)","Fourier: f(h,p,q)","Inverse Fourier"
     do i=1,nx-(nx-2)
        do j=1,ny-(ny-2)
           do k=1,nz-(nz-2)
              write(*,*)i,j,k,fxyz(i,j,k),fhpq(i,j,k),fxyz_cpu(i,j,k)/(nx*ny*nz)
           end do
        end do
     end do
  end if

  deallocate(fxyz)
  deallocate(fhpq)
  deallocate(fxyz_cpu)

  call bye_fft_cpu()
  !shutting down the timers
  call stop_Alltimers_cpu()

end program testing_fft_3D_D2Z_cpu




