! ============================================================================
! Name        : testing_fft_2D_S2C_cpu.f90
! Author      : Frederic Bonnet
! Version     :
! Date        : 3rd of July 2015
! Description : tests the functionality of the 2D FFT library and calls
!             :
! ============================================================================
!
program testing_fft_2D_S2C_cpu

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
  real(sp) :: a,b,c,n
  real(sp) :: delta
  !2D stuff
  integer  :: nx,ny
  !real(sp)
  real(sp),allocatable :: x(:),y(:)
  real(sp),allocatable :: fxy(:,:)
  real(sp),allocatable :: fxy_cpu(:,:)
  !complex(sp)
  complex(sp) :: ui,uf
  complex(sp),allocatable :: fhp(:,:)
  !counters
  integer  :: i,j

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
  ui = cmplx(-10.0,-10.0 )
  uf = cmplx( 10.0, 10.0 )
  a = 2
  b = 2
  c = 1
  n = 1
  nx = 4500
  ny = 4500

!  call t_fft_cpu%new_fft_2D_SC_cpu(nx,ny)

  write(*,*)'                                                               '
  write(*,*)'************** GPU 2D FFTW S2C -> C2S *************************'
  write(*,*)'                                                               '

  !2D complex function

  allocate(fxy(nx,ny))
  allocate(x(nx))
  allocate(y(ny))

  x = 0.0
  y = 0.0
  fxy = 0.0
  
  call atom_like_2D_single(ui,uf,nx,ny,a,b,c,x,y,fxy,n)

  allocate(fhp(nx,(ny/2+1)))
  allocate(fxy_cpu(nx,ny))

  fhp = 0.0
  call start_timer_cpu("S2C_gpu")
  call t_fft_cpu%gather_fft_2D_S2C_cpu(nx,ny,fxy,fhp)
  call stop_timer_cpu("S2C_gpu")

  fxy_cpu = 0.0
  call start_timer_cpu("C2S_gpu")
  call t_fft_cpu%gather_fft_2D_C2S_cpu(nx,ny,fhp,fxy_cpu)
  call stop_timer_cpu("C2S_gpu")
  
  if ( nx >= 5 ) then
     write(*,*) "the first 3x(i,j) = 9 entries of:"
     write(*,'(2x,a,2x,a,5(15x,a))')"i","j","x","y","data:f(x,y)(D)","Fourier: f(h,p)(Z)","Inverse Fourier(D)"
     do i=1,nx-(nx-3)
        do j=1,ny-(ny-3)
           write(*,*)i,j,x(i),y(j),fxy(i,j),fhp(i,j),fxy_cpu(i,j)/(nx*ny)
        end do
     end do
  end if

  !freeing the ressources
  deallocate(fhp)
  deallocate(fxy_cpu)

  deallocate(x)
  deallocate(y)
  deallocate(fxy)

  call bye_fft_cpu()
  !shutting down the timers
  call stop_Alltimers_cpu()

end program testing_fft_2D_S2C_cpu




