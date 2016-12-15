! ============================================================================
! Name        : testing_fft_1D_Z2Z_cpu.f90
! Author      : Frederic Bonnet
! Version     :
! Date        : 17th of June 2015
! Description : tests the functionality of the 1D FFT library and calls
!             :
! ============================================================================
!
program testing_fft_1D_Z2Z_cpu

  use simple_defs
  use simple_cuda_defs
  use greeting_version
  use simple_timing
  use fft123D_cpu
  use simple_testfunction
  use simple_systemQuery_cpu

  implicit none
#define devptr_t integer*8

  type(systemQuery_cpu)           :: sysQ
  type(systemDetails)             :: hstD
  type(fft_cpu)                   :: t_fft_cpu
  integer                         :: err
  !local variables
  integer                         :: nCPU_cores
  complex(dp) :: ui,uf
  integer  :: nx
  real(dp) :: a,b,c,n
  real(dp) :: delta
  complex(dp),allocatable :: u(:),fu(:)
  complex(dp),allocatable :: fh(:)
  complex(dp),allocatable :: fu_cpu(:)
  !counters
  integer  :: i

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

  !real function 
  !test function details
  nx = 100000
  a = 2
  b = 2
  c = 1
  n = 1

!  call t_fft_cpu%new_fft_1D_ZZ_cpu(nx)

  allocate(u(nx))
  allocate(fu(nx))
  allocate(fh(nx))
  allocate(fu_cpu(nx))

  !complex function
  ui = cmplx(-10.0, (4.0d0 * atan(1.0d0)/2.0d0)**(1/b) )
  uf = cmplx( 10.0, (4.0d0 * atan(1.0d0)/2.0d0)**(1/b) )
  fu = atom_like_complex(ui,uf,nx,a,b,c,u,n)
!  call t_fft_cpu%set_Zdata_1D(nx,u,fu)

  write(*,*)'                                                               '
  write(*,*)'************** CPU fftw 1D Z2Z -> Z2Z *************************'
  write(*,*)'                                                               '

  fh = 0.0d0
  call t_fft_cpu%gather_fft_1D_Z2Z_cpu(nx,fu,fh,FFTW_FORWARD)
  fu_cpu = 0.0d0
  call t_fft_cpu%gather_fft_1D_Z2Z_cpu(nx,fh,fu_cpu,FFTW_BACKWARD)

  if ( nx >= 5 ) then
     write(*,*) "the first 5 entries of:"
     write(*,*) "u         data:f(u)     Fourier: f(h)      Inverse Fourier"
     do i=1,nx-(nx-5)
        write(*,*)u(i), fu(i), fh(i), fu_cpu(i)/nx
     end do
  end if

  call bye_fft_cpu()
  !shutting down the timers
  call stop_Alltimers_cpu()

  !freeing the ressources
  deallocate(u)
  deallocate(fu)
  deallocate(fh)
  deallocate(fu_cpu)

end program testing_fft_1D_Z2Z_cpu




