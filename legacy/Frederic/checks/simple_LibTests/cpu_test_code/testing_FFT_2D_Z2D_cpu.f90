! ============================================================================
! Name        : testing_FFT_2D_Z2D_cpu.f90
! Author      : Frederic Bonnet
! Version     :
! Date        : 3rd of July 2015
! Description : tests the functionality of the 2D FFT library and calls
!             :
! ============================================================================
!
program testing_fft_2D_Z2D_cpu

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
  integer  :: nx,ny
  !real(dp)
  real(dp),allocatable :: fhp(:,:)
  !complex(dp)
  complex(dp) :: ui,uf,vi,vf
  complex(dp),allocatable :: u(:),v(:),fuv(:,:)
  complex(dp),allocatable :: fuv_cpu(:,:)
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
  ui = cmplx(-10.0d0,-10.0d0 )
  uf = cmplx( 10.0d0, 10.0d0 )
  a = 2
  b = 2
  c = 1
  n = 1
  nx = 4500
  ny = 4500

!  call t_fft_cpu%new_fft_2D_DZ_cpu(nx,ny)

  write(*,*)'                                                               '
  write(*,*)'************** CPU 2D FFTW Z2D -> D2Z *************************'
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
  allocate(fuv_cpu(nx,ny))

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
  call t_fft_cpu%gather_fft_2D_Z2D_cpu(nx,ny,fuv,fhp)

  write(*,*)" !! Caution f(u,v)(Z) has been modified !! "

  fuv_cpu = 0.0d0
  call t_fft_cpu%gather_fft_2D_D2Z_cpu(nx,ny,fhp,fuv_cpu)
  
  if ( nx >= 5 ) then
     write(*,*) "the first 3x(i,j) = 9 entries of:"
     write(*,'(2x,a,2x,a,5(15x,a))')"i","j","u","v","data:f(u,v)(Z)","Fourier: f(h,p)(D)","Inverse Fourier(Z)"
     do i=1,nx-(nx-3)
        do j=1,ny-(ny-3)
           write(*,*)i,j,u(i),v(j),fuv(i,j),fhp(i,j),fuv_cpu(i,j)/(nx*ny)
        end do
     end do
  end if

  !freeing the ressources
  deallocate(u)
  deallocate(v)
!  deallocate(fuv)
  deallocate(fhp)
!  deallocate(fuv_cpu)

  call bye_fft_cpu()
  !shutting down the timers
  call stop_Alltimers_cpu()

end program testing_fft_2D_Z2D_cpu




