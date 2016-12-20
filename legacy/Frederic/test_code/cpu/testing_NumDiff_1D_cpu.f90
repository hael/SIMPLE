! ============================================================================
! Name        : testing_NumDiff_1D_cpu.f90
! Author      : Frederic Bonnet
! Version     :
! Date        : 3rd of July 2015
! Description : testing the numerical differentiation method
!             :
! ============================================================================
!
program testing_NumDiff_1D_cpu

  use simple_defs
  use simple_cuda_defs
  use greeting_version
  use simple_timing
  use simple_testfunction
  use simple_math

  implicit none

  integer                         :: err
  !local variables

  real(sp) :: a,b,c,n
  real(sp) :: delta
  !2D stuff
  integer  :: nx,ny
  !real(sp)
  real(sp) :: xi,xf
  real(sp) :: deriv,deriv_err
  !1D functions
  real(sp),allocatable :: x(:),fx(:)
  real(sp),allocatable :: fx_int(:)
  real(sp),allocatable :: dfx(:)
  real(sp),allocatable :: dfx_err(:)

  !counters
  integer  :: i,j
!  real(sp) :: func

  interface
     function func(point) result ( val )
       use simple_defs
       implicit none
       real(sp), intent(in) :: point 
       real(sp) :: val
     end function func
     function func_dfx(point) result (val)
       use simple_defs
       implicit none
       real(sp), intent(in) :: point 
       real(sp) :: val
     end function func_dfx
  end interface

  !start of the execution commands
  !start of the greeting message
  call timestamp()
  call start_Alltimers_cpu()

  !test function details
  a = 2
  b = 2
  c = 1
  n = 1
  nx = 1500
  ny = 1500

  write(*,*)'                                                               '
  write(*,*)'************** CPU 1D differentiation R ***********************'
  write(*,*)'                                                               '

  !1D real function

  xi = -10.0
  xf = 10.0

  allocate(x(nx))
  allocate(fx(nx))
  allocate(dfx(nx))
  allocate(dfx_err(nx))
  allocate(fx_int(nx))

  fx = atom_like_real(xi,xf,nx,a,b,c,x,n)

  delta = (xf - xi) / (nx-1)
  do i = 1,nx
     call numderiv(func,x(i),delta,deriv_err,deriv)
     dfx(i) = deriv
     dfx_err(i) = err
  end do

  if ( nx >= 5 ) then
     write(*,*) "the first 3x(i,j) = 9 entries of:"
     write(*,'(2x,a,2(10x,a),6x,a,3x,a)')"i","x(i)","data:f(x)","diff_num: f(x)","diff_symbolic f(x)"
     do i=1,nx-(nx-3)
        write(*,'(2x,i1,4(2x,f15.8))')i,x(i),fx(i),dfx(i),func_dfx(x(i))
     end do
  end if

!  do i=1,nx
!     write(2000,*)x(i),fx(i)
!     write(3000,*)x(i),dfx(i)
!     write(4000,*)x(i),fx_int(i)
!  end do

  !freeing the ressources
  deallocate(x)
  deallocate(fx)
  deallocate(dfx)
  deallocate(dfx_err)
  deallocate(fx_int)

  !shutting down the timers
  call stop_Alltimers_cpu()

end program testing_NumDiff_1D_cpu

function func_dfx(point) result (val)
  use simple_defs
  use simple_testfunction
  implicit none
  real(sp), intent(in) :: point 
  real(sp) :: val

  !this fucntion call is a check for a special case when
  !a = 2  b = 2  c = 1  n = 1

  val  = atom_like_real_ptwise_symbolicDiff(point)

end function func_dfx

function func(point) result ( val )
  use simple_defs
  use simple_testfunction
  implicit none
  real(sp), intent(in) :: point 
  real(sp) :: val
  !local variables
  real(sp) :: a,b,c,n
  real(sp) :: xi,xf

  a = 2
  b = 2
  c = 1
  n = 1

  val  = atom_like_real_ptwise(a,b,c,point,n)

end function func


