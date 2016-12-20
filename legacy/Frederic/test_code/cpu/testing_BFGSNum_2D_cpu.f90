! ============================================================================
! Name        : testing_BFGSNum_cpu
! Author      : Frederic Bonnet
! Version     :
! Date        : 3rd of July 2015
! Description : tests the BFGS minimizer
!             :
! ============================================================================
!
program testing_BFGSNum_cpu

  use simple_defs
  use simple_cuda_defs
  use greeting_version
  use simple_timing
  use simple_testfunction
  use simple_math
  !simple module to get the spec going and working
  use simple_params,only:params
  use simple_opt_spec
  use simple_bfgs_opt

  implicit none

  integer, parameter :: D = 2 !dimension of the function 2D:={x,y} \in R
  integer            :: err
  !local variables

  real(sp) :: a,b,c,n
  real(sp) :: delta
  !2D stuff
  integer  :: nx,ny
  !real(sp)
  real(sp) :: xi,xf
  real(sp) :: deriv,deriv_err
  !2D functions
  real(sp) :: deltax,deltay
  real(sp) :: x_D(D)
  real(sp) :: delta_D(D)
  real(sp) :: grad_D(D)
  real(sp) :: grad_symboli_D(D)
  real(sp),allocatable :: x(:),y(:),fxy(:,:)
  real(sp),allocatable :: grad_fxy(:,:)
  !complex(dp)
  complex(sp) :: ui,uf
  !bfgs stuff
  real(sp) :: lowest_cost
  !object stuff
  type(opt_spec) :: porispec
  type(bfgs_opt) :: bfgs_zer
  
  !counters
  integer  :: i,j
!  real(sp) :: func

  interface
     function func(point_d,D) result ( val )
       use simple_defs
       implicit none
       integer, intent(in) :: D
       real(sp), intent(in) :: point_d(D)
       real(sp) :: val
     end function func
     function gradfxy_symbolic(point_d,D) result (val)
       use simple_defs
       use simple_testfunction
       implicit none
       integer, intent(in) :: D
       real(sp), intent(inout) :: point_d(D) 
       real(sp) :: val(D)
     end function gradfxy_symbolic
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
  nx = 500
  ny = 500

  write(*,*)'                                                               '
  write(*,*)'************** CPU 2D Gradient R ******************************'
  write(*,*)'                                                               '

  ui = cmplx(-10.0d0,-10.0d0 )
  uf = cmplx( 10.0d0, 10.0d0 )

  allocate(x(nx))
  allocate(y(ny))
  allocate(fxy(nx,ny))
  allocate(grad_fxy(D*nx,ny))

  call atom_like_2D_single(ui,uf,nx,ny,a,b,c,x,y,fxy,n)

  deltax = ( real(uf) - real(ui) ) / (nx-1)
  deltay = ( imag(uf) - imag(ui) ) / (ny-1)
  delta_D(1) = deltax
  delta_D(2) = deltay
 
  do i=1,nx
     do j=1,ny
        x_D(1) = x(i)
        x_D(2) = y(j)
        grad_D = numgrad( func, x_D, D, delta_D )
        grad_fxy(i,j) = grad_D(1)
        grad_fxy(nx+i,j) = grad_D(2)
     end do
  end do

  if ( nx >= 5 ) then
     write(*,*) "the first 3x(i,j) = 9 entries of:"
     write(*,'(2x,a,2x,a,15x,a,15x,a,15x,a,15x,a,15x,a)') &
          "i","j","x","y",                                &
          "data:f(x,y)(R)",                               &
          "Gradient: f(x,y)(R)",                          &
          "Symbolic gradient: f(x,y)(R)"
     write(*,'(74x,a,6x,a,7x,a,7x,a)')                    &
          "pd_x f(x,y)(R)","pd_y f(x,y)(R)",              &
          "pd_x f(x,y)(R)","pd_y f(x,y)(R)"
     do i=1,nx-(nx-3)
        do j=1,ny-(ny-3)
           x_D(1) = x(i)
           x_D(2) = y(j)
            grad_symboli_D= gradfxy_symbolic(x_D,D)
           write(*,'(2x,i1,2x,i1,7(5x,f15.8))') &
                i,j,x(i),y(j),fxy(i,j),grad_fxy(i,j),grad_fxy(nx+i,j), &
                grad_symboli_D(1),grad_symboli_D(2)
        end do
     end do
  end if

!  do i=1,nx
!     do j=1,ny
!        write(2000,*) x(i),y(j),fxy(i,j),grad_fxy(i,j),grad_fxy(nx+i,j)
!     end do
!  end do

  write(*,*)'                                                               '
  write(*,*)'************** CPU 2D BFGS Hessian Symbolic R *****************'
  write(*,*)'                                                               '

  call porispec%specify('bfgs',nx)
  call porispec%set_costfun(func)
  call porispec%set_gcostfun(gradfxy_symbolic)

  write(*,*) "from the opt_spec object: ndim =",porispec%ndim

  call bfgs_zer%new(porispec)
  call bfgs_zer%minimize(porispec,lowest_cost)

  write(*,*) "lowest_cost: ",lowest_cost

  !realeasing ressources
  deallocate(x)
  deallocate(y)
  deallocate(fxy)
  deallocate(grad_fxy)

  !shutting down the timers
  call stop_Alltimers_cpu()

end program testing_BFGSNum_cpu

function gradfxy_symbolic(point_d,D) result (val)
  use simple_defs
  use simple_testfunction
  implicit none
  integer, intent(in) :: D
  real(sp), intent(inout) :: point_d(D) 
  real(sp) :: val(D)
  !local variables
  real(sp) :: a,b,c,n

  !this fucntion call is a check for a special case when
  a = 2
  b = 2
  c = 1
  n = 1

  val = atom_like_2D_single_ptwise_symbolicDiff(a,b,c,n,point_d,D)

end function gradfxy_symbolic

function func(point_d,D) result ( val )
  use simple_defs
  use simple_testfunction
  implicit none
  integer, intent(in) :: D
  real(sp), intent(in) :: point_d(D)
  real(sp) :: val
  !local variables
  real(sp) :: a,b,c,n
  real(sp) :: xi,xf

  a = 2
  b = 2
  c = 1
  n = 1

  val  = atom_like_2D_single_ptwise(a,b,c,point_d(1),point_d(2),n)

end function func


