program simple_test_complex_grad
include 'simple_lib.f08'
use simple_optimizer,         only: optimizer
use simple_opt_factory,       only: opt_factory
use simple_opt_spec,          only: opt_spec
implicit none
class(optimizer), pointer   :: opt_ptr=>null()      ! the generic optimizer object
integer,          parameter :: NDIM=2, NRESTARTS=1
real,             parameter :: C = 0.3 * PI, T = 0.2 * PI
type(opt_factory) :: ofac                           ! the optimization factory object
type(opt_spec)    :: spec                           ! the optimizer specification object
character(len=8)  :: str_opts                       ! string descriptors for the NOPTS optimizers
real              :: lowest_cost, truth_xy(2), xy(2), lims(2,2), carg
complex           :: A, B
A         = cmplx(ran3(), ran3())
truth_xy  = [ran3() * 4. - 2., ran3() * 4. - 2.]
B         = A * cmplx(cos(C * truth_xy(1)), sin(T * truth_xy(2)))
xy        = [ran3() - 0.5, ran3() - 0.5]
str_opts  = 'lbfgsb'
lims(1,1) = -5.
lims(1,2) =  5.
lims(2,1) = -5.
lims(2,2) =  5.
call spec%specify(str_opts, NDIM, limits=lims, nrestarts=NRESTARTS, factr  = 1.0d+5, pgtol = 1.0d-7)
call spec%set_costfun(costfct)                                      ! set pointer to costfun
call spec%set_gcostfun(gradfct)                                     ! set pointer to gradient of costfun
call ofac%new(spec, opt_ptr)                                        ! generate optimizer object with the factory
spec%x = xy
call opt_ptr%minimize(spec, opt_ptr, lowest_cost)                   ! minimize the test function
print *, 'CASE: quite stable complex-variable cost function'
print *, 'starting xy = ', xy
print *, 'truth    xy = ', truth_xy
print *, 'searched xy = ', spec%x
! current shift searching cost function
carg      = C * truth_xy(1) + T * truth_xy(2)
B         = A * cmplx(cos(carg), sin(carg))
call spec%set_costfun(costfct_2)                                      ! set pointer to costfun
call spec%set_gcostfun(gradfct_2)                                     ! set pointer to gradient of costfun
spec%x = xy
call opt_ptr%minimize(spec, opt_ptr, lowest_cost)                   ! minimize the test function
print *, '-------------------------------------------------'
print *, 'CASE: unstable complex-variable cost function'
print *, 'starting xy = ', xy
print *, 'truth    xy = ', truth_xy
print *, 'searched xy = ', spec%x

contains

    function costfct( fun_self, x, d ) result( r )
        class(*), intent(inout) :: fun_self
        integer,  intent(in)    :: d
        real,     intent(in)    :: x(d)
        real    :: r
        complex :: diff
        diff = A * cmplx(cos(C * x(1)), sin(T * x(2))) - B
        r    = real(diff * conjg(diff))
    end function

    subroutine gradfct( fun_self, x, grad, d )
        class(*), intent(inout) :: fun_self
        integer,  intent(in)    :: d
        real,     intent(inout) :: x(d)
        real,     intent(out)   :: grad(d)
        complex :: diff
        real    :: r
        diff    = A * cmplx(cos(C * x(1)), sin(T * x(2))) - B
        r       = real(diff * conjg(diff))
        grad(1) = - 2. * C * sin(C*x(1)) *  real(A * conjg(diff))
        grad(2) = - 2. * T * cos(T*x(2)) * aimag(A * conjg(diff))
    end subroutine

    function costfct_2( fun_self, x, d ) result( r )
        class(*), intent(inout) :: fun_self
        integer,  intent(in)    :: d
        real,     intent(in)    :: x(d)
        real    :: r, carg
        complex :: diff
        carg = C * x(1) + T * x(2)
        diff = A * cmplx(cos(carg), sin(carg)) - B
        r    = real(diff * conjg(diff))
    end function

    subroutine gradfct_2( fun_self, x, grad, d )
        class(*), intent(inout) :: fun_self
        integer,  intent(in)    :: d
        real,     intent(inout) :: x(d)
        real,     intent(out)   :: grad(d)
        complex :: diff
        real    :: r, carg, tmp
        carg    = C * x(1) + T * x(2)
        diff    = A * cmplx(cos(carg), sin(carg)) - B
        r       = real(diff * conjg(diff))
        tmp     = -2. * (real(A * conjg(diff)) * sin(carg) + aimag(A * conjg(diff)) * cos(carg))
        grad(1) = C * tmp
        grad(2) = T * tmp
    end subroutine

end program simple_test_complex_grad
