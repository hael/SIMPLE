program simple_test_seg_regression
include 'simple_lib.f08'
use simple_optimizer,   only: optimizer
use simple_opt_factory, only: opt_factory
use simple_opt_spec,    only: opt_spec
implicit none
integer,          parameter :: N = 100, SEG_IND = 30
real,             parameter :: A1 = 0.2, B1 = 1., A2 = 0.8, B2 = -17
integer,          parameter :: NDIM = 1, NRESTARTS = 1
class(optimizer), pointer   :: opt_ptr=>null()      ! the generic optimizer object
type(opt_factory) :: ofac                           ! the optimization factory object
type(opt_spec)    :: spec                           ! the optimizer specification object
character(len=8)  :: str_opts                       ! string descriptors for the NOPTS optimizers
real              :: lims(2), lowest_cost, data_ori(N)
integer           :: i
! data generation
data_ori = 0.
do i = 1, SEG_IND
    data_ori(i) = A1 * real(i) + B1 + (ran3() - 0.5)
enddo
do i = SEG_IND+1, N
    data_ori(i) = A2 * real(i) + B2 + (ran3() - 0.5)
enddo
! setting up the optimization
str_opts = 'de'
lims(1)  =  1.
lims(2)  =  N
call spec%specify(str_opts, ndim=NDIM, limits=lims, nrestarts=NRESTARTS)
call spec%set_costfun(costfun)
call ofac%new(spec, opt_ptr)
spec%x = (1. + real(N))/2.
call opt_ptr%minimize(spec, opt_ptr, lowest_cost)

write(*, *) lowest_cost, spec%x

call opt_ptr%kill
deallocate(opt_ptr)

contains
    function costfun( self, vec, D ) result( cost )
        class(*), intent(inout) :: self
        integer,  intent(in)    :: D
        real,     intent(in)    :: vec(D)
        real                    :: cost
        cost = (vec(1) - 7.)**2
    end function costfun
end program simple_test_seg_regression
