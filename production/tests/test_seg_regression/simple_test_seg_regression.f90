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
integer           :: i, inds(N)
! data generation
data_ori = 0.
inds     = (/(i,i=1,N)/)
do i = 1, SEG_IND
    data_ori(i) = A1 * real(i) + B1 + (ran3() - 0.5)
enddo
do i = SEG_IND+1, N
    data_ori(i) = A2 * real(i) + B2 + (ran3() - 0.5)
enddo
! setting up the optimization
str_opts = 'de'
lims(1)  = 1
lims(2)  = real(N)
call spec%specify(str_opts, ndim=NDIM, limits=lims, limits_init=lims, nrestarts=NRESTARTS)
call spec%set_costfun(costfun)
call ofac%new(spec, opt_ptr)
! spec%x = (1. + real(N))/2.
spec%x = 30.
call opt_ptr%minimize(spec, opt_ptr, lowest_cost)

print *, lowest_cost, spec%x
call opt_ptr%kill
deallocate(opt_ptr)

contains
    function costfun( self, vec, D ) result( cost )
        class(*), intent(inout)  :: self
        integer,  intent(in)     :: D
        real,     intent(in)     :: vec(D)
        real(dp) :: mean_x, mean_y, a, b, dcost, denom
        real     :: cost
        integer  :: mid
        mid    = floor(vec(1))
        ! linear fitting (y = a + bx ) from 1 to mid
        mean_x = sum(real(inds(    1:mid), dp)) / real(mid, dp)
        mean_y = sum(     data_ori(1:mid)     ) / real(mid, dp)
        denom  = sum((real(inds(1:mid), dp) - mean_x)**2)
        if( denom > TINY )then
            b  = sum( (real(inds(1:mid), dp) - mean_x) * (data_ori(1:mid) - mean_y) ) / denom
        else
            b  = 0.
        endif
        a      = mean_y - b * mean_x
        dcost  = sum( ((a + b * real(inds(1:mid), dp)) -  data_ori(1:mid))**2 ) / real(mid, dp)
        ! linear fitting from mid+1 to N
        if( mid < N )then
            mean_x = sum(real(inds(    mid+1:N), dp)) / real(N-mid, dp)
            mean_y = sum(     data_ori(mid+1:N)     ) / real(N-mid, dp)
            denom  = sum((real(inds(mid+1:N), dp) - mean_x)**2)
            if( denom > TINY )then
                b  = sum( (real(inds(mid+1:N), dp) - mean_x) * (data_ori(mid+1:N) - mean_y) ) / denom
            else
                b  = 0.
            endif
            a      = mean_y - b * mean_x
            dcost  = dcost + sum( ((a + b * real(inds(mid+1:N), dp)) -  data_ori(mid+1:N))**2 ) / real(N-mid, dp)
        endif
        cost = real(dcost)
    end function costfun
end program simple_test_seg_regression
