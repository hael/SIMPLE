program simple_test_order_corr
use simple_core_module_api
implicit none
type(oris)           :: os
type(ori)            :: o
integer, allocatable :: order(:)
integer              :: i
call seed_rnd
call os%new(11, is_ptcl=.false.)
call os%set_all2single('state',1.)
do i=1,11
    call os%set(i, 'corr', ran3())
end do
call os%set(7, 'state', 0.)
order = os%order()
do i=1,11
    call os%get_ori(order(i), o)
    call o%print_ori()
end do
call o%kill
if( size(order) == 11 )then
    print *,'PASSED'
else
    print *,'FAILED'
endif
end program simple_test_order_corr
