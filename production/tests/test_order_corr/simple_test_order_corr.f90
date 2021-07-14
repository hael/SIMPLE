program simple_test_order_corr
use simple_oris, only: oris
use simple_ori,  only: ori
use simple_rnd,  only: ran3, seed_rnd
implicit none
type(oris) :: os
type(ori)  :: o
integer, allocatable :: order(:)
integer :: i
call seed_rnd
call os%new(10, is_ptcl=.false.)
do i=1,10
    call os%set(i, 'corr', ran3())
end do
order = os%order_corr()
do i=1,10
    call os%get_ori(order(i), o)
    call o%print_ori()
end do
call o%kill
end program simple_test_order_corr
