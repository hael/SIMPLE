program simple_test_tseries_neigh
include 'simple_lib.f08'
implicit none

type(oris) :: os
integer :: i, inds(20)
integer, allocatable :: ptcls2neigh(:,:)

call os%new(20, is_ptcl=.false.)
do i=1,20,5
    inds(i:) = i
end do
do i=1,20
    call os%set(i, 'class', real(inds(i)))
    print *, inds(i)
end do

call os%get_tseries_neighs(2, ptcls2neigh)

do i=1,20
    print *, i, inds(i), ptcls2neigh(i,:)
end do



end program simple_test_tseries_neigh
