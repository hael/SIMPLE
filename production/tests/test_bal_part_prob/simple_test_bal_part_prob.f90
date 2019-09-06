program simple_test_bal_part_prob
include 'simple_lib.f08'
use simple_map_reduce
implicit none
integer, parameter :: NPTCLS=20000, NPARTS=5
real,    parameter :: MAXNPEAKS=80.
real    :: npeaks(NPTCLS), npeaks_copy(NPTCLS), sums(NPARTS)
integer :: i, kassgn(NPTCLS)
call seed_rnd
do i=1,NPTCLS
    npeaks(i) = nint(ran3()*MAXNPEAKS)
end do
npeaks_copy = npeaks
call approx_balanced_partitioning_prob(npeaks, NPTCLS, NPARTS, kassgn)
! regenerate the sums
print *, 'presumably approximatively balanced sums'
do i=1,NPARTS
    sums(i) = sum(npeaks_copy, mask=kassgn .eq. i)
    print *, i, sums(i)
end do
end program simple_test_bal_part_prob
