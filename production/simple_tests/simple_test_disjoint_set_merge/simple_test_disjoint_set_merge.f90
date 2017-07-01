program simple_test_disjoint_set_merge
use simple_combinatorics
implicit none

integer :: nnmat(5,3), i
integer, allocatable :: disjoint(:)

nnmat(1,:) = [1,2,3]
nnmat(2,:) = [2,3,4]
nnmat(3,:) = [6,7,8]
nnmat(4,:) = [3,4,5]
nnmat(5,:) = [9,10,11] ! correct answer: 1,2,3,4,5

disjoint = merge_into_disjoint_set( 5, 3, nnmat, [1,2,4] )

do i=1,size(disjoint)
	print *, disjoint(i)
end do


end program simple_test_disjoint_set_merge
