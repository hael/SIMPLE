program simple_test_stringmatch
include 'simple_lib.f08'
implicit none
character(len=*), parameter :: projname = 'state_2'
character(len=*), parameter :: stkname  = '2_state.mrc'
character(len=*), parameter :: inds_str = ' 1,3,  5,7,15  '
integer,        allocatable :: inds(:)
integer :: i 
inds = list_of_ints2arr(inds_str)
print *, 'found # integer numbers: ', size(inds)
do i = 1, size(inds)
    print *, inds(i)
end do
end program
