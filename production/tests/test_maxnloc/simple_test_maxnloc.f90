program simple_test_maxnloc
use simple_math
use simple_ran_tabu
implicit none
integer, parameter :: NNRS = 1000, NSEL = 10, NTST=100000
real    :: arr(NNRS), arr_copy(NNRS)
integer :: indxarr(NNRS), i, loc(NSEL)
type(ran_tabu) :: rt
rt = ran_tabu(NNRS)
do i = 1, NNRS
    arr(i) = real(i)
end do
print *, 'testing maxnnloc'
call rt%shuffle(arr)
loc = maxnloc(arr, NSEL)
arr_copy = arr
indxarr = (/(i,i=1,NNRS)/)
call hpsort(arr, indxarr)
call reverse(indxarr)
do i=1,NSEL
    print *, i, arr_copy(indxarr(i)), arr_copy(loc(i))
end do
print *, ''
print *, 'testing minnnloc'
call rt%shuffle(arr)
loc = minnloc(arr, NSEL)
arr_copy = arr
indxarr = (/(i,i=1,NNRS)/)
call hpsort(arr, indxarr)
do i=1,NSEL
    print *, i, arr_copy(indxarr(i)), arr_copy(loc(i))
end do
end program simple_test_maxnloc
