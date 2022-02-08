program simple_test_stringmatch
include 'simple_lib.f08'
implicit none


character(len=*), parameter :: projname = 'state_2'
character(len=*), parameter :: stkname  = '2_state.mrc'


print *, str_has_substr(projname, stkname)



end program
