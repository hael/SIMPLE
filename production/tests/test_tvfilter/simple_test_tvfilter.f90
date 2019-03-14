program simple_test_tvfilter
include 'simple_lib.f08'
use simple_tvfilter, only: test_tvfilter
implicit none
call test_tvfilter( [512,400,1], 1., 5. )
end program simple_test_tvfilter
