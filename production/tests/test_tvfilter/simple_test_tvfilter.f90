program simple_test_tvfilter
include 'simple_lib.f08'
use simple_tvfilter, only: test_tvfilter, test_tvfilter_3d
implicit none
write (*,*) '2D:'
call test_tvfilter( [512,400,1], 1., 5. )
write (*,*) '2D: finished'
write (*,*) '3D:'
call test_tvfilter_3d( [100,100,100], 1., 5. )
write (*,*) '3D: finished'
end program simple_test_tvfilter
