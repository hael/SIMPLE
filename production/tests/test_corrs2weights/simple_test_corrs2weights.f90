program simple_test_corrs2weights
include 'simple_lib.f08'
implicit none
real    :: corrs(12), weights(12)
integer :: i
corrs(1) = -1.
corrs(2) = 0.0
corrs(3) = 0.005
corrs(4) = 0.1
corrs(5) = 0.2
corrs(6) = 0.3
corrs(7) = 0.4
corrs(8) = 0.5
corrs(9) = 0.51
corrs(10) = 0.52
corrs(11) = 0.53
corrs(12) = 0.6
weights = corrs2weights(corrs, CORRW_CRIT)
do i=1,size(corrs)
    print *, 'corr/weight: ', corrs(i), weights(i)
end do
end program simple_test_corrs2weights
