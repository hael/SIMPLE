program simple_test_rank_weights
use gnufor2, only: plot
use simple_stat
implicit none
real    :: ranks(200), weights(200)
integer :: i
do i=1,100
    ranks(i) = real(i)
end do
call rank_sum_weights(200, weights)
! call plot(ranks, weights)
call rank_inverse_weights(200, weights)
! call plot(ranks, weights)
call rank_centroid_weights(200, weights)
! call plot(ranks, weights)
call rank_exponent_weights(200, 10.0, weights)
call plot(ranks, weights)
end program
