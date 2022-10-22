program simple_test_rank_weights
use gnufor2, only: plot
include 'simple_lib.f08'
implicit none
real    :: ranks(200), weights(200)
integer :: i
do i=1,200
    ranks(i) = real(i)
end do
call rank_sum_weights(200, weights)
! call plot(ranks, weights)
call rank_inverse_weights(200, weights)
! call plot(ranks, weights)
call rank_centroid_weights(200, weights)
! call plot(ranks, weights)
call rank_exponent_weights(200, 10.0, weights)
! call plot(ranks, weights)
do i=1,100
    weights(i) = 101.0 - real(i)
end do
do i=101,200
    weights(i) = real(i) - 99.5
end do
call plot(ranks, weights)
call conv2rank_weights(200, weights, RANK_SUM_CRIT)
call plot(ranks, weights)
call conv2rank_weights(200, weights, RANK_CEN_CRIT)
! call plot(ranks, weights)
call conv2rank_weights(200, weights, RANK_EXP_CRIT, p=2.0)
! call plot(ranks, weights)
call conv2rank_weights(200, weights, RANK_INV_CRIT)
! call plot(ranks, weights)
end program
