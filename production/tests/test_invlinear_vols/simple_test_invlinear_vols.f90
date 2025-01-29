program simple_test_invlinear_vols
include 'simple_lib.f08'
use simple_opt_filter, only: uni_inv_linear
implicit none
integer, parameter :: NSTATES = 3, NP = 10, LOW_IND = 1, HIGH_IND = 4
integer :: ivol, j
real    :: truths(NP,NSTATES), vols(NP,NSTATES), probs(NSTATES,NSTATES), truths_inv(NP,NSTATES)
truths(:,1) = [5.,   4.,  3., 4., 5., 6., 7., 8.,  9., 10.]
truths(:,2) = [1.1,  1.9, 3., 4., 3., 5., 7., 8., 11., 14.]
truths(:,3) = [1.05, 2.1, 3., 4., 5., 6., 7., 8.,  9., 10.]
probs(1,:)  = [0.8,  0.1,  0.1]
probs(2,:)  = [0.12, 0.7,  0.18]
probs(3,:)  = [0.04, 0.06, 0.9]
vols        = 0.
do ivol = 1, NSTATES
    do j = 1, NSTATES
        vols(:,ivol) = vols(:,ivol) + probs(ivol,j) * truths(:,j)
    enddo
    print *, vols(:,ivol)
enddo
call uni_inv_linear(NSTATES, HIGH_IND - LOW_IND + 1, vols(LOW_IND:HIGH_IND,:), truths_inv(LOW_IND:HIGH_IND,:), verbose=.true.)
end program simple_test_invlinear_vols
