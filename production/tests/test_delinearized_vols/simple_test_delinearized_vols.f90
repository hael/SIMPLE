program simple_test_delinearized_vols
include 'simple_lib.f08'
implicit none
integer :: ivol, j, errflg
real    :: truths(4,3), vols(4,3), probs(3,3), probs_inv(3,3), truths_inv(4,3)
truths(:,1) = [1., 2., 3., 4.]
truths(:,2) = [.5, 2., 3., 6.]
truths(:,3) = [.8, 2., 3., 5.]
probs(1,:)  = [0.8, 0.1,  0.1]
probs(2,:)  = [0.7, 0.12, 0.18]
probs(3,:)  = [0.9, 0.06, 0.04]
vols        = 0.
do ivol = 1, 3
    do j = 1, 3
        vols(:,ivol) = vols(:,ivol) + probs(ivol,j) * truths(:,j)
    enddo
enddo
call matinv(probs, probs_inv, 3, errflg)
! recovering
truths_inv = 0.
do ivol = 1, 3
    do j = 1, 3
        truths_inv(:,ivol) = truths_inv(:,ivol) + probs_inv(ivol,j) * vols(:,j)
    enddo
    print *, truths_inv(:,ivol)
enddo
end program simple_test_delinearized_vols
