program simple_test_invlinear_vols
include 'simple_lib.f08'
use simple_opt_filter, only: uni_delinear
use simple_pca_svd,    only: pca_svd
implicit none
integer, parameter :: NSTATES = 3, NP = 10, LOW_IND = 1, HIGH_IND = 4, NC = 2, MAXPCAITS = 10, NINDS = HIGH_IND - LOW_IND + 1
type(pca_svd)      :: pca_obj
integer :: ivol, j
real    :: truths(NP,NSTATES), vols(NP,NSTATES), probs(NSTATES,NSTATES), truths_inv(NP,NSTATES), avg(NSTATES),&
        &data_cen(NSTATES, NINDS), tmpvec(NSTATES), data_pca(NSTATES, NINDS)
truths(:,1) = [5.,   4.,  3., 4., 5., 6., 7., 8.,  9., 10.]
truths(:,2) = [1.1,  1.9, 3., 4., 3., 5., 7., 8., 11., 14.]
truths(:,3) = [1.05, 2.1, 3., 4., 5., 6., 7., 8.,  9., 10.]
probs(1,:)  = [0.8,  0.1,  0.1]
probs(2,:)  = [0.12, 0.7,  0.18]
probs(3,:)  = [0.04, 0.06, 0.9]
vols        = 0.
print *, 'current data'
do ivol = 1, NSTATES
    do j = 1, NSTATES
        vols(:,ivol) = vols(:,ivol) + probs(ivol,j) * truths(:,j)
    enddo
    print *, vols(:,ivol)
enddo
print *, '-----'
call uni_delinear(NSTATES, HIGH_IND - LOW_IND + 1, vols(LOW_IND:HIGH_IND,:), truths_inv(LOW_IND:HIGH_IND,:), verbose=.true.)
print *, '-----'
print *, 'reconstructed data'
do ivol = 1, NSTATES
    print *, truths_inv(:,ivol)
enddo
! PCA
data_cen = transpose(vols(LOW_IND:HIGH_IND,:))
avg      = sum(data_cen, dim=2) / real(NSTATES)
do j = 1, NSTATES
    data_cen(j,:) = data_cen(j,:) - avg(j)
enddo
call pca_obj%new(NINDS, NSTATES, NC)
call pca_obj%master(data_cen, MAXPCAITS)
do j = 1, NINDS
    call pca_obj%generate(j, avg, tmpvec)
    data_pca(:,j) = tmpvec
end do
print *, 'Pre-imaged data using PCA:'
do j = 1, NSTATES
    print *, data_pca(j,:)
enddo
print *, sum(data_pca, dim=1)/real(NSTATES)
end program simple_test_invlinear_vols
