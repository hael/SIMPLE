program simple_test_ppca_kpca
include 'simple_lib.f08'
use simple_ppca_inmem, only: ppca_inmem
use simple_pca,        only: pca
implicit none
integer, parameter :: NP = 3, NS = 4, NC = 2, MAXPCAITS = 15
type(ppca_inmem)   :: prob_pca
type(pca)          :: pca_obj
integer :: j
real    :: data_ori(NP, NS), avg(NP), tmpvec(NP), data_pca(NP, NS), E_zn(NC, NS), eig_vecs(NP, NS), eig_vals(NS), tmp(NS, NS)
data_ori(1,:) = [ 1, 2, 3, 4]
data_ori(2,:) = [ 3, 1, 5, 8]
data_ori(3,:) = [-1, 0, 4, 10]
print *, data_ori(1,:)
print *, data_ori(2,:)
avg           = sum(data_ori, dim=2) / real(NS)
data_ori(1,:) = data_ori(1,:) - avg(1)
data_ori(2,:) = data_ori(2,:) - avg(2)
data_ori(3,:) = data_ori(3,:) - avg(3)
call prob_pca%new(NS, NP, NC)
call prob_pca%master(data_ori, MAXPCAITS)
!$omp parallel do private(j,tmpvec) default(shared) proc_bind(close) schedule(static)
do j = 1, NS
    call prob_pca%generate(j, avg, tmpvec)
    data_pca(:,j) = tmpvec
end do
!$omp end parallel do
print *, data_pca(1,:)
print *, data_pca(2,:)
do j = 1, NS
    E_zn(:,j) = prob_pca%get_feat(j)
enddo
do j = 1, NC
    print *, E_zn(j, :)
    print *, prob_pca%get_var(j)
enddo
! SVD test
eig_vecs = data_ori
call svdcmp(eig_vecs,eig_vals,tmp)
print *, eig_vals**2 / 3.
print *, eig_vecs(:,1)
! PCA test
call pca_obj%new(NS, NP, NC)
call pca_obj%master(data_ori)
end program simple_test_ppca_kpca
