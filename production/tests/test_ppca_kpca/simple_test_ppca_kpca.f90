program simple_test_ppca_kpca
use simple_ppca_inmem, only: ppca_inmem
implicit none
integer, parameter :: NP = 2, NS = 4, NC = 1, MAXPCAITS = 25
type(ppca_inmem)   :: prob_pca
integer :: j
real    :: data_ori(NP, NS), avg(NP), tmpvec(NP), data_pca(NP, NS), W(NP, NC), E_zn(NC, NS)
data_ori(1,:) = [1, 2, 3, 4]
data_ori(2,:) = [3, 1, 5, 8]
print *, data_ori(1,:)
print *, data_ori(2,:)
avg           = sum(data_ori, dim=2) / real(NS)
data_ori(1,:) = data_ori(1,:) - avg(1)
data_ori(2,:) = data_ori(2,:) - avg(2)
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
call prob_pca%get_W(W)
print *, W
call prob_pca%get_E_zn(E_zn)
print *, E_zn
end program simple_test_ppca_kpca
