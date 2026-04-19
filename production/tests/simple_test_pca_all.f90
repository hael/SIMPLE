program simple_test_pca_all
use simple_core_module_api
!$ use omp_lib
use simple_ppca,       only: ppca
use simple_mppca,      only: mppca
use simple_pca_svd,    only: pca_svd
use simple_kpca_svd,   only: kpca_svd
use simple_cmdline,    only: cmdline
use simple_parameters, only: parameters
implicit none
integer, parameter :: NP = 3, NS = 4, NC = 3, MAXPCAITS = 15
type(ppca)         :: prob_pca
type(mppca)        :: mix_prob_pca
type(pca_svd)      :: pca_obj
type(kpca_svd)     :: kpca_obj
type(kpca_svd)     :: kpca_nystrom_obj
type(kpca_svd)     :: kpca_rbf_obj
type(kpca_svd)     :: kpca_rbf_nystrom_obj
type(parameters)   :: params
type(cmdline)      :: cline
    integer :: j, nthr_test
real    :: data_ori(NP, NS), avg(NP), tmpvec(NP), data_pca(NP, NS), E_zn(NC, NS), data_cen(NP, NS)
call params%new(cline)
!$ nthr_test = min(4, max(2, omp_get_max_threads()))
!$ call omp_set_num_threads(nthr_test)
!$ params%nthr = nthr_test
!$ if( nthr_test > 1 ) print *, 'Testing kPCA with OpenMP threads:', nthr_test
data_ori(1,:) = [ 1, 2, 3, 4]
data_ori(2,:) = [ 3, 1, 5, 8]
data_ori(3,:) = [-1, 0, 4, 10]
! data_ori(4,:) = [ 0, 0, 7, 10]
! data_ori(5,:) = [-2, 0, 1, 10]
print *, 'Original data:'
print *, data_ori(1,:)
print *, data_ori(2,:)
print *, data_ori(3,:)
avg = sum(data_ori, dim=2) / real(NS)
do j = 1, NP
     data_cen(j,:) = data_ori(j,:) - avg(j)
enddo
call prob_pca%new(NS, NP, NC)
call prob_pca%master(data_cen, MAXPCAITS)
print *, 'PPCA eigenvalues:'
!$omp parallel do private(j,tmpvec) default(shared) proc_bind(close) schedule(static)
do j = 1, NS
    call prob_pca%generate(j, avg, tmpvec)
    data_pca(:,j) = tmpvec
end do
!$omp end parallel do
print *, 'Pre-imaged data using PPCA:'
do j = 1, NP
    print *, data_pca(j,:)
enddo
print *, 'Feature vecs using PPCA:'
do j = 1, NS
    E_zn(:,j) = prob_pca%get_feat(j)
enddo
do j = 1, NC
    print *, E_zn(j, :)
enddo
print *, '---------------------------------------------------'
call mix_prob_pca%new(NS, NP, NC)
call mix_prob_pca%set_params(2, params%nthr)
call mix_prob_pca%master(data_cen, MAXPCAITS)
!$omp parallel do private(j,tmpvec) default(shared) proc_bind(close) schedule(static)
do j = 1, NS
    call mix_prob_pca%generate(j, avg, tmpvec)
    data_pca(:,j) = tmpvec
end do
!$omp end parallel do
print *, 'Pre-imaged data using mPPCA:'
do j = 1, NP
    print *, data_pca(j,:)
enddo
print *, 'Feature vecs using mPPCA:'
do j = 1, NS
    E_zn(:,j) = mix_prob_pca%get_feat(j)
enddo
do j = 1, NC
    print *, E_zn(j, :)
enddo
print *, '---------------------------------------------------'
! PCA test
print *, 'PCA eigenvalues/eigenvectors:'
call pca_obj%new(NS, NP, NC)
call pca_obj%master(data_cen)
print *, 'Feature vecs using PCA:'
do j = 1, NS
    E_zn(:,j) = pca_obj%get_feat(j)
enddo
do j = 1, NC
    print *, E_zn(j, :)
enddo
!$omp parallel do private(j,tmpvec) default(shared) proc_bind(close) schedule(static)
do j = 1, NS
    call pca_obj%generate(j, avg, tmpvec)
    data_pca(:,j) = tmpvec
end do
!$omp end parallel do
print *, 'Pre-imaged data using PCA:'
do j = 1, NP
    print *, data_pca(j,:)
enddo
print *, '---------------------------------------------------'
    ! kPCA test
call kpca_obj%new(NS, NP, NC)
    call kpca_obj%set_params(params%nthr, params%kpca_ker, 'exact', kpca_rbf_gamma=params%kpca_rbf_gamma)
call kpca_obj%master(data_cen)
!$omp parallel do private(j,tmpvec) default(shared) proc_bind(close) schedule(static)
do j = 1, NS
    call kpca_obj%generate(j, avg, tmpvec)
    data_pca(:,j) = tmpvec
end do
!$omp end parallel do
print *, 'Pre-imaged data using kPCA:'
do j = 1, NP
    print *, data_pca(j,:)
enddo
    print *, '---------------------------------------------------'
    ! Nyström kPCA test
    call kpca_nystrom_obj%new(NS, NP, NC)
    call kpca_nystrom_obj%set_params(params%nthr, params%kpca_ker, 'nystrom', NS, params%kpca_rbf_gamma, &
        kpca_cosine_weight_power=params%kpca_cosine_weight_power)
call kpca_nystrom_obj%master(data_cen)
!$omp parallel do private(j,tmpvec) default(shared) proc_bind(close) schedule(static)
do j = 1, NS
    call kpca_nystrom_obj%generate(j, avg, tmpvec)
    data_pca(:,j) = tmpvec
end do
!$omp end parallel do
print *, 'Pre-imaged data using Nyström kPCA:'
do j = 1, NP
    print *, data_pca(j,:)
enddo
print *, '---------------------------------------------------'
! RBF kPCA smoke test
call kpca_rbf_obj%new(NS, NP, NC)
    call kpca_rbf_obj%set_params(params%nthr, 'rbf', 'exact', kpca_rbf_gamma=0.)
call kpca_rbf_obj%master(data_cen)
!$omp parallel do private(j,tmpvec) default(shared) proc_bind(close) schedule(static)
do j = 1, NS
    call kpca_rbf_obj%generate(j, avg, tmpvec)
    data_pca(:,j) = tmpvec
end do
!$omp end parallel do
print *, 'Pre-imaged data using RBF kPCA:'
do j = 1, NP
    print *, data_pca(j,:)
enddo
print *, '---------------------------------------------------'
! RBF Nyström kPCA smoke test
call kpca_rbf_nystrom_obj%new(NS, NP, NC)
    call kpca_rbf_nystrom_obj%set_params(params%nthr, 'rbf', 'nystrom', NS, 0.)
call kpca_rbf_nystrom_obj%master(data_cen)
!$omp parallel do private(j,tmpvec) default(shared) proc_bind(close) schedule(static)
do j = 1, NS
    call kpca_rbf_nystrom_obj%generate(j, avg, tmpvec)
    data_pca(:,j) = tmpvec
end do
!$omp end parallel do
print *, 'Pre-imaged data using RBF Nyström kPCA:'
do j = 1, NP
    print *, data_pca(j,:)
enddo
end program simple_test_pca_all
