program simple_test_pca_imgvar
include 'simple_lib.f08'
use simple_ppca_inmem, only: ppca_inmem
use simple_pca_svd,    only: pca_svd
use simple_kpca_svd,   only: kpca_svd
use simple_cmdline,    only: cmdline
use simple_parameters, only: parameters
implicit none
integer, parameter :: NX = 5, NY = 5, NP = NX*NY, NC = 2, MAXPCAITS = 15
type(ppca_inmem)   :: prob_pca
integer :: i, j, cnt
real    :: imgs(NX, NY, NC), flat_img(NP), dist_x(NP), dist_y(NP)
real    :: data_ori(NP, NP), avg(NP), E_zn(NC, NP), data_cen(NP, NP), tmpvec(NP)
real    :: lists(4), result

lists  = [2., 1., 3., 4.]
result = 0.
do i = 1, size(lists)
    do j = 1, size(lists)
        result = result + lists(i) * lists(j)
    enddo
enddo
print *, result

imgs(1,:,1) = [ 0., 2., 3., 4., 0.]
imgs(2,:,1) = [ 0., 1., 5., 0., 0.]
imgs(3,:,1) = [ 0., 0., 0., 0., 0.]
imgs(4,:,1) = [ 0., 0., 0., 0., 0.]
imgs(5,:,1) = [ 0., 0., 0., 0., 0.]

imgs(1,:,2) = [ 0.,  0.,  0.,  0., 0.]
imgs(2,:,2) = [ 0.,  0.,  0.,  0., 0.]
imgs(3,:,2) = [ 0.,  0.,  0.,  0., 0.]
imgs(4,:,2) = [ 0.,  0., -5., -9., 0.]
imgs(5,:,2) = [ 0.,-11., -1., -2., 0.]

cnt = 1
do i = 1, NY
    do j = 1, NX
        flat_img(cnt) = imgs(i,j,1) + imgs(i,j,2)
        dist_x(cnt)   = real(i)
        dist_y(cnt)   = real(j)
        cnt           = cnt + 1
    enddo
enddo

data_ori = 0.
do i = 1, NP
    do j = 1, NP
        if( i == j ) cycle
        data_ori(i,j) = abs(flat_img(i) - flat_img(j)) / abs(real(i - j))
    enddo
enddo
avg = sum(data_ori, dim=2) / real(NP)
do j = 1, NP
    data_cen(j,:) = data_ori(j,:) - avg(j)
enddo
call prob_pca%new(NP, NP, NC)
call prob_pca%master(data_cen, MAXPCAITS)
!$omp parallel do private(j,tmpvec) default(shared) proc_bind(close) schedule(static)
do j = 1, NP
    call prob_pca%generate(j, avg, tmpvec)
end do
!$omp end parallel do
print *, 'Feature vecs using PPCA:'
do j = 1, NP
    E_zn(:,j) = prob_pca%get_feat(j)
enddo
do j = 1, 1
    print *, E_zn(j, 1:5)
    print *, E_zn(j, 6:10)
    print *, E_zn(j,11:15)
    print *, E_zn(j,16:20)
    print *, E_zn(j,21:25)
enddo
do j = 1, 1
    print *, imgs(1,:,1)
    print *, imgs(2,:,1)
    print *, imgs(3,:,1)
    print *, imgs(4,:,1)
    print *, imgs(5,:,1)
enddo

end program simple_test_pca_imgvar
