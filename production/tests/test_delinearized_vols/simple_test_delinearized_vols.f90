program simple_test_delinearized_vols
include 'simple_lib.f08'
implicit none
integer :: ivol, ivar, j, errflg, var_inds(3), dist_ind
real    :: truths(4,3), vols(4,3), probs(3,3), probs_inv(3,3), truths_inv(4,3), avg_vol(4), var(3), probs_dist(3,3), var_sorted(3)
truths(:,1) = [1., 2., 3., 4.]
truths(:,2) = [.5, 2., 3., 6.]
truths(:,3) = [.8, 2., 3., 5.]
probs(1,:)  = [0.8,  0.1,  0.1]
probs(2,:)  = [0.12, 0.7,  0.18]
probs(3,:)  = [0.04, 0.06, 0.9]
vols        = 0.
do ivol = 1, 3
    do j = 1, 3
        vols(:,ivol) = vols(:,ivol) + probs(ivol,j) * truths(:,j)
    enddo
    print *, vols(:,ivol)
enddo
avg_vol = sum(vols, dim=2)/3.
print *, avg_vol
do ivol = 1, 3
    var(ivol) = sum((vols(:,ivol) - avg_vol)**2)
enddo
var_inds   = [1, 2, 3]
var_sorted = var
call hpsort(var_sorted, var_inds)
print *, 'var_sorted = ', var_sorted
print *, 'var_inds   = ', var_inds
! distribution guess
do ivol = 1, 3
    dist_ind = 0
    probs_dist(ivol,ivol) = exp(-(dist_ind - 0.)**2/2./var(ivol))/sqrt(2. * PI * var(ivol))
    dist_ind = dist_ind + 1
    do ivar = 1, 3
        if( var_inds(ivar) == ivol ) cycle
        probs_dist(ivol,var_inds(ivar)) = exp(-(dist_ind - 0.)**2/2./var(ivol))/sqrt(2. * PI * var(ivol))
        dist_ind = dist_ind + 1
    enddo
    probs_dist(ivol,:) = probs_dist(ivol,:) / sum(probs_dist(ivol,:))
    print *, 'probs', ivol, ' = ', probs_dist(ivol,:)
enddo
call matinv(probs_dist, probs_inv, 3, errflg)
! recovering
truths_inv = 0.
do ivol = 1, 3
    do j = 1, 3
        truths_inv(:,ivol) = truths_inv(:,ivol) + probs_inv(ivol,j) * vols(:,j)
    enddo
    print *, truths_inv(:,ivol)
enddo
end program simple_test_delinearized_vols
