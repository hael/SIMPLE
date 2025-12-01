program simple_test_eigh
include 'simple_lib.f08'
integer, parameter :: N = 5, N_EIGS = 3, N_L = 15000, N_EIGS_L = 12
integer            :: i
real               :: mat(N, N), eigvals(N_EIGS), eigvecs(N,N_EIGS), tmp(N,N), eig_vals(N), mat_ori(N, N)
mat(:,1) = [ 0.67,-0.20, 0.19,-1.06, 0.46]
mat(:,2) = [-0.20, 3.82,-0.13, 1.06,-0.48]
mat(:,3) = [ 0.19,-0.13, 3.27, 0.11, 1.10]
mat(:,4) = [-1.06, 1.06, 0.11, 5.86,-0.98]
mat(:,5) = [ 0.46,-0.48, 1.10,-0.98, 3.54]
mat_ori  = mat
print *, 'SVDCMP result:'
call svdcmp(mat, eig_vals, tmp)
print *, eig_vals
print *, 'EIGH result:'
call eigh( N, mat_ori, N_EIGS, eigvals, eigvecs )
print *, 'Selected eigenvalues', eigvals
print *, 'Selected eigenvectors (stored columnwise)'
do i = 1, N
    print *, eigvecs(i,:)
enddo
call test_eigh(N_L, N_EIGS_L)
end program simple_test_eigh
