program simple_test_eigh
use simple_linalg
implicit none
integer, parameter :: N = 5, N_EIGS = 3, N_L = 15000, N_EIGS_L = 12
type dense_matvec_ctx
    real, pointer :: mat(:,:) => null()
end type dense_matvec_ctx
integer            :: i
integer            :: sparse_info
type(dense_matvec_ctx) :: sparse_ctx
real               :: mat(N, N), eigvals(N_EIGS), eigvecs(N,N_EIGS), tmp(N,N), eig_vals(N), mat_ori(N, N)
real               :: sparse_vals(N_EIGS), sparse_vecs(N,N_EIGS)
real, target       :: mat_sparse(N,N)
mat(:,1) = [ 0.67,-0.20, 0.19,-1.06, 0.46]
mat(:,2) = [-0.20, 3.82,-0.13, 1.06,-0.48]
mat(:,3) = [ 0.19,-0.13, 3.27, 0.11, 1.10]
mat(:,4) = [-1.06, 1.06, 0.11, 5.86,-0.98]
mat(:,5) = [ 0.46,-0.48, 1.10,-0.98, 3.54]
mat_ori  = mat
mat_sparse = mat
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
print *, 'SPARSE_EIGH result:'
sparse_ctx%mat => mat_sparse
call sparse_eigh(test_sparse_matvec, sparse_ctx, N, N_EIGS, sparse_vals, sparse_vecs, tol=1.e-6, max_basis=N, info=sparse_info)
print *, 'Selected sparse eigenvalues', sparse_vals
if( sparse_info /= 0 ) stop 'sparse_eigh returned nonzero info'
if( maxval(abs(sparse_vals - eigvals)) > 1.e-4 ) stop 'sparse_eigh eigenvalue mismatch'
call test_eigh(N_L, N_EIGS_L)

contains

subroutine test_sparse_matvec(ctx, x, y)
    class(*), intent(in)  :: ctx
    real,     intent(in)  :: x(:)
    real,     intent(out) :: y(:)
    select type(ctx)
    type is(dense_matvec_ctx)
        y = matmul(ctx%mat, x)
    class default
        stop 'unsupported sparse_eigh test context'
    end select
end subroutine test_sparse_matvec

end program simple_test_eigh
