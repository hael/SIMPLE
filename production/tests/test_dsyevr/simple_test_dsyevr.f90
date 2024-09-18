program simple_test_dsyevr
include 'simple_lib.f08'
use simple_dsyevr, only: s_dsyevr
integer, parameter :: N = 5, NSELECT = 3
integer, parameter :: LDA = N, LDZ = N
integer, parameter :: LWMAX = 1000
integer  :: info, lwork, liwork, il, iu, m, isuppz(N), iwork(LWMAX), i
real(dp) :: abstol, vl, vu, mat(LDA, N), w(N), z(LDZ, NSELECT), work(LWMAX)
mat(:,1) = [ 0.67, 0.00, 0.00, 0.00, 0.00]
mat(:,2) = [-0.20, 3.82, 0.00, 0.00, 0.00]
mat(:,3) = [ 0.19,-0.13, 3.27, 0.00, 0.00]
mat(:,4) = [-1.06, 1.06, 0.11, 5.86, 0.00]
mat(:,5) = [ 0.46,-0.48, 1.10,-0.98, 3.54]
print *, 'DSYEVR Example Program Results'
abstol = -1.0
il     = 1
iu     = NSELECT
! Query the optimal workspace.
lwork  = -1
liwork = -1
call s_dsyevr( 'Vectors', 'Indices', 'Upper', N, mat, LDA, vl, vu, il, iu, abstol,&
            & m, w, z, LDZ, isuppz, work, lwork, iwork, liwork, info )
lwork  = min( LWMAX, int( work( 1 ) ) )
liwork = min( LWMAX, iwork( 1 ) )
! Solve eigenproblem.
call s_dsyevr( 'Vectors', 'Indices', 'Upper', N, mat, LDA, vl, vu, il, iu, abstol,&
            & m, w, z, LDZ, isuppz, work, lwork, iwork, liwork, info )
! Check for convergence.
if( info > 0 )then
    print *, 'The algorithm failed to compute eigenvalues.'
    STOP
endif
print *, 'The total number of eigenvalues found:', m
print *, 'Selected eigenvalues', w
print *, 'Selected eigenvectors (stored columnwise)'
do i = 1, N
    print *, z(i,:)
enddo
end program simple_test_dsyevr