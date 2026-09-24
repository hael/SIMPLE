!@descr: unit test routines for simple_linalg (eigensolvers, SVD, inversion, least-squares fits, vector helpers)
! Symmetric eigensolvers (dense, sparse, Jacobi), SVD, inversion, least-squares fits and the small vector
! helpers, pinned on matrices with independently computed answers.
module simple_linalg_tester
use simple_test_utils    ! assertions etc.
use simple_defs          ! sp, dp, PI
use simple_string_utils, only: int2str
use simple_linalg
implicit none
private
public :: run_all_linalg_tests

integer, parameter :: N = 5
real,    parameter :: TOL = 1.0e-4 ! single-precision LAPACK on a well-conditioned 5x5
! the LAPACK ssyevr example matrix (symmetric positive definite); reference values from numpy.linalg
real, parameter :: EIGVALS_ASC(N) = [0.43302180, 2.14494666, 3.36808674, 4.27915302, 6.93479178]
real, parameter :: INV_ROW1(N)    = [ 2.22312154, -0.01178349, -0.08938773,  0.37958447, -0.15761941]
real, parameter :: INV_ROW2(N)    = [-0.01178349,  0.27781983,  0.00522982, -0.04844234,  0.02416600]
real, parameter :: INV_ROW3(N)    = [-0.08938773,  0.00522982,  0.34880701, -0.04165620, -0.10759382]
real, parameter :: INV_ROW4(N)    = [ 0.37958447, -0.04844234, -0.04165620,  0.25340427,  0.02720249]
real, parameter :: INV_ROW5(N)    = [-0.15761941,  0.02416600, -0.10759382,  0.02720249,  0.34720798]
real, parameter :: TRACE_REF      = 17.16

! context handed through sparse_eigh's class(*) argument: a dense matrix behind the matvec callback
type :: dense_ctx
    real :: mat(N,N)
end type dense_ctx

contains

    subroutine run_all_linalg_tests()
        write(*,'(A)') '**** running all linalg tests ****'
        call test_eigh_largest()
        call test_eigh_smallest()
        call test_sparse_eigh()
        call test_svdcmp()
        call test_matinv()
        call test_jacobi_eigsrt()
        call test_svdfit_polynomial()
        call test_svd_multifit_plane()
        call test_fit_straight_line()
        call test_fit_straight_line_recovery()
        call test_plane_fits()
        call test_vector_helpers()
        call test_gemm_tn()
    end subroutine run_all_linalg_tests

    !---------------- fixtures ----------------

    function example_matrix() result( mat )
        real :: mat(N,N)
        mat(:,1) = [ 0.67,-0.20, 0.19,-1.06, 0.46]
        mat(:,2) = [-0.20, 3.82,-0.13, 1.06,-0.48]
        mat(:,3) = [ 0.19,-0.13, 3.27, 0.11, 1.10]
        mat(:,4) = [-1.06, 1.06, 0.11, 5.86,-0.98]
        mat(:,5) = [ 0.46,-0.48, 1.10,-0.98, 3.54]
    end function example_matrix

    function example_inverse() result( inv )
        real :: inv(N,N)
        inv(1,:) = INV_ROW1
        inv(2,:) = INV_ROW2
        inv(3,:) = INV_ROW3
        inv(4,:) = INV_ROW4
        inv(5,:) = INV_ROW5
    end function example_inverse

    subroutine dense_matvec( ctx, x, y )
        class(*), intent(in)  :: ctx
        real,     intent(in)  :: x(:)
        real,     intent(out) :: y(:)
        select type(ctx)
            type is(dense_ctx)
                y = matmul(ctx%mat, x)
            class default
                y = 0.
        end select
    end subroutine dense_matvec

    ! eigenpairs of the example matrix: unit columns, mutually orthogonal, A v = lambda v
    subroutine assert_eigenpairs( mat, eigvals, eigvecs, label )
        real,             intent(in) :: mat(N,N), eigvals(:), eigvecs(:,:)
        character(len=*), intent(in) :: label
        real    :: resid(N)
        integer :: i, j
        do i = 1,size(eigvals)
            call assert_real(1.0, sqrt(sum(eigvecs(:,i)**2)), TOL, label//': eigenvector '//int2str(i)//' has unit norm')
            do j = i+1,size(eigvals)
                call assert_real(0.0, dot_product(eigvecs(:,i), eigvecs(:,j)), TOL, &
                    &label//': eigenvectors '//int2str(i)//' and '//int2str(j)//' are orthogonal')
            enddo
            resid = matmul(mat, eigvecs(:,i)) - eigvals(i) * eigvecs(:,i)
            call assert_true(maxval(abs(resid)) < TOL * max(1.0, abs(eigvals(i))), &
                &label//': A v = lambda v for eigenpair '//int2str(i))
        enddo
    end subroutine assert_eigenpairs

    !---------------- eigensolvers ----------------

    subroutine test_eigh_largest()
        integer, parameter :: NEIGS = 3
        real :: mat(N,N), work(N,N), eigvals(NEIGS), eigvecs(N,NEIGS)
        integer :: i
        write(*,'(A)') 'test_eigh_largest'
        mat  = example_matrix()
        work = mat
        call eigh(N, work, NEIGS, eigvals, eigvecs)
        do i = 1,NEIGS
            call assert_real(EIGVALS_ASC(N-NEIGS+i), eigvals(i), TOL, 'eigh: largest eigenvalues in ascending order, '//int2str(i))
        enddo
        call assert_eigenpairs(mat, eigvals, eigvecs, 'eigh largest')
    end subroutine test_eigh_largest

    subroutine test_eigh_smallest()
        integer, parameter :: NEIGS = 3
        real :: mat(N,N), work(N,N), eigvals(NEIGS), eigvecs(N,NEIGS)
        integer :: i
        write(*,'(A)') 'test_eigh_smallest'
        mat  = example_matrix()
        work = mat
        call eigh(N, work, NEIGS, eigvals, eigvecs, smallest=.true.)
        do i = 1,NEIGS
            call assert_real(EIGVALS_ASC(i), eigvals(i), TOL, 'eigh: smallest eigenvalues in ascending order, '//int2str(i))
        enddo
        call assert_eigenpairs(mat, eigvals, eigvecs, 'eigh smallest')
    end subroutine test_eigh_smallest

    ! the ARPACK path agrees with the dense one when the operator is the same matrix
    subroutine test_sparse_eigh()
        integer, parameter :: NEIGS = 3
        type(dense_ctx) :: ctx
        real    :: eigvals(NEIGS), eigvecs(N,NEIGS)
        integer :: info, i
        write(*,'(A)') 'test_sparse_eigh'
        ctx%mat = example_matrix()
        info    = -1
        call sparse_eigh(dense_matvec, ctx, N, NEIGS, eigvals, eigvecs, tol=1.e-6, max_basis=N, info=info)
        call assert_int(0, info, 'sparse_eigh: info is zero')
        do i = 1,NEIGS
            call assert_real(EIGVALS_ASC(N-NEIGS+i), eigvals(i), TOL, 'sparse_eigh: largest eigenvalues in ascending order, '//int2str(i))
        enddo
        call assert_eigenpairs(ctx%mat, eigvals, eigvecs, 'sparse_eigh')
    end subroutine test_sparse_eigh

    subroutine test_jacobi_eigsrt()
        real    :: a(N,N), mat(N,N), d(N), v(N,N)
        integer :: nrot, i
        write(*,'(A)') 'test_jacobi_eigsrt'
        mat  = example_matrix()
        a    = mat
        nrot = -1
        call jacobi(a, N, N, d, v, nrot)
        call assert_int(0, nrot, 'jacobi: the ssyev-backed replacement reports no rotations (nrot = 0)')
        call eigsrt(d, v, N, N)
        do i = 1,N
            call assert_real(EIGVALS_ASC(N+1-i), d(i), TOL, 'jacobi + eigsrt: eigenvalue '//int2str(i)//' in descending order')
        enddo
        call assert_eigenpairs(mat, d, v, 'jacobi')
    end subroutine test_jacobi_eigsrt

    !---------------- SVD and inversion ----------------

    subroutine test_svdcmp()
        integer, parameter :: M = 4, K = 3
        real    :: a(N,N), u(N,N), w(N), v(N,N), rec(N,N), wsorted(N)
        real    :: b(M,K), ub(M,K), wb(K), vb(K,K), recb(M,K)
        integer :: i, j
        write(*,'(A)') 'test_svdcmp'
        a = example_matrix()
        u = a
        call svdcmp(u, w, v)
        ! SPD matrix: the singular values are the eigenvalues
        wsorted = w
        call sort_desc(wsorted)
        do i = 1,N
            call assert_real(EIGVALS_ASC(N+1-i), wsorted(i), TOL, 'svdcmp: singular value '//int2str(i)//' of an SPD matrix is its eigenvalue')
        enddo
        do j = 1,N
            rec(:,j) = matmul(u, w * v(j,:))
        enddo
        call assert_true(maxval(abs(rec - a)) < TOL, 'svdcmp: U diag(w) V^T reconstructs the square matrix')
        call assert_true(maxval(abs(matmul(transpose(v), v) - identity(N))) < TOL, 'svdcmp: V is orthogonal')
        ! a rectangular matrix
        do j = 1,K
            do i = 1,M
                b(i,j) = real(i) + 0.5 * real(j) * real(i - 2) + 0.1 * real(j)**2
            enddo
        enddo
        ub = b
        call svdcmp(ub, wb, vb)
        do j = 1,K
            recb(:,j) = matmul(ub, wb * vb(j,:))
        enddo
        call assert_true(maxval(abs(recb - b)) < TOL, 'svdcmp: reconstruction of a 4x3 matrix')
        call assert_true(all(wb >= 0.), 'svdcmp: singular values are non-negative')
    end subroutine test_svdcmp

    subroutine test_matinv()
        real    :: a(N,N), inv(N,N), sing(3,3), inv3(3,3)
        integer :: errflg
        write(*,'(A)') 'test_matinv'
        a = example_matrix()
        call matinv(a, inv, N, errflg)
        call assert_int(0, errflg, 'matinv: no error on a regular matrix')
        call assert_true(maxval(abs(inv - example_inverse())) < TOL, 'matinv: inverse equals the reference')
        call assert_true(maxval(abs(matmul(a, inv) - identity(N))) < TOL, 'matinv: A A^-1 = I')
        sing = 1.0 ! rank one
        errflg = 0
        call matinv(sing, inv3, 3, errflg)
        call assert_int(1, errflg, 'matinv: a singular matrix is flagged')
    end subroutine test_matinv

    !---------------- least-squares fits ----------------

    subroutine test_svdfit_polynomial()
        integer, parameter :: NPTS = 20, MA = 3
        real    :: x(NPTS), y(NPTS), sig(NPTS), a(MA), v(MA,MA), w(MA), chisq
        integer :: i
        write(*,'(A)') 'test_svdfit_polynomial'
        do i = 1,NPTS
            x(i) = 0.1 * real(i - 1)
            y(i) = 1.0 + 2.0 * x(i) + 3.0 * x(i)**2
        enddo
        sig = 1.0
        call svdfit(x, y, sig, a, v, w, chisq, poly3)
        call assert_real(1.0, a(1), 1.e-3, 'svdfit: constant coefficient')
        call assert_real(2.0, a(2), 1.e-3, 'svdfit: linear coefficient')
        call assert_real(3.0, a(3), 1.e-3, 'svdfit: quadratic coefficient')
        call assert_true(chisq < 1.e-4, 'svdfit: exact data gives chi^2 ~ 0')
        ! noisy data: the fit passes through the noise with a chi^2 of the noise size
        y = y + 0.01 * [( real(mod(3*i, 7) - 3), i=1,NPTS )]
        call svdfit(x, y, sig, a, v, w, chisq, poly3)
        call assert_real(2.0, a(2), 5.e-2, 'svdfit: linear coefficient survives small noise')
        call assert_real(7.3e-3, chisq, 1.e-3, 'svdfit: chi^2 of the noisy fit is the residual noise (numpy lstsq: 7.32e-3)')
    end subroutine test_svdfit_polynomial

    subroutine test_svd_multifit_plane()
        integer, parameter :: NPTS = 12, MA = 3
        real    :: x(2,NPTS), y(NPTS), sig(NPTS), a(MA), v(MA,MA), w(MA), chisq
        integer :: i
        write(*,'(A)') 'test_svd_multifit_plane'
        do i = 1,NPTS
            x(1,i) = real(mod(i, 4))
            x(2,i) = real(i / 4)
            y(i)   = 0.5 + 1.5 * x(1,i) - 2.0 * x(2,i)
        enddo
        sig = 1.0
        call svd_multifit(x, y, sig, a, v, w, chisq, affine2)
        call assert_real( 0.5, a(1), 1.e-3, 'svd_multifit: intercept')
        call assert_real( 1.5, a(2), 1.e-3, 'svd_multifit: first slope')
        call assert_real(-2.0, a(3), 1.e-3, 'svd_multifit: second slope')
        call assert_true(chisq < 1.e-4, 'svd_multifit: exact data gives chi^2 ~ 0')
    end subroutine test_svd_multifit_plane

    subroutine test_fit_straight_line()
        integer, parameter :: NPTS = 10
        real    :: datavec(NPTS,2), slope, intercept, corr
        integer :: i
        write(*,'(A)') 'test_fit_straight_line'
        do i = 1,NPTS
            datavec(i,1) = real(i)
            datavec(i,2) = 2.0 * real(i) + 1.0
        enddo
        call fit_straight_line(NPTS, datavec, slope, intercept, corr)
        call assert_real(2.0, slope,     1.e-5, 'fit_straight_line: slope of an exact line')
        call assert_real(1.0, intercept, 1.e-5, 'fit_straight_line: intercept of an exact line')
        call assert_real(1.0, corr,      1.e-5, 'fit_straight_line: corr (r squared) of an exact line is 1')
        ! symmetric perturbation: same slope and intercept, r squared below 1
        datavec(2,2) = datavec(2,2) + 0.5
        datavec(9,2) = datavec(9,2) - 0.5
        call fit_straight_line(NPTS, datavec, slope, intercept, corr)
        call assert_real(0.99890, corr,   1.e-4, 'fit_straight_line: r squared of the perturbed line (161.5^2 / (82.5 * 316.5))')
        call assert_real(1.95758, slope,  1.e-4, 'fit_straight_line: slope of the perturbed line (161.5 / 82.5)')
        call assert_real(1.23333, intercept, 1.e-4, 'fit_straight_line: intercept of the perturbed line (12 - slope * 5.5)')
    end subroutine test_fit_straight_line

    !> 35 exact lines on 100 points of x in [-1, 1) (single precision, as a caller hands them in),
    !! slopes from -5 to 5 through near-flat and flat, intercepts from -10 to 10: slope and intercept
    !! come back within 1e-5 (a float32 emulation of the fit gives 6e-8 at worst). r squared is not
    !! asserted: for a near-flat line it is 0/0 in single precision (below 0.9999 for |slope| under
    !! about 5e-6 |intercept|), and the one production caller (guinier_bfac) uses only the slope.
    !! Replaces the unit_numerics sub-suite `straight-line fit`, which drew 10 000 random lines and
    !! required r squared >= 0.9999, so it failed for about one seed in twenty.
    subroutine test_fit_straight_line_recovery()
        integer, parameter :: NPTS = 100
        real,    parameter :: SLOPES(7)     = [-5., -0.5, -1.e-5, 0., 1.e-5, 0.5, 5.]
        real,    parameter :: INTERCEPTS(5) = [-10., -1., 0., 3., 10.]
        real    :: datavec(NPTS,2), slope, intercept, corr, x, err_slope, err_intercept
        integer :: i, j, k
        write(*,'(A)') 'test_fit_straight_line_recovery'
        err_slope     = 0.
        err_intercept = 0.
        do i = 1,size(SLOPES)
            do j = 1,size(INTERCEPTS)
                x = -1.
                do k = 1,NPTS
                    datavec(k,1) = x
                    datavec(k,2) = SLOPES(i) * x + INTERCEPTS(j)
                    x = x + 0.02
                enddo
                call fit_straight_line(NPTS, datavec, slope, intercept, corr)
                err_slope     = max(err_slope,     abs(slope     - SLOPES(i)))
                err_intercept = max(err_intercept, abs(intercept - INTERCEPTS(j)))
            enddo
        enddo
        call assert_real(0., err_slope,     1.e-5, 'fit_straight_line: slope of 35 exact lines, steep to flat')
        call assert_real(0., err_intercept, 1.e-5, 'fit_straight_line: intercept of 35 exact lines, steep to flat')
    end subroutine test_fit_straight_line_recovery

    subroutine test_plane_fits()
        integer, parameter :: NPTS = 9
        real    :: xyz(NPTS,3), points(3,NPTS), a, b, c, sol(3), line(4,3)
        logical :: err
        integer :: i
        write(*,'(A)') 'test_plane_fits'
        do i = 1,NPTS
            xyz(i,1) = real(mod(i-1, 3)) - 1.0
            xyz(i,2) = real((i-1) / 3) - 1.0
            xyz(i,3) = 0.5 * xyz(i,1) - 0.25 * xyz(i,2) + 2.0
        enddo
        call fit_lsq_plane(NPTS, xyz, a, b, c, err)
        call assert_false(err, 'fit_lsq_plane: nine points on a plane fit without error')
        call assert_real( 0.50, a, 1.e-4, 'fit_lsq_plane: x slope')
        call assert_real(-0.25, b, 1.e-4, 'fit_lsq_plane: y slope')
        call assert_real( 2.00, c, 1.e-4, 'fit_lsq_plane: offset')
        points = transpose(xyz)
        sol = plane_from_points(points)
        call assert_real( 0.50, sol(1), 1.e-4, 'plane_from_points: x slope')
        call assert_real(-0.25, sol(2), 1.e-4, 'plane_from_points: y slope')
        call assert_real( 2.00, sol(3), 1.e-4, 'plane_from_points: offset')
        ! collinear points do not define a plane
        do i = 1,4
            line(i,:) = [real(i), 2.0*real(i), 3.0*real(i)]
        enddo
        call fit_lsq_plane(4, line, a, b, c, err)
        call assert_true(err, 'fit_lsq_plane: collinear points are flagged')
    end subroutine test_plane_fits

    !---------------- vector helpers ----------------

    subroutine test_vector_helpers()
        real :: v(3), w(3), m22a(2,2), m22b(2,2), mat(N,N)
        write(*,'(A)') 'test_vector_helpers'
        call assert_real(5.0, euclid([0.,0.], [3.,4.]),  1.e-6, 'euclid: 3-4-5 in 1D arrays')
        m22a = 0.; m22b = reshape([3.,0.,0.,4.], [2,2])
        call assert_real(5.0, euclid(m22a, m22b),       1.e-6, 'euclid: 3-4-5 in 2D arrays')
        call assert_real(5.0, real(euclid([0.d0,0.d0], [3.d0,4.d0])), 1.e-6, 'euclid: double precision')
        call assert_real(5.0, hyp(3.,4.),               1.e-6, 'hyp: two reals')
        call assert_real(3.0, hyp(1.,2.,2.),            1.e-6, 'hyp: three reals')
        call assert_real(5.0, hyp(3,4),                 1.e-6, 'hyp: two integers')
        call assert_real(3.0, hyp(1,2,2),               1.e-6, 'hyp: three integers')
        call assert_real(5.0, arg([3.,4.]),             1.e-6, 'arg: vector length')
        call assert_real(5.0, norm_2([3.,4.]),          1.e-6, 'norm_2: vector length (returned 0 through Accelerate snrm2 until 2026-09-23)')
        call assert_real(0.0, myacos(1.5),              1.e-6, 'myacos: clamps above 1')
        call assert_real(PI,  myacos(-2.0),             1.e-6, 'myacos: clamps below -1')
        call assert_real(PI/3., myacos(0.5),            1.e-6, 'myacos: plain acos inside [-1,1]')
        v = [1.,0.,0.]; w = [0.,1.,0.]
        call assert_real(PI/2., vector_angle_norm(v, w), 1.e-6, 'vector_angle_norm: orthogonal unit vectors')
        call assert_real(0.0,   vector_angle_norm(v, v), 1.e-6, 'vector_angle_norm: a vector with itself')
        mat = example_matrix()
        call assert_real(TRACE_REF, trace(mat),         1.e-5, 'trace of the example matrix')
        call assert_int(3, ang2vox(2.6, 1.3),                  'ang2vox: int(ang/smpd) + 1')
        call assert_real(2.6, vox2ang(3, 1.3),          1.e-6, 'vox2ang: (vox - 1) smpd')
        call assert_real(180., rad2deg(PI),             1.e-4, 'rad2deg')
        call assert_real(PI,   deg2rad(180.),           1.e-6, 'deg2rad')
    end subroutine test_vector_helpers

    subroutine test_gemm_tn()
        real :: a(3,2), b(3,4), c(2,4)
        integer :: i, j
        write(*,'(A)') 'test_gemm_tn'
        do j = 1,2
            do i = 1,3
                a(i,j) = real(i) + 10. * real(j)
            enddo
        enddo
        do j = 1,4
            do i = 1,3
                b(i,j) = real(i * j) - 2.5
            enddo
        enddo
        call gemm_tn(a, b, c)
        call assert_true(maxval(abs(c - matmul(transpose(a), b))) < 1.e-4, 'gemm_tn: C = A^T B')
    end subroutine test_gemm_tn

    !---------------- helpers ----------------

    function identity( n ) result( eye )
        integer, intent(in) :: n
        real    :: eye(n,n)
        integer :: i
        eye = 0.
        do i = 1,n
            eye(i,i) = 1.
        enddo
    end function identity

    subroutine sort_desc( arr )
        real, intent(inout) :: arr(:)
        real    :: tmp
        integer :: i, j
        do i = 1,size(arr)-1
            do j = i+1,size(arr)
                if( arr(j) > arr(i) )then
                    tmp = arr(i); arr(i) = arr(j); arr(j) = tmp
                endif
            enddo
        enddo
    end subroutine sort_desc

    function poly3( x, n ) result( f )
        real,    intent(in) :: x
        integer, intent(in) :: n
        real :: f(n)
        f(1) = 1.0
        if( n >= 2 ) f(2) = x
        if( n >= 3 ) f(3) = x * x
    end function poly3

    function affine2( x, n ) result( f )
        real,    intent(in) :: x(:)
        integer, intent(in) :: n
        real :: f(n)
        f(1) = 1.0
        if( n >= 2 ) f(2) = x(1)
        if( n >= 3 ) f(3) = x(2)
    end function affine2

end module simple_linalg_tester
