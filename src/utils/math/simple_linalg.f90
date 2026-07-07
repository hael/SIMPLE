!@descr: SIMPLE linear algebra helpers with BLAS/LAPACK-backed kernels.
module simple_linalg
use simple_defs
use simple_is_check_assert
implicit none
private

public :: arg, eigsrt, eigh, euclid, fit_lsq_plane, fit_straight_line
public :: hyp, jacobi, l1dist, matinv, myacos, norm_2, outerprod
public :: plane_from_points, projz, pythag, rad2deg, deg2rad
public :: same_energy_euclid, svbksb, svdcmp, svdfit, svd_multifit
public :: svd_solve, normal_solve, qr_solve, svdvar, test_eigh
public :: trace, vabs, vector_angle_norm, vox2ang, ang2vox
public :: sparse_eigh
public :: hermitian_eigh, hermitian_invert, hermitian_solve

abstract interface
    subroutine sparse_matvec_sp_proc(ctx, x, y)
        class(*), intent(in)  :: ctx
        real,     intent(in)  :: x(:)
        real,     intent(out) :: y(:)
    end subroutine sparse_matvec_sp_proc
end interface

interface eigsrt
    module procedure eigsrt_sp, eigsrt_dp
end interface eigsrt

interface eigh
    module procedure eigh_sp
end interface eigh

interface hermitian_eigh
    module procedure hermitian_eigh_z
end interface hermitian_eigh

interface hermitian_invert
    module procedure hermitian_invert_dp
    module procedure hermitian_invert_z
end interface hermitian_invert

interface hermitian_solve
    module procedure hermitian_solve_dp
    module procedure hermitian_solve_z
end interface hermitian_solve

interface euclid
    module procedure euclid_sp_1, euclid_sp_2, euclid_dp
end interface euclid

interface hyp
    module procedure hyp_1, hyp_2, hyp_3, hyp_4
end interface hyp

interface jacobi
    module procedure jacobi_sp, jacobi_dp
end interface jacobi

interface l1dist
    module procedure l1dist_sp, l1dist_dp
end interface l1dist

interface matinv
    module procedure matinv_sp, matinv_dp
end interface matinv

interface myacos
    module procedure myacos_sp, myacos_dp
end interface myacos

interface norm_2
    module procedure norm_2_sp, norm_2_dp
end interface norm_2

interface outerprod
    module procedure outerprod_r, outerprod_d
end interface outerprod

interface pythag
    module procedure pythag_sp, pythag_dp
end interface pythag

interface rad2deg
    module procedure rad2deg_1, rad2deg_2
end interface rad2deg

interface deg2rad
    module procedure deg2rad_sp, deg2rad_dp
end interface deg2rad

interface svbksb
    module procedure svbksb_sp, svbksb_dp
end interface svbksb

interface svdcmp
    module procedure svdcmp_sp, svdcmp_dp
end interface svdcmp

interface svdfit
    module procedure svdfit_sp, svdfit_dp
end interface svdfit

interface svd_multifit
    module procedure svd_multifit_sp, svd_multifit_dp
end interface svd_multifit

interface vabs
    module procedure vabs_sp, vabs_dp
end interface vabs

interface
    function snrm2(n, x, incx) result(res)
        integer(kind=4), intent(in) :: n, incx
        real(kind=4), intent(in) :: x(*)
        real(kind=4) :: res
    end function snrm2

    function dnrm2(n, x, incx) result(res)
        integer(kind=4), intent(in) :: n, incx
        real(kind=8), intent(in) :: x(*)
        real(kind=8) :: res
    end function dnrm2

    subroutine sgelsy(m, n, nrhs, a, lda, b, ldb, jpvt, rcond, rank, work, lwork, info)
        integer(kind=4), intent(in) :: m, n, nrhs, lda, ldb, lwork
        integer(kind=4), intent(inout) :: jpvt(*)
        integer(kind=4), intent(out) :: rank, info
        real(kind=4), intent(inout) :: a(lda,*), b(ldb,*), work(*)
        real(kind=4), intent(in) :: rcond
    end subroutine sgelsy

    subroutine dgelsy(m, n, nrhs, a, lda, b, ldb, jpvt, rcond, rank, work, lwork, info)
        integer(kind=4), intent(in) :: m, n, nrhs, lda, ldb, lwork
        integer(kind=4), intent(inout) :: jpvt(*)
        integer(kind=4), intent(out) :: rank, info
        real(kind=8), intent(inout) :: a(lda,*), b(ldb,*), work(*)
        real(kind=8), intent(in) :: rcond
    end subroutine dgelsy

    subroutine dgelss(m, n, nrhs, a, lda, b, ldb, s, rcond, rank, work, lwork, info)
        integer(kind=4), intent(in) :: m, n, nrhs, lda, ldb, lwork
        integer(kind=4), intent(out) :: rank, info
        real(kind=8), intent(inout) :: a(lda,*), b(ldb,*), work(*)
        real(kind=8), intent(out) :: s(*)
        real(kind=8), intent(in) :: rcond
    end subroutine dgelss

    subroutine dposv(uplo, n, nrhs, a, lda, b, ldb, info)
        character(len=1), intent(in) :: uplo
        integer(kind=4), intent(in) :: n, nrhs, lda, ldb
        integer(kind=4), intent(out) :: info
        real(kind=8), intent(inout) :: a(lda,*), b(ldb,*)
    end subroutine dposv

    subroutine zposv(uplo, n, nrhs, a, lda, b, ldb, info)
        character(len=1), intent(in) :: uplo
        integer(kind=4), intent(in) :: n, nrhs, lda, ldb
        integer(kind=4), intent(out) :: info
        complex(kind=8), intent(inout) :: a(lda,*), b(ldb,*)
    end subroutine zposv

    subroutine sgesvd(jobu, jobvt, m, n, a, lda, s, u, ldu, vt, ldvt, work, lwork, info)
        character(len=1), intent(in) :: jobu, jobvt
        integer(kind=4), intent(in) :: m, n, lda, ldu, ldvt, lwork
        integer(kind=4), intent(out) :: info
        real(kind=4), intent(inout) :: a(lda,*), work(*)
        real(kind=4), intent(out) :: s(*), u(ldu,*), vt(ldvt,*)
    end subroutine sgesvd

    subroutine dgesvd(jobu, jobvt, m, n, a, lda, s, u, ldu, vt, ldvt, work, lwork, info)
        character(len=1), intent(in) :: jobu, jobvt
        integer(kind=4), intent(in) :: m, n, lda, ldu, ldvt, lwork
        integer(kind=4), intent(out) :: info
        real(kind=8), intent(inout) :: a(lda,*), work(*)
        real(kind=8), intent(out) :: s(*), u(ldu,*), vt(ldvt,*)
    end subroutine dgesvd

    subroutine sgetrf(m, n, a, lda, ipiv, info)
        integer(kind=4), intent(in) :: m, n, lda
        integer(kind=4), intent(out) :: ipiv(*), info
        real(kind=4), intent(inout) :: a(lda,*)
    end subroutine sgetrf

    subroutine dgetrf(m, n, a, lda, ipiv, info)
        integer(kind=4), intent(in) :: m, n, lda
        integer(kind=4), intent(out) :: ipiv(*), info
        real(kind=8), intent(inout) :: a(lda,*)
    end subroutine dgetrf

    subroutine sgetri(n, a, lda, ipiv, work, lwork, info)
        integer(kind=4), intent(in) :: n, lda, lwork
        integer(kind=4), intent(in) :: ipiv(*)
        integer(kind=4), intent(out) :: info
        real(kind=4), intent(inout) :: a(lda,*), work(*)
    end subroutine sgetri

    subroutine dgetri(n, a, lda, ipiv, work, lwork, info)
        integer(kind=4), intent(in) :: n, lda, lwork
        integer(kind=4), intent(in) :: ipiv(*)
        integer(kind=4), intent(out) :: info
        real(kind=8), intent(inout) :: a(lda,*), work(*)
    end subroutine dgetri

    subroutine ssyev(jobz, uplo, n, a, lda, w, work, lwork, info)
        character(len=1), intent(in) :: jobz, uplo
        integer(kind=4), intent(in) :: n, lda, lwork
        integer(kind=4), intent(out) :: info
        real(kind=4), intent(inout) :: a(lda,*), work(*)
        real(kind=4), intent(out) :: w(*)
    end subroutine ssyev

    subroutine dsyev(jobz, uplo, n, a, lda, w, work, lwork, info)
        character(len=1), intent(in) :: jobz, uplo
        integer(kind=4), intent(in) :: n, lda, lwork
        integer(kind=4), intent(out) :: info
        real(kind=8), intent(inout) :: a(lda,*), work(*)
        real(kind=8), intent(out) :: w(*)
    end subroutine dsyev

    subroutine zheev(jobz, uplo, n, a, lda, w, work, lwork, rwork, info)
        character(len=1), intent(in) :: jobz, uplo
        integer(kind=4), intent(in) :: n, lda, lwork
        integer(kind=4), intent(out) :: info
        complex(kind=8), intent(inout) :: a(lda,*), work(*)
        real(kind=8), intent(out) :: w(*), rwork(*)
    end subroutine zheev

    subroutine ssyevr(jobz, range, uplo, n, a, lda, vl, vu, il, iu, abstol, &
        m, w, z, ldz, isuppz, work, lwork, iwork, liwork, info)
        character(len=1), intent(in) :: jobz, range, uplo
        integer(kind=4), intent(in) :: n, il, iu, ldz, lda, lwork, liwork
        integer(kind=4), intent(out) :: m, isuppz(*), iwork(*), info
        real(kind=4), intent(in) :: abstol, vl, vu
        real(kind=4), intent(inout) :: a(lda,*)
        real(kind=4), intent(out) :: work(*), z(ldz,*), w(*)
    end subroutine ssyevr

    subroutine ssaupd(ido, bmat, n, which, nev, tol, resid, ncv, v, ldv, &
        iparam, ipntr, workd, workl, lworkl, info)
        character(len=1), intent(in) :: bmat
        character(len=2), intent(in) :: which
        integer(kind=4), intent(in) :: n, nev, ncv, ldv, lworkl
        integer(kind=4), intent(inout) :: ido, info
        integer(kind=4), intent(inout) :: iparam(*), ipntr(*)
        real(kind=4), intent(in) :: tol
        real(kind=4), intent(inout) :: resid(*), v(ldv,*), workd(*), workl(*)
    end subroutine ssaupd

    subroutine sseupd(rvec, howmny, select, d, z, ldz, sigma, bmat, n, which, &
        nev, tol, resid, ncv, v, ldv, iparam, ipntr, workd, workl, lworkl, info)
        logical, intent(in) :: rvec
        character(len=1), intent(in) :: howmny, bmat
        character(len=2), intent(in) :: which
        integer(kind=4), intent(in) :: n, nev, ncv, ldz, ldv, lworkl
        integer(kind=4), intent(inout) :: iparam(*), ipntr(*)
        integer(kind=4), intent(out) :: info
        logical, intent(inout) :: select(*)
        real(kind=4), intent(in) :: tol, sigma
        real(kind=4), intent(inout) :: resid(*), v(ldv,*), workd(*), workl(*)
        real(kind=4), intent(out) :: d(*), z(ldz,*)
    end subroutine sseupd
end interface

contains

pure function arg(vec) result(length)
    real(kind=4), intent(in) :: vec(:)
    real(kind=4) :: length
    length = sqrt(sum(vec * vec))
end function arg

subroutine eigsrt_sp(d, v, n, np)
    integer, intent(in) :: n, np
    real(kind=4), intent(inout) :: d(np), v(np,np)
    integer :: i, j, k
    real(kind=4) :: p
    do i = 1, n - 1
        k = i
        p = d(i)
        do j = i + 1, n
            if(d(j) > p)then
                k = j
                p = d(j)
            endif
        enddo
        if(k /= i)then
            d(k) = d(i)
            d(i) = p
            do j = 1, n
                p = v(j,i)
                v(j,i) = v(j,k)
                v(j,k) = p
            enddo
        endif
    enddo
end subroutine eigsrt_sp

subroutine eigsrt_dp(d, v, n, np)
    integer, intent(in) :: n, np
    real(kind=8), intent(inout) :: d(np), v(np,np)
    integer :: i, j, k
    real(kind=8) :: p
    do i = 1, n - 1
        k = i
        p = d(i)
        do j = i + 1, n
            if(d(j) > p)then
                k = j
                p = d(j)
            endif
        enddo
        if(k /= i)then
            d(k) = d(i)
            d(i) = p
            do j = 1, n
                p = v(j,i)
                v(j,i) = v(j,k)
                v(j,k) = p
            enddo
        endif
    enddo
end subroutine eigsrt_dp

subroutine fit_lsq_plane(n, xyz, a, b, c, err)
    integer, intent(in) :: n
    real(kind=4), intent(in) :: xyz(n,3)
    real(kind=4), intent(out) :: a, b, c
    logical, intent(out) :: err
    real(kind=8) :: sx, sy, sz, sxx, syy, sxy, sxz, syz, denom, rn
    err = .false.
    a = 0.0
    b = 0.0
    c = 0.0
    rn = real(n, dp)
    sx = sum(real(xyz(:,1), dp))
    sy = sum(real(xyz(:,2), dp))
    sz = sum(real(xyz(:,3), dp))
    sxx = sum(real(xyz(:,1), dp)**2)
    syy = sum(real(xyz(:,2), dp)**2)
    sxy = sum(real(xyz(:,1), dp) * real(xyz(:,2), dp))
    sxz = sum(real(xyz(:,1), dp) * real(xyz(:,3), dp))
    syz = sum(real(xyz(:,2), dp) * real(xyz(:,3), dp))
    denom = sx*sx*syy - 2.d0*sxy*sx*sy + sxx*sy*sy + rn*(sxy*sxy - sxx*syy)
    if(abs(denom) < 1.d-10)then
        err = .true.
        return
    endif
    a = real((sy*sy*sxz - syy*rn*sxz + sxy*rn*syz + sx*syy*sz - sy*(sx*syz + sxy*sz)) / denom, sp)
    b = real((sxy*rn*sxz + sx*sx*syz - sxx*rn*syz + sxx*sy*sz - sx*(sy*sxz + sxy*sz)) / denom, sp)
    c = real((sx*syy*sxz - sxy*sy*sxz - sxy*sx*syz + sxx*sy*syz + sz*(sxy*sxy - sxx*syy)) / denom, sp)
end subroutine fit_lsq_plane

subroutine fit_straight_line(n, datavec, slope, intercept, corr)
    integer, intent(in) :: n
    real(kind=4), intent(in) :: datavec(n,2)
    real(kind=4), intent(out) :: slope, intercept, corr
    real(kind=8) :: ave_x, ave_y, ss_xx, ss_yy, ss_xy, x, y, dn
    integer :: i
    ave_x = 0.d0
    ave_y = 0.d0
    ss_xx = 0.d0
    ss_yy = 0.d0
    ss_xy = 0.d0
    do i = 1, n
        x = datavec(i,1)
        y = datavec(i,2)
        ave_x = ave_x + x
        ave_y = ave_y + y
        ss_xx = ss_xx + x*x
        ss_yy = ss_yy + y*y
        ss_xy = ss_xy + x*y
    enddo
    dn = dble(n)
    ave_x = ave_x / dn
    ave_y = ave_y / dn
    ss_xx = ss_xx - dn*ave_x*ave_x
    ss_yy = ss_yy - dn*ave_y*ave_y
    ss_xy = ss_xy - dn*ave_x*ave_y
    slope = real(ss_xy / ss_xx, sp)
    intercept = real(ave_y - slope*ave_x, sp)
    corr = real((ss_xy*ss_xy) / (ss_xx*ss_yy), sp)
end subroutine fit_straight_line

subroutine jacobi_sp(a, n, np, d, v, nrot)
    integer, intent(in) :: n, np
    real(kind=4), intent(inout) :: a(np,np), v(np,np), d(np)
    integer, intent(inout) :: nrot
    integer :: info, lwork
    real(kind=4), allocatable :: work(:), amat(:,:)
    real(kind=4) :: work_query(1)
    allocate(amat(n,n))
    amat = a(1:n,1:n)
    lwork = -1
    call ssyev('V', 'U', n, amat, n, d, work_query, lwork, info)
    if(info /= 0) call lapack_stop('JACOBI workspace query', 'SSYEV', info)
    lwork = max(1, int(work_query(1)))
    allocate(work(lwork))
    call ssyev('V', 'U', n, amat, n, d, work, lwork, info)
    if(info /= 0) call lapack_stop('JACOBI', 'SSYEV', info)
    v = 0.0_sp
    v(1:n,1:n) = amat
    if(np > n) d(n+1:np) = 0.0_sp
    nrot = 0
end subroutine jacobi_sp

subroutine jacobi_dp(a, n, np, d, v, nrot)
    integer, intent(in) :: n, np
    real(kind=8), intent(inout) :: a(np,np), v(np,np), d(np)
    integer, intent(inout) :: nrot
    integer :: info, lwork
    real(kind=8), allocatable :: work(:), amat(:,:)
    real(kind=8) :: work_query(1)
    allocate(amat(n,n))
    amat = a(1:n,1:n)
    lwork = -1
    call dsyev('V', 'U', n, amat, n, d, work_query, lwork, info)
    if(info /= 0) call lapack_stop('JACOBI workspace query', 'DSYEV', info)
    lwork = max(1, int(work_query(1)))
    allocate(work(lwork))
    call dsyev('V', 'U', n, amat, n, d, work, lwork, info)
    if(info /= 0) call lapack_stop('JACOBI', 'DSYEV', info)
    v = 0.0_dp
    v(1:n,1:n) = amat
    if(np > n) d(n+1:np) = 0.0_dp
    nrot = 0
end subroutine jacobi_dp

subroutine matinv_sp(matrix, inverse, n, errflg)
    integer, intent(in) :: n
    real(kind=4), intent(in) :: matrix(n,n)
    real(kind=4), intent(out) :: inverse(n,n)
    integer, intent(out) :: errflg
    integer :: info, lwork
    integer, allocatable :: ipiv(:)
    real(kind=4), allocatable :: work(:)
    real(kind=4) :: work_query(1)
    errflg = 0
    inverse = matrix
    allocate(ipiv(n))
    call sgetrf(n, n, inverse, n, ipiv, info)
    if(info /= 0)then
        errflg = 1
        return
    endif
    lwork = -1
    call sgetri(n, inverse, n, ipiv, work_query, lwork, info)
    if(info /= 0)then
        errflg = 1
        return
    endif
    lwork = max(1, int(work_query(1)))
    allocate(work(lwork))
    call sgetri(n, inverse, n, ipiv, work, lwork, info)
    if(info /= 0) errflg = 1
end subroutine matinv_sp

subroutine matinv_dp(matrix, inverse, n, errflg)
    integer, intent(in) :: n
    real(kind=8), intent(in) :: matrix(n,n)
    real(kind=8), intent(out) :: inverse(n,n)
    integer, intent(out) :: errflg
    integer :: info, lwork
    integer, allocatable :: ipiv(:)
    real(kind=8), allocatable :: work(:)
    real(kind=8) :: work_query(1)
    errflg = 0
    inverse = matrix
    allocate(ipiv(n))
    call dgetrf(n, n, inverse, n, ipiv, info)
    if(info /= 0)then
        errflg = 1
        return
    endif
    lwork = -1
    call dgetri(n, inverse, n, ipiv, work_query, lwork, info)
    if(info /= 0)then
        errflg = 1
        return
    endif
    lwork = max(1, int(work_query(1)))
    allocate(work(lwork))
    call dgetri(n, inverse, n, ipiv, work, lwork, info)
    if(info /= 0) errflg = 1
end subroutine matinv_dp

function norm_2_sp(v) result(r)
    real(kind=4), intent(in) :: v(:)
    real(kind=4) :: r
    r = snrm2(size(v), v, 1)
end function norm_2_sp

function norm_2_dp(v) result(r)
    real(kind=8), intent(in) :: v(:)
    real(kind=8) :: r
    r = dnrm2(size(v), v, 1)
end function norm_2_dp

function outerprod_r(a, b)
    real(kind=4), intent(in) :: a(:), b(:)
    real(kind=4) :: outerprod_r(size(a),size(b))
    outerprod_r = spread(a, dim=2, ncopies=size(b)) * spread(b, dim=1, ncopies=size(a))
end function outerprod_r

function outerprod_d(a, b)
    real(kind=8), intent(in) :: a(:), b(:)
    real(kind=8) :: outerprod_d(size(a),size(b))
    outerprod_d = spread(a, dim=2, ncopies=size(b)) * spread(b, dim=1, ncopies=size(a))
end function outerprod_d

function plane_from_points(points) result(sol)
    real(kind=4), intent(inout) :: points(:,:)
    real(kind=4) :: sol(3)
    real(kind=4) :: m(size(points, dim=2),3), b(size(points, dim=2))
    real(kind=4) :: prod(3,3), prod_inv(3,3), prod1(3,size(points, dim=2))
    integer :: errflg, p, npts
    sol = 0.0_sp
    if(size(points, dim=1) /= 3)then
        write(logfhandle,*) 'Need to input points in 3D!; plane_from_points'
        return
    endif
    if(size(points, dim=2) < 3)then
        write(logfhandle,*) 'Not enough input points for fitting!; plane_from_points'
        return
    endif
    npts = size(points, dim=2)
    do p = 1, npts
        m(p,1) = points(1,p)
        m(p,2) = points(2,p)
        m(p,3) = 1.0_sp
        b(p) = points(3,p)
    enddo
    prod = matmul(transpose(m), m)
    call matinv(prod, prod_inv, 3, errflg)
    if(errflg /= 0)then
        write(logfhandle,*) 'Couldn t find inverse matrix! ;plane_from_points'
        stop
    endif
    prod1 = matmul(prod_inv, transpose(m))
    sol = matmul(prod1, b)
end function plane_from_points

subroutine projz(vec3, vec2)
    real(kind=4), intent(in) :: vec3(3)
    real(kind=4), intent(out) :: vec2(2)
    vec2(1) = vec3(1)
    vec2(2) = vec3(2)
end subroutine projz

subroutine svbksb_sp(u, w, v, b, x)
    real(kind=4), intent(in) :: u(:,:), v(:,:), w(:), b(:)
    real(kind=4), intent(out) :: x(:)
    real(kind=4) :: tmp(size(x))
    integer :: j
    tmp = 0.0_sp
    do j = 1, size(x)
        if(w(j) /= 0.0_sp) tmp(j) = dot_product(u(:,j), b) / w(j)
    enddo
    x = matmul(v, tmp)
end subroutine svbksb_sp

subroutine svbksb_dp(u, w, v, b, x)
    real(kind=8), intent(in) :: u(:,:), v(:,:), w(:), b(:)
    real(kind=8), intent(out) :: x(:)
    real(kind=8) :: tmp(size(x))
    integer :: j
    tmp = 0.0_dp
    do j = 1, size(x)
        if(w(j) /= 0.0_dp) tmp(j) = dot_product(u(:,j), b) / w(j)
    enddo
    x = matmul(v, tmp)
end subroutine svbksb_dp

subroutine svdcmp_sp(a, w, v)
    real(kind=4), intent(inout) :: a(:,:)
    real(kind=4), intent(out) :: w(:), v(:,:)
    integer :: m, n, minmn, info, lwork
    real(kind=4), allocatable :: acopy(:,:), u(:,:), vt(:,:), work(:)
    real(kind=4) :: work_query(1)
    m = size(a, 1)
    n = size(a, 2)
    minmn = min(m, n)
    allocate(acopy(m,n), u(m,minmn), vt(n,n))
    acopy = a
    lwork = -1
    call sgesvd('S', 'A', m, n, acopy, m, w, u, m, vt, n, work_query, lwork, info)
    if(info /= 0) call lapack_stop('SVDCMP workspace query', 'SGESVD', info)
    lwork = max(1, int(work_query(1)))
    allocate(work(lwork))
    acopy = a
    call sgesvd('S', 'A', m, n, acopy, m, w, u, m, vt, n, work, lwork, info)
    if(info /= 0) call lapack_stop('SVDCMP', 'SGESVD', info)
    a = 0.0_sp
    a(:,1:minmn) = u(:,1:minmn)
    v = transpose(vt)
    if(size(w) > minmn) w(minmn+1:) = 0.0_sp
end subroutine svdcmp_sp

subroutine svdcmp_dp(a, w, v)
    real(kind=8), intent(inout) :: a(:,:)
    real(kind=8), intent(out) :: w(:), v(:,:)
    integer :: m, n, minmn, info, lwork
    real(kind=8), allocatable :: acopy(:,:), u(:,:), vt(:,:), work(:)
    real(kind=8) :: work_query(1)
    m = size(a, 1)
    n = size(a, 2)
    minmn = min(m, n)
    allocate(acopy(m,n), u(m,minmn), vt(n,n))
    acopy = a
    lwork = -1
    call dgesvd('S', 'A', m, n, acopy, m, w, u, m, vt, n, work_query, lwork, info)
    if(info /= 0) call lapack_stop('SVDCMP workspace query', 'DGESVD', info)
    lwork = max(1, int(work_query(1)))
    allocate(work(lwork))
    acopy = a
    call dgesvd('S', 'A', m, n, acopy, m, w, u, m, vt, n, work, lwork, info)
    if(info /= 0) call lapack_stop('SVDCMP', 'DGESVD', info)
    a = 0.0_dp
    a(:,1:minmn) = u(:,1:minmn)
    v = transpose(vt)
    if(size(w) > minmn) w(minmn+1:) = 0.0_dp
end subroutine svdcmp_dp

subroutine svdfit_sp(x, y, sig, a, v, w, chisq, funcs)
    real(kind=4), intent(in) :: x(:), y(:), sig(:)
    real(kind=4), intent(out) :: a(:), v(:,:), w(:), chisq
    interface
        function funcs(x, n)
            real(kind=4), intent(in) :: x
            integer, intent(in) :: n
            real(kind=4) :: funcs(n)
        end function funcs
    end interface
    integer :: i, ma, n
    real(kind=4) :: b(size(x)), sigi(size(x)), u(size(x),size(a)), usav(size(x),size(a))
    real(kind=4), parameter :: tol = 1.0e-5_sp
    n = assert_eq(size(x), size(y), size(sig), 'svdfit_sp: n')
    ma = assert_eq(size(a), size(v,1), size(v,2), size(w), 'svdfit_sp: ma')
    sigi = 1.0_sp / sig
    b = y * sigi
    do i = 1, n
        usav(i,:) = funcs(x(i), ma)
    enddo
    u = usav * spread(sigi, dim=2, ncopies=ma)
    usav = u
    call svdcmp(u, w, v)
    where(w < tol * maxval(w)) w = 0.0_sp
    call svbksb(u, w, v, b, a)
    chisq = vabs(matmul(usav, a) - b)**2
end subroutine svdfit_sp

subroutine svdfit_dp(x, y, sig, a, v, w, chisq, funcs)
    real(kind=8), intent(in) :: x(:), y(:), sig(:)
    real(kind=8), intent(out) :: a(:), v(:,:), w(:), chisq
    interface
        function funcs(x, n)
            real(kind=8), intent(in) :: x
            integer, intent(in) :: n
            real(kind=8) :: funcs(n)
        end function funcs
    end interface
    integer :: i, ma, n
    real(kind=8) :: b(size(x)), sigi(size(x)), u(size(x),size(a)), usav(size(x),size(a))
    real(kind=8), parameter :: tol = 1.0e-14_dp
    n = assert_eq(size(x), size(y), size(sig), 'svdfit_dp: n')
    ma = assert_eq(size(a), size(v,1), size(v,2), size(w), 'svdfit_dp: ma')
    sigi = 1.0_dp / sig
    b = y * sigi
    do i = 1, n
        usav(i,:) = funcs(x(i), ma)
    enddo
    u = usav * spread(sigi, dim=2, ncopies=ma)
    usav = u
    call svdcmp(u, w, v)
    where(w < tol * maxval(w)) w = 0.0_dp
    call svbksb(u, w, v, b, a)
    chisq = vabs(matmul(usav, a) - b)**2
end subroutine svdfit_dp

subroutine svd_multifit_sp(x, y, sig, a, v, w, chisq, funcs)
    real(kind=4), intent(in) :: x(:,:), y(:), sig(:)
    real(kind=4), intent(out) :: a(:), v(:,:), w(:), chisq
    interface
        function funcs(x, n)
            real(kind=4), intent(in) :: x(:)
            integer, intent(in) :: n
            real(kind=4) :: funcs(n)
        end function funcs
    end interface
    integer :: i, ma, n
    real(kind=4) :: b(size(sig)), sigi(size(sig)), u(size(x,dim=2),size(a))
    real(kind=4) :: usav(size(x,dim=2),size(a))
    real(kind=4), parameter :: tol = 1.0e-5_sp
    n = assert_eq(size(x, dim=2), size(y), size(sig), 'svd_multifit_sp: n')
    ma = assert_eq(size(a), size(v,1), size(v,2), size(w), 'svd_multifit_sp: ma')
    sigi = 1.0_sp / sig
    b = y * sigi
    do i = 1, n
        usav(i,:) = funcs(x(:,i), ma)
    enddo
    u = usav * spread(sigi, dim=2, ncopies=ma)
    usav = u
    call svdcmp(u, w, v)
    where(w < tol * maxval(w)) w = 0.0_sp
    call svbksb(u, w, v, b, a)
    chisq = vabs(matmul(usav, a) - b)**2
end subroutine svd_multifit_sp

subroutine svd_multifit_dp(x, y, sig, a, v, w, chisq, funcs)
    real(kind=8), intent(in) :: x(:,:), y(:), sig(:)
    real(kind=8), intent(out) :: a(:), v(:,:), w(:), chisq
    interface
        function funcs(x, n)
            real(kind=8), intent(in) :: x(:)
            integer, intent(in) :: n
            real(kind=8) :: funcs(n)
        end function funcs
    end interface
    integer :: i, ma, n
    real(kind=8) :: b(size(sig)), sigi(size(sig)), u(size(x,dim=2),size(a))
    real(kind=8) :: usav(size(x,dim=2),size(a))
    real(kind=8), parameter :: tol = 1.0e-14_dp
    n = assert_eq(size(x, dim=2), size(y), size(sig), 'svd_multifit_dp: n')
    ma = assert_eq(size(a), size(v,1), size(v,2), size(w), 'svd_multifit_dp: ma')
    sigi = 1.0_dp / sig
    b = y * sigi
    do i = 1, n
        usav(i,:) = funcs(x(:,i), ma)
    enddo
    u = usav * spread(sigi, dim=2, ncopies=ma)
    usav = u
    call svdcmp(u, w, v)
    where(w < tol * maxval(w)) w = 0.0_dp
    call svbksb(u, w, v, b, a)
    chisq = vabs(matmul(usav, a) - b)**2
end subroutine svd_multifit_dp

subroutine svdvar(v, w, cvm)
    real(kind=4), intent(in) :: v(:,:), w(:)
    real(kind=4), intent(out) :: cvm(:,:)
    integer :: ma
    real(kind=4) :: wti(size(w))
    ma = assert_eq((/size(v,1), size(v,2), size(w), size(cvm,1), size(cvm,2)/), 'svdvar')
    where(is_equal(w, 0.0))
        wti = 0.0_sp
    elsewhere
        wti = 1.0_sp / (w * w)
    endwhere
    cvm = v * spread(wti, dim=1, ncopies=ma)
    cvm = matmul(cvm, transpose(v))
end subroutine svdvar

pure function trace(mat) result(tr)
    real(kind=4), intent(in) :: mat(:,:)
    real(kind=4) :: tr
    integer :: i
    tr = 0.0_sp
    do i = 1, min(size(mat,1), size(mat,2))
        tr = tr + mat(i,i)
    enddo
end function trace

function vabs_sp(v)
    real(kind=4), intent(in) :: v(:)
    real(kind=4) :: vabs_sp
    vabs_sp = snrm2(size(v), v, 1)
end function vabs_sp

function vabs_dp(v)
    real(kind=8), intent(in) :: v(:)
    real(kind=8) :: vabs_dp
    vabs_dp = dnrm2(size(v), v, 1)
end function vabs_dp

pure function vector_angle_norm(v, w) result(eulerdist)
    real(kind=4), intent(in) :: v(3), w(3)
    real(kind=4) :: eulerdist
    eulerdist = acos(max(-1.0_sp, min(1.0_sp, dot_product(v, w))))
end function vector_angle_norm

pure function myacos_sp(arg) result(r)
    real(kind=4), intent(in) :: arg
    real(kind=4) :: r, x, y
    x = min(1.0_sp, abs(arg))
    y = sign(x, arg)
    r = acos(y)
end function myacos_sp

pure function myacos_dp(arg) result(r)
    real(kind=8), intent(in) :: arg
    real(kind=8) :: r, x, y
    x = min(1.0_dp, abs(arg))
    y = sign(x, arg)
    r = acos(y)
end function myacos_dp

elemental function deg2rad_sp(deg) result(rad)
    real(kind=4), intent(in) :: deg
    real(kind=4) :: rad
    rad = (deg / 180.0_sp) * PI
end function deg2rad_sp

elemental function deg2rad_dp(deg) result(rad)
    real(kind=8), intent(in) :: deg
    real(kind=8) :: rad
    rad = (deg / 180.0_dp) * DPI
end function deg2rad_dp

elemental function rad2deg_1(rad) result(deg)
    real(kind=4), intent(in) :: rad
    real(kind=4) :: deg
    deg = (rad / PI) * 180.0_sp
end function rad2deg_1

elemental function rad2deg_2(rad) result(deg)
    real(kind=8), intent(in) :: rad
    real(kind=8) :: deg
    deg = (rad / DPI) * 180.0_dp
end function rad2deg_2

elemental function ang2vox(ang, smpd) result(vox)
    real(kind=4), intent(in) :: ang, smpd
    integer :: vox
    vox = int(ang / smpd) + 1
end function ang2vox

elemental function vox2ang(vox, smpd) result(ang)
    integer, intent(in) :: vox
    real(kind=4), intent(in) :: smpd
    real(kind=4) :: ang
    ang = (vox - 1.0_sp) * smpd
end function vox2ang

real(kind=4) pure function hyp_1(x1, x2)
    real(kind=4), intent(in) :: x1, x2
    hyp_1 = sqrt(x1*x1 + x2*x2)
end function hyp_1

real(kind=4) pure function hyp_2(x1, x2, x3)
    real(kind=4), intent(in) :: x1, x2, x3
    hyp_2 = sqrt(x1*x1 + x2*x2 + x3*x3)
end function hyp_2

real(kind=4) pure function hyp_3(x1, x2)
    integer, intent(in) :: x1, x2
    hyp_3 = sqrt(real(x1*x1 + x2*x2, sp))
end function hyp_3

real(kind=4) pure function hyp_4(x1, x2, x3)
    integer, intent(in) :: x1, x2, x3
    hyp_4 = sqrt(real(x1*x1 + x2*x2 + x3*x3, sp))
end function hyp_4

pure function euclid_sp_1(vec1, vec2) result(dist)
    real(kind=4), intent(in) :: vec1(:), vec2(:)
    real(kind=4) :: dist
    dist = sqrt(sum((vec1 - vec2)**2))
end function euclid_sp_1

pure function euclid_sp_2(vec1, vec2) result(dist)
    real(kind=4), intent(in) :: vec1(:,:), vec2(:,:)
    real(kind=4) :: dist
    dist = sqrt(sum((vec1 - vec2)**2))
end function euclid_sp_2

pure function euclid_dp(vec1, vec2) result(dist)
    real(kind=8), intent(in) :: vec1(:), vec2(:)
    real(kind=8) :: dist
    dist = sqrt(sum((vec1 - vec2)**2))
end function euclid_dp

pure function l1dist_sp(vec1, vec2) result(dist)
    real(kind=4), intent(in) :: vec1(:), vec2(:)
    real(kind=4) :: dist
    dist = sum(abs(vec1 - vec2))
end function l1dist_sp

pure function l1dist_dp(vec1, vec2) result(dist)
    real(kind=8), intent(in) :: vec1(:), vec2(:)
    real(kind=8) :: dist
    dist = sum(abs(vec1 - vec2))
end function l1dist_dp

function same_energy_euclid(vec1, vec2) result(dist)
    real(kind=4), intent(in) :: vec1(:), vec2(:)
    real(kind=4) :: avg1, avg2, dist
    real(kind=4), allocatable :: diff1(:), diff2(:)
    integer :: sz1, sz2
    sz1 = size(vec1)
    sz2 = size(vec2)
    allocate(diff1(sz1), diff2(sz2))
    avg1 = sum(vec1) / real(sz1, sp)
    avg2 = sum(vec2) / real(sz2, sp)
    diff1 = vec1 - avg1
    diff2 = vec2 - avg2
    dist = euclid(diff1, diff2)
end function same_energy_euclid

function pythag_sp(a, b)
    real(kind=4), intent(in) :: a, b
    real(kind=4) :: pythag_sp
    pythag_sp = sqrt(a*a + b*b)
end function pythag_sp

function pythag_dp(a, b)
    real(kind=8), intent(in) :: a, b
    real(kind=8) :: pythag_dp
    pythag_dp = sqrt(a*a + b*b)
end function pythag_dp

subroutine eigh_sp(n, mat, neigs, eigvals, eigvecs, smallest)
    integer, intent(in) :: n, neigs
    real(kind=4), intent(inout) :: mat(n,n)
    real(kind=4), intent(out) :: eigvals(neigs), eigvecs(n,neigs)
    logical, optional, intent(in) :: smallest
    integer, allocatable :: iwork(:), isuppz(:)
    real(kind=4), allocatable :: work(:)
    integer :: info, il, iu, m, lwork, liwork
    real(kind=4) :: vl, vu
    logical :: l_smallest
    l_smallest = .false.
    vl = 0.0_sp
    vu = 0.0_sp
    if(present(smallest)) l_smallest = smallest
    if(l_smallest)then
        il = 1
        iu = neigs
    else
        il = n - neigs + 1
        iu = n
    endif
    lwork = max(1, 26 * n)
    liwork = max(1, 10 * n)
    allocate(work(lwork), iwork(liwork), isuppz(2 * max(1, neigs)))
    call ssyevr('V', 'I', 'L', n, mat, n, vl, vu, il, iu, -1.0_sp, &
        m, eigvals, eigvecs, n, isuppz, work, lwork, iwork, liwork, info)
    if(info /= 0) call lapack_stop('EIGH', 'SSYEVR', info)
end subroutine eigh_sp

subroutine normal_solve(m, n, a, b, x, flag)
    integer(kind=4), intent(in) :: m, n
    real(kind=8), intent(in) :: a(m,n), b(m)
    real(kind=8), intent(out) :: x(n)
    integer(kind=4), intent(out) :: flag
    integer(kind=4) :: info
    real(kind=8), allocatable :: ata(:,:), atb(:,:)
    flag = 0
    if(m < n .or. m <= 0 .or. n <= 0)then
        flag = 1
        return
    endif
    allocate(ata(n,n), atb(n,1))
    ata = matmul(transpose(a), a)
    atb(:,1) = matmul(transpose(a), b)
    call dposv('U', n, 1, ata, n, atb, n, info)
    if(info /= 0)then
        flag = 1
        return
    endif
    x = atb(:,1)
end subroutine normal_solve

subroutine qr_solve(m, n, a, b, x)
    integer(kind=4), intent(in) :: m, n
    real(kind=8), intent(in) :: a(m,n), b(m)
    real(kind=8), intent(out) :: x(n)
    integer(kind=4), allocatable :: jpvt(:)
    integer(kind=4) :: info, lda, ldb, lwork, nrhs, rank
    real(kind=8), allocatable :: a_qr(:,:), rhs(:,:), work(:)
    real(kind=8) :: rcond, work_query(1)
    if(m <= 0 .or. n <= 0) call lapack_stop('QR_SOLVE', 'DGELSY', -1)
    lda = max(1, m)
    ldb = max(1, m, n)
    nrhs = 1
    rcond = epsilon(1.0_dp)
    allocate(a_qr(lda,n), rhs(ldb,nrhs), jpvt(n))
    a_qr = a
    rhs = 0.0_dp
    rhs(1:m,1) = b
    jpvt = 0
    lwork = -1
    call dgelsy(m, n, nrhs, a_qr, lda, rhs, ldb, jpvt, rcond, rank, work_query, lwork, info)
    if(info /= 0) call lapack_stop('QR_SOLVE workspace query', 'DGELSY', info)
    lwork = max(1, int(work_query(1)))
    allocate(work(lwork))
    a_qr = a
    rhs = 0.0_dp
    rhs(1:m,1) = b
    jpvt = 0
    call dgelsy(m, n, nrhs, a_qr, lda, rhs, ldb, jpvt, rcond, rank, work, lwork, info)
    if(info /= 0) call lapack_stop('QR_SOLVE', 'DGELSY', info)
    x = rhs(1:n,1)
end subroutine qr_solve

subroutine svd_solve(m, n, a, b, x)
    integer(kind=4), intent(in) :: m, n
    real(kind=8), intent(in) :: a(m,n), b(m)
    real(kind=8), intent(out) :: x(n)
    integer(kind=4) :: info, lda, ldb, lwork, nrhs, rank
    real(kind=8), allocatable :: a_copy(:,:), rhs(:,:), s(:), work(:)
    real(kind=8) :: rcond, work_query(1)
    if(m <= 0 .or. n <= 0) call lapack_stop('SVD_SOLVE', 'DGELSS', -1)
    lda = max(1, m)
    ldb = max(1, m, n)
    nrhs = 1
    rcond = epsilon(1.0_dp)
    allocate(a_copy(lda,n), rhs(ldb,nrhs), s(min(m,n)))
    a_copy = a
    rhs = 0.0_dp
    rhs(1:m,1) = b
    lwork = -1
    call dgelss(m, n, nrhs, a_copy, lda, rhs, ldb, s, rcond, rank, work_query, lwork, info)
    if(info /= 0) call lapack_stop('SVD_SOLVE workspace query', 'DGELSS', info)
    lwork = max(1, int(work_query(1)))
    allocate(work(lwork))
    a_copy = a
    rhs = 0.0_dp
    rhs(1:m,1) = b
    call dgelss(m, n, nrhs, a_copy, lda, rhs, ldb, s, rcond, rank, work, lwork, info)
    if(info /= 0) call lapack_stop('SVD_SOLVE', 'DGELSS', info)
    x = rhs(1:n,1)
end subroutine svd_solve

subroutine hermitian_solve_dp(a, b, x, flag)
    real(kind=8), intent(in)  :: a(:,:), b(:)
    real(kind=8), intent(out) :: x(:)
    integer(kind=4), optional, intent(out) :: flag
    real(kind=8), allocatable :: acopy(:,:), rhs(:,:)
    integer(kind=4) :: info, n
    if(present(flag)) flag = 0
    n = size(b)
    x = 0.0_dp
    if(size(a,1) /= n .or. size(a,2) /= n .or. size(x) /= n .or. n <= 0)then
        if(present(flag))then
            flag = 1
            return
        endif
        call lapack_stop('HERMITIAN_SOLVE', 'DPOSV', -1)
    endif
    allocate(acopy(n,n), rhs(n,1))
    acopy = 0.5_dp * (a + transpose(a))
    rhs(:,1) = b
    call dposv('U', n, 1, acopy, n, rhs, n, info)
    if(info /= 0)then
        if(present(flag))then
            flag = info
            return
        endif
        call lapack_stop('HERMITIAN_SOLVE', 'DPOSV', info)
    endif
    x = rhs(:,1)
end subroutine hermitian_solve_dp

subroutine hermitian_invert_dp(a, ainv, flag)
    real(kind=8), intent(in)  :: a(:,:)
    real(kind=8), intent(out) :: ainv(:,:)
    integer(kind=4), optional, intent(out) :: flag
    real(kind=8), allocatable :: acopy(:,:), rhs(:,:)
    integer(kind=4) :: info, n, i
    if(present(flag)) flag = 0
    n = size(a,1)
    ainv = 0.0_dp
    if(size(a,2) /= n .or. size(ainv,1) /= n .or. size(ainv,2) /= n .or. n <= 0)then
        if(present(flag))then
            flag = 1
            return
        endif
        call lapack_stop('HERMITIAN_INVERT', 'DPOSV', -1)
    endif
    allocate(acopy(n,n), rhs(n,n), source=0.0_dp)
    acopy = 0.5_dp * (a + transpose(a))
    do i = 1,n
        rhs(i,i) = 1.0_dp
    end do
    call dposv('U', n, n, acopy, n, rhs, n, info)
    if(info /= 0)then
        if(present(flag))then
            flag = info
            return
        endif
        call lapack_stop('HERMITIAN_INVERT', 'DPOSV', info)
    endif
    ainv = 0.5_dp * (rhs + transpose(rhs))
end subroutine hermitian_invert_dp

subroutine hermitian_solve_z(a, b, x, flag)
    complex(kind=8), intent(in)  :: a(:,:), b(:)
    complex(kind=8), intent(out) :: x(:)
    integer(kind=4), optional, intent(out) :: flag
    complex(kind=8), allocatable :: acopy(:,:), rhs(:,:)
    integer(kind=4) :: info, n, i
    if(present(flag)) flag = 0
    n = size(b)
    x = cmplx(0.0_dp, 0.0_dp, kind=8)
    if(size(a,1) /= n .or. size(a,2) /= n .or. size(x) /= n .or. n <= 0)then
        if(present(flag))then
            flag = 1
            return
        endif
        call lapack_stop('HERMITIAN_SOLVE', 'ZPOSV', -1)
    endif
    allocate(acopy(n,n), rhs(n,1))
    acopy = 0.5_dp * (a + transpose(conjg(a)))
    do i = 1,n
        acopy(i,i) = cmplx(real(acopy(i,i), dp), 0.0_dp, kind=8)
    end do
    rhs(:,1) = b
    call zposv('U', n, 1, acopy, n, rhs, n, info)
    if(info /= 0)then
        if(present(flag))then
            flag = info
            return
        endif
        call lapack_stop('HERMITIAN_SOLVE', 'ZPOSV', info)
    endif
    x = rhs(:,1)
end subroutine hermitian_solve_z

subroutine hermitian_invert_z(a, ainv, flag)
    complex(kind=8), intent(in)  :: a(:,:)
    complex(kind=8), intent(out) :: ainv(:,:)
    integer(kind=4), optional, intent(out) :: flag
    complex(kind=8), allocatable :: acopy(:,:), rhs(:,:)
    integer(kind=4) :: info, n, i
    if(present(flag)) flag = 0
    n = size(a,1)
    ainv = cmplx(0.0_dp, 0.0_dp, kind=8)
    if(size(a,2) /= n .or. any(shape(ainv) /= [n,n]) .or. n <= 0)then
        if(present(flag))then
            flag = 1
            return
        endif
        call lapack_stop('HERMITIAN_INVERT', 'ZPOSV', -1)
    endif
    allocate(acopy(n,n), rhs(n,n), source=cmplx(0.0_dp, 0.0_dp, kind=8))
    acopy = 0.5_dp * (a + transpose(conjg(a)))
    do i = 1,n
        acopy(i,i) = cmplx(real(acopy(i,i), dp), 0.0_dp, kind=8)
        rhs(i,i) = cmplx(1.0_dp, 0.0_dp, kind=8)
    end do
    call zposv('U', n, n, acopy, n, rhs, n, info)
    if(info /= 0)then
        if(present(flag))then
            flag = info
            return
        endif
        call lapack_stop('HERMITIAN_INVERT', 'ZPOSV', info)
    endif
    ainv = 0.5_dp * (rhs + transpose(conjg(rhs)))
    do i = 1,n
        ainv(i,i) = cmplx(real(ainv(i,i), dp), 0.0_dp, kind=8)
    end do
end subroutine hermitian_invert_z

subroutine hermitian_eigh_z(a, eigvals, eigvecs, flag)
    complex(kind=8), intent(in)  :: a(:,:)
    real(kind=8),    intent(out) :: eigvals(:)
    complex(kind=8), intent(out) :: eigvecs(:,:)
    integer(kind=4), optional, intent(out) :: flag
    complex(kind=8), allocatable :: acopy(:,:), work(:), vec_tmp(:)
    complex(kind=8) :: work_query(1)
    real(kind=8), allocatable :: rwork(:)
    real(kind=8) :: eval_tmp
    integer(kind=4) :: info, lwork, n, i, j
    if(present(flag)) flag = 0
    n = size(a,1)
    eigvals = 0.0_dp
    eigvecs = cmplx(0.0_dp, 0.0_dp, kind=8)
    if(size(a,2) /= n .or. size(eigvals) < n .or. size(eigvecs,1) < n .or. size(eigvecs,2) < n .or. n <= 0)then
        if(present(flag))then
            flag = 1
            return
        endif
        call lapack_stop('HERMITIAN_EIGH', 'ZHEEV', -1)
    endif
    allocate(acopy(n,n), rwork(max(1, 3*n - 2)))
    acopy = 0.5_dp * (a + transpose(conjg(a)))
    do i = 1,n
        acopy(i,i) = cmplx(real(acopy(i,i), dp), 0.0_dp, kind=8)
    end do
    lwork = -1
    call zheev('V', 'U', n, acopy, n, eigvals, work_query, lwork, rwork, info)
    if(info /= 0)then
        if(present(flag))then
            flag = info
            return
        endif
        call lapack_stop('HERMITIAN_EIGH workspace query', 'ZHEEV', info)
    endif
    lwork = max(1, int(real(work_query(1), dp)))
    allocate(work(lwork))
    call zheev('V', 'U', n, acopy, n, eigvals, work, lwork, rwork, info)
    if(info /= 0)then
        if(present(flag))then
            flag = info
            return
        endif
        call lapack_stop('HERMITIAN_EIGH', 'ZHEEV', info)
    endif
    eigvecs(1:n,1:n) = acopy
    allocate(vec_tmp(n))
    do i = 1,n / 2
        j = n - i + 1
        eval_tmp   = eigvals(i)
        eigvals(i) = eigvals(j)
        eigvals(j) = eval_tmp
        vec_tmp = eigvecs(1:n,i)
        eigvecs(1:n,i) = eigvecs(1:n,j)
        eigvecs(1:n,j) = vec_tmp
    end do
end subroutine hermitian_eigh_z

subroutine sparse_eigh(matvec, ctx, n, neigs, eigvals, eigvecs, tol, max_basis, info)
    procedure(sparse_matvec_sp_proc) :: matvec
    class(*),          intent(in)  :: ctx
    integer,           intent(in)  :: n, neigs
    real,              intent(out) :: eigvals(neigs)
    real,              intent(out) :: eigvecs(n,neigs)
    real,    optional, intent(in)  :: tol
    integer, optional, intent(in)  :: max_basis
    integer, optional, intent(out) :: info
    character(len=1) :: bmat
    character(len=2) :: which
    logical, allocatable :: select(:)
    real, allocatable :: resid(:), v(:,:), workd(:), workl(:)
    real, allocatable :: dense_mat(:,:), basis_vec(:), prod_vec(:), tmp_vec(:)
    real :: sigma, tol_loc, tmp_val
    integer :: ido, arpack_info, ssaupd_info, ncv, lworkl, i, j, k, nvalid
    integer :: iparam(11), ipntr(11)
    if( present(info) ) info = 0
    eigvals = 0.0_sp
    eigvecs = 0.0_sp
    if( n < 1 ) STOP 'sparse_eigh: empty matrix'
    if( neigs < 1 .or. neigs > n ) STOP 'sparse_eigh: invalid neigs'
    tol_loc = 1.e-5
    if( present(tol) ) tol_loc = tol

    if( neigs == n )then
        allocate(dense_mat(n,n), basis_vec(n), prod_vec(n), source=0.0_sp)
        do i = 1,n
            basis_vec = 0.0_sp
            basis_vec(i) = 1.0_sp
            call matvec(ctx, basis_vec, prod_vec)
            dense_mat(:,i) = prod_vec
        end do
        dense_mat = 0.5_sp * (dense_mat + transpose(dense_mat))
        call eigh(n, dense_mat, neigs, eigvals, eigvecs)
        deallocate(dense_mat, basis_vec, prod_vec)
        return
    endif

    ncv = min(n, max(80, 4 * neigs + 40))
    if( present(max_basis) ) ncv = min(n, max(neigs + 1, max_basis))
    ncv = max(neigs + 1, ncv)
    lworkl = ncv * (ncv + 8)
    allocate(resid(n), v(n,ncv), workd(3*n), workl(lworkl), select(ncv))
    resid = 0.0_sp
    v = 0.0_sp
    workd = 0.0_sp
    workl = 0.0_sp
    select = .false.
    do i = 1,n
        resid(i) = sin(real(37*i,sp)) + 0.25_sp * cos(real(17*i,sp))
    end do

    bmat = 'I'
    which = 'LA'
    ido = 0
    arpack_info = 1
    iparam = 0
    ipntr = 0
    iparam(1) = 1
    iparam(3) = max(300, 10 * ncv)
    iparam(7) = 1

    do
        call ssaupd(ido, bmat, n, which, neigs, tol_loc, resid, ncv, v, n, &
            iparam, ipntr, workd, workl, lworkl, arpack_info)
        if( ido == -1 .or. ido == 1 )then
            call matvec(ctx, workd(ipntr(1):ipntr(1)+n-1), workd(ipntr(2):ipntr(2)+n-1))
        elseif( ido == 99 )then
            exit
        else
            arpack_info = -9999
            exit
        endif
    end do

    ssaupd_info = arpack_info
    if( ssaupd_info < 0 )then
        if( present(info) )then
            info = ssaupd_info
            deallocate(resid, v, workd, workl, select)
            return
        endif
        call arpack_stop('SPARSE_EIGH', 'SSAUPD', ssaupd_info)
    endif

    sigma = 0.0_sp
    arpack_info = 0
    call sseupd(.true., 'A', select, eigvals, eigvecs, n, sigma, bmat, n, which, &
        neigs, tol_loc, resid, ncv, v, n, iparam, ipntr, workd, workl, lworkl, arpack_info)
    if( arpack_info /= 0 )then
        if( present(info) )then
            info = arpack_info
            deallocate(resid, v, workd, workl, select)
            return
        endif
        call arpack_stop('SPARSE_EIGH', 'SSEUPD', arpack_info)
    endif

    nvalid = min(neigs, max(0, iparam(5)))
    if( nvalid == 0 ) nvalid = neigs
    allocate(tmp_vec(n))
    do i = 1,nvalid - 1
        k = i
        do j = i + 1,nvalid
            if( eigvals(j) < eigvals(k) ) k = j
        end do
        if( k /= i )then
            tmp_val = eigvals(i)
            eigvals(i) = eigvals(k)
            eigvals(k) = tmp_val
            tmp_vec = eigvecs(:,i)
            eigvecs(:,i) = eigvecs(:,k)
            eigvecs(:,k) = tmp_vec
        endif
    end do
    deallocate(tmp_vec)

    if( present(info) )then
        info = ssaupd_info
    endif
    deallocate(resid, v, workd, workl, select)
end subroutine sparse_eigh

subroutine test_eigh(n, n_eigs)
    integer, intent(in) :: n, n_eigs
    real(kind=4) :: a(n,n), eigvals(n_eigs), eigvecs(n,n_eigs)
    call random_number(a)
    a = 0.5_sp * (a + transpose(a))
    call eigh(n, a, n_eigs, eigvals, eigvecs)
end subroutine test_eigh

subroutine lapack_stop(caller, routine, info)
    character(len=*), intent(in) :: caller, routine
    integer(kind=4), intent(in) :: info
    write(*,'(a)') ' '
    write(*,'(a)') trim(caller) // ' - Failure!'
    write(*,'(a)') '  LAPACK routine ' // trim(routine) // ' returned a nonzero INFO value.'
    write(*,'(a,i8)') '  INFO = ', info
    stop
end subroutine lapack_stop

subroutine arpack_stop(caller, routine, info)
    character(len=*), intent(in) :: caller, routine
    integer(kind=4), intent(in) :: info
    write(*,'(a)') ' '
    write(*,'(a)') trim(caller) // ' - Failure!'
    write(*,'(a)') '  ARPACK routine ' // trim(routine) // ' returned a nonzero INFO value.'
    write(*,'(a,i8)') '  INFO = ', info
    stop
end subroutine arpack_stop

end module simple_linalg
