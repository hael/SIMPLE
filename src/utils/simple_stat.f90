! statistics utility functions
module simple_stat
!$ use omp_lib
!$ use omp_lib_kinds
use simple_defs
use simple_error, only: simple_exception
use simple_srch_sort_loc
use simple_is_check_assert
implicit none

public :: moment, moment_serial, pearsn, normalize, normalize_sigm, normalize_minmax, stdev
public :: corrs2weights, corr2distweight, analyze_smat, dev_from_dmat, mad, mad_gau, robust_sigma_thres, z_scores
public :: median, median_nocopy, robust_z_scores, robust_normalization, pearsn_serial_8, kstwo
public :: rank_sum_weights, rank_inverse_weights, rank_centroid_weights, rank_exponent_weights
public :: conv2rank_weights, calc_stats, pearsn_serial, norm_corr, norm_corr_8
private
#include "simple_local_flags.inc"

interface moment
    module procedure moment_1
    module procedure moment_2
    module procedure moment_3
    module procedure moment_4
end interface

interface pearsn
    module procedure pearsn_1
    module procedure pearsn_2
    module procedure pearsn_3
end interface

interface normalize
    module procedure normalize_1
    module procedure normalize_2
    module procedure normalize_3
    module procedure normalize_4
end interface

interface normalize_sigm
    module procedure normalize_sigm_1
    module procedure normalize_sigm_2
    module procedure normalize_sigm_3
end interface

real, parameter :: NNET_CONST = exp(1.)-1.

contains

    subroutine moment_1( data, ave, sdev, var, err )
        real,    intent(out) :: ave, sdev, var
        logical, intent(out) :: err
        real,    intent(in)  :: data(:)
        integer  :: n, i
        real     :: nr, dev
        real(dp) :: ep_dp, var_dp
        err  = .false.
        n    = size(data,1)
        nr   = real(n)
        ave  = 0.
        sdev = 0.
        var  = 0.
        if( n < 2 )then
            THROW_WARN('ERROR: n must be at least 2')
            return
        endif
        ! calc average
        ave = sum(real(data,dp))/real(n,dp)
        ! calc sum of devs and sum of devs squared
        ep_dp  = 0.d0
        var_dp = 0.d0
        !$omp parallel do default(shared) private(i,dev) schedule(static)&
        !$omp reduction(+:ep_dp,var_dp) proc_bind(close)
        do i=1,n
            dev    = data(i) - ave
            ep_dp  = ep_dp   + real(dev,dp)
            var_dp = var_dp  + real(dev*dev,dp)
        end do
        !$omp end parallel do
        var  = real((var_dp - ep_dp**2./nr)/(nr-1.)) ! corrected two-pass formula
        sdev = 0.
        if( var > 0. ) sdev = sqrt(var)
        if( abs(var) < TINY )then
            if( var < 0. )then
                err  = .true.
                ave  = 0.
            endif
            var  = 0.
            sdev = 0.
        endif
    end subroutine moment_1

    subroutine moment_2( data, ave, sdev, var, err )
        real,    intent(out) :: ave, sdev, var
        logical, intent(out) :: err
        real,    intent(in)  :: data(:,:)
        integer  :: nx, ny, n, i, j
        real     :: nr, dev
        real(dp) :: ep_dp, var_dp
        err  = .false.
        nx   = size(data,1)
        ny   = size(data,2)
        n    = nx*ny
        nr   = real(n)
        ave  = 0.
        sdev = 0.
        var  = 0.
        if( n < 2 )then
            THROW_WARN('ERROR: n must be at least 2')
            return
        endif
        ! calc average
        ave = sum(real(data,dp))/real(n,dp)
        ! calc sum of devs and sum of devs squared
        ep_dp = 0._dp
        var_dp = 0._dp
        !$omp parallel do default(shared) private(i,j,dev) schedule(static)&
        !$omp reduction(+:ep_dp,var_dp) proc_bind(close) collapse(2)
        do i=1,nx
            do j=1,ny
                dev = data(i,j)-ave
                ep_dp = ep_dp+real(dev,dp)
                var_dp = var_dp+real(dev*dev,dp)
            end do
        end do
        !$omp end parallel do
        var  = real((var_dp-ep_dp**2./nr)/(nr-1.)) ! corrected two-pass formula
        sdev = 0.
        if( var > 0. ) sdev = sqrt(var)
        if( abs(var) < TINY )then
            if( var < 0. )then
                err  = .true.
                ave  = 0.
            endif
            var  = 0.
            sdev = 0.
        endif
    end subroutine moment_2

    subroutine moment_3( data, ave, sdev, var, err )
        real,    intent(out) :: ave, sdev, var
        logical, intent(out) :: err
        real,    intent(in)  :: data(:,:,:)
        integer  :: nx, ny, nz, n, i, j, k
        real     :: nr, dev
        real(dp) :: ep_dp, var_dp
        err  = .false.
        nx   = size(data,1)
        ny   = size(data,2)
        nz   = size(data,3)
        n    = nx*ny*nz
        nr   = real(n)
        ave  = 0.
        sdev = 0.
        var  = 0.
        if( n < 2 )then
            THROW_WARN('ERROR: n must be at least 2')
            return
        endif
        ave = sum(real(data,dp))/real(n,dp)
        ! calc sum of devs and sum of devs squared
        ep_dp = 0._dp
        var_dp = 0._dp
        !$omp parallel do default(shared) private(i,j,k,dev) schedule(static)&
        !$omp reduction(+:ep_dp,var_dp) proc_bind(close) collapse(3)
        do i=1,nx
            do j=1,ny
                do k=1,nz
                    dev = data(i,j,k)-ave
                    ep_dp = ep_dp+real(dev,dp)
                    var_dp = var_dp+real(dev*dev,dp)
                end do
            end do
        end do
        !$omp end parallel do
        var = real((var_dp-ep_dp**2./nr)/(nr-1.)) ! corrected two-pass formula
        sdev = 0.
        if( var > 0. ) sdev = sqrt(var)
        if( abs(var) < TINY )then
            if( var < 0. )then
                err  = .true.
                ave  = 0.
            endif
            var  = 0.
            sdev = 0.
        endif
    end subroutine moment_3

    subroutine moment_4( data, ave, sdev, var, err, mask )
        real,    intent(out) :: ave, sdev, var
        logical, intent(out) :: err
        real,    intent(in)  :: data(:)
        logical, intent(in)  :: mask(:)
        integer  :: n, i, sz
        real     :: nr, dev, abs_var
        real(dp) :: ep_dp, var_dp
        err  = .false.
        sz   = size(data)
        if( sz /= size(mask) ) THROW_HARD('mask does not conform with data; moment_4')
        n    = count(mask)
        nr   = real(n)
        ave  = 0.
        sdev = 0.
        var  = 0.
        if( n < 2 )then
            THROW_WARN('ERROR: n must be at least 2')
            return
        endif
        ! calc average
        ave = sum(real(data,dp), mask=mask)/real(n,dp)
        ! calc sum of devs and sum of devs squared
        ep_dp  = 0._dp
        var_dp = 0._dp
        !$omp parallel do default(shared) private(i,dev) schedule(static)&
        !$omp reduction(+:ep_dp,var_dp) proc_bind(close)
        do i=1,sz
            if( mask(i) )then
                dev = data(i) - ave
                ep_dp  = ep_dp + real(dev,dp)
                var_dp = var_dp + real(dev*dev,dp)
            endif
        end do
        !$omp end parallel do
        var = real((var_dp-ep_dp**2./nr)/(nr-1.)) ! corrected two-pass formula
        sdev = 0.
        if( var > 0. ) sdev = sqrt(var)
        if( abs(var) < TINY )then
            if( var < 0. )then
                err  = .true.
                ave  = 0.
            endif
            var  = 0.
            sdev = 0.
        endif
    end subroutine moment_4

    subroutine moment_serial( data, ave, sdev, var, err, mask )
        real,    intent(out) :: ave, sdev, var
        logical, intent(out) :: err
        real,    intent(in)  :: data(:)
        logical, intent(in)  :: mask(:)
        integer  :: n, i, sz
        real     :: nr, dev, abs_var
        real(dp) :: ep_dp, var_dp
        err  = .false.
        sz   = size(data)
        if( sz /= size(mask) ) THROW_HARD('mask does not conform with data; moment_4')
        n    = count(mask)
        nr   = real(n)
        ave  = 0.
        sdev = 0.
        var  = 0.
        if( n < 2 )then
            THROW_WARN('ERROR: n must be at least 2')
            return
        endif
        ! calc average
        ave = sum(real(data,dp), mask=mask)/real(n,dp)
        ! calc sum of devs and sum of devs squared
        ep_dp  = 0._dp
        var_dp = 0._dp
        do i=1,sz
            if( mask(i) )then
                dev = data(i) - ave
                ep_dp  = ep_dp + real(dev,dp)
                var_dp = var_dp + real(dev*dev,dp)
            endif
        end do
        var = real((var_dp-ep_dp**2./nr)/(nr-1.)) ! corrected two-pass formula
        sdev = 0.
        if( var > 0. ) sdev = sqrt(var)
        if( abs(var) < TINY )then
            if( var < 0. )then
                err  = .true.
                ave  = 0.
            endif
            var  = 0.
            sdev = 0.
        endif
    end subroutine moment_serial

    function pearsn_1( x, y ) result( r )
        real, intent(in) :: x(:),y(:)
        real    :: r,ax,ay,sxx,syy,sxy,xt,yt
        integer :: j, n
        n = size(x)
        if( size(y) /= n ) THROW_HARD('arrays not equal size in pearsn_1')
        ax  = sum(x)/real(n)
        ay  = sum(y)/real(n)
        sxx = 0.
        syy = 0.
        sxy = 0.
        !$omp parallel do default(shared) private(j,xt,yt) &
        !$omp reduction(+:sxx,syy,sxy) schedule(static) proc_bind(close)
        do j=1,n
            xt  = x(j)-ax
            yt  = y(j)-ay
            sxx = sxx+xt**2
            syy = syy+yt**2
            sxy = sxy+xt*yt
        end do
        !$omp end parallel do
        r = max(-1.,min(1.,sxy/sqrt(sxx*syy)))
    end function pearsn_1

    function pearsn_serial( x, y ) result( r )
        real, intent(in) :: x(:),y(:)
        real    :: r,ax,ay,sxx,syy,sxy,xt,yt
        integer :: j, n
        n = size(x)
        if( size(y) /= n ) THROW_HARD('pearsn_serial')
        ax  = sum(x) / real(n)
        ay  = sum(y) / real(n)
        sxx = 0.
        syy = 0.
        sxy = 0.
        do j=1,n
            xt  = x(j) - ax
            yt  = y(j) - ay
            sxx = sxx + xt**2
            syy = syy + yt**2
            sxy = sxy + xt*yt
        end do
        r = sxy / sqrt(sxx * syy)
    end function pearsn_serial

    function pearsn_serial_8( n, x, y ) result( r )
        integer,  intent(in) :: n
        real(dp), intent(in) :: x(:),y(:)
        real(dp) :: ax,ay,sxx,syy,sxy,xt,yt,prod,dn
        real     :: r
        integer  :: j
        dn  = dble(n)
        ax  = sum(x)/dn
        ay  = sum(y)/dn
        sxx = 0.d0
        syy = 0.d0
        sxy = 0.d0
        do j=1,n
            xt  = x(j) - ax
            yt  = y(j) - ay
            sxx = sxx + xt * xt
            syy = syy + yt * yt
            sxy = sxy + xt * yt
        end do
        prod = sxx * syy
        r    = 0.
        if( prod > 0.d0 ) r = real(max(-1.d0,min(1.d0,sxy/sqrt(prod))),kind=4)
    end function pearsn_serial_8

    function pearsn_2( x, y ) result( r )
        real, intent(in) :: x(:,:),y(:,:)
        real    :: r,ax,ay,sxx,syy,sxy,xt,yt
        integer :: i, j, nx, ny
        nx = size(x,1)
        ny = size(x,2)
        if( size(y,1) /= nx .or. size(y,2) /= ny ) THROW_HARD('arrays not equal size in pearsn_2')
        ax  = sum(x)/real(nx*ny)
        ay  = sum(y)/real(nx*ny)
        sxx = 0.
        syy = 0.
        sxy = 0.
        !$omp parallel do default(shared) private(i,j,xt,yt) collapse(2)&
        !$omp reduction(+:sxx,syy,sxy) schedule(static) proc_bind(close)
        do i=1,nx
            do j=1,ny
                xt  = x(i,j)-ax
                yt  = y(i,j)-ay
                sxx = sxx+xt**2
                syy = syy+yt**2
                sxy = sxy+xt*yt
            end do
        end do
        !$omp end parallel do
        r = max(-1.,min(1.,sxy/sqrt(sxx*syy)))
    end function pearsn_2

    function pearsn_3( x, y ) result( r )
        real, intent(in) :: x(:,:,:),y(:,:,:)
        real    :: r,ax,ay,sxx,syy,sxy,xt,yt
        integer :: i, j, k, nx, ny, nz
        nx = size(x,1)
        ny = size(x,2)
        nz = size(x,3)
        if( size(y,1) /= nx .or. size(y,2) /= ny .or. size(y,3) /= nz )&
        THROW_HARD('arrays not equal size, in pearsn_3')
        ax  = sum(x)/real(nx*ny*nz)
        ay  = sum(y)/real(nx*ny*nz)
        sxx = 0.
        syy = 0.
        sxy = 0.
        !$omp parallel do default(shared) private(i,j,k,xt,yt) collapse(2)&
        !$omp reduction(+:sxx,syy,sxy) schedule(static) proc_bind(close)
        do i=1,nx
            do j=1,ny
                do k=1,nz
                    xt  = x(i,j,k)-ax
                    yt  = y(i,j,k)-ay
                    sxx = sxx+xt**2
                    syy = syy+yt**2
                    sxy = sxy+xt*yt
                end do
            end do
        end do
        !$omp end parallel do
        r = max(-1.,min(1.,sxy/sqrt(sxx*syy)))
    end function pearsn_3

    subroutine normalize_1( arr, err )
        real, intent(inout)  :: arr(:)
        logical, intent(out) :: err
        real :: ave, sdev, var
        call moment_1( arr, ave, sdev, var, err )
        if( err ) return
        arr = (arr-ave)/sdev
    end subroutine normalize_1

    subroutine normalize_2( arr, err )
        real,    intent(inout)  :: arr(:,:)
        logical, intent(out) :: err
        real :: ave, sdev, var
        call moment_2( arr, ave, sdev, var, err )
        if( err ) return
        arr = (arr-ave)/sdev
    end subroutine normalize_2

    subroutine normalize_3( arr, err )
        real, intent(inout)  :: arr(:,:,:)
        logical, intent(out) :: err
        real :: ave, sdev, var
        call moment_3( arr, ave, sdev, var, err )
        if( err ) return
        arr = (arr-ave)/sdev
    end subroutine normalize_3

    subroutine normalize_4( arr, err, mask )
        real, intent(inout)  :: arr(:)
        logical, intent(out) :: err
        logical, intent(in)  :: mask(:)
        real :: ave, sdev, var
        call moment_4( arr, ave, sdev, var, err, mask )
        if( err ) return
        where( mask ) arr = (arr-ave)/sdev
    end subroutine normalize_4

    subroutine normalize_sigm_1( arr )
        real, intent(inout) :: arr(:)
        real                :: smin, smax, delta
        if( size(arr) == 1 )then
            arr(1) = max(0., min(arr(1), 1.))
            return
        endif
        ! find minmax
        smin  = minval(arr)
        smax  = maxval(arr)
        delta = smax-smin
        if( delta > TINY )then
            ! create [0,1]-normalized vector
            !$omp parallel workshare default(shared)
            arr = (exp((arr-smin)/delta)-1.)/NNET_CONST
            !$omp end parallel workshare
        else
            THROW_WARN('normalize_sigm_1, division with zero, no normalisation done')
        endif
    end subroutine normalize_sigm_1

    subroutine normalize_sigm_2( arr )
        real, intent(inout) :: arr(:,:)
        real                :: smin, smax, delta
        ! find minmax
        smin  = minval(arr)
        smax  = maxval(arr)
        delta = smax-smin
        if( delta > TINY )then
            ! create [0,1]-normalized vector
            !$omp parallel workshare default(shared)
            arr = (exp((arr-smin)/delta)-1.)/NNET_CONST
            !$omp end parallel workshare
        else
            THROW_WARN('normalize_sigm_2, division with zero, no normalisation done')
        endif
    end subroutine normalize_sigm_2

    subroutine normalize_sigm_3( arr )
        real, intent(inout) :: arr(:,:,:)
        real                :: smin, smax, delta
        ! find minmax
        smin  = minval(arr)
        smax  = maxval(arr)
        delta = smax-smin
        if( delta > TINY )then
            ! create [0,1]-normalized vector
            !$omp parallel workshare default(shared)
            arr = (exp((arr-smin)/delta)-1.)/NNET_CONST
            !$omp end parallel workshare
        else
            THROW_WARN('normalize_sigm_3, division with zero, no normalisation done')
        endif
    end subroutine normalize_sigm_3

    subroutine normalize_minmax( arr )
        real, intent(inout) :: arr(:)
        real                :: smin, smax, delta
        smin  = minval(arr)
        smax  = maxval(arr)
        delta = smax - smin
        arr   = (arr - smin)/delta
    end subroutine normalize_minmax

    function stdev (a)
        real, intent(in) :: a(:)
        real    :: stdev, avg, SumSQR, var, n
        n= real(size(a))
        if(n < 2) return
        avg = sum(a)/n
        SumSQR = sum(sqrt(a))/n
        var = (SumSQR - avg*avg)/(n-1)
        stdev = sqrt(var)
    end function stdev

    !>  \brief  is for calculating variable statistics
    subroutine calc_stats( arr, statvars, mask )
        real,               intent(in)  :: arr(:)
        type(stats_struct), intent(out) :: statvars
        logical, optional,  intent(in)  :: mask(:)
        real, allocatable :: tmparr(:)
        real    :: var
        logical :: err
        if( present(mask) )then
            call moment(arr, statvars%avg, statvars%sdev, var, err, mask)
            tmparr        = pack(arr, mask=mask)
            statvars%med  = median_nocopy(tmparr)
            statvars%minv = minval(arr, mask)
            statvars%maxv = maxval(arr, mask)
            deallocate(tmparr)
        else
            call moment(arr, statvars%avg, statvars%sdev, var, err)
            statvars%med  = median(arr)
            statvars%minv = minval(arr)
            statvars%maxv = maxval(arr)
        endif
    end subroutine calc_stats

    elemental pure function norm_corr( xprod, sqprod1, sqprod2 ) result( cc_norm )
        real, intent(in) :: xprod, sqprod1, sqprod2
        real :: eps, sqrt_denom, cc_norm
        sqrt_denom = sqrt(sqprod1 * sqprod2)
        eps = epsilon(xprod)
        if( xprod < eps .and. sqprod1 < eps .and. sqprod2 < eps )then
            cc_norm = 1.
        elseif( sqrt_denom < eps )then
            cc_norm = 0.
        else
            cc_norm = xprod / sqrt_denom
        endif
    end function norm_corr

    elemental pure function norm_corr_8( xprod, sqprod1, sqprod2 ) result( cc_norm )
        real(dp), intent(in) :: xprod, sqprod1, sqprod2
        real(dp) :: eps, sqrt_denom, cc_norm
        sqrt_denom = sqrt(sqprod1 * sqprod2)
        eps = epsilon(xprod)
        if( xprod < eps .and. sqprod1 < eps .and. sqprod2 < eps )then
            cc_norm = 1.
        elseif( sqrt_denom < eps )then
            cc_norm = 0.
        else
            cc_norm = xprod / sqrt_denom
        endif
    end function norm_corr_8

    function corrs2weights( corrs, crit, p, norm_sigm ) result( weights )
        real,                           intent(in) :: corrs(:) !< correlation input
        integer(kind=kind(ENUM_WCRIT)), intent(in) :: crit
        real,                 optional, intent(in) :: p
        logical,              optional, intent(in) :: norm_sigm
        real, allocatable :: weights(:), corrs_copy(:)
        real, parameter   :: THRESHOLD=1.5
        real    :: maxminratio, corrmax, corrmin, minw
        integer :: ncorrs
        logical :: nnorm_sigm
        nnorm_sigm = .true.
        if( present(norm_sigm) ) nnorm_sigm = norm_sigm
        ncorrs = size(corrs)
        allocate(weights(ncorrs), source=0.)
        allocate(corrs_copy(ncorrs), source=corrs)
        corrmax = maxval(corrs_copy)
        if( corrmax <= 0. )then
            ! weighting does not make sense, put them all to 0.
            return
        endif
        ! remove negatives to prevent corrs around zero to recieve any weight power
        where( corrs_copy <= TINY ) corrs_copy = 0.
        ! correlation-based weights
        select case(crit)
            case(CORRW_CRIT)
                if( nnorm_sigm )then
                    corrmin     = minval(corrs_copy, mask=corrs_copy > TINY)
                    maxminratio = corrmax / corrmin
                    if( maxminratio >= THRESHOLD )then
                       ! min/max normalise the correlations
                       call normalize_sigm(corrs_copy)
                    endif
                endif
                ! calculate the exponential of the negative distances
                where( corrs_copy > TINY )
                    weights = exp(-(1. - corrs_copy))
                else where
                    weights = 0.
                end where
                ! normalize weights
                weights = weights / sum(weights)
            case(CORRW_ZSCORE_CRIT)
                weights = z_scores(corrs_copy)
                minw    = minval(weights)
                weights = weights + abs(minw)
                where( corrs_copy <= 0. ) weights = 0.
            case(RANK_SUM_CRIT,RANK_CEN_CRIT,RANK_EXP_CRIT,RANK_INV_CRIT)
                ! convert to rank-based weights
                weights = corrs_copy
                call conv2rank_weights(ncorrs, weights, crit, p)
            case DEFAULT
                THROW_HARD('Unsupported conversion criterion; corrs2weights')
        end select
    end function corrs2weights

    elemental real(dp) function corr2distweight( cc, npix, tau )
        real,    intent(in) :: cc   ! correlation
        integer, intent(in) :: npix ! total number of carteasian pixels
        real,    intent(in) :: tau  !< scale factor
        corr2distweight = 2.d0 * (1.d0-real(cc,dp)) * real(npix,dp) / real(tau,dp)
    end function corr2distweight

    !>   Kolmogorov-Smirnov test to deduce equivalence or non-equivalence
    !>  between two distributions.
    !          The routine returns the K-S statistic d, and the significance
    !          level prob for the null hypothesis that the data sets are drawn
    !          from the same distribution. Small values for prob show that the
    !          cumulative distribution function of data1 is significantly
    !          different from that of data2. The input arrays are modified
    !          (sorted)
    subroutine kstwo( data1, n1, data2, n2, d, prob )
        integer, intent(in)    :: n1, n2
        real,    intent(inout) :: data1(n1), data2(n2), d, prob !< significance
        integer                :: j1, j2
        real                   :: d1, d2, dt, en1, en2, en, fn1, fn2
        call hpsort(data1)
        call hpsort(data2)
        en1 = n1
        en2 = n2
        j1  = 1
        j2  = 1
        fn1 = 0.
        fn2 = 0.
        d   = 0.
        do while(j1.le.n1.and.j2.le.n2)
            d1 = data1(j1)
            d2 = data2(j2)
            if(d1.le.d2)then
                fn1 = j1/en1
                j1  = j1+1
            endif
            if(d2.le.d1)then
                fn2 = j2/en2
                j2  = j2+1
            endif
            dt = abs(fn2-fn1)
            if(dt.gt.d) d = dt
        end do
        en   = sqrt(en1*en2/(en1+en2))
        prob = probks((en+0.12+0.11/en)*d) ! significance

        contains

            !< Calculate K-S significance test
            function probks( alam ) result( p )
                real, intent(in) :: alam
                real :: p, a2, fac, term, termbf
                real, parameter :: EPS1=0.001, EPS2=1.e-8
                integer :: j
                a2  = -2.*alam**2
                fac = 2.
                p = 0.
                termbf = 0. ! prev term in sum
                do j=1,100
                    term = fac*exp(a2*j**2)
                    p = p+term
                    if(abs(term).le.EPS1*termbf.or.abs(term).le.EPS2*p) return
                    fac = -fac ! alternate signs in sum
                    termbf = abs(term)
                end do
                p=1. ! fail to converge
            end function probks

    end subroutine kstwo

    subroutine analyze_smat( s, symmetrize, smin, smax )
        real, intent(inout) :: s(:,:)          !< similarity matrix
        logical, intent(in) :: symmetrize      !< force diag symmetry
        real, intent(out)   :: smin, smax
        integer             :: i, j, n, npairs
        if( size(s,1) .ne. size(s,2) )then
            THROW_HARD('not a similarity matrix')
        endif
        n      = size(s,1)
        npairs = (n*(n-1))/2
        smin   = s(1,2)
        smax   = smin
        do i=1,n-1
            do j=i+1,n
                if( symmetrize ) s(j,i)  = s(i,j)
                if( s(i,j) < smin ) smin = s(i,j)
                if( s(i,j) > smax ) smax = s(i,j)
            end do
        end do
    end subroutine analyze_smat

    !>  measure of cluster spread (for use in one-class clustering)
    !!  (1) estimate median of cluster (member most similar to all others)
    !!  (2) calculate the median of the similarities between the median and all others
    !!  suggested exclusion based on 2 * ddev (sigma) criterion
    subroutine dev_from_smat( smat, i_median, sdev )
        real,    intent(in)  :: smat(:,:)
        integer, intent(out) :: i_median
        real,    intent(out) :: sdev
        real, allocatable :: sims(:)
        integer :: loc(1), i, j, n
        n = size(smat,1)
        if( n /= size(smat,2) ) THROW_HARD('symmetric similarity matrix assumed; stat :: dev_from_smat')
        allocate(sims(n))
        do i=1,n
            sims(i) = 0.0
            do j=1,n
                if( i /= j )then
                    sims(i) = sims(i) + smat(i,j)
                endif
            end do
        end do
        loc      = maxloc(sims)
        i_median = loc(1)
        sdev     = median(smat(i_median,:))
    end subroutine dev_from_smat

    !>  measure of cluster spread (for use in one-class clustering)
    !!  (1) estimate median of cluster (member most similar to all others)
    !!  (2) calculate the median of the distances between the median and all others
    !!  suggested exclusion based on 2 * ddev (sigma) criterion
    subroutine dev_from_dmat( dmat, i_median, ddev )
        real,    intent(in)  :: dmat(:,:)
        integer, intent(out) :: i_median
        real,    intent(out) :: ddev
        real, allocatable :: dists(:)
        integer :: loc(1), i, j, n
        n = size(dmat,1)
        if( n /= size(dmat,2) ) THROW_HARD('symmetric distance matrix assumed; dev_from_dmat')
        allocate(dists(n))
        do i=1,n
            dists(i) = 0.0
            do j=1,n
                if( i /= j )then
                    dists(i) = dists(i) + dmat(i,j)
                endif
            end do
        end do
        loc      = minloc(dists)
        i_median = loc(1)
        ddev     = median(dmat(i_median,:))
    end subroutine dev_from_dmat

    ! ROBUST STATISTICS

    !>   for calculating the median
    function median( arr ) result( val )
        real, intent(in)  :: arr(:)
        real              :: copy(size(arr))
        real    :: val, val1, val2
        integer :: n, pos1, pos2
        n = size(arr)
        if( is_even(n) )then
            pos1 = n/2
            pos2 = pos1+1
        else
            pos1 = nint(real(n)/2.)
            pos2 = pos1
        endif
        copy = arr
        if( pos1 == pos2 )then
            val  = selec(pos1,n,copy)
        else
            val1 = selec(pos1,n,copy)
            val2 = selec(pos2,n,copy)
            val  = (val1+val2)/2.
        endif
    end function median

    !>   for calculating the median
    function median_nocopy( arr ) result( val )
        real, intent(inout) :: arr(:)
        real    :: val, val1, val2
        integer :: n, pos1, pos2
        n = size(arr)
        if( is_even(n) )then
            pos1 = n/2
            pos2 = pos1+1
        else
            pos1 = nint(real(n)/2.)
            pos2 = pos1
        endif
        if( pos1 == pos2 )then
            val  = selec(pos1,n,arr)
        else
            val1 = selec(pos1,n,arr)
            val2 = selec(pos2,n,arr)
            val  = (val1+val2)/2.
        endif
    end function median_nocopy

    ! median absolute deviation
    ! calculated as the median of absolute deviations of the data points
    real function mad( x, med )
        real, intent(in) :: x(:) ! data points
        real, intent(in) :: med  ! median of data points
        real, allocatable :: absdevs(:)
        allocate(absdevs(size(x)), source=abs(x - med))
        mad = median_nocopy(absdevs)
    end function mad

    ! median absolute deviation, assuming underlying Gaussian distribution
    real function mad_gau( x, med )
        real, intent(in) :: x(:) ! data points
        real, intent(in) :: med  ! median of data points
        mad_gau = 1.4826 * mad(x, med)
    end function mad_gau

    function robust_sigma_thres( x, nsig ) result( t )
        real, intent(in)    :: x(:) ! data points
        real, intent(in)    :: nsig ! # sigmas, can be +/-
        real, allocatable   :: absdevs(:)
        real :: med, mad, t
        med = median(x)
        allocate(absdevs(size(x)), source=abs(x - med))
        mad = median_nocopy(absdevs)
        t = med + nsig * mad
    end function robust_sigma_thres

    ! the Z-score calculates the number of standard deviations a data point is away from the mean
    function z_scores( x, mask ) result( zscores )
        real,              intent(in) :: x(:) ! data points
        logical, optional, intent(in) :: mask(:)
        real, allocatable :: zscores(:)
        logical :: err
        allocate(zscores(size(x)), source=x)
        if( present(mask) )then
            call normalize(zscores, err, mask)
            if( err )then
                where( mask ) zscores = 1.0
            endif
        else
            call normalize(zscores, err)
            if( err ) zscores = 1.0
        endif
    end function z_scores

    ! the Z-score calculates the number of standard deviations a data point is away from the mean
    ! the robust Z-score does so in the presence of outliers
    function robust_z_scores( x ) result( z_scores )
        real, intent(in)  :: x(:) ! data points
        real, allocatable :: z_scores(:)
        real :: med
        med = median(x)
        allocate(z_scores(size(x)), source=(x - med) / mad_gau(x, med))
    end function robust_z_scores

    subroutine robust_normalization( x )
        real, intent(inout) :: x(:) ! data points
        real :: med, dev
        med = median(x)
        dev = mad_gau(x, med)
        x   = (x - med) / dev
    end subroutine robust_normalization

    subroutine expweights_from_scores( n, scores, weights )
        integer, intent(in)    :: n
        real,    intent(in)    :: scores(n)
        real,    intent(inout) :: weights(n)
        real(dp) :: best_score, sumw,  diffs(n), sigma, cnt
        integer  :: i, loc(1)
        loc        = maxloc(scores)
        best_score = real(scores(loc(1)),kind=dp)
        sigma      = 0.d0
        cnt        = 0.d0
        do i = 1,n
            ! the argument to the exponential function is always negative
            diffs(i) = real(scores(i),kind=dp) - best_score
            sigma    = sigma + diffs(i) * diffs(i)
            cnt      = cnt   + 1.d0
        end do
        sigma = sqrt(sigma / (cnt - 1.d0))
        sumw  = 0.d0
        do i = 1,n
            ! for the best score the exponent will be 1 and < 1 for all others
            if( i == loc(1) )then
                sumw  = sumw + 1.d0
            else if( diffs(i) < 0.d0 )then
                sumw  = sumw + exp(diffs(i) / sigma)
            endif
        end do
        do i = 1,n
            weights(i) = max(0.,min(1.,real(diffs(i) / sumw)))
        end do
    end subroutine expweights_from_scores

    ! RANK STUFF

    subroutine conv2rank_weights( n, weights, crit, p )
        integer,                        intent(in)    :: n
        real,                           intent(inout) :: weights(n)
        integer(kind=kind(ENUM_WCRIT)), intent(in)    :: crit
        real, optional,                 intent(in)    :: p
        real,    allocatable :: weights_tmp(:)
        integer, allocatable :: ranks(:)
        logical :: mask(n)
        integer :: n_nonzero, i, cnt
        mask = weights > TINY
        n_nonzero = count(mask)
        if( n_nonzero == 1 )then
            where( mask )
                weights = 1.
            elsewhere
                weights = 0.
            endwhere
            return
        else if( n_nonzero < 1 )then
            weights = 0.
            return
        endif
        ! produce ranking
        allocate(ranks(n_nonzero),weights_tmp(n_nonzero))
        cnt = 0
        do i = 1,n
            if( mask(i) )then
                cnt              = cnt + 1
                ranks(cnt)       = i
                weights_tmp(cnt) = weights(i)
            endif
        enddo
        call hpsort(weights_tmp, ranks) ! largest last
        call reverse(ranks)             ! largest first
        call reverse(weights_tmp)
        ! calculate weights from ranks
        select case(crit)
            case(RANK_SUM_CRIT)
                call rank_sum_weights(n_nonzero, weights_tmp)
            case(RANK_CEN_CRIT)
                call rank_centroid_weights(n_nonzero, weights_tmp)
            case(RANK_EXP_CRIT)
                if( present(p) )then
                    call rank_exponent_weights(n_nonzero, p, weights_tmp)
                else
                    THROW_HARD('need exponent (p) for rank exponent weight calculation; conv2rank_weights')
                endif
            case(RANK_INV_CRIT)
                call rank_inverse_weights(n_nonzero, weights_tmp)
            case DEFAULT
                THROW_HARD('unsupported rank ordering criteria weighting method; conv2rank_weights')
        end select
        weights = 0.
        do i=1,n_nonzero
            weights(ranks(i)) = weights_tmp(i)
        enddo
        ! re-normalize
        weights = weights / sum(weights_tmp)
        deallocate(weights_tmp, ranks)
    end subroutine conv2rank_weights

    subroutine rank_sum_weights( n, weights )
        integer, intent(in)  :: n
        real,    intent(out) :: weights(n)
        real    :: denom, rn, rn_plus_1
        integer :: irank
        rn        = real(n)
        rn_plus_1 = rn + 1.0
        denom     = rn * rn_plus_1
        do irank=1,n
            weights(irank) = (2.0 * (rn_plus_1 - real(irank))) / denom
        end do
    end subroutine rank_sum_weights

    ! ROC weights
    subroutine rank_centroid_weights( n, weights )
        integer, intent(in)  :: n
        real,    intent(out) :: weights(n)
        real    :: inv_ranks(n), rn
        integer :: irank
        rn = real(n)
        do irank=1,n
            inv_ranks(irank) = 1.0 / real(irank)
        end do
        do irank=1,n
            weights(irank) = (1.0 / rn) * sum(inv_ranks(irank:))
        end do
    end subroutine rank_centroid_weights

    subroutine rank_exponent_weights( n, p, weights )
        integer, intent(in)  :: n
        real,    intent(in)  :: p
        real,    intent(out) :: weights(n)
        real    :: rn
        integer :: irank
        rn = real(n)
        do irank=1,n
            weights(irank) = (rn - real(irank) + 1.0)**p
        end do
        weights = weights / sum(weights)
    end subroutine rank_exponent_weights

    subroutine rank_inverse_weights( n, weights )
        integer, intent(in)  :: n
        real,    intent(out) :: weights(n)
        integer :: irank
        do irank=1,n
            weights(irank) = 1.0 / real(irank)
        end do
        weights = weights / sum(weights)
    end subroutine rank_inverse_weights

end module simple_stat
