! statistics utility functions
module simple_stat
use simple_defs   ! singleton
use simple_error, only: allocchk, simple_exception
use simple_math,  only: hpsort, median, median_nocopy
implicit none

public :: moment, pearsn, normalize, normalize_sigm, normalize_minmax
public :: corrs2weights, analyze_smat, dev_from_dmat, mad, mad_gau, z_scores
public :: robust_z_scores, robust_normalization, pearsn_serial_8, kstwo
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

contains

    ! MOMENTS/NORMALIZATION

    !>    given a 1D real array of data, this routine returns its mean: _ave_,
    !!          standard deviation: _sdev_, and variance: _var_
    !! \param data input data
    !! \param ave,sdev,var  Geometric average,   standard deviation  and variance
    subroutine moment_1( data, ave, sdev, var, err )
        !$ use omp_lib
        !$ use omp_lib_kinds
        real,    intent(out) :: ave, sdev, var
        logical, intent(out) :: err            !< error status
        real,    intent(in)  :: data(:)
        integer              :: n, i
        real                 :: ep, nr, dev
        err = .false.
        n   = size(data,1)
        nr  = real(n)
        if( n <= 1 ) then
            write(logfhandle,*) 'ERROR: n must be at least 2'
            write(logfhandle,*) 'In: moment_1, module: simple_stat.f90'
            stop
        endif
        ! calc average
        ave = sum(data)/nr
        ! calc sum of devs and sum of devs squared
        ep = 0.
        var = 0.
        !$omp parallel do default(shared) private(i,dev) schedule(static)&
        !$omp reduction(+:ep,var) proc_bind(close)
        do i=1,n
            dev = data(i)-ave
            ep = ep+dev
            var = var+dev*dev
        end do
        !$omp end parallel do
        var = (var-ep**2./nr)/(nr-1.) ! corrected two-pass formula
        sdev = 0.
        if( var > 0. ) sdev = sqrt(var)
        if( abs(var) < TINY )then
            err  = .true.
            ave  = 0.
            sdev = 0.
            var  = 0.
        endif
    end subroutine moment_1

    !>    given a 2D real array of data, this routine returns its mean: _ave_,
    !!          standard deviation: _sdev_, and variance: _var_
    !! \param data input data
    !! \param ave,sdev,var  Geometric mean, standard deviation and variance
    subroutine moment_2( data, ave, sdev, var, err )
        !$ use omp_lib
        !$ use omp_lib_kinds
        real, intent(out)    :: ave, sdev, var
        logical, intent(out) :: err              !< error status
        real, intent(in)     :: data(:,:)        !< input data
        integer              :: nx, ny, n, i, j
        real                 :: ep, nr, dev
        err = .false.
        nx = size(data,1)
        ny = size(data,2)
        n  = nx*ny
        nr = real(n)
        if( n <= 1 ) then
            write(logfhandle,*) 'ERROR: n must be at least 2'
            write(logfhandle,*) 'In: moment_2, module: simple_stat.f90'
            stop
        endif
        ! calc average
        ave = sum(data)/nr
        ! calc sum of devs and sum of devs squared
        ep = 0.
        var = 0.
        !$omp parallel do default(shared) private(i,j,dev) schedule(static)&
        !$omp reduction(+:ep,var) proc_bind(close) collapse(2)
        do i=1,nx
            do j=1,ny
                dev = data(i,j)-ave
                ep = ep+dev
                var = var+dev*dev
            end do
        end do
        !$omp end parallel do
        var = (var-ep**2./nr)/(nr-1.) ! corrected two-pass formula
        sdev = 0.
        if( var > 0. ) sdev = sqrt(var)
        if( abs(var) < TINY )then
            err  = .true.
            ave  = 0.
            sdev = 0.
            var  = 0.
        endif
    end subroutine moment_2

    !>    given a 3D real array of data, this routine returns its mean: _ave_,
    !!          standard deviation: _sdev_, and variance: _var_
    !! \param data input data
    !! \param ave,sdev,var  Geometric average,   standard deviation  and variance
    subroutine moment_3( data, ave, sdev, var, err )
        !$ use omp_lib
        !$ use omp_lib_kinds
        real, intent(out)    :: ave, sdev, var
        logical, intent(out) :: err               !< error status
        real, intent(in)     :: data(:,:,:)       !< input data
        integer              :: nx, ny, nz, n, i, j, k
        real                 :: ep, nr, dev
        err = .false.
        nx = size(data,1)
        ny = size(data,2)
        nz = size(data,3)
        n  = nx*ny*nz
        nr = real(n)
        if( n <= 1 ) then
            write(logfhandle,*) 'ERROR: n must be at least 2'
            write(logfhandle,*) 'In: moment_3, module: simple_stat.f90'
            stop
        endif
        ave = sum(data)/nr
        ! calc sum of devs and sum of devs squared
        ep = 0.
        var = 0.
        !$omp parallel do default(shared) private(i,j,k,dev) schedule(static)&
        !$omp reduction(+:ep,var) proc_bind(close) collapse(3)
        do i=1,nx
            do j=1,ny
                do k=1,nz
                    dev = data(i,j,k)-ave
                    ep = ep+dev
                    var = var+dev*dev
                end do
            end do
        end do
        !$omp end parallel do
        var = (var-ep**2./nr)/(nr-1.) ! corrected two-pass formula
        sdev = 0.
        if( var > 0. ) sdev = sqrt(var)
        if( abs(var) < TINY )then
            err  = .true.
            ave  = 0.
            sdev = 0.
            var  = 0.
        endif
    end subroutine moment_3

    !>    given a 1D real array of data, this routine returns its mean: _ave_,
    !!          standard deviation: _sdev_, and variance: _var_
    !! \param data input data
    !! \param ave,sdev,var  Geometric average,   standard deviation  and variance
    subroutine moment_4( data, ave, sdev, var, err, mask )
        !$ use omp_lib
        !$ use omp_lib_kinds
        real,    intent(out) :: ave, sdev, var
        logical, intent(out) :: err !< error status
        real,    intent(in)  :: data(:)
        logical, intent(in)  :: mask(:)
        integer :: n, i, sz
        real    :: ep, nr, dev
        err = .false.
        sz  = size(data)
        if( sz /= size(mask) ) THROW_HARD('mask does not conform with data; moment_4')
        n   = count(mask)
        nr  = real(n)
        if( n <= 1 ) THROW_HARD('n must be at least 2; moment_4')
        ! calc average
        ave = sum(data, mask=mask)/nr
        ! calc sum of devs and sum of devs squared
        ep  = 0.
        var = 0.
        !$omp parallel do default(shared) private(i,dev) schedule(static)&
        !$omp reduction(+:ep,var) proc_bind(close)
        do i=1,sz
            if( mask(i) )then
                dev = data(i) - ave
                ep  = ep + dev
                var = var + dev*dev
            endif
        end do
        !$omp end parallel do
        var = (var-ep**2./nr)/(nr-1.) ! corrected two-pass formula
        sdev = 0.
        if( var > 0. ) sdev = sqrt(var)
        if( abs(var) < TINY )then
            err  = .true.
            ave  = 0.
            sdev = 0.
            var  = 0.
        endif
    end subroutine moment_4

    !>    is for statistical normalization of an array
    subroutine normalize_1( arr, err )
        real, intent(inout)  :: arr(:)           !< input data
        logical, intent(out) :: err              !< error status
        real :: ave, sdev, var
        call moment_1( arr, ave, sdev, var, err )
        if( err ) return
        arr = (arr-ave)/sdev ! array op
    end subroutine normalize_1

    !>    is for statistical normalization of a 2D matrix
    subroutine normalize_2( arr, err )
        real,    intent(inout)  :: arr(:,:)      !< input data
        logical, intent(out) :: err              !< error status
        real :: ave, sdev, var
        call moment_2( arr, ave, sdev, var, err )
        if( err ) return
        arr = (arr-ave)/sdev ! array op
    end subroutine normalize_2

    !>    is for statistical normalization of a 3D matrix
    subroutine normalize_3( arr, err )
        real, intent(inout)  :: arr(:,:,:)       !< input data
        logical, intent(out) :: err              !< error status
        real :: ave, sdev, var
        call moment_3( arr, ave, sdev, var, err )
        if( err ) return
        arr = (arr-ave)/sdev ! array op
    end subroutine normalize_3

    !>    is for statistical normalization of an array
    subroutine normalize_4( arr, err, mask )
        real, intent(inout)  :: arr(:)           !< input data
        logical, intent(out) :: err              !< error status
        logical, intent(in)  :: mask(:)          !< logical mask
        real :: ave, sdev, var
        call moment_4( arr, ave, sdev, var, err, mask )
        if( err ) return
        where( mask ) arr = (arr-ave)/sdev ! array op
    end subroutine normalize_4

    !>    is for sigmoid normalisation [0,1]
    subroutine normalize_sigm_1( arr )
        real, intent(inout) :: arr(:)
        real                :: smin, smax, delta
        real, parameter     :: NNET_CONST = exp(1.)-1.
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

    !>    is for sigmoid normalisation [0,1]
    subroutine normalize_sigm_2( arr )
        real, intent(inout) :: arr(:,:) !< input data
        real                :: smin, smax, delta
        real, parameter     :: NNET_CONST = exp(1.)-1.
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

    !>    is for sigmoid normalisation [0,1]
    subroutine normalize_sigm_3( arr )
        real, intent(inout) :: arr(:,:,:) !< input data
        real                :: smin, smax, delta
        real, parameter     :: NNET_CONST = exp(1.)-1.
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
        if( size(arr) == 1 )then
            arr(1) = max(0., min(arr(1), 1.))
            return
        endif
         ! find minmax
        smin  = minval(arr)
        smax  = maxval(arr)
        delta = smax - smin
        if( delta > TINY )then
            !$omp parallel workshare default(shared)
            arr = (arr - smin)/delta
            !$omp end parallel workshare
        else
            THROW_WARN('normalize_minmax, division with zero, no normalisation done')
        endif
    end subroutine normalize_minmax

    ! CORRELATION

    !>    calculates Pearson's correlation coefficient
    !! \param x input reference array
    !! \param y input test array
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

    !>    calculates Pearson's correlation coefficient
    !! \param x input reference array
    !! \param y input test array
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

    !>    calculates Pearson's correlation coefficient
    !! \param x input reference array
    !! \param y input test array
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

    !>    calculates Pearson's correlation coefficient
    !! \param x input reference array
    !! \param y input test array
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

    function corrs2weights( corrs ) result( weights )
        real, intent(in)  :: corrs(:) !< correlation input
        real, allocatable :: weights(:), corrs_copy(:)
        real, parameter   :: THRESHOLD=1.5
        real    :: maxminratio, corrmax, corrmin
        integer :: ncorrs
        ncorrs = size(corrs)
        allocate(weights(ncorrs), stat=alloc_stat)
        if(alloc_stat /= 0) call allocchk("In: corrs2weights; simple_stat 1" , alloc_stat)
        weights = 0.
        allocate(corrs_copy(ncorrs), source=corrs,stat=alloc_stat)
        if(alloc_stat /= 0) call allocchk("In: corrs2weights; simple_stat 2" , alloc_stat)
        corrmax = maxval(corrs_copy)
        if( corrmax <= 0. )then
            ! weighting does not make sense, put them all to 1/ncorrs
            weights = 1. / real(ncorrs)
            return
        endif
        ! remove negatives to prevent corrs around zero to recieve any weight power
        where( corrs_copy <= 0. ) corrs_copy = 0.
        corrmin     = minval(corrs_copy, mask=corrs_copy > TINY)
        maxminratio = corrmax / corrmin
        if( maxminratio >= THRESHOLD )then
           ! min/max normalise the correlations
           call normalize_sigm(corrs_copy)
        endif
        ! calculate the exponential of the negative distances
        where( corrs_copy > TINY )
            weights = exp(-(1. - corrs_copy))
        else where
            weights = 0.
        end where
        ! normalize weights
        weights = weights / sum(weights)
    end function corrs2weights

    ! integer STUFF

    !>   Kolmogorov-Smirnov test to deduce equivalence or non-equivalence
    !>  between two distributions.
    !          The routine returns the K-S statistic d, and the significance
    !          level prob for the null hypothesis that the data sets are drawn
    !          from the same distribution. Small values for prob show that the
    !          cumulative distribution function of data1 is significantly
    !          different from that of data2. The input arrays are modified
    !          (sorted)
    subroutine kstwo( data1, n1, data2, n2, d, prob )
        integer, intent(in) :: n1, n2
        real, intent(inout) :: data1(n1), data2(n2), d, prob  !< significance
        integer             :: j1, j2
        real                :: d1, d2, dt, en1, en2, en, fn1, fn2
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
        en = sqrt(en1*en2/(en1+en2))
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

    !>    4 statistical analysis of similarity matrix, returns smin,smax
    !! \param smin minimum similarity
    !! \param smax maximum similarity
    subroutine analyze_smat( s, symmetrize, smin, smax )
        real, intent(inout) :: s(:,:)          !< similarity matrix
        logical, intent(in) :: symmetrize      !< force diag symmetry
        real, intent(out)   :: smin, smax
        integer             :: i, j, n, npairs
        if( size(s,1) .ne. size(s,2) )then
            THROW_HARD('not a similarity matrix; analyze_smat')
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
    ! subroutine dev_from_smat( smat, i_median, sdev )
    !     real,    intent(in)  :: smat(:,:)
    !     integer, intent(out) :: i_median
    !     real,    intent(out) :: sdev
    !     real, allocatable :: sims(:)
    !     integer :: loc(1), i, j, n
    !     n = size(smat,1)
    !     if( n /= size(smat,2) ) THROW_HARD('symmetric similarity matrix assumed; stat :: dev_from_smat')
    !     allocate(sims(n))
    !     do i=1,n
    !         sims(i) = 0.0
    !         do j=1,n
    !             if( i /= j )then
    !                 sims(i) = sims(i) + smat(i,j)
    !             endif
    !         end do
    !     end do
    !     loc      = maxloc(sims)
    !     i_median = loc(1)
    !     sdev     = median(smat(i_median,:))
    ! end subroutine dev_from_smat

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

    ! the Z-score calculates the number of standard deviations a data point is away from the mean
    function z_scores( x, mask ) result( zscores )
        real,              intent(in) :: x(:) ! data points
        logical, optional, intent(in) :: mask(:)
        real, allocatable :: zscores(:)
        logical :: err
        allocate(zscores(size(x)), source=x)
        if( present(mask) )then
            call normalize(zscores, err, mask)
        else
            call normalize(zscores, err)
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

end module simple_stat
