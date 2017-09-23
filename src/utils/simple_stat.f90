! statistics utility functions
module simple_stat
    use simple_defs ! singleton
    use simple_syslib
use simple_math, only: hpsort
implicit none

private :: moment_1, moment_2, moment_3, normalize_1, normalize_2, normalize_3

interface moment
    module procedure moment_1
    module procedure moment_2
    module procedure moment_3
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
end interface

interface normalize_sigm
    module procedure normalize_sigm_1
    module procedure normalize_sigm_2
    module procedure normalize_sigm_3
end interface

interface rank_transform
    module procedure rank_transform_1
    module procedure rank_transform_2
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
            write(*,*) 'ERROR: n must be at least 2'
            write(*,*) 'In: moment_1, module: simple_stat.f90'
            stop
        endif
        ! calc average
        ave = sum(data)/nr
        ! calc sum of devs and sum of devs squared
        ep = 0.
        var = 0.
        !$omp parallel do default(shared) private(i,dev) schedule(auto) reduction(+:ep,var)
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
            write(*,*) 'ERROR: n must be at least 2'
            write(*,*) 'In: moment_2, module: simple_stat.f90'
            stop
        endif
        ! calc average
        ave = sum(data)/nr
        ! calc sum of devs and sum of devs squared
        ep = 0.
        var = 0.
        !$omp parallel do default(shared) private(i,j,dev) schedule(auto) reduction(+:ep,var)
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
            write(*,*) 'ERROR: n must be at least 2'
            write(*,*) 'In: moment_3, module: simple_stat.f90'
            stop
        endif
        ave = sum(data)/nr
        ! calc sum of devs and sum of devs squared
        ep = 0.
        var = 0.
        !$omp parallel do default(shared) private(i,j,k,dev) schedule(auto) reduction(+:ep,var)
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

    !>    is for statistical normalization of an array
    subroutine normalize_1( arr, err )
        real, intent(inout)  :: arr(:)           !< input data
        real                 :: ave, sdev, var   !< temp stats
        logical, intent(out) :: err              !< error status
        call moment_1( arr, ave, sdev, var, err )
        if( err ) return
        arr = (arr-ave)/sdev ! array op
    end subroutine normalize_1

    !>    is for statistical normalization of a 2D matrix
    subroutine normalize_2( arr, err )
        real, intent(inout)  :: arr(:,:)         !< input data
        real                 :: ave, sdev, var   !< temp stats
        logical, intent(out) :: err              !< error status
        call moment_2( arr, ave, sdev, var, err )
        if( err ) return
        arr = (arr-ave)/sdev ! array op
    end subroutine normalize_2

    !>    is for statistical normalization of a 3D matrix
    subroutine normalize_3( arr, err )
        real, intent(inout)  :: arr(:,:,:)       !< input data
        real                 :: ave, sdev, var   !< temp stats
        logical, intent(out) :: err              !< error status
        call moment_3( arr, ave, sdev, var, err )
        if( err ) return
        arr = (arr-ave)/sdev ! array op
    end subroutine normalize_3

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
            write(*,'(a)') 'WARNING! stat :: normalize_sigm_1, division with zero'
            write(*,'(a)') 'no normalisation done'
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
            write(*,'(a)') 'WARNING! stat :: normalize_sigm_2, division with zero'
            write(*,'(a)') 'no normalisation done'
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
            write(*,'(a)') 'WARNING! stat :: normalize_sigm_2, division with zero'
            write(*,'(a)') 'no normalisation done'
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
            write(*,'(a)') 'WARNING! stat :: normalize_minmax, division with zero'
            write(*,'(a)') 'no normalisation done'
        endif
    end subroutine normalize_minmax

    !>    calculates the devation around point
    !! \param data input data
    !! \param sdev standard deviation
    subroutine deviation( data, point, sdev, var, err )
        !$ use omp_lib
        !$ use omp_lib_kinds
        real, intent(out)    :: sdev, var       !< var variance
        logical, intent(out) :: err             !< error status
        real, intent(in)     :: data(:), point  !< input data and deviation point
        integer              :: n, i
        real                 :: ep, nr, dev
        err = .false.
        n   = size(data,1)
        nr  = n
        ! calc sum of devs and sum of devs squared
        ep = 0.
        var = 0.
        !$omp parallel do default(shared) private(i,dev) schedule(auto) reduction(+:ep,var)
        do i=1,n
            dev = data(i)-point
            ep = ep+dev
            var = var+dev*dev
        end do
        !$omp end parallel do
        var = (var-ep**2./nr)/(nr-1.) ! corrected two-pass formula
        sdev = sqrt(var)
        if( abs(var) < TINY ) err = .true.
    end subroutine deviation

    !>    given a 2D real array of data, this routine returns its windowed mean
    !! \param data input data
    !! \param winsz window size
    !! \param ave output array
    subroutine mean_2D( data, winsz, ave, err )
        !$ use omp_lib
        !$ use omp_lib_kinds
        real, intent(out), allocatable    :: ave(:,:)
        integer, intent(in) :: winsz
        logical, intent(out) :: err              !< error status
        real, intent(in)     :: data(:,:)        !< input data
        integer              :: nx, ny, i, j, h, k, px, py, tmpcnt
        real                 :: n, tmpave!, nrpatch
        err = .false.
        nx = size(data,1)
        ny = size(data,2)
        allocate(ave(nx,ny))
        n  = nx*ny
        if( n <= 1 ) then
            write(*,*) 'ERROR: n must be at least 2'
            write(*,*) 'In: mean_2D, module: simple_stat.f90'
            stop
        endif
        if( winsz > nx/2 .or. winsz > ny/2  ) then
            write(*,*) 'ERROR: winsz must be smaller than half image dimensions'
            write(*,*) 'In: mean_2D, module: simple_stat.f90'
            stop
        endif
        ! calc average
        !$omp parallel do default(shared) private(i,j) schedule(auto)
        do i=1,nx
            do j=1,ny
                tmpave=0.
                tmpcnt=0
                do px=-winsz,winsz
                    ! reflect on x- boundary
                    h= i+px
                    h = merge(i - px,  h,  (h < 1).or.(h > nx)) 
                    do py=-winsz,winsz
                        k= j + py
                        k = merge(j - py, k, (k < 1).or.(k > ny))
                        ! calc avg of window
                        tmpave = tmpave + data(h,k)
                        tmpcnt = tmpcnt+1
                    end do
                end do
                ! calc avg of window
                if( tmpcnt < 1) err = .true.
                ave(i,j) = tmpave / REAL(tmpcnt)

            end do
        end do
        !$omp end parallel do
    
    end subroutine mean_2D

    !> stdev_2D Standard deviation 2D filter
    !! \param data 2D array
    !! \param winsz window size
    !! \param ave Average (mean) 2D filter
    !! \param stdev Output std dev filter
    !! \param var Optional output variance 2D array
    !! \param err error flag
    !!
    subroutine stdev_2D( data, winsz, ave, stdev, var, err)
        !$ use omp_lib
        !$ use omp_lib_kinds
        integer, intent(in) :: winsz
        real, intent(in)     :: data(:,:), ave(:,:)        !< input data
        real, intent(out), allocatable  :: stdev(:,:),var(:,:)
        logical, intent(out) :: err              !< error status
        integer              :: nx, ny, i, j, h, k,px,py, tmpcnt
        real                 :: n, tmpdev, nr, tmpvar, ep
        err = .false.
        nx = size(data,1)
        ny = size(data,2)
        allocate(stdev(nx,ny),var(nx,ny))
        n  = nx*ny
        if( n <= 1 ) then
            write(*,*) 'ERROR: n must be at least 2'
            write(*,*) 'In: mean_2D, module: simple_stat.f90'
            stop
        endif
        if( winsz > nx/2 .or. winsz > ny/2  ) then
            write(*,*) 'ERROR: winsz must be smaller than half image dimensions'
            write(*,*) 'In: mean_2D, module: simple_stat.f90'
            stop
        endif
        ! calc average
        !$omp parallel do default(shared) private(i,j) schedule(auto)
        do i=1,nx
            do j=1,ny
                ep=0.;tmpdev=0.; tmpvar=0.
                tmpcnt=0
                do px=-winsz,winsz
                    ! reflect on x- boundary
                    h= i+px
                    h = merge(i - px,  h,  (h < 1).or.(h > nx)) 
                    do py=-winsz,winsz
                        k= j + py
                        k = merge(j - py, k, (k < 1).or.(k > ny))
                        tmpdev = data(h,k)-ave(h,k)
                        ep = ep+tmpdev
                        tmpvar = tmpvar+tmpdev*tmpdev
                        tmpcnt = tmpcnt+1
                    end do
                end do
                ! calc avg of window
                nr = REAL(tmpcnt)
                var(i,j) = (tmpvar-( (ep**2.) / nr)) / REAL(nr -1.)
                stdev(i,j) = sqrt(var(i,j)); if( abs(var(i,j)) < TINY ) err = .true.
            end do
        end do
        !$omp end parallel do
    
    end subroutine stdev_2D
    
    ! CORRELATION

    !>    calculates Pearson's correlation coefficient
    !! \param x input reference array
    !! \param y input test array
    function pearsn_1( x, y ) result( r )
        real, intent(in) :: x(:),y(:)
        real    :: r,ax,ay,sxx,syy,sxy,xt,yt
        integer :: j, n
        n = size(x)
        if( size(y) /= n ) stop 'Arrays not equal size, in pearsn_1, module: simple_stat'
        ax  = sum(x)/real(n)
        ay  = sum(y)/real(n)
        sxx = 0.
        syy = 0.
        sxy = 0.
        !$omp parallel do default(shared) private(j,xt,yt) &
        !$omp reduction(+:sxx,syy,sxy) schedule(auto)
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
    function pearsn_2( x, y ) result( r )
        real, intent(in) :: x(:,:),y(:,:)
        real    :: r,ax,ay,sxx,syy,sxy,xt,yt
        integer :: i, j, nx, ny
        nx = size(x,1)
        ny = size(x,2)
        if( size(y,1) /= nx .or. size(y,2) /= ny ) stop 'Arrays not equal size, in pearsn_2, module: simple_stat'
        ax  = sum(x)/real(nx*ny)
        ay  = sum(y)/real(nx*ny)
        sxx = 0.
        syy = 0.
        sxy = 0.
        !$omp parallel do default(shared) private(i,j,xt,yt) &
        !$omp reduction(+:sxx,syy,sxy) schedule(auto)
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
        real    :: r,ax,ay,az,sxx,syy,sxy,xt,yt
        integer :: i, j, k, nx, ny, nz
        nx = size(x,1)
        ny = size(x,2)
        nz = size(x,3)
        if( size(y,1) /= nx .or. size(y,2) /= ny .or. size(y,3) /= nz )&
        stop 'Arrays not equal size, in pearsn_3, module: simple_stat'
        ax  = sum(x)/real(nx*ny*nz)
        ay  = sum(y)/real(nx*ny*nz)
        sxx = 0.
        syy = 0.
        sxy = 0.
        !$omp parallel do default(shared) private(i,j,k,xt,yt) &
        !$omp reduction(+:sxx,syy,sxy) schedule(auto)
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

    !>    calculates the Pearson correlation for pre-normalized data
    !! \param x input reference array
    !! \param y input test array
    function pearsn_prenorm( x, y ) result( r )
        real    :: x(:),y(:)
        real    :: r,sxx,syy,sxy,den
        integer :: j, n
        n = size(x)
        if( size(y) /= n ) stop 'Arrays not equal size, in pearsn_prenorm, module: simple_stat'
        sxx = 0.
        syy = 0.
        sxy = 0.
        !$omp parallel do default(shared) private(j) &
        !$omp reduction(+:sxx,syy,sxy) schedule(auto)
        do j=1,n
            sxx = sxx+x(j)**2
            syy = syy+y(j)**2
            sxy = sxy+x(j)*y(j)
        end do
        !$omp end parallel do
        den = sxx*syy
        if( den > 0. )then
            r = sxy/sqrt(den)
        else
            r = 0.
        endif
    end function pearsn_prenorm

    function corrs2weights( corrs ) result( weights )
        real, intent(in)  :: corrs(:) !< correlation input
        real, allocatable :: weights(:), corrs_copy(:), expnegdists(:)
        real, parameter   :: THRESHOLD=1.5
        real    :: maxminratio, normfac, corrmax, corrmin
        integer :: icorr, ncorrs, alloc_stat
        ncorrs = size(corrs)
        allocate(weights(ncorrs), corrs_copy(ncorrs), expnegdists(ncorrs),stat=alloc_stat)
        call alloc_errchk("In: corrs2weights; simple_stat", alloc_stat )
        weights     = 0.
        corrs_copy  = corrs
        expnegdists = 0.
        corrmax     = maxval(corrs_copy)
        if( corrmax < 0. )then
            ! weighting does not make sense, put them all to 1/ncorrs
            weights = 1./real(ncorrs)
            return
        endif
        corrmin = minval(corrs_copy, mask=corrs_copy > TINY)
        maxminratio = corrmax/corrmin
        if( maxminratio >= THRESHOLD )then
           ! min/max normalise the correlations
           call normalize_sigm(corrs_copy)
        endif
        ! calculate the exponential of the negative distances
        ! so that when diff==0 the weights are maximum and when
        ! diff==corrmax the weights are minimum
        do icorr=1,ncorrs
            if( corrs_copy(icorr) >= 0. )then
                expnegdists(icorr) = exp(-(1.-corrs_copy(icorr)))
            else
                expnegdists(icorr) = 0.
            endif
        end do
        ! calculate weight normalisation factor
        normfac = sum(expnegdists)
        ! normalise and set weights where corrs==0 to 0
        do icorr=1,ncorrs
            if( corrs_copy(icorr) > TINY )then
                weights(icorr) = expnegdists(icorr)/normfac
            else
                weights(icorr) = 0.
            endif
        end do
    end function corrs2weights

    ! INTEGER STUFF

    !>    is for rank transformation of an array
    subroutine rank_transform_1( arr )
        real, intent(inout)  :: arr(:)  !< array to be modified
        integer              :: j, n, alloc_stat
        integer, allocatable :: order(:)
        real, allocatable    :: vals(:)
        n = size(arr)
        allocate( vals(n), order(n), stat=alloc_stat )
        call alloc_errchk("In: rank_transform_1; simple_stat", alloc_stat )
        do j=1,n
            order(j) = j
            vals(j)  = arr(j)
        end do
        call hpsort(n, vals, order)
        do j=1,n
            arr(order(j)) = real(j)
        end do
        deallocate(vals, order, stat=alloc_stat )
        call alloc_errchk("In: rank_transform_1; simple_stat", alloc_stat )
    end subroutine rank_transform_1

    !>    is for rank transformation of a 2D matrix
    subroutine rank_transform_2( mat )
        real, intent(inout)   :: mat(:,:)  !< matrix to be modified
        integer, allocatable  :: order(:), indices(:,:)
        real, allocatable     :: vals(:)
        integer               :: n, alloc_stat, i, j, cnt, nx, ny
        nx = size(mat,1)
        ny = size(mat,2)
        n = nx*ny
        allocate( vals(n), order(n), indices(n,2), stat=alloc_stat )
        call alloc_errchk("In: rank_transform_2; simple_stat", alloc_stat )
        cnt = 0
        do i=1,nx
            do j=1,ny
                cnt            = cnt+1
                indices(cnt,:) = [i,j]
                order(cnt)     = cnt
                vals(cnt)      = mat(i,j)
            end do
        end do
        call hpsort(n, vals, order)
        ! convert to rank
        do j=1,n
            mat(indices(order(j),1),indices(order(j),2)) = real(j)
        end do
        deallocate(vals, order, indices, stat=alloc_stat )
         call alloc_errchk("In: rank_transform_2; simple_stat", alloc_stat )
    end subroutine rank_transform_2

    !>    Spearman rank correlation
    !<  \param pi1 test data 1 and \param pi2 reference data
    function spear( n, pi1, pi2 ) result( corr )
        integer, intent(in) :: n
        real, intent(in)    :: pi1(n), pi2(n)  !<  pi1 test data 1  pi2 test data 2
        real                :: corr, sqsum, rn
        integer             :: k
        rn = real(n)
        sqsum = 0.
        do k=1,n
            sqsum = sqsum+(pi1(k)-pi2(k))**2.
        end do
        corr = 1.-(6.*sqsum)/(rn**3.-rn)
    end function spear

    !>   Kolmogorov-Smirnov test to deduce equivalence or non-equivalence
    !>  between two distributions.
    !!          The routine returns the K-S statistic d, and the significance
    !!          level prob for the null hypothesis that the data sets are drawn
    !!          from the same distribution. Small values for prob show that the
    !!          cumulative distribution function of data1 is significantly
    !!          different from that of data2. The input arrays are modified
    !!          (sorted)
    !! \param  data1,data2 distribution arrays
    !! \param  n1,n2 size of distribution arrays
    !! \param  d K-S statistic
    subroutine kstwo( data1, n1, data2, n2, d, prob )
        integer, intent(in) :: n1, n2
        real, intent(inout) :: data1(n1), data2(n2), d, prob  !< significance
        integer             :: j1, j2
        real                :: d1, d2, dt, en1, en2, en, fn1, fn2
        call hpsort( n1, data1 )
        call hpsort( n2, data2 )
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
            stop 'not a similarity matrix; analyze_smat; simple_stat'
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
        use simple_math, only: median
        real,    intent(in)  :: smat(:,:)
        integer, intent(out) :: i_median
        real,    intent(out) :: sdev
        real, allocatable :: sims(:)
        integer :: loc(1), i, j, n
        n = size(smat,1)
        if( n /= size(smat,2) ) stop 'symmetric similarity matrix assumed; stat :: median_dev_from_smat'
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
        use simple_math, only: median
        real,    intent(in)  :: dmat(:,:)
        integer, intent(out) :: i_median
        real,    intent(out) :: ddev
        real, allocatable :: dists(:)
        integer :: loc(1), i, j, n
        n = size(dmat,1)
        if( n /= size(dmat,2) ) stop 'symmetric distance matrix assumed; stat :: median_dev_from_dmat'
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

    ! SPECIAL FUNCTIONS

    !>    is the factorial function
    recursive function factorial( n )result( f )
        integer, intent(in) :: n !< factorial arg
        integer :: f
        if( n < 0 )then
            stop 'Negative factorial in simple_stat%factorial'
        else if( n==0 )then
            f = 1
        else
            f = n * factorial( n-1 )
        endif
    end function factorial

    !>    is the binomial coefficient
    function bin_coeff( n, k )result( val )
        integer, intent(in) :: n, k
        integer :: val
        if( k>n )stop 'Inconsistent input in simple_stat%bin_coeff 1'
        if( k<0 )stop 'Inconsistent input in simple_stat%bin_coeff 2'
        if( n<0 )stop 'Inconsistent input in simple_stat%bin_coeff 3'
        val = factorial(n) / (factorial(k) * factorial(n-k))
    end function bin_coeff

    !>    generates a primitive histogram given an array of real data
    function get_hist( arr, nbins )result( h )
        real, intent(in)     :: arr(:) !< input array
        integer, intent(in)  :: nbins  !< num histogram bins
        real                 :: binwidth, minv, maxv
        integer, allocatable :: h(:)
        integer              :: bin, i, n, alloc_stat
        if( nbins<2 )stop 'Invalib number of bins in simple_stat%get_hist'
        n        = size(arr,1)
        minv     = minval(arr)
        maxv     = maxval(arr)
        binwidth = ( maxv - minv ) / real( nbins )
        allocate( h(nbins) , stat=alloc_stat )
        call alloc_errchk('In: simple_stat; get_hist', alloc_stat)
        h = 0
        do i=1,n
            bin = nint((arr(i)-minv)/binwidth)   ! int(1.+(arr(i)-minv)/binwidth)
            if( bin < 1 )     bin = 1            ! check for underflows
            if( bin > nbins ) bin = nbins        ! check for overflows
            h(bin) = h(bin)+1
        end do
    end function get_hist

    !>    generates a primitive joint histogram given two arrays of real data
    !! \param x input reference array
    !! \param y input test array
    function get_jointhist( x, y, nbins )result( h )
        real,    intent(in)  :: x(:), y(:)
        integer, intent(in)  :: nbins              !< num histogram bins
        integer, allocatable :: h(:,:)             !< output histogram
        real                 :: binwidth1, minv1
        real                 :: binwidth2, minv2
        integer              :: i, n, bin1, bin2, alloc_stat
        if( nbins<2 )stop 'Invalib number of bins in simple_stat%get_hist'
        ! first array
        n         = size(x,1)
        if( n/=size(y,1) )stop 'Invalid dimensions in simple_stat%get_joint_hist'
        minv1     = minval(x)
        binwidth1 = ( maxval(x)-minv1 ) / real( nbins )
        ! second array
        minv2     = minval(y)
        binwidth2 = ( maxval(y)-minv2 ) / real( nbins )
        ! Joint
        allocate( h(nbins,nbins) , stat=alloc_stat )
        call alloc_errchk('In: simple_stat; get_jointhist', alloc_stat) 
        h = 0
        do i=1,n
            bin1 = bin( x(i), minv1, binwidth1 )
            bin2 = bin( y(i), minv2, binwidth2 )
            h( bin1, bin2 ) = h( bin1, bin2 ) + 1
        end do

        contains
            function bin( v, minv, width )result( ind )
                real, intent(in) :: v,minv,width
                integer :: ind
                ind = nint( (v-minv)/width ) ! int(1.+(arr(i)-minv)/binwidth)
                ind = max(1,ind)
                ind = min(ind,nbins)
            end function bin

    end function get_jointhist

    !>    In information theory, the Hamming distance between two strings of equal length
    !!          is the number of positions at which the corresponding symbols are different.
    !! \param x input reference array
    !! \param y input test array
    function hamming_dist( x, y ) result( dist )
        integer, intent(in) :: x(:), y(:)
        real :: dist
        if( size(x)/=size(y) )stop 'Invalid dimensions in simple_stat :: hamming_dist'
        dist = real(count( x /= y ))
    end function hamming_dist

    !> distance metric based on the likelihood ratio
    !! could be useful for Fourier transforms
    real function likelihood_ratio( a, b )
        real, intent(in) :: a(:), b(:)
        if( size(a) /= size(b) ) stop 'ERROR, noncongruent arrays; stat :: likelihood_ratio '
        likelihood_ratio = sum(2.0 * log(a + b) - log(a) - log(b))
    end function likelihood_ratio

    !>   is the Normalized Mutual Information (in bits)
    !! Mutual Information \f$ I(X;Y)=\sum_{y\in Y}\sum_{x\in X}p(x,y)\log{\left({\frac{p(x,y)}{p(x)\,p(y)}}\right)} \f$
    !!
    !! Normalised MI also known as Information Quality Ratio (IQR)
    !! \f[ IQR(X,Y)=E\left[I(X;Y)\right]={\frac{I(X;Y)}{\mathrm{H}(X,Y)}}=
    !!   {\frac{\sum_{x\in X}\sum_{y\in Y}p(x,y)\log {p(x)p(y)}}
    !!   {\sum _{x\in X}\sum _{y\in Y}p(x,y)\log {p(x,y)}}-1}
    !! \f]
    !! \see https://en.wikipedia.org/wiki/Mutual_information#Normalized_variants
    !! \param x input reference array
    !! \param y input test array
    function nmi( x, y, nbins )result( val )
        real,    intent(in)  :: x(:), y(:)
        integer, intent(in)  :: nbins         !< num histogram bins
        real,    allocatable :: rh(:,:), pxs(:), pys(:)
        real    :: mi, val, ex, ey, pxy, px, py, logtwo
        integer :: i, j, n
        if( nbins<2 )stop 'Invalid number of bins in simple_stat%nmi'
        ! Init
        logtwo  = log(2.)
        n       = size(x)
        if( n/=size(y) )stop 'Invalid dimensions in simple_stat%nmi'
        mi = 0.
        ex = 0.
        ey = 0.
        allocate( rh(nbins,nbins), pxs(nbins), pys(nbins) )
        rh = real( get_jointhist( x, y, nbins ) ) / real( n )
        ! marginal entropies
        do i=1,nbins
            px     = sum( rh(i,:) )
            py     = sum( rh(:,i) )
            pxs(i) = px
            pys(i) = py
            if( px>0. )ex = ex - px * log(px) / logtwo
            if( py>0. )ey = ey - py * log(py) / logtwo
        enddo
        ! mutual information
        do i=1,nbins
            px = pxs(i)
            if( px<=0. )cycle
            do j=1,nbins
                py  = pys(j)
                if( px<=0. )cycle
                pxy = rh(i,j)
                if( pxy<=0. )cycle
                mi = mi + pxy * log( pxy/ (px*py) ) / logtwo
            enddo
        enddo
        if( ( ex /= 0.0) .and. (ey /= 0.0) )then
            val = mi / sqrt(ex*ey)
        else
            val = 0.0
        endif
        val = max(0.,val)
    end function nmi

    !>    calculates rand index
    !! \param x input reference array
    !! \param y input test array
    function rand_index( x, y ) result( rind )
        integer, intent(in) :: x(:),y(:)
        integer  :: i,j, n, a,b,cd
        real     :: rind
        n = size(x)
        if( n /= size(y) )stop 'Inconsistent dimensions in simple_stat%rand_index'
        a  = 0 ! common co-occuring pairs
        b  = 0 ! common not co-occuring pairs
        cd = 0 ! the rest
        do i=1,n
            do j=i+1,n
                if ( x(i)==x(j) )then
                    if( y(i)==y(j) )then
                        a  =  a + 1
                    else
                        cd = cd + 1
                    endif
                else
                    if( y(i)/=y(j) )then
                        b  =  b + 1
                    else
                        cd = cd + 1
                    endif
                endif
            enddo
        enddo
        rind = real(a+b) / real(a+b+cd)
        if( rind>1.0001 )stop 'Error in simple_stat%rand_index 1'
        if( rind<0. )stop 'Error in simple_stat%rand_index 2'
    end function rand_index

    !>    calculates Jaccard index
    !! \f$  J(A,B) = \frac{|A \cap B|}{|A \cup B|} = \frac{|A \cap B|}{|A| + |B| - |A \cap B|} \f$
    !! \param x input reference array
    !! \param y input test array
    function jaccard_index( x, y )result( jind )
        integer, intent(in) :: x(:), y(:)
        real :: jind
        integer  :: i,j, n, a,cd
        n = size(x)
        if( n /= size(y) )stop 'Inconsistent dimensions in simple_stat%jaccard_index'
        a  = 0 ! co-occuring pairs
        cd = 0 ! the rest and excluding common not co-occuring pairs (b in rand index)
        do i=1,n
            do j=i+1,n
                if ( x(i)==x(j) )then
                    if( y(i)==y(j) )then
                        a  =  a + 1
                    else
                        cd = cd + 1
                    endif
                else
                    if( y(i)==y(j) )cd = cd + 1
                endif
            enddo
        enddo
        jind = real(a) / real(a+cd)
        if( jind>1.0001 )stop 'Error in simple_stat%jaccard_index 1'
        if( jind<0. )stop 'Error in simple_stat%jaccard_index 2'
    end function jaccard_index

    !>    calculates Jaccard distance
    !! \f$  d_J(A,B) = 1 - J(A,B) = \frac{ |A \cup B| - |A \cap B| }{ |A \cup B| } \f$
    !! \param x input reference array
    !! \param y input test array
    function jaccard_dist( x, y )result( jdist )
        integer, intent(in) :: x(:), y(:)
        real :: jdist
        if( size(x) /= size(y) )stop 'Inconsistent dimensions in simple_stat%jaccard_dist'
        jdist = 1. - jaccard_index( x,y )
    end function jaccard_dist

end module simple_stat
