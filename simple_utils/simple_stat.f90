module simple_stat
use simple_defs ! singleton
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

    ! VISUALIZATION

    !>  \brief  prints a primitive histogram given an array of real data
    subroutine plot_hist( arr, bins, sc )
        real, intent(in)    :: arr(:)
        integer, intent(in) :: bins, sc
        real                :: a(bins), minv, maxv, mean, std, var, binwidth
        integer             :: h(bins), bin, i, j, n
        logical             :: err
        write(*,*) ">>> GENERATING HISTOGRAM <<<"
        n = size(arr,1)
        call moment( arr, mean, std, var, err )
        minv = minval(arr)
        maxv = maxval(arr)
        write(*,*) "Number of values = ", n
        write(*,*) "Minimum value = ", minv
        write(*,*) "Maximum value = ", maxv
        write(*,*) "Mean value = ", mean
        write(*,*) "Standard deviation = ", std
        write(*,*) "Variance = ", var
        h = 0
        binwidth = (maxv-minv)/real(bins)
        do i=1,n 
            bin = nint((arr(i)-minv)/binwidth)  ! int(1.+(arr(i)-minv)/binwidth)
            if( bin < 1 )    bin = 1            ! check for underflows
            if( bin > bins ) bin = bins         ! check for overflows
            h(bin) = h(bin)+1
            a(bin) = a(bin)+arr(i)
        end do
        do i=1,bins
            a(i) = a(i)/real(h(i))
            write(*,*) a(i),h(i),"|",('#', j=1,nint(real(h(i))/real(sc)))
        end do
    end subroutine plot_hist

    ! MOMENTS/NORMALIZATION

    !>  \brief  given a 1D real array of data, this routine returns its mean: _ave_,
    !!          standard deviation: _sdev_, and variance: _var_
    subroutine moment_1( data, ave, sdev, var, err )
        !$ use omp_lib
        !$ use omp_lib_kinds
        real,    intent(out) :: ave, sdev, var
        logical, intent(out) :: err
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
    
    !>  \brief  given a 2D real array of data, this routine returns its mean: _ave_,
    !!          standard deviation: _sdev_, and variance: _var_
    subroutine moment_2( data, ave, sdev, var, err )
        !$ use omp_lib
        !$ use omp_lib_kinds
        real, intent(out)    :: ave, sdev, var
        logical, intent(out) :: err
        real, intent(in)     :: data(:,:)
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
    
    !>  \brief  given a 3D real array of data, this routine returns its mean: _ave_,
    !!          standard deviation: _sdev_, and variance: _var_
    subroutine moment_3( data, ave, sdev, var, err )
        !$ use omp_lib
        !$ use omp_lib_kinds
        real, intent(out)    :: ave, sdev, var
        logical, intent(out) :: err
        real, intent(in)     :: data(:,:,:)
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
    
    !>  \brief  is for statistical normalization of an array
    subroutine normalize_1( arr, err )
        real, intent(inout)  :: arr(:)
        real                 :: ave, sdev, var
        logical, intent(out) :: err
        call moment_1( arr, ave, sdev, var, err )
        if( err ) return
        arr = (arr-ave)/sdev ! array op    
    end subroutine normalize_1
    
    !>  \brief  is for statistical normalization of a 2D matrix
    subroutine normalize_2( arr, err )
        real, intent(inout)  :: arr(:,:)
        real                 :: ave, sdev, var
        logical, intent(out) :: err
        call moment_2( arr, ave, sdev, var, err )
        if( err ) return
        arr = (arr-ave)/sdev ! array op    
    end subroutine normalize_2
    
    !>  \brief  is for statistical normalization of a 3D matrix
    subroutine normalize_3( arr, err )
        real, intent(inout)  :: arr(:,:,:)
        real                 :: ave, sdev, var
        logical, intent(out) :: err
        call moment_3( arr, ave, sdev, var, err )
        if( err ) return
        arr = (arr-ave)/sdev ! array op    
    end subroutine normalize_3
    
    !>  \brief  is for sigmoid normalisation [0,1]
    subroutine normalize_sigm_1( arr ) 
        real, intent(inout) :: arr(:)
        real                :: smin, smax, delta
        real, parameter     :: nnet_const = exp(1.)-1.
        ! find minmax
        smin  = minval(arr)
        smax  = maxval(arr)
        delta = smax-smin
        ! create [0,1]-normalized vector
        !$omp parallel workshare default(shared)
        arr = (exp((arr-smin)/delta)-1.)/nnet_const
        !$omp end parallel workshare
    end subroutine normalize_sigm_1
    
    !>  \brief  is for sigmoid normalisation [0,1]
    subroutine normalize_sigm_2( arr ) 
        real, intent(inout) :: arr(:,:)
        real                :: smin, smax, delta
        real, parameter     :: nnet_const = exp(1.)-1.
        ! find minmax
        smin  = minval(arr)
        smax  = maxval(arr)
        delta = smax-smin
        ! create [0,1]-normalized vector
        !$omp parallel workshare default(shared)
        arr = (exp((arr-smin)/delta)-1.)/nnet_const
        !$omp end parallel workshare
    end subroutine normalize_sigm_2
    
    !>  \brief  is for sigmoid normalisation [0,1]
    subroutine normalize_sigm_3( arr ) 
        real, intent(inout) :: arr(:,:,:)
        real                :: smin, smax, delta
        real, parameter     :: nnet_const = exp(1.)-1.
        ! find minmax
        smin  = minval(arr)
        smax  = maxval(arr)
        delta = smax-smin
        ! create [0,1]-normalized vector
        !$omp parallel workshare default(shared)
        arr = (exp((arr-smin)/delta)-1.)/nnet_const
        !$omp end parallel workshare
    end subroutine normalize_sigm_3
    
    !>  \brief  calculates the devation around point
    subroutine deviation( data, point, sdev, var, err )
        !$ use omp_lib
        !$ use omp_lib_kinds
        real, intent(out)    :: sdev, var
        logical, intent(out) :: err
        real, intent(in)     :: data(:), point
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
    
    ! CORRELATION
    
    !>  \brief  calculates Pearson's correlation coefficient
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
    
    !>  \brief  calculates Pearson's correlation coefficient
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
    
    !>  \brief  calculates Pearson's correlation coefficient
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
   
    !>  \brief  calculates the Pearson correlation for pre-normalized data
    function pearsn_prenorm( x, y ) result( r )
        real    :: x(:),y(:),r,sxx,syy,sxy,den
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
        real, intent(in)  :: corrs(:)
        real, allocatable :: weights(:), corrs_copy(:), expnegdists(:)
        real, parameter   :: THRESHOLD=1.5
        real    :: maxminratio, normfac, corrmax, corrmin
        integer :: icorr, ncorrs
        ncorrs = size(corrs)
        allocate(weights(ncorrs), corrs_copy(ncorrs), expnegdists(ncorrs))
        weights     = 0. 
        corrs_copy  = corrs
        expnegdists = 0. 
        corrmax     = maxval(corrs_copy)
        if( corrmax < 0. )then 
            ! weighting does not make sense, put them all to 1/ncorrs
            weights = 1./real(ncorrs)
            return
        endif
        corrmin = minval(corrs_copy)
        maxminratio = 2.0
        if( corrmin > 0. ) maxminratio = corrmax/corrmin
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
            if( corrs_copy(icorr) >= 0. )then
                weights(icorr) = expnegdists(icorr)/normfac
            else
                weights(icorr) = 0.
            endif
        end do
        deallocate(corrs_copy)
    end function corrs2weights

    ! INTEGER STUFF
    
    !>  \brief  is for rank transformation of an array
    subroutine rank_transform_1( arr )
        use simple_jiffys, only: alloc_err
        real, intent(inout)  :: arr(:)
        integer              :: j, n, alloc_stat
        integer, allocatable :: order(:)
        real, allocatable    :: vals(:)
        n = size(arr)
        allocate( vals(n), order(n), stat=alloc_stat )
        call alloc_err("In: rank_transform_1; simple_stat", alloc_stat )
        do j=1,n
            order(j) = j
            vals(j)  = arr(j)
        end do  
        call hpsort(n, vals, order)
        do j=1,n
            arr(order(j)) = real(j)
        end do
        deallocate(vals, order)
    end subroutine rank_transform_1
    
    !>  \brief  is for rank transformation of a 2D matrix
    subroutine rank_transform_2( mat )
        use simple_jiffys, only: alloc_err
        real, intent(inout)   :: mat(:,:)
        integer, allocatable  :: order(:), indices(:,:)
        real, allocatable     :: vals(:)
        integer               :: n, alloc_stat, i, j, cnt, nx, ny
        nx = size(mat,1)
        ny = size(mat,2)
        n = nx*ny
        allocate( vals(n), order(n), indices(n,2), stat=alloc_stat )
        call alloc_err("In: rank_transform_2; simple_stat", alloc_stat )
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
        deallocate(vals, order, indices)
    end subroutine rank_transform_2
    
    !>  \brief  Spearman rank correlation
    function spear( n, pi1, pi2 ) result( corr )
        integer, intent(in) :: n
        real, intent(in)    :: pi1(n), pi2(n)
        real                :: corr, sqsum, rn
        integer             :: k
        rn = real(n)
        sqsum = 0. 
        do k=1,n
            sqsum = sqsum+(pi1(k)-pi2(k))**2.
        end do
        corr = 1.-(6.*sqsum)/(rn**3.-rn)
    end function spear

    !>  \brief  Kolmogorov-Smirnov test to deduce equivalence or non-equivalence between two distributions.
    !!          The routine returns the K-S statistic d, and the significance level prob for the null hypothesis
    !!          that the data sets are drawn from the same distribution. Small values for prob show that the cumulative
    !!          distribution function of data1 is significantly different from that of data2. The input arrays are
    !!          modified (sorted)
    subroutine kstwo( data1, n1, data2, n2, d, prob )
        integer, intent(in) :: n1, n2
        real, intent(inout) :: data1(n1), data2(n2), d, prob
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
    
    !>  \brief  4 statistical analysis of similarity matrix
    subroutine analyze_smat( s, symmetrize, smin, smax )
        real, intent(inout) :: s(:,:)
        logical, intent(in) :: symmetrize
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

    ! SPECIAL FUNCTIONS

    !>  \brief  is the factorial function
    recursive function factorial( n )result( f )
        integer, intent(in) :: n
        integer :: f
        if( n < 0 )then
            stop 'Negative factorial in simple_stat%factorial'
        else if( n==0 )then
            f = 1
        else
            f = n * factorial( n-1 )
        endif
    end function factorial

    !>  \brief  is the binomial coefficient
    function bin_coeff( n, k )result( val )
        integer, intent(in) :: n, k
        integer :: val
        if( k>n )stop 'Inconsistent input in simple_stat%bin_coeff 1'
        if( k<0 )stop 'Inconsistent input in simple_stat%bin_coeff 2'
        if( n<0 )stop 'Inconsistent input in simple_stat%bin_coeff 3'
        val = factorial(n) / (factorial(k) * factorial(n-k))
    end function bin_coeff

    !>  \brief  generates a primitive histogram given an array of real data
    function get_hist( arr, nbins )result( h )
        real, intent(in)     :: arr(:)
        integer, intent(in)  :: nbins
        real                 :: binwidth, minv, maxv
        integer, allocatable :: h(:)
        integer              :: bin, i, n
        if( nbins<2 )stop 'Invalib number of bins in simple_stat%get_hist'
        n        = size(arr,1)
        minv     = minval(arr)
        maxv     = maxval(arr)
        binwidth = ( maxv - minv ) / real( nbins )
        allocate( h(nbins) )
        h = 0
        do i=1,n 
            bin = nint((arr(i)-minv)/binwidth)   ! int(1.+(arr(i)-minv)/binwidth)
            if( bin < 1 )     bin = 1            ! check for underflows
            if( bin > nbins ) bin = nbins        ! check for overflows
            h(bin) = h(bin)+1
        end do
    end function get_hist

    !>  \brief  generates a primitive joint histogram given two arrays of real data
    function get_jointhist( arr1, arr2, nbins )result( h )
        real,    intent(in)  :: arr1(:), arr2(:)
        integer, intent(in)  :: nbins
        integer, allocatable :: h(:,:)
        real                 :: binwidth1, minv1
        real                 :: binwidth2, minv2
        integer              :: i, n, bin1, bin2
        if( nbins<2 )stop 'Invalib number of bins in simple_stat%get_hist'
        ! first array
        n         = size(arr1,1)
        if( n/=size(arr2,1) )stop 'Invalid dimensions in simple_stat%get_joint_hist'
        minv1     = minval(arr1)
        binwidth1 = ( maxval(arr1)-minv1 ) / real( nbins )
        ! second array
        minv2     = minval(arr2)
        binwidth2 = ( maxval(arr2)-minv2 ) / real( nbins )
        ! Joint
        allocate( h(nbins,nbins) )
        h = 0
        do i=1,n
            bin1 = bin( arr1(i), minv1, binwidth1 ) 
            bin2 = bin( arr2(i), minv2, binwidth2 ) 
            h( bin1, bin2 ) = h( bin1, bin2 ) + 1
        end do

        contains
            function bin( x, minv, width )result( ind )
                real, intent(in) :: x,minv,width
                integer :: ind
                ind = nint( (x-minv)/width ) ! int(1.+(arr(i)-minv)/binwidth)
                ind = max(1,ind)
                ind = min(ind,nbins)
            end function bin

    end function get_jointhist

    !>   \biref In information theory, the Hamming distance between two strings of equal length 
    !!          is the number of positions at which the corresponding symbols are different. 
    function hamming_dist( x, y ) result( dist )
        integer, intent(in) :: x(:), y(:)
        real :: dist
        if( size(x)/=size(y) )stop 'Invalid dimensions in simple_stat :: hamming_dist'
        dist = real(count( x /= y ))
    end function hamming_dist

    !>  \brief is the Normalized Mutual Information (in bits)
    function nmi( x, y, nbins )result( val )
        real,    intent(in)  :: x(:), y(:)
        integer, intent(in)  :: nbins
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
        if( ex/=0. .and. ey/=0. )then
            val = mi / sqrt(ex*ey)
        else
            val = 0.
        endif
        val = max(0.,val)
    end function nmi

    !>  \brief  calculates rand index
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

    !>  \brief  calculates Jaccard index
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

    !>  \brief  calculates Jaccard distance
    function jaccard_dist( x, y )result( jdist )
        integer, intent(in) :: x(:), y(:)
        real :: jdist
        if( size(x) /= size(y) )stop 'Inconsistent dimensions in simple_stat%jaccard_dist'
        jdist = 1. - jaccard_index( x,y )
    end function jaccard_dist

end module simple_stat
