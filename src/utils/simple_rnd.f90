! random number generation module

module simple_rnd
use simple_defs ! singleton
use simple_math
use simple_error,  only: allocchk, simple_exception
use simple_syslib, only: get_process_id
implicit none

public :: seed_rnd, ran3, ran3arr, randn, multinomal, gasdev, irnd_uni, irnd_uni_pair
public :: irnd_gasdev, rnd_4dim_sphere_pnt, shcloc, mnorm_smp, r8po_fa
private
#include "simple_local_flags.inc"

interface ran3arr
    module procedure ran3arr_1
    module procedure ran3arr_2
end interface

interface gasdev
    module procedure gasdev_1
    module procedure gasdev_2
    module procedure gasdev_3
end interface

interface randn
    module procedure randn_1
    module procedure randn_2
end interface

integer(long), save :: idum

contains

    !>  \brief  random seed
    !>  solution from https://stackoverflow.com/questions/34797938/random-number-generator-in-pgi-fortran-not-so-random
    !>  the old version was replaced because the bug described in the thread was observed with GCC and PGI
    subroutine seed_rnd
        use iso_fortran_env, only: int64
        integer, allocatable :: seed(:)
        integer              :: i, n, istat, dt(8), pid
        integer(int64)       :: t
        integer, parameter   :: un=703
        call random_seed(size = n)
        allocate(seed(n))
        ! First try if the OS provides a random number generator
        open(unit=un, file="/dev/urandom", access="stream", &
            form="unformatted", action="read", status="old", iostat=istat)
        if (istat == 0) then
            read(un) seed
            close(un)
        else
            ! The PID is
            ! useful in case one launches multiple instances of the same
            ! program in parallel.
            call system_clock(t)
            if (t == 0) then
                call date_and_time(values=dt)
                t = (dt(1) - 1970) * 365_int64 * 24 * 60 * 60 * 1000 &
                    + dt(2) * 31_int64 * 24 * 60 * 60 * 1000 &
                    + dt(3) * 24_int64 * 60 * 60 * 1000 &
                    + dt(5) * 60 * 60 * 1000 &
                    + dt(6) * 60 * 1000 + dt(7) * 1000 &
                    + dt(8)
            end if
            pid = get_process_id()
            t = ieor( t, int(pid, kind(t)) )
            do i = 1, n
                seed(i) = lcg(t)
            end do
        end if
        call random_seed(put=seed)

    contains

        function lcg(s) result(res)
            integer                       :: res
            integer(int64), intent(inout) :: s
            if (s == 0) then
                s = 104729
            else
                s = mod(s, 4294967296_int64)
            end if
            s = mod(s * 279470273_int64, 4294967291_int64)
            res = int(mod(s, int(huge(0), 8)), kind(0))
        end function lcg

    end subroutine seed_rnd

    !>  \brief  wrapper for the intrinsic Fortran random number generator
    function ran3( ) result( harvest )
        real :: harvest
        call random_number(harvest)
    end function ran3

    !>  \brief  returns a random matrix
    subroutine ran3arr_1( harvest )
        real, intent(inout) :: harvest(:)
        call random_number(harvest)
    end subroutine ran3arr_1

     !>  \brief  returns a random matrix
    subroutine ran3arr_2( harvest )
        real, intent(inout) :: harvest(:,:)
        call random_number(harvest)
    end subroutine ran3arr_2

    !>  \brief  returns a random array [-1,1]
    !! \param n dimension of required rand array
    function randn_1( n ) result( a )
        integer, intent(in) :: n
        real, allocatable   :: a(:)
        allocate( a(n) )
        call random_number(a)
        a = a * 2. - 1.
    end function randn_1

    !>  \brief  returns a random matrix [-1,1]
    !! \param n1,n2 dimensions of required rand array
    function randn_2( n1, n2 ) result( a )
        integer, intent(in) :: n1, n2
        real, allocatable   :: a(:,:)
        allocate( a(n1,n2) )
        call random_number(a)
        a = a * 2. - 1.
    end function randn_2

    !>  \brief  generates a multinomal 1-of-K random number according to the
    !!          distribution in pvec
    function multinomal( pvec ) result( which )
        use simple_math, only: hpsort
        real,     intent(in) :: pvec(:) !< probabilities
        real,    allocatable :: pvec_sorted(:)
        integer, allocatable :: inds(:)
        integer :: i, n, which
        real    :: rnd, bound
        n = size(pvec)
        allocate(pvec_sorted(n), source=pvec)
        inds = (/(i,i=1,n)/)
        call hpsort(pvec_sorted, inds)
        rnd = ran3()
        do which=n,1,-1
            bound = sum(pvec_sorted(which:n))
            if( rnd <= bound )exit
        enddo
        which = inds(which)
        deallocate(pvec_sorted,inds)
    end function multinomal

    !>  \brief  random number generator yielding normal distribution
    !!          with zero mean and unit variance (from NR)
    function gasdev_1( ) result( harvest )
        real                   :: v1=0., v2=0., r, fac, harvest
        real,save              :: gset
        integer,save           :: iset=0
        if( idum < 0 ) iset = 0
        if( iset == 0 )then ! reinitialize
            r = 99.
            do while( r >= 1. .or. r < TINY )
                v1 = 2.*ran3( )-1.
                v2 = 2.*ran3( )-1.
                r  = v1**2.+v2**2.
            end do
            fac = sqrt(-2.*log(r)/r)
            ! Now make Box-Muller transformation to get the two normal deviates.
            ! Return one and save the other for the next call.
            gset    = v1*fac
            iset    = 1
            harvest = v2*fac
        else
            iset    = 0
            harvest = gset
        endif
    end function gasdev_1

    !>  \brief  random number generator yielding normal distribution
    !!          with given _mean_ and _stdev_ (from NR). standard _mean_
    !!          and _stdev_ values are 0 and 1, respectively
    function gasdev_2( mean, stdev ) result( harvest )
        real, intent(in) :: mean, stdev
        real             :: harvest
        harvest = stdev * gasdev( ) + mean
    end function gasdev_2

    !>  \brief  random number generator yielding normal distribution with
    !!          given mean and stdev (from NR). Added acceptance-rejection
    !!          according to limits (for constrained CE optimization)
    function gasdev_3( mean, stdev, limits ) result( harvest )
        real, intent(in)   :: mean, stdev
        real, intent(in)   :: limits(2)
        real               :: harvest
        integer, parameter :: maxits=500
        integer            :: cnt
        harvest = limits(2) + 99.
        cnt = 0
        do while( harvest < limits(1) .or. harvest > limits(2) )
            harvest = gasdev( mean, stdev )
            cnt = cnt + 1
            if( cnt == maxits )then
                THROW_WARN('gasdev_3 exceeding maxits')
                return
             endif
        end do
    end function gasdev_3

    !>  \brief  generate a uniformly distributed random integer [_1_,_NP_]
    function irnd_uni( NP ) result( irnd )
        integer, intent(in) :: NP
        integer             :: irnd
        irnd = 1
        if( NP == 0 )then
            write(logfhandle,*) 'Uniform random integer must be generated from a non-empty set!'
            write(logfhandle,*) 'In: irnd_uni, module: simple_rnd'
            stop
        else if( NP == 1 )then
            irnd = 1
        else
            irnd = ceiling(ran3()*real(NP))
            irnd = max(1,irnd)
            irnd = min(NP,irnd)
        endif
    end function irnd_uni

    !>  \brief  generate a pair of disjoint uniformly distributed random integers [_1_,_NP_]
    function irnd_uni_pair( NP ) result( rp )
        integer, intent(in) :: NP
        integer             :: rp(2)
        rp(1) = irnd_uni( NP )
        rp(2) = irnd_uni( NP )
        do while( rp(1) == rp(2) )
            rp(2) = irnd_uni( NP )
        end do
    end function irnd_uni_pair

    !>  \brief  generates a normally distributed random integer [_1_,_NP_]
    function irnd_gasdev( mean, stdev, NP )result( irnd )
        real, intent(in)    :: mean, stdev
        integer, intent(in) :: NP
        integer             :: irnd
        real                :: limits(2)
        irnd = 1
        if( NP == 0 )then
            write(logfhandle,*) 'Gaussian random integer must be generated from a non-empty (.ne. 0) set!'
            write(logfhandle,*) 'In: irnd_gasdev, module: simple_rnd'
            stop
        else if( NP == 1 )then
            irnd = 1
        else
            limits(1) = 1.
            limits(2) = real(NP)
            irnd      = max(1,nint(gasdev( mean, stdev, limits )))
        endif
    end function irnd_gasdev

    !>  \brief  pick a random point on the surface of the 4-dimensional sphere
    function rnd_4dim_sphere_pnt( ) result( rsph )
        real :: u0, u1, u2, u3, rsph(4), sca
        ! generate 4 uniform random deviates in [-1,1]
        u0      = 2.0*ran3()-1.0
        u1      = 2.0*ran3()-1.0
        u2      = 2.0*ran3()-1.0
        u3      = 2.0*ran3()-1.0
        ! the first two dimensions are deviates u0 & u1
        rsph(1) = u0
        rsph(2) = u1
        ! the next two dimensions are the remaining deviates scaled by the below
        sca     = sqrt((1.0-u0*u0-u1*u1)/(u2*u2+u3*u3))
        rsph(3) = u2*sca
        rsph(4) = u3*sca
    end function rnd_4dim_sphere_pnt

    !>  \brief  pick a random improving correlation, given the previous
    function shcloc( ncorrs, corrs, corr_prev ) result( this )
        integer, intent(in)  :: ncorrs
        real, intent(in)     :: corrs(ncorrs)
        real, intent(in)     :: corr_prev
        real    :: diffs(ncorrs)
        logical :: avail(ncorrs)
        integer :: this
        this  = 0
        diffs = corrs-corr_prev
        avail = .true.
        where( diffs < 0. ) avail = .false. ! diffs < 0. means corr_prev > all(corrs)
        if( .not. any(avail) )then          ! if no one is avalable, we return
            return
        else if( all(avail) )then
            this = irnd_uni(ncorrs)         ! if all are available, we pick a random one
            return
        else                                ! if only some are available, we pick a random of the ones available
            this = irnd_uni(ncorrs)
            do while( .not. avail(this) )
                this = irnd_uni(ncorrs)
            end do
        endif
    end function shcloc

    !>  \brief  mnorm_smp samples a multivariate normal distribution.
    !!          The multivariate normal distribution for the M dimensional vector X has the form:
    !!          pdf(X) = (2*pi*det(A))**(-M/2) * exp(-0.5*(X-MU)'*inverse(A)*(X-MU))
    !!          where MU is the mean vector, and A is a positive definite symmetric
    !!          matrix called the variance-covariance matrix. M=the dimension of the space.
    !!          N=the number of points. Input, real A(M,M), the variance-covariance
    !!          matrix.  A must be positive definite symmetric. Input, real MU(M), the mean vector.
    !!          Output, real X(M), the points.
    function mnorm_smp( cov, m, means ) result( x )
        integer, intent(in) :: m
        real, intent(in)    :: cov(m,m), means(m)
        integer             :: info, i
        real                :: r(m,m), x(m), xtmp(1,m)
        ! Compute the upper triangular Cholesky factor R of the variance-covariance matrix.
        r = cov
        call r8po_fa ( m, r, info )
        if ( info /= 0 ) then
            write ( *, '(a)' ) 'mnorm_smp - Fatal error!'
            write ( *, '(a)' ) 'The variance-covariance matrix is not positive definite symmetric'
            stop
        end if
        ! Samples of the 1D normal distribution with mean 0 and variance 1.
        do i=1,m
            x(i) = gasdev()
        end do
        ! Compute R' * X.
        xtmp(1,:) = x
        xtmp = matmul(xtmp,r)
        x = xtmp(1,:)+means
    end function mnorm_smp

    !>  \brief  R8PO_FA factors an R8PO matrix. The R8PO storage format is used for a symmetric
    !!          positive definite matrix and its inverse.  (The Cholesky factor of an R8PO matrix is an
    !!          upper triangular matrix, so it will be in R8GE storage format.) Only the diagonal and
    !!          upper triangle of the square array are used. This same storage scheme is used when the
    !!          matrix is factored by R8PO_FA, or inverted by R8PO_INVERSE.  For clarity, the lower triangle
    !!          is set to zero. R8PO storage is used by LINPACK and LAPACK. The positive definite symmetric
    !!          matrix A has a Cholesky factorization of the form:
    !!
    !!          A = R' * R
    !!
    !!          where R is an upper triangular matrix with positive elements on its diagonal. This routine
    !!          overwrites the matrix A with its factor R
    !!          Reference: Jack Dongarra, Jim Bunch, Cleve Moler, Pete Stewart,
    !!          LINPACK User's Guide, SIAM, 1979, ISBN13: 978-0-898711-72-1, LC: QA214.L56.
    subroutine r8po_fa( n, a, info )
       integer, intent(in)  :: n
       real, intent(inout)  :: a(n,n)
       integer, intent(out) :: info
       integer :: i, j, k
       real :: s
       do j = 1, n
           do k = 1, j - 1
               a(k,j) = ( a(k,j) - sum ( a(1:k-1,k) * a(1:k-1,j) ) ) / a(k,k)
           end do
           s = a(j,j) - sum ( a(1:j-1,j)**2 )
           if ( s <= 0.0D+00 ) then
               info = j
               return
           end if
           a(j,j) = sqrt(s)
       end do
       info = 0
       ! Since the Cholesky factor is stored in R8GE format, be sure to
       ! zero out the lower triangle
       do i = 1, n
           do j = 1, i-1
               a(i,j) = 0.0D+00
           end do
       end do
    end subroutine r8po_fa

end module simple_rnd
