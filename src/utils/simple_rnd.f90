! random number generation module
module simple_rnd
use simple_defs ! singleton
use simple_syslib, only: alloc_errchk
!use simple_jiffys
implicit none

private :: idum, r8po_fa
public

interface gasdev
    module procedure gasdev_1
    module procedure gasdev_2
    module procedure gasdev_3
end interface

interface randn
    module procedure randn_1
    module procedure randn_2
end interface

integer(long), save   :: idum
#include "simple_lib.f08"
contains

    !>  \brief  random seed
    subroutine seed_rnd
        integer :: i, n, clock
        integer, dimension(:), allocatable :: seed
        call random_seed(size = n)
        allocate(seed(n))
        call system_clock(count=clock)
        seed = clock + 37 * (/ (i - 1, i = 1, n) /)
        call random_seed(put = seed)
        deallocate(seed)
    end subroutine seed_rnd

    !>  \brief  wrapper for the intrinsic Fortran random number generator
    function ran3( ) result( harvest )
        real :: harvest
        call random_number(harvest)
    end function ran3

    !>  \brief  returns a random matrix
    subroutine ran3arr( harvest )
        real, intent(inout) :: harvest(:)
        call random_number(harvest)
    end subroutine ran3arr

    !>  \brief  returns a random array [-1,1]
    !! \param n dimension of required rand array
    function randn_1( n ) result( a )
        integer, intent(in) :: n
        real, allocatable   :: a(:)
        integer             :: i
        allocate( a(n), stat=alloc_stat )
        if(alloc_stat /= 0) allocchk("In: randn_1; simple_rnd")
        do i=1,n
            a(i) = -1.+2.*ran3()
        end do
    end function randn_1

    !>  \brief  returns a random matrix [-1,1]
    !! \param n1,n2 dimensions of required rand array
    function randn_2( n1, n2 ) result( a )
        integer, intent(in) :: n1, n2
        real, allocatable   :: a(:,:)
        integer             :: i, j
        allocate( a(n1,n2), stat=alloc_stat )
        if(alloc_stat /= 0) allocchk("In: randn_2; simple_rnd")
        do i=1,n1
            do j=1,n2
                a(i,j) = -1.+2.*ran3()
            end do
        end do
    end function randn_2

    !>  \brief  is for generating a Bernoulli random number _b_,
    !!          which is a random number that is 1.0 with probablility
    !!          _p_ and 0.0 otherwise. Therefore, for a uniform random
    !!          number _u_ drawn between zero and one, we take _b_ = _1_._0_
    !!          when _u_ <= _p_ and _b_ = _0_._0_ otherwise
    function bran( p ) result( b )
        real, intent(in)       :: p !< probablility
        real                   :: b
        ! generate the Bernoulli random number
        if( ran3() <= p ) then
            b = 1. ! success
        else
            b = 0. ! failure
        endif
    end function bran

    !>  \brief  generates a multinomal 1-of-K random number according to the
    !!          distribution in pvec
    function multinomal( pvec ) result( which )
        use simple_math, only: hpsort
        real,     intent(in) :: pvec(:) !< probabilities
        real,    allocatable :: pvec_sorted(:)
        integer, allocatable :: inds(:)
        integer :: i, n, which !, irnd
        real    :: rnd, bound
        n = size(pvec)
        allocate(pvec_sorted(n), source=pvec)
        inds = (/(i,i=1,n)/)
        call hpsort(n, pvec_sorted, inds)
        if( sum(pvec_sorted) >= 1.001 )then
            stop 'probability distribution does not sum up to 1.; multinomal; simple_rnd;'
        endif
        rnd = ran3()
        do which=1,n
            bound = sum(pvec_sorted(1:which))
            if( rnd <= bound )exit
        enddo
        if( which > n ) which = n ! to deal with numerical instability
        which = inds(which)
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
        harvest = stdev*gasdev( )+mean
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
        harvest = limits(2)+99.
        cnt = 0
        do while( harvest < limits(1) .or. harvest > limits(2) )
            harvest = gasdev( mean, stdev )
            cnt = cnt+1
            if( cnt == maxits )then
                write(*,*) 'WARNING! gasdev_3 exceeding maxits; simple_rnd'
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
            write(*,*) 'Uniform random integer must be generated from a non-empty set!'
            write(*,*) 'In: irnd_uni, module: simple_rnd'
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
            write(*,*) 'Gaussian random integer must be generated from a non-empty (.ne. 0) set!'
            write(*,*) 'In: irnd_gasdev, module: simple_rnd'
            stop
        else if( NP == 1 )then
            irnd = 1
        else
            limits(1) = 1.
            limits(2) = real(NP)
            irnd      = max(1,nint(gasdev( mean, stdev, limits )))
        endif
    end function irnd_gasdev

    !>  \brief  generates an array of random integers [_1_,_NP_]
    subroutine ran_iarr( iarr, NP )
        integer, intent(in)  :: NP
        integer, intent(out) :: iarr(:)
        integer :: i
        do i=1,size( iarr )
            iarr(i) = irnd_uni(NP)
        end do
    end subroutine ran_iarr

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

end module simple_rnd
