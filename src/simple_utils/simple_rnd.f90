!> Simple random number generation module
!
! simple_rnd contains routines for generation of random numbers. The code is distributed
! with the hope that it will be useful, but _WITHOUT_ _ANY_ _WARRANTY_. Redistribution
! or modification is regulated by the GNU General Public License. *Author:* Hans Elmlund, 2009-05-12.
!
!==Changes are documented below
!
!* incorporated in the _SIMPLE_ library, HE 2009-06-25
!
module simple_rnd
use simple_defs ! singleton
use simple_jiffys, only: alloc_err
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

contains

    !>  \brief  set idum to any negative value to initialize or reinitialize the sequence
    subroutine seed_rnd
        ! real :: rrnd
        call init_random_seed
        ! call random_number(rrnd)
        ! idum = -nint(rrnd*(real(3000000)))

        contains

            subroutine init_random_seed
                integer :: i, n, clock
                integer, dimension(:), allocatable :: seed
                call random_seed(size = n)
                allocate(seed(n))
                call system_clock(count=clock)
                seed = clock + 37 * (/ (i - 1, i = 1, n) /)
                call random_seed(put = seed)
                deallocate(seed)
            end subroutine

    end subroutine seed_rnd

    !! NEW: the intrinsic Fortran random number generator (because ran3() is buggy in OpenMP sections)
    !! OLD:
    !>  \brief  returns a random deviate between 0.0 and 1.0. The algorithm is developed
    !!          by Donald Knuth, fetched from numerical recepies. If suspecting that the
    !!          randomness is not sufficient using this routine, implement ran4 from NR or
    !!          go back to the standard generator
    function ran3( ) result( harvest )
        ! integer(long), parameter :: mbig=1000000000, mseed=161803398, mz=0
        ! real, parameter          :: fac=1.e-9
        ! integer(long), save      :: ma(55),inext,inextp,iff=0
        ! integer(long)            :: mj,mk,i,ii,k
        real                     :: harvest
        ! if( idum < 0 .or. iff == 0 )then ! Initialization
        !     iff    = 1
        !     mj     = mseed-iabs(idum) ! Initialize ma(55) using the seed idum and the large number mseed
        !     mj     = mod(mj,mbig)
        !     ma(55) = mj
        !     mk     = 1
        !     do i=1,54                 ! Now, initialize the rest of the table
        !         ii     = mod(21*i,55) ! in a slightly random order
        !         ma(ii) = mk           ! with numbers that are not especially random
        !         mk     = mj-mk
        !         if( mk < mz ) mk = mk+mbig
        !         mj     = ma(ii)
        !     end do
        !     do k=1,4          ! randomize by "warming up the generator"
        !         do i=1,55
        !               ma(i) = ma(i)-ma(1+mod(i+30,55))
        !               if( ma(i) < mz ) ma(i) = ma(i)+mbig
        !         end do
        !     end do
        !     inext  = 0  ! prepare indices for our first generated number
        !     inextp = 31 ! constant 31 is special, see Knuth
        !     idum   = 1
        ! endif
        ! ! Here is were it starts, except on initialization
        ! inext  = inext+1
        ! if( inext == 56 ) inext = 1
        ! inextp = inextp+1
        ! if( inextp == 56 ) inextp = 1
        ! mj = ma(inext)-ma(inextp)
        ! if( mj < mz ) mj = mj+mbig
        ! ma(inext) = mj
        ! harvest = real(mj)*fac
        call random_number(harvest)
    end function ran3

    !>  \brief  returns a random matrix
    subroutine ran3arr( harvest )
        real, intent(inout) :: harvest(:)
        !integer :: i
        call random_number(harvest)
        ! do i=1,size(harvest)
        !     harvest(i) = ran3()
        ! end do
    end subroutine ran3arr

    !>  \brief  returns a random array [-1,1]
    !! \param n dimension of required rand array
    function randn_1( n ) result( a )
        integer, intent(in) :: n
        real, allocatable   :: a(:)
        integer             :: alloc_stat, i
        allocate( a(n), stat=alloc_stat )
        call alloc_err("In: randn_1; simple_rnd", alloc_stat)
        do i=1,n
            a(i) = -1.+2.*ran3()
        end do
    end function randn_1

    !>  \brief  returns a random matrix [-1,1]
    !! \param n1,n2 dimensions of required rand array
    function randn_2( n1, n2 ) result( a )
        integer, intent(in) :: n1, n2
        real, allocatable   :: a(:,:)
        integer             :: alloc_stat, i, j
        allocate( a(n1,n2), stat=alloc_stat )
        call alloc_err("In: randn_2; simple_rnd", alloc_stat)
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
        !real                :: rrnd
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
        use simple_jiffys, only: alloc_err
        integer, intent(in)  :: ncorrs
        real, intent(in)     :: corrs(ncorrs)
        real, intent(in)     :: corr_prev
        real    :: diffs(ncorrs)
        logical :: avail(ncorrs)
        integer :: this !, alloc_stat
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
