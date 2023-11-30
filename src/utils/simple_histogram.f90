module simple_histogram
include 'simple_lib.f08'
use simple_image, only: image
implicit none

public :: histogram
private
#include "simple_local_flags.inc"

type :: histogram
    private
    real,    allocatable :: x(:)
    integer, allocatable :: counts(:)
    integer              :: nbins = 0
    integer              :: ntot  = 0
    real                 :: dx=0., xrange=0.
    logical              :: exists = .false.

  contains
    ! Initialization
    procedure, private :: new_1, new_2
    generic            :: new => new_1, new_2
    procedure          :: reset
    procedure          :: quantize
    ! Getters
    procedure          :: get
    procedure          :: get_bin
    ! Arithmetics
    procedure          :: add
    ! Descriptors
    procedure          :: npeaks
    procedure, private :: moments
    procedure, private :: sided_moments
    procedure          :: mean
    procedure          :: variance
    procedure          :: skewness
    procedure          :: kurtosis
    procedure          :: hmode
    procedure          :: entropy
    ! Operations
    procedure          :: smooth
    ! Comparison
    procedure          :: TVD
    procedure          :: KLD
    procedure          :: HD
    ! Destructor
    procedure          :: kill

end type histogram

contains

    subroutine new_1( self, nbins )
        class(histogram), intent(inout) :: self
        integer,             intent(in) :: nbins
        logical :: alloc
        alloc= .true.
        if( self%exists )then
            if( self%nbins == nbins )then
                call self%reset
                alloc = .false.
            else
                call self%kill
            endif
        endif
        if( alloc )then
            if( nbins < 2 ) THROW_HARD('Invalid NBINS value: '//int2str(nbins))
            self%nbins = nbins
            allocate(self%x(self%nbins+1),self%counts(self%nbins))
            self%x      = 0.
            self%xrange = 0.
            self%counts = 0
            self%dx     = 0.
            self%ntot   = 0
        endif
        self%exists = .true.
    end subroutine new_1

    subroutine new_2( self, img, nbins, grad, minmax, radius )
        class(histogram),  intent(inout) :: self
        class(image),      intent(in)    :: img
        integer,           intent(in)    :: nbins
        logical, optional, intent(in)    :: grad
        real,    optional, intent(in)    :: minmax(2), radius
        call self%new_1(nbins)
        call self%quantize(img, grad, minmax, radius)
    end subroutine new_2

    pure subroutine reset( self)
        class(histogram), intent(inout) :: self
        if( self%exists )then
            self%x      = 0.
            self%xrange = 0.
            self%counts = 0
            self%ntot   = 0
            self%nbins  = 0
            self%dx     = 0.
        endif
    end subroutine reset

    subroutine quantize( self, img, grad, minmax, radius )
        class(histogram),  intent(inout) :: self
        class(image),      intent(in)    :: img
        logical, optional, intent(in)    :: grad
        real,    optional, intent(in)    :: minmax(2), radius
        real, pointer :: prmat(:,:,:)
        real    :: minmax_here(2), diff, radsq, dsq
        integer :: dims(3), center(3), i,j,bin,djsq
        if( .not.self%exists ) THROW_HARD('Object has not been instanciated!')
        if( img%is_ft() ) THROW_HARD('Real space only!')
        dims = img%get_ldim()
        if( dims(3) /= 1 ) THROW_HARD('2D images only!')
        if( present(minmax) ) minmax_here = minmax
        radsq = real(maxval(dims)**2)
        if( present(radius) ) radsq = radius**2
        call img%get_rmat_ptr(prmat)
        if( present(grad) )then
            if( grad )then
                ! gradients magitude, img destroyed on output
                prmat(2:dims(1)-1,2:dims(2)-1,1) = &
                &(prmat(3:dims(1),2:dims(2)-1,1) - prmat(1:dims(1)-2,2:dims(2)-1,1))**2 +&
                &(prmat(2:dims(1)-1,3:dims(2),1) - prmat(2:dims(1)-1,1:dims(2)-2,1))**2
                prmat(2:dims(1)-1,2:dims(2)-1,1) = sqrt(prmat(2:dims(1)-1,2:dims(2)-1,1)) / 2.
                prmat(1,:,1)       = 0.
                prmat(dims(1),:,1) = 0.
                prmat(:,1,1)       = 0.
                prmat(:,dims(2),1) = 0.
                if( .not.present(minmax) )then
                    minmax_here    = img%minmax(radius=radius)
                    minmax_here(1) = 0.
                endif
            else
                minmax_here = img%minmax(radius=radius)
            endif
        else
            minmax_here = img%minmax(radius=radius)
        endif
        diff = minmax_here(2)-minmax_here(1)
        if( diff < TINY ) THROW_HARD('Invalid bounds!')
        self%xrange = diff
        self%dx = self%xrange / real(self%nbins)
        self%x(1) = minmax_here(1)
        do i = 2,self%nbins
            self%x(i) = self%x(1) + real(i-1)*self%dx
        enddo
        self%x(self%nbins+1) = minmax_here(2)
        center = dims/2 + 1
        self%counts = 0
        do j = 1,dims(2)
            djsq = (j-center(2))**2
            do i = 1,dims(1)
                dsq  = real(djsq + (i-center(1))**2)
                if( dsq > radsq ) cycle
                bin = self%get_bin(prmat(i,j,1))
                if( bin /= 0 ) self%counts(bin) = self%counts(bin)+1
            enddo
        enddo
        self%ntot = sum(self%counts)
        nullify(prmat)
    end subroutine quantize

    !> Getters

    elemental real function get( self , i, prob )
        class(histogram),  intent(in) :: self
        integer,           intent(in) :: i
        logical, optional, intent(in) :: prob
        if( (i < 1) .or. (i > self%nbins) )then
            get = 0.
            return
        endif
        if( present(prob) )then
            if( prob )then
                get = real(self%counts(i)) / real(self%ntot)
            else
                get = real(self%counts(i))
            endif
        else
            get = real(self%counts(i))
        endif
    end function get

    elemental integer function get_bin( self, val )
        class(histogram), intent(in) :: self
        real,             intent(in) :: val
        real :: v
        get_bin = 0
        v = val - self%x(1)
        if( v > -TINY .and. v < self%xrange+TINY )then
            v = real(self%nbins) * v / self%xrange
            get_bin = min(self%nbins,max(1, floor(v)+1))
        endif
    end function get_bin

    !> Arithmetics

    subroutine add( self, other )
        class(histogram), intent(inout) :: self, other
        if( self%nbins /= other%nbins ) THROW_HARD('Invalid dimensions')
        self%counts = self%counts + other%counts
        self%ntot   = self%ntot + other%ntot
    end subroutine add

    !> Calculators

    pure integer function npeaks( self, include_lims )
        class(histogram),  intent(in) :: self
        logical, optional, intent(in) :: include_lims
        integer :: i
        npeaks = 0
        do i = 2,self%nbins-1
            if( (self%counts(i) > self%counts(i-1)) .and. (self%counts(i) > self%counts(i+1)) )then
                npeaks = npeaks+1
            endif
        enddo
        if( present(include_lims) )then
            if(include_lims)then
                if((self%counts(1) > self%counts(2)))                     npeaks = npeaks+1
                if((self%counts(self%nbins) > self%counts(self%nbins+1))) npeaks = npeaks+1
            endif
        endif
    end function npeaks

    real function moments( self, order, mean )
        class(histogram), intent(in) :: self
        integer,          intent(in) :: order
        real,             intent(in) :: mean
        real :: v
        if( order<1 .or. order>4 )then
            THROW_HARD('Unsupported moment order!')
        endif
        v = self%dx/2. - mean
        moments = sum( real(self%counts) * (self%x(:self%nbins)+v)**order) / real(self%ntot)
    end function moments

    subroutine sided_moments( self, order, val, ms )
        class(histogram), intent(in)  :: self
        integer,          intent(in)  :: order
        real,             intent(in)  :: val
        real,             intent(out) :: ms(2)
        real    :: v, w
        integer :: i, n1,n2
        if( order<2 .or. order>4 )then
            THROW_HARD('Unsupported moment order!')
        endif
        v  = self%dx/2. - val
        ms = 0.
        n1 = 0
        n2 = 0
        do i = 1,self%nbins
            w = real(self%counts(i)) * (self%x(i)+v)**order
            if( (self%x(i)+self%dx/2.) < v )then
                ms(1) = ms(1) + w
                n1 = n1 + 1
            else
                ms(2) = ms(2) + w
                n2 = n2 + 1
            endif
        enddo
        if( n1 > 0 ) ms(1) = ms(1) / real(n1) 
        if( n2 > 0 ) ms(2) = ms(2) / real(n2) 
    end subroutine sided_moments

    real function mean( self )
        class(histogram), intent(in) :: self
        mean = self%moments(1,0.)
    end function mean

    real function variance( self, mean )
        class(histogram), intent(in) :: self
        real, optional,   intent(in) :: mean
        real :: m
        if( present(mean) )then
            m = mean
        else
            m = self%mean()
        endif
        variance = self%moments(2,m)
    end function variance

    ! Fisher's skewness, third moment
    real function skewness( self, mean )
        class(histogram), intent(in) :: self
        real,   optional, intent(in) :: mean
        real :: m
        if( present(mean) )then
            m = mean
        else
            m = self%mean()
        endif
        skewness = self%moments(3,m)
    end function skewness

    real function kurtosis( self, mean )
        class(histogram), intent(in) :: self
        real,   optional, intent(in) :: mean
        real :: m
        if( present(mean) )then
            m = mean
        else
            m = self%mean()
        endif
        kurtosis = self%moments(4,m)
    end function kurtosis

    pure real function hmode( self )
        class(histogram), intent(in) :: self
        integer :: bin
        bin   = maxloc(self%counts,dim=1)
        hmode = self%x(bin) + self%dx/2.
    end function hmode

    pure real function entropy( self, bounded )
        class(histogram),  intent(in) :: self
        logical, optional, intent(in) :: bounded
        real    :: p
        integer :: i
        entropy = 0.
        do i =1,self%nbins
            if( self%counts(i) == 0 ) cycle
            p = min(1., max(TINY, real(self%counts(i))/real(self%ntot)))
            entropy = entropy - p * log(p) ! in nats
        enddo
        if( present(bounded) )then
            if( bounded ) entropy = entropy / log(real(self%nbins))
        endif
    end function entropy

    ! Operations

    ! convolute with window [decay, 1., decay], borders assumed 0.
    subroutine smooth( self, w )
        class(histogram), intent(inout) :: self
        real,             intent(in)    :: w
        real    :: vec(self%nbins)
        vec                 = real(self%counts)
        vec(2:self%nbins)   = vec(2:self%nbins)   + w * real(self%counts(1:self%nbins-1))
        vec(1:self%nbins-1) = vec(1:self%nbins-1) + w * real(self%counts(2:self%nbins))
        self%counts = nint(vec / (1.+2.*w))
        self%ntot   = sum(self%counts)
    end subroutine smooth
    
    ! Comparison

    ! Total Variation Distance
    real function TVD( self, other )
        class(histogram), intent(in) :: self, other
        if( self%nbins /= other%nbins ) THROW_HARD('Invalid dimensions')
        TVD = 0.5 * sum(abs(real(self%counts)/real(self%ntot) - real(other%counts)/real(other%ntot)))
    end function TVD

    ! K-L Divergence
    real function KLD( self, other )
        class(histogram), intent(in) :: self, other
        real    :: pself, pother
        integer :: i
        if( self%nbins /= other%nbins ) THROW_HARD('Invalid dimensions')
        KLD = 0.
        do i =1,self%nbins
            pself = real(self%counts(i)) / real(self%ntot)
            if( pself < TINY ) cycle
            pother = real(other%counts(i))  / real(other%ntot)
            if( pother < TINY ) cycle
            KLD = KLD + pself *log(pself/pother)
        enddo
    end function KLD

    ! Hellinger Distance
    real function HD( self, other )
        class(histogram), intent(in) :: self, other
        if( self%nbins /= other%nbins ) THROW_HARD('Invalid dimensions')
        HD = sum( (sqrt(real(self%counts)/real(self%ntot)) - sqrt(real(other%counts)/real(other%ntot)))**2 )
        HD = sqrt(HD/2.)
    end function HD

    ! Destructor
    pure subroutine kill( self )
        class(histogram), intent(inout) :: self
        if( self%exists )then
            deallocate(self%x,self%counts)
            self%ntot   = 0
            self%dx     = 0.
            self%xrange = 0.
            self%nbins  = 0
            self%exists = .false.
        endif
    end subroutine kill

end module simple_histogram