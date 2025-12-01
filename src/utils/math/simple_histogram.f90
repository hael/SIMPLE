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
    real,    allocatable :: counts(:)
    integer              :: nbins = 0
    real                 :: ntot  = 0.
    real                 :: dx=0., xrange=0.
    logical              :: exists = .false.
  contains
    ! Initialization
    procedure, private :: new_1, new_2, new_3, new_4, new_5
    generic            :: new => new_1, new_2, new_3, new_4, new_5
    procedure          :: reset
    procedure          :: zero
    procedure          :: quantize
    ! Getters
    procedure          :: get, get_x, get_bin, get_nbins
    ! Arithmetics
    procedure          :: add
    procedure          :: div
    procedure          :: update
    ! Descriptors
    procedure          :: npeaks
    procedure          :: find_hill, find_next_valley
    procedure, private :: moments
    procedure          :: mean
    procedure          :: variance
    procedure          :: skew
    procedure          :: hmode
    procedure          :: entropy
    ! Operations
    procedure          :: topdf
    procedure          :: smooth
    procedure          :: plot
    ! Comparison
    procedure          :: TVD
    procedure          :: KLD
    procedure          :: JSD
    procedure          :: HD
    ! Destructor
    procedure          :: kill

end type histogram

contains

    subroutine new_1( self, nbins )
        class(histogram), intent(inout) :: self
        integer,          intent(in)    :: nbins
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
        self%nbins = nbins
        if( alloc )then
            if( nbins < 2 ) THROW_HARD('Invalid NBINS value: '//int2str(nbins))
            allocate(self%x(self%nbins+1),self%counts(self%nbins))
            self%x      = 0.
            self%xrange = 0.
            self%counts = 0.
            self%dx     = 0.
            self%ntot   = 0.
        endif
        self%exists = .true.
    end subroutine new_1

    subroutine new_2( self, img, nbins, minmax, radius )
        class(histogram),  intent(inout) :: self
        class(image),      intent(in)    :: img
        integer,           intent(in)    :: nbins
        real,    optional, intent(in)    :: minmax(2), radius
        call self%new_1(nbins)
        call self%quantize(img, minmax, radius)
    end subroutine new_2

    subroutine new_3( self, other, copy )
        class(histogram),  intent(inout) :: self
        type(histogram),   intent(in)    :: other
        logical, optional, intent(in)    :: copy
        call self%new_1(other%nbins)
        self%x      = other%x
        self%xrange = other%xrange
        self%dx     = other%dx
        if( present(copy) )then
            if( copy )then
                self%counts = other%counts
                self%ntot   = other%ntot
            endif
        endif
    end subroutine new_3

    subroutine new_4( self, nbins, rvec )
        class(histogram),  intent(inout) :: self
        integer,           intent(in)    :: nbins
        real, allocatable, intent(in)    :: rvec(:)
        real    :: minmax(2), diff
        integer :: i,n,bin
        if( .not.allocated(rvec) ) THROW_HARD('Input vector not allocated')
        call self%new_1(nbins)
        n         = size(rvec)
        minmax(1) = minval(rvec)
        minmax(2) = maxval(rvec)
        diff      = minmax(2) - minmax(1)
        !if( diff < TINY ) THROW_HARD('Invalid bounds!')
        if( diff < TINY ) then
            minmax(1) = minmax(1) * 0.9
            minmax(2) = minmax(2) * 1.1
            diff      = minmax(2) - minmax(1)
        end if
        self%xrange = diff
        self%dx     = self%xrange / real(self%nbins)
        self%x(1)   = minmax(1)
        do i = 2,self%nbins
            self%x(i) = self%x(1) + real(i-1)*self%dx
        enddo
        self%x(self%nbins+1) = minmax(2)
        ! quantization
        self%counts = 0.
        do i = 1,n
            bin = self%get_bin(rvec(i))
            bin = min(self%nbins,max(1,bin))
            self%counts(bin) = self%counts(bin) + 1.
        enddo
        self%ntot = sum(self%counts)
    end subroutine new_4

    subroutine new_5( self, rvec )
        class(histogram),  intent(inout) :: self
        real, allocatable, intent(in)    :: rvec(:) ! are the desired bins centers
        real    :: minmax(2), diff
        integer :: i,nbins
        if( .not.allocated(rvec) ) THROW_HARD('Input vector not allocated')
        nbins = size(rvec)
        call self%new_1(nbins)
        minmax(1) = minval(rvec)
        minmax(2) = maxval(rvec)
        diff      = minmax(2) - minmax(1)
        if( diff < TINY ) THROW_HARD('Invalid bounds!')
        minmax(1) = minmax(1) - diff / real(self%nbins-1) / 2.0
        minmax(2) = minmax(2) + diff / real(self%nbins-1) / 2.0
        diff      = minmax(2) - minmax(1)
        self%xrange = diff
        self%dx     = self%xrange / real(self%nbins)
        self%x(1)   = minmax(1)
        do i = 2,self%nbins
            self%x(i) = self%x(1) + real(i-1)*self%dx
        enddo
        self%x(self%nbins+1) = minmax(2)
        call self%zero
    end subroutine new_5

    pure subroutine reset( self)
        class(histogram), intent(inout) :: self
        if( self%exists )then
            self%x      = 0.
            self%xrange = 0.
            self%nbins  = 0
            self%dx     = 0.
            call self%zero
        endif
    end subroutine reset

    pure subroutine zero( self)
        class(histogram), intent(inout) :: self
        if( self%exists )then
            self%counts = 0.
            self%ntot   = 0.
        endif
    end subroutine zero

    subroutine quantize( self, img, minmax, radius )
        class(histogram),  intent(inout) :: self
        class(image),      intent(in)    :: img
        real,    optional, intent(in)    :: minmax(2), radius
        real, pointer :: prmat(:,:,:)
        real          :: minmax_here(2), diff, radsq, dsq
        integer       :: dims(3), center(3), i,j,bin,djsq
        if( .not.self%exists ) THROW_HARD('Object has not been instanciated!')
        if( img%is_ft() ) THROW_HARD('Real space only!')
        dims = img%get_ldim()
        if( dims(3) /= 1 ) THROW_HARD('2D images only!')
        if( present(radius) )then
            radsq = radius**2
        else
            radsq = real(maxval(dims)**2)
        endif
        if( present(minmax) )then
            minmax_here = minmax
        else
            minmax_here = img%minmax(radius=radius)
        endif
        diff = minmax_here(2)-minmax_here(1)
        if( diff < TINY ) THROW_HARD('Invalid bounds!')
        self%xrange = diff
        self%dx     = self%xrange / real(self%nbins)
        self%x(1)   = minmax_here(1)
        do i = 2,self%nbins
            self%x(i) = self%x(1) + real(i-1)*self%dx
        enddo
        self%x(self%nbins+1) = minmax_here(2)
        center = dims/2 + 1
        self%counts = 0.
        call img%get_rmat_ptr(prmat)
        do j = 1,dims(2)
            djsq = (j-center(2))**2
            do i = 1,dims(1)
                dsq  = real(djsq + (i-center(1))**2)
                if( dsq > radsq ) cycle
                bin = self%get_bin(prmat(i,j,1))
                bin = min(self%nbins,max(1,bin))
                self%counts(bin) = self%counts(bin)+1.
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
        get = self%counts(i)
        if( present(prob) )then
            if( prob ) get = get / self%ntot
        endif
    end function get

    real function get_x( self , bin )
        class(histogram),  intent(in) :: self
        integer,           intent(in) :: bin
        if( (bin < 1) .or. (bin > self%nbins) ) THROW_HARD('Invalid index!')
        get_x = self%x(bin) + self%dx/2.
    end function get_x

    elemental integer function get_bin( self, val )
        class(histogram), intent(in) :: self
        real,             intent(in) :: val
        real :: v
        get_bin = 0
        v = val - self%x(1)
        if( v < 0. )then
            get_bin = 0
        else if( v > self%xrange )then
            get_bin = self%nbins+1
        else
            v = real(self%nbins) * v / self%xrange
            get_bin = min(self%nbins,max(1, ceiling(v)))
        endif
    end function get_bin

    elemental integer function get_nbins( self )
        class(histogram), intent(in) :: self
        get_nbins = self%nbins
    end function get_nbins

    !> Arithmetics

    subroutine add( self, other )
        class(histogram), intent(inout) :: self
        type(histogram),  intent(in)    :: other
        if( self%nbins /= other%nbins )then
            THROW_HARD('Invalid dimensions: '//int2str(self%nbins)//' vs. '//int2str(other%nbins))
        endif
        self%counts = self%counts + other%counts
        self%ntot   = self%ntot   + other%ntot
    end subroutine add

    subroutine div( self, v )
        class(histogram), intent(inout) :: self
        real,             intent(in)    :: v
        if( abs(v) > TINY )then
            self%counts = self%counts / v
            self%ntot    = sum(self%counts)
        endif
    end subroutine div

    subroutine update( self, v )
        class(histogram), intent(inout) :: self
        real,             intent(in)    :: v
        integer :: bin
        bin = self%get_bin(v)
        if( bin < 1 )          return
        if( bin > self%nbins ) return
        self%counts(bin) = self%counts(bin) + 1.
        self%ntot        = sum(self%counts)
    end subroutine update

    !> Calculators

    pure integer function npeaks( self, include_lims )
        class(histogram),  intent(in) :: self
        logical, optional, intent(in) :: include_lims
        integer :: i
        npeaks = 0
        do i = 2,self%nbins-1
            ! if( self%counts(i) == self%counts(i-1) ) cycle
            if( (self%counts(i) > self%counts(i-1)) .and. (self%counts(i) > self%counts(i+1)) )then
                npeaks = npeaks+1
            endif
        enddo
        if( present(include_lims) )then
            if(include_lims)then
                if((self%counts(1) > self%counts(2)))                     npeaks = npeaks+1
                if((self%counts(self%nbins) > self%counts(self%nbins-1))) npeaks = npeaks+1
            endif
        endif
    end function npeaks

    integer function find_hill( self, npeak )
        class(histogram), intent(in) :: self
        integer,          intent(in) :: npeak ! Nth peak
        real    :: prev_count
        integer :: i, n
        find_hill  = 0
        n          = 0
        prev_count = self%counts(1)
        do i = 2,self%nbins-1
            prev_count = min(prev_count,self%counts(i-1))
            if( (self%counts(i) > prev_count) .and. (self%counts(i) > self%counts(i+1)) )then
                n = n+1
                if( n == npeak )then
                    find_hill = i
                    exit
                endif
                prev_count = self%counts(i)
            endif
        enddo   
    end function find_hill

    pure integer function find_next_valley( self, bin )
        class(histogram), intent(in) :: self
        integer,          intent(in) :: bin
        integer :: i
        find_next_valley = 0
        do i = bin+1,self%nbins-1
            if( self%counts(i) < self%counts(i+1) )then
                find_next_valley = i
                exit
            endif
        enddo   
    end function find_next_valley

    real function moments( self, order, mean )
        class(histogram), intent(in) :: self
        integer,          intent(in) :: order
        real,             intent(in) :: mean
        real :: v
        if( order<1 .or. order>3 )then
            THROW_HARD('Unsupported moment order!')
        endif
        v = self%dx/2. - mean
        moments = sum(self%counts * (self%x(:self%nbins)+v)**order) / self%ntot
    end function moments

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

    ! adjusted Fisher-Pearson coefficient of skewness
    real function skew( self )
        class(histogram), intent(in) :: self
        real :: m, var
        m    = self%mean()
        var  = self%variance(m)
        skew = self%moments(3,m) / var**1.5
        skew = skew * sqrt(self%ntot*(self%ntot-1.)) / (self%ntot-2.)
    end function skew

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
            if( self%counts(i) < TINY ) cycle
            p = min(1., max(TINY, self%counts(i)/self%ntot))
            entropy = entropy - p * log(p) ! in nats
        enddo
        if( present(bounded) )then
            if( bounded ) entropy = entropy / log(real(self%nbins))
        endif
    end function entropy

    ! Operations

    subroutine topdf( self )
        class(histogram), intent(inout) :: self
        if( .not.self%exists ) THROW_HARD('Instance dos not exist!')
        if( self%ntot < TINY ) THROW_HARD('Empty histogram')
        self%counts = self%counts / self%ntot
    end subroutine topdf

    ! convolute with window [w, 1., w], borders assumed 0.
    subroutine smooth( self, w, nits, npeaks_target, npeaks_out )
        class(histogram), intent(inout) :: self
        real,             intent(in)    :: w
        integer,          intent(in)    :: nits, npeaks_target
        integer,          intent(inout) :: npeaks_out
        real    :: vec(0:self%nbins+1), vec2(0:self%nbins+1)
        integer :: i,it
        npeaks_out = 0
        vec(0)            = 0.
        vec(1:self%nbins) = self%counts(1:self%nbins)
        vec(self%nbins+1) = 0.
        do it = 1,nits
            ! smooth
            vec2(0)            = 0.
            vec2(1:self%nbins) = vec(1:self%nbins) + w * (vec(0:self%nbins-1) + vec(2:self%nbins+1))
            vec2(self%nbins+1) = 0.
            vec2(1:self%nbins) = vec2(1:self%nbins) / (1.+2.*w)
            ! rescale
            vec2(1:self%nbins) = vec2(1:self%nbins) * (self%ntot/ sum(vec2(1:self%nbins)))
            npeaks_out = 0
            do i = 1,self%nbins
                if( (vec2(i) > vec2(i-1)) .and. (vec2(i) > vec2(i+1)) )then
                    ! counting peaks and ignoring spurious ones
                    if( vec2(i) > 0.5 ) npeaks_out = npeaks_out+1
                endif
            enddo
            vec(0:self%nbins+1) = vec2(0:self%nbins+1)
            if( npeaks_out <= npeaks_target )exit
        enddo
        where( vec(1:self%nbins) < 0.5 ) vec(1:self%nbins) = 0.
        self%counts = vec(1:self%nbins)
        self%ntot   = sum(self%counts)
    end subroutine smooth

    subroutine plot( self, fname, abscissa )
        use CPlot2D_wrapper_module
        class(histogram),           intent(in) :: self
        character(len=*),           intent(in) :: fname
        character(len=*), optional, intent(in) :: abscissa
        type(string)                 :: title
        type(CPlot2D_type)           :: fig
        type(CDataSet_type)          :: dataSet
        character(len=LONGSTRLEN)    :: ps2pdf_cmd, fname_pdf, fname_eps
        integer                      :: i,iostat
        if( .not.self%exists ) THROW_HARD('Object not instantiated')
        fname_eps  = trim(fname)//'.eps'
        fname_pdf  = trim(fname)//'.pdf'
        call CPlot2D__new(fig, trim(fname)//C_NULL_CHAR)
        call CPlot2D__SetXAxisSize(fig, 400.d0)
        call CPlot2D__SetYAxisSize(fig, 400.d0)
        call CPlot2D__SetDrawLegend(fig, C_FALSE)
        call CPlot2D__SetFlipY(fig, C_FALSE)
        call CDataSet__new(dataSet)
        call CDataSet__SetDrawMarker(dataSet, C_FALSE)
        call CDataSet__SetDatasetColor(dataSet, 0.d0,0.d0,1.d0)
        call CDataSet_addpoint(dataSet, self%x(1), 0.)
        do i = 1,self%nbins
            call CDataSet_addpoint(dataSet, self%x(i)+self%dx/2., self%counts(i))
        end do
        call CDataSet_addpoint(dataSet, self%x(self%nbins)+self%dx, 0.)
        call CPlot2D__AddDataSet(fig, dataset)
        call CDataSet__delete(dataset)
        if( present(abscissa) )then
            title = trim(abscissa)//C_NULL_CHAR
        else
            title = 'X'//C_NULL_CHAR
        endif
        call CPlot2D__SetXAxisTitle(fig, title%to_char())
        title = 'Counts'//C_NULL_CHAR
        call CPlot2D__SetYAxisTitle(fig, title%to_char())
        call CPlot2D__OutputPostScriptPlot(fig, trim(fname_eps)//C_NULL_CHAR)
        call CPlot2D__delete(fig)
        ! conversion to PDF
        ps2pdf_cmd = 'gs -q -sDEVICE=pdfwrite -dNOPAUSE -dBATCH -dSAFER -dDEVICEWIDTHPOINTS=600 -dDEVICEHEIGHTPOINTS=600 -sOutputFile='&
            &//trim(fname_pdf)//' '//trim(fname_eps)
        call exec_cmdline(ps2pdf_cmd, suppress_errors=.true., exitstat=iostat)
        if( iostat == 0 ) call del_file(fname_eps)
    end subroutine plot
    
    ! Comparison

    ! Total Variation Distance
    real function TVD( self, other )
        class(histogram), intent(in) :: self, other
        if( self%nbins /= other%nbins ) THROW_HARD('Invalid dimensions')
        TVD = 0.5 * sum(abs(self%counts/self%ntot - other%counts/other%ntot))
    end function TVD

    ! K-L Divergence
    real function KLD( self, other )
        class(histogram), intent(in) :: self, other
        real    :: pself, pother
        integer :: i
        if( self%nbins /= other%nbins ) THROW_HARD('Invalid dimensions')
        KLD = 0.
        do i =1,self%nbins
            pself = self%counts(i) / self%ntot
            if( pself < TINY ) cycle
            pother = other%counts(i)  / other%ntot
            if( pother < TINY ) cycle
            KLD = KLD + pself *log(pself/pother)
        enddo
    end function KLD

    ! Jensen-Shannon Divergence, symmetrized KLD
    real function JSD( self, other )
        class(histogram), intent(in) :: self, other
        JSD = 0.5 *(self%KLD(other) + other%KLD(self))
    end function JSD

    ! Hellinger Distance
    real function HD( self, other )
        class(histogram), intent(in) :: self, other
        if( self%nbins /= other%nbins ) THROW_HARD('Invalid dimensions')
        HD = sum( (sqrt(self%counts/self%ntot) - sqrt(other%counts/other%ntot))**2 )
        HD = sqrt(HD/2.)
    end function HD

    ! Destructor
    elemental subroutine kill( self )
        class(histogram), intent(inout) :: self
        if( self%exists )then
            deallocate(self%x,self%counts)
            self%ntot   = 0.
            self%dx     = 0.
            self%xrange = 0.
            self%nbins  = 0.
            self%exists = .false.
        endif
    end subroutine kill

end module simple_histogram