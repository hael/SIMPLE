!>  \brief  SIMPLE comlin_corrcalc class
module simple_comlin_corrcalc
use simple_defs    ! singleton
use simple_jiffys, only: alloc_err
implicit none

public :: comlin_corrcalc
private

! CLASS PARAMETERS/VARIABLES
complex, parameter :: zero=cmplx(0.,0.) !< just a complex zero
logical, parameter :: debug=.false.     !< debug indicator

type :: comlin_corrcalc
    private
    integer              :: nptcls     = 1       !< the total number of particles in partition (logically indexded [fromp,top])
    integer              :: nrots      = 0       !< number of in-plane rotations for one pft (determined by radius of molecule)
    integer              :: ring2      = 0       !< radius of molecule
    integer              :: ldim(3)    = 0       !< logical dimensions of original cartesian image
    integer              :: kfromto(2) = 0       !< Fourier index range
    real,    allocatable :: sqsums(:,:)          !< memoized square sums for the correlation calculations
    real,    allocatable :: angtab(:)            !< table of in-plane angles (in degrees)
    real,    allocatable :: polar(:,:)           !< table of polar coordinates (in Cartesian coordinates)
    complex, allocatable :: pfts(:,:,:)          !< 3D complex matrix of polar sections
    logical              :: existence=.false.    !< to indicate existence
  contains
    ! CONSTRUCTOR
    procedure :: new
    ! GETTERS/SETTERS
    procedure :: get_nptcls
    procedure :: get_nrots
    procedure :: get_ring2
    procedure :: get_ldim
    procedure :: get_kfromto
    procedure :: get_coord
    procedure :: set_pft
    ! CHECKUPS
    procedure :: exists
    ! MEMOIZERS
    procedure :: memoize_sqsums
    ! CALCULATORS
    procedure :: corr
    ! DESTRUCTOR
    procedure :: kill
end type

contains
    
    ! CONSTRUCTORS
    
    !>  \brief  is a constructor
    subroutine new( self, nptcls, ldim, kfromto, ring2 )
        use simple_math, only: rad2deg, is_even, round2even
        class(comlin_corrcalc), intent(inout) :: self
        integer, intent(in)                   :: nptcls, ldim(3), kfromto(2), ring2
        integer :: alloc_stat, irot, k
        logical :: even_dims, test(3)
        real    :: ang
        ! kill possibly pre-existing object
        call self%kill
        ! error check
        if( kfromto(2)-kfromto(1) <= 2 )then
            write(*,*) 'kfromto: ', kfromto(1), kfromto(2)
            stop 'resolution range too narrow; new; simple_comlin_corrcalc'
        endif
        if( ring2 < 1 )then
            write(*,*) 'ring2: ', ring2
            stop 'ring2 must be > 0; new; simple_comlin_corrcalc'
        endif
        if( nptcls < 1 )then
            write(*,*) 'nptcls: ', nptcls
            stop 'nptcls (# of particle sections) must be > 0; new; simple_comlin_corrcalc'
        endif
        if( any(ldim == 0) )then
            write(*,*) 'ldim: ', ldim
            stop 'ldim is not conforming (is zero); new; simple_comlin'
        endif
        if( ldim(3) > 1 )then
            write(*,*) 'ldim: ', ldim
            stop '3D comlins are not yet supported; new; simple_comlin_corrcalc'
        endif
        test = .false.
        test(1) = is_even(ldim(1))
        test(2) = is_even(ldim(2))
        test(3) = ldim(3) == 1
        even_dims = all(test)
        if( .not. even_dims )then
            write(*,*) 'ldim: ', ldim
            stop 'only even logical dims supported; new; simple_comlin_corrcalc'
        endif
        ! set constants
        self%nptcls  = nptcls                          !< the total number of particles in partition (logically indexded [fromp,top])
        self%ring2   = ring2                           !< radius of molecule
        self%nrots   = round2even(twopi*real(ring2))/2 !< number of in-plane rotations for one pft  (determined by radius of molecule)
        self%ldim    = ldim                            !< logical dimensions of original cartesian image
        self%kfromto = kfromto                         !< Fourier index range
        ! generate polar coordinates
        allocate( self%polar(2*self%nrots,self%kfromto(1):self%kfromto(2)), self%angtab(self%nrots), stat=alloc_stat)
        call alloc_err('polar coordinate arrays; new; simple_comlin_corrcalc', alloc_stat)
        ang = twopi/real(2*self%nrots)
        do irot=1,self%nrots
            self%angtab(irot) = (irot-1)*ang
            do k=self%kfromto(1),self%kfromto(2)
                self%polar(irot,k)            = cos(self%angtab(irot))*real(k) ! x-coordinate
                self%polar(irot+self%nrots,k) = sin(self%angtab(irot))*real(k) ! y-coordinate
            end do
            self%angtab(irot) = rad2deg(self%angtab(irot)) ! angle (in degrees)
        end do
        ! allocate pfts and sqsums
        allocate( self%pfts(self%nptcls,self%nrots,self%kfromto(1):self%kfromto(2)),&
                  self%sqsums(self%nptcls,self%nrots), stat=alloc_stat)
        call alloc_err('pfts and sqsums; new; simple_comlin_corrcalc', alloc_stat)
        self%pfts      = zero
        self%sqsums    = 0.
        self%existence = .true.
    end subroutine
    
    !>  \brief  for getting the number of particles
    pure function get_nptcls( self ) result( nptcls )
        class(comlin_corrcalc), intent(in) :: self
        integer :: nptcls
        nptcls = self%nptcls
    end function
    
    !>  \brief  for getting the number of in-plane rotations
    pure function get_nrots( self ) result( nrots )
        class(comlin_corrcalc), intent(in) :: self
        integer :: nrots
        nrots = self%nrots
    end function
    
    !>  \brief  for getting the particle radius (ring2)
    function get_ring2( self ) result( ring2 )
        class(comlin_corrcalc), intent(in) :: self
        integer :: ring2
        ring2 = self%ring2
    end function
    
    !>  \brief  for getting the logical dimension of the original
    !!          Cartesian image
    function get_ldim( self ) result( ldim )
        class(comlin_corrcalc), intent(in) :: self
        integer :: ldim(3)
        ldim = self%ldim
    end function
    
    !>  \brief  for getting the Fourier index range (hp/lp)
    function get_kfromto( self ) result( lim )
        class(comlin_corrcalc), intent(inout) :: self
        integer :: lim(2)
        lim = self%kfromto
    end function

    !>  \brief returns polar coordinate for rotation rot
    !!         and Fourier index k
    function get_coord( self, rot, k ) result( xy )
        class(comlin_corrcalc), intent(in) :: self
        integer, intent(in)                :: rot, k
        real :: xy(2)
        xy(1) = self%polar(rot,k)
        xy(2) = self%polar(self%nrots+rot,k)
    end function
    
    !>  \brief  sets particle pft iptcl
    subroutine set_pft( self, iptcl, pft )
        class(comlin_corrcalc), intent(inout) :: self
        integer, intent(in) :: iptcl
        complex, intent(in) :: pft(:,:)
        self%pfts(iptcl,:,:) = pft
        ! calculate the square sum required for correlation calculation
        call self%memoize_sqsums(iptcl)
    end subroutine
    
    !>  \brief  for printing info about the object
    subroutine print( self )
        class(comlin_corrcalc), intent(in) :: self
        write(*,*) "number of particles:                            ", self%nptcls
        write(*,*) "number of rotations:                            ", self%nrots
        write(*,*) "radius of molecule:                             ", self%ring2
        write(*,*) "logical dimensions of original Cartesian image: ", self%ldim
        write(*,*) "high-pass limit Fourier index                   ", self%kfromto(1)
        write(*,*) "low-pass limit Fourier index                    ", self%kfromto(2)
    end subroutine
    
    ! CHECKUPS

    !>  \brief  checks for existence
    function exists( self ) result( yes )
        class(comlin_corrcalc), intent(in) :: self
        logical :: yes
        yes = self%existence
    end function
    
    ! MEMOIZERS

    !>  \brief  is for memoization of the complex square sums reqired for correlation calculations
    subroutine memoize_sqsums( self, iptcl )
        !$ use omp_lib
        !$ use omp_lib_kinds
        use simple_math, only: csq
        class(comlin_corrcalc), intent(inout) :: self
        integer, intent(in) :: iptcl
        integer :: irot
        !$omp parallel do default(shared) schedule(auto) private(irot)
        do irot=1,self%nrots
            self%sqsums(iptcl,irot) = sum(csq(self%pfts(iptcl,irot,:)))
        end do 
        !$omp end parallel do 
    end subroutine

    ! CALCULATOR
    
    !>  \brief  for calculating the maximum common line correlation between particles iptcl and jptcl
    function corr( self, iptcl, jptcl ) result( cc )
        !$ use omp_lib
        !$ use omp_lib_kinds
        class(comlin_corrcalc), intent(inout) :: self          !< instance
        integer, intent(in)                   :: iptcl, jptcl  !< particle indices
        real    :: corrs(self%nrots), cc
        integer :: irot
        !$omp parallel do default(shared) schedule(auto) private(irot)
        do irot=1,self%nrots
            corrs(irot) = sum(real(self%pfts(iptcl,irot,:)*conjg(self%pfts(jptcl,irot,:))))/&
            sqrt(self%sqsums(iptcl,irot)*self%sqsums(jptcl,irot))
        end do
        !$omp end parallel do
        cc = maxval(corrs)
    end function

    ! DESTRUCTOR

    !>  \brief  is a destructor
    subroutine kill( self )
        class(comlin_corrcalc), intent(inout) :: self
        if( self%existence )then
            deallocate(self%sqsums, self%angtab, self%polar, self%pfts)
            self%existence = .false.
        endif
    end subroutine
    
end module simple_comlin_corrcalc
