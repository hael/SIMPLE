!>  \brief  SIMPLE polarft_corrcalc class
module simple_polarft_corrcalc
use simple_defs ! singleton
use simple_jiffys, only: alloc_err
implicit none

public :: polarft_corrcalc, test_polarft_corrcalc, test_polarft_corrcalc_inds
private

! CLASS PARAMETERS/VARIABLES
complex, parameter :: zero=cmplx(0.,0.) !< just a complex zero
logical, parameter :: debug=.false.     !< debug indicator

type :: polarft_corrcalc
    private
    integer              :: fromp     = 1      !< from particle index (in parallel execution)
    integer              :: top       = 1      !< to particle index (in parallel execution)
    integer              :: nptcls    = 1      !< the total number of particles in partition (logically indexded [fromp,top])
    integer              :: nrefs     = 1      !< the number of references (logically indexded [1,nrefs])
    integer              :: npfts     = 2      !< the total number of polarfts in the structure
    integer              :: nrots     = 0      !< number of in-plane rotations for one pft  (determined by radius of molecule)
    integer              :: nrots_ref = 0      !< total number of reference section in-plane rotations
    integer              :: nrots_tot = 0      !< total number of in-plane rotations
    integer              :: ring2     = 0      !< radius of molecule
    integer              :: refsz     = 0      !< size of reference (number of vectors used for matching)
    integer              :: ptclsz    = 0      !< size of particle (2*nrots)
    integer              :: ldim(3)   = 0      !< logical dimensions of original cartesian image
    integer              :: khp       = 2      !< high-pass limit Fourier index
    integer              :: klp       = 11     !< low-pass limit Fourier index
    real, allocatable    :: sqsums(:,:)        !< memoized square sums for the correlation calcualtions
    real, allocatable    :: angtab(:)          !< table of in-plane angles (in degrees)
    real, allocatable    :: argtransf(:,:)     !< argument transfer constants for shifting the references
    real, allocatable    :: polar(:,:)         !< table of polar coordinates (in Cartesian coordinates)
    real, allocatable    :: carte(:,:)         !< table of cartesian coordinates 
    !Partial derivative variables
    real, allocatable    :: grad_polar_psi(:,:)!< table of Gradient polar coordinates wrt psi in x0
    real, allocatable    :: grad_polar_k(:,:)  !< table of Gradient polar coordinates wrt k in x0
    real, allocatable    :: grad_argtransf_x0(:,:) !< argument transfer constants for shifting the references
    real, allocatable    :: grad_argtransf_y0(:,:) !< argument transfer constants for shifting the references
    !pft varaibles
    complex, allocatable :: pft(:,:)           !< first coord is angle (1:nrots), second is ring
    logical              :: existence=.false.  !< to indicate existence
  contains
    ! CONSTRUCTOR
    procedure :: new
    ! GETTERS/SETTERS
    procedure :: get_ptcl_range
    procedure :: get_nptcls
    procedure :: get_nrefs
    procedure :: get_nrots
    procedure :: get_ring2
    procedure :: get_klp
    procedure :: get_khp
    procedure :: get_ldim
    procedure :: get_rot
    procedure :: get_roind
    procedure :: get_coord
    procedure :: get_ref_dims
    procedure :: get_refsh_dims
    procedure :: get_ptcl_dims
    procedure :: get_ptcl_pft
    procedure :: get_ref_pft
    procedure :: set_fcomp
    procedure :: set_ref_pft
    procedure :: set_ptcl_pft
    procedure :: print
    ! PRIVATE ROUTINES FOR DEALING WITH MATRIX INDICES
    procedure, private :: rlim      ! in-plane rotational range for reference with logical index i .in. [1,nrefs]
    procedure, private :: rlimsh    ! in-plane rotational range for shifted reference with logical index i .in. [1,nrefs]
    procedure, private :: plim      ! in-plane rotational range for particle with logical index j .in. [fromp,top] and rotation r .in. [1,nrots]
    procedure, private :: pindsqsum ! 2D index for memoized sqsums of particle with logical index j .in. [fromp,top]
    ! CHECKUPS
    procedure :: exists
    procedure, private :: same_dims
    generic :: operator(.eqdims.) => same_dims
    ! MEMOIZERS
    procedure :: memoize_ref_sqsums
    procedure, private :: memoize_ptcl_sqsums
    ! MODIFIERS
    procedure, private :: grad_polar_x0
    procedure, private :: grad_polar_y0
    procedure, private :: shift
    procedure, private :: pdiff_shift_polar_x0
    procedure, private :: pdiff_shift_polar_y0
    procedure, private :: pdiff_shift_psi_x0
    procedure, private :: pdiff_shift_psi_y0
    ! CALCULATORS
    procedure :: grid_srch
    procedure :: gencorrs
    procedure :: corr
    procedure :: pdiff_corr_polar_x0  !partial derivatives for the corr in x0
    procedure :: pdiff_corr_polar_y0  !partial derivatives for the corr in y0
    procedure :: pdiff_corr_psi_x0    !partial derivatives for the corr in x0
    procedure :: pdiff_corr_psi_y0    !partial derivatives for the corr in y0
    procedure, private :: calc_corr
    procedure, private :: calc_corr_shifted
    ! DESTRUCTOR
    procedure :: kill
end type

interface polarft_corrcalc
    module procedure constructor
end interface

contains
    
    ! CONSTRUCTORS
    
    !>  \brief  is a constructor
    function constructor( nrefs, prange, ldim, kfromto, ring2 ) result( self )
        integer, intent(in)    :: nrefs, prange(2), ldim(3), kfromto(2), ring2
        type(polarft_corrcalc) :: self
        call self%new(nrefs, prange, ldim, kfromto, ring2)
    end function

    !>  \brief  is a constructor
    subroutine new( self, nrefs, prange, ldim, kfromto, ring2 )
        use simple_math, only: rad2deg, is_even, round2even
        class(polarft_corrcalc), intent(inout) :: self
        integer, intent(in) :: nrefs, prange(2), ldim(3), kfromto(2), ring2
        integer :: alloc_stat, r, i, j
        logical :: even_dims, test(3)
        real    :: ang
        ! kill possibly pre-existing object
        call self%kill
        ! error check
        if( kfromto(2)-kfromto(1) <= 2 )then
            write(*,*) 'fromto:', kfromto(1), kfromto(2)
            stop 'resolution range too narrow; new; simple_polarft_corrcalc'
        endif
        if( ring2 < 1 )      stop &
             'ring2 must be > 0; new; simple_polarft_corrcalc'
        if( prange(2)-prange(1)+1 < 1 ) stop &
        'nptcls (# of particles) must be > 0; new; simple_polarft_corrcalc'
        if( nrefs < 1 )      stop &
             'nrefs (# of refrc sectns) must be > 0;new;simple_polarft_corrcalc'
        if( any(ldim == 0) ) stop &
             'ldim is not conforming (is zero); new; simple_polarft'
        if( ldim(3) > 1 )    stop &
             '3D polarfts are not yet supported; new; simple_polarft_corrcalc'
        test = .false.
        test(1) = is_even(ldim(1))
        test(2) = is_even(ldim(2))
        test(3) = ldim(3) == 1
        even_dims = all(test)
        ! set constants
        if( .not. even_dims ) stop &
             'only even logical dims supported; new; simple_polarft_corrcalc'
        self%fromp      = prange(1)                     ! from particle index (in parallel execution)
        self%top        = prange(2)                     ! to particle index (in parallel execution)
        self%nptcls     = self%top-self%fromp+1         ! the total number of particles in partition (logically indexded [fromp,top])
        self%nrefs      = nrefs                         ! the number of references (logically indexded [1,nrefs])
        self%npfts      = self%nptcls+self%nrefs        ! the total number of polarfts in the structure
        self%ring2      = ring2                         ! radius of molecule
        self%ldim       = ldim                          ! logical dimensions of original Cartesian image
        self%khp        = kfromto(1)                    ! high-pass limit Fourier index
        self%klp        = kfromto(2)                    ! low-pass limit Fourier index
        self%nrots      = round2even(twopi*real(ring2)) ! number of in-plane rotations for one pft (determined by radius of molecule)
        self%refsz      = self%nrots/2                  ! size (number of rotations) of a referece pft
        self%ptclsz     = self%nrots*2                  ! size (number of rotations) of a particle pft
        self%nrots_ref  = self%nrefs*self%nrots*2       ! total number of reference section in-plane rotations
        self%nrots_tot  = self%nrots_ref+self%nptcls*self%ptclsz ! total number of in-plane rotations in polarft_corrcalc object             
        ! generate polar coordinates
        allocate( self%polar(2*self%nrots,self%khp:self%klp), &
             self%angtab(self%nrots), stat=alloc_stat)
        call alloc_err('when trying to allocate polar coordinate arrays; new; simple_polarft_corrcalc', alloc_stat)
        allocate( self%grad_polar_psi(2*self%nrots,self%khp:self%klp), &
             stat=alloc_stat) !alloctng grad polar wrt psi
        call alloc_err('when trying to allocate polar coordinate arrays; new; simple_polarft_corrcalc', alloc_stat)
        !TODO: remove the memory allocation wrt k after test
        allocate( self%grad_polar_k(2*self%nrots,self%khp:self%klp), &
             stat=alloc_stat) !alloctng grad polar wrt k
        call alloc_err('when trying to allocate polar coordinate arrays; new; simple_polarft_corrcalc', alloc_stat)
        ang = twopi/real(self%nrots)
        do i=1,self%nrots
           self%angtab(i) = (i-1)*ang
           do j=self%khp,self%klp
              self%polar(i,j)            = cos(self%angtab(i))*real(j) ! x-coordinate
              self%polar(i+self%nrots,j) = sin(self%angtab(i))*real(j) ! y-coordinate
           end do
           self%angtab(i) = rad2deg(self%angtab(i)) ! angle (in degrees)
        end do
        ! generate the argument transfer constants for shifting reference polarfts 
        allocate( self%argtransf(self%nrots,self%khp:self%klp), stat=alloc_stat)
        call alloc_err('when trying to allocate argument transfer array; new; simple_polarft_corrcalc', alloc_stat)

        self%argtransf(:self%nrots/2,:)   = &
            self%polar(:self%nrots/2,:)   * &
            (PI/real(self%ldim(1)/2))    ! x-part
        self%argtransf(self%nrots/2+1:,:) = &
            self%polar(self%nrots+1:self%nrots+self%nrots/2,:) * &
            (PI/real(self%ldim(2)/2))    ! y-part

        !generating the gradient polar coordinates in both x0 and y0

        call self%grad_polar_x0()
        call self%grad_polar_y0()

        allocate( self%grad_argtransf_x0(self%nrots,self%khp:self%klp), &
             stat=alloc_stat)
        call alloc_err('when trying to allocate argument transfer array; new; simple_polarft_corrcalc', alloc_stat)
        allocate( self%grad_argtransf_y0(self%nrots,self%khp:self%klp),   &
             stat=alloc_stat)
        call alloc_err('when trying to allocate argument transfer array; new; simple_polarft_corrcalc', alloc_stat)
! NEEDED TO COMMENT THIS OUT BEACUSE I GOT AN OUT OF BOUNDS BUG
! At line 214 of file src/simple/simple_polarft_corrcalc.f90
! Fortran runtime error: Array bound mismatch for dimension 1 of array 'self' (171/170)
! bash-3.2$ simple_prime2 stk=projs.spi vol1=startvol.spi smpd=2 msk=55 lp=20 oritab=rndoris_discrete.txt
!         self%grad_argtransf_x0(:self%nrots/2,:)   = &
!            self%grad_polar_psi(:self%nrots/2,:)   * &
!            (PI/real(self%ldim(1)/2.))               ! x-part for x0
!         !TODO: remove the memory allocation wrt k after test
!         self%grad_argtransf_x0(self%nrots/2+1:,:) = &
!              self%grad_polar_k(:self%nrots/2,:)   * &
!              (PI/real(self%ldim(1)/2.))             ! y-part for x0
!
!         self%grad_argtransf_y0(self%nrots/2:,:)   = &
!            self%grad_polar_psi(self%nrots+1:self%nrots+self%nrots/2,:) * &
!            (PI/real(self%ldim(2)/2.))               ! x-part for y0
!         !TODO: remove the memory allocation wrt k after test
!         self%grad_argtransf_y0(self%nrots/2+1:,:) = &
!              self%grad_polar_k(self%nrots+1:self%nrots+self%nrots/2,:) * &
!              (PI/real(self%ldim(2)/2.))             ! y-part for y0
  
        ! allocate polarft and sqsums
        allocate( self%pft(self%nrots_tot,self%khp:self%klp), &
                  self%sqsums(self%npfts,4), stat=alloc_stat)
        call alloc_err('when trying to allocate polarfts & sqsums; new; simple_polarft_corrcalc', alloc_stat)
        ! for references the first half 1:nrots/2 represents the unshifted
        ! polarft whereas the second half nrots/2+1:nrots represents the 
        ! shifted version. For reference polarfts we need four square sums 
        ! self%sqsums(1) represents the sum for the unshifted part (1:nrots/2),
        ! self%sqsums(2) represents the sum for the shifted part 
        ! (nrots/2+1:nrots), In particle pfts, the first half 1:nrots is
        ! identical to the second half  nrots+1:2*nrots. There is only one
        ! sqsums. The particle pfts are copied to remove the time-consuming
        ! and difficult to handle on GPU index modifications
        self%pft       = zero
        self%sqsums    = 0.
        self%existence = .true.
    end subroutine
    
    ! GETTERS/SETTERS
    
    !>  \brief  for getting the logical particle range
    function get_ptcl_range( self ) result( lim )
        class(polarft_corrcalc), intent(inout) :: self
        integer :: lim(2)
        lim = [self%fromp,self%top]
    end function

    !>  \brief  for getting the number of particles
    function get_nptcls( self ) result( nptcls )
        class(polarft_corrcalc), intent(in) :: self
        integer :: nptcls
        nptcls = self%nptcls
    end function

    !>  \brief  for getting the number of references
    function get_nrefs( self ) result( nrefs )
        class(polarft_corrcalc), intent(in) :: self
        integer :: nrefs
        nrefs = self%nrefs
    end function

    !>  \brief  for getting the number of in-plane rotations
    function get_nrots( self ) result( nrots )
        class(polarft_corrcalc), intent(in) :: self
        integer :: nrots
        nrots = self%nrots
    end function

    !>  \brief  for getting the particle radius (ring2)
    function get_ring2( self ) result( ring2 )
        class(polarft_corrcalc), intent(in) :: self
        integer :: ring2
        ring2 = self%ring2
    end function
    
    !>  \brief  for getting the Fourier index resolution range
    function get_khp( self ) result( khp )
        class(polarft_corrcalc), intent(in) :: self
        integer :: khp
        khp = self%khp
    end function
    
    !>  \brief  for getting the Fourier index resolution range
    function get_klp( self ) result( klp )
        class(polarft_corrcalc), intent(in) :: self
        integer :: klp
        klp = self%klp
    end function

    !>  \brief  for getting the Fourier index resolution range
    function get_resrange( self ) result( kfromto )
        class(polarft_corrcalc), intent(in) :: self
        integer :: kfromto(2)
        kfromto = [self%khp,self%klp]
    end function

    !>  \brief  for getting the logical dimension of the original
    !!          Cartesian image
    function get_ldim( self ) result( ldim )
        class(polarft_corrcalc), intent(in) :: self
        integer :: ldim(3)
        ldim = self%ldim
    end function

    !>  \brief is for getting the continuous in-plane rotation
    !!         corresponding to in-plane rotation index roind
    function get_rot( self, ind ) result( rot )
        class(polarft_corrcalc), intent(in) :: self
        integer, intent(in) :: ind
        real :: rot
        if( ind < 1 .or. ind > self%nrots )then
            stop 'i is out of range; get_rot; simple_polarft'
        endif
        rot = self%angtab(ind)
    end function

    !>  \brief is for getting the discrete in-plane rotational
    !!         index corresponding to continuous rotation rot
    function get_roind( self, rot ) result( ind )
        class(polarft_corrcalc), intent(in) :: self
        real, intent(in) :: rot
        integer :: ind, i, alloc_stat, loc(1)
        real    :: dists(self%nrots)
        do i=1,self%nrots
            dists(i) = sqrt((self%angtab(i)-rot)**2.)
        end do
        loc = minloc(dists)
        ind = loc(1)
    end function

    !>  \brief returns polar coordinate for rotation rot
    !!         and Fourier index find
    function get_coord( self, rot, find ) result( xy )
        class(polarft_corrcalc), intent(in) :: self
        integer, intent(in)        :: rot, find
        real                       :: xy(2)
        xy(1) = self%polar(rot,find)
        xy(2) = self%polar(self%nrots+rot,find)
    end function

    !>  \brief  returns the dimensions of reference i
    function get_ref_dims( self, i ) result( dims )
        class(polarft_corrcalc), intent(in) :: self
        integer, intent(in) :: i
        integer :: dims(4)
        dims(1:2) = self%rlim(i)
        dims(3:4) = [self%khp,self%klp]
    end function

    !>  \brief  returns the dimensions of reference i
    function get_refsh_dims( self, i ) result( dims )
        class(polarft_corrcalc), intent(in) :: self
        integer, intent(in) :: i
        integer :: dims(4)
        dims(1:2) = self%rlimsh(i)
        dims(3:4) = [self%khp,self%klp]
    end function

    !>  \brief  returns the dimensions of reference i
    function get_ptcl_dims( self, i, r ) result( dims )
        class(polarft_corrcalc), intent(in) :: self
        integer, intent(in) :: i, r
        integer :: dims(4)
        dims(1:2) = self%plim(i,r)
        dims(3:4) = [self%khp,self%klp]
    end function
    
    !>  \brief  returns polar Fourier transform of particle j
    function get_ptcl_pft( self, j, r ) result( pft )
        class(polarft_corrcalc), intent(in) :: self
        integer, intent(in)  :: j, r
        complex, allocatable :: pft(:,:)
        integer :: lims(2), alloc_stat
        allocate(pft(self%refsz,self%khp:self%klp), stat=alloc_stat)
        call alloc_err("In: get_ptcl_pft; simple_polarft_corrcalc", alloc_stat)
        lims = self%plim(j,r)
        pft = self%pft(lims(1):lims(2),:)
    end function
    
    !>  \brief  returns polar Fourier transform of reference i
    function get_ref_pft( self, i ) result( pft )
        class(polarft_corrcalc), intent(in) :: self
        integer, intent(in)  :: i
        complex, allocatable :: pft(:,:)
        integer :: lims(2), alloc_stat
        allocate(pft(self%refsz,self%khp:self%klp), stat=alloc_stat)
        call alloc_err("In: get_ref_pft; simple_polarft_corrcalc", alloc_stat)
        lims = self%rlim(i)
        pft = self%pft(lims(1):lims(2),:)
    end function

    !>  \brief  sets Fourier component (i,j)
    subroutine set_fcomp( self, j, k, comp )
        class(polarft_corrcalc), intent(inout) :: self
        integer, intent(in)                    :: j, k
        complex, intent(in)                    :: comp
        self%pft(j,k) = comp
    end subroutine

    !>  \brief  sets reference pft i
    subroutine set_ref_pft( self, i, pft )
        class(polarft_corrcalc), intent(inout) :: self
        integer, intent(in) :: i
        complex, intent(in) :: pft(:,:)
        integer :: lim(2)
        lim = self%rlim(i)
        self%pft(lim(1):lim(2),:) = pft(1:self%refsz,:)
        ! calculate the square sums required for correlation calculation
        call self%memoize_ref_sqsums(i)
    end subroutine

    !>  \brief  sets particle pft j
    subroutine set_ptcl_pft( self, j, pft )
        class(polarft_corrcalc), intent(inout) :: self
        integer, intent(in) :: j
        complex, intent(in) :: pft(:,:)
        integer :: lim1(2), lim2(2)
        lim1    = self%plim(j,1)
        lim1(2) = lim1(1)+self%nrots-1
        lim2(1) = lim1(2)+1
        lim2(2) = lim2(1)+self%nrots-1
        self%pft(lim1(1):lim1(2),:) = pft
        self%pft(lim2(1):lim2(2),:) = pft
        ! calculate the square sums required for correlation calculation
        call self%memoize_ptcl_sqsums(j)
    end subroutine
    
    !>  \brief  for printing info about the object
    subroutine print( self )
        class(polarft_corrcalc), intent(in) :: self
        print *, "from particle index: ", self%fromp
        print *, "to particle index", self%top
        print *, "total number of particles in partition", self%nptcls
        print *, "number of references", self%nrefs
        print *, "total number of polarfts in the structure", self%npfts
        print *, "radius of molecule", self%ring2
        print *, "logical dimensions of original Cartesian image", self%ldim
        print *, "high-pass limit Fourier index", self%khp
        print *, "low-pass limit Fourier index", self%klp
        print *, "number of in-plane rotations for one pft", self%nrots
        print *, "number of rotations of a referece pft", self%refsz
        print *, "number of rotations of a particle pft", self%ptclsz
        print *, "total number of reference section in-plane rotations", self%nrots_ref
        print *, "total number of in-plane rotations in polarft_corrcalc object", self%nrots_tot
    end subroutine

    ! PRIVATE ROUTINES FOR DEALING WITH MATRIX INDICES

    !>  \brief  in-plane rotational range for reference
    !!          with logical index i .in. [1,nrefs]
    function rlim( self, i ) result( lim )
        class(polarft_corrcalc), intent(in) :: self
        integer, intent(in) :: i
        integer :: lim(2)
        lim(1) = (i-1)*self%nrots+1
        lim(2) = lim(1)+self%refsz-1
    end function

    !>  \brief  in-plane rotational range for shifted reference
    !!          with logical index i .in. [1,nrefs]
    function rlimsh( self, i ) result( lim )
         class(polarft_corrcalc), intent(in) :: self
        integer, intent(in) :: i
        integer :: lim(2)
        lim(1) = (i-1)*self%nrots+self%refsz+1
        lim(2) = lim(1)+self%refsz-1
    end function

    !>  \brief  in-plane rotational range for particle with
    !!          logical index j .in. [fromp,top] and rotation r .in. [1,nrots]
    function plim( self, j, r ) result( lim )
        class(polarft_corrcalc), intent(in) :: self
        integer, intent(in) :: j, r
        integer :: lim(2)
        lim(1) = self%nrots_ref+(j-self%fromp)*self%ptclsz+r
        lim(2) = lim(1)+self%refsz-1
    end function

    !>  \brief  2D index for memoized sqsums of particle
    !!          with logical index j .in. [fromp,top]
    function pindsqsum( self, j ) result( ind )
         class(polarft_corrcalc), intent(in) :: self
        integer, intent(in) :: j
        integer :: ind
        ind = self%nrefs+j-self%fromp+1
    end function

    ! CHECKUPS

    !>  \brief  checks for existence
    function exists( self ) result( yes )
        class(polarft_corrcalc), intent(in) :: self
        logical :: yes
        yes = self%existence
    end function

    !>  \brief  checks for same dimensions, overloaded as (.eqdims.)
    function same_dims( self1, self2 ) result( yep )
        class(polarft_corrcalc), intent(in) :: self1
        class(polarft_corrcalc), intent(in) :: self2
        logical :: yep, test(4)
        test = .false.
        test(1) = self1%nrots == self2%nrots
        test(2) = self1%khp == self2%khp
        test(3) = self1%klp == self2%klp
        test(4) = all(self1%ldim == self2%ldim)
        yep = all(test)
    end function

    ! MEMOIZERS

    !>  \brief  is for memoization of the complex square sums reqired for correlation calculation
    subroutine memoize_ref_sqsums( self, i )
        use simple_math, only: csq
        class(polarft_corrcalc), intent(inout) :: self
        integer, intent(in) :: i
        integer :: lim(2)
        !!!!!!!!!!!<<<GPU KERNEL START>>>!!!!!!!!!!!
        lim = self%rlim(i)
        self%sqsums(i,1) = sum(csq(self%pft(lim(1):lim(2),:)))
        !!!!!!!!!!!<<<GPU KERNEL END>>>!!!!!!!!!!!
    end subroutine

    !>  \brief  is for memoization of the complex square sums reqired for correlation calculation
    subroutine memoize_ptcl_sqsums( self, j )
        use simple_math, only: csq
        use simple_stat, only: moment
        class(polarft_corrcalc), intent(inout) :: self
        integer, intent(in) :: j
        integer :: lim(2), ind
        !!!!!!!!!!!<<<GPU KERNEL START>>>!!!!!!!!!!!
        ind = self%pindsqsum(j)
        lim = self%plim(j,1)
        self%sqsums(ind,1) = sum(csq(self%pft(lim(1):lim(2),:)))
        !!!!!!!!!!!<<<GPU KERNEL END>>>!!!!!!!!!!!
    end subroutine

    ! MODIFIERS

    !TODO: remove the memory allocation wrt k after test
    !> \brief Partial derivatives for the shifted transfer matrix in direction 1
    subroutine grad_polar_x0(self)
      use simple_math, only: rad2deg
      class(polarft_corrcalc), intent(inout) :: self !< instance
      !local variables
      integer :: i,j
      real :: ang

      ang = twopi/real(self%nrots)
      do i=1,self%nrots
         self%angtab(i) = (i-1)*ang
         do j=self%khp,self%klp
            self%grad_polar_psi(i,j) = - sin(self%angtab(i))*real(j) ! x-coordinate
            self%grad_polar_k(i,j)   =   cos(self%angtab(i))         ! x-coordinate
         end do
         self%angtab(i) = rad2deg(self%angtab(i)) ! angle (in degrees)
      end do

      return
    end subroutine grad_polar_x0

    !TODO: remove the memory allocation wrt k after test
    !> \brief Partial derivatives for the shifted transfer matrix in direction 2
    subroutine grad_polar_y0(self)
      use simple_math, only: rad2deg
      class(polarft_corrcalc), intent(inout) :: self !< instance
      !local variables
      integer :: i,j
      real :: ang
      ang = twopi/real(self%nrots)
      do i=1,self%nrots
         self%angtab(i) = (i-1)*ang
         do j=self%khp,self%klp
            self%grad_polar_psi(i+self%nrots,j) = cos(self%angtab(i))*real(j) ! y-coordinate
            self%grad_polar_k(i+self%nrots,j)   = sin(self%angtab(i))         ! y-coordinate
         end do
         self%angtab(i) = rad2deg(self%angtab(i)) ! angle (in degrees)
      end do
      return
    end subroutine grad_polar_y0

    !>  \brief  is for gradient shifting in x0 a reference 
    subroutine pdiff_shift_psi_x0(self, i, shvec)
      use simple_math, only: csq
      class(polarft_corrcalc), intent(inout) :: self !< instance
      integer, intent(in) :: i                       !< reference to shift
      real, intent(in)    :: shvec(2)                !< origin shift vector
      real    :: grad_argmat_x0(self%nrots/2,self%khp:self%klp)
      complex :: grad_shmat_x0(self%nrots/2,self%khp:self%klp)
      integer :: j, lim(2), limsh(2)

      ! generate the argument matrix from memoized components in grad_argtransf
      grad_argmat_x0 = self%grad_argtransf_x0(:self%nrots/2,:) * shvec(1)
      ! + self%grad_argtransf_x0(self%nrots/2+1:,:) * shvec(2)

      ! generate the complex shift transformation matrix
      grad_shmat_x0 = cmplx(cos(grad_argmat_x0),sin(grad_argmat_x0))
      ! get the range for the i:th reference...
      lim = self%rlim(i)
      !...and its shifted version
      limsh = self%rlimsh(i)
      ! shift
      self%pft(limsh(1):limsh(2),:) = self%pft(lim(1):lim(2),:) * &
                                      grad_shmat_x0
      ! memoize the grad shifted square sums in psi
      self%sqsums(i,2) = sum(csq(self%pft(limsh(1):limsh(2),:))) ! self%sqsums(i,2) because shifted ref

      return
    end subroutine pdiff_shift_psi_x0

    !>  \brief  is for gradient shifting in y0 a reference 
    subroutine pdiff_shift_psi_y0(self, i, shvec)
      use simple_math, only: csq
      class(polarft_corrcalc), intent(inout) :: self !< instance
      integer, intent(in) :: i                       !< reference to shift
      real, intent(in)    :: shvec(2)                !< origin shift vector
      real    :: grad_argmat_y0(self%nrots/2,self%khp:self%klp)
      complex :: grad_shmat_y0(self%nrots/2,self%khp:self%klp)
      integer :: j, lim(2), limsh(2)

      ! generate the argument matrix from memoized components in grad_argtransf
      grad_argmat_y0 = self%grad_argtransf_y0(:self%nrots/2,:) * shvec(1)
      ! + self%grad_argtransf_y0(self%nrots/2+1:,:) * shvec(2)

      ! generate the complex shift transformation matrix
      grad_shmat_y0 = cmplx(cos(grad_argmat_y0),sin(grad_argmat_y0))
      ! get the range for the i:th reference...
      lim = self%rlim(i)
      !...and its shifted version
      limsh = self%rlimsh(i)
      ! shift
      self%pft(limsh(1):limsh(2),:) = self%pft(lim(1):lim(2),:) * &
                                      grad_shmat_y0
      ! memoize the grad shifted square sums in k
      self%sqsums(i,2) = sum(csq(self%pft(limsh(1):limsh(2),:))) ! self%sqsums(i,2) because shifted ref

      return
    end subroutine pdiff_shift_psi_y0

    !>  \brief  is for gradient shifting in x0 a reference 
    subroutine pdiff_shift_polar_x0(self, i, shvec)
      use simple_math, only: csq
      class(polarft_corrcalc), intent(inout) :: self !< instance
      integer, intent(in) :: i                       !< reference to shift
      real, intent(in)    :: shvec(2)                !< origin shift vector
      real    :: argmat(self%nrots/2,self%khp:self%klp)
      complex :: grad_shmat_x0(self%nrots/2,self%khp:self%klp)
      integer :: j, lim(2), limsh(2)

      ! generate the argument matrix from memoized components in argtransf
      argmat = self%argtransf(:self%nrots/2,:)   * shvec(1) + &
               self%argtransf(self%nrots/2+1:,:) * shvec(2)
      ! generate the complex shift transformation matrix
      grad_shmat_x0 = cmplx( - self%argtransf(:self%nrots/2,:) * sin(argmat),  &
                               self%argtransf(:self%nrots/2,:) * cos(argmat)   )
      ! get the range for the i:th reference...
      lim = self%rlim(i)
      !...and its shifted version
      limsh = self%rlimsh(i)
      ! shift
      self%pft(limsh(1):limsh(2),:) = self%pft(lim(1):lim(2),:) * &
                                      grad_shmat_x0
      ! memoize the grad shifted square sums in psi
      self%sqsums(i,2) = sum(csq(self%pft(limsh(1):limsh(2),:))) ! self%sqsums(i,2) because shifted ref

      return
    end subroutine pdiff_shift_polar_x0

    !>  \brief is for gradient shifting in y0 a reference 
    subroutine pdiff_shift_polar_y0(self, i, shvec)
      use simple_math, only: csq
      class(polarft_corrcalc), intent(inout) :: self !< instance
      integer, intent(in) :: i                       !< reference to shift
      real, intent(in)    :: shvec(2)                !< origin shift vector
      real    :: argmat(self%nrots/2,self%khp:self%klp)
      complex :: grad_shmat_y0(self%nrots/2,self%khp:self%klp)
      integer :: j, lim(2), limsh(2)

      ! generate the argument matrix from memoized components in argtransf
      argmat = self%argtransf(:self%nrots/2,:)   * shvec(1) + &
               self%argtransf(self%nrots/2+1:,:) * shvec(2)
      ! generate the complex shift transformation matrix
      grad_shmat_y0 = cmplx( - self%argtransf(self%nrots/2+1:,:) * sin(argmat),  &
                               self%argtransf(self%nrots/2+1:,:) * cos(argmat)   )
      ! get the range for the i:th reference...
      lim = self%rlim(i)
      !...and its shifted version
      limsh = self%rlimsh(i)
      ! shift
      self%pft(limsh(1):limsh(2),:) = self%pft(lim(1):lim(2),:) * &
                                      grad_shmat_y0
      ! memoize the grad shifted square sums in k
      self%sqsums(i,2) = sum(csq(self%pft(limsh(1):limsh(2),:))) ! self%sqsums(i,2) because shifted ref

      return
    end subroutine pdiff_shift_polar_y0

    !>  \brief  is for shifting a reference
    subroutine shift( self, i, shvec )
      use simple_math, only: csq
      class(polarft_corrcalc), intent(inout) :: self !< instance
      integer, intent(in) :: i                       !< reference to shift
      real, intent(in)    :: shvec(2)                !< origin shift vector
      real    :: argmat(self%nrots/2,self%khp:self%klp)
      complex :: shmat(self%nrots/2,self%khp:self%klp)
      integer :: j, lim(2), limsh(2)
      ! get the range for the i:th reference...
      lim = self%rlim(i)
      !...and its shifted version
      limsh = self%rlimsh(i)
      ! generate the argument matrix from memoized components in argtransf
      argmat = self%argtransf(:self%nrots/2,:)*shvec(1)+self%argtransf(self%nrots/2+1:,:)*shvec(2)
      ! generate the complex shift transformation matrix
      shmat = cmplx(cos(argmat),sin(argmat))
      ! shift
      self%pft(limsh(1):limsh(2),:) = self%pft(lim(1):lim(2),:)*shmat
      ! memoize the shifted square sums
      self%sqsums(i,2) = sum(csq(self%pft(limsh(1):limsh(2),:))) ! self%sqsums(i,2) because shifted ref
    end subroutine shift

    ! CALCULATORS

    !>  \brief  for generating all in-plane correlations between reference i and particle j
    
    !!!! THE SHIFTING FUNCTIONALITY IS NEVER USED SO IN THE NEW IMPLEMENTATION IT SHOULD BE OMITTED FOR CLARITY !!!
    
    function gencorrs( self, i, j, shvec ) result( cc )
        !$ use omp_lib
        !$ use omp_lib_kinds
        class(polarft_corrcalc), intent(inout) :: self     !< instance
        integer, intent(in)                    :: i, j     !< i=reference, j=particle
        real, intent(in), optional             :: shvec(2) !< origin shift vector
        real    :: cc(self%nrots)
        integer :: r
        if( present(shvec) )then
            !!!!!!!!!!!<<<GPU KERNEL START>>>!!!!!!!!!!!
            call self%shift(i, shvec)
            !$omp parallel do default(shared) private(r)
            do r=1,self%nrots
                cc(r) = self%calc_corr_shifted(i, j, r)
            end do
            !$omp end parallel do
            !!!!!!!!!!!<<<GPU KERNEL END>>>!!!!!!!!!!!
        else
            !!!!!!!!!!!<<<GPU KERNEL START>>>!!!!!!!!!!!
            !$omp parallel do default(shared) private(r)
            do r=1,self%nrots
                cc(r) = self%corr(i, j, r)
            end do
            !$omp end parallel do
            !!!!!!!!!!!<<<GPU KERNEL END>>>!!!!!!!!!!!
        endif
    end function

    !>  \brief  for calculating the correlation between reference i and particle j in rotation r
    function pdiff_corr_psi_x0( self, i, j, r, shvec ) result( pdiff_cc_psi_x0 )
      class(polarft_corrcalc), intent(inout) :: self     !< instance
      integer, intent(in)                    :: i, j, r  !< i=reference, j=particle, r=rotation
      real, intent(in), optional             :: shvec(2) !< origin shift vector
      real :: pdiff_cc_psi_x0
      integer :: limr(2), limp(2), indsq
      if( present(shvec) )then
         !!!!!!!!!!!<<<GPU KERNEL START>>>!!!!!!!!!!!
         call self%pdiff_shift_psi_x0(i, shvec)
         pdiff_cc_psi_x0 = self%calc_corr_shifted(i, j, r)
         !!!!!!!!!!!<<<GPU KERNEL END>>>!!!!!!!!!!!
      else
         !!!!!!!!!!!<<<GPU KERNEL START>>>!!!!!!!!!!!
         pdiff_cc_psi_x0 = self%calc_corr(i, j, r)
         !!!!!!!!!!!<<<GPU KERNEL END>>>!!!!!!!!!!!
      endif
    end function pdiff_corr_psi_x0

    !>  \brief  for calculating the correlation between reference i and particle j in rotation r
    function pdiff_corr_psi_y0( self, i, j, r, shvec ) result( pdiff_cc_psi_y0 )
      class(polarft_corrcalc), intent(inout) :: self     !< instance
      integer, intent(in)                    :: i, j, r  !< i=reference, j=particle, r=rotation
      real, intent(in), optional             :: shvec(2) !< origin shift vector
      real :: pdiff_cc_psi_y0
      integer :: limr(2), limp(2), indsq
      if( present(shvec) )then
         !!!!!!!!!!!<<<GPU KERNEL START>>>!!!!!!!!!!!
         call self%pdiff_shift_psi_y0(i, shvec)
         pdiff_cc_psi_y0 = self%calc_corr_shifted(i, j, r)
         !!!!!!!!!!!<<<GPU KERNEL END>>>!!!!!!!!!!!
      else
         !!!!!!!!!!!<<<GPU KERNEL START>>>!!!!!!!!!!!
         pdiff_cc_psi_y0 = self%calc_corr(i, j, r)
         !!!!!!!!!!!<<<GPU KERNEL END>>>!!!!!!!!!!!
      endif
    end function pdiff_corr_psi_y0

    !>  \brief  for calculating the correlation between reference i and particle j in rotation r
    function pdiff_corr_polar_x0( self, i, j, r, shvec ) result( pdiff_cc_x0 )
      class(polarft_corrcalc), intent(inout) :: self     !< instance
      integer, intent(in)                    :: i, j, r  !< i=reference, j=particle, r=rotation
      real, intent(in), optional             :: shvec(2) !< origin shift vector
      real :: pdiff_cc_x0
      integer :: limr(2), limp(2), indsq
      if( present(shvec) )then
         !!!!!!!!!!!<<<GPU KERNEL START>>>!!!!!!!!!!!
         call self%pdiff_shift_polar_x0(i, shvec)
         pdiff_cc_x0 = self%calc_corr_shifted(i, j, r)
         !!!!!!!!!!!<<<GPU KERNEL END>>>!!!!!!!!!!!
      else
         !!!!!!!!!!!<<<GPU KERNEL START>>>!!!!!!!!!!!
         pdiff_cc_x0 = self%calc_corr(i, j, r)
         !!!!!!!!!!!<<<GPU KERNEL END>>>!!!!!!!!!!!
      endif
    end function pdiff_corr_polar_x0

    !>  \brief  for calculating the correlation between reference i and particle j in rotation r
    function pdiff_corr_polar_y0( self, i, j, r, shvec ) result( pdiff_cc_y0 )
      class(polarft_corrcalc), intent(inout) :: self     !< instance
      integer, intent(in)                    :: i, j, r  !< i=reference, j=particle, r=rotation
      real, intent(in), optional             :: shvec(2) !< origin shift vector
      real :: pdiff_cc_y0
      integer :: limr(2), limp(2), indsq
      if( present(shvec) )then
         !!!!!!!!!!!<<<GPU KERNEL START>>>!!!!!!!!!!!
         call self%pdiff_shift_polar_y0(i, shvec)
         pdiff_cc_y0 = self%calc_corr_shifted(i, j, r)
         !!!!!!!!!!!<<<GPU KERNEL END>>>!!!!!!!!!!!
      else
         !!!!!!!!!!!<<<GPU KERNEL START>>>!!!!!!!!!!!
         pdiff_cc_y0 = self%calc_corr(i, j, r)
         !!!!!!!!!!!<<<GPU KERNEL END>>>!!!!!!!!!!!
      endif
    end function pdiff_corr_polar_y0

    !>  \brief  for calculating the correlation between reference i and particle j in rotation r
    
    !!!! THE SHIFTING FUNCTIONALITY IS USED BY SIMPLE_PFTCC_SHSRCH BUT THE SHIFTING IS DONE ON THE FLY    !!!
    !!!! BY THE SIMPLEX OPTIMIZER SO THERE ARE NO BIG MATRICES TO DEAL WITH THAT ARE SUITABLE FOR THE GPU !!!
    
    function corr( self, i, j, r, shvec ) result( cc )
        class(polarft_corrcalc), intent(inout) :: self     !< instance
        integer, intent(in)                    :: i, j, r  !< i=reference, j=particle, r=rotation
        real, intent(in), optional             :: shvec(2) !< origin shift vector
        real :: cc
        integer :: limr(2), limp(2), indsq
        if( present(shvec) )then
            !!!!!!!!!!!<<<GPU KERNEL START>>>!!!!!!!!!!!
            call self%shift(i, shvec)
            cc = self%calc_corr_shifted(i, j, r)
            !!!!!!!!!!!<<<GPU KERNEL END>>>!!!!!!!!!!!
        else
            !!!!!!!!!!!<<<GPU KERNEL START>>>!!!!!!!!!!!
            cc = self%calc_corr(i, j, r)
            !!!!!!!!!!!<<<GPU KERNEL END>>>!!!!!!!!!!!
        endif
    end function

    !>  \brief  for calculating the correlation between reference i and particle j
    function calc_corr( self, i, j, r ) result( cc )
        class(polarft_corrcalc), intent(inout) :: self     !< instance
        integer, intent(in)                    :: i, j, r  !< i=reference, j=particle, r=rotation
        real    :: cc
        integer :: limr(2), limp(2), indsq
        limr  = self%rlim(i)
        limp  = self%plim(j,r)
        indsq = self%pindsqsum(j)
        cc = sum(real(self%pft(limr(1):limr(2),:)*conjg(self%pft(limp(1):limp(2),:))))
        cc = cc/sqrt(self%sqsums(i,1)*self%sqsums(indsq,1)) ! self%sqsums(i,1) because unshifted ref
    end function

    !>  \brief  for calculating the correlation between shifted reference i and particle j
    function calc_corr_shifted( self, i, j, r ) result( cc )
        class(polarft_corrcalc), intent(inout) :: self     !< instance
        integer, intent(in)                    :: i, j, r  !< i=reference, j=particle, r=rotation
        real    :: cc
        integer :: limr(2), limp(2), indsq
        limr  = self%rlimsh(i)
        limp  = self%plim(j,r)
        indsq = self%pindsqsum(j)
        cc = sum(real(self%pft(limr(1):limr(2),:)*conjg(self%pft(limp(1):limp(2),:))))
        cc = cc/sqrt(self%sqsums(i,2)*self%sqsums(indsq,1)) ! self%sqsums(i,2) because shifted ref
    end function

    !> \brief  exhaustive search of rotation and shift parameters. This routine puts
    !!         particle j in register with reference i
    
    !!!! SHOULD NOT RE-IMPLEMENT THE GRID_SRCH OPTION IN THE RE-IMPLEMENTATION, REMOVE THE EXHAUSTIVE MODE    !!!!
    !!!! FROM PRIME. IT IS TOO SLOW AND THE SAME FUNCTIONALITY CAN BE OBTAINED BY SIMPLY NOT GIVING AN ORITAB !!!!
    
    subroutine grid_srch( self, i, j, trs, r, shvec, cc )
        use simple_math, only: rotmat2d
        class(polarft_corrcalc), intent(inout) :: self     !< instance
        integer, intent(in)                    :: i, j     !< i=reference, j=particle
        real, intent(in)                       :: trs      !< origin shift half interval size [-trs,trs]
        integer, intent(out)                   :: r        !< in-plane rotation index output
        real, intent(out)                      :: shvec(2) !< origin shift output
        real, intent(out)                      :: cc       !< cross correlation output
        real    :: corrs(self%nrots), shvec_tst(2), rxsh, rysh, corr, m(2,2)
        integer :: xsh, ysh, sh, loc(1), lims(2,2)
        ! initialize correlation and shift vector
        shvec = 0.
        cc    = -1.
        if( abs(trs) < 0.5 )then
            ! initialize correlation
            shvec = 0.
            cc    = -1.
            ! rotation grid search
            corrs = self%gencorrs(i, j)
            loc   = maxloc(corrs)
            if( corrs(loc(1)) > cc )then
                r  = loc(1)
                cc = corrs(loc(1))
            endif 
        else
            r = 1
            ! translation grid search
            sh = nint(trs)
            do xsh=-sh,sh
                do ysh=-sh,sh
                    shvec_tst = [real(xsh),real(ysh)]
                    ! rotation grid search
                    corrs = self%gencorrs(i, j, shvec_tst)
                    loc   = maxloc(corrs)
                    if( corrs(loc(1)) >= cc )then
                        r     = loc(1)
                        shvec = shvec_tst
                        cc    = corrs(loc(1))
                    endif 
                end do
            end do
            ! origin shift refinement
            lims(1,1) = shvec(1)-0.75
            lims(1,2) = shvec(1)+0.75
            lims(2,1) = shvec(2)-0.75
            lims(2,2) = shvec(2)+0.75
            rxsh = lims(1,1)
            do while(rxsh <= lims(1,2))
                rysh = lims(2,1)
                do while(rysh <= lims(2,2))
                    shvec_tst = [rxsh,rysh]
                    corr = self%corr(i, j, r, shvec_tst)
                    if( corr > cc )then
                        shvec = shvec_tst
                        cc    = corr
                    endif
                    rysh = rysh+0.25
                end do
                rxsh = rxsh+0.25
           end do
           m = rotmat2d(self%get_rot(r)) 
           shvec = matmul(shvec,m)
       endif 
    end subroutine

    ! DESTRUCTOR

    !>  \brief  is a destructor
    subroutine kill( self )
        class(polarft_corrcalc), intent(inout) :: self
        if( self%existence )then
            deallocate(self%pft, self%sqsums, self%angtab, self%polar, self%argtransf)

            deallocate(self%grad_polar_psi) !< dealloc the gradient wrt psi in x0
            deallocate(self%grad_polar_k)   !< dealloc the gradient wrt k in x0
            deallocate(self%grad_argtransf_x0) !< dealloc the gradient wrt k in x0
            deallocate(self%grad_argtransf_y0) !< dealloc the gradient wrt k in x0

            self%existence = .false.
        endif
    end subroutine
    
    ! UNIT TESTS
    
    subroutine test_polarft_corrcalc_inds
        type(polarft_corrcalc) :: pftcc
        integer :: i, j, rlim(2), rlimsh(2)
        integer :: plim(2), pindsqsum(2)
        print *, 'assuming two references & particle range of [3,5]'
        call pftcc%new(2, [3,5], [100,100,1], [2,20], 44)
        call pftcc%print
        print *, '******************'
        do i=1,pftcc%nrefs
            print *, 'reference: ', i
            rlim      = pftcc%rlim(i)
            rlimsh    = pftcc%rlimsh(i)
            print *, 'unshifted reference range: ', rlim
            if( rlim(2)-rlim(1)+1 .ne. pftcc%refsz )then
                stop 'reference  has incorrect range width'
            endif
            if( rlimsh(1) .ne. rlim(2)+1 ) stop 'ERROR 1'
            print *, 'shifted reference range: ', rlimsh
            if( rlimsh(2)-rlimsh(1)+1 .ne. pftcc%refsz )then
                stop 'reference  has incorrect range width'
            endif
        end do
        do j=3,5
            print *, 'particle: ', j
            plim = pftcc%plim(j, 1)
            print *, 'rotind = 1 particle range: ', plim
            if( plim(1) <= pftcc%nrots_ref )then
                stop 'particle is eating into the reference section, out of bound'
            endif
            if( j==3 .and. plim(1) .ne. pftcc%nrots_ref+1 ) stop 'ERROR 4'
            if( plim(2)-plim(1)+1 .ne. pftcc%refsz )then
                stop 'particle (rot=1)  has incorrect range width'
            endif
            pindsqsum = pftcc%pindsqsum(j)
            print *, 'rotind = 1 sqsum indices: ', pindsqsum
            plim = pftcc%plim(j, pftcc%nrots)
            print *, 'rotind = nrots particle range: ', plim
            if( plim(2)-plim(1)+1 .ne. pftcc%refsz )then
                stop 'particle (rot=nrots)  has incorrect range width'
            endif
            if( plim(2) > pftcc%nrots_tot )then
                stop 'particle range is out of bound (righthand side)'
            endif
            pindsqsum = pftcc%pindsqsum(j)
            print *, 'rotind = nrots sqsum indices: ', pindsqsum
        end do
    end subroutine
    
    !>  \brief  is a unit test for polarft
    subroutine test_polarft_corrcalc
        use simple_polarft,   only: polarft
        use simple_image,     only: image
        type(polarft)          :: pft1, pft2
        type(polarft_corrcalc) :: pftcc
        integer, parameter     :: NITER=3000
        real                   :: corr1, corr2, diff
        integer                :: i, j, nrots, r
        real, allocatable      :: corrs1(:), corrs2(:)
        complex, allocatable   :: pft(:,:)
        write(*,'(a)') '**info(simple_polarft_corrcalc_unit_test): testing all functionality'
        call pft1%new([2,11], 45, ptcl=.false.) ! implemented as reference
        call pft2%new([2,11], 45, ptcl=.true.)  ! implemented as particle
        call pft1%set_ldim([100,100,1])
        call pft2%set_ldim([100,100,1])
        nrots = pft1%get_nrots()
        allocate(corrs1(nrots), corrs2(nrots))
        ! make a conforming polarft_corrcalc object of 3 refs and 5 ptcls
        call pftcc%new(3, [1,5], [100,100,1], [2,11], 45)
        write(*,'(a)') '**info(simple_polarft_corrcalc_unit_test, part 1):&
        comparing the newly implemented correlation with the original one'
        do i=1,3
            do j=1,5
                call pft1%noise_fill
                call pft2%noise_fill
                call pft1%memoize_sqsums
                call pft2%memoize_sqsums
                ! pft1 is ref
                pft = pft1%get_pft()
                call pftcc%set_ref_pft(i, pft)
                deallocate(pft)
                ! pft2 is ptcl
                pft = pft2%get_pft()
                call pftcc%set_ptcl_pft(j, pft)
                deallocate(pft)
                ! compare the two correlation functions
                diff = 0.
                do r=1,nrots
                    corr1 = pft1%corr(pft2, r)  ! the original correlation
                    corr2 = pftcc%corr(i, j, r) ! the newly implemented one
                    diff = diff+abs(corr1-corr2)
                end do
                diff = diff/real(nrots)
                if( diff > 1.5e-2 )then
                    print *, 'diff = ', diff
                    stop 'the new correlation does not conform with the old; test_polarft_corrcalc'
                endif
            end do
        end do
        write(*,'(a)') '**info(simple_polarft_corrcalc_unit_test, part 2): testing the gencorrs function'
        do i=1,3
            do j=1,5
                call pft1%noise_fill
                call pft2%noise_fill
                call pft1%memoize_sqsums
                call pft2%memoize_sqsums
                ! pft1 is ref
                pft = pft1%get_pft()
                call pftcc%set_ref_pft(i, pft)
                deallocate(pft)
                ! pft2 is ptcl
                pft = pft2%get_pft()
                call pftcc%set_ptcl_pft(j, pft)
                deallocate(pft)
                ! compare the two correlation functions
                diff = 0.
                call pft1%gencorrs(pft2, corrs1)
                corrs2 = pftcc%gencorrs(i,j)
                diff = sum(abs(corrs1-corrs2))/real(nrots*3*5)
                if( diff > 1.5e-2 )then
                    print *, 'diff = ', diff
                    stop 'the new correlation does not conform with the old; test_polarft_corrcalc'
                endif
            end do
        end do
        deallocate(corrs1, corrs2)
        write(*,'(a)') 'SIMPLE_POLARFT_CORRCALC_UNIT_TEST COMPLETED SUCCESSFULLY ;-)'
    end subroutine

end module simple_polarft_corrcalc
