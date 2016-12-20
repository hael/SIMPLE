!>  \brief  SIMPLE polarft_corrcalc class
module simple_polarft_corrcalc
use simple_defs      ! singleton
use simple_params,   only: params
use simple_ran_tabu, only: ran_tabu
use simple_jiffys
implicit none

public :: polarft_corrcalc, test_polarft_corrcalc
private

! CLASS PARAMETERS/VARIABLES
complex, parameter :: zero=cmplx(0.,0.) !< just a complex zero

type :: polarft_corrcalc
    private
    integer                          :: pfromto(2) = 1        !< from/to particle indices (in parallel execution)
    integer                          :: nptcls     = 1        !< the total number of particles in partition (logically indexded [fromp,top])
    integer                          :: nrefs      = 1        !< the number of references (logically indexded [1,nrefs])
    integer                          :: nrots      = 0        !< number of in-plane rotations for one pft (determined by radius of molecule)
    integer                          :: ring2      = 0        !< radius of molecule
    integer                          :: refsz      = 0        !< size of reference (nrots/2) (number of vectors used for matching)
    integer                          :: ptclsz     = 0        !< size of particle (2*nrots)
    integer                          :: ldim(3)    = 0        !< logical dimensions of original cartesian image
    integer                          :: kfromto(2) = 0        !< Fourier index range
    integer                          :: nk         = 0        !< number of resolution elements (Fourier components)
    real,        allocatable         :: sqsums_refs(:)        !< memoized square sums for the correlation calculations
    real,        allocatable         :: sqsums_ptcls(:)       !< memoized square sums for the correlation calculations
    real,        allocatable         :: angtab(:)             !< table of in-plane angles (in degrees)
    real,        allocatable         :: argtransf(:,:)        !< argument transfer constants for shifting the references
    real,        allocatable         :: polar(:,:)            !< table of polar coordinates (in Cartesian coordinates)
    complex(sp), allocatable         :: pfts_refs(:,:,:)      !< 3D complex matrix of polar reference sections (nrefs,refsz,nk)
    complex(sp), allocatable         :: pfts_refs_ctf(:,:,:)  !< 3D complex matrix of polar reference sections with CTF applied
    complex(sp), allocatable         :: pfts_ptcls(:,:,:)     !< 3D complex matrix of particle sections
    real,        allocatable         :: pfts_ptcls_ctf(:,:,:) !< real valued CTF values for implementation of the CTF on GPU 
    logical                          :: xfel=.false.          !< to indicate whether we process xfel patterns or not
    logical                          :: dim_expanded=.false.  !< to indicate whether dim has been expanded or not
    logical                          :: existence=.false.     !< to indicate existence
  contains
    ! CONSTRUCTOR
    procedure          :: new
    ! INDEX MAPPING ROUTINES
    procedure, private :: k_ind
    procedure, private :: ptcl_ind
    ! SETTERS
    procedure          :: set_ref_pft
    procedure          :: set_ptcl_pft
    procedure          :: set_ref_fcomp
    procedure          :: set_ptcl_fcomp
    procedure          :: cp_ptcls2refs
    ! GETTERS
    procedure          :: get_pfromto
    procedure          :: get_nptcls
    procedure          :: get_nrefs
    procedure          :: get_nrots
    procedure          :: get_ring2
    procedure          :: get_refsz
    procedure          :: get_ptclsz
    procedure          :: get_ldim
    procedure          :: get_kfromto
    procedure          :: get_rot
    procedure          :: get_roind
    procedure          :: get_coord
    procedure          :: get_ptcl_pft
    procedure          :: get_ref_pft
    procedure          :: exists
    ! PRINTERS/VISUALISERS
    procedure          :: print
    procedure          :: vis_ptcl
    procedure          :: vis_ref    
    ! MEMOIZERS
    procedure          :: memoize_sqsum_ref
    procedure          :: memoize_sqsum_ref_ctf
    procedure, private :: memoize_sqsum_ptcl
    ! MODIFIERS
    procedure          :: apply_ctf
    procedure          :: xfel_subtract_shell_mean
    ! GPU PREP ROUTINES
    procedure          :: expand_dim
    procedure          :: prep_ctf4gpu
    ! CALCULATORS
    procedure, private :: create_polar_ctfmat
    procedure          :: gencorrs_all_cpu
    procedure, private :: gencorrs_all_tester_1
    procedure, private :: gencorrs_all_tester_2
    generic            :: gencorrs_all_tester => gencorrs_all_tester_1, gencorrs_all_tester_2
    procedure          :: gencorrs
    procedure          :: rosrch_shc
    procedure, private :: corr_1
    procedure, private :: corr_2
    generic            :: corr => corr_1, corr_2
    ! DESTRUCTOR
    procedure          :: kill
end type polarft_corrcalc

contains
    
    ! CONSTRUCTORS
    
    !>  \brief  is a constructor
    subroutine new( self, nrefs, pfromto, ldim, kfromto, ring2, ctfflag, isxfel )
        use simple_math, only: rad2deg, is_even, round2even
        class(polarft_corrcalc),    intent(inout) :: self
        integer,                    intent(in)    :: nrefs, pfromto(2), ldim(3), kfromto(2), ring2
        character(len=*),           intent(in)    :: ctfflag
        character(len=*), optional, intent(in)    :: isxfel
        integer     :: alloc_stat, irot, k, k_ind, err!, idev
        logical     :: even_dims, test(3)
        real        :: ang
        ! function call
        ! integer :: get_dev_count_c
        ! kill possibly pre-existing object
        call self%kill
        ! error check
        if( kfromto(2)-kfromto(1) <= 2 )then
            write(*,*) 'kfromto: ', kfromto(1), kfromto(2)
            stop 'resolution range too narrow; new; simple_polarft_corrcalc'
        endif
        if( ring2 < 1 )then
            write(*,*) 'ring2: ', ring2
            stop 'ring2 must be > 0; new; simple_polarft_corrcalc'
        endif
        if( pfromto(2)-pfromto(1)+1 < 1 )then
            write(*,*) 'pfromto: ', pfromto(1), pfromto(2)
            stop 'nptcls (# of particles) must be > 0; new; simple_polarft_corrcalc'
        endif
        if( nrefs < 1 )then
            write(*,*) 'nrefs: ', nrefs
            stop 'nrefs (# of reference sections) must be > 0; new; simple_polarft_corrcalc'
        endif
        if( any(ldim == 0) )then
            write(*,*) 'ldim: ', ldim
            stop 'ldim is not conforming (is zero); new; simple_polarft'
        endif
        if( ldim(3) > 1 )then
            write(*,*) 'ldim: ', ldim
            stop '3D polarfts are not yet supported; new; simple_polarft_corrcalc'
        endif
        test = .false.
        test(1) = is_even(ldim(1))
        test(2) = is_even(ldim(2))
        test(3) = ldim(3) == 1
        even_dims = all(test)
        if( .not. even_dims )then
            write(*,*) 'ldim: ', ldim
            stop 'only even logical dims supported; new; simple_polarft_corrcalc'
        endif
        self%xfel = .false.
        if( present(isxfel) )then
            if( isxfel .eq. 'yes' ) self%xfel = .true.
        end if
        ! set constants
        self%pfromto = pfromto                       !< from/to particle indices (in parallel execution)
        self%nptcls  = pfromto(2)-pfromto(1)+1       !< the total number of particles in partition (logically indexded [fromp,top])
        self%nrefs   = nrefs                         !< the number of references (logically indexded [1,nrefs])
        self%ring2   = ring2                         !< radius of molecule
        self%nrots   = round2even(twopi*real(ring2)) !< number of in-plane rotations for one pft  (determined by radius of molecule)
        self%refsz   = self%nrots/2                  !< size of reference (nrots/2) (number of vectors used for matching)
        self%ptclsz  = self%nrots*2                  !< size of particle (2*nrots)
        self%ldim    = ldim                          !< logical dimensions of original cartesian image
        self%kfromto = kfromto                       !< Fourier index range
        self%nk      = kfromto(2)-kfromto(1)+1       !< number of resolution elements (Fourier components)
        ! generate polar coordinates
        allocate( self%polar(self%ptclsz,self%nk), self%angtab(self%nrots), stat=alloc_stat)
        call alloc_err('polar coordinate arrays; new; simple_polarft_corrcalc', alloc_stat)
        ang = twopi/real(self%nrots)
        do irot=1,self%nrots
            self%angtab(irot) = (irot-1)*ang
            do k=self%kfromto(1),self%kfromto(2)
                k_ind = self%k_ind(k) 
                self%polar(irot,k_ind)            = cos(self%angtab(irot))*real(k) ! x-coordinate
                self%polar(irot+self%nrots,k_ind) = sin(self%angtab(irot))*real(k) ! y-coordinate
            end do
            self%angtab(irot) = rad2deg(self%angtab(irot)) ! angle (in degrees)
        end do
        ! generate the argument transfer constants for shifting reference polarfts
        allocate( self%argtransf(self%nrots,self%nk), stat=alloc_stat)
        call alloc_err('shift argument transfer array; new; simple_polarft_corrcalc', alloc_stat)
        self%argtransf(:self%refsz,:)   = &
            self%polar(:self%refsz,:)   * &
            (PI/real(self%ldim(1)/2))    ! x-part
        self%argtransf(self%refsz+1:,:) = &
            self%polar(self%nrots+1:self%nrots+self%refsz,:) * &
            (PI/real(self%ldim(2)/2))    ! y-part
        ! allocate polarfts and sqsums
        allocate(   self%pfts_refs(self%nrefs,self%refsz,self%nk),&
                    self%pfts_ptcls(self%nptcls,self%ptclsz,self%nk),&
                    self%sqsums_refs(self%nrefs),&
                    self%sqsums_ptcls(self%nptcls), stat=alloc_stat)
        call alloc_err('polarfts and sqsums; new; simple_polarft_corrcalc', alloc_stat)
        self%pfts_refs    = zero
        self%pfts_ptcls   = zero
        self%sqsums_refs  = 0.
        self%sqsums_ptcls = 0.
        ! pfts_refs_ctf if needed
        if( ctfflag .ne. 'no' )then
            allocate(self%pfts_refs_ctf(self%nrefs,self%refsz,self%nk), stat=alloc_stat)
            call alloc_err('pfts_refs_ctf; new; simple_polarft_corrcalc', alloc_stat)
            self%pfts_refs_ctf = zero
            ! we will allocate the congruent CTF matrix for GPU execution when we expand the dimensions
        endif
        self%dim_expanded = .false.
        self%existence    = .true.
    end subroutine new

    ! INDEX MAPPING ROUTINES

    !>  \brief  for mapping physical Fourier index to implemented
    function k_ind( self, k_in ) result( k_out )
        class(polarft_corrcalc), intent(in) :: self
        integer,                 intent(in) :: k_in
        integer :: k_out
        k_out = k_in-self%kfromto(1)+1
    end function k_ind
    
    !>  \brief  for mapping physical particle index to implemented
    function ptcl_ind( self, iptcl_in ) result( iptcl_out )
        class(polarft_corrcalc), intent(in) :: self
        integer,                 intent(in) :: iptcl_in
        integer :: iptcl_out
        iptcl_out = iptcl_in-self%pfromto(1)+1
    end function ptcl_ind

    ! SETTERS

    !>  \brief  sets reference pft iref
    subroutine set_ref_pft( self, iref, pft )
        class(polarft_corrcalc), intent(inout) :: self
        integer,                 intent(in)    :: iref
        complex,                 intent(in)    :: pft(:,:)
        self%pfts_refs(iref,:,:) = pft(:self%refsz,:)
        ! calculate the square sum required for correlation calculation
        call self%memoize_sqsum_ref(iref)
    end subroutine set_ref_pft

    !>  \brief  sets particle pft iptcl
    subroutine set_ptcl_pft( self, iptcl, pft )
        class(polarft_corrcalc), intent(inout) :: self
        integer,                 intent(in)    :: iptcl
        complex,                 intent(in)    :: pft(:,:)
        integer :: ptcl_ind
        ptcl_ind = self%ptcl_ind(iptcl)
        self%pfts_ptcls(ptcl_ind,:self%nrots,:)   = pft
        self%pfts_ptcls(ptcl_ind,self%nrots+1:,:) = pft ! because rot dim is expanded
        ! calculate the square sum required for correlation calculation
        call self%memoize_sqsum_ptcl(iptcl)
    end subroutine set_ptcl_pft

    !>  \brief  sets a reference Fourier component
    subroutine set_ref_fcomp( self, iref, irot, k, comp )
        class(polarft_corrcalc), intent(inout) :: self
        integer,                 intent(in)    :: iref, irot, k
        complex,                 intent(in)    :: comp
        integer :: k_ind
        k_ind = self%k_ind(k)
        self%pfts_refs(iref,irot,k_ind) = comp
    end subroutine set_ref_fcomp
    
    !>  \brief  sets a particle Fourier component
    subroutine set_ptcl_fcomp( self, iptcl, irot, k, comp )
        class(polarft_corrcalc), intent(inout) :: self
        integer,                 intent(in)    :: iptcl, irot, k
        complex,                 intent(in)    :: comp
        integer                                :: ptcl_ind, k_ind
        ptcl_ind = self%ptcl_ind(iptcl)
        k_ind    = self%k_ind(k)
        self%pfts_ptcls(ptcl_ind,irot,k_ind) = comp
        self%pfts_ptcls(ptcl_ind,irot+self%nrots,k_ind) = comp ! because rot dim is expanded
    end subroutine set_ptcl_fcomp

    !>  \brief  copies the particles to the references
    subroutine cp_ptcls2refs( self )
        class(polarft_corrcalc), intent(inout) :: self
        if( self%nrefs .eq. self%nptcls )then
            self%pfts_refs(:,:,:) = self%pfts_ptcls(:,:self%refsz,:)
            self%sqsums_refs = self%sqsums_ptcls           
        else
            stop 'pfts_refs and pfts_ptcls not congruent (nrefs .ne. nptcls)'
        endif
    end subroutine cp_ptcls2refs
    
    ! GETTERS
    
    !>  \brief  for getting the logical particle range
    function get_pfromto( self ) result( lim )
        class(polarft_corrcalc), intent(inout) :: self
        integer :: lim(2)
        lim = self%pfromto
    end function get_pfromto
    
    !>  \brief  for getting the number of particles
    pure function get_nptcls( self ) result( nptcls )
        class(polarft_corrcalc), intent(in) :: self
        integer :: nptcls
        nptcls = self%nptcls
    end function get_nptcls
    
    !>  \brief  for getting the number of references
    pure function get_nrefs( self ) result( nrefs )
        class(polarft_corrcalc), intent(in) :: self
        integer :: nrefs
        nrefs = self%nrefs
    end function get_nrefs
    
    !>  \brief  for getting the number of in-plane rotations
    pure function get_nrots( self ) result( nrots )
        class(polarft_corrcalc), intent(in) :: self
        integer :: nrots
        nrots = self%nrots
    end function get_nrots
    
    !>  \brief  for getting the particle radius (ring2)
    function get_ring2( self ) result( ring2 )
        class(polarft_corrcalc), intent(in) :: self
        integer :: ring2
        ring2 = self%ring2
    end function get_ring2
    
    !>  \brief  for getting the number of reference rotations (size of second dim of self%pfts_refs)
    function get_refsz( self ) result( refsz )
        class(polarft_corrcalc), intent(in) :: self
        integer :: refsz
        refsz = self%refsz
    end function get_refsz
    
    !>  \brief  for getting the number of particle rotations (size of second dim of self%pfts_ptcls)
    function get_ptclsz( self ) result( ptclsz )
        class(polarft_corrcalc), intent(in) :: self
        integer :: ptclsz
        ptclsz = self%ptclsz
    end function get_ptclsz
    
    !>  \brief  for getting the logical dimension of the original
    !!          Cartesian image
    function get_ldim( self ) result( ldim )
        class(polarft_corrcalc), intent(in) :: self
        integer :: ldim(3)
        ldim = self%ldim
    end function get_ldim
    
    !>  \brief  for getting the Fourier index range (hp/lp)
    function get_kfromto( self ) result( lim )
        class(polarft_corrcalc), intent(inout) :: self
        integer :: lim(2)
        lim = self%kfromto
    end function get_kfromto
    
    !>  \brief is for getting the continuous in-plane rotation
    !!         corresponding to in-plane rotation index roind
    function get_rot( self, roind ) result( rot )
        class(polarft_corrcalc), intent(in) :: self
        integer,                 intent(in) :: roind
        real :: rot
        if( roind < 1 .or. roind > self%nrots )then
            stop 'roind is out of range; get_rot; simple_polarft_corrcalc'
        endif
        rot = self%angtab(roind)
    end function get_rot

    !>  \brief is for getting the discrete in-plane rotational
    !!         index corresponding to continuous rotation rot
    function get_roind( self, rot ) result( ind )
        class(polarft_corrcalc), intent(in) :: self
        real,                    intent(in) :: rot
        integer :: ind, irot, alloc_stat, loc(1)
        real    :: dists(self%nrots)
        do irot=1,self%nrots
            dists(irot) = sqrt((self%angtab(irot)-rot)**2.)
        end do
        loc = minloc(dists)
        ind = loc(1)
    end function get_roind

    !>  \brief returns polar coordinate for rotation rot
    !!         and Fourier index k
    function get_coord( self, rot, k ) result( xy )
        class(polarft_corrcalc), intent(in) :: self
        integer,                 intent(in) :: rot, k
        real    :: xy(2)
        integer :: k_ind
        k_ind = self%k_ind(k)    
        xy(1) = self%polar(rot,k_ind)
        xy(2) = self%polar(self%nrots+rot,k_ind)
    end function get_coord
    
    !>  \brief  returns polar Fourier transform of particle iptcl in rotation irot
    function get_ptcl_pft( self, iptcl, irot ) result( pft )
        class(polarft_corrcalc), intent(in) :: self
        integer,                 intent(in) :: iptcl, irot
        complex, allocatable :: pft(:,:)
        integer :: alloc_stat, ptcl_ind
        ptcl_ind = self%ptcl_ind(iptcl)
        allocate(pft(self%refsz,self%kfromto(1):self%kfromto(2)),&
        source=self%pfts_ptcls(ptcl_ind,irot:irot+self%refsz-1,:), stat=alloc_stat)
        call alloc_err("In: get_ptcl_pft; simple_polarft_corrcalc", alloc_stat)
    end function get_ptcl_pft
    
    !>  \brief  returns polar Fourier transform of reference iref
    function get_ref_pft( self, iref ) result( pft )
        class(polarft_corrcalc), intent(in) :: self
        integer,                 intent(in) :: iref
        complex, allocatable :: pft(:,:)
        integer :: alloc_stat
        allocate(pft(self%refsz,self%kfromto(1):self%kfromto(2)),&
        source=self%pfts_refs(iref,:,:), stat=alloc_stat)
        call alloc_err("In: get_ref_pft; simple_polarft_corrcalc", alloc_stat)
    end function get_ref_pft

    !>  \brief  checks for existence
    function exists( self ) result( yes )
        class(polarft_corrcalc), intent(in) :: self
        logical :: yes
        yes = self%existence
    end function exists

    ! PRINTERS/VISUALISERS

    !>  \brief  is for plotting a particle polar FT
    subroutine vis_ptcl( self, iptcl )
        use gnufor2
        class(polarft_corrcalc), intent(in) :: self
        integer,                 intent(in) :: iptcl
        integer :: ptcl_ind
        ptcl_ind = self%ptcl_ind(iptcl)
        call gnufor_image(real(self%pfts_ptcls(ptcl_ind,:self%refsz,:)),  palette='gray')
        call gnufor_image(aimag(self%pfts_ptcls(ptcl_ind,:self%refsz,:)), palette='gray')
    end subroutine vis_ptcl
    
    !>  \brief  is for plotting a particle polar FT
    subroutine vis_ref( self, iref )
        use gnufor2
        class(polarft_corrcalc), intent(in) :: self
        integer,                 intent(in) :: iref
        call gnufor_image(real(self%pfts_refs(iref,:,:)),  palette='gray')
        call gnufor_image(aimag(self%pfts_refs(iref,:,:)), palette='gray')
    end subroutine vis_ref
      
    !>  \brief  for printing info about the object
    subroutine print( self )
      class(polarft_corrcalc), intent(in) :: self
      write(*,*) "from/to particle indices              (self%pfromto): ", self%pfromto
      write(*,*) "total n particles in partition         (self%nptcls): ", self%nptcls
      write(*,*) "number of references                    (self%nrefs): ", self%nrefs
      write(*,*) "number of rotations                     (self%nrots): ", self%nrots
      write(*,*) "number of nk                               (self%nk): ", self%nk
      write(*,*) "radius of molecule                      (self%ring2): ", self%ring2
      write(*,*) "nr of rots for ref (2nd dim of pftmat)  (self%refsz): ", self%refsz
      write(*,*) "n rots for ptcl (2nd dim of pftmat)    (self%ptclsz): ", self%ptclsz 
      write(*,*) "logical dim. of original Cartesian image (self%ldim): ", self%ldim
      write(*,*) "high-pass limit Fourier index      (self%kfromto(1)): ", self%kfromto(1)
      write(*,*) "low-pass limit Fourier index       (self%kfromto(2)): ", self%kfromto(2)
      return
    end subroutine print
   
    ! MEMOIZERS

    !>  \brief  is for memoization of the complex square sums reqired for correlation calculation
    subroutine memoize_sqsum_ref( self, iref )
        use simple_math, only: csq
        class(polarft_corrcalc), intent(inout) :: self
        integer,                 intent(in)    :: iref
        if( self%xfel )then
            self%sqsums_refs(iref) = sum(real(self%pfts_refs(iref,:,:))**2.)
        else
            self%sqsums_refs(iref) = sum(csq(self%pfts_refs(iref,:,:)))
        endif
    end subroutine memoize_sqsum_ref
    
    !>  \brief  is for memoization of the complex square sums reqired for correlation calculation
    subroutine memoize_sqsum_ref_ctf( self, iref )
        use simple_math, only: csq
        class(polarft_corrcalc), intent(inout) :: self
        integer,                 intent(in)    :: iref
        self%sqsums_refs(iref) = sum(csq(self%pfts_refs_ctf(iref,:,:)))
    end subroutine memoize_sqsum_ref_ctf

    !>  \brief  is for memoization of the complex square sums reqired for correlation calculation
    subroutine memoize_sqsum_ptcl( self, iptcl )
        use simple_math, only: csq
        class(polarft_corrcalc), intent(inout) :: self
        integer,                 intent(in)    :: iptcl
        integer :: ptcl_ind
        ptcl_ind = self%ptcl_ind(iptcl)
        if( self%xfel )then
            self%sqsums_ptcls(ptcl_ind) = sum(real(self%pfts_ptcls(ptcl_ind,:self%refsz,:self%nk))**2.)
        else
            self%sqsums_ptcls(ptcl_ind) = sum(csq(self%pfts_ptcls(ptcl_ind,:self%refsz,:self%nk)))
        endif
    end subroutine memoize_sqsum_ptcl
    
    ! MODIFIERS

    !>  \brief  is for applying CTF to references and updating the memoized ref sqsums
    subroutine apply_ctf( self, smpd, tfun, dfx, dfy, angast, ref )
        !$ use omp_lib
        !$ use omp_lib_kinds
        use simple_ctf,   only: ctf
        use simple_image, only: image
        class(polarft_corrcalc), intent(inout) :: self
        real,                    intent(in)    :: smpd
        class(ctf),              intent(inout) :: tfun
        real,                    intent(in)    :: dfx
        real,    optional,       intent(in)    :: dfy, angast
        integer, optional,       intent(in)    :: ref
        real, allocatable :: ctfmat(:,:)
        integer           :: iref
        type(image)       :: img
        ! create the congruent polar matrix of real CTF values
        ctfmat = self%create_polar_ctfmat(tfun, dfx, dfy, angast, self%refsz)
        ! multiply the references with the CTF
        if( present(ref) )then
            ! apply to single reference
            if( ref<1 .or. ref>self%nrefs )stop 'reference index out of bounds; simple_polarft_corrcalc::apply_ctf'
            iref = ref
            self%pfts_refs_ctf(iref,:,:) = self%pfts_refs(iref,:,:)*ctfmat
            call self%memoize_sqsum_ref_ctf(iref)
        else
            ! apply to all references
            !$omp parallel do default(shared) schedule(auto) private(iref)
            do iref=1,self%nrefs
                self%pfts_refs_ctf(iref,:,:) = self%pfts_refs(iref,:,:)*ctfmat
                call self%memoize_sqsum_ref_ctf(iref)
            end do
            !$omp end parallel do
        endif
        deallocate(ctfmat)
    end subroutine apply_ctf

    !>  \brief  is for preparing for XFEL pattern corr calc
    subroutine xfel_subtract_shell_mean( self )
        class(polarft_corrcalc), intent(inout) :: self
        real, allocatable    :: ptcls_mean_tmp(:,:,:)
        real, allocatable    :: refs_mean_tmp(:,:)
        integer :: iptcl, iref, irot, k, k_ind
        allocate( ptcls_mean_tmp(2*self%nptcls,self%ptclsz,self%nk), refs_mean_tmp(self%nrefs,self%nk))
        ! calculate the mean of each reference at each k shell
        do iref=1,self%nrefs
            do k=self%kfromto(1),self%kfromto(2)
                k_ind = self%k_ind(k)
                refs_mean_tmp(iref,k_ind) = sum(self%pfts_refs(iref,:, k_ind))/self%refsz
            end do
        end do
        ! calculate the mean of each reference at each k shell
        do iref=1,self%nrefs
            do irot=1,self%refsz
                do k=self%kfromto(1),self%kfromto(2)
                    k_ind = self%k_ind(k)
                    self%pfts_refs(iref,irot, k_ind) = &
                    self%pfts_refs(iref,irot, k_ind) - refs_mean_tmp(iref,k_ind) 
                end do
            end do
        end do
        ! calculate the mean of each particle at each k shell at each in plane rotation
        do iptcl=self%pfromto(1),self%pfromto(2)
            do k=self%kfromto(1),self%kfromto(2)
                k_ind = self%k_ind(k)
                ptcls_mean_tmp(iptcl,1,k_ind) = sum(self%pfts_ptcls(iptcl,1:self%refsz-1,k_ind))/self%refsz
            end do
        end do
        ! subtract the mean of each particle at each k shell at each in plane rotation
        do iptcl=self%pfromto(1),self%pfromto(2)
            do irot=1,self%ptclsz
                do k=self%kfromto(1),self%kfromto(2)
                    k_ind = self%k_ind(k)
                    self%pfts_ptcls(iptcl,irot, k_ind) = &
                    self%pfts_ptcls(iptcl,irot, k_ind) - ptcls_mean_tmp(iptcl,1,k_ind)
                end do
            end do
        end do
        deallocate( ptcls_mean_tmp, refs_mean_tmp )
    end subroutine xfel_subtract_shell_mean

    ! GPU PREP ROUTINES

    !>  \brief  is for expanding the projection direction dimension to allow efficient
    !!          transfer to the GPU via subsectioning
    subroutine expand_dim( self )
        !$ use omp_lib
        !$ use omp_lib_kinds
        class(polarft_corrcalc), intent(inout) :: self
        complex, allocatable :: pfts_ptcls_tmp(:,:,:)
        real, allocatable    :: sqsums_ptcls_tmp(:)
        integer :: alloc_stat
        if( size(self%pfts_ptcls,dim=1) == 2*self%nptcls )then
            ! the arrays are already expanded, all we need to do is to copy data over
            !$omp parallel workshare
            self%pfts_ptcls(self%nptcls+1:,:,:) = self%pfts_ptcls(:self%nptcls,:,:)
            self%sqsums_ptcls(self%nptcls+1:)   = self%sqsums_ptcls(:self%nptcls)
            !$omp end parallel workshare
        else if( size(self%pfts_ptcls,dim=1) == self%nptcls )then
            allocate( pfts_ptcls_tmp(self%nptcls,self%ptclsz,self%nk), source=self%pfts_ptcls, stat=alloc_stat )
            call alloc_err("In: expand_dim, simple_polarft_corrcalc, pfts_ptcls_tmp", alloc_stat)
            allocate( sqsums_ptcls_tmp(self%nptcls), source=self%sqsums_ptcls, stat=alloc_stat )
            call alloc_err("In: expand_dim, simple_polarft_corrcalc, sqsums_ptcls_tmp", alloc_stat)
            deallocate(self%pfts_ptcls)
            deallocate(self%sqsums_ptcls)
            allocate( self%pfts_ptcls(2*self%nptcls,self%ptclsz,self%nk), stat=alloc_stat )
            call alloc_err("In: expand_dim, simple_polarft_corrcalc, self%pfts_ptcls", alloc_stat)
            allocate( self%sqsums_ptcls(2*self%nptcls), stat=alloc_stat )
            call alloc_err("In: expand_dim, simple_polarft_corrcalc, 4", alloc_stat)
            !$omp parallel workshare
            self%pfts_ptcls(:self%nptcls,:,:)   = pfts_ptcls_tmp
            self%pfts_ptcls(self%nptcls+1:,:,:) = pfts_ptcls_tmp
            self%sqsums_ptcls(:self%nptcls)     = sqsums_ptcls_tmp
            self%sqsums_ptcls(self%nptcls+1:)   = sqsums_ptcls_tmp
            !$omp end parallel workshare
        else
            stop 'nonconforming dimension (1) of self%pfts_ptcls; expand_dim; simple_polarft_corrcalc'
        endif
        self%dim_expanded = .true.
    end subroutine expand_dim

    !>  \brief  prepare the real CTF matrix for GPU-based corr calc
    subroutine prep_ctf4gpu( self, smpd, a, tfun )
        !$ use omp_lib
        !$ use omp_lib_kinds
        use simple_ctf,   only: ctf
        use simple_oris,  only: oris
        use simple_image, only: image
        class(polarft_corrcalc), intent(inout) :: self
        real,                    intent(in)    :: smpd
        class(oris),             intent(inout) :: a
        class(ctf),              intent(inout) :: tfun
        type(image)       :: img
        real, allocatable :: ctfparams(:,:), ctfmat(:,:)
        integer           :: alloc_stat, iptcl
        ! create template image
        call img%new(self%ldim, smpd)
        ! get the CTF parameters
        ctfparams = a%get_ctfparams(self%pfromto)
        ! allocate the array that will carry the real CTF values
        if( allocated(self%pfts_ptcls_ctf) ) deallocate(self%pfts_ptcls_ctf)
        allocate( self%pfts_ptcls_ctf(2*self%nptcls,self%ptclsz,self%nk), stat=alloc_stat )
        call alloc_err("In: prep_ctf4gpu, simple_polarft_corrcalc, self%pfts_ptcls_ctf", alloc_stat)
        !$omp parallel default(shared) private(iptcl,ctfmat)
        !$omp do schedule(auto)
        do iptcl=1,self%nptcls
            ! create the congruent polar matrix of real CTF values
            ctfmat = self%create_polar_ctfmat(tfun, ctfparams(iptcl,1),&
                     ctfparams(iptcl,2), ctfparams(iptcl,3), self%nrots)
            ! set it
            self%pfts_ptcls_ctf(iptcl,:self%nrots,:)   = ctfmat
            ! expand in in-plane rotation dimension
            self%pfts_ptcls_ctf(iptcl,self%nrots+1:,:) = ctfmat
        end do
        !$omp end do nowait
        !$omp workshare
        ! expand in projection direction dimension
        self%pfts_ptcls_ctf(self%nptcls+1:,:,:) = self%pfts_ptcls_ctf(:self%nptcls,:,:)
        !$omp end workshare nowait
        !$omp end parallel
        deallocate(ctfmat)
        call img%kill
    end subroutine prep_ctf4gpu

    ! CALCULATORS

    !>  \brief  is for generating a matrix of CTF values
    function create_polar_ctfmat( self, tfun, dfx, dfy, angast, endrot ) result( ctfmat )
        !$ use omp_lib
        !$ use omp_lib_kinds
        use simple_ctf, only: ctf
        class(polarft_corrcalc), intent(inout) :: self
        class(ctf),              intent(inout) :: tfun
        real,                    intent(in)    :: dfx
        real, optional,          intent(in)    :: dfy, angast
        integer,                 intent(in)    :: endrot
        real, allocatable :: ctfmat(:,:)
        real              :: inv_ldim(3),hinv,kinv,spaFreqSq,ddfy,aangast,ang
        integer           :: irot,k,k_ind
        allocate( ctfmat(endrot,self%nk) )
        ddfy = dfx
        if( present(dfy) ) ddfy = dfy
        aangast = 0.
        if( present(angast) ) aangast = angast    
        inv_ldim = 1./real(self%ldim)
        !$omp parallel do default(shared) private(irot,k,k_ind,hinv,kinv,spaFreqSq,ang) schedule(auto)
        do irot=1,endrot
            do k=self%kfromto(1),self%kfromto(2)
                k_ind              = self%k_ind(k)
                hinv               = self%polar(irot,k_ind)*inv_ldim(1)
                kinv               = self%polar(irot+self%nrots,k_ind)*inv_ldim(2)
                spaFreqSq          = hinv*hinv+kinv*kinv
                ang                = atan2(self%polar(irot+self%nrots,k_ind),self%polar(irot,k_ind))
                ctfmat(irot,k_ind) = tfun%eval(spaFreqSq,dfx,ddfy,aangast,ang)
            end do
        end do
        !$omp end parallel do
    end function create_polar_ctfmat

    !>  \brief  is for generating all rotational correlations for all nptcls particles and all nrefs
    !!          references (nptcls=nrefs) (assumes that the first and second dimensions of the particle
    !!          matrix pfts_ptcls have been expanded and that self%nptcls=self%nrefs)
    subroutine gencorrs_all_cpu( self, corrmat3dout )
        class(polarft_corrcalc), intent(inout) :: self !< instance
        real,                    intent(out)   :: corrmat3dout(self%nptcls,self%nptcls,self%nrots)
        real    :: hadamard_prod(self%nptcls,self%refsz,self%nk)
        integer :: iref, iptcl, irot, d1lims(2), d2lims(2)
        if( .not. self%dim_expanded )&
        stop 'expand dimension with expand_dim before calling simple_polarft_corrcalc :: gencorrs_all_cpu'
        do iref=1,self%nptcls
            d1lims(1) = iref
            d1lims(2) = iref + self%nptcls - 1
            do irot=1,self%nrots
                d2lims(1) = irot 
                d2lims(2) = irot + self%refsz - 1
                ! calculate the Hadamard product
                !$omp parallel default(shared) private(iptcl)
                !$omp workshare
                hadamard_prod(:,:,:) = &
                real(self%pfts_refs(:,:,:)*conjg(self%pfts_ptcls(d1lims(1):d1lims(2),d2lims(1):d2lims(2),:)))
                !$omp end workshare
                ! sum over the rotational and resolution dimensions
                !$omp do schedule(auto)
                do iptcl=1,self%nptcls
                    corrmat3dout(iptcl,iref,irot) = sum(hadamard_prod(iptcl,:,:))
                end do
                !$omp end do
                ! normalise correlations
                !$omp workshare
                corrmat3dout(:,iref,irot) =&
                corrmat3dout(:,iref,irot)/&
                sqrt(self%sqsums_refs(:)*self%sqsums_ptcls(d1lims(1):d1lims(2)))
                !$omp end workshare nowait
                !$omp end parallel
            end do
        end do
    end subroutine gencorrs_all_cpu
    
    !>  \brief  tester routine for generating all rotational correlations for all nptcls particles and all nrefs references
    subroutine gencorrs_all_tester_1( self, corrmat2dout, inplmat2dout )
        class(polarft_corrcalc), intent(inout) :: self               !< instance
        real,    intent(out) :: corrmat2dout(self%nptcls,self%nrefs) !< output correlation matrix
        integer, intent(out) :: inplmat2dout(self%nptcls,self%nrefs) !< output inplane rot index matrix
        real    :: cc(self%nrots)
        integer :: iref, iptcl, loc(1)
        do iptcl=self%pfromto(1),self%pfromto(2)
            do iref=1,self%nrefs
                cc = self%gencorrs(iref, iptcl)
                loc = maxloc(cc)
                inplmat2dout(iptcl,iref) = loc(1)
                corrmat2dout(iptcl,iref) = cc(loc(1))
            end do
        end do
    end subroutine gencorrs_all_tester_1
    
    !>  \brief  tester routine for generating all rotational correlations for all nptcls particles and all nrefs references
    subroutine gencorrs_all_tester_2( self, corrmat3dout )
        class(polarft_corrcalc), intent(inout) :: self                       !< instance
        real, intent(out) :: corrmat3dout(self%nptcls,self%nrefs,self%nrots) !< output correlation matrix
        integer :: iref, iptcl
        do iptcl=self%pfromto(1),self%pfromto(2)
            do iref=1,self%nrefs
                corrmat3dout(iptcl,iref,:) = self%gencorrs(iref, iptcl)
            end do
        end do
    end subroutine gencorrs_all_tester_2

    !>  \brief  is for generating rotational correlations
    function gencorrs( self, iref, iptcl ) result( cc )
        !$ use omp_lib
        !$ use omp_lib_kinds
        class(polarft_corrcalc), intent(inout) :: self        !< instance
        integer, intent(in)                    :: iref, iptcl !< ref & ptcl indices
        real    :: cc(self%nrots)
        integer :: irot
        !$omp parallel do default(shared) private(irot)
        do irot=1,self%nrots
            cc(irot) = self%corr_1(iref, iptcl, irot)
        end do
        !$omp end parallel do
    end function gencorrs

    !>  \brief  for calculating the correlation between reference iref and particle iptcl in rotation irot
    function corr_1( self, iref, iptcl, irot ) result( cc )
        class(polarft_corrcalc), intent(inout) :: self              !< instance
        integer, intent(in)                    :: iref, iptcl, irot !< reference, particle, rotation
        real    :: cc, ptcl_sqsum_tmp, ref_sqsum_tmp, ptcl_mean_tmp, ref_mean_tmp
        integer :: ptcl_ind, kcount, k, kind
        ptcl_ind = self%ptcl_ind(iptcl)
        if( allocated(self%pfts_refs_ctf) )then
           cc = sum( real(self%pfts_refs_ctf(iref,:,:) * &
                     conjg(  self%pfts_ptcls(ptcl_ind,irot:irot+self%refsz-1,:))))
        else if(self%xfel)then
            !
            ! Previous XFEL correlation
            !
            !cc = sum(real(self%pfts_refs(iref,:,:))*real(self%pfts_ptcls(ptcl_ind,irot:irot+self%refsz-1,:)))
            !
            ! (Corrected) Normalisation terms for standard correlation
            !
            !ptcl_sqsum_tmp = sum(real(self%pfts_ptcls(ptcl_ind,irot:irot+self%refsz-1,:)) &
            !*real(self%pfts_ptcls(ptcl_ind,irot:irot+self%refsz-1,:)))
            !ref_sqsum_tmp = sum(real(self%pfts_refs(iref,:,:))*real(self%pfts_refs(iref,:,:)))
            !print *, 'unnorm cc:', cc
            !
            cc = 0.
            kcount = 0
            ! Shell-by-shell xfel correlation
            do k=self%kfromto(1),self%kfromto(2)
               kind = self%k_ind(k)
               ptcl_sqsum_tmp = sum(real(self%pfts_ptcls(ptcl_ind,irot:irot+self%refsz-1,kind)) &
                    *real(self%pfts_ptcls(ptcl_ind,irot:irot+self%refsz-1,kind)))
               ref_sqsum_tmp = sum(real(self%pfts_refs(iref,:,kind))*real(self%pfts_refs(iref,:,kind)))
               cc = cc + sum(real(self%pfts_refs(iref,:,kind)) &
                    *real(self%pfts_ptcls(ptcl_ind,irot:irot+self%refsz-1,kind))) &
                    / sqrt(ref_sqsum_tmp * ptcl_sqsum_tmp)
               kcount = kcount + 1
            end do
            cc = cc / kcount
        else
           cc = sum( real(  self%pfts_refs(iref,:,:) * &
                     conjg(self%pfts_ptcls(ptcl_ind,irot:irot+self%refsz-1,:))))
        endif
        if( self%xfel) then
            ! 
            ! Previous prime correlation (xfel version)
            !
            ! cc = cc/sqrt(ref_sqsum_tmp * ptcl_sqsum_tmp)
        else 
            if( self%sqsums_refs(iref)<TINY .or. self%sqsums_ptcls(ptcl_ind)<TINY )then
                cc = 0.
            else
                cc = cc/sqrt(self%sqsums_refs(iref)*self%sqsums_ptcls(ptcl_ind))
            endif
        endif
    end function corr_1
    
    !>  \brief  for calculating the on-fly shifted correlation between reference iref and particle iptcl in rotation irot
    function corr_2( self, iref, iptcl, irot, shvec ) result( cc )
        !$ use omp_lib
        !$ use omp_lib_kinds
        use simple_math, only: csq
        class(polarft_corrcalc), intent(inout) :: self !< instance
        integer, intent(in) :: iref, iptcl, irot       !< reference, particle, rotation
        real, intent(in)    :: shvec(2)                !< origin shift vector
        real    :: argmat(self%refsz,self%nk), sqsum_ref_sh, cc
        complex :: shmat(self%refsz,self%nk)
        complex :: pft_ref_sh(self%refsz,self%nk)
        integer :: ptcl_ind
        ptcl_ind = self%ptcl_ind(iptcl)
        if( allocated(self%pfts_refs_ctf) )then
            !$omp parallel workshare 
            ! generate the argument matrix from memoized components in argtransf
            argmat = self%argtransf(:self%refsz,:)*shvec(1)+self%argtransf(self%refsz+1:,:)*shvec(2)
            ! generate the complex shift transformation matrix
            shmat = cmplx(cos(argmat),sin(argmat))
            ! shift
            pft_ref_sh = self%pfts_refs_ctf(iref,:,:)*shmat
            ! calculate correlation precursors
            argmat = real(pft_ref_sh*conjg(self%pfts_ptcls(ptcl_ind,irot:irot+self%refsz-1,:)))
            !$omp end parallel workshare
            cc = sum(argmat)
            sqsum_ref_sh = sum(csq(pft_ref_sh))
        else        
            !$omp parallel workshare
            ! generate the argument matrix from memoized components in argtransf
            argmat = self%argtransf(:self%refsz,:)*shvec(1)+self%argtransf(self%refsz+1:,:)*shvec(2)
            ! generate the complex shift transformation matrix
            shmat = cmplx(cos(argmat),sin(argmat))
            ! shift
            pft_ref_sh = self%pfts_refs(iref,:,:)*shmat
            ! calculate correlation precursors
            argmat = real(pft_ref_sh*conjg(self%pfts_ptcls(ptcl_ind,irot:irot+self%refsz-1,:)))
            !$omp end parallel workshare
            cc = sum(argmat)
            sqsum_ref_sh = sum(csq(pft_ref_sh))
        endif
        cc = cc/sqrt(sqsum_ref_sh*self%sqsums_ptcls(ptcl_ind))
    end function corr_2

    !>  \brief  is for in-plane rotation search using stochastic hill climbing (shc)
    !!          we exit as soon as a better one has been identified
    subroutine rosrch_shc( self, iref, iptcl, ccinout, inplinout )
        class(polarft_corrcalc), intent(inout) :: self        !< instance
        integer, intent(in)                    :: iref, iptcl !< ref & ptcl indices
        real, intent(inout)                    :: ccinout     !< in/out because inval serves as treshold
        integer, intent(inout)                 :: inplinout   !< in/out because we may want to preserve input value
        type(ran_tabu) :: rt          
        integer        :: srch_order(self%nrots), irot
        real           :: cc
        ! make random rotation index order
        rt = ran_tabu(self%nrots)
        call rt%ne_ran_iarr(srch_order)
        do irot=1,self%nrots
            ! calculate the cross correlation coefficient
            cc = self%corr_1(iref, iptcl, srch_order(irot))
            if( cc > ccinout )then
                ! update the parameters
                ccinout   = cc
                inplinout = srch_order(irot)
                call rt%kill
                return
            endif
        end do
        call rt%kill
    end subroutine rosrch_shc
    
    ! DESTRUCTOR

    !>  \brief  is a destructor
    subroutine kill( self )
        class(polarft_corrcalc), intent(inout) :: self
        if( self%existence )then
            deallocate( self%sqsums_refs,  &
                        self%sqsums_ptcls, &
                        self%angtab,       &
                        self%argtransf,    &
                        self%polar,        &
                        self%pfts_refs,    &
                        self%pfts_ptcls    )
            if( allocated(self%pfts_refs_ctf)  ) deallocate(self%pfts_refs_ctf)
            if( allocated(self%pfts_ptcls_ctf) ) deallocate(self%pfts_ptcls_ctf)
            self%existence = .false.
        endif
    end subroutine kill

    ! UNIT TEST

    !>  \brief  is a unit test for polarft
    subroutine test_polarft_corrcalc
        use simple_polarft, only: polarft
        use simple_image,   only: image
        use simple_stat,    only: pearsn
        use simple_math,    only: euclid
        type(polarft)          :: pfts1(3,5), pfts2(3,5)
        type(polarft_corrcalc) :: pftcc
        integer, parameter     :: NITER=3000
        real                   :: corr1, corr2, diff, lscape_dist, lscape_corr
        integer                :: i, j, nrots, r
        real, allocatable      :: corrs1(:), corrs2(:)
        complex, allocatable   :: pft(:,:)
        write(*,'(a)') '**info(simple_polarft_corrcalc_unit_test): testing all functionality'
        do i=1,3
            do j=1,5
                call pfts1(i,j)%new([2,11], 45, ptcl=.false.) ! implemented as reference
                call pfts2(i,j)%new([2,11], 45, ptcl=.true.)  ! implemented as particle
                call pfts1(i,j)%set_ldim([100,100,1])
                call pfts2(i,j)%set_ldim([100,100,1])
            end do
        end do
        nrots = pfts1(1,1)%get_nrots()
        allocate(corrs1(nrots), corrs2(nrots))
        ! make a conforming polarft_corrcalc object of 3 refs and 5 ptcls
        call pftcc%new(3, [1,5], [100,100,1], [2,11], 45, 'no')
        do i=1,3
            do j=1,5
                call pfts1(i,j)%noise_fill
                call pfts2(i,j)%noise_fill
                call pfts1(i,j)%memoize_sqsums
                call pfts2(i,j)%memoize_sqsums
                ! pfts1 is ref
                pft = pfts1(i,j)%get_pft()
                call pftcc%set_ref_pft(i, pft)
                deallocate(pft)
                ! pfts2 is ptcl
                pft = pfts2(i,j)%get_pft()
                call pftcc%set_ptcl_pft(j, pft)
                deallocate(pft)
            end do
        end do
        write(*,'(a)') '**info(simple_polarft_corrcalc_unit_test, part 1):&
        &comparing the gencorrs correlations with the ones generated here'
        do i=1,3
            do j=1,5
                corrs1 = pftcc%gencorrs(i,j)
                do r=1,nrots
                    corrs2(r) = pftcc%corr(i, j, r)
                end do
                lscape_corr = pearsn(corrs1,corrs2)
                lscape_dist = euclid(corrs1,corrs2)
                if( lscape_dist > 1e-6 .or. lscape_corr < 0.9999 )then
                    write(*,*) 'landscape correlation: ', lscape_corr
                    write(*,*) 'landscape dist: ', lscape_dist
                    stop 'corr and gencorrs not consistent in pftcc'
                endif
            end do
        end do
        write(*,'(a)') '**info(simple_polarft_corrcalc_unit_test, part 2):&
        &testing the gencorrs function'
        do i=1,3
            do j=1,5
                call pfts1(i,j)%noise_fill
                call pfts2(i,j)%noise_fill
                call pfts1(i,j)%memoize_sqsums
                call pfts2(i,j)%memoize_sqsums
                ! pfts1 is ref
                pft = pfts1(i,j)%get_pft()
                call pftcc%set_ref_pft(i, pft)
                deallocate(pft)
                ! pfts2(i,j) is ptcl
                pft = pfts2(i,j)%get_pft()
                call pftcc%set_ptcl_pft(j, pft)
                deallocate(pft)
                ! compare the two correlation functions
                diff = 0.
                call pfts1(i,j)%gencorrs(pfts2(i,j), corrs1)
                corrs2 = pftcc%gencorrs(i,j)   
                diff = sum(abs(corrs1-corrs2))/real(nrots*3*5)
                if( diff > 1.5e-2 )then
                    print *, 'diff = ', diff
                    stop 'the new gencorrs function does not conform with the old; test_polarft_corrcalc'
                endif
            end do
        end do
        do i=1,3
            do j=1,5
                call pfts1(i,j)%noise_fill
                call pfts2(i,j)%noise_fill
                call pfts1(i,j)%memoize_sqsums
                call pfts2(i,j)%memoize_sqsums
                ! pfts1 is ref
                pft = pfts1(i,j)%get_pft()
                call pftcc%set_ref_pft(i, pft)
                deallocate(pft)
                ! pfts2(i,j) is ptcl
                pft = pfts2(i,j)%get_pft()
                call pftcc%set_ptcl_pft(j, pft)
                deallocate(pft)
            end do
        end do
        write(*,'(a)') '**info(simple_polarft_corrcalc_unit_test, part 3):&
        &testing the gencorrs function with expand_dim'
        call pftcc%expand_dim
        do i=1,3
            do j=1,5
                ! compare the two correlation functions
                diff = 0.
                call pfts1(i,j)%gencorrs(pfts2(i,j), corrs1)
                corrs2 = pftcc%gencorrs(i,j)   
                diff = sum(abs(corrs1-corrs2))/real(nrots*3*5)
                if( diff > 1.5e-2 )then
                    print *, 'diff = ', diff
                    stop 'the new  gencorrs function with expand_dim does not&
                    &conform with the old; test_polarft_corrcalc'
                endif
            end do
        end do
        write(*,'(a)') '**info(simple_polarft_corrcalc_unit_test, part 4):&
        &testing the corr function with expand_dim'
        do i=1,3
            do j=1,5
                ! compare the two correlation functions
                diff = 0.
                do r=1,nrots
                    corrs1(r) = pfts1(i,j)%corr(pfts2(i,j), r)
                    corrs2(r) = pftcc%corr(i, j, r) ! the newly implemented one
                end do
                diff = sum(abs(corrs1-corrs2))/real(nrots*3*5)
                if( diff > 1.5e-2 )then
                    print *, 'diff = ', diff
                    stop 'the new correlation function does not conform with the old; test_polarft_corrcalc'
                endif
            end do
        end do
        deallocate(corrs1, corrs2)
        write(*,'(a)') 'SIMPLE_POLARFT_CORRCALC_UNIT_TEST COMPLETED SUCCESSFULLY ;-)'
    end subroutine test_polarft_corrcalc
    
end module simple_polarft_corrcalc
