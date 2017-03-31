!>  \brief  SIMPLE polarft_corrcalc class
module simple_polarft_corrcalc
use simple_defs      ! use all in there
use simple_params,   only: params
use simple_ran_tabu, only: ran_tabu
use simple_jiffys,   only: alloc_err
implicit none

public :: polarft_corrcalc
private

! CLASS PARAMETERS/VARIABLES
complex, parameter :: zero=cmplx(0.,0.) !< just a complex zero
logical, parameter :: DEBUG = .true.

type :: polarft_corrcalc
    private
    integer                  :: pfromto(2) = 1         !< from/to particle indices (in parallel execution)
    integer                  :: nptcls     = 1         !< the total number of particles in partition (logically indexded [fromp,top])
    integer                  :: nrefs      = 1         !< the number of references (logically indexded [1,nrefs])
    integer                  :: nrots      = 0         !< number of in-plane rotations for one pft (determined by radius of molecule)
    integer                  :: ring2      = 0         !< radius of molecule
    integer                  :: refsz      = 0         !< size of reference (nrots/2) (number of vectors used for matching)
    integer                  :: ptclsz     = 0         !< size of particle (2*nrots)
    integer                  :: winsz      = 0         !< size of moving window in correlation cacluations
    integer                  :: ldim(3)    = 0         !< logical dimensions of original cartesian image
    integer                  :: kfromto(2) = 0         !< Fourier index range
    integer                  :: nthr       = 0         !< nr of CPU threads
    real,        allocatable :: sqsums_refs(:)         !< memoized square sums for the correlation calculations
    real,        allocatable :: sqsums_ptcls(:)        !< memoized square sums for the correlation calculations
    real,        allocatable :: angtab(:)              !< table of in-plane angles (in degrees)
    real,        allocatable :: argtransf(:,:)         !< argument transfer constants for shifting the references
    real,        allocatable :: polar(:,:)             !< table of polar coordinates (in Cartesian coordinates)
    real,        allocatable :: ctfmats(:,:,:)         !< expandd set of CTF matrices (for efficient parallel exec)
    complex(sp), allocatable :: pfts_refs(:,:,:)       !< 3D complex matrix of polar reference sections (nrefs,refsz,nk)
    complex(sp), allocatable :: pfts_refs_ctf(:,:,:)   !< 3D complex matrix of polar reference sections with CTF applied
    complex(sp), allocatable :: pfts_ptcls(:,:,:)      !< 3D complex matrix of particle sections
    logical                  :: with_ctf     = .false. !< CTF flag
    logical                  :: xfel         = .false. !< to indicate whether we process xfel patterns or not
    logical                  :: dim_expanded = .false. !< to indicate whether dim has been expanded or not
    logical                  :: existence    = .false. !< to indicate existence
  contains
    ! CONSTRUCTOR
    procedure          :: new
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
    procedure          :: get_pdim
    procedure          :: get_rot
    procedure          :: get_roind
    procedure          :: get_win_roind
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
    ! I/O
    procedure          :: write_pfts_ptcls
    procedure          :: read_pfts_ptcls
    ! MODIFIERS
    procedure, private :: apply_ctf_1
    procedure, private :: apply_ctf_2
    generic            :: apply_ctf => apply_ctf_1, apply_ctf_2
    procedure          :: apply_ctf_single
    procedure          :: xfel_subtract_shell_mean
    ! CALCULATORS
    procedure          :: create_polar_ctfmat
    procedure          :: create_polar_ctfmats
    procedure, private :: gencorrs_all_cpu_1
    procedure, private :: gencorrs_all_cpu_2
    generic            :: gencorrs_all_cpu => gencorrs_all_cpu_1, gencorrs_all_cpu_2
    procedure          :: gencorrs_serial
    procedure          :: gencorrs
    procedure          :: genfrc
    procedure          :: corrs
    procedure, private :: corr_1
    procedure, private :: corr_2
    generic            :: corr => corr_1, corr_2
    ! DESTRUCTOR
    procedure          :: kill
end type polarft_corrcalc

contains
    
    ! CONSTRUCTORS
    
    !>  \brief  is a constructor
    subroutine new( self, nrefs, pfromto, ldim, kfromto, ring2, nthr, ctfflag, isxfel )
        use simple_math, only: rad2deg, is_even, round2even
        class(polarft_corrcalc),    intent(inout) :: self
        integer,                    intent(in)    :: nrefs, pfromto(2), ldim(3), kfromto(2), ring2, nthr
        character(len=*),           intent(in)    :: ctfflag
        character(len=*), optional, intent(in)    :: isxfel
        integer :: alloc_stat, irot, k, err
        logical :: even_dims, test(3)
        real    :: ang
        ! kill possibly pre-existing object
        call self%kill
        ! error check
        if( kfromto(2) - kfromto(1) <= 2 )then
            write(*,*) 'kfromto: ', kfromto(1), kfromto(2)
            stop 'resolution range too narrow; new; simple_polarft_corrcalc'
        endif
        if( ring2 < 1 )then
            write(*,*) 'ring2: ', ring2
            stop 'ring2 must be > 0; new; simple_polarft_corrcalc'
        endif
        if( pfromto(2) - pfromto(1) + 1 < 1 )then
            write(*,*) 'pfromto: ', pfromto(1), pfromto(2)
            stop 'nptcls (# of particles) must be > 0; new; simple_polarft_corrcalc'
        endif
        if( nrefs < 1 )then
            write(*,*) 'nrefs: ', nrefs
            stop 'nrefs (# of reference sections) must be > 0; new; simple_polarft_corrcalc'
        endif
        if( any(ldim == 0) )then
            write(*,*) 'ldim: ', ldim
            stop 'ldim is not conforming (is zero); new; simple_polarft_corrcalc'
        endif
        if( ldim(3) > 1 )then
            write(*,*) 'ldim: ', ldim
            stop '3D polarfts are not yet supported; new; simple_polarft_corrcalc'
        endif
        test    = .false.
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
        self%pfromto = pfromto                         !< from/to particle indices (in parallel execution)
        self%nptcls  = pfromto(2) - pfromto(1) + 1     !< the total number of particles in partition (logically indexded [fromp,top])
        self%nrefs   = nrefs                           !< the number of references (logically indexded [1,nrefs])
        self%ring2   = ring2                           !< radius of molecule
        self%nrots   = round2even(twopi * real(ring2)) !< number of in-plane rotations for one pft  (determined by radius of molecule)
        self%refsz   = self%nrots / 2                  !< size of reference (nrots/2) (number of vectors used for matching)
        self%winsz   = self%refsz - 1                  !< size of moving window in correlation cacluations
        self%ptclsz  = self%nrots * 2                  !< size of particle (2*nrots)
        self%ldim    = ldim                            !< logical dimensions of original cartesian image
        self%kfromto = kfromto                         !< Fourier index range
        self%nthr    = nthr                            !< number of CPU threads
        ! generate polar coordinates
        allocate( self%polar(self%ptclsz,self%kfromto(1):self%kfromto(2)), self%angtab(self%nrots), stat=alloc_stat)
        call alloc_err('polar coordinate arrays; new; simple_polarft_corrcalc', alloc_stat)
        ang = twopi/real(self%nrots)
        do irot=1,self%nrots
            self%angtab(irot) = real(irot-1)*ang
            do k=self%kfromto(1),self%kfromto(2)
                self%polar(irot,k)            = cos(self%angtab(irot))*real(k) ! x-coordinate
                self%polar(irot+self%nrots,k) = sin(self%angtab(irot))*real(k) ! y-coordinate
            end do
            self%angtab(irot) = rad2deg(self%angtab(irot)) ! angle (in degrees)
        end do
        ! generate the argument transfer constants for shifting reference polarfts
        allocate( self%argtransf(self%nrots,self%kfromto(1):self%kfromto(2)), stat=alloc_stat)
        call alloc_err('shift argument transfer array; new; simple_polarft_corrcalc', alloc_stat)
        self%argtransf(:self%refsz,:)   = &
            self%polar(:self%refsz,:)   * &
            (PI/real(self%ldim(1)/2))    ! x-part
        self%argtransf(self%refsz+1:,:) = &
            self%polar(self%nrots+1:self%nrots+self%refsz,:) * &
            (PI/real(self%ldim(2)/2))    ! y-part
        ! allocate polarfts and sqsums
        allocate(   self%pfts_refs(self%nrefs,self%refsz,self%kfromto(1):self%kfromto(2)),&
                    self%pfts_ptcls(self%pfromto(1):self%pfromto(2),self%ptclsz,self%kfromto(1):self%kfromto(2)),&
                    self%sqsums_refs(self%nrefs),&
                    self%sqsums_ptcls(self%pfromto(1):self%pfromto(2)),&
                    self%pfts_refs_ctf(self%nrefs,self%refsz,self%kfromto(1):self%kfromto(2)), stat=alloc_stat)
        call alloc_err('polarfts and sqsums; new; simple_polarft_corrcalc', alloc_stat)
        self%pfts_refs     = zero
        self%pfts_ptcls    = zero
        self%sqsums_refs   = 0.
        self%sqsums_ptcls  = 0.
        self%pfts_refs_ctf = zero
        ! pfts_refs_ctf if needed
        self%with_ctf = .false.
        if( ctfflag .ne. 'no' ) self%with_ctf = .true.
        self%dim_expanded = .false.
        self%existence    = .true.
    end subroutine new

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
        self%pfts_ptcls(iptcl,:self%nrots,:)   = pft
        self%pfts_ptcls(iptcl,self%nrots+1:,:) = pft ! because rot dim is expanded
        ! calculate the square sum required for correlation calculation
        call self%memoize_sqsum_ptcl(iptcl)
    end subroutine set_ptcl_pft

    !>  \brief  sets a reference Fourier component
    subroutine set_ref_fcomp( self, iref, irot, k, comp )
        class(polarft_corrcalc), intent(inout) :: self
        integer,                 intent(in)    :: iref, irot, k
        complex,                 intent(in)    :: comp
        self%pfts_refs(iref,irot,k) = comp
    end subroutine set_ref_fcomp
    
    !>  \brief  sets a particle Fourier component
    subroutine set_ptcl_fcomp( self, iptcl, irot, k, comp )
        class(polarft_corrcalc), intent(inout) :: self
        integer,                 intent(in)    :: iptcl, irot, k
        complex,                 intent(in)    :: comp
        self%pfts_ptcls(iptcl,irot,k) = comp
        self%pfts_ptcls(iptcl,irot+self%nrots,k) = comp ! because rot dim is expanded
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

    !>  \brief  for getting the dimensions of the reference polar FT
    function get_pdim( self, isptcl ) result( pdim )
        class(polarft_corrcalc), intent(in) :: self
        logical,                 intent(in) :: isptcl
        integer :: pdim(3)
        if( isptcl )then
            pdim = [self%nrots,self%kfromto(1),self%kfromto(2)]
        else
            pdim = [self%refsz,self%kfromto(1),self%kfromto(2)]
        endif
    end function get_pdim
    
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
        integer :: ind, irot, loc(1)
        real    :: dists(self%nrots)
        dists = abs(self%angtab-rot)
        where(dists>180.)dists = 360.-dists
        loc = minloc(dists)
        ind = loc(1)
    end function get_roind

    !>  \brief is for getting the discrete in-plane rotational
    !!         indices within a window of +-winsz degrees of ang
    !!         For use together with gencorrs
    function get_win_roind( self, ang, winsz )result( roind_vec )
        use simple_math, only: rad2deg
        class(polarft_corrcalc), intent(in) :: self
        real,                    intent(in) :: ang, winsz
        integer, allocatable :: roind_vec(:)
        real    :: dist(self%nrots)
        integer :: i, irot, nrots, alloc_stat
        if(ang>360. .or. ang<TINY)stop 'input angle outside of the conventional range; simple_polarft_corrcalc::get_win_roind'
        if(winsz<0. .or. winsz>180.)stop 'invalid window size; simple_polarft_corrcalc::get_win_roind'
        if(winsz < 360./real(self%nrots))stop 'too small window size; simple_polarft_corrcalc::get_win_roind'
        i    = self%get_roind( ang )
        dist = abs(self%angtab(i) - self%angtab)
        where( dist>180. )dist = 360.-dist
        nrots = count(dist <= winsz)
        allocate( roind_vec(nrots), stat=alloc_stat )
        call alloc_err("In: get_win_roind; simple_polarft_corrcalc", alloc_stat)
        irot = 0
        do i = 1,self%nrots
            if( dist(i)<=winsz )then
                irot = irot+1
                roind_vec(irot) = i
            endif
        enddo
    end function get_win_roind

    !>  \brief returns polar coordinate for rotation rot
    !!         and Fourier index k
    function get_coord( self, rot, k ) result( xy )
        class(polarft_corrcalc), intent(in) :: self
        integer,                 intent(in) :: rot, k
        real :: xy(2)
        xy(1) = self%polar(rot,k)
        xy(2) = self%polar(self%nrots+rot,k)
    end function get_coord
    
    !>  \brief  returns polar Fourier transform of particle iptcl in rotation irot
    function get_ptcl_pft( self, iptcl, irot ) result( pft )
        class(polarft_corrcalc), intent(in) :: self
        integer,                 intent(in) :: iptcl, irot
        complex, allocatable :: pft(:,:)
        integer :: alloc_stat
        allocate(pft(self%refsz,self%kfromto(1):self%kfromto(2)),&
        source=self%pfts_ptcls(iptcl,irot:irot+self%winsz,:), stat=alloc_stat)
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
        call gnufor_image(real(self%pfts_ptcls(iptcl,:self%refsz,:)),  palette='gray')
        call gnufor_image(aimag(self%pfts_ptcls(iptcl,:self%refsz,:)), palette='gray')
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
        write(*,*) "radius of molecule                      (self%ring2): ", self%ring2
        write(*,*) "nr of rots for ref (2nd dim of pftmat)  (self%refsz): ", self%refsz
        write(*,*) "n rots for ptcl (2nd dim of pftmat)    (self%ptclsz): ", self%ptclsz 
        write(*,*) "logical dim. of original Cartesian image (self%ldim): ", self%ldim
        write(*,*) "high-pass limit Fourier index      (self%kfromto(1)): ", self%kfromto(1)
        write(*,*) "low-pass limit Fourier index       (self%kfromto(2)): ", self%kfromto(2)
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
    subroutine memoize_sqsum_ref_ctf( self, iref, irot )
        use simple_math, only: csq
        class(polarft_corrcalc), intent(inout) :: self
        integer,                 intent(in)    :: iref
        integer,       optional, intent(in)    :: irot
        if( present(irot) )then
            self%sqsums_refs(iref) = sum(csq(self%pfts_refs_ctf(iref,irot:irot+self%winsz,:)))
        else
            self%sqsums_refs(iref) = sum(csq(self%pfts_refs_ctf(iref,:,:)))
        endif
    end subroutine memoize_sqsum_ref_ctf

    !>  \brief  is for memoization of the complex square sums reqired for correlation calculation
    subroutine memoize_sqsum_ptcl( self, iptcl )
        use simple_math, only: csq
        class(polarft_corrcalc), intent(inout) :: self
        integer,                 intent(in)    :: iptcl
        if( self%xfel )then
            self%sqsums_ptcls(iptcl) = sum(real(self%pfts_ptcls(iptcl,:self%refsz,:))**2.)
        else
            self%sqsums_ptcls(iptcl) = sum(csq(self%pfts_ptcls(iptcl,:self%refsz,:)))
        endif
    end subroutine memoize_sqsum_ptcl

    ! I/O

    !>  \brief  is for writing particle pfts to file
    subroutine write_pfts_ptcls( self, fname )
        use simple_filehandling, only: get_fileunit
        class(polarft_corrcalc), intent(in) :: self
        character(len=*),        intent(in) :: fname
        integer :: funit, io_stat
        funit = get_fileunit()
        open(unit=funit, status='REPLACE', action='WRITE', file=trim(fname), access='STREAM')
        write(unit=funit,pos=1,iostat=io_stat) self%pfts_ptcls
        ! Check if the write was successful
        if( io_stat .ne. 0 )then
            write(*,'(a,i0,2a)') '**ERROR(simple_polarft_corrcalc): I/O error ', io_stat, ' when writing file: ', trim(fname)
            stop 'I/O error; write_pfts_ptcls'
        endif
        close(funit)
    end subroutine write_pfts_ptcls

    !>  \brief  is for reading particle pfts from file
    subroutine read_pfts_ptcls( self, fname )
        use simple_filehandling, only: get_fileunit
        class(polarft_corrcalc), intent(inout) :: self
        character(len=*),        intent(in)    :: fname
        integer :: funit, io_stat, iptcl
        funit = get_fileunit()
        open(unit=funit, status='OLD', action='READ', file=trim(fname), access='STREAM')
        read(unit=funit,pos=1,iostat=io_stat) self%pfts_ptcls
        ! Check if the read was successful
        if( io_stat .ne. 0 )then
            write(*,'(a,i0,2a)') '**ERROR(simple_polarft_corrcalc): I/O error ', io_stat, ' when reading file: ', trim(fname)
            stop 'I/O error; read_pfts_ptcls'
        endif
        close(funit)
        ! memoize sqsum_ptcls
        !$omp parallel do schedule(auto) default(shared) private(iptcl)
        do iptcl=self%pfromto(1),self%pfromto(2)
            call self%memoize_sqsum_ptcl(iptcl)
        end do
        !$omp end parallel do
    end subroutine read_pfts_ptcls
    
    ! MODIFIERS

    !>  \brief  is for applying CTF to references and updating the memoized ref sqsums
    subroutine apply_ctf_1( self, tfun, dfx, dfy, angast, refvec, ctfmat )
        !$ use omp_lib
        !$ use omp_lib_kinds
        use simple_ctf,   only: ctf
        class(polarft_corrcalc), intent(inout) :: self
        class(ctf),              intent(inout) :: tfun
        real,                    intent(in)    :: dfx
        real,    optional,       intent(in)    :: dfy, angast
        integer, optional,       intent(in)    :: refvec(2)
        real,    optional,       intent(in)    :: ctfmat(:,:)
        real, allocatable :: ctfmat_here(:,:)
        integer           :: iref, ref_start, ref_end
        if( present(ctfmat) )then
            ctfmat_here = ctfmat
        else
            ! create the congruent polar matrix of real CTF values
            ctfmat_here = self%create_polar_ctfmat(tfun, dfx, dfy, angast, self%refsz)
        endif
        ! multiply the references with the CTF
        if( present(refvec) )then
            ! slice of references
            if( any(refvec<1) .or. any(refvec>self%nrefs) .or. refvec(1)>refvec(2) )then
                stop 'invalid reference indices; simple_polarft_corrcalc::apply_ctf:'
            endif
            ref_start = refvec(1)
            ref_end   = refvec(2)
        else
            ! all references
            ref_start = 1
            ref_end   = self%nrefs
        endif
        !$omp parallel do default(shared) schedule(auto) private(iref)
        do iref=ref_start,ref_end
            self%pfts_refs_ctf(iref,:,:) = self%pfts_refs(iref,:,:)*ctfmat_here
            call self%memoize_sqsum_ref_ctf(iref)
        end do
        !$omp end parallel do
        if( allocated(ctfmat_here) )deallocate(ctfmat_here)
    end subroutine apply_ctf_1

    !>  \brief  is for applying CTF to references and updating the memoized ref sqsums
    subroutine apply_ctf_2( self, iptcl, refvec )
        !$ use omp_lib
        !$ use omp_lib_kinds
        class(polarft_corrcalc), intent(inout) :: self
        integer,                 intent(in)    :: iptcl
        integer, optional,       intent(in)    :: refvec(2)
        integer :: iref, ref_start, ref_end
        if( self%with_ctf )then
            ! multiply the references with the CTF
            if( present(refvec) )then
                ! slice of references
                if( any(refvec<1) .or. any(refvec>self%nrefs) .or. refvec(1)>refvec(2) )then
                    stop 'invalid reference indices; simple_polarft_corrcalc::apply_ctf_2'
                endif
                ref_start = refvec(1)
                ref_end   = refvec(2)
            else
                ! all references
                ref_start = 1
                ref_end   = self%nrefs
            endif
            !$omp parallel do default(shared) schedule(auto) private(iref)
            do iref=ref_start,ref_end
                self%pfts_refs_ctf(iref,:,:) = self%pfts_refs(iref,:,:)*self%ctfmats(iptcl,:,:)
                call self%memoize_sqsum_ref_ctf(iref)
            end do
            !$omp end parallel do
        endif
    end subroutine apply_ctf_2

    !>  \brief  is for applying CTF to references and updating the memoized ref sqsums
    subroutine apply_ctf_single( self, iptcl, iref )
        class(polarft_corrcalc), intent(inout) :: self
        integer,                 intent(in)    :: iptcl, iref
        self%pfts_refs_ctf(iref,:,:) = self%pfts_refs(iref,:,:)*self%ctfmats(iptcl,:,:)
        call self%memoize_sqsum_ref_ctf(iref)
    end subroutine apply_ctf_single

    !>  \brief  is for preparing for XFEL pattern corr calc
    subroutine xfel_subtract_shell_mean( self )
        class(polarft_corrcalc), intent(inout) :: self
        real, allocatable    :: ptcls_mean_tmp(:,:,:)
        real, allocatable    :: refs_mean_tmp(:,:)
        integer :: iptcl, iref, irot, k
        allocate( ptcls_mean_tmp(2*self%nptcls,self%ptclsz,self%kfromto(1):self%kfromto(2)),&
        refs_mean_tmp(self%nrefs,self%kfromto(1):self%kfromto(2)))
        ! calculate the mean of each reference at each k shell
        do iref=1,self%nrefs
            do k=self%kfromto(1),self%kfromto(2)
                refs_mean_tmp(iref,k) = sum(self%pfts_refs(iref,:,k))/self%refsz
            end do
        end do
        ! calculate the mean of each reference at each k shell
        do iref=1,self%nrefs
            do irot=1,self%refsz
                do k=self%kfromto(1),self%kfromto(2)
                    self%pfts_refs(iref,irot,k) = &
                    self%pfts_refs(iref,irot,k) - refs_mean_tmp(iref,k) 
                end do
            end do
        end do
        ! calculate the mean of each particle at each k shell at each in plane rotation
        do iptcl=self%pfromto(1),self%pfromto(2)
            do k=self%kfromto(1),self%kfromto(2)
                ptcls_mean_tmp(iptcl,1,k) = sum(self%pfts_ptcls(iptcl,1:self%winsz,k))/self%refsz
            end do
        end do
        ! subtract the mean of each particle at each k shell at each in plane rotation
        do iptcl=self%pfromto(1),self%pfromto(2)
            do irot=1,self%ptclsz
                do k=self%kfromto(1),self%kfromto(2)
                    self%pfts_ptcls(iptcl,irot,k) = &
                    self%pfts_ptcls(iptcl,irot,k) - ptcls_mean_tmp(iptcl,1,k)
                end do
            end do
        end do
        deallocate( ptcls_mean_tmp, refs_mean_tmp )
    end subroutine xfel_subtract_shell_mean

    ! CALCULATORS

    !>  \brief  is for generating a matrix of CTF values
    function create_polar_ctfmat( self, tfun, dfx, dfy, angast, endrot ) result( ctfmat )
        !$ use omp_lib
        !$ use omp_lib_kinds
        use simple_ctf, only: ctf
        class(polarft_corrcalc), intent(inout) :: self
        class(ctf),              intent(inout) :: tfun
        real,                    intent(in)    :: dfx, dfy, angast
        integer,                 intent(in)    :: endrot
        real, allocatable :: ctfmat(:,:)
        real              :: inv_ldim(3),hinv,kinv,spaFreqSq,ang
        integer           :: irot,k
        allocate( ctfmat(endrot,self%kfromto(1):self%kfromto(2)) )
        inv_ldim = 1./real(self%ldim)
        !$omp parallel do collapse(2) default(shared) private(irot,k,hinv,kinv,spaFreqSq,ang) schedule(auto)
        do irot=1,endrot
            do k=self%kfromto(1),self%kfromto(2)
                hinv           = self%polar(irot,k)*inv_ldim(1)
                kinv           = self%polar(irot+self%nrots,k)*inv_ldim(2)
                spaFreqSq      = hinv*hinv+kinv*kinv
                ang            = atan2(self%polar(irot+self%nrots,k),self%polar(irot,k))
                ctfmat(irot,k) = tfun%eval(spaFreqSq,dfx,dfy,angast,ang)
            end do
        end do
        !$omp end parallel do
    end function create_polar_ctfmat

    !>  \brief  is for generating all matrices of CTF values
    subroutine create_polar_ctfmats( self, smpd, a )
        use simple_ctf,  only: ctf
        use simple_oris, only: oris
        class(polarft_corrcalc), intent(inout) :: self
        real,                    intent(in)    :: smpd
        class(oris),             intent(inout) :: a
        type(ctf) :: tfun
        integer   :: iptcl,alloc_stat 
        real      :: kv,cs,fraca,dfx,dfy,angast
        logical   :: astig
        astig = a%isthere('dfy')
        if( allocated(self%ctfmats) ) deallocate(self%ctfmats)
        allocate(self%ctfmats(self%pfromto(1):self%pfromto(2),self%refsz,self%kfromto(1):self%kfromto(2)), stat=alloc_stat)
        call alloc_err("In: simple_polarft_corrcalc :: create_polar_ctfmats, 2", alloc_stat)
        do iptcl=self%pfromto(1),self%pfromto(2)
            kv     = a%get(iptcl, 'kv'   )
            cs     = a%get(iptcl, 'cs'   )
            fraca  = a%get(iptcl, 'fraca')
            dfx    = a%get(iptcl, 'dfx'  )
            dfy    = dfx
            angast = 0.            
            if( astig )then
                dfy    = a%get(iptcl, 'dfy'   )
                angast = a%get(iptcl, 'angast')
            endif
            tfun = ctf(smpd, kv, cs, fraca)
            self%ctfmats(iptcl,:,:) = self%create_polar_ctfmat(tfun, dfx, dfy, angast, self%refsz)
        end do
    end subroutine create_polar_ctfmats

    !>  \brief  routine for generating all rotational correlations
    subroutine gencorrs_all_cpu_1( self, corrmat3dout )
        !$ use omp_lib
        !$ use omp_lib_kinds
        class(polarft_corrcalc), intent(inout) :: self
        real,                    intent(out)   :: corrmat3dout(self%pfromto(1):self%pfromto(2),self%nrefs,self%nrots)
        integer :: iptcl, iref, nptcls
        !$omp parallel default(shared) private(iref)
        do iptcl=self%pfromto(1),self%pfromto(2)
            ! tried to parallelize this one level up, which doesn't work because 
            ! then we would need one CTF modulated reference array per thread
            ! as we would otherwise get a race condition because diferent threads
            ! try to write to the same memory location
            !$omp do schedule(auto)
            do iref=1,self%nrefs
                if( self%with_ctf ) call self%apply_ctf_single(iptcl, iref)
                corrmat3dout(iptcl,iref,:) = self%gencorrs_serial(iref,iptcl)
            end do
            !$omp end do
        end do
        !$omp end parallel
    end subroutine gencorrs_all_cpu_1

    !>  \brief  routine for generating all rotational correlations
    !!          in a per-particle neighborhood structure
    !!          need to remeber to create this structure for 3D implementation
    subroutine gencorrs_all_cpu_2( self, nspace, nnn, nnmat, prevprojs, corrmat3dout )
        !$ use omp_lib
        !$ use omp_lib_kinds
        use simple_math, only: csq
        class(polarft_corrcalc), intent(inout) :: self
        integer,                 intent(in)    :: nspace, nnn, nnmat(nspace,nnn)
        integer,                 intent(in)    :: prevprojs(self%pfromto(1):self%pfromto(2))
        real,                    intent(out)   :: corrmat3dout(self%pfromto(1):self%pfromto(2),nnn,self%nrots)
        integer :: iptcl, iiref, iref
        !$omp parallel default(shared) private(iiref,iref)
        do iptcl=self%pfromto(1),self%pfromto(2)
            ! tried to parallelize this one level up, which doesn't work because 
            ! then we would need one CTF modulated reference array per thread
            ! as we would otherwise get a race condition because diferent threads
            ! try to write to the same memory location
            !$omp do schedule(auto)
            do iiref=1,nnn
                iref = nnmat(prevprojs(iptcl),iiref)
                if( self%with_ctf ) call self%apply_ctf_single(iptcl, iref)
                corrmat3dout(iptcl,iiref,:) = self%gencorrs_serial(iref,iptcl)
            end do
            !$omp end do
        end do
        !$omp end parallel
    end subroutine gencorrs_all_cpu_2

    !>  \brief  is for generating rotational correlations
    function gencorrs_serial( self, iref, iptcl, roind_vec ) result( cc )
        class(polarft_corrcalc), intent(inout) :: self        !< instance
        integer,                 intent(in)    :: iref, iptcl !< ref & ptcl indices
        integer,       optional, intent(in)    :: roind_vec(:)
        real      :: cc(self%nrots)
        integer   :: irot, i, nrots
        ! all correlations
        do irot=1,self%nrots
            cc(irot) = self%corr_1(iref, iptcl, irot)
        end do
    end function gencorrs_serial

    !>  \brief  is for generating rotational correlations
    function gencorrs( self, iref, iptcl, roind_vec ) result( cc )
        !$ use omp_lib
        !$ use omp_lib_kinds
        class(polarft_corrcalc), intent(inout) :: self        !< instance
        integer,                 intent(in)    :: iref, iptcl !< ref & ptcl indices
        integer,       optional, intent(in)    :: roind_vec(:)
        real      :: cc(self%nrots)
        integer   :: irot, i, nrots
        if( present(roind_vec) )then
            ! calculates only corrs for rotational indices provided in roind_vec
            ! see get_win_roind. returns -1.0 when not calculated
            if( any(roind_vec<=0) .or. any(roind_vec>self%nrots) )&
                &stop'index out of range; simple_polarft_corrcalc::gencorrs'
            cc    = -1.
            nrots = size(roind_vec)
            !$omp parallel do schedule(auto) default(shared) private(i,irot)
            do i=1,nrots
                irot = roind_vec(i)
                cc(irot) = self%corr_1(iref, iptcl, irot)
            end do
            !$omp end parallel do 
        else
            ! all correlations
            !$omp parallel do schedule(auto) default(shared) private(irot)
            do irot=1,self%nrots
                cc(irot) = self%corr_1(iref, iptcl, irot)
            end do
            !$omp end parallel do 
        endif
    end function gencorrs

    !>  \brief  is for generating rotational correlations
    function corrs( self, refvec, nrefs_in, iptcl, irot) result( cc )
        class(polarft_corrcalc), intent(inout) :: self                  !< instance
        integer,                 intent(in)    :: refvec(nrefs_in)      !< ref & ptcl indices
        integer,                 intent(in)    :: irot, iptcl, nrefs_in !< ref & ptcl indices
        real    :: cc(nrefs_in)
        integer :: i, iref
        cc = -1.
        do i=1,nrefs_in
            iref = refvec(i)
            if(iref>0 .and. iref<=self%nrefs )cc(i) = self%corr_1(iref, iptcl, irot)
        end do
    end function corrs

    !>  \brief  is for generating resolution dependent correlations
    function genfrc( self, iref, iptcl, irot ) result( frc )
        !$ use omp_lib
        !$ use omp_lib_kinds
        use simple_math, only: csq
        class(polarft_corrcalc), target, intent(inout) :: self              !< instance
        integer,                         intent(in)    :: iref, iptcl, irot !< reference, particle, rotation
        real, allocatable :: frc(:)
        complex, pointer  :: refs(:,:,:) => null()
        real    :: sumsqref, sumsqptcl
        integer :: k
        allocate( frc(self%kfromto(1):self%kfromto(2)) )
        if( allocated(self%pfts_refs_ctf) )then
            refs => self%pfts_refs_ctf
        else
            refs => self%pfts_refs
        endif
        !$omp parallel do default(shared) private(k,sumsqref,sumsqptcl) schedule(auto)
        do k=self%kfromto(1),self%kfromto(2)
            frc(k)    = sum(refs(iref,:,k)*conjg(self%pfts_ptcls(iptcl,irot:irot+self%winsz,k)))
            sumsqref  = sum(csq(refs(iref,:,k)))
            sumsqptcl = sum(csq(self%pfts_ptcls(iptcl,:self%refsz,k)))
            frc(k)    = frc(k)/sqrt(sumsqref*sumsqptcl)
        end do
        !$omp end parallel do
    end function genfrc

    !>  \brief  for calculating the correlation between reference iref and particle iptcl in rotation irot
    function corr_1( self, iref, iptcl, irot ) result( cc )
        class(polarft_corrcalc), intent(inout) :: self              !< instance
        integer,                 intent(in)    :: iref, iptcl, irot !< reference, particle, rotation
        real    :: cc
        if( self%sqsums_refs(iref) < TINY .or. self%sqsums_ptcls(iptcl) < TINY )then
            cc = 0.
            return
        endif
        if( self%with_ctf )then
            cc = sum(real(self%pfts_refs_ctf(iref,:,:) * conjg(self%pfts_ptcls(iptcl,irot:irot+self%winsz,:))))
        else
            cc = sum(real(self%pfts_refs(iref,:,:) * conjg(self%pfts_ptcls(iptcl,irot:irot+self%winsz,:))))
        endif
        cc = cc/sqrt(self%sqsums_refs(iref)*self%sqsums_ptcls(iptcl))
    end function corr_1

    !>  \brief  for calculating the on-fly shifted correlation between reference iref and particle iptcl in rotation irot
    function corr_2( self, iref, iptcl, irot, shvec ) result( cc )
        !$ use omp_lib
        !$ use omp_lib_kinds
        use simple_math, only: csq
        class(polarft_corrcalc), intent(inout) :: self              !< instance
        integer,                 intent(in)    :: iref, iptcl, irot !< reference, particle, rotation
        real,                    intent(in)    :: shvec(2)          !< origin shift vector
        real    :: argmat(self%refsz,self%kfromto(1):self%kfromto(2)), sqsum_ref_sh, cc
        complex :: shmat(self%refsz,self%kfromto(1):self%kfromto(2))
        complex :: pft_ref_sh(self%refsz,self%kfromto(1):self%kfromto(2))
        if( allocated(self%ctfmats) )then
            ! generate the argument matrix from memoized components in argtransf
            argmat = self%argtransf(:self%refsz,:) * shvec(1) + self%argtransf(self%refsz+1:,:) * shvec(2)
            ! generate the complex shift transformation matrix
            shmat = cmplx(cos(argmat),sin(argmat))
            ! shift
            pft_ref_sh  = (self%pfts_refs(iref,:,:) * self%ctfmats(iptcl,:,:)) * shmat
            ! calculate correlation precursors
            argmat = real(pft_ref_sh * conjg(self%pfts_ptcls(iptcl,irot:irot+self%winsz,:)))
            cc = sum(argmat)
            sqsum_ref_sh = sum(csq(pft_ref_sh))
        else if( allocated(self%pfts_refs_ctf) )then
            ! generate the argument matrix from memoized components in argtransf
            argmat = self%argtransf(:self%refsz,:) * shvec(1)+self%argtransf(self%refsz+1:,:) * shvec(2)
            ! generate the complex shift transformation matrix
            shmat = cmplx(cos(argmat),sin(argmat))
            ! shift
            pft_ref_sh = self%pfts_refs_ctf(iref,:,:) * shmat
            ! calculate correlation precursors
            argmat = real(pft_ref_sh * conjg(self%pfts_ptcls(iptcl,irot:irot+self%winsz,:)))
            cc = sum(argmat)
            sqsum_ref_sh = sum(csq(pft_ref_sh))
        else
            ! generate the argument matrix from memoized components in argtransf
            argmat = self%argtransf(:self%refsz,:) * shvec(1)+self%argtransf(self%refsz+1:,:) * shvec(2)
            ! generate the complex shift transformation matrix
            shmat = cmplx(cos(argmat),sin(argmat))
            ! shift
            pft_ref_sh = self%pfts_refs(iref,:,:) * shmat
            ! calculate correlation precursors
            argmat = real(pft_ref_sh * conjg(self%pfts_ptcls(iptcl,irot:irot+self%winsz,:)))
            cc = sum(argmat)
            sqsum_ref_sh = sum(csq(pft_ref_sh))
        endif
        cc = cc/sqrt(sqsum_ref_sh*self%sqsums_ptcls(iptcl))
    end function corr_2
    
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
            self%existence = .false.
        endif
    end subroutine kill

    !>  \brief  is a unit test for polarft
    ! subroutine test_polarft_corrcalc
    !     use simple_image,   only: image
    !     use simple_stat,    only: pearsn
    !     use simple_math,    only: euclid
    !     type(polarft_corrcalc) :: pftcc
    !     integer, parameter     :: NITER=3000
    !     real                   :: diff, lscape_dist, lscape_corr
    !     integer                :: i, j, nrots, r
    !     real, allocatable      :: corrs1(:), corrs2(:)
    !     complex, allocatable   :: pft(:,:)
    !     write(*,'(a)') '**info(simple_polarft_corrcalc_unit_test): testing all functionality'
    !     do i=1,3
    !         do j=1,5
    !             call pfts1(i,j)%new([2,11], 45, ptcl=.false.) ! implemented as reference
    !             call pfts2(i,j)%new([2,11], 45, ptcl=.true.)  ! implemented as particle
    !             call pfts1(i,j)%set_ldim([100,100,1])
    !             call pfts2(i,j)%set_ldim([100,100,1])
    !         end do
    !     end do
    !     nrots = pfts1(1,1)%get_nrots()
    !     allocate(corrs1(nrots), corrs2(nrots))
    !     ! make a conforming polarft_corrcalc object of 3 refs and 5 ptcls
    !     call pftcc%new(3, [1,5], [100,100,1], [2,11], 45, 'no')
    !     do i=1,3
    !         do j=1,5
    !             call pfts1(i,j)%noise_fill
    !             call pfts2(i,j)%noise_fill
    !             call pfts1(i,j)%memoize_sqsums
    !             call pfts2(i,j)%memoize_sqsums
    !             ! pfts1 is ref
    !             pft = pfts1(i,j)%get_pft()
    !             call pftcc%set_ref_pft(i, pft)
    !             deallocate(pft)
    !             ! pfts2 is ptcl
    !             pft = pfts2(i,j)%get_pft()
    !             call pftcc%set_ptcl_pft(j, pft)
    !             deallocate(pft)
    !         end do
    !     end do
    !     write(*,'(a)') '**info(simple_polarft_corrcalc_unit_test, part 1):&
    !     &comparing the gencorrs correlations with the ones generated here'
    !     do i=1,3
    !         do j=1,5
    !             corrs1 = pftcc%gencorrs(i,j)
    !             do r=1,nrots
    !                 corrs2(r) = pftcc%corr(i, j, r)
    !             end do
    !             lscape_corr = pearsn(corrs1,corrs2)
    !             lscape_dist = euclid(corrs1,corrs2)
    !             if( lscape_dist > 1e-6 .or. lscape_corr < 0.9999 )then
    !                 write(*,*) 'landscape correlation: ', lscape_corr
    !                 write(*,*) 'landscape dist: ', lscape_dist
    !                 stop 'corr and gencorrs not consistent in pftcc'
    !             endif
    !         end do
    !     end do
    !     write(*,'(a)') '**info(simple_polarft_corrcalc_unit_test, part 2):&
    !     &testing the gencorrs function'
    !     do i=1,3
    !         do j=1,5
    !             call pfts1(i,j)%noise_fill
    !             call pfts2(i,j)%noise_fill
    !             call pfts1(i,j)%memoize_sqsums
    !             call pfts2(i,j)%memoize_sqsums
    !             ! pfts1 is ref
    !             pft = pfts1(i,j)%get_pft()
    !             call pftcc%set_ref_pft(i, pft)
    !             deallocate(pft)
    !             ! pfts2(i,j) is ptcl
    !             pft = pfts2(i,j)%get_pft()
    !             call pftcc%set_ptcl_pft(j, pft)
    !             deallocate(pft)
    !             ! compare the two correlation functions
    !             diff = 0.
    !             call pfts1(i,j)%gencorrs(pfts2(i,j), corrs1)
    !             corrs2 = pftcc%gencorrs(i,j)   
    !             diff = sum(abs(corrs1-corrs2))/real(nrots*3*5)
    !             if( diff > 1.5e-2 )then
    !                 print *, 'diff = ', diff
    !                 stop 'the new gencorrs function does not conform with the old; test_polarft_corrcalc'
    !             endif
    !         end do
    !     end do
    !     do i=1,3
    !         do j=1,5
    !             call pfts1(i,j)%noise_fill
    !             call pfts2(i,j)%noise_fill
    !             call pfts1(i,j)%memoize_sqsums
    !             call pfts2(i,j)%memoize_sqsums
    !             ! pfts1 is ref
    !             pft = pfts1(i,j)%get_pft()
    !             call pftcc%set_ref_pft(i, pft)
    !             deallocate(pft)
    !             ! pfts2(i,j) is ptcl
    !             pft = pfts2(i,j)%get_pft()
    !             call pftcc%set_ptcl_pft(j, pft)
    !             deallocate(pft)
    !         end do
    !     end do
    !     write(*,'(a)') 'SIMPLE_POLARFT_CORRCALC_UNIT_TEST COMPLETED SUCCESSFULLY ;-)'
    ! end subroutine test_polarft_corrcalc

    !<<< THIS ROUTINE SHOULD BE SIMD INSTRUCTED >>>
    !>  \brief  for calculating the correlation between reference iref and particle iptcl in rotation irot
    ! function corr_1( self, iref, iptcl, irot ) result( cc )
    !     use simple_math, only: csq
    !     class(polarft_corrcalc), intent(inout) :: self              !< instance
    !     integer,                 intent(in)    :: iref, iptcl, irot !< reference, particle, rotation
    !     real    :: cc, ptcl_sqsum_tmp, ref_sqsum_tmp
    !     integer :: kcount, k
    !     if( self%sqsums_refs(iref) < TINY .or. self%sqsums_ptcls(iptcl) < TINY )then
    !         cc = 0.
    !         return
    !     endif
        ! if( allocated(self%pfts_refs_ctf) )then
            ! cc = sum(real(self%pfts_refs_ctf(iref,:,:) * conjg(self%pfts_ptcls(iptcl,irot:irot+self%refsz-1,:))))
        ! else if( self%xfel )then
        !     cc = 0.
        !     kcount = 0
        !     ! Shell-by-shell xfel correlation
        !     do k=self%kfromto(1),self%kfromto(2)
        !        ptcl_sqsum_tmp = sum(real(self%pfts_ptcls(iptcl,irot:irot+self%refsz-1,k)) &
        !             *real(self%pfts_ptcls(iptcl,irot:irot+self%refsz-1,k)))
        !        ref_sqsum_tmp = sum(real(self%pfts_refs(iref,:,k))*real(self%pfts_refs(iref,:,k)))
        !        cc = cc + sum(real(self%pfts_refs(iref,:,k)) &
        !             *real(self%pfts_ptcls(iptcl,irot:irot+self%refsz-1,k))) &
        !             / sqrt(ref_sqsum_tmp * ptcl_sqsum_tmp)
        !        kcount = kcount + 1
        !     end do
        !     cc = cc / kcount
        ! else
           ! cc = sum(real(self%pfts_refs(iref,:,:) * conjg(self%pfts_ptcls(iptcl,irot:irot+self%refsz-1,:))))
        ! endif
        ! if( self%xfel) then
        ! else 
            ! cc = cc/sqrt(self%sqsums_refs(iref)*self%sqsums_ptcls(iptcl))
        ! endif
    ! end function corr_1

    ! RETIRED OLD GPU-STYLE ROUTINES WITH CSHIFT BS

    !>  \brief  is for generating all rotational correlations for all nptcls particles and all nrefs
    !!          references (nptcls=nrefs) (assumes that the first and second dimensions of the particle
    !!          matrix pfts_ptcls have been expanded and that self%nptcls=self%nrefs)
    ! subroutine gencorrs_all_cpu( self, corrmat3dout )
    !     class(polarft_corrcalc), intent(inout) :: self !< instance
    !     real,                    intent(out)   :: corrmat3dout(self%nptcls,self%nptcls,self%nrots)
    !     real    :: hadamard_prod(self%nptcls,self%refsz,self%nk)
    !     integer :: iref, iptcl, irot, d1lims(2), d2lims(2)
    !     if( .not. self%dim_expanded )&
    !     stop 'expand dimension with expand_dim before calling simple_polarft_corrcalc :: gencorrs_all_cpu'
    !     do iref=1,self%nptcls
    !         d1lims(1) = iref
    !         d1lims(2) = iref + self%nptcls - 1
    !         do irot=1,self%nrots
    !             d2lims(1) = irot 
    !             d2lims(2) = irot + self%refsz - 1
    !             ! calculate the Hadamard product
    !             !$omp parallel default(shared) private(iptcl)
    !             !$omp workshare
    !             hadamard_prod(:,:,:) = &
    !             real(self%pfts_refs(:,:,:)*conjg(self%pfts_ptcls(d1lims(1):d1lims(2),d2lims(1):d2lims(2),:)))
    !             !$omp end workshare
    !             ! sum over the rotational and resolution dimensions
    !             !$omp do schedule(auto)
    !             do iptcl=1,self%nptcls
    !                 corrmat3dout(iptcl,iref,irot) = sum(hadamard_prod(iptcl,:,:))
    !             end do
    !             !$omp end do
    !             ! normalise correlations
    !             !$omp workshare
    !             corrmat3dout(:,iref,irot) =&
    !             corrmat3dout(:,iref,irot)/&
    !             sqrt(self%sqsums_refs(:)*self%sqsums_ptcls(d1lims(1):d1lims(2)))
    !             !$omp end workshare nowait
    !             !$omp end parallel
    !         end do
    !     end do
    ! end subroutine gencorrs_all_cpu

    !>  \brief  prepare the real CTF matrix for GPU-based corr calc
    ! subroutine prep_ctf4gpu( self, smpd, a )
    !     !$ use omp_lib
    !     !$ use omp_lib_kinds
    !     use simple_ctf,   only: ctf
    !     use simple_oris,  only: oris
    !     use simple_image, only: image
    !     class(polarft_corrcalc), intent(inout) :: self
    !     real,                    intent(in)    :: smpd
    !     class(oris),             intent(inout) :: a
    !     type(image)       :: img
    !     type(ctf)         :: tfun
    !     real, allocatable :: ctfparams(:,:), ctfmat(:,:)
    !     integer           :: alloc_stat, iptcl
    !     ! create template image
    !     call img%new(self%ldim, smpd)
    !     ! get the CTF parameters
    !     ctfparams = a%get_ctfparams(self%pfromto)
    !     ! allocate the array that will carry the real CTF values
    !     if( allocated(self%pfts_ptcls_ctf) ) deallocate(self%pfts_ptcls_ctf)
    !     allocate( self%pfts_ptcls_ctf(2*self%nptcls,self%ptclsz,self%nk), stat=alloc_stat )
    !     call alloc_err("In: prep_ctf4gpu, simple_polarft_corrcalc, self%pfts_ptcls_ctf", alloc_stat)
    !     !$omp parallel default(shared) private(iptcl,tfun,ctfmat)
    !     !$omp do schedule(auto)
    !     do iptcl=1,self%nptcls
    !         ! we here need to create the CTF object online as kV/cs/fraca are now per-particle params
    !         tfun = ctf(smpd, ctfparams(iptcl,1), ctfparams(iptcl,2), ctfparams(iptcl,3))
    !         ! create the congruent polar matrix of real CTF values
    !         ctfmat = self%create_polar_ctfmat(tfun, ctfparams(iptcl,4),&
    !                  ctfparams(iptcl,5), ctfparams(iptcl,6), self%nrots)
    !         ! set it
    !         self%pfts_ptcls_ctf(iptcl,:self%nrots,:)   = ctfmat
    !         ! expand in in-plane rotation dimension
    !         self%pfts_ptcls_ctf(iptcl,self%nrots+1:,:) = ctfmat
    !     end do
    !     !$omp end do nowait
    !     !$omp workshare
    !     ! expand in projection direction dimension
    !     self%pfts_ptcls_ctf(self%nptcls+1:,:,:) = self%pfts_ptcls_ctf(:self%nptcls,:,:)
    !     !$omp end workshare nowait
    !     !$omp end parallel
    !     deallocate(ctfmat)
    !     call img%kill
    ! end subroutine prep_ctf4gpu

    !>  \brief  is for expanding the projection direction dimension to allow efficient
    !!          transfer to the GPU via subsectioning
    ! subroutine expand_dim( self )
    !     !$ use omp_lib
    !     !$ use omp_lib_kinds
    !     class(polarft_corrcalc), intent(inout) :: self
    !     complex, allocatable :: pfts_ptcls_tmp(:,:,:)
    !     real, allocatable    :: sqsums_ptcls_tmp(:)
    !     integer :: alloc_stat
    !     if( size(self%pfts_ptcls,dim=1) == 2*self%nptcls )then
    !         ! the arrays are already expanded, all we need to do is to copy data over
    !         !$omp parallel workshare
    !         self%pfts_ptcls(self%nptcls+1:,:,:) = self%pfts_ptcls(:self%nptcls,:,:)
    !         self%sqsums_ptcls(self%nptcls+1:)   = self%sqsums_ptcls(:self%nptcls)
    !         !$omp end parallel workshare
    !     else if( size(self%pfts_ptcls,dim=1) == self%nptcls )then
    !         allocate( pfts_ptcls_tmp(self%nptcls,self%ptclsz,self%nk), source=self%pfts_ptcls, stat=alloc_stat )
    !         call alloc_err("In: expand_dim, simple_polarft_corrcalc, pfts_ptcls_tmp", alloc_stat)
    !         allocate( sqsums_ptcls_tmp(self%nptcls), source=self%sqsums_ptcls, stat=alloc_stat )
    !         call alloc_err("In: expand_dim, simple_polarft_corrcalc, sqsums_ptcls_tmp", alloc_stat)
    !         deallocate(self%pfts_ptcls)
    !         deallocate(self%sqsums_ptcls)
    !         allocate( self%pfts_ptcls(2*self%nptcls,self%ptclsz,self%nk), stat=alloc_stat )
    !         call alloc_err("In: expand_dim, simple_polarft_corrcalc, self%pfts_ptcls", alloc_stat)
    !         allocate( self%sqsums_ptcls(2*self%nptcls), stat=alloc_stat )
    !         call alloc_err("In: expand_dim, simple_polarft_corrcalc, 4", alloc_stat)
    !         !$omp parallel workshare
    !         self%pfts_ptcls(:self%nptcls,:,:)   = pfts_ptcls_tmp
    !         self%pfts_ptcls(self%nptcls+1:,:,:) = pfts_ptcls_tmp
    !         self%sqsums_ptcls(:self%nptcls)     = sqsums_ptcls_tmp
    !         self%sqsums_ptcls(self%nptcls+1:)   = sqsums_ptcls_tmp
    !         !$omp end parallel workshare
    !     else
    !         stop 'nonconforming dimension (1) of self%pfts_ptcls; expand_dim; simple_polarft_corrcalc'
    !     endif
    !     self%dim_expanded = .true.
    ! end subroutine expand_dim
 
end module simple_polarft_corrcalc
