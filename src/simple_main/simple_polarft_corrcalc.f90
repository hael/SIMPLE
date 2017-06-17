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
complex(sp), parameter :: zero=cmplx(0.,0.) !< just a complex zero
logical,     parameter :: DEBUG = .true.

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
    real(sp),    allocatable :: sqsums_ptcls(:)        !< memoized square sums for the correlation calculations
    real(sp),    allocatable :: angtab(:)              !< table of in-plane angles (in degrees)
    real(sp),    allocatable :: argtransf(:,:)         !< argument transfer constants for shifting the references
    real(sp),    allocatable :: polar(:,:)             !< table of polar coordinates (in Cartesian coordinates)
    real(sp),    allocatable :: ctfmats(:,:,:)         !< expandd set of CTF matrices (for efficient parallel exec)
    complex(sp), allocatable :: pfts_refs(:,:,:)       !< 3D complex matrix of polar reference sections (nrefs,refsz,nk)
    complex(sp), allocatable :: pfts_ptcls(:,:,:)      !< 3D complex matrix of particle sections
    logical                  :: with_ctf     = .false. !< CTF flag
    logical                  :: xfel         = .false. !< to indicate whether we process xfel patterns or not
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
    procedure          :: cp_ptcl2ref
    procedure          :: zero_ref
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
    ! MEMOIZER
    procedure, private :: memoize_sqsum_ptcl
    ! I/O
    procedure          :: write_pfts_ptcls
    procedure          :: read_pfts_ptcls
    ! MODIFIERS
    procedure, private :: prep_ref4corr
    procedure          :: xfel_subtract_shell_mean
    ! CALCULATORS
    procedure, private :: create_polar_ctfmat
    procedure          :: create_polar_ctfmats
    procedure          :: gencorrs
    procedure          :: genfrc
    procedure, private :: corr_1
    procedure, private :: corr_2
    generic            :: corr => corr_1, corr_2
    procedure          :: euclid
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
        integer  :: alloc_stat, irot, k
        logical  :: even_dims, test(3)
        real(sp) :: ang
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
                    self%sqsums_ptcls(self%pfromto(1):self%pfromto(2)), stat=alloc_stat)
        call alloc_err('polarfts and sqsums; new; simple_polarft_corrcalc', alloc_stat)
        self%pfts_refs    = zero
        self%pfts_ptcls   = zero
        self%sqsums_ptcls = 0.
        self%with_ctf     = .false.
        if( ctfflag .ne. 'no' ) self%with_ctf = .true.
        self%existence    = .true.
    end subroutine new

    ! SETTERS

    !>  \brief  sets reference pft iref
    subroutine set_ref_pft( self, iref, pft )
        class(polarft_corrcalc), intent(inout) :: self
        integer,                 intent(in)    :: iref
        complex(sp),             intent(in)    :: pft(:,:)
        self%pfts_refs(iref,:,:) = pft(:self%refsz,:)
    end subroutine set_ref_pft

    !>  \brief  sets particle pft iptcl
    subroutine set_ptcl_pft( self, iptcl, pft )
        class(polarft_corrcalc), intent(inout) :: self
        integer,                 intent(in)    :: iptcl
        complex(sp),             intent(in)    :: pft(:,:)
        self%pfts_ptcls(iptcl,:self%nrots,:)   = pft
        self%pfts_ptcls(iptcl,self%nrots+1:,:) = pft ! because rot dim is expanded
        ! calculate the square sum required for correlation calculation
        call self%memoize_sqsum_ptcl(iptcl)
    end subroutine set_ptcl_pft

    !>  \brief  sets a reference Fourier component
    subroutine set_ref_fcomp( self, iref, irot, k, comp )
        class(polarft_corrcalc), intent(inout) :: self
        integer,                 intent(in)    :: iref, irot, k
        complex(sp),             intent(in)    :: comp
        self%pfts_refs(iref,irot,k) = comp
    end subroutine set_ref_fcomp
    
    !>  \brief  sets a particle Fourier component
    subroutine set_ptcl_fcomp( self, iptcl, irot, k, comp )
        class(polarft_corrcalc), intent(inout) :: self
        integer,                 intent(in)    :: iptcl, irot, k
        complex(sp),             intent(in)    :: comp
        self%pfts_ptcls(iptcl,irot,k) = comp
        self%pfts_ptcls(iptcl,irot+self%nrots,k) = comp ! because rot dim is expanded
    end subroutine set_ptcl_fcomp

    !>  \brief  copies the particles to the references
    subroutine cp_ptcls2refs( self )
        class(polarft_corrcalc), intent(inout) :: self
        if( self%nrefs .eq. self%nptcls )then
            self%pfts_refs(:,:,:) = self%pfts_ptcls(:,:self%refsz,:)         
        else
            stop 'pfts_refs and pfts_ptcls not congruent (nrefs .ne. nptcls)'
        endif
    end subroutine cp_ptcls2refs

    !>  \brief  copies the particles to the references
    subroutine cp_ptcl2ref( self, iptcl, iref, irot )
        class(polarft_corrcalc), intent(inout) :: self
        integer,                 intent(in)    :: iptcl, iref
        integer, optional,       intent(in)    :: irot
        if( present(irot) )then
            self%pfts_refs(iref,:,:) = self%pfts_ptcls(iptcl,irot:irot+self%winsz,:)
        else
            self%pfts_refs(iref,:,:) = self%pfts_ptcls(iptcl,:self%refsz,:)
        endif
    end subroutine cp_ptcl2ref

    !>  \brief  zeroes the iref reference
    subroutine zero_ref( self, iref )
        class(polarft_corrcalc), intent(inout) :: self
        integer,                 intent(in)    :: iref
        self%pfts_refs(iref,:,:) = zero
    end subroutine zero_ref
    
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
        real(sp) :: rot
        if( roind < 1 .or. roind > self%nrots )then
            stop 'roind is out of range; get_rot; simple_polarft_corrcalc'
        endif
        rot = self%angtab(roind)
    end function get_rot

    !>  \brief is for getting the discrete in-plane rotational
    !!         index corresponding to continuous rotation rot
    function get_roind( self, rot ) result( ind )
        class(polarft_corrcalc), intent(in) :: self
        real(sp),                intent(in) :: rot
        real(sp) :: dists(self%nrots)
        integer  :: ind, loc(1)
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
        real(sp),                intent(in) :: ang, winsz
        integer, allocatable :: roind_vec(:)
        real(sp) :: dist(self%nrots)
        integer  :: i, irot, nrots, alloc_stat
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
        real(sp) :: xy(2)
        xy(1) = self%polar(rot,k)
        xy(2) = self%polar(self%nrots+rot,k)
    end function get_coord
    
    !>  \brief  returns polar Fourier transform of particle iptcl in rotation irot
    function get_ptcl_pft( self, iptcl, irot ) result( pft )
        class(polarft_corrcalc), intent(in) :: self
        integer,                 intent(in) :: iptcl, irot
        complex(sp), allocatable :: pft(:,:)
        integer :: alloc_stat
        allocate(pft(self%refsz,self%kfromto(1):self%kfromto(2)),&
        source=self%pfts_ptcls(iptcl,irot:irot+self%winsz,:), stat=alloc_stat)
        call alloc_err("In: get_ptcl_pft; simple_polarft_corrcalc", alloc_stat)
    end function get_ptcl_pft
    
    !>  \brief  returns polar Fourier transform of reference iref
    function get_ref_pft( self, iref ) result( pft )
        class(polarft_corrcalc), intent(in) :: self
        integer,                 intent(in) :: iref
        complex(sp), allocatable :: pft(:,:)
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

    !>  \brief  is for memoization of the complex square sums required for correlation calculation
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
        !$ use omp_lib
        !$ use omp_lib_kinds
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
        !$omp parallel do schedule(static) default(shared) private(iptcl) proc_bind(close)
        do iptcl=self%pfromto(1),self%pfromto(2)
            call self%memoize_sqsum_ptcl(iptcl)
        end do
        !$omp end parallel do
    end subroutine read_pfts_ptcls
    
    ! MODIFIERS

    subroutine prep_ref4corr( self, iptcl, iref, pft_ref, sqsum_ref )
        use simple_math, only: csq
        class(polarft_corrcalc), intent(inout) :: self
        integer,                 intent(in)    :: iptcl, iref
        complex(sp),             intent(out)   :: pft_ref(self%refsz,self%kfromto(1):self%kfromto(2))
        real(sp),                intent(out)   :: sqsum_ref
        if( self%with_ctf )then
            pft_ref = self%pfts_refs(iref,:,:) * self%ctfmats(iptcl,:,:)
        else
            pft_ref = self%pfts_refs(iref,:,:)
        endif
        sqsum_ref = sum(csq(pft_ref))
    end subroutine prep_ref4corr

    !>  \brief  is for preparing for XFEL pattern corr calc
    subroutine xfel_subtract_shell_mean( self )
        class(polarft_corrcalc), intent(inout) :: self
        real(sp), allocatable :: ptcls_mean_tmp(:,:,:)
        real(sp), allocatable :: refs_mean_tmp(:,:)
        integer :: iptcl, iref, irot, k
        allocate( ptcls_mean_tmp(2*self%nptcls,self%ptclsz,self%kfromto(1):self%kfromto(2)),&
        refs_mean_tmp(self%nrefs,self%kfromto(1):self%kfromto(2)))
        ! calculate the mean of each reference at each k shell
        do iref=1,self%nrefs
            do k=self%kfromto(1),self%kfromto(2)
                refs_mean_tmp(iref,k) = sum(real(self%pfts_refs(iref,:,k)))/real(self%refsz)
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
                ptcls_mean_tmp(iptcl,1,k) = sum(real(self%pfts_ptcls(iptcl,1:self%winsz,k)))/real(self%refsz)
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
        real(sp),                intent(in)    :: dfx, dfy, angast
        integer,                 intent(in)    :: endrot
        real(sp), allocatable :: ctfmat(:,:)
        real(sp)              :: inv_ldim(3),hinv,kinv,spaFreqSq,ang
        integer               :: irot,k
        allocate( ctfmat(endrot,self%kfromto(1):self%kfromto(2)) )
        inv_ldim = 1./real(self%ldim)
        !$omp parallel do collapse(2) default(shared) private(irot,k,hinv,kinv,spaFreqSq,ang)&
        !$omp schedule(static) proc_bind(close)
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
        real(sp),                intent(in)    :: smpd
        class(oris),             intent(inout) :: a
        type(ctf) :: tfun
        integer   :: iptcl,alloc_stat 
        real(sp)  :: kv,cs,fraca,dfx,dfy,angast
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

    !>  \brief  is for generating rotational correlations
    function gencorrs( self, iref, iptcl, roind_vec ) result( cc )
        use simple_math, only: csq
        class(polarft_corrcalc), intent(inout) :: self        !< instance
        integer,                 intent(in)    :: iref, iptcl !< ref & ptcl indices
        integer,       optional, intent(in)    :: roind_vec(:)
        complex(sp) :: pft_ref(self%refsz,self%kfromto(1):self%kfromto(2))
        real(sp)    :: cc(self%nrots), sqsum_ref
        integer     :: irot, i
        call self%prep_ref4corr(iptcl, iref, pft_ref, sqsum_ref)
        if( present(roind_vec) )then
            ! calculates only corrs for rotational indices provided in roind_vec
            ! see get_win_roind. returns -1. when not calculated
            if( any(roind_vec<=0) .or. any(roind_vec>self%nrots) )&
                &stop 'index out of range; simple_polarft_corrcalc::gencorrs'
            cc = -1.
            do i = 1, size(roind_vec)
                irot = roind_vec(i)
                cc(irot) = sum(real( pft_ref * conjg(self%pfts_ptcls(iptcl,irot:irot+self%winsz,:)) ))
                cc(irot) = cc(irot) / sqrt(sqsum_ref * self%sqsums_ptcls(iptcl))
            end do
        else
            ! all rotations
            ! numerator
            do irot = 1, self%nrots
                cc(irot) = sum(real( pft_ref * conjg(self%pfts_ptcls(iptcl,irot:irot+self%winsz,:)) ))
            end do
            ! denominator
            cc = cc / sqrt(sqsum_ref * self%sqsums_ptcls(iptcl))
        endif
    end function gencorrs

    !>  \brief  is for generating resolution dependent correlations
    function genfrc( self, iref, iptcl, irot ) result( frc )
        use simple_math, only: csq
        class(polarft_corrcalc), target, intent(inout) :: self              !< instance
        integer,                         intent(in)    :: iref, iptcl, irot !< reference, particle, rotation
        real(sp), allocatable :: frc(:)
        complex(sp) :: pft_ref_ctf(self%refsz,self%kfromto(1):self%kfromto(2))        
        real(sp)    :: sumsqref, sumsqptcl
        integer     :: k
        allocate( frc(self%kfromto(1):self%kfromto(2)) )
        if( self%with_ctf )then
            ! multiply reference with CTF
            pft_ref_ctf = self%pfts_refs(iref,:,:) * self%ctfmats(iptcl,:,:)
            ! calc FRC
            do k=self%kfromto(1),self%kfromto(2)
                frc(k)    = sum(real(pft_ref_ctf(:,k)*conjg(self%pfts_ptcls(iptcl,irot:irot+self%winsz,k))))
                sumsqref  = sum(csq(pft_ref_ctf(:,k)))
                sumsqptcl = sum(csq(self%pfts_ptcls(iptcl,:self%refsz,k)))
                if( sumsqref < TINY .or. sumsqptcl < TINY )then
                    frc(k) = 0.
                else
                    frc(k) = frc(k)/sqrt(sumsqref*sumsqptcl)
                endif
            end do
        else
            ! calc FRC
            do k=self%kfromto(1),self%kfromto(2)
                frc(k)    = sum(real(self%pfts_refs(iref,:,k)*conjg(self%pfts_ptcls(iptcl,irot:irot+self%winsz,k))))
                sumsqref  = sum(csq(self%pfts_refs(iref,:,k)))
                sumsqptcl = sum(csq(self%pfts_ptcls(iptcl,:self%refsz,k)))
                if( sumsqref < TINY .or. sumsqptcl < TINY )then
                    frc(k) = 0.
                else
                    frc(k) = frc(k)/sqrt(sumsqref*sumsqptcl)
                endif
            end do
        endif
    end function genfrc

    !>  \brief  for calculating the correlation between reference iref and particle iptcl in rotation irot
    function corr_1( self, iref, iptcl, irot ) result( cc )
        class(polarft_corrcalc), intent(inout) :: self              !< instance
        integer,                 intent(in)    :: iref, iptcl, irot !< reference, particle, rotation
        real(sp)    :: cc, sqsum_ref
        complex(sp) :: pft_ref(self%refsz,self%kfromto(1):self%kfromto(2))
        cc = 0.
        ! floating point check
        if( self%sqsums_ptcls(iptcl) < TINY ) return
        call self%prep_ref4corr(iptcl, iref, pft_ref, sqsum_ref)
        ! floating point check
        if( sqsum_ref < TINY ) return
        ! numerator
        cc = sum(real( pft_ref * conjg(self%pfts_ptcls(iptcl,irot:irot+self%winsz,:)) ))
        ! denominator
        cc = cc / sqrt(sqsum_ref * self%sqsums_ptcls(iptcl))
        ! check
        ! if( cc >= 1. ) print *,'cc out of range', iref, iptcl, cc
    end function corr_1

    !>  \brief  for calculating the on-fly shifted correlation between reference iref and particle iptcl in rotation irot
    function corr_2( self, iref, iptcl, irot, shvec ) result( cc )
        use simple_math, only: csq
        class(polarft_corrcalc), intent(inout) :: self              !< instance
        integer,                 intent(in)    :: iref, iptcl, irot !< reference, particle, rotation
        real(sp),                intent(in)    :: shvec(2)          !< origin shift vector
        real(sp)    :: sqsum_ref_sh, cc
        real(sp)    :: argmat(self%refsz,self%kfromto(1):self%kfromto(2))
        complex(sp) :: pft_ref_sh(self%refsz,self%kfromto(1):self%kfromto(2))
        complex(sp) :: shmat(self%refsz,self%kfromto(1):self%kfromto(2))
        cc = 0.
        ! floating point check
        if( self%sqsums_ptcls(iptcl) < TINY ) return
        ! generate the argument matrix from memoized components in argtransf
        argmat = self%argtransf(:self%refsz,:) * shvec(1) + self%argtransf(self%refsz+1:,:) * shvec(2)
        ! generate the complex shift transformation matrix
        shmat = cmplx(cos(argmat),sin(argmat))
        ! shift
        if( self%with_ctf)then
            pft_ref_sh = (self%pfts_refs(iref,:,:) * self%ctfmats(iptcl,:,:)) * shmat
        else
            pft_ref_sh = self%pfts_refs(iref,:,:) * shmat
        endif
        ! calculate correlation precursors
        sqsum_ref_sh = sum(csq(pft_ref_sh))
        ! floating point check
        if( sqsum_ref_sh < TINY  ) return
        argmat = real(pft_ref_sh * conjg(self%pfts_ptcls(iptcl,irot:irot+self%winsz,:)))
        cc     = sum(argmat)
        ! finalize cross-correlation
        cc = cc/sqrt(sqsum_ref_sh*self%sqsums_ptcls(iptcl))
    end function corr_2

    !>  \brief  for calculating the Euclidean distance between reference & particle
    !!          This ought to be a better choice for the orientation weights
    function euclid( self, iref, iptcl, irot ) result( dist )
        use gnufor2
        class(polarft_corrcalc), intent(inout) :: self
        integer,                 intent(in)    :: iref, iptcl, irot
        integer     :: k, nk
        real        :: pow_ref, pow_ptcl, norm, dist, sqsum_ref
        complex(sp) :: pft_ref_norm(self%refsz,self%kfromto(1):self%kfromto(2))
        complex(sp) :: pft_ptcl_norm(self%refsz,self%kfromto(1):self%kfromto(2))
        call self%prep_ref4corr(iptcl, iref, pft_ref_norm, sqsum_ref)
        do k=self%kfromto(1),self%kfromto(2)
            pow_ref  = sum(real(pft_ref_norm(:,k))*&
                &conjg(pft_ref_norm(:,k)))/real(self%refsz)
            pow_ptcl = sum(real(self%pfts_ptcls(iptcl,irot:irot+self%winsz,k))*&
                &conjg(self%pfts_ptcls(iptcl,irot:irot+self%winsz,k)))/real(self%refsz)
            if( pow_ref > TINY )then
                norm = sqrt(pow_ref)
                pft_ref_norm(:,k) = pft_ref_norm(:,k)/norm
            else
                pft_ref_norm(:,k) = cmplx(0.,0.)
            endif
            if( pow_ptcl > TINY )then
                norm = sqrt(pow_ptcl)
                pft_ptcl_norm(:,k) = self%pfts_ptcls(iptcl,irot:irot+self%winsz,k)/norm
            else
                pft_ptcl_norm(:,k) = cmplx(0.,0.)
            endif
        end do
        nk   = self%refsz * (self%kfromto(2) - self%kfromto(1) + 1)
        dist = sum(cabs(pft_ref_norm - pft_ptcl_norm)**2.0)/(nk * self%refsz)
    end function euclid

    ! DESTRUCTOR

    !>  \brief  is a destructor
    subroutine kill( self )
        class(polarft_corrcalc), intent(inout) :: self
        if( self%existence )then
            deallocate( self%sqsums_ptcls, &
                        self%angtab,       &
                        self%argtransf,    &
                        self%polar,        &
                        self%pfts_refs,    &
                        self%pfts_ptcls    )
            self%existence = .false.
        endif
    end subroutine kill
 
end module simple_polarft_corrcalc
