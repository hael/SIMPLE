! for calculation of band-pass limited cross-correlation of polar Fourier transforms
module simple_polarft_corrcalc
use simple_defs      ! use all in there
use simple_params,   only: params
use simple_ran_tabu, only: ran_tabu
use simple_syslib,   only: alloc_errchk, simple_stop
implicit none

public :: polarft_corrcalc
private

! CLASS PARAMETERS/VARIABLES
complex(sp), parameter :: zero=cmplx(0.,0.) !< just a complex zero

type :: polarft_corrcalc
    private
    integer                  :: pfromto(2) = 1         !< from/to particle indices (in parallel execution)
    integer                  :: nptcls     = 1         !< the total number of particles in partition (logically indexded [fromp,top])
    integer                  :: nrefs      = 1         !< the number of references (logically indexded [1,nrefs])
    integer                  :: nrots      = 0         !< number of in-plane rotations for one pft (determined by radius of molecule)
    integer                  :: ring2      = 0         !< radius of molecule
    integer                  :: refsz      = 0         !< size of reference (nrots/2) (number of vectors used for matching)
    integer                  :: ptclsz     = 0         !< size of particle (2*nrots)
    integer                  :: winsz      = 0         !< size of moving window in correlation calculations
    integer                  :: ldim(3)    = 0         !< logical dimensions of original cartesian image
    integer                  :: kfromto(2) = 0         !< Fourier index range
    real(sp)                 :: smpd       = 0.        !< sampling distance
    real(sp),    allocatable :: sqsums_ptcls(:)        !< memoized square sums for the correlation calculations
    real(sp),    allocatable :: angtab(:)              !< table of in-plane angles (in degrees)
    real(sp),    allocatable :: argtransf(:,:)         !< argument transfer constants for shifting the references
    real(sp),    allocatable :: polar(:,:)             !< table of polar coordinates (in Cartesian coordinates)
    real(sp),    allocatable :: ctfmats(:,:,:)         !< expandd set of CTF matrices (for efficient parallel exec)
    complex(sp), allocatable :: pfts_refs(:,:,:)       !< 3D complex matrix of polar reference sections (nrefs,refsz,nk)
    complex(sp), allocatable :: pfts_ptcls(:,:,:)      !< 3D complex matrix of particle sections
    logical                  :: with_ctf     = .false. !< CTF flag
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
    ! CALCULATORS
    procedure, private :: create_polar_ctfmat
    procedure          :: create_polar_ctfmats
    procedure, private :: gencorrs_1
    procedure, private :: gencorrs_2
    generic            :: gencorrs => gencorrs_1, gencorrs_2
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
    !! \param pfromto  from/to particle indices (in parallel execution)
    !! \param nrefs    the number of references (logically indexded [1,nrefs])
    !! \param ring2    radius of molecule
    !! \param ldim     logical dimensions of original cartesian image                              
    !! \param smpd     sampling distance                                                           
    !! \param kfromto  Fourier index range
    subroutine new( self, nrefs, pfromto, ldim, smpd, kfromto, ring2, ctfflag )
        use simple_math, only: rad2deg, is_even, round2even
        class(polarft_corrcalc), intent(inout) :: self    !< this object
        integer,                 intent(in)    :: nrefs, pfromto(2), ldim(3), kfromto(2), ring2
        real,                    intent(in)    :: smpd
        character(len=*),        intent(in)    :: ctfflag !< are we using ctf
        integer  :: alloc_stat, irot, k
        logical  :: even_dims, test(3)
        real(sp) :: ang
        ! kill possibly pre-existing object
        call self%kill
        ! error check
        if( kfromto(2) - kfromto(1) <= 2 )then
            write(*,*) 'kfromto: ', kfromto(1), kfromto(2)
            call simple_stop( 'resolution range too narrow; new; simple_polarft_corrcalc')
        endif
        if( ring2 < 1 )then
            write(*,*) 'ring2: ', ring2
            call simple_stop ( 'ring2 must be > 0; new; simple_polarft_corrcalc')
        endif
        if( pfromto(2) - pfromto(1) + 1 < 1 )then
            write(*,*) 'pfromto: ', pfromto(1), pfromto(2)
            call simple_stop ('nptcls (# of particles) must be > 0; new; simple_polarft_corrcalc')
        endif
        if( nrefs < 1 )then
            write(*,*) 'nrefs: ', nrefs
            call simple_stop ('nrefs (# of reference sections) must be > 0; new; simple_polarft_corrcalc')
        endif
        if( any(ldim == 0) )then
            write(*,*) 'ldim: ', ldim
            call simple_stop ('ldim is not conforming (is zero); new; simple_polarft_corrcalc')
        endif
        if( ldim(3) > 1 )then
            write(*,*) 'ldim: ', ldim
            call simple_stop ('3D polarfts are not yet supported; new; simple_polarft_corrcalc')
        endif
        test    = .false.
        test(1) = is_even(ldim(1))
        test(2) = is_even(ldim(2))
        test(3) = ldim(3) == 1
        even_dims = all(test)
        if( .not. even_dims )then
            write(*,*) 'ldim: ', ldim
            call simple_stop ('only even logical dims supported; new; simple_polarft_corrcalc')
        endif
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
        self%smpd    = smpd                            !< sampling distance
        self%kfromto = kfromto                         !< Fourier index range
        ! generate polar coordinates
        allocate( self%polar(self%ptclsz,self%kfromto(1):self%kfromto(2)), self%angtab(self%nrots), stat=alloc_stat)
        call alloc_errchk('polar coordinate arrays; new; simple_polarft_corrcalc', alloc_stat)
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
        call alloc_errchk('shift argument transfer array; new; simple_polarft_corrcalc', alloc_stat)
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
        call alloc_errchk('polarfts and sqsums; new; simple_polarft_corrcalc', alloc_stat)
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
        class(polarft_corrcalc), intent(inout) :: self     !< this object
        integer,                 intent(in)    :: iref     !< reference index
        complex(sp),             intent(in)    :: pft(:,:) !< reference pft
        self%pfts_refs(iref,:,:) = pft(:self%refsz,:)
    end subroutine set_ref_pft

    !>  \brief  sets particle pft iptcl
    subroutine set_ptcl_pft( self, iptcl, pft )
        class(polarft_corrcalc), intent(inout) :: self  !< this object
        integer,                 intent(in)    :: iptcl  !< particle index
        complex(sp),             intent(in)    :: pft(:,:) !< particle's pft
        self%pfts_ptcls(iptcl,:self%nrots,:)   = pft
        self%pfts_ptcls(iptcl,self%nrots+1:,:) = pft ! because rot dim is expanded
        ! calculate the square sum required for correlation calculation
        call self%memoize_sqsum_ptcl(iptcl)
    end subroutine set_ptcl_pft

    !>  \brief set_ref_fcomp sets a reference Fourier component
    !! \param iref reference index 
    !! \param irot rotation index
    !! \param k  index (third dim ptfs_refs)
    !! \param comp Fourier component
    !!
    subroutine set_ref_fcomp( self, iref, irot, k, comp )
        class(polarft_corrcalc), intent(inout) :: self  !< this object
        integer,                 intent(in)    :: iref, irot, k
        complex(sp),             intent(in)    :: comp
        self%pfts_refs(iref,irot,k) = comp
    end subroutine set_ref_fcomp
    
    !>  \brief  sets a particle Fourier component
    !! \param iptcl particle index 
    !! \param irot rotation index
    !! \param k  index (third dim ptfs_ptcls)
    !! \param comp Fourier component
    !!
    subroutine set_ptcl_fcomp( self, iptcl, irot, k, comp )
        class(polarft_corrcalc), intent(inout) :: self  !< this object
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
            call simple_stop ('pfts_refs and pfts_ptcls not congruent (nrefs .ne. nptcls)')
        endif
    end subroutine cp_ptcls2refs

    !>  \brief  copies the particles to the references
    !! \param iref reference index
    !! \param iptcl particle index 
    !! \param irot rotation index
    !!
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
     !! \param iref reference index
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
        integer,                 intent(in) :: roind !<  in-plane rotation index
        real(sp) :: rot
        if( roind < 1 .or. roind > self%nrots )then
            call simple_stop( 'roind is out of range; get_rot; simple_polarft_corrcalc')
        endif
        rot = self%angtab(roind)
    end function get_rot

    !>  \brief is for getting the discrete in-plane rotational
    !!         index corresponding to continuous rotation rot
    function get_roind( self, rot ) result( ind )
        class(polarft_corrcalc), intent(in) :: self
        real(sp),                intent(in) :: rot !<  continuous rotation 
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
    !! \param ang angle
    !! \param winsz window size
    !! \return roind_vec in-plane rotational indices
    !!
    function get_win_roind( self, ang, winsz )result( roind_vec )
        use simple_math, only: rad2deg
        class(polarft_corrcalc), intent(in) :: self
        real(sp),                intent(in) :: ang, winsz
        integer, allocatable :: roind_vec(:)  
        real(sp) :: dist(self%nrots)
        integer  :: i, irot, nrots, alloc_stat
        if(ang>360. .or. ang<TINY)call simple_stop ('input angle outside of the conventional range; simple_polarft_corrcalc::get_win_roind')
        if(winsz<0. .or. winsz>180.)call simple_stop ('invalid window size; simple_polarft_corrcalc::get_win_roind')
        if(winsz < 360./real(self%nrots))call simple_stop ('too small window size; simple_polarft_corrcalc::get_win_roind')
        i    = self%get_roind( ang )
        dist = abs(self%angtab(i) - self%angtab)
        where( dist>180. )dist = 360.-dist
        nrots = count(dist <= winsz)
        allocate( roind_vec(nrots), stat=alloc_stat )
        call alloc_errchk("In: get_win_roind; simple_polarft_corrcalc", alloc_stat)
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
    !! \param irot rotation index
    !! \param iptcl particle index 
    function get_ptcl_pft( self, iptcl, irot ) result( pft )
        class(polarft_corrcalc), intent(in) :: self
        integer,                 intent(in) :: iptcl, irot
        complex(sp), allocatable :: pft(:,:)
        integer :: alloc_stat
        allocate(pft(self%refsz,self%kfromto(1):self%kfromto(2)),&
        source=self%pfts_ptcls(iptcl,irot:irot+self%winsz,:), stat=alloc_stat)
        call alloc_errchk("In: get_ptcl_pft; simple_polarft_corrcalc", alloc_stat)
    end function get_ptcl_pft
    
    !>  \brief  returns polar Fourier transform of reference iref
    !! \param iref reference index
    function get_ref_pft( self, iref ) result( pft )
        class(polarft_corrcalc), intent(in) :: self
        integer,                 intent(in) :: iref
        complex(sp), allocatable :: pft(:,:)
        integer :: alloc_stat
        allocate(pft(self%refsz,self%kfromto(1):self%kfromto(2)),&
        source=self%pfts_refs(iref,:,:), stat=alloc_stat)
        call alloc_errchk("In: get_ref_pft; simple_polarft_corrcalc", alloc_stat)
    end function get_ref_pft

    !>  \brief  checks for existence
    function exists( self ) result( yes )
        class(polarft_corrcalc), intent(in) :: self
        logical :: yes
        yes = self%existence
    end function exists

    ! PRINTERS/VISUALISERS

    !>  \brief  is for plotting a particle polar FT
    !! \param iptcl particle index 
    subroutine vis_ptcl( self, iptcl )
        use gnufor2
        class(polarft_corrcalc), intent(in) :: self
        integer,                 intent(in) :: iptcl
        call gnufor_image(real(self%pfts_ptcls(iptcl,:self%refsz,:)),  palette='gray')
        call gnufor_image(aimag(self%pfts_ptcls(iptcl,:self%refsz,:)), palette='gray')
    end subroutine vis_ptcl
    
    !>  \brief  is for plotting a particle polar FT
    !! \param iref reference index
   
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
    !! \param iptcl particle index 
    subroutine memoize_sqsum_ptcl( self, iptcl )
        use simple_math, only: csq
        class(polarft_corrcalc), intent(inout) :: self
        integer,                 intent(in)    :: iptcl
        self%sqsums_ptcls(iptcl) = sum(csq(self%pfts_ptcls(iptcl,:self%refsz,:)))
    end subroutine memoize_sqsum_ptcl

    ! I/O

    !>  \brief  is for writing particle pfts to file
    subroutine write_pfts_ptcls( self, fname )
        use simple_fileio      
        class(polarft_corrcalc), intent(in) :: self
        character(len=*),        intent(in) :: fname !< output filename
        integer :: funit, io_stat
        if(.not.fopen(funit, status='REPLACE', action='WRITE', file=trim(fname), access='STREAM', iostat=io_stat))&
             call fileio_errmsg('polarft_corrcalc write_ptfs_ptcls  ', io_stat)
        write(unit=funit,pos=1,iostat=io_stat) self%pfts_ptcls
        ! Check if the write was successful
        if( io_stat .ne. 0 )then
            call fileio_errmsg('**ERROR(simple_polarft_corrcalc): I/O error when writing file: '// trim(fname), io_stat)
            call simple_stop ('I/O error; write_pfts_ptcls')
        endif
         if(.not.fclose(funit, iostat=io_stat))&
             call fileio_errmsg('polarft_corrcalc write_ptfs_ptcls  ', io_stat)
    end subroutine write_pfts_ptcls

    !>  \brief  is for reading particle pfts from file
    subroutine read_pfts_ptcls( self, fname )
        !$ use omp_lib
        !$ use omp_lib_kinds
        use simple_fileio 
        class(polarft_corrcalc), intent(inout) :: self
        character(len=*),        intent(in)    :: fname !< input filename
        integer :: funit, io_stat, iptcl
        
         if(.not.fopen(funit, status='OLD', action='READ', file=trim(fname), access='STREAM', iostat=io_stat))&
             call fileio_errmsg('polarft_corrcalc read_ptfs_ptcls fopen ', io_stat)
        read(unit=funit,pos=1,iostat=io_stat) self%pfts_ptcls
        ! Check if the read was successful
        if( io_stat .ne. 0 )then
            call fileio_errmsg('**ERROR(simple_polarft_corrcalc): I/O error when reading file: '// trim(fname), io_stat)
        endif
         if(.not.fclose(funit, iostat=io_stat))&
             call fileio_errmsg('polarft_corrcalc read_ptfs_ptcls fclose ', io_stat)
        ! memoize sqsum_ptcls
        !$omp parallel do schedule(static) default(shared) private(iptcl) proc_bind(close)
        do iptcl=self%pfromto(1),self%pfromto(2)
            call self%memoize_sqsum_ptcl(iptcl)
        end do
        !$omp end parallel do
    end subroutine read_pfts_ptcls
    
    ! MODIFIERS
    !> prep_ref4corr
    !! \param iref reference index
    !! \param iptcl particle index 
    !! \param pft_ref references
    !! \param sqsum_ref squared sum references
    !! \param kstop end point
    !!
    subroutine prep_ref4corr( self, iptcl, iref, pft_ref, sqsum_ref, kstop )
        use simple_math, only: csq
        class(polarft_corrcalc), intent(inout) :: self
        integer,                 intent(in)    :: iptcl, iref
        complex(sp),             intent(out)   :: pft_ref(self%refsz,self%kfromto(1):self%kfromto(2))
        real(sp),                intent(out)   :: sqsum_ref
        integer, optional,       intent(in)    :: kstop
        if( self%with_ctf )then
            pft_ref = self%pfts_refs(iref,:,:) * self%ctfmats(iptcl,:,:)
        else
            pft_ref = self%pfts_refs(iref,:,:)
        endif
        if( present(kstop) )then
            sqsum_ref = sum(csq(pft_ref(:,self%kfromto(1):kstop)))
        else
            sqsum_ref = sum(csq(pft_ref))
        endif
    end subroutine prep_ref4corr

    ! CALCULATORS

    !>  \brief create_polar_ctfmat  is for generating a matrix of CTF values
    !! \param tfun transfer function object
    !! \param dfx,dfy resolution along Fourier axes 
    !! \param angast astigmatic angle
    !! \param endrot number of rotations
    !! \return ctfmat matrix with CTF values
    !!
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
    subroutine create_polar_ctfmats( self, a )
        use simple_ctf,  only: ctf
        use simple_oris, only: oris
        class(polarft_corrcalc), intent(inout) :: self
        class(oris),             intent(inout) :: a !< oris object
        type(ctf) :: tfun
        integer   :: iptcl,alloc_stat 
        real(sp)  :: kv,cs,fraca,dfx,dfy,angast
        logical   :: astig
        astig = a%isthere('dfy')
        if( allocated(self%ctfmats) ) deallocate(self%ctfmats)
        allocate(self%ctfmats(self%pfromto(1):self%pfromto(2),self%refsz,self%kfromto(1):self%kfromto(2)), stat=alloc_stat)
        call alloc_errchk("In: simple_polarft_corrcalc :: create_polar_ctfmats, 2", alloc_stat)
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
            tfun = ctf(self%smpd, kv, cs, fraca)
            self%ctfmats(iptcl,:,:) = self%create_polar_ctfmat(tfun, dfx, dfy, angast, self%refsz)
        end do
    end subroutine create_polar_ctfmats

    !>  \brief gencorrs is for generating rotational correlations
    !! \param iref reference index
    !! \param iptcl particle index
    !! \param roind_vec vector of rotational indices
    function gencorrs_1( self, iref, iptcl, roind_vec ) result( cc )
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
    end function gencorrs_1

    !>  \brief  is for generating rotational correlations
    !! \param iref reference index
    !! \param iptcl particle index
    !! \param roind_vec vector of rotational indices
    function gencorrs_2( self, iref, iptcl, kstop, roind_vec ) result( cc )
    use simple_math, only: csq
        class(polarft_corrcalc), intent(inout) :: self        !< instance
        integer,                 intent(in)    :: iref, iptcl
        integer,                 intent(in)    :: kstop    !< last frequency
        integer,       optional, intent(in)    :: roind_vec(:)
        complex(sp) :: pft_ref(self%refsz,self%kfromto(1):self%kfromto(2))
        real(sp)    :: cc(self%nrots), sqsum_ref, sqsum_ptcl
        integer     :: irot, i
        call self%prep_ref4corr(iptcl, iref, pft_ref, sqsum_ref, kstop)
        sqsum_ptcl = sum(csq(self%pfts_ptcls(iptcl,:self%refsz,self%kfromto(1):kstop)))
        if( present(roind_vec) )then
            ! calculates only corrs for rotational indices provided in roind_vec
            ! see get_win_roind. returns -1. when not calculated
            if( any(roind_vec<=0) .or. any(roind_vec>self%nrots) )&
                &call simple_stop ('index out of range; simple_polarft_corrcalc::gencorrs')
            cc = -1.
            do i = 1, size(roind_vec)
                irot = roind_vec(i)
                cc(irot) = sum(real( pft_ref(:,self%kfromto(1):kstop) * &
                    conjg(self%pfts_ptcls(iptcl,irot:irot+self%winsz,self%kfromto(1):kstop)) ))      
                cc(irot) = cc(irot) / sqrt(sqsum_ref * sqsum_ptcl)
            end do
        else
            ! all rotations
            ! numerator
            do irot = 1, self%nrots
                cc(irot) = sum(real( pft_ref(:,self%kfromto(1):kstop) * &
                    conjg(self%pfts_ptcls(iptcl,irot:irot+self%winsz,self%kfromto(1):kstop)) ))
            end do
            ! denominator
            cc = cc / sqrt(sqsum_ref * sqsum_ptcl)
        endif
    end function gencorrs_2

    !>  \brief  is for generating resolution dependent correlations
    !! \param iref reference index
    !! \param iptcl particle index
    function genfrc( self, iref, iptcl, irot ) result( frc )
        use simple_math, only: csq
        class(polarft_corrcalc), target, intent(inout) :: self              !< instance
        integer,                         intent(in)    :: iref, iptcl, irot !< rotation index
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
    !! \param iref reference index
    !! \param iptcl particle index
    !! \param irot rotational 
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
        ! finalize cross-correlation
        cc = cc / sqrt(sqsum_ref * self%sqsums_ptcls(iptcl))
    end function corr_1

    !>  \brief  for calculating the on-fly shifted correlation between reference iref and particle iptcl in rotation irot
    !! \param iref reference index
    !! \param iptcl particle index
    !! \param irot rotational index
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

    !>  \brief  for calculating the normalized Euclidean distance between reference & particle
    !! \param iref reference index
    !! \param iptcl particle index
    !! \param irot rotational indices
    function euclid( self, iref, iptcl, irot, bfac ) result( dist )
        class(polarft_corrcalc), intent(inout) :: self
        integer,                 intent(in)    :: iref, iptcl, irot
        real, optional,          intent(in)    :: bfac  !< B-factor
        real        :: dist, sqsum_ref, rn
        integer     :: k
        complex(sp) :: pft_ref(self%refsz,self%kfromto(1):self%kfromto(2))
        complex(sp) :: pft_ptcl(self%refsz,self%kfromto(1):self%kfromto(2))
        ! prep ptcl PFT
        if( present(bfac) )then
             ! apply B-factor
            do k=self%kfromto(1),self%kfromto(2)
                pft_ptcl(:,k) = self%pfts_ptcls(iptcl,irot:irot+self%winsz,k) * eval_bfac(k)
            end do
        else
            ! just copy
            pft_ptcl = self%pfts_ptcls(iptcl,irot:irot+self%winsz,:)
        endif
        ! prep ref PFT
        call self%prep_ref4corr(iptcl, iref, pft_ref, sqsum_ref)
        dist = sum(cabs(pft_ref - pft_ptcl)**2.0)/real(self%refsz * (self%kfromto(2) - self%kfromto(1) + 1))

        contains

            real function eval_bfac( k )
                integer, intent(in) :: k
                real :: res
                res = real(k)/(real(self%ldim(1))*self%smpd) ! assuming square dimensions
                eval_bfac = max(0.,exp(-(bfac/4.)*res*res))
            end function eval_bfac

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
