!>  \brief  SIMPLE polarft_corrcalc class
module simple_polarft_corrcalc
use, intrinsic :: iso_c_binding
use simple_defs ! singleton
use simple_cuda_defs
use simple_jiffys, only: alloc_err
use simple_ran_tabu, only: ran_tabu
use simple_timing
implicit none

public :: polarft_corrcalc, test_polarft_corrcalc
private

! CLASS PARAMETERS/VARIABLES
complex, parameter :: zero=cmplx(0.,0.) !< just a complex zero
logical, parameter :: debug=.false.     !< debug indicator

type :: polarft_corrcalc
    private
    type(polar_corr_calc):: s_polar              !< Data structure used in the kernel
    integer              :: pfromto(2) = 1       !< from/to particle indices (in parallel execution)
    integer              :: nptcls     = 1       !< the total number of particles in partition (logically indexded [fromp,top])
    integer              :: nrefs      = 1       !< the number of references (logically indexded [1,nrefs])
    integer              :: nrots      = 0       !< number of in-plane rotations for one pft (determined by radius of molecule)
    integer              :: ring2      = 0       !< radius of molecule
    integer              :: refsz      = 0       !< size of reference (nrots/2) (number of vectors used for matching)
    integer              :: ptclsz     = 0       !< size of particle (2*nrots)
    integer              :: ldim(3)    = 0       !< logical dimensions of original cartesian image
    integer              :: kfromto(2) = 0       !< Fourier index range
    integer              :: nk         = 0       !< number of resolution elements (Fourier components)
    real,    allocatable :: sqsums_refs(:)       !< memoized square sums for the correlation calculations
    real,    allocatable :: sqsums_ptcls(:)      !< memoized square sums for the correlation calculations
    real,    allocatable :: angtab(:)            !< table of in-plane angles (in degrees)
    real,    allocatable :: argtransf(:,:)       !< argument transfer constants for shifting the references
    real,    allocatable :: polar(:,:)           !< table of polar coordinates (in Cartesian coordinates)
    complex(sp), allocatable :: pfts_refs(:,:,:)     !< 3D complex matrix of polar reference sections (nrefs,refsz,nk)
    complex, allocatable :: pfts_refs_ctf(:,:,:) !< 3D complex matrix of polar reference sections with abs(CTF) applied
    complex, allocatable :: pfts_ptcls(:,:,:)    !< 3D complex matrix of particle sections (phase-flipped)
    logical              :: xfel=.false.         !< to indicate whether we process xfel patterns or not
    logical              :: existence=.false.    !< to indicate existence

  contains
    ! CONSTRUCTOR
    procedure :: new
    ! INDEX MAPPING ROUTINES
    procedure, private :: k_ind
    procedure, private :: ptcl_ind
    ! SETTERS/GETTERS
    procedure :: get_pfromto
    procedure :: get_nptcls
    procedure :: get_nrefs
    procedure :: get_nrots
    procedure :: get_ring2
    procedure :: get_refsz
    procedure :: get_ptclsz
    procedure :: get_ldim
    procedure :: get_kfromto
    procedure :: get_rot
    procedure :: get_roind
    procedure :: get_coord
    procedure :: get_ptcl_pft
    procedure :: get_ref_pft
    procedure :: set_ref_fcomp
    procedure :: set_ptcl_fcomp
    procedure :: set_ref_pft
    procedure :: set_ptcl_pft
    ! PRINTERS
    procedure :: print
    ! CHECKUPS
    procedure :: exists
    ! VISUALISATORS
    procedure :: vis_ptcl
    procedure :: vis_ref
    ! MEMOIZERS
    procedure :: memoize_sqsum_ref
    procedure :: memoize_sqsum_ref_ctf
    procedure, private :: memoize_sqsum_ptcl
    ! MODIFIERS
    procedure :: expand_dim
    ! CTF MODIFIERS
    procedure :: apply_ctf
    procedure :: copy_refs2ctf
    ! CALCULATORS
    procedure :: gencorrs_all
    procedure :: gencorrs_all_tester
    procedure :: gencorrs
    procedure :: rosrch_shc
    procedure, private :: corr_1
    procedure, private :: corr_2
    generic :: corr => corr_1, corr_2
    ! DESTRUCTOR
    procedure :: kill
end type

contains

    ! CONSTRUCTORS

    !>  \brief  is a constructor
    subroutine new( self, nrefs, pfromto, ldim, kfromto, ring2, ctfflag, isxfel )
        use simple_math, only: rad2deg, is_even, round2even
        class(polarft_corrcalc), intent(inout) :: self
        integer, intent(in)                    :: nrefs, pfromto(2), ldim(3), kfromto(2), ring2
        character(len=*), intent(in)           :: ctfflag
        character(len=*), intent(in), optional :: isxfel
        integer :: alloc_stat, irot, k, k_ind
        logical :: even_dims, test(3)
        real    :: ang
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

        !!!! IFDEF CUDA WE NEED TO CHECK NPTCLS == NREFS !!!!!!

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
        allocate(     self%pfts_refs(self%nrefs,self%refsz,self%nk),   &
                     self%pfts_ptcls(self%nptcls,self%ptclsz,self%nk), &
                    self%sqsums_refs(self%nrefs),                      &
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
        endif
        self%existence = .true.
<<<<<<< Updated upstream
    end subroutine new
   
=======
    end subroutine

    !>  \brief  for mapping physical Fourier index to implemented
    function k_ind( self, k_in ) result( k_out )
        class(polarft_corrcalc), intent(in) :: self
        integer, intent(in)                 :: k_in
        integer :: k_out
        k_out = k_in-self%kfromto(1)+1
    end function

    !>  \brief  for mapping physical particle index to implemented
    function ptcl_ind( self, iptcl_in ) result( iptcl_out )
        class(polarft_corrcalc), intent(in) :: self
        integer, intent(in)                 :: iptcl_in
        integer :: iptcl_out
        iptcl_out = iptcl_in-self%pfromto(1)+1
    end function

    !>  \brief  for getting the logical particle range
    function get_pfromto( self ) result( lim )
        class(polarft_corrcalc), intent(inout) :: self
        integer :: lim(2)
        lim = self%pfromto
    end function

    !>  \brief  for getting the number of particles
    pure function get_nptcls( self ) result( nptcls )
        class(polarft_corrcalc), intent(in) :: self
        integer :: nptcls
        nptcls = self%nptcls
    end function

    !>  \brief  for getting the number of references
    pure function get_nrefs( self ) result( nrefs )
        class(polarft_corrcalc), intent(in) :: self
        integer :: nrefs
        nrefs = self%nrefs
    end function

    !>  \brief  for getting the number of in-plane rotations
    pure function get_nrots( self ) result( nrots )
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

    !>  \brief  for getting the number of reference rotations (size of second dim of self%pfts_refs)
    function get_refsz( self ) result( refsz )
        class(polarft_corrcalc), intent(in) :: self
        integer :: refsz
        refsz = self%refsz
    end function

    !>  \brief  for getting the number of particle rotations (size of second dim of self%pfts_ptcls)
    function get_ptclsz( self ) result( ptclsz )
        class(polarft_corrcalc), intent(in) :: self
        integer :: ptclsz
        ptclsz = self%ptclsz
    end function

    !>  \brief  for getting the logical dimension of the original
    !!          Cartesian image
    function get_ldim( self ) result( ldim )
        class(polarft_corrcalc), intent(in) :: self
        integer :: ldim(3)
        ldim = self%ldim
    end function

    !>  \brief  for getting the Fourier index range (hp/lp)
    function get_kfromto( self ) result( lim )
        class(polarft_corrcalc), intent(inout) :: self
        integer :: lim(2)
        lim = self%kfromto
    end function

    !>  \brief is for getting the continuous in-plane rotation
    !!         corresponding to in-plane rotation index roind
    function get_rot( self, roind ) result( rot )
        class(polarft_corrcalc), intent(in) :: self
        integer, intent(in) :: roind
        real :: rot
        if( roind < 1 .or. roind > self%nrots )then
            stop 'roind is out of range; get_rot; simple_polarft_corrcalc'
        endif
        rot = self%angtab(roind)
    end function

    !>  \brief is for getting the discrete in-plane rotational
    !!         index corresponding to continuous rotation rot
    function get_roind( self, rot ) result( ind )
        class(polarft_corrcalc), intent(in) :: self
        real, intent(in) :: rot
        integer :: ind, irot, alloc_stat, loc(1)
        real    :: dists(self%nrots)
        do irot=1,self%nrots
            dists(irot) = sqrt((self%angtab(irot)-rot)**2.)
        end do
        loc = minloc(dists)
        ind = loc(1)
    end function

    !>  \brief returns polar coordinate for rotation rot
    !!         and Fourier index k
    function get_coord( self, rot, k ) result( xy )
        class(polarft_corrcalc), intent(in) :: self
        integer, intent(in)                 :: rot, k
        real    :: xy(2)
        integer :: k_ind
        k_ind = self%k_ind(k)
        xy(1) = self%polar(rot,k_ind)
        xy(2) = self%polar(self%nrots+rot,k_ind)
    end function

    !>  \brief  returns polar Fourier transform of particle iptcl in rotation irot
    function get_ptcl_pft( self, iptcl, irot ) result( pft )
        class(polarft_corrcalc), intent(in) :: self
        integer, intent(in)  :: iptcl, irot
        complex, allocatable :: pft(:,:)
        integer :: alloc_stat, ptcl_ind
        ptcl_ind = self%ptcl_ind(iptcl)
        allocate(pft(self%refsz,self%kfromto(1):self%kfromto(2)),&
        source=self%pfts_ptcls(ptcl_ind,irot:irot+self%refsz-1,:), stat=alloc_stat)
        call alloc_err("In: get_ptcl_pft; simple_polarft_corrcalc", alloc_stat)
    end function

    !>  \brief  returns polar Fourier transform of reference iref
    function get_ref_pft( self, iref ) result( pft )
        class(polarft_corrcalc), intent(in) :: self
        integer, intent(in)  :: iref
        complex, allocatable :: pft(:,:)
        integer :: alloc_stat
        allocate(pft(self%refsz,self%kfromto(1):self%kfromto(2)),&
        source=self%pfts_refs(iref,:,:), stat=alloc_stat)
        call alloc_err("In: get_ref_pft; simple_polarft_corrcalc", alloc_stat)
    end function

    !>  \brief  sets a reference Fourier component
    subroutine set_ref_fcomp( self, iref, irot, k, comp )
        class(polarft_corrcalc), intent(inout) :: self
        integer, intent(in)                    :: iref, irot, k
        complex, intent(in)                    :: comp
        integer :: k_ind
        k_ind = self%k_ind(k)
        self%pfts_refs(iref,irot,k_ind) = comp
    end subroutine

    !>  \brief  sets a particle Fourier component
    subroutine set_ptcl_fcomp( self, iptcl, irot, k, comp )
        class(polarft_corrcalc), intent(inout) :: self
        integer, intent(in)                    :: iptcl, irot, k
        complex, intent(in)                    :: comp
        integer                                :: ptcl_ind, k_ind
        ptcl_ind = self%ptcl_ind(iptcl)
        k_ind    = self%k_ind(k)
        self%pfts_ptcls(ptcl_ind,irot,k_ind) = comp
        self%pfts_ptcls(ptcl_ind,irot+self%nrots,k_ind) = comp ! because rot dim is expanded
    end subroutine

    !>  \brief  sets reference pft iref
    subroutine set_ref_pft( self, iref, pft )
        class(polarft_corrcalc), intent(inout) :: self
        integer, intent(in) :: iref
        complex, intent(in) :: pft(:,:)
        self%pfts_refs(iref,:,:) = pft(:self%refsz,:)
        ! calculate the square sum required for correlation calculation
        call self%memoize_sqsum_ref(iref)
    end subroutine

    !>  \brief  sets particle pft iptcl
    subroutine set_ptcl_pft( self, iptcl, pft )
        class(polarft_corrcalc), intent(inout) :: self
        integer, intent(in) :: iptcl
        complex, intent(in) :: pft(:,:)
        integer :: ptcl_ind
        ptcl_ind = self%ptcl_ind(iptcl)
        self%pfts_ptcls(ptcl_ind,:self%nrots,:)   = pft
        self%pfts_ptcls(ptcl_ind,self%nrots+1:,:) = pft ! because rot dim is expanded
        ! calculate the square sum required for correlation calculation
        call self%memoize_sqsum_ptcl(iptcl)
    end subroutine

    !>  \brief  for printing info about the object
    subroutine print( self )
        class(polarft_corrcalc), intent(in) :: self
        write(*,*) "from/to particle indices: ",                       self%pfromto
        write(*,*) "total number of particles in partition: ",         self%nptcls
        write(*,*) "number of references: ",                           self%nrefs
        write(*,*) "number of rotations: ",                            self%nrots
        write(*,*) "radius of molecule: ",                             self%ring2
        write(*,*) "nr of rots for ref (2nd dim of pftmat): ",         self%refsz
        write(*,*) "nr of rots for ptcl (2nd dim of pftmat): ",        self%ptclsz
        write(*,*) "logical dimensions of original Cartesian image: ", self%ldim
        write(*,*) "high-pass limit Fourier index",                    self%kfromto(1)
        write(*,*) "low-pass limit Fourier index",                     self%kfromto(2)
    end subroutine

    ! I/O

    !>  \brief  for writing the references generated to disk
!     subroutine write_refs( self, fname )
!         use simple_jiffys, only: get_fileunit
!         class(polarft_corrcalc), intent(in) :: self
!         character(len=*), intent(in)        :: fname
!         integer :: recsz, iref, funit
!         inquire( iolength=recsz ) self%pfts_refs(1,:,:)
!         funit = get_fileunit()
!         open(unit=funit, action='write', status='replace', file=fname,&
!         access='direct', form='unformatted', recl=recsz)
!         do iref=1,self%nrefs
!             write(funit,rec=iref) self%pfts_refs(iref,:,:)
!         end do
!         close(funit)
!     end subroutine
!
!     !>  \brief  for reading references from disk
!     subroutine read_refs( self, fname )
!         use simple_jiffys, only: get_fileunit
!         class(polarft_corrcalc), intent(inout) :: self
!         character(len=*), intent(in)           :: fname
!         integer :: recsz, iref, funit
!         inquire( iolength=recsz ) self%pfts_refs(1,:,:)
!         funit = get_fileunit()
!         open(unit=funit,  action='read', status='old', file=fname,&
!         access='direct', form='unformatted', recl=recsz)
!         do iref=1,self%nrefs
!             read(funit,rec=iref) self%pfts_refs(iref,:,:)
!         end do
!         close(funit)
!     end subroutine

    ! CHECKUPS

    !>  \brief  checks for existence
    function exists( self ) result( yes )
        class(polarft_corrcalc), intent(in) :: self
        logical :: yes
        yes = self%existence
    end function

    !>  \brief  is for plotting a particle polar FT
    subroutine vis_ptcl( self, iptcl )
        use gnufor2
        class(polarft_corrcalc), intent(in) :: self
        integer, intent(in) :: iptcl
        integer :: ptcl_ind
        ptcl_ind = self%ptcl_ind(iptcl)
        call gnufor_image(real(self%pfts_ptcls(ptcl_ind,:self%refsz,:)),  palette='gray')
        call gnufor_image(aimag(self%pfts_ptcls(ptcl_ind,:self%refsz,:)), palette='gray')
    end subroutine

    !>  \brief  is for plotting a particle polar FT
    subroutine vis_ref( self, iref )
        use gnufor2
        class(polarft_corrcalc), intent(in) :: self
        integer, intent(in) :: iref
        if( allocated(self%pfts_refs_ctf) )then
            call gnufor_image(real(self%pfts_refs_ctf(iref,:,:)),  palette='gray')
            call gnufor_image(aimag(self%pfts_refs_ctf(iref,:,:)), palette='gray')
        else
            call gnufor_image(real(self%pfts_refs(iref,:,:)),  palette='gray')
            call gnufor_image(aimag(self%pfts_refs(iref,:,:)), palette='gray')
        endif
    end subroutine

>>>>>>> Stashed changes
    ! MEMOIZERS

    !>  \brief  is for memoization of the complex square sums reqired for correlation calculation
    subroutine memoize_sqsum_ref( self, iref )
        use simple_math, only: csq
        class(polarft_corrcalc), intent(inout) :: self
        integer, intent(in) :: iref
        self%sqsums_refs(iref) = sum(csq(self%pfts_refs(iref,:,:)))
    end subroutine

    !>  \brief  is for memoization of the complex square sums reqired for correlation calculation
    subroutine memoize_sqsum_ref_ctf( self, iref )
        use simple_math, only: csq
        class(polarft_corrcalc), intent(inout) :: self
        integer, intent(in) :: iref
        if( self%xfel )then
            self%sqsums_refs(iref) = sum(real(self%pfts_refs_ctf(iref,:,:))**2.)
        else
            self%sqsums_refs(iref) = sum(csq(self%pfts_refs_ctf(iref,:,:)))
        endif
    end subroutine

    !>  \brief  is for memoization of the complex square sums reqired for correlation calculation
    subroutine memoize_sqsum_ptcl( self, iptcl )
        use simple_math, only: csq
        class(polarft_corrcalc), intent(inout) :: self
        integer, intent(in) :: iptcl
        integer :: ptcl_ind
        ptcl_ind = self%ptcl_ind(iptcl)
        if( self%xfel )then
            self%sqsums_ptcls(ptcl_ind) = sum(real(self%pfts_ptcls(ptcl_ind,:self%refsz,:))**2.)
        else
            self%sqsums_ptcls(ptcl_ind) = sum(csq(self%pfts_ptcls(ptcl_ind,:self%refsz,:)))
        endif
    end subroutine
<<<<<<< Updated upstream
    
    ! MODIFIERS
   
    !>  \brief  is for expanding the projection direction dimension to allow efficient
    !!          transfer to the GPU via subsectioning
    subroutine expand_dim( self )
      !$ use omp_lib
      !$ use omp_lib_kinds
      use simple_jiffys, only: alloc_err
      class(polarft_corrcalc), intent(inout) :: self
      complex, allocatable :: pfts_ptcls_tmp(:,:,:)
      real, allocatable    :: sqsums_ptcls_tmp(:)
      integer :: alloc_stat ! shapearr(3), 
      if( size(self%pfts_ptcls,dim=1) == 2*self%nptcls )then
         ! the arrays are already expanded, all we need to do is to copy data over
         !$omp parallel workshare default(shared)
         self%pfts_ptcls(self%nptcls+1:,:,:) = self%pfts_ptcls(:self%nptcls,:,:)
         self%sqsums_ptcls(self%nptcls+1:)   = self%sqsums_ptcls(:self%nptcls)
         !$omp end parallel workshare
      else if( size(self%pfts_ptcls,dim=1) == self%nptcls )then

         allocate( pfts_ptcls_tmp(self%nptcls,self%ptclsz,self%nk), source=self%pfts_ptcls, stat=alloc_stat )

         call alloc_err("In: expand_dim, simple_polarft_corrcalc, 1", alloc_stat)

         allocate( sqsums_ptcls_tmp(self%nptcls), source=self%sqsums_ptcls, stat=alloc_stat )

         call alloc_err("In: expand_dim, simple_polarft_corrcalc, 2", alloc_stat)
         deallocate(self%pfts_ptcls)
         deallocate(self%sqsums_ptcls)
            
         allocate( self%pfts_ptcls(2*self%nptcls,self%ptclsz,self%nk), stat=alloc_stat )
         call alloc_err("In: expand_dim, simple_polarft_corrcalc, 3", alloc_stat)

         allocate( self%sqsums_ptcls(2*self%nptcls), stat=alloc_stat )
         call alloc_err("In: expand_dim, simple_polarft_corrcalc, 4", alloc_stat)
         !$omp parallel workshare default(shared)
         self%pfts_ptcls(:self%nptcls,:,:)   = pfts_ptcls_tmp
         self%pfts_ptcls(self%nptcls+1:,:,:) = pfts_ptcls_tmp
         self%sqsums_ptcls(:self%nptcls)     = sqsums_ptcls_tmp
         self%sqsums_ptcls(self%nptcls+1:)   = sqsums_ptcls_tmp
         !$omp end parallel workshare
      else
         stop 'nonconforming dimension (1) of self%pfts_ptcls; expand_dim; simple_polarft_corrcalc'
      endif
      return
    end subroutine expand_dim

    ! CALCULATORS

    !>  \brief  is for generating all rotational correlations for all nptcls particles and all nrefs
    !           references (nptcls=nrefs) (assumes that the first and second dimensions of the particle
    !           matrix pfts_ptcls have been expanded)
    subroutine gencorrs_all( self, corrmat2dout, inplmat2dout )
=======

    ! MODIFIERES

    !>  \brief  is for generating a matrix of abs(CTF) values
    subroutine apply_ctf( self, tfun, dfx, dfy, angast )
        !$ use omp_lib
        !$ use omp_lib_kinds
        use simple_ctf, only: ctf
        class(polarft_corrcalc), intent(inout) :: self
        class(ctf), intent(inout)              :: tfun
        real, intent(in)                       :: dfx
        real, intent(in), optional             :: dfy, angast
        real :: absctfmat(self%refsz,self%nk)
        real :: inv_ldim(3),hinv,kinv,spaFreqSq,ddfy,aangast,ang
        integer :: irot, k, k_ind, iref
        ddfy = dfx
        if( present(dfy) ) ddfy = dfy
        aangast = 0.
        if( present(angast) ) aangast = angast
        inv_ldim = 1./real(self%ldim)
        ! CREATE CONFORMING MATRIX OF ABS(CTF) VALUES
        !$omp parallel default(shared) private(irot,k,k_ind,hinv,kinv,spaFreqSq,ang,iref)
        !$omp do schedule(auto)
        do irot=1,self%refsz
            do k=self%kfromto(1),self%kfromto(2)
                k_ind = self%k_ind(k)
                hinv   = self%polar(irot,k_ind)*inv_ldim(1)
                kinv   = self%polar(irot+self%nrots,k_ind)*inv_ldim(2)
                spaFreqSq = hinv*hinv+kinv*kinv
                ang    = atan2(self%polar(irot+self%nrots,k_ind),self%polar(irot,k_ind))
                absctfmat(irot,k_ind) = abs(tfun%eval(spaFreqSq,dfx,ddfy,aangast,ang))
            end do
        end do
        !$omp end do
        ! MULTIPLY THE REFERENCES WITH THE CTF
        !$omp do schedule(auto)
        do iref=1,self%nrefs
            self%pfts_refs_ctf(iref,:,:) = self%pfts_refs(iref,:,:)*absctfmat
            call self%memoize_sqsum_ref_ctf(iref)
        end do
        !$omp end do nowait
        !$omp end parallel
    end subroutine

    !>  \brief  is for copying the references to the CTF location
    subroutine copy_refs2ctf( self )
        class(polarft_corrcalc), intent(inout) :: self
        self%pfts_refs_ctf = self%pfts_refs
    end subroutine

    !>  \brief  is for expanding the projection direction dimension to allow efficient
    !!          transfer to the GPU via matrix subsectioning
    subroutine expand_dim( self )
>>>>>>> Stashed changes
        !$ use omp_lib
        !$ use omp_lib_kinds
        use simple_jiffys, only: progress
        use simple_jiffys, only: alloc_err
<<<<<<< Updated upstream
        class(polarft_corrcalc), intent(inout) :: self                         !< instance
=======
        class(polarft_corrcalc), intent(inout) :: self
        complex, allocatable :: pfts_ptcls_tmp(:,:,:)
        real, allocatable    :: sqsums_ptcls_tmp(:)
        integer :: alloc_stat ! shapearr(3),
        if( size(self%pfts_ptcls,dim=1) == 2*self%nptcls )then
            ! the arrays are already expanded, all we need to do is to copy data over
            !$omp parallel workshare default(shared)
            self%pfts_ptcls(self%nptcls+1:,:,:) = self%pfts_ptcls(:self%nptcls,:,:)
            self%sqsums_ptcls(self%nptcls+1:)   = self%sqsums_ptcls(:self%nptcls)
            !$omp end parallel workshare
        else if( size(self%pfts_ptcls,dim=1) == self%nptcls )then
            allocate( pfts_ptcls_tmp(self%nptcls,self%ptclsz,self%nk), source=self%pfts_ptcls, stat=alloc_stat )
            call alloc_err("In: expand_dim, simple_polarft_corrcalc, 1", alloc_stat)
            allocate( sqsums_ptcls_tmp(self%nptcls), source=self%sqsums_ptcls, stat=alloc_stat )
            call alloc_err("In: expand_dim, simple_polarft_corrcalc, 2", alloc_stat)
            deallocate(self%pfts_ptcls)
            deallocate(self%sqsums_ptcls)
            allocate( self%pfts_ptcls(2*self%nptcls,self%ptclsz,self%nk), stat=alloc_stat )
            call alloc_err("In: expand_dim, simple_polarft_corrcalc, 3", alloc_stat)
            allocate( self%sqsums_ptcls(2*self%nptcls), stat=alloc_stat )
            call alloc_err("In: expand_dim, simple_polarft_corrcalc, 4", alloc_stat)
            !$omp parallel workshare default(shared)
            self%pfts_ptcls(:self%nptcls,:,:)   = pfts_ptcls_tmp
            self%pfts_ptcls(self%nptcls+1:,:,:) = pfts_ptcls_tmp
            self%sqsums_ptcls(:self%nptcls)     = sqsums_ptcls_tmp
            self%sqsums_ptcls(self%nptcls+1:)   = sqsums_ptcls_tmp
            !$omp end parallel workshare
        else
            stop 'nonconforming dimension (1) of self%pfts_ptcls; expand_dim; simple_polarft_corrcalc'
        endif
    end subroutine
>>>>>>> Stashed changes

        real, intent(out)             :: corrmat2dout(self%nptcls,self%nptcls) !< output correlation matrix
        integer, intent(out)          :: inplmat2dout(self%nptcls,self%nptcls) !< output inplane rot index matrix

<<<<<<< Updated upstream
        real    ::     corrmat2d(self%nptcls,self%nptcls)
        real    ::     corrmat3d(self%nptcls,self%nptcls,self%nrots)
=======
    !>  \brief  is for generating all rotational correlations for all nptcls particles and all nrefs references (nptcls=nrefs)
    !!          (assumes that the first and second dimensions of the particle matrix pfts_ptcls have been expanded)
    subroutine gencorrs_all( self, corrmat2dout, inplmat2dout )
        !$ use omp_lib
        !$ use omp_lib_kinds
        use simple_jiffys, only: progress
        class(polarft_corrcalc), intent(inout) :: self                                  !< instance
        real, intent(out)                      :: corrmat2dout(self%nptcls,self%nptcls) !< output correlation matrix
        integer, intent(out)                   :: inplmat2dout(self%nptcls,self%nptcls) !< output inplane rot index matrix
        real    :: corrmat3d(self%nptcls,self%nptcls,self%nrots)
        real    :: corrmat2d(self%nptcls,self%nptcls)
>>>>>>> Stashed changes
        real    :: hadamard_prod(self%nptcls,self%refsz,self%nk)
        real(sp), allocatable         :: sumb_vec(:)
        complex(sp), allocatable      :: temp(:,:,:)

        integer :: ptcl_ind, iptcl, iref, irot, d1lim(2), d2lim(2), indices(self%nptcls)
        integer :: inplmat2d(self%nptcls,self%nptcls)

        !kernel variables
        integer                       :: err
        integer                       :: lda,ldb,ldc
        real(sp)                      :: alpha
        integer :: get_polarft_corr_gpu_c_

        !timming variables
        double precision              :: elps_L
        double precision,dimension(2) :: st_L_r, et_L_r

        !temporray index variables
        integer :: i,j,k
        
        if( allocated(self%pfts_refs_ctf) )then
            stop 'CTF modulation of the references not possible with simple_polarft_corrcalc :: gencorrs_all'
        endif

        allocate(sumb_vec(self%nptcls))
        allocate(temp(self%nptcls,self%refsz,self%nk))
        
        corrmat3d     = 0.
        corrmat2d     = 0.
        hadamard_prod = 0.
        sumb_vec      = 0.0
        temp          = 0.0
        inplmat2d     = 0
        self%s_polar%r_polar         = R_POLAR_D
        self%s_polar%sumasq_polar    = SUMASQ_POLAR_D
        self%s_polar%sumbsq_polar    = SUMBSQ_POLAR_D
        self%s_polar%ikrnl           = KERNEL_D
        self%s_polar%threadsPerBlock = THREADSPERBLOCK_D
        self%s_polar%nx              = NX_3D_D
        self%s_polar%ny              = NY_3D_D
        self%s_polar%nz              = NZ_3D_D
        alpha         = 1.0
        lda           = self%nptcls
        ldb           = self%refsz
        ldc           = self%nptcls
        
        do ptcl_ind=1,self%nptcls
            d1lim(1) = ptcl_ind
            d1lim(2) = ptcl_ind+self%nptcls-1
            do irot =1,self%nrots
                d2lim(1) = irot
                d2lim(2) = irot+self%refsz-1
#if defined (CUDA) && defined (MAGMA)
<<<<<<< Updated upstream
=======
                ! (1) do the Hadamard product on the GPU

                !$omp parallel default(shared) private(iptcl)
                !$omp workshare
                hadamard_prod = real(self%pfts_refs(:,:,:)*&
                conjg(self%pfts_ptcls(d1lim(1):d1lim(2),d2lim(1):d2lim(2),:)))
                !$omp end workshare
>>>>>>> Stashed changes

! GPU starts here
                temp = 0.0
                temp = self%pfts_ptcls(d1lim(1):d1lim(2),d2lim(1):d2lim(2),:)

                err = get_polarft_corr_gpu_c_(self%s_polar,"X","N",              &
                                              sumb_vec,                          &
                                              self%pfts_refs,                    &
                                              temp,                              &
                                              self%nptcls,self%refsz,self%nk,    &
                                              lda,ldb,ldc,alpha)
! GPU ends here

                !$omp parallel default(shared) private(iptcl)
                
                !$omp workshare
                corrmat3d(:,ptcl_ind,irot) = sumb_vec(:)/&
                sqrt(self%sqsums_refs(:)*self%sqsums_ptcls(d1lim(1):d1lim(2)))
                !$omp end workshare nowait

                !$omp end parallel
<<<<<<< Updated upstream
#else                                
                write(*,*) "#elseif defined (CUDA) && defined (MAGMA)"
=======
#else
>>>>>>> Stashed changes
                !$omp parallel default(shared) private(iptcl)
                !$omp workshare
                hadamard_prod = real(self%pfts_refs(:,:,:)*&
                conjg(self%pfts_ptcls(d1lim(1):d1lim(2),d2lim(1):d2lim(2),:)))
                !$omp end workshare

                ! correlation normalization

                !$omp do schedule(auto)
                do iptcl=1,self%nptcls
                    corrmat3d(iptcl,ptcl_ind,irot) = sum(hadamard_prod(iptcl,:,:))
                end do
                !$omp end do
                
                !$omp workshare
                corrmat3d(:,ptcl_ind,irot) = corrmat3d(:,ptcl_ind,irot)/&
                sqrt(self%sqsums_refs(:)*self%sqsums_ptcls(d1lim(1):d1lim(2)))
                !$omp end workshare nowait
                !$omp end parallel
#endif
            end do
        end do
        
        deallocate(sumb_vec)
        deallocate(temp)

        !$omp parallel default(shared) private(iptcl)
        !$omp workshare
        corrmat2d = maxval(corrmat3d, dim=3)
        inplmat2d = maxloc(corrmat3d, dim=3)
        !$omp end workshare
        !$omp do schedule(auto)
        do iptcl=1,self%nptcls
            indices(iptcl) = iptcl
        end do
        !$omp end do nowait

#if defined (BENCH)
        call gettimeofday_c(st_L_r)
#endif
        !$omp end parallel
        do iref=1,self%nptcls
<<<<<<< Updated upstream
           if( iref /= 1 )then
              indices = cshift(indices, shift=1)
           endif
           !$omp parallel do schedule(auto) private(iptcl)
           do iptcl=1,self%nptcls
              corrmat2dout(indices(iptcl),iref) = corrmat2d(iref,iptcl)
              inplmat2dout(indices(iptcl),iref) = inplmat2d(iref,iptcl)
           end do
           !$omp end parallel do
        end do
#if defined (BENCH)
        call gettimeofday_c(et_L_r)
        call elapsed_time_c(st_L_r,et_L_r,elps_L)
        write(*,*)"Matrix mapping"
        write(*,*)"---------------------------------------------------"
        write(*,'(f15.8)') elps_L
#endif

    end subroutine gencorrs_all
    
=======
            if( iref /= 1 )then
                indices = cshift(indices, shift=1)
            endif
            !$omp parallel do schedule(auto) private(iptcl)
            do iptcl=1,self%nptcls
                corrmat2dout(indices(iptcl),iref) = corrmat2d(iref,iptcl)
                inplmat2dout(indices(iptcl),iref) = inplmat2d(iref,iptcl)
            end do
            !$omp end parallel do
        end do
    end subroutine

>>>>>>> Stashed changes
    !>  \brief  tester routine for generating all rotational correlations for all nptcls particles and all nrefs references
    subroutine gencorrs_all_tester( self, corrmat2dout, inplmat2dout )
        use simple_jiffys, only: progress
        class(polarft_corrcalc), intent(inout) :: self               !< instance
        real, intent(out)    :: corrmat2dout(self%nptcls,self%nrefs) !< output correlation matrix
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
    end subroutine

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
    end function

<<<<<<< Updated upstream
=======
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
    end subroutine

>>>>>>> Stashed changes
    !>  \brief  for calculating the correlation between reference iref and particle iptcl in rotation irot
    function corr_1( self, iref, iptcl, irot ) result( cc )
        class(polarft_corrcalc), intent(inout) :: self              !< instance
        integer, intent(in)                    :: iref, iptcl, irot !< reference, particle, rotation
        real    :: cc
        integer :: ptcl_ind
        ptcl_ind = self%ptcl_ind(iptcl)
        if( allocated(self%pfts_refs_ctf) )then
<<<<<<< Updated upstream
           cc = sum( real(self%pfts_refs_ctf(iref,:,:) * &
                     conjg(  self%pfts_ptcls(ptcl_ind,irot:irot+self%refsz-1,:))))
=======
            cc = sum(real(self%pfts_refs_ctf(iref,:,:)*conjg(self%pfts_ptcls(ptcl_ind,irot:irot+self%refsz-1,:))))
        else if(self%xfel)then
            cc = sum(real(self%pfts_refs(iref,:,:)*real(self%pfts_ptcls(ptcl_ind,irot:irot+self%refsz-1,:))))
!           print *, 'unnorm cc:', cc
>>>>>>> Stashed changes
        else
           cc = sum( real(  self%pfts_refs(iref,:,:) * &
                     conjg(self%pfts_ptcls(ptcl_ind,irot:irot+self%refsz-1,:))))
        endif
!        print *, 'sqsums iref:',     self%sqsums_refs(iref)
!        print *, 'sqsums_ptcl_ind:', self%sqsums_ptcls(ptcl_ind)
        cc = cc/sqrt(self%sqsums_refs(iref)*self%sqsums_ptcls(ptcl_ind))
<<<<<<< Updated upstream
    end function corr_1
    
=======
!        print *, 'norm cc:', cc
    end function

>>>>>>> Stashed changes
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
            !$omp parallel workshare default(shared)
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
            !$omp parallel workshare default(shared)
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
    end subroutine

    ! CTF MODIFIERS

    !>  \brief  is for generating a matrix of abs(CTF) values
    subroutine apply_ctf( self, tfun, dfx, dfy, angast )
        !$ use omp_lib
        !$ use omp_lib_kinds
        use simple_ctf, only: ctf
        class(polarft_corrcalc), intent(inout) :: self
        class(ctf), intent(inout)              :: tfun
        real, intent(in)                       :: dfx
        real, intent(in), optional             :: dfy, angast
        real :: absctfmat(self%refsz,self%nk)
        real :: inv_ldim(3),hinv,kinv,spaFreqSq,ddfy,aangast,ang
        integer :: irot, k, k_ind, iref
        ddfy = dfx
        if( present(dfy) ) ddfy = dfy
        aangast = 0.
        if( present(angast) ) aangast = angast       
        inv_ldim = 1./real(self%ldim)
        ! CREATE CONFORMING MATRIX OF ABS(CTF) VALUES
        !$omp parallel default(shared) private(irot,k,k_ind,hinv,kinv,spaFreqSq,ang,iref)
        !$omp do schedule(auto)
        do irot=1,self%refsz
            do k=self%kfromto(1),self%kfromto(2)
                k_ind = self%k_ind(k)
                hinv   = self%polar(irot,k_ind)*inv_ldim(1)
                kinv   = self%polar(irot+self%nrots,k_ind)*inv_ldim(2)
                spaFreqSq = hinv*hinv+kinv*kinv
                ang    = atan2(self%polar(irot+self%nrots,k_ind),self%polar(irot,k_ind))
                absctfmat(irot,k_ind) = abs(tfun%eval(spaFreqSq,dfx,ddfy,aangast,ang))
            end do
        end do
        !$omp end do
        ! MULTIPLY THE REFERENCES WITH THE CTF
        !$omp do schedule(auto)
        do iref=1,self%nrefs
            self%pfts_refs_ctf(iref,:,:) = self%pfts_refs(iref,:,:)*absctfmat
            call self%memoize_sqsum_ref_ctf(iref)
        end do
        !$omp end do nowait
        !$omp end parallel
    end subroutine
    
    !>  \brief  is for copying the references to the CTF location
    subroutine copy_refs2ctf( self )
        class(polarft_corrcalc), intent(inout) :: self
        self%pfts_refs_ctf = self%pfts_refs
    end subroutine
    
    ! SETTERS AND GETTERS

    !>  \brief  for mapping physical Fourier index to implemented
    function k_ind( self, k_in ) result( k_out )
        class(polarft_corrcalc), intent(in) :: self
        integer, intent(in)                 :: k_in
        integer :: k_out
        k_out = k_in-self%kfromto(1)+1
    end function
    
    !>  \brief  for mapping physical particle index to implemented
    function ptcl_ind( self, iptcl_in ) result( iptcl_out )
        class(polarft_corrcalc), intent(in) :: self
        integer, intent(in)                 :: iptcl_in
        integer :: iptcl_out
        iptcl_out = iptcl_in-self%pfromto(1)+1
    end function ptcl_ind
    
    !>  \brief  for getting the logical particle range
    function get_pfromto( self ) result( lim )
        class(polarft_corrcalc), intent(inout) :: self
        integer :: lim(2)
        lim = self%pfromto
    end function
    
    !>  \brief  for getting the number of particles
    pure function get_nptcls( self ) result( nptcls )
        class(polarft_corrcalc), intent(in) :: self
        integer :: nptcls
        nptcls = self%nptcls
    end function
    
    !>  \brief  for getting the number of references
    pure function get_nrefs( self ) result( nrefs )
        class(polarft_corrcalc), intent(in) :: self
        integer :: nrefs
        nrefs = self%nrefs
    end function
    
    !>  \brief  for getting the number of in-plane rotations
    pure function get_nrots( self ) result( nrots )
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
    
    !>  \brief  for getting the number of reference rotations (size of second dim of self%pfts_refs)
    function get_refsz( self ) result( refsz )
        class(polarft_corrcalc), intent(in) :: self
        integer :: refsz
        refsz = self%refsz
    end function
    
    !>  \brief  for getting the number of particle rotations (size of second dim of self%pfts_ptcls)
    function get_ptclsz( self ) result( ptclsz )
        class(polarft_corrcalc), intent(in) :: self
        integer :: ptclsz
        ptclsz = self%ptclsz
    end function
    
    !>  \brief  for getting the logical dimension of the original
    !!          Cartesian image
    function get_ldim( self ) result( ldim )
        class(polarft_corrcalc), intent(in) :: self
        integer :: ldim(3)
        ldim = self%ldim
    end function
    
    !>  \brief  for getting the Fourier index range (hp/lp)
    function get_kfromto( self ) result( lim )
        class(polarft_corrcalc), intent(inout) :: self
        integer :: lim(2)
        lim = self%kfromto
    end function
    
    !>  \brief is for getting the continuous in-plane rotation
    !!         corresponding to in-plane rotation index roind
    function get_rot( self, roind ) result( rot )
        class(polarft_corrcalc), intent(in) :: self
        integer, intent(in) :: roind
        real :: rot
        if( roind < 1 .or. roind > self%nrots )then
            stop 'roind is out of range; get_rot; simple_polarft_corrcalc'
        endif
        rot = self%angtab(roind)
    end function

    !>  \brief is for getting the discrete in-plane rotational
    !!         index corresponding to continuous rotation rot
    function get_roind( self, rot ) result( ind )
        class(polarft_corrcalc), intent(in) :: self
        real, intent(in) :: rot
        integer :: ind, irot, alloc_stat, loc(1)
        real    :: dists(self%nrots)
        do irot=1,self%nrots
            dists(irot) = sqrt((self%angtab(irot)-rot)**2.)
        end do
        loc = minloc(dists)
        ind = loc(1)
    end function

    !>  \brief returns polar coordinate for rotation rot
    !!         and Fourier index k
    function get_coord( self, rot, k ) result( xy )
        class(polarft_corrcalc), intent(in) :: self
        integer, intent(in)                 :: rot, k
        real    :: xy(2)
        integer :: k_ind
        k_ind = self%k_ind(k)    
        xy(1) = self%polar(rot,k_ind)
        xy(2) = self%polar(self%nrots+rot,k_ind)
    end function
    
    !>  \brief  returns polar Fourier transform of particle iptcl in rotation irot
    function get_ptcl_pft( self, iptcl, irot ) result( pft )
        class(polarft_corrcalc), intent(in) :: self
        integer, intent(in)  :: iptcl, irot
        complex, allocatable :: pft(:,:)
        integer :: alloc_stat, ptcl_ind
        ptcl_ind = self%ptcl_ind(iptcl)
        allocate(pft(self%refsz,self%kfromto(1):self%kfromto(2)),&
        source=self%pfts_ptcls(ptcl_ind,irot:irot+self%refsz-1,:), stat=alloc_stat)
        call alloc_err("In: get_ptcl_pft; simple_polarft_corrcalc", alloc_stat)
    end function
    
    !>  \brief  returns polar Fourier transform of reference iref
    function get_ref_pft( self, iref ) result( pft )
        class(polarft_corrcalc), intent(in) :: self
        integer, intent(in)  :: iref
        complex, allocatable :: pft(:,:)
        integer :: alloc_stat
        allocate(pft(self%refsz,self%kfromto(1):self%kfromto(2)),&
        source=self%pfts_refs(iref,:,:), stat=alloc_stat)
        call alloc_err("In: get_ref_pft; simple_polarft_corrcalc", alloc_stat)
    end function
    
    !>  \brief  sets a reference Fourier component
    subroutine set_ref_fcomp( self, iref, irot, k, comp )
        class(polarft_corrcalc), intent(inout) :: self
        integer, intent(in)                    :: iref, irot, k
        complex, intent(in)                    :: comp
        integer :: k_ind
        k_ind = self%k_ind(k)
        self%pfts_refs(iref,irot,k_ind) = comp
    end subroutine
    
    !>  \brief  sets a particle Fourier component
    subroutine set_ptcl_fcomp( self, iptcl, irot, k, comp )
        class(polarft_corrcalc), intent(inout) :: self
        integer, intent(in)                    :: iptcl, irot, k
        complex, intent(in)                    :: comp
        integer                                :: ptcl_ind, k_ind
        ptcl_ind = self%ptcl_ind(iptcl)
        k_ind    = self%k_ind(k)
        self%pfts_ptcls(ptcl_ind,irot,k_ind) = comp
        self%pfts_ptcls(ptcl_ind,irot+self%nrots,k_ind) = comp ! because rot dim is expanded
    end subroutine
    
    !>  \brief  sets reference pft iref
    subroutine set_ref_pft( self, iref, pft )
        class(polarft_corrcalc), intent(inout) :: self
        integer, intent(in) :: iref
        complex, intent(in) :: pft(:,:)
        self%pfts_refs(iref,:,:) = pft(:self%refsz,:)
        ! calculate the square sum required for correlation calculation
        call self%memoize_sqsum_ref(iref)
    end subroutine
    
    !>  \brief  sets particle pft iptcl
    subroutine set_ptcl_pft( self, iptcl, pft )
        class(polarft_corrcalc), intent(inout) :: self
        integer, intent(in) :: iptcl
        complex, intent(in) :: pft(:,:)
        integer :: ptcl_ind
        ptcl_ind = self%ptcl_ind(iptcl)
        self%pfts_ptcls(ptcl_ind,:self%nrots,:)   = pft
        self%pfts_ptcls(ptcl_ind,self%nrots+1:,:) = pft ! because rot dim is expanded
        ! calculate the square sum required for correlation calculation
        call self%memoize_sqsum_ptcl(iptcl)
    end subroutine set_ptcl_pft

    ! VISUALISATORS
    
    !>  \brief  is for plotting a particle polar FT
    subroutine vis_ptcl( self, iptcl )
        use gnufor2
        class(polarft_corrcalc), intent(in) :: self
        integer, intent(in) :: iptcl
        integer :: ptcl_ind
        ptcl_ind = self%ptcl_ind(iptcl)
        call gnufor_image(real(self%pfts_ptcls(ptcl_ind,:self%refsz,:)),  palette='gray')
        call gnufor_image(aimag(self%pfts_ptcls(ptcl_ind,:self%refsz,:)), palette='gray')
    end subroutine
    
    !>  \brief  is for plotting a particle polar FT
    subroutine vis_ref( self, iref )
        use gnufor2
        class(polarft_corrcalc), intent(in) :: self
        integer, intent(in) :: iref
        if( allocated(self%pfts_refs_ctf) )then
            call gnufor_image(real(self%pfts_refs_ctf(iref,:,:)),  palette='gray')
            call gnufor_image(aimag(self%pfts_refs_ctf(iref,:,:)), palette='gray')
        else
            call gnufor_image(real(self%pfts_refs(iref,:,:)),  palette='gray')
            call gnufor_image(aimag(self%pfts_refs(iref,:,:)), palette='gray')
        endif
    end subroutine

<<<<<<< Updated upstream
    ! PRINTERS
      
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

    ! CHECKUPS

    !>  \brief  checks for existence
    function exists( self ) result( yes )
        class(polarft_corrcalc), intent(in) :: self
        logical :: yes
        yes = self%existence
    end function
    
    ! DESTRUCTOR

    !>  \brief  is a destructor
    subroutine kill( self )
      class(polarft_corrcalc), intent(inout) :: self
      if( self%existence )then
         deallocate(self%sqsums_refs, self%sqsums_ptcls, self%angtab,&
              self%argtransf, self%polar, self%pfts_refs, self%pfts_ptcls )
         if( allocated(self%pfts_refs_ctf) ) deallocate(self%pfts_refs_ctf)
         self%existence = .false.
      endif
    end subroutine kill
    
    ! TESTERS

=======
>>>>>>> Stashed changes
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
        comparing the gencorrs correlations with the ones generated here'
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
        testing the gencorrs function'
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
        testing the gencorrs function with expand_dim'
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
                    conform with the old; test_polarft_corrcalc'
                endif
            end do
        end do
        write(*,'(a)') '**info(simple_polarft_corrcalc_unit_test, part 4):&
        testing the corr function with expand_dim'
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
    end subroutine

end module simple_polarft_corrcalc
