!>  \brief  SIMPLE prime2Dcalc class
module simple_prime2Dcalc
use simple_defs
use simple_params,   only: params
use simple_ran_tabu, only: ran_tabu
use simple_jiffys,   only: alloc_err
implicit none

public :: prime2Dcalc, test_prime2Dcalc
private

! CLASS PARAMETERS/VARIABLES
complex, parameter :: zero=cmplx(0.,0.) !< just a complex zero

type :: prime2Dcalc
    private
    integer                  :: nptcls     = 1        !< the total number of particles in partition (logically indexded [1,nptcls])
    integer                  :: nrefs      = 1        !< the number of references (logically indexded [1,nrefs])
    integer                  :: nrots      = 0        !< number of in-plane rotations for one pft (determined by radius of molecule)
    integer                  :: ring2      = 0        !< radius of molecule
    integer                  :: refsz      = 0        !< size of reference (nrots/2) (number of vectors used for matching)
    integer                  :: ptclsz     = 0        !< size of particle (2*nrots)
    integer                  :: ldim(3)    = 0        !< logical dimensions of original cartesian image
    integer                  :: kfromto(2) = 0        !< Fourier index range
    integer                  :: nk         = 0        !< number of resolution elements (Fourier components)
    real,        allocatable :: sqsums_refs(:)        !< memoized square sums for the correlation calculations
    real,        allocatable :: sqsums_ptcls(:)       !< memoized square sums for the correlation calculations
    real,        allocatable :: angtab(:)             !< table of in-plane angles (in degrees)
    real,        allocatable :: argtransf(:,:)        !< argument transfer constants for shifting the references
    real,        allocatable :: polar(:,:)            !< table of polar coordinates (in Cartesian coordinates)
    complex(sp), allocatable :: pfts_refs(:,:,:)      !< 3D complex matrix of polar reference sections (nrefs,refsz,nk)
    complex(sp), allocatable :: pfts_refs_ctf(:,:,:)  !< 3D complex matrix of polar reference sections with CTF applied
    complex(sp), allocatable :: pfts_sums(:,:,:)      !< 3D complex matrix of class sums
    real,        allocatable :: ctfsq_sums(:,:,:)     !< 3D complex matrix of CTF**2.0 terms
    complex(sp), allocatable :: pfts_ptcls(:,:,:)     !< 3D complex matrix of particle sections
    logical                  :: existence=.false.     !< to indicate existence
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

end type prime2Dcalc

contains


    !>  \brief  is a constructor
    subroutine new( self, nrefs, nptcls, ldim, kfromto, ring2, ctfflag )
        use simple_math, only: rad2deg, is_even, round2even
        class(prime2Dcalc), intent(inout) :: self
        integer,            intent(in)    :: nrefs, pfromto(2), ldim(3), kfromto(2), ring2
        character(len=*),   intent(in)    :: ctfflag
        integer :: alloc_stat, irot, k, k_ind, err
        logical :: even_dims, test(3)
        real    :: ang
        ! kill possibly pre-existing object
        call self%kill
        ! error check
        if( kfromto(2) - kfromto(1) <= 2 )then
            write(*,*) 'kfromto: ', kfromto(1), kfromto(2)
            stop 'resolution range too narrow; new; simple_prime2Dcalc'
        endif
        if( ring2 < 1 )then
            write(*,*) 'ring2: ', ring2
            stop 'ring2 must be > 0; new; simple_prime2Dcalc'
        endif
        if( pfromto(2) - pfromto(1) + 1 < 1 )then
            write(*,*) 'pfromto: ', pfromto(1), pfromto(2)
            stop 'nptcls (# of particles) must be > 0; new; simple_prime2Dcalc'
        endif
        if( nrefs < 1 )then
            write(*,*) 'nrefs: ', nrefs
            stop 'nrefs (# of reference sections) must be > 0; new; simple_prime2Dcalc'
        endif
        if( any(ldim == 0) )then
            write(*,*) 'ldim: ', ldim
            stop 'ldim is not conforming (is zero); new; simple_prime2Dcalc'
        endif
        if( ldim(3) > 1 )then
            write(*,*) 'ldim: ', ldim
            stop '3D polarfts are not yet supported; new; simple_prime2Dcalc'
        endif
        test    = .false.
        test(1) = is_even(ldim(1))
        test(2) = is_even(ldim(2))
        test(3) = ldim(3) == 1
        even_dims = all(test)
        if( .not. even_dims )then
            write(*,*) 'ldim: ', ldim
            stop 'only even logical dims supported; new; simple_prime2Dcalc'
        endif
        self%nptcls  = nptcls                        !< the total number of particles in partition (logically indexded [1,nptcls])
        self%nrefs   = nrefs                         !< the number of references (logically indexded [1,nrefs])
        self%ring2   = ring2                         !< radius of molecule
        self%nrots   = round2even(twopi*real(ring2)) !< number of in-plane rotations for one pft (determined by radius of molecule)
        self%refsz   = self%nrots/2                  !< size of reference (nrots/2) (number of vectors used for matching)
        self%ptclsz  = self%nrots*2                  !< size of particle (2*nrots)
        self%ldim    = ldim                          !< logical dimensions of original cartesian image
        self%kfromto = kfromto                       !< Fourier index range
        self%nk      = kfromto(2)-kfromto(1)+1       !< number of resolution elements (Fourier components)
        ! generate polar coordinates
        allocate( self%polar(self%ptclsz,self%nk), self%angtab(self%nrots), stat=alloc_stat)
        call alloc_err('polar coordinate arrays; new; simple_prime2Dcalc', alloc_stat)
        ang = twopi/real(self%nrots)
        do irot=1,self%nrots
            self%angtab(irot) = real(irot-1)*ang
            do k=self%kfromto(1),self%kfromto(2)
                k_ind = k-self%kfromto(1)+1
                self%polar(irot,k_ind)            = cos(self%angtab(irot))*real(k) ! x-coordinate
                self%polar(irot+self%nrots,k_ind) = sin(self%angtab(irot))*real(k) ! y-coordinate
            end do
            self%angtab(irot) = rad2deg(self%angtab(irot)) ! angle (in degrees)
        end do
        ! generate the argument transfer constants for shifting reference polarfts
        allocate( self%argtransf(self%nrots,self%nk), stat=alloc_stat)
        call alloc_err('shift argument transfer array; new; simple_prime2Dcalc', alloc_stat)
        self%argtransf(:self%refsz,:)   = &
            self%polar(:self%refsz,:)   * &
            (PI/real(self%ldim(1)/2))    ! x-part
        self%argtransf(self%refsz+1:,:) = &
            self%polar(self%nrots+1:self%nrots+self%refsz,:) * &
            (PI/real(self%ldim(2)/2))    ! y-part
        ! allocate polarfts, sums, and sqsums
        allocate(   self%pfts_refs(self%nrefs,self%refsz,self%nk),&
                    self%pfts_ptcls(self%nptcls,self%ptclsz,self%nk),&
                    self%pfts_sums(self%nrefs,self%refsz,self%nk),&
                    self%ctfsq_sums(self%nrefs,self%refsz,self%nk),&
                    self%sqsums_refs(self%nrefs),&
                    self%sqsums_ptcls(self%nptcls), stat=alloc_stat)
        call alloc_err('polarfts, sums, and sqsums; new; simple_prime2Dcalc', alloc_stat)
        self%pfts_refs    = zero
        self%pfts_ptcls   = zero
        self%pfts_sums    = zero
        self%ctfsq_sums   = 0.
        self%sqsums_refs  = 0.
        self%sqsums_ptcls = 0.
        ! pfts_refs_ctf if needed
        if( ctfflag .ne. 'no' )then
            allocate(self%pfts_refs_ctf(self%nrefs,self%refsz,self%nk), stat=alloc_stat)
            call alloc_err('pfts_refs_ctf; new; simple_prime2Dcalc', alloc_stat)
            self%pfts_refs_ctf = zero
        endif
        self%existence = .true.
    end subroutine new



end module simple_prime2Dcalc