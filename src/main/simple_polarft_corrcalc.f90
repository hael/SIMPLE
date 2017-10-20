! for calculation of band-pass limited cross-correlation of polar Fourier transforms
module simple_polarft_corrcalc
#include "simple_lib.f08"
use simple_params,   only: params
use simple_ran_tabu, only: ran_tabu
use simple_fftw3
!$ use omp_lib
!$ use omp_lib_kinds
implicit none

public :: polarft_corrcalc
private
#include "simple_local_flags.inc"

! CLASS PARAMETERS/VARIABLES
complex(sp), parameter :: zero=cmplx(0.,0.) !< just a complex zero

! the fftw_arrs data structures are needed for thread-safe FFTW exec. Letting OpenMP copy out the per-threads
! arrays leads to bugs because of inconsistency between data in memory and the fftw_plan
type fftw_carr
    type(c_ptr)                            :: p_re             !< pointer for C-style allocation
    type(c_ptr)                            :: p_im             !< -"-
    real(kind=c_float),            pointer :: re(:) => null()  !< corresponding Fortran pointers
    complex(kind=c_float_complex), pointer :: im(:) => null()  !< -"-
end type fftw_carr

type fftw_carr_fft
    type(c_ptr)                            :: p_re             !< pointer for C-style allocation
    type(c_ptr)                            :: p_im             !< -"-
    complex(kind=c_float_complex), pointer :: re(:) => null()  !< corresponding Fortran pointers
    complex(kind=c_float_complex), pointer :: im(:) => null()  !< -"-
end type fftw_carr_fft

type fftw_arrs
    type(c_ptr)                            :: p_ref_re                  !< pointer for C-style allocation
    type(c_ptr)                            :: p_ref_im                  !< -"-
    type(c_ptr)                            :: p_ref_fft_re              !< -"-
    type(c_ptr)                            :: p_ref_fft_im              !< -"-
    type(c_ptr)                            :: p_product_fft             !< -"-
    type(c_ptr)                            :: p_backtransf              !< -"-
    real(kind=c_float),            pointer :: ref_re(:)      => null()  !< corresponding Fortran pointers
    complex(kind=c_float_complex), pointer :: ref_im(:)      => null()  !< -"-
    complex(kind=c_float_complex), pointer :: ref_fft_re(:)  => null()  !< -"-
    complex(kind=c_float_complex), pointer :: ref_fft_im(:)  => null()  !< -"-
    complex(kind=c_float_complex), pointer :: product_fft(:) => null()  !< -"-
    real(kind=c_float),            pointer :: backtransf(:)  => null()  !< -"-
end type fftw_arrs

type :: polarft_corrcalc
  private
    integer                          :: pfromto(2) = 1        !< from/to particle indices (in parallel execution)
    integer                          :: nptcls     = 1        !< the total number of particles in partition (logically indexded [fromp,top])
    integer                          :: nrefs      = 1        !< the number of references (logically indexded [1,nrefs])
    integer                          :: nrots      = 0        !< number of in-plane rotations for one pft (determined by radius of molecule)
    integer                          :: ring2      = 0        !< radius of molecule
    integer                          :: pftsz      = 0        !< size of reference and particle pft (nrots/2)
    integer                          :: nk         = 0        !< # resolution elements in the band-pass limited PFTs
    integer                          :: nthr       = 0        !< # OpenMP threads
    integer                          :: winsz      = 0        !< size of moving window in correlation calculations
    integer                          :: ldim(3)    = 0        !< logical dimensions of original cartesian image
    integer                          :: kfromto(2) = 0        !< Fourier index range
    real(sp)                         :: smpd       = 0.       !< sampling distance
    real(sp),            allocatable :: sqsums_ptcls(:)       !< memoized square sums for the correlation calculations
    real(sp),            allocatable :: angtab(:)             !< table of in-plane angles (in degrees)
    real(sp),            allocatable :: argtransf(:,:)        !< argument transfer constants for shifting the references
    real(sp),            allocatable :: polar(:,:)            !< table of polar coordinates (in Cartesian coordinates)
    complex(sp),         allocatable :: pfts_refs_even(:,:,:) !< 3D complex matrix of polar reference sections (nrefs,pftsz,nk), even
    complex(sp),         allocatable :: pfts_refs_odd(:,:,:)  !< -"-, odd
    complex(sp),         allocatable :: pfts_ptcls(:,:,:)     !< 3D complex matrix of particle sections
    complex(sp),         allocatable :: fft_factors(:)        !< phase factors for accelerated gencorrs routines
    type(fftw_arrs),     allocatable :: fftdat(:)             !< arrays for accelerated gencorrs routines
    type(fftw_carr_fft), allocatable :: fftdat_ptcls(:,:)     !< for memoization of particle  FFTs in accelerated gencorrs routines
    type(fftw_carr_fft), allocatable :: fftdat_refs_even(:,:) !< for memoization of reference FFTs in accelerated gencorrs routines, even
    type(fftw_carr_fft), allocatable :: fftdat_refs_odd(:,:)  !< -"-, odd
    logical,             allocatable :: iseven(:)             !< eo assignment for gold-standard FSC
    logical                          :: phaseplate            !< images obtained with the Volta
    type(c_ptr)                      :: plan_fwd_1            !< FFTW plans for gencorrs
    type(c_ptr)                      :: plan_fwd_2            !< -"-
    type(c_ptr)                      :: plan_bwd              !< -"-
    logical                          :: existence = .false.   !< to indicate existence
  contains
    ! CONSTRUCTOR
    procedure          :: new
    ! SETTERS
    procedure          :: set_ref_pft
    procedure          :: set_ptcl_pft
    procedure          :: set_ref_fcomp
    procedure          :: set_ptcl_fcomp
    procedure          :: zero_ref
    procedure          :: cp_even2odd_ref
    ! GETTERS
    procedure          :: get_pfromto
    procedure          :: get_nptcls
    procedure          :: get_nrefs
    procedure          :: get_nrots
    procedure          :: get_ring2
    procedure          :: get_pftsz
    procedure          :: get_ldim
    procedure          :: get_kfromto
    procedure          :: get_pdim
    procedure          :: get_rot
    procedure          :: get_rots_for_applic
    procedure          :: get_roind
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
    procedure          :: memoize_ffts
    ! CALCULATORS
    procedure, private :: create_polar_ctfmat
    procedure          :: apply_ctf_to_ptcls
    procedure, private :: calc_corrs_over_k
    procedure, private :: calc_corrs_over_k_wrefmem
    procedure, private :: calc_k_corrs_wrefmem
    procedure, private :: gencorrs_1
    procedure, private :: gencorrs_2
    procedure, private :: gencorrs_3
    generic            :: gencorrs => gencorrs_1, gencorrs_2, gencorrs_3
    procedure          :: genfrc
    ! DESTRUCTOR
    procedure          :: kill
end type polarft_corrcalc

contains

    ! CONSTRUCTORS

    !>  \brief  is a constructor
    subroutine new( self, nrefs, p, eoarr )
        use simple_math,   only: rad2deg, is_even, round2even
        use simple_params, only: params
        class(polarft_corrcalc), intent(inout) :: self
        integer,                 intent(in)    :: nrefs
        class(params),           intent(inout) :: p
        integer, optional,       intent(in)    :: eoarr(p%fromp:p%top)
        integer  :: alloc_stat, irot, k, ithr, iptcl, ik, iref
        logical  :: even_dims, test(2)
        real(sp) :: ang
        integer(kind=c_int) :: wsdm_ret
        ! kill possibly pre-existing object
        call self%kill
        ! error check
        if( p%kfromto(2) - p%kfromto(1) <= 2 )then
            write(*,*) 'p%kfromto: ', p%kfromto(1), p%kfromto(2)
            call simple_stop( 'resolution range too narrow; new; simple_polarft_corrcalc')
        endif
        if( p%ring2 < 1 )then
            write(*,*) 'p%ring2: ', p%ring2
            call simple_stop ( 'p%ring2 must be > 0; new; simple_polarft_corrcalc')
        endif
        if( p%top - p%fromp + 1 < 1 )then
            write(*,*) 'pfromto: ', p%fromp, p%top
            call simple_stop ('nptcls (# of particles) must be > 0; new; simple_polarft_corrcalc')
        endif
        if( nrefs < 1 )then
            write(*,*) 'nrefs: ', nrefs
            call simple_stop ('nrefs (# of reference sections) must be > 0; new; simple_polarft_corrcalc')
        endif
        self%ldim = [p%boxmatch,p%boxmatch,1] !< logical dimensions of original cartesian image
        test    = .false.
        test(1) = is_even(self%ldim(1))
        test(2) = is_even(self%ldim(2))
        even_dims = all(test)
        if( .not. even_dims )then
            write(*,*) 'self%ldim: ', self%ldim
            call simple_stop ('only even logical dims supported; new; simple_polarft_corrcalc')
        endif
        ! set constants
        self%pfromto    = [p%fromp,p%top]                       !< from/to particle indices (in parallel execution)
        self%nptcls     = p%top - p%fromp + 1                   !< the total number of particles in partition (logically indexded [fromp,top])
        self%nrefs      = nrefs                                 !< the number of references (logically indexded [1,nrefs])
        self%ring2      = p%ring2                               !< radius of molecule
        self%nrots      = round2even(twopi * real(p%ring2))     !< number of in-plane rotations for one pft  (determined by radius of molecule)
        self%pftsz      = self%nrots / 2                        !< size of reference (nrots/2) (number of vectors used for matching)
        self%smpd       = p%smpd                                !< sampling distance
        self%kfromto    = p%kfromto                             !< Fourier index range
        self%nk         = self%kfromto(2) - self%kfromto(1) + 1 !< # resolution elements
        self%nthr       = p%nthr                                !< # OpenMP threads
        self%phaseplate = p%tfplan%l_phaseplate                 !< images obtained with the Volta
        ! generate polar coordinates & eo assignment
        allocate( self%polar(2*self%nrots,self%kfromto(1):self%kfromto(2)),&
                 &self%angtab(self%nrots), self%iseven(self%pfromto(1):self%pfromto(2)), stat=alloc_stat)
        allocchk('polar coordinate arrays; new; simple_polarft_corrcalc')
        ang = twopi/real(self%nrots)
        do irot=1,self%nrots
            self%angtab(irot) = real(irot-1)*ang
            do k=self%kfromto(1),self%kfromto(2)
                self%polar(irot,k)            = cos(self%angtab(irot))*real(k) ! x-coordinate
                self%polar(irot+self%nrots,k) = sin(self%angtab(irot))*real(k) ! y-coordinate
            end do
            self%angtab(irot) = rad2deg(self%angtab(irot)) ! angle (in degrees)
        end do
        ! eo assignment
        if( present(eoarr) )then
            if( all(eoarr == - 1) )then
                self%iseven = .true.
            else
                where(eoarr == 0)
                    self%iseven = .true.
                else where
                    self%iseven = .false.
                end where
            endif
        else
            self%iseven = .true.
        endif
        ! generate the argument transfer constants for shifting reference polarfts
        allocate( self%argtransf(self%nrots,self%kfromto(1):self%kfromto(2)), stat=alloc_stat)
        allocchk('shift argument transfer array; new; simple_polarft_corrcalc')
        self%argtransf(:self%pftsz,:)   = &
            self%polar(:self%pftsz,:)   * &
            (PI/real(self%ldim(1)/2))    ! x-part
        self%argtransf(self%pftsz+1:,:) = &
            self%polar(self%nrots+1:self%nrots+self%pftsz,:) * &
            (PI/real(self%ldim(2)/2))    ! y-part
        ! allocate others
        allocate( self%pfts_refs_even(self%nrefs,self%pftsz,self%kfromto(1):self%kfromto(2)),&
                 &self%pfts_refs_odd(self%nrefs,self%pftsz,self%kfromto(1):self%kfromto(2)),&
                 &self%pfts_ptcls(self%pfromto(1):self%pfromto(2),self%pftsz,self%kfromto(1):self%kfromto(2)),&
                 &self%sqsums_ptcls(self%pfromto(1):self%pfromto(2)), self%fftdat(self%nthr),&
                 &self%fftdat_ptcls(self%pfromto(1):self%pfromto(2),self%kfromto(1):self%kfromto(2)),&
                 &self%fftdat_refs_even(1:self%nrefs,self%kfromto(1):self%kfromto(2)),&
                 &self%fftdat_refs_odd(1:self%nrefs,self%kfromto(1):self%kfromto(2)), stat=alloc_stat)
        allocchk('polarfts and sqsums; new; simple_polarft_corrcalc')
        self%pfts_refs_even = zero
        self%pfts_refs_odd  = zero
        self%pfts_ptcls     = zero
        self%sqsums_ptcls   = 0.
        ! thread-safe c-style allocatables for gencorrs
        do ithr=1,self%nthr
            self%fftdat(ithr)%p_ref_re      = fftwf_alloc_real   (int(self%pftsz, c_size_t))
            self%fftdat(ithr)%p_ref_im      = fftwf_alloc_complex(int(self%pftsz, c_size_t))
            self%fftdat(ithr)%p_ref_fft_re  = fftwf_alloc_complex(int(self%pftsz, c_size_t))
            self%fftdat(ithr)%p_ref_fft_im  = fftwf_alloc_complex(int(self%pftsz, c_size_t))
            self%fftdat(ithr)%p_product_fft = fftwf_alloc_complex(int(self%nrots, c_size_t))
            self%fftdat(ithr)%p_backtransf  = fftwf_alloc_real   (int(self%nrots, c_size_t))
            call c_f_pointer(self%fftdat(ithr)%p_ref_re,      self%fftdat(ithr)%ref_re,      [self%pftsz])
            call c_f_pointer(self%fftdat(ithr)%p_ref_im,      self%fftdat(ithr)%ref_im,      [self%pftsz])
            call c_f_pointer(self%fftdat(ithr)%p_ref_fft_re,  self%fftdat(ithr)%ref_fft_re,  [self%pftsz])
            call c_f_pointer(self%fftdat(ithr)%p_ref_fft_im,  self%fftdat(ithr)%ref_fft_im,  [self%pftsz])
            call c_f_pointer(self%fftdat(ithr)%p_product_fft, self%fftdat(ithr)%product_fft, [self%nrots])
            call c_f_pointer(self%fftdat(ithr)%p_backtransf,  self%fftdat(ithr)%backtransf,  [self%nrots])
        end do
        ! thread-safe c-style allocatables for gencorrs, particle memoization
        do iptcl = self%pfromto(1),self%pfromto(2)
            do ik = self%kfromto(1),self%kfromto(2)
                self%fftdat_ptcls(iptcl,ik)%p_re = fftwf_alloc_complex(int(self%pftsz, c_size_t))
                self%fftdat_ptcls(iptcl,ik)%p_im = fftwf_alloc_complex(int(self%pftsz, c_size_t))
                call c_f_pointer(self%fftdat_ptcls(iptcl,ik)%p_re, self%fftdat_ptcls(iptcl,ik)%re, [self%pftsz])
                call c_f_pointer(self%fftdat_ptcls(iptcl,ik)%p_im, self%fftdat_ptcls(iptcl,ik)%im, [self%pftsz])
            end do
        end do
        ! thread-safe c-style allocatables for gencorrs, reference memoization
        do iref = 1,self%nrefs
            do ik = self%kfromto(1),self%kfromto(2)
                self%fftdat_refs_even(iref,ik)%p_re = fftwf_alloc_complex(int(self%pftsz, c_size_t))
                self%fftdat_refs_even(iref,ik)%p_im = fftwf_alloc_complex(int(self%pftsz, c_size_t))
                call c_f_pointer(self%fftdat_refs_even(iref,ik)%p_re, self%fftdat_refs_even(iref,ik)%re, [self%pftsz])
                call c_f_pointer(self%fftdat_refs_even(iref,ik)%p_im, self%fftdat_refs_even(iref,ik)%im, [self%pftsz])
                self%fftdat_refs_odd(iref,ik)%p_re = fftwf_alloc_complex(int(self%pftsz, c_size_t))
                self%fftdat_refs_odd(iref,ik)%p_im = fftwf_alloc_complex(int(self%pftsz, c_size_t))
                call c_f_pointer(self%fftdat_refs_odd(iref,ik)%p_re, self%fftdat_refs_odd(iref,ik)%re, [self%pftsz])
                call c_f_pointer(self%fftdat_refs_odd(iref,ik)%p_im, self%fftdat_refs_odd(iref,ik)%im, [self%pftsz])
            end do
        end do
        ! FFTW plans
        wsdm_ret = fftw_import_wisdom_from_filename(WISDOM_FNAME)
        self%plan_fwd_1 = fftwf_plan_dft_r2c_1d(self%pftsz, self%fftdat(1)%ref_re, &
             self%fftdat(1)%ref_fft_re, FFTW_PATIENT)
        self%plan_fwd_2 = fftwf_plan_dft_1d    (self%pftsz, self%fftdat(1)%ref_im, &
             self%fftdat(1)%ref_fft_im, FFTW_FORWARD, FFTW_PATIENT)
        self%plan_bwd   = fftwf_plan_dft_c2r_1d(self%nrots, self%fftdat(1)%product_fft, &
             self%fftdat(1)%backtransf, FFTW_PATIENT)
        wsdm_ret = fftw_export_wisdom_to_filename(WISDOM_FNAME)
        if (wsdm_ret == 0) then
           write (*, *) 'Error: could not write FFTW3 wisdom file! Check permissions.'
        end if
        ! factors for expansion of phase terms
        allocate(self%fft_factors(self%pftsz))
        do irot = 1,self%pftsz
            self%fft_factors(irot) = exp(-(0.,1.)*PI*real(irot-1)/real(self%pftsz))
        end do
        ! flag existence
        self%existence    = .true.
    end subroutine new

    ! SETTERS

    !>  \brief  sets reference pft iref
    subroutine set_ref_pft( self, iref, pft, iseven )
        class(polarft_corrcalc), intent(inout) :: self     !< this object
        integer,                 intent(in)    :: iref     !< reference index
        complex(sp),             intent(in)    :: pft(:,:) !< reference pft
        logical,                 intent(in)    :: iseven   !< logical eo-flag
        if( iseven )then
            self%pfts_refs_even(iref,:,:) = pft
        else
            self%pfts_refs_odd(iref,:,:)  = pft
        endif
    end subroutine set_ref_pft

    !>  \brief  sets particle pft iptcl
    subroutine set_ptcl_pft( self, iptcl, pft )
        class(polarft_corrcalc), intent(inout) :: self     !< this object
        integer,                 intent(in)    :: iptcl    !< particle index
        complex(sp),             intent(in)    :: pft(:,:) !< particle's pft
        self%pfts_ptcls(iptcl,:,:) = pft
        ! calculate the square sum required for correlation calculation
        call self%memoize_sqsum_ptcl(iptcl)
    end subroutine set_ptcl_pft

    !>  \brief set_ref_fcomp sets a reference Fourier component
    !! \param iref reference index
    !! \param irot rotation index
    !! \param k  index (third dim ptfs_refs)
    !! \param comp Fourier component
    subroutine set_ref_fcomp( self, iref, irot, k, comp, iseven )
        class(polarft_corrcalc), intent(inout) :: self
        integer,                 intent(in)    :: iref, irot, k
        complex(sp),             intent(in)    :: comp
        logical,                 intent(in)    :: iseven
        if( iseven )then
            self%pfts_refs_even(iref,irot,k) = comp
        else
            self%pfts_refs_odd(iref,irot,k)  = comp
        endif
    end subroutine set_ref_fcomp

    !>  \brief  sets a particle Fourier component
    !! \param iptcl particle index
    !! \param irot rotation index
    !! \param k  index (third dim ptfs_ptcls)
    !! \param comp Fourier component
    subroutine set_ptcl_fcomp( self, iptcl, irot, k, comp )
        class(polarft_corrcalc), intent(inout) :: self
        integer,                 intent(in)    :: iptcl, irot, k
        complex(sp),             intent(in)    :: comp
        self%pfts_ptcls(iptcl,irot,k) = comp
    end subroutine set_ptcl_fcomp

    !>  \brief  zeroes the iref reference
     !! \param iref reference index
    subroutine zero_ref( self, iref )
        class(polarft_corrcalc), intent(inout) :: self
        integer,                 intent(in)    :: iref
        self%pfts_refs_even(iref,:,:) = zero
        self%pfts_refs_odd(iref,:,:)  = zero
    end subroutine zero_ref

    subroutine cp_even2odd_ref( self, iref )
        class(polarft_corrcalc), intent(inout) :: self
        integer,                 intent(in)    :: iref
        self%pfts_refs_odd(iref,:,:) = self%pfts_refs_even(iref,:,:)
    end subroutine cp_even2odd_ref

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

    !>  \brief  for getting the number of reference rotations (size of second dim of self%pfts_refs_even)
    function get_pftsz( self ) result( pftsz )
        class(polarft_corrcalc), intent(in) :: self
        integer :: pftsz
        pftsz = self%pftsz
    end function get_pftsz

    !>  \brief  for getting the logical dimension of the original Cartesian image
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
    function get_pdim( self ) result( pdim )
        class(polarft_corrcalc), intent(in) :: self
        integer :: pdim(3)
        pdim = [self%pftsz,self%kfromto(1),self%kfromto(2)]
    end function get_pdim

    !>  \brief is for getting the continuous in-plane rotation
    !!         corresponding to in-plane rotation index roind
    function get_rot( self, roind ) result( rot )
        class(polarft_corrcalc), intent(in) :: self
        integer,                 intent(in) :: roind !<  in-plane rotation index
        real(sp) :: rot
        if( roind < 1 .or. roind > self%nrots )then
            print *, 'roind: ', roind
            print *, 'nrots: ', self%nrots
            stop 'roind is out of range; get_rot; simple_polarft_corrcalc'
        endif
        rot = self%angtab(roind)
    end function get_rot

    !>  \brief is for getting the (sign inversed) rotations for application in
    !!         classaverager/reconstructor
    function get_rots_for_applic( self ) result( inplrots )
        class(polarft_corrcalc), intent(in) :: self
        real, allocatable :: inplrots(:)
        allocate(inplrots(self%nrots), source=360. - self%angtab)
    end function get_rots_for_applic

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
        allocate(pft(self%pftsz,self%kfromto(1):self%kfromto(2)),&
        source=self%pfts_ptcls(iptcl,:,:), stat=alloc_stat)
        allocchk("In: get_ptcl_pft; simple_polarft_corrcalc")
    end function get_ptcl_pft

    !>  \brief  returns polar Fourier transform of reference iref
    !! \param iref reference index
    function get_ref_pft( self, iref, iseven ) result( pft )
        class(polarft_corrcalc), intent(in) :: self
        integer,                 intent(in) :: iref
        logical,                 intent(in) :: iseven
        complex(sp), allocatable :: pft(:,:)
        if( iseven )then
            allocate(pft(self%pftsz,self%kfromto(1):self%kfromto(2)),&
            source=self%pfts_refs_even(iref,:,:), stat=alloc_stat)
        else
            allocate(pft(self%pftsz,self%kfromto(1):self%kfromto(2)),&
            source=self%pfts_refs_odd(iref,:,:), stat=alloc_stat)
        endif
        allocchk("In: get_ref_pft; simple_polarft_corrcalc")
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
        call gnufor_image(real(self%pfts_ptcls(iptcl,:,:)),  palette='gray')
        call gnufor_image(aimag(self%pfts_ptcls(iptcl,:,:)), palette='gray')
    end subroutine vis_ptcl

    !>  \brief  is for plotting a particle polar FT
    !! \param iref reference index
    subroutine vis_ref( self, iref, iseven )
        use gnufor2
        class(polarft_corrcalc), intent(in) :: self
        integer,                 intent(in) :: iref
        logical,                 intent(in) :: iseven
        if( iseven )then
            call gnufor_image(real(self%pfts_refs_even(iref,:,:)),  palette='gray')
            call gnufor_image(aimag(self%pfts_refs_even(iref,:,:)), palette='gray')
        else
            call gnufor_image(real(self%pfts_refs_odd(iref,:,:)),  palette='gray')
            call gnufor_image(aimag(self%pfts_refs_odd(iref,:,:)), palette='gray')
        endif
    end subroutine vis_ref

    !>  \brief  for printing info about the object
    subroutine print( self )
        class(polarft_corrcalc), intent(in) :: self
        write(*,*) "from/to particle indices              (self%pfromto): ", self%pfromto
        write(*,*) "total n particles in partition         (self%nptcls): ", self%nptcls
        write(*,*) "number of references                    (self%nrefs): ", self%nrefs
        write(*,*) "number of rotations                     (self%nrots): ", self%nrots
        write(*,*) "radius of molecule                      (self%ring2): ", self%ring2
        write(*,*) "size of pft                             (self%pftsz): ", self%pftsz
        write(*,*) "logical dim. of original Cartesian image (self%ldim): ", self%ldim
        write(*,*) "high-pass limit Fourier index      (self%kfromto(1)): ", self%kfromto(1)
        write(*,*) "low-pass limit Fourier index       (self%kfromto(2)): ", self%kfromto(2)
    end subroutine print

    ! MEMOIZERS

    !>  \brief  is for memoization of the complex square sums required for correlation calculation
    !! \param iptcl particle index
    subroutine memoize_sqsum_ptcl( self, iptcl )
        class(polarft_corrcalc), intent(inout) :: self
        integer,                 intent(in)    :: iptcl
        self%sqsums_ptcls(iptcl) = sum(csq(self%pfts_ptcls(iptcl,:,:)))
    end subroutine memoize_sqsum_ptcl

    !>  \brief  is for memoization of the ffts used in gencorrs
    subroutine memoize_ffts( self )
        class(polarft_corrcalc), intent(inout) :: self
        type(fftw_carr) :: carray(self%nthr)
        integer         :: iref, iptcl, ik, ithr
        ! allocate local memory in a thread-safe manner
        do ithr = 1,self%nthr
            carray(ithr)%p_re = fftwf_alloc_real(int(self%pftsz, c_size_t))
            call c_f_pointer(carray(ithr)%p_re, carray(ithr)%re, [self%pftsz])
            carray(ithr)%p_im = fftwf_alloc_complex(int(self%pftsz, c_size_t))
            call c_f_pointer(carray(ithr)%p_im, carray(ithr)%im, [self%pftsz])
        end do
        ! memoize reference FFTs in parallel
        !$omp parallel default(shared) private(iref,iptcl,ik,ithr) proc_bind(close)
        !$omp do collapse(2) schedule(static)
        do iref = 1,self%nrefs
            do ik = self%kfromto(1),self%kfromto(2)
                ! get thread index
                ithr = omp_get_thread_num() + 1
                ! copy even reference pfts
                carray(ithr)%re = real(self%pfts_refs_even(iref,:,ik))
                carray(ithr)%im = aimag(self%pfts_refs_even(iref,:,ik)) * self%fft_factors
                ! FFT
                call fftwf_execute_dft_r2c(self%plan_fwd_1, carray(ithr)%re, self%fftdat_refs_even(iref,ik)%re)
                call fftwf_execute_dft    (self%plan_fwd_2, carray(ithr)%im, self%fftdat_refs_even(iref,ik)%im)
                ! copy odd reference pfts
                carray(ithr)%re = real(self%pfts_refs_odd(iref,:,ik))
                carray(ithr)%im = aimag(self%pfts_refs_odd(iref,:,ik)) * self%fft_factors
                ! FFT
                call fftwf_execute_dft_r2c(self%plan_fwd_1, carray(ithr)%re, self%fftdat_refs_odd(iref,ik)%re)
                call fftwf_execute_dft    (self%plan_fwd_2, carray(ithr)%im, self%fftdat_refs_odd(iref,ik)%im)
            end do
        end do
        !$omp end do nowait
        ! memoize particle FFTs in parallel
        !$omp do collapse(2) schedule(static)
        do iptcl = self%pfromto(1),self%pfromto(2)
            do ik = self%kfromto(1),self%kfromto(2)
                ! get thread index
                ithr = omp_get_thread_num() + 1
                ! copy particle pfts
                carray(ithr)%re = real(self%pfts_ptcls(iptcl,:,ik))
                carray(ithr)%im = aimag(self%pfts_ptcls(iptcl,:,ik)) * self%fft_factors
                ! FFT
                call fftwf_execute_dft_r2c(self%plan_fwd_1, carray(ithr)%re, self%fftdat_ptcls(iptcl,ik)%re)
                call fftwf_execute_dft    (self%plan_fwd_2, carray(ithr)%im, self%fftdat_ptcls(iptcl,ik)%im)
            end do
        end do
        !$omp end do
        !$omp end parallel
        ! free memory
        do ithr = 1,self%nthr
            call fftwf_free(carray(ithr)%p_re)
            call fftwf_free(carray(ithr)%p_im)
        end do
    end subroutine memoize_ffts

    ! CALCULATORS

    !>  \brief create_polar_ctfmat  is for generating a matrix of CTF values
    !! \param tfun transfer function object
    !! \param dfx,dfy resolution along Fourier axes
    !! \param angast astigmatic angle (degrees)
    !! \param add_phshift additional phase shift (radians) introduced by the Volta
    !! \param endrot number of rotations
    !! \return ctfmat matrix with CTF values
    function create_polar_ctfmat( self, tfun, dfx, dfy, angast, add_phshift, endrot ) result( ctfmat )
        !$ use omp_lib
        !$ use omp_lib_kinds
        use simple_ctf, only: ctf
        class(polarft_corrcalc), intent(inout) :: self
        class(ctf),              intent(inout) :: tfun
        real(sp),                intent(in)    :: dfx, dfy, angast, add_phshift
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
                if( self%phaseplate )then
                    ctfmat(irot,k) = tfun%eval(spaFreqSq,dfx,dfy,angast,ang,add_phshift)
                else
                    ctfmat(irot,k) = tfun%eval(spaFreqSq,dfx,dfy,angast,ang)
                endif
            end do
        end do
        !$omp end parallel do
    end function create_polar_ctfmat

    subroutine apply_ctf_to_ptcls( self, a )
        use simple_ctf,  only: ctf
        use simple_oris, only: oris
        class(polarft_corrcalc), intent(inout) :: self
        class(oris),             intent(inout) :: a
        type(ctf)             :: tfun
        integer               :: iptcl
        real(sp)              :: kv,cs,fraca,dfx,dfy,angast,phshift
        logical               :: astig
        real(sp), allocatable :: ctfmat(:,:)
        astig = a%isthere('dfy')
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
            phshift = 0.
            if( self%phaseplate ) phshift = a%get(iptcl, 'phshift')
            tfun   = ctf(self%smpd, kv, cs, fraca)
            ctfmat = self%create_polar_ctfmat(tfun, dfx, dfy, angast, phshift, self%pftsz)
            self%pfts_ptcls(iptcl,:,:) = self%pfts_ptcls(iptcl,:,:) * ctfmat
        end do
    end subroutine apply_ctf_to_ptcls

    subroutine calc_corrs_over_k( self, pft_ref, iptcl, kstop, corrs_over_k )
        class(polarft_corrcalc), intent(inout) :: self
        integer,                 intent(in)    :: iptcl, kstop
        complex(sp),             intent(in)    :: pft_ref(1:self%pftsz,self%kfromto(1):kstop)
        real,                    intent(out)   :: corrs_over_k(self%nrots)
        integer :: ithr, ik
        ! get thread index
        ithr = omp_get_thread_num() + 1
        ! sum up correlations over k-rings
        corrs_over_k = 0.
        do ik = self%kfromto(1),kstop
            ! move reference into Fourier Fourier space (particles are memoized)
            self%fftdat(ithr)%ref_re(:) =  real(pft_ref(:,ik))
            self%fftdat(ithr)%ref_im(:) = aimag(pft_ref(:,ik)) * self%fft_factors
            call fftwf_execute_dft_r2c(self%plan_fwd_1, self%fftdat(ithr)%ref_re, self%fftdat(ithr)%ref_fft_re)
            call fftwf_execute_dft    (self%plan_fwd_2, self%fftdat(ithr)%ref_im, self%fftdat(ithr)%ref_fft_im)
            ! correlate
            self%fftdat(ithr)%ref_fft_re = self%fftdat(ithr)%ref_fft_re * conjg(self%fftdat_ptcls(iptcl,ik)%re)
            self%fftdat(ithr)%ref_fft_im = self%fftdat(ithr)%ref_fft_im * conjg(self%fftdat_ptcls(iptcl,ik)%im)
            self%fftdat(ithr)%product_fft(1:1+2*int(self%pftsz/2):2) = &
                4. * self%fftdat(ithr)%ref_fft_re(1:1+int(self%pftsz/2))
            self%fftdat(ithr)%product_fft(2:2+2*int(self%pftsz/2):2) = &
                4. * self%fftdat(ithr)%ref_fft_im(1:int(self%pftsz/2)+1)
            ! back transform
            call fftwf_execute_dft_c2r(self%plan_bwd, self%fftdat(ithr)%product_fft, self%fftdat(ithr)%backtransf)
            ! accumulate corrs
            corrs_over_k = corrs_over_k + self%fftdat(ithr)%backtransf
        end do
        ! fftw3 routines are not properly normalized, hence division by self%nrots * 2
        corrs_over_k = corrs_over_k  / real(self%nrots * 2)
        ! corrs_over_k needs to be reordered
        corrs_over_k = corrs_over_k(self%nrots:1:-1) ! step 1 is reversing
        corrs_over_k = cshift(corrs_over_k, -1)      ! step 2 is circular shift by 1
    end subroutine calc_corrs_over_k

    subroutine calc_corrs_over_k_wrefmem( self, iref, iptcl, kstop, corrs_over_k )
        class(polarft_corrcalc), intent(inout) :: self
        integer,                 intent(in)    :: iref, iptcl, kstop
        real,                    intent(out)   :: corrs_over_k(self%nrots)
        integer :: ithr, ik
        ! get thread index
        ithr = omp_get_thread_num() + 1
        ! sum up correlations over k-rings
        corrs_over_k = 0.
        do ik = self%kfromto(1),kstop
            ! correlate FFTs
            if( self%iseven(iptcl) )then
                self%fftdat(ithr)%ref_fft_re = self%fftdat_refs_even(iref,ik)%re * conjg(self%fftdat_ptcls(iptcl,ik)%re)
                self%fftdat(ithr)%ref_fft_im = self%fftdat_refs_even(iref,ik)%im * conjg(self%fftdat_ptcls(iptcl,ik)%im)
            else
                self%fftdat(ithr)%ref_fft_re = self%fftdat_refs_odd(iref,ik)%re * conjg(self%fftdat_ptcls(iptcl,ik)%re)
                self%fftdat(ithr)%ref_fft_im = self%fftdat_refs_odd(iref,ik)%im * conjg(self%fftdat_ptcls(iptcl,ik)%im)
            endif
            self%fftdat(ithr)%product_fft(1:1+2*int(self%pftsz/2):2) = &
                4. * self%fftdat(ithr)%ref_fft_re(1:1+int(self%pftsz/2))
            self%fftdat(ithr)%product_fft(2:2+2*int(self%pftsz/2):2) = &
                4. * self%fftdat(ithr)%ref_fft_im(1:int(self%pftsz/2)+1)
            ! back transform
            call fftwf_execute_dft_c2r(self%plan_bwd, self%fftdat(ithr)%product_fft, self%fftdat(ithr)%backtransf)
            ! accumulate corrs
            corrs_over_k = corrs_over_k + self%fftdat(ithr)%backtransf
        end do
        ! fftw3 routines are not properly normalized, hence division by self%nrots * 2
        corrs_over_k = corrs_over_k  / real(self%nrots * 2)
        ! corrs_over_k needs to be reordered
        corrs_over_k = corrs_over_k(self%nrots:1:-1) ! step 1 is reversing
        corrs_over_k = cshift(corrs_over_k, -1)      ! step 2 is circular shift by 1
    end subroutine calc_corrs_over_k_wrefmem

    subroutine calc_k_corrs_wrefmem( self, iref, iptcl, k, kcorrs )
        class(polarft_corrcalc), intent(inout) :: self
        integer,                 intent(in)    :: iref, iptcl, k
        real,                    intent(out)   :: kcorrs(self%nrots)
        integer :: ithr
        ! get thread index
        ithr = omp_get_thread_num() + 1
        ! correlate FFTs
        if( self%iseven(iptcl) )then
            self%fftdat(ithr)%ref_fft_re = self%fftdat_refs_even(iref,k)%re * conjg(self%fftdat_ptcls(iptcl,k)%re)
            self%fftdat(ithr)%ref_fft_im = self%fftdat_refs_even(iref,k)%im * conjg(self%fftdat_ptcls(iptcl,k)%im)
        else
            self%fftdat(ithr)%ref_fft_re = self%fftdat_refs_odd(iref,k)%re * conjg(self%fftdat_ptcls(iptcl,k)%re)
            self%fftdat(ithr)%ref_fft_im = self%fftdat_refs_odd(iref,k)%im * conjg(self%fftdat_ptcls(iptcl,k)%im)
        endif
        self%fftdat(ithr)%product_fft(1:1+2*int(self%pftsz/2):2) = &
            4. * self%fftdat(ithr)%ref_fft_re(1:1+int(self%pftsz/2))
        self%fftdat(ithr)%product_fft(2:2+2*int(self%pftsz/2):2) = &
            4. * self%fftdat(ithr)%ref_fft_im(1:int(self%pftsz/2)+1)
        ! back transform
        call fftwf_execute_dft_c2r(self%plan_bwd, self%fftdat(ithr)%product_fft, self%fftdat(ithr)%backtransf)
        ! fftw3 routines are not properly normalized, hence division by self%nrots * 2
        kcorrs = self%fftdat(ithr)%backtransf / real(self%nrots * 2)
        ! kcorrs needs to be reordered
        kcorrs = kcorrs(self%nrots:1:-1) ! step 1 is reversing
        kcorrs = cshift(kcorrs, -1)      ! step 2 is circular shift by 1
    end subroutine calc_k_corrs_wrefmem

    function gencorrs_1( self, iref, iptcl ) result( cc )
        use simple_math, only: csq
        class(polarft_corrcalc), intent(inout) :: self
        integer,                 intent(in)    :: iref, iptcl
        real(sp) :: cc(self%nrots), sqsum_ref
        real     :: corrs_over_k(self%nrots)
        call self%calc_corrs_over_k_wrefmem(iref, iptcl, self%kfromto(2), corrs_over_k)
        if( self%iseven(iptcl) )then
            sqsum_ref = sum(csq(self%pfts_refs_even(iref,:,:)))
        else
            sqsum_ref = sum(csq(self%pfts_refs_odd(iref,:,:)))
        endif
        cc = corrs_over_k / sqrt(sqsum_ref * self%sqsums_ptcls(iptcl))
    end function gencorrs_1

    function gencorrs_2( self, iref, iptcl, kstop ) result( cc )
        use simple_math, only: csq
        class(polarft_corrcalc), intent(inout) :: self
        integer,                 intent(in)    :: iref, iptcl, kstop
        real(sp)    :: cc(self%nrots), sqsum_ref, sqsum_ptcl
        real        :: corrs_over_k(self%nrots)
        call self%calc_corrs_over_k_wrefmem(iref, iptcl, kstop, corrs_over_k)
        sqsum_ptcl = sum(csq(self%pfts_ptcls(iptcl, :, self%kfromto(1):kstop)))
        if( self%iseven(iptcl) )then
            sqsum_ref = sum(csq(self%pfts_refs_even(iref,:,self%kfromto(1):kstop)))
        else
            sqsum_ref = sum(csq(self%pfts_refs_odd(iref,:,self%kfromto(1):kstop)))
        endif
        cc = corrs_over_k / sqrt(sqsum_ref * sqsum_ptcl)
    end function gencorrs_2

    function gencorrs_3( self, iref, iptcl, shvec ) result( cc )
        use simple_math, only: csq
        class(polarft_corrcalc), intent(inout) :: self
        integer,                 intent(in)    :: iref, iptcl
        real(sp),                intent(in)    :: shvec(2)
        complex(sp) :: pft_ref(self%pftsz,self%kfromto(1):self%kfromto(2))
        complex(sp) :: shmat(self%pftsz,self%kfromto(1):self%kfromto(2))
        real(sp)    :: corrs_over_k(self%nrots), argmat(self%pftsz,self%kfromto(1):self%kfromto(2))
        real(sp)    :: cc(self%nrots), sqsum_ref
        ! generate the argument matrix from memoized components in argtransf
        argmat    = self%argtransf(:self%pftsz,:) * shvec(1) + self%argtransf(self%pftsz+1:,:) * shvec(2)
        ! generate the complex shift transformation matrix
        shmat     = cmplx(cos(argmat),sin(argmat))
        ! shift
        if( self%iseven(iptcl) )then
            pft_ref = self%pfts_refs_even(iref,:,:) * shmat
        else
            pft_ref = self%pfts_refs_odd(iref,:,:) * shmat
        endif
        sqsum_ref = sum(csq(pft_ref))
        ! correlate
        call self%calc_corrs_over_k(pft_ref, iptcl, self%kfromto(2), corrs_over_k)
        cc = corrs_over_k  / sqrt(sqsum_ref * self%sqsums_ptcls(iptcl))
    end function gencorrs_3

    !>  \brief  is for generating resolution dependent correlations
    function genfrc( self, iref, iptcl, irot ) result( frc )
        use simple_math, only: csq
        class(polarft_corrcalc), target, intent(inout) :: self
        integer,                         intent(in)    :: iref, iptcl, irot
        real(sp), allocatable :: frc(:)
        real(sp)              :: kcorrs(self%nrots), sumsqref, sumsqptcl
        integer               :: k
        allocate( frc(self%kfromto(1):self%kfromto(2)) )
        ! calc k-corrs and norms
        do k=self%kfromto(1),self%kfromto(2)
            call self%calc_k_corrs_wrefmem(iref, iptcl, k, kcorrs)
            sumsqptcl = sum(csq(self%pfts_ptcls(iptcl,:,k)))
            if( self%iseven(iptcl) )then
                sumsqref = sum(csq(self%pfts_refs_even(iref,:,k)))
            else
                sumsqref = sum(csq(self%pfts_refs_odd(iref,:,k)))
            endif
            frc(k) = kcorrs(irot) / sqrt(sumsqref * sumsqptcl)
        end do
    end function genfrc

    ! DESTRUCTOR

    !>  \brief  is a destructor
    subroutine kill( self )
        class(polarft_corrcalc), intent(inout) :: self
        integer :: ithr, iptcl, ik, iref
        if( self%existence )then
            do ithr=1,self%nthr
                call fftwf_free(self%fftdat(ithr)%p_ref_re)
                call fftwf_free(self%fftdat(ithr)%p_ref_im)
                call fftwf_free(self%fftdat(ithr)%p_ref_fft_re)
                call fftwf_free(self%fftdat(ithr)%p_ref_fft_im)
                call fftwf_free(self%fftdat(ithr)%p_product_fft)
                call fftwf_free(self%fftdat(ithr)%p_backtransf)
            end do
            do iptcl = self%pfromto(1),self%pfromto(2)
                do ik = self%kfromto(1),self%kfromto(2)
                    call fftwf_free(self%fftdat_ptcls(iptcl,ik)%p_re)
                    call fftwf_free(self%fftdat_ptcls(iptcl,ik)%p_im)
                end do
            end do
            do iref = 1,self%nrefs
                do ik = self%kfromto(1),self%kfromto(2)
                    call fftwf_free(self%fftdat_refs_even(iref,ik)%p_re)
                    call fftwf_free(self%fftdat_refs_even(iref,ik)%p_im)
                    call fftwf_free(self%fftdat_refs_odd(iref,ik)%p_re)
                    call fftwf_free(self%fftdat_refs_odd(iref,ik)%p_im)
                end do
            end do
            deallocate( self%sqsums_ptcls, self%angtab, self%argtransf,&
                &self%polar, self%pfts_refs_even, self%pfts_refs_odd, self%pfts_ptcls,&
                &self%fft_factors, self%fftdat, self%fftdat_ptcls, self%fftdat_refs_even,&
                &self%fftdat_refs_odd, self%iseven)
            call fftwf_destroy_plan(self%plan_bwd)
            call fftwf_destroy_plan(self%plan_fwd_1)
            call fftwf_destroy_plan(self%plan_fwd_2)
            self%existence = .false.
        endif
    end subroutine kill

end module simple_polarft_corrcalc
