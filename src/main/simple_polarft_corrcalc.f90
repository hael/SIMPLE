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
integer,     parameter :: FFTW_USE_WISDOM=16

! the fftw_arrs data structures are needed for thread-safe FFTW exec. Letting OpenMP copy out the per-threads
! arrays leads to bugs because of inconsistency between data in memory and the fftw_plan
type fftw_carr
    type(c_ptr)                            :: p_re                      !< pointer for C-style allocation
    type(c_ptr)                            :: p_im                      !< -"-
    real(kind=c_float),            pointer :: re(:) => null()           !< corresponding Fortran pointers
    complex(kind=c_float_complex), pointer :: im(:) => null()           !< -"-
end type fftw_carr

type fftw_carr_fft
    type(c_ptr)                            :: p_re                      !< pointer for C-style allocation
    type(c_ptr)                            :: p_im                      !< -"-bfac_norms
    complex(kind=c_float_complex), pointer :: re(:) => null()           !< corresponding Fortran pointers
    complex(kind=c_float_complex), pointer :: im(:) => null()           !< -"-
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

type heap_vars
    complex(sp), pointer :: pft_ref(:,:)        => null()
    complex(sp), pointer :: pft_ref_tmp(:,:)    => null()
    complex(sp), pointer :: pft_ref_tmp1(:,:)   => null()
    complex(sp), pointer :: pft_ref_tmp2(:,:)   => null()
    real,        pointer :: corrs_over_k(:)     => null()
    real(sp),    pointer :: argmat(:,:)         => null()
    complex(sp), pointer :: shmat(:,:)          => null()
    real,        pointer :: kcorrs(:)           => null()
    complex(dp), pointer :: pft_ref_8(:,:)      => null()
    complex(dp), pointer :: pft_ref_tmp_8(:,:)  => null()
    complex(dp), pointer :: pft_ref_tmp1_8(:,:) => null()
    complex(dp), pointer :: pft_ref_tmp2_8(:,:) => null()
    complex(dp), pointer :: shmat_8(:,:)        => null()
    real(dp),    pointer :: argmat_8(:,:)       => null()
end type heap_vars

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
    integer,             allocatable :: pinds(:)              !< index array (to reduce memory when frac_update < 1)
    real(sp),            allocatable :: ptcl_bfac_weights(:,:)!< B-factor per particle array for weighting of the correlation
    real(sp),            allocatable :: ptcl_bfac_norms(:)    !< normalisation constants for B-factor weighted ccres
    real(sp),            allocatable :: sqsums_ptcls(:)       !< memoized square sums for the correlation calculations
    real(sp),            allocatable :: angtab(:)             !< table of in-plane angles (in degrees)
    real(sp),            allocatable :: argtransf(:,:)        !< argument transfer constants for shifting the references
    real(sp),            allocatable :: polar(:,:)            !< table of polar coordinates (in Cartesian coordinates)
    real(sp),            allocatable :: ctfmats(:,:,:)        !< expand set of CTF matrices (for efficient parallel exec)
    complex(sp),         allocatable :: pfts_refs_even(:,:,:) !< 3D complex matrix of polar reference sections (nrefs,pftsz,nk), even
    complex(sp),         allocatable :: pfts_refs_odd(:,:,:)  !< -"-, odd
    complex(sp),         allocatable :: pfts_ptcls(:,:,:)     !< 3D complex matrix of particle sections
    complex(sp),         allocatable :: fft_factors(:)        !< phase factors for accelerated gencorrs routines
    type(fftw_arrs),     allocatable :: fftdat(:)             !< arrays for accelerated gencorrs routines
    type(fftw_carr_fft), allocatable :: fftdat_ptcls(:,:)     !< for memoization of particle  FFTs in accelerated gencorrs routines
    logical,             allocatable :: iseven(:)             !< eo assignment for gold-standard FSC
    logical                          :: phaseplate            !< images obtained with the Volta
    type(c_ptr)                      :: plan_fwd_1            !< FFTW plans for gencorrs
    type(c_ptr)                      :: plan_fwd_2            !< -"-
    type(c_ptr)                      :: plan_bwd              !< -"-
    logical                          :: with_ctf    = .false. !< CTF flag
    logical                          :: existence   = .false. !< to indicate existence
    logical                          :: l_cc_objfun = .true.  !< objective function(cc|ccres)
    logical                          :: l_cc_bfac   = .false. !< flag for B-factor weighting of the correlation
    type(heap_vars),     allocatable :: heap_vars(:)          !< allocated fields to save stack allocation in subroutines and functions
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
    procedure          :: cp_even_ref2ptcl
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
    procedure          :: objfun_is_ccres
    procedure          :: exists
    ! PRINTERS/VISUALISERS
    procedure          :: print
    procedure          :: vis_ptcl
    procedure          :: vis_ref
    ! MEMOIZER
    procedure, private :: memoize_sqsum_ptcl
    procedure          :: memoize_ffts
    procedure          :: memoize_bfac
    procedure          :: memoize_bfacs
    ! CALCULATORS
    procedure, private :: create_polar_ctfmat
    procedure          :: create_polar_ctfmats
    procedure, private :: prep_ref4corr
    procedure, private :: calc_corrs_over_k
    procedure, private :: calc_k_corrs
    procedure, private :: calc_corr_for_rot
    procedure, private :: calc_corr_for_rot_8
    procedure, private :: calc_corrk_for_rot
    procedure, private :: calc_corrk_for_rot_8
    procedure, private :: gencorrs_cc_1
    procedure, private :: gencorrs_cc_2
    procedure, private :: gencorrs_cc_3
    procedure, private :: gencorrs_resnorm_1
    procedure, private :: gencorrs_resnorm_2
    procedure, private :: gencorrs_resnorm_3
    procedure, private :: gencorrs_1
    procedure, private :: gencorrs_2
    procedure, private :: gencorrs_3
    generic            :: gencorrs => gencorrs_1, gencorrs_2, gencorrs_3
    procedure          :: gencorrs_cc_grad
    procedure          :: gencorrs_cc_grad_only
    procedure          :: gencorr_for_rot
    procedure          :: gencorr_for_rot_8
    procedure          :: gencorr_grad_for_rot
    procedure          :: gencorr_grad_only_for_rot
    procedure          :: gencorr_grad_for_rot_8
    procedure          :: gencorr_grad_only_for_rot_8
    procedure          :: gencorr_cc_for_rot
    procedure          :: gencorr_cc_for_rot_8
    procedure          :: gencorr_cc_grad_for_rot
    procedure          :: gencorr_cc_grad_only_for_rot
    procedure          :: gencorr_cc_grad_for_rot_8
    procedure          :: gencorr_cc_grad_only_for_rot_8
    procedure          :: gencorr_resnorm_for_rot
    procedure          :: gencorr_resnorm_for_rot_8
    procedure          :: gencorr_resnorm_grad_for_rot
    procedure          :: gencorr_resnorm_grad_only_for_rot
    procedure          :: gencorr_resnorm_grad_for_rot_8
    procedure          :: gencorr_resnorm_grad_only_for_rot_8
    procedure, private :: genfrc
    procedure          :: calc_frc
    procedure          :: specscore
    procedure          :: fit_bfac
    ! DESTRUCTOR
    procedure          :: kill
end type polarft_corrcalc

contains

    ! CONSTRUCTORS

    !>  \brief  is a constructor
    subroutine new( self, nrefs, p, ptcl_mask, eoarr )
        use simple_math,   only: rad2deg, is_even, round2even
        use simple_params, only: params
        class(polarft_corrcalc), intent(inout) :: self
        integer,                 intent(in)    :: nrefs
        class(params),           intent(in)    :: p
        logical, optional,       intent(in)    :: ptcl_mask(p%fromp:p%top)
        integer, optional,       intent(in)    :: eoarr(p%fromp:p%top)
        character(kind=c_char, len=:), allocatable :: fft_wisdoms_fname ! FFTW wisdoms (per part or suffer I/O lag)
        integer             :: alloc_stat, irot, k, ithr, i, ik, cnt
        logical             :: even_dims, test(2)
        real(sp)            :: ang, res
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
        self%pfromto     = [p%fromp,p%top]                       !< from/to particle indices (in parallel execution)
        if( present(ptcl_mask) )then
            self%nptcls  = count(ptcl_mask)                      !< the total number of particles in partition
        else
            self%nptcls  = p%top - p%fromp + 1                   !< the total number of particles in partition
        endif
        self%nrefs       = nrefs                                 !< the number of references (logically indexded [1,nrefs])
        self%ring2       = p%ring2                               !< radius of molecule
        self%nrots       = round2even(twopi * real(p%ring2))     !< number of in-plane rotations for one pft  (determined by radius of molecule)
        self%pftsz       = self%nrots / 2                        !< size of reference (nrots/2) (number of vectors used for matching)
        self%smpd        = p%smpd                                !< sampling distance
        self%kfromto     = p%kfromto                             !< Fourier index range
        self%nk          = self%kfromto(2) - self%kfromto(1) + 1 !< # resolution elements
        self%nthr        = p%nthr                                !< # OpenMP threads
        self%phaseplate  = p%tfplan%l_phaseplate                 !< images obtained with the Volta
        ! self%l_cc_bfac   = p%l_cc_bfac                           !< B-factor weighting of the correlation
        ! take care of objective function flags
        select case(trim(p%objfun))
            case('cc')
                self%l_cc_objfun   = .true.
            case('ccres')
                self%l_cc_objfun = .false.
                self%l_cc_bfac   = .true.
                allocate(self%ptcl_bfac_weights(1:self%nptcls, self%kfromto(1):self%kfromto(2)), self%ptcl_bfac_norms(1:self%nptcls))
                self%ptcl_bfac_weights = 1.0
                self%ptcl_bfac_norms   = real(self%nk)
            case DEFAULT
                write(*,*) 'unsupported objective function: ', trim(p%objfun)
                stop 'ABORTING, simple_polarft_corrcalc :: new'
        end select
        ! generate polar coordinates & eo assignment
        allocate( self%polar(2*self%nrots,self%kfromto(1):self%kfromto(2)),&
                 &self%angtab(self%nrots), self%iseven(1:self%nptcls), stat=alloc_stat)
        allocchk('polar coordinate arrays; new; simple_polarft_corrcalc, 1')
        ang = twopi/real(self%nrots)
        do irot=1,self%nrots
            self%angtab(irot) = real(irot-1)*ang
            do k=self%kfromto(1),self%kfromto(2)
                self%polar(irot,k)            = cos(self%angtab(irot))*real(k) ! x-coordinate
                self%polar(irot+self%nrots,k) = sin(self%angtab(irot))*real(k) ! y-coordinate
            end do
            self%angtab(irot) = rad2deg(self%angtab(irot)) ! angle (in degrees)
        end do
        ! index translation table
        allocate( self%pinds(p%fromp:p%top), source=0, stat=alloc_stat)
        allocchk('polar coordinate arrays; new; simple_polarft_corrcalc, 2')
        if( present(ptcl_mask) )then
            cnt = 0
            do i=p%fromp,p%top
                if( ptcl_mask(i) )then
                    cnt = cnt + 1
                    self%pinds(i) = cnt
                endif
            end do
        else
            self%pinds = (/(i,i=1,self%nptcls)/)
        endif
        ! eo assignment
        if( present(eoarr) )then
            if( all(eoarr == - 1) )then
                self%iseven = .true.
            else
                do i=p%fromp,p%top
                    if( self%pinds(i) > 0 )then
                        if( eoarr(i) == 0 )then
                            self%iseven(self%pinds(i)) = .true.
                        else
                            self%iseven(self%pinds(i)) = .false.
                        endif
                    endif
                end do
            endif
        else
            self%iseven = .true.
        endif
        ! generate the argument transfer constants for shifting reference polarfts
        allocate( self%argtransf(self%nrots,self%kfromto(1):self%kfromto(2)), stat=alloc_stat)
        allocchk('shift argument transfer array; new; simple_polarft_corrcalc')
        self%argtransf(:self%pftsz,:)   = &
            self%polar(:self%pftsz,:)   * &
            (PI/real(self%ldim(1) / 2))    ! x-part
        self%argtransf(self%pftsz + 1:,:) = &
            self%polar(self%nrots + 1:self%nrots+self%pftsz,:) * &
            (PI/real(self%ldim(2) / 2))    ! y-part
        ! allocate others
        allocate(self%pfts_refs_even(self%nrefs,self%pftsz,self%kfromto(1):self%kfromto(2)),&
                 &self%pfts_refs_odd(self%nrefs,self%pftsz,self%kfromto(1):self%kfromto(2)),&
                 &self%pfts_ptcls(1:self%nptcls,self%pftsz,self%kfromto(1):self%kfromto(2)),&
                 &self%sqsums_ptcls(1:self%nptcls),self%fftdat(self%nthr),&
                 &self%fftdat_ptcls(1:self%nptcls,self%kfromto(1):self%kfromto(2)),&
                 &self%heap_vars(self%nthr),stat=alloc_stat)
        do ithr=1,self%nthr
            allocate(self%heap_vars(ithr)%pft_ref(self%pftsz,self%kfromto(1):self%kfromto(2)),&
                &self%heap_vars(ithr)%pft_ref_tmp(self%pftsz,self%kfromto(1):self%kfromto(2)),&
                &self%heap_vars(ithr)%pft_ref_tmp1(self%pftsz,self%kfromto(1):self%kfromto(2)),&
                &self%heap_vars(ithr)%pft_ref_tmp2(self%pftsz,self%kfromto(1):self%kfromto(2)),&
                &self%heap_vars(ithr)%corrs_over_k(self%nrots),&
                &self%heap_vars(ithr)%argmat(self%pftsz,self%kfromto(1):self%kfromto(2)),&
                &self%heap_vars(ithr)%shmat(self%pftsz,self%kfromto(1):self%kfromto(2)),&
                &self%heap_vars(ithr)%kcorrs(self%nrots),&
                &self%heap_vars(ithr)%pft_ref_8(self%pftsz,self%kfromto(1):self%kfromto(2)),&
                &self%heap_vars(ithr)%pft_ref_tmp_8(self%pftsz,self%kfromto(1):self%kfromto(2)),&
                &self%heap_vars(ithr)%pft_ref_tmp1_8(self%pftsz,self%kfromto(1):self%kfromto(2)),&
                &self%heap_vars(ithr)%pft_ref_tmp2_8(self%pftsz,self%kfromto(1):self%kfromto(2)),&
                &self%heap_vars(ithr)%shmat_8(self%pftsz,self%kfromto(1):self%kfromto(2)),&
                &self%heap_vars(ithr)%argmat_8(self%pftsz,self%kfromto(1):self%kfromto(2)),&
                &stat=alloc_stat)
        end do
        allocchk('polarfts and sqsums; new; simple_polarft_corrcalc')
        self%pfts_refs_even = zero
        self%pfts_refs_odd  = zero
        self%pfts_ptcls     = zero
        self%sqsums_ptcls   = 0.
        ! set CTF flag
        self%with_ctf = .false.
        if( p%ctf .ne. 'no' ) self%with_ctf = .true.
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
        do i = 1,self%nptcls
            do ik = self%kfromto(1),self%kfromto(2)
                self%fftdat_ptcls(i,ik)%p_re = fftwf_alloc_complex(int(self%pftsz, c_size_t))
                self%fftdat_ptcls(i,ik)%p_im = fftwf_alloc_complex(int(self%pftsz, c_size_t))
                call c_f_pointer(self%fftdat_ptcls(i,ik)%p_re, self%fftdat_ptcls(i,ik)%re, [self%pftsz])
                call c_f_pointer(self%fftdat_ptcls(i,ik)%p_im, self%fftdat_ptcls(i,ik)%im, [self%pftsz])
            end do
        end do
        ! FFTW3 wisdoms file
        if( p%l_distr_exec )then
            allocate(fft_wisdoms_fname, source='fft_wisdoms_part'//int2str_pad(p%part,p%numlen)//'.dat'//c_null_char)
        else
            allocate(fft_wisdoms_fname, source='fft_wisdoms.dat'//c_null_char)
        endif
        ! FFTW plans
        wsdm_ret = fftw_import_wisdom_from_filename(fft_wisdoms_fname)
        self%plan_fwd_1 = fftwf_plan_dft_r2c_1d(self%pftsz, self%fftdat(1)%ref_re, &
             self%fftdat(1)%ref_fft_re, ior(FFTW_PATIENT, FFTW_USE_WISDOM))
        self%plan_fwd_2 = fftwf_plan_dft_1d(self%pftsz, self%fftdat(1)%ref_im, &
             self%fftdat(1)%ref_fft_im, FFTW_FORWARD, ior(FFTW_PATIENT, FFTW_USE_WISDOM))
        self%plan_bwd   = fftwf_plan_dft_c2r_1d(self%nrots, self%fftdat(1)%product_fft, &
             self%fftdat(1)%backtransf, ior(FFTW_PATIENT, FFTW_USE_WISDOM))
        wsdm_ret = fftw_export_wisdom_to_filename(fft_wisdoms_fname)
        deallocate(fft_wisdoms_fname)
        if (wsdm_ret == 0) then
           write (*, *) 'Error: could not write FFTW3 wisdom file! Check permissions.'
        end if
        ! factors for expansion of phase terms
        allocate(self%fft_factors(self%pftsz))
        do irot = 1,self%pftsz
            self%fft_factors(irot) = exp(-(0.,1.) * PI * real(irot - 1) / real(self%pftsz))
        end do
        ! flag existence
        self%existence = .true.
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
        self%pfts_ptcls(self%pinds(iptcl),:,:) = pft
        ! calculate the square sum required for correlation calculation
        call self%memoize_sqsum_ptcl(self%pinds(iptcl))
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
        self%pfts_ptcls(self%pinds(iptcl),irot,k) = comp
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

    subroutine cp_even_ref2ptcl( self, iref, iptcl )
        class(polarft_corrcalc), intent(inout) :: self
        integer,                 intent(in)    :: iref, iptcl
        self%pfts_ptcls(self%pinds(iptcl),:,:) = self%pfts_refs_even(iref,:,:)
        call self%memoize_sqsum_ptcl(self%pinds(iptcl))
    end subroutine cp_even_ref2ptcl

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
        integer,                 intent(in) :: roind !< in-plane rotation index
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

    !>  \brief  returns polar Fourier transform of particle iptcl
    !! \param iptcl particle index
    function get_ptcl_pft( self, iptcl) result( pft )
        class(polarft_corrcalc), intent(in) :: self
        integer,                 intent(in) :: iptcl
        complex(sp), allocatable :: pft(:,:)
        allocate(pft(self%pftsz,self%kfromto(1):self%kfromto(2)),&
        source=self%pfts_ptcls(self%pinds(iptcl),:,:), stat=alloc_stat)
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

    function objfun_is_ccres( self ) result( is )
        class(polarft_corrcalc), intent(in) :: self
        logical :: is
        is = .not. self%l_cc_objfun
    end function objfun_is_ccres

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
        call gnufor_image(real(self%pfts_ptcls(self%pinds(iptcl),:,:)),  palette='gray')
        call gnufor_image(aimag(self%pfts_ptcls(self%pinds(iptcl),:,:)), palette='gray')
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
    subroutine memoize_sqsum_ptcl( self, i )
        class(polarft_corrcalc), intent(inout) :: self
        integer,                 intent(in)    :: i
        self%sqsums_ptcls(i) = sum(csq(self%pfts_ptcls(i,:,:)))
    end subroutine memoize_sqsum_ptcl

    !>  \brief  is for memoization of the ffts used in gencorrs
    subroutine memoize_ffts( self )
        class(polarft_corrcalc), intent(inout) :: self
        type(fftw_carr) :: carray(self%nthr)
        integer         :: i, ik, ithr
        ! allocate local memory in a thread-safe manner
        do ithr = 1,self%nthr
            carray(ithr)%p_re = fftwf_alloc_real(int(self%pftsz, c_size_t))
            call c_f_pointer(carray(ithr)%p_re, carray(ithr)%re, [self%pftsz])
            carray(ithr)%p_im = fftwf_alloc_complex(int(self%pftsz, c_size_t))
            call c_f_pointer(carray(ithr)%p_im, carray(ithr)%im, [self%pftsz])
        end do
        ! memoize particle FFTs in parallel
        !$omp parallel do default(shared) private(i,ik,ithr) proc_bind(close) collapse(2) schedule(static)
        do i = 1, self%nptcls
            do ik = self%kfromto(1),self%kfromto(2)
                ! get thread index
                ithr = omp_get_thread_num() + 1
                ! copy particle pfts
                carray(ithr)%re = real(self%pfts_ptcls(i,:,ik))
                carray(ithr)%im = aimag(self%pfts_ptcls(i,:,ik)) * self%fft_factors
                ! FFT
                call fftwf_execute_dft_r2c(self%plan_fwd_1, carray(ithr)%re, self%fftdat_ptcls(i,ik)%re)
                call fftwf_execute_dft    (self%plan_fwd_2, carray(ithr)%im, self%fftdat_ptcls(i,ik)%im)
            end do
        end do
        !$omp end parallel do
        ! free memory
        do ithr = 1,self%nthr
            call fftwf_free(carray(ithr)%p_re)
            call fftwf_free(carray(ithr)%p_im)
        end do
    end subroutine memoize_ffts

    subroutine memoize_bfac( self, iptcl, bfac )
        class(polarft_corrcalc), intent(inout) :: self
        integer,                 intent(in)    :: iptcl
        real,                    intent(in)    :: bfac
        integer :: k, i
        real    :: res
        if( self%l_cc_objfun )then
            ! nothing to do
        else
            i = self%pinds(iptcl)
            do k=self%kfromto(1),self%kfromto(2)
                res = real(k) / (real(self%ldim(1)) * self%smpd) ! assuming square dimensions
                self%ptcl_bfac_weights(i,k) = max(0., exp(-(bfac/4.) * res * res))
            end do
            self%ptcl_bfac_norms(i) = sum(self%ptcl_bfac_weights(i,:))
        endif
    end subroutine memoize_bfac

    !>  \brief  is for memoization of the ffts used in gencorrs
    subroutine memoize_bfacs( self, a )
        use simple_oris, only: oris
        class(polarft_corrcalc), intent(inout) :: self
        class(oris),             intent(inout) :: a
        integer  :: iptcl, k
        real(sp) :: resarrsq(self%kfromto(1):self%kfromto(2)), bfac, tmp(self%kfromto(1):self%kfromto(2))
        if( self%l_cc_objfun )then
            !nothing to do
        else
            ! pre-calc squared spatial frequencies
            do k=self%kfromto(1),self%kfromto(2)
                resarrsq(k) = real(k) / (real(self%ldim(1)) * self%smpd) ! assuming square dimensions
                resarrsq(k) = resarrsq(k) * resarrsq(k)
            end do
            ! pre-calc B-factor weight arrays
            !$omp parallel do default(shared) private(iptcl,bfac,tmp) proc_bind(close) schedule(static)
            do iptcl=self%pfromto(1),self%pfromto(2)
                if( self%pinds(iptcl) > 0 )then
                    bfac = a%get(iptcl, 'bfac')
                    tmp  = exp(-(bfac/4.) * resarrsq(:))
                    where( tmp > 0. )
                        self%ptcl_bfac_weights(self%pinds(iptcl),:) = tmp
                    else where
                        self%ptcl_bfac_weights(self%pinds(iptcl),:) = 0.
                    end where
                    self%ptcl_bfac_norms(self%pinds(iptcl)) = sum(self%ptcl_bfac_weights(self%pinds(iptcl),:))
                endif
            end do
            !$omp end parallel do
        endif
    end subroutine memoize_bfacs

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

    !>  \brief  is for generating all matrices of CTF values
    subroutine create_polar_ctfmats( self, a )
        use simple_ctf,  only: ctf
        use simple_oris, only: oris
        class(polarft_corrcalc), intent(inout) :: self
        class(oris),             intent(inout) :: a
        type(ctf) :: tfun
        integer   :: iptcl
        real(sp)  :: kv,cs,fraca,dfx,dfy,angast,phshift
        logical   :: astig
        astig = a%isthere('dfy')
        if( allocated(self%ctfmats) ) deallocate(self%ctfmats)
        allocate(self%ctfmats(1:self%nptcls,self%pftsz,self%kfromto(1):self%kfromto(2)), stat=alloc_stat)
        allocchk("In: simple_polarft_corrcalc :: create_polar_ctfmats, 2")
        do iptcl=self%pfromto(1),self%pfromto(2)
            if( self%pinds(iptcl) > 0 )then
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
                tfun = ctf(self%smpd, kv, cs, fraca)
                self%ctfmats(self%pinds(iptcl),:,:) = self%create_polar_ctfmat(tfun, dfx, dfy, angast, phshift, self%pftsz)
            endif
        end do
    end subroutine create_polar_ctfmats

    subroutine prep_ref4corr( self, iref, i, pft_ref, sqsum_ref, kstop )
         use simple_estimate_ssnr, only: fsc2optlp
        class(polarft_corrcalc), intent(inout) :: self
        integer,                 intent(in)    :: iref, i
        complex(sp),             intent(out)   :: pft_ref(self%pftsz,self%kfromto(1):kstop)
        real(sp),                intent(out)   :: sqsum_ref
        integer,                 intent(in)    :: kstop
        ! copy
        if( self%iseven(i) )then
            pft_ref = self%pfts_refs_even(iref,:,:)
        else
            pft_ref = self%pfts_refs_odd(iref,:,:)
        endif
        ! multiply with CTF
        if( self%with_ctf ) pft_ref = pft_ref * self%ctfmats(i,:,:)
        ! for corr normalisation
        sqsum_ref = sum(csq(pft_ref(:,self%kfromto(1):kstop)))
    end subroutine prep_ref4corr

    subroutine calc_corrs_over_k( self, pft_ref, i, kstop, corrs_over_k )
        class(polarft_corrcalc), intent(inout) :: self
        integer,                 intent(in)    :: i, kstop
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
            self%fftdat(ithr)%ref_fft_re = self%fftdat(ithr)%ref_fft_re * conjg(self%fftdat_ptcls(i,ik)%re)
            self%fftdat(ithr)%ref_fft_im = self%fftdat(ithr)%ref_fft_im * conjg(self%fftdat_ptcls(i,ik)%im)
            self%fftdat(ithr)%product_fft(1:1 + 2 * int(self%pftsz / 2):2) = &
                4. * self%fftdat(ithr)%ref_fft_re(1:1+int(self%pftsz/2))
            self%fftdat(ithr)%product_fft(2:2 + 2 * int(self%pftsz / 2):2) = &
                4. * self%fftdat(ithr)%ref_fft_im(1:int(self%pftsz / 2) + 1)
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

    subroutine calc_k_corrs( self, pft_ref, i, k, kcorrs )
        class(polarft_corrcalc), intent(inout) :: self
        complex(sp),             intent(in)    :: pft_ref(1:self%pftsz,self%kfromto(1):self%kfromto(2))
        integer,                 intent(in)    :: i, k
        real,                    intent(out)   :: kcorrs(self%nrots)
        integer :: ithr
        ! get thread index
        ithr = omp_get_thread_num() + 1
        ! move reference into Fourier Fourier space (particles are memoized)
        self%fftdat(ithr)%ref_re(:) =  real(pft_ref(:,k))
        self%fftdat(ithr)%ref_im(:) = aimag(pft_ref(:,k)) * self%fft_factors
        call fftwf_execute_dft_r2c(self%plan_fwd_1, self%fftdat(ithr)%ref_re, self%fftdat(ithr)%ref_fft_re)
        call fftwf_execute_dft    (self%plan_fwd_2, self%fftdat(ithr)%ref_im, self%fftdat(ithr)%ref_fft_im)
        ! correlate FFTs
        self%fftdat(ithr)%ref_fft_re = self%fftdat(ithr)%ref_fft_re * conjg(self%fftdat_ptcls(i,k)%re)
        self%fftdat(ithr)%ref_fft_im = self%fftdat(ithr)%ref_fft_im * conjg(self%fftdat_ptcls(i,k)%im)
        self%fftdat(ithr)%product_fft(1:1 + 2 * int(self%pftsz / 2):2) = &
            4. * self%fftdat(ithr)%ref_fft_re(1:1 + int(self%pftsz / 2))
        self%fftdat(ithr)%product_fft(2:2 + 2 * int(self%pftsz / 2):2) = &
            4. * self%fftdat(ithr)%ref_fft_im(1:int(self%pftsz / 2) + 1)
        ! back transform
        call fftwf_execute_dft_c2r(self%plan_bwd, self%fftdat(ithr)%product_fft, self%fftdat(ithr)%backtransf)
        ! fftw3 routines are not properly normalized, hence division by self%nrots * 2
        kcorrs = self%fftdat(ithr)%backtransf / real(self%nrots * 2)
        ! kcorrs needs to be reordered
        kcorrs = kcorrs(self%nrots:1:-1) ! step 1 is reversing
        kcorrs = cshift(kcorrs, -1)      ! step 2 is circular shift by 1
    end subroutine calc_k_corrs

    function calc_corr_for_rot(self, pft_ref, i, kstop, irot) result(corr)
        class(polarft_corrcalc), intent(inout) :: self
        integer,                 intent(in)    :: i, kstop
        complex(sp),             intent(in)    :: pft_ref(1:self%pftsz,self%kfromto(1):kstop)
        integer,                 intent(in)    :: irot
        integer :: rot
        real    :: corr
        complex :: tmp
        corr = 0.
        tmp  = 0.
        if( irot >= self%pftsz + 1 )then
            rot = irot - self%pftsz
        else
            rot = irot
        end if
        if( irot == 1 )then
            tmp = sum( pft_ref(:,:) * conjg(self%pfts_ptcls(i,:,:)))
        else if( irot <= self%pftsz )then
            tmp =       sum( pft_ref(               1:self%pftsz-rot+1,:) * conjg(self%pfts_ptcls(i,rot:self%pftsz,:)))
            tmp = tmp + sum( pft_ref(self%pftsz-rot+2:self%pftsz,      :) *       self%pfts_ptcls(i,  1:rot-1,     :))
        else if( irot == self%pftsz + 1 )then
            tmp = sum( pft_ref(:,:) * self%pfts_ptcls(i,:,:) )
        else
            tmp =       sum( pft_ref(1:self%pftsz-rot+1,:)          *        self%pfts_ptcls(i,rot:self%pftsz,:))
            tmp = tmp + sum( pft_ref(self%pftsz-rot+2:self%pftsz,:) * conjg( self%pfts_ptcls(i,1:rot-1,:) ))
        end if
        corr = real(tmp)
    end function calc_corr_for_rot

    function calc_corr_for_rot_8(self, pft_ref, i, kstop, irot) result(corr)
        class(polarft_corrcalc), intent(inout) :: self
        integer,                 intent(in)    :: i, kstop
        complex(dp),             intent(in)    :: pft_ref(1:self%pftsz,self%kfromto(1):kstop)
        integer,                 intent(in)    :: irot
        integer     :: rot
        real(dp)    :: corr
        complex(dp) :: tmp
        corr = 0.
        tmp = 0.
        if (irot >= self%pftsz + 1) then
            rot = irot - self%pftsz
        else
            rot = irot
        end if
        if (irot == 1) then
            tmp = sum( pft_ref(:,:) * conjg(self%pfts_ptcls(i,:,:)))
        else if (irot <= self%pftsz) then
            tmp =       sum( pft_ref(               1:self%pftsz-rot+1,:) * conjg(self%pfts_ptcls(i,rot:self%pftsz,:)))
            tmp = tmp + sum( pft_ref(self%pftsz-rot+2:self%pftsz,      :) *       self%pfts_ptcls(i,  1:rot-1,     :))
        else if (irot == self%pftsz + 1) then
            tmp = sum( pft_ref(:,:) * self%pfts_ptcls(i,:,:) )
        else
            tmp =       sum( pft_ref(1:self%pftsz-rot+1,:)          *        self%pfts_ptcls(i,rot:self%pftsz,:))
            tmp = tmp + sum( pft_ref(self%pftsz-rot+2:self%pftsz,:) * conjg( self%pfts_ptcls(i,1:rot-1,:) ))
        end if
        corr = real(tmp)
    end function calc_corr_for_rot_8

    function calc_corrk_for_rot(self, pft_ref, i, kstop, k, irot) result(corr)
        class(polarft_corrcalc), intent(inout) :: self
        integer,                 intent(in)    :: i, kstop, k
        complex(sp),             intent(in)    :: pft_ref(1:self%pftsz,self%kfromto(1):kstop)
        integer,                 intent(in)    :: irot
        integer :: rot
        real    :: corr
        complex :: tmp
        corr = 0.
        tmp  = 0.
        if( irot >= self%pftsz + 1 )then
            rot = irot - self%pftsz
        else
            rot = irot
        end if
        if( irot == 1 )then
            tmp = sum( pft_ref(:,k) * conjg(self%pfts_ptcls(i,:,k)))
        else if( irot <= self%pftsz )then
            tmp =       sum( pft_ref(               1:self%pftsz-rot+1,k) * conjg(self%pfts_ptcls(i,rot:self%pftsz,k)))
            tmp = tmp + sum( pft_ref(self%pftsz-rot+2:self%pftsz,      k) *       self%pfts_ptcls(i,  1:rot-1,     k))
        else if( irot == self%pftsz + 1 )then
            tmp = sum( pft_ref(:,k) * self%pfts_ptcls(i,:,k) )
        else
            tmp =       sum( pft_ref(1:self%pftsz-rot+1,k)          *        self%pfts_ptcls(i,rot:self%pftsz,k))
            tmp = tmp + sum( pft_ref(self%pftsz-rot+2:self%pftsz,k) * conjg( self%pfts_ptcls(i,1:rot-1,k) ))
        end if
        corr = real(tmp)
    end function calc_corrk_for_rot

    function calc_corrk_for_rot_8(self, pft_ref, i, kstop, k, irot) result(corr)
        class(polarft_corrcalc), intent(inout) :: self
        integer,                 intent(in)    :: i, kstop, k
        complex(dp),             intent(in)    :: pft_ref(1:self%pftsz,self%kfromto(1):kstop)
        integer,                 intent(in)    :: irot
        integer     :: rot
        real(dp)    :: corr
        complex(dp) :: tmp
        corr = 0.
        tmp = 0.
        if (irot >= self%pftsz + 1) then
            rot = irot - self%pftsz
        else
            rot = irot
        end if
        if (irot == 1) then
            tmp = sum( pft_ref(:,k) * conjg(self%pfts_ptcls(i,:,k)))
        else if (irot <= self%pftsz) then
            tmp =       sum( pft_ref(               1:self%pftsz-rot+1,k) * conjg(self%pfts_ptcls(i,rot:self%pftsz,k)))
            tmp = tmp + sum( pft_ref(self%pftsz-rot+2:self%pftsz,      k) *       self%pfts_ptcls(i,  1:rot-1,     k))
        else if (irot == self%pftsz + 1) then
            tmp = sum( pft_ref(:,k) * self%pfts_ptcls(i,:,k) )
        else
            tmp =       sum( pft_ref(1:self%pftsz-rot+1,k)          *        self%pfts_ptcls(i,rot:self%pftsz,k))
            tmp = tmp + sum( pft_ref(self%pftsz-rot+2:self%pftsz,k) * conjg( self%pfts_ptcls(i,1:rot-1,k) ))
        end if
        corr = real(tmp)
    end function calc_corrk_for_rot_8

    !>  \brief  is for generating resolution dependent correlations
    subroutine genfrc( self, iref, i, irot, frc )
        use simple_math, only: csq
        class(polarft_corrcalc),  intent(inout) :: self
        integer,                  intent(in)    :: iref, i, irot
        real(sp),                 intent(out)   :: frc(self%kfromto(1):self%kfromto(2))
        complex(sp), pointer :: pft_ref(:,:)
        real(sp),    pointer :: kcorrs(:)
        real(sp) :: sumsqref, sumsqptcl, sqsum_ref
        integer  :: k, ithr
        ithr     =  omp_get_thread_num() + 1
        pft_ref  => self%heap_vars(ithr)%pft_ref
        kcorrs   => self%heap_vars(ithr)%kcorrs
        call self%prep_ref4corr(iref, i, pft_ref, sqsum_ref, self%kfromto(2))
        do k=self%kfromto(1),self%kfromto(2)
            call self%calc_k_corrs(pft_ref, i, k, kcorrs)
            sumsqptcl = sum(csq(self%pfts_ptcls(i,:,k)))
            sumsqref  = sum(csq(pft_ref(:,k)))
            frc(k)    = kcorrs(irot) / sqrt(sumsqref * sumsqptcl)
        end do
    end subroutine genfrc

    !>  \brief  is for generating resolution dependent correlations with shift
    subroutine calc_frc( self, iref, iptcl, irot, shvec, frc )
        use simple_math, only: csq
        class(polarft_corrcalc),  intent(inout) :: self
        integer,                  intent(in)    :: iref, iptcl, irot
        real(sp),                 intent(in)    :: shvec(2)
        real(sp),                 intent(out)   :: frc(self%kfromto(1):self%kfromto(2))
        complex(sp), pointer :: pft_ref(:,:), shmat(:,:)
        real(sp),    pointer :: kcorrs(:), argmat(:,:)
        real(sp) :: sumsqref, sumsqptcl, sqsum_ref
        integer  :: k, ithr, i
        i       = self%pinds(iptcl)
        ithr    =  omp_get_thread_num() + 1
        pft_ref => self%heap_vars(ithr)%pft_ref
        kcorrs  => self%heap_vars(ithr)%kcorrs
        shmat   => self%heap_vars(ithr)%shmat
        argmat  => self%heap_vars(ithr)%argmat
        argmat  =  self%argtransf(:self%pftsz,:) * shvec(1) + self%argtransf(self%pftsz + 1:,:) * shvec(2)
        shmat   =  cmplx(cos(argmat),sin(argmat))
        if( self%with_ctf )then
            if( self%iseven(i) )then
                pft_ref = (self%pfts_refs_even(iref,:,:) * self%ctfmats(i,:,:)) * shmat
            else
                pft_ref = (self%pfts_refs_odd(iref,:,:)  * self%ctfmats(i,:,:)) * shmat
            endif
        else
            if( self%iseven(i) )then
                pft_ref = self%pfts_refs_even(iref,:,:) * shmat
            else
                pft_ref = self%pfts_refs_odd(iref,:,:)  * shmat
            endif
        endif
        call self%prep_ref4corr(iref, i, pft_ref, sqsum_ref, self%kfromto(2))
        do k=self%kfromto(1),self%kfromto(2)
            call self%calc_k_corrs(pft_ref, i, k, kcorrs)
            sumsqptcl = sum(csq(self%pfts_ptcls(i,:,k)))
            sumsqref  = sum(csq(pft_ref(:,k)))
            frc(k)    = kcorrs(irot) / sqrt(sumsqref * sumsqptcl)
        end do
    end subroutine calc_frc

    subroutine gencorrs_cc_1( self, iref, iptcl, cc )
        use simple_math, only: csq
        class(polarft_corrcalc), intent(inout) :: self
        integer,                 intent(in)    :: iref, iptcl
        real(sp),                intent(out)   :: cc(self%nrots)
        complex(sp), pointer :: pft_ref(:,:)
        real(sp),    pointer :: corrs_over_k(:)
        real(sp) :: sqsum_ref
        integer  :: ithr
        ithr         =  omp_get_thread_num() + 1
        pft_ref      => self%heap_vars(ithr)%pft_ref
        corrs_over_k => self%heap_vars(ithr)%corrs_over_k
        call self%prep_ref4corr(iref, self%pinds(iptcl), pft_ref, sqsum_ref, self%kfromto(2))
        call self%calc_corrs_over_k(pft_ref, self%pinds(iptcl), self%kfromto(2), corrs_over_k)
        cc = corrs_over_k / sqrt(sqsum_ref * self%sqsums_ptcls(self%pinds(iptcl)))
    end subroutine gencorrs_cc_1

    subroutine gencorrs_cc_2( self, iref, iptcl, kstop, cc )
        use simple_math, only: csq
        class(polarft_corrcalc), intent(inout) :: self
        integer,                 intent(in)    :: iref, iptcl, kstop
        real,                    intent(out)   :: cc(self%nrots)
        complex(sp), pointer :: pft_ref(:,:)
        real(sp),    pointer :: corrs_over_k(:)
        real(sp) :: sqsum_ref, sqsum_ptcl
        integer  :: ithr
        ithr         =  omp_get_thread_num() + 1
        pft_ref      => self%heap_vars(ithr)%pft_ref
        corrs_over_k => self%heap_vars(ithr)%corrs_over_k
        call self%prep_ref4corr(iref, self%pinds(iptcl), pft_ref, sqsum_ref, kstop)
        call self%calc_corrs_over_k(pft_ref, self%pinds(iptcl), kstop, corrs_over_k)
        sqsum_ptcl = sum(csq(self%pfts_ptcls(self%pinds(iptcl), :, self%kfromto(1):kstop)))
        cc = corrs_over_k / sqrt(sqsum_ref * sqsum_ptcl)
    end subroutine gencorrs_cc_2

    subroutine gencorrs_cc_3( self, iref, iptcl, shvec, cc )
        use simple_math, only: csq
        class(polarft_corrcalc), intent(inout) :: self
        integer,                 intent(in)    :: iref, iptcl
        real(sp),                intent(in)    :: shvec(2)
        real(sp),                intent(out)   :: cc(self%nrots)
        complex(sp), pointer :: pft_ref(:,:), shmat(:,:)
        real(sp),    pointer :: corrs_over_k(:), argmat(:,:)
        real(sp) :: sqsum_ref
        integer  :: ithr
        ithr         = omp_get_thread_num() + 1
        pft_ref      => self%heap_vars(ithr)%pft_ref
        shmat        => self%heap_vars(ithr)%shmat
        corrs_over_k => self%heap_vars(ithr)%corrs_over_k
        argmat       => self%heap_vars(ithr)%argmat
        argmat = self%argtransf(:self%pftsz,:) * shvec(1) + self%argtransf(self%pftsz + 1:,:) * shvec(2)
        shmat  = cmplx(cos(argmat),sin(argmat))
        if( self%with_ctf )then
            if( self%iseven(self%pinds(iptcl)) )then
                pft_ref = (self%pfts_refs_even(iref,:,:) * self%ctfmats(self%pinds(iptcl),:,:)) * shmat
            else
                pft_ref = (self%pfts_refs_odd(iref,:,:)  * self%ctfmats(self%pinds(iptcl),:,:)) * shmat
            endif
        else
            if( self%iseven(self%pinds(iptcl)) )then
                pft_ref = self%pfts_refs_even(iref,:,:) * shmat
            else
                pft_ref = self%pfts_refs_odd(iref,:,:)  * shmat
            endif
        endif
        sqsum_ref = sum(csq(pft_ref))
        call self%calc_corrs_over_k(pft_ref, self%pinds(iptcl), self%kfromto(2), corrs_over_k)
        cc = corrs_over_k  / sqrt(sqsum_ref * self%sqsums_ptcls(self%pinds(iptcl)))
    end subroutine gencorrs_cc_3

    subroutine gencorrs_resnorm_1( self, iref, iptcl, cc )
        class(polarft_corrcalc), intent(inout) :: self
        integer,                 intent(in)    :: iref, iptcl
        real(sp),                intent(out)   :: cc(self%nrots)
        complex(sp), pointer :: pft_ref(:,:)
        real(sp),    pointer :: kcorrs(:)
        real(sp) :: sumsqref, sumsqptcl
        integer  :: k, ithr
        ithr    =  omp_get_thread_num() + 1
        pft_ref => self%heap_vars(ithr)%pft_ref
        kcorrs  => self%heap_vars(ithr)%kcorrs
        call self%prep_ref4corr(iref, self%pinds(iptcl), pft_ref, sumsqref, self%kfromto(2))
        cc(:) = 0.
        if( self%l_cc_bfac )then
            ! with B-factor weighting
            do k=self%kfromto(1),self%kfromto(2)
                call self%calc_k_corrs(pft_ref, self%pinds(iptcl), k, kcorrs)
                sumsqptcl = sum(csq(self%pfts_ptcls(self%pinds(iptcl),:,k)))
                sumsqref  = sum(csq(pft_ref(:,k)))
                ! B-factor weighted correlation
                cc(:) = cc(:) + (kcorrs(:) * self%ptcl_bfac_weights(self%pinds(iptcl),k)) / sqrt(sumsqref * sumsqptcl)
            end do
            cc(:) = cc(:) / self%ptcl_bfac_norms(self%pinds(iptcl))
        else
            ! without B-factor weighting
            do k=self%kfromto(1),self%kfromto(2)
                call self%calc_k_corrs(pft_ref, self%pinds(iptcl), k, kcorrs)
                sumsqptcl = sum(csq(self%pfts_ptcls(self%pinds(iptcl),:,k)))
                sumsqref  = sum(csq(pft_ref(:,k)))
                ! all rotational correlations in one shell: kcorrs(:) / sqrt(sumsqref * sumsqptcl)
                ! sum over shells
                cc(:) = cc(:) + kcorrs(:) / sqrt(sumsqref * sumsqptcl)
            end do
            cc(:) = cc(:) / real(self%nk)
        endif
    end subroutine gencorrs_resnorm_1

    subroutine gencorrs_resnorm_2( self, iref, iptcl, kstop, cc )
        class(polarft_corrcalc), intent(inout) :: self
        integer,                 intent(in)    :: iref, iptcl, kstop
        real(sp),                intent(out)   :: cc(self%nrots)
        complex(sp), pointer :: pft_ref(:,:)
        real(sp),    pointer :: kcorrs(:)
        real(sp) :: sumsqref, sumsqptcl
        integer  :: k, ithr
        ithr    =  omp_get_thread_num() + 1
        pft_ref => self%heap_vars(ithr)%pft_ref
        kcorrs  => self%heap_vars(ithr)%kcorrs
        call self%prep_ref4corr(iref, self%pinds(iptcl), pft_ref, sumsqref, kstop)
        cc(:) = 0.
        if( self%l_cc_bfac )then
            ! with B-factor weighting
            do k=self%kfromto(1),kstop
                call self%calc_k_corrs(pft_ref, self%pinds(iptcl), k, kcorrs)
                sumsqptcl = sum(csq(self%pfts_ptcls(self%pinds(iptcl),:,k)))
                sumsqref  = sum(csq(pft_ref(:,k)))
                ! B-factor weighted correlation
                cc(:) = cc(:) + (kcorrs(:) * self%ptcl_bfac_weights(self%pinds(iptcl),k)) / sqrt(sumsqref * sumsqptcl)
            end do
            cc(:) = cc(:) / sum(self%ptcl_bfac_weights(self%pinds(iptcl),self%kfromto(1):kstop))
        else
            ! without B-factor weighting
            do k=self%kfromto(1),kstop
                call self%calc_k_corrs(pft_ref, self%pinds(iptcl), k, kcorrs)
                sumsqptcl = sum(csq(self%pfts_ptcls(self%pinds(iptcl),:,k)))
                sumsqref  = sum(csq(pft_ref(:,k)))
                ! all rotational correlations in one shell: kcorrs(:) / sqrt(sumsqref * sumsqptcl)
                ! sum over shells
                cc(:) = cc(:) + kcorrs(:) / sqrt(sumsqref * sumsqptcl)
            end do
            cc(:) = cc(:) / real(self%kfromto(2) - kstop + 1)
        endif
    end subroutine gencorrs_resnorm_2

    subroutine gencorrs_resnorm_3( self, iref, iptcl, shvec, cc )
        use simple_math, only: csq
        class(polarft_corrcalc), intent(inout) :: self
        integer,                 intent(in)    :: iref, iptcl
        real(sp),                intent(in)    :: shvec(2)
        real(sp),                intent(out)   :: cc(self%nrots)
        complex(sp), pointer :: pft_ref(:,:), shmat(:,:)
        real(sp),    pointer :: kcorrs(:), argmat(:,:)
        real(sp) :: sumsqptcl, sumsqref
        integer  :: ithr, k
        ithr    = omp_get_thread_num() + 1
        pft_ref => self%heap_vars(ithr)%pft_ref
        shmat   => self%heap_vars(ithr)%shmat
        kcorrs  => self%heap_vars(ithr)%kcorrs
        argmat  => self%heap_vars(ithr)%argmat
        argmat = self%argtransf(:self%pftsz,:) * shvec(1) + self%argtransf(self%pftsz + 1:,:) * shvec(2)
        shmat  = cmplx(cos(argmat),sin(argmat))
        if( self%with_ctf )then
            if( self%iseven(self%pinds(iptcl)) )then
                pft_ref = (self%pfts_refs_even(iref,:,:) * self%ctfmats(self%pinds(iptcl),:,:)) * shmat
            else
                pft_ref = (self%pfts_refs_odd(iref,:,:)  * self%ctfmats(self%pinds(iptcl),:,:)) * shmat
            endif
        else
            if( self%iseven(self%pinds(iptcl)) )then
                pft_ref = self%pfts_refs_even(iref,:,:) * shmat
            else
                pft_ref = self%pfts_refs_odd(iref,:,:)  * shmat
            endif
        endif
        cc(:) = 0.
        if( self%l_cc_bfac )then
            ! with B-factor weighting
            do k=self%kfromto(1),self%kfromto(2)
                call self%calc_k_corrs(pft_ref, self%pinds(iptcl), k, kcorrs)
                sumsqptcl = sum(csq(self%pfts_ptcls(self%pinds(iptcl),:,k)))
                sumsqref  = sum(csq(pft_ref(:,k)))
                ! B-factor weighted correlation
                cc(:) = cc(:) + (kcorrs(:) * self%ptcl_bfac_weights(self%pinds(iptcl),k)) / sqrt(sumsqref * sumsqptcl)
            end do
            cc(:) = cc(:) / self%ptcl_bfac_norms(self%pinds(iptcl))
        else
            ! without B-factor weighting
            do k=self%kfromto(1),self%kfromto(2)
                call self%calc_k_corrs(pft_ref, self%pinds(iptcl), k, kcorrs)
                sumsqptcl = sum(csq(self%pfts_ptcls(self%pinds(iptcl),:,k)))
                sumsqref  = sum(csq(pft_ref(:,k)))
                ! all rotational correlations in one shell: kcorrs(:) / sqrt(sumsqref * sumsqptcl)
                ! sum over shells
                cc(:) = cc(:) + kcorrs(:) / sqrt(sumsqref * sumsqptcl)
            end do
            cc(:) = cc(:) / real(self%nk)
        endif
    end subroutine gencorrs_resnorm_3

    subroutine gencorrs_1( self, iref, iptcl, cc )
        class(polarft_corrcalc), intent(inout) :: self
        integer,                 intent(in)    :: iref, iptcl
        real(sp),                intent(out)   :: cc(self%nrots)
        if( self%l_cc_objfun )then
            call self%gencorrs_cc_1(iref, iptcl, cc)
        else
            call self%gencorrs_resnorm_1(iref, iptcl, cc)
        endif
    end subroutine gencorrs_1

    subroutine gencorrs_2( self, iref, iptcl, kstop, cc )
        class(polarft_corrcalc), intent(inout) :: self
        integer,                 intent(in)    :: iref, iptcl, kstop
        real,                    intent(out)   :: cc(self%nrots)
        if( self%l_cc_objfun )then
            call self%gencorrs_cc_2(iref, iptcl, kstop, cc)
        else
            call self%gencorrs_resnorm_2(iref, iptcl, kstop, cc)
        endif
    end subroutine gencorrs_2

    subroutine gencorrs_3( self, iref, iptcl, shvec, cc )
        class(polarft_corrcalc), intent(inout) :: self
        integer,                 intent(in)    :: iref, iptcl
        real(sp),                intent(in)    :: shvec(2)
        real(sp),                intent(out)   :: cc(self%nrots)
        if( self%l_cc_objfun )then
            call self%gencorrs_cc_3(iref, iptcl, shvec, cc)
        else
            call self%gencorrs_resnorm_3(iref, iptcl, shvec, cc)
        endif
    end subroutine gencorrs_3

    !< brief  generates correlation for one specific rotation angle
    function gencorr_for_rot( self, iref, iptcl, shvec, irot ) result( cc )
        class(polarft_corrcalc), intent(inout) :: self
        integer,                 intent(in)    :: iref, iptcl
        real(sp),                intent(in)    :: shvec(2)
        integer,                 intent(in)    :: irot
        real(sp)                               :: cc
        if( self%l_cc_objfun )then
            cc = self%gencorr_cc_for_rot( iref, iptcl, shvec, irot )
        else
            cc = self%gencorr_resnorm_for_rot( iref, iptcl, shvec, irot )
        endif
    end function gencorr_for_rot

    !< brief  generates correlation for one specific rotation angle, double precision
    function gencorr_for_rot_8( self, iref, iptcl, shvec, irot ) result( cc )
        class(polarft_corrcalc), intent(inout) :: self
        integer,                 intent(in)    :: iref, iptcl
        real(dp),                intent(in)    :: shvec(2)
        integer,                 intent(in)    :: irot
        real(dp)                               :: cc
        if( self%l_cc_objfun )then
            cc = self%gencorr_cc_for_rot_8( iref, iptcl, shvec, irot )
        else
            cc = self%gencorr_resnorm_for_rot_8( iref, iptcl, shvec, irot )
        endif
    end function gencorr_for_rot_8

    !< brief  generates correlation for one specific rotation angle
    function gencorr_cc_for_rot( self, iref, iptcl, shvec, irot ) result( cc )
        use simple_math, only: csq
        class(polarft_corrcalc), intent(inout) :: self
        integer,                 intent(in)    :: iref, iptcl
        real(sp),                intent(in)    :: shvec(2)
        integer,                 intent(in)    :: irot
        complex(sp), pointer :: pft_ref(:,:), shmat(:,:)
        real(sp),    pointer :: argmat(:,:)
        real(sp) :: cc, corr, sqsum_ref
        integer  :: ithr
        ithr    =  omp_get_thread_num() + 1
        pft_ref => self%heap_vars(ithr)%pft_ref
        shmat   => self%heap_vars(ithr)%shmat
        argmat  => self%heap_vars(ithr)%argmat
        argmat  =  self%argtransf(:self%pftsz,:) * shvec(1) + self%argtransf(self%pftsz + 1:,:) * shvec(2)
        shmat   =  cmplx(cos(argmat),sin(argmat))
        if( self%with_ctf )then
            if( self%iseven(self%pinds(iptcl)) )then
                pft_ref = (self%pfts_refs_even(iref,:,:) * self%ctfmats(self%pinds(iptcl),:,:)) * shmat
            else
                pft_ref = (self%pfts_refs_odd(iref,:,:)  * self%ctfmats(self%pinds(iptcl),:,:)) * shmat
            endif
        else
            if( self%iseven(self%pinds(iptcl)) )then
                pft_ref = self%pfts_refs_even(iref,:,:) * shmat
            else
                pft_ref = self%pfts_refs_odd(iref,:,:)  * shmat
            endif
        endif
        sqsum_ref = sum(csq(pft_ref))
        corr      = self%calc_corr_for_rot(pft_ref, self%pinds(iptcl), self%kfromto(2), irot)
        cc        = corr  / sqrt(sqsum_ref * self%sqsums_ptcls(self%pinds(iptcl)))
    end function gencorr_cc_for_rot

    !< brief  generates correlation for one specific rotation angle, double precision
    function gencorr_cc_for_rot_8( self, iref, iptcl, shvec, irot ) result( cc )
        use simple_math, only: csq
        class(polarft_corrcalc), intent(inout) :: self
        integer,                 intent(in)    :: iref, iptcl
        real(dp),                intent(in)    :: shvec(2)
        integer,                 intent(in)    :: irot
        complex(dp), pointer :: pft_ref(:,:), shmat(:,:)
        real(dp),    pointer :: argmat(:,:)
        real(dp) :: cc, corr, sqsum_ref
        integer  :: ithr
        ithr    =  omp_get_thread_num() + 1
        pft_ref => self%heap_vars(ithr)%pft_ref_8
        shmat   => self%heap_vars(ithr)%shmat_8
        argmat  => self%heap_vars(ithr)%argmat_8
        argmat  =  self%argtransf(:self%pftsz,:) * shvec(1) + self%argtransf(self%pftsz + 1:,:) * shvec(2)
        shmat   =  cmplx(cos(argmat),sin(argmat),dp)
        if( self%with_ctf )then
            if( self%iseven(self%pinds(iptcl)) )then
                pft_ref = (self%pfts_refs_even(iref,:,:) * self%ctfmats(self%pinds(iptcl),:,:)) * shmat
            else
                pft_ref = (self%pfts_refs_odd(iref,:,:)  * self%ctfmats(self%pinds(iptcl),:,:)) * shmat
            endif
        else
            if( self%iseven(self%pinds(iptcl)) )then
                pft_ref = self%pfts_refs_even(iref,:,:) * shmat
            else
                pft_ref = self%pfts_refs_odd(iref,:,:)  * shmat
            endif
        endif
        sqsum_ref = sum(csq(pft_ref))
        corr      = self%calc_corr_for_rot_8(pft_ref, self%pinds(iptcl), self%kfromto(2), irot)
        cc        = corr  / sqrt(sqsum_ref * self%sqsums_ptcls(self%pinds(iptcl)))
    end function gencorr_cc_for_rot_8

    !< brief  generates correlation for one specific rotation angle
    function gencorr_resnorm_for_rot( self, iref, iptcl, shvec, irot ) result( cc )
        use simple_math, only: csq
        class(polarft_corrcalc), intent(inout) :: self
        integer,                 intent(in)    :: iref, iptcl
        real(sp),                intent(in)    :: shvec(2)
        integer,                 intent(in)    :: irot
        complex(sp), pointer :: pft_ref(:,:), shmat(:,:)
        real(sp),    pointer :: argmat(:,:)
        real(sp) :: cc, corrk, sqsumk_ref, sqsumk_ptcl
        integer  :: ithr, k
        ithr    =  omp_get_thread_num() + 1
        pft_ref => self%heap_vars(ithr)%pft_ref
        shmat   => self%heap_vars(ithr)%shmat
        argmat  => self%heap_vars(ithr)%argmat
        argmat  =  self%argtransf(:self%pftsz,:) * shvec(1) + self%argtransf(self%pftsz + 1:,:) * shvec(2)
        shmat   =  cmplx(cos(argmat),sin(argmat))
        if( self%with_ctf )then
            if( self%iseven(self%pinds(iptcl)) )then
                pft_ref = (self%pfts_refs_even(iref,:,:) * self%ctfmats(self%pinds(iptcl),:,:)) * shmat
            else
                pft_ref = (self%pfts_refs_odd(iref,:,:)  * self%ctfmats(self%pinds(iptcl),:,:)) * shmat
            endif
        else
            if( self%iseven(self%pinds(iptcl)) )then
                pft_ref = self%pfts_refs_even(iref,:,:) * shmat
            else
                pft_ref = self%pfts_refs_odd(iref,:,:)  * shmat
            endif
        endif
        cc = 0.0
        if( self%l_cc_bfac )then
            ! with B-factor weighting
            do k = self%kfromto(1),self%kfromto(2)
                sqsumk_ref  = sum(csq(pft_ref(:,k)))
                sqsumk_ptcl = sum(csq(self%pfts_ptcls(self%pinds(iptcl),:,k)))
                corrk       = self%calc_corrk_for_rot(pft_ref, self%pinds(iptcl), self%kfromto(2), k, irot)
                cc          = cc + (corrk * self%ptcl_bfac_weights(self%pinds(iptcl),k)) / sqrt(sqsumk_ref * sqsumk_ptcl)
            end do
            cc = cc / self%ptcl_bfac_norms(self%pinds(iptcl))
        else
            ! without B-factor weighting
            do k = self%kfromto(1),self%kfromto(2)
                sqsumk_ref  = sum(csq(pft_ref(:,k)))
                sqsumk_ptcl = sum(csq(self%pfts_ptcls(self%pinds(iptcl),:,k)))
                corrk       = self%calc_corrk_for_rot(pft_ref, self%pinds(iptcl), self%kfromto(2), k, irot)
                cc          = cc + corrk / sqrt(sqsumk_ref * sqsumk_ptcl)
            end do
            cc = cc / real(self%nk)
        endif
    end function gencorr_resnorm_for_rot

    !< brief  generates correlation for one specific rotation angle, double precision
    function gencorr_resnorm_for_rot_8( self, iref, iptcl, shvec, irot ) result( cc )
        use simple_math, only: csq
        class(polarft_corrcalc), intent(inout) :: self
        integer,                 intent(in)    :: iref, iptcl
        real(dp),                intent(in)    :: shvec(2)
        integer,                 intent(in)    :: irot
        complex(dp), pointer :: pft_ref(:,:), shmat(:,:)
        real(dp),    pointer :: argmat(:,:)
        real(dp) :: cc, corrk, sqsumk_ref, sqsumk_ptcl
        integer  :: ithr, k
        ithr    =  omp_get_thread_num() + 1
        pft_ref => self%heap_vars(ithr)%pft_ref_8
        shmat   => self%heap_vars(ithr)%shmat_8
        argmat  => self%heap_vars(ithr)%argmat_8
        argmat  =  self%argtransf(:self%pftsz,:) * shvec(1) + self%argtransf(self%pftsz + 1:,:) * shvec(2)
        shmat   =  cmplx(cos(argmat),sin(argmat),dp)
        if( self%with_ctf )then
            if( self%iseven(self%pinds(iptcl)) )then
                pft_ref = (self%pfts_refs_even(iref,:,:) * self%ctfmats(self%pinds(iptcl),:,:)) * shmat
            else
                pft_ref = (self%pfts_refs_odd(iref,:,:)  * self%ctfmats(self%pinds(iptcl),:,:)) * shmat
            endif
        else
            if( self%iseven(self%pinds(iptcl)) )then
                pft_ref = self%pfts_refs_even(iref,:,:) * shmat
            else
                pft_ref = self%pfts_refs_odd(iref,:,:)  * shmat
            endif
        endif
        cc = 0.0_dp
        if( self%l_cc_bfac )then
            ! with B-factor weighting
            do k = self%kfromto(1),self%kfromto(2)
                sqsumk_ref  = sum(csq(pft_ref(:,k)))
                sqsumk_ptcl = sum(csq(self%pfts_ptcls(self%pinds(iptcl),:,k)))
                corrk       = self%calc_corrk_for_rot_8(pft_ref, self%pinds(iptcl), self%kfromto(2), k, irot)
                cc          = cc + (corrk * real(self%ptcl_bfac_weights(self%pinds(iptcl),k), kind=dp)) / sqrt(sqsumk_ref * sqsumk_ptcl)
            end do
            cc = cc / real(self%ptcl_bfac_norms(self%pinds(iptcl)), kind=dp)
        else
            ! without B-factor weighting
            do k = self%kfromto(1),self%kfromto(2)
                sqsumk_ref = sum(csq(pft_ref(:,k)))
                sqsumk_ptcl = sum(csq(self%pfts_ptcls(self%pinds(iptcl),:,k)))
                corrk      = self%calc_corrk_for_rot_8(pft_ref, self%pinds(iptcl), self%kfromto(2), k, irot)
                cc         = cc + corrk / sqrt(sqsumk_ref * sqsumk_ptcl)
            end do
            cc = cc / real(self%nk, kind=dp)
        endif
    end function gencorr_resnorm_for_rot_8

    !< brief  calculates correlations and gradient for origin shift
    subroutine gencorrs_cc_grad( self, iref, iptcl, shvec, cc, grad )
        use simple_math, only: csq
        class(polarft_corrcalc), intent(inout) :: self
        integer,                 intent(in)    :: iref, iptcl
        real(sp),                intent(in)    :: shvec(2)
        real(sp),                intent(out)   :: cc(self%nrots)
        real(sp),                intent(out)   :: grad(2, self%nrots)
        complex(sp), pointer :: pft_ref(:,:), pft_ref_tmp(:,:), shmat(:,:)
        real(sp),    pointer :: corrs_over_k(:), argmat(:,:)
        real(sp) :: denom
        integer  :: ithr
        ithr         =  omp_get_thread_num() + 1
        pft_ref      => self%heap_vars(ithr)%pft_ref
        pft_ref_tmp  => self%heap_vars(ithr)%pft_ref_tmp
        shmat        => self%heap_vars(ithr)%shmat
        corrs_over_k => self%heap_vars(ithr)%corrs_over_k
        argmat       => self%heap_vars(ithr)%argmat
        argmat       =  self%argtransf(:self%pftsz,:) * shvec(1) + self%argtransf(self%pftsz + 1:,:) * shvec(2)
        shmat        =  cmplx(cos(argmat),sin(argmat))
        if( self%with_ctf )then
            if( self%iseven(self%pinds(iptcl)) )then
                pft_ref = (self%pfts_refs_even(iref,:,:) * self%ctfmats(self%pinds(iptcl),:,:)) * shmat
            else
                pft_ref = (self%pfts_refs_odd(iref,:,:)  * self%ctfmats(self%pinds(iptcl),:,:)) * shmat
            endif
        else
            if( self%iseven(self%pinds(iptcl)) )then
                pft_ref = self%pfts_refs_even(iref,:,:) * shmat
            else
                pft_ref = self%pfts_refs_odd(iref,:,:)  * shmat
            endif
        endif
        denom       = sqrt(sum(csq(pft_ref)) * self%sqsums_ptcls(self%pinds(iptcl)))
        call self%calc_corrs_over_k(pft_ref, self%pinds(iptcl), self%kfromto(2), corrs_over_k)
        cc          = corrs_over_k / denom
        pft_ref_tmp = pft_ref * (0., 1.) * self%argtransf(:self%pftsz,:)
        call self%calc_corrs_over_k(pft_ref_tmp, self%pinds(iptcl), self%kfromto(2), corrs_over_k)
        grad(1,:)   = corrs_over_k / denom
        pft_ref_tmp = pft_ref * (0., 1.) * self%argtransf(self%pftsz + 1:,:)
        call self%calc_corrs_over_k(pft_ref_tmp, self%pinds(iptcl), self%kfromto(2), corrs_over_k)
        grad(2,:)   = corrs_over_k / denom
    end subroutine gencorrs_cc_grad

    !< brief  calculates only gradient for correlations
    subroutine gencorrs_cc_grad_only( self, iref, iptcl, shvec, grad )
        use simple_math, only: csq
        class(polarft_corrcalc), intent(inout) :: self
        integer,                 intent(in)    :: iref, iptcl
        real(sp),                intent(in)    :: shvec(2)
        real(sp),                intent(out)   :: grad(2, self%nrots)
        complex(sp), pointer :: pft_ref(:,:), pft_ref_tmp(:,:), shmat(:,:)
        real(sp),    pointer :: corrs_over_k(:), argmat(:,:)
        real(sp) :: denom
        integer  :: ithr
        ithr         =  omp_get_thread_num() + 1
        pft_ref      => self%heap_vars(ithr)%pft_ref
        pft_ref_tmp  => self%heap_vars(ithr)%pft_ref_tmp
        shmat        => self%heap_vars(ithr)%shmat
        corrs_over_k => self%heap_vars(ithr)%corrs_over_k
        argmat       => self%heap_vars(ithr)%argmat
        argmat       =  self%argtransf(:self%pftsz,:) * shvec(1) + self%argtransf(self%pftsz + 1:,:) * shvec(2)
        shmat        =  cmplx(cos(argmat),sin(argmat))
        if( self%with_ctf )then
            if( self%iseven(self%pinds(iptcl)) )then
                pft_ref = (self%pfts_refs_even(iref,:,:) * self%ctfmats(self%pinds(iptcl),:,:)) * shmat
            else
                pft_ref = (self%pfts_refs_odd(iref,:,:)  * self%ctfmats(self%pinds(iptcl),:,:)) * shmat
            endif
        else
            if( self%iseven(self%pinds(iptcl)) )then
                pft_ref = self%pfts_refs_even(iref,:,:) * shmat
            else
                pft_ref = self%pfts_refs_odd(iref,:,:)  * shmat
            endif
        endif
        denom       = sqrt(sum(csq(pft_ref)) * self%sqsums_ptcls(self%pinds(iptcl)))
        pft_ref_tmp = pft_ref * (0., 1.) * self%argtransf(:self%pftsz,:)
        call self%calc_corrs_over_k(pft_ref_tmp, self%pinds(iptcl), self%kfromto(2), corrs_over_k)
        grad(1,:)   = corrs_over_k / denom
        pft_ref_tmp = pft_ref * (0., 1.) * self%argtransf(self%pftsz + 1:,:)
        call self%calc_corrs_over_k(pft_ref_tmp, self%pinds(iptcl), self%kfromto(2), corrs_over_k)
        grad(2,:)   = corrs_over_k / denom
    end subroutine gencorrs_cc_grad_only

    !< brief  calculates correlation and gradient for origin shift, for one specific rotation angle
    subroutine gencorr_grad_for_rot( self, iref, iptcl, shvec, irot, f, grad )
        class(polarft_corrcalc), intent(inout) :: self
        integer,                 intent(in)    :: iref, iptcl
        real(sp),                intent(in)    :: shvec(2)
        integer,                 intent(in)    :: irot
        real(sp),                intent(out)   :: f, grad(2)
        if( self%l_cc_objfun )then
            call self%gencorr_cc_grad_for_rot( iref, iptcl, shvec, irot, f, grad )
        else
            call self%gencorr_resnorm_grad_for_rot( iref, iptcl, shvec, irot, f, grad )
        endif
    end subroutine gencorr_grad_for_rot

    !< brief  calculates correlation and gradient for origin shift, for one specific rotation angle, double precision
    subroutine gencorr_grad_for_rot_8( self, iref, iptcl, shvec, irot, f, grad )
        use simple_math, only: csq
        class(polarft_corrcalc), intent(inout) :: self
        integer,                 intent(in)    :: iref, iptcl
        real(dp),                intent(in)    :: shvec(2)
        integer,                 intent(in)    :: irot
        real(dp),                intent(out)   :: f, grad(2)
        if( self%l_cc_objfun )then
            call self%gencorr_cc_grad_for_rot_8( iref, iptcl, shvec, irot, f, grad )
        else
            call self%gencorr_resnorm_grad_for_rot_8( iref, iptcl, shvec, irot, f, grad )
        endif
    end subroutine gencorr_grad_for_rot_8

    !< brief  calculates correlation and gradient for origin shift, for one specific rotation angle
    subroutine gencorr_cc_grad_for_rot( self, iref, iptcl, shvec, irot, f, grad )
        use simple_math, only: csq
        class(polarft_corrcalc), intent(inout) :: self
        integer,                 intent(in)    :: iref, iptcl
        real(sp),                intent(in)    :: shvec(2)
        integer,                 intent(in)    :: irot
        real(sp),                intent(out)   :: f, grad(2)
        complex(sp), pointer :: pft_ref(:,:), pft_ref_tmp(:,:), shmat(:,:)
        real(sp),    pointer :: argmat(:,:)
        real(sp) :: denom, corr
        integer  :: ithr
        ithr        =  omp_get_thread_num() + 1
        pft_ref     => self%heap_vars(ithr)%pft_ref
        pft_ref_tmp => self%heap_vars(ithr)%pft_ref_tmp
        shmat       => self%heap_vars(ithr)%shmat
        argmat      => self%heap_vars(ithr)%argmat
        argmat      =  self%argtransf(:self%pftsz,:) * shvec(1) + self%argtransf(self%pftsz + 1:,:) * shvec(2)
        shmat       =  cmplx(cos(argmat),sin(argmat))
        if( self%with_ctf )then
            if( self%iseven(self%pinds(iptcl)) )then
                pft_ref = (self%pfts_refs_even(iref,:,:) * self%ctfmats(self%pinds(iptcl),:,:)) * shmat
            else
                pft_ref = (self%pfts_refs_odd(iref,:,:)  * self%ctfmats(self%pinds(iptcl),:,:)) * shmat
            endif
        else
            if( self%iseven(self%pinds(iptcl)) )then
                pft_ref = self%pfts_refs_even(iref,:,:) * shmat
            else
                pft_ref = self%pfts_refs_odd(iref,:,:)  * shmat
            endif
        endif
        denom       = sqrt(sum(csq(pft_ref)) * self%sqsums_ptcls(self%pinds(iptcl)))
        corr        = self%calc_corr_for_rot(pft_ref, self%pinds(iptcl), self%kfromto(2), irot)
        f           = corr  / denom
        pft_ref_tmp = pft_ref * (0., 1.) * self%argtransf(:self%pftsz,:)
        corr        = self%calc_corr_for_rot(pft_ref_tmp, self%pinds(iptcl), self%kfromto(2), irot)
        grad(1)     = corr / denom
        pft_ref_tmp = pft_ref * (0., 1.) * self%argtransf(self%pftsz + 1:,:)
        corr        = self%calc_corr_for_rot(pft_ref_tmp, self%pinds(iptcl), self%kfromto(2), irot)
        grad(2)     = corr / denom
    end subroutine gencorr_cc_grad_for_rot

    !< brief  calculates correlation and gradient for origin shift, for one specific rotation angle, double precision
    subroutine gencorr_cc_grad_for_rot_8( self, iref, iptcl, shvec, irot, f, grad )
        use simple_math, only: csq
        class(polarft_corrcalc), intent(inout) :: self
        integer,                 intent(in)    :: iref, iptcl
        real(dp),                intent(in)    :: shvec(2)
        integer,                 intent(in)    :: irot
        real(dp),                intent(out)   :: f, grad(2)
        complex(dp), pointer :: pft_ref(:,:), pft_ref_tmp(:,:), shmat(:,:)
        real(dp),    pointer :: argmat(:,:)
        real(dp) :: corr, denom
        integer  :: ithr
        ithr        =  omp_get_thread_num() + 1
        pft_ref     => self%heap_vars(ithr)%pft_ref_8
        pft_ref_tmp => self%heap_vars(ithr)%pft_ref_tmp_8
        shmat       => self%heap_vars(ithr)%shmat_8
        argmat      => self%heap_vars(ithr)%argmat_8
        argmat      =  self%argtransf(:self%pftsz,:) * shvec(1) + self%argtransf(self%pftsz + 1:,:) * shvec(2)
        shmat       =  cmplx(cos(argmat),sin(argmat),dp)
        if( self%with_ctf )then
            if( self%iseven(self%pinds(iptcl)) )then
                pft_ref = (self%pfts_refs_even(iref,:,:) * self%ctfmats(self%pinds(iptcl),:,:)) * shmat
            else
                pft_ref = (self%pfts_refs_odd(iref,:,:)  * self%ctfmats(self%pinds(iptcl),:,:)) * shmat
            endif
        else
            if( self%iseven(self%pinds(iptcl)) )then
                pft_ref = self%pfts_refs_even(iref,:,:) * shmat
            else
                pft_ref = self%pfts_refs_odd(iref,:,:)  * shmat
            endif
        endif
        denom       = sqrt(sum(csq(pft_ref)) * self%sqsums_ptcls(self%pinds(iptcl)))
        corr        = self%calc_corr_for_rot_8(pft_ref, self%pinds(iptcl), self%kfromto(2), irot)
        f           = corr  / denom
        pft_ref_tmp = pft_ref * (0., 1.) * self%argtransf(:self%pftsz,:)
        corr        = self%calc_corr_for_rot_8(pft_ref_tmp, self%pinds(iptcl), self%kfromto(2), irot)
        grad(1)     = corr / denom
        pft_ref_tmp = pft_ref * (0., 1.) * self%argtransf(self%pftsz + 1:,:)
        corr        = self%calc_corr_for_rot_8(pft_ref_tmp, self%pinds(iptcl), self%kfromto(2), irot)
        grad(2)     = corr / denom
    end subroutine gencorr_cc_grad_for_rot_8

    !< brief  calculates correlation and gradient for origin shift, for one specific rotation angle
    subroutine gencorr_resnorm_grad_for_rot( self, iref, iptcl, shvec, irot, f, grad )
        use simple_math, only: csq
        class(polarft_corrcalc), intent(inout) :: self
        integer,                 intent(in)    :: iref, iptcl
        real(sp),                intent(in)    :: shvec(2)
        integer,                 intent(in)    :: irot
        real(sp),                intent(out)   :: f, grad(2)
        complex(sp), pointer :: pft_ref(:,:), pft_ref_tmp1(:,:), pft_ref_tmp2(:,:), shmat(:,:)
        real(sp),    pointer :: argmat(:,:)
        real(sp) :: sqsumk_ref, corrk, sqsumk_ptcl, denom
        integer  :: ithr, k
        ithr         =  omp_get_thread_num() + 1
        pft_ref      => self%heap_vars(ithr)%pft_ref
        pft_ref_tmp1 => self%heap_vars(ithr)%pft_ref_tmp1
        pft_ref_tmp2 => self%heap_vars(ithr)%pft_ref_tmp2
        shmat        => self%heap_vars(ithr)%shmat
        argmat       => self%heap_vars(ithr)%argmat
        argmat       =  self%argtransf(:self%pftsz,:) * shvec(1) + self%argtransf(self%pftsz + 1:,:) * shvec(2)
        shmat        =  cmplx(cos(argmat),sin(argmat))
        if( self%with_ctf )then
            if( self%iseven(self%pinds(iptcl)) )then
                pft_ref = (self%pfts_refs_even(iref,:,:) * self%ctfmats(self%pinds(iptcl),:,:)) * shmat
            else
                pft_ref = (self%pfts_refs_odd(iref,:,:)  * self%ctfmats(self%pinds(iptcl),:,:)) * shmat
            endif
        else
            if( self%iseven(self%pinds(iptcl)) )then
                pft_ref = self%pfts_refs_even(iref,:,:) * shmat
            else
                pft_ref = self%pfts_refs_odd(iref,:,:)  * shmat
            endif
        endif
        f    = 0.0
        grad = 0.0
        pft_ref_tmp1 = pft_ref * (0., 1.) * self%argtransf(:self%pftsz,:)
        pft_ref_tmp2 = pft_ref * (0., 1.) * self%argtransf(self%pftsz + 1:,:)
        if( self%l_cc_bfac )then
            do k = self%kfromto(1), self%kfromto(2)
                sqsumk_ref  = sum(csq(pft_ref(:,k)))
                sqsumk_ptcl = sum(csq(self%pfts_ptcls(self%pinds(iptcl),:,k)))
                denom       = sqrt(sqsumk_ref * sqsumk_ptcl)
                corrk       = self%calc_corrk_for_rot(pft_ref, self%pinds(iptcl), self%kfromto(2), k, irot)
                f           = f +       (corrk * self%ptcl_bfac_weights(self%pinds(iptcl),k)) / denom
                corrk       = self%calc_corrk_for_rot(pft_ref_tmp1, self%pinds(iptcl), self%kfromto(2), k, irot)
                grad(1)     = grad(1) + (corrk * self%ptcl_bfac_weights(self%pinds(iptcl),k)) / denom
                corrk       = self%calc_corrk_for_rot(pft_ref_tmp2, self%pinds(iptcl), self%kfromto(2), k, irot)
                grad(2)     = grad(2) + (corrk * self%ptcl_bfac_weights(self%pinds(iptcl),k)) / denom
            end do
            f = f             / real(self%ptcl_bfac_norms(self%pinds(iptcl)))
            grad(:) = grad(:) / real(self%ptcl_bfac_norms(self%pinds(iptcl)))
        else
            do k = self%kfromto(1), self%kfromto(2)
                sqsumk_ref  = sum(csq(pft_ref(:,k)))
                sqsumk_ptcl = sum(csq(self%pfts_ptcls(self%pinds(iptcl),:,k)))
                denom       = sqrt(sqsumk_ref * sqsumk_ptcl)
                corrk       = self%calc_corrk_for_rot(pft_ref, self%pinds(iptcl), self%kfromto(2), k, irot)
                f           = f + corrk       / denom
                corrk       = self%calc_corrk_for_rot(pft_ref_tmp1, self%pinds(iptcl), self%kfromto(2), k, irot)
                grad(1)     = grad(1) + corrk / denom
                corrk       = self%calc_corrk_for_rot(pft_ref_tmp2, self%pinds(iptcl), self%kfromto(2), k, irot)
                grad(2)     = grad(2) + corrk / denom
            end do
            f = f             / real(self%nk)
            grad(:) = grad(:) / real(self%nk)
        endif
    end subroutine gencorr_resnorm_grad_for_rot

    !< brief  calculates correlation and gradient for origin shift, for one specific rotation angle, double precision
    subroutine gencorr_resnorm_grad_for_rot_8( self, iref, iptcl, shvec, irot, f, grad )
        use simple_math, only: csq
        class(polarft_corrcalc), intent(inout) :: self
        integer,                 intent(in)    :: iref, iptcl
        real(dp),                intent(in)    :: shvec(2)
        integer,                 intent(in)    :: irot
        real(dp),                intent(out)   :: f, grad(2)
        complex(dp), pointer :: pft_ref(:,:), pft_ref_tmp1(:,:), pft_ref_tmp2(:,:), shmat(:,:)
        real(dp),    pointer :: argmat(:,:)
        real(dp) :: sqsumk_ref, corrk, sqsumk_ptcl, denom
        integer  :: ithr, k
        ithr         =  omp_get_thread_num() + 1
        pft_ref      => self%heap_vars(ithr)%pft_ref_8
        pft_ref_tmp1 => self%heap_vars(ithr)%pft_ref_tmp1_8
        pft_ref_tmp2 => self%heap_vars(ithr)%pft_ref_tmp2_8
        shmat        => self%heap_vars(ithr)%shmat_8
        argmat       => self%heap_vars(ithr)%argmat_8
        argmat       =  self%argtransf(:self%pftsz,:) * shvec(1) + self%argtransf(self%pftsz + 1:,:) * shvec(2)
        shmat        =  cmplx(cos(argmat),sin(argmat),dp)
        if( self%with_ctf )then
            if( self%iseven(self%pinds(iptcl)) )then
                pft_ref = (self%pfts_refs_even(iref,:,:) * self%ctfmats(self%pinds(iptcl),:,:)) * shmat
            else
                pft_ref = (self%pfts_refs_odd(iref,:,:)  * self%ctfmats(self%pinds(iptcl),:,:)) * shmat
            endif
        else
            if( self%iseven(self%pinds(iptcl)) )then
                pft_ref = self%pfts_refs_even(iref,:,:) * shmat
            else
                pft_ref = self%pfts_refs_odd(iref,:,:)  * shmat
            endif
        endif
        f    = 0.0_dp
        grad = 0.0_dp
        pft_ref_tmp1 = pft_ref * (0., 1.) * self%argtransf(:self%pftsz,:)
        pft_ref_tmp2 = pft_ref * (0., 1.) * self%argtransf(self%pftsz + 1:,:)
        if( self%l_cc_bfac )then
            do k = self%kfromto(1), self%kfromto(2)
                sqsumk_ref  = sum(csq(pft_ref(:,k)))
                sqsumk_ptcl = sum(csq(self%pfts_ptcls(self%pinds(iptcl),:,k)))
                denom       = sqrt(sqsumk_ref * sqsumk_ptcl)
                corrk       = self%calc_corrk_for_rot_8(pft_ref, self%pinds(iptcl), self%kfromto(2), k, irot)
                f           = f +       (corrk * real(self%ptcl_bfac_weights(self%pinds(iptcl),k), kind=dp)) / denom
                corrk       = self%calc_corrk_for_rot_8(pft_ref_tmp1, self%pinds(iptcl), self%kfromto(2), k, irot)
                grad(1)     = grad(1) + (corrk * real(self%ptcl_bfac_weights(self%pinds(iptcl),k), kind=dp)) / denom
                corrk       = self%calc_corrk_for_rot_8(pft_ref_tmp2, self%pinds(iptcl), self%kfromto(2), k, irot)
                grad(2)     = grad(2) + (corrk * real(self%ptcl_bfac_weights(self%pinds(iptcl),k), kind=dp)) / denom
            end do
            f       = f       / real(self%ptcl_bfac_norms(self%pinds(iptcl)), kind=dp)
            grad(:) = grad(:) / real(self%ptcl_bfac_norms(self%pinds(iptcl)), kind=dp)
        else
            do k = self%kfromto(1), self%kfromto(2)
                sqsumk_ref  = sum(csq(pft_ref(:,k)))
                sqsumk_ptcl = sum(csq(self%pfts_ptcls(self%pinds(iptcl),:,k)))
                denom       = sqrt(sqsumk_ref * sqsumk_ptcl)
                corrk      = self%calc_corrk_for_rot_8(pft_ref, self%pinds(iptcl), self%kfromto(2), k, irot)
                f          = f + corrk       / denom
                corrk      = self%calc_corrk_for_rot_8(pft_ref_tmp1, self%pinds(iptcl), self%kfromto(2), k, irot)
                grad(1)    = grad(1) + corrk / denom
                corrk      = self%calc_corrk_for_rot_8(pft_ref_tmp2, self%pinds(iptcl), self%kfromto(2), k, irot)
                grad(2)    = grad(2) + corrk / denom
            end do
            f       = f       / real(self%nk, kind=dp)
            grad(:) = grad(:) / real(self%nk, kind=dp)
        endif
    end subroutine gencorr_resnorm_grad_for_rot_8

    !< brief  calculates only gradient for correlation, for one specific rotation angle
    subroutine gencorr_grad_only_for_rot( self, iref, iptcl, shvec, irot, grad )
        use simple_math, only: csq
        class(polarft_corrcalc), intent(inout) :: self
        integer,                 intent(in)    :: iref, iptcl
        real(sp),                intent(in)    :: shvec(2)
        integer,                 intent(in)    :: irot
        real(sp),                intent(out)   :: grad(2)
        if( self%l_cc_objfun )then
            call self%gencorr_cc_grad_only_for_rot( iref, iptcl, shvec, irot, grad )
        else
            call self%gencorr_resnorm_grad_only_for_rot( iref, iptcl, shvec, irot, grad )
        endif
    end subroutine gencorr_grad_only_for_rot

    !< brief  calculates only gradient for correlation, for one specific rotation angle, double precision
    subroutine gencorr_grad_only_for_rot_8( self, iref, iptcl, shvec, irot, grad )
        use simple_math, only: csq
        class(polarft_corrcalc), intent(inout) :: self
        integer,                 intent(in)    :: iref, iptcl
        real(dp),                intent(in)    :: shvec(2)
        integer,                 intent(in)    :: irot
        real(dp),                intent(out)   :: grad(2)
        if( self%l_cc_objfun )then
            call self%gencorr_cc_grad_only_for_rot_8( iref, iptcl, shvec, irot, grad )
        else
            call self%gencorr_resnorm_grad_only_for_rot_8( iref, iptcl, shvec, irot, grad )
        endif
    end subroutine gencorr_grad_only_for_rot_8

    !< brief  calculates only gradient for correlation, for one specific rotation angle
    subroutine gencorr_cc_grad_only_for_rot( self, iref, iptcl, shvec, irot, grad )
        use simple_math, only: csq
        class(polarft_corrcalc), intent(inout) :: self
        integer,                 intent(in)    :: iref, iptcl
        real(sp),                intent(in)    :: shvec(2)
        integer,                 intent(in)    :: irot
        real(sp),                intent(out)   :: grad(2)
        complex(sp), pointer :: pft_ref(:,:), pft_ref_tmp(:,:), shmat(:,:)
        real(sp),    pointer :: argmat(:,:)
        real(sp) :: denom, corr
        integer  :: ithr
        ithr        =  omp_get_thread_num() + 1
        pft_ref     => self%heap_vars(ithr)%pft_ref
        pft_ref_tmp => self%heap_vars(ithr)%pft_ref_tmp
        shmat       => self%heap_vars(ithr)%shmat
        argmat      => self%heap_vars(ithr)%argmat
        argmat      =  self%argtransf(:self%pftsz,:) * shvec(1) + self%argtransf(self%pftsz + 1:,:) * shvec(2)
        shmat       =  cmplx(cos(argmat),sin(argmat))
        if( self%with_ctf )then
            if( self%iseven(self%pinds(iptcl)) )then
                pft_ref = (self%pfts_refs_even(iref,:,:) * self%ctfmats(self%pinds(iptcl),:,:)) * shmat
            else
                pft_ref = (self%pfts_refs_odd(iref,:,:)  * self%ctfmats(self%pinds(iptcl),:,:)) * shmat
            endif
        else
            if( self%iseven(self%pinds(iptcl)) )then
                pft_ref = self%pfts_refs_even(iref,:,:) * shmat
            else
                pft_ref = self%pfts_refs_odd(iref,:,:)  * shmat
            endif
        endif
        denom       = sqrt(sum(csq(pft_ref)) * self%sqsums_ptcls(self%pinds(iptcl)))
        pft_ref_tmp = pft_ref * (0., 1.) * self%argtransf(:self%pftsz,:)
        corr        = self%calc_corr_for_rot(pft_ref_tmp, self%pinds(iptcl), self%kfromto(2), irot)
        grad(1)     = corr / denom
        pft_ref_tmp = pft_ref * (0., 1.) * self%argtransf(self%pftsz + 1:,:)
        corr        = self%calc_corr_for_rot(pft_ref_tmp, self%pinds(iptcl), self%kfromto(2), irot)
        grad(2)     = corr / denom
    end subroutine gencorr_cc_grad_only_for_rot

    !< brief  calculates only gradient for correlation, for one specific rotation angle, double precision
    subroutine gencorr_cc_grad_only_for_rot_8( self, iref, iptcl, shvec, irot, grad )
        use simple_math, only: csq
        class(polarft_corrcalc), intent(inout) :: self
        integer,                 intent(in)    :: iref, iptcl
        real(dp),                intent(in)    :: shvec(2)
        integer,                 intent(in)    :: irot
        real(dp),                intent(out)   :: grad(2)
        complex(dp), pointer :: pft_ref(:,:), pft_ref_tmp(:,:), shmat(:,:)
        real(dp),    pointer :: argmat(:,:)
        real(dp) :: denom, corr
        integer  :: ithr
        ithr        =  omp_get_thread_num() + 1
        pft_ref     => self%heap_vars(ithr)%pft_ref_8
        pft_ref_tmp => self%heap_vars(ithr)%pft_ref_tmp_8
        shmat       => self%heap_vars(ithr)%shmat_8
        argmat      => self%heap_vars(ithr)%argmat_8
        argmat      =  self%argtransf(:self%pftsz,:) * shvec(1) + self%argtransf(self%pftsz + 1:,:) * shvec(2)
        shmat       =  cmplx(cos(argmat),sin(argmat),dp)
        if( self%with_ctf )then
            if( self%iseven(self%pinds(iptcl)) )then
                pft_ref = (self%pfts_refs_even(iref,:,:) * self%ctfmats(self%pinds(iptcl),:,:)) * shmat
            else
                pft_ref = (self%pfts_refs_odd(iref,:,:)  * self%ctfmats(self%pinds(iptcl),:,:)) * shmat
            endif
        else
            if( self%iseven(self%pinds(iptcl)) )then
                pft_ref = self%pfts_refs_even(iref,:,:) * shmat
            else
                pft_ref = self%pfts_refs_odd(iref,:,:)  * shmat
            endif
        endif
        denom       = sqrt(sum(csq(pft_ref)) * self%sqsums_ptcls(self%pinds(iptcl)))
        pft_ref_tmp = pft_ref * (0., 1.) * self%argtransf(:self%pftsz,:)
        corr        = self%calc_corr_for_rot_8(pft_ref_tmp, self%pinds(iptcl), self%kfromto(2), irot)
        grad(1)     = corr / denom
        pft_ref_tmp = pft_ref * (0., 1.) * self%argtransf(self%pftsz + 1:,:)
        corr        = self%calc_corr_for_rot_8(pft_ref_tmp, self%pinds(iptcl), self%kfromto(2), irot)
        grad(2)     = corr / denom
    end subroutine gencorr_cc_grad_only_for_rot_8

    !< brief  calculates only gradient for correlation, for one specific rotation angle
    subroutine gencorr_resnorm_grad_only_for_rot( self, iref, iptcl, shvec, irot, grad )
        use simple_math, only: csq
        class(polarft_corrcalc), intent(inout) :: self
        integer,                 intent(in)    :: iref, iptcl
        real(sp),                intent(in)    :: shvec(2)
        integer,                 intent(in)    :: irot
        real(sp),                intent(out)   :: grad(2)
        complex(sp), pointer :: pft_ref(:,:), pft_ref_tmp1(:,:), pft_ref_tmp2(:,:), shmat(:,:)
        real(sp),    pointer :: argmat(:,:)
        real(sp) :: sqsumk_ref, corrk, sqsumk_ptcl, denom
        integer  :: ithr, k
        ithr         =  omp_get_thread_num() + 1
        pft_ref      => self%heap_vars(ithr)%pft_ref
        pft_ref_tmp1 => self%heap_vars(ithr)%pft_ref_tmp1
        pft_ref_tmp2 => self%heap_vars(ithr)%pft_ref_tmp2
        shmat        => self%heap_vars(ithr)%shmat
        argmat       => self%heap_vars(ithr)%argmat
        argmat       =  self%argtransf(:self%pftsz,:) * shvec(1) + self%argtransf(self%pftsz + 1:,:) * shvec(2)
        shmat        =  cmplx(cos(argmat),sin(argmat))
        if( self%with_ctf )then
            if( self%iseven(self%pinds(iptcl)) )then
                pft_ref = (self%pfts_refs_even(iref,:,:) * self%ctfmats(self%pinds(iptcl),:,:)) * shmat
            else
                pft_ref = (self%pfts_refs_odd(iref,:,:)  * self%ctfmats(self%pinds(iptcl),:,:)) * shmat
            endif
        else
            if( self%iseven(self%pinds(iptcl)) )then
                pft_ref = self%pfts_refs_even(iref,:,:) * shmat
            else
                pft_ref = self%pfts_refs_odd(iref,:,:)  * shmat
            endif
        endif
        pft_ref_tmp1 = pft_ref * (0., 1.) * self%argtransf(:self%pftsz,:)
        pft_ref_tmp2 = pft_ref * (0., 1.) * self%argtransf(self%pftsz + 1:,:)
        grad = 0.0
        if( self%l_cc_bfac )then
            do k = self%kfromto(1), self%kfromto(1)
                sqsumk_ref  = sum(csq(pft_ref(:,k)))
                sqsumk_ptcl = sum(csq(self%pfts_ptcls(self%pinds(iptcl),:,k)))
                denom       = sqrt(sqsumk_ref * sqsumk_ptcl)
                corrk       = self%calc_corrk_for_rot(pft_ref_tmp1, self%pinds(iptcl), self%kfromto(2), k, irot)
                grad(1)     = grad(1) + (corrk * self%ptcl_bfac_weights(self%pinds(iptcl),k)) / denom
                corrk       = self%calc_corrk_for_rot(pft_ref_tmp2, self%pinds(iptcl), self%kfromto(2), k, irot)
                grad(2)     = grad(2) + (corrk * self%ptcl_bfac_weights(self%pinds(iptcl),k)) / denom
            end do
            grad = grad / self%ptcl_bfac_norms(self%pinds(iptcl))
        else
            do k = self%kfromto(1), self%kfromto(1)
                sqsumk_ref  = sum(csq(pft_ref(:,k)))
                sqsumk_ptcl = sum(csq(self%pfts_ptcls(self%pinds(iptcl),:,k)))
                denom       = sqrt(sqsumk_ref * sqsumk_ptcl)
                corrk       = self%calc_corrk_for_rot(pft_ref_tmp1, self%pinds(iptcl), self%kfromto(2), k, irot)
                grad(1)     = grad(1) + corrk / denom
                corrk       = self%calc_corrk_for_rot(pft_ref_tmp2, self%pinds(iptcl), self%kfromto(2), k, irot)
                grad(2)     = grad(2) + corrk / denom
            end do
            grad = grad / real(self%nk)
        endif
    end subroutine gencorr_resnorm_grad_only_for_rot

    !< brief  calculates only gradient for correlation, for one specific rotation angle
    subroutine gencorr_resnorm_grad_only_for_rot_8( self, iref, iptcl, shvec, irot, grad )
        use simple_math, only: csq
        class(polarft_corrcalc), intent(inout) :: self
        integer,                 intent(in)    :: iref, iptcl
        real(dp),                intent(in)    :: shvec(2)
        integer,                 intent(in)    :: irot
        real(dp),                intent(out)   :: grad(2)
        complex(dp), pointer :: pft_ref(:,:), pft_ref_tmp1(:,:), pft_ref_tmp2(:,:), shmat(:,:)
        real(dp),    pointer :: argmat(:,:)
        real(dp) :: sqsumk_ref, corrk, sqsumk_ptcl, denom
        integer  :: ithr, k
        ithr         =  omp_get_thread_num() + 1
        pft_ref      => self%heap_vars(ithr)%pft_ref_8
        pft_ref_tmp1 => self%heap_vars(ithr)%pft_ref_tmp1_8
        pft_ref_tmp2 => self%heap_vars(ithr)%pft_ref_tmp2_8
        shmat        => self%heap_vars(ithr)%shmat_8
        argmat       => self%heap_vars(ithr)%argmat_8
        argmat       =  self%argtransf(:self%pftsz,:) * shvec(1) + self%argtransf(self%pftsz + 1:,:) * shvec(2)
        shmat        =  cmplx(cos(argmat),sin(argmat),dp)
        if( self%with_ctf )then
            if( self%iseven(self%pinds(iptcl)) )then
                pft_ref = (self%pfts_refs_even(iref,:,:) * self%ctfmats(self%pinds(iptcl),:,:)) * shmat
            else
                pft_ref = (self%pfts_refs_odd(iref,:,:)  * self%ctfmats(self%pinds(iptcl),:,:)) * shmat
            endif
        else
            if( self%iseven(self%pinds(iptcl)) )then
                pft_ref = self%pfts_refs_even(iref,:,:) * shmat
            else
                pft_ref = self%pfts_refs_odd(iref,:,:)  * shmat
            endif
        endif
        pft_ref_tmp1 = pft_ref * (0., 1.) * self%argtransf(:self%pftsz,:)
        pft_ref_tmp2 = pft_ref * (0., 1.) * self%argtransf(self%pftsz + 1:,:)
        grad = 0.0_dp
        if( self%l_cc_bfac )then
            do k = self%kfromto(1), self%kfromto(1)
                sqsumk_ref  = sum(csq(pft_ref(:,k)))
                sqsumk_ptcl = sum(csq(self%pfts_ptcls(self%pinds(iptcl),:,k)))
                denom       = sqrt(sqsumk_ref * sqsumk_ptcl)
                corrk       = self%calc_corrk_for_rot_8(pft_ref_tmp1, self%pinds(iptcl), self%kfromto(2), k, irot)
                grad(1)     = grad(1) + (corrk * real(self%ptcl_bfac_weights(self%pinds(iptcl),k), kind=dp)) / denom
                corrk       = self%calc_corrk_for_rot_8(pft_ref_tmp2, self%pinds(iptcl), self%kfromto(2), k, irot)
                grad(2)     = grad(2) + (corrk * real(self%ptcl_bfac_weights(self%pinds(iptcl),k), kind=dp)) / denom
            end do
            grad = grad / real(self%ptcl_bfac_norms(self%pinds(iptcl)), kind=dp)
        else
            do k = self%kfromto(1), self%kfromto(1)
                sqsumk_ref  = sum(csq(pft_ref(:,k)))
                sqsumk_ptcl = sum(csq(self%pfts_ptcls(self%pinds(iptcl),:,k)))
                denom       = sqrt(sqsumk_ref * sqsumk_ptcl)
                corrk       = self%calc_corrk_for_rot_8(pft_ref_tmp1, self%pinds(iptcl), self%kfromto(2), k, irot)
                grad(1)     = grad(1) + corrk / denom
                corrk       = self%calc_corrk_for_rot_8(pft_ref_tmp2, self%pinds(iptcl), self%kfromto(2), k, irot)
                grad(2)     = grad(2) + corrk / denom
            end do
            grad = grad / real(self%nk, kind=dp)
        endif
    end subroutine gencorr_resnorm_grad_only_for_rot_8

    !>  \brief  is for generating resolution dependent correlations
    real function specscore( self, iref, iptcl, irot )
        class(polarft_corrcalc), intent(inout) :: self
        integer,                 intent(in)    :: iref, iptcl, irot
        real :: frc(self%kfromto(1):self%kfromto(2))
        call self%genfrc(iref, self%pinds(iptcl), irot, frc)
        specscore = max(0.,median_nocopy(frc))
    end function specscore

    !>  \brief  is for fitting a bfactor to ptcl vs. ref FRC
    real function fit_bfac( self, iref, iptcl, irot, shvec )
        ! Fitting to Y = A * exp( -B/(4s2) )
        ! from mathworld.wolfram.com/LeastSquaresFittingExponential.html
        class(polarft_corrcalc),  intent(inout) :: self
        integer,                  intent(in)    :: iref, iptcl, irot
        real(sp),                 intent(in)    :: shvec(2)
        real(sp) :: denom, frc(self%kfromto(1):self%kfromto(2)), X(self%kfromto(1):self%kfromto(2))
        real(sp) :: sumxfrc, sumfrc
        logical  :: peakmsk(self%kfromto(1):self%kfromto(2))
        integer  :: i
        call self%calc_frc(iref, iptcl, irot, shvec, frc)
        call peakfinder_inplace(frc, peakmsk)
        where(frc <= TINY) peakmsk = .false.
        if( count(peakmsk) < 3 )then
            where(frc > TINY) peakmsk = .true.
            if( count(peakmsk) < 3 )then
                fit_bfac = 1000.
                return
            endif
        endif
        X = self%smpd * real(self%ldim(1)) / real((/(i,i=self%kfromto(1),self%kfromto(2))/))
        X = -1. / (4.*X*X)
        ! all points are weighted by Y
        sumfrc   = sum(frc,mask=peakmsk)
        sumxfrc  = sum(X*frc,mask=peakmsk)
        denom    = sumfrc * sum(X*X*frc,mask=peakmsk) - sumxfrc**2.
        fit_bfac = sumfrc * sum(X*frc*log(frc), mask=peakmsk)
        fit_bfac = fit_bfac - sumxfrc * sum(frc*log(frc), mask=peakmsk)
        fit_bfac = fit_bfac / denom
    end function fit_bfac

    ! DESTRUCTOR

    !>  \brief  is a destructor
    subroutine kill( self )
        class(polarft_corrcalc), intent(inout) :: self
        integer :: ithr, i, ik
        if( self%existence )then
            do ithr=1,self%nthr
                call fftwf_free(self%fftdat(ithr)%p_ref_re)
                call fftwf_free(self%fftdat(ithr)%p_ref_im)
                call fftwf_free(self%fftdat(ithr)%p_ref_fft_re)
                call fftwf_free(self%fftdat(ithr)%p_ref_fft_im)
                call fftwf_free(self%fftdat(ithr)%p_product_fft)
                call fftwf_free(self%fftdat(ithr)%p_backtransf)
                deallocate(self%heap_vars(ithr)%pft_ref,self%heap_vars(ithr)%pft_ref_tmp,&
                    &self%heap_vars(ithr)%pft_ref_tmp1, self%heap_vars(ithr)%pft_ref_tmp2,&
                    &self%heap_vars(ithr)%corrs_over_k,self%heap_vars(ithr)%argmat,&
                    &self%heap_vars(ithr)%shmat,self%heap_vars(ithr)%kcorrs,&
                    &self%heap_vars(ithr)%pft_ref_8,self%heap_vars(ithr)%pft_ref_tmp_8,&
                    &self%heap_vars(ithr)%pft_ref_tmp1_8, self%heap_vars(ithr)%pft_ref_tmp2_8,&
                    &self%heap_vars(ithr)%shmat_8,self%heap_vars(ithr)%argmat_8)
            end do
            do i = 1, self%nptcls
                do ik = self%kfromto(1),self%kfromto(2)
                    call fftwf_free(self%fftdat_ptcls(i,ik)%p_re)
                    call fftwf_free(self%fftdat_ptcls(i,ik)%p_im)
                end do
            end do
            if( allocated(self%ptcl_bfac_weights) ) deallocate(self%ptcl_bfac_weights)
            if( allocated(self%ptcl_bfac_norms)   ) deallocate(self%ptcl_bfac_norms)
            deallocate( self%sqsums_ptcls, self%angtab, self%argtransf,&
                &self%polar, self%pfts_refs_even, self%pfts_refs_odd, self%pfts_ptcls,&
                &self%fft_factors, self%fftdat, self%fftdat_ptcls,&
                &self%iseven, self%pinds)
            call fftwf_destroy_plan(self%plan_bwd)
            call fftwf_destroy_plan(self%plan_fwd_1)
            call fftwf_destroy_plan(self%plan_fwd_2)
            self%existence = .false.
        endif
    end subroutine kill

end module simple_polarft_corrcalc
