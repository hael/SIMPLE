! for calculation of band-pass limited cross-correlation of polar Fourier transforms
module simple_polarft_corrcalc
!$ use omp_lib
!$ use omp_lib_kinds
include 'simple_lib.f08'
use simple_parameters, only: params_glob
use simple_ori,        only: geodesic_frobdev
implicit none

public :: polarft_corrcalc, pftcc_glob
private
#include "simple_local_flags.inc"

type fftw_cvec
    type(c_ptr)                            :: p
    complex(kind=c_float_complex), pointer :: c(:) => null()
end type fftw_cvec

type fftw_rvec
    type(c_ptr)                 :: p
    real(kind=c_float), pointer :: r(:) => null()
end type fftw_rvec

type fftw_drvec
    type(c_ptr)                  :: p
    real(kind=c_double), pointer :: r(:) => null()
end type fftw_drvec

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
    type(c_ptr)                            :: p_im                      !< -"-
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
    complex(sp), pointer :: pft_ref(:,:)       => null()
    complex(sp), pointer :: pft_ref_tmp(:,:)   => null()
    complex(sp), pointer :: pft_dref(:,:,:)    => null()
    real,        pointer :: corrs_over_k(:)    => null()
    real(dp),    pointer :: argvec(:)          => null()
    complex(sp), pointer :: shmat(:,:)         => null()
    real(dp),    pointer :: kcorrs(:)          => null()
    complex(dp), pointer :: pft_ref_8(:,:)     => null()
    complex(dp), pointer :: pft_ref_tmp_8(:,:) => null()
    complex(dp), pointer :: pft_dref_8(:,:,:)  => null()
    complex(dp), pointer :: shvec(:)           => null()
    complex(dp), pointer :: shmat_8(:,:)       => null()
    real(dp),    pointer :: argmat_8(:,:)      => null()
    real(dp),    pointer :: fdf_y_8(:)         => null()
    real(dp),    pointer :: fdf_T1_8(:,:)      => null()
    real(dp),    pointer :: fdf_T2_8(:,:)      => null()
end type heap_vars

type :: polarft_corrcalc
    ! private
    integer                          :: nptcls     = 1              !< the total number of particles in partition (logically indexded [fromp,top])
    integer                          :: nrefs      = 1              !< the number of references (logically indexded [1,nrefs])
    integer                          :: nrots      = 0              !< number of in-plane rotations for one pft (determined by radius of molecule)
    integer                          :: pftsz      = 0              !< size of reference and particle pft (nrots/2)
    integer                          :: pfromto(2) = 0              !< particle index range
    integer                          :: ldim(3)    = 0              !< logical dimensions of original cartesian image
    integer                          :: kfromto(2)                  !< band-pass Fourier index limits
    integer                          :: nk                          !< number of shells used durring alignement
    integer,             allocatable :: pinds(:)                    !< index array (to reduce memory when frac_update < 1)
    real                             :: delta                       !< voxel size in the frequency domain
    real,                allocatable :: npix_per_shell(:)           !< number of (cartesian) pixels per shell
    real(dp),            allocatable :: sqsums_ptcls(:)             !< memoized square sums for the correlation calculations (taken from kfromto(1):kfromto(2))
    real(dp),            allocatable :: wsqsums_ptcls(:)            !< memoized square sums weighted by sigmas^2 (taken from kfromto(1):kfromto(2))
    real(sp),            allocatable :: angtab(:)                   !< table of in-plane angles (in degrees)
    real(dp),            allocatable :: argtransf(:,:)              !< argument transfer constants for shifting the references
    real(sp),            allocatable :: polar(:,:)                  !< table of polar coordinates (in Cartesian coordinates)
    real(sp),            allocatable :: ctfmats(:,:,:)              !< expand set of CTF matrices (for efficient parallel exec)
    real(dp),            allocatable :: argtransf_shellone(:)       !< one dimensional argument transfer constants (shell k=1) for shifting the references
    real(dp),            allocatable :: refs_reg(:,:,:)             !< -"-, reference reg terms
    real(dp),            allocatable :: regs_denom(:,:,:)           !< -"-
    complex(sp),         allocatable :: pfts_refs_even(:,:,:)       !< 3D complex matrix of polar reference sections (nrefs,pftsz,nk), even
    complex(sp),         allocatable :: pfts_refs_odd(:,:,:)        !< -"-, odd
    complex(sp),         allocatable :: pfts_drefs_even(:,:,:,:)    !< derivatives w.r.t. orientation angles of 3D complex matrices
    complex(sp),         allocatable :: pfts_drefs_odd(:,:,:,:)     !< derivatives w.r.t. orientation angles of 3D complex matrices
    complex(sp),         allocatable :: pfts_ptcls(:,:,:)           !< 3D complex matrix of particle sections
    complex(sp),         allocatable :: fft_factors(:)              !< phase factors for accelerated gencorrs routines
    type(fftw_arrs),     allocatable :: fftdat(:)                   !< arrays for accelerated gencorrs routines
    type(fftw_carr_fft), allocatable :: fftdat_ptcls(:,:)           !< for memoization of particle  FFTs in accelerated gencorrs routines
    type(fftw_carr),     allocatable :: fft_carray(:)               !< for on-the-fly memoization of particle  FFTs
    logical,             allocatable :: iseven(:)                   !< eo assignment for gold-standard FSC
    real,                pointer     :: sigma2_noise(:,:) => null() !< for euclidean distances
    type(c_ptr)                      :: plan_fwd_1                  !< FFTW plans for gencorrs
    type(c_ptr)                      :: plan_fwd_2                  !< -"-
    type(c_ptr)                      :: plan_bwd                    !< -"-
    logical                          :: l_filt_set   = .false.      !< to indicate whether filter is set
    logical                          :: with_ctf     = .false.      !< CTF flag
    logical                          :: existence    = .false.      !< to indicate existence
    type(heap_vars),     allocatable :: heap_vars(:)                !< allocated fields to save stack allocation in subroutines and functions
    ! for testing
    type(c_ptr)                   :: plan_fwd1, plan_bwd1
    type(fftw_cvec),  allocatable :: ft_ptcl_ctf(:,:),  ft_ctf2(:,:)
    type(fftw_cvec),  allocatable :: ft_ref_even(:,:),  ft_ref_odd(:,:)
    type(fftw_cvec),  allocatable :: ft_ref2_even(:,:), ft_ref2_odd(:,:)
    type(fftw_cvec),  allocatable :: cvec1(:), cvec2(:)
    type(fftw_rvec),  allocatable :: rvec1(:)
    type(fftw_drvec), allocatable :: drvec(:)

    contains
    ! CONSTRUCTOR
    procedure          :: new
    ! SETTERS
    procedure          :: reallocate_ptcls
    procedure          :: set_ref_pft
    procedure          :: set_ptcl_pft
    procedure          :: set_ref_fcomp
    procedure          :: set_dref_fcomp
    procedure          :: set_ptcl_fcomp
    procedure          :: cp_even2odd_ref
    procedure          :: cp_odd2even_ref
    procedure          :: cp_even_ref2ptcl
    procedure          :: cp_refs
    procedure          :: swap_ptclsevenodd
    procedure          :: set_eo
    procedure          :: set_eos
    procedure          :: assign_sigma2_noise
    procedure          :: update_sigma
    ! GETTERS
    procedure          :: get_nrots
    procedure          :: get_pdim
    procedure          :: get_pftsz
    procedure          :: get_rot
    procedure          :: get_roind
    procedure          :: get_coord
    procedure          :: get_ref_pft
    procedure          :: get_nrefs
    procedure          :: exists
    procedure          :: ptcl_iseven
    procedure          :: get_nptcls
    procedure          :: assign_pinds
    procedure          :: get_npix
    procedure          :: get_work_pft_ptr
    ! PRINTERS/VISUALISERS
    procedure          :: print
    procedure          :: vis_ptcl
    procedure          :: vis_ref
    procedure          :: polar2cartesian
    ! MODIFIERS
    procedure          :: shift_ptcl
    ! MEMOIZER
    procedure          :: memoize_sqsum_ptcl
    procedure, private :: memoize_fft
    procedure          :: ref_reg_cc
    procedure          :: ref_reg_cc_neigh
    procedure          :: ref_reg_cc_dev
    procedure          :: regularize_refs
    procedure          :: reset_regs
    procedure          :: memoize_ffts
    procedure, private :: setup_npix_per_shell
    ! CALCULATORS
    procedure          :: create_polar_absctfmats, calc_polar_ctf
    procedure, private :: prep_ref4corr_sp
    procedure, private :: prep_ref4corr_dp
    generic            :: prep_ref4corr => prep_ref4corr_sp, prep_ref4corr_dp
    procedure, private :: gen_shmat
    procedure, private :: gen_shmat_8
    procedure          :: calc_corr_rot_shift
    procedure, private :: calc_k_corrs
    procedure, private :: calc_corr_for_rot
    procedure, private :: calc_corr_for_rot_8
    procedure, private :: calc_T1_T2_for_rot_8
    procedure, private :: calc_euclid_for_rot
    procedure, private :: calc_euclid_for_rot_8
    procedure, private :: calc_prob_for_rot
    procedure, private :: calc_prob_for_rot_8
    procedure, private :: calc_corrk_for_rot
    procedure, private :: calc_corrk_for_rot_8
    procedure, private :: calc_euclidk_for_rot
    procedure, private :: calc_euclidk_for_rot_8
    procedure          :: genmaxcorr_comlin
    procedure, private :: gencorrs_cc
    procedure, private :: gencorrs_euclid
    procedure, private :: gencorrs_prob
    procedure, private :: gencorrs_1
    procedure, private :: gencorrs_2
    generic            :: gencorrs => gencorrs_1, gencorrs_2
    procedure          :: gencorr_for_rot_8
    procedure          :: gencorr_grad_for_rot_8
    procedure          :: gencorr_grad_only_for_rot_8
    procedure          :: gencorr_cc_for_rot_8
    procedure          :: gencorr_cont_grad_cc_for_rot_8
    procedure          :: gencorr_cont_cc_for_rot_8
    procedure          :: gencorr_cont_shift_grad_cc_for_rot_8
    procedure          :: gencorr_cc_grad_for_rot_8
    procedure          :: gencorr_cc_grad_only_for_rot_8
    procedure          :: gencorr_euclid_for_rot_8
    procedure          :: gencorr_prob_for_rot_8
    procedure          :: gencorr_cont_grad_euclid_for_rot_8
    procedure          :: gencorr_cont_shift_grad_euclid_for_rot_8
    procedure          :: gencorr_euclid_grad_for_rot_8
    procedure          :: gencorr_euclid_grad_only_for_rot_8
    procedure          :: gencorr_prob_grad_for_rot_8
    procedure          :: gencorr_prob_grad_only_for_rot_8
    procedure          :: gencorr_sigma_contrib
    procedure, private :: calc_frc, calc_raw_frc, calc_pspec
    procedure, private :: rotate_polar_real, rotate_polar_complex
    generic            :: rotate_polar => rotate_polar_real, rotate_polar_complex
    procedure, private :: rotate_ref
    procedure, private :: specscore_1, specscore_2
    generic            :: specscore => specscore_1, specscore_2
    procedure, private :: weight_ref_ptcl_sp, weight_ref_ptcl_dp
    generic,   private :: weight_ref_ptcl => weight_ref_ptcl_sp, weight_ref_ptcl_dp
    procedure, private :: deweight_ref_ptcl_sp, deweight_ref_ptcl_dp
    generic,   private :: deweight_ref_ptcl => deweight_ref_ptcl_sp, deweight_ref_ptcl_dp
    ! DESTRUCTOR
    procedure          :: kill
    ! testing
    procedure          :: memoize_ptcls, memoize_refs, kill_memoized_ptcls, kill_memoized_refs
    procedure          :: allocate_ptcls_memoization, allocate_refs_memoization
    procedure, private :: gencorrs_dev_1, gencorrs_dev_2
    generic            :: gencorrs_dev => gencorrs_dev_1, gencorrs_dev_2
    procedure, private :: gencorrs_cc_dev, gencorrs_euclid_dev, gencorrs_prob_dev
    procedure, private :: gencorrs_shifted_cc_dev, gencorrs_shifted_euclid_dev, gencorrs_shifted_prob_dev
    procedure          :: gencorr_for_rot_8_dev
    procedure, private :: gencorr_cc_for_rot_8_dev, gencorr_euclid_for_rot_8_dev, gencorr_prob_for_rot_8_dev
    procedure          :: gencorr_grad_for_rot_8_dev
    procedure, private :: gencorr_cc_grad_for_rot_8_dev, gencorr_euclid_grad_for_rot_8_dev, gencorr_prob_grad_for_rot_8_dev
    procedure          :: gencorr_grad_only_for_rot_8_dev
    procedure, private :: gencorr_cc_grad_only_for_rot_8_dev
    procedure          :: gencorr_sigma_contrib_dev
    procedure          :: calc_frc_dev
    procedure, private :: gencorr_cont_cc_for_rot_8_dev, gencorr_cont_shift_grad_cc_for_rot_8_dev, gencorr_cont_grad_cc_for_rot_8_dev
end type polarft_corrcalc

! CLASS PARAMETERS/VARIABLES
complex(sp), parameter           :: zero            = cmplx(0.,0.) !< just a complex zero
integer,     parameter           :: FFTW_USE_WISDOM = 16
class(polarft_corrcalc), pointer :: pftcc_glob => null()

contains

    ! CONSTRUCTORS

    subroutine new( self, nrefs, pfromto, kfromto, ptcl_mask, eoarr )
        class(polarft_corrcalc), target, intent(inout) :: self
        integer,                         intent(in)    :: nrefs
        integer,                         intent(in)    :: pfromto(2), kfromto(2)
        logical, optional,               intent(in)    :: ptcl_mask(pfromto(1):pfromto(2))
        integer, optional,               intent(in)    :: eoarr(pfromto(1):pfromto(2))
        character(kind=c_char, len=:), allocatable :: fft_wisdoms_fname ! FFTW wisdoms (per part or suffer I/O lag)
        real(sp), allocatable :: polar_here(:)
        real(dp)              :: A(2)
        real(sp)              :: ang
        integer(kind=c_int)   :: wsdm_ret
        integer               :: local_stat,irot, k, ithr, i, ik, cnt
        logical               :: even_dims, test(2)
        ! kill possibly pre-existing object
        call self%kill
        ! set particle index range
        self%pfromto = pfromto
        ! set band-pass Fourier index limits
        self%kfromto = kfromto
        self%nk      = self%kfromto(2) - self%kfromto(1) + 1
        ! error check
        if( self%pfromto(2) - self%pfromto(1) + 1 < 1 )then
            write(logfhandle,*) 'pfromto: ', self%pfromto(1), self%pfromto(2)
            THROW_HARD ('nptcls (# of particles) must be > 0; new')
        endif
        if( nrefs < 1 )then
            write(logfhandle,*) 'nrefs: ', nrefs
            THROW_HARD ('nrefs (# of reference sections) must be > 0; new')
        endif
        self%ldim = [params_glob%box,params_glob%box,1] !< logical dimensions of original cartesian image
        test      = .false.
        test(1)   = is_even(self%ldim(1))
        test(2)   = is_even(self%ldim(2))
        even_dims = all(test)
        if( .not. even_dims )then
            write(logfhandle,*) 'self%ldim: ', self%ldim
            THROW_HARD ('only even logical dims supported; new')
        endif
        ! set constants
        if( present(ptcl_mask) )then
            self%nptcls  = count(ptcl_mask)                      !< the total number of particles in partition
        else
            self%nptcls  = self%pfromto(2) - self%pfromto(1) + 1 !< the total number of particles in partition
        endif
        self%nrefs = nrefs                                   !< the number of references (logically indexded [1,nrefs])
        self%pftsz = magic_pftsz(nint(params_glob%msk_crop)) !< size of reference (number of vectors used for matching,determined by radius of molecule)
        self%nrots = 2 * self%pftsz                          !< number of in-plane rotations for one pft  (pftsz*2)
        ! generate polar coordinates
        allocate( self%polar(2*self%nrots,self%kfromto(1):self%kfromto(2)),&
                    &self%angtab(self%nrots), self%iseven(1:self%nptcls), polar_here(2*self%nrots))
        ang = twopi/real(self%nrots)
        do irot=1,self%nrots
            self%angtab(irot) = real(irot-1)*ang
            ! cycling over non-redundant logical dimensions
            do k=self%kfromto(1),self%kfromto(2)
                self%polar(irot,k)            =  sin(self%angtab(irot))*real(k) ! x-coordinate
                self%polar(irot+self%nrots,k) = -cos(self%angtab(irot))*real(k) ! y-coordinate
            end do
            ! for k = 1
            polar_here(irot)            =  sin(real(self%angtab(irot)))
            polar_here(irot+self%nrots) = -cos(real(self%angtab(irot)))
            ! angle (in degrees) from now
            self%angtab(irot) = rad2deg(self%angtab(irot))
        end do
        ! index translation table
        allocate( self%pinds(self%pfromto(1):self%pfromto(2)), source=0 )
        if( present(ptcl_mask) )then
            cnt = 0
            do i=self%pfromto(1),self%pfromto(2)
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
                do i=self%pfromto(1),self%pfromto(2)
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
        allocate( self%argtransf(self%nrots,self%kfromto(1):self%kfromto(2)),&
            &self%argtransf_shellone(self%nrots) )
        A = DPI / real(self%ldim(1:2)/2,dp) ! argument transfer matrix normalization constant
        ! shell = 1
        self%argtransf_shellone(:self%pftsz  ) = real(polar_here(:self%pftsz),dp)                        * A(1) ! x-part
        self%argtransf_shellone(self%pftsz+1:) = real(polar_here(self%nrots+1:self%nrots+self%pftsz),dp) * A(2) ! y-part
        ! all shells in resolution range
        self%argtransf(:self%pftsz,:)     = real(self%polar(:self%pftsz,:),dp)                          * A(1)  ! x-part
        self%argtransf(self%pftsz + 1:,:) = real(self%polar(self%nrots + 1:self%nrots+self%pftsz,:),dp) * A(2)  ! y-part
        ! allocate others
        allocate(self%pfts_refs_even(self%pftsz,self%kfromto(1):self%kfromto(2),self%nrefs),&
                    &self%pfts_refs_odd(self%pftsz,self%kfromto(1):self%kfromto(2),self%nrefs),&
                    &self%pfts_drefs_even(self%pftsz,self%kfromto(1):self%kfromto(2),3,params_glob%nthr),&
                    &self%pfts_drefs_odd (self%pftsz,self%kfromto(1):self%kfromto(2),3,params_glob%nthr),&
                    &self%pfts_ptcls(self%pftsz,self%kfromto(1):self%kfromto(2),1:self%nptcls),&
                    &self%sqsums_ptcls(1:self%nptcls),self%wsqsums_ptcls(1:self%nptcls),self%fftdat(params_glob%nthr),self%fft_carray(params_glob%nthr),&
                    &self%fftdat_ptcls(self%kfromto(1):self%kfromto(2),1:self%nptcls),self%regs_denom(self%pftsz,self%kfromto(1):self%kfromto(2),self%nrefs),&
                    &self%heap_vars(params_glob%nthr),self%refs_reg(self%pftsz,self%kfromto(1):self%kfromto(2),self%nrefs))
        local_stat=0
        do ithr=1,params_glob%nthr
            allocate(self%heap_vars(ithr)%pft_ref(self%pftsz,self%kfromto(1):self%kfromto(2)),&
                &self%heap_vars(ithr)%pft_ref_tmp(self%pftsz,self%kfromto(1):self%kfromto(2)),&
                &self%heap_vars(ithr)%pft_dref(self%pftsz,self%kfromto(1):self%kfromto(2),3),&
                &self%heap_vars(ithr)%corrs_over_k(self%nrots),&
                &self%heap_vars(ithr)%argvec(self%pftsz),&
                &self%heap_vars(ithr)%shvec(self%pftsz),&
                &self%heap_vars(ithr)%shmat(self%pftsz,self%kfromto(1):self%kfromto(2)),&
                &self%heap_vars(ithr)%kcorrs(self%nrots),&
                &self%heap_vars(ithr)%pft_ref_8(self%pftsz,self%kfromto(1):self%kfromto(2)),&
                &self%heap_vars(ithr)%pft_ref_tmp_8(self%pftsz,self%kfromto(1):self%kfromto(2)),&
                &self%heap_vars(ithr)%pft_dref_8(self%pftsz,self%kfromto(1):self%kfromto(2),3),&
                &self%heap_vars(ithr)%shmat_8(self%pftsz,self%kfromto(1):self%kfromto(2)),&
                &self%heap_vars(ithr)%argmat_8(self%pftsz,self%kfromto(1):self%kfromto(2)),&
                &self%heap_vars(ithr)%fdf_y_8(self%kfromto(1):self%kfromto(2)),&
                &self%heap_vars(ithr)%fdf_T1_8(self%kfromto(1):self%kfromto(2),3),&
                &self%heap_vars(ithr)%fdf_T2_8(self%kfromto(1):self%kfromto(2),3))
        end do
        self%pfts_refs_even = zero
        self%pfts_refs_odd  = zero
        self%pfts_ptcls     = zero
        self%sqsums_ptcls   = 0.d0
        self%wsqsums_ptcls  = 0.d0
        self%refs_reg       = 0.d0
        self%regs_denom     = 0.d0
        ! set CTF flag
        self%with_ctf = .false.
        if( params_glob%ctf .ne. 'no' ) self%with_ctf = .true.
        ! thread-safe c-style allocatables for gencorrs
        do ithr=1,params_glob%nthr
            self%fftdat(ithr)%p_ref_re       = fftwf_alloc_real   (int(self%pftsz, c_size_t))
            self%fftdat(ithr)%p_ref_im       = fftwf_alloc_complex(int(self%pftsz, c_size_t))
            self%fftdat(ithr)%p_ref_fft_re   = fftwf_alloc_complex(int(self%pftsz, c_size_t))
            self%fftdat(ithr)%p_ref_fft_im   = fftwf_alloc_complex(int(self%pftsz, c_size_t))
            self%fftdat(ithr)%p_product_fft  = fftwf_alloc_complex(int(self%nrots, c_size_t))
            self%fftdat(ithr)%p_backtransf   = fftwf_alloc_real   (int(self%nrots, c_size_t))
            call c_f_pointer(self%fftdat(ithr)%p_ref_re,        self%fftdat(ithr)%ref_re,        [self%pftsz])
            call c_f_pointer(self%fftdat(ithr)%p_ref_im,        self%fftdat(ithr)%ref_im,        [self%pftsz])
            call c_f_pointer(self%fftdat(ithr)%p_ref_fft_re,    self%fftdat(ithr)%ref_fft_re,    [self%pftsz])
            call c_f_pointer(self%fftdat(ithr)%p_ref_fft_im,    self%fftdat(ithr)%ref_fft_im,    [self%pftsz])
            call c_f_pointer(self%fftdat(ithr)%p_product_fft,   self%fftdat(ithr)%product_fft,   [self%nrots])
            call c_f_pointer(self%fftdat(ithr)%p_backtransf,    self%fftdat(ithr)%backtransf,    [self%nrots])
            ! thread-safe c-style allocatables for on-the-fly particle memoization
            self%fft_carray(ithr)%p_re = fftwf_alloc_real(int(self%pftsz, c_size_t))
            self%fft_carray(ithr)%p_im = fftwf_alloc_complex(int(self%pftsz, c_size_t))
            call c_f_pointer(self%fft_carray(ithr)%p_re, self%fft_carray(ithr)%re, [self%pftsz])
            call c_f_pointer(self%fft_carray(ithr)%p_im, self%fft_carray(ithr)%im, [self%pftsz])
        end do
        ! thread-safe c-style allocatables for gencorrs, particle memoization
        do i = 1,self%nptcls
            do ik = self%kfromto(1),self%kfromto(2)
                self%fftdat_ptcls(ik,i)%p_re = fftwf_alloc_complex(int(self%pftsz, c_size_t))
                self%fftdat_ptcls(ik,i)%p_im = fftwf_alloc_complex(int(self%pftsz, c_size_t))
                call c_f_pointer(self%fftdat_ptcls(ik,i)%p_re, self%fftdat_ptcls(ik,i)%re, [self%pftsz])
                call c_f_pointer(self%fftdat_ptcls(ik,i)%p_im, self%fftdat_ptcls(ik,i)%im, [self%pftsz])
            end do
        end do
        ! FFTW3 wisdoms file
        if( params_glob%l_distr_exec )then
            allocate(fft_wisdoms_fname, source='fft_wisdoms_part'//int2str_pad(params_glob%part,params_glob%numlen)//'.dat'//c_null_char)
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
        ! setup npix_per_shell
        call self%setup_npix_per_shell
        ! testing
        if( L_CTFROTDEV ) call self%allocate_ptcls_memoization
        ! flag existence
        self%existence = .true.
        ! set pointer to global instance
        pftcc_glob => self
        ! voxel-size in the frequency domain
        self%delta = 1. / real(params_glob%box) / params_glob%smpd
    end subroutine new

    ! SETTERS

    subroutine reallocate_ptcls( self, nptcls, pinds )
        class(polarft_corrcalc), intent(inout) :: self
        integer,                 intent(in)    :: nptcls
        integer,                 intent(in)    :: pinds(nptcls)
        integer :: i,iptcl,ik
        self%pfromto(1) = minval(pinds)
        self%pfromto(2) = maxval(pinds)
        if( allocated(self%pinds) ) deallocate(self%pinds)
        if( self%nptcls == nptcls )then
            ! just need to update particles indexing
        else
            ! re-index & reallocate
            self%nptcls = nptcls
            if( allocated(self%sqsums_ptcls) ) deallocate(self%sqsums_ptcls)
            if( allocated(self%wsqsums_ptcls)) deallocate(self%wsqsums_ptcls)
            if( allocated(self%iseven) )       deallocate(self%iseven)
            if( allocated(self%pfts_ptcls) )   deallocate(self%pfts_ptcls)
            if( allocated(self%fftdat_ptcls) )then
                do i = 1, size(self%fftdat_ptcls,dim=2)
                    do ik = self%kfromto(1),self%kfromto(2)
                        call fftwf_free(self%fftdat_ptcls(ik,i)%p_re)
                        call fftwf_free(self%fftdat_ptcls(ik,i)%p_im)
                    end do
                end do
                deallocate(self%fftdat_ptcls)
            endif
            allocate( self%pfts_ptcls(self%pftsz,self%kfromto(1):self%kfromto(2),1:self%nptcls),&
                        &self%sqsums_ptcls(1:self%nptcls),self%wsqsums_ptcls(1:self%nptcls),self%iseven(1:self%nptcls),&
                        &self%fftdat_ptcls(self%kfromto(1):self%kfromto(2), 1:self%nptcls) )
            do i = 1,self%nptcls
                do ik = self%kfromto(1),self%kfromto(2)
                    self%fftdat_ptcls(ik,i)%p_re = fftwf_alloc_complex(int(self%pftsz, c_size_t))
                    self%fftdat_ptcls(ik,i)%p_im = fftwf_alloc_complex(int(self%pftsz, c_size_t))
                    call c_f_pointer(self%fftdat_ptcls(ik,i)%p_re, self%fftdat_ptcls(ik,i)%re, [self%pftsz])
                    call c_f_pointer(self%fftdat_ptcls(ik,i)%p_im, self%fftdat_ptcls(ik,i)%im, [self%pftsz])
                end do
            end do
            if( L_CTFROTDEV )then
                call self%kill_memoized_ptcls
                call self%allocate_ptcls_memoization
            endif
        endif
        self%pfts_ptcls    = zero
        self%sqsums_ptcls  = 0.d0
        self%wsqsums_ptcls = 0.d0
        self%iseven        = .true.
        allocate(self%pinds(self%pfromto(1):self%pfromto(2)), source=0)
        do i = 1,self%nptcls
            iptcl = pinds(i)
            self%pinds( iptcl ) = i
        enddo
    end subroutine reallocate_ptcls

    subroutine set_ref_pft( self, iref, pft, iseven )
        class(polarft_corrcalc), intent(inout) :: self     !< this object
        integer,                 intent(in)    :: iref     !< reference index
        complex(sp),             intent(in)    :: pft(:,:) !< reference pft
        logical,                 intent(in)    :: iseven   !< logical eo-flag
        if( iseven )then
            self%pfts_refs_even(:,:,iref) = pft
        else
            self%pfts_refs_odd(:,:,iref)  = pft
        endif
    end subroutine set_ref_pft

    subroutine set_ptcl_pft( self, iptcl, pft )
        class(polarft_corrcalc), intent(inout) :: self     !< this object
        integer,                 intent(in)    :: iptcl    !< particle index
        complex(sp),             intent(in)    :: pft(:,:) !< particle's pft
        self%pfts_ptcls(:,:,self%pinds(iptcl)) = pft
        call self%memoize_sqsum_ptcl(iptcl)
    end subroutine set_ptcl_pft

    subroutine set_ref_fcomp( self, iref, irot, k, comp, iseven )
        class(polarft_corrcalc), intent(inout) :: self
        integer,                 intent(in)    :: iref, irot, k
        complex(sp),             intent(in)    :: comp
        logical,                 intent(in)    :: iseven
        if( iseven )then
            self%pfts_refs_even(irot,k,iref) = comp
        else
            self%pfts_refs_odd(irot,k,iref)  = comp
        endif
    end subroutine set_ref_fcomp

    subroutine set_dref_fcomp( self, iref, irot, k, dcomp, iseven )
        class(polarft_corrcalc), intent(inout) :: self
        integer,                 intent(in)    :: iref, irot, k
        complex(sp),             intent(in)    :: dcomp(3)
        logical,                 intent(in)    :: iseven
        if( iseven )then
            self%pfts_drefs_even(irot,k,:,iref) = dcomp
        else
            self%pfts_drefs_odd(irot,k,:,iref)  = dcomp
        endif
    end subroutine set_dref_fcomp

    subroutine set_ptcl_fcomp( self, iptcl, irot, k, comp )
        class(polarft_corrcalc), intent(inout) :: self
        integer,                 intent(in)    :: iptcl, irot, k
        complex(sp),             intent(in)    :: comp
        self%pfts_ptcls(irot,k,self%pinds(iptcl)) = comp
    end subroutine set_ptcl_fcomp

    subroutine cp_even2odd_ref( self, iref )
        class(polarft_corrcalc), intent(inout) :: self
        integer,                 intent(in)    :: iref
        self%pfts_refs_odd(:,:,iref) = self%pfts_refs_even(:,:,iref)
    end subroutine cp_even2odd_ref

    subroutine cp_odd2even_ref( self, iref )
        class(polarft_corrcalc), intent(inout) :: self
        integer,                 intent(in)    :: iref
        self%pfts_refs_even(:,:,iref) = self%pfts_refs_odd(:,:,iref)
    end subroutine cp_odd2even_ref

    subroutine cp_refs( self, self2 )
        class(polarft_corrcalc), intent(inout) :: self, self2
        self%pfts_refs_odd  = self2%pfts_refs_odd
        self%pfts_refs_even = self2%pfts_refs_even
    end subroutine cp_refs

    subroutine cp_even_ref2ptcl( self, iref, iptcl )
        class(polarft_corrcalc), intent(inout) :: self
        integer,                 intent(in)    :: iref, iptcl
        self%pfts_ptcls(:,:,self%pinds(iptcl)) = self%pfts_refs_even(:,:,iref)
        call self%memoize_sqsum_ptcl(self%pinds(iptcl))
    end subroutine cp_even_ref2ptcl

    subroutine swap_ptclsevenodd( self )
        class(polarft_corrcalc), intent(inout) :: self
        self%iseven = .not.self%iseven
    end subroutine swap_ptclsevenodd

    subroutine set_eo( self, iptcl, is_even )
        class(polarft_corrcalc), intent(inout) :: self
        integer,                 intent(in)    :: iptcl
        logical,                 intent(in)    :: is_even
        self%iseven(self%pinds(iptcl)) = is_even
    end subroutine set_eo

    subroutine set_eos( self, eoarr )
        class(polarft_corrcalc), intent(inout) :: self
        integer,                 intent(in)    :: eoarr(self%nptcls)
        integer :: i
        if( all(eoarr == - 1) )then
            self%iseven = .true.
        else
            do i=1,self%nptcls
                if( eoarr(i) == 0 )then
                    self%iseven(i) = .true.
                else
                    self%iseven(i) = .false.
                endif
            end do
        endif
    end subroutine set_eos

    subroutine assign_sigma2_noise( self, sigma2_noise )
        class(polarft_corrcalc),      intent(inout) :: self
        real,    allocatable, target, intent(inout) :: sigma2_noise(:,:)
        self%sigma2_noise => sigma2_noise
    end subroutine assign_sigma2_noise

    ! GETTERS

    !>  \brief  for getting the number of in-plane rotations
    pure function get_nrots( self ) result( nrots )
        class(polarft_corrcalc), intent(in) :: self
        integer :: nrots
        nrots = self%nrots
    end function get_nrots

    !>  \brief  for getting the dimensions of the reference polar FT
    pure function get_pdim( self ) result( pdim )
        class(polarft_corrcalc), intent(in) :: self
        integer :: pdim(3)
        pdim = [self%pftsz,self%kfromto(1),self%kfromto(2)]
    end function get_pdim

    ! !>  \brief  for getting the dimension of the reference polar FT
    pure integer function get_pftsz( self )
        class(polarft_corrcalc), intent(in) :: self
        get_pftsz = self%pftsz
    end function get_pftsz

    !>  \brief is for getting the continuous in-plane rotation
    !!         corresponding to in-plane rotation index roind
    function get_rot( self, roind ) result( rot )
        class(polarft_corrcalc), intent(in) :: self
        integer,                 intent(in) :: roind !< in-plane rotation index
        real(sp) :: rot
        if( roind < 1 .or. roind > self%nrots )then
            write(logfhandle,*) 'roind: ', roind
            write(logfhandle,*) 'nrots: ', self%nrots
            THROW_HARD('roind is out of range; get_rot')
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

    !>  \brief returns polar coordinate for rotation rot
    !!         and Fourier index k
    function get_coord( self, rot, k ) result( xy )
        class(polarft_corrcalc), intent(in) :: self
        integer,                 intent(in) :: rot, k
        real(sp) :: xy(2)
        xy(1) = self%polar(rot,k)
        xy(2) = self%polar(self%nrots+rot,k)
    end function get_coord

    !>  \brief  returns polar Fourier transform of reference iref
    function get_ref_pft( self, iref, iseven ) result( pft )
        class(polarft_corrcalc), intent(in) :: self
        integer,                 intent(in) :: iref
        logical,                 intent(in) :: iseven
        complex(sp), allocatable :: pft(:,:)
        if( iseven )then
            allocate(pft(self%pftsz,self%kfromto(1):self%kfromto(2)),&
            source=self%pfts_refs_even(:,:,iref))
        else
            allocate(pft(self%pftsz,self%kfromto(1):self%kfromto(2)),&
            source=self%pfts_refs_odd(:,:,iref))
        endif
    end function get_ref_pft

    integer function get_nrefs( self )
        class(polarft_corrcalc), intent(in) :: self
        get_nrefs = self%nrefs
    end function get_nrefs

    logical function exists( self )
        class(polarft_corrcalc), intent(in) :: self
        exists = self%existence
    end function exists

    logical function ptcl_iseven( self, iptcl )
        class(polarft_corrcalc), intent(in) :: self
        integer,                 intent(in) :: iptcl
        ptcl_iseven = self%iseven(self%pinds(iptcl))
    end function ptcl_iseven

    integer function get_nptcls( self )
        class(polarft_corrcalc), intent(in) :: self
        get_nptcls = self%nptcls
    end function get_nptcls

    subroutine assign_pinds( self, pinds )
        class(polarft_corrcalc), intent(inout) :: self
        integer, allocatable,    intent(out)   :: pinds(:)
        pinds = self%pinds
    end subroutine assign_pinds

    integer function get_npix( self )
        class(polarft_corrcalc), intent(in) :: self
        get_npix = sum(nint(self%npix_per_shell(self%kfromto(1):self%kfromto(2))))
    end function get_npix

    ! returns pointer to temporary pft according to current thread
    subroutine get_work_pft_ptr( self, ptr )
        class(polarft_corrcalc), intent(in) :: self
        complex(sp),   pointer, intent(out) :: ptr(:,:)
        integer :: ithr
        ithr = omp_get_thread_num()+1
        ptr => self%heap_vars(ithr)%pft_ref_tmp
    end subroutine get_work_pft_ptr

    ! PRINTERS/VISUALISERS

    subroutine vis_ptcl( self, iptcl )
        use gnufor2
        class(polarft_corrcalc), intent(in) :: self
        integer,                 intent(in) :: iptcl
        call gnufor_image( real(self%pfts_ptcls(:,:,self%pinds(iptcl))), palette='gray')
        call gnufor_image(aimag(self%pfts_ptcls(:,:,self%pinds(iptcl))), palette='gray')
    end subroutine vis_ptcl

    subroutine vis_ref( self, iref, iseven )
        use gnufor2
        class(polarft_corrcalc), intent(in) :: self
        integer,                 intent(in) :: iref
        logical,                 intent(in) :: iseven
        if( iseven )then
            call gnufor_image( real(self%pfts_refs_even(:,:,iref)), palette='gray')
            call gnufor_image(aimag(self%pfts_refs_even(:,:,iref)), palette='gray')
        else
            call gnufor_image( real(self%pfts_refs_odd(:,:,iref)), palette='gray')
            call gnufor_image(aimag(self%pfts_refs_odd(:,:,iref)), palette='gray')
        endif
    end subroutine vis_ref

    subroutine polar2cartesian( self, i, isref, cmat, box )
        class(polarft_corrcalc), intent(in)    :: self
        integer,                 intent(in)    :: i
        logical,                 intent(in)    :: isref
        complex,    allocatable, intent(inout) :: cmat(:,:)
        integer,                 intent(out)   :: box
        integer, allocatable :: norm(:,:)
        complex :: comp
        integer :: k,c,irot,physh,physk
        if( allocated(cmat) ) deallocate(cmat)
        box = 2*self%kfromto(2)
        c   = box/2+1
        allocate(cmat(box/2+1,box),source=cmplx(0.0,0.0))
        allocate(norm(box/2+1,box),source=0)
        do irot=1,self%pftsz
            do k=self%kfromto(1),self%kfromto(2)
                ! Nearest-neighbour interpolation
                physh = nint(self%polar(irot,k)) + 1
                physk = nint(self%polar(irot+self%nrots,k)) + c
                if( physk > box ) cycle
                if( isref )then
                    comp = self%pfts_refs_even(irot,k,i)
                else
                    comp = self%pfts_ptcls(irot,k,i)
                endif
                cmat(physh,physk) = cmat(physh,physk) + comp
                norm(physh,physk) = norm(physh,physk) + 1
            end do
        end do
        ! normalization
        where(norm>0)
            cmat = cmat / real(norm)
        end where
        ! irot = self%pftsz+1, eg. angle=180.
        do k = 1,box/2-1
            cmat(1,k+c) = conjg(cmat(1,c-k))
        enddo
        ! arbitrary magnitude
        cmat(1,c) = (100.0,0.0)
    end subroutine polar2cartesian

    subroutine print( self )
        class(polarft_corrcalc), intent(in) :: self
        write(logfhandle,*) "total n particles in partition         (self%nptcls): ", self%nptcls
        write(logfhandle,*) "number of references                    (self%nrefs): ", self%nrefs
        write(logfhandle,*) "number of rotations                     (self%nrots): ", self%nrots
        write(logfhandle,*) "size of pft                             (self%pftsz): ", self%pftsz
        write(logfhandle,*) "logical dim. of original Cartesian image (self%ldim): ", self%ldim
    end subroutine print

    ! MODIFIERS

    subroutine shift_ptcl( self, iptcl, shvec)
        class(polarft_corrcalc),  intent(inout) :: self
        integer,                  intent(in)    :: iptcl
        real(sp),                 intent(in)    :: shvec(2)
        complex(sp), pointer :: shmat(:,:)
        integer  :: ithr, i
        ithr  = omp_get_thread_num() + 1
        i     = self%pinds(iptcl)
        shmat => self%heap_vars(ithr)%shmat
        call self%gen_shmat(ithr, shvec, shmat)
        self%pfts_ptcls(:,:,i) = self%pfts_ptcls(:,:,i) * shmat
    end subroutine shift_ptcl

    ! MEMOIZERS

    subroutine memoize_sqsum_ptcl( self, iptcl )
        class(polarft_corrcalc), intent(inout) :: self
        integer,                 intent(in)    :: iptcl
        real(dp) :: sumsqk
        integer  :: i, ik
        logical  :: l_sigma
        i       = self%pinds(iptcl)
        l_sigma = associated(self%sigma2_noise)
        self%sqsums_ptcls(i) = 0.d0
        if( l_sigma ) self%wsqsums_ptcls(i) = 0.d0
        do ik = self%kfromto(1),self%kfromto(2)
            sumsqk               = real(ik,dp) * real(sum(csq_fast(self%pfts_ptcls(:,ik,i))),dp)
            self%sqsums_ptcls(i) = self%sqsums_ptcls(i) + sumsqk
            if( l_sigma ) self%wsqsums_ptcls(i) = self%wsqsums_ptcls(i) + sumsqk / real(self%sigma2_noise(ik,iptcl),dp)
        enddo
    end subroutine memoize_sqsum_ptcl

    ! memoize all particles ffts
    subroutine memoize_ffts( self )
        class(polarft_corrcalc), intent(inout) :: self
        integer :: i
        ! memoize particle FFTs in parallel
        !$omp parallel do default(shared) private(i) proc_bind(close) schedule(static)
        do i=1,self%nptcls
            call self%memoize_fft(i)
        end do
        !$omp end parallel do
    end subroutine memoize_ffts

    ! memoize particle fft, serial only, private use
    subroutine memoize_fft( self, i )
        class(polarft_corrcalc), intent(inout) :: self
        integer  :: i, ik, ithr
        ithr = omp_get_thread_num() + 1
        ! memoize particle FFTs
        do ik = self%kfromto(1),self%kfromto(2)
            ! copy particle pfts
            self%fft_carray(ithr)%re = real(self%pfts_ptcls(:,ik,i))
            self%fft_carray(ithr)%im = aimag(self%pfts_ptcls(:,ik,i)) * self%fft_factors
            ! FFT
            call fftwf_execute_dft_r2c(self%plan_fwd_1, self%fft_carray(ithr)%re, self%fftdat_ptcls(ik,i)%re)
            call fftwf_execute_dft    (self%plan_fwd_2, self%fft_carray(ithr)%im, self%fftdat_ptcls(ik,i)%im)
        end do
    end subroutine memoize_fft

    ! accumulating reference reg terms for each batch of particles, with cc-based global objfunc
    subroutine ref_reg_cc( self, eulspace, ptcl_eulspace, glob_pinds )
        use simple_oris
        class(polarft_corrcalc), intent(inout) :: self
        type(oris),              intent(in)    :: eulspace
        type(oris),              intent(in)    :: ptcl_eulspace
        integer,                 intent(in)    :: glob_pinds(self%nptcls)
        integer  :: i, iref, iptcl, loc
        real     :: inpl_corrs(self%nrots), ptcl_ref_dist, ptcl_ctf(self%pftsz,self%kfromto(1):self%kfromto(2),self%nptcls)
        real     :: euls(3), euls_ref(3)
        real(dp) :: ctf_rot(self%pftsz,self%kfromto(1):self%kfromto(2)), ptcl_ctf_rot(self%pftsz,self%kfromto(1):self%kfromto(2))
        ptcl_ctf = real(self%pfts_ptcls * self%ctfmats)
        !$omp parallel do collapse(2) default(shared) private(i,iref,euls_ref,euls,ptcl_ref_dist,iptcl,inpl_corrs,loc,ptcl_ctf_rot, ctf_rot) proc_bind(close) schedule(static)
        do iref = 1, self%nrefs
            do i = 1, self%nptcls
                iptcl    = glob_pinds(i)
                euls_ref = eulspace%get_euler(iref)
                euls     = ptcl_eulspace%get_euler(iptcl)
                ! projection direction distance, euler_dist could be used instead
                euls_ref(3)   = 0.
                euls(3)       = 0.
                ptcl_ref_dist = geodesic_frobdev(euls_ref,euls)
                ! find best irot for this pair of iref, iptcl
                call self%gencorrs( iref, iptcl, inpl_corrs )
                loc = maxloc(inpl_corrs, dim=1)
                if( inpl_corrs(loc) < TINY ) cycle
                ! distance & correlation weighing
                ptcl_ref_dist = inpl_corrs(loc) / ( 1. + ptcl_ref_dist )
                ! computing the reg terms as the gradients w.r.t 2D references of the probability
                loc = (self%nrots+1)-(loc-1)
                if( loc > self%nrots ) loc = loc - self%nrots
                call self%rotate_polar(    ptcl_ctf(:,:,i), ptcl_ctf_rot, loc)
                call self%rotate_polar(self%ctfmats(:,:,i),      ctf_rot, loc)
                self%refs_reg(  :,:,iref) = self%refs_reg(  :,:,iref) + ptcl_ctf_rot * real(ptcl_ref_dist, dp)
                self%regs_denom(:,:,iref) = self%regs_denom(:,:,iref) +   ctf_rot**2 * real(ptcl_ref_dist, dp)
            enddo
        enddo
        !$omp end parallel do
    end subroutine ref_reg_cc

    subroutine ref_reg_cc_neigh( self, eulspace, ptcl_eulspace, glob_pinds )
        use simple_oris
        class(polarft_corrcalc), intent(inout) :: self
        type(oris),              intent(in)    :: eulspace
        type(oris),              intent(in)    :: ptcl_eulspace
        integer,                 intent(in)    :: glob_pinds(self%nptcls)
        integer  :: i, iref, iptcl, loc
        real     :: inpl_corrs(self%nrots), ptcl_ref_dist, ptcl_ctf(self%pftsz,self%kfromto(1):self%kfromto(2),self%nptcls)
        real     :: euls(3), euls_ref(3), theta
        real(dp) :: ctf_rot(self%pftsz,self%kfromto(1):self%kfromto(2)), ptcl_ctf_rot(self%pftsz,self%kfromto(1):self%kfromto(2))
        ptcl_ctf = real(self%pfts_ptcls * self%ctfmats)
        !$omp parallel do collapse(2) default(shared) private(i,iref,euls_ref,euls,ptcl_ref_dist,iptcl,inpl_corrs,loc,ptcl_ctf_rot,ctf_rot,theta) proc_bind(close) schedule(static)
        do iref = 1, self%nrefs
            do i = 1, self%nptcls
                iptcl    = glob_pinds(i)
                euls_ref =      eulspace%get_euler(iref)  * pi / 180.
                euls     = ptcl_eulspace%get_euler(iptcl) * pi / 180.
                theta    = acos(cos(euls_ref(2))*cos(euls(2)) + sin(euls_ref(2))*sin(euls(2))*cos(euls_ref(1) - euls(1)))
                if( theta <= params_glob%arc_thres*pi/180.  .and. theta >= 0. )then
                    ! find best irot for this pair of iref, iptcl
                    call self%gencorrs( iref, iptcl, inpl_corrs )
                    loc = maxloc(inpl_corrs, dim=1)
                    if( inpl_corrs(loc) < TINY ) cycle
                    ! distance & correlation weighing
                    ptcl_ref_dist = inpl_corrs(loc) / ( 1. + theta )
                    ! computing the reg terms as the gradients w.r.t 2D references of the probability
                    loc = (self%nrots+1)-(loc-1)
                    if( loc > self%nrots ) loc = loc - self%nrots
                    call self%rotate_polar(    ptcl_ctf(:,:,i), ptcl_ctf_rot, loc)
                    call self%rotate_polar(self%ctfmats(:,:,i),      ctf_rot, loc)
                    self%refs_reg(  :,:,iref) = self%refs_reg(  :,:,iref) + ptcl_ctf_rot * real(ptcl_ref_dist, dp)
                    self%regs_denom(:,:,iref) = self%regs_denom(:,:,iref) +   ctf_rot**2 * real(ptcl_ref_dist, dp)
                endif
            enddo
        enddo
        !$omp end parallel do
    end subroutine ref_reg_cc_neigh

    subroutine ref_reg_cc_dev( self, eulspace, ptcl_eulspace, glob_pinds )
        use simple_oris
        class(polarft_corrcalc), intent(inout) :: self
        type(oris),              intent(in)    :: eulspace
        type(oris),              intent(in)    :: ptcl_eulspace
        integer,                 intent(in)    :: glob_pinds(self%nptcls)
        integer  :: i, iref, iptcl, loc, k, k_int
        real     :: inpl_corrs(self%nrots), ptcl_ref_dist, ptcl_ctf(self%pftsz,self%kfromto(1):self%kfromto(2),self%nptcls)
        real     :: euls(3), euls_ref(3), cos_theta, k_proj, theta
        real(dp) :: ctf_rot(self%pftsz,self%kfromto(1):self%kfromto(2)), ptcl_ctf_rot(self%pftsz,self%kfromto(1):self%kfromto(2))
        ptcl_ctf = real(self%pfts_ptcls * self%ctfmats)
        !$omp parallel do collapse(2) default(shared) private(i,iref,euls_ref,euls,ptcl_ref_dist,iptcl,inpl_corrs,loc,ptcl_ctf_rot,ctf_rot,cos_theta,theta,k_proj,k_int,k) proc_bind(close) schedule(static)
        do iref = 1, self%nrefs
            do i = 1, self%nptcls
                iptcl     = glob_pinds(i)
                euls_ref  =      eulspace%get_euler(iref)  * pi / 180.
                euls      = ptcl_eulspace%get_euler(iptcl) * pi / 180.
                cos_theta = cos(euls_ref(2))*cos(euls(2)) + sin(euls_ref(2))*sin(euls(2))*cos(euls_ref(1) - euls(1))
                theta     = acos(cos_theta)
                if( theta < params_glob%arc_thres*pi/180. .and. theta >= 0. .and. &
                 &(self%kfromto(2)*cos_theta >= self%kfromto(1) .or. self%kfromto(1)*cos_theta <= self%kfromto(2)) )then
                    ! find best irot for this pair of iref, iptcl
                    call self%gencorrs( iref, iptcl, inpl_corrs )
                    loc = maxloc(inpl_corrs, dim=1)
                    if( inpl_corrs(loc) < TINY ) cycle
                    ! distance & correlation weighing
                    ptcl_ref_dist = 1.
                    ! computing the reg terms as the gradients w.r.t 2D references of the probability
                    loc = (self%nrots+1)-(loc-1)
                    if( loc > self%nrots ) loc = loc - self%nrots
                    call self%rotate_polar(    ptcl_ctf(:,:,i), ptcl_ctf_rot, loc)
                    call self%rotate_polar(self%ctfmats(:,:,i),      ctf_rot, loc)
                    do k = self%kfromto(1),self%kfromto(2)
                        k_proj = real(k) * cos_theta
                        if( k_proj > self%kfromto(1) .and. k_proj < self%kfromto(2) )then
                            k_int = floor(k_proj)
                            self%refs_reg(  :,k_int,iref) = self%refs_reg(  :,k_int,iref) + ptcl_ctf_rot(:,k)    * real(ptcl_ref_dist, dp) * (k_int + 1. - k_proj)
                            self%regs_denom(:,k_int,iref) = self%regs_denom(:,k_int,iref) +      ctf_rot(:,k)**2 * real(ptcl_ref_dist, dp) * (k_int + 1. - k_proj)
                            k_int = ceiling(k_proj)
                            self%refs_reg(  :,k_int,iref) = self%refs_reg(  :,k_int,iref) + ptcl_ctf_rot(:,k)    * real(ptcl_ref_dist, dp) * (k_proj - k_int + 1)
                            self%regs_denom(:,k_int,iref) = self%regs_denom(:,k_int,iref) +      ctf_rot(:,k)**2 * real(ptcl_ref_dist, dp) * (k_proj - k_int + 1)
                        elseif( int(k_proj) == self%kfromto(1) )then
                            k_int = self%kfromto(1)
                            self%refs_reg(  :,k_int,iref) = self%refs_reg(  :,k_int,iref) + ptcl_ctf_rot(:,k)    * real(ptcl_ref_dist, dp)
                            self%regs_denom(:,k_int,iref) = self%regs_denom(:,k_int,iref) +      ctf_rot(:,k)**2 * real(ptcl_ref_dist, dp)
                        elseif( int(k_proj) == self%kfromto(2) )then
                            k_int = self%kfromto(2)
                            self%refs_reg(  :,k_int,iref) = self%refs_reg(  :,k_int,iref) + ptcl_ctf_rot(:,k)    * real(ptcl_ref_dist, dp)
                            self%regs_denom(:,k_int,iref) = self%regs_denom(:,k_int,iref) +      ctf_rot(:,k)**2 * real(ptcl_ref_dist, dp)
                        endif
                    enddo
                endif
            enddo
        enddo
        !$omp end parallel do
    end subroutine ref_reg_cc_dev

    subroutine regularize_refs( self )
        class(polarft_corrcalc), intent(inout) :: self
        integer  :: iref, k
        real(dp) :: prob_cc_odd(self%nrefs), prob_cc_even(self%nrefs)
        !$omp parallel do default(shared) private(k) proc_bind(close) schedule(static)
        do k = self%kfromto(1),self%kfromto(2)
            where( abs(self%regs_denom(:,k,:)) < TINY )
                self%refs_reg(:,k,:) = real(k, dp) * self%refs_reg(:,k,:)
            elsewhere
                self%refs_reg(:,k,:) = real(k, dp) * self%refs_reg(:,k,:) / self%regs_denom(:,k,:)
            endwhere
        enddo
        !$omp end parallel do
        if( trim(params_glob%reg_mode) == 'globdev' .or.  trim(params_glob%reg_mode) == 'neigh' )then
            prob_cc_odd  = 0._dp
            prob_cc_even = 0._dp
            !$omp parallel do default(shared) private(iref, k) proc_bind(close) schedule(static)
            do iref = 1, self%nrefs
                do k = self%kfromto(1),self%kfromto(2)
                    prob_cc_odd( iref) = prob_cc_odd( iref) + real(k, dp) * sum( self%refs_reg(:,k,iref) * real(self%pfts_refs_odd( :,k,iref), dp) )
                    prob_cc_even(iref) = prob_cc_even(iref) + real(k, dp) * sum( self%refs_reg(:,k,iref) * real(self%pfts_refs_even(:,k,iref), dp) )
                enddo
                prob_cc_odd( iref) = 1._dp + max(0._dp, prob_cc_odd( iref))
                prob_cc_even(iref) = 1._dp + max(0._dp, prob_cc_even(iref))
            enddo
            !$omp end parallel do
            !$omp parallel do default(shared) private(iref) proc_bind(close) schedule(static)
            do iref = 1, self%nrefs
                self%pfts_refs_even(:,:,iref) = self%pfts_refs_even(:,:,iref) + params_glob%eps * real(self%refs_reg(:,:,iref) / prob_cc_even(iref))
                self%pfts_refs_odd( :,:,iref) = self%pfts_refs_odd( :,:,iref) + params_glob%eps * real(self%refs_reg(:,:,iref) / prob_cc_odd( iref))
            enddo
            !$omp end parallel do
        else
            !$omp parallel do default(shared) private(iref) proc_bind(close) schedule(static)
            do iref = 1, self%nrefs
                self%pfts_refs_even(:,:,iref) = (1. - params_glob%eps) * self%pfts_refs_even(:,:,iref) + params_glob%eps * real(self%refs_reg(:,:,iref))
                self%pfts_refs_odd( :,:,iref) = (1. - params_glob%eps) * self%pfts_refs_odd( :,:,iref) + params_glob%eps * real(self%refs_reg(:,:,iref))
            enddo
            !$omp end parallel do
        endif
        
    end subroutine regularize_refs

    subroutine reset_regs( self )
        class(polarft_corrcalc), intent(inout) :: self
        self%refs_reg   = 0._dp
        self%regs_denom = 0._dp
    end subroutine reset_regs

    subroutine rotate_polar_real( self, ptcl_ctf, ptcl_ctf_rot, irot )
        class(polarft_corrcalc), intent(inout) :: self
        real(sp),                intent(in)    :: ptcl_ctf(    self%pftsz,self%kfromto(1):self%kfromto(2))
        real(dp),                intent(inout) :: ptcl_ctf_rot(self%pftsz,self%kfromto(1):self%kfromto(2))
        integer,                 intent(in)    :: irot
        integer :: rot
        if( irot >= self%pftsz + 1 )then
            rot = irot - self%pftsz
        else
            rot = irot
        end if
        ! just need the realpart
        if( irot == 1 .or. irot == self%pftsz + 1 )then
            ptcl_ctf_rot = ptcl_ctf
        else
            ptcl_ctf_rot(  1:rot-1    , :) = ptcl_ctf(self%pftsz-rot+2:self%pftsz      ,:)
            ptcl_ctf_rot(rot:self%pftsz,:) = ptcl_ctf(               1:self%pftsz-rot+1,:)
        end if
    end subroutine rotate_polar_real

    subroutine rotate_polar_complex( self, ptcl_ctf, ptcl_ctf_rot, irot )
        class(polarft_corrcalc), intent(inout) :: self
        complex(dp),             intent(in)    :: ptcl_ctf(    self%pftsz,self%kfromto(1):self%kfromto(2))
        complex(dp),             intent(inout) :: ptcl_ctf_rot(self%pftsz,self%kfromto(1):self%kfromto(2))
        integer,                 intent(in)    :: irot
        integer :: rot
        if( irot >= self%pftsz + 1 )then
            rot = irot - self%pftsz
        else
            rot = irot
        end if
        if( irot == 1 )then
            ptcl_ctf_rot = ptcl_ctf
        else if( irot <= self%pftsz )then
            ptcl_ctf_rot(rot:self%pftsz,:) =       ptcl_ctf(               1:self%pftsz-rot+1,:)
            ptcl_ctf_rot(  1:rot-1     ,:) = conjg(ptcl_ctf(self%pftsz-rot+2:self%pftsz      ,:))
        else if( irot == self%pftsz + 1 )then
            ptcl_ctf_rot = conjg(ptcl_ctf)
        else
            ptcl_ctf_rot(rot:self%pftsz,:) = conjg(ptcl_ctf(               1:self%pftsz-rot+1,:))
            ptcl_ctf_rot(  1:rot-1     ,:) =       ptcl_ctf(self%pftsz-rot+2:self%pftsz      ,:)
        end if
    end subroutine rotate_polar_complex

    ! Brings the reference in register with particle
    subroutine rotate_ref( self, ref, irot, ref_rot)
        class(polarft_corrcalc), intent(inout) :: self
        complex(dp),             intent(in)    :: ref(1:self%pftsz,self%kfromto(1):self%kfromto(2))
        integer,                 intent(in)    :: irot
        complex(dp),             intent(out)   :: ref_rot(1:self%pftsz,self%kfromto(1):self%kfromto(2))
        integer :: rot, jrot
        do jrot = 1,self%pftsz
            rot = jrot - (irot - 1) ! reverse rotation
            if( rot < 1 ) rot = rot + self%nrots
            if( rot > self%pftsz )then
                ref_rot(jrot,:) = conjg(ref(rot-self%pftsz,:))
            else
                ref_rot(jrot,:) = ref(rot,:)
            endif
        enddo
    end subroutine rotate_ref

    subroutine calc_polar_ctf( self, iptcl, smpd, kv, cs, fraca, dfx, dfy, angast )
        use simple_ctf,        only: ctf
        class(polarft_corrcalc),   intent(inout) :: self
        integer,                   intent(in)    :: iptcl
        real,                      intent(in)    :: smpd, kv, cs, fraca, dfx, dfy, angast
        type(ctf)       :: tfun
        real(sp)        :: spaFreqSq_mat(self%pftsz,self%kfromto(1):self%kfromto(2))
        real(sp)        :: ang_mat(self%pftsz,self%kfromto(1):self%kfromto(2))
        real(sp)        :: inv_ldim(3),hinv,kinv
        integer         :: i,irot,k
        if( .not.allocated(self%ctfmats) )then
            allocate(self%ctfmats(self%pftsz,self%kfromto(1):self%kfromto(2),1:self%nptcls), source=1.)
        endif
        ! if(.not. self%with_ctf ) return
        inv_ldim = 1./real(self%ldim)
        !$omp parallel do default(shared) private(irot,k,hinv,kinv) schedule(static) proc_bind(close)
        do irot=1,self%pftsz
            do k=self%kfromto(1),self%kfromto(2)
                hinv = self%polar(irot,k) * inv_ldim(1)
                kinv = self%polar(irot+self%nrots,k) * inv_ldim(2)
                spaFreqSq_mat(irot,k) = hinv*hinv+kinv*kinv
                ang_mat(irot,k)       = atan2(self%polar(irot+self%nrots,k),self%polar(irot,k))
            end do
        end do
        !$omp end parallel do
        i = self%pinds(iptcl)
        if( i > 0 )then
            tfun   = ctf(smpd, kv, cs, fraca)
            call tfun%init(dfx, dfy, angast)
            self%ctfmats(:,:,i) = tfun%eval(spaFreqSq_mat(:,:), ang_mat(:,:), 0.0, .not.params_glob%l_wiener_part)
        endif
    end subroutine calc_polar_ctf

    subroutine setup_npix_per_shell( self )
        class(polarft_corrcalc), intent(inout) :: self
        integer :: h,k,sh
        if( allocated(self%npix_per_shell) ) deallocate(self%npix_per_shell)
        allocate(self%npix_per_shell(self%kfromto(1):self%kfromto(2)),source=0.0)
        do h = 0,self%kfromto(2)
            do k = -self%kfromto(2),self%kfromto(2)
                if( (h==0) .and. (k>0) ) cycle
                sh = nint(sqrt(real(h**2+k**2)))
                if( sh < self%kfromto(1) ) cycle
                if( sh > self%kfromto(2) ) cycle
                self%npix_per_shell(sh) = self%npix_per_shell(sh) + 1.0
            end do
        end do
    end subroutine setup_npix_per_shell

    ! CALCULATORS

    subroutine create_polar_absctfmats( self, spproj, oritype, pfromto )
        use simple_ctf,        only: ctf
        use simple_sp_project, only: sp_project
        class(polarft_corrcalc),   intent(inout) :: self
        class(sp_project), target, intent(inout) :: spproj
        character(len=*),          intent(in)    :: oritype
        integer, optional,         intent(in)    :: pfromto(2)
        type(ctfparams) :: ctfparms(nthr_glob)
        type(ctf)       :: tfuns(nthr_glob)
        real(sp)        :: spaFreqSq_mat(self%pftsz,self%kfromto(1):self%kfromto(2))
        real(sp)        :: ang_mat(self%pftsz,self%kfromto(1):self%kfromto(2)), hinv,kinv
        integer         :: i,irot,k,iptcl,ithr,ppfromto(2),ctfmatind
        logical         :: present_pfromto
        present_pfromto = present(pfromto)
        ppfromto = self%pfromto
        if( present_pfromto ) ppfromto = pfromto
        if( allocated(self%ctfmats) ) deallocate(self%ctfmats)
        allocate(self%ctfmats(self%pftsz,self%kfromto(1):self%kfromto(2),1:self%nptcls), source=1.)
        if(.not. self%with_ctf ) return
        !$omp parallel do default(shared) private(irot,k,hinv,kinv) schedule(static) proc_bind(close)
        do irot=1,self%pftsz
            do k=self%kfromto(1),self%kfromto(2)
                hinv = self%polar(irot,k) / self%ldim(1)
                kinv = self%polar(irot+self%nrots,k) / self%ldim(2)
                spaFreqSq_mat(irot,k) = hinv*hinv+kinv*kinv
                ang_mat(irot,k)       = atan2(self%polar(irot+self%nrots,k),self%polar(irot,k))
            end do
        end do
        !$omp end parallel do
        if( params_glob%l_wiener_part )then
            ! taking into account CTF is intact before limit
            !$omp parallel do default(shared) private(i,iptcl,ctfmatind,ithr) schedule(static) proc_bind(close)
            do i=ppfromto(1),ppfromto(2)
                if( .not. present_pfromto )then
                    iptcl     = i
                    ctfmatind = i
                else
                    iptcl     = i
                    ctfmatind = i - ppfromto(1) + 1
                endif
                if( self%pinds(iptcl) > 0 )then
                    ithr           = omp_get_thread_num() + 1
                    ctfparms(ithr) = spproj%get_ctfparams(trim(oritype), iptcl)
                    tfuns(ithr)    = ctf(ctfparms(ithr)%smpd, ctfparms(ithr)%kv, ctfparms(ithr)%cs, ctfparms(ithr)%fraca)
                    call tfuns(ithr)%init(ctfparms(ithr)%dfx, ctfparms(ithr)%dfy, ctfparms(ithr)%angast)
                    if( ctfparms(ithr)%l_phaseplate )then
                        self%ctfmats(:,:,self%pinds(ctfmatind)) = abs(tfuns(ithr)%eval(spaFreqSq_mat(:,:), ang_mat(:,:), ctfparms(ithr)%phshift, .false. ))
                    else
                        self%ctfmats(:,:,self%pinds(ctfmatind)) = abs(tfuns(ithr)%eval(spaFreqSq_mat(:,:), ang_mat(:,:), 0.0,                    .false.))
                    endif
                endif
            end do
            !$omp end parallel do
        else
            !$omp parallel do default(shared) private(i,iptcl,ctfmatind,ithr) schedule(static) proc_bind(close)
            do i=ppfromto(1),ppfromto(2)
                if( .not. present_pfromto )then
                    iptcl     = i
                    ctfmatind = i
                else
                    iptcl     = i
                    ctfmatind = i - ppfromto(1) + 1
                endif
                if( self%pinds(iptcl) > 0 )then
                    ithr           = omp_get_thread_num() + 1
                    ctfparms(ithr) = spproj%get_ctfparams(trim(oritype), iptcl)
                    tfuns(ithr)    = ctf(ctfparms(ithr)%smpd, ctfparms(ithr)%kv, ctfparms(ithr)%cs, ctfparms(ithr)%fraca)
                    call tfuns(ithr)%init(ctfparms(ithr)%dfx, ctfparms(ithr)%dfy, ctfparms(ithr)%angast)
                    if( ctfparms(ithr)%l_phaseplate )then
                        self%ctfmats(:,:,self%pinds(ctfmatind)) = abs(tfuns(ithr)%eval(spaFreqSq_mat(:,:), ang_mat(:,:), ctfparms(ithr)%phshift) )
                    else
                        self%ctfmats(:,:,self%pinds(ctfmatind)) = abs(tfuns(ithr)%eval(spaFreqSq_mat(:,:), ang_mat(:,:)))
                    endif
                endif
            end do
            !$omp end parallel do
        endif
    end subroutine create_polar_absctfmats

    subroutine prep_ref4corr_sp( self, iref, iptcl, pft_ref, ithr)
        class(polarft_corrcalc), intent(inout) :: self
        integer,                 intent(in)    :: iref, iptcl
        complex(sp), pointer,    intent(out)   :: pft_ref(:, :)
        integer,                 intent(out)   :: ithr
        integer :: i
        i       =  self%pinds(iptcl)
        ithr    =  omp_get_thread_num() + 1
        pft_ref => self%heap_vars(ithr)%pft_ref
        ! copy
        if( self%iseven(i) )then
            pft_ref = self%pfts_refs_even(:,:,iref)
        else
            pft_ref = self%pfts_refs_odd(:,:,iref)
        endif
        ! multiply with CTF
        if( self%with_ctf ) pft_ref = pft_ref * self%ctfmats(:,:,i)
    end subroutine prep_ref4corr_sp

    subroutine prep_ref4corr_dp( self, iref, iptcl, pft_ref_8, ithr)
        class(polarft_corrcalc), intent(inout) :: self
        integer,                 intent(in)    :: iref, iptcl
        complex(dp), pointer,    intent(out)   :: pft_ref_8(:, :)
        integer,                 intent(out)   :: ithr
        integer :: i
        i         =  self%pinds(iptcl)
        ithr      =  omp_get_thread_num() + 1
        pft_ref_8 => self%heap_vars(ithr)%pft_ref_8
        ! copy
        if( self%iseven(i) )then
            pft_ref_8 = self%pfts_refs_even(:,:,iref)
        else
            pft_ref_8 = self%pfts_refs_odd(:,:,iref)
        endif
        ! multiply with CTF
        if( self%with_ctf ) pft_ref_8 = pft_ref_8 * self%ctfmats(:,:,i)
    end subroutine prep_ref4corr_dp

    !>  Generate polar shift matrix by means of de Moivre's formula, double precision
    subroutine gen_shmat_8( self, ithr, shift_8 , shmat_8 )
        class(polarft_corrcalc), intent(inout) :: self
        integer,                 intent(in)    :: ithr
        real(dp),                intent(in)    :: shift_8(2)
        complex(dp),    pointer, intent(inout) :: shmat_8(:,:)
        integer     :: k
        ! first shell, analytic
        self%heap_vars(ithr)%argvec = self%argtransf(:self%pftsz,  self%kfromto(1)) * shift_8(1) +&
                                    & self%argtransf(self%pftsz+1:,self%kfromto(1)) * shift_8(2)
        shmat_8(:,self%kfromto(1))  = dcmplx(dcos(self%heap_vars(ithr)%argvec), dsin(self%heap_vars(ithr)%argvec))
        ! one shell to the next
        self%heap_vars(ithr)%argvec = self%argtransf_shellone(:self%pftsz)   * shift_8(1) +&
                                    & self%argtransf_shellone(self%pftsz+1:) * shift_8(2)
        self%heap_vars(ithr)%shvec  = dcmplx(dcos(self%heap_vars(ithr)%argvec), dsin(self%heap_vars(ithr)%argvec))
        ! remaining shells, cos(kx)+isin(kx) = (cos(x)+isin(x))**k-1 * (cos(x)+isin(x))
        do k = self%kfromto(1)+1,self%kfromto(2)
            shmat_8(:,k) = shmat_8(:,k-1) * self%heap_vars(ithr)%shvec
        enddo
        ! alternative to:
        ! argmat  => self%heap_vars(ithr)%argmat_8
        ! argmat  =  self%argtransf(:self%pftsz,:)*shvec(1) + self%argtransf(self%pftsz + 1:,:)*shvec(2)
        ! shmat   =  cmplx(cos(argmat),sin(argmat),dp)
    end subroutine gen_shmat_8

    !>  Generate shift matrix following de Moivre's formula, single precision
    subroutine gen_shmat( self, ithr, shift, shmat )
        class(polarft_corrcalc), intent(inout) :: self
        integer,                 intent(in)    :: ithr
        real(sp),                intent(in)    :: shift(2)
        complex(sp),    pointer, intent(inout) :: shmat(:,:)
        call self%gen_shmat_8(ithr, real(shift,dp), self%heap_vars(ithr)%shmat_8)
        shmat = cmplx(self%heap_vars(ithr)%shmat_8)
    end subroutine gen_shmat

    ! Benchmarck for correlation calculation
    ! Is not FFT-accelerated, does not rely on memoization, for reference only
    real function calc_corr_rot_shift( self, iref, iptcl, shvec, irot)
        class(polarft_corrcalc), intent(inout) :: self
        integer,                 intent(in)    :: iref, iptcl
        real(sp),                intent(in)    :: shvec(2)
        integer,                 intent(in)    :: irot
        complex(dp), pointer :: pft_ref(:,:), shmat(:,:), pft_rot_ref(:,:)
        real(dp)    :: sqsumref, sqsumptcl, num
        integer     :: i, k, ithr!, rot
        i    = self%pinds(iptcl)
        ithr = omp_get_thread_num() + 1
        pft_ref     => self%heap_vars(ithr)%pft_ref_8
        pft_rot_ref => self%heap_vars(ithr)%pft_ref_tmp_8
        shmat       => self%heap_vars(ithr)%shmat_8
        if( self%iseven(i) )then
            pft_ref = self%pfts_refs_even(:,:,iref)
        else
            pft_ref = self%pfts_refs_odd(:,:,iref)
        endif
        call self%gen_shmat_8(ithr, real(shvec,dp),shmat)
        pft_ref = pft_ref * shmat
        call self%rotate_ref(pft_ref, irot, pft_rot_ref)
        if( self%with_ctf ) pft_rot_ref = pft_rot_ref * self%ctfmats(:,:,i)
        sqsumref  = 0.d0
        sqsumptcl = 0.d0
        num       = 0.d0
        do k = self%kfromto(1),self%kfromto(2)
            sqsumptcl = sqsumptcl + real(k,dp) * real(sum(self%pfts_ptcls(:,k,i) * conjg(self%pfts_ptcls(:,k,i))),dp)
            sqsumref  = sqsumref  + real(k,dp) * real(sum(pft_rot_ref(:,k)       * conjg(pft_rot_ref(:,k))) ,dp)
            num       = num       + real(k,dp) * real(sum(pft_rot_ref(:,k)       * conjg(self%pfts_ptcls(:,k,i))),dp)
        enddo
        calc_corr_rot_shift = real(num/sqrt(sqsumref*sqsumptcl))
    end function calc_corr_rot_shift

    subroutine calc_k_corrs( self, pft_ref, iptcl, k, kcorrs )
        class(polarft_corrcalc), intent(inout) :: self
        complex(sp),             intent(in)    :: pft_ref(1:self%pftsz,self%kfromto(1):self%kfromto(2))
        integer,                 intent(in)    :: iptcl, k
        real(dp),                intent(out)   :: kcorrs(self%nrots)
        integer :: ithr, i
        i = self%pinds(iptcl)
        ! get thread index
        ithr = omp_get_thread_num() + 1
        ! move reference into Fourier Fourier space (particles are memoized)
        self%fftdat(ithr)%ref_re(:) =  real(pft_ref(:,k))
        self%fftdat(ithr)%ref_im(:) = aimag(pft_ref(:,k)) * self%fft_factors
        call fftwf_execute_dft_r2c(self%plan_fwd_1, self%fftdat(ithr)%ref_re, self%fftdat(ithr)%ref_fft_re)
        call fftwf_execute_dft    (self%plan_fwd_2, self%fftdat(ithr)%ref_im, self%fftdat(ithr)%ref_fft_im)
        ! correlate FFTs
        self%fftdat(ithr)%ref_fft_re = conjg(self%fftdat(ithr)%ref_fft_re) * self%fftdat_ptcls(k,i)%re
        self%fftdat(ithr)%ref_fft_im = conjg(self%fftdat(ithr)%ref_fft_im) * self%fftdat_ptcls(k,i)%im
        self%fftdat(ithr)%product_fft(1:1 + 2 * int(self%pftsz / 2):2) = &
            4. * self%fftdat(ithr)%ref_fft_re(1:1 + int(self%pftsz / 2))
        self%fftdat(ithr)%product_fft(2:2 + 2 * int(self%pftsz / 2):2) = &
            4. * self%fftdat(ithr)%ref_fft_im(1:int(self%pftsz / 2) + 1)
        ! back transform
        call fftwf_execute_dft_c2r(self%plan_bwd, self%fftdat(ithr)%product_fft, self%fftdat(ithr)%backtransf)
        ! fftw3 routines are not properly normalized, hence division by self%nrots * 2
        kcorrs = real(self%fftdat(ithr)%backtransf / real(self%nrots * 2), dp)
    end subroutine calc_k_corrs

    function calc_corr_for_rot( self, pft_ref, iptcl, irot )result( corr )
        class(polarft_corrcalc), intent(inout) :: self
        integer,                 intent(in)    :: iptcl, irot
        complex(sp),             intent(in)    :: pft_ref(1:self%pftsz,self%kfromto(1):self%kfromto(2))
        integer :: rot, i
        real    :: corr
        complex :: tmp
        i    = self%pinds(iptcl)
        corr = 0.
        tmp  = 0.
        if( irot >= self%pftsz + 1 )then
            rot = irot - self%pftsz
        else
            rot = irot
        end if
        if( irot == 1 )then
            tmp = sum( pft_ref(:,self%kfromto(1):self%kfromto(2)) * conjg(self%pfts_ptcls(:,self%kfromto(1):self%kfromto(2),i)))
        else if( irot <= self%pftsz )then
            tmp =       sum( pft_ref(               1:self%pftsz-rot+1,self%kfromto(1):self%kfromto(2)) * conjg(self%pfts_ptcls(rot:self%pftsz,self%kfromto(1):self%kfromto(2),i)))
            tmp = tmp + sum( pft_ref(self%pftsz-rot+2:self%pftsz,      self%kfromto(1):self%kfromto(2)) *       self%pfts_ptcls(  1:rot-1,     self%kfromto(1):self%kfromto(2),i))
        else if( irot == self%pftsz + 1 )then
            tmp = sum( pft_ref(:,self%kfromto(1):self%kfromto(2)) * self%pfts_ptcls(:,self%kfromto(1):self%kfromto(2),i) )
        else
            tmp =       sum( pft_ref(               1:self%pftsz-rot+1,self%kfromto(1):self%kfromto(2)) *        self%pfts_ptcls(rot:self%pftsz,self%kfromto(1):self%kfromto(2),i))
            tmp = tmp + sum( pft_ref(self%pftsz-rot+2:self%pftsz,      self%kfromto(1):self%kfromto(2)) * conjg( self%pfts_ptcls(  1:rot-1,     self%kfromto(1):self%kfromto(2),i)))
        end if
        corr = real(tmp)
    end function calc_corr_for_rot

    function calc_corr_for_rot_8( self, pft_ref, iptcl, irot )result( corr )
        class(polarft_corrcalc), intent(inout) :: self
        integer,                 intent(in)    :: iptcl, irot
        complex(dp),             intent(in)    :: pft_ref(1:self%pftsz,self%kfromto(1):self%kfromto(2))
        integer     :: rot, i
        real(dp)    :: corr
        complex(dp) :: tmp
        i   = self%pinds(iptcl)
        tmp = 0.
        if( irot >= self%pftsz + 1 )then
            rot = irot - self%pftsz
        else
            rot = irot
        end if
        if( irot == 1 )then
            tmp = sum( pft_ref(:,self%kfromto(1):self%kfromto(2)) * conjg(self%pfts_ptcls(:,self%kfromto(1):self%kfromto(2),i)))
        else if( irot <= self%pftsz )then
            tmp =       sum( pft_ref(               1:self%pftsz-rot+1,self%kfromto(1):self%kfromto(2)) * conjg(self%pfts_ptcls(rot:self%pftsz,self%kfromto(1):self%kfromto(2),i)))
            tmp = tmp + sum( pft_ref(self%pftsz-rot+2:self%pftsz,      self%kfromto(1):self%kfromto(2)) *       self%pfts_ptcls(  1:rot-1,     self%kfromto(1):self%kfromto(2),i))
        else if( irot == self%pftsz + 1 ) then
            tmp = sum( pft_ref(:,self%kfromto(1):self%kfromto(2)) * self%pfts_ptcls(:,self%kfromto(1):self%kfromto(2),i) )
        else
            tmp =       sum( pft_ref(               1:self%pftsz-rot+1,self%kfromto(1):self%kfromto(2)) *        self%pfts_ptcls(rot:self%pftsz,self%kfromto(1):self%kfromto(2),i))
            tmp = tmp + sum( pft_ref(self%pftsz-rot+2:self%pftsz,      self%kfromto(1):self%kfromto(2)) * conjg( self%pfts_ptcls(  1:rot-1,     self%kfromto(1):self%kfromto(2),i)))
        end if
        corr = real(tmp, kind=dp)
    end function calc_corr_for_rot_8

    !<  \brief  compute the terms T1, T2 necessary for finding the derivative of the correlations, double precision
    subroutine calc_T1_T2_for_rot_8( self, pft_ref, pft_dref, iptcl, irot, nderivs, T1, T2)
        class(polarft_corrcalc), intent(inout) :: self
        integer,                 intent(in)    :: iptcl, irot, nderivs
        complex(dp),             intent(in)    :: pft_ref( 1:self%pftsz,self%kfromto(1):self%kfromto(2)), &
                                                  pft_dref(1:self%pftsz,self%kfromto(1):self%kfromto(2),nderivs)
        real(dp),                intent(out)   :: T1(self%kfromto(1):self%kfromto(2),nderivs), T2(self%kfromto(1):self%kfromto(2),nderivs)
        integer :: k, rot, j, i
        i = self%pinds(iptcl)
        if( irot >= self%pftsz + 1 )then
            rot = irot - self%pftsz
        else
            rot = irot
        end if
        do j = 1, nderivs
            do k = self%kfromto(1),self%kfromto(2)
                if( irot == 1 ) then
                    T1(k,j) = real(sum( pft_dref(:,k,j) * conjg(self%pfts_ptcls(:,k,i))), kind=dp)
                else if( irot <= self%pftsz ) then
                    T1(k,j) =           real(sum( pft_dref(               1:self%pftsz-rot+1,k,j) * conjg(self%pfts_ptcls(rot:self%pftsz,k,i))), kind=dp)
                    T1(k,j) = T1(k,j) + real(sum( pft_dref(self%pftsz-rot+2:self%pftsz,      k,j) *       self%pfts_ptcls(  1:rot-1,     k,i)),  kind=dp)
                else if( irot == self%pftsz + 1 ) then
                    T1(k,j) = real(sum( pft_dref(:,k,j) * self%pfts_ptcls(:,k,i) ), kind=dp)
                else
                    T1(k,j) =           real(sum( pft_dref(               1:self%pftsz-rot+1,k,j) *       self%pfts_ptcls(rot:self%pftsz,k,i)),   kind=dp)
                    T1(k,j) = T1(k,j) + real(sum( pft_dref(self%pftsz-rot+2:self%pftsz,      k,j) * conjg(self%pfts_ptcls(  1:rot-1,     k,i)) ), kind=dp)
                end if
                T2(k,j) = real(sum( pft_dref(:,k,j) * conjg(pft_ref(:,k))), kind=dp)
            end do
        end do
    end subroutine calc_T1_T2_for_rot_8

    function calc_euclid_for_rot( self, pft_ref, iptcl, irot ) result( euclid )
        class(polarft_corrcalc), intent(inout) :: self
        integer,                 intent(in)    :: iptcl, irot
        complex(sp),             intent(in)    :: pft_ref(1:self%pftsz,self%kfromto(1):self%kfromto(2))
        integer  :: rot, k, i
        real(sp) :: euclid, tmp
        i = self%pinds(iptcl)
        if( irot >= self%pftsz + 1 )then
            rot = irot - self%pftsz
        else
            rot = irot
        end if
        euclid = 0.
        do k = self%kfromto(1), self%kfromto(2)
            if( irot == 1 )then
                tmp =       sum(csq_fast(pft_ref(:,k) - self%pfts_ptcls(:,k,i)))
            else if( irot <= self%pftsz )then
                tmp =       sum(csq_fast(pft_ref(               1:self%pftsz-rot+1, k) -       self%pfts_ptcls(rot:self%pftsz,k,i)))
                tmp = tmp + sum(csq_fast(pft_ref(self%pftsz-rot+2:self%pftsz,       k) - conjg(self%pfts_ptcls(  1:rot-1,     k,i))))
            else if( irot == self%pftsz + 1 )then
                tmp = sum(csq_fast(pft_ref(:,k) - conjg(self%pfts_ptcls(:,k,i))))
            else
                tmp =       sum(csq_fast(pft_ref(               1:self%pftsz-rot+1, k) - conjg(self%pfts_ptcls(rot:self%pftsz,k,i))))
                tmp = tmp + sum(csq_fast(pft_ref(self%pftsz-rot+2:self%pftsz,       k) -       self%pfts_ptcls(  1:rot-1,     k,i)))
            end if
            euclid = euclid + tmp
        end do
    end function calc_euclid_for_rot

    function calc_prob_for_rot( self, pft_ref, iptcl, irot ) result( euclid_prob )
        class(polarft_corrcalc), intent(inout) :: self
        integer,                 intent(in)    :: iptcl, irot
        complex(sp),             intent(in)    :: pft_ref(1:self%pftsz,self%kfromto(1):self%kfromto(2))
        integer  :: rot, k, i
        real(sp) :: euclid_prob, tmp, denom
        i = self%pinds(iptcl)
        if( irot >= self%pftsz + 1 )then
            rot = irot - self%pftsz
        else
            rot = irot
        end if
        euclid_prob = 0.
        denom       = sum(csq_fast(self%pfts_ptcls(:, self%kfromto(1):self%kfromto(2),i)))
        do k = self%kfromto(1), self%kfromto(2)
            if( irot == 1 )then
                tmp =       sum(csq_fast(pft_ref(:,k) - self%pfts_ptcls(:,k,i)))
            else if( irot <= self%pftsz )then
                tmp =       sum(csq_fast(pft_ref(               1:self%pftsz-rot+1, k) -       self%pfts_ptcls(rot:self%pftsz,k,i)))
                tmp = tmp + sum(csq_fast(pft_ref(self%pftsz-rot+2:self%pftsz,       k) - conjg(self%pfts_ptcls(  1:rot-1,     k,i))))
            else if( irot == self%pftsz + 1 )then
                tmp = sum(csq_fast(pft_ref(:,k) - conjg(self%pfts_ptcls(:,k,i))))
            else
                tmp =       sum(csq_fast(pft_ref(               1:self%pftsz-rot+1, k) - conjg(self%pfts_ptcls(rot:self%pftsz,k,i))))
                tmp = tmp + sum(csq_fast(pft_ref(self%pftsz-rot+2:self%pftsz,       k) -       self%pfts_ptcls(  1:rot-1,     k,i)))
            end if
            euclid_prob = euclid_prob + exp( -tmp/denom )
        end do
        euclid_prob = euclid_prob / real(self%nk)
    end function calc_prob_for_rot

    function calc_euclid_for_rot_8( self, pft_ref, iptcl, irot ) result( euclid )
        class(polarft_corrcalc), intent(inout) :: self
        integer,                 intent(in)    :: iptcl, irot
        complex(dp),             intent(in)    :: pft_ref(1:self%pftsz,self%kfromto(1):self%kfromto(2))
        integer  :: rot, k, i
        real(dp) :: euclid, tmp
        i = self%pinds(iptcl)
        if( irot >= self%pftsz + 1 )then
            rot = irot - self%pftsz
        else
            rot = irot
        end if
        euclid = 0.d0
        do k = self%kfromto(1), self%kfromto(2)
            if( irot == 1 )then
                tmp =       sum(csq_fast(pft_ref(:,k) - self%pfts_ptcls(:,k,i)))
            else if( irot <= self%pftsz )then
                tmp =       sum(csq_fast(pft_ref(               1:self%pftsz-rot+1, k) -       self%pfts_ptcls(rot:self%pftsz,k,i)))
                tmp = tmp + sum(csq_fast(pft_ref(self%pftsz-rot+2:self%pftsz,       k) - conjg(self%pfts_ptcls(  1:rot-1,     k,i))))
            else if( irot == self%pftsz + 1 )then
                tmp = sum(csq_fast(pft_ref(:,k) - conjg(self%pfts_ptcls(:,k,i))))
            else
                tmp =       sum(csq_fast(pft_ref(               1:self%pftsz-rot+1, k) - conjg(self%pfts_ptcls(rot:self%pftsz,k,i))))
                tmp = tmp + sum(csq_fast(pft_ref(self%pftsz-rot+2:self%pftsz,       k) -       self%pfts_ptcls(  1:rot-1,     k,i)))
            end if
            euclid = euclid + tmp
        end do
    end function calc_euclid_for_rot_8

    function calc_prob_for_rot_8( self, pft_ref, iptcl, irot ) result( euclid_prob )
        class(polarft_corrcalc), intent(inout) :: self
        integer,                 intent(in)    :: iptcl, irot
        complex(dp),             intent(in)    :: pft_ref(1:self%pftsz,self%kfromto(1):self%kfromto(2))
        integer  :: rot, k, i
        real(dp) :: euclid_prob, tmp, denom
        i = self%pinds(iptcl)
        if( irot >= self%pftsz + 1 )then
            rot = irot - self%pftsz
        else
            rot = irot
        end if
        euclid_prob = 0.d0
        denom       = sum(real(csq_fast(self%pfts_ptcls(:, self%kfromto(1):self%kfromto(2),i)), dp))
        do k = self%kfromto(1), self%kfromto(2)
            if( irot == 1 )then
                tmp =       sum(csq_fast(pft_ref(:,k) - self%pfts_ptcls(:,k,i)))
            else if( irot <= self%pftsz )then
                tmp =       sum(csq_fast(pft_ref(               1:self%pftsz-rot+1, k) -       self%pfts_ptcls(rot:self%pftsz,k,i)))
                tmp = tmp + sum(csq_fast(pft_ref(self%pftsz-rot+2:self%pftsz,       k) - conjg(self%pfts_ptcls(  1:rot-1,     k,i))))
            else if( irot == self%pftsz + 1 )then
                tmp = sum(csq_fast(pft_ref(:,k) - conjg(self%pfts_ptcls(:,k,i))))
            else
                tmp =       sum(csq_fast(pft_ref(               1:self%pftsz-rot+1, k) - conjg(self%pfts_ptcls(rot:self%pftsz,k,i))))
                tmp = tmp + sum(csq_fast(pft_ref(self%pftsz-rot+2:self%pftsz,       k) -       self%pfts_ptcls(  1:rot-1,     k,i)))
            end if
            euclid_prob = euclid_prob + exp( -tmp/denom )
        end do
        euclid_prob = euclid_prob / real(self%nk)
    end function calc_prob_for_rot_8

    function calc_corrk_for_rot( self, pft_ref, iptcl, k, irot ) result( corr )
        class(polarft_corrcalc), intent(inout) :: self
        integer,                 intent(in)    :: iptcl, k, irot
        complex,                 intent(in)    :: pft_ref(1:self%pftsz,self%kfromto(1):self%kfromto(2))
        integer :: rot, i
        real    :: corr
        complex :: tmp
        i    = self%pinds(iptcl)
        corr = 0.
        tmp  = 0.
        if (irot >= self%pftsz + 1) then
            rot = irot - self%pftsz
        else
            rot = irot
        end if
        if (irot == 1) then
            tmp = sum( pft_ref(:,k) * conjg(self%pfts_ptcls(:,k,i)))
        else if (irot <= self%pftsz) then
            tmp =       sum( pft_ref(               1:self%pftsz-rot+1,k) * conjg(self%pfts_ptcls(rot:self%pftsz,k,i)))
            tmp = tmp + sum( pft_ref(self%pftsz-rot+2:self%pftsz,      k) *       self%pfts_ptcls(  1:rot-1,     k,i))
        else if (irot == self%pftsz + 1) then
            tmp = sum( pft_ref(:,k) * self%pfts_ptcls(:,k,i) )
        else
            tmp =       sum( pft_ref(               1:self%pftsz-rot+1,k) *        self%pfts_ptcls(rot:self%pftsz,k,i))
            tmp = tmp + sum( pft_ref(self%pftsz-rot+2:self%pftsz,      k) * conjg( self%pfts_ptcls(  1:rot-1,     k,i) ))
        end if
        corr = real(tmp)
    end function calc_corrk_for_rot

    function calc_corrk_for_rot_8( self, pft_ref, iptcl, k, irot ) result( corr )
        class(polarft_corrcalc), intent(inout) :: self
        integer,                 intent(in)    :: iptcl, k, irot
        complex(dp),             intent(in)    :: pft_ref(1:self%pftsz,self%kfromto(1):self%kfromto(2))
        integer     :: rot, i
        real(dp)    :: corr
        complex(dp) :: tmp
        i    = self%pinds(iptcl)
        corr = 0.d0
        tmp  = (0.d0, 0.d0)
        if (irot >= self%pftsz + 1) then
            rot = irot - self%pftsz
        else
            rot = irot
        end if
        if (irot == 1) then
            tmp = sum( pft_ref(:,k) * conjg(self%pfts_ptcls(:,k,i)))
        else if (irot <= self%pftsz) then
            tmp =       sum( pft_ref(               1:self%pftsz-rot+1,k) * conjg(self%pfts_ptcls(rot:self%pftsz,k,i)))
            tmp = tmp + sum( pft_ref(self%pftsz-rot+2:self%pftsz,      k) *       self%pfts_ptcls(  1:rot-1,     k,i))
        else if (irot == self%pftsz + 1) then
            tmp = sum( pft_ref(:,k) * self%pfts_ptcls(:,k,i) )
        else
            tmp =       sum( pft_ref(               1:self%pftsz-rot+1,k) *        self%pfts_ptcls(rot:self%pftsz,k,i))
            tmp = tmp + sum( pft_ref(self%pftsz-rot+2:self%pftsz,      k) * conjg( self%pfts_ptcls(  1:rot-1,     k,i) ))
        end if
        corr = real(tmp,dp)
    end function calc_corrk_for_rot_8

    real(sp) function calc_euclidk_for_rot( self, pft_ref, iptcl, k, irot )
        class(polarft_corrcalc), intent(inout) :: self
        integer,                 intent(in)    :: iptcl, irot, k
        complex(sp),             intent(in)    :: pft_ref(1:self%pftsz,self%kfromto(1):self%kfromto(2))
        real    :: euclid
        integer :: rot, i
        i = self%pinds(iptcl)
        if( irot >= self%pftsz + 1 )then
            rot = irot - self%pftsz
        else
            rot = irot
        end if
        if( irot == 1 )then
            euclid = sum(csq_fast(pft_ref(:,k) - self%pfts_ptcls(:,k,i)))
        else if( irot <= self%pftsz )then
            euclid =          sum(csq_fast(pft_ref(               1:self%pftsz-rot+1,k) -       self%pfts_ptcls(rot:self%pftsz,k,i)))
            euclid = euclid + sum(csq_fast(pft_ref(self%pftsz-rot+2:self%pftsz,      k) - conjg(self%pfts_ptcls(  1:rot-1,     k,i))))
        else if( irot == self%pftsz + 1 )then
            euclid = sum(csq_fast(pft_ref(:,k) - conjg(self%pfts_ptcls(:,k,i))))
        else
            euclid =          sum(csq_fast(pft_ref(               1:self%pftsz-rot+1,k) - conjg(self%pfts_ptcls(rot:self%pftsz,k,i))))
            euclid = euclid + sum(csq_fast(pft_ref(self%pftsz-rot+2:self%pftsz,      k) -       self%pfts_ptcls(  1:rot-1,     k,i)))
        end if
        calc_euclidk_for_rot = euclid
    end function calc_euclidk_for_rot

    real(dp) function calc_euclidk_for_rot_8( self, pft_ref, iptcl, k, irot )
        class(polarft_corrcalc), intent(inout) :: self
        integer,                 intent(in)    :: iptcl, irot, k
        complex(dp),             intent(in)    :: pft_ref(1:self%pftsz,self%kfromto(1):self%kfromto(2))
        real(dp) :: euclid
        integer  :: rot, i
        i = self%pinds(iptcl)
        if( irot >= self%pftsz + 1 )then
            rot = irot - self%pftsz
        else
            rot = irot
        end if
        if( irot == 1 )then
            euclid = sum(csq_fast(pft_ref(:,k) - self%pfts_ptcls(:,k,i)))
        else if( irot <= self%pftsz )then
            euclid =          sum(csq_fast(pft_ref(               1:self%pftsz-rot+1,k) -       self%pfts_ptcls(rot:self%pftsz,k,i)))
            euclid = euclid + sum(csq_fast(pft_ref(self%pftsz-rot+2:self%pftsz,      k) - conjg(self%pfts_ptcls(  1:rot-1,     k,i))))
        else if( irot == self%pftsz + 1 )then
            euclid = sum(csq_fast(pft_ref(:,k) - conjg(self%pfts_ptcls(:,k,i))))
        else
            euclid =          sum(csq_fast(pft_ref(               1:self%pftsz-rot+1,k) - conjg(self%pfts_ptcls(rot:self%pftsz,k,i))))
            euclid = euclid + sum(csq_fast(pft_ref(self%pftsz-rot+2:self%pftsz,      k) -       self%pfts_ptcls(  1:rot-1,     k,i)))
        end if
        calc_euclidk_for_rot_8 = euclid
    end function calc_euclidk_for_rot_8

    subroutine calc_frc( self, iref, iptcl, irot, shvec, frc )
        class(polarft_corrcalc),  intent(inout) :: self
        integer,                  intent(in)    :: iref, iptcl, irot
        real(sp),                 intent(in)    :: shvec(2)
        real(sp),                 intent(out)   :: frc(self%kfromto(1):self%kfromto(2))
        complex(sp), pointer :: pft_ref(:,:), shmat(:,:)
        real(dp),    pointer :: kcorrs(:)
        real(sp) :: sumsqref, sumsqptcl, denom
        integer  :: k, ithr, i
        if( L_CTFROTDEV )then
            call self%calc_frc_dev(iref, iptcl, irot, shvec, frc)
            return
        endif
        i = self%pinds(iptcl)
        call self%prep_ref4corr(iref, iptcl, pft_ref, ithr)
        kcorrs  => self%heap_vars(ithr)%kcorrs
        shmat   => self%heap_vars(ithr)%shmat
        call self%gen_shmat(ithr, shvec, shmat)
        pft_ref =  pft_ref * shmat
        do k=self%kfromto(1),self%kfromto(2)
            call self%calc_k_corrs(pft_ref, iptcl, k, kcorrs)
            sumsqptcl = sum(csq_fast(self%pfts_ptcls(:,k,i)))
            sumsqref  = sum(csq_fast(pft_ref(:,k)))
            denom     = sqrt(sumsqptcl * sumsqref)
            if( denom < 1.e-12 )then
                frc(k) = 0.
            else
                frc(k) = real(kcorrs(irot)) / denom
            endif
        end do
    end subroutine calc_frc

    ! Calculates frc between two PFTs, rotation, shift & ctf are not factored in
    subroutine calc_raw_frc( self, pft1, pft2, frc )
        class(polarft_corrcalc), intent(inout) :: self
        complex(dp),             intent(in)    :: pft1(self%pftsz,self%kfromto(1):self%kfromto(2))
        complex(dp),             intent(in)    :: pft2(self%pftsz,self%kfromto(1):self%kfromto(2))
        real,                    intent(out)   :: frc(self%kfromto(1):self%kfromto(2))
        real(dp) :: num, denom
        integer  :: k
        do k = self%kfromto(1),self%kfromto(2)
            num   = real(sum(pft1(:,k)*conjg(pft2(:,k))),dp)
            denom = real(sum(pft1(:,k)*conjg(pft1(:,k))),dp) * real(sum(pft2(:,k)*conjg(pft2(:,k))),dp)
            if( denom > DTINY )then
                frc(k) = real(num / dsqrt(denom))
            else
                frc(k) = 0.0
            endif
        end do
    end subroutine calc_raw_frc

    ! Calculates normalized PFT power spectrum
    subroutine calc_pspec( self, pft, pspec )
        class(polarft_corrcalc), intent(inout) :: self
        complex(dp),             intent(in)    :: pft(self%pftsz,self%kfromto(1):self%kfromto(2))
        real,                    intent(out)   :: pspec(self%kfromto(1):self%kfromto(2))
        integer :: k
        do k = self%kfromto(1),self%kfromto(2)
            pspec(k) = real( real(sum(pft(:,k)*conjg(pft(:,k))),dp) / real(self%pftsz,dp) )
        end do
    end subroutine calc_pspec

    function genmaxcorr_comlin( self, ieven, jeven ) result( cc_max )
        class(polarft_corrcalc), intent(inout) :: self
        integer,                 intent(in)    :: ieven, jeven
        complex(sp), pointer :: pft_ref_i(:,:), pft_ref_j(:,:)
        real     :: sqsum_i, sqsum_j, cc, cc_max
        integer  :: ithr, i, j
        ithr      =       omp_get_thread_num() + 1
        pft_ref_i =>      self%heap_vars(ithr)%pft_ref
        pft_ref_j =>      self%heap_vars(ithr)%pft_ref_tmp
        pft_ref_i =       self%pfts_refs_even(:,:,ieven)
        pft_ref_j = conjg(self%pfts_refs_even(:,:,jeven))
        ! no CTF to worry about since this is intended for class avgs
        cc_max = -1.
        do i = 1, self%pftsz
            sqsum_i = sum(csq_fast(pft_ref_i(i,self%kfromto(1):self%kfromto(2))))
            do j = 1, self%pftsz
                sqsum_j = sum(csq_fast(pft_ref_j(j,self%kfromto(1):self%kfromto(2))))
                cc      = sum(    real(pft_ref_i(i,self%kfromto(1):self%kfromto(2)) *&
                                       pft_ref_j(j,self%kfromto(1):self%kfromto(2))) ) / sqrt(sqsum_i * sqsum_j)
                if( cc > cc_max ) cc_max = cc
            end do
        end do
    end function genmaxcorr_comlin

    subroutine gencorrs_1( self, iref, iptcl, cc )
        class(polarft_corrcalc), intent(inout) :: self
        integer,                 intent(in)    :: iref, iptcl
        real(sp),                intent(out)   :: cc(self%nrots)
        complex(sp), pointer :: pft_ref(:,:)
        integer :: ithr
        if( L_CTFROTDEV )then
            call self%gencorrs_dev_1(iref, iptcl, cc)
            return
        endif
        call self%prep_ref4corr(iref, iptcl, pft_ref, ithr)
        select case(params_glob%cc_objfun)
            case(OBJFUN_CC)
                call self%gencorrs_cc(    pft_ref, iptcl, ithr, iref, self%heap_vars(ithr)%kcorrs)
                cc = real(self%heap_vars(ithr)%kcorrs)
            case(OBJFUN_EUCLID)
                call self%gencorrs_euclid(pft_ref, self%heap_vars(ithr)%kcorrs, iptcl, cc)
            case(OBJFUN_PROB)
                call self%gencorrs_prob(  pft_ref, self%heap_vars(ithr)%kcorrs, iptcl, cc)
        end select
    end subroutine gencorrs_1

    subroutine gencorrs_2( self, iref, iptcl, shvec, cc )
        class(polarft_corrcalc), intent(inout) :: self
        integer,                 intent(in)    :: iref, iptcl
        real(sp),                intent(in)    :: shvec(2)
        real(sp),                intent(out)   :: cc(self%nrots)
        complex(sp), pointer :: pft_ref(:,:), shmat(:,:)
        integer              :: ithr
        if( L_CTFROTDEV )then
            call self%gencorrs_dev_2(iref, iptcl, shvec, cc)
            return
        endif

        call self%prep_ref4corr(iref, iptcl, pft_ref, ithr)
        shmat => self%heap_vars(ithr)%shmat
        call self%gen_shmat(ithr, shvec, shmat)
        pft_ref = pft_ref * shmat
        select case(params_glob%cc_objfun)
            case(OBJFUN_CC)
                call self%gencorrs_cc(    pft_ref, iptcl, ithr, iref, self%heap_vars(ithr)%kcorrs)
                cc = real(self%heap_vars(ithr)%kcorrs)
            case(OBJFUN_EUCLID)
                call self%gencorrs_euclid(pft_ref, self%heap_vars(ithr)%kcorrs, iptcl, cc)
            case(OBJFUN_PROB)
                call self%gencorrs_prob(  pft_ref, self%heap_vars(ithr)%kcorrs, iptcl, cc)
        end select
    end subroutine gencorrs_2

    subroutine gencorrs_cc( self, pft_ref, iptcl, ithr, iref, cc )
        class(polarft_corrcalc), intent(inout) :: self
        complex(sp), pointer,    intent(in)    :: pft_ref(:,:)
        integer,                 intent(in)    :: iptcl, ithr, iref
        real(dp),                intent(out)   :: cc(self%nrots)
        real(dp) :: sqsum_ref, sqsum_ptcl
        integer  :: ik, i
        i          = self%pinds(iptcl)
        sqsum_ref  = 0._dp
        sqsum_ptcl = 0._dp
        cc         = 0._dp
        ! sum up correlations over k-rings
        do ik = self%kfromto(1),self%kfromto(2)
            sqsum_ref  = sqsum_ref  + real(ik, dp) * sum(real(csq_fast(pft_ref(:,ik)), dp))
            sqsum_ptcl = sqsum_ptcl + real(ik, dp) * sum(real(csq_fast(self%pfts_ptcls(:,ik,i)), dp))
            ! move reference into Fourier Fourier space (particles are memoized)
            self%fftdat(ithr)%ref_re(:) =  real(pft_ref(:,ik))
            self%fftdat(ithr)%ref_im(:) = aimag(pft_ref(:,ik)) * self%fft_factors
            call fftwf_execute_dft_r2c(self%plan_fwd_1, self%fftdat(ithr)%ref_re, self%fftdat(ithr)%ref_fft_re)
            call fftwf_execute_dft    (self%plan_fwd_2, self%fftdat(ithr)%ref_im, self%fftdat(ithr)%ref_fft_im)
            ! correlate
            self%fftdat(ithr)%ref_fft_re = conjg(self%fftdat(ithr)%ref_fft_re) * self%fftdat_ptcls(ik,i)%re
            self%fftdat(ithr)%ref_fft_im = conjg(self%fftdat(ithr)%ref_fft_im) * self%fftdat_ptcls(ik,i)%im
            self%fftdat(ithr)%product_fft(1:1 + 2 * int(self%pftsz / 2):2) = &
                4. * self%fftdat(ithr)%ref_fft_re(1:1+int(self%pftsz/2))
            self%fftdat(ithr)%product_fft(2:2 + 2 * int(self%pftsz / 2):2) = &
                4. * self%fftdat(ithr)%ref_fft_im(1:int(self%pftsz / 2) + 1)
            ! back transform
            call fftwf_execute_dft_c2r(self%plan_bwd, self%fftdat(ithr)%product_fft, self%fftdat(ithr)%backtransf)
            ! accumulate corrs
            cc = cc + real(ik, dp) * real(self%fftdat(ithr)%backtransf, dp)
        end do
        ! fftw3 routines are not properly normalized, hence division by self%nrots * 2
        cc = cc / real(self%nrots * 2, dp) / dsqrt(sqsum_ref * sqsum_ptcl)
    end subroutine gencorrs_cc

    subroutine gencorrs_euclid( self, pft_ref, keuclids, iptcl, euclids )
        class(polarft_corrcalc), intent(inout) :: self
        complex(sp), pointer,    intent(inout) :: pft_ref(:,:)
        real(dp),    pointer,    intent(inout) :: keuclids(:)
        integer,                 intent(in)    :: iptcl
        real(sp),                intent(out)   :: euclids(self%nrots)
        real(dp) :: denom, sumsqref, sumsqptcl
        integer  :: k, i
        i = self%pinds(iptcl)
        call self%weight_ref_ptcl(pft_ref, iptcl)
        euclids(:) = 0.
        denom      = sum(real(csq_fast(self%pfts_ptcls(:, self%kfromto(1):self%kfromto(2),i)), dp))
        do k=self%kfromto(1),self%kfromto(2)
            call self%calc_k_corrs(pft_ref, iptcl, k, keuclids)
            sumsqptcl = sum(real(csq_fast(self%pfts_ptcls(:,k,i)), dp))
            sumsqref  = sum(real(csq_fast(pft_ref(:,k)), dp))
            euclids   = euclids + real(sumsqptcl + sumsqref - 2. * keuclids(:))
        end do
        euclids = real(dexp( - euclids/denom ))
        call self%deweight_ref_ptcl(pft_ref, iptcl)
    end subroutine gencorrs_euclid

    subroutine gencorrs_prob( self, pft_ref, keuclids, iptcl, euclids )
        class(polarft_corrcalc), intent(inout) :: self
        complex(sp), pointer,    intent(inout) :: pft_ref(:,:)
        real(dp),    pointer,    intent(inout) :: keuclids(:)
        integer,                 intent(in)    :: iptcl
        real(sp),                intent(out)   :: euclids(self%nrots)
        real(dp) :: denom, sumsqref, sumsqptcl
        integer  :: k, i
        i = self%pinds(iptcl)
        call self%weight_ref_ptcl(pft_ref, iptcl)
        euclids(:) = 0.
        denom      = sum(real(csq_fast(self%pfts_ptcls(:, self%kfromto(1):self%kfromto(2),i)), dp))
        do k=self%kfromto(1),self%kfromto(2)
            call self%calc_k_corrs(pft_ref, iptcl, k, keuclids)
            sumsqptcl = sum(real(csq_fast(self%pfts_ptcls(:,k,i)), dp))
            sumsqref  = sum(real(csq_fast(pft_ref(:,k)), dp))
            euclids   = euclids + real(dexp( -(sumsqptcl + sumsqref - 2. * keuclids(:))/denom ))
        end do
        euclids = euclids/real(self%kfromto(2) - self%kfromto(1) + 1)
        call self%deweight_ref_ptcl(pft_ref, iptcl)
    end subroutine gencorrs_prob

    real(dp) function gencorr_for_rot_8( self, iref, iptcl, shvec, irot )
        class(polarft_corrcalc), intent(inout) :: self
        integer,                 intent(in)    :: iref, iptcl
        real(dp),                intent(in)    :: shvec(2)
        integer,                 intent(in)    :: irot
        complex(dp), pointer :: pft_ref_8(:,:), shmat_8(:,:)
        integer              :: ithr
        if( L_CTFROTDEV )then
            gencorr_for_rot_8 = self%gencorr_for_rot_8_dev(iref, iptcl, shvec, irot)
            return
        endif
        call self%prep_ref4corr(iref, iptcl, pft_ref_8, ithr)
        shmat_8 => self%heap_vars(ithr)%shmat_8
        call self%gen_shmat_8(ithr, shvec, shmat_8)
        pft_ref_8 = pft_ref_8 * shmat_8
        gencorr_for_rot_8 = 0._dp
        select case(params_glob%cc_objfun)
            case(OBJFUN_CC)
                gencorr_for_rot_8 = self%gencorr_cc_for_rot_8(     pft_ref_8, iptcl, irot, iref )
            case(OBJFUN_EUCLID)
                gencorr_for_rot_8 = self%gencorr_euclid_for_rot_8( pft_ref_8, iptcl, irot )
            case(OBJFUN_PROB)
                gencorr_for_rot_8 = self%gencorr_prob_for_rot_8(   pft_ref_8, iptcl, irot )
        end select
    end function gencorr_for_rot_8

    function gencorr_cc_for_rot_8( self, pft_ref_8, iptcl, irot, iref ) result( cc )
        class(polarft_corrcalc), intent(inout) :: self
        complex(dp), pointer,    intent(inout) :: pft_ref_8(:,:)
        integer,                 intent(in)    :: iptcl, irot, iref
        real(dp) :: cc, sqsum_ref, sqsum_ptcl
        integer  :: ik, i
        i          = self%pinds(iptcl)
        sqsum_ref  = 0._dp
        sqsum_ptcl = 0._dp
        cc         = 0._dp
        do ik = self%kfromto(1),self%kfromto(2)
            sqsum_ref  = sqsum_ref  + real(ik,kind=dp) * sum(real(csq_fast(pft_ref_8(:,ik)), dp))
            sqsum_ptcl = sqsum_ptcl + real(ik,kind=dp) * sum(real(csq_fast(self%pfts_ptcls(:,ik,i)), dp))
            cc         = cc         + real(ik,kind=dp) * self%calc_corrk_for_rot_8(pft_ref_8, iptcl, ik, irot)
        end do
        cc = cc / dsqrt(sqsum_ref * sqsum_ptcl)
    end function gencorr_cc_for_rot_8

    real(dp) function gencorr_euclid_for_rot_8( self, pft_ref_8, iptcl, irot )
        class(polarft_corrcalc), intent(inout) :: self
        complex(dp), pointer,    intent(inout) :: pft_ref_8(:,:)
        integer,                 intent(in)    :: iptcl, irot
        real(dp) :: denom
        integer  :: i
        i = self%pinds(iptcl)
        call self%weight_ref_ptcl(pft_ref_8, iptcl)
        denom                    = sum(real(csq_fast(self%pfts_ptcls(:, self%kfromto(1):self%kfromto(2),i)), dp))
        gencorr_euclid_for_rot_8 = dexp( - self%calc_euclid_for_rot_8(pft_ref_8, iptcl, irot)/denom )
        call self%deweight_ref_ptcl(pft_ref_8, iptcl)
    end function gencorr_euclid_for_rot_8

    real(dp) function gencorr_prob_for_rot_8( self, pft_ref_8, iptcl, irot )
        class(polarft_corrcalc), intent(inout) :: self
        complex(dp), pointer,    intent(inout) :: pft_ref_8(:,:)
        integer,                 intent(in)    :: iptcl, irot
        call self%weight_ref_ptcl(pft_ref_8, iptcl)
        gencorr_prob_for_rot_8 = self%calc_prob_for_rot_8(pft_ref_8, iptcl, irot)
        call self%deweight_ref_ptcl(pft_ref_8, iptcl)
    end function gencorr_prob_for_rot_8

    subroutine gencorr_grad_for_rot_8( self, iref, iptcl, shvec, irot, f, grad )
        class(polarft_corrcalc), intent(inout) :: self
        integer,                 intent(in)    :: iref, iptcl
        real(dp),                intent(in)    :: shvec(2)
        integer,                 intent(in)    :: irot
        real(dp),                intent(out)   :: f, grad(2)
        complex(dp), pointer :: pft_ref_8(:,:), shmat_8(:,:), pft_ref_tmp(:,:)
        integer              :: ithr
        if( L_CTFROTDEV )then
            call self%gencorr_grad_for_rot_8_dev(iref, iptcl, shvec, irot, f, grad)
            return
        endif
        call self%prep_ref4corr(iref, iptcl, pft_ref_8, ithr)
        pft_ref_tmp => self%heap_vars(ithr)%pft_ref_tmp_8
        shmat_8     => self%heap_vars(ithr)%shmat_8
        call self%gen_shmat_8(ithr, shvec, shmat_8)
        pft_ref_8 = pft_ref_8 * shmat_8
        select case(params_glob%cc_objfun)
            case(OBJFUN_CC)
                call self%gencorr_cc_grad_for_rot_8(     pft_ref_8, pft_ref_tmp, iptcl, irot, iref, f, grad )
            case(OBJFUN_EUCLID)
                call self%gencorr_euclid_grad_for_rot_8( pft_ref_8, pft_ref_tmp, iptcl, irot, f, grad )
            case(OBJFUN_PROB)
                call self%gencorr_prob_grad_for_rot_8(   pft_ref_8, pft_ref_tmp, iptcl, irot, f, grad )
        end select
    end subroutine gencorr_grad_for_rot_8

    subroutine gencorr_cc_grad_for_rot_8( self, pft_ref, pft_ref_tmp, iptcl, irot, iref, f, grad )
        class(polarft_corrcalc), intent(inout) :: self
        complex(dp), pointer,    intent(inout) :: pft_ref(:,:), pft_ref_tmp(:,:)
        integer,                 intent(in)    :: iptcl, irot, iref
        real(dp),                intent(out)   :: f, grad(2)
        real(dp) :: sqsum_ref, sqsum_ptcl
        integer  :: ik, i
        ! use jacobian resolution weights
        i           = self%pinds(iptcl)
        sqsum_ref   = 0._dp
        sqsum_ptcl  = 0._dp
        f           = 0._dp
        grad        = 0._dp
        pft_ref_tmp = pft_ref * (0.d0, 1.d0) * self%argtransf(:self%pftsz,:)
        do ik = self%kfromto(1),self%kfromto(2)
            sqsum_ref  = sqsum_ref  + real(ik,kind=dp) * sum(real(csq_fast(pft_ref(:,ik)), dp))
            sqsum_ptcl = sqsum_ptcl + real(ik,kind=dp) * sum(real(csq_fast(self%pfts_ptcls(:,ik,i)), dp))
            f          = f          + real(ik,kind=dp) * self%calc_corrk_for_rot_8(pft_ref,     iptcl, ik, irot)
            grad(1)    = grad(1)    + real(ik,kind=dp) * self%calc_corrk_for_rot_8(pft_ref_tmp, iptcl, ik, irot)
        end do
        pft_ref_tmp = pft_ref * (0.d0, 1.d0) * self%argtransf(self%pftsz + 1:,:)
        do ik = self%kfromto(1),self%kfromto(2)
            grad(2) = grad(2) + real(ik,kind=dp) * self%calc_corrk_for_rot_8(pft_ref_tmp, iptcl, ik, irot)
        end do
        f    = f    / dsqrt(sqsum_ref*sqsum_ptcl)
        grad = grad / dsqrt(sqsum_ref*sqsum_ptcl)
    end subroutine gencorr_cc_grad_for_rot_8

    subroutine gencorr_euclid_grad_for_rot_8( self, pft_ref, pft_ref_tmp, iptcl, irot, f, grad )
        class(polarft_corrcalc), intent(inout) :: self
        complex(dp), pointer,    intent(inout) :: pft_ref(:,:), pft_ref_tmp(:,:)
        integer,                 intent(in)    :: iptcl, irot
        real(dp),                intent(out)   :: f, grad(2)
        real(dp) :: diffsq, denom
        integer  :: k, i
        i = self%pinds(iptcl)
        call self%weight_ref_ptcl(pft_ref, iptcl)
        denom       = sum(real(csq_fast(self%pfts_ptcls(:, self%kfromto(1):self%kfromto(2),i)), dp))
        f           = self%calc_euclid_for_rot_8(pft_ref, iptcl, irot)
        grad        = 0._dp
        pft_ref_tmp = pft_ref * (0.d0, 1.d0) * self%argtransf(:self%pftsz,:)
        do k = self%kfromto(1), self%kfromto(2)
            diffsq  = real(sum(pft_ref_tmp(:,k)*conjg(pft_ref(:,k)))) - self%calc_corrk_for_rot_8(pft_ref_tmp, iptcl, k, irot)
            grad(1) = grad(1) + diffsq
        end do
        pft_ref_tmp = pft_ref * (0.d0, 1.d0) * self%argtransf(self%pftsz + 1:,:)
        do k = self%kfromto(1), self%kfromto(2)
            diffsq  = real(sum(pft_ref_tmp(:,k)*conjg(pft_ref(:,k)))) - self%calc_corrk_for_rot_8(pft_ref_tmp, iptcl, k, irot)
            grad(2) = grad(2) + diffsq
        end do
        f    = dexp( -f/denom )
        grad = -f * 2._dp * grad/denom
        call self%deweight_ref_ptcl(pft_ref, iptcl)
    end subroutine gencorr_euclid_grad_for_rot_8

    subroutine gencorr_prob_grad_for_rot_8( self, pft_ref, pft_ref_tmp, iptcl, irot, f, grad )
        class(polarft_corrcalc), intent(inout) :: self
        complex(dp), pointer,    intent(inout) :: pft_ref(:,:), pft_ref_tmp(:,:)
        integer,                 intent(in)    :: iptcl, irot
        real(dp),                intent(out)   :: f, grad(2)
        real(dp) :: diffsq, denom, gradsq
        integer  :: k, i
        i = self%pinds(iptcl)
        call self%weight_ref_ptcl(pft_ref, iptcl)
        denom       = sum(real(csq_fast(self%pfts_ptcls(:, self%kfromto(1):self%kfromto(2),i)), dp))
        f           = self%calc_prob_for_rot_8(pft_ref, iptcl, irot)
        grad        = 0._dp
        pft_ref_tmp = pft_ref * (0.d0, 1.d0) * self%argtransf(:self%pftsz,:)
        do k = self%kfromto(1), self%kfromto(2)
            diffsq  = self%calc_euclidk_for_rot_8(pft_ref, iptcl, k, irot)
            gradsq  = real(sum(pft_ref_tmp(:,k)*conjg(pft_ref(:,k)))) - self%calc_corrk_for_rot_8(pft_ref_tmp, iptcl, k, irot)
            grad(1) = grad(1) - dexp(-diffsq/denom) * 2._dp * gradsq/denom
        end do
        pft_ref_tmp = pft_ref * (0.d0, 1.d0) * self%argtransf(self%pftsz + 1:,:)
        do k = self%kfromto(1), self%kfromto(2)
            diffsq  = self%calc_euclidk_for_rot_8(pft_ref, iptcl, k, irot)
            gradsq  = real(sum(pft_ref_tmp(:,k)*conjg(pft_ref(:,k)))) - self%calc_corrk_for_rot_8(pft_ref_tmp, iptcl, k, irot)
            grad(2) = grad(2) - dexp(-diffsq/denom) * 2._dp * gradsq/denom
        end do
        grad = grad / real(self%nk)
        call self%deweight_ref_ptcl(pft_ref, iptcl)
    end subroutine gencorr_prob_grad_for_rot_8

    subroutine gencorr_grad_only_for_rot_8( self, iref, iptcl, shvec, irot, grad )
        class(polarft_corrcalc), intent(inout) :: self
        integer,                 intent(in)    :: iref, iptcl
        real(dp),                intent(in)    :: shvec(2)
        integer,                 intent(in)    :: irot
        real(dp),                intent(out)   :: grad(2)
        complex(dp), pointer :: pft_ref_8(:,:), shmat_8(:,:), pft_ref_tmp(:,:)
        integer              :: ithr, i
        if( L_CTFROTDEV )then
            call self%gencorr_grad_only_for_rot_8_dev(iref, iptcl, shvec, irot, grad)
            return
        endif
        call self%prep_ref4corr(iref, iptcl, pft_ref_8, ithr)
        pft_ref_tmp => self%heap_vars(ithr)%pft_ref_tmp_8
        shmat_8     => self%heap_vars(ithr)%shmat_8
        call self%gen_shmat_8(ithr, shvec, shmat_8)
        pft_ref_8 = pft_ref_8 * shmat_8
        select case(params_glob%cc_objfun)
            case(OBJFUN_CC)
                call self%gencorr_cc_grad_only_for_rot_8(     pft_ref_8, pft_ref_tmp, iptcl, i, irot, iref, grad )
            case(OBJFUN_EUCLID)
                call self%gencorr_euclid_grad_only_for_rot_8( pft_ref_8, pft_ref_tmp, iptcl, i, irot, grad )
            case(OBJFUN_PROB)
                call self%gencorr_prob_grad_only_for_rot_8(   pft_ref_8, pft_ref_tmp, iptcl, i, irot, grad )
        end select
    end subroutine gencorr_grad_only_for_rot_8

    subroutine gencorr_cc_grad_only_for_rot_8( self, pft_ref, pft_ref_tmp, iptcl, i, irot, iref, grad )
        class(polarft_corrcalc), intent(inout) :: self
        complex(dp), pointer,    intent(inout) :: pft_ref(:,:), pft_ref_tmp(:,:)
        integer,                 intent(in)    :: iptcl, i, irot, iref
        real(dp),                intent(out)   :: grad(2)
        integer  :: ik
        real(dp) :: sqsum_ref, sqsum_ptcl
        sqsum_ref  = 0._dp
        sqsum_ptcl = 0._dp
        grad       = 0._dp
        pft_ref_tmp = pft_ref * (0.d0, 1.d0) * self%argtransf(:self%pftsz,:)
        do ik = self%kfromto(1),self%kfromto(2)
            sqsum_ref  = sqsum_ref  + real(ik,kind=dp) * sum(real(csq_fast(pft_ref(:,ik)), dp))
            sqsum_ptcl = sqsum_ptcl + real(ik,kind=dp) * sum(real(csq_fast(self%pfts_ptcls(:,ik,i)), dp))
            grad(1)    = grad(1)    + real(ik,kind=dp) * self%calc_corrk_for_rot_8(pft_ref_tmp, iptcl, ik, irot)
        end do
        pft_ref_tmp = pft_ref * (0.d0, 1.d0) * self%argtransf(self%pftsz + 1:,:)
        do ik = self%kfromto(1),self%kfromto(2)
            grad(2) = grad(2) + real(ik,kind=dp) * self%calc_corrk_for_rot_8(pft_ref_tmp, iptcl, ik, irot)
        end do
        grad = grad / dsqrt(sqsum_ref*sqsum_ptcl)
    end subroutine gencorr_cc_grad_only_for_rot_8

    subroutine gencorr_euclid_grad_only_for_rot_8( self, pft_ref, pft_ref_tmp, iptcl, i, irot, grad )
        class(polarft_corrcalc), intent(inout) :: self
        complex(dp), pointer,    intent(inout) :: pft_ref(:,:), pft_ref_tmp(:,:)
        integer,                 intent(in)    :: iptcl, i, irot
        real(dp),                intent(out)   :: grad(2)
        real(dp) :: f
        call self%gencorr_euclid_grad_for_rot_8(pft_ref, pft_ref_tmp, iptcl, irot, f, grad)
    end subroutine gencorr_euclid_grad_only_for_rot_8

    subroutine gencorr_prob_grad_only_for_rot_8( self, pft_ref, pft_ref_tmp, iptcl, i, irot, grad )
        class(polarft_corrcalc), intent(inout) :: self
        complex(dp), pointer,    intent(inout) :: pft_ref(:,:), pft_ref_tmp(:,:)
        integer,                 intent(in)    :: iptcl, i, irot
        real(dp),                intent(out)   :: grad(2)
        real(dp) :: f
        call self%gencorr_prob_grad_for_rot_8(pft_ref, pft_ref_tmp, iptcl, irot, f, grad)
    end subroutine gencorr_prob_grad_only_for_rot_8

    function gencorr_cont_grad_cc_for_rot_8( self, iref, iptcl, shvec, irot, dcc ) result( cc )
        class(polarft_corrcalc), intent(inout) :: self
        integer,                 intent(in)    :: iref, iptcl
        real(dp),                intent(in)    :: shvec(2)
        integer,                 intent(in)    :: irot
        real(dp),                intent(out)   :: dcc(3)
        complex(dp), pointer :: pft_ref(:,:), shmat(:,:), pft_dref(:,:,:)
        real(dp),    pointer :: fdf_y(:), fdf_T1(:,:), fdf_T2(:,:)
        real(dp) :: cc, sqsum_ref, denom
        integer  :: ithr, j, k
        ithr     =  omp_get_thread_num() + 1
        pft_ref  => self%heap_vars(ithr)%pft_ref_8
        pft_dref => self%heap_vars(ithr)%pft_dref_8
        shmat    => self%heap_vars(ithr)%shmat_8
        fdf_y    => self%heap_vars(ithr)%fdf_y_8
        fdf_T1   => self%heap_vars(ithr)%fdf_T1_8
        fdf_T2   => self%heap_vars(ithr)%fdf_T2_8
        call self%gen_shmat_8(ithr, shvec, shmat)
        if( self%iseven(self%pinds(iptcl)) )then
            pft_ref  = self%pfts_refs_even(:,:,iref)
            pft_dref = self%pfts_drefs_even(:,:,:,iref)
        else
            pft_ref  = self%pfts_refs_odd(:,:,iref)
            pft_dref = self%pfts_drefs_odd(:,:,:,iref)
        endif
        if( self%with_ctf )then
            pft_ref = (pft_ref * self%ctfmats(:,:,self%pinds(iptcl))) * shmat
            do j = 1,3
                pft_dref(:,:,j) = (pft_dref(:,:,j) * self%ctfmats(:,:,self%pinds(iptcl))) * shmat
            end do
        else
            pft_ref = pft_ref * shmat
            do j = 1,3
                pft_dref(:,:,j) = pft_dref(:,:,j) * shmat
            end do
        endif
        sqsum_ref = sum(csq_fast(pft_ref(:,self%kfromto(1):self%kfromto(2))))
        denom     = sqrt(sqsum_ref * self%sqsums_ptcls(self%pinds(iptcl)))
        do k = self%kfromto(1),self%kfromto(2)
            fdf_y(k) = self%calc_corrk_for_rot_8(pft_ref, iptcl, k, irot)
        end do
        call self%calc_T1_T2_for_rot_8( pft_ref, pft_dref, iref, irot, 3, fdf_T1, fdf_T2)
        cc = sum(fdf_y) / denom
        do j = 1,3
            dcc(j) = ( sum(fdf_T1(:,j)) - sum(fdf_y) * sum(fdf_T2(:,j)) / sqsum_ref ) / denom
        end do
    end function gencorr_cont_grad_cc_for_rot_8

    real(dp) function gencorr_cont_cc_for_rot_8( self, iref, iptcl, shvec, irot )
        class(polarft_corrcalc), intent(inout) :: self
        integer,                 intent(in)    :: iref, iptcl
        real(dp),                intent(in)    :: shvec(2)
        integer,                 intent(in)    :: irot
        complex(dp), pointer :: pft_ref(:,:), shmat(:,:)
        real(dp),    pointer :: fdf_y(:)
        real(dp) :: sqsum_ref, denom
        integer  :: ithr, k
        ithr        =  omp_get_thread_num() + 1
        pft_ref     => self%heap_vars(ithr)%pft_ref_8
        shmat       => self%heap_vars(ithr)%shmat_8
        fdf_y       => self%heap_vars(ithr)%fdf_y_8
        call self%gen_shmat_8(ithr, shvec, shmat)
        if( self%iseven(self%pinds(iptcl)) )then
            pft_ref  = self%pfts_refs_even(:,:,iref)
        else
            pft_ref  = self%pfts_refs_odd(:,:,iref)
        endif
        if( self%with_ctf )then
            pft_ref = (pft_ref * self%ctfmats(:,:,self%pinds(iptcl))) * shmat
        else
            pft_ref = pft_ref * shmat
        endif
        sqsum_ref = sum(csq_fast(pft_ref(:,self%kfromto(1):self%kfromto(2))))
        denom     = sqrt(sqsum_ref * self%sqsums_ptcls(self%pinds(iptcl)))
        do k = self%kfromto(1),self%kfromto(2)
            fdf_y(k) = self%calc_corrk_for_rot_8(pft_ref, iptcl, k, irot)
        end do
        gencorr_cont_cc_for_rot_8 = sum(fdf_y) / denom
    end function gencorr_cont_cc_for_rot_8

    subroutine gencorr_cont_shift_grad_cc_for_rot_8( self, iref, iptcl, shvec, irot, f, grad )
        class(polarft_corrcalc), intent(inout) :: self
        integer,                 intent(in)    :: iref, iptcl
        real(dp),                intent(in)    :: shvec(2)
        integer,                 intent(in)    :: irot
        real(dp),                intent(out)   :: f
        real(dp),                intent(out)   :: grad(5) ! 3 orientation angles, 2 shifts
        complex(dp), pointer :: pft_ref(:,:), pft_ref_tmp(:,:), shmat(:,:), pft_dref(:,:,:)
        real(dp),    pointer :: fdf_y(:), fdf_T1(:,:), fdf_T2(:,:)
        real(dp) :: sqsum_ref, denom, corr
        integer  :: ithr, j, k
        ithr        =  omp_get_thread_num() + 1
        pft_ref     => self%heap_vars(ithr)%pft_ref_8
        pft_ref_tmp => self%heap_vars(ithr)%pft_ref_tmp_8
        pft_dref    => self%heap_vars(ithr)%pft_dref_8
        shmat       => self%heap_vars(ithr)%shmat_8
        fdf_y       => self%heap_vars(ithr)%fdf_y_8
        fdf_T1      => self%heap_vars(ithr)%fdf_T1_8
        fdf_T2      => self%heap_vars(ithr)%fdf_T2_8
        call self%gen_shmat_8(ithr, shvec, shmat)
        if( self%iseven(self%pinds(iptcl)) )then
            pft_ref  = self%pfts_refs_even(:,:,iref)
            pft_dref = self%pfts_drefs_even(:,:,:,iref)
        else
            pft_ref  = self%pfts_refs_odd(:,:,iref)
            pft_dref = self%pfts_drefs_odd(:,:,:,iref)
        endif
        if( self%with_ctf )then
            pft_ref = (pft_ref * self%ctfmats(:,:,self%pinds(iptcl))) * shmat
            do j = 1,3
                pft_dref(:,:,j) = (pft_dref(:,:,j) * self%ctfmats(:,:,self%pinds(iptcl))) * shmat
            end do
        else
            pft_ref = pft_ref * shmat
            do j = 1,3
                pft_dref(:,:,j) = pft_dref(:,:,j) * shmat
            end do
        endif
        sqsum_ref = sum(csq_fast(pft_ref(:,self%kfromto(1):self%kfromto(2))))
        denom     = sqrt(sqsum_ref * self%sqsums_ptcls(self%pinds(iptcl)))
        do k = self%kfromto(1),self%kfromto(2)
            fdf_y(k) = self%calc_corrk_for_rot_8(pft_ref, iptcl, k, irot)
        end do
        call self%calc_T1_T2_for_rot_8( pft_ref, pft_dref, iref, irot, 3, fdf_T1, fdf_T2)
        f = sum(fdf_y) / denom
        do j = 1,3
            grad(j) = ( sum(fdf_T1(:,j)) - sum(fdf_y) * sum(fdf_T2(:,j)) / sqsum_ref ) / denom
        end do
        pft_ref_tmp = pft_ref * (0.d0, 1.d0) * self%argtransf(:self%pftsz,:)
        corr        = self%calc_corr_for_rot_8(pft_ref_tmp, iptcl, irot)
        grad(4)     = corr / denom
        pft_ref_tmp = pft_ref * (0.d0, 1.d0) * self%argtransf(self%pftsz + 1:,:)
        corr        = self%calc_corr_for_rot_8(pft_ref_tmp, iptcl, irot)
        grad(5)     = corr / denom
    end subroutine gencorr_cont_shift_grad_cc_for_rot_8

    function gencorr_cont_grad_euclid_for_rot_8( self, iref, iptcl, shvec, irot, dcc ) result( cc )
        class(polarft_corrcalc), intent(inout) :: self
        integer,                 intent(in)    :: iref, iptcl
        real(dp),                intent(in)    :: shvec(2)
        integer,                 intent(in)    :: irot
        real(dp),                intent(out)   :: dcc(3)
        real(dp) :: cc
        cc = 0._dp
        ! TODO: implement me
    end function gencorr_cont_grad_euclid_for_rot_8

    subroutine gencorr_cont_shift_grad_euclid_for_rot_8( self, iref, iptcl, shvec, irot, f, grad )
        class(polarft_corrcalc), intent(inout) :: self
        integer,                 intent(in)    :: iref, iptcl
        real(dp),                intent(in)    :: shvec(2)
        integer,                 intent(in)    :: irot
        real(dp),                intent(out)   :: f
        real(dp),                intent(out)   :: grad(5) ! 3 orientation angles, 2 shifts
        ! TODO: implement me
    end subroutine gencorr_cont_shift_grad_euclid_for_rot_8

    subroutine gencorr_sigma_contrib( self, iref, iptcl, shvec, irot, sigma_contrib)
        class(polarft_corrcalc), intent(inout) :: self
        integer,                 intent(in)    :: iref, iptcl
        real(sp),                intent(in)    :: shvec(2)
        integer,                 intent(in)    :: irot
        real(sp),                intent(out)   :: sigma_contrib(self%kfromto(1):self%kfromto(2))
        complex(sp), pointer :: pft_ref(:,:), shmat(:,:)
        integer  :: ithr, k
        if( L_CTFROTDEV )then
            call self%gencorr_sigma_contrib_dev(iref, iptcl, shvec, irot, sigma_contrib)
            return
        endif
        call self%prep_ref4corr(iref, iptcl, pft_ref, ithr)
        shmat => self%heap_vars(ithr)%shmat
        call self%gen_shmat(ithr, shvec, shmat)
        pft_ref = pft_ref * shmat
        do k = self%kfromto(1), self%kfromto(2)
            sigma_contrib(k) = 0.5 * self%calc_euclidk_for_rot(pft_ref, iptcl, k, irot) / real(self%pftsz)
        end do
    end subroutine gencorr_sigma_contrib

    !< updating sigma for this particle/reference pair
    subroutine update_sigma( self, iref, iptcl, shvec, irot)
        class(polarft_corrcalc), intent(inout) :: self
        integer,                 intent(in)    :: iref, iptcl
        real(sp),                intent(in)    :: shvec(2)
        integer,                 intent(in)    :: irot
        call self%gencorr_sigma_contrib( iref, iptcl, shvec, irot, self%sigma2_noise(self%kfromto(1):self%kfromto(2), iptcl))
    end subroutine update_sigma

    real function specscore_1( self, iref, iptcl, irot )
        class(polarft_corrcalc), intent(inout) :: self
        integer,                 intent(in)    :: iref, iptcl, irot
        real :: frc(self%kfromto(1):self%kfromto(2))
        call self%calc_frc(iref, iptcl, irot, [0.0,0.0], frc )
        specscore_1 = max(0.,median_nocopy(frc))
    end function specscore_1

    real function specscore_2( self, iref, iptcl, irot, shvec )
        class(polarft_corrcalc), intent(inout) :: self
        integer,                 intent(in)    :: iref, iptcl, irot
        real,                    intent(in)    :: shvec(2)
        real :: frc(self%kfromto(1):self%kfromto(2))
        call self%calc_frc(iref, iptcl, irot, shvec, frc )
        specscore_2 = max(0.,median_nocopy(frc))
    end function specscore_2

    subroutine weight_ref_ptcl_sp( self, pft_ref, iptcl )
        class(polarft_corrcalc), intent(inout) :: self
        complex(sp),    pointer, intent(inout) :: pft_ref(:,:)
        integer,                 intent(in)    :: iptcl
        integer  :: k, i
        real(sp) :: w
        i = self%pinds(iptcl)
        do k=self%kfromto(1),self%kfromto(2)
            w                         = sqrt( self%delta * k / self%sigma2_noise(k,iptcl) )
            pft_ref(:,k)              = w * pft_ref(:,k)
            self%pfts_ptcls(:,k,i)    = w * self%pfts_ptcls(:,k,i)
            self%fftdat_ptcls(k,i)%re = w * self%fftdat_ptcls(k,i)%re
            self%fftdat_ptcls(k,i)%im = w * self%fftdat_ptcls(k,i)%im
        end do
    end subroutine weight_ref_ptcl_sp

    subroutine weight_ref_ptcl_dp( self, pft_ref, iptcl )
        class(polarft_corrcalc), intent(inout) :: self
        complex(dp),    pointer, intent(inout) :: pft_ref(:,:)
        integer,                 intent(in)    :: iptcl
        integer  :: k, i
        real(dp) :: w
        i = self%pinds(iptcl)
        do k=self%kfromto(1),self%kfromto(2)
            w                         = dsqrt( real(self%delta, dp) * k / self%sigma2_noise(k,iptcl) )
            pft_ref(:,k)              =      w  * pft_ref(:,k)
            self%pfts_ptcls(:,k,i)    = real(w) * self%pfts_ptcls(:,k,i)
            self%fftdat_ptcls(k,i)%re = real(w) * self%fftdat_ptcls(k,i)%re
            self%fftdat_ptcls(k,i)%im = real(w) * self%fftdat_ptcls(k,i)%im
        end do
    end subroutine weight_ref_ptcl_dp

    subroutine deweight_ref_ptcl_sp( self, pft_ref, iptcl )
        class(polarft_corrcalc), intent(inout) :: self
        complex(sp),    pointer, intent(inout) :: pft_ref(:,:)
        integer,                 intent(in)    :: iptcl
        integer  :: k, i
        real(sp) :: w
        i = self%pinds(iptcl)
        do k=self%kfromto(1),self%kfromto(2)
            w                         = sqrt( self%delta * k / self%sigma2_noise(k,iptcl) )
            pft_ref(:,k)              = pft_ref(:,k)              / w
            self%pfts_ptcls(:,k,i)    = self%pfts_ptcls(:,k,i)    / w
            self%fftdat_ptcls(k,i)%re = self%fftdat_ptcls(k,i)%re / w
            self%fftdat_ptcls(k,i)%im = self%fftdat_ptcls(k,i)%im / w
        end do
    end subroutine deweight_ref_ptcl_sp

    subroutine deweight_ref_ptcl_dp( self, pft_ref, iptcl )
        class(polarft_corrcalc), intent(inout) :: self
        complex(dp),    pointer, intent(inout) :: pft_ref(:,:)
        integer,                 intent(in)    :: iptcl
        integer  :: k, i
        real(dp) :: w
        i = self%pinds(iptcl)
        do k=self%kfromto(1),self%kfromto(2)
            w                         = dsqrt( real(self%delta, dp) * k / self%sigma2_noise(k,iptcl) )
            pft_ref(:,k)              = pft_ref(:,k)              / w
            self%pfts_ptcls(:,k,i)    = self%pfts_ptcls(:,k,i)    / real(w)
            self%fftdat_ptcls(k,i)%re = self%fftdat_ptcls(k,i)%re / real(w)
            self%fftdat_ptcls(k,i)%im = self%fftdat_ptcls(k,i)%im / real(w)
        end do
    end subroutine deweight_ref_ptcl_dp

    ! DESTRUCTOR

    subroutine kill( self )
        class(polarft_corrcalc), intent(inout) :: self
        integer :: ithr, i, ik
        if( self%existence )then
            do ithr=1,params_glob%nthr
                call fftwf_free(self%fftdat(ithr)%p_ref_re)
                call fftwf_free(self%fftdat(ithr)%p_ref_im)
                call fftwf_free(self%fftdat(ithr)%p_ref_fft_re)
                call fftwf_free(self%fftdat(ithr)%p_ref_fft_im)
                call fftwf_free(self%fftdat(ithr)%p_product_fft)
                call fftwf_free(self%fftdat(ithr)%p_backtransf)
                call fftwf_free(self%fft_carray(ithr)%p_re)
                call fftwf_free(self%fft_carray(ithr)%p_im)
                deallocate(self%heap_vars(ithr)%pft_ref,self%heap_vars(ithr)%pft_ref_tmp,&
                    &self%heap_vars(ithr)%pft_dref,&
                    &self%heap_vars(ithr)%argvec, self%heap_vars(ithr)%shvec,&
                    &self%heap_vars(ithr)%corrs_over_k,&
                    &self%heap_vars(ithr)%shmat,self%heap_vars(ithr)%kcorrs,&
                    &self%heap_vars(ithr)%pft_ref_8,self%heap_vars(ithr)%pft_ref_tmp_8,&
                    &self%heap_vars(ithr)%pft_dref_8,&
                    &self%heap_vars(ithr)%shmat_8,self%heap_vars(ithr)%argmat_8,&
                    &self%heap_vars(ithr)%fdf_y_8,self%heap_vars(ithr)%fdf_T1_8,&
                    &self%heap_vars(ithr)%fdf_T2_8)
            end do
            do i = 1, self%nptcls
                do ik = self%kfromto(1),self%kfromto(2)
                    call fftwf_free(self%fftdat_ptcls(ik,i)%p_re)
                    call fftwf_free(self%fftdat_ptcls(ik,i)%p_im)
                end do
            end do
            if( allocated(self%ctfmats)        ) deallocate(self%ctfmats)
            if( allocated(self%npix_per_shell) ) deallocate(self%npix_per_shell)
            deallocate( self%sqsums_ptcls, self%wsqsums_ptcls, self%angtab, self%argtransf,&
                &self%polar, self%pfts_refs_even, self%pfts_refs_odd, self%pfts_drefs_even, self%pfts_drefs_odd,&
                self%pfts_ptcls, self%fft_factors, self%fftdat, self%fftdat_ptcls, self%fft_carray,&
                &self%iseven, self%pinds, self%heap_vars, self%argtransf_shellone,self%refs_reg,self%regs_denom)
            call fftwf_destroy_plan(self%plan_bwd)
            call fftwf_destroy_plan(self%plan_fwd_1)
            call fftwf_destroy_plan(self%plan_fwd_2)
            if( L_CTFROTDEV )then
                call self%kill_memoized_ptcls
                call self%kill_memoized_refs
            endif
            nullify(self%sigma2_noise, pftcc_glob)
            self%l_filt_set   = .false.
            self%existence    = .false.
        endif
    end subroutine kill

    ! For testing

    subroutine memoize_ptcls( self )
        class(polarft_corrcalc), intent(inout) :: self
        type(fftw_cvec), allocatable :: cvecs1(:), cvecs2(:)
        type(c_ptr)                  :: fwd_plan, r2c_plan
        integer :: ithr,i,k
        ! allocations
        allocate(cvecs1(nthr_glob),cvecs2(nthr_glob))
        do ithr = 1,nthr_glob
            cvecs1(ithr)%p = fftwf_alloc_complex(int(self%nrots, c_size_t))
            cvecs2(ithr)%p = fftwf_alloc_complex(int(self%nrots, c_size_t))
            call c_f_pointer(cvecs1(ithr)%p, cvecs1(ithr)%c, [self%nrots])
            call c_f_pointer(cvecs2(ithr)%p, cvecs2(ithr)%c, [self%nrots])
        enddo
        fwd_plan = fftwf_plan_dft_1d(    self%nrots, cvecs1(1)%c,     cvecs2(1)%c, FFTW_FORWARD, ior(FFTW_PATIENT, FFTW_USE_WISDOM))
        r2c_plan = fftwf_plan_dft_r2c_1d(self%nrots, self%rvec1(1)%r, self%cvec1(1)%c,           ior(FFTW_PATIENT, FFTW_USE_WISDOM))
        ! memoization
        !$omp parallel do collapse(2) private(i,k,ithr) default(shared) proc_bind(close) schedule(static)
        do i = 1,self%nptcls
            do k = self%kfromto(1),self%kfromto(2)
                ithr = omp_get_thread_num() + 1
                ! FT(X.CTF)
                if( self%with_ctf )then
                    cvecs1(ithr)%c(1:self%pftsz) = self%pfts_ptcls(:,k,i) * self%ctfmats(:,k,i)
                else
                    cvecs1(ithr)%c(1:self%pftsz) = self%pfts_ptcls(:,k,i)
                endif
                cvecs1(ithr)%c(self%pftsz+1:self%nrots) = conjg(cvecs1(ithr)%c(1:self%pftsz))
                call fftwf_execute_dft(fwd_plan, cvecs1(ithr)%c, cvecs2(ithr)%c)
                self%ft_ptcl_ctf(k,i)%c(1:self%pftsz+1) = cvecs2(ithr)%c(1:self%pftsz+1)
                ! FT(CTF2)
                if( self%with_ctf )then
                    self%rvec1(ithr)%r(1:self%pftsz)            = self%ctfmats(:,k,i)*self%ctfmats(:,k,i)
                    self%rvec1(ithr)%r(self%pftsz+1:self%nrots) = self%rvec1(ithr)%r(1:self%pftsz)
                else
                    self%rvec1(ithr)%r = 1.0
                endif
                call fftwf_execute_dft_r2c(r2c_plan, self%rvec1(ithr)%r, self%cvec1(ithr)%c)
                self%ft_ctf2(k,i)%c(1:self%pftsz+1) = self%cvec1(ithr)%c(1:self%pftsz+1)
            enddo
        enddo
        !$omp end parallel do
        ! clean-up
        do ithr = 1,nthr_glob
            call fftwf_free(cvecs1(ithr)%p)
            call fftwf_free(cvecs2(ithr)%p)
        enddo
        deallocate(cvecs1,cvecs2)
        call fftwf_destroy_plan(fwd_plan)
        call fftwf_destroy_plan(r2c_plan)
    end subroutine memoize_ptcls

    subroutine memoize_refs( self )
        class(polarft_corrcalc), intent(inout) :: self
        type(fftw_cvec), allocatable :: cvecs1(:)
        type(c_ptr)                  :: r2c_plan
        integer :: k, ithr, iref
        ! allocations
        call self%allocate_refs_memoization
        allocate(cvecs1(nthr_glob))
        r2c_plan = fftwf_plan_dft_r2c_1d(self%nrots, self%rvec1(1)%r, self%cvec1(1)%c, ior(FFTW_PATIENT, FFTW_USE_WISDOM))
        do ithr = 1,nthr_glob
            cvecs1(ithr)%p = fftwf_alloc_complex(int(self%nrots, c_size_t))
            call c_f_pointer(cvecs1(ithr)%p, cvecs1(ithr)%c, [self%nrots])
        enddo
        ! memoization
        !$omp parallel do collapse(2) private(iref,k,ithr) default(shared) proc_bind(close) schedule(static)
        do iref = 1,self%nrefs
            do k = self%kfromto(1),self%kfromto(2)
                ithr = omp_get_thread_num() + 1
                ! FT(REFeven)*
                self%cvec2(ithr)%c(           1:self%pftsz) = self%pfts_refs_even(:,k,iref)
                self%cvec2(ithr)%c(self%pftsz+1:self%nrots) = conjg(self%cvec2(ithr)%c(1:self%pftsz))
                call fftwf_execute_dft(self%plan_fwd1, self%cvec2(ithr)%c, self%cvec2(ithr)%c)
                self%ft_ref_even(k,iref)%c = conjg(self%cvec2(ithr)%c(1:self%pftsz+1))
                ! FT(REFodd)*
                self%cvec2(ithr)%c(           1:self%pftsz) = self%pfts_refs_odd(:,k,iref)
                self%cvec2(ithr)%c(self%pftsz+1:self%nrots) = conjg(self%cvec2(ithr)%c(1:self%pftsz))
                call fftwf_execute_dft(self%plan_fwd1, self%cvec2(ithr)%c, self%cvec2(ithr)%c)
                self%ft_ref_odd(k,iref)%c = conjg(self%cvec2(ithr)%c(1:self%pftsz+1))
                ! FT(REF2even)*
                self%rvec1(ithr)%r(           1:self%pftsz) = csq_fast(self%pfts_refs_even(:,k,iref))
                self%rvec1(ithr)%r(self%pftsz+1:self%nrots) = self%rvec1(ithr)%r(1:self%pftsz)
                call fftwf_execute_dft_r2c(r2c_plan, self%rvec1(ithr)%r, self%cvec1(ithr)%c)
                self%ft_ref2_even(k,iref)%c = conjg(self%cvec1(ithr)%c(1:self%pftsz+1))
                ! FT(REF2odd)*
                self%rvec1(ithr)%r(           1:self%pftsz) = csq_fast(self%pfts_refs_odd(:,k,iref))
                self%rvec1(ithr)%r(self%pftsz+1:self%nrots) = self%rvec1(ithr)%r(1:self%pftsz)
                call fftwf_execute_dft_r2c(r2c_plan, self%rvec1(ithr)%r, self%cvec1(ithr)%c)
                self%ft_ref2_odd(k,iref)%c = conjg(self%cvec1(ithr)%c(1:self%pftsz+1))
            enddo
        enddo
        !$omp end parallel do
        ! clean-up
        do ithr = 1,nthr_glob
            call fftwf_free(cvecs1(ithr)%p)
        enddo
        deallocate(cvecs1)
        call fftwf_destroy_plan(r2c_plan)
    end subroutine memoize_refs

    subroutine kill_memoized_ptcls( self )
        class(polarft_corrcalc), intent(inout) :: self
        integer :: i,j
        do i = 1,size(self%ft_ptcl_ctf,dim=1)
            do j = 1,size(self%ft_ptcl_ctf,dim=2)
                call fftwf_free(self%ft_ptcl_ctf(i,j)%p)
                call fftwf_free(self%ft_ctf2(i,j)%p)
            enddo
        enddo
        deallocate(self%ft_ptcl_ctf,self%ft_ctf2)
    end subroutine kill_memoized_ptcls

    subroutine kill_memoized_refs( self )
        class(polarft_corrcalc), intent(inout) :: self
        integer :: i,j
        do i = 1,size(self%ft_ref_even,dim=1)
            do j = 1,size(self%ft_ref_even,dim=2)
                call fftwf_free(self%ft_ref_even(i,j)%p)
                call fftwf_free(self%ft_ref_odd(i,j)%p)
                call fftwf_free(self%ft_ref2_even(i,j)%p)
                call fftwf_free(self%ft_ref2_odd(i,j)%p)
            enddo
        enddo
        do i = 1,size(self%cvec1,dim=1)
            call fftwf_free(self%cvec1(i)%p)
            call fftwf_free(self%cvec2(i)%p)
            call fftw_free(self%drvec(i)%p)
        enddo
        deallocate(self%ft_ref_even,self%ft_ref_odd,&
        &self%ft_ref2_even,self%ft_ref2_odd,self%cvec1,self%cvec2,self%drvec)
        call fftwf_destroy_plan(self%plan_fwd1)
        call fftwf_destroy_plan(self%plan_bwd1)
    end subroutine kill_memoized_refs

    subroutine allocate_ptcls_memoization( self )
        class(polarft_corrcalc), intent(inout) :: self
        integer :: i, k
        allocate(self%ft_ptcl_ctf(self%kfromto(1):self%kfromto(2),self%nptcls),&
                &self%ft_ctf2(self%kfromto(1):self%kfromto(2),self%nptcls))
        do i = 1,self%nptcls
            do k = self%kfromto(1),self%kfromto(2)
                self%ft_ptcl_ctf(k,i)%p = fftwf_alloc_complex(int(self%pftsz+1, c_size_t))
                self%ft_ctf2(k,i)%p     = fftwf_alloc_complex(int(self%pftsz+1, c_size_t))
                call c_f_pointer(self%ft_ptcl_ctf(k,i)%p, self%ft_ptcl_ctf(k,i)%c, [self%pftsz+1])
                call c_f_pointer(self%ft_ctf2(    k,i)%p, self%ft_ctf2(    k,i)%c, [self%pftsz+1])
            enddo
        enddo
    end subroutine allocate_ptcls_memoization

    subroutine allocate_refs_memoization( self )
        class(polarft_corrcalc), intent(inout) :: self
        integer :: ik, ithr, iref
        allocate(self%ft_ref_even( self%kfromto(1):self%kfromto(2),self%nrefs),&
        &self%ft_ref_odd(  self%kfromto(1):self%kfromto(2),self%nrefs),&
        &self%ft_ref2_even(self%kfromto(1):self%kfromto(2),self%nrefs),&
        &self%ft_ref2_odd( self%kfromto(1):self%kfromto(2),self%nrefs),&
        &self%rvec1(nthr_glob), self%cvec1(nthr_glob),self%cvec2(nthr_glob),&
        &self%drvec(nthr_glob))
        ! convenience objects
        do ithr = 1,nthr_glob
            self%cvec1(ithr)%p = fftwf_alloc_complex(int(self%pftsz+1, c_size_t))
            self%cvec2(ithr)%p = fftwf_alloc_complex(int(self%nrots, c_size_t))
            call c_f_pointer(self%cvec1(ithr)%p, self%cvec1(ithr)%c, [self%pftsz+1])
            call c_f_pointer(self%cvec1(ithr)%p, self%rvec1(ithr)%r, [self%nrots+2])
            call c_f_pointer(self%cvec2(ithr)%p, self%cvec2(ithr)%c, [self%nrots])
            self%drvec(ithr)%p = fftw_alloc_real(int(self%nrots, c_size_t))
            call c_f_pointer(self%drvec(ithr)%p, self%drvec(ithr)%r, [self%nrots])
        enddo
        ! references
        do iref = 1,self%nrefs
            do ik = self%kfromto(1),self%kfromto(2)
                self%ft_ref_even( ik,iref)%p = fftwf_alloc_complex(int(self%pftsz+1,c_size_t))
                self%ft_ref_odd(  ik,iref)%p = fftwf_alloc_complex(int(self%pftsz+1,c_size_t))
                self%ft_ref2_even(ik,iref)%p = fftwf_alloc_complex(int(self%pftsz+1,c_size_t))
                self%ft_ref2_odd( ik,iref)%p = fftwf_alloc_complex(int(self%pftsz+1,c_size_t))
                call c_f_pointer(self%ft_ref_even( ik,iref)%p, self%ft_ref_even( ik,iref)%c, [self%pftsz+1])
                call c_f_pointer(self%ft_ref_odd(  ik,iref)%p, self%ft_ref_odd(  ik,iref)%c, [self%pftsz+1])
                call c_f_pointer(self%ft_ref2_even(ik,iref)%p, self%ft_ref2_even(ik,iref)%c, [self%pftsz+1])
                call c_f_pointer(self%ft_ref2_odd( ik,iref)%p, self%ft_ref2_odd( ik,iref)%c, [self%pftsz+1])
            enddo
        enddo
        self%plan_fwd1 = fftwf_plan_dft_1d(    self%nrots, self%cvec2(1)%c, self%cvec2(1)%c, FFTW_FORWARD, ior(FFTW_PATIENT, FFTW_USE_WISDOM))
        self%plan_bwd1 = fftwf_plan_dft_c2r_1d(self%nrots, self%cvec1(1)%c, self%rvec1(1)%r,               ior(FFTW_PATIENT, FFTW_USE_WISDOM))
    end subroutine allocate_refs_memoization

    subroutine gencorrs_dev_1( self, iref, iptcl, cc )
        class(polarft_corrcalc), intent(inout) :: self
        integer,                 intent(in)    :: iref, iptcl
        real(sp),                intent(out)   :: cc(self%nrots)
        select case(params_glob%cc_objfun)
            case(OBJFUN_CC)
                call self%gencorrs_cc_dev(iptcl, iref, cc)
            case(OBJFUN_EUCLID)
                call self%gencorrs_euclid_dev(iptcl, iref, cc)
            case(OBJFUN_PROB)
                call self%gencorrs_prob_dev(iptcl, iref, cc)
        end select
    end subroutine gencorrs_dev_1

    subroutine gencorrs_dev_2( self, iref, iptcl, shift, cc )
        class(polarft_corrcalc), intent(inout) :: self
        integer,                 intent(in)    :: iref, iptcl
        real(sp),                intent(in)    :: shift(2)
        real(sp),                intent(out)   :: cc(self%nrots)
        complex(sp), pointer :: pft_ref(:,:), shmat(:,:)
        integer :: i, ithr
        ithr    = omp_get_thread_num() + 1
        i       = self%pinds(iptcl)
        shmat   => self%heap_vars(ithr)%shmat
        pft_ref => self%heap_vars(ithr)%pft_ref
        call self%gen_shmat(ithr, shift, shmat)
        if( self%iseven(i) )then
            pft_ref = shmat * self%pfts_refs_even(:,:,iref)
        else
            pft_ref = shmat * self%pfts_refs_odd(:,:,iref)
        endif
        select case(params_glob%cc_objfun)
            case(OBJFUN_CC)
                call self%gencorrs_shifted_cc_dev(pft_ref, iptcl, iref, shift, cc)
            case(OBJFUN_EUCLID)
                call self%gencorrs_shifted_euclid_dev(pft_ref, iptcl, iref, shift, cc)
            case(OBJFUN_PROB)
                call self%gencorrs_shifted_prob_dev(pft_ref, iptcl, iref, shift, cc)
        end select
    end subroutine gencorrs_dev_2

    subroutine gencorrs_cc_dev( self, iptcl, iref, corrs)
        class(polarft_corrcalc), intent(inout) :: self
        integer,                 intent(in)    :: iptcl, iref
        real(sp),                intent(out)   :: corrs(self%nrots)
        integer  :: k, i, ithr
        logical  :: even
        ithr = omp_get_thread_num() + 1
        i    = self%pinds(iptcl)
        even = self%iseven(i)
        self%heap_vars(ithr)%kcorrs = 0.d0
        self%drvec(ithr)%r          = 0.d0
        do k = self%kfromto(1),self%kfromto(2)
            ! FT(CTF2) x FT(REF2)*)
            if( even )then
                self%cvec1(ithr)%c(1:self%pftsz+1) = self%ft_ctf2(k,i)%c(1:self%pftsz+1) * self%ft_ref2_even(k,iref)%c(1:self%pftsz+1)
            else
                self%cvec1(ithr)%c(1:self%pftsz+1) = self%ft_ctf2(k,i)%c(1:self%pftsz+1) * self%ft_ref2_odd(k,iref)%c(1:self%pftsz+1)
            endif
            ! IFFT(FT(CTF2) x FT(REF2)*)
            call fftwf_execute_dft_c2r(self%plan_bwd1, self%cvec1(ithr)%c, self%rvec1(ithr)%r)
            self%drvec(ithr)%r(1:self%nrots) = self%drvec(ithr)%r(1:self%nrots) + real(k,dp) * real(self%rvec1(ithr)%r(1:self%nrots),dp)
            ! FT(X.CTF) x FT(REF)*
            if( even )then
                self%cvec1(ithr)%c(1:self%pftsz+1) = self%ft_ptcl_ctf(k,i)%c(1:self%pftsz+1) * self%ft_ref_even(k,iref)%c(1:self%pftsz+1)
            else
                self%cvec1(ithr)%c(1:self%pftsz+1) = self%ft_ptcl_ctf(k,i)%c(1:self%pftsz+1) * self%ft_ref_odd(k,iref)%c(1:self%pftsz+1)
            endif
            ! IFFT( FT(X.CTF) x FT(REF)* )
            call fftwf_execute_dft_c2r(self%plan_bwd1, self%cvec1(ithr)%c, self%rvec1(ithr)%r)
            self%heap_vars(ithr)%kcorrs(1:self%nrots) = self%heap_vars(ithr)%kcorrs(1:self%nrots) + real(k,dp) * real(self%rvec1(ithr)%r(1:self%nrots),dp)
        end do
        self%drvec(ithr)%r(1:self%nrots) = self%drvec(ithr)%r(1:self%nrots) * (self%sqsums_ptcls(i) * real(2*self%nrots,dp))
        corrs = real(self%heap_vars(ithr)%kcorrs(1:self%nrots) / dsqrt(self%drvec(ithr)%r(1:self%nrots)))
    end subroutine gencorrs_cc_dev

    subroutine gencorrs_shifted_cc_dev( self, pft_ref, iptcl, iref, shift, corrs)
        class(polarft_corrcalc), intent(inout) :: self
        complex(sp),             intent(in)    :: pft_ref(1:self%pftsz,self%kfromto(1):self%kfromto(2))
        integer,                 intent(in)    :: iptcl, iref
        real(sp),                intent(in)    :: shift(2)
        real(sp),                intent(out)   :: corrs(self%nrots)
        integer  :: k, i, ithr
        logical  :: even
        ithr = omp_get_thread_num() + 1
        i    = self%pinds(iptcl)
        even = self%iseven(i)
        self%heap_vars(ithr)%kcorrs = 0.d0
        self%drvec(ithr)%r          = 0.d0
        do k = self%kfromto(1),self%kfromto(2)
            ! FT(CTF2) x FT(REF2)), REF2 is shift invariant
            if( even )then
                self%cvec1(ithr)%c = self%ft_ctf2(k,i)%c * self%ft_ref2_even(k,iref)%c
            else
                self%cvec1(ithr)%c = self%ft_ctf2(k,i)%c * self%ft_ref2_odd(k,iref)%c
            endif
            ! IFFT(FT(CTF2) x FT(REF2))
            call fftwf_execute_dft_c2r(self%plan_bwd1, self%cvec1(ithr)%c, self%rvec1(ithr)%r)
            self%drvec(ithr)%r = self%drvec(ithr)%r + real(k,dp) * real(self%rvec1(ithr)%r(1:self%nrots),dp)
            ! FT(S.REF), shifted reference
            self%cvec2(ithr)%c(1:self%pftsz)            = pft_ref(:,k)
            self%cvec2(ithr)%c(self%pftsz+1:self%nrots) = conjg(pft_ref(:,k))
            call fftwf_execute_dft(self%plan_fwd1, self%cvec2(ithr)%c, self%cvec2(ithr)%c)
            ! FT(X.CTF) x FT(S.REF)*
            self%cvec1(ithr)%c = self%ft_ptcl_ctf(k,i)%c * conjg(self%cvec2(ithr)%c(1:self%pftsz+1))
            ! IFFT(FT(X.CTF) x FT(S.REF)*)
            call fftwf_execute_dft_c2r(self%plan_bwd1, self%cvec1(ithr)%c, self%rvec1(ithr)%r)
            self%heap_vars(ithr)%kcorrs = self%heap_vars(ithr)%kcorrs + real(k,dp) * real(self%rvec1(ithr)%r(1:self%nrots),dp)
        end do
        self%drvec(ithr)%r = self%drvec(ithr)%r * real(self%sqsums_ptcls(i) * real(2*self%nrots),dp)
        corrs = real(self%heap_vars(ithr)%kcorrs / dsqrt(self%drvec(ithr)%r))
    end subroutine gencorrs_shifted_cc_dev

    subroutine gencorrs_euclid_dev( self, iptcl, iref, euclids )
        class(polarft_corrcalc), intent(inout) :: self
        integer,                 intent(in)    :: iptcl, iref
        real(sp),                intent(out)   :: euclids(self%nrots)
        real(dp) :: w, sumsqptcl
        integer  :: k, i, ithr
        logical  :: even
        ithr = omp_get_thread_num() + 1
        i    = self%pinds(iptcl)
        even = self%iseven(i)
        self%heap_vars(ithr)%kcorrs = 0.d0
        do k = self%kfromto(1),self%kfromto(2)
            w         = real(k,dp) / real(self%sigma2_noise(k,iptcl),dp)
            sumsqptcl = sum(real(csq_fast(self%pfts_ptcls(:,k,i)), dp))
            ! FT(CTF2) x FT(REF2)*)
            if( even )then
                self%cvec1(ithr)%c = self%ft_ctf2(k,i)%c * self%ft_ref2_even(k,iref)%c
            else
                self%cvec1(ithr)%c = self%ft_ctf2(k,i)%c * self%ft_ref2_odd(k,iref)%c
            endif
            ! IFFT(FT(CTF2) x FT(REF2)*)
            call fftwf_execute_dft_c2r(self%plan_bwd1, self%cvec1(ithr)%c, self%rvec1(ithr)%r)
            self%drvec(ithr)%r = real(self%rvec1(ithr)%r(1:self%nrots),dp) ! fftw unnormalized |CTF.REF|2
            ! FT(X.CTF) x FT(REF)*
            if( even )then
                self%cvec1(ithr)%c = self%ft_ptcl_ctf(k,i)%c * self%ft_ref_even(k,iref)%c
            else
                self%cvec1(ithr)%c = self%ft_ptcl_ctf(k,i)%c * self%ft_ref_odd(k,iref)%c
            endif
            ! X.CTF.REF = IFFT( FT(X.CTF) x FT(REF)* )
            call fftwf_execute_dft_c2r(self%plan_bwd1, self%cvec1(ithr)%c, self%rvec1(ithr)%r)
            ! k/sig2 x ( |CTF.REF|2 - 2.X.CTF.REF ), fftw normalized
            self%drvec(ithr)%r = (w / real(2*self%nrots,dp)) * (self%drvec(ithr)%r - 2.d0*real(self%rvec1(ithr)%r(1:self%nrots),dp))
            ! k/sig2 x ( |X|2 + |CTF.REF|2 - 2.X.CTF.REF )
            self%heap_vars(ithr)%kcorrs = self%heap_vars(ithr)%kcorrs + w * sumsqptcl + self%drvec(ithr)%r
        end do
        euclids = real( dexp( -self%heap_vars(ithr)%kcorrs / self%wsqsums_ptcls(i) ) )
    end subroutine gencorrs_euclid_dev

    subroutine gencorrs_shifted_euclid_dev( self, pft_ref, iptcl, iref, shift, euclids )
        class(polarft_corrcalc), intent(inout) :: self
        complex(sp),             intent(in)    :: pft_ref(1:self%pftsz,self%kfromto(1):self%kfromto(2))
        integer,                 intent(in)    :: iptcl, iref
        real(sp),                intent(in)    :: shift(2)
        real(sp),                intent(out)   :: euclids(self%nrots)
        real(dp) :: w, sumsqptcl
        integer  :: k, i, ithr
        logical  :: even
        ithr = omp_get_thread_num() + 1
        i    = self%pinds(iptcl)
        even = self%iseven(i)
        self%heap_vars(ithr)%kcorrs = 0.d0
        do k = self%kfromto(1),self%kfromto(2)
            w         = real(k,dp) / real(self%sigma2_noise(k,iptcl),dp)
            sumsqptcl = sum(real(csq_fast(self%pfts_ptcls(:,k,i)), dp))
            ! FT(CTF2) x FT(REF2)*)
            if( even )then
                self%cvec1(ithr)%c = self%ft_ctf2(k,i)%c * self%ft_ref2_even(k,iref)%c
            else
                self%cvec1(ithr)%c = self%ft_ctf2(k,i)%c * self%ft_ref2_odd(k,iref)%c
            endif
            ! IFFT(FT(CTF2) x FT(REF2)*)
            call fftwf_execute_dft_c2r(self%plan_bwd1, self%cvec1(ithr)%c, self%rvec1(ithr)%r)
            self%drvec(ithr)%r = real(self%rvec1(ithr)%r(1:self%nrots),dp) ! fftw unnormalized |CTF.REF|2
            ! FT(S.REF), shifted reference
            self%cvec2(ithr)%c(1:self%pftsz)            = pft_ref(:,k)
            self%cvec2(ithr)%c(self%pftsz+1:self%nrots) = conjg(pft_ref(:,k))
            call fftwf_execute_dft(self%plan_fwd1, self%cvec2(ithr)%c, self%cvec2(ithr)%c)
            ! FT(X.CTF) x FT(REF)*
            self%cvec1(ithr)%c = self%ft_ptcl_ctf(k,i)%c * conjg(self%cvec2(ithr)%c(1:self%pftsz+1))
            ! X.CTF.REF = IFFT( FT(X.CTF) x FT(REF)* )
            call fftwf_execute_dft_c2r(self%plan_bwd1, self%cvec1(ithr)%c, self%rvec1(ithr)%r)
            ! k/sig2 x ( |CTF.REF|2 - 2X.CTF.REF ), fftw normalized
            self%drvec(ithr)%r = (w / real(2*self%nrots,dp)) * (self%drvec(ithr)%r - 2.d0*real(self%rvec1(ithr)%r(1:self%nrots),dp))
            ! k/sig2 x ( |X|2 + |CTF.REF|2 - 2X.CTF.REF )
            self%heap_vars(ithr)%kcorrs = self%heap_vars(ithr)%kcorrs + w * sumsqptcl + self%drvec(ithr)%r
        end do
        euclids = real( dexp( -self%heap_vars(ithr)%kcorrs / self%wsqsums_ptcls(i) ) )
    end subroutine gencorrs_shifted_euclid_dev

    subroutine gencorrs_prob_dev( self, iptcl, iref, prob )
        class(polarft_corrcalc), intent(inout) :: self
        integer,                 intent(in)    :: iptcl, iref
        real(sp),                intent(out)   :: prob(self%nrots)
        real(dp) :: w, sumsqptcl, denom
        integer  :: k, i, ithr
        logical  :: even
        ithr  = omp_get_thread_num() + 1
        i     = self%pinds(iptcl)
        even  = self%iseven(i)
        denom =  self%wsqsums_ptcls(i)
        self%heap_vars(ithr)%kcorrs = 0.d0
        do k = self%kfromto(1),self%kfromto(2)
            w         = real(k,dp) / real(self%sigma2_noise(k,iptcl),dp)
            sumsqptcl = sum(real(csq_fast(self%pfts_ptcls(:,k,i)), dp))
            ! FT(CTF2) x FT(REF2)*)
            if( even )then
                self%cvec1(ithr)%c = self%ft_ctf2(k,i)%c * self%ft_ref2_even(k,iref)%c
            else
                self%cvec1(ithr)%c = self%ft_ctf2(k,i)%c * self%ft_ref2_odd(k,iref)%c
            endif
            ! IFFT(FT(CTF2) x FT(REF2)*)
            call fftwf_execute_dft_c2r(self%plan_bwd1, self%cvec1(ithr)%c, self%rvec1(ithr)%r)
            self%drvec(ithr)%r = real(self%rvec1(ithr)%r(1:self%nrots),dp) ! fftw unnormalized |CTF.REF|2
            ! FT(X.CTF) x FT(REF)*
            if( even )then
                self%cvec1(ithr)%c = self%ft_ptcl_ctf(k,i)%c * self%ft_ref_even(k,iref)%c
            else
                self%cvec1(ithr)%c = self%ft_ptcl_ctf(k,i)%c * self%ft_ref_odd(k,iref)%c
            endif
            ! X.CTF.REF = IFFT( FT(X.CTF) x FT(REF)* )
            call fftwf_execute_dft_c2r(self%plan_bwd1, self%cvec1(ithr)%c, self%rvec1(ithr)%r)
            ! k/sig2 x ( |CTF.REF|2 - 2X.CTF.REF), fftw normalized
            self%drvec(ithr)%r = (w / real(2*self%nrots,dp)) * (self%drvec(ithr)%r - 2.d0*self%rvec1(ithr)%r(1:self%nrots))
            ! exp( -k/sig2 x ( |X|2 + |CTF.REF|2 - 2X.CTF.REF )  / SUM(exp( -k/sig2 x |X|2) )
            self%heap_vars(ithr)%kcorrs = self%heap_vars(ithr)%kcorrs + dexp( -(w * sumsqptcl + self%drvec(ithr)%r) / denom)
        end do
        prob = real( self%heap_vars(ithr)%kcorrs / real(self%nk,dp) )
    end subroutine gencorrs_prob_dev

    subroutine gencorrs_shifted_prob_dev( self, pft_ref, iptcl, iref, shift, prob )
        class(polarft_corrcalc), intent(inout) :: self
        complex(sp),             intent(in)    :: pft_ref(1:self%pftsz,self%kfromto(1):self%kfromto(2))
        integer,                 intent(in)    :: iptcl, iref
        real(sp),                intent(in)    :: shift(2)
        real(sp),                intent(out)   :: prob(self%nrots)
        real(dp) :: w, sumsqptcl, denom
        integer  :: k, i, ithr
        logical  :: even
        ithr  = omp_get_thread_num() + 1
        i     = self%pinds(iptcl)
        even  = self%iseven(i)
        denom = self%wsqsums_ptcls(i)
        self%heap_vars(ithr)%kcorrs = 0.d0
        do k = self%kfromto(1),self%kfromto(2)
            w         = real(k,dp) / real(self%sigma2_noise(k,iptcl),dp)
            sumsqptcl = sum(real(csq_fast(self%pfts_ptcls(:,k,i)), dp))
            ! FT(CTF2) x FT(REF2)*)
            if( even )then
                self%cvec1(ithr)%c = self%ft_ctf2(k,i)%c * self%ft_ref2_even(k,iref)%c
            else
                self%cvec1(ithr)%c = self%ft_ctf2(k,i)%c * self%ft_ref2_odd(k,iref)%c
            endif
            ! IFFT(FT(CTF2) x FT(REF2)*)
            call fftwf_execute_dft_c2r(self%plan_bwd1, self%cvec1(ithr)%c, self%rvec1(ithr)%r)
            self%drvec(ithr)%r = real(self%rvec1(ithr)%r(1:self%nrots),dp) ! fftw unnormalized |CTF.REF|2
            ! FT(S.REF), shifted reference
            self%cvec2(ithr)%c(1:self%pftsz)            = pft_ref(:,k)
            self%cvec2(ithr)%c(self%pftsz+1:self%nrots) = conjg(pft_ref(:,k))
            call fftwf_execute_dft(self%plan_fwd1, self%cvec2(ithr)%c, self%cvec2(ithr)%c)
            ! FT(X.CTF) x FT(REF)*
            self%cvec1(ithr)%c = self%ft_ptcl_ctf(k,i)%c * conjg(self%cvec2(ithr)%c(1:self%pftsz+1))
            ! X.CTF.REF = IFFT( FT(X.CTF) x FT(REF)* )
            call fftwf_execute_dft_c2r(self%plan_bwd1, self%cvec1(ithr)%c, self%rvec1(ithr)%r)
            ! k/sig2 x ( |CTF.REF|2 - 2X.CTF.REF ), fftw normalized
            self%drvec(ithr)%r = (w / real(2*self%nrots,dp)) * (self%drvec(ithr)%r - 2.d0*self%rvec1(ithr)%r(1:self%nrots))
            ! exp( -k/sig2 x ( |X|2 + |CTF.REF|2 - 2X.CTF.REF )  / SUM(exp( -k/sig2 x |X|2) )
            self%heap_vars(ithr)%kcorrs = self%heap_vars(ithr)%kcorrs + dexp( -(w * sumsqptcl + self%drvec(ithr)%r) / denom )
        end do
        prob = real( self%heap_vars(ithr)%kcorrs / real(self%nk,dp) )
    end subroutine gencorrs_shifted_prob_dev

    real(dp) function gencorr_for_rot_8_dev( self, iref, iptcl, shvec, irot )
        class(polarft_corrcalc), intent(inout) :: self
        integer,                 intent(in)    :: iref, iptcl
        real(dp),                intent(in)    :: shvec(2)
        integer,                 intent(in)    :: irot
        complex(dp), pointer :: pft_ref_8(:,:), pft_ref_tmp_8(:,:), shmat_8(:,:)
        integer              :: ithr, i
        i    =  self%pinds(iptcl)
        ithr = omp_get_thread_num() + 1
        pft_ref_8     => self%heap_vars(ithr)%pft_ref_8
        pft_ref_tmp_8 => self%heap_vars(ithr)%pft_ref_tmp_8
        shmat_8       => self%heap_vars(ithr)%shmat_8
        if( self%iseven(i) )then
            pft_ref_8 = self%pfts_refs_even(:,:,iref)
        else
            pft_ref_8 = self%pfts_refs_odd(:,:,iref)
        endif
        ! shift
        call self%gen_shmat_8(ithr, shvec, shmat_8)
        pft_ref_8 = pft_ref_8 * shmat_8
        ! rotation
        call self%rotate_ref(pft_ref_8, irot, pft_ref_tmp_8)
        ! ctf
        if( self%with_ctf ) pft_ref_tmp_8 = pft_ref_tmp_8 * self%ctfmats(:,:,i)
        gencorr_for_rot_8_dev = 0.d0
        select case(params_glob%cc_objfun)
            case(OBJFUN_CC)
                gencorr_for_rot_8_dev = self%gencorr_cc_for_rot_8_dev(pft_ref_tmp_8, i)
            case(OBJFUN_EUCLID)
                gencorr_for_rot_8_dev = self%gencorr_euclid_for_rot_8_dev(pft_ref_tmp_8, iptcl)
            case(OBJFUN_PROB)
                gencorr_for_rot_8_dev = self%gencorr_prob_for_rot_8_dev(pft_ref_tmp_8, iptcl)
        end select
    end function gencorr_for_rot_8_dev

    real(dp) function gencorr_cc_for_rot_8_dev( self, pft_ref, i )
        class(polarft_corrcalc), intent(inout) :: self
        complex(dp), pointer,    intent(inout) :: pft_ref(:,:)
        integer,                 intent(in)    :: i
        real(dp) :: sqsum_ref, denom
        integer  :: k
        sqsum_ref  = 0.d0
        gencorr_cc_for_rot_8_dev = 0.d0
        do k = self%kfromto(1),self%kfromto(2)
            sqsum_ref                = sqsum_ref +                real(k,kind=dp) * sum(real(pft_ref(:,k) * conjg(pft_ref(:,k)),dp))
            gencorr_cc_for_rot_8_dev = gencorr_cc_for_rot_8_dev + real(k,kind=dp) * sum(real(pft_ref(:,k) * conjg(self%pfts_ptcls(:,k,i)),dp))
        end do
        denom = sqsum_ref * self%sqsums_ptcls(i)
        if( denom > DTINY )then
            gencorr_cc_for_rot_8_dev = gencorr_cc_for_rot_8_dev / dsqrt(denom)
        else
            gencorr_cc_for_rot_8_dev = 0.d0
        endif
    end function gencorr_cc_for_rot_8_dev

    real(dp) function gencorr_euclid_for_rot_8_dev( self, pft_ref, iptcl )
        class(polarft_corrcalc), intent(inout) :: self
        complex(dp), pointer,    intent(inout) :: pft_ref(:,:)
        integer,                 intent(in)    :: iptcl
        integer  :: i,k
        i       = self%pinds(iptcl)
        pft_ref = pft_ref - self%pfts_ptcls(:,:,i)
        gencorr_euclid_for_rot_8_dev  = 0.d0
        do k = self%kfromto(1),self%kfromto(2)
            gencorr_euclid_for_rot_8_dev = gencorr_euclid_for_rot_8_dev +&
            &(real(k,dp) / self%sigma2_noise(k,iptcl)) * sum(real(csq_fast(pft_ref(:,k)),dp))
        end do
        gencorr_euclid_for_rot_8_dev = dexp( -gencorr_euclid_for_rot_8_dev / self%wsqsums_ptcls(i) )
    end function gencorr_euclid_for_rot_8_dev

    real(dp) function gencorr_prob_for_rot_8_dev( self, pft_ref, iptcl )
        class(polarft_corrcalc), intent(inout) :: self
        complex(dp), pointer,    intent(inout) :: pft_ref(:,:)
        integer,                 intent(in)    :: iptcl
        real(dp) :: diffsq, denom
        integer  :: i,k
        i       = self%pinds(iptcl)
        denom   = self%wsqsums_ptcls(i)
        pft_ref = pft_ref - self%pfts_ptcls(:,:,i)
        gencorr_prob_for_rot_8_dev = 0.d0
        do k = self%kfromto(1),self%kfromto(2)
            diffsq = (real(k,dp) / self%sigma2_noise(k,iptcl)) * sum(real(csq_fast(pft_ref(:,k)),dp))
            gencorr_prob_for_rot_8_dev   = gencorr_prob_for_rot_8_dev + dexp(-diffsq / denom)
        end do
        gencorr_prob_for_rot_8_dev = gencorr_prob_for_rot_8_dev / real(self%nk,dp)
    end function gencorr_prob_for_rot_8_dev

    subroutine gencorr_grad_for_rot_8_dev( self, iref, iptcl, shvec, irot, f, grad )
        class(polarft_corrcalc), intent(inout) :: self
        integer,                 intent(in)    :: iref, iptcl
        real(dp),                intent(in)    :: shvec(2)
        integer,                 intent(in)    :: irot
        real(dp),                intent(out)   :: f, grad(2)
        complex(dp), pointer :: pft_ref_8(:,:), shmat_8(:,:), pft_ref_tmp_8(:,:)
        integer              :: ithr, i
        i    =  self%pinds(iptcl)
        ithr = omp_get_thread_num() + 1
        pft_ref_8     => self%heap_vars(ithr)%pft_ref_8
        pft_ref_tmp_8 => self%heap_vars(ithr)%pft_ref_tmp_8
        shmat_8       => self%heap_vars(ithr)%shmat_8
        if( self%iseven(i) )then
            pft_ref_8 = self%pfts_refs_even(:,:,iref)
        else
            pft_ref_8 = self%pfts_refs_odd(:,:,iref)
        endif
        call self%gen_shmat_8(ithr, shvec, shmat_8)
        pft_ref_8 = pft_ref_8 * shmat_8
        select case(params_glob%cc_objfun)
            case(OBJFUN_CC)
                call self%gencorr_cc_grad_for_rot_8_dev(    pft_ref_8, pft_ref_tmp_8, iptcl, irot, f, grad)
            case(OBJFUN_EUCLID)
                call self%gencorr_euclid_grad_for_rot_8_dev(pft_ref_8, pft_ref_tmp_8, iptcl, irot, f, grad)
            case(OBJFUN_PROB)
                call self%gencorr_prob_grad_for_rot_8_dev(  pft_ref_8, pft_ref_tmp_8, iptcl, irot, f, grad)
        end select
    end subroutine gencorr_grad_for_rot_8_dev

    subroutine gencorr_cc_grad_for_rot_8_dev( self, pft_ref, pft_ref_tmp, iptcl, irot, f, grad )
        class(polarft_corrcalc), intent(inout) :: self
        complex(dp), pointer,    intent(inout) :: pft_ref(:,:), pft_ref_tmp(:,:)
        integer,                 intent(in)    :: iptcl, irot
        real(dp),                intent(out)   :: f, grad(2)
        real(dp) :: sqsum_ref, denom
        integer  :: k, i
        i           = self%pinds(iptcl)
        sqsum_ref   = 0.d0
        f           = 0.d0
        grad        = 0.d0
        if( self%with_ctf )then
            call self%rotate_ref(pft_ref, irot, pft_ref_tmp)
            do k = self%kfromto(1),self%kfromto(2)
                sqsum_ref = sqsum_ref + real(k,kind=dp) * sum(real(self%ctfmats(:,k,i)*self%ctfmats(:,k,i) * pft_ref_tmp(:,k) * conjg(pft_ref_tmp(:,k)),dp))
                f         = f         + real(k,kind=dp) * sum(real(self%ctfmats(:,k,i)                     * pft_ref_tmp(:,k) * conjg(self%pfts_ptcls(:,k,i)),dp))
            enddo
            call self%rotate_ref(pft_ref * dcmplx(0.d0,self%argtransf(:self%pftsz,:)), irot, pft_ref_tmp)
            do k = self%kfromto(1),self%kfromto(2)
                grad(1) = grad(1) + real(k,kind=dp) * sum(real(self%ctfmats(:,k,i) * pft_ref_tmp(:,k) * conjg(self%pfts_ptcls(:,k,i)),dp))
            enddo
            call self%rotate_ref(pft_ref * dcmplx(0.d0,self%argtransf(self%pftsz+1:,:)), irot, pft_ref_tmp)
            do k = self%kfromto(1),self%kfromto(2)
                grad(2) = grad(2) + real(k,kind=dp) * sum(real(self%ctfmats(:,k,i) * pft_ref_tmp(:,k) * conjg(self%pfts_ptcls(:,k,i)),dp))
            end do
        else
            call self%rotate_ref(pft_ref, irot, pft_ref_tmp)
            do k = self%kfromto(1),self%kfromto(2)
                sqsum_ref = sqsum_ref + real(k,kind=dp) * sum(real(pft_ref_tmp(:,k) * conjg(pft_ref_tmp(:,k)),dp))
                f         = f         + real(k,kind=dp) * sum(real(pft_ref_tmp(:,k) * conjg(self%pfts_ptcls(:,k,i)),dp))
            enddo
            call self%rotate_ref(pft_ref * dcmplx(0.d0,self%argtransf(:self%pftsz,:)), irot, pft_ref_tmp)
            do k = self%kfromto(1),self%kfromto(2)
                grad(1) = grad(1) + real(k,kind=dp) * sum(real(pft_ref_tmp(:,k) * conjg(self%pfts_ptcls(:,k,i)),dp))
            enddo
            call self%rotate_ref(pft_ref * dcmplx(0.d0,self%argtransf(self%pftsz+1:,:)), irot, pft_ref_tmp)
            do k = self%kfromto(1),self%kfromto(2)
                grad(2) = grad(2) + real(k,kind=dp) * sum(real(pft_ref_tmp(:,k) * conjg(self%pfts_ptcls(:,k,i)),dp))
            end do
        endif
        denom = dsqrt(sqsum_ref * self%sqsums_ptcls(i))
        f     = f    / denom
        grad  = grad / denom
    end subroutine gencorr_cc_grad_for_rot_8_dev

    subroutine gencorr_euclid_grad_for_rot_8_dev( self, pft_ref, pft_ref_tmp, iptcl, irot, f, grad )
        class(polarft_corrcalc), intent(inout) :: self
        complex(dp), pointer,    intent(inout) :: pft_ref(:,:), pft_ref_tmp(:,:)
        integer,                 intent(in)    :: iptcl, irot
        real(dp),                intent(out)   :: f, grad(2)
        complex(dp), pointer :: pft_diff(:,:)
        real(dp) :: denom, w
        integer  :: k, i, ithr
        ithr        = omp_get_thread_num() + 1
        i           = self%pinds(iptcl)
        f           = 0.d0
        grad        = 0.d0
        denom       = self%wsqsums_ptcls(i)
        pft_diff => self%heap_vars(ithr)%shmat_8
        call self%rotate_ref(pft_ref, irot, pft_ref_tmp)
        if( self%with_ctf ) pft_ref_tmp = pft_ref_tmp * self%ctfmats(:,:,i)
        pft_diff = pft_ref_tmp - self%pfts_ptcls(:,:,i) ! Ref(shift + rotation + CTF) - Ptcl
        call self%rotate_ref(pft_ref * dcmplx(0.d0,self%argtransf(:self%pftsz,:)), irot, pft_ref_tmp)
        if( self%with_ctf ) pft_ref_tmp = pft_ref_tmp * self%ctfmats(:,:,i)
        do k = self%kfromto(1),self%kfromto(2)
            w       = real(k,dp) / real(self%sigma2_noise(k,iptcl))
            f       = f + w * sum(real(csq_fast(pft_diff(:,k)),dp))
            grad(1) = grad(1) + w * real(sum(pft_ref_tmp(:,k) * conjg(pft_diff(:,k))),dp)
        end do
        call self%rotate_ref(pft_ref * dcmplx(0.d0,self%argtransf(self%pftsz+1:,:)), irot, pft_ref_tmp)
        if( self%with_ctf ) pft_ref_tmp = pft_ref_tmp * self%ctfmats(:,:,i)
        do k = self%kfromto(1),self%kfromto(2)
            w      = real(k,dp) / real(self%sigma2_noise(k,iptcl))
            grad(2) = grad(2) + w * real(sum(pft_ref_tmp(:,k) * conjg(pft_diff(:,k))),dp)
        end do
        f    = dexp( -f / denom )
        grad = -f * 2.d0 * grad / denom
    end subroutine gencorr_euclid_grad_for_rot_8_dev

    subroutine gencorr_prob_grad_for_rot_8_dev( self, pft_ref, pft_ref_tmp, iptcl, irot, f, grad )
        class(polarft_corrcalc), intent(inout) :: self
        complex(dp), pointer,    intent(inout) :: pft_ref(:,:), pft_ref_tmp(:,:)
        integer,                 intent(in)    :: iptcl, irot
        real(dp),                intent(out)   :: f, grad(2)
        complex(dp), pointer :: pft_diff(:,:)
        real(dp) :: expdiffsq, gradsq, denom, w
        integer  :: k, i, ithr
        ithr     = omp_get_thread_num() + 1
        i        = self%pinds(iptcl)
        f        = 0.d0
        grad     = 0.d0
        denom    = self%wsqsums_ptcls(i)
        pft_diff => self%heap_vars(ithr)%shmat_8
        call self%rotate_ref(pft_ref, irot, pft_ref_tmp)
        if( self%with_ctf ) pft_ref_tmp = pft_ref_tmp * self%ctfmats(:,:,i)
        pft_diff = pft_ref_tmp - self%pfts_ptcls(:,:,i) ! Ref(shift + rotation + CTF) - Ptcl
        call self%rotate_ref(pft_ref * dcmplx(0.d0,self%argtransf(:self%pftsz,:)), irot, pft_ref_tmp)
        if( self%with_ctf ) pft_ref_tmp = pft_ref_tmp * self%ctfmats(:,:,i)
        do k = self%kfromto(1),self%kfromto(2)
            w         = real(k,dp) / real(self%sigma2_noise(k,iptcl))
            expdiffsq =  dexp(-w * sum(real(csq_fast(pft_diff(:,k)),dp)) / denom)
            f         = f + expdiffsq
            gradsq    = real(sum(pft_ref_tmp(:,k) * conjg(pft_diff(:,k))),dp)
            grad(1)   = grad(1) - expdiffsq * 2.d0 * w*gradsq/denom
        end do
        call self%rotate_ref(pft_ref * dcmplx(0.d0,self%argtransf(self%pftsz+1:,:)), irot, pft_ref_tmp)
        if( self%with_ctf ) pft_ref_tmp = pft_ref_tmp * self%ctfmats(:,:,i)
        do k = self%kfromto(1),self%kfromto(2)
            w         = real(k,dp) / real(self%sigma2_noise(k,iptcl))
            expdiffsq =  dexp(-w * sum(real(csq_fast(pft_diff(:,k)),dp)) / denom)
            gradsq    = real(sum(pft_ref_tmp(:,k) * conjg(pft_diff(:,k))),dp)
            grad(2)   = grad(2) - expdiffsq * 2.d0 * w*gradsq/denom
        end do
        f    = f    / real(self%nk,dp)
        grad = grad / real(self%nk,dp)
    end subroutine gencorr_prob_grad_for_rot_8_dev

    subroutine gencorr_grad_only_for_rot_8_dev( self, iref, iptcl, shvec, irot, grad )
        class(polarft_corrcalc), intent(inout) :: self
        integer,                 intent(in)    :: iref, iptcl
        real(dp),                intent(in)    :: shvec(2)
        integer,                 intent(in)    :: irot
        real(dp),                intent(out)   :: grad(2)
        complex(dp), pointer :: pft_ref_8(:,:), shmat_8(:,:), pft_ref_tmp_8(:,:)
        real(dp) :: f
        integer  :: ithr, i
        i    =  self%pinds(iptcl)
        ithr = omp_get_thread_num() + 1
        pft_ref_8     => self%heap_vars(ithr)%pft_ref_8
        pft_ref_tmp_8 => self%heap_vars(ithr)%pft_ref_tmp_8
        shmat_8       => self%heap_vars(ithr)%shmat_8
        if( self%iseven(i) )then
            pft_ref_8 = self%pfts_refs_even(:,:,iref)
        else
            pft_ref_8 = self%pfts_refs_odd(:,:,iref)
        endif
        call self%gen_shmat_8(ithr, shvec, shmat_8)
        pft_ref_8 = pft_ref_8 * shmat_8
        select case(params_glob%cc_objfun)
            case(OBJFUN_CC)
                call self%gencorr_cc_grad_only_for_rot_8_dev(pft_ref_8, pft_ref_tmp_8, i, irot, grad)
            case(OBJFUN_EUCLID)
                call self%gencorr_euclid_grad_for_rot_8_dev(pft_ref_8, pft_ref_tmp_8, iptcl, irot, f, grad)
            case(OBJFUN_PROB)
                call self%gencorr_prob_grad_for_rot_8_dev(  pft_ref_8, pft_ref_tmp_8, iptcl, irot, f, grad)
        end select
    end subroutine gencorr_grad_only_for_rot_8_dev

    subroutine gencorr_cc_grad_only_for_rot_8_dev( self, pft_ref, pft_ref_tmp, i, irot, grad )
        class(polarft_corrcalc), intent(inout) :: self
        complex(dp), pointer,    intent(inout) :: pft_ref(:,:), pft_ref_tmp(:,:)
        integer,                 intent(in)    :: i, irot
        real(dp),                intent(out)   :: grad(2)
        real(dp) :: sqsum_ref
        integer  :: k
        sqsum_ref = 0.d0
        grad      = 0.d0
        if( self%with_ctf )then
            call self%rotate_ref(pft_ref, irot, pft_ref_tmp)
            do k = self%kfromto(1),self%kfromto(2)
                sqsum_ref = sqsum_ref + real(k,kind=dp) * sum(real(self%ctfmats(:,k,i)*self%ctfmats(:,k,i) * pft_ref_tmp(:,k) * conjg(pft_ref_tmp(:,k)),dp))
            enddo
            call self%rotate_ref(pft_ref * dcmplx(0.d0,self%argtransf(:self%pftsz,:)), irot, pft_ref_tmp)
            do k = self%kfromto(1),self%kfromto(2)
                grad(1) = grad(1) + real(k,kind=dp) * sum(real(self%ctfmats(:,k,i) * pft_ref_tmp(:,k) * conjg(self%pfts_ptcls(:,k,i)),dp))
            enddo
            call self%rotate_ref(pft_ref * dcmplx(0.d0,self%argtransf(self%pftsz+1:,:)), irot, pft_ref_tmp)
            do k = self%kfromto(1),self%kfromto(2)
                grad(2) = grad(2) + real(k,kind=dp) * sum(real(self%ctfmats(:,k,i) * pft_ref_tmp(:,k) * conjg(self%pfts_ptcls(:,k,i)),dp))
            end do
        else
            call self%rotate_ref(pft_ref, irot, pft_ref_tmp)
            do k = self%kfromto(1),self%kfromto(2)
                sqsum_ref = sqsum_ref + real(k,kind=dp) * sum(real(pft_ref_tmp(:,k) * conjg(pft_ref_tmp(:,k)),dp))
            enddo
            call self%rotate_ref(pft_ref * dcmplx(0.d0,self%argtransf(:self%pftsz,:)), irot, pft_ref_tmp)
            do k = self%kfromto(1),self%kfromto(2)
                grad(1) = grad(1) + real(k,kind=dp) * sum(real(pft_ref_tmp(:,k) * conjg(self%pfts_ptcls(:,k,i)),dp))
            enddo
            call self%rotate_ref(pft_ref * dcmplx(0.d0,self%argtransf(self%pftsz+1:,:)), irot, pft_ref_tmp)
            do k = self%kfromto(1),self%kfromto(2)
                grad(2) = grad(2) + real(k,kind=dp) * sum(real(pft_ref_tmp(:,k) * conjg(self%pfts_ptcls(:,k,i)),dp))
            end do
        endif
        grad = grad / dsqrt(sqsum_ref * self%sqsums_ptcls(i))
    end subroutine gencorr_cc_grad_only_for_rot_8_dev

    subroutine gencorr_sigma_contrib_dev( self, iref, iptcl, shvec, irot, sigma_contrib)
        class(polarft_corrcalc), intent(inout) :: self
        integer,                 intent(in)    :: iref, iptcl
        real(sp),                intent(in)    :: shvec(2)
        integer,                 intent(in)    :: irot
        real(sp),                intent(out)   :: sigma_contrib(self%kfromto(1):self%kfromto(2))
        complex(dp), pointer :: pft_ref_8(:,:), shmat_8(:,:), pft_ref_tmp_8(:,:)
        integer :: i,ithr
        i    =  self%pinds(iptcl)
        ithr = omp_get_thread_num() + 1
        pft_ref_8     => self%heap_vars(ithr)%pft_ref_8
        pft_ref_tmp_8 => self%heap_vars(ithr)%pft_ref_tmp_8
        shmat_8       => self%heap_vars(ithr)%shmat_8
        ! e/o
        if( self%iseven(i) )then
            pft_ref_8 = self%pfts_refs_even(:,:,iref)
        else
            pft_ref_8 = self%pfts_refs_odd(:,:,iref)
        endif
        ! shift
        call self%gen_shmat_8(ithr, real(shvec,dp), shmat_8)
        pft_ref_8 = pft_ref_8 * shmat_8
        ! rotation
        call self%rotate_ref(pft_ref_8, irot, pft_ref_tmp_8)
        ! ctf
        if( self%with_ctf ) pft_ref_tmp_8 = pft_ref_tmp_8 * real(self%ctfmats(:,:,i),dp)
        ! difference
        pft_ref_tmp_8 = pft_ref_tmp_8 - self%pfts_ptcls(:,:,i)
        ! sigma2
        sigma_contrib = real(sum(real(pft_ref_tmp_8 * conjg(pft_ref_tmp_8),dp), dim=1) / (2.d0*real(self%pftsz,dp)))
    end subroutine gencorr_sigma_contrib_dev

    subroutine calc_frc_dev( self, iref, iptcl, irot, shvec, frc )
        class(polarft_corrcalc),  intent(inout) :: self
        integer,                  intent(in)    :: iref, iptcl, irot
        real(sp),                 intent(in)    :: shvec(2)
        real(sp),                 intent(out)   :: frc(self%kfromto(1):self%kfromto(2))
        complex(dp), pointer :: pft_ref(:,:), shmat(:,:), pft_rot_ref(:,:)
        real(dp) :: sumsqref, sumsqptcl, denom, num
        integer  :: k, ithr, i
        i    =  self%pinds(iptcl)
        ithr = omp_get_thread_num() + 1
        pft_ref     => self%heap_vars(ithr)%pft_ref_8
        pft_rot_ref => self%heap_vars(ithr)%pft_ref_tmp_8
        shmat       => self%heap_vars(ithr)%shmat_8
        if( self%iseven(i) )then
            pft_ref = self%pfts_refs_even(:,:,iref)
        else
            pft_ref = self%pfts_refs_odd(:,:,iref)
        endif
        call self%gen_shmat_8(ithr, real(shvec,dp), shmat)
        pft_ref = pft_ref * shmat
        call self%rotate_ref(pft_ref, irot, pft_rot_ref)
        if( self%with_ctf ) pft_rot_ref = pft_rot_ref * real(self%ctfmats(:,:,i),dp)
        do k = self%kfromto(1),self%kfromto(2)
            num       = real(sum(pft_rot_ref(:,k)       * conjg(self%pfts_ptcls(:,k,i))),dp)
            sumsqptcl = real(sum(self%pfts_ptcls(:,k,i) * conjg(self%pfts_ptcls(:,k,i))),dp)
            sumsqref  = real(sum(pft_rot_ref(:,k)       * conjg(pft_rot_ref(:,k))),dp)
            denom     = sumsqptcl * sumsqref
            if( denom < 1.d-16 )then
                frc(k) = 0.
            else
                frc(k) = real(num / denom)
            endif
        end do
    end subroutine calc_frc_dev

    function gencorr_cont_cc_for_rot_8_dev( self, iref, iptcl, shvec, irot )result(cc)
        class(polarft_corrcalc), intent(inout) :: self
        integer,                 intent(in)    :: iref, iptcl
        real(dp),                intent(in)    :: shvec(2)
        integer,                 intent(in)    :: irot
        complex(dp), pointer :: pft_ref_8(:,:), pft_ref_tmp_8(:,:), shmat_8(:,:)
        real(dp)             :: sqsum_ref, cc
        integer              :: ithr, i, k
        i    =  self%pinds(iptcl)
        ithr = omp_get_thread_num() + 1
        pft_ref_8     => self%heap_vars(ithr)%pft_ref_8
        pft_ref_tmp_8 => self%heap_vars(ithr)%pft_ref_tmp_8
        shmat_8       => self%heap_vars(ithr)%shmat_8
        if( self%iseven(i) )then
            pft_ref_8 = self%pfts_refs_even(:,:,iref)
        else
            pft_ref_8 = self%pfts_refs_odd(:,:,iref)
        endif
        ! shift
        call self%gen_shmat_8(ithr, shvec, shmat_8)
        pft_ref_8 = pft_ref_8 * shmat_8
        ! rotation
        call self%rotate_ref(pft_ref_8, irot, pft_ref_tmp_8)
        ! ctf
        if( self%with_ctf ) pft_ref_tmp_8 = pft_ref_tmp_8 * self%ctfmats(:,:,i)
        ! correlation
        sqsum_ref = 0.0
        cc        = 0.d0
        do k = self%kfromto(1),self%kfromto(2)
            sqsum_ref = sqsum_ref + real(k,kind=dp) * sum(real(pft_ref_tmp_8(:,k) * conjg(pft_ref_tmp_8(:,k)),dp))
            cc        = cc        + real(k,kind=dp) * sum(real(pft_ref_tmp_8(:,k) * conjg(self%pfts_ptcls(:,k,i)),dp))
        end do
        cc = cc / dsqrt(sqsum_ref * self%sqsums_ptcls(i))
    end function gencorr_cont_cc_for_rot_8_dev

    function gencorr_cont_grad_cc_for_rot_8_dev( self, iref, iptcl, shvec, irot, dcc ) result( cc )
        class(polarft_corrcalc), intent(inout) :: self
        integer,                 intent(in)    :: iref, iptcl
        real(dp),                intent(in)    :: shvec(2)
        integer,                 intent(in)    :: irot
        real(dp),                intent(out)   :: dcc(3)
        complex(dp), pointer :: pft_ref_8(:,:), pft_ref_tmp_8(:,:), shmat_8(:,:), pft_dref_8(:,:,:)
        real(dp) :: T1(3), T2(3), cc, sqsum_ref, denom, num
        integer  :: ithr, j, k, i
        i    =  self%pinds(iptcl)
        ithr = omp_get_thread_num() + 1
        pft_ref_8     => self%heap_vars(ithr)%pft_ref_8
        pft_ref_tmp_8 => self%heap_vars(ithr)%pft_ref_tmp_8
        pft_dref_8    => self%heap_vars(ithr)%pft_dref_8
        shmat_8       => self%heap_vars(ithr)%shmat_8
        ! e/o
        if( self%iseven(i) )then
            pft_ref_8  = self%pfts_refs_even(:,:,iref)
            pft_dref_8 = self%pfts_drefs_even(:,:,:,iref)
        else
            pft_ref_8  = self%pfts_refs_odd(:,:,iref)
            pft_dref_8 = self%pfts_drefs_odd(:,:,:,iref)
        endif
        ! shift
        call self%gen_shmat_8(ithr, shvec, shmat_8)
        pft_ref_8 = pft_ref_8 * shmat_8
        do j = 1,3
            pft_dref_8(:,:,j) = pft_dref_8(:,:,j) * shmat_8
        end do
        ! rotation
        call self%rotate_ref(pft_ref_8, irot, pft_ref_tmp_8)
        pft_ref_8 = pft_ref_tmp_8
        do j = 1,3
            call self%rotate_ref(pft_dref_8(:,:,j), irot, pft_ref_tmp_8)
            pft_dref_8(:,:,j) = pft_ref_tmp_8
        end do
        ! ctf
        if( self%with_ctf )then
            pft_ref_8 = pft_ref_8 * self%ctfmats(:,:,i)
            do j = 1,3
                pft_dref_8(:,:,j) = pft_dref_8(:,:,j) * self%ctfmats(:,:,i)
            end do
        endif
        ! correlation & derivatives
        sqsum_ref = 0.d0
        num       = 0.d0
        T1        = 0.d0
        T2        = 0.d0
        do k = self%kfromto(1),self%kfromto(2)
            sqsum_ref  = sqsum_ref + real(k,kind=dp) * sum(real(pft_ref_8(:,k) * conjg(pft_ref_8(:,k)),dp))
            num        = num       + real(k,kind=dp) * sum(real(pft_ref_8(:,k) * conjg(self%pfts_ptcls(:,k,i)),dp))
            do j = 1,3
                T1(j) = T1(j) + real(k,kind=dp) * real(sum(pft_dref_8(:,k,j) * conjg(self%pfts_ptcls(:,k,i))),dp)
                T2(j) = T2(j) + real(k,kind=dp) * real(sum(pft_dref_8(:,k,j) * conjg(pft_ref_8(:,k))),dp)
            enddo
        enddo
        denom = sqrt(sqsum_ref * self%sqsums_ptcls(i))
        cc    = num / denom
        dcc   = (T1 - num * T2 / sqsum_ref) / denom
    end function gencorr_cont_grad_cc_for_rot_8_dev

    subroutine gencorr_cont_shift_grad_cc_for_rot_8_dev( self, iref, iptcl, shvec, irot, f, grad )
        class(polarft_corrcalc), intent(inout) :: self
        integer,                 intent(in)    :: iref, iptcl
        real(dp),                intent(in)    :: shvec(2)
        integer,                 intent(in)    :: irot
        real(dp),                intent(out)   :: f
        real(dp),                intent(out)   :: grad(5) ! 3 orientation angles, 2 shifts
        complex(dp), pointer :: pft_ref_8(:,:), pft_ref_tmp_8(:,:), shmat_8(:,:), pft_dref_8(:,:,:)
        real(dp) :: T1(3), T2(3), sqsum_ref, denom, num
        integer  :: ithr, j, k, i
        i    =  self%pinds(iptcl)
        ithr = omp_get_thread_num() + 1
        pft_ref_8     => self%heap_vars(ithr)%pft_ref_8
        pft_ref_tmp_8 => self%heap_vars(ithr)%pft_ref_tmp_8
        pft_dref_8    => self%heap_vars(ithr)%pft_dref_8
        shmat_8       => self%heap_vars(ithr)%shmat_8
        ! e/o
        if( self%iseven(i) )then
            pft_ref_8  = self%pfts_refs_even(:,:,iref)
            pft_dref_8 = self%pfts_drefs_even(:,:,:,iref)
        else
            pft_ref_8  = self%pfts_refs_odd(:,:,iref)
            pft_dref_8 = self%pfts_drefs_odd(:,:,:,iref)
        endif
        ! shift
        call self%gen_shmat_8(ithr, shvec, shmat_8)
        pft_ref_8 = pft_ref_8 * shmat_8
        do j = 1,3
            pft_dref_8(:,:,j) = pft_dref_8(:,:,j) * shmat_8
        end do
        ! rotation
        call self%rotate_ref(pft_ref_8, irot, pft_ref_tmp_8)
        pft_ref_8 = pft_ref_tmp_8
        do j = 1,3
            call self%rotate_ref(pft_dref_8(:,:,j), irot, pft_ref_tmp_8)
            pft_dref_8(:,:,j) = pft_ref_tmp_8
        end do
        ! ctf
        if( self%with_ctf )then
            pft_ref_8 = pft_ref_8 * self%ctfmats(:,:,i)
            do j = 1,3
                pft_dref_8(:,:,j) = pft_dref_8(:,:,j) * self%ctfmats(:,:,i)
            end do
        endif
        ! correlation & orientation derivatives
        sqsum_ref = 0.d0
        num       = 0.d0
        T1        = 0.d0
        T2        = 0.d0
        grad      = 0.d0
        do k = self%kfromto(1),self%kfromto(2)
            sqsum_ref  = sqsum_ref + real(k,kind=dp) * sum(real(pft_ref_8(:,k) * conjg(pft_ref_8(:,k)),dp))
            num        = num       + real(k,kind=dp) * sum(real(pft_ref_8(:,k) * conjg(self%pfts_ptcls(:,k,i)),dp))
            do j = 1,3
                T1(j) = T1(j) + real(k,kind=dp) * real(sum(pft_dref_8(:,k,j) * conjg(self%pfts_ptcls(:,k,i))),dp)
                T2(j) = T2(j) + real(k,kind=dp) * real(sum(pft_dref_8(:,k,j) * conjg(pft_ref_8(:,k))),dp)
            enddo
        enddo
        denom     = sqrt(sqsum_ref * self%sqsums_ptcls(i))
        f         = num / denom
        grad(1:3) = (T1 - num * T2 / sqsum_ref) / denom
        ! shift derivatives
        pft_ref_tmp_8 = pft_ref_8 * (0.d0, 1.d0) * self%argtransf(:self%pftsz,:)
        num = 0.d0
        do k = self%kfromto(1),self%kfromto(2)
            num = num + real(k,kind=dp) * sum(real(pft_ref_tmp_8(:,k) * conjg(self%pfts_ptcls(:,k,i)),dp))
        enddo
        grad(4) = num / denom
        pft_ref_tmp_8 = pft_ref_8 * (0.d0, 1.d0) * self%argtransf(self%pftsz+1:,:)
        num = 0.d0
        do k = self%kfromto(1),self%kfromto(2)
            num = num + real(k,kind=dp) * sum(real(pft_ref_tmp_8(:,k) * conjg(self%pfts_ptcls(:,k,i)),dp))
        enddo
        grad(5) = num / denom
    end subroutine gencorr_cont_shift_grad_cc_for_rot_8_dev

end module simple_polarft_corrcalc
