! for calculation of band-pass limited cross-correlation of polar Fourier transforms
module simple_polarft_corrcalc
!$ use omp_lib
!$ use omp_lib_kinds
include 'simple_lib.f08'
use simple_fftw3
use simple_parameters, only: params_glob
implicit none

public :: polarft_corrcalc, pftcc_glob
private
#include "simple_local_flags.inc"

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
    complex(sp), pointer :: pft_dref(:,:,:)     => null()
    real,        pointer :: corrs_over_k(:)     => null()
    real(sp),    pointer :: argmat(:,:)         => null()
    complex(sp), pointer :: shmat(:,:)          => null()
    real,        pointer :: kcorrs(:)           => null()
    complex(dp), pointer :: pft_ref_8(:,:)      => null()
    complex(dp), pointer :: pft_ref_tmp_8(:,:)  => null()
    complex(dp), pointer :: pft_ref_tmp1_8(:,:) => null()
    complex(dp), pointer :: pft_ref_tmp2_8(:,:) => null()
    complex(dp), pointer :: pft_dref_8(:,:,:)   => null()
    complex(dp), pointer :: shmat_8(:,:)        => null()
    real(dp),    pointer :: argmat_8(:,:)       => null()
    real(dp),    pointer :: fdf_y_8(:)          => null()
    real(dp),    pointer :: fdf_T1_8(:,:)       => null()
    real(dp),    pointer :: fdf_T2_8(:,:)       => null()
end type heap_vars

type :: polarft_corrcalc
    private
    integer                          :: hpind_fsc  = 0        !< high-pass index FSC
    integer                          :: nptcls     = 1        !< the total number of particles in partition (logically indexded [fromp,top])
    integer                          :: nrefs      = 1        !< the number of references (logically indexded [1,nrefs])
    integer                          :: nrots      = 0        !< number of in-plane rotations for one pft (determined by radius of molecule)
    integer                          :: pftsz      = 0        !< size of reference and particle pft (nrots/2)
    integer                          :: nthr       = 0        !< # OpenMP threads
    integer                          :: pfromto(2) = 0        !< particle index range
    integer                          :: winsz      = 0        !< size of moving window in correlation calculations
    integer                          :: ldim(3)    = 0        !< logical dimensions of original cartesian image
    integer,             allocatable :: pinds(:)              !< index array (to reduce memory when frac_update < 1)
    real(sp),            allocatable :: ptcl_bfac_weights(:,:)!< B-factor per particle array for weighting of the correlation
    real(sp),            allocatable :: ptcl_bfac_norms(:)    !< normalisation constants for B-factor weighted ccres
    real(sp),            allocatable :: inv_resarrsq(:)       !< memoized -1./(4*res^2) for B-factors calculation
    real(sp),            allocatable :: sqsums_ptcls(:)       !< memoized square sums for the correlation calculations
    real(sp),            allocatable :: angtab(:)             !< table of in-plane angles (in degrees)
    real(sp),            allocatable :: argtransf(:,:)        !< argument transfer constants for shifting the references
    real(sp),            allocatable :: polar(:,:)            !< table of polar coordinates (in Cartesian coordinates)
    real(sp),            allocatable :: ctfmats(:,:,:)        !< expand set of CTF matrices (for efficient parallel exec)
    complex(sp),         allocatable :: pfts_refs_even(:,:,:) !< 3D complex matrix of polar reference sections (nrefs,pftsz,nk), even
    complex(sp),         allocatable :: pfts_refs_odd(:,:,:)  !< -"-, odd
    complex(sp),         allocatable :: pfts_drefs_even(:,:,:,:) !< derivatives w.r.t. orientation angles of 3D complex matrices
    complex(sp),         allocatable :: pfts_drefs_odd(:,:,:,:)  !< derivatives w.r.t. orientation angles of 3D complex matrices
    complex(sp),         allocatable :: pfts_ptcls(:,:,:)     !< 3D complex matrix of particle sections
    complex(sp),         allocatable :: fft_factors(:)        !< phase factors for accelerated gencorrs routines
    type(fftw_arrs),     allocatable :: fftdat(:)             !< arrays for accelerated gencorrs routines
    type(fftw_carr_fft), allocatable :: fftdat_ptcls(:,:)     !< for memoization of particle  FFTs in accelerated gencorrs routines
    type(fftw_carr),     allocatable :: fft_carray(:)         !< for on-the-fly memoization of particle  FFTs
    logical,             allocatable :: iseven(:)             !< eo assignment for gold-standard FSC
    real,                pointer     :: sigma2_noise(:,:)    => null() !< for euclidean distances
    logical,             pointer     :: sigma2_exists_msk(:) => null() !< mask indicating if sigma2_noise exists for particle
    type(c_ptr)                      :: plan_fwd_1            !< FFTW plans for gencorrs
    type(c_ptr)                      :: plan_fwd_2            !< -"-
    type(c_ptr)                      :: plan_bwd              !< -"-
    logical                          :: with_ctf    = .false. !< CTF flag
    logical                          :: existence   = .false. !< to indicate existence
    type(heap_vars),     allocatable :: heap_vars(:)          !< allocated fields to save stack allocation in subroutines and functions
  contains
    ! CONSTRUCTOR
    procedure          :: new
    ! SETTERS
    procedure          :: set_ref_pft
    procedure          :: set_ptcl_pft
    procedure          :: set_ref_fcomp
    procedure          :: set_dref_fcomp
    procedure          :: set_ptcl_fcomp
    procedure          :: cp_even2odd_ref
    procedure          :: cp_even_ref2ptcl
    procedure          :: cp_refs
    procedure          :: swap_ptclsevenodd
    procedure          :: set_eos
    procedure          :: assign_sigma2_noise
    ! GETTERS
    procedure          :: get_nrots
    procedure          :: get_pdim
    procedure          :: get_rot
    procedure          :: get_roind
    procedure          :: get_coord
    procedure          :: get_ref_pft
    procedure          :: get_nrefs
    procedure          :: exists
    procedure          :: ptcl_iseven
    procedure          :: get_nptcls
    procedure          :: assign_pinds
    procedure          :: is_euclid
    ! PRINTERS/VISUALISERS
    procedure          :: print
    procedure          :: vis_ptcl
    procedure          :: vis_ref
    ! MEMOIZER
    procedure, private :: memoize_sqsum_ptcl
    procedure, private :: memoize_fft
    procedure          :: memoize_ptcl_fft
    procedure          :: memoize_ffts
    procedure          :: memoize_bfac
    ! CALCULATORS
    procedure          :: create_polar_absctfmats
    procedure, private :: prep_ref4corr
    procedure, private :: calc_corrs_over_k
    procedure, private :: calc_k_corrs
    procedure, private :: calc_corr_for_rot
    procedure, private :: calc_corr_for_rot_8
    procedure, private :: calc_T1_T2_for_rot_8
    procedure, private :: calc_euclid_for_rot
    procedure, private :: calc_euclid_for_rot_8
    procedure, private :: calc_corrk_for_rot_8
    procedure, private :: calc_euclidk_for_rot
    procedure, private :: gencorrs_cc_1
    procedure, private :: gencorrs_cc_2
    procedure          :: gencorrs_cc_valid
    procedure, private :: gencorrs_resnorm_1
    procedure, private :: gencorrs_resnorm_2
    procedure, private :: gencorrs_euclid_1
    procedure, private :: gencorrs_euclid_2
    procedure, private :: gencorrs_1
    procedure, private :: gencorrs_2
    generic            :: gencorrs => gencorrs_1, gencorrs_2
    procedure          :: gencorr_for_rot_8
    procedure          :: gencorr_grad_for_rot_8
    procedure          :: gencorr_grad_only_for_rot_8
    procedure          :: gencorr_cc_for_rot
    procedure          :: gencorr_cc_for_rot_8
    procedure          :: gencorr_cont_grad_cc_for_rot_8
    procedure          :: gencorr_cont_shift_grad_cc_for_rot_8
    procedure          :: gencorr_cc_grad_for_rot_8
    procedure          :: gencorr_cc_grad_only_for_rot_8
    procedure          :: gencorr_resnorm_for_rot_8
    procedure          :: gencorr_resnorm_grad_for_rot_8
    procedure          :: gencorr_cont_grad_resnorm_for_rot_8
    procedure          :: gencorr_cont_shift_grad_resnorm_for_rot_8
    procedure          :: gencorr_resnorm_grad_only_for_rot_8
    procedure          :: gencorr_euclid_for_rot
    procedure          :: gencorr_euclid_for_rot_8
    procedure          :: gencorr_cont_grad_euclid_for_rot_8
    procedure          :: gencorr_cont_shift_grad_euclid_for_rot_8
    procedure          :: gencorr_euclid_grad_for_rot_8
    procedure          :: gencorr_euclid_grad_only_for_rot_8
    procedure          :: gencorr_sigma_contrib
    procedure, private :: genfrc
    procedure          :: calc_frc
    procedure          :: specscore
    procedure          :: fit_bfac
    procedure          :: calc_roinv_corrmat
    ! DESTRUCTOR
    procedure          :: kill
end type polarft_corrcalc

! CLASS PARAMETERS/VARIABLES
complex(sp), parameter           :: zero=cmplx(0.,0.) !< just a complex zero
integer,     parameter           :: FFTW_USE_WISDOM=16
class(polarft_corrcalc), pointer :: pftcc_glob

contains

    ! CONSTRUCTORS

    !>  \brief  is a constructor
    subroutine new( self, nrefs, pfromto, ptcl_mask, eoarr )
        class(polarft_corrcalc), target, intent(inout) :: self
        integer,                         intent(in)    :: nrefs
        integer,                         intent(in)    :: pfromto(2)
        logical, optional,               intent(in)    :: ptcl_mask(pfromto(1):pfromto(2))
        integer, optional,               intent(in)    :: eoarr(pfromto(1):pfromto(2))
        character(kind=c_char, len=:), allocatable :: fft_wisdoms_fname ! FFTW wisdoms (per part or suffer I/O lag)
        integer             :: local_stat,irot, k, ithr, i, ik, cnt
        logical             :: even_dims, test(2)
        real(sp)            :: ang
        integer(kind=c_int) :: wsdm_ret
        ! kill possibly pre-existing object
        call self%kill
        ! set particle index range
        self%pfromto = pfromto
        ! error check
        if( self%pfromto(2) - self%pfromto(1) + 1 < 1 )then
            write(*,*) 'pfromto: ', self%pfromto(1), self%pfromto(2)
            call simple_stop ('nptcls (# of particles) must be > 0; new; simple_polarft_corrcalc')
        endif
        if( nrefs < 1 )then
            write(*,*) 'nrefs: ', nrefs
            call simple_stop ('nrefs (# of reference sections) must be > 0; new; simple_polarft_corrcalc')
        endif
        self%ldim = [params_glob%boxmatch,params_glob%boxmatch,1] !< logical dimensions of original cartesian image
        test      = .false.
        test(1)   = is_even(self%ldim(1))
        test(2)   = is_even(self%ldim(2))
        even_dims = all(test)
        if( .not. even_dims )then
            write(*,*) 'self%ldim: ', self%ldim
            call simple_stop ('only even logical dims supported; new; simple_polarft_corrcalc')
        endif
        ! set constants
        if( present(ptcl_mask) )then
            self%nptcls  = count(ptcl_mask)                      !< the total number of particles in partition
        else
            self%nptcls  = self%pfromto(2) - self%pfromto(1) + 1                   !< the total number of particles in partition
        endif
        self%nrefs       = nrefs                                 !< the number of references (logically indexded [1,nrefs])
        self%nrots       = round2even(twopi * real(params_glob%ring2))     !< number of in-plane rotations for one pft  (determined by radius of molecule)
        self%pftsz       = self%nrots / 2                        !< size of reference (nrots/2) (number of vectors used for matching)
        ! take care of objective function flags
        allocate(self%inv_resarrsq(params_glob%kfromto(1):params_glob%kstop), stat=alloc_stat)
        if(alloc_stat/=0)call allocchk("Inv_resarrsq failed allocation; simple_polarft_corrcalc ; new")
        self%inv_resarrsq = params_glob%smpd * real(self%ldim(1)) / real((/(k,k=params_glob%kfromto(1),params_glob%kstop)/))
        self%inv_resarrsq = -1. / (4.*self%inv_resarrsq*self%inv_resarrsq)
        select case(params_glob%cc_objfun)
            case (OBJFUN_CC, OBJFUN_EUCLID)
            case (OBJFUN_RES)
                allocate(self%ptcl_bfac_weights(params_glob%kfromto(1):params_glob%kstop, 1:self%nptcls),&
                    &self%ptcl_bfac_norms(1:self%nptcls))
                self%ptcl_bfac_weights = 1.0
                self%ptcl_bfac_norms   = real(params_glob%kstop - params_glob%kfromto(1) + 1) !< # resolution elements in the band-pass limited PFTs
            case DEFAULT
                write(*,*) 'unsupported objective function: ', trim(params_glob%objfun)
                stop 'ABORTING, simple_polarft_corrcalc :: new'
        end select
        ! generate polar coordinates & eo assignment
        allocate( self%polar(2*self%nrots,params_glob%kfromto(1):params_glob%kfromto(2)),&
                 &self%angtab(self%nrots), self%iseven(1:self%nptcls), stat=alloc_stat)
        if(alloc_stat/=0)call allocchk('polar coordinate arrays; new; simple_polarft_corrcalc, 1')
        ang = twopi/real(self%nrots)
        do irot=1,self%nrots
            self%angtab(irot) = real(irot-1)*ang
            do k=params_glob%kfromto(1),params_glob%kfromto(2)
                self%polar(irot,k)            = cos(self%angtab(irot))*real(k) ! x-coordinate
                self%polar(irot+self%nrots,k) = sin(self%angtab(irot))*real(k) ! y-coordinate
            end do
            self%angtab(irot) = rad2deg(self%angtab(irot)) ! angle (in degrees)
        end do
        ! index translation table
        allocate( self%pinds(self%pfromto(1):self%pfromto(2)), source=0, stat=alloc_stat)
        if(alloc_stat/=0)call allocchk('polar coordinate arrays; new; simple_polarft_corrcalc, 2')
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
        allocate( self%argtransf(self%nrots,params_glob%kfromto(1):params_glob%kfromto(2)), stat=alloc_stat)
        if(alloc_stat.ne.0)call allocchk('shift argument transfer array; new; simple_polarft_corrcalc')
        self%argtransf(:self%pftsz,:)   = &
            self%polar(:self%pftsz,:)   * &
            (PI/real(self%ldim(1) / 2))    ! x-part
        self%argtransf(self%pftsz + 1:,:) = &
            self%polar(self%nrots + 1:self%nrots+self%pftsz,:) * &
            (PI/real(self%ldim(2) / 2))    ! y-part
        ! allocate others
        allocate(self%pfts_refs_even(self%pftsz,params_glob%kfromto(1):params_glob%kfromto(2),self%nrefs),&
                 &self%pfts_refs_odd(self%pftsz,params_glob%kfromto(1):params_glob%kfromto(2),self%nrefs),&
                 &self%pfts_drefs_even(self%pftsz,params_glob%kfromto(1):params_glob%kstop,3,params_glob%nthr),&
                 &self%pfts_drefs_odd (self%pftsz,params_glob%kfromto(1):params_glob%kstop,3,params_glob%nthr),&
                 &self%pfts_ptcls(self%pftsz,params_glob%kfromto(1):params_glob%kfromto(2),1:self%nptcls),&
                 &self%sqsums_ptcls(1:self%nptcls),self%fftdat(params_glob%nthr),self%fft_carray(params_glob%nthr),&
                 &self%fftdat_ptcls(1:self%nptcls,params_glob%kfromto(1):params_glob%kstop),&
                 &self%heap_vars(params_glob%nthr),stat=alloc_stat)
        if(alloc_stat.ne.0)call allocchk('shared arrays; new; simple_polarft_corrcalc')
        local_stat=0
        do ithr=1,params_glob%nthr
            allocate(self%heap_vars(ithr)%pft_ref(self%pftsz,params_glob%kfromto(1):params_glob%kfromto(2)),&
                &self%heap_vars(ithr)%pft_ref_tmp(self%pftsz,params_glob%kfromto(1):params_glob%kfromto(2)),&
                &self%heap_vars(ithr)%pft_ref_tmp1(self%pftsz,params_glob%kfromto(1):params_glob%kfromto(2)),&
                &self%heap_vars(ithr)%pft_ref_tmp2(self%pftsz,params_glob%kfromto(1):params_glob%kfromto(2)),&
                &self%heap_vars(ithr)%pft_dref(self%pftsz,params_glob%kfromto(1):params_glob%kfromto(2),3),&
                &self%heap_vars(ithr)%corrs_over_k(self%nrots),&
                &self%heap_vars(ithr)%argmat(self%pftsz,params_glob%kfromto(1):params_glob%kfromto(2)),&
                &self%heap_vars(ithr)%shmat(self%pftsz,params_glob%kfromto(1):params_glob%kfromto(2)),&
                &self%heap_vars(ithr)%kcorrs(self%nrots),&
                &self%heap_vars(ithr)%pft_ref_8(self%pftsz,params_glob%kfromto(1):params_glob%kfromto(2)),&
                &self%heap_vars(ithr)%pft_ref_tmp_8(self%pftsz,params_glob%kfromto(1):params_glob%kfromto(2)),&
                &self%heap_vars(ithr)%pft_ref_tmp1_8(self%pftsz,params_glob%kfromto(1):params_glob%kfromto(2)),&
                &self%heap_vars(ithr)%pft_ref_tmp2_8(self%pftsz,params_glob%kfromto(1):params_glob%kfromto(2)),&
                &self%heap_vars(ithr)%pft_dref_8(self%pftsz,params_glob%kfromto(1):params_glob%kfromto(2),3),&
                &self%heap_vars(ithr)%shmat_8(self%pftsz,params_glob%kfromto(1):params_glob%kfromto(2)),&
                &self%heap_vars(ithr)%argmat_8(self%pftsz,params_glob%kfromto(1):params_glob%kfromto(2)),&
                &self%heap_vars(ithr)%fdf_y_8(params_glob%kfromto(1):params_glob%kfromto(2)),&
                &self%heap_vars(ithr)%fdf_T1_8(params_glob%kfromto(1):params_glob%kfromto(2),3),&
                &self%heap_vars(ithr)%fdf_T2_8(params_glob%kfromto(1):params_glob%kfromto(2),3),&
                &stat=alloc_stat)
            if(alloc_stat.ne.0) exit
        end do
        if(alloc_stat.ne.0)call allocchk('polarfts and sqsums; new; simple_polarft_corrcalc')
        self%pfts_refs_even = zero
        self%pfts_refs_odd  = zero
        self%pfts_ptcls     = zero
        self%sqsums_ptcls   = 0.
        ! set CTF flag
        self%with_ctf = .false.
        if( params_glob%ctf .ne. 'no' ) self%with_ctf = .true.
        ! thread-safe c-style allocatables for gencorrs
        do ithr=1,params_glob%nthr
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
            ! thread-safe c-style allocatables for on-the-fly particle memoization
            self%fft_carray(ithr)%p_re = fftwf_alloc_real(int(self%pftsz, c_size_t))
            self%fft_carray(ithr)%p_im = fftwf_alloc_complex(int(self%pftsz, c_size_t))
            call c_f_pointer(self%fft_carray(ithr)%p_re, self%fft_carray(ithr)%re, [self%pftsz])
            call c_f_pointer(self%fft_carray(ithr)%p_im, self%fft_carray(ithr)%im, [self%pftsz])
        end do
        ! thread-safe c-style allocatables for gencorrs, particle memoization
        do i = 1,self%nptcls
            do ik = params_glob%kfromto(1),params_glob%kstop
                self%fftdat_ptcls(i,ik)%p_re = fftwf_alloc_complex(int(self%pftsz, c_size_t))
                self%fftdat_ptcls(i,ik)%p_im = fftwf_alloc_complex(int(self%pftsz, c_size_t))
                call c_f_pointer(self%fftdat_ptcls(i,ik)%p_re, self%fftdat_ptcls(i,ik)%re, [self%pftsz])
                call c_f_pointer(self%fftdat_ptcls(i,ik)%p_im, self%fftdat_ptcls(i,ik)%im, [self%pftsz])
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
        ! flag existence
        self%existence = .true.
        ! set pointer to global instance
        pftcc_glob => self
    end subroutine new

    ! SETTERS

    !>  \brief  sets reference pft iref
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
        ! calculate the square sum required for correlation calculation
        call self%memoize_sqsum_ptcl(self%pinds(iptcl))
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

    subroutine zero_ref( self, iref )
        class(polarft_corrcalc), intent(inout) :: self
        integer,                 intent(in)    :: iref
        self%pfts_refs_even(:,:,iref) = zero
        self%pfts_refs_odd(:,:,iref)  = zero
    end subroutine zero_ref

    subroutine cp_even2odd_ref( self, iref )
        class(polarft_corrcalc), intent(inout) :: self
        integer,                 intent(in)    :: iref
        self%pfts_refs_odd(:,:,iref) = self%pfts_refs_even(:,:,iref)
    end subroutine cp_even2odd_ref

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

    subroutine assign_sigma2_noise( self, sigma2_noise, sigma2_exists_msk )
        class(polarft_corrcalc),      intent(inout) :: self
        real,    allocatable, target, intent(inout) :: sigma2_noise(:,:)
        logical, allocatable, target, intent(inout) :: sigma2_exists_msk(:)
        self%sigma2_noise      => sigma2_noise
        self%sigma2_exists_msk => sigma2_exists_msk
    end subroutine assign_sigma2_noise

    ! GETTERS

    !>  \brief  for getting the number of in-plane rotations
    pure function get_nrots( self ) result( nrots )
        class(polarft_corrcalc), intent(in) :: self
        integer :: nrots
        nrots = self%nrots
    end function get_nrots

    ! !>  \brief  for getting the dimensions of the reference polar FT
    function get_pdim( self ) result( pdim )
        class(polarft_corrcalc), intent(in) :: self
        integer :: pdim(3)
        pdim = [self%pftsz,params_glob%kfromto(1),params_glob%kfromto(2)]
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
        allocate(pft(self%pftsz,params_glob%kfromto(1):params_glob%kfromto(2)),&
        source=self%pfts_ptcls(:,:,self%pinds(iptcl)), stat=alloc_stat)
        if(alloc_stat.ne.0)call allocchk("In: get_ptcl_pft; simple_polarft_corrcalc",alloc_stat)
    end function get_ptcl_pft

    !>  \brief  returns polar Fourier transform of reference iref
    !! \param iref reference index
    function get_ref_pft( self, iref, iseven ) result( pft )
        class(polarft_corrcalc), intent(in) :: self
        integer,                 intent(in) :: iref
        logical,                 intent(in) :: iseven
        complex(sp), allocatable :: pft(:,:)
        if( iseven )then
            allocate(pft(self%pftsz,params_glob%kfromto(1):params_glob%kfromto(2)),&
            source=self%pfts_refs_even(:,:,iref), stat=alloc_stat)
        else
            allocate(pft(self%pftsz,params_glob%kfromto(1):params_glob%kfromto(2)),&
            source=self%pfts_refs_odd(:,:,iref), stat=alloc_stat)
        endif
        if(alloc_stat.ne.0)call allocchk("In: get_ref_pft; simple_polarft_corrcalc")
    end function get_ref_pft

    !>  \brief  returns # references
    integer function get_nrefs( self )
        class(polarft_corrcalc), intent(in) :: self
        get_nrefs = self%nrefs
    end function get_nrefs

    !>  \brief  checks for existence
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

    function is_euclid( self, iptcl ) result( res )
        class(polarft_corrcalc), intent(in) :: self
        integer,                 intent(in) :: iptcl
        logical                             :: res
        if (params_glob%cc_objfun .ne. OBJFUN_EUCLID) then
            res = .false.
        else
            if (self%pinds(iptcl) .eq. 0) then
                write (*,*) 'error in simple_polarft_corrcalc: sigma2_exists. iptcl = ', iptcl
                stop 'error in simple_polarft_corrcalc: is_euclid. iptcl wrong! '
            else
                res = self%sigma2_exists_msk(self%pinds(iptcl))
            end if
        end if
    end function is_euclid

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

    !>  \brief  for printing info about the object
    subroutine print( self )
        class(polarft_corrcalc), intent(in) :: self
        write(*,*) "total n particles in partition         (self%nptcls): ", self%nptcls
        write(*,*) "number of references                    (self%nrefs): ", self%nrefs
        write(*,*) "number of rotations                     (self%nrots): ", self%nrots
        write(*,*) "size of pft                             (self%pftsz): ", self%pftsz
        write(*,*) "logical dim. of original Cartesian image (self%ldim): ", self%ldim
    end subroutine print

    ! MEMOIZERS

    subroutine memoize_sqsum_ptcl( self, i )
        class(polarft_corrcalc), intent(inout) :: self
        integer,                 intent(in)    :: i
        self%sqsums_ptcls(i) = sum(csq(self%pfts_ptcls(:,:,i)))
    end subroutine memoize_sqsum_ptcl

    ! ! memoize particles ffts
    ! subroutine memoize_ffts( self )
    !     class(polarft_corrcalc), intent(inout) :: self
    !     type(fftw_carr) :: carray(params_glob%nthr)
    !     integer         :: i, ik, ithr
    !     ! memoize particle FFTs in parallel
    !     !$omp parallel do default(shared) private(i,ik,ithr) proc_bind(close) collapse(2) schedule(static)
    !     do i = 1, self%nptcls
    !         do ik = params_glob%kfromto(1),params_glob%kstop
    !             ! get thread index
    !             ithr = omp_get_thread_num() + 1
    !             ! copy particle pfts
    !             carray(ithr)%re = real(self%pfts_ptcls(:,ik,i))
    !             carray(ithr)%im = aimag(self%pfts_ptcls(:,ik,i)) * self%fft_factors
    !             ! FFT
    !             call fftwf_execute_dft_r2c(self%plan_fwd_1, carray(ithr)%re, self%fftdat_ptcls(i,ik)%re)
    !             call fftwf_execute_dft    (self%plan_fwd_2, carray(ithr)%im, self%fftdat_ptcls(i,ik)%im)
    !         end do
    !     end do
    !     !$omp end parallel do
    ! end subroutine memoize_ffts

    ! memoize all particles ffts
    subroutine memoize_ffts( self )
        class(polarft_corrcalc), intent(inout) :: self
        integer :: i
        ! memoize particle FFTs in parallel
        !$omp parallel do default(shared) private(i) proc_bind(close) schedule(static)
        do i = 1, self%nptcls
            call self%memoize_fft(i)
        end do
        !$omp end parallel do
    end subroutine memoize_ffts

    ! memoize single particle fft
    subroutine memoize_ptcl_fft( self, iptcl )
        class(polarft_corrcalc), intent(inout) :: self
        integer,                 intent(in)    :: iptcl
        integer :: i
        i = self%pinds(iptcl)
        call self%memoize_fft(i)
    end subroutine memoize_ptcl_fft

    ! memoize particle fft, serial only
    subroutine memoize_fft( self, i )
        class(polarft_corrcalc), intent(inout) :: self
        integer  :: i, ik, ithr
        ithr = omp_get_thread_num() + 1
        ! memoize particle FFTs
        do ik = params_glob%kfromto(1),params_glob%kstop
            ! copy particle pfts
            self%fft_carray(ithr)%re = real(self%pfts_ptcls(:,ik,i))
            self%fft_carray(ithr)%im = aimag(self%pfts_ptcls(:,ik,i)) * self%fft_factors
            ! FFT
            call fftwf_execute_dft_r2c(self%plan_fwd_1, self%fft_carray(ithr)%re, self%fftdat_ptcls(i,ik)%re)
            call fftwf_execute_dft    (self%plan_fwd_2, self%fft_carray(ithr)%im, self%fftdat_ptcls(i,ik)%im)
        end do
    end subroutine memoize_fft

    subroutine memoize_bfac( self, iptcl, bfac )
        class(polarft_corrcalc), intent(inout) :: self
        integer,                 intent(in)    :: iptcl
        real,                    intent(in)    :: bfac
        if( params_glob%cc_objfun == OBJFUN_CC )then
            ! nothing to do
        else if ( params_glob%cc_objfun == OBJFUN_RES ) then
            self%ptcl_bfac_weights(:,self%pinds(iptcl)) = exp( bfac * self%inv_resarrsq(:) ) ! exp( -bfac/(4.*res^2) )
            where( self%ptcl_bfac_weights(:,self%pinds(iptcl)) < TINY) self%ptcl_bfac_weights(:,self%pinds(iptcl)) = 0.
            self%ptcl_bfac_norms(self%pinds(iptcl)) = sum(self%ptcl_bfac_weights(:,self%pinds(iptcl)))
        else
            ! nothing to do
        end if
    end subroutine memoize_bfac

    ! CALCULATORS

    subroutine create_polar_absctfmats( self, spproj, oritype, pfromto )
        !$ use omp_lib
        !$ use omp_lib_kinds
        use simple_ctf,        only: ctf
        use simple_sp_project, only: sp_project
        class(polarft_corrcalc),   intent(inout) :: self
        class(sp_project), target, intent(inout) :: spproj
        character(len=*),          intent(in)    :: oritype
        integer, optional,         intent(in)    :: pfromto(2)
        type(ctfparams) :: ctfparms(nthr_glob)
        type(ctf)       :: tfuns(nthr_glob)
        real(sp)        :: inv_ldim(3),hinv,kinv,spaFreqSq,ang
        integer         :: i,irot,k,iptcl,ithr,ppfromto(2),ctfmatind
        logical         :: present_pfromto
        present_pfromto = present(pfromto)
        ppfromto = self%pfromto
        if( present_pfromto ) ppfromto = pfromto
        if( allocated(self%ctfmats) ) deallocate(self%ctfmats)
        allocate(self%ctfmats(self%pftsz,params_glob%kfromto(1):params_glob%kfromto(2),1:self%nptcls), stat=alloc_stat)
        if(alloc_stat.ne.0)call allocchk("In: simple_polarft_corrcalc :: create_polar_ctfmats, 2",alloc_stat)
        inv_ldim = 1./real(self%ldim)
        !$omp parallel do default(shared) private(i,iptcl,ctfmatind,ithr,irot,k,hinv,kinv,spaFreqSq,ang)&
        !$omp schedule(static) proc_bind(close)
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
                do irot=1,self%pftsz
                    do k=params_glob%kfromto(1),params_glob%kfromto(2)
                        hinv           = self%polar(irot,k)*inv_ldim(1)
                        kinv           = self%polar(irot+self%nrots,k)*inv_ldim(2)
                        spaFreqSq      = hinv*hinv+kinv*kinv
                        ang            = atan2(self%polar(irot+self%nrots,k),self%polar(irot,k))
                        if( ctfparms(ithr)%l_phaseplate )then
                            self%ctfmats(irot,k,self%pinds(ctfmatind)) = abs( tfuns(ithr)%eval(spaFreqSq, ang, ctfparms(ithr)%phshift) )
                        else
                            self%ctfmats(irot,k,self%pinds(ctfmatind)) = abs( tfuns(ithr)%eval(spaFreqSq, ang) )
                        endif
                    end do
                end do
            endif
        end do
        !$omp end parallel do
    end subroutine create_polar_absctfmats

    subroutine prep_ref4corr( self, iref, i, pft_ref, sqsum_ref, kstop )
        class(polarft_corrcalc), intent(inout) :: self
        integer,                 intent(in)    :: iref, i
        integer,                 intent(in)    :: kstop
        complex(sp),             intent(out)   :: pft_ref(self%pftsz,params_glob%kfromto(1):params_glob%kfromto(2))
        real(sp),                intent(out)   :: sqsum_ref
        ! copy
        if( self%iseven(i) )then
            pft_ref = self%pfts_refs_even(:,:,iref)
        else
            pft_ref = self%pfts_refs_odd(:,:,iref)
        endif
        ! multiply with CTF
        if( self%with_ctf ) pft_ref = pft_ref * self%ctfmats(:,:,i)
        ! for corr normalisation
        sqsum_ref = sum(csq(pft_ref(:,params_glob%kfromto(1):kstop)))
    end subroutine prep_ref4corr

    subroutine calc_corrs_over_k( self, pft_ref, i, kfromto, corrs_over_k )
        class(polarft_corrcalc), intent(inout) :: self
        integer,                 intent(in)    :: i, kfromto(2)
        complex(sp),             intent(in)    :: pft_ref(1:self%pftsz,params_glob%kfromto(1):params_glob%kfromto(2))
        real,                    intent(out)   :: corrs_over_k(self%nrots)
        integer :: ithr, ik
        ! get thread index
        ithr = omp_get_thread_num() + 1
        ! sum up correlations over k-rings
        corrs_over_k = 0.
        do ik = kfromto(1),kfromto(2)
            ! move reference into Fourier Fourier space (particles are memoized)
            self%fftdat(ithr)%ref_re(:) =  real(pft_ref(:,ik))
            self%fftdat(ithr)%ref_im(:) = aimag(pft_ref(:,ik)) * self%fft_factors
            call fftwf_execute_dft_r2c(self%plan_fwd_1, self%fftdat(ithr)%ref_re, self%fftdat(ithr)%ref_fft_re)
            call fftwf_execute_dft    (self%plan_fwd_2, self%fftdat(ithr)%ref_im, self%fftdat(ithr)%ref_fft_im)
            ! correlate
            self%fftdat(ithr)%ref_fft_re = conjg(self%fftdat(ithr)%ref_fft_re) * self%fftdat_ptcls(i,ik)%re
            self%fftdat(ithr)%ref_fft_im = conjg(self%fftdat(ithr)%ref_fft_im) * self%fftdat_ptcls(i,ik)%im
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
    end subroutine calc_corrs_over_k

    subroutine calc_k_corrs( self, pft_ref, i, k, kcorrs )
        class(polarft_corrcalc), intent(inout) :: self
        complex(sp),             intent(in)    :: pft_ref(1:self%pftsz,params_glob%kfromto(1):params_glob%kfromto(2))
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
        self%fftdat(ithr)%ref_fft_re = conjg(self%fftdat(ithr)%ref_fft_re) * self%fftdat_ptcls(i,k)%re
        self%fftdat(ithr)%ref_fft_im = conjg(self%fftdat(ithr)%ref_fft_im) * self%fftdat_ptcls(i,k)%im
        self%fftdat(ithr)%product_fft(1:1 + 2 * int(self%pftsz / 2):2) = &
            4. * self%fftdat(ithr)%ref_fft_re(1:1 + int(self%pftsz / 2))
        self%fftdat(ithr)%product_fft(2:2 + 2 * int(self%pftsz / 2):2) = &
            4. * self%fftdat(ithr)%ref_fft_im(1:int(self%pftsz / 2) + 1)
        ! back transform
        call fftwf_execute_dft_c2r(self%plan_bwd, self%fftdat(ithr)%product_fft, self%fftdat(ithr)%backtransf)
        ! fftw3 routines are not properly normalized, hence division by self%nrots * 2
        kcorrs = self%fftdat(ithr)%backtransf / real(self%nrots * 2)
    end subroutine calc_k_corrs

    function calc_corr_for_rot(self, pft_ref, i, irot) result(corr)
        class(polarft_corrcalc), intent(inout) :: self
        integer,                 intent(in)    :: i, irot
        complex(sp),             intent(in)    :: pft_ref(1:self%pftsz,params_glob%kfromto(1):params_glob%kfromto(2))
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
            tmp = sum( pft_ref(:,:) * conjg(self%pfts_ptcls(:,:,i)))
        else if( irot <= self%pftsz )then
            tmp =       sum( pft_ref(               1:self%pftsz-rot+1,:) * conjg(self%pfts_ptcls(rot:self%pftsz,:,i)))
            tmp = tmp + sum( pft_ref(self%pftsz-rot+2:self%pftsz,      :) *       self%pfts_ptcls(  1:rot-1,     :,i))
        else if( irot == self%pftsz + 1 )then
            tmp = sum( pft_ref(:,:) * self%pfts_ptcls(:,:,i) )
        else
            tmp =       sum( pft_ref(1:self%pftsz-rot+1,:)          *        self%pfts_ptcls(rot:self%pftsz,:,i))
            tmp = tmp + sum( pft_ref(self%pftsz-rot+2:self%pftsz,:) * conjg( self%pfts_ptcls( 1:rot-1,      :,i)))
        end if
        corr = real(tmp)
    end function calc_corr_for_rot

    function calc_corr_for_rot_8( self, pft_ref, i, irot ) result( corr )
        class(polarft_corrcalc), intent(inout) :: self
        integer,                 intent(in)    :: i, irot
        complex(dp),             intent(in)    :: pft_ref(1:self%pftsz,params_glob%kfromto(1):params_glob%kfromto(2))
        integer     :: rot
        real(dp)    :: corr
        complex(sp) :: tmp
        tmp = 0.
        if( irot >= self%pftsz + 1 )then
            rot = irot - self%pftsz
        else
            rot = irot
        end if
        if( irot == 1 )then
            tmp = sum( pft_ref(:,:) * conjg(self%pfts_ptcls(:,:,i)))
        else if( irot <= self%pftsz )then
            tmp =       sum( pft_ref(               1:self%pftsz-rot+1,:) * conjg(self%pfts_ptcls(rot:self%pftsz,:,i)))
            tmp = tmp + sum( pft_ref(self%pftsz-rot+2:self%pftsz,      :) *       self%pfts_ptcls(  1:rot-1,     :,i))
        else if( irot == self%pftsz + 1 ) then
            tmp = sum( pft_ref(:,:) * self%pfts_ptcls(:,:,i) )
        else
            tmp =       sum( pft_ref(1:self%pftsz-rot+1,:)          *        self%pfts_ptcls(rot:self%pftsz,:,i))
            tmp = tmp + sum( pft_ref(self%pftsz-rot+2:self%pftsz,:) * conjg( self%pfts_ptcls( 1:rot-1,      :,i)))
        end if
        corr = real(tmp, kind=dp)
    end function calc_corr_for_rot_8

    !<  \brief  compute the terms T1, T2 necessary for finding the derivative of the correlations, double precision
    subroutine calc_T1_T2_for_rot_8( self, pft_ref, pft_dref, i, irot, nderivs, T1_k, T2_k)
        class(polarft_corrcalc), intent(inout) :: self
        integer,                 intent(in)    :: i, irot, nderivs
        complex(dp),             intent(in)    :: pft_ref(1:self%pftsz,params_glob%kfromto(1):params_glob%kfromto(2)), &
            pft_dref(1:self%pftsz,params_glob%kfromto(1):params_glob%kstop,nderivs)
        real(dp),                intent(out)   :: T1_k(params_glob%kfromto(1):params_glob%kstop,nderivs), T2_k(params_glob%kfromto(1):params_glob%kstop,nderivs)
        integer                                :: k, rot, j
        complex(dp)                            :: tmp_T1(params_glob%kfromto(1):params_glob%kstop,nderivs)
        complex(dp)                            :: tmp_T2(params_glob%kfromto(1):params_glob%kstop,nderivs)
        if( irot >= self%pftsz + 1 )then
            rot = irot - self%pftsz
        else
            rot = irot
        end if
        do j = 1, nderivs
            do k = params_glob%kfromto(1),params_glob%kstop
                if( irot == 1 ) then
                    tmp_T1(k,j) = sum( pft_dref(:,k,j) * conjg(self%pfts_ptcls(:,k,i)))
                else if( irot <= self%pftsz ) then
                    tmp_T1(k,j) =               sum( pft_dref(               1:self%pftsz-rot+1,k,j) * conjg(self%pfts_ptcls(rot:self%pftsz,k,i)))
                    tmp_T1(k,j) = tmp_T1(k,j) + sum( pft_dref(self%pftsz-rot+2:self%pftsz,      k,j) *       self%pfts_ptcls(  1:rot-1,     k,i))
                else if( irot == self%pftsz + 1 ) then
                    tmp_T1(k,j) = sum( pft_dref(:,k,j) * self%pfts_ptcls(:,k,i) )
                else
                    tmp_T1(k,j) =               sum( pft_dref(               1:self%pftsz-rot+1,k,j) *       self%pfts_ptcls(rot:self%pftsz,k,i))
                    tmp_T1(k,j) = tmp_T1(k,j) + sum( pft_dref(self%pftsz-rot+2:self%pftsz,      k,j) * conjg(self%pfts_ptcls(  1:rot-1,     k,i)) )
                end if
                tmp_T2(k,j) = sum( pft_dref(:,k,j) * conjg(pft_ref(:,k)))
            end do
        end do
        T1_k = real(tmp_T1, kind=dp)
        T2_k = real(tmp_T2, kind=dp)
    end subroutine calc_T1_T2_for_rot_8

    function calc_euclid_for_rot( self, pft_ref, i, irot ) result( euclid )
        class(polarft_corrcalc), intent(inout) :: self
        integer,                 intent(in)    :: i, irot
        complex(sp),             intent(in)    :: pft_ref(1:self%pftsz,params_glob%kfromto(1):params_glob%kfromto(2))
        integer     :: rot, k
        real(sp)    :: euclid, euclidtmp
        complex(sp) :: tmp
        if( irot >= self%pftsz + 1 )then
            rot = irot - self%pftsz
        else
            rot = irot
        end if
        euclid = 1.
        do k = params_glob%kfromto(1), params_glob%kstop
            if( irot == 1 )then
                tmp =       sum(csq(pft_ref(:,k) - self%pfts_ptcls(:,k,i)))
            else if( irot <= self%pftsz )then
                tmp =       sum(csq(pft_ref(1:self%pftsz-rot+1,k) - self%pfts_ptcls(rot:self%pftsz,k,i)))
                tmp = tmp + sum(csq(pft_ref(self%pftsz-rot+2:self%pftsz,k) - conjg(self%pfts_ptcls(1:rot-1,k,i))))
            else if( irot == self%pftsz + 1 )then
                tmp = sum(csq(pft_ref(:,k) - conjg(self%pfts_ptcls(:,k,i))))
            else
                tmp =       sum(csq(pft_ref(1:self%pftsz-rot+1,k) - conjg(self%pfts_ptcls(rot:self%pftsz,k,i))))
                tmp = tmp + sum(csq(pft_ref(self%pftsz-rot+2:self%pftsz,k) - self%pfts_ptcls(1:rot-1,k,i)))
            end if
            euclid = euclid - tmp / ( 2. * self%sigma2_noise(k, i))
        end do
    end function calc_euclid_for_rot

    function calc_euclid_for_rot_8( self, pft_ref, i, irot ) result( euclid )
        class(polarft_corrcalc), intent(inout) :: self
        integer,                 intent(in)    :: i, irot
        complex(dp),             intent(in)    :: pft_ref(1:self%pftsz,params_glob%kfromto(1):params_glob%kfromto(2))
        integer     :: rot, k
        real(dp)    :: euclid
        complex(dp) :: tmp
        if( irot >= self%pftsz + 1 )then
            rot = irot - self%pftsz
        else
            rot = irot
        end if
        euclid = 1._dp
        do k = params_glob%kfromto(1), params_glob%kstop
        if( irot == 1 )then
                tmp =       sum(csq(pft_ref(:,k) - self%pfts_ptcls(:,k,i)))
        else if( irot <= self%pftsz )then
                tmp =       sum(csq(pft_ref(1:self%pftsz-rot+1,k) - self%pfts_ptcls(rot:self%pftsz,k,i)))
                tmp = tmp + sum(csq(pft_ref(self%pftsz-rot+2:self%pftsz,k) - conjg(self%pfts_ptcls(1:rot-1,k,i))))
        else if( irot == self%pftsz + 1 )then
                tmp = sum(csq(pft_ref(:,k) - conjg(self%pfts_ptcls(:,k,i))))
        else
                tmp =       sum(csq(pft_ref(1:self%pftsz-rot+1,k) - conjg(self%pfts_ptcls(rot:self%pftsz,k,i))))
                tmp = tmp + sum(csq(pft_ref(self%pftsz-rot+2:self%pftsz,k) - self%pfts_ptcls(1:rot-1,k,i)))
        end if
            euclid = euclid - tmp / ( 2._dp * self%sigma2_noise(k, i))
        end do
    end function calc_euclid_for_rot_8

    function calc_corrk_for_rot_8( self, pft_ref, i, k, irot ) result( corr )
        class(polarft_corrcalc), intent(inout) :: self
        integer,                 intent(in)    :: i, k, irot
        complex(dp),             intent(in)    :: pft_ref(1:self%pftsz,params_glob%kfromto(1):params_glob%kfromto(2))
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
            tmp = sum( pft_ref(:,k) * conjg(self%pfts_ptcls(:,k,i)))
        else if (irot <= self%pftsz) then
            tmp =       sum( pft_ref(               1:self%pftsz-rot+1,k) * conjg(self%pfts_ptcls(rot:self%pftsz,k,i)))
            tmp = tmp + sum( pft_ref(self%pftsz-rot+2:self%pftsz,      k) *       self%pfts_ptcls(  1:rot-1,     k,i))
        else if (irot == self%pftsz + 1) then
            tmp = sum( pft_ref(:,k) * self%pfts_ptcls(:,k,i) )
        else
            tmp =       sum( pft_ref(1:self%pftsz-rot+1,k)          *        self%pfts_ptcls(rot:self%pftsz,k,i))
            tmp = tmp + sum( pft_ref(self%pftsz-rot+2:self%pftsz,k) * conjg( self%pfts_ptcls(1:rot-1,k,i) ))
        end if
        corr = real(tmp)
    end function calc_corrk_for_rot_8

    function calc_euclidk_for_rot( self, pft_ref, i, k, irot ) result( euclid )
        class(polarft_corrcalc), intent(inout) :: self
        integer,                 intent(in)    :: i, irot, k
        complex(sp),             intent(in)    :: pft_ref(1:self%pftsz,params_glob%kfromto(1):params_glob%kfromto(2))
        integer     :: rot
        real(sp)    :: euclid
        complex(sp) :: tmp
        euclid = 0.
        tmp = 0.
        if( irot >= self%pftsz + 1 )then
            rot = irot - self%pftsz
        else
            rot = irot
        end if
        if( irot == 1 )then
            tmp = sum(csq(pft_ref(:,k) - self%pfts_ptcls(:,k,i)))
        else if( irot <= self%pftsz )then
            tmp =       sum(csq(pft_ref(               1:self%pftsz-rot+1,k) -       self%pfts_ptcls(rot:self%pftsz,k,i)))
            tmp = tmp + sum(csq(pft_ref(self%pftsz-rot+2:self%pftsz,      k) - conjg(self%pfts_ptcls(  1:rot-1,     k,i))))
        else if( irot == self%pftsz + 1 )then
            tmp = sum(csq(pft_ref(:,k) - conjg(self%pfts_ptcls(:,k,i))))
        else
            tmp =       sum(csq(pft_ref(               1:self%pftsz-rot+1,k) - conjg(self%pfts_ptcls(rot:self%pftsz,k,i))))
            tmp = tmp + sum(csq(pft_ref(self%pftsz-rot+2:self%pftsz,      k) -       self%pfts_ptcls(  1:rot-1,     k,i)))
        end if
        euclid = real(tmp)
    end function calc_euclidk_for_rot

    subroutine genfrc( self, iref, i, irot, frc )
        class(polarft_corrcalc),  intent(inout) :: self
        integer,                  intent(in)    :: iref, i, irot
        real(sp),                 intent(out)   :: frc(params_glob%kfromto(1):params_glob%kstop)
        complex(sp), pointer :: pft_ref(:,:)
        real(sp),    pointer :: kcorrs(:)
        real(sp) :: sumsqref, sumsqptcl, sqsum_ref
        integer  :: k, ithr
        ithr     =  omp_get_thread_num() + 1
        pft_ref  => self%heap_vars(ithr)%pft_ref
        kcorrs   => self%heap_vars(ithr)%kcorrs
        call self%prep_ref4corr(iref, i, pft_ref, sqsum_ref, params_glob%kstop)
        do k=params_glob%kfromto(1),params_glob%kstop
            call self%calc_k_corrs(pft_ref, i, k, kcorrs)
            sumsqptcl = sum(csq(self%pfts_ptcls(:,k,i)))
            sumsqref  = sum(csq(pft_ref(:,k)))
            frc(k)    = kcorrs(irot) / sqrt(sumsqref * sumsqptcl)
        end do
    end subroutine genfrc

    subroutine calc_frc( self, iref, iptcl, irot, shvec, frc )
        class(polarft_corrcalc),  intent(inout) :: self
        integer,                  intent(in)    :: iref, iptcl, irot
        real(sp),                 intent(in)    :: shvec(2)
        real(sp),                 intent(out)   :: frc(params_glob%kfromto(1):params_glob%kstop)
        complex(sp), pointer :: pft_ref(:,:), shmat(:,:)
        real(sp),    pointer :: kcorrs(:), argmat(:,:)
        real(sp) :: sumsqref, sumsqptcl
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
                pft_ref = (self%pfts_refs_even(:,:,iref) * self%ctfmats(:,:,i)) * shmat
            else
                pft_ref = (self%pfts_refs_odd(:,:,iref)  * self%ctfmats(:,:,i)) * shmat
            endif
        else
            if( self%iseven(i) )then
                pft_ref = self%pfts_refs_even(:,:,iref) * shmat
            else
                pft_ref = self%pfts_refs_odd(:,:,iref)  * shmat
            endif
        endif
        do k=params_glob%kfromto(1),params_glob%kstop
            call self%calc_k_corrs(pft_ref, i, k, kcorrs)
            sumsqptcl = sum(csq(self%pfts_ptcls(:,k,i)))
            sumsqref  = sum(csq(pft_ref(:,k)))
            frc(k)    = kcorrs(irot) / sqrt(sumsqref * sumsqptcl)
        end do
    end subroutine calc_frc

    subroutine gencorrs_cc_1( self, iref, iptcl, cc )
        class(polarft_corrcalc), intent(inout) :: self
        integer,                 intent(in)    :: iref, iptcl
        real,                    intent(out)   :: cc(self%nrots)
        complex(sp), pointer :: pft_ref(:,:)
        real(sp),    pointer :: corrs_over_k(:)
        real(sp) :: sqsum_ref, sqsum_ptcl
        integer  :: ithr
        ithr         =  omp_get_thread_num() + 1
        pft_ref      => self%heap_vars(ithr)%pft_ref
        corrs_over_k => self%heap_vars(ithr)%corrs_over_k
        call self%prep_ref4corr(iref, self%pinds(iptcl), pft_ref, sqsum_ref, params_glob%kstop)
        call self%calc_corrs_over_k(pft_ref, self%pinds(iptcl), [params_glob%kfromto(1),params_glob%kstop], corrs_over_k)
        sqsum_ptcl = sum(csq(self%pfts_ptcls(:, params_glob%kfromto(1):params_glob%kstop, self%pinds(iptcl))))
        cc = corrs_over_k / sqrt(sqsum_ref * sqsum_ptcl)
    end subroutine gencorrs_cc_1

    subroutine gencorrs_cc_2( self, iref, iptcl, shvec, cc )
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
                pft_ref = (self%pfts_refs_even(:,:,iref) * self%ctfmats(:,:,self%pinds(iptcl))) * shmat
            else
                pft_ref = (self%pfts_refs_odd (:,:,iref) * self%ctfmats(:,:,self%pinds(iptcl))) * shmat
            endif
        else
            if( self%iseven(self%pinds(iptcl)) )then
                pft_ref = self%pfts_refs_even(:,:,iref) * shmat
            else
                pft_ref = self%pfts_refs_odd (:,:,iref) * shmat
            endif
        endif
        sqsum_ref = sum(csq(pft_ref))
        call self%calc_corrs_over_k(pft_ref, self%pinds(iptcl), [params_glob%kfromto(1),params_glob%kstop], corrs_over_k)
        cc = corrs_over_k  / sqrt(sqsum_ref * self%sqsums_ptcls(self%pinds(iptcl)))
    end subroutine gencorrs_cc_2

    subroutine gencorrs_cc_valid( self, iref, iptcl, kfromto, cc )
        class(polarft_corrcalc), intent(inout) :: self
        integer,                 intent(in)    :: iref, iptcl, kfromto(2)
        real,                    intent(out)   :: cc(self%nrots)
        complex(sp), pointer :: pft_ref(:,:)
        real(sp),    pointer :: corrs_over_k(:)
        real(sp) :: sqsum_ref, sqsum_ptcl
        integer  :: ithr
        ithr         =  omp_get_thread_num() + 1
        pft_ref      => self%heap_vars(ithr)%pft_ref
        corrs_over_k => self%heap_vars(ithr)%corrs_over_k
        ! copy the opposite reference as compared to search
        if( .not. self%iseven(self%pinds(iptcl)) )then
            pft_ref = self%pfts_refs_even(:,:,iref)
        else
            pft_ref = self%pfts_refs_odd(:,:,iref)
        endif
        ! multiply with CTF
        if( self%with_ctf ) pft_ref = pft_ref * self%ctfmats(:,:,self%pinds(iptcl))
        ! for corr normalisation
        sqsum_ref = sum(csq(pft_ref(:,kfromto(1):kfromto(2))))
        call self%calc_corrs_over_k(pft_ref, self%pinds(iptcl), kfromto, corrs_over_k)
        sqsum_ptcl = sum(csq(self%pfts_ptcls(:,kfromto(1):kfromto(2),self%pinds(iptcl))))
        cc = corrs_over_k / sqrt(sqsum_ref * sqsum_ptcl)
    end subroutine gencorrs_cc_valid

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
        call self%prep_ref4corr(iref, self%pinds(iptcl), pft_ref, sumsqref, params_glob%kstop)
        cc(:) = 0.
        do k=params_glob%kfromto(1),params_glob%kstop
            call self%calc_k_corrs(pft_ref, self%pinds(iptcl), k, kcorrs)
            sumsqptcl = sum(csq(self%pfts_ptcls(:,k,self%pinds(iptcl))))
            sumsqref  = sum(csq(pft_ref(:,k)))
            cc(:) = cc(:) + (kcorrs(:) * self%ptcl_bfac_weights(k,self%pinds(iptcl))) / sqrt(sumsqref * sumsqptcl)
        end do
        cc(:) = cc(:) / sum(self%ptcl_bfac_weights(params_glob%kfromto(1):params_glob%kstop,self%pinds(iptcl)))
    end subroutine gencorrs_resnorm_1

    subroutine gencorrs_resnorm_2( self, iref, iptcl, shvec, cc )
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
                pft_ref = (self%pfts_refs_even(:,:,iref) * self%ctfmats(:,:,self%pinds(iptcl))) * shmat
            else
                pft_ref = (self%pfts_refs_odd (:,:,iref) * self%ctfmats(:,:,self%pinds(iptcl))) * shmat
            endif
        else
            if( self%iseven(self%pinds(iptcl)) )then
                pft_ref = self%pfts_refs_even(:,:,iref) * shmat
            else
                pft_ref = self%pfts_refs_odd (:,:,iref) * shmat
            endif
        endif
        cc(:) = 0.
        do k=params_glob%kfromto(1),params_glob%kstop
            call self%calc_k_corrs(pft_ref, self%pinds(iptcl), k, kcorrs)
            sumsqptcl = sum(csq(self%pfts_ptcls(:,k,self%pinds(iptcl))))
            sumsqref  = sum(csq(pft_ref(:,k)))
            cc(:) = cc(:) + (kcorrs(:) * self%ptcl_bfac_weights(k,self%pinds(iptcl))) / sqrt(sumsqref * sumsqptcl)
        end do
        cc(:) = cc(:) / self%ptcl_bfac_norms(self%pinds(iptcl))
    end subroutine gencorrs_resnorm_2

    subroutine gencorrs_euclid_1( self, iref, iptcl, euclids )
        class(polarft_corrcalc), intent(inout) :: self
        integer,                 intent(in)    :: iref, iptcl
        real(sp),                intent(out)   :: euclids(self%nrots)
        complex(sp), pointer :: pft_ref(:,:)
        real(sp),    pointer :: keuclids(:)
        real(sp) :: sumsqref, sumsqptcl
        integer  :: k, ithr
        ithr           =  omp_get_thread_num() + 1
        pft_ref        => self%heap_vars(ithr)%pft_ref
        keuclids       => self%heap_vars(ithr)%kcorrs !can be reused
        call self%prep_ref4corr(iref, self%pinds(iptcl), pft_ref, sumsqref, params_glob%kstop)
        euclids(:) = 0.
        do k=params_glob%kfromto(1),params_glob%kstop
            call self%calc_k_corrs(pft_ref, self%pinds(iptcl), k, keuclids)
            sumsqptcl = sum(csq(self%pfts_ptcls(:,k,self%pinds(iptcl))))
            sumsqref  = sum(csq(pft_ref(:,k)))
            euclids(:) = euclids(:) + (2. * keuclids(:) - sumsqptcl - sumsqref ) / (2. * self%sigma2_noise(k, self%pinds(iptcl)))
        end do
    end subroutine gencorrs_euclid_1

    subroutine gencorrs_euclid_2( self, iref, iptcl, shvec, euclids )
        class(polarft_corrcalc), intent(inout) :: self
        integer,                 intent(in)    :: iref, iptcl
        real(sp),                intent(in)    :: shvec(2)
        real(sp),                intent(out)   :: euclids(self%nrots)
        complex(sp), pointer :: pft_ref(:,:), shmat(:,:)
        real(sp),    pointer :: keuclids(:), argmat(:,:)
        real(sp) :: sumsqptcl, sumsqref
        integer  :: ithr, k
        ithr           = omp_get_thread_num() + 1
        pft_ref        => self%heap_vars(ithr)%pft_ref
        shmat          => self%heap_vars(ithr)%shmat
        keuclids       => self%heap_vars(ithr)%kcorrs !can be reused
        argmat         => self%heap_vars(ithr)%argmat
        argmat = self%argtransf(:self%pftsz,:) * shvec(1) + self%argtransf(self%pftsz + 1:,:) * shvec(2)
        shmat  = cmplx(cos(argmat),sin(argmat))
        if( self%with_ctf )then
            if( self%iseven(self%pinds(iptcl)) )then
                pft_ref = (self%pfts_refs_even(:,:,iref) * self%ctfmats(:,:,self%pinds(iptcl))) * shmat
            else
                pft_ref = (self%pfts_refs_odd (:,:,iref) * self%ctfmats(:,:,self%pinds(iptcl))) * shmat
            endif
        else
            if( self%iseven(self%pinds(iptcl)) )then
                pft_ref = self%pfts_refs_even(:,:,iref) * shmat
            else
                pft_ref = self%pfts_refs_odd (:,:,iref) * shmat
            endif
        endif
        euclids(:) = 0.
        do k=params_glob%kfromto(1),params_glob%kstop
            call self%calc_k_corrs(pft_ref, self%pinds(iptcl), k, keuclids)
            sumsqptcl = sum(csq(self%pfts_ptcls(:,k,self%pinds(iptcl))))
            sumsqref  = sum(csq(pft_ref(:,k)))
            euclids(:) = euclids(:) + (2. * keuclids(:) - sumsqptcl - sumsqref ) / (2. * self%sigma2_noise(k, self%pinds(iptcl)))
        end do
    end subroutine gencorrs_euclid_2

    subroutine gencorrs_1( self, iref, iptcl, cc )
        class(polarft_corrcalc), intent(inout) :: self
        integer,                 intent(in)    :: iref, iptcl
        real(sp),                intent(out)   :: cc(self%nrots)
        if     ( params_glob%cc_objfun == OBJFUN_CC )then
            call self%gencorrs_cc_1(iref, iptcl, cc)
        else if( params_glob%cc_objfun == OBJFUN_RES )then
            call self%gencorrs_resnorm_1(iref, iptcl, cc)
        else if( params_glob%cc_objfun == OBJFUN_EUCLID )then
            if (self%sigma2_exists_msk(self%pinds(iptcl))) then
                call self%gencorrs_euclid_1(iref, iptcl, cc)
            else
                call self%gencorrs_cc_1(iref, iptcl, cc)
            end if
        end if
    end subroutine gencorrs_1

    subroutine gencorrs_2( self, iref, iptcl, shvec, cc )
        class(polarft_corrcalc), intent(inout) :: self
        integer,                 intent(in)    :: iref, iptcl
        real(sp),                intent(in)    :: shvec(2)
        real(sp),                intent(out)   :: cc(self%nrots)
        if     ( params_glob%cc_objfun == OBJFUN_CC )then
            call self%gencorrs_cc_2(iref, iptcl, shvec, cc)
        else if( params_glob%cc_objfun == OBJFUN_RES )then
            call self%gencorrs_resnorm_2(iref, iptcl, shvec, cc)
        else if( params_glob%cc_objfun == OBJFUN_EUCLID )then
            if (self%sigma2_exists_msk(self%pinds(iptcl)))then
                call self%gencorrs_euclid_2(iref, iptcl, shvec, cc)
            else
                call self%gencorrs_cc_2(iref, iptcl, shvec, cc)
            end if
        end if
    end subroutine gencorrs_2

    function gencorr_for_rot_8( self, iref, iptcl, shvec, irot ) result( cc )
        class(polarft_corrcalc), intent(inout) :: self
        integer,                 intent(in)    :: iref, iptcl
        real(dp),                intent(in)    :: shvec(2)
        integer,                 intent(in)    :: irot
        real(dp)                               :: cc
        if     ( params_glob%cc_objfun == OBJFUN_CC )then
            cc = self%gencorr_cc_for_rot_8( iref, iptcl, shvec, irot )
        else if( params_glob%cc_objfun == OBJFUN_RES )then
            cc = self%gencorr_resnorm_for_rot_8( iref, iptcl, shvec, irot )
        else if( params_glob%cc_objfun == OBJFUN_EUCLID )then
            if (self%sigma2_exists_msk(self%pinds(iptcl)))then
                cc = self%gencorr_euclid_for_rot_8( iref, iptcl, shvec, irot )
            else
                cc = self%gencorr_cc_for_rot_8( iref, iptcl, shvec, irot )
            end if
        end if
    end function gencorr_for_rot_8

    function gencorr_cc_for_rot( self, iref, iptcl, shvec, irot ) result( cc )
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
                pft_ref = (self%pfts_refs_even(:,:,iref) * self%ctfmats(:,:,self%pinds(iptcl))) * shmat
            else
                pft_ref = (self%pfts_refs_odd (:,:,iref) * self%ctfmats(:,:,self%pinds(iptcl))) * shmat
            endif
        else
            if( self%iseven(self%pinds(iptcl)) )then
                pft_ref = self%pfts_refs_even(:,:,iref) * shmat
            else
                pft_ref = self%pfts_refs_odd (:,:,iref) * shmat
            endif
        endif
        sqsum_ref = sum(csq(pft_ref))
        corr      = self%calc_corr_for_rot(pft_ref, self%pinds(iptcl), irot)
        cc        = corr  / sqrt(sqsum_ref * self%sqsums_ptcls(self%pinds(iptcl)))
    end function gencorr_cc_for_rot

    function gencorr_cc_for_rot_8( self, iref, iptcl, shvec, irot ) result( cc )
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
                pft_ref = (self%pfts_refs_even(:,:,iref) * self%ctfmats(:,:,self%pinds(iptcl))) * shmat
            else
                pft_ref = (self%pfts_refs_odd(:,:,iref)  * self%ctfmats(:,:,self%pinds(iptcl))) * shmat
            endif
        else
            if( self%iseven(self%pinds(iptcl)) )then
                pft_ref = self%pfts_refs_even(:,:,iref) * shmat
            else
                pft_ref = self%pfts_refs_odd(:,:,iref)  * shmat
            endif
        endif
        sqsum_ref = sum(csq(pft_ref))
        corr      = self%calc_corr_for_rot_8(pft_ref, self%pinds(iptcl), irot)
        cc        = corr  / sqrt(sqsum_ref * self%sqsums_ptcls(self%pinds(iptcl)))
    end function gencorr_cc_for_rot_8

    function gencorr_cont_grad_cc_for_rot_8( self, iref, iptcl, shvec, irot, dcc ) result( cc )
        class(polarft_corrcalc), intent(inout) :: self
        integer,                 intent(in)    :: iref, iptcl
        real(dp),                intent(in)    :: shvec(2)
        integer,                 intent(in)    :: irot
        real(dp),                intent(out)   :: dcc(3)
        complex(dp), pointer :: pft_ref(:,:), shmat(:,:), pft_dref(:,:,:)
        real(dp),    pointer :: argmat(:,:)
        real(dp),    pointer :: fdf_y(:), fdf_T1(:,:), fdf_T2(:,:)
        real(dp) :: cc, sqsum_ref, denom
        integer  :: ithr, j, k
        ithr     =  omp_get_thread_num() + 1
        pft_ref  => self%heap_vars(ithr)%pft_ref_8
        pft_dref => self%heap_vars(ithr)%pft_dref_8
        shmat    => self%heap_vars(ithr)%shmat_8
        argmat   => self%heap_vars(ithr)%argmat_8
        fdf_y    => self%heap_vars(ithr)%fdf_y_8
        fdf_T1   => self%heap_vars(ithr)%fdf_T1_8
        fdf_T2   => self%heap_vars(ithr)%fdf_T2_8
        argmat   =  self%argtransf(:self%pftsz,:) * shvec(1) + self%argtransf(self%pftsz + 1:,:) * shvec(2)
        shmat    =  cmplx(cos(argmat),sin(argmat),dp)
        if( self%with_ctf ) then
            if( self%iseven(self%pinds(iptcl)) ) then
                pft_ref  = (self%pfts_refs_even(:,:,iref)    * self%ctfmats(:,:,self%pinds(iptcl))) * shmat
                pft_dref = self%pfts_drefs_even(:,:,:,iref)
            else
                pft_ref = (self%pfts_refs_odd(:,:,iref)  * self%ctfmats(:,:,self%pinds(iptcl))) * shmat
                pft_dref = self%pfts_drefs_odd(:,:,:,iref)
            endif
            do j = 1,3
                pft_dref(:,:,j) = (pft_dref(:,:,j) * self%ctfmats(:,:,self%pinds(iptcl))) * shmat
            end do
        else
            if( self%iseven(self%pinds(iptcl)) ) then
                pft_ref  = self%pfts_refs_even(:,:,iref) * shmat
                pft_dref = self%pfts_drefs_even(:,:,:,iref)
            else
                pft_ref  = self%pfts_refs_odd(:,:,iref)  * shmat
                pft_dref = self%pfts_drefs_odd(:,:,:,iref)
            endif
            do j = 1,3
                pft_dref(:,:,j) = pft_dref(:,:,j) * shmat
            end do
        endif
        sqsum_ref = sum(csq(pft_ref))
        denom    = sqrt(sqsum_ref * self%sqsums_ptcls(self%pinds(iptcl)))
        do k = params_glob%kfromto(1),params_glob%kstop
            fdf_y(k) = self%calc_corrk_for_rot_8(pft_ref, self%pinds(iptcl), k, irot)
        end do
        call self%calc_T1_T2_for_rot_8( pft_ref, pft_dref, iref, irot, 3, fdf_T1, fdf_T2)
        cc = sum(fdf_y) / denom
        do j = 1,3
            dcc(j) = ( sum(fdf_T1(:,j)) - sum(fdf_y) * sum(fdf_T2(:,j)) / sqsum_ref ) / denom
        end do
    end function gencorr_cont_grad_cc_for_rot_8

    subroutine gencorr_cont_shift_grad_cc_for_rot_8( self, iref, iptcl, shvec, irot, f, grad )
        class(polarft_corrcalc), intent(inout) :: self
        integer,                 intent(in)    :: iref, iptcl
        real(dp),                intent(in)    :: shvec(2)
        integer,                 intent(in)    :: irot
        real(dp),                intent(out)   :: f
        real(dp),                intent(out)   :: grad(5) ! 3 orientation angles, 2 shifts
        complex(dp), pointer :: pft_ref(:,:), pft_ref_tmp(:,:), shmat(:,:), pft_dref(:,:,:)
        real(dp),    pointer :: argmat(:,:)
        real(dp),    pointer :: fdf_y(:), fdf_T1(:,:), fdf_T2(:,:)
        real(dp) :: sqsum_ref, denom, corr
        integer  :: ithr, j, k
        ithr        =  omp_get_thread_num() + 1
        pft_ref     => self%heap_vars(ithr)%pft_ref_8
        pft_ref_tmp => self%heap_vars(ithr)%pft_ref_tmp_8
        pft_dref    => self%heap_vars(ithr)%pft_dref_8
        shmat       => self%heap_vars(ithr)%shmat_8
        argmat      => self%heap_vars(ithr)%argmat_8
        fdf_y       => self%heap_vars(ithr)%fdf_y_8
        fdf_T1      => self%heap_vars(ithr)%fdf_T1_8
        fdf_T2      => self%heap_vars(ithr)%fdf_T2_8
        argmat      =  self%argtransf(:self%pftsz,:) * shvec(1) + self%argtransf(self%pftsz + 1:,:) * shvec(2)
        shmat       =  cmplx(cos(argmat),sin(argmat),dp)
        if( self%with_ctf ) then
            if( self%iseven(self%pinds(iptcl)) ) then
                pft_ref  = (self%pfts_refs_even(:,:,iref)    * self%ctfmats(:,:,self%pinds(iptcl))) * shmat
                pft_dref = self%pfts_drefs_even(:,:,:,iref)
            else
                pft_ref = (self%pfts_refs_odd(:,:,iref)  * self%ctfmats(:,:,self%pinds(iptcl))) * shmat
                pft_dref = self%pfts_drefs_odd(:,:,:,iref)
            endif
            do j = 1,3
                pft_dref(:,:,j) = (pft_dref(:,:,j) * self%ctfmats(:,:,self%pinds(iptcl))) * shmat
            end do
        else
            if( self%iseven(self%pinds(iptcl)) ) then
                pft_ref  = self%pfts_refs_even(:,:,iref) * shmat
                pft_dref = self%pfts_drefs_even(:,:,:,iref)
            else
                pft_ref  = self%pfts_refs_odd(:,:,iref)  * shmat
                pft_dref = self%pfts_drefs_odd(:,:,:,iref)
            endif
            do j = 1,3
                pft_dref(:,:,j) = pft_dref(:,:,j) * shmat
            end do
        endif
        sqsum_ref = sum(csq(pft_ref))
        denom    = sqrt(sqsum_ref * self%sqsums_ptcls(self%pinds(iptcl)))
        do k = params_glob%kfromto(1),params_glob%kstop
            fdf_y(k) = self%calc_corrk_for_rot_8(pft_ref, self%pinds(iptcl), k, irot)
        end do
        call self%calc_T1_T2_for_rot_8( pft_ref, pft_dref, iref, irot, 3, fdf_T1, fdf_T2)
        f = sum(fdf_y) / denom
        do j = 1,3
            grad(j) = ( sum(fdf_T1(:,j)) - sum(fdf_y) * sum(fdf_T2(:,j)) / sqsum_ref ) / denom
        end do
        pft_ref_tmp = pft_ref * (0., 1.) * self%argtransf(:self%pftsz,:)
        corr        = self%calc_corr_for_rot_8(pft_ref_tmp, self%pinds(iptcl), irot)
        grad(4)     = corr / denom
        pft_ref_tmp = pft_ref * (0., 1.) * self%argtransf(self%pftsz + 1:,:)
        corr        = self%calc_corr_for_rot_8(pft_ref_tmp, self%pinds(iptcl), irot)
        grad(5)     = corr / denom
    end subroutine gencorr_cont_shift_grad_cc_for_rot_8

    function gencorr_resnorm_for_rot_8( self, iref, iptcl, shvec, irot ) result( cc )
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
                pft_ref = (self%pfts_refs_even(:,:,iref) * self%ctfmats(:,:,self%pinds(iptcl))) * shmat
            else
                pft_ref = (self%pfts_refs_odd(:,:,iref)  * self%ctfmats(:,:,self%pinds(iptcl))) * shmat
            endif
        else
            if( self%iseven(self%pinds(iptcl)) )then
                pft_ref = self%pfts_refs_even(:,:,iref) * shmat
            else
                pft_ref = self%pfts_refs_odd(:,:,iref)  * shmat
            endif
        endif
        cc = 0.0_dp
        ! with B-factor weighting
        do k = params_glob%kfromto(1),params_glob%kstop
            sqsumk_ref  = sum(csq(pft_ref(:,k)))
            sqsumk_ptcl = sum(csq(self%pfts_ptcls(:,k,self%pinds(iptcl))))
            corrk       = self%calc_corrk_for_rot_8(pft_ref, self%pinds(iptcl), k, irot)
            cc          = cc + (corrk * real(self%ptcl_bfac_weights(k,self%pinds(iptcl)), kind=dp)) / sqrt(sqsumk_ref * sqsumk_ptcl)
        end do
        cc = cc / real(self%ptcl_bfac_norms(self%pinds(iptcl)), kind=dp)
    end function gencorr_resnorm_for_rot_8

    subroutine gencorr_grad_for_rot_8( self, iref, iptcl, shvec, irot, f, grad )
        class(polarft_corrcalc), intent(inout) :: self
        integer,                 intent(in)    :: iref, iptcl
        real(dp),                intent(in)    :: shvec(2)
        integer,                 intent(in)    :: irot
        real(dp),                intent(out)   :: f, grad(2)
        if     ( params_glob%cc_objfun == OBJFUN_CC )then
            call self%gencorr_cc_grad_for_rot_8( iref, iptcl, shvec, irot, f, grad )
        else if( params_glob%cc_objfun == OBJFUN_RES)then
            call self%gencorr_resnorm_grad_for_rot_8( iref, iptcl, shvec, irot, f, grad )
        else if( params_glob%cc_objfun == OBJFUN_EUCLID)then
            if (self%sigma2_exists_msk(self%pinds(iptcl)))then
                call self%gencorr_euclid_grad_for_rot_8( iref, iptcl, shvec, irot, f, grad )
            else
                call self%gencorr_cc_grad_for_rot_8( iref, iptcl, shvec, irot, f, grad )
            end if
        end if
    end subroutine gencorr_grad_for_rot_8

    subroutine gencorr_cc_grad_for_rot_8( self, iref, iptcl, shvec, irot, f, grad )
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
                pft_ref = (self%pfts_refs_even(:,:,iref) * self%ctfmats(:,:,self%pinds(iptcl))) * shmat
            else
                pft_ref = (self%pfts_refs_odd(:,:,iref)  * self%ctfmats(:,:,self%pinds(iptcl))) * shmat
            endif
        else
            if( self%iseven(self%pinds(iptcl)) )then
                pft_ref = self%pfts_refs_even(:,:,iref) * shmat
            else
                pft_ref = self%pfts_refs_odd(:,:,iref)  * shmat
            endif
        endif
        denom       = sqrt(sum(csq(pft_ref)) * self%sqsums_ptcls(self%pinds(iptcl)))
        corr        = self%calc_corr_for_rot_8(pft_ref, self%pinds(iptcl), irot)
        f           = corr  / denom
        pft_ref_tmp = pft_ref * (0., 1.) * self%argtransf(:self%pftsz,:)
        corr        = self%calc_corr_for_rot_8(pft_ref_tmp, self%pinds(iptcl), irot)
        grad(1)     = corr / denom
        pft_ref_tmp = pft_ref * (0., 1.) * self%argtransf(self%pftsz + 1:,:)
        corr        = self%calc_corr_for_rot_8(pft_ref_tmp, self%pinds(iptcl), irot)
        grad(2)     = corr / denom
    end subroutine gencorr_cc_grad_for_rot_8

    subroutine gencorr_resnorm_grad_for_rot_8( self, iref, iptcl, shvec, irot, f, grad )
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
                pft_ref = (self%pfts_refs_even(:,:,iref) * self%ctfmats(:,:,self%pinds(iptcl))) * shmat
            else
                pft_ref = (self%pfts_refs_odd(:,:,iref)  * self%ctfmats(:,:,self%pinds(iptcl))) * shmat
            endif
        else
            if( self%iseven(self%pinds(iptcl)) )then
                pft_ref = self%pfts_refs_even(:,:,iref) * shmat
            else
                pft_ref = self%pfts_refs_odd (:,:,iref)  * shmat
            endif
        endif
        f    = 0.0_dp
        grad = 0.0_dp
        pft_ref_tmp1 = pft_ref * (0., 1.) * self%argtransf(:self%pftsz,    :)
        pft_ref_tmp2 = pft_ref * (0., 1.) * self%argtransf(self%pftsz + 1:,:)
        do k = params_glob%kfromto(1), params_glob%kstop
            sqsumk_ref  = sum(csq(pft_ref(:,k)))
            sqsumk_ptcl = sum(csq(self%pfts_ptcls(:,k,self%pinds(iptcl))))
            denom       = sqrt(sqsumk_ref * sqsumk_ptcl)
            corrk       = self%calc_corrk_for_rot_8(pft_ref, self%pinds(iptcl), k, irot)
            f           = f +       (corrk * real(self%ptcl_bfac_weights(k,self%pinds(iptcl)), kind=dp)) / denom
            corrk       = self%calc_corrk_for_rot_8(pft_ref_tmp1, self%pinds(iptcl), k, irot)
            grad(1)     = grad(1) + (corrk * real(self%ptcl_bfac_weights(k,self%pinds(iptcl)), kind=dp)) / denom
            corrk       = self%calc_corrk_for_rot_8(pft_ref_tmp2, self%pinds(iptcl), k, irot)
            grad(2)     = grad(2) + (corrk * real(self%ptcl_bfac_weights(k,self%pinds(iptcl)), kind=dp)) / denom
        end do
        f       = f       / real(self%ptcl_bfac_norms(self%pinds(iptcl)), kind=dp)
        grad(:) = grad(:) / real(self%ptcl_bfac_norms(self%pinds(iptcl)), kind=dp)
    end subroutine gencorr_resnorm_grad_for_rot_8

    function gencorr_cont_grad_resnorm_for_rot_8( self, iref, iptcl, shvec, irot, dcc ) result( cc )
        class(polarft_corrcalc), intent(inout) :: self
        integer,                 intent(in)    :: iref, iptcl
        real(dp),                intent(in)    :: shvec(2)
        integer,                 intent(in)    :: irot
        real(dp),                intent(out)   :: dcc(3)
        complex(dp), pointer :: pft_ref(:,:), shmat(:,:), pft_dref(:,:,:)
        real(dp),    pointer :: argmat(:,:)
        real(dp),    pointer :: fdf_T1(:,:), fdf_T2(:,:)
        real(dp) :: cc, sqsumk_ref, sqsumk_ptcl, fdf_yk
        real(dp) :: bfac_weight, denomk
        integer  :: ithr, k, j
        ithr    =  omp_get_thread_num() + 1
        pft_ref => self%heap_vars(ithr)%pft_ref_8
        pft_dref => self%heap_vars(ithr)%pft_dref_8
        shmat   => self%heap_vars(ithr)%shmat_8
        argmat  => self%heap_vars(ithr)%argmat_8
        fdf_T1   => self%heap_vars(ithr)%fdf_T1_8
        fdf_T2   => self%heap_vars(ithr)%fdf_T2_8
        argmat  =  self%argtransf(:self%pftsz,:) * shvec(1) + self%argtransf(self%pftsz + 1:,:) * shvec(2)
        shmat   =  cmplx(cos(argmat),sin(argmat),dp)
        if( self%with_ctf ) then
            if( self%iseven(self%pinds(iptcl)) ) then
                pft_ref = (self%pfts_refs_even(:,:,iref) * self%ctfmats(:,:,self%pinds(iptcl))) * shmat
                pft_dref = self%pfts_drefs_even(:,:,:,iref)
            else
                pft_ref = (self%pfts_refs_odd(:,:,iref)  * self%ctfmats(:,:,self%pinds(iptcl))) * shmat
                pft_dref = self%pfts_drefs_odd(:,:,:,iref)
            endif
            do j = 1,3
                pft_dref(:,:,j) = (pft_dref(:,:,j) * self%ctfmats(:,:,self%pinds(iptcl))) * shmat
            end do
        else
            if( self%iseven(self%pinds(iptcl)) ) then
                pft_ref = self%pfts_refs_even(:,:,iref) * shmat
                pft_dref = self%pfts_drefs_even(:,:,:,iref)
            else
                pft_ref = self%pfts_refs_odd(:,:,iref)  * shmat
                pft_dref = self%pfts_drefs_odd(:,:,:,iref)
            endif
            do j = 1,3
                pft_dref(:,:,j) = pft_dref(:,:,j) * shmat
            end do
        endif
        call self%calc_T1_T2_for_rot_8( pft_ref, pft_dref, iref, irot, 3, fdf_T1, fdf_T2)
        cc  = 0.0_dp
        dcc = 0.0_dp
        ! with B-factor weighting
        do k = params_glob%kfromto(1),params_glob%kstop
            sqsumk_ref  = sum(csq(pft_ref(:,k)))
            sqsumk_ptcl = sum(csq(self%pfts_ptcls(:,k,self%pinds(iptcl))))
            denomk      = sqrt(sqsumk_ref * sqsumk_ptcl)
            fdf_yk       = self%calc_corrk_for_rot_8(pft_ref, self%pinds(iptcl), k, irot)
            bfac_weight = real(self%ptcl_bfac_weights(k,self%pinds(iptcl)), kind=dp)
            cc          = cc + (fdf_yk * bfac_weight) / denomk
            do j = 1,3
                dcc(j) = dcc(j) + bfac_weight * &
                    &( fdf_T1(k,j) - fdf_yk * fdf_T2(k,j) / sqsumk_ref) / denomk
            end do
        end do
        cc  = cc  / real(self%ptcl_bfac_norms(self%pinds(iptcl)), kind=dp)
        dcc = dcc / real(self%ptcl_bfac_norms(self%pinds(iptcl)), kind=dp)
    end function gencorr_cont_grad_resnorm_for_rot_8

    subroutine gencorr_cont_shift_grad_resnorm_for_rot_8( self, iref, iptcl, shvec, irot, f, grad)
        class(polarft_corrcalc), intent(inout) :: self
        integer,                 intent(in)    :: iref, iptcl
        real(dp),                intent(in)    :: shvec(2)
        integer,                 intent(in)    :: irot
        real(dp),                intent(out)   :: f
        real(dp),                intent(out)   :: grad(5) ! 3 orientation angles, 2 shifts
        complex(dp), pointer :: pft_ref(:,:), pft_ref_tmp1(:,:), pft_ref_tmp2(:,:), shmat(:,:), pft_dref(:,:,:)
        real(dp),    pointer :: argmat(:,:)
        real(dp),    pointer :: fdf_T1(:,:), fdf_T2(:,:)
        real(dp) :: sqsumk_ref, sqsumk_ptcl, fdf_yk
        real(dp) :: bfac_weight, denomk
        integer  :: ithr, k, j
        ithr    =  omp_get_thread_num() + 1
        pft_ref => self%heap_vars(ithr)%pft_ref_8
        pft_dref => self%heap_vars(ithr)%pft_dref_8
        shmat   => self%heap_vars(ithr)%shmat_8
        argmat  => self%heap_vars(ithr)%argmat_8
        fdf_T1   => self%heap_vars(ithr)%fdf_T1_8
        fdf_T2   => self%heap_vars(ithr)%fdf_T2_8
        argmat  =  self%argtransf(:self%pftsz,:) * shvec(1) + self%argtransf(self%pftsz + 1:,:) * shvec(2)
        shmat   =  cmplx(cos(argmat),sin(argmat),dp)
        if( self%with_ctf ) then
            if( self%iseven(self%pinds(iptcl)) ) then
                pft_ref = (self%pfts_refs_even(:,:,iref) * self%ctfmats(:,:,self%pinds(iptcl))) * shmat
                pft_dref = self%pfts_drefs_even(:,:,:,iref)
            else
                pft_ref = (self%pfts_refs_odd(:,:,iref)  * self%ctfmats(:,:,self%pinds(iptcl))) * shmat
                pft_dref = self%pfts_drefs_odd(:,:,:,iref)
            endif
            do j = 1,3
                pft_dref(:,:,j) = (pft_dref(:,:,j) * self%ctfmats(:,:,self%pinds(iptcl))) * shmat
            end do
        else
            if( self%iseven(self%pinds(iptcl)) ) then
                pft_ref = self%pfts_refs_even(:,:,iref) * shmat
                pft_dref = self%pfts_drefs_even(:,:,:,iref)
            else
                pft_ref = self%pfts_refs_odd(:,:,iref)  * shmat
                pft_dref = self%pfts_drefs_odd(:,:,:,iref)
            endif
            do j = 1,3
                pft_dref(:,:,j) = pft_dref(:,:,j) * shmat
            end do
        endif
        call self%calc_T1_T2_for_rot_8( pft_ref, pft_dref, iref, irot, 3, fdf_T1, fdf_T2)
        f    = 0.0_dp
        grad = 0.0_dp
        pft_ref_tmp1 = pft_ref * (0., 1.) * self%argtransf(:self%pftsz,    :)
        pft_ref_tmp2 = pft_ref * (0., 1.) * self%argtransf(self%pftsz + 1:,:)
        ! with B-factor weighting
        do k = params_glob%kfromto(1),params_glob%kstop
            sqsumk_ref  = sum(csq(pft_ref(:,k)))
            sqsumk_ptcl = sum(csq(self%pfts_ptcls(:,k,self%pinds(iptcl))))
            denomk      = sqrt(sqsumk_ref * sqsumk_ptcl)
            fdf_yk      = self%calc_corrk_for_rot_8(pft_ref, self%pinds(iptcl), k, irot)
            bfac_weight = real(self%ptcl_bfac_weights(k,self%pinds(iptcl)), kind=dp)
            f           = f + (fdf_yk * bfac_weight) / denomk
            do j = 1,3
                grad(j) = grad(j) + bfac_weight * &
                    &( fdf_T1(k,j) - fdf_yk * fdf_T2(k,j) / sqsumk_ref) / denomk
            end do
            fdf_yk      = self%calc_corrk_for_rot_8(pft_ref_tmp1, self%pinds(iptcl), k, irot)
            grad(4)     = grad(1) + (fdf_yk * real(self%ptcl_bfac_weights(k,self%pinds(iptcl)), kind=dp)) / denomk
            fdf_yk      = self%calc_corrk_for_rot_8(pft_ref_tmp2, self%pinds(iptcl), k, irot)
            grad(5)     = grad(2) + (fdf_yk * real(self%ptcl_bfac_weights(k,self%pinds(iptcl)), kind=dp)) / denomk
        end do
        f       = f       / real(self%ptcl_bfac_norms(self%pinds(iptcl)), kind=dp)
        grad(:) = grad(:) / real(self%ptcl_bfac_norms(self%pinds(iptcl)), kind=dp)
    end subroutine gencorr_cont_shift_grad_resnorm_for_rot_8

    subroutine gencorr_grad_only_for_rot_8( self, iref, iptcl, shvec, irot, grad )
        class(polarft_corrcalc), intent(inout) :: self
        integer,                 intent(in)    :: iref, iptcl
        real(dp),                intent(in)    :: shvec(2)
        integer,                 intent(in)    :: irot
        real(dp),                intent(out)   :: grad(2)
        if     ( params_glob%cc_objfun == OBJFUN_CC )then
            call self%gencorr_cc_grad_only_for_rot_8( iref, iptcl, shvec, irot, grad )
        else if( params_glob%cc_objfun == OBJFUN_RES )then
            call self%gencorr_resnorm_grad_only_for_rot_8( iref, iptcl, shvec, irot, grad )
        else if( params_glob%cc_objfun == OBJFUN_EUCLID )then
            if (self%sigma2_exists_msk(self%pinds(iptcl))) then
                call self%gencorr_euclid_grad_only_for_rot_8( iref, iptcl, shvec, irot, grad )
            else
                call self%gencorr_cc_grad_only_for_rot_8( iref, iptcl, shvec, irot, grad )
            end if
        end if
    end subroutine gencorr_grad_only_for_rot_8

    subroutine gencorr_cc_grad_only_for_rot_8( self, iref, iptcl, shvec, irot, grad )
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
                pft_ref = (self%pfts_refs_even(:,:,iref) * self%ctfmats(:,:,self%pinds(iptcl))) * shmat
            else
                pft_ref = (self%pfts_refs_odd(:,:,iref)  * self%ctfmats(:,:,self%pinds(iptcl))) * shmat
            endif
        else
            if( self%iseven(self%pinds(iptcl)) )then
                pft_ref = self%pfts_refs_even(:,:,iref) * shmat
            else
                pft_ref = self%pfts_refs_odd(:,:,iref)  * shmat
            endif
        endif
        denom       = sqrt(sum(csq(pft_ref)) * self%sqsums_ptcls(self%pinds(iptcl)))
        pft_ref_tmp = pft_ref * (0., 1.) * self%argtransf(:self%pftsz,:)
        corr        = self%calc_corr_for_rot_8(pft_ref_tmp, self%pinds(iptcl), irot)
        grad(1)     = corr / denom
        pft_ref_tmp = pft_ref * (0., 1.) * self%argtransf(self%pftsz + 1:,:)
        corr        = self%calc_corr_for_rot_8(pft_ref_tmp, self%pinds(iptcl), irot)
        grad(2)     = corr / denom
    end subroutine gencorr_cc_grad_only_for_rot_8

    subroutine gencorr_resnorm_grad_only_for_rot_8( self, iref, iptcl, shvec, irot, grad )
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
                pft_ref = (self%pfts_refs_even(:,:,iref) * self%ctfmats(:,:,self%pinds(iptcl))) * shmat
            else
                pft_ref = (self%pfts_refs_odd(:,:,iref)  * self%ctfmats(:,:,self%pinds(iptcl))) * shmat
            endif
        else
            if( self%iseven(self%pinds(iptcl)) )then
                pft_ref = self%pfts_refs_even(:,:,iref) * shmat
            else
                pft_ref = self%pfts_refs_odd(:,:,iref)  * shmat
            endif
        endif
        pft_ref_tmp1 = pft_ref * (0., 1.) * self%argtransf(:self%pftsz,:)
        pft_ref_tmp2 = pft_ref * (0., 1.) * self%argtransf(self%pftsz + 1:,:)
        grad = 0.0_dp
        do k = params_glob%kfromto(1), params_glob%kstop
            sqsumk_ref  = sum(csq(pft_ref(:,k)))
            sqsumk_ptcl = sum(csq(self%pfts_ptcls(:,k,self%pinds(iptcl))))
            denom       = sqrt(sqsumk_ref * sqsumk_ptcl)
            corrk       = self%calc_corrk_for_rot_8(pft_ref_tmp1, self%pinds(iptcl), k, irot)
            grad(1)     = grad(1) + (corrk * real(self%ptcl_bfac_weights(k,self%pinds(iptcl)), kind=dp)) / denom
            corrk       = self%calc_corrk_for_rot_8(pft_ref_tmp2, self%pinds(iptcl), k, irot)
            grad(2)     = grad(2) + (corrk * real(self%ptcl_bfac_weights(k,self%pinds(iptcl)), kind=dp)) / denom
        end do
        grad = grad / real(self%ptcl_bfac_norms(self%pinds(iptcl)), kind=dp)
    end subroutine gencorr_resnorm_grad_only_for_rot_8

    function gencorr_euclid_for_rot( self, iref, iptcl, shvec, irot ) result( cc )
        class(polarft_corrcalc), intent(inout) :: self
        integer,                 intent(in)    :: iref, iptcl
        real(sp),                intent(in)    :: shvec(2)
        integer,                 intent(in)    :: irot
        complex(sp), pointer :: pft_ref(:,:), shmat(:,:)
        real(sp),    pointer :: argmat(:,:)
        real(sp) :: cc
        integer  :: ithr
        ithr    =  omp_get_thread_num() + 1
        pft_ref => self%heap_vars(ithr)%pft_ref
        shmat   => self%heap_vars(ithr)%shmat
        argmat  => self%heap_vars(ithr)%argmat
        argmat  =  self%argtransf(:self%pftsz,:) * shvec(1) + self%argtransf(self%pftsz + 1:,:) * shvec(2)
        shmat   =  cmplx(cos(argmat),sin(argmat))
        if( self%with_ctf )then
            if( self%iseven(self%pinds(iptcl)) )then
                pft_ref = (self%pfts_refs_even(:,:,iref) * self%ctfmats(:,:,self%pinds(iptcl))) * shmat
            else
                pft_ref = (self%pfts_refs_odd (:,:,iref) * self%ctfmats(:,:,self%pinds(iptcl))) * shmat
            endif
        else
            if( self%iseven(self%pinds(iptcl)) )then
                pft_ref = self%pfts_refs_even(:,:,iref) * shmat
            else
                pft_ref = self%pfts_refs_odd (:,:,iref) * shmat
            endif
        endif
        cc = self%calc_euclid_for_rot(pft_ref, self%pinds(iptcl), irot)
    end function gencorr_euclid_for_rot

    function gencorr_euclid_for_rot_8( self, iref, iptcl, shvec, irot ) result( cc )
        class(polarft_corrcalc), intent(inout) :: self
        integer,                 intent(in)    :: iref, iptcl
        real(dp),                intent(in)    :: shvec(2)
        integer,                 intent(in)    :: irot
        complex(dp), pointer :: pft_ref(:,:), shmat(:,:)
        real(dp),    pointer :: argmat(:,:)
        real(dp) :: cc
        integer  :: ithr
        ithr    =  omp_get_thread_num() + 1
        pft_ref => self%heap_vars(ithr)%pft_ref_8
        shmat   => self%heap_vars(ithr)%shmat_8
        argmat  => self%heap_vars(ithr)%argmat_8
        argmat  =  self%argtransf(:self%pftsz,:) * shvec(1) + self%argtransf(self%pftsz + 1:,:) * shvec(2)
        shmat   =  cmplx(cos(argmat),sin(argmat),dp)
        if( self%with_ctf )then
            if( self%iseven(self%pinds(iptcl)) )then
                pft_ref = (self%pfts_refs_even(:,:,iref) * self%ctfmats(:,:,self%pinds(iptcl))) * shmat
            else
                pft_ref = (self%pfts_refs_odd(:,:,iref)  * self%ctfmats(:,:,self%pinds(iptcl))) * shmat
            endif
        else
            if( self%iseven(self%pinds(iptcl)) )then
                pft_ref = self%pfts_refs_even(:,:,iref) * shmat
            else
                pft_ref = self%pfts_refs_odd(:,:,iref)  * shmat
            endif
        endif
        cc = self%calc_euclid_for_rot_8(pft_ref, self%pinds(iptcl), irot)
    end function gencorr_euclid_for_rot_8

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

    subroutine gencorr_euclid_grad_for_rot_8( self, iref, iptcl, shvec, irot, f, grad )
        class(polarft_corrcalc), intent(inout) :: self
        integer,                 intent(in)    :: iref, iptcl
        real(dp),                intent(in)    :: shvec(2)
        integer,                 intent(in)    :: irot
        real(dp),                intent(out)   :: f, grad(2)
        complex(dp), pointer :: pft_ref(:,:), pft_ref_tmp(:,:), shmat(:,:)
        real(dp),    pointer :: argmat(:,:)
        real(dp) :: euclid, denom, corr
        integer  :: ithr, k
        ithr        =  omp_get_thread_num() + 1
        pft_ref     => self%heap_vars(ithr)%pft_ref_8
        pft_ref_tmp => self%heap_vars(ithr)%pft_ref_tmp_8
        shmat       => self%heap_vars(ithr)%shmat_8
        argmat      => self%heap_vars(ithr)%argmat_8
        argmat      =  self%argtransf(:self%pftsz,:) * shvec(1) + self%argtransf(self%pftsz + 1:,:) * shvec(2)
        shmat       =  cmplx(cos(argmat),sin(argmat),dp)
        if( self%with_ctf )then
            if( self%iseven(self%pinds(iptcl)) )then
                pft_ref = (self%pfts_refs_even(:,:,iref) * self%ctfmats(:,:,self%pinds(iptcl))) * shmat
            else
                pft_ref = (self%pfts_refs_odd(:,:,iref)  * self%ctfmats(:,:,self%pinds(iptcl))) * shmat
            endif
        else
            if( self%iseven(self%pinds(iptcl)) )then
                pft_ref = self%pfts_refs_even(:,:,iref) * shmat
            else
                pft_ref = self%pfts_refs_odd(:,:,iref)  * shmat
            endif
        endif
        euclid      = self%calc_euclid_for_rot_8(pft_ref, self%pinds(iptcl), irot)
        f           = euclid
        grad(1)     = 0._dp
        pft_ref_tmp = pft_ref * (0., 1.) * self%argtransf(:self%pftsz,:)
        do k = params_glob%kfromto(1), params_glob%kstop
            grad(1) = grad(1) + &
                self%calc_corrk_for_rot_8(pft_ref_tmp, self%pinds(iptcl), k, irot) &
                / self%sigma2_noise(k, self%pinds(iptcl))
            grad(1) = grad(1) - real(sum(pft_ref_tmp(:,k)*conjg(pft_ref(:,k)) )) &
                / self%sigma2_noise(k, self%pinds(iptcl))
        end do
        grad(2)     = 0._dp
        pft_ref_tmp = pft_ref * (0., 1.) * self%argtransf(self%pftsz + 1:,:)
        do k = params_glob%kfromto(1), params_glob%kstop
            grad(2) = grad(2) + &
                self%calc_corrk_for_rot_8(pft_ref_tmp, self%pinds(iptcl), k, irot) &
                / self%sigma2_noise(k, self%pinds(iptcl))
            grad(2) = grad(2) - real(sum(pft_ref_tmp(:,k) * conjg(pft_ref(:,k)) )) &
                / self%sigma2_noise(k, self%pinds(iptcl))
        end do
    end subroutine gencorr_euclid_grad_for_rot_8

    subroutine gencorr_euclid_grad_only_for_rot_8( self, iref, iptcl, shvec, irot, grad )
        class(polarft_corrcalc), intent(inout) :: self
        integer,                 intent(in)    :: iref, iptcl
        real(dp),                intent(in)    :: shvec(2)
        integer,                 intent(in)    :: irot
        real(dp),                intent(out)   :: grad(2)
        complex(dp), pointer :: pft_ref(:,:), pft_ref_tmp(:,:), shmat(:,:)
        real(dp),    pointer :: argmat(:,:)
        real(dp) :: euclid, corr, denom
        integer  :: ithr, k
        ithr        =  omp_get_thread_num() + 1
        pft_ref     => self%heap_vars(ithr)%pft_ref_8
        pft_ref_tmp => self%heap_vars(ithr)%pft_ref_tmp_8
        shmat       => self%heap_vars(ithr)%shmat_8
        argmat      => self%heap_vars(ithr)%argmat_8
        argmat      =  self%argtransf(:self%pftsz,:) * shvec(1) + self%argtransf(self%pftsz + 1:,:) * shvec(2)
        shmat       =  cmplx(cos(argmat),sin(argmat),dp)
        if( self%with_ctf )then
            if( self%iseven(self%pinds(iptcl)) )then
                pft_ref = (self%pfts_refs_even(:,:,iref) * self%ctfmats(:,:,self%pinds(iptcl))) * shmat
            else
                pft_ref = (self%pfts_refs_odd(:,:,iref)  * self%ctfmats(:,:,self%pinds(iptcl))) * shmat
            endif
        else
            if( self%iseven(self%pinds(iptcl)) )then
                pft_ref = self%pfts_refs_even(:,:,iref) * shmat
            else
                pft_ref = self%pfts_refs_odd(:,:,iref)  * shmat
            endif
        endif
        grad(1)     = 0._dp
        pft_ref_tmp = pft_ref * (0., 1.) * self%argtransf(:self%pftsz,:)
        do k = params_glob%kfromto(1), params_glob%kstop
            grad(1) = grad(1) + &
                self%calc_corrk_for_rot_8(pft_ref_tmp, self%pinds(iptcl), k, irot) &
                / self%sigma2_noise(k, self%pinds(iptcl))
            grad(1) = grad(1) - real(sum(pft_ref_tmp(:,k)*conjg(pft_ref(:,k)) )) &
                / self%sigma2_noise(k, self%pinds(iptcl))
        end do
        grad(2)     = 0._dp
        pft_ref_tmp = pft_ref * (0., 1.) * self%argtransf(self%pftsz + 1:,:)
        do k = params_glob%kfromto(1), params_glob%kstop
            grad(2) = grad(2) + &
                self%calc_corrk_for_rot_8(pft_ref_tmp, self%pinds(iptcl), k, irot) &
                / self%sigma2_noise(k, self%pinds(iptcl))
            grad(2) = grad(2) - real(sum(pft_ref_tmp(:,k) * conjg(pft_ref(:,k)) )) &
                / self%sigma2_noise(k, self%pinds(iptcl))
        end do
    end subroutine gencorr_euclid_grad_only_for_rot_8

    subroutine gencorr_sigma_contrib( self, iref, iptcl, shvec, irot, sigma_contrib )
        class(polarft_corrcalc), intent(inout) :: self
        integer,                 intent(in)    :: iref, iptcl
        real(sp),                intent(in)    :: shvec(2)
        integer,                 intent(in)    :: irot
        real(sp),                intent(out)   :: sigma_contrib(params_glob%kfromto(1):params_glob%kfromto(2))
        complex(sp), pointer :: pft_ref(:,:), shmat(:,:)
        real(sp),    pointer :: argmat(:,:)
        real(sp) :: cc, corr
        integer  :: ithr, ik
        ithr    =  omp_get_thread_num() + 1
        pft_ref => self%heap_vars(ithr)%pft_ref
        shmat   => self%heap_vars(ithr)%shmat
        argmat  => self%heap_vars(ithr)%argmat
        argmat  =  self%argtransf(:self%pftsz,:) * shvec(1) + self%argtransf(self%pftsz + 1:,:) * shvec(2)
        shmat   =  cmplx(cos(argmat),sin(argmat))
        if( self%with_ctf )then
            if( self%iseven(self%pinds(iptcl)) )then
                pft_ref = (self%pfts_refs_even(:,:,iref) * self%ctfmats(:,:,self%pinds(iptcl))) * shmat
            else
                pft_ref = (self%pfts_refs_odd (:,:,iref) * self%ctfmats(:,:,self%pinds(iptcl))) * shmat
            endif
        else
            if( self%iseven(self%pinds(iptcl)) )then
                pft_ref = self%pfts_refs_even(:,:,iref) * shmat
            else
                pft_ref = self%pfts_refs_odd (:,:,iref) * shmat
            endif
        endif
        do ik = params_glob%kfromto(1), params_glob%kfromto(2)
            sigma_contrib(ik) = self%calc_euclidk_for_rot(pft_ref, self%pinds(iptcl), ik, irot)
        end do
        sigma_contrib = sigma_contrib !/ real(self%nrots)
    end subroutine gencorr_sigma_contrib

    real function specscore( self, iref, iptcl, irot )
        class(polarft_corrcalc), intent(inout) :: self
        integer,                 intent(in)    :: iref, iptcl, irot
        real :: frc(params_glob%kfromto(1):params_glob%kstop)
        call self%genfrc(iref, self%pinds(iptcl), irot, frc)
        specscore = max(0.,median_nocopy(frc))
    end function specscore

    real function fit_bfac( self, iref, iptcl, irot, shvec )
        ! Fitting to Y = A * exp( -B/(4s2) )
        ! from mathworld.wolfram.com/LeastSquaresFittingExponential.html
        class(polarft_corrcalc),  intent(inout) :: self
        integer,                  intent(in)    :: iref, iptcl, irot
        real(sp),                 intent(in)    :: shvec(2)
        real(sp) :: denom, frc(params_glob%kfromto(1):params_glob%kstop), logfrc(params_glob%kfromto(1):params_glob%kstop)
        real(sp) :: sumxfrc, sumfrc
        logical  :: peakmsk(params_glob%kfromto(1):params_glob%kstop)
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
        ! all points are weighted by Y
        sumfrc   = sum(frc, mask=peakmsk)
        sumxfrc  = sum(self%inv_resarrsq * frc, mask=peakmsk)
        where( frc > TINY )
            logfrc = log(frc)
        else where
            logfrc = 0.
        end where
        denom    = sumfrc * sum(self%inv_resarrsq * self%inv_resarrsq * frc, mask=peakmsk) - sumxfrc*sumxfrc
        fit_bfac = sumfrc * sum(self%inv_resarrsq * frc * logfrc, mask=peakmsk)
        fit_bfac = fit_bfac - sumxfrc * sum(frc*logfrc, mask=peakmsk)
        fit_bfac = fit_bfac / denom
        fit_bfac = max(0., fit_bfac)
    end function fit_bfac

    function calc_roinv_corrmat( self) result(corrmat )
        ! assumes particles/references in pftcc identical sets
        class(polarft_corrcalc), intent(inout) :: self
        real, allocatable   :: corrmat(:,:)
        real    :: corrs(self%nrots)
        integer :: iref, iptcl, loc(1)
        if( self%nptcls /= self%nrefs ) stop 'nptcls == nrefs in pftcc required; simple_corrmat :: calc_roinv_corrmat'
        if( allocated(corrmat) ) deallocate(corrmat)
        allocate(corrmat(self%nptcls,self%nptcls), stat=alloc_stat)
        if(alloc_stat/=0)call allocchk('In: calc_roinv_corrmat; simple_corrmat')
        corrmat = 1.
        !$omp parallel do default(shared) schedule(guided) private(iref,iptcl,corrs,loc) proc_bind(close)
        do iref=1,self%nptcls - 1
            do iptcl=iref + 1,self%nptcls
                ! rotational corr
                call self%gencorrs(iref, iptcl, corrs)
                loc  = maxloc(corrs)
                corrmat(iref,iptcl) = corrs(loc(1))
                ! symmetrize
                corrmat(iptcl,iref) = corrmat(iref,iptcl)
            end do
        end do
        !$omp end parallel do
    end function calc_roinv_corrmat

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
                    &self%heap_vars(ithr)%pft_ref_tmp1, self%heap_vars(ithr)%pft_ref_tmp2,&
                    &self%heap_vars(ithr)%pft_dref,&
                    &self%heap_vars(ithr)%corrs_over_k,self%heap_vars(ithr)%argmat,&
                    &self%heap_vars(ithr)%shmat,self%heap_vars(ithr)%kcorrs,&
                    &self%heap_vars(ithr)%pft_ref_8,self%heap_vars(ithr)%pft_ref_tmp_8,&
                    &self%heap_vars(ithr)%pft_ref_tmp1_8, self%heap_vars(ithr)%pft_ref_tmp2_8,&
                    &self%heap_vars(ithr)%pft_dref_8,&
                    &self%heap_vars(ithr)%shmat_8,self%heap_vars(ithr)%argmat_8,&
                    &self%heap_vars(ithr)%fdf_y_8,self%heap_vars(ithr)%fdf_T1_8,&
                    &self%heap_vars(ithr)%fdf_T2_8)
            end do
            do i = 1, self%nptcls
                do ik = params_glob%kfromto(1),params_glob%kstop
                    call fftwf_free(self%fftdat_ptcls(i,ik)%p_re)
                    call fftwf_free(self%fftdat_ptcls(i,ik)%p_im)
                end do
            end do
            if( allocated(self%ptcl_bfac_weights) ) deallocate(self%ptcl_bfac_weights)
            if( allocated(self%ptcl_bfac_norms)   ) deallocate(self%ptcl_bfac_norms)
            if( allocated(self%inv_resarrsq)      ) deallocate(self%inv_resarrsq)
            deallocate( self%sqsums_ptcls, self%angtab, self%argtransf,&
                &self%polar, self%pfts_refs_even, self%pfts_refs_odd, self%pfts_drefs_even, self%pfts_drefs_odd,&
                self%pfts_ptcls, self%fft_factors, self%fftdat, self%fftdat_ptcls,&
                &self%iseven, self%pinds)
            call fftwf_destroy_plan(self%plan_bwd)
            call fftwf_destroy_plan(self%plan_fwd_1)
            call fftwf_destroy_plan(self%plan_fwd_2)
            self%sigma2_noise      => null()
            self%sigma2_exists_msk => null()
            self%existence = .false.
        endif
    end subroutine kill

end module simple_polarft_corrcalc
