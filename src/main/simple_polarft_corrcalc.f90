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
    type(c_ptr)                            :: p_im                      !< -"-
    complex(kind=c_float_complex), pointer :: re(:) => null()           !< corresponding Fortran pointers
    complex(kind=c_float_complex), pointer :: im(:) => null()           !< -"-
end type fftw_carr_fft

type fftw_arrs
    type(c_ptr)                            :: p_ref_re                  !< pointer for C-style allocation
    type(c_ptr)                            :: p_ref_im                  !< -"-
    type(c_ptr)                            :: p_ref_fft_re              !< -"-
    type(c_ptr)                            :: p_ref_fft_im              !< -"-
    type(c_ptr)                            :: p_ref_fft_re_2            !< -"-
    type(c_ptr)                            :: p_ref_fft_im_2            !< -"-
    type(c_ptr)                            :: p_product_fft             !< -"-
    type(c_ptr)                            :: p_backtransf              !< -"-
    real(kind=c_float),            pointer :: ref_re(:)       => null()  !< corresponding Fortran pointers
    complex(kind=c_float_complex), pointer :: ref_im(:)       => null()  !< -"-
    complex(kind=c_float_complex), pointer :: ref_fft_re(:)   => null()  !< -"-
    complex(kind=c_float_complex), pointer :: ref_fft_im(:)   => null()  !< -"-
    complex(kind=c_float_complex), pointer :: ref_fft_re_2(:) => null()  !< -"-
    complex(kind=c_float_complex), pointer :: ref_fft_im_2(:) => null()  !< -"-
    complex(kind=c_float_complex), pointer :: product_fft(:)  => null()  !< -"-
    real(kind=c_float),            pointer :: backtransf(:)   => null()  !< -"-
end type fftw_arrs

type heap_vars
    complex(sp), pointer :: pft_ref(:,:)        => null()
    complex(sp), pointer :: pft_ref_tmp(:,:)    => null()
    complex(sp), pointer :: pft_ref_tmp1(:,:)   => null()
    complex(sp), pointer :: pft_ref_tmp2(:,:)   => null()
    complex(sp), pointer :: pft_dref(:,:,:)     => null()
    real,        pointer :: corrs_over_k(:)     => null()
    real,        pointer :: corrs_mir_over_k(:) => null()
    real(sp),    pointer :: argmat(:,:)         => null()
    complex(sp), pointer :: shmat(:,:)          => null()
    real,        pointer :: kcorrs(:)           => null()
    real,        pointer :: kcorrs_mir(:)       => null()
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
    !private
    integer                          :: nptcls     = 1        !< the total number of particles in partition (logically indexded [fromp,top])
    integer                          :: nrefs      = 1        !< the number of references (logically indexded [1,nrefs])
    integer                          :: nrots      = 0        !< number of in-plane rotations for one pft (determined by radius of molecule)
    integer                          :: pftsz      = 0        !< size of reference and particle pft (nrots/2)
    integer                          :: nthr       = 0        !< # OpenMP threads
    integer                          :: pfromto(2) = 0        !< particle index range
    integer                          :: winsz      = 0        !< size of moving window in correlation calculations
    integer                          :: ldim(3)    = 0        !< logical dimensions of original cartesian image
    integer,             allocatable :: pinds(:)              !< index array (to reduce memory when frac_update < 1)
    real(sp),            allocatable :: sqsums_ptcls(:)       !< memoized square sums for the correlation calculations (taken from kfromto(1):kstop)
    real(sp),            allocatable :: angtab(:)             !< table of in-plane angles (in degrees)
    real(sp),            allocatable :: argtransf(:,:)        !< argument transfer constants for shifting the references
    real(sp),            allocatable :: polar(:,:)            !< table of polar coordinates (in Cartesian coordinates)
    real(sp),            allocatable :: ctfmats(:,:,:)        !< expand set of CTF matrices (for efficient parallel exec)
    real(sp),            allocatable :: ref_optlp(:,:)        !< references optimal filter
    real(sp),            allocatable :: pssnr_filt(:,:)       !< filter for particle ssnr
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
    type(c_ptr)                      :: plan_fwd_1            !< FFTW plans for gencorrs
    type(c_ptr)                      :: plan_fwd_2            !< -"-
    type(c_ptr)                      :: plan_bwd              !< -"-
    logical                          :: l_clsfrcs   = .false. !< CLS2D/3DRefs flag
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
    procedure          :: set_ref_optlp
    procedure          :: set_pssnr_filt
    procedure          :: cp_even2odd_ref
    procedure          :: cp_even_ref2ptcl
    procedure          :: cp_refs
    procedure          :: swap_ptclsevenodd
    procedure          :: set_eos
    procedure          :: assign_sigma2_noise
    ! GETTERS
    procedure          :: get_nrots
    procedure          :: get_pdim
    procedure          :: get_pftsz
    procedure          :: get_box
    procedure          :: get_rot
    procedure          :: get_roind
    procedure          :: get_coord
    procedure          :: get_ref_pft
    procedure          :: get_nrefs
    procedure          :: exists
    procedure          :: ptcl_iseven
    procedure          :: get_nptcls
    procedure          :: assign_pinds
    ! PRINTERS/VISUALISERS
    procedure          :: print
    procedure          :: vis_ptcl
    procedure          :: vis_ref
    ! MODIFIERS
    procedure, private :: shellnorm_and_filter_ref
    procedure, private :: shellnorm_and_filter_ref_8
    procedure, private :: shellnorm_and_filter_ref_dref_8
    ! MEMOIZER
    procedure, private :: memoize_sqsum_ptcl
    procedure, private :: memoize_fft
    procedure          :: memoize_ffts
    ! CALCULATORS
    procedure          :: create_polar_absctfmats
    procedure, private :: prep_ref4corr
    procedure          :: prep_matchfilt
    procedure, private :: calc_corrs_over_k_1
    procedure, private :: calc_corrs_over_k_2
    generic            :: calc_corrs_over_k => calc_corrs_over_k_1, calc_corrs_over_k_2
    procedure, private :: calc_k_corrs_1
    procedure, private :: calc_k_corrs_2
    generic            :: calc_k_corrs => calc_k_corrs_1, calc_k_corrs_2
    procedure, private :: calc_corr_for_rot
    procedure, private :: calc_corr_for_rot_8
    procedure, private :: calc_T1_T2_for_rot_8
    procedure, private :: calc_euclid_for_rot
    procedure, private :: calc_euclid_for_rot_8
    procedure, private :: calc_corrk_for_rot
    procedure, private :: calc_corrk_for_rot_8
    procedure, private :: calc_euclidk_for_rot
    procedure, private :: gencorrs_cc_1
    procedure, private :: gencorrs_cc_2
    procedure, private :: gencorrs_cc_mir
    procedure, private :: gencorrs_euclid_1
    procedure, private :: gencorrs_euclid_2
    procedure, private :: gencorrs_euclid_mir
    procedure, private :: gencorrs_1
    procedure, private :: gencorrs_2
    generic            :: gencorrs => gencorrs_1, gencorrs_2
    procedure          :: gencorrs_mir
    procedure          :: gencorr_for_rot_8
    procedure          :: gencorr_grad_for_rot_8
    procedure          :: gencorr_grad_only_for_rot_8
    procedure          :: gencorr_cc_for_rot
    procedure          :: gencorr_cc_for_rot_8
    procedure          :: gencorr_cont_grad_cc_for_rot_8
    procedure          :: gencorr_cont_shift_grad_cc_for_rot_8
    procedure          :: gencorr_cc_grad_for_rot_8
    procedure          :: gencorr_cc_grad_only_for_rot_8
    procedure          :: gencorr_euclid_for_rot
    procedure          :: gencorr_euclid_for_rot_8
    procedure          :: gencorr_cont_grad_euclid_for_rot_8
    procedure          :: gencorr_cont_shift_grad_euclid_for_rot_8
    procedure          :: gencorr_euclid_grad_for_rot_8
    procedure          :: gencorr_euclid_grad_only_for_rot_8
    procedure          :: gencorr_sigma_contrib
    procedure, private :: genfrc
    procedure, private :: calc_frc
    procedure, private :: specscore_1
    procedure, private :: specscore_2
    generic            :: specscore => specscore_1, specscore_2
    procedure          :: calc_roinv_corrmat
    ! DESTRUCTOR
    procedure          :: kill
end type polarft_corrcalc

! CLASS PARAMETERS/VARIABLES
complex(sp), parameter           :: zero            = cmplx(0.,0.) !< just a complex zero
integer,     parameter           :: FFTW_USE_WISDOM = 16
class(polarft_corrcalc), pointer :: pftcc_glob

contains

    ! CONSTRUCTORS

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
            write(logfhandle,*) 'pfromto: ', self%pfromto(1), self%pfromto(2)
            THROW_HARD ('nptcls (# of particles) must be > 0; new')
        endif
        if( nrefs < 1 )then
            write(logfhandle,*) 'nrefs: ', nrefs
            THROW_HARD ('nrefs (# of reference sections) must be > 0; new')
        endif
        self%ldim = [params_glob%boxmatch,params_glob%boxmatch,1] !< logical dimensions of original cartesian image
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
            self%nptcls  = count(ptcl_mask)                            !< the total number of particles in partition
        else
            self%nptcls  = self%pfromto(2) - self%pfromto(1) + 1       !< the total number of particles in partition
        endif
        self%nrefs       = nrefs                                       !< the number of references (logically indexded [1,nrefs])
        self%nrots       = round2even(twopi * real(params_glob%ring2)) !< number of in-plane rotations for one pft  (determined by radius of molecule)
        self%pftsz       = self%nrots / 2                              !< size of reference (nrots/2) (number of vectors used for matching)
        ! allocate optimal low-pass filter if matched filter is on
        if( params_glob%l_match_filt )then
            if( params_glob%l_pssnr )then
                allocate(self%pssnr_filt(params_glob%kfromto(1):params_glob%kstop,2),source=1.)
            else
                allocate(self%ref_optlp(params_glob%kfromto(1):params_glob%kstop,self%nrefs),source=1.)
            endif
        endif
        ! generate polar coordinates & eo assignment
        allocate( self%polar(2*self%nrots,params_glob%kfromto(1):params_glob%kfromto(2)),&
                 &self%angtab(self%nrots), self%iseven(1:self%nptcls), stat=alloc_stat)
        if(alloc_stat/=0)call allocchk('polar coordinate arrays; new; simple_polarft_corrcalc, 1')
        ang = twopi/real(self%nrots)
        do irot=1,self%nrots
            self%angtab(irot) = real(irot-1)*ang
            ! cycling over non-redundant logical dimensions
            do k=params_glob%kfromto(1),params_glob%kfromto(2)
                self%polar(irot,k)            =  sin(self%angtab(irot))*real(k) ! x-coordinate
                self%polar(irot+self%nrots,k) = -cos(self%angtab(irot))*real(k) ! y-coordinate
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
                 &self%fftdat_ptcls(1:self%nptcls,params_glob%kfromto(1):params_glob%kfromto(2)),&
                 &self%heap_vars(params_glob%nthr), stat=alloc_stat)
        if(alloc_stat.ne.0)call allocchk('shared arrays; new; simple_polarft_corrcalc')
        local_stat=0
        do ithr=1,params_glob%nthr
            allocate(self%heap_vars(ithr)%pft_ref(self%pftsz,params_glob%kfromto(1):params_glob%kfromto(2)),&
                &self%heap_vars(ithr)%pft_ref_tmp(self%pftsz,params_glob%kfromto(1):params_glob%kfromto(2)),&
                &self%heap_vars(ithr)%pft_ref_tmp1(self%pftsz,params_glob%kfromto(1):params_glob%kfromto(2)),&
                &self%heap_vars(ithr)%pft_ref_tmp2(self%pftsz,params_glob%kfromto(1):params_glob%kfromto(2)),&
                &self%heap_vars(ithr)%pft_dref(self%pftsz,params_glob%kfromto(1):params_glob%kfromto(2),3),&
                &self%heap_vars(ithr)%corrs_over_k(self%nrots),&
                &self%heap_vars(ithr)%corrs_mir_over_k(self%nrots),&
                &self%heap_vars(ithr)%argmat(self%pftsz,params_glob%kfromto(1):params_glob%kfromto(2)),&
                &self%heap_vars(ithr)%shmat(self%pftsz,params_glob%kfromto(1):params_glob%kfromto(2)),&
                &self%heap_vars(ithr)%kcorrs(self%nrots),&
                &self%heap_vars(ithr)%kcorrs_mir(self%nrots),&
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
            self%fftdat(ithr)%p_ref_re       = fftwf_alloc_real   (int(self%pftsz, c_size_t))
            self%fftdat(ithr)%p_ref_im       = fftwf_alloc_complex(int(self%pftsz, c_size_t))
            self%fftdat(ithr)%p_ref_fft_re   = fftwf_alloc_complex(int(self%pftsz, c_size_t))
            self%fftdat(ithr)%p_ref_fft_im   = fftwf_alloc_complex(int(self%pftsz, c_size_t))
            self%fftdat(ithr)%p_ref_fft_re_2 = fftwf_alloc_complex(int(self%pftsz, c_size_t))
            self%fftdat(ithr)%p_ref_fft_im_2 = fftwf_alloc_complex(int(self%pftsz, c_size_t))
            self%fftdat(ithr)%p_product_fft  = fftwf_alloc_complex(int(self%nrots, c_size_t))
            self%fftdat(ithr)%p_backtransf   = fftwf_alloc_real   (int(self%nrots, c_size_t))
            call c_f_pointer(self%fftdat(ithr)%p_ref_re,        self%fftdat(ithr)%ref_re,        [self%pftsz])
            call c_f_pointer(self%fftdat(ithr)%p_ref_im,        self%fftdat(ithr)%ref_im,        [self%pftsz])
            call c_f_pointer(self%fftdat(ithr)%p_ref_fft_re,    self%fftdat(ithr)%ref_fft_re,    [self%pftsz])
            call c_f_pointer(self%fftdat(ithr)%p_ref_fft_im,    self%fftdat(ithr)%ref_fft_im,    [self%pftsz])
            call c_f_pointer(self%fftdat(ithr)%p_ref_fft_re_2,  self%fftdat(ithr)%ref_fft_re_2,  [self%pftsz])
            call c_f_pointer(self%fftdat(ithr)%p_ref_fft_im_2,  self%fftdat(ithr)%ref_fft_im_2,  [self%pftsz])
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
            do ik = params_glob%kfromto(1),params_glob%kfromto(2)
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
        ! 2Dclass/3Drefs mapping on/off
        self%l_clsfrcs = params_glob%clsfrcs.eq.'yes'
        ! flag existence
        self%existence = .true.
        ! set pointer to global instance
        pftcc_glob => self
    end subroutine new

    ! SETTERS

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

    subroutine assign_sigma2_noise( self, sigma2_noise )
        class(polarft_corrcalc),      intent(inout) :: self
        real,    allocatable, target, intent(inout) :: sigma2_noise(:,:)
        self%sigma2_noise      => sigma2_noise
    end subroutine assign_sigma2_noise

    subroutine set_ref_optlp( self, iref, optlp )
        class(polarft_corrcalc), intent(inout) :: self
        integer,                 intent(in)    :: iref
        real,                    intent(in)    :: optlp(params_glob%kfromto(1):params_glob%kstop)
        self%ref_optlp(:,iref) = optlp(:)
    end subroutine set_ref_optlp

    subroutine set_pssnr_filt( self, filt, iseven )
        class(polarft_corrcalc), intent(inout) :: self
        real,                    intent(in)    :: filt(params_glob%kfromto(1):params_glob%kstop)
        logical,                 intent(in)    :: iseven
        if( iseven )then
            ! even is 1
            self%pssnr_filt(:,1) = filt(:)
        else
            self%pssnr_filt(:,2) = filt(:)
        endif
    end subroutine set_pssnr_filt

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
        pdim = [self%pftsz,params_glob%kfromto(1),params_glob%kfromto(2)]
    end function get_pdim

    ! !>  \brief  for getting the dimension of the reference polar FT
    pure integer function get_pftsz( self )
        class(polarft_corrcalc), intent(in) :: self
        get_pftsz = self%pftsz
    end function get_pftsz

    pure function get_box( self ) result( box )
        class(polarft_corrcalc), intent(in) :: self
        integer :: box
        box = self%ldim(1)
    end function get_box

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

    !>  \brief  returns polar Fourier transform of particle iptcl
    function get_ptcl_pft( self, iptcl) result( pft )
        class(polarft_corrcalc), intent(in) :: self
        integer,                 intent(in) :: iptcl
        complex(sp), allocatable :: pft(:,:)
        allocate(pft(self%pftsz,params_glob%kfromto(1):params_glob%kfromto(2)),&
        source=self%pfts_ptcls(:,:,self%pinds(iptcl)), stat=alloc_stat)
        if(alloc_stat.ne.0)call allocchk("In: get_ptcl_pft; simple_polarft_corrcalc",alloc_stat)
    end function get_ptcl_pft

    !>  \brief  returns polar Fourier transform of reference iref
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

    subroutine print( self )
        class(polarft_corrcalc), intent(in) :: self
        write(logfhandle,*) "total n particles in partition         (self%nptcls): ", self%nptcls
        write(logfhandle,*) "number of references                    (self%nrefs): ", self%nrefs
        write(logfhandle,*) "number of rotations                     (self%nrots): ", self%nrots
        write(logfhandle,*) "size of pft                             (self%pftsz): ", self%pftsz
        write(logfhandle,*) "logical dim. of original Cartesian image (self%ldim): ", self%ldim
    end subroutine print

    ! MODIFIERS

    subroutine shellnorm_and_filter_ref( self, iptcl, iref, pft )
        class(polarft_corrcalc), intent(in)    :: self
        integer,                 intent(in)    :: iptcl, iref
        complex(sp),             intent(inout) :: pft(self%pftsz,params_glob%kfromto(1):params_glob%kfromto(2))
        real    :: pw
        integer :: k, j
        if( params_glob%l_pssnr )then
            j = merge(1, 2, self%iseven(self%pinds(iptcl)))
            do k=params_glob%kfromto(1),params_glob%kstop
                pw = sum(csq(pft(:,k))) / real(self%pftsz)
                if( pw > 0.000001 )then
                    pft(:,k) = pft(:,k) * sqrt(self%pssnr_filt(k,j) / pw)
                else
                    pft(:,k) = pft(:,k) * sqrt(self%pssnr_filt(k,j))
                endif
            enddo
        else
            do k=params_glob%kfromto(1),params_glob%kstop
                pw       = sum(csq(pft(:,k))) / real(self%pftsz)
                if( pw > 0.000001 )then
                    pft(:,k) = pft(:,k) * (self%ref_optlp(k,iref) / sqrt(pw))
                else
                    pft(:,k) = pft(:,k) * self%ref_optlp(k,iref)
                endif
            enddo
        endif
    end subroutine shellnorm_and_filter_ref

    subroutine shellnorm_and_filter_ref_8( self, iptcl, iref, pft )
        class(polarft_corrcalc), intent(in)    :: self
        integer,                 intent(in)    :: iptcl, iref
        complex(dp),             intent(inout) :: pft(self%pftsz,params_glob%kfromto(1):params_glob%kfromto(2))
        real(dp) :: pw
        integer  :: j, k
        if( params_glob%l_pssnr )then
            j = merge(1, 2, self%iseven(self%pinds(iptcl)))
            do k=params_glob%kfromto(1),params_glob%kstop
                pw = sum(csq(pft(:,k))) / real(self%pftsz,kind=dp)
                if( pw > 0.000001d0 )then
                    pft(:,k) = pft(:,k) * dsqrt(real(self%pssnr_filt(k,j),kind=dp) / pw)
                else
                    pft(:,k) = pft(:,k) * dsqrt(real(self%pssnr_filt(k,j),kind=dp))
                endif
            enddo
        else
            do k=params_glob%kfromto(1),params_glob%kstop
                pw = sum(csq(pft(:,k))) / real(self%pftsz,kind=dp)
                if( pw > 0.000001d0 )then
                    pft(:,k) = pft(:,k) * (real(self%ref_optlp(k,iref),kind=dp) / dsqrt(pw))
                else
                    pft(:,k) = pft(:,k) * real(self%ref_optlp(k,iref),kind=dp)
                endif
            enddo
        endif
    end subroutine shellnorm_and_filter_ref_8

    subroutine shellnorm_and_filter_ref_dref_8( self, iptcl, iref, pft, dpft )
        class(polarft_corrcalc), intent(in)    :: self
        integer,                 intent(in)    :: iptcl, iref
        complex(dp),             intent(inout) :: pft(self%pftsz,params_glob%kfromto(1):params_glob%kfromto(2))
        complex(dp),             intent(inout) :: dpft(self%pftsz,params_glob%kfromto(1):params_glob%kfromto(2),3)
        real(dp) :: w, pw
        integer  :: k, j
        if( params_glob%l_pssnr )then
            j = merge(1, 2, self%iseven(self%pinds(iptcl)))
            do k=params_glob%kfromto(1),params_glob%kstop
                pw = sum(csq(pft(:,k))) / real(self%pftsz,kind=dp)
                if( pw > 0.000001d0 )then
                    w  = dsqrt( real(self%pssnr_filt(k,j),kind=dp) / pw)
                else
                    w  = dsqrt(real(self%pssnr_filt(k,j),kind=dp))
                endif
                pft(:,k)    = w * pft(:,k)
                dpft(:,k,:) = w * dpft(:,k,:)
            enddo
        else
            do k=params_glob%kfromto(1),params_glob%kstop
                pw = sum(csq(pft(:,k))) / real(self%pftsz,kind=dp)
                if( pw > 0.000001d0 )then
                    w  = real(self%ref_optlp(k,iref),kind=dp) / dsqrt(pw)
                else
                    w  = real(self%ref_optlp(k,iref),kind=dp)
                endif
                pft(:,k)    = w * pft(:,k)
                dpft(:,k,:) = w * dpft(:,k,:)
            enddo
        endif
    end subroutine shellnorm_and_filter_ref_dref_8

    ! MEMOIZERS

    subroutine memoize_sqsum_ptcl( self, i )
        class(polarft_corrcalc), intent(inout) :: self
        integer,                 intent(in)    :: i
        self%sqsums_ptcls(i) = sum(csq(self%pfts_ptcls(:,params_glob%kfromto(1):params_glob%kstop,i)))
    end subroutine memoize_sqsum_ptcl

    ! memoize all particles ffts
    subroutine memoize_ffts( self )
        class(polarft_corrcalc), intent(inout) :: self
        integer :: i
        if( params_glob%l_match_filt )then
            ! because with match_filt=yes the particle filtering is done on the fly
            return
        endif
        ! memoize particle FFTs in parallel
        !$omp parallel do default(shared) private(i) proc_bind(close) schedule(static)
        do i=1,self%nptcls
            call self%memoize_fft(i)
        end do
        !$omp end parallel do
    end subroutine memoize_ffts

    ! memoize particle fft, serial only
    subroutine memoize_fft( self, i )
        class(polarft_corrcalc), intent(inout) :: self
        integer  :: i, ik, ithr
        ithr = omp_get_thread_num() + 1
        ! memoize particle FFTs
        do ik = params_glob%kfromto(1),params_glob%kfromto(2)
            ! copy particle pfts
            self%fft_carray(ithr)%re = real(self%pfts_ptcls(:,ik,i))
            self%fft_carray(ithr)%im = aimag(self%pfts_ptcls(:,ik,i)) * self%fft_factors
            ! FFT
            call fftwf_execute_dft_r2c(self%plan_fwd_1, self%fft_carray(ithr)%re, self%fftdat_ptcls(i,ik)%re)
            call fftwf_execute_dft    (self%plan_fwd_2, self%fft_carray(ithr)%im, self%fftdat_ptcls(i,ik)%im)
        end do
    end subroutine memoize_fft

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

    subroutine prep_ref4corr( self, iref, iptcl, pft_ref, sqsum_ref, kstop )
        class(polarft_corrcalc), intent(inout) :: self
        integer,                 intent(in)    :: iref, iptcl
        integer,                 intent(in)    :: kstop
        complex(sp),             intent(out)   :: pft_ref(self%pftsz,params_glob%kfromto(1):params_glob%kfromto(2))
        real(sp),                intent(out)   :: sqsum_ref
        integer :: i
        i = self%pinds(iptcl)
        ! copy
        if( self%iseven(i) )then
            pft_ref = self%pfts_refs_even(:,:,iref)
        else
            pft_ref = self%pfts_refs_odd(:,:,iref)
        endif
        ! shell normalization and filtering
        if( params_glob%l_match_filt )then
            if( self%l_clsfrcs )then
                call self%shellnorm_and_filter_ref(iptcl, iptcl, pft_ref)
            else
                call self%shellnorm_and_filter_ref(iptcl, iref, pft_ref)
            endif
        endif
        ! multiply with CTF
        if( self%with_ctf ) pft_ref = pft_ref * self%ctfmats(:,:,i)
        ! for corr normalisation
        sqsum_ref = sum(csq(pft_ref(:,params_glob%kfromto(1):kstop)))
    end subroutine prep_ref4corr

    subroutine prep_matchfilt( self, iptcl, iref, irot)
        class(polarft_corrcalc), intent(inout) :: self
        integer,                 intent(in)    :: iptcl, iref, irot
        complex(sp), pointer :: pft_ref(:,:)
        complex(sp) :: pft_ptcl(self%pftsz,params_glob%kfromto(1):params_glob%kfromto(2))
        real        :: pw_diff(params_glob%kfromto(1):params_glob%kfromto(2))
        real        :: pw_diff_fit(params_glob%kfromto(1):params_glob%kfromto(2))
        real        :: pw_ptcl, pw_ref, w, ssnr
        integer     :: i, j, k, ithr, rot
        ! particle is assumed phase-flipped, reference untouched
        i = self%pinds(iptcl)
        if( iref == 0 .or. .not. params_glob%l_match_filt )then
            ! just memoize particle
            call self%memoize_fft(i)
            return
        endif
        ! init
        ithr     =  omp_get_thread_num() + 1
        pft_ref  => self%heap_vars(ithr)%pft_ref
        if( self%iseven(i) )then
            pft_ref = self%pfts_refs_even(:,:,iref)
        else
            pft_ref = self%pfts_refs_odd(:,:,iref)
        endif
        pft_ptcl = self%pfts_ptcls(:,:,i)
        if( params_glob%l_pssnr )then
            j = merge(1, 2, self%iseven(i))
            do k=params_glob%kfromto(1),params_glob%kstop
                ! particle power spectrum
                pw_ptcl = sum(csq(pft_ptcl(:,k))) / real(self%pftsz)
                ! shell normalization
                pft_ptcl(:,k) = pft_ptcl(:,k) / sqrt(pw_ptcl)
                ! pssnr filter
                pft_ptcl(:,k) = pft_ptcl(:,k) * sqrt(1. + self%pssnr_filt(k,j) * self%ctfmats(:,k,i)**2.)
            enddo
        else
            ! CTF
            if( self%with_ctf )then
                ! particle is phase-flipped
                ! reference: x|CTF|
                pft_ref = pft_ref * self%ctfmats(:,:,i)
            endif
            ! match spectrums
            do k=params_glob%kfromto(1),params_glob%kstop
                ! particle
                pw_ptcl = sum(csq(pft_ptcl(:,k))) / real(self%pftsz)
                if( pw_ptcl > 0.000001 )then
                    pft_ptcl(:,k) = pft_ptcl(:,k) * self%ref_optlp(k,iref) / sqrt(pw_ptcl)
                else
                    pft_ptcl(:,k) = pft_ptcl(:,k) * self%ref_optlp(k,iref)
                endif
                ! reference
                pw_ref = sum(csq(pft_ref(:,k))) / real(self%pftsz)
                if( pw_ref > 0.000001 )then
                    pft_ref(:,k) = pft_ref(:,k) * self%ref_optlp(k,iref) / sqrt(pw_ref)
                else
                    pft_ref(:,k) = pft_ref(:,k) * self%ref_optlp(k,iref)
                endif
            enddo
            ! power spectrum difference
            rot = merge(irot - self%pftsz, irot, irot >= self%pftsz + 1)
            do k = params_glob%kfromto(1), params_glob%kstop
                if( irot == 1 )then
                    pw_diff(k) = sum(csq(pft_ref(:,k) - pft_ptcl(:,k)))
                else if( irot <= self%pftsz )then
                    pw_diff(k) =              sum(csq(pft_ref(1:self%pftsz-rot+1,k)          - pft_ptcl(rot:self%pftsz,k)))
                    pw_diff(k) = pw_diff(k) + sum(csq(pft_ref(self%pftsz-rot+2:self%pftsz,k) - conjg(pft_ptcl(1:rot-1,k))))
                else if( irot == self%pftsz + 1 )then
                    pw_diff(k) = sum(csq(pft_ref(:,k) - conjg(pft_ptcl(:,k))))
                else
                    pw_diff(k) =              sum(csq(pft_ref(1:self%pftsz-rot+1,k)          - conjg(pft_ptcl(rot:self%pftsz,k))))
                    pw_diff(k) = pw_diff(k) + sum(csq(pft_ref(self%pftsz-rot+2:self%pftsz,k) - pft_ptcl(1:rot-1,k)))
                end if
                pw_diff(k) = pw_diff(k) / real(self%pftsz)
            end do
            ! fitting of difference
            pw_diff_fit = pw_diff
            call SavitzkyGolay_filter(params_glob%kstop-params_glob%kfromto(1)+1, pw_diff_fit)
            where( pw_diff_fit <= 0.) pw_diff_fit = pw_diff                 ! takes care of range
            pw_diff_fit(params_glob%kfromto(1):params_glob%kfromto(1)+2) =& ! to avoid left border effect
                &pw_diff(params_glob%kfromto(1):params_glob%kfromto(1)+2)
            ! shell normalize & filter particle
            do k=params_glob%kfromto(1),params_glob%kstop
                ! particle
                pw_ptcl = sum(csq(pft_ptcl(:,k))) / real(self%pftsz)
                w = 1.
                if( pw_ptcl > 0.000001 ) w  = 1. / sqrt(pw_ptcl)
                ! reference
                pw_ref = sum(csq(pft_ref(:,k))) / real(self%pftsz)
                ! optimal filter
                if( pw_diff_fit(k) > 0.000001 )then
                    ssnr = pw_ref / pw_diff_fit(k)
                    w    = w * sqrt(ssnr / (ssnr + 1.))
                else
                    ! ssnr -> inf
                endif
                ! filter
                pft_ptcl(:,k) = w * pft_ptcl(:,k)
            enddo
        endif
        ! update particle pft & re-memoize
        call self%set_ptcl_pft(iptcl, pft_ptcl)
        call self%memoize_fft(i)
    end subroutine prep_matchfilt

    subroutine calc_corrs_over_k_1( self, pft_ref, i, kfromto, corrs_over_k)
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
    end subroutine calc_corrs_over_k_1

    subroutine calc_corrs_over_k_2( self, pft_ref, i, kfromto, corrs_over_k, corrs_mir_over_k )
        class(polarft_corrcalc), intent(inout) :: self
        integer,                 intent(in)    :: i, kfromto(2)
        complex(sp),             intent(in)    :: pft_ref(1:self%pftsz,params_glob%kfromto(1):params_glob%kfromto(2))
        real,                    intent(out)   :: corrs_over_k(self%nrots)
        real,                    intent(out)   :: corrs_mir_over_k(self%nrots)
        integer :: ithr, ik
        ! get thread index
        ithr = omp_get_thread_num() + 1
        ! sum up correlations over k-rings
        corrs_over_k     = 0.
        corrs_mir_over_k = 0.
        do ik = kfromto(1),kfromto(2)
            ! move reference into Fourier Fourier space (particles are memoized)
            self%fftdat(ithr)%ref_re(:) =  real(pft_ref(:,ik))
            self%fftdat(ithr)%ref_im(:) = aimag(pft_ref(:,ik)) * self%fft_factors
            call fftwf_execute_dft_r2c(self%plan_fwd_1, self%fftdat(ithr)%ref_re, self%fftdat(ithr)%ref_fft_re)
            call fftwf_execute_dft    (self%plan_fwd_2, self%fftdat(ithr)%ref_im, self%fftdat(ithr)%ref_fft_im)
            ! correlate
            self%fftdat(ithr)%ref_fft_re_2 =       self%fftdat(ithr)%ref_fft_re  * self%fftdat_ptcls(i,ik)%re
            self%fftdat(ithr)%ref_fft_im_2 =       self%fftdat(ithr)%ref_fft_im  * self%fftdat_ptcls(i,ik)%im
            self%fftdat(ithr)%ref_fft_re   = conjg(self%fftdat(ithr)%ref_fft_re) * self%fftdat_ptcls(i,ik)%re
            self%fftdat(ithr)%ref_fft_im   = conjg(self%fftdat(ithr)%ref_fft_im) * self%fftdat_ptcls(i,ik)%im
            self%fftdat(ithr)%product_fft(1:1 + 2 * int(self%pftsz / 2):2) = &
                4. * self%fftdat(ithr)%ref_fft_re(1:1+int(self%pftsz/2))
            self%fftdat(ithr)%product_fft(2:2 + 2 * int(self%pftsz / 2):2) = &
                4. * self%fftdat(ithr)%ref_fft_im(1:int(self%pftsz / 2) + 1)
            ! back transform
            call fftwf_execute_dft_c2r(self%plan_bwd, self%fftdat(ithr)%product_fft, self%fftdat(ithr)%backtransf)
            ! accumulate corrs
            corrs_over_k = corrs_over_k + self%fftdat(ithr)%backtransf
            ! TODO: use fftw_plan_many_dft_c2r interface
            ! correlate
            self%fftdat(ithr)%product_fft(1:1 + 2 * int(self%pftsz / 2):2) = &
                4. * self%fftdat(ithr)%ref_fft_re(1:1+int(self%pftsz/2))
            self%fftdat(ithr)%product_fft(2:2 + 2 * int(self%pftsz / 2):2) = &
                4. * self%fftdat(ithr)%ref_fft_im(1:int(self%pftsz / 2) + 1)
            ! back transform
            call fftwf_execute_dft_c2r(self%plan_bwd, self%fftdat(ithr)%product_fft, self%fftdat(ithr)%backtransf)
            ! accumulate corrs
            corrs_mir_over_k = corrs_mir_over_k + self%fftdat(ithr)%backtransf
        end do
        ! fftw3 routines are not properly normalized, hence division by self%nrots * 2
        corrs_over_k     = corrs_over_k     / real(self%nrots * 2)
        corrs_mir_over_k = corrs_mir_over_k / real(self%nrots * 2)
    end subroutine calc_corrs_over_k_2

    subroutine calc_k_corrs_1( self, pft_ref, i, k, kcorrs )
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
    end subroutine calc_k_corrs_1

    subroutine calc_k_corrs_2( self, pft_ref, i, k, kcorrs, kcorrs_mir )
        class(polarft_corrcalc), intent(inout) :: self
        complex(sp),             intent(in)    :: pft_ref(1:self%pftsz,params_glob%kfromto(1):params_glob%kfromto(2))
        integer,                 intent(in)    :: i, k
        real,                    intent(out)   :: kcorrs(self%nrots)
        real,                    intent(out)   :: kcorrs_mir(self%nrots)
        integer :: ithr
        ! get thread index
        ithr = omp_get_thread_num() + 1
        ! move reference into Fourier Fourier space (particles are memoized)
        self%fftdat(ithr)%ref_re(:) =  real(pft_ref(:,k))
        self%fftdat(ithr)%ref_im(:) = aimag(pft_ref(:,k)) * self%fft_factors
        call fftwf_execute_dft_r2c(self%plan_fwd_1, self%fftdat(ithr)%ref_re, self%fftdat(ithr)%ref_fft_re)
        call fftwf_execute_dft    (self%plan_fwd_2, self%fftdat(ithr)%ref_im, self%fftdat(ithr)%ref_fft_im)
        ! correlate FFTs
        self%fftdat(ithr)%ref_fft_re_2 =       self%fftdat(ithr)%ref_fft_re  * self%fftdat_ptcls(i,k)%re
        self%fftdat(ithr)%ref_fft_im_2 =       self%fftdat(ithr)%ref_fft_im  * self%fftdat_ptcls(i,k)%im
        self%fftdat(ithr)%ref_fft_re   = conjg(self%fftdat(ithr)%ref_fft_re) * self%fftdat_ptcls(i,k)%re
        self%fftdat(ithr)%ref_fft_im   = conjg(self%fftdat(ithr)%ref_fft_im) * self%fftdat_ptcls(i,k)%im
        self%fftdat(ithr)%product_fft(1:1 + 2 * int(self%pftsz / 2):2) = &
            4. * self%fftdat(ithr)%ref_fft_re(1:1 + int(self%pftsz / 2))
        self%fftdat(ithr)%product_fft(2:2 + 2 * int(self%pftsz / 2):2) = &
            4. * self%fftdat(ithr)%ref_fft_im(1:int(self%pftsz / 2) + 1)
        ! back transform
        call fftwf_execute_dft_c2r(self%plan_bwd, self%fftdat(ithr)%product_fft, self%fftdat(ithr)%backtransf)
        ! fftw3 routines are not properly normalized, hence division by self%nrots * 2
        kcorrs = self%fftdat(ithr)%backtransf / real(self%nrots * 2)
        ! correlate FFTs
        self%fftdat(ithr)%product_fft(1:1 + 2 * int(self%pftsz / 2):2) = &
            4. * self%fftdat(ithr)%ref_fft_re_2(1:1 + int(self%pftsz / 2))
        self%fftdat(ithr)%product_fft(2:2 + 2 * int(self%pftsz / 2):2) = &
            4. * self%fftdat(ithr)%ref_fft_im_2(1:int(self%pftsz / 2) + 1)
        ! back transform
        call fftwf_execute_dft_c2r(self%plan_bwd, self%fftdat(ithr)%product_fft, self%fftdat(ithr)%backtransf)
        ! fftw3 routines are not properly normalized, hence division by self%nrots * 2
        kcorrs_mir = self%fftdat(ithr)%backtransf / real(self%nrots * 2)
        ! TODO: use fftw_plan_many_dft_c2r interface
    end subroutine calc_k_corrs_2

    function calc_corr_for_rot( self, pft_ref, i, irot )result( corr )
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
            tmp = sum( pft_ref(:,params_glob%kfromto(1):params_glob%kstop) * conjg(self%pfts_ptcls(:,params_glob%kfromto(1):params_glob%kstop,i)))
        else if( irot <= self%pftsz )then
            tmp =       sum( pft_ref(               1:self%pftsz-rot+1,params_glob%kfromto(1):params_glob%kstop) * conjg(self%pfts_ptcls(rot:self%pftsz,params_glob%kfromto(1):params_glob%kstop,i)))
            tmp = tmp + sum( pft_ref(self%pftsz-rot+2:self%pftsz,      params_glob%kfromto(1):params_glob%kstop) *       self%pfts_ptcls(  1:rot-1,     params_glob%kfromto(1):params_glob%kstop,i))
        else if( irot == self%pftsz + 1 )then
            tmp = sum( pft_ref(:,params_glob%kfromto(1):params_glob%kstop) * self%pfts_ptcls(:,params_glob%kfromto(1):params_glob%kstop,i) )
        else
            tmp =       sum( pft_ref(1:self%pftsz-rot+1,params_glob%kfromto(1):params_glob%kstop)          *        self%pfts_ptcls(rot:self%pftsz,params_glob%kfromto(1):params_glob%kstop,i))
            tmp = tmp + sum( pft_ref(self%pftsz-rot+2:self%pftsz,params_glob%kfromto(1):params_glob%kstop) * conjg( self%pfts_ptcls( 1:rot-1,      params_glob%kfromto(1):params_glob%kstop,i)))
        end if
        corr = real(tmp)
    end function calc_corr_for_rot

    function calc_corr_for_rot_8( self, pft_ref, i, irot )result( corr )
        class(polarft_corrcalc), intent(inout) :: self
        integer,                 intent(in)    :: i, irot
        complex(dp),             intent(in)    :: pft_ref(1:self%pftsz,params_glob%kfromto(1):params_glob%kfromto(2))
        integer     :: rot
        real(dp)    :: corr
        complex(dp) :: tmp
        tmp = 0.
        if( irot >= self%pftsz + 1 )then
            rot = irot - self%pftsz
        else
            rot = irot
        end if
        if( irot == 1 )then
            tmp = sum( pft_ref(:,params_glob%kfromto(1):params_glob%kstop) * conjg(self%pfts_ptcls(:,params_glob%kfromto(1):params_glob%kstop,i)))
        else if( irot <= self%pftsz )then
            tmp =       sum( pft_ref(               1:self%pftsz-rot+1,params_glob%kfromto(1):params_glob%kstop) * conjg(self%pfts_ptcls(rot:self%pftsz,params_glob%kfromto(1):params_glob%kstop,i)))
            tmp = tmp + sum( pft_ref(self%pftsz-rot+2:self%pftsz,      params_glob%kfromto(1):params_glob%kstop) *       self%pfts_ptcls(  1:rot-1,     params_glob%kfromto(1):params_glob%kstop,i))
        else if( irot == self%pftsz + 1 ) then
            tmp = sum( pft_ref(:,params_glob%kfromto(1):params_glob%kstop) * self%pfts_ptcls(:,params_glob%kfromto(1):params_glob%kstop,i) )
        else
            tmp =       sum( pft_ref(               1:self%pftsz-rot+1,params_glob%kfromto(1):params_glob%kstop) *        self%pfts_ptcls(rot:self%pftsz,params_glob%kfromto(1):params_glob%kstop,i))
            tmp = tmp + sum( pft_ref(self%pftsz-rot+2:self%pftsz,      params_glob%kfromto(1):params_glob%kstop) * conjg( self%pfts_ptcls( 1:rot-1,      params_glob%kfromto(1):params_glob%kstop,i)))
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
        real(sp)    :: euclid, tmp
        if( irot >= self%pftsz + 1 )then
            rot = irot - self%pftsz
        else
            rot = irot
        end if
        euclid = 0.
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
        integer  :: rot, k
        real(dp) :: euclid, tmp
        if( irot >= self%pftsz + 1 )then
            rot = irot - self%pftsz
        else
            rot = irot
        end if
        euclid = 0.d0
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
            euclid = euclid - tmp / ( 2.d0 * self%sigma2_noise(k, i))
        end do
    end function calc_euclid_for_rot_8

    function calc_corrk_for_rot( self, pft_ref, i, k, irot ) result( corr )
        class(polarft_corrcalc), intent(inout) :: self
        integer,                 intent(in)    :: i, k, irot
        complex,                 intent(in)    :: pft_ref(1:self%pftsz,params_glob%kfromto(1):params_glob%kfromto(2))
        integer :: rot
        real    :: corr
        complex :: tmp
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
            tmp =       sum( pft_ref(1:self%pftsz-rot+1,k)          *        self%pfts_ptcls(rot:self%pftsz,k,i))
            tmp = tmp + sum( pft_ref(self%pftsz-rot+2:self%pftsz,k) * conjg( self%pfts_ptcls(1:rot-1,k,i) ))
        end if
        corr = real(tmp)
    end function calc_corrk_for_rot

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

    real(sp) function calc_euclidk_for_rot( self, pft_ref, i, k, irot )
        class(polarft_corrcalc), intent(inout) :: self
        integer,                 intent(in)    :: i, irot, k
        complex(sp),             intent(in)    :: pft_ref(1:self%pftsz,params_glob%kfromto(1):params_glob%kfromto(2))
        real    :: euclid
        integer :: rot
        if( irot >= self%pftsz + 1 )then
            rot = irot - self%pftsz
        else
            rot = irot
        end if
        if( irot == 1 )then
            euclid = sum(csq(pft_ref(:,k) - self%pfts_ptcls(:,k,i)))
        else if( irot <= self%pftsz )then
            euclid =          sum(csq(pft_ref(               1:self%pftsz-rot+1,k) -       self%pfts_ptcls(rot:self%pftsz,k,i)))
            euclid = euclid + sum(csq(pft_ref(self%pftsz-rot+2:self%pftsz,      k) - conjg(self%pfts_ptcls(  1:rot-1,     k,i))))
        else if( irot == self%pftsz + 1 )then
            euclid = sum(csq(pft_ref(:,k) - conjg(self%pfts_ptcls(:,k,i))))
        else
            euclid =          sum(csq(pft_ref(               1:self%pftsz-rot+1,k) - conjg(self%pfts_ptcls(rot:self%pftsz,k,i))))
            euclid = euclid + sum(csq(pft_ref(self%pftsz-rot+2:self%pftsz,      k) -       self%pfts_ptcls(  1:rot-1,     k,i)))
        end if
    end function calc_euclidk_for_rot

    subroutine genfrc( self, iref, iptcl, irot, frc )
        class(polarft_corrcalc),  intent(inout) :: self
        integer,                  intent(in)    :: iref, iptcl, irot
        real(sp),                 intent(out)   :: frc(params_glob%kfromto(1):params_glob%kstop)
        complex(sp), pointer :: pft_ref(:,:)
        real(sp),    pointer :: kcorrs(:)
        real(sp) :: sumsqref, sumsqptcl, sqsum_ref
        integer  :: k, ithr, i
        i       =  self%pinds(iptcl)
        ithr    =  omp_get_thread_num() + 1
        pft_ref => self%heap_vars(ithr)%pft_ref
        kcorrs  => self%heap_vars(ithr)%kcorrs
        call self%prep_ref4corr(iref, iptcl, pft_ref, sqsum_ref, params_glob%kstop)
        do k=params_glob%kfromto(1),params_glob%kstop
            call self%calc_k_corrs(pft_ref, i, k, kcorrs)
            sumsqptcl = sum(csq(self%pfts_ptcls(:,k,i)))
            sumsqref  = sum(csq(pft_ref(:,k)))
            frc(k)    = kcorrs(irot) / sqrt(sumsqptcl * sumsqref)
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
        i       =  self%pinds(iptcl)
        ithr    =  omp_get_thread_num() + 1
        pft_ref => self%heap_vars(ithr)%pft_ref
        kcorrs  => self%heap_vars(ithr)%kcorrs
        shmat   => self%heap_vars(ithr)%shmat
        argmat  => self%heap_vars(ithr)%argmat
        argmat  =  self%argtransf(:self%pftsz,:) * shvec(1) + self%argtransf(self%pftsz + 1:,:) * shvec(2)
        shmat   =  cmplx(cos(argmat),sin(argmat))
        if( self%iseven(i) )then
            pft_ref = self%pfts_refs_even(:,:,iref)
        else
            pft_ref = self%pfts_refs_odd(:,:,iref)
        endif
        if( params_glob%l_match_filt )then
            if( self%l_clsfrcs )then
                call self%shellnorm_and_filter_ref(iptcl, iptcl, pft_ref)
            else
                call self%shellnorm_and_filter_ref(iptcl, iref, pft_ref)
            endif
        endif
        if( self%with_ctf )then
            pft_ref = (pft_ref * self%ctfmats(:,:,i)) * shmat
        else
            pft_ref = pft_ref * shmat
        endif
        do k=params_glob%kfromto(1),params_glob%kstop
            call self%calc_k_corrs(pft_ref, i, k, kcorrs)
            sumsqptcl = sum(csq(self%pfts_ptcls(:,k,i)))
            sumsqref  = sum(csq(pft_ref(:,k)))
            frc(k)    = kcorrs(irot) / sqrt(sumsqptcl * sumsqref)
        end do
    end subroutine calc_frc

    subroutine gencorrs_cc_1( self, iref, iptcl, cc )
        class(polarft_corrcalc), intent(inout) :: self
        integer,                 intent(in)    :: iref, iptcl
        real,                    intent(out)   :: cc(self%nrots)
        complex(sp), pointer :: pft_ref(:,:)
        real(sp),    pointer :: corrs_over_k(:)
        real(sp) :: sqsum_ref
        integer  :: ithr
        ithr         =  omp_get_thread_num() + 1
        pft_ref      => self%heap_vars(ithr)%pft_ref
        corrs_over_k => self%heap_vars(ithr)%corrs_over_k
        call self%prep_ref4corr(iref, iptcl, pft_ref, sqsum_ref, params_glob%kstop)
        call self%calc_corrs_over_k(pft_ref, self%pinds(iptcl), [params_glob%kfromto(1),params_glob%kstop], corrs_over_k)
        cc = corrs_over_k / sqrt(sqsum_ref * self%sqsums_ptcls(self%pinds(iptcl)))
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
        ithr         =  omp_get_thread_num() + 1
        pft_ref      => self%heap_vars(ithr)%pft_ref
        shmat        => self%heap_vars(ithr)%shmat
        corrs_over_k => self%heap_vars(ithr)%corrs_over_k
        argmat       => self%heap_vars(ithr)%argmat
        argmat = self%argtransf(:self%pftsz,:) * shvec(1) + self%argtransf(self%pftsz + 1:,:) * shvec(2)
        shmat  = cmplx(cos(argmat),sin(argmat))
        if( self%iseven(self%pinds(iptcl)) )then
            pft_ref = self%pfts_refs_even(:,:,iref)
        else
            pft_ref = self%pfts_refs_odd(:,:,iref)
        endif
        if( params_glob%l_match_filt )then
            if( self%l_clsfrcs )then
                call self%shellnorm_and_filter_ref(iptcl, iptcl, pft_ref)
            else
                call self%shellnorm_and_filter_ref(iptcl, iref, pft_ref)
            endif
        endif
        if( self%with_ctf )then
            pft_ref = (pft_ref * self%ctfmats(:,:,self%pinds(iptcl))) * shmat
        else
            pft_ref = pft_ref * shmat
        endif
        sqsum_ref = sum(csq(pft_ref))
        call self%calc_corrs_over_k(pft_ref, self%pinds(iptcl), [params_glob%kfromto(1),params_glob%kstop], corrs_over_k)
        cc = corrs_over_k  / sqrt(sqsum_ref * self%sqsums_ptcls(self%pinds(iptcl)))
    end subroutine gencorrs_cc_2

    subroutine gencorrs_cc_mir( self, iref, iptcl, cc, cc_mir )
        class(polarft_corrcalc), intent(inout) :: self
        integer,                 intent(in)    :: iref, iptcl
        real,                    intent(out)   :: cc(self%nrots)
        real,                    intent(out)   :: cc_mir(self%nrots)
        complex(sp), pointer :: pft_ref(:,:)
        real(sp),    pointer :: corrs_over_k(:)
        real(sp),    pointer :: corrs_mir_over_k(:)
        real(sp) :: sqsum_ref
        integer  :: ithr
        ithr             =  omp_get_thread_num() + 1
        pft_ref          => self%heap_vars(ithr)%pft_ref
        corrs_over_k     => self%heap_vars(ithr)%corrs_over_k
        corrs_mir_over_k => self%heap_vars(ithr)%corrs_mir_over_k
        call self%prep_ref4corr(iref, iptcl, pft_ref, sqsum_ref, params_glob%kstop)
        call self%calc_corrs_over_k(pft_ref, self%pinds(iptcl), [params_glob%kfromto(1),params_glob%kstop], corrs_over_k, &
            corrs_mir_over_k)
        cc     = corrs_over_k     / sqrt(sqsum_ref * self%sqsums_ptcls(self%pinds(iptcl)))
        cc_mir = corrs_mir_over_k / sqrt(sqsum_ref * self%sqsums_ptcls(self%pinds(iptcl)))
    end subroutine gencorrs_cc_mir

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
        keuclids       => self%heap_vars(ithr)%kcorrs ! can be reused
        call self%prep_ref4corr(iref, iptcl, pft_ref, sumsqref, params_glob%kstop)
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
        ithr     =  omp_get_thread_num() + 1
        pft_ref  => self%heap_vars(ithr)%pft_ref
        shmat    => self%heap_vars(ithr)%shmat
        keuclids => self%heap_vars(ithr)%kcorrs ! can be reused
        argmat   => self%heap_vars(ithr)%argmat
        argmat   =  self%argtransf(:self%pftsz,:) * shvec(1) + self%argtransf(self%pftsz + 1:,:) * shvec(2)
        shmat    =  cmplx(cos(argmat),sin(argmat))
        if( self%iseven(self%pinds(iptcl)) )then
            pft_ref = self%pfts_refs_even(:,:,iref)
        else
            pft_ref = self%pfts_refs_odd(:,:,iref)
        endif
        if( self%with_ctf )then
            pft_ref = (pft_ref * self%ctfmats(:,:,self%pinds(iptcl))) * shmat
        else
            pft_ref = pft_ref * shmat
        endif
        euclids(:) = 0.
        do k=params_glob%kfromto(1),params_glob%kstop
            call self%calc_k_corrs(pft_ref, self%pinds(iptcl), k, keuclids)
            sumsqptcl = sum(csq(self%pfts_ptcls(:,k,self%pinds(iptcl))))
            sumsqref  = sum(csq(pft_ref(:,k)))
            euclids(:) = euclids(:) + (2. * keuclids(:) - sumsqptcl - sumsqref ) / (2. * self%sigma2_noise(k, self%pinds(iptcl)))
        end do
    end subroutine gencorrs_euclid_2

    subroutine gencorrs_euclid_mir( self, iref, iptcl, euclids, euclids_mir )
        class(polarft_corrcalc), intent(inout) :: self
        integer,                 intent(in)    :: iref, iptcl
        real(sp),                intent(out)   :: euclids(self%nrots)
        real(sp),                intent(out)   :: euclids_mir(self%nrots)
        complex(sp), pointer :: pft_ref(:,:)
        real(sp),    pointer :: keuclids(:)
        real(sp),    pointer :: keuclids_mir(:)
        real(sp) :: sumsqref, sumsqptcl
        integer  :: k, ithr
        ithr         =  omp_get_thread_num() + 1
        pft_ref      => self%heap_vars(ithr)%pft_ref
        keuclids     => self%heap_vars(ithr)%kcorrs     ! can be reused
        keuclids_mir => self%heap_vars(ithr)%kcorrs_mir ! can be reused
        call self%prep_ref4corr(iref, iptcl, pft_ref, sumsqref, params_glob%kstop)
        euclids(:) = 0.
        do k=params_glob%kfromto(1),params_glob%kstop
            call self%calc_k_corrs(pft_ref, self%pinds(iptcl), k, keuclids, keuclids_mir)
            sumsqptcl = sum(csq(self%pfts_ptcls(:,k,self%pinds(iptcl))))
            sumsqref  = sum(csq(pft_ref(:,k)))
            euclids(:)     = euclids(:) + (2. * keuclids(:)    - sumsqptcl - sumsqref ) / &
                (2. * self%sigma2_noise(k, self%pinds(iptcl)))
            euclids_mir(:) = euclids(:) + (2. * keuclids_mir(:) - sumsqptcl - sumsqref ) / &
                (2. * self%sigma2_noise(k, self%pinds(iptcl)))
        end do
    end subroutine gencorrs_euclid_mir

    subroutine gencorrs_1( self, iref, iptcl, cc )
        class(polarft_corrcalc), intent(inout) :: self
        integer,                 intent(in)    :: iref, iptcl
        real(sp),                intent(out)   :: cc(self%nrots)
        select case(params_glob%cc_objfun)
            case(OBJFUN_CC)
                call self%gencorrs_cc_1(iref, iptcl, cc)
            case(OBJFUN_EUCLID)
                call self%gencorrs_euclid_1(iref, iptcl, cc)
        end select
    end subroutine gencorrs_1

    subroutine gencorrs_2( self, iref, iptcl, shvec, cc )
        class(polarft_corrcalc), intent(inout) :: self
        integer,                 intent(in)    :: iref, iptcl
        real(sp),                intent(in)    :: shvec(2)
        real(sp),                intent(out)   :: cc(self%nrots)
        select case(params_glob%cc_objfun)
            case(OBJFUN_CC)
                call self%gencorrs_cc_2(iref, iptcl, shvec, cc)
            case(OBJFUN_EUCLID)
                call self%gencorrs_euclid_2(iref, iptcl, shvec, cc)
        end select
    end subroutine gencorrs_2

    subroutine gencorrs_mir( self, iref, iptcl, cc, cc_mir )
        class(polarft_corrcalc), intent(inout) :: self
        integer,                 intent(in)    :: iref, iptcl
        real(sp),                intent(out)   :: cc(self%nrots)
        real(sp),                intent(out)   :: cc_mir(self%nrots)
        select case(params_glob%cc_objfun)
            case(OBJFUN_CC)
                call self%gencorrs_cc_mir(iref, iptcl, cc, cc_mir)
            case(OBJFUN_EUCLID)
                call self%gencorrs_euclid_mir(iref, iptcl, cc, cc_mir)
        end select
    end subroutine gencorrs_mir

    function gencorr_for_rot_8( self, iref, iptcl, shvec, irot ) result( cc )
        class(polarft_corrcalc), intent(inout) :: self
        integer,                 intent(in)    :: iref, iptcl
        real(dp),                intent(in)    :: shvec(2)
        integer,                 intent(in)    :: irot
        real(dp)                               :: cc
        select case(params_glob%cc_objfun)
            case(OBJFUN_CC)
                cc = self%gencorr_cc_for_rot_8( iref, iptcl, shvec, irot )
            case(OBJFUN_EUCLID)
                cc = self%gencorr_euclid_for_rot_8( iref, iptcl, shvec, irot )
        end select
    end function gencorr_for_rot_8

    function gencorr_cc_for_rot( self, iref, iptcl, shvec, irot ) result( cc )
        class(polarft_corrcalc), intent(inout) :: self
        integer,                 intent(in)    :: iref, iptcl
        real(sp),                intent(in)    :: shvec(2)
        integer,                 intent(in)    :: irot
        complex(sp), pointer :: pft_ref(:,:), shmat(:,:)
        real(sp),    pointer :: argmat(:,:)
        real(sp) :: cc, corr, sqsum_ref, sqsum_ptcl
        integer  :: ithr, ik
        ithr    =  omp_get_thread_num() + 1
        pft_ref => self%heap_vars(ithr)%pft_ref
        shmat   => self%heap_vars(ithr)%shmat
        argmat  => self%heap_vars(ithr)%argmat
        argmat  =  self%argtransf(:self%pftsz,:) * shvec(1) + self%argtransf(self%pftsz + 1:,:) * shvec(2)
        shmat   =  cmplx(cos(argmat),sin(argmat))
        if( self%iseven(self%pinds(iptcl)) )then
            pft_ref = self%pfts_refs_even(:,:,iref)
        else
            pft_ref = self%pfts_refs_odd(:,:,iref)
        endif
        if( params_glob%l_match_filt )then
            if( self%l_clsfrcs )then
                call self%shellnorm_and_filter_ref(iptcl, iptcl, pft_ref)
            else
                call self%shellnorm_and_filter_ref(iptcl, iref, pft_ref)
            endif
        endif
        if( self%with_ctf )then
            pft_ref = (pft_ref * self%ctfmats(:,:,self%pinds(iptcl))) * shmat
        else
            pft_ref = pft_ref * shmat
        endif
        if( .not.params_glob%l_match_filt )then
            sqsum_ref  = 0.
            sqsum_ptcl = 0.
            corr       = 0.
            do ik = params_glob%kfromto(1),params_glob%kstop
                sqsum_ref  = sqsum_ref  + real(ik) * sum(csq(pft_ref(:,ik)))
                sqsum_ptcl = sqsum_ptcl + real(ik) * sum(csq(self%pfts_ptcls(:,ik,self%pinds(iptcl))))
                corr       = corr + &
                    real(ik) * self%calc_corrk_for_rot(pft_ref, self%pinds(iptcl), ik, irot)
            end do
            cc = corr / sqrt(sqsum_ref * sqsum_ptcl)
        else
            sqsum_ref = sum(csq(pft_ref(:,params_glob%kfromto(1):params_glob%kstop)))
            corr      = self%calc_corr_for_rot(pft_ref, self%pinds(iptcl), irot)
            cc        = corr  / sqrt(sqsum_ref * self%sqsums_ptcls(self%pinds(iptcl)))
        endif
    end function gencorr_cc_for_rot

    function gencorr_cc_for_rot_8( self, iref, iptcl, shvec, irot ) result( cc )
        class(polarft_corrcalc), intent(inout) :: self
        integer,                 intent(in)    :: iref, iptcl
        real(dp),                intent(in)    :: shvec(2)
        integer,                 intent(in)    :: irot
        complex(dp), pointer :: pft_ref(:,:), shmat(:,:)
        real(dp),    pointer :: argmat(:,:)
        real(dp) :: cc, corr, sqsum_ref, sqsum_ptcl
        integer  :: ithr, ik
        ithr    =  omp_get_thread_num() + 1
        pft_ref => self%heap_vars(ithr)%pft_ref_8
        shmat   => self%heap_vars(ithr)%shmat_8
        argmat  => self%heap_vars(ithr)%argmat_8
        argmat  =  self%argtransf(:self%pftsz,:) * shvec(1) + self%argtransf(self%pftsz + 1:,:) * shvec(2)
        shmat   =  cmplx(cos(argmat),sin(argmat),dp)
        if( self%iseven(self%pinds(iptcl)) )then
            pft_ref = self%pfts_refs_even(:,:,iref)
        else
            pft_ref = self%pfts_refs_odd(:,:,iref)
        endif
        if( params_glob%l_match_filt )then
            if( self%l_clsfrcs )then
                call self%shellnorm_and_filter_ref_8(iptcl, iptcl, pft_ref)
            else
                call self%shellnorm_and_filter_ref_8(iptcl, iref, pft_ref)
            endif
        endif
        if( self%with_ctf )then
            pft_ref = (pft_ref * self%ctfmats(:,:,self%pinds(iptcl))) * shmat
        else
            pft_ref = pft_ref * shmat
        endif
        if( .not.params_glob%l_match_filt )then
            sqsum_ref  = 0._dp
            sqsum_ptcl = 0._dp
            corr       = 0._dp
            do ik = params_glob%kfromto(1),params_glob%kstop
                sqsum_ref  = sqsum_ref  + real(ik,kind=dp) * sum(csq(pft_ref(:,ik)))
                sqsum_ptcl = sqsum_ptcl + real(ik,kind=dp) * sum(csq(self%pfts_ptcls(:,ik,self%pinds(iptcl))))
                corr       = corr + &
                    real(ik,kind=dp) * self%calc_corrk_for_rot_8(pft_ref, self%pinds(iptcl), ik, irot)
            end do
            cc = corr / sqrt(sqsum_ref * sqsum_ptcl)
        else
            sqsum_ref = sum(csq(pft_ref(:,params_glob%kfromto(1):params_glob%kstop)))
            corr      = self%calc_corr_for_rot_8(pft_ref, self%pinds(iptcl), irot)
            cc        = corr  / sqrt(sqsum_ref * self%sqsums_ptcls(self%pinds(iptcl)))
        endif
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
        if( self%iseven(self%pinds(iptcl)) )then
            pft_ref  = self%pfts_refs_even(:,:,iref)
            pft_dref = self%pfts_drefs_even(:,:,:,iref)
        else
            pft_ref  = self%pfts_refs_odd(:,:,iref)
            pft_dref = self%pfts_drefs_odd(:,:,:,iref)
        endif
        if( params_glob%l_match_filt )then
            if( self%l_clsfrcs )then
                call self%shellnorm_and_filter_ref_8(iptcl, iptcl, pft_ref)
            else
                call self%shellnorm_and_filter_ref_8(iptcl, iref, pft_ref)
            endif
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
        sqsum_ref = sum(csq(pft_ref(:,params_glob%kfromto(1):params_glob%kstop)))
        denom     = sqrt(sqsum_ref * self%sqsums_ptcls(self%pinds(iptcl)))
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
        if( self%iseven(self%pinds(iptcl)) )then
            pft_ref  = self%pfts_refs_even(:,:,iref)
            pft_dref = self%pfts_drefs_even(:,:,:,iref)
        else
            pft_ref  = self%pfts_refs_odd(:,:,iref)
            pft_dref = self%pfts_drefs_odd(:,:,:,iref)
        endif
        if( params_glob%l_match_filt )then
            if( self%l_clsfrcs )then
                call self%shellnorm_and_filter_ref_8(iptcl, iptcl, pft_ref)
            else
                call self%shellnorm_and_filter_ref_8(iptcl, iref, pft_ref)
            endif
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
        sqsum_ref = sum(csq(pft_ref(:,params_glob%kfromto(1):params_glob%kstop)))
        denom     = sqrt(sqsum_ref * self%sqsums_ptcls(self%pinds(iptcl)))
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

    subroutine gencorr_grad_for_rot_8( self, iref, iptcl, shvec, irot, f, grad )
        class(polarft_corrcalc), intent(inout) :: self
        integer,                 intent(in)    :: iref, iptcl
        real(dp),                intent(in)    :: shvec(2)
        integer,                 intent(in)    :: irot
        real(dp),                intent(out)   :: f, grad(2)
        select case(params_glob%cc_objfun)
            case(OBJFUN_CC)
                call self%gencorr_cc_grad_for_rot_8( iref, iptcl, shvec, irot, f, grad )
            case(OBJFUN_EUCLID)
                call self%gencorr_euclid_grad_for_rot_8( iref, iptcl, shvec, irot, f, grad )
        end select
    end subroutine gencorr_grad_for_rot_8

    subroutine gencorr_cc_grad_for_rot_8( self, iref, iptcl, shvec, irot, f, grad )
        class(polarft_corrcalc), intent(inout) :: self
        integer,                 intent(in)    :: iref, iptcl
        real(dp),                intent(in)    :: shvec(2)
        integer,                 intent(in)    :: irot
        real(dp),                intent(out)   :: f, grad(2)
        complex(dp), pointer :: pft_ref(:,:), pft_ref_tmp(:,:), shmat(:,:)
        real(dp),    pointer :: argmat(:,:)
        real(dp) :: corr, denom, sqsum_ref, sqsum_ptcl
        integer  :: ithr, ik
        ithr        =  omp_get_thread_num() + 1
        pft_ref     => self%heap_vars(ithr)%pft_ref_8
        pft_ref_tmp => self%heap_vars(ithr)%pft_ref_tmp_8
        shmat       => self%heap_vars(ithr)%shmat_8
        argmat      => self%heap_vars(ithr)%argmat_8
        argmat      =  self%argtransf(:self%pftsz,:) * shvec(1) + self%argtransf(self%pftsz + 1:,:) * shvec(2)
        shmat       =  cmplx(cos(argmat),sin(argmat),dp)
        if( self%iseven(self%pinds(iptcl)) )then
            pft_ref = self%pfts_refs_even(:,:,iref)
        else
            pft_ref = self%pfts_refs_odd(:,:,iref)
        endif
        if( params_glob%l_match_filt )then
            if( self%l_clsfrcs )then
                call self%shellnorm_and_filter_ref_8(iptcl, iptcl, pft_ref)
            else
                call self%shellnorm_and_filter_ref_8(iptcl, iref, pft_ref)
            endif
        endif
        if( self%with_ctf )then
            pft_ref = (pft_ref * self%ctfmats(:,:,self%pinds(iptcl))) * shmat
        else
            pft_ref = pft_ref * shmat
        endif
        if( .not.params_glob%l_match_filt )then
            sqsum_ref  = 0._dp
            sqsum_ptcl = 0._dp
            corr       = 0._dp
            grad(1)    = 0._dp
            grad(2)    = 0._dp
            do ik = params_glob%kfromto(1),params_glob%kstop
                sqsum_ref  = sqsum_ref  + real(ik,kind=dp) * sum(csq(pft_ref(:,ik)))
                sqsum_ptcl = sqsum_ptcl + real(ik,kind=dp) * sum(csq(self%pfts_ptcls(:,ik,self%pinds(iptcl))))
                corr       = corr + &
                    real(ik,kind=dp) * self%calc_corrk_for_rot_8(pft_ref, self%pinds(iptcl), ik, irot)
            end do
            pft_ref_tmp = pft_ref * (0., 1.) * self%argtransf(:self%pftsz,:)
            do ik = params_glob%kfromto(1),params_glob%kstop
                grad(1) = grad(1) + &
                    real(ik,kind=dp) * self%calc_corrk_for_rot_8(pft_ref_tmp, self%pinds(iptcl), ik, irot)
            end do
            pft_ref_tmp = pft_ref * (0., 1.) * self%argtransf(self%pftsz + 1:,:)
            do ik = params_glob%kfromto(1),params_glob%kstop
                grad(2) = grad(2) + &
                    real(ik,kind=dp) * self%calc_corrk_for_rot_8(pft_ref_tmp, self%pinds(iptcl), ik, irot)
            end do
            f    = corr / sqrt(sqsum_ref*sqsum_ptcl)
            grad = grad / sqrt(sqsum_ref*sqsum_ptcl)
        else
            denom       = sqrt(sum(csq(pft_ref(:,params_glob%kfromto(1):params_glob%kstop))) * self%sqsums_ptcls(self%pinds(iptcl)))
            corr        = self%calc_corr_for_rot_8(pft_ref, self%pinds(iptcl), irot)
            f           = corr  / denom
            pft_ref_tmp = pft_ref * (0., 1.) * self%argtransf(:self%pftsz,:)
            corr        = self%calc_corr_for_rot_8(pft_ref_tmp, self%pinds(iptcl), irot)
            grad(1)     = corr / denom
            pft_ref_tmp = pft_ref * (0., 1.) * self%argtransf(self%pftsz + 1:,:)
            corr        = self%calc_corr_for_rot_8(pft_ref_tmp, self%pinds(iptcl), irot)
            grad(2)     = corr / denom
        endif
    end subroutine gencorr_cc_grad_for_rot_8

    subroutine gencorr_grad_only_for_rot_8( self, iref, iptcl, shvec, irot, grad )
        class(polarft_corrcalc), intent(inout) :: self
        integer,                 intent(in)    :: iref, iptcl
        real(dp),                intent(in)    :: shvec(2)
        integer,                 intent(in)    :: irot
        real(dp),                intent(out)   :: grad(2)
        select case(params_glob%cc_objfun)
            case(OBJFUN_CC)
                call self%gencorr_cc_grad_only_for_rot_8( iref, iptcl, shvec, irot, grad )
            case(OBJFUN_EUCLID)
                call self%gencorr_euclid_grad_only_for_rot_8( iref, iptcl, shvec, irot, grad )
        end select
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
        integer  :: ithr, ik
        real(dp) :: sqsum_ref, sqsum_ptcl
        ithr        =  omp_get_thread_num() + 1
        pft_ref     => self%heap_vars(ithr)%pft_ref_8
        pft_ref_tmp => self%heap_vars(ithr)%pft_ref_tmp_8
        shmat       => self%heap_vars(ithr)%shmat_8
        argmat      => self%heap_vars(ithr)%argmat_8
        argmat      =  self%argtransf(:self%pftsz,:) * shvec(1) + self%argtransf(self%pftsz + 1:,:) * shvec(2)
        shmat       =  cmplx(cos(argmat),sin(argmat),dp)
        if( self%iseven(self%pinds(iptcl)) )then
            pft_ref = self%pfts_refs_even(:,:,iref)
        else
            pft_ref = self%pfts_refs_odd(:,:,iref)
        endif
        if( params_glob%l_match_filt )then
            if( self%l_clsfrcs )then
                call self%shellnorm_and_filter_ref_8(iptcl, iptcl, pft_ref)
            else
                call self%shellnorm_and_filter_ref_8(iptcl, iref, pft_ref)
            endif
        endif
        if( self%with_ctf )then
            pft_ref = (pft_ref * self%ctfmats(:,:,self%pinds(iptcl))) * shmat
        else
            pft_ref = pft_ref * shmat
        endif
        if( .not.params_glob%l_match_filt )then
            sqsum_ref  = 0._dp
            sqsum_ptcl = 0._dp
            grad(1)    = 0._dp
            grad(2)    = 0._dp
            do ik = params_glob%kfromto(1),params_glob%kstop
                sqsum_ref  = sqsum_ref  + real(ik,kind=dp) * sum(csq(pft_ref(:,ik)))
                sqsum_ptcl = sqsum_ptcl + real(ik,kind=dp) * sum(csq(self%pfts_ptcls(:,ik,self%pinds(iptcl))))
            end do
            pft_ref_tmp = pft_ref * (0., 1.) * self%argtransf(:self%pftsz,:)
            do ik = params_glob%kfromto(1),params_glob%kstop
                grad(1) = grad(1) + &
                    real(ik,kind=dp) * self%calc_corrk_for_rot_8(pft_ref_tmp, self%pinds(iptcl), ik, irot)
            end do
            pft_ref_tmp = pft_ref * (0., 1.) * self%argtransf(self%pftsz + 1:,:)
            do ik = params_glob%kfromto(1),params_glob%kstop
                grad(2) = grad(2) + &
                    real(ik,kind=dp) * self%calc_corrk_for_rot_8(pft_ref_tmp, self%pinds(iptcl), ik, irot)
            end do
            grad = grad / sqrt(sqsum_ref*sqsum_ptcl)
        else
            denom       = sqrt(sum(csq(pft_ref(:,params_glob%kfromto(1):params_glob%kstop))) * self%sqsums_ptcls(self%pinds(iptcl)))
            pft_ref_tmp = pft_ref * (0., 1.) * self%argtransf(:self%pftsz,:)
            corr        = self%calc_corr_for_rot_8(pft_ref_tmp, self%pinds(iptcl), irot)
            grad(1)     = corr / denom
            pft_ref_tmp = pft_ref * (0., 1.) * self%argtransf(self%pftsz + 1:,:)
            corr        = self%calc_corr_for_rot_8(pft_ref_tmp, self%pinds(iptcl), irot)
            grad(2)     = corr / denom
        endif
    end subroutine gencorr_cc_grad_only_for_rot_8

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
        if( self%iseven(self%pinds(iptcl)) )then
            pft_ref = self%pfts_refs_even(:,:,iref)
        else
            pft_ref = self%pfts_refs_odd (:,:,iref)
        endif
        if( self%with_ctf )then
            pft_ref = (pft_ref * self%ctfmats(:,:,self%pinds(iptcl))) * shmat
        else
            pft_ref = pft_ref * shmat
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
        if( self%iseven(self%pinds(iptcl)) )then
            pft_ref = self%pfts_refs_even(:,:,iref)
        else
            pft_ref = self%pfts_refs_odd(:,:,iref)
        endif
        if( self%with_ctf )then
            pft_ref = (pft_ref * self%ctfmats(:,:,self%pinds(iptcl))) * shmat
        else
            pft_ref = pft_ref * shmat
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
        integer  :: ithr, k
        ithr        =  omp_get_thread_num() + 1
        pft_ref     => self%heap_vars(ithr)%pft_ref_8
        pft_ref_tmp => self%heap_vars(ithr)%pft_ref_tmp_8
        shmat       => self%heap_vars(ithr)%shmat_8
        argmat      => self%heap_vars(ithr)%argmat_8
        argmat      =  self%argtransf(:self%pftsz,:) * shvec(1) + self%argtransf(self%pftsz + 1:,:) * shvec(2)
        shmat       =  cmplx(cos(argmat),sin(argmat),dp)
        if( self%iseven(self%pinds(iptcl)) )then
            pft_ref = self%pfts_refs_even(:,:,iref)
        else
            pft_ref = self%pfts_refs_odd(:,:,iref)
        endif
        if( self%with_ctf )then
            pft_ref = (pft_ref * self%ctfmats(:,:,self%pinds(iptcl))) * shmat
        else
            pft_ref = pft_ref * shmat
        endif
        f           = self%calc_euclid_for_rot_8(pft_ref, self%pinds(iptcl), irot)
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
        integer  :: ithr, k
        ithr        =  omp_get_thread_num() + 1
        pft_ref     => self%heap_vars(ithr)%pft_ref_8
        pft_ref_tmp => self%heap_vars(ithr)%pft_ref_tmp_8
        shmat       => self%heap_vars(ithr)%shmat_8
        argmat      => self%heap_vars(ithr)%argmat_8
        argmat      =  self%argtransf(:self%pftsz,:) * shvec(1) + self%argtransf(self%pftsz + 1:,:) * shvec(2)
        shmat       =  cmplx(cos(argmat),sin(argmat),dp)
        if( self%iseven(self%pinds(iptcl)) )then
            pft_ref = self%pfts_refs_even(:,:,iref)
        else
            pft_ref = self%pfts_refs_odd(:,:,iref)
        endif
        if( self%with_ctf )then
            pft_ref = (pft_ref * self%ctfmats(:,:,self%pinds(iptcl))) * shmat
        else
            pft_ref = pft_ref * shmat
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
        integer  :: i, ithr, k
        i       =  self%pinds(iptcl)
        ithr    =  omp_get_thread_num() + 1
        pft_ref => self%heap_vars(ithr)%pft_ref
        shmat   => self%heap_vars(ithr)%shmat
        argmat  => self%heap_vars(ithr)%argmat
        argmat  =  self%argtransf(:self%pftsz,:) * shvec(1) + self%argtransf(self%pftsz + 1:,:) * shvec(2)
        shmat   =  cmplx(cos(argmat),sin(argmat))
        if( self%iseven(i) )then
            pft_ref = self%pfts_refs_even(:,:,iref)
        else
            pft_ref = self%pfts_refs_odd (:,:,iref)
        endif
        if( self%with_ctf )then
            pft_ref = (pft_ref * self%ctfmats(:,:,i)) * shmat
        else
            pft_ref = pft_ref * shmat
        endif
        do k = params_glob%kfromto(1), params_glob%kfromto(2)
            sigma_contrib(k) = self%calc_euclidk_for_rot(pft_ref, i, k, irot)
        end do
        sigma_contrib = 0.5 * sigma_contrib / real(self%pftsz)
    end subroutine gencorr_sigma_contrib

    real function specscore_1( self, iref, iptcl, irot )
        class(polarft_corrcalc), intent(inout) :: self
        integer,                 intent(in)    :: iref, iptcl, irot
        real :: frc(params_glob%kfromto(1):params_glob%kstop)
        call self%genfrc(iref, iptcl, irot, frc)
        specscore_1 = max(0.,min(real(count(frc>0.143)) / real(size(frc)),1.))
    end function specscore_1

    real function specscore_2( self, iref, iptcl, irot, shvec )
        class(polarft_corrcalc), intent(inout) :: self
        integer,                 intent(in)    :: iref, iptcl, irot
        real,                    intent(in)    :: shvec(2)
        real :: frc(params_glob%kfromto(1):params_glob%kstop)
        call self%calc_frc(iref, iptcl, irot, shvec, frc )
        specscore_2 = max(0.,min(real(count(frc>0.143)) / real(size(frc)),1.))
    end function specscore_2

    function calc_roinv_corrmat( self) result(corrmat )
        ! assumes particles/references in pftcc identical sets
        class(polarft_corrcalc), intent(inout) :: self
        real, allocatable   :: corrmat(:,:)
        real    :: corrs(self%nrots)
        integer :: iref, iptcl, loc(1)
        if( self%nptcls /= self%nrefs ) THROW_HARD('nptcls == nrefs in pftcc required; calc_roinv_corrmat')
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
                call fftwf_free(self%fftdat(ithr)%p_ref_fft_re_2)
                call fftwf_free(self%fftdat(ithr)%p_ref_fft_im_2)
                call fftwf_free(self%fftdat(ithr)%p_product_fft)
                call fftwf_free(self%fftdat(ithr)%p_backtransf)
                call fftwf_free(self%fft_carray(ithr)%p_re)
                call fftwf_free(self%fft_carray(ithr)%p_im)
                deallocate(self%heap_vars(ithr)%pft_ref,self%heap_vars(ithr)%pft_ref_tmp,&
                    &self%heap_vars(ithr)%pft_ref_tmp1, self%heap_vars(ithr)%pft_ref_tmp2,&
                    &self%heap_vars(ithr)%pft_dref,&
                    &self%heap_vars(ithr)%corrs_over_k,self%heap_vars(ithr)%argmat,&
                    &self%heap_vars(ithr)%corrs_mir_over_k,&
                    &self%heap_vars(ithr)%shmat,self%heap_vars(ithr)%kcorrs,&
                    &self%heap_vars(ithr)%kcorrs_mir,&
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
            if( allocated(self%ctfmats)    ) deallocate(self%ctfmats)
            if( allocated(self%ref_optlp)  ) deallocate(self%ref_optlp)
            if( allocated(self%pssnr_filt) ) deallocate(self%pssnr_filt)
            deallocate( self%sqsums_ptcls, self%angtab, self%argtransf,&
                &self%polar, self%pfts_refs_even, self%pfts_refs_odd, self%pfts_drefs_even, self%pfts_drefs_odd,&
                self%pfts_ptcls, self%fft_factors, self%fftdat, self%fftdat_ptcls, self%fft_carray,&
                &self%iseven, self%pinds, self%heap_vars)
            call fftwf_destroy_plan(self%plan_bwd)
            call fftwf_destroy_plan(self%plan_fwd_1)
            call fftwf_destroy_plan(self%plan_fwd_2)
            self%sigma2_noise      => null()
            self%existence = .false.
        endif
    end subroutine kill

end module simple_polarft_corrcalc
