! for calculation of band-pass limited cross-correlation of polar Fourier transforms
module simple_polarft_corrcalc
!$ use omp_lib
!$ use omp_lib_kinds
include 'simple_lib.f08'
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
    real,        pointer :: kcorrs(:)          => null()
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
    private
    integer                          :: nptcls     = 1              !< the total number of particles in partition (logically indexded [fromp,top])
    integer                          :: nrefs      = 1              !< the number of references (logically indexded [1,nrefs])
    integer                          :: nrots      = 0              !< number of in-plane rotations for one pft (determined by radius of molecule)
    integer                          :: pftsz      = 0              !< size of reference and particle pft (nrots/2)
    integer                          :: pfromto(2) = 0              !< particle index range
    integer                          :: ldim(3)    = 0              !< logical dimensions of original cartesian image
    integer                          :: kfromto(2)                  !< band-pass Fourier index limits
    integer,             allocatable :: pinds(:)                    !< index array (to reduce memory when frac_update < 1)
    real,                allocatable :: npix_per_shell(:)           !< number of (cartesian) pixels per shell
    real(sp),            allocatable :: sqsums_ptcls(:)             !< memoized square sums for the correlation calculations (taken from kfromto(1):kfromto(2))
    real(sp),            allocatable :: angtab(:)                   !< table of in-plane angles (in degrees)
    real(dp),            allocatable :: argtransf(:,:)              !< argument transfer constants for shifting the references
    real(sp),            allocatable :: polar(:,:)                  !< table of polar coordinates (in Cartesian coordinates)
    real(sp),            allocatable :: ctfmats(:,:,:)              !< expand set of CTF matrices (for efficient parallel exec)
    real(dp),            allocatable :: argtransf_shellone(:)       !< one dimensional argument transfer constants (shell k=1) for shifting the references
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
    ! MODIFIERS
    procedure          :: shift_ptcl
    ! MEMOIZER
    procedure, private :: memoize_sqsum_ptcl
    procedure, private :: memoize_fft
    procedure          :: memoize_ffts
    procedure, private :: setup_npix_per_shell
    ! CALCULATORS
    procedure          :: create_polar_absctfmats
    procedure, private :: prep_ref4corr_sp
    procedure, private :: prep_ref4corr_dp
    generic            :: prep_ref4corr => prep_ref4corr_sp, prep_ref4corr_dp
    procedure, private :: gen_shmat
    procedure, private :: gen_shmat_8
    procedure, private :: calc_corrs_over_k
    procedure, private :: calc_k_corrs
    procedure, private :: calc_corr_for_rot
    procedure, private :: calc_corr_for_rot_8
    procedure, private :: calc_T1_T2_for_rot_8
    procedure, private :: calc_euclid_for_rot
    procedure, private :: calc_euclid_for_rot_8
    procedure, private :: calc_corrk_for_rot
    procedure, private :: calc_corrk_for_rot_8
    procedure, private :: calc_euclidk_for_rot
    procedure, private :: gencorrs_cc
    procedure          :: genmaxcorr_comlin
    procedure, private :: gencorrs_euclid
    procedure, private :: gencorrs_1
    procedure, private :: gencorrs_2
    generic            :: gencorrs => gencorrs_1, gencorrs_2
    procedure          :: gencorr_for_rot
    procedure          :: gencorr_for_rot_8
    procedure          :: gencorr_grad_for_rot_8
    procedure          :: gencorr_grad_only_for_rot_8
    procedure          :: gencorr_cc_for_rot
    procedure          :: gencorr_cc_for_rot_8
    procedure          :: gencorr_cont_grad_cc_for_rot_8
    procedure          :: gencorr_cont_cc_for_rot_8
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
    procedure, private :: specscore_1, specscore_2
    generic            :: specscore => specscore_1, specscore_2
    procedure, private :: weight_ref_ptcl_sp, weight_ref_ptcl_dp
    generic,   private :: weight_ref_ptcl => weight_ref_ptcl_sp, weight_ref_ptcl_dp
    procedure, private :: deweight_ref_ptcl_sp, deweight_ref_ptcl_dp
    generic,   private :: deweight_ref_ptcl => deweight_ref_ptcl_sp, deweight_ref_ptcl_dp
    ! DESTRUCTOR
    procedure          :: kill
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
        real(sp)              :: ang, freq
        integer(kind=c_int)   :: wsdm_ret
        integer               :: local_stat,irot, k, ithr, i, ik, cnt
        logical               :: even_dims, test(2)
        ! kill possibly pre-existing object
        call self%kill
        ! set particle index range
        self%pfromto = pfromto
        ! set band-pass Fourier indexc limits
        self%kfromto = kfromto
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
                    &self%sqsums_ptcls(1:self%nptcls),self%fftdat(params_glob%nthr),self%fft_carray(params_glob%nthr),&
                    &self%fftdat_ptcls(1:self%nptcls,self%kfromto(1):self%kfromto(2)),&
                    &self%heap_vars(params_glob%nthr) )
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
        ! setup npix_per_shell
        call self%setup_npix_per_shell
        ! flag existence
        self%existence = .true.
        ! set pointer to global instance
        pftcc_glob => self
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
            if( allocated(self%iseven) )       deallocate(self%iseven)
            if( allocated(self%pfts_ptcls) )   deallocate(self%pfts_ptcls)
            if( allocated(self%fftdat_ptcls) )then
                do i = 1, size(self%fftdat_ptcls,dim=1)
                    do ik = self%kfromto(1),self%kfromto(2)
                        call fftwf_free(self%fftdat_ptcls(i,ik)%p_re)
                        call fftwf_free(self%fftdat_ptcls(i,ik)%p_im)
                    end do
                end do
                deallocate(self%fftdat_ptcls)
            endif
            allocate( self%pfts_ptcls(self%pftsz,self%kfromto(1):self%kfromto(2),1:self%nptcls),&
                        &self%sqsums_ptcls(1:self%nptcls),self%iseven(1:self%nptcls),&
                        &self%fftdat_ptcls(1:self%nptcls,self%kfromto(1):self%kfromto(2)) )
            do i = 1,self%nptcls
                do ik = self%kfromto(1),self%kfromto(2)
                    self%fftdat_ptcls(i,ik)%p_re = fftwf_alloc_complex(int(self%pftsz, c_size_t))
                    self%fftdat_ptcls(i,ik)%p_im = fftwf_alloc_complex(int(self%pftsz, c_size_t))
                    call c_f_pointer(self%fftdat_ptcls(i,ik)%p_re, self%fftdat_ptcls(i,ik)%re, [self%pftsz])
                    call c_f_pointer(self%fftdat_ptcls(i,ik)%p_im, self%fftdat_ptcls(i,ik)%im, [self%pftsz])
                end do
            end do
            endif
            self%pfts_ptcls   = zero
            self%sqsums_ptcls = 0.
            self%iseven       = .true.
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

    !>  \brief  returns polar Fourier transform of particle iptcl
    function get_ptcl_pft( self, iptcl) result( pft )
        class(polarft_corrcalc), intent(in) :: self
        integer,                 intent(in) :: iptcl
        complex(sp), allocatable :: pft(:,:)
        allocate(pft(self%pftsz,self%kfromto(1):self%kfromto(2)),&
        source=self%pfts_ptcls(:,:,self%pinds(iptcl)))
    end function get_ptcl_pft

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
        ithr  =  omp_get_thread_num() + 1
        i     = self%pinds(iptcl)
        shmat => self%heap_vars(ithr)%shmat
        call self%gen_shmat(ithr, shvec, shmat)
        self%pfts_ptcls(:,:,i) = self%pfts_ptcls(:,:,i) * shmat
    end subroutine shift_ptcl

    ! MEMOIZERS

    subroutine memoize_sqsum_ptcl( self, i )
        class(polarft_corrcalc), intent(inout) :: self
        integer,                 intent(in)    :: i
        self%sqsums_ptcls(i) = sum(csq_fast(self%pfts_ptcls(:,self%kfromto(1):self%kfromto(2),i)))
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
            call fftwf_execute_dft_r2c(self%plan_fwd_1, self%fft_carray(ithr)%re, self%fftdat_ptcls(i,ik)%re)
            call fftwf_execute_dft    (self%plan_fwd_2, self%fft_carray(ithr)%im, self%fftdat_ptcls(i,ik)%im)
        end do
    end subroutine memoize_fft

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

    subroutine prep_ref4corr_sp( self, iref, iptcl, pft_ref, i, ithr)
        class(polarft_corrcalc), intent(inout) :: self
        integer,                 intent(in)    :: iref, iptcl
        complex(sp), pointer,    intent(out)   :: pft_ref(:, :)
        integer,                 intent(out)   :: i, ithr
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

    subroutine prep_ref4corr_dp( self, iref, iptcl, pft_ref_8, i, ithr)
        class(polarft_corrcalc), intent(inout) :: self
        integer,                 intent(in)    :: iref, iptcl
        complex(dp), pointer,    intent(out)   :: pft_ref_8(:, :)
        integer,                 intent(out)   :: i, ithr
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
        self%heap_vars(ithr)%argvec       = self%argtransf(:self%pftsz,  self%kfromto(1)) * shift_8(1) +&
                                            & self%argtransf(self%pftsz+1:,self%kfromto(1)) * shift_8(2)
        shmat_8(:,self%kfromto(1)) = dcmplx(dcos(self%heap_vars(ithr)%argvec), dsin(self%heap_vars(ithr)%argvec))
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

    subroutine calc_corrs_over_k( self, pft_ref, i, corrs_over_k)
        class(polarft_corrcalc), intent(in)  :: self
        integer,                 intent(in)  :: i
        complex(sp),             intent(in)  :: pft_ref(1:self%pftsz,self%kfromto(1):self%kfromto(2))
        real,                    intent(out) :: corrs_over_k(self%nrots)
        integer :: ithr, ik
        ! get thread index
        ithr = omp_get_thread_num() + 1
        ! sum up correlations over k-rings
        corrs_over_k = 0.
        do ik = self%kfromto(1),self%kfromto(2)
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

    function calc_corr_for_rot( self, pft_ref, i, irot )result( corr )
        class(polarft_corrcalc), intent(inout) :: self
        integer,                 intent(in)    :: i, irot
        complex(sp),             intent(in)    :: pft_ref(1:self%pftsz,self%kfromto(1):self%kfromto(2))
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
            tmp = sum( pft_ref(:,self%kfromto(1):self%kfromto(2)) * conjg(self%pfts_ptcls(:,self%kfromto(1):self%kfromto(2),i)))
        else if( irot <= self%pftsz )then
            tmp =       sum( pft_ref(               1:self%pftsz-rot+1,self%kfromto(1):self%kfromto(2)) * conjg(self%pfts_ptcls(rot:self%pftsz,self%kfromto(1):self%kfromto(2),i)))
            tmp = tmp + sum( pft_ref(self%pftsz-rot+2:self%pftsz,      self%kfromto(1):self%kfromto(2)) *       self%pfts_ptcls(  1:rot-1,     self%kfromto(1):self%kfromto(2),i))
        else if( irot == self%pftsz + 1 )then
            tmp = sum( pft_ref(:,self%kfromto(1):self%kfromto(2)) * self%pfts_ptcls(:,self%kfromto(1):self%kfromto(2),i) )
        else
            tmp =       sum( pft_ref(1:self%pftsz-rot+1,self%kfromto(1):self%kfromto(2))          *        self%pfts_ptcls(rot:self%pftsz,self%kfromto(1):self%kfromto(2),i))
            tmp = tmp + sum( pft_ref(self%pftsz-rot+2:self%pftsz,self%kfromto(1):self%kfromto(2)) * conjg( self%pfts_ptcls( 1:rot-1,      self%kfromto(1):self%kfromto(2),i)))
        end if
        corr = real(tmp)
    end function calc_corr_for_rot

    function calc_corr_for_rot_8( self, pft_ref, i, irot )result( corr )
        class(polarft_corrcalc), intent(inout) :: self
        integer,                 intent(in)    :: i, irot
        complex(dp),             intent(in)    :: pft_ref(1:self%pftsz,self%kfromto(1):self%kfromto(2))
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
            tmp = sum( pft_ref(:,self%kfromto(1):self%kfromto(2)) * conjg(self%pfts_ptcls(:,self%kfromto(1):self%kfromto(2),i)))
        else if( irot <= self%pftsz )then
            tmp =       sum( pft_ref(               1:self%pftsz-rot+1,self%kfromto(1):self%kfromto(2)) * conjg(self%pfts_ptcls(rot:self%pftsz,self%kfromto(1):self%kfromto(2),i)))
            tmp = tmp + sum( pft_ref(self%pftsz-rot+2:self%pftsz,      self%kfromto(1):self%kfromto(2)) *       self%pfts_ptcls(  1:rot-1,     self%kfromto(1):self%kfromto(2),i))
        else if( irot == self%pftsz + 1 ) then
            tmp = sum( pft_ref(:,self%kfromto(1):self%kfromto(2)) * self%pfts_ptcls(:,self%kfromto(1):self%kfromto(2),i) )
        else
            tmp =       sum( pft_ref(               1:self%pftsz-rot+1,self%kfromto(1):self%kfromto(2)) *        self%pfts_ptcls(rot:self%pftsz,self%kfromto(1):self%kfromto(2),i))
            tmp = tmp + sum( pft_ref(self%pftsz-rot+2:self%pftsz,      self%kfromto(1):self%kfromto(2)) * conjg( self%pfts_ptcls( 1:rot-1,      self%kfromto(1):self%kfromto(2),i)))
        end if
        corr = real(tmp, kind=dp)
    end function calc_corr_for_rot_8

    !<  \brief  compute the terms T1, T2 necessary for finding the derivative of the correlations, double precision
    subroutine calc_T1_T2_for_rot_8( self, pft_ref, pft_dref, i, irot, nderivs, T1, T2)
        class(polarft_corrcalc), intent(inout) :: self
        integer,                 intent(in)    :: i, irot, nderivs
        complex(dp),             intent(in)    :: pft_ref( 1:self%pftsz,self%kfromto(1):self%kfromto(2)), &
                                                  pft_dref(1:self%pftsz,self%kfromto(1):self%kfromto(2),nderivs)
        real(dp),                intent(out)   :: T1(self%kfromto(1):self%kfromto(2),nderivs), T2(self%kfromto(1):self%kfromto(2),nderivs)
        integer :: k, rot, j
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
                    T1(k,j) = T1(k,j) + real(sum( pft_dref(self%pftsz-rot+2:self%pftsz,      k,j) *       self%pfts_ptcls(  1:rot-1,     k,i)), kind=dp)
                else if( irot == self%pftsz + 1 ) then
                    T1(k,j) = real(sum( pft_dref(:,k,j) * self%pfts_ptcls(:,k,i) ), kind=dp)
                else
                    T1(k,j) =           real(sum( pft_dref(               1:self%pftsz-rot+1,k,j) *       self%pfts_ptcls(rot:self%pftsz,k,i)), kind=dp)
                    T1(k,j) = T1(k,j) + real(sum( pft_dref(self%pftsz-rot+2:self%pftsz,      k,j) * conjg(self%pfts_ptcls(  1:rot-1,     k,i)) ), kind=dp)
                end if
                T2(k,j) = real(sum( pft_dref(:,k,j) * conjg(pft_ref(:,k))), kind=dp)
            end do
        end do
    end subroutine calc_T1_T2_for_rot_8

    function calc_euclid_for_rot( self, pft_ref, iptcl, i, irot ) result( euclid )
        class(polarft_corrcalc), intent(inout) :: self
        integer,                 intent(in)    :: iptcl, i, irot
        complex(sp),             intent(in)    :: pft_ref(1:self%pftsz,self%kfromto(1):self%kfromto(2))
        integer     :: rot, k
        real(sp)    :: euclid, tmp
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
                tmp =       sum(csq_fast(pft_ref(1:self%pftsz-rot+1,k) - self%pfts_ptcls(rot:self%pftsz,k,i)))
                tmp = tmp + sum(csq_fast(pft_ref(self%pftsz-rot+2:self%pftsz,k) - conjg(self%pfts_ptcls(1:rot-1,k,i))))
            else if( irot == self%pftsz + 1 )then
                tmp = sum(csq_fast(pft_ref(:,k) - conjg(self%pfts_ptcls(:,k,i))))
            else
                tmp =       sum(csq_fast(pft_ref(1:self%pftsz-rot+1,k) - conjg(self%pfts_ptcls(rot:self%pftsz,k,i))))
                tmp = tmp + sum(csq_fast(pft_ref(self%pftsz-rot+2:self%pftsz,k) - self%pfts_ptcls(1:rot-1,k,i)))
            end if
            euclid = euclid + tmp
        end do
    end function calc_euclid_for_rot

    function calc_euclid_for_rot_8( self, pft_ref, iptcl, i, irot ) result( euclid )
        class(polarft_corrcalc), intent(inout) :: self
        integer,                 intent(in)    :: iptcl, i, irot
        complex(dp),             intent(in)    :: pft_ref(1:self%pftsz,self%kfromto(1):self%kfromto(2))
        integer  :: rot, k
        real(dp) :: euclid, tmp
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
                tmp =       sum(csq_fast(pft_ref(1:self%pftsz-rot+1,k) - self%pfts_ptcls(rot:self%pftsz,k,i)))
                tmp = tmp + sum(csq_fast(pft_ref(self%pftsz-rot+2:self%pftsz,k) - conjg(self%pfts_ptcls(1:rot-1,k,i))))
            else if( irot == self%pftsz + 1 )then
                tmp = sum(csq_fast(pft_ref(:,k) - conjg(self%pfts_ptcls(:,k,i))))
            else
                tmp =       sum(csq_fast(pft_ref(1:self%pftsz-rot+1,k) - conjg(self%pfts_ptcls(rot:self%pftsz,k,i))))
                tmp = tmp + sum(csq_fast(pft_ref(self%pftsz-rot+2:self%pftsz,k) - self%pfts_ptcls(1:rot-1,k,i)))
            end if
            euclid = euclid + tmp
        end do
    end function calc_euclid_for_rot_8

    function calc_corrk_for_rot( self, pft_ref, i, k, irot ) result( corr )
        class(polarft_corrcalc), intent(inout) :: self
        integer,                 intent(in)    :: i, k, irot
        complex,                 intent(in)    :: pft_ref(1:self%pftsz,self%kfromto(1):self%kfromto(2))
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
        complex(dp),             intent(in)    :: pft_ref(1:self%pftsz,self%kfromto(1):self%kfromto(2))
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
        complex(sp),             intent(in)    :: pft_ref(1:self%pftsz,self%kfromto(1):self%kfromto(2))
        real    :: euclid
        integer :: rot
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

    subroutine genfrc( self, iref, iptcl, irot, frc )
        class(polarft_corrcalc),  intent(inout) :: self
        integer,                  intent(in)    :: iref, iptcl, irot
        real(sp),                 intent(out)   :: frc(self%kfromto(1):self%kfromto(2))
        complex(sp), pointer :: pft_ref(:,:)
        real(sp),    pointer :: kcorrs(:)
        real(sp) :: sumsqref, sumsqptcl, denom
        integer  :: k, ithr, i
        call self%prep_ref4corr(iref, iptcl, pft_ref, i, ithr)
        kcorrs => self%heap_vars(ithr)%kcorrs
        do k=self%kfromto(1),self%kfromto(2)
            call self%calc_k_corrs(pft_ref, i, k, kcorrs)
            sumsqptcl = sum(csq_fast(self%pfts_ptcls(:,k,i)))
            sumsqref  = sum(csq_fast(pft_ref(:,k)))
            denom     = sqrt(sumsqptcl * sumsqref)
            if( denom < 1.e-12 )then
                frc(k) = 0.
            else
                frc(k) = kcorrs(irot) / denom
            endif
        end do
    end subroutine genfrc

    subroutine calc_frc( self, iref, iptcl, irot, shvec, frc )
        class(polarft_corrcalc),  intent(inout) :: self
        integer,                  intent(in)    :: iref, iptcl, irot
        real(sp),                 intent(in)    :: shvec(2)
        real(sp),                 intent(out)   :: frc(self%kfromto(1):self%kfromto(2))
        complex(sp), pointer :: pft_ref(:,:), shmat(:,:)
        real(sp),    pointer :: kcorrs(:)
        real(sp) :: sumsqref, sumsqptcl, denom
        integer  :: k, ithr, i
        call self%prep_ref4corr(iref, iptcl, pft_ref, i, ithr)
        kcorrs  => self%heap_vars(ithr)%kcorrs
        shmat   => self%heap_vars(ithr)%shmat
        call self%gen_shmat(ithr, shvec, shmat)
        pft_ref =  pft_ref * shmat
        do k=self%kfromto(1),self%kfromto(2)
            call self%calc_k_corrs(pft_ref, i, k, kcorrs)
            sumsqptcl = sum(csq_fast(self%pfts_ptcls(:,k,i)))
            sumsqref  = sum(csq_fast(pft_ref(:,k)))
            denom     = sqrt(sumsqptcl * sumsqref)
            if( denom < 1.e-12 )then
                frc(k) = 0.
            else
                frc(k) = kcorrs(irot) / denom
            endif
        end do
    end subroutine calc_frc

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
                cc      = sum( real(pft_ref_i(i,self%kfromto(1):self%kfromto(2)) *&
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
        integer :: i, ithr
        call self%prep_ref4corr(iref, iptcl, pft_ref, i, ithr)
        select case(params_glob%cc_objfun)
            case(OBJFUN_CC)
                call self%gencorrs_cc(    pft_ref, iptcl, i, cc)
            case(OBJFUN_EUCLID)
                call self%gencorrs_euclid(pft_ref, self%heap_vars(ithr)%kcorrs, iptcl, i, cc)
        end select
    end subroutine gencorrs_1

    subroutine gencorrs_2( self, iref, iptcl, shvec, cc )
        class(polarft_corrcalc), intent(inout) :: self
        integer,                 intent(in)    :: iref, iptcl
        real(sp),                intent(in)    :: shvec(2)
        real(sp),                intent(out)   :: cc(self%nrots)
        complex(sp), pointer :: pft_ref(:,:), shmat(:,:)
        integer              :: ithr, i
        call self%prep_ref4corr(iref, iptcl, pft_ref, i, ithr)
        shmat => self%heap_vars(ithr)%shmat
        call self%gen_shmat(ithr, shvec, shmat)
        pft_ref = pft_ref * shmat
        select case(params_glob%cc_objfun)
            case(OBJFUN_CC)
                call self%gencorrs_cc(    pft_ref, iptcl, i, cc)
            case(OBJFUN_EUCLID)
                call self%gencorrs_euclid(pft_ref, self%heap_vars(ithr)%kcorrs, iptcl, i, cc)
        end select
    end subroutine gencorrs_2

    subroutine gencorrs_cc( self, pft_ref, iptcl, i, cc )
        class(polarft_corrcalc), intent(inout) :: self
        complex(sp), pointer,    intent(in)    :: pft_ref(:,:)
        integer,                 intent(in)    :: iptcl, i
        real,                    intent(out)   :: cc(self%nrots)
        real(sp) :: sqsum_ref
        sqsum_ref = sum(csq_fast(pft_ref(:, self%kfromto(1):self%kfromto(2))))
        call self%calc_corrs_over_k(pft_ref, i, cc)
        cc = cc / sqrt(sqsum_ref * self%sqsums_ptcls(i))
    end subroutine gencorrs_cc

    subroutine gencorrs_euclid( self, pft_ref, keuclids, iptcl, i, euclids )
        class(polarft_corrcalc), intent(inout) :: self
        complex(sp), pointer,    intent(inout) :: pft_ref(:,:)
        real(sp),    pointer,    intent(inout) :: keuclids(:)
        integer,                 intent(in)    :: iptcl, i
        real(sp),                intent(out)   :: euclids(self%nrots)
        real(sp) :: sumsqref, sumsqptcl
        real(dp) :: denom
        integer  :: k
        call self%weight_ref_ptcl(pft_ref, i, iptcl)
        euclids(:) = 0.
        denom      = sum(real(csq_fast(self%pfts_ptcls(:, self%kfromto(1):self%kfromto(2),i)), dp))
        do k=self%kfromto(1),self%kfromto(2)
            call self%calc_k_corrs(pft_ref, i, k, keuclids)
            sumsqptcl = sum(csq_fast(self%pfts_ptcls(:,k,i)))
            sumsqref  = sum(csq_fast(pft_ref(:,k)))
            euclids   = euclids + (sumsqptcl + sumsqref - 2. * keuclids(:))
        end do
        euclids = exp( - euclids/denom )
        call self%deweight_ref_ptcl(pft_ref, i, iptcl)
    end subroutine gencorrs_euclid

    real function gencorr_for_rot( self, iref, iptcl, shvec, irot )
        class(polarft_corrcalc), intent(inout) :: self
        integer,                 intent(in)    :: iref, iptcl
        real    ,                intent(in)    :: shvec(2)
        integer,                 intent(in)    :: irot
        complex(sp), pointer :: pft_ref(:,:), shmat(:,:)
        integer              :: ithr, i
        call self%prep_ref4corr(iref, iptcl, pft_ref, i, ithr)
        shmat => self%heap_vars(ithr)%shmat
        call self%gen_shmat(ithr, shvec, shmat)
        pft_ref = pft_ref * shmat
        select case(params_glob%cc_objfun)
            case(OBJFUN_CC)
                gencorr_for_rot = self%gencorr_cc_for_rot(     pft_ref, iptcl, i, irot )
            case(OBJFUN_EUCLID)
                gencorr_for_rot = self%gencorr_euclid_for_rot( pft_ref, iptcl, i, irot)
        end select
    end function gencorr_for_rot

    function gencorr_cc_for_rot( self, pft_ref, iptcl, i, irot ) result( cc )
        class(polarft_corrcalc), intent(inout) :: self
        complex(sp), pointer,    intent(inout) :: pft_ref(:,:)
        integer,                 intent(in)    :: iptcl, i, irot
        real(sp) :: cc, sqsum_ref, sqsum_ptcl
        integer  :: ik
        sqsum_ref  = 0.
        sqsum_ptcl = 0.
        cc         = 0.
        do ik = self%kfromto(1),self%kfromto(2)
            sqsum_ref  = sqsum_ref  + real(ik) * sum(csq_fast(pft_ref(:,ik)))
            sqsum_ptcl = sqsum_ptcl + real(ik) * sum(csq_fast(self%pfts_ptcls(:,ik,i)))
            cc         = cc         + real(ik) * self%calc_corrk_for_rot(pft_ref, i, ik, irot)
        end do
        cc = cc / sqrt(sqsum_ref * sqsum_ptcl)
    end function gencorr_cc_for_rot

    real(sp) function gencorr_euclid_for_rot( self, pft_ref, iptcl, i, irot )
        class(polarft_corrcalc), intent(inout) :: self
        complex(sp), pointer,    intent(inout) :: pft_ref(:,:)
        integer,                 intent(in)    :: iptcl, i, irot
        real(dp) :: denom
        call self%weight_ref_ptcl(pft_ref, i, iptcl)
        denom                  = sum(real(csq_fast(self%pfts_ptcls(:, self%kfromto(1):self%kfromto(2),i)), dp))
        gencorr_euclid_for_rot = exp( - self%calc_euclid_for_rot(pft_ref, iptcl, i, irot)/denom )
        call self%deweight_ref_ptcl(pft_ref, i, iptcl)
    end function gencorr_euclid_for_rot

    real(dp) function gencorr_for_rot_8( self, iref, iptcl, shvec, irot )
        class(polarft_corrcalc), intent(inout) :: self
        integer,                 intent(in)    :: iref, iptcl
        real(dp),                intent(in)    :: shvec(2)
        integer,                 intent(in)    :: irot
        complex(dp), pointer :: pft_ref_8(:,:), shmat_8(:,:)
        integer              :: ithr, i
        call self%prep_ref4corr(iref, iptcl, pft_ref_8, i, ithr)
        shmat_8 => self%heap_vars(ithr)%shmat_8
        call self%gen_shmat_8(ithr, shvec, shmat_8)
        pft_ref_8 = pft_ref_8 * shmat_8
        select case(params_glob%cc_objfun)
            case(OBJFUN_CC)
                gencorr_for_rot_8 = self%gencorr_cc_for_rot_8(     pft_ref_8, iptcl, i, irot )
            case(OBJFUN_EUCLID)
                gencorr_for_rot_8 = self%gencorr_euclid_for_rot_8( pft_ref_8, iptcl, i, irot )
        end select
    end function gencorr_for_rot_8

    function gencorr_cc_for_rot_8( self, pft_ref_8, iptcl, i, irot ) result( cc )
        class(polarft_corrcalc), intent(inout) :: self
        complex(dp), pointer,    intent(inout) :: pft_ref_8(:,:)
        integer,                 intent(in)    :: iptcl, i, irot
        real(dp) :: cc, sqsum_ref, sqsum_ptcl
        integer  :: ik
        sqsum_ref  = 0._dp
        sqsum_ptcl = 0._dp
        cc         = 0._dp
        do ik = self%kfromto(1),self%kfromto(2)
            sqsum_ref  = sqsum_ref  + real(ik,kind=dp) * sum(csq_fast(pft_ref_8(:,ik)))
            sqsum_ptcl = sqsum_ptcl + real(ik,kind=dp) * sum(csq_fast(self%pfts_ptcls(:,ik,i)))
            cc         = cc         + real(ik,kind=dp) * self%calc_corrk_for_rot_8(pft_ref_8, i, ik, irot)
        end do
        cc = cc / sqrt(sqsum_ref * sqsum_ptcl)
    end function gencorr_cc_for_rot_8

    real(dp) function gencorr_euclid_for_rot_8( self, pft_ref_8, iptcl, i, irot )
        class(polarft_corrcalc), intent(inout) :: self
        complex(dp), pointer,    intent(inout) :: pft_ref_8(:,:)
        integer,                 intent(in)    :: iptcl, i, irot
        real(dp) :: denom
        call self%weight_ref_ptcl(pft_ref_8, i, iptcl)
        denom                    = sum(real(csq_fast(self%pfts_ptcls(:, self%kfromto(1):self%kfromto(2),i)), dp))
        gencorr_euclid_for_rot_8 = dexp( - self%calc_euclid_for_rot_8(pft_ref_8, iptcl, i, irot)/denom )
        call self%deweight_ref_ptcl(pft_ref_8, i, iptcl)
    end function gencorr_euclid_for_rot_8

    subroutine gencorr_grad_for_rot_8( self, iref, iptcl, shvec, irot, f, grad )
        class(polarft_corrcalc), intent(inout) :: self
        integer,                 intent(in)    :: iref, iptcl
        real(dp),                intent(in)    :: shvec(2)
        integer,                 intent(in)    :: irot
        real(dp),                intent(out)   :: f, grad(2)
        complex(dp), pointer :: pft_ref_8(:,:), shmat_8(:,:), pft_ref_tmp(:,:)
        integer              :: ithr, i
        call self%prep_ref4corr(iref, iptcl, pft_ref_8, i, ithr)
        pft_ref_tmp => self%heap_vars(ithr)%pft_ref_tmp_8
        shmat_8     => self%heap_vars(ithr)%shmat_8
        call self%gen_shmat_8(ithr, shvec, shmat_8)
        pft_ref_8 = pft_ref_8 * shmat_8
        select case(params_glob%cc_objfun)
            case(OBJFUN_CC)
                call self%gencorr_cc_grad_for_rot_8(     pft_ref_8, pft_ref_tmp, iptcl, i, irot, f, grad )
            case(OBJFUN_EUCLID)
                call self%gencorr_euclid_grad_for_rot_8( pft_ref_8, pft_ref_tmp, iptcl, i, irot, f, grad )
        end select
    end subroutine gencorr_grad_for_rot_8

    subroutine gencorr_cc_grad_for_rot_8( self, pft_ref, pft_ref_tmp, iptcl, i, irot, f, grad )
        class(polarft_corrcalc), intent(inout) :: self
        complex(dp), pointer,    intent(inout) :: pft_ref(:,:), pft_ref_tmp(:,:)
        integer,                 intent(in)    :: iptcl, i, irot
        real(dp),                intent(out)   :: f, grad(2)
        real(dp) :: denom, sqsum_ref, sqsum_ptcl
        integer  :: ik
        ! use jacobian resolution weights
        sqsum_ref   = 0._dp
        sqsum_ptcl  = 0._dp
        f           = 0._dp
        grad        = 0._dp
        pft_ref_tmp = pft_ref * (0.d0, 1.d0) * self%argtransf(:self%pftsz,:)
        do ik = self%kfromto(1),self%kfromto(2)
            sqsum_ref  = sqsum_ref  + real(ik,kind=dp) * sum(csq_fast(pft_ref(:,ik)))
            sqsum_ptcl = sqsum_ptcl + real(ik,kind=dp) * sum(csq_fast(self%pfts_ptcls(:,ik,i)))
            f          = f          + real(ik,kind=dp) * self%calc_corrk_for_rot_8(pft_ref,     i, ik, irot)
            grad(1)    = grad(1)    + real(ik,kind=dp) * self%calc_corrk_for_rot_8(pft_ref_tmp, i, ik, irot)
        end do
        pft_ref_tmp = pft_ref * (0.d0, 1.d0) * self%argtransf(self%pftsz + 1:,:)
        do ik = self%kfromto(1),self%kfromto(2)
            grad(2) = grad(2) + real(ik,kind=dp) * self%calc_corrk_for_rot_8(pft_ref_tmp, i, ik, irot)
        end do
        f    = f    / sqrt(sqsum_ref*sqsum_ptcl)
        grad = grad / sqrt(sqsum_ref*sqsum_ptcl)
    end subroutine gencorr_cc_grad_for_rot_8

    subroutine gencorr_euclid_grad_for_rot_8( self, pft_ref, pft_ref_tmp, iptcl, i, irot, f, grad )
        class(polarft_corrcalc), intent(inout) :: self
        complex(dp), pointer,    intent(inout) :: pft_ref(:,:), pft_ref_tmp(:,:)
        integer,                 intent(in)    :: iptcl, i, irot
        real(dp),                intent(out)   :: f, grad(2)
        real(dp) :: diffsq, denom
        integer  :: k
        call self%weight_ref_ptcl(pft_ref, i, iptcl)
        denom       = sum(real(csq_fast(self%pfts_ptcls(:, self%kfromto(1):self%kfromto(2),i)), dp))
        f           = self%calc_euclid_for_rot_8(pft_ref, iptcl, i, irot)
        grad        = 0._dp
        pft_ref_tmp = pft_ref * (0.d0, 1.d0) * self%argtransf(:self%pftsz,:)
        do k = self%kfromto(1), self%kfromto(2)
            diffsq  = real(sum(pft_ref_tmp(:,k)*conjg(pft_ref(:,k)))) - self%calc_corrk_for_rot_8(pft_ref_tmp, i, k, irot)
            grad(1) = grad(1) + diffsq
        end do
        pft_ref_tmp = pft_ref * (0.d0, 1.d0) * self%argtransf(self%pftsz + 1:,:)
        do k = self%kfromto(1), self%kfromto(2)
            diffsq  = real(sum(pft_ref_tmp(:,k)*conjg(pft_ref(:,k)))) - self%calc_corrk_for_rot_8(pft_ref_tmp, i, k, irot)
            grad(2) = grad(2) + diffsq
        end do
        f    = dexp( -f/denom )
        grad = -f * 2._dp * grad/denom
        call self%deweight_ref_ptcl(pft_ref, i, iptcl)
    end subroutine gencorr_euclid_grad_for_rot_8

    subroutine gencorr_grad_only_for_rot_8( self, iref, iptcl, shvec, irot, grad )
        class(polarft_corrcalc), intent(inout) :: self
        integer,                 intent(in)    :: iref, iptcl
        real(dp),                intent(in)    :: shvec(2)
        integer,                 intent(in)    :: irot
        real(dp),                intent(out)   :: grad(2)
        complex(dp), pointer :: pft_ref_8(:,:), shmat_8(:,:), pft_ref_tmp(:,:)
        integer              :: ithr, i
        call self%prep_ref4corr(iref, iptcl, pft_ref_8, i, ithr)
        pft_ref_tmp => self%heap_vars(ithr)%pft_ref_tmp_8
        shmat_8     => self%heap_vars(ithr)%shmat_8
        call self%gen_shmat_8(ithr, shvec, shmat_8)
        pft_ref_8 = pft_ref_8 * shmat_8
        select case(params_glob%cc_objfun)
            case(OBJFUN_CC)
                call self%gencorr_cc_grad_only_for_rot_8(     pft_ref_8, pft_ref_tmp, iptcl, i, irot, grad )
            case(OBJFUN_EUCLID)
                call self%gencorr_euclid_grad_only_for_rot_8( pft_ref_8, pft_ref_tmp, iptcl, i, irot, grad )
        end select
    end subroutine gencorr_grad_only_for_rot_8

    subroutine gencorr_cc_grad_only_for_rot_8( self, pft_ref, pft_ref_tmp, iptcl, i, irot, grad )
        class(polarft_corrcalc), intent(inout) :: self
        complex(dp), pointer,    intent(inout) :: pft_ref(:,:), pft_ref_tmp(:,:)
        integer,                 intent(in)    :: iptcl, i, irot
        real(dp),                intent(out)   :: grad(2)
        real(dp) :: denom, corr
        integer  :: ik
        real(dp) :: sqsum_ref, sqsum_ptcl
        sqsum_ref  = 0._dp
        sqsum_ptcl = 0._dp
        grad       = 0._dp
        pft_ref_tmp = pft_ref * (0.d0, 1.d0) * self%argtransf(:self%pftsz,:)
        do ik = self%kfromto(1),self%kfromto(2)
            sqsum_ref  = sqsum_ref  + real(ik,kind=dp) * sum(csq_fast(pft_ref(:,ik)))
            sqsum_ptcl = sqsum_ptcl + real(ik,kind=dp) * sum(csq_fast(self%pfts_ptcls(:,ik,i)))
            grad(1)    = grad(1)    + real(ik,kind=dp) * self%calc_corrk_for_rot_8(pft_ref_tmp, i, ik, irot)
        end do
        pft_ref_tmp = pft_ref * (0.d0, 1.d0) * self%argtransf(self%pftsz + 1:,:)
        do ik = self%kfromto(1),self%kfromto(2)
            grad(2) = grad(2) + real(ik,kind=dp) * self%calc_corrk_for_rot_8(pft_ref_tmp, i, ik, irot)
        end do
        grad = grad / sqrt(sqsum_ref*sqsum_ptcl)
    end subroutine gencorr_cc_grad_only_for_rot_8

    subroutine gencorr_euclid_grad_only_for_rot_8( self, pft_ref, pft_ref_tmp, iptcl, i, irot, grad )
        class(polarft_corrcalc), intent(inout) :: self
        complex(dp), pointer,    intent(inout) :: pft_ref(:,:), pft_ref_tmp(:,:)
        integer,                 intent(in)    :: iptcl, i, irot
        real(dp),                intent(out)   :: grad(2)
        real(dp) :: f
        call self%gencorr_euclid_grad_for_rot_8(pft_ref, pft_ref_tmp, iptcl, i, irot, f, grad)
    end subroutine gencorr_euclid_grad_only_for_rot_8

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
            fdf_y(k) = self%calc_corrk_for_rot_8(pft_ref, self%pinds(iptcl), k, irot)
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
            fdf_y(k) = self%calc_corrk_for_rot_8(pft_ref, self%pinds(iptcl), k, irot)
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
            fdf_y(k) = self%calc_corrk_for_rot_8(pft_ref, self%pinds(iptcl), k, irot)
        end do
        call self%calc_T1_T2_for_rot_8( pft_ref, pft_dref, iref, irot, 3, fdf_T1, fdf_T2)
        f = sum(fdf_y) / denom
        do j = 1,3
            grad(j) = ( sum(fdf_T1(:,j)) - sum(fdf_y) * sum(fdf_T2(:,j)) / sqsum_ref ) / denom
        end do
        pft_ref_tmp = pft_ref * (0.d0, 1.d0) * self%argtransf(:self%pftsz,:)
        corr        = self%calc_corr_for_rot_8(pft_ref_tmp, self%pinds(iptcl), irot)
        grad(4)     = corr / denom
        pft_ref_tmp = pft_ref * (0.d0, 1.d0) * self%argtransf(self%pftsz + 1:,:)
        corr        = self%calc_corr_for_rot_8(pft_ref_tmp, self%pinds(iptcl), irot)
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
        integer  :: i, ithr, k
        call self%prep_ref4corr(iref, iptcl, pft_ref, i, ithr)
        shmat => self%heap_vars(ithr)%shmat
        call self%gen_shmat(ithr, shvec, shmat)
        pft_ref = pft_ref * shmat
        do k = self%kfromto(1), self%kfromto(2)
            sigma_contrib(k) = 0.5 * self%calc_euclidk_for_rot(pft_ref, i, k, irot) / real(self%pftsz)
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
        call self%genfrc(iref, iptcl, irot, frc)
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

    subroutine weight_ref_ptcl_sp( self, pft_ref, i, iptcl )
        class(polarft_corrcalc), intent(inout) :: self
        complex(sp),    pointer, intent(inout) :: pft_ref(:,:)
        integer,                 intent(in)    :: i, iptcl
        integer  :: k
        real(sp) :: w
        do k=self%kfromto(1),self%kfromto(2)
            w                         = sqrt(self%npix_per_shell(k) / real(self%pftsz) / (2. * self%sigma2_noise(k,iptcl)))
            pft_ref(:,k)              = w * pft_ref(:,k)
            self%pfts_ptcls(:,k,i)    = w * self%pfts_ptcls(:,k,i)
            self%fftdat_ptcls(i,k)%re = w * self%fftdat_ptcls(i,k)%re
            self%fftdat_ptcls(i,k)%im = w * self%fftdat_ptcls(i,k)%im
        end do
    end subroutine weight_ref_ptcl_sp

    subroutine weight_ref_ptcl_dp( self, pft_ref, i, iptcl )
        class(polarft_corrcalc), intent(inout) :: self
        complex(dp),    pointer, intent(inout) :: pft_ref(:,:)
        integer,                 intent(in)    :: i, iptcl
        integer  :: k
        real(dp) :: w
        do k=self%kfromto(1),self%kfromto(2)
            w                         = dsqrt(self%npix_per_shell(k) / real(self%pftsz, dp) / (2.0_dp * self%sigma2_noise(k,iptcl)))
            pft_ref(:,k)              = w * pft_ref(:,k)
            self%pfts_ptcls(:,k,i)    = w * self%pfts_ptcls(:,k,i)
            self%fftdat_ptcls(i,k)%re = w * self%fftdat_ptcls(i,k)%re
            self%fftdat_ptcls(i,k)%im = w * self%fftdat_ptcls(i,k)%im
        end do
    end subroutine weight_ref_ptcl_dp

    subroutine deweight_ref_ptcl_sp( self, pft_ref, i, iptcl )
        class(polarft_corrcalc), intent(inout) :: self
        complex(sp),    pointer, intent(inout) :: pft_ref(:,:)
        integer,                 intent(in)    :: i, iptcl
        integer  :: k
        real(sp) :: w
        do k=self%kfromto(1),self%kfromto(2)
            w                         = sqrt(self%npix_per_shell(k) / real(self%pftsz) / (2. * self%sigma2_noise(k,iptcl)))
            pft_ref(:,k)              = pft_ref(:,k)              / w
            self%pfts_ptcls(:,k,i)    = self%pfts_ptcls(:,k,i)    / w
            self%fftdat_ptcls(i,k)%re = self%fftdat_ptcls(i,k)%re / w
            self%fftdat_ptcls(i,k)%im = self%fftdat_ptcls(i,k)%im / w
        end do
    end subroutine deweight_ref_ptcl_sp

    subroutine deweight_ref_ptcl_dp( self, pft_ref, i, iptcl )
        class(polarft_corrcalc), intent(inout) :: self
        complex(dp),    pointer, intent(inout) :: pft_ref(:,:)
        integer,                 intent(in)    :: i, iptcl
        integer  :: k
        real(dp) :: w
        do k=self%kfromto(1),self%kfromto(2)
            w                         = dsqrt(self%npix_per_shell(k) / real(self%pftsz, dp) / (2.0_dp * self%sigma2_noise(k,iptcl)))
            pft_ref(:,k)              = pft_ref(:,k)              / w
            self%pfts_ptcls(:,k,i)    = self%pfts_ptcls(:,k,i)    / w
            self%fftdat_ptcls(i,k)%re = self%fftdat_ptcls(i,k)%re / w
            self%fftdat_ptcls(i,k)%im = self%fftdat_ptcls(i,k)%im / w
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
                    call fftwf_free(self%fftdat_ptcls(i,ik)%p_re)
                    call fftwf_free(self%fftdat_ptcls(i,ik)%p_im)
                end do
            end do
            if( allocated(self%ctfmats)        ) deallocate(self%ctfmats)
            if( allocated(self%npix_per_shell) ) deallocate(self%npix_per_shell)
            deallocate( self%sqsums_ptcls, self%angtab, self%argtransf,&
                &self%polar, self%pfts_refs_even, self%pfts_refs_odd, self%pfts_drefs_even, self%pfts_drefs_odd,&
                self%pfts_ptcls, self%fft_factors, self%fftdat, self%fftdat_ptcls, self%fft_carray,&
                &self%iseven, self%pinds, self%heap_vars, self%argtransf_shellone)
            call fftwf_destroy_plan(self%plan_bwd)
            call fftwf_destroy_plan(self%plan_fwd_1)
            call fftwf_destroy_plan(self%plan_fwd_2)
            nullify(self%sigma2_noise, pftcc_glob)
            self%l_filt_set   = .false.
            self%existence    = .false.
        endif
    end subroutine kill

end module simple_polarft_corrcalc
