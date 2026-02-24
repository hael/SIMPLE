!@descr: polarft class complete interface
module simple_polarft_calc
use simple_pftc_api
implicit none

public :: polarft_calc, pftc_glob, polaft_dims_from_file_header
private
#include "simple_local_flags.inc"

type fftw_drvec
    type(c_ptr)                  :: p
    real(kind=c_double), pointer :: r(:) => null()
end type fftw_drvec

type fftw_cmat
    type(c_ptr)                            :: p
    complex(kind=c_float_complex), pointer :: c(:,:) => null()
end type fftw_cmat

type fftw_crmat
    type(c_ptr)                            :: p
    complex(kind=c_float_complex), pointer :: c(:,:) => null()  ! (pftsz+1, nk)
    real(kind=c_float),            pointer :: r(:,:) => null()  ! (nrots+2, nk) view on same memory
end type fftw_crmat

type heap_vars
    complex(sp), pointer :: pft_ref(:,:)       => null()
    complex(sp), pointer :: pft_ref_tmp(:,:)   => null()
    real(dp),    pointer :: argvec(:)          => null()
    complex(sp), pointer :: shmat(:,:)         => null()
    real(dp),    pointer :: kcorrs(:)          => null()
    complex(dp), pointer :: pft_ref_8(:,:)     => null()
    complex(dp), pointer :: pft_ref_tmp_8(:,:) => null()
    complex(dp), pointer :: pft_dref_8(:,:,:)  => null()
    complex(dp), pointer :: shvec(:)           => null()
    complex(dp), pointer :: shmat_8(:,:)       => null()
    real(dp),    pointer :: pft_r1_8(:,:)      => null()
    real(sp),    pointer :: pft_r(:,:)         => null()
    ! cache arrays for batched computations
    real(dp),    pointer :: w_weights(:)       => null()  ! k/sigma2 weights
    real(dp),    pointer :: sumsq_cache(:)     => null()  ! particle sumsq per k
end type heap_vars

type :: polarft_calc
    private
    class(parameters), pointer :: p_ptr => null()        !< pointer to parameters object (for access to runtime parameters)
    integer                  :: nptcls     = 1           !< the total number of particles in partition (logically indexded [fromp,top])
    integer                  :: nrefs      = 1           !< the number of references (logically indexded [1,nrefs])
    integer                  :: ncls       = 0
    integer                  :: nrots      = 0           !< number of in-plane rotations for one pft (determined by radius of molecule)
    integer                  :: pftsz      = 0           !< size of reference and particle pft (nrots/2)
    integer                  :: pfromto(2) = 0           !< particle index range
    integer                  :: ldim(3)    = 0           !< logical dimensions of original cartesian image
    integer                  :: kfromto(2)               !< band-pass Fourier index limits
    integer                  :: nk                       !< number of shells used during alignment
    real                     :: dang                     !< angular increment
    integer,     allocatable :: pinds(:)                 !< index array (to reduce memory when frac_update < 1)
    real,        allocatable :: npix_per_shell(:)        !< number of (cartesian) pixels per shell
    real(dp),    allocatable :: sqsums_ptcls(:)          !< memoized square sums for the correlation calculations (taken from kfromto(1):kfromto(2))
    real(dp),    allocatable :: ksqsums_ptcls(:)         !< memoized k-weighted square sums for the correlation calculations (taken from kfromto(1):kfromto(2))
    real(dp),    allocatable :: wsqsums_ptcls(:)         !< memoized square sums weighted by k and  sigmas^2 (taken from kfromto(1):kfromto(2))
    real(sp),    allocatable :: angtab(:)                !< table of in-plane angles (in degrees)
    real(dp),    allocatable :: argtransf(:,:)           !< argument transfer constants for shifting the references
    real(sp),    allocatable :: polar(:,:)               !< table of polar coordinates (in Cartesian coordinates)
    real(sp),    allocatable :: ctfmats(:,:,:)           !< expand set of CTF matrices (for efficient parallel exec)
    real(dp),    allocatable :: argtransf_shellone(:)    !< one dimensional argument transfer constants (shell k=1) for shifting the references
    complex(sp), allocatable :: pfts_refs_even(:,:,:)    !< 3D complex matrix of polar reference sections (pftsz,nk,nrefs), even
    complex(sp), allocatable :: pfts_refs_odd(:,:,:)     !< -"-, odd
    complex(sp), allocatable :: pfts_drefs_even(:,:,:,:) !< derivatives w.r.t. orientation angles of 3D complex matrices
    complex(sp), allocatable :: pfts_drefs_odd(:,:,:,:)  !< derivatives w.r.t. orientation angles of 3D complex matrices
    complex(sp), allocatable :: pfts_ptcls(:,:,:)        !< 3D complex matrix of particle sections
    ! FFTW plans
    ! batched FFT plans for vectors of length nrots (nk transforms)
    type(c_ptr) :: plan_fwd1_many, plan_bwd1_many, plan_mem_r2c_many
    ! Memoized terms
    complex(kind=c_float_complex), allocatable :: ft_ptcl_ctf(:,:,:)                      !< Fourier Transform of particle times CTF
    complex(kind=c_float_complex), allocatable :: ft_absptcl_ctf(:,:,:)                   !< Fourier Transform of (particle times CTF)**2
    complex(kind=c_float_complex), allocatable :: ft_ctf2(:,:,:)                          !< Fourier Transform of CTF squared modulus
    complex(kind=c_float_complex), allocatable :: ft_ref_even(:,:,:),  ft_ref_odd(:,:,:)  !< Fourier Transform of even/odd references
    complex(kind=c_float_complex), allocatable :: ft_ref2_even(:,:,:), ft_ref2_odd(:,:,:) !< Fourier Transform of even/odd references squared modulus
    ! buffer: (nrots,nk) holding cvec2 for all k
    type(fftw_cmat),  allocatable :: cmat2_many(:)  ! per thread
    ! batched backward (c2r) plan: nk transforms of length nrots, in-place
    type(fftw_crmat), allocatable :: crmat1_many(:) ! per thread
    ! Convenience vectors, thread memoization
    type(heap_vars),  allocatable :: heap_vars(:)
    type(fftw_drvec), allocatable :: drvec(:)
    ! for the subset of polarft_ops submodules
    complex(dp),      allocatable :: pfts_even(:,:,:), pfts_odd(:,:,:), pfts_merg(:,:,:)
    real(dp),         allocatable :: ctf2_even(:,:,:), ctf2_odd(:,:,:)
    integer,          allocatable :: prev_eo_pops(:,:), eo_pops(:,:)
    logical                       :: l_comlin   = .false.
    ! Others
    logical, allocatable :: iseven(:)                   !< eo assignment for gold-standard FSC
    real,    pointer     :: sigma2_noise(:,:) => null() !< for euclidean distances
    logical              :: with_ctf  = .false.         !< CTF flag
    logical              :: existence = .false.         !< to indicate existence
  contains
    ! ===== CORE (new, kill, setters, getters, pointer helpers): simple_polarft_core.f90
    procedure          :: new
    procedure          :: kill
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
    procedure          :: set_with_ctf
    procedure          :: assign_sigma2_noise
    ! ===== GETTERS + POINTER ACCESSORS: simple_polarft_access.f90
    procedure          :: get_nrots
    procedure          :: get_pdim
    procedure          :: get_kfromto
    procedure          :: get_pftsz
    procedure          :: get_rot
    procedure          :: get_roind
    procedure          :: get_roind_fast
    procedure          :: get_dang
    procedure          :: get_coord
    procedure          :: get_ref_pft
    procedure          :: get_nrefs
    procedure          :: exists
    procedure          :: ptcl_iseven
    procedure          :: get_nptcls
    procedure          :: get_pinds
    procedure          :: get_npix
    procedure          :: is_with_ctf
    procedure          :: allocate_pft
    procedure, private :: get_work_pft_ptr
    procedure, private :: get_work_rpft_ptr
    procedure, private :: get_work_rpft8_ptr
    ! ===== VIS (print, vis_ptcl, vis_ref, read/write): simple_polarft_vis.f90
    procedure          :: vis_ptcl
    procedure          :: vis_ref
    procedure          :: print
    ! ===== CTF: simple_polarft_ctf.f90
    procedure          :: create_polar_absctfmats
    ! ===== GEOM: simple_polarft_geom.f90
    procedure, private :: gen_shmat
    procedure, private :: gen_shmat_8
    procedure          :: gen_clin_weights
    procedure          :: shift_ptcl
    procedure          :: shift_ref
    procedure          :: mirror_ref_pft
    procedure, private :: rotate_pft_1, rotate_pft_2, rotate_pft_3, rotate_pft_4
    generic            :: rotate_pft => rotate_pft_1, rotate_pft_2, rotate_pft_3, rotate_pft_4
    procedure, private :: rotate_iref_1, rotate_iref_2
    generic            :: rotate_iref => rotate_iref_1, rotate_iref_2
    ! ===== MEMO: simple_polarft_memo.f90
    procedure          :: memoize_sqsum_ptcl
    procedure, private :: setup_npix_per_shell
    procedure          :: memoize_ptcls, memoize_refs
    procedure, private :: kill_memoized_ptcls, kill_memoized_refs
    procedure          :: allocate_ptcls_memoization, allocate_refs_memoization
    ! ===== CORR: simple_polarft_corr.f90
    procedure          :: calc_corr_rot_shift
    procedure          :: calc_frc
    procedure          :: gen_objfun_vals
    procedure, private :: gen_corrs
    procedure, private :: gen_euclids
    procedure, private :: gen_corr_for_rot_8_1, gen_corr_for_rot_8_2
    generic            :: gen_corr_for_rot_8 => gen_corr_for_rot_8_1, gen_corr_for_rot_8_2
    procedure          :: gen_corr_grad_for_rot_8
    procedure          :: gen_corr_grad_only_for_rot_8
    procedure, private :: gen_corr_cc_for_rot_8
    procedure, private :: gen_corr_cc_grad_for_rot_8
    procedure, private :: gen_corr_cc_grad_only_for_rot_8
    procedure, private :: gen_euclid_for_rot_8
    procedure, private :: gen_euclid_grad_for_rot_8
    procedure          :: gen_sigma_contrib

    ! ===== STATE: simple_polarft_ops_state.f90
    procedure          :: polar_cavger_new
    procedure          :: polar_cavger_zero_pft_refs
    procedure          :: polar_cavger_set_ref_pft
    procedure          :: polar_cavger_calc_pops
    procedure          :: polar_cavger_update_sums
    procedure          :: polar_cavger_kill
    procedure          :: center_3Dpolar_refs
    ! ===== RESTORE: simple_polarft_ops_restore.f90
    procedure          :: polar_cavger_merge_eos_and_norm2D
    procedure          :: polar_cavger_merge_eos_and_norm
    procedure          :: polar_cavger_calc_and_write_frcs_and_eoavg
    procedure          :: polar_prep2Dref
    procedure          :: polar_cavger_gen2Dclassdoc
    procedure          :: polar_filterrefs
    procedure, private :: calc_comlin_contrib
    procedure, private :: mirror_ctf2
    procedure, private :: mirror_pft
    procedure, private :: safe_norm
    procedure, private :: get_line
    procedure          :: polar_cavger_calc_frc
    ! ===== I/O: simple_polarft_ops_io.f90
    procedure          :: polar_cavger_refs2cartesian
    procedure          :: polar_cavger_read
    procedure          :: polar_cavger_write
    procedure          :: polar_cavger_writeall
    procedure          :: polar_cavger_writeall_pftcrefs
    procedure          :: polar_cavger_write_cartrefs
    procedure          :: polar_cavger_read_all
    procedure          :: polar_cavger_readwrite_partial_sums
    procedure          :: polar_cavger_assemble_sums_from_parts
    procedure, private :: write_pft_array_local
    procedure, private :: write_pft_array
    procedure, private :: write_ctf2_array_local
    procedure, private :: write_ctf2_array
    procedure, private :: open_pft_array_for_read
    procedure, private :: transfer_pft_array_buffer
    procedure, private :: read_pft_array
    procedure, private :: open_ctf2_array_for_read
    procedure, private :: transfer_ctf2_array_buffer
end type polarft_calc

interface

    ! ===== CORE (new, kill, setters, getters, pointer helpers) =====

   module subroutine new(self, params, nrefs, pfromto, kfromto, eoarr)
        use simple_parameters, only: parameters
        class(polarft_calc), target, intent(inout) :: self
        class(parameters),   target, intent(in)    :: params
        integer,                     intent(in)    :: nrefs
        integer,                     intent(in)    :: pfromto(2), kfromto(2)
        integer, optional,           intent(in)    :: eoarr(pfromto(1):pfromto(2))
    end subroutine new

    module subroutine kill(self)
        class(polarft_calc), intent(inout) :: self
    end subroutine kill

    module subroutine reallocate_ptcls(self, nptcls, pinds)
        class(polarft_calc), intent(inout) :: self
        integer,             intent(in)    :: nptcls
        integer,             intent(in)    :: pinds(nptcls)
    end subroutine reallocate_ptcls

    module subroutine set_ref_pft(self, iref, pft, iseven)
        class(polarft_calc), intent(inout) :: self
        integer,             intent(in)    :: iref
        complex(sp),         intent(in)    :: pft(self%pftsz,self%kfromto(1):self%kfromto(2))
        logical,             intent(in)    :: iseven
    end subroutine set_ref_pft

    module subroutine set_ptcl_pft(self, iptcl, pft)
        class(polarft_calc), intent(inout) :: self
        integer,             intent(in)    :: iptcl
        complex(sp),         intent(in)    :: pft(self%pftsz,self%kfromto(1):self%kfromto(2))
    end subroutine set_ptcl_pft

    module subroutine set_ref_fcomp(self, iref, irot, k, comp, iseven)
        class(polarft_calc), intent(inout) :: self
        integer,             intent(in)    :: iref, irot, k
        complex(sp),         intent(in)    :: comp
        logical,             intent(in)    :: iseven
    end subroutine set_ref_fcomp

    module subroutine set_dref_fcomp(self, iref, irot, k, dcomp, iseven)
        class(polarft_calc), intent(inout) :: self
        integer,             intent(in)    :: iref, irot, k
        complex(sp),         intent(in)    :: dcomp(3)
        logical,             intent(in)    :: iseven
    end subroutine set_dref_fcomp

    module subroutine set_ptcl_fcomp(self, iptcl, irot, k, comp)
        class(polarft_calc), intent(inout) :: self
        integer,             intent(in)    :: iptcl, irot, k
        complex(sp),         intent(in)    :: comp
    end subroutine set_ptcl_fcomp

    module subroutine cp_even2odd_ref(self, iref)
        class(polarft_calc), intent(inout) :: self
        integer,             intent(in)    :: iref
    end subroutine cp_even2odd_ref

    module subroutine cp_odd2even_ref(self, iref)
        class(polarft_calc), intent(inout) :: self
        integer,             intent(in)    :: iref
    end subroutine cp_odd2even_ref

    module subroutine cp_refs(self, self2)
        class(polarft_calc), intent(inout) :: self, self2
    end subroutine cp_refs

    module subroutine cp_even_ref2ptcl(self, iref, iptcl)
        class(polarft_calc), intent(inout) :: self
        integer,             intent(in)    :: iref, iptcl
    end subroutine cp_even_ref2ptcl

    module subroutine swap_ptclsevenodd(self)
        class(polarft_calc), intent(inout) :: self
    end subroutine swap_ptclsevenodd

    module subroutine set_eo(self, iptcl, is_even)
        class(polarft_calc), intent(inout) :: self
        integer,             intent(in)    :: iptcl
        logical,             intent(in)    :: is_even
    end subroutine set_eo

    module subroutine set_eos(self, eoarr)
        class(polarft_calc), intent(inout) :: self
        integer,             intent(in)    :: eoarr(self%nptcls)
    end subroutine set_eos

    module subroutine set_with_ctf(self, l_wctf)
        class(polarft_calc), intent(inout) :: self
        logical,             intent(in)    :: l_wctf
    end subroutine set_with_ctf

    module subroutine assign_sigma2_noise(self, sigma2_noise)
        class(polarft_calc),       intent(inout) :: self
        real, allocatable, target, intent(inout) :: sigma2_noise(:,:)
    end subroutine assign_sigma2_noise

    ! ===== GETTERS + POINTER ACCESSORS =====

    module pure function get_nrots(self) result(nrots)
        class(polarft_calc), intent(in) :: self
        integer :: nrots
    end function get_nrots

    module pure function get_pdim(self) result(pdim)
        class(polarft_calc), intent(in) :: self
        integer :: pdim(3)
    end function get_pdim

    module pure function get_kfromto(self) result(kfromto)
        class(polarft_calc), intent(in) :: self
        integer :: kfromto(2)
    end function get_kfromto

    module pure integer function get_pftsz(self)
        class(polarft_calc), intent(in) :: self
    end function get_pftsz

    module function get_rot(self, roind) result(rot)
        class(polarft_calc), intent(in) :: self
        integer,             intent(in) :: roind
        real(sp) :: rot
    end function get_rot

    module function get_roind(self, rot) result(ind)
        class(polarft_calc), intent(in) :: self
        real(sp),            intent(in) :: rot
        integer :: ind
    end function get_roind

    module pure integer function get_roind_fast(self, psi)
        class(polarft_calc), intent(in) :: self
        real,                intent(in) :: psi
    end function get_roind_fast

    module pure real function get_dang(self)
        class(polarft_calc), intent(in) :: self
    end function get_dang

    module pure function get_coord(self, rot, k) result(xy)
        class(polarft_calc), intent(in) :: self
        integer,             intent(in) :: rot, k
        real(sp) :: xy(2)
    end function get_coord

    module subroutine get_ref_pft(self, iref, iseven, pft)
        class(polarft_calc), intent(in)    :: self
        integer,             intent(in)    :: iref
        logical,             intent(in)    :: iseven
        complex(sp),         intent(inout) :: pft(self%pftsz,self%kfromto(1):self%kfromto(2))
    end subroutine get_ref_pft

    module pure integer function get_nrefs(self)
        class(polarft_calc), intent(in) :: self
    end function get_nrefs

    module logical function exists(self)
        class(polarft_calc), intent(in) :: self
    end function exists

    module logical function ptcl_iseven(self, iptcl)
        class(polarft_calc), intent(in) :: self
        integer, intent(in) :: iptcl
    end function ptcl_iseven

    module integer function get_nptcls(self)
        class(polarft_calc), intent(in) :: self
    end function get_nptcls

    module subroutine get_pinds(self, pinds)
        class(polarft_calc),  intent(in)  :: self
        integer, allocatable, intent(out) :: pinds(:)
    end subroutine get_pinds

    module integer function get_npix(self)
        class(polarft_calc), intent(in) :: self
    end function get_npix

    module pure logical function is_with_ctf( self )
        class(polarft_calc), intent(in) :: self
    end function is_with_ctf

    module function allocate_pft( self ) result( pft )
        class(polarft_calc),  intent(in)  :: self
        complex(sp), allocatable :: pft(:,:)
    end function allocate_pft

    module subroutine get_work_pft_ptr( self, ptr )
        class(polarft_calc),  intent(in)  :: self
        complex(sp), pointer, intent(out) :: ptr(:,:)
    end subroutine get_work_pft_ptr

    module subroutine get_work_rpft_ptr( self, ptr )
        class(polarft_calc), intent(in)  :: self
        real(sp), pointer,   intent(out) :: ptr(:,:)
    end subroutine get_work_rpft_ptr

    module subroutine get_work_rpft8_ptr( self, ptr )
        class(polarft_calc), intent(in)  :: self
        real(dp), pointer,   intent(out) :: ptr(:,:)
    end subroutine get_work_rpft8_ptr

    ! ===== VIS (print, vis_ptcl, vis_ref, read/write) =====

    module subroutine vis_ptcl(self, iptcl)
        class(polarft_calc), intent(in) :: self
        integer,             intent(in) :: iptcl
    end subroutine vis_ptcl

    module subroutine vis_ref(self, iref, iseven)
        class(polarft_calc), intent(in) :: self
        integer,             intent(in) :: iref
        logical,             intent(in) :: iseven
    end subroutine vis_ref

    module subroutine print(self)
        class(polarft_calc), intent(in) :: self
    end subroutine print

    ! ===== CTF =====

    module subroutine create_polar_absctfmats(self, spproj, oritype, pfromto)
        class(polarft_calc),       intent(inout) :: self
        class(sp_project), target, intent(inout) :: spproj
        character(len=*),          intent(in)    :: oritype
        integer, optional,         intent(in)    :: pfromto(2)
    end subroutine create_polar_absctfmats

    ! ===== GEOM (rotate_pft, rotate_iref, shift matrices) =====

    module subroutine gen_shmat(self, ithr, shift, shmat)
        class(polarft_calc),  intent(inout) :: self
        integer,              intent(in)    :: ithr
        real(sp),             intent(in)    :: shift(2)
        complex(sp), pointer, intent(inout) :: shmat(:,:)
    end subroutine gen_shmat

    module subroutine gen_shmat_8(self, ithr, shift_8, shmat_8)
        class(polarft_calc),  intent(inout) :: self
        integer,              intent(in)    :: ithr
        real(dp),             intent(in)    :: shift_8(2)
        complex(dp), pointer, intent(inout) :: shmat_8(:,:)
    end subroutine gen_shmat_8

    module subroutine gen_clin_weights( self, psi, lrot, rrot, lw, rw )
        class(polarft_calc), intent(in)    :: self
        real,                intent(in)    :: psi
        integer,             intent(inout) :: lrot
        integer,             intent(inout) :: rrot
        real(dp),            intent(out)   :: lw(self%kfromto(1):self%kfromto(2))
        real(dp),            intent(out)   :: rw(self%kfromto(1):self%kfromto(2))
    end subroutine gen_clin_weights

    module subroutine shift_ptcl( self, iptcl, shvec)
        class(polarft_calc), intent(inout) :: self
        integer,             intent(in)    :: iptcl
        real(sp),            intent(in)    :: shvec(2)
    end subroutine shift_ptcl

    module subroutine shift_ref( self, iref, shvec)
        class(polarft_calc), intent(inout) :: self
        integer,             intent(in)    :: iref
        real(sp),            intent(in)    :: shvec(2)
    end subroutine shift_ref

    module subroutine mirror_ref_pft( self, iref )
        class(polarft_calc), target, intent(in) :: self
        integer,                     intent(in) :: iref
    end subroutine mirror_ref_pft

    module subroutine rotate_pft_1(self, pft, irot, pft_rot)
        class(polarft_calc), intent(in)  :: self
        complex(dp),         intent(in)  :: pft(self%pftsz,self%kfromto(1):self%kfromto(2))
        integer,             intent(in)  :: irot
        complex(dp),         intent(out) :: pft_rot(self%pftsz,self%kfromto(1):self%kfromto(2))
    end subroutine rotate_pft_1

    module subroutine rotate_pft_2(self, pft, irot, pft_rot)
        class(polarft_calc), intent(in)  :: self
        complex(sp),         intent(in)  :: pft(self%pftsz,self%kfromto(1):self%kfromto(2))
        integer,             intent(in)  :: irot
        complex(sp),         intent(out) :: pft_rot(self%pftsz,self%kfromto(1):self%kfromto(2))
    end subroutine rotate_pft_2

    module subroutine rotate_pft_3(self, pft, irot, pft_rot)
        class(polarft_calc), intent(in)  :: self
        real(sp),            intent(in)  :: pft(self%pftsz,self%kfromto(1):self%kfromto(2))
        integer,             intent(in)  :: irot
        real(sp),            intent(out) :: pft_rot(self%pftsz,self%kfromto(1):self%kfromto(2))
    end subroutine rotate_pft_3

    module subroutine rotate_pft_4(self, pft, irot, pft_rot)
        class(polarft_calc), intent(in)  :: self
        real(dp),            intent(in)  :: pft(self%pftsz,self%kfromto(1):self%kfromto(2))
        integer,             intent(in)  :: irot
        real(dp),            intent(out) :: pft_rot(self%pftsz,self%kfromto(1):self%kfromto(2))
    end subroutine rotate_pft_4

    module subroutine rotate_iref_1(self, iref, irot, sh)
        class(polarft_calc), intent(inout) :: self
        integer,             intent(in)    :: iref, irot
        real,                intent(in)    :: sh(2)
    end subroutine rotate_iref_1

    module subroutine rotate_iref_2(self, pft_ref, irot, sh, pft_ref_out)
        class(polarft_calc),  intent(inout) :: self
        complex(dp), pointer, intent(in)    :: pft_ref(:,:)
        integer,              intent(in)    :: irot
        real,                 intent(in)    :: sh(2)
        complex(dp), pointer, intent(out)   :: pft_ref_out(:,:)
    end subroutine rotate_iref_2

    ! ===== MEMO  =====

    module subroutine setup_npix_per_shell(self)
        class(polarft_calc), intent(inout) :: self
    end subroutine setup_npix_per_shell

    module subroutine memoize_sqsum_ptcl(self, iptcl)
        class(polarft_calc), intent(inout) :: self
        integer,             intent(in) :: iptcl
    end subroutine memoize_sqsum_ptcl

    module subroutine memoize_ptcls(self)
        class(polarft_calc), intent(inout) :: self
    end subroutine memoize_ptcls

    module subroutine memoize_refs(self)
        class(polarft_calc), intent(inout) :: self
    end subroutine memoize_refs

    module subroutine allocate_ptcls_memoization(self)
        class(polarft_calc), intent(inout) :: self
    end subroutine allocate_ptcls_memoization

    module subroutine allocate_refs_memoization(self)
        class(polarft_calc), intent(inout) :: self
    end subroutine allocate_refs_memoization

    module subroutine kill_memoized_ptcls(self)
        class(polarft_calc), intent(inout) :: self
    end subroutine kill_memoized_ptcls

    module subroutine kill_memoized_refs(self)
        class(polarft_calc), intent(inout) :: self
    end subroutine kill_memoized_refs

    ! ===== CORR  =====

    module function calc_corr_rot_shift(self, iref, iptcl, shvec, irot, kweight) result(val)
        class(polarft_calc), intent(inout) :: self
        integer,             intent(in)    :: iref, iptcl, irot
        real(sp),            intent(in)    :: shvec(2)
        logical, optional,   intent(in)    :: kweight
        real :: val
    end function calc_corr_rot_shift

    module subroutine calc_frc(self, iref, iptcl, irot, shvec, frc)
        class(polarft_calc), intent(inout) :: self
        integer,             intent(in) :: iref, iptcl, irot
        real(sp),            intent(in) :: shvec(2)
        real(sp),            intent(out) :: frc(self%kfromto(1):self%kfromto(2))
    end subroutine calc_frc

     module subroutine gen_objfun_vals(self, iref, iptcl, shift, vals)
        class(polarft_calc), intent(inout) :: self
        integer,             intent(in)    :: iref, iptcl
        real(sp),            intent(in)    :: shift(2)
        real(sp),            intent(out)   :: vals(self%nrots)
    end subroutine gen_objfun_vals

    module subroutine gen_corrs(self, iref, iptcl, shift, cc)
        class(polarft_calc), target, intent(inout) :: self
        integer,                     intent(in)    :: iref, iptcl
        real(sp),                    intent(in)    :: shift(2)
        real(sp),                    intent(out)   :: cc(self%nrots)
    end subroutine gen_corrs

    module subroutine gen_euclids(self, iref, iptcl, shift, euclids)
        class(polarft_calc), target, intent(inout) :: self
        integer,                     intent(in)    :: iref, iptcl
        real,                        intent(in)    :: shift(2)
        real(sp),                    intent(out)   :: euclids(self%nrots)
    end subroutine gen_euclids

    module function gen_corr_for_rot_8_1(self, iref, iptcl, irot) result(val)
        class(polarft_calc), target, intent(inout) :: self
        integer,                     intent(in)    :: iref, iptcl, irot
        real(dp) :: val
    end function gen_corr_for_rot_8_1

    module function gen_corr_for_rot_8_2(self, iref, iptcl, shvec, irot) result(val)
        class(polarft_calc), target, intent(inout) :: self
        integer,                     intent(in)    :: iref, iptcl, irot
        real(dp),                    intent(in)    :: shvec(2)
        real(dp) :: val
    end function gen_corr_for_rot_8_2

    module function gen_corr_cc_for_rot_8(self, pft_ref, i) result(val)
        class(polarft_calc),  intent(inout) :: self
        complex(dp), pointer, intent(inout) :: pft_ref(:,:)
        integer,              intent(in)    :: i
        real(dp) :: val
    end function gen_corr_cc_for_rot_8

    module function gen_euclid_for_rot_8(self, pft_ref, iptcl) result(val)
        class(polarft_calc),  intent(inout) :: self
        complex(dp), pointer, intent(inout) :: pft_ref(:,:)
        integer,              intent(in)    :: iptcl
        real(dp) :: val
    end function gen_euclid_for_rot_8

    module subroutine gen_corr_grad_for_rot_8(self, iref, iptcl, shvec, irot, f, grad)
        class(polarft_calc), target, intent(inout) :: self
        integer,                     intent(in)    :: iref, iptcl, irot
        real(dp),                    intent(in)    :: shvec(2)
        real(dp),                    intent(out)   :: f, grad(2)
    end subroutine gen_corr_grad_for_rot_8

    module subroutine gen_corr_cc_grad_for_rot_8(self, pft_ref, pft_ref_tmp, iptcl, irot, f, grad)
        class(polarft_calc),  intent(inout) :: self
        complex(dp), pointer, intent(inout) :: pft_ref(:,:), pft_ref_tmp(:,:)
        integer,              intent(in)    :: iptcl, irot
        real(dp),             intent(out)   :: f, grad(2)
    end subroutine gen_corr_cc_grad_for_rot_8

    module subroutine gen_euclid_grad_for_rot_8(self, pft_ref, pft_ref_tmp, iptcl, irot, f, grad)
        class(polarft_calc), target, intent(inout) :: self
        complex(dp), pointer,        intent(inout) :: pft_ref(:,:), pft_ref_tmp(:,:)
        integer,                     intent(in)    :: iptcl, irot
        real(dp),                    intent(out)   :: f, grad(2)
    end subroutine gen_euclid_grad_for_rot_8

    module subroutine gen_corr_grad_only_for_rot_8(self, iref, iptcl, shvec, irot, grad)
        class(polarft_calc), target, intent(inout) :: self
        integer,                     intent(in)    :: iref, iptcl, irot
        real(dp),                    intent(in)    :: shvec(2)
        real(dp),                    intent(out)   :: grad(2)
    end subroutine gen_corr_grad_only_for_rot_8

    module subroutine gen_corr_cc_grad_only_for_rot_8(self, pft_ref, pft_ref_tmp, i, irot, grad)
        class(polarft_calc),  intent(inout) :: self
        complex(dp), pointer, intent(inout) :: pft_ref(:,:), pft_ref_tmp(:,:)
        integer,              intent(in)    :: i, irot
        real(dp),             intent(out)   :: grad(2)
    end subroutine gen_corr_cc_grad_only_for_rot_8

    module subroutine gen_sigma_contrib(self, iref, iptcl, shvec, irot, sigma_contrib)
        class(polarft_calc), target, intent(inout) :: self
        integer,                     intent(in)    :: iref, iptcl, irot
        real(sp),                    intent(in)    :: shvec(2)
        real(sp), optional,          intent(out)   :: sigma_contrib(self%kfromto(1):self%kfromto(2))
    end subroutine gen_sigma_contrib

    ! The below routines were previoulsy implemented as an independent set of submodules.
    ! This design decision was incorrect and forced the use of pointers to access encapsulated data.
    ! This is now implemented using a decorator design pattern, where we preserve the original interface
    ! while adding the necessary functionality with direct access to encapsulated data. Submodules in
    ! modern Fortran essentially captures inheritance by composition and type extension in one hit, and 
    ! it also enables trivial implementation of the decorator design pattern for extending derived types
    ! while preserving their interface. The only drawback is that we expose any decorating 
    ! data in the original class, so if we absolutely need data integrity at this level, then we use
    ! type extension or composition with an object that is properly encapsulated.

    ! ===== STATE: simple_polarft_ops_state.f90

    module subroutine polar_cavger_new( self, l_comlin, nrefs )
        class(polarft_calc), intent(inout) :: self
        logical,             intent(in)    :: l_comlin
        integer,   optional, intent(in)    :: nrefs
    end subroutine polar_cavger_new

    module subroutine polar_cavger_zero_pft_refs( self )
        class(polarft_calc), intent(inout) :: self
    end  subroutine polar_cavger_zero_pft_refs

    module subroutine polar_cavger_set_ref_pft( self, icls, which )
        class(polarft_calc), intent(inout) :: self
        integer,             intent(in)    :: icls
        character(len=*),    intent(in)    :: which
    end subroutine polar_cavger_set_ref_pft

    module subroutine polar_cavger_calc_pops( self, spproj )
        class(polarft_calc),       intent(inout) :: self
        class(sp_project), target, intent(in)    :: spproj
    end subroutine polar_cavger_calc_pops

    module subroutine polar_cavger_update_sums( self, nptcls, pinds, spproj, incr_shifts, is3D )
        class(polarft_calc),         intent(inout) :: self
        integer,                     intent(in)    :: nptcls, pinds(nptcls)
        class(sp_project),           intent(inout) :: spproj
        real,              optional, intent(in)    :: incr_shifts(2,nptcls)
        logical,           optional, intent(in)    :: is3d
    end subroutine polar_cavger_update_sums

    module subroutine polar_cavger_kill( self )
        class(polarft_calc), intent(inout) :: self
    end subroutine polar_cavger_kill

    module subroutine center_3Dpolar_refs( self, algndoc, algnrefs )
        class(polarft_calc), intent(inout) :: self
        class(oris),         intent(inout) :: algndoc
        class(oris),         intent(in)    :: algnrefs
    end subroutine center_3Dpolar_refs

    ! ===== RESTORE: simple_polarft_ops_restore.f90

    module subroutine polar_cavger_merge_eos_and_norm2D( self )
        class(polarft_calc), intent(inout) :: self
    end subroutine polar_cavger_merge_eos_and_norm2D

    module subroutine polar_cavger_merge_eos_and_norm( self, reforis, symop, cl_weight )
        class(polarft_calc),  intent(inout) :: self
        type(oris),           intent(in)    :: reforis
        type(sym),            intent(in)    :: symop
        real,       optional, intent(in)    :: cl_weight
    end subroutine polar_cavger_merge_eos_and_norm

    module subroutine polar_cavger_calc_and_write_frcs_and_eoavg( self, clsfrcs, update_frac, fname, cline )
        use simple_cmdline, only: cmdline
        use simple_class_frcs, only: class_frcs
        class(polarft_calc), intent(inout) :: self
        class(class_frcs),   intent(inout) :: clsfrcs
        real,                intent(in)    :: update_frac
        class(string),       intent(in)    :: fname
        type(cmdline),       intent(in)    :: cline
    end subroutine polar_cavger_calc_and_write_frcs_and_eoavg

    module subroutine polar_prep2Dref( self, clsfrcs, icls, gaufilt )
        use simple_class_frcs, only: class_frcs
        class(polarft_calc), intent(inout) :: self
        class(class_frcs),   intent(inout) :: clsfrcs
        integer,             intent(in)    :: icls
        logical,             intent(in)    :: gaufilt
    end subroutine polar_prep2Dref

    module subroutine polar_cavger_gen2Dclassdoc( self, spproj, clsfrcs )
        use simple_class_frcs, only: class_frcs
        class(polarft_calc),       intent(in)    :: self
        class(sp_project), target, intent(inout) :: spproj
        class(class_frcs),         intent(inout) :: clsfrcs
    end subroutine polar_cavger_gen2Dclassdoc

    module subroutine polar_filterrefs( self, icls, filter )
        class(polarft_calc), intent(inout) :: self
        integer,             intent(in)    :: icls
        real,                intent(in)    :: filter(:)
    end subroutine polar_filterrefs

    module subroutine calc_comlin_contrib( self, ref_space, symop, pfts_cl_even, pfts_cl_odd, ctf2_cl_even, ctf2_cl_odd )
        class(polarft_calc), intent(in)    :: self
        type(oris),          intent(in)    :: ref_space
        type(sym),           intent(in)    :: symop
        complex(kind=dp),    intent(inout) :: pfts_cl_even(self%pftsz,self%kfromto(1):self%kfromto(2),self%ncls)
        complex(kind=dp),    intent(inout) :: pfts_cl_odd(self%pftsz,self%kfromto(1):self%kfromto(2),self%ncls)
        real(kind=dp),       intent(inout) :: ctf2_cl_even(self%pftsz,self%kfromto(1):self%kfromto(2),self%ncls)
        real(kind=dp),       intent(inout) :: ctf2_cl_odd(self%pftsz,self%kfromto(1):self%kfromto(2),self%ncls)
    end subroutine calc_comlin_contrib

    module pure subroutine mirror_ctf2( self, ctf2in, ctf2out )
        class(polarft_calc), intent(in)    :: self
        real(dp),            intent(in)    :: ctf2in(self%pftsz,self%kfromto(1):self%kfromto(2))
        real(dp),            intent(inout) :: ctf2out(self%pftsz,self%kfromto(1):self%kfromto(2))
    end subroutine mirror_ctf2

    module pure subroutine mirror_pft( self, pftin, pftout )
        class(polarft_calc), intent(in)    :: self
        complex(dp),         intent(in)    :: pftin(self%pftsz,self%kfromto(1):self%kfromto(2))
        complex(dp),         intent(inout) :: pftout(self%pftsz,self%kfromto(1):self%kfromto(2))
    end subroutine mirror_pft

    module pure subroutine safe_norm( self, Mnum, Mdenom, Mout )
        class(polarft_calc), intent(in)    :: self
        complex(dp),         intent(in)    :: Mnum(self%pftsz,self%kfromto(1):self%kfromto(2))
        real(dp),            intent(inout) :: Mdenom(self%pftsz,self%kfromto(1):self%kfromto(2))
        complex(dp),         intent(inout) :: Mout(self%pftsz,self%kfromto(1):self%kfromto(2))
    end subroutine safe_norm
    
    module pure subroutine get_line( self, ref, rot, even, pftline, ctf2line )
        class(polarft_calc), intent(in)  :: self
        integer,             intent(in)  :: ref, rot
        logical,             intent(in)  :: even
        complex(dp),         intent(out) :: pftline(self%kfromto(1):self%kfromto(2))
        real(dp),            intent(out) :: ctf2line(self%kfromto(1):self%kfromto(2))
    end subroutine get_line

    module subroutine polar_cavger_calc_frc( self, pft1, pft2, n, frc )
        class(polarft_calc), intent(in)    :: self
        complex(dp),         intent(in)    :: pft1(self%pftsz,self%kfromto(1):self%kfromto(2)), pft2(self%pftsz,self%kfromto(1):self%kfromto(2))
        integer,             intent(in)    :: n
        real(sp),            intent(inout) :: frc(1:n)
    end subroutine polar_cavger_calc_frc

    ! ===== I/O: simple_polarft_ops_io.f90

    module subroutine polar_cavger_refs2cartesian( self, cavgs, which, pfts_in )
        class(polarft_calc),     intent(in)    :: self
        type(image),             intent(inout) :: cavgs(self%ncls)
        character(len=*),        intent(in)    :: which
        complex(dp),   optional, intent(in)    :: pfts_in(1:self%pftsz,self%kfromto(1):self%kfromto(2),1:self%ncls)
    end subroutine polar_cavger_refs2cartesian

    module subroutine polar_cavger_read( self, fname, which )
        class(polarft_calc), intent(inout) :: self
        class(string),       intent(in)    :: fname
        character(len=*),    intent(in)    :: which
    end subroutine polar_cavger_read

    module subroutine polar_cavger_write( self, fname, which )
        class(polarft_calc), intent(in) :: self
        class(string),       intent(in) :: fname
        character(len=*),    intent(in) :: which
    end subroutine polar_cavger_write

    module subroutine polar_cavger_writeall( self, tmpl_fname )
        class(polarft_calc), intent(inout) :: self
        class(string),       intent(in)    :: tmpl_fname
    end subroutine polar_cavger_writeall

    module subroutine polar_cavger_writeall_pftcrefs( self, tmpl_fname )
        class(polarft_calc), intent(inout) :: self
        class(string),       intent(in)    :: tmpl_fname
    end subroutine polar_cavger_writeall_pftcrefs

    module subroutine polar_cavger_write_cartrefs( self, tmpl_fname, which )
        class(polarft_calc), intent(in) :: self
        class(string),       intent(in) :: tmpl_fname
        character(len=*),    intent(in) :: which
    end subroutine polar_cavger_write_cartrefs

    module subroutine polar_cavger_read_all( self, fname )
        class(polarft_calc), intent(inout) :: self
        class(string),       intent(in)    :: fname
    end subroutine polar_cavger_read_all

    module subroutine polar_cavger_readwrite_partial_sums( self, which )
        class(polarft_calc), intent(inout) :: self
        character(len=*),    intent(in)    :: which
    end subroutine polar_cavger_readwrite_partial_sums

    module subroutine polar_cavger_assemble_sums_from_parts( self, reforis, symop, clin_anneal )
        class(polarft_calc),  intent(inout) :: self
        type(oris), optional, intent(in)    :: reforis
        type(sym),  optional, intent(in)    :: symop
        real,       optional, intent(in)    :: clin_anneal
    end subroutine polar_cavger_assemble_sums_from_parts

    module subroutine write_pft_array_local( self, funit, array )
        class(polarft_calc), intent(in) :: self
        integer,             intent(in) :: funit
        complex(dp),         intent(in) :: array(self%pftsz,self%kfromto(1):self%kfromto(2),self%ncls)
    end subroutine write_pft_array_local

    module subroutine write_pft_array( self, array, fname )
        class(polarft_calc), intent(in) :: self
        complex(dp),         intent(in) :: array(self%pftsz,self%kfromto(1):self%kfromto(2),self%ncls)
        class(string),       intent(in) :: fname
    end subroutine write_pft_array

    module subroutine write_ctf2_array_local( self, funit, array )
        class(polarft_calc), intent(in) :: self
        integer,             intent(in) :: funit
        real(dp),            intent(in) :: array(self%pftsz,self%kfromto(1):self%kfromto(2),self%ncls)
    end subroutine write_ctf2_array_local

    module subroutine write_ctf2_array( self, array, fname )
        class(polarft_calc), intent(in) :: self
        real(dp),            intent(in) :: array(self%pftsz,self%kfromto(1):self%kfromto(2),self%ncls)
        class(string),       intent(in) :: fname
    end subroutine write_ctf2_array

    module subroutine open_pft_array_for_read( self, fname, array, funit, dims, buffer )
        class(polarft_calc),      intent(in)    :: self
        class(string),            intent(in)    :: fname
        complex(dp), allocatable, intent(inout) :: array(:,:,:)
        integer,                  intent(out)   :: funit, dims(4)
        complex(sp), allocatable, intent(inout) :: buffer(:,:,:)
    end subroutine open_pft_array_for_read

    module subroutine transfer_pft_array_buffer( self, array, funit, dims, buffer )
        class(polarft_calc), intent(in)    :: self
        complex(dp),         intent(inout) :: array(self%pftsz,self%kfromto(1):self%kfromto(2),self%ncls)
        integer,             intent(in)    :: funit, dims(4)
        complex(sp),         intent(inout) :: buffer(dims(1),dims(2):dims(3),dims(4))
    end subroutine transfer_pft_array_buffer
    
    module subroutine read_pft_array( self, fname, array)
        class(polarft_calc),      intent(inout) :: self
        class(string),            intent(in)    :: fname
        complex(dp), allocatable, intent(inout) :: array(:,:,:)
    end subroutine read_pft_array

    module subroutine open_ctf2_array_for_read( self, fname, array, funit, dims, buffer )
        class(polarft_calc),   intent(in)    :: self
        class(string),         intent(in)    :: fname
        real(dp), allocatable, intent(inout) :: array(:,:,:)
        integer,               intent(out)   :: funit, dims(4)
        real(sp), allocatable, intent(inout) :: buffer(:,:,:)
    end subroutine open_ctf2_array_for_read

    module subroutine transfer_ctf2_array_buffer( self, array, funit, dims, buffer )
        class(polarft_calc), intent(in) :: self
        real(dp), intent(inout) :: array(self%pftsz,self%kfromto(1):self%kfromto(2),self%ncls)
        integer,  intent(in)    :: funit, dims(4)
        real(sp), intent(inout) :: buffer(dims(1),dims(2):dims(3),dims(4))
    end subroutine transfer_ctf2_array_buffer

end interface

! CLASS PARAMETERS/VARIABLES
complex(sp), parameter       :: zero            = cmplx(0.,0.) !< just a complex zero
integer,     parameter       :: FFTW_USE_WISDOM = 16
class(polarft_calc), pointer :: pftc_glob       => null()

contains

    ! PUBLIC UTILITIES

    subroutine polaft_dims_from_file_header( fname, pftsz_here, kfromto_here, ncls_here )
        class(string), intent(in)    :: fname
        integer,       intent(inout) :: pftsz_here, kfromto_here(2), ncls_here
        integer :: dims(4), funit, io_stat
        if( .not.file_exists(fname) ) THROW_HARD(fname%to_char()//' does not exist')
        call fopen(funit, fname, access='STREAM', action='READ', status='OLD', iostat=io_stat)
        call fileiochk('dims_from_header; fopen failed: '//fname%to_char(), io_stat)
        read(unit=funit,pos=1) dims
        call fclose(funit)
        pftsz_here   = dims(1)
        kfromto_here = dims(2:3)
        ncls_here    = dims(4)
    end subroutine polaft_dims_from_file_header

end module simple_polarft_calc
