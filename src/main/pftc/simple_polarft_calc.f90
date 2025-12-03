! for calculation of band-pass limited cross-correlation of polar Fourier transforms

! For consideration:

! 1. Memoization + thread pools should be isolated and reusable
! Right now heap_vars has many duplicated arrays. Many can be:
! allocated once, reused across calls
! replaced with a scratch allocator
! stored in a thread-local pool to avoid repeated malloc/free overhead
! Refactoring isolates memo routines → easy to replace with better memory management.

! 2. FFT plan creation can move into allocate_refs_memoization / constructor
! Instead of creating FFTW plans multiple times or destroying them manually, encapsulate plan lifetimes cleanly:
! create plan on first call
! reuse across all operations (FFT reuse = major speedup)
! destroy only in kill
! A dedicated FFT submodule makes this easy.

! 4. Correlation routines are highly redundant
! You have dozens of near-identical kernels:
! gen_corrs_cc
! gen_corrs_shifted_cc
! gen_corrs_weighted_cc
! gen_corrs_shifted_weighted_cc
! gen_corr_for_rot_8_*
! calc_corr_rot_shift
! many more
! These can be unified into:
! shared template functions
! macro-generated variants (if you want)
! or a strategy pattern (function pointers)
! After refactor, these functions live together → easier to unify performance-critical kernels.

! 5. Opportunity for GPU offload
! If/when you consider GPU acceleration, the new submodules naturally isolate the FFT and kernel operations.
! The correlation kernels in particular are highly parallelizable.

module simple_polarft_calc
!$ use omp_lib
!$ use omp_lib_kinds
include 'simple_lib.f08'
use simple_sp_project, only: sp_project
use simple_parameters, only: params_glob
implicit none

public :: polarft_calc, pftc_glob
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
end type heap_vars

type :: polarft_calc
    private
    integer                  :: nptcls     = 1           !< the total number of particles in partition (logically indexded [fromp,top])
    integer                  :: nrefs      = 1           !< the number of references (logically indexded [1,nrefs])
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
    type(c_ptr)              :: plan_fwd1, plan_bwd1
    type(c_ptr)              :: plan_mem_r2c
    ! Memoized terms
    complex(kind=c_float_complex), allocatable :: ft_ptcl_ctf(:,:,:)                      !< Fourier Transform of particle times CTF
    complex(kind=c_float_complex), allocatable :: ft_absptcl_ctf(:,:,:)                   !< Fourier Transform of (particle times CTF)**2
    complex(kind=c_float_complex), allocatable :: ft_ctf2(:,:,:)                          !< Fourier Transform of CTF squared modulus
    complex(kind=c_float_complex), allocatable :: ft_ref_even(:,:,:),  ft_ref_odd(:,:,:)  !< Fourier Transform of even/odd references
    complex(kind=c_float_complex), allocatable :: ft_ref2_even(:,:,:), ft_ref2_odd(:,:,:) !< Fourier Transform of even/odd references squared modulus
    ! Convenience vectors, thread memoization
    type(heap_vars),               allocatable :: heap_vars(:)
    type(fftw_cvec),               allocatable :: cvec1(:), cvec2(:)
    type(fftw_rvec),               allocatable :: rvec1(:)
    type(fftw_drvec),              allocatable :: drvec(:)
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
    procedure          :: get_work_pft_ptr
    procedure          :: get_work_rpft_ptr
    procedure          :: get_work_rpft8_ptr
    procedure          :: get_ptcls_ptr
    procedure          :: get_ctfmats_ptr
    procedure          :: get_refs_ptr
    ! ===== VIS (print, vis_ptcl, vis_ref, read/write): simple_polarft_vis.f90
    procedure          :: vis_ptcl
    procedure          :: vis_ref
    procedure          :: print
    ! ===== CTF: simple_polarft_ctf.f90
    procedure          :: create_polar_absctfmats, calc_polar_ctf
    ! ===== GEOM: simple_polarft_geom.f90
    procedure          :: gen_shmat
    procedure, private :: gen_shmat_8
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
    procedure          :: calc_corr_rot_shift, calc_magcorr_rot
    procedure          :: gen_corrs_mag, gen_corrs_mag_cc
    procedure          :: gen_corrs_weighted_cc, gen_corrs_shifted_weighted_cc
    procedure          :: gen_corrs_cc,          gen_corrs_shifted_cc
    procedure, private :: gen_corrs_1, gen_corrs_2
    generic            :: gen_corrs => gen_corrs_1, gen_corrs_2
    procedure, private :: gen_corr_for_rot_8_1, gen_corr_for_rot_8_2
    generic            :: gen_corr_for_rot_8 => gen_corr_for_rot_8_1, gen_corr_for_rot_8_2
    procedure          :: gen_corr_grad_for_rot_8
    procedure          :: gen_corr_grad_only_for_rot_8
    procedure          :: gen_corr_cc_for_rot_8
    procedure          :: gen_corr_cc_grad_for_rot_8
    procedure          :: gen_corr_cc_grad_only_for_rot_8
    procedure          :: gen_euclids
    procedure          :: gen_euclids_shifted
    procedure          :: gen_euclid_for_rot_8
    procedure          :: gen_euclid_grad_for_rot_8
    procedure          :: gen_sigma_contrib
    procedure          :: calc_frc
    procedure          :: bidirectional_shift_search  
end type polarft_calc

interface

    ! ===== CORE (new, kill, setters, getters, pointer helpers) =====

    module subroutine new(self, nrefs, pfromto, kfromto, eoarr)
        class(polarft_calc), target, intent(inout) :: self
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
        complex(sp),         intent(in)    :: pft(:,:)
        logical,             intent(in)    :: iseven
    end subroutine set_ref_pft

    module subroutine set_ptcl_pft(self, iptcl, pft)
        class(polarft_calc), intent(inout) :: self
        integer,             intent(in)    :: iptcl
        complex(sp),         intent(in)    :: pft(:,:)
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

    module subroutine get_ptcls_ptr( self, ptr )
        class(polarft_calc), target,  intent(in)  :: self
        complex(sp),         pointer, intent(out) :: ptr(:,:,:)
    end subroutine get_ptcls_ptr

    module subroutine get_ctfmats_ptr( self, ptr )
        class(polarft_calc), target,  intent(in)  :: self
        real(sp),            pointer, intent(out) :: ptr(:,:,:)
    end subroutine get_ctfmats_ptr

    module subroutine get_refs_ptr( self, ptre, ptro )
        class(polarft_calc), target,  intent(in)  :: self
        complex(sp),         pointer, intent(out) :: ptre(:,:,:), ptro(:,:,:)
    end subroutine get_refs_ptr

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

    ! module subroutine read_ptcls_ctfs(self, ibatch, success)
    !     class(polarft_calc), intent(inout) :: self
    !     integer,             intent(in)    :: ibatch
    !     logical,             intent(inout) :: success
    ! end subroutine read_ptcls_ctfs

    ! module subroutine write_ptcls_ctfs(self, ibatch)
    !     class(polarft_calc), intent(in) :: self
    !     integer,             intent(in) :: ibatch
    ! end subroutine write_ptcls_ctfs

    module subroutine print(self)
        class(polarft_calc), intent(in) :: self
    end subroutine print

    ! ===== CTF =====

    module subroutine calc_polar_ctf(self, iptcl, smpd, kv, cs, fraca, dfx, dfy, angast)
        class(polarft_calc), intent(inout) :: self
        integer,             intent(in)    :: iptcl
        real,                intent(in)    :: smpd, kv, cs, fraca, dfx, dfy, angast
    end subroutine calc_polar_ctf

    module subroutine create_polar_absctfmats(self, spproj, oritype, pfromto)
        class(polarft_calc),       intent(inout) :: self
        class(sp_project), target, intent(inout) :: spproj
        character(len=*),          intent(in) :: oritype
        integer, optional,         intent(in) :: pfromto(2)
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
        integer,              intent(in) :: ithr
        real(dp),             intent(in) :: shift_8(2)
        complex(dp), pointer, intent(inout) :: shmat_8(:,:)
    end subroutine gen_shmat_8

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
        class(polarft_calc), intent(in) :: self
        complex(dp),         intent(in) :: pft(:,:)
        integer,             intent(in) :: irot
        complex(dp),         intent(out) :: pft_rot(:,:)
    end subroutine rotate_pft_1

    module subroutine rotate_pft_2(self, pft, irot, pft_rot)
        class(polarft_calc), intent(in) :: self
        complex(sp),         intent(in) :: pft(:,:)
        integer,             intent(in) :: irot
        complex(sp),         intent(out) :: pft_rot(:,:)
    end subroutine rotate_pft_2

    module subroutine rotate_pft_3(self, pft, irot, pft_rot)
        class(polarft_calc), intent(in) :: self
        real(sp),            intent(in) :: pft(:,:)
        integer,             intent(in) :: irot
        real(sp),            intent(out) :: pft_rot(:,:)
    end subroutine rotate_pft_3

    module subroutine rotate_pft_4(self, pft, irot, pft_rot)
        class(polarft_calc), intent(in)  :: self
        real(dp),            intent(in)  :: pft(:,:)
        integer,             intent(in)  :: irot
        real(dp),            intent(out) :: pft_rot(:,:)
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
        integer,             intent(in) :: iref, iptcl, irot
        real(sp),            intent(in) :: shvec(2)
        logical, optional,   intent(in) :: kweight
        real :: val
    end function calc_corr_rot_shift

    module function calc_magcorr_rot(self, iref, iptcl, irot, kweight) result(val)
        class(polarft_calc), intent(inout) :: self
        integer,             intent(in) :: iref, iptcl, irot
        logical, optional,   intent(in) :: kweight
        real :: val
    end function calc_magcorr_rot

    module subroutine gen_corrs_mag_cc( self, iref, iptcl, ccs, kweight )
        class(polarft_calc), intent(inout) :: self
        integer,             intent(in)    :: iref, iptcl
        real,                intent(inout) :: ccs(self%pftsz)
        logical,   optional, intent(in)    :: kweight
    end subroutine gen_corrs_mag_cc

    module subroutine gen_corrs_mag(self, iref, iptcl, ccs, kweight)
        class(polarft_calc), intent(inout) :: self
        integer,             intent(in) :: iref, iptcl
        real,                intent(inout) :: ccs(self%pftsz)
        logical, optional,   intent(in) :: kweight
    end subroutine gen_corrs_mag

    module subroutine calc_frc(self, iref, iptcl, irot, shvec, frc)
        class(polarft_calc), intent(inout) :: self
        integer,             intent(in) :: iref, iptcl, irot
        real(sp),            intent(in) :: shvec(2)
        real(sp),            intent(out) :: frc(self%kfromto(1):self%kfromto(2))
    end subroutine calc_frc

    module subroutine gen_corrs_1(self, iref, iptcl, cc, kweight)
        class(polarft_calc), intent(inout) :: self
        integer,             intent(in) :: iref, iptcl
        real(sp),            intent(out) :: cc(self%nrots)
        logical, optional,   intent(in) :: kweight
    end subroutine gen_corrs_1

    module subroutine gen_corrs_2(self, iref, iptcl, shift, cc, kweight)
        class(polarft_calc), intent(inout) :: self
        integer,             intent(in) :: iref, iptcl
        real(sp),            intent(in) :: shift(2)
        real(sp),            intent(out) :: cc(self%nrots)
        logical, optional,   intent(in) :: kweight
    end subroutine gen_corrs_2

    module subroutine gen_corrs_cc(self, iptcl, iref, corrs)
        class(polarft_calc), intent(inout) :: self
        integer,             intent(in) :: iptcl, iref
        real(sp),            intent(out) :: corrs(self%nrots)
    end subroutine gen_corrs_cc

    module subroutine gen_corrs_shifted_cc(self, pft_ref, iptcl, iref, corrs)
        class(polarft_calc), intent(inout) :: self
        complex(sp),         intent(in) :: pft_ref(1:self%pftsz,self%kfromto(1):self%kfromto(2))
        integer,             intent(in) :: iptcl, iref
        real(sp),            intent(out) :: corrs(self%nrots)
    end subroutine gen_corrs_shifted_cc

    module subroutine gen_corrs_weighted_cc(self, iptcl, iref, corrs)
        class(polarft_calc), intent(inout) :: self
        integer,             intent(in) :: iptcl, iref
        real(sp),            intent(out) :: corrs(self%nrots)
    end subroutine gen_corrs_weighted_cc

    module subroutine gen_corrs_shifted_weighted_cc(self, pft_ref, iptcl, iref, corrs)
        class(polarft_calc), intent(inout) :: self
        complex(sp),         intent(in) :: pft_ref(1:self%pftsz,self%kfromto(1):self%kfromto(2))
        integer,             intent(in) :: iptcl, iref
        real(sp),            intent(out) :: corrs(self%nrots)
    end subroutine gen_corrs_shifted_weighted_cc

    module subroutine gen_euclids(self, iptcl, iref, euclids)
        class(polarft_calc), intent(inout) :: self
        integer,             intent(in) :: iptcl, iref
        real(sp),            intent(out) :: euclids(self%nrots)
    end subroutine gen_euclids

    module subroutine gen_euclids_shifted(self, pft_ref, iptcl, iref, euclids)
        class(polarft_calc),  intent(inout) :: self
        complex(sp), pointer, intent(in) :: pft_ref(:,:)
        integer,              intent(in) :: iptcl, iref
        real(sp),             intent(out) :: euclids(self%nrots)
    end subroutine gen_euclids_shifted

    module subroutine bidirectional_shift_search(self, iref, iptcl, irot, hn, shifts, grid1, grid2)
        class(polarft_calc), intent(inout) :: self
        integer,             intent(in) :: iref, iptcl, irot, hn
        real,                intent(in) :: shifts(-hn:hn)
        real,                intent(out) :: grid1(-hn:hn,-hn:hn), grid2(-hn:hn,-hn:hn)
    end subroutine bidirectional_shift_search

    module function gen_corr_for_rot_8_1(self, iref, iptcl, irot) result(val)
        class(polarft_calc), intent(inout) :: self
        integer,             intent(in) :: iref, iptcl, irot
        real(dp) :: val
    end function gen_corr_for_rot_8_1

    module function gen_corr_for_rot_8_2(self, iref, iptcl, shvec, irot) result(val)
        class(polarft_calc), intent(inout) :: self
        integer,             intent(in) :: iref, iptcl, irot
        real(dp),            intent(in) :: shvec(2)
        real(dp) :: val
    end function gen_corr_for_rot_8_2

    module function gen_corr_cc_for_rot_8(self, pft_ref, i) result(val)
        class(polarft_calc),  intent(inout) :: self
        complex(dp), pointer, intent(inout) :: pft_ref(:,:)
        integer,              intent(in) :: i
        real(dp) :: val
    end function gen_corr_cc_for_rot_8

    module function gen_euclid_for_rot_8(self, pft_ref, iptcl) result(val)
        class(polarft_calc),  intent(inout) :: self
        complex(dp), pointer, intent(inout) :: pft_ref(:,:)
        integer,              intent(in) :: iptcl
        real(dp) :: val
    end function gen_euclid_for_rot_8

    module subroutine gen_corr_grad_for_rot_8(self, iref, iptcl, shvec, irot, f, grad)
        class(polarft_calc), intent(inout) :: self
        integer,             intent(in) :: iref, iptcl, irot
        real(dp),            intent(in) :: shvec(2)
        real(dp),            intent(out) :: f, grad(2)
    end subroutine gen_corr_grad_for_rot_8

    module subroutine gen_corr_cc_grad_for_rot_8(self, pft_ref, pft_ref_tmp, iptcl, irot, f, grad)
        class(polarft_calc),  intent(inout) :: self
        complex(dp), pointer, intent(inout) :: pft_ref(:,:), pft_ref_tmp(:,:)
        integer,              intent(in) :: iptcl, irot
        real(dp),             intent(out) :: f, grad(2)
    end subroutine gen_corr_cc_grad_for_rot_8

    module subroutine gen_euclid_grad_for_rot_8(self, pft_ref, pft_ref_tmp, iptcl, irot, f, grad)
        class(polarft_calc),  intent(inout) :: self
        complex(dp), pointer, intent(inout) :: pft_ref(:,:), pft_ref_tmp(:,:)
        integer,              intent(in) :: iptcl, irot
        real(dp),             intent(out) :: f, grad(2)
    end subroutine gen_euclid_grad_for_rot_8

    module subroutine gen_corr_grad_only_for_rot_8(self, iref, iptcl, shvec, irot, grad)
        class(polarft_calc), intent(inout) :: self
        integer,             intent(in) :: iref, iptcl, irot
        real(dp),            intent(in) :: shvec(2)
        real(dp),            intent(out) :: grad(2)
    end subroutine gen_corr_grad_only_for_rot_8

    module subroutine gen_corr_cc_grad_only_for_rot_8(self, pft_ref, pft_ref_tmp, i, irot, grad)
        class(polarft_calc),  intent(inout) :: self
        complex(dp), pointer, intent(inout) :: pft_ref(:,:), pft_ref_tmp(:,:)
        integer,              intent(in) :: i, irot
        real(dp),             intent(out) :: grad(2)
    end subroutine gen_corr_cc_grad_only_for_rot_8

    module subroutine gen_sigma_contrib(self, iref, iptcl, shvec, irot, sigma_contrib)
        class(polarft_calc), intent(inout) :: self
        integer,             intent(in) :: iref, iptcl, irot
        real(sp),            intent(in) :: shvec(2)
        real(sp), optional,  intent(out) :: sigma_contrib(self%kfromto(1):self%kfromto(2))
    end subroutine gen_sigma_contrib

end interface

! CLASS PARAMETERS/VARIABLES
complex(sp), parameter       :: zero            = cmplx(0.,0.) !< just a complex zero
integer,     parameter       :: FFTW_USE_WISDOM = 16
class(polarft_calc), pointer :: pftc_glob       => null()

end module simple_polarft_calc
