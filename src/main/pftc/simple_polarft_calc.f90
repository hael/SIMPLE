!@descr: polarft class complete interface
module simple_polarft_calc
use simple_pftc_api
implicit none

public :: polarft_calc, polarft_dims_from_file_header, polarft_estimate_lplim3D
public :: vol_pad2ref_pfts, vol_pad2ref_pfts_write_range
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
    complex(sp), pointer :: pft_ref(:,:)        => null()
    real(dp),    pointer :: argvec(:)           => null()
    complex(sp), pointer :: shmat(:,:)          => null()
    real(dp),    pointer :: kcorrs(:)           => null()
    complex(dp), pointer :: pft_ref_8(:,:)      => null()
    complex(dp), pointer :: pft_ref_tmp_8(:,:)  => null()
    complex(dp), pointer :: pft_ref_tmp2_8(:,:) => null()
    complex(dp), pointer :: shvec(:)            => null()
    complex(dp), pointer :: shmat_8(:,:)        => null()
    real(dp),    pointer :: pft_r1_8(:,:)       => null()
    real(sp),    pointer :: pft_r(:,:)          => null()
    ! cache arrays for batched computations
    real(dp),    pointer :: w_weights(:)        => null()  ! k/sigma2 weights
    real(dp),    pointer :: sumsq_cache(:)      => null()  ! particle sumsq per k
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
    integer                  :: interpklim               !< highest shell index used for interpolation
    real                     :: dang                     !< angular increment
    integer,     allocatable :: pinds(:)                 !< index array (to reduce memory when frac_update < 1)
    real(dp),    allocatable :: sqsums_ptcls(:)          !< memoized square sums for the correlation calculations (taken from kfromto(1):kfromto(2))
    real(dp),    allocatable :: ksqsums_ptcls(:)         !< memoized k-weighted square sums for the correlation calculations (taken from kfromto(1):kfromto(2))
    real(dp),    allocatable :: wsqsums_ptcls(:)         !< memoized square sums weighted by k and  sigmas^2 (taken from kfromto(1):kfromto(2))
    real(sp),    allocatable :: angtab(:)                !< table of in-plane angles (in degrees)
    real(dp),    allocatable :: argtransf(:,:)           !< argument transfer constants for shifting the references
    real(sp),    allocatable :: polar(:,:,:)             !< table of polar coordinates (2,k,irot): polar(1,k,irot)=h, polar(2,k,irot)=kc
    real(sp),    allocatable :: ctfmats(:,:,:)           !< expand set of CTF matrices (for efficient parallel exec)
    real(dp),    allocatable :: argtransf_shellone(:)    !< one dimensional argument transfer constants (shell k=1) for shifting the references
    complex(sp), allocatable :: pfts_refs_even(:,:,:)    !< 3D complex matrix of polar reference sections (pftsz,nk,nrefs), even
    complex(sp), allocatable :: pfts_refs_odd(:,:,:)     !< -"-, odd
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
    procedure          :: set_ptcl_fcomp
    procedure          :: cp_even2odd_ref
    procedure          :: cp_odd2even_ref
    procedure          :: cp_even_ref2ptcl
    procedure          :: cp_refs
    procedure          :: swap_ptclsevenodd
    procedure          :: set_eo
    procedure          :: set_with_ctf
    procedure          :: assign_sigma2_noise
    ! ===== GETTERS + POINTER ACCESSORS: simple_polarft_access.f90
    procedure          :: get_nrots
    procedure          :: get_pdim_interp
    procedure          :: get_pdim_srch
    procedure          :: get_kfromto
    procedure          :: get_pftsz
    procedure          :: get_rot
    procedure          :: get_roind
    procedure          :: get_roind_fast
    procedure          :: get_dang
    procedure          :: get_coord
    procedure          :: get_ref_pft
    procedure          :: get_ptcl_pft
    procedure          :: get_ptcl_line
    procedure          :: get_nrefs
    procedure          :: exists
    procedure          :: ptcl_iseven
    procedure          :: get_nptcls
    procedure          :: get_pinds
    procedure          :: is_with_ctf
    procedure          :: allocate_pft
    procedure          :: allocate_ptcl_pft
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
    procedure, private :: rotate_ref_8
    procedure, private :: rotate_ptcl
    procedure, private :: rotate_ctf
    ! ===== MEMO: simple_polarft_memo.f90
    procedure          :: memoize_sqsum_ptcl
    procedure          :: memoize_ptcls, memoize_refs
    procedure, private :: kill_memo_ptcls, kill_memo_refs
    procedure, private :: allocate_memo_workspace, kill_memo_workspace
    procedure, private :: alloc_memo_ptcls, alloc_memo_refs
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
    procedure          :: gen_objfun_vals_mirr_vals
    procedure, private :: gen_corrs_mirr_corrs, gen_euclids_mirr_euclids
    ! ===== STATE: simple_polarft_ops_state.f90
    procedure          :: polar_cavger_new
    procedure          :: polar_cavger_zero_pft_refs
    procedure          :: polar_cavger_set_ref_pft
    procedure          :: polar_cavger_calc_pops
    procedure          :: polar_cavger_update_sums
    procedure          :: polar_cavger_insert_ptcls_direct
    procedure          :: polar_cavger_insert_ptcls_obsfield
    procedure          :: polar_cavger_kill
    procedure          :: center_3Dpolar_refs
    ! ===== RESTORE: simple_polarft_ops_restore.f90
    procedure          :: polar_filterrefs
    procedure, private :: polar_cavger_calc_frc
    procedure          :: polar_cavger_merge_eos_and_norm
    procedure          :: polar_cavger_merge_eos_and_norm_direct
    procedure, private :: finalize_trail_rec
    procedure, private :: mirror_slices
    procedure, private :: calc_fsc
    procedure, private :: add_invtausq2rho
    procedure, private :: restore_references
    procedure, private :: calc_comlin_contrib
    procedure, private :: mirror_ctf2
    procedure, private :: mirror_pft
    procedure, private :: shell_floor_norm
    ! ===== I/O: simple_polarft_ops_io.f90
    procedure          :: polar_cavger_refs2cartesian
    procedure          :: polar_cavger_read
    procedure          :: polar_cavger_write
    procedure          :: polar_cavger_writeall
    procedure          :: polar_cavger_write_eo_pftcrefs
    procedure          :: write_ptcl_pft_range
    procedure          :: polar_cavger_read_all
    procedure          :: polar_cavger_readwrite_partial_sums
    procedure          :: polar_cavger_assemble_sums_from_parts
    procedure, private :: write_pft_array_local
    procedure, private :: write_pft_array
    procedure, private :: write_ctf2_array_local
    procedure, private :: write_ctf2_array
    procedure, private :: get_pft_array_dims
    procedure, private :: open_pft_array_for_read
    procedure, private :: transfer_pft_array_buffer
    procedure, private :: read_pft_array
    procedure, private :: read_any_pft_array
    procedure, private :: open_ctf2_array_for_read
    procedure, private :: transfer_ctf2_array_buffer
    procedure          :: assemble_projected_refs_from_parts
    procedure, private :: read_ref_pfts_range
end type polarft_calc

interface

    ! ===== CORE (new, kill, setters, getters, pointer helpers) =====

   module subroutine new(self, params, nrefs, pfromto, kfromto)
        use simple_parameters, only: parameters
        class(polarft_calc), target, intent(inout) :: self
        class(parameters),   target, intent(in)    :: params
        integer,                     intent(in)    :: nrefs
        integer,                     intent(in)    :: pfromto(2), kfromto(2)
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
        complex(sp),         intent(in)    :: pft(self%pftsz,self%kfromto(1):self%interpklim)
        logical,             intent(in)    :: iseven
    end subroutine set_ref_pft

    module subroutine set_ptcl_pft(self, iptcl, pft)
        class(polarft_calc), intent(inout) :: self
        integer,             intent(in)    :: iptcl
        complex(sp),         intent(in)    :: pft(self%pftsz,self%kfromto(1):self%interpklim)
    end subroutine set_ptcl_pft

    module pure subroutine set_ref_fcomp(self, iref, irot, k, comp, iseven)
        class(polarft_calc), intent(inout) :: self
        integer,             intent(in)    :: iref, irot, k
        complex(sp),         intent(in)    :: comp
        logical,             intent(in)    :: iseven
    end subroutine set_ref_fcomp

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

    module subroutine set_with_ctf(self, l_wctf)
        class(polarft_calc), intent(inout) :: self
        logical,             intent(in)    :: l_wctf
    end subroutine set_with_ctf

    module subroutine assign_sigma2_noise(self, sigma2_noise)
        class(polarft_calc),       intent(inout) :: self
        real, allocatable, target, intent(inout) :: sigma2_noise(:,:)
    end subroutine assign_sigma2_noise

    module subroutine vol_pad2ref_pfts(self, vol_pad, eulspace, state, iseven, mask)
        use simple_projector, only: projector
        class(polarft_calc), intent(inout) :: self
        class(projector),    intent(in)    :: vol_pad
        class(oris),         intent(inout) :: eulspace
        integer,             intent(in)    :: state
        logical,             intent(in)    :: iseven
        logical,             intent(in)    :: mask(:)
    end subroutine vol_pad2ref_pfts

    module subroutine vol_pad2ref_pfts_write_range(self, vol_pad, eulspace, state, iproj_from, iproj_to, mask, tmpl_fname)
        use simple_projector, only: projector
        class(polarft_calc), intent(inout) :: self
        class(projector),    intent(in)    :: vol_pad
        class(oris),         intent(inout) :: eulspace
        integer,             intent(in)    :: state, iproj_from, iproj_to
        logical,             intent(in)    :: mask(:)
        class(string),       intent(in)    :: tmpl_fname
    end subroutine vol_pad2ref_pfts_write_range

    ! ===== GETTERS + POINTER ACCESSORS =====

    module pure function get_nrots(self) result(nrots)
        class(polarft_calc), intent(in) :: self
        integer :: nrots
    end function get_nrots

    module pure function get_pdim_srch(self) result(pdim)
        class(polarft_calc), intent(in) :: self
        integer :: pdim(3)
    end function get_pdim_srch

    module pure function get_pdim_interp(self) result(dims)
        class(polarft_calc), intent(in) :: self
        integer :: dims(3)
    end function get_pdim_interp

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

    module subroutine get_ptcl_pft(self, iptcl, pft)
        class(polarft_calc), intent(in)    :: self
        integer,             intent(in)    :: iptcl
        complex(sp),         intent(inout) :: pft(self%pftsz,self%kfromto(1):self%kfromto(2))
    end subroutine get_ptcl_pft

    module subroutine get_ptcl_line(self, iptcl, irot, line)
        class(polarft_calc), intent(in)    :: self
        integer,             intent(in)    :: iptcl, irot
        complex(sp),         intent(inout) :: line(self%kfromto(1):self%kfromto(2))
    end subroutine get_ptcl_line

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

    module pure logical function is_with_ctf( self )
        class(polarft_calc), intent(in) :: self
    end function is_with_ctf

    module function allocate_pft( self ) result( pft )
        class(polarft_calc),  intent(in)  :: self
        complex(sp), allocatable :: pft(:,:)
    end function allocate_pft

    module function allocate_ptcl_pft( self ) result( pft )
        class(polarft_calc),  intent(in)  :: self
        complex(sp), allocatable :: pft(:,:)
    end function allocate_ptcl_pft

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

    ! ===== GEOM (rotation & shift of matrices) =====

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

    module pure subroutine gen_clin_weights( self, psi, lrot, rrot, lw, rw )
        class(polarft_calc), intent(in)    :: self
        real(dp),            intent(in)    :: psi
        integer,             intent(inout) :: lrot
        integer,             intent(inout) :: rrot
        real(dp),            intent(out)   :: lw(self%kfromto(1):self%interpklim)
        real(dp),            intent(out)   :: rw(self%kfromto(1):self%interpklim)
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

    module pure subroutine mirror_ref_pft( self, pft, pftm )
        class(polarft_calc), intent(in)  :: self
        complex(sp),         intent(in)  :: pft(1:self%pftsz, self%kfromto(1):self%interpklim)
        complex(sp),         intent(out) :: pftm(1:self%pftsz, self%kfromto(1):self%interpklim)
    end subroutine mirror_ref_pft

    module subroutine rotate_ref_8(self, pft, irot, pft_rot)
        class(polarft_calc), intent(in)  :: self
        complex(dp),         intent(in)  :: pft(self%pftsz,self%kfromto(1):self%kfromto(2))
        integer,             intent(in)  :: irot
        complex(dp),         intent(out) :: pft_rot(self%pftsz,self%kfromto(1):self%kfromto(2))
    end subroutine rotate_ref_8

    module subroutine rotate_ptcl(self, i, irot, pft_rot)
        class(polarft_calc), intent(in)  :: self
        integer,             intent(in)  :: i, irot
        complex(sp),         intent(out) :: pft_rot(self%pftsz,self%kfromto(1):self%interpklim)
    end subroutine rotate_ptcl

    module subroutine rotate_ctf(self, i, irot, ctf_rot)
        class(polarft_calc), intent(in)  :: self
        integer,             intent(in)  :: i, irot
        real(sp),            intent(out) :: ctf_rot(self%pftsz,self%kfromto(1):self%interpklim)
    end subroutine rotate_ctf

    ! ===== MEMO  =====

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

    module subroutine alloc_memo_ptcls(self)
        class(polarft_calc), intent(inout) :: self
    end subroutine alloc_memo_ptcls

    module subroutine alloc_memo_refs(self)
        class(polarft_calc), intent(inout) :: self
    end subroutine alloc_memo_refs

    module subroutine allocate_memo_workspace(self)
        class(polarft_calc), intent(inout) :: self
    end subroutine allocate_memo_workspace

    module subroutine kill_memo_ptcls(self)
        class(polarft_calc), intent(inout) :: self
    end subroutine kill_memo_ptcls

    module subroutine kill_memo_refs(self)
        class(polarft_calc), intent(inout) :: self
    end subroutine kill_memo_refs

    module subroutine kill_memo_workspace(self)
        class(polarft_calc), intent(inout) :: self
    end subroutine kill_memo_workspace

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

    module subroutine gen_objfun_vals_mirr_vals( self, iref, iptcl, vals, mvals )
        class(polarft_calc), intent(inout) :: self
        integer,             intent(in)    :: iref, iptcl
        real(sp),            intent(out)   :: vals(self%nrots), mvals(self%nrots)
    end subroutine gen_objfun_vals_mirr_vals

    module subroutine gen_corrs_mirr_corrs( self, iref, iptcl, ccs, mccs )
        class(polarft_calc), target, intent(inout) :: self
        integer,                     intent(in)    :: iref, iptcl
        real(sp),                    intent(out)   :: ccs(self%nrots), mccs(self%nrots)
    end subroutine gen_corrs_mirr_corrs

    module subroutine gen_euclids_mirr_euclids( self, iref, iptcl, euclids, meuclids )
        class(polarft_calc), target, intent(inout) :: self
        integer,                     intent(in)    :: iref, iptcl
        real(sp),                    intent(out)   :: euclids(self%nrots), meuclids(self%nrots)
    end subroutine gen_euclids_mirr_euclids

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

    module subroutine polar_cavger_insert_ptcls_direct( self, eulspace, ptcl_field, symop, nptcls, pinds, fpls )
        use simple_math_ft, only: fplane_get_cmplx, fplane_get_ctfsq
        class(polarft_calc),        intent(inout) :: self
        class(oris),                intent(in)    :: eulspace
        class(oris),  pointer,      intent(inout) :: ptcl_field
        class(sym),                 intent(in)    :: symop
        integer,                    intent(in)    :: nptcls, pinds(nptcls)
        class(fplane_type), target, intent(inout) :: fpls(nptcls)
    end subroutine polar_cavger_insert_ptcls_direct

    module subroutine polar_cavger_insert_ptcls_obsfield( self, eulspace, ptcl_field, symop, nptcls, pinds, fpls )
        class(polarft_calc),        intent(inout) :: self
        class(oris),                intent(inout) :: eulspace
        class(oris),  pointer,      intent(inout) :: ptcl_field
        class(sym),                 intent(inout) :: symop
        integer,                    intent(in)    :: nptcls, pinds(nptcls)
        class(fplane_type), target, intent(inout) :: fpls(nptcls)
    end subroutine polar_cavger_insert_ptcls_obsfield

    module subroutine polar_cavger_kill( self )
        class(polarft_calc), intent(inout) :: self
    end subroutine polar_cavger_kill

    module subroutine center_3Dpolar_refs( self, algndoc, algnrefs )
        class(polarft_calc), intent(inout) :: self
        class(oris),         intent(inout) :: algndoc
        class(oris),         intent(in)    :: algnrefs
    end subroutine center_3Dpolar_refs

    ! ===== RESTORE: simple_polarft_ops_restore.f90

    module subroutine polar_filterrefs( self, icls, filter )
        class(polarft_calc), intent(inout) :: self
        integer,             intent(in)    :: icls
        real,                intent(in)    :: filter(:)
    end subroutine polar_filterrefs

    module subroutine polar_cavger_calc_frc( self, pft1, pft2, n, frc )
        class(polarft_calc), intent(in)    :: self
        complex(dp),         intent(in)    :: pft1(self%pftsz,self%kfromto(1):self%interpklim), pft2(self%pftsz,self%kfromto(1):self%interpklim)
        integer,             intent(in)    :: n
        real(sp),            intent(inout) :: frc(1:n)
    end subroutine polar_cavger_calc_frc

    module subroutine polar_cavger_merge_eos_and_norm( self, reforis, symop, cline, update_frac )
        class(polarft_calc), intent(inout) :: self
        type(oris),          intent(in)    :: reforis
        type(sym),           intent(in)    :: symop
        type(cmdline),       intent(in)    :: cline
        real,                intent(in)    :: update_frac
    end subroutine polar_cavger_merge_eos_and_norm

    module subroutine polar_cavger_merge_eos_and_norm_direct( self, reforis, cline, update_frac )
        class(polarft_calc), intent(inout) :: self
        type(oris),          intent(in)    :: reforis
        type(cmdline),       intent(in)    :: cline
        real,                intent(in)    :: update_frac
    end subroutine polar_cavger_merge_eos_and_norm_direct

    module subroutine finalize_trail_rec( self, ufrac_trec, prev_even, prev_odd )
        class(polarft_calc),      intent(inout) :: self
        real(dp),                 intent(in)    :: ufrac_trec
        complex(dp), allocatable, intent(inout) :: prev_even(:,:,:)
        complex(dp), allocatable, intent(inout) :: prev_odd(:,:,:)
    end subroutine finalize_trail_rec

    module subroutine mirror_slices( self, ref_space )
        class(polarft_calc), intent(inout) :: self
        type(oris),          intent(in)    :: ref_space
    end subroutine mirror_slices

    module subroutine calc_fsc( self, pfts_even, pfts_odd, ctf2_even, ctf2_odd, fsc, ufrac_trec, prev_even, prev_odd )
        class(polarft_calc),   intent(in)  :: self
        complex(dp),           intent(in)  :: pfts_even(self%pftsz,self%kfromto(1):self%interpklim,self%ncls)
        complex(dp),           intent(in)  :: pfts_odd(self%pftsz,self%kfromto(1):self%interpklim,self%ncls)
        real(dp),              intent(in)  :: ctf2_even(self%pftsz,self%kfromto(1):self%interpklim,self%ncls)
        real(dp),              intent(in)  :: ctf2_odd(self%pftsz,self%kfromto(1):self%interpklim,self%ncls)
        real(dp),              intent(out) :: fsc(self%kfromto(1):self%interpklim)
        real(dp),    optional, intent(in)  :: ufrac_trec
        complex(dp), optional, intent(in)  :: prev_even(:,:,:)
        complex(dp), optional, intent(in)  :: prev_odd(:,:,:)
    end subroutine calc_fsc

    module subroutine add_invtausq2rho( self, ctf2_even, ctf2_odd, fsc )
        class(polarft_calc), intent(inout) :: self
        real(dp),            intent(inout) :: ctf2_even(self%pftsz,self%kfromto(1):self%interpklim,self%ncls)
        real(dp),            intent(inout) :: ctf2_odd(self%pftsz,self%kfromto(1):self%interpklim,self%ncls)
        real(dp),            intent(in)    :: fsc(self%kfromto(1):self%interpklim)
    end subroutine add_invtausq2rho

    module subroutine restore_references( self, reforis, pfts_even, pfts_odd, ctf2_even, ctf2_odd )
        class(polarft_calc), intent(inout) :: self
        type(oris),          intent(in)    :: reforis
        complex(dp),         intent(in)    :: pfts_even(self%pftsz,self%kfromto(1):self%interpklim,self%ncls)
        complex(dp),         intent(in)    :: pfts_odd(self%pftsz,self%kfromto(1):self%interpklim,self%ncls)
        real(dp),            intent(in)    :: ctf2_even(self%pftsz,self%kfromto(1):self%interpklim,self%ncls)
        real(dp),            intent(in)    :: ctf2_odd(self%pftsz,self%kfromto(1):self%interpklim,self%ncls)
    end subroutine restore_references

    module subroutine calc_comlin_contrib( self, ref_space, symop, pfts_cl_even, pfts_cl_odd, ctf2_cl_even, ctf2_cl_odd )
        class(polarft_calc), intent(in)    :: self
        type(oris),          intent(in)    :: ref_space
        type(sym),           intent(in)    :: symop
        complex(kind=dp),    intent(inout) :: pfts_cl_even(self%pftsz,self%kfromto(1):self%interpklim,self%ncls)
        complex(kind=dp),    intent(inout) :: pfts_cl_odd(self%pftsz,self%kfromto(1):self%interpklim,self%ncls)
        real(kind=dp),       intent(inout) :: ctf2_cl_even(self%pftsz,self%kfromto(1):self%interpklim,self%ncls)
        real(kind=dp),       intent(inout) :: ctf2_cl_odd(self%pftsz,self%kfromto(1):self%interpklim,self%ncls)
    end subroutine calc_comlin_contrib

    module pure subroutine mirror_ctf2( self, ctf2in, ctf2out )
        class(polarft_calc), intent(in)    :: self
        real(dp),            intent(in)    :: ctf2in(self%pftsz,self%kfromto(1):self%interpklim)
        real(dp),            intent(inout) :: ctf2out(self%pftsz,self%kfromto(1):self%interpklim)
    end subroutine mirror_ctf2

    module pure subroutine mirror_pft( self, pftin, pftout )
        class(polarft_calc), intent(in)    :: self
        complex(dp),         intent(in)    :: pftin(self%pftsz,self%kfromto(1):self%interpklim)
        complex(dp),         intent(inout) :: pftout(self%pftsz,self%kfromto(1):self%interpklim)
    end subroutine mirror_pft

    module pure subroutine shell_floor_norm( self, Mnum, Mdenom, Mout )
        class(polarft_calc), intent(in)    :: self
        complex(dp),         intent(in)    :: Mnum(self%pftsz,self%kfromto(1):self%interpklim)
        real(dp),            intent(inout) :: Mdenom(self%pftsz,self%kfromto(1):self%interpklim)
        complex(dp),         intent(inout) :: Mout(self%pftsz,self%kfromto(1):self%interpklim)
    end subroutine shell_floor_norm

    ! ===== I/O: simple_polarft_ops_io.f90

    module subroutine polar_cavger_refs2cartesian( self, cavgs, which )
        class(polarft_calc),     intent(in)    :: self
        type(image),             intent(inout) :: cavgs(self%ncls)
        character(len=*),        intent(in)    :: which
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

    module subroutine polar_cavger_write_eo_pftcrefs( self, tmpl_fname )
        class(polarft_calc), intent(inout) :: self
        class(string),       intent(in)    :: tmpl_fname
    end subroutine polar_cavger_write_eo_pftcrefs

    module subroutine write_ptcl_pft_range( self, fname, nptcls_total, iwrite_from, iwrite_to )
        class(polarft_calc), intent(in) :: self
        class(string),       intent(in) :: fname
        integer,             intent(in) :: nptcls_total
        integer,             intent(in) :: iwrite_from, iwrite_to
    end subroutine write_ptcl_pft_range

    module subroutine polar_cavger_read_all( self, fname )
        class(polarft_calc), intent(inout) :: self
        class(string),       intent(in)    :: fname
    end subroutine polar_cavger_read_all

    module subroutine polar_cavger_readwrite_partial_sums( self, which )
        class(polarft_calc), intent(inout) :: self
        character(len=*),    intent(in)    :: which
    end subroutine polar_cavger_readwrite_partial_sums

    module subroutine polar_cavger_assemble_sums_from_parts( self )
        class(polarft_calc),  intent(inout) :: self
    end subroutine polar_cavger_assemble_sums_from_parts

    module subroutine write_pft_array_local( self, funit, array )
        class(polarft_calc), intent(in) :: self
        integer,             intent(in) :: funit
        complex(dp),         intent(in) :: array(self%pftsz,self%kfromto(1):self%interpklim,self%ncls)
    end subroutine write_pft_array_local

    module subroutine write_pft_array( self, array, fname )
        class(polarft_calc), intent(in) :: self
        complex(dp),         intent(in) :: array(self%pftsz,self%kfromto(1):self%interpklim,self%ncls)
        class(string),       intent(in) :: fname
    end subroutine write_pft_array

    module subroutine write_ctf2_array_local( self, funit, array )
        class(polarft_calc), intent(in) :: self
        integer,             intent(in) :: funit
        real(dp),            intent(in) :: array(self%pftsz,self%kfromto(1):self%interpklim,self%ncls)
    end subroutine write_ctf2_array_local

    module subroutine write_ctf2_array( self, array, fname )
        class(polarft_calc), intent(in) :: self
        real(dp),            intent(in) :: array(self%pftsz,self%kfromto(1):self%interpklim,self%ncls)
        class(string),       intent(in) :: fname
    end subroutine write_ctf2_array

    module subroutine get_pft_array_dims( self, fname, pftsz, kfromto, nrefs )
        class(polarft_calc),      intent(in)  :: self
        class(string),            intent(in)  :: fname
        integer,                  intent(out) :: pftsz, kfromto(2), nrefs
    end subroutine get_pft_array_dims

    module subroutine open_pft_array_for_read( self, fname, array, funit, dims, buffer )
        class(polarft_calc),      intent(in)    :: self
        class(string),            intent(in)    :: fname
        complex(dp), allocatable, intent(inout) :: array(:,:,:)
        integer,                  intent(out)   :: funit, dims(4)
        complex(sp), allocatable, intent(inout) :: buffer(:,:,:)
    end subroutine open_pft_array_for_read

    module subroutine transfer_pft_array_buffer( self, array, funit, dims, buffer )
        class(polarft_calc), intent(in)    :: self
        complex(dp),         intent(inout) :: array(self%pftsz,self%kfromto(1):self%interpklim,self%ncls)
        integer,             intent(in)    :: funit, dims(4)
        complex(sp),         intent(inout) :: buffer(dims(1),dims(2):dims(3),dims(4))
    end subroutine transfer_pft_array_buffer
    
    module subroutine read_pft_array( self, fname, array)
        class(polarft_calc),      intent(in)    :: self
        class(string),            intent(in)    :: fname
        complex(dp), allocatable, intent(inout) :: array(:,:,:)
    end subroutine read_pft_array

    module subroutine read_any_pft_array( self, fname, array )
        class(polarft_calc),      intent(in)    :: self
        class(string),            intent(in)    :: fname
        complex(sp), allocatable, intent(inout) :: array(:,:,:)
    end subroutine read_any_pft_array

    module subroutine open_ctf2_array_for_read( self, fname, array, funit, dims, buffer )
        class(polarft_calc),   intent(in)    :: self
        class(string),         intent(in)    :: fname
        real(dp), allocatable, intent(inout) :: array(:,:,:)
        integer,               intent(out)   :: funit, dims(4)
        real(sp), allocatable, intent(inout) :: buffer(:,:,:)
    end subroutine open_ctf2_array_for_read

    module subroutine transfer_ctf2_array_buffer( self, array, funit, dims, buffer )
        class(polarft_calc), intent(in) :: self
        real(dp), intent(inout) :: array(self%pftsz,self%kfromto(1):self%interpklim,self%ncls)
        integer,  intent(in)    :: funit, dims(4)
        real(sp), intent(inout) :: buffer(dims(1),dims(2):dims(3),dims(4))
    end subroutine transfer_ctf2_array_buffer

    module subroutine assemble_projected_refs_from_parts( self, nparts, numlen )
        class(polarft_calc), intent(inout) :: self
        integer,             intent(in)    :: nparts, numlen
    end subroutine assemble_projected_refs_from_parts

    module subroutine read_ref_pfts_range( self, fname, iseven, iref_from, iref_to )
        class(polarft_calc), intent(inout) :: self
        class(string),       intent(in)    :: fname
        logical,             intent(in)    :: iseven
        integer,             intent(in)    :: iref_from, iref_to
    end subroutine read_ref_pfts_range

end interface

! CLASS PARAMETERS/VARIABLES
complex(sp), parameter       :: zero            = cmplx(0.,0.) !< just a complex zero
integer,     parameter       :: FFTW_USE_WISDOM = 16

contains

    ! PUBLIC UTILITIES

    ! To obtain dimensions of a PFT array from file
    subroutine polarft_dims_from_file_header( fname, pftsz_here, kfromto_here, ncls_here )
        class(string), intent(in)    :: fname
        integer,       intent(inout) :: pftsz_here, kfromto_here(2), ncls_here
        integer :: dims(4), funit, io_stat
        if( .not.file_exists(fname) ) THROW_HARD(fname%to_char()//' does not exist')
        call fopen(funit, fname, access='STREAM', action='READ', status='OLD', iostat=io_stat)
        call fileiochk('dims_from_header; fopen failed: '//fname%to_char(), io_stat)
        read(unit=funit,pos=1) dims
        call fclose(funit)
        pftsz_here   = dims(1)
        ! On-disk contract: dims(2:3) store [kfromto(1), interpklim] at write time.
        kfromto_here = dims(2:3)
        ncls_here    = dims(4)
    end subroutine polarft_dims_from_file_header

    ! To estimate resolution limit
    subroutine polarft_estimate_lplim3D( box, smpd, lpstart, lpstop, lpopt )
        use simple_butterworth, only: butterworth_filter
        integer,      intent(in)  :: box
        real,         intent(in)  :: smpd, lpstart, lpstop
        real,         intent(out) :: lpopt
        type(polarft_calc)        :: pftc ! local  & not instanciated
        type(string)              :: str_even, str_odd
        complex(dp),  allocatable :: even(:,:,:), odd(:,:,:), diff(:,:)
        real(sp),     allocatable :: bwfilter(:)
        real(dp),     allocatable :: costs(:)
        real(dp) :: cost, best_cost, wk
        integer  :: kreq(2), ksearch(2), kavail_on_disk(2), best_k, k, kk, pftsz_on_disk, ncls_on_disk, kfallback
        logical, parameter :: DEBUG = .false.
        601 format(A,1X,F12.3)
        602 format(A,1X,I8)
        603 format(A,1X,2I8)
        604 format(A,1X,ES16.8)
        if( lpstart < lpstop )then
            THROW_HARD('Invalid low-pass range ordering in polarft_estimate_lplim3D: lpstart must be >= lpstop')
        endif
        ! read current references
        str_even = POLAR_REFS_FBODY//'_even'//BIN_EXT
        str_odd  = POLAR_REFS_FBODY//'_odd'//BIN_EXT
        if( .not.file_exists(str_even) .or. .not.file_exists(str_odd) )then
            if( DEBUG )then
                write(logfhandle,'(A)') '>>> DEBUG: polarft_estimate_lplim3D'
                write(logfhandle,'(A)') '    Missing even/odd PFT files, returning lpstart limit'
            endif
            ! no update to existing limit
            lpopt = calc_lowpass_lim(calc_fourier_index(lpstart, box, smpd), box, smpd)
            return
        endif
        ! Read on-disk PFT dimensions first.
        ! Note: kavail_on_disk(2) is the stored interpklim (not the requested search kend)
        call polarft_dims_from_file_header(str_even, pftsz_on_disk, kavail_on_disk, ncls_on_disk)
        pftc%pftsz      = pftsz_on_disk
        pftc%kfromto    = kavail_on_disk
        pftc%interpklim = kavail_on_disk(2)
        pftc%ncls       = ncls_on_disk
        ! Determine requested search range from lpstart/lpstop
        kreq(1) = calc_fourier_index(lpstart, box, smpd)
        kreq(2) = calc_fourier_index(lpstop, box, smpd)
        if( kreq(1) > kreq(2) )then
            THROW_HARD('Invalid requested k-range in polarft_estimate_lplim3D: kstart must be <= kstop')
        endif
        if( DEBUG )then
            write(logfhandle,'(A)') '>>> DEBUG: polarft_estimate_lplim3D'
            write(logfhandle,602) '    box                 =', box
            write(logfhandle,601) '    smpd                =', smpd
            write(logfhandle,601) '    lpstart [A]         =', lpstart
            write(logfhandle,601) '    lpstop  [A]         =', lpstop
            write(logfhandle,603) '    requested ksearch   =', kreq(1), kreq(2)
            write(logfhandle,603) '    k available on disk =', kavail_on_disk(1), kavail_on_disk(2)
            write(logfhandle,602) '    pftc%pftsz          =', pftc%pftsz
            write(logfhandle,602) '    pftc%ncls           =', pftc%ncls
        endif
        if( kreq(2) < kavail_on_disk(1) .or. kreq(1) > kavail_on_disk(2) )then
            ! No overlap between requested search interval and on-disk available k-range.
            if( kreq(1) > kavail_on_disk(2) )then
                kfallback = kavail_on_disk(2)
            else
                kfallback = kavail_on_disk(1)
            endif
            if( DEBUG )then
                write(logfhandle,'(A)') '    No overlap between requested range and disk range, using nearest available shell'
                write(logfhandle,602) '    fallback k          =', kfallback
            endif
            lpopt = calc_lowpass_lim(kfallback, box, smpd)
            return
        endif
        ! Clamp requested interval to what is available on disk
        ksearch(1) = max(kreq(1), kavail_on_disk(1))
        ksearch(2) = min(kreq(2), kavail_on_disk(2))
        if( DEBUG )then
            if( (kreq(1) /= ksearch(1)) .or. (kreq(2) /= ksearch(2)) )then
                write(logfhandle,603) '    capped ksearch      =', ksearch(1), ksearch(2)
            endif
            write(logfhandle,603) '    optimization ksearch=', ksearch(1), ksearch(2)
        endif
        if( ksearch(1) == ksearch(2) )then
            if( DEBUG ) write(logfhandle,'(A)') '    Degenerate search range, returning current cutoff'
            ! nothing to search, return limit from current ksearch(1)
            lpopt = calc_lowpass_lim(ksearch(1), box, smpd)
            return
        endif
        ! Read PFT arrays
        call pftc%read_pft_array(str_even, even)
        call pftc%read_pft_array(str_odd,  odd)
        ! Optimization: find best resolution cutoff
        if( DEBUG )then
            write(logfhandle,'(A)') '    --- optimization trace (k, cost, is_new_best) ---'
        endif
        ! allocations (bwfilter and diff are allocated per OMP thread inside the parallel loop)
        allocate(costs(ksearch(1):ksearch(2)))
        !$omp parallel do default(shared) schedule(dynamic,1) &
        !$omp private(k, cost, bwfilter, diff, kk, wk)
        do k = ksearch(1), ksearch(2) ! search range
            ! Per-thread private workspace
            allocate(bwfilter(1:pftc%kfromto(2)), diff(pftc%pftsz,pftc%ncls))
            cost = 0.d0
            call butterworth_filter(k, bwfilter)
            do kk = pftc%kfromto(1),pftc%kfromto(2) ! resolution range available
                ! diff = even - bw(odd)
                diff = even(:,kk,:) - real(bwfilter(kk),dp)*odd(:,kk,:)
                ! Weighted sum(|diff|2) with shell volume and normalization factor.
                wk = real(box,dp)**3 * (4.d0/3.d0)*DPI * (real(kk,dp)**3 - real(kk-1,dp)**3)
                wk = wk / real(2*pftc%ncls*pftc%pftsz,dp)
                cost = cost + wk * sum(real(diff*conjg(diff),dp))
            enddo
            costs(k) = cost
            deallocate(bwfilter, diff)
        enddo
        !$omp end parallel do
        ! Sequential reduction: find k that achieves the minimum
        best_cost = huge(best_cost)
        best_k    = ksearch(1)
        do k = ksearch(1), ksearch(2)
            if( DEBUG ) write(logfhandle,'(A,I8,A,ES16.8,A,L1)') '    k=', k, ' cost=', costs(k), ' new_best=', costs(k) <= best_cost
            if( costs(k) <= best_cost )then
                best_cost = costs(k)
                best_k    = k
            endif
        enddo
        deallocate(even, odd, costs)
        ! Solution
        lpopt = calc_lowpass_lim(best_k, box, smpd)
        if( DEBUG )then
            write(logfhandle,602) '    best_k=', best_k
            write(logfhandle,604) '    best_cost=', best_cost
            write(logfhandle,601) '    lpopt [A]=', lpopt
        endif
    end subroutine polarft_estimate_lplim3D

end module simple_polarft_calc
