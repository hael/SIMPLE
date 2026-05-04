!@descr: polarft class complete interface
module simple_polarft_calc
use simple_pftc_api
implicit none

public :: polarft_calc
public :: vol_pad2ref_pfts
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

type fftw_crvec
    type(c_ptr)                            :: p
    complex(kind=c_float_complex), pointer :: c(:) => null()  ! (pftsz+1)
    real(kind=c_float),            pointer :: r(:) => null()  ! (nrots+2) view on same memory
end type fftw_crvec

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
    type(c_ptr) :: plan_fwd1_many, plan_bwd1_many, plan_bwd1_single, plan_mem_r2c_many
    ! Memoized terms
    complex(kind=c_float_complex), allocatable :: ft_ptcl_ctf(:,:,:)                      !< Fourier Transform of particle times CTF
    complex(kind=c_float_complex), allocatable :: ft_absptcl_ctf(:,:,:)                   !< Fourier Transform of (particle times CTF)**2
    complex(kind=c_float_complex), allocatable :: ft_ctf2(:,:,:)                          !< Fourier Transform of CTF squared modulus
    complex(kind=c_float_complex), allocatable :: ft_ref_even(:,:,:),  ft_ref_odd(:,:,:)  !< Fourier Transform of even/odd references
    complex(kind=c_float_complex), allocatable :: ft_ref2_even(:,:,:), ft_ref2_odd(:,:,:) !< Fourier Transform of even/odd references squared modulus
    ! buffer: (nrots,nk) holding cvec2 for all k
    type(fftw_cmat),               allocatable :: cmat2_many(:)  ! per thread
    ! batched backward (c2r) plan: nk transforms of length nrots, in-place
    type(fftw_crmat),              allocatable :: crmat1_many(:) ! per thread
    ! single backward (c2r) plan: one k-collapsed transform of length nrots
    type(fftw_crvec),              allocatable :: crvec1(:)      ! per thread
    ! Convenience vectors, thread memoization
    type(heap_vars),               allocatable :: heap_vars(:)
    type(fftw_drvec),              allocatable :: drvec(:)
    ! Others
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
    ! ===== CTF: simple_polarft_ctf.f90
    procedure          :: create_polar_absctfmats
    ! ===== GEOM: simple_polarft_geom.f90
    procedure, private :: gen_shmat, gen_shmat4aln
    procedure, private :: gen_shmat_8, gen_shmat4aln_8
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
    procedure, private :: gen_corr_cc_for_rot_8_1, gen_corr_cc_for_rot_8_2
    procedure, private :: gen_euclid_for_rot_8_1, gen_euclid_for_rot_8_2
    procedure          :: gen_corr_grad_for_rot_8
    procedure, private :: gen_corr_cc_grad_for_rot_8
    procedure, private :: gen_euclid_grad_for_rot_8
    procedure          :: gen_corr_grad_only_for_rot_8
    procedure          :: gen_sigma_contrib
    procedure          :: gen_objfun_vals_mirr_vals
    procedure, private :: gen_corrs_mirr_corrs, gen_euclids_mirr_euclids
    ! ===== I/O: simple_polarft_ops_io.f90
    procedure          :: write_ptcl_pft_range
    procedure          :: write_ref_pfts
    procedure          :: read_ref_pfts
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

    module subroutine vol_pad2ref_pfts(self, vol_pad, eulspace, state, iseven)
        use simple_projector, only: projector
        class(polarft_calc), intent(inout) :: self
        class(projector),    intent(in)    :: vol_pad
        class(oris),         intent(inout) :: eulspace
        integer,             intent(in)    :: state
        logical,             intent(in)    :: iseven
    end subroutine vol_pad2ref_pfts

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

    module subroutine gen_shmat4aln(self, ithr, shift, shmat)
        class(polarft_calc),  intent(inout) :: self
        integer,              intent(in)    :: ithr
        real(sp),             intent(in)    :: shift(2)
        complex(sp), pointer, intent(inout) :: shmat(:,:)
    end subroutine gen_shmat4aln

    module subroutine gen_shmat_8(self, ithr, shift_8, shmat_8)
        class(polarft_calc),  intent(inout) :: self
        integer,              intent(in)    :: ithr
        real(dp),             intent(in)    :: shift_8(2)
        complex(dp), pointer, intent(inout) :: shmat_8(:,:)
    end subroutine gen_shmat_8

    module subroutine gen_shmat4aln_8(self, ithr, shift_8, shmat_8)
        class(polarft_calc),  intent(inout) :: self
        integer,              intent(in)    :: ithr
        real(dp),             intent(in)    :: shift_8(2)
        complex(dp), pointer, intent(inout) :: shmat_8(:,:)
    end subroutine gen_shmat4aln_8

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

    module subroutine memoize_refs( self, eulspace )
        class(polarft_calc),           intent(inout) :: self
        class(oris),         optional, intent(in)    :: eulspace
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

    module real(dp) function gen_corr_for_rot_8_1( self, iref, iptcl, irot )
        class(polarft_calc), target, intent(inout) :: self
        integer,                     intent(in)    :: iref, iptcl, irot
    end function gen_corr_for_rot_8_1

    module real(dp) function gen_corr_for_rot_8_2( self, iref, iptcl, shvec, irot )
        class(polarft_calc), target, intent(inout) :: self
        integer,                     intent(in)    :: iref, iptcl
        real(dp),                    intent(in)    :: shvec(2)
        integer,                     intent(in)    :: irot
    end function gen_corr_for_rot_8_2

    module real(dp) function gen_corr_cc_for_rot_8_1( self, pft_ref, i, irot )
        class(polarft_calc),  intent(inout) :: self
        complex(dp), pointer, intent(inout) :: pft_ref(:,:)
        integer,              intent(in)    :: i, irot
    end function gen_corr_cc_for_rot_8_1

    module real(dp) function gen_corr_cc_for_rot_8_2( self, pft_ref, i, shvec, irot )
        class(polarft_calc), target, intent(inout) :: self
        complex(dp),        pointer, intent(inout) :: pft_ref(:,:)
        integer,                     intent(in)    :: i, irot
        real(dp),                    intent(in)    :: shvec(2)
    end function gen_corr_cc_for_rot_8_2

    module real(dp) function gen_euclid_for_rot_8_1( self, pft_ref, iptcl, irot )
        class(polarft_calc), target, intent(inout) :: self
        complex(dp),        pointer, intent(inout) :: pft_ref(:,:)
        integer,                     intent(in)    :: iptcl, irot
    end function gen_euclid_for_rot_8_1

    module real(dp) function gen_euclid_for_rot_8_2( self, pft_ref, iptcl, shvec, irot )
        class(polarft_calc), target, intent(inout) :: self
        complex(dp),        pointer, intent(inout) :: pft_ref(:,:)
        integer,                     intent(in)    :: iptcl, irot
        real(dp),                    intent(in)    :: shvec(2)
    end function gen_euclid_for_rot_8_2

    module subroutine gen_corr_grad_for_rot_8(self, iref, iptcl, shvec, irot, f, grad)
        class(polarft_calc), target, intent(inout) :: self
        integer,                     intent(in)    :: iref, iptcl, irot
        real(dp),                    intent(in)    :: shvec(2)
        real(dp),                    intent(out)   :: f, grad(2)
    end subroutine gen_corr_grad_for_rot_8

    module subroutine gen_corr_cc_grad_for_rot_8( self, pft_ref, i, shvec, irot, f, grad)
        class(polarft_calc), target, intent(inout) :: self
        complex(dp),        pointer, intent(inout) :: pft_ref(:,:)
        integer,                     intent(in)    :: i, irot
        real(dp),                    intent(in)    :: shvec(2)
        real(dp),                    intent(out)   :: f, grad(2)
    end subroutine gen_corr_cc_grad_for_rot_8

    module subroutine gen_euclid_grad_for_rot_8(self, pft_ref, iptcl, shvec, irot, f, grad)
        class(polarft_calc), target, intent(inout) :: self
        complex(dp), pointer,        intent(inout) :: pft_ref(:,:)
        integer,                     intent(in)    :: iptcl, irot
        real(dp),                    intent(in)    :: shvec(2)
        real(dp),                    intent(out)   :: f, grad(2)
    end subroutine gen_euclid_grad_for_rot_8

    module subroutine gen_corr_grad_only_for_rot_8( self, iref, iptcl, shvec, irot, grad )
        class(polarft_calc), target, intent(inout) :: self
        integer,                     intent(in)    :: iref, iptcl
        real(dp),                    intent(in)    :: shvec(2)
        integer,                     intent(in)    :: irot
        real(dp),                    intent(out)   :: grad(2)
    end subroutine gen_corr_grad_only_for_rot_8

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

    ! ===== I/O: simple_polarft_ops_io.f90

    module subroutine write_ptcl_pft_range( self, fname, nptcls_total, iwrite_from, iwrite_to )
        class(polarft_calc), intent(in) :: self
        class(string),       intent(in) :: fname
        integer,             intent(in) :: nptcls_total
        integer,             intent(in) :: iwrite_from, iwrite_to
    end subroutine write_ptcl_pft_range

    module subroutine write_ref_pfts( self, fname, iseven )
        class(polarft_calc), intent(in) :: self
        class(string),       intent(in) :: fname
        logical,             intent(in) :: iseven
    end subroutine write_ref_pfts

    module subroutine read_ref_pfts( self, fname, iseven )
        class(polarft_calc), intent(inout) :: self
        class(string),       intent(in)    :: fname
        logical,             intent(in)    :: iseven
    end subroutine read_ref_pfts

end interface

! CLASS PARAMETERS/VARIABLES
complex(sp), parameter       :: zero            = cmplx(0.,0.) !< just a complex zero
integer,     parameter       :: FFTW_USE_WISDOM = 16

contains

end module simple_polarft_calc
