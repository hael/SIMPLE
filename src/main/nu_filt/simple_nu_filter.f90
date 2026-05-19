!@descr: volume-domain nonuniform filtering of even/odd volumes
!
! Static low-pass bank sequence:
!    call setup_nu_dmats(vol_even, vol_odd, l_mask, [real ::])
!    call optimize_nu_cutoff_finds()
!    call nu_filter_vols(vol_even_filt, vol_odd_filt)
!    call cleanup_nu_filter()
!
! Iterative high-resolution refinement can challenge one Fourier shell at a
! time after the static bank has been optimized. Each challenger allocates and
! evaluates only one extra candidate; callers may loop while the shell is
! accepted. The default acceptance gate for extend_nu_filter_highres_shell_next
! is the refinement threshold.
!    call setup_nu_dmats(vol_even, vol_odd, l_mask, [real ::])
!    call optimize_nu_cutoff_finds()
!    do
!        call extend_nu_filter_highres_shell_next(vol_even, vol_odd, stats=ext_stats)
!        if( .not.ext_stats%attempted    ) exit
!        if( .not.ext_stats%applied      ) exit
!        if( .not.ext_stats%promote_next ) exit
!    end do
!    call nu_filter_vols(vol_even_filt, vol_odd_filt)
!    call cleanup_nu_filter()
!
! Postprocess workflows can instead challenge Fourier-shell candidates
! sequentially with the loop helper, then apply the resulting local filter map
! to a merged map using resolution-dependent B factors followed by a shared
! global antialiasing filter:
!    call setup_nu_dmats(vol_even, vol_odd, l_mask, [real ::])
!    call optimize_nu_cutoff_finds()
!    call extend_nu_filter_highres_shells(vol_even, vol_odd, accept_pct=0.)
!    call nu_postprocess_vol(vol, vol_lp, vol_pproc, fsc0143, bfac, optlp)
!
! Auxiliary candidate pairs supplied through setup_nu_dmats compete with the
! base low-pass bank during voxelwise optimization.
!
module simple_nu_filter
use simple_core_module_api
use simple_image, only: image
use simple_butterworth
use simple_tent_smooth, only: tent_smooth_3d
use simple_neighs,      only: neigh_8_3D
implicit none

public :: setup_nu_dmats, optimize_nu_cutoff_finds, nu_filter_vols, nu_filter_vol, nu_postprocess_vol, &
          cleanup_nu_filter, pack_filtmap_lowpass_limits,&
          calc_filtmap_lowpass_stats, print_nu_filtmap_lowpass_stats, calc_filtmap_lowpass_histogram,&
          print_filtmap_lowpass_histogram, extend_nu_filter_highres_shell_next, extend_nu_filter_highres_shells,&
          analyze_filtmap_neighbor_continuity, nu_highres_extension_stats, get_nu_filter_bank_finest_lp
private
#include "simple_local_flags.inc"

real,             parameter   :: lowpass_limits(8) = [20.,15.,12.,10.,8.,6.,5.,4.]
! Minimum finest-frontier fraction of the NU mask required before testing a
! finer shell. Zero means challenge whenever at least one frontier voxel exists.
real,             parameter   :: NU_HIGHRES_EXTENSION_THRESHOLD_PCT = 0.
real,             parameter   :: NU_REFINE_EXTENSION_ACCEPT_PCT     = 5.
real,             parameter   :: NU_POSTPROCESS_EXTENSION_ACCEPT_PCT = 0.
! NU postprocess B-factor model tuning:
! - ALPHA controls positive-B damping for bins worse than the global FSC
!   resolution. Increase it if low-resolution regions remain too sharp/noisy;
!   decrease it if those regions become too blurred.
! - BETA controls the hard negative-B sharpening ceiling for bins better than
!   the global FSC resolution. Increase it if high-resolution regions remain
!   too soft; decrease it if they look noisy or over-sharpened.
! - RATIO_MAX controls how much the squared resolution ratio can drive the
!   interpolation before the hard B-factor clamp. Increase it for stronger
!   separation between bins; decrease it for more conservative variation.
real,             parameter   :: NU_POSTPROCESS_BFAC_ALPHA           = 0.75
real,             parameter   :: NU_POSTPROCESS_BFAC_BETA            = 2.
real,             parameter   :: NU_POSTPROCESS_BFAC_RATIO_MIN       = 0.
real,             parameter   :: NU_POSTPROCESS_BFAC_RATIO_MAX       = 2.
! Physical half-width of the tent regularization kernel. The smoother consumes
! this as an integer pixel radius, so the full tent base spans 2*radius + 1
! voxels along each axis; 8 A at 1 A/px gives radius=8 and a 17-voxel base.
real,             parameter   :: WINSZ_TENT_ANGSTROM       = 8.
integer,          parameter   :: DISCONT_STEP_THRESH       = 1
integer,          parameter   :: NU_LABEL_SMOOTH_MAXITS    = 3
integer,          parameter   :: NU_LABEL_SMOOTH_STEP_TOL  = 1
integer,          parameter   :: NU_LABEL_SMOOTH_NNEIGH    = 26
integer,          parameter   :: NU_LABEL_SMOOTH_NCOLORS   = 8
real,             parameter   :: NU_LABEL_SMOOTH_BETA_FRAC = 2.0
real,             parameter   :: NU_LABEL_SMOOTH_QUAD_FRAC = 1.0
real,             parameter   :: NU_LABEL_SMOOTH_TIE_EPS   = 1.e-6
character(len=*), parameter   :: NU_FILTER_CACHE_EVEN      = 'nu_filter_cache_even'
character(len=*), parameter   :: NU_FILTER_CACHE_ODD       = 'nu_filter_cache_odd'
real,             allocatable :: dmats_mask(:,:)
real,             allocatable :: bwfilters(:,:)
real,             allocatable :: candidate_coords(:)
integer,          allocatable :: filtmap(:,:,:)
integer,          allocatable :: srcmap(:,:,:)
integer,          allocatable :: cutoff_finds(:)
real,             allocatable :: dmat_finest_cached(:,:,:)
logical,          allocatable :: nu_lmask(:,:,:)
integer,          allocatable :: nu_mask_index(:,:,:)
integer,          allocatable :: nu_mask_vox(:,:)
real,             allocatable :: nu_smooth_norm(:,:,:)
type(image),      allocatable :: aux_even_bank(:), aux_odd_bank(:)
integer :: ldim(3), box
integer :: n_nu_mask = 0
integer :: winsz_tent
real    :: smpd

type :: nu_highres_extension_stats
    logical :: attempted    = .false.
    logical :: applied      = .false.
    logical :: promote_next = .false.
    real    :: new_limit = 0.
    integer :: n_mask     = 0
    integer :: n_tested   = 0
    integer :: n_extended = 0
    real    :: pct_tested_mask     = 0.
    real    :: pct_extended_tested = 0.
end type nu_highres_extension_stats


interface

    ! In submodule: simple_nu_filter_state.f90
    module real function cutoff_find_to_lowpass_limit( icut )
        integer, intent(in) :: icut
    end function cutoff_find_to_lowpass_limit

    module subroutine init_nu_filter( vol_even, vol_odd, n_highres_steps )
        class(image), intent(in) :: vol_even, vol_odd
        integer, optional, intent(in) :: n_highres_steps
    end subroutine init_nu_filter

    module function filtered_vol_fname( cache_prefix, cutoff_find ) result( fname )
        class(string), intent(in) :: cache_prefix
        integer,       intent(in) :: cutoff_find
        type(string) :: fname
    end function filtered_vol_fname

    module logical function filtered_vols_cached( cache_prefix )
        class(string), intent(in) :: cache_prefix
    end function filtered_vols_cached

    module subroutine delete_cached_filtered_vols( cache_prefix )
        class(string), intent(in) :: cache_prefix
    end subroutine delete_cached_filtered_vols

    module subroutine cleanup_nu_filter
    end subroutine cleanup_nu_filter

    module subroutine cleanup_aux_bank
    end subroutine cleanup_aux_bank

    module subroutine validate_aux_volumes( aux_even, aux_odd )
        type(image), intent(in) :: aux_even(:), aux_odd(:)
    end subroutine validate_aux_volumes

    module subroutine stash_aux_volumes( aux_even, aux_odd )
        type(image), intent(in) :: aux_even(:), aux_odd(:)
    end subroutine stash_aux_volumes

    module subroutine setup_nu_mask_index
    end subroutine setup_nu_mask_index

    module subroutine smooth_nu_objective( dmat, tmp )
        real, intent(inout) :: dmat(:,:,:)
        real, intent(inout) :: tmp(:,:,:)
    end subroutine smooth_nu_objective

    module subroutine pack_nu_dmat_candidate( dmat_full, icand )
        real,    intent(in) :: dmat_full(:,:,:)
        integer, intent(in) :: icand
    end subroutine pack_nu_dmat_candidate

    module subroutine unpack_nu_dmat_candidate( icand, dmat_full )
        integer, intent(in)  :: icand
        real,    intent(out) :: dmat_full(:,:,:)
    end subroutine unpack_nu_dmat_candidate

    module subroutine cache_filtered_vols( vol_even, vol_odd )
        class(image), intent(in) :: vol_even, vol_odd
    end subroutine cache_filtered_vols

    module subroutine generate_single_filtered_pair( vol_even, vol_odd, cutoff_find, even_cache_fname, odd_cache_fname )
        class(image),  intent(in) :: vol_even, vol_odd
        integer,       intent(in) :: cutoff_find
        class(string), intent(in) :: even_cache_fname, odd_cache_fname
    end subroutine generate_single_filtered_pair

    module subroutine delete_cached_filtered_pair( cutoff_find )
        integer, intent(in) :: cutoff_find
    end subroutine delete_cached_filtered_pair

    ! In submodule: simple_nu_filter_bank.f90
    module subroutine setup_nu_dmats( vol_even, vol_odd, l_mask, aux_resolutions, aux_even, aux_odd, n_highres_steps )
        class(image),          intent(in) :: vol_even, vol_odd
        logical,               intent(in) :: l_mask(:,:,:)
        real,                  intent(in) :: aux_resolutions(:)
        type(image), optional, intent(in) :: aux_even(:), aux_odd(:)
        integer,     optional, intent(in) :: n_highres_steps
    end subroutine setup_nu_dmats

    module subroutine setup_nu_candidate_coords( n_candidates, aux_resolutions )
        integer, intent(in) :: n_candidates
        real,    intent(in) :: aux_resolutions(:)
    end subroutine setup_nu_candidate_coords

    module real function lowpass_limit_to_candidate_coord( lp_angstrom )
        real, intent(in) :: lp_angstrom
    end function lowpass_limit_to_candidate_coord

    module real function get_nu_filter_bank_finest_lp()
    end function get_nu_filter_bank_finest_lp

    module subroutine optimize_nu_cutoff_finds()
    end subroutine optimize_nu_cutoff_finds

    module subroutine candidate_map_to_filt_and_src( candmap, n_base )
        integer, intent(in) :: candmap(:,:,:), n_base
    end subroutine candidate_map_to_filt_and_src

    module subroutine log_nu_candidate_selection_counts( candmap, n_base, stage )
        integer,          intent(in) :: candmap(:,:,:), n_base
        character(len=*), intent(in) :: stage
    end subroutine log_nu_candidate_selection_counts

    module subroutine log_nu_aux_unary_margin_stats( n_base )
        integer, intent(in) :: n_base
    end subroutine log_nu_aux_unary_margin_stats

    module subroutine log_nu_candidate_coords
    end subroutine log_nu_candidate_coords

    module real function nu_candidate_coord_for_label( ilabel )
        integer, intent(in) :: ilabel
    end function nu_candidate_coord_for_label

    module integer function nu_effective_base_label_for_candidate( icand, n_base )
        integer, intent(in) :: icand, n_base
    end function nu_effective_base_label_for_candidate

    ! In submodule: simple_nu_filter_potts.f90
    module subroutine refine_nu_candidate_map_ordered_labels( candmap, n_candidates )
        integer, intent(inout) :: candmap(:,:,:)
        integer, intent(in)    :: n_candidates
    end subroutine refine_nu_candidate_map_ordered_labels

    module real function estimate_nu_label_smooth_beta( n_candidates )
        integer, intent(in) :: n_candidates
    end function estimate_nu_label_smooth_beta

    module real function nu_label_smooth_neighborhood_cost( icand, candmap, neigh, nsz )
        integer, intent(in) :: icand, candmap(:,:,:), neigh(3,NU_LABEL_SMOOTH_NNEIGH), nsz
    end function nu_label_smooth_neighborhood_cost

    module real function nu_label_smooth_pair_cost( icand, jcand )
        integer, intent(in) :: icand, jcand
    end function nu_label_smooth_pair_cost

    module real function nu_label_smooth_coord_pair_cost( icoord, jcoord )
        real, intent(in) :: icoord, jcoord
    end function nu_label_smooth_coord_pair_cost

    module integer function nu_label_smooth_color( i, j, k )
        integer, intent(in) :: i, j, k
    end function nu_label_smooth_color

    module logical function nu_label_smooth_is_better( e, best_e )
        real, intent(in) :: e, best_e
    end function nu_label_smooth_is_better

    module real function calc_nu_label_smooth_site_energy( candmap, beta )
        integer, intent(in) :: candmap(:,:,:)
        real,    intent(in) :: beta
    end function calc_nu_label_smooth_site_energy

    ! In submodule: simple_nu_filter_extend.f90
    module subroutine extend_nu_filter_highres( vol_even, vol_odd, threshold_pct, accept_pct, new_limit, stats )
        class(image), intent(in) :: vol_even, vol_odd
        real,         intent(in) :: threshold_pct   ! e.g. 10.0
        real,         intent(in) :: accept_pct      ! minimum selected frontier percentage
        real,         intent(in) :: new_limit        ! Angstrom limit for the proposed shell
        type(nu_highres_extension_stats), optional, intent(out) :: stats
    end subroutine extend_nu_filter_highres

    module subroutine extend_nu_filter_highres_shell_next( vol_even, vol_odd, stats, accept_pct )
        class(image), intent(in) :: vol_even, vol_odd
        type(nu_highres_extension_stats), optional, intent(out) :: stats
        real, optional, intent(in) :: accept_pct
    end subroutine extend_nu_filter_highres_shell_next

    module subroutine extend_nu_filter_highres_shells( vol_even, vol_odd, nsteps, accept_pct )
        class(image), intent(in) :: vol_even, vol_odd
        integer, optional, intent(out) :: nsteps
        real, optional, intent(in) :: accept_pct
    end subroutine extend_nu_filter_highres_shells

    module subroutine init_nu_highres_extension_selection( extend_mask, dmat_old, dmat_new, extend_to_new, n_extended )
        logical, intent(in)    :: extend_mask(:,:,:)
        real,    intent(in)    :: dmat_old(:,:,:), dmat_new(:,:,:)
        logical, intent(inout) :: extend_to_new(:,:,:)
        integer, intent(out)   :: n_extended
    end subroutine init_nu_highres_extension_selection

    module subroutine refine_nu_highres_extension_selection( extend_mask, dmat_old, dmat_new, extend_to_new, old_label, n_extended )
        logical, intent(in)    :: extend_mask(:,:,:)
        real,    intent(in)    :: dmat_old(:,:,:), dmat_new(:,:,:)
        logical, intent(inout) :: extend_to_new(:,:,:)
        integer, intent(in)    :: old_label
        integer, intent(out)   :: n_extended
    end subroutine refine_nu_highres_extension_selection

    module real function estimate_nu_highres_extension_beta( extend_mask, dmat_old, dmat_new )
        logical, intent(in) :: extend_mask(:,:,:)
        real,    intent(in) :: dmat_old(:,:,:), dmat_new(:,:,:)
    end function estimate_nu_highres_extension_beta

    module real function nu_highres_extension_neighborhood_cost( icoord, extend_mask, extend_to_new, old_label, &
        &new_coord, neigh, nsz )
        real,    intent(in) :: icoord, new_coord
        logical, intent(in) :: extend_mask(:,:,:), extend_to_new(:,:,:)
        integer, intent(in) :: old_label, neigh(3,NU_LABEL_SMOOTH_NNEIGH), nsz
    end function nu_highres_extension_neighborhood_cost

    module real function nu_highres_extension_current_coord( i, j, k, extend_mask, extend_to_new, old_label, new_coord )
        integer, intent(in) :: i, j, k, old_label
        logical, intent(in) :: extend_mask(:,:,:), extend_to_new(:,:,:)
        real,    intent(in) :: new_coord
    end function nu_highres_extension_current_coord

    module subroutine apply_nu_highres_extension_selection( extend_to_new, new_label )
        logical, intent(in) :: extend_to_new(:,:,:)
        integer, intent(in) :: new_label
    end subroutine apply_nu_highres_extension_selection

    module subroutine append_nu_highres_candidate_coord( old_n_base, new_coord )
        integer, intent(in) :: old_n_base
        real,    intent(in) :: new_coord
    end subroutine append_nu_highres_candidate_coord

    ! In submodule: simple_nu_filter_apply.f90
    module subroutine nu_filter_vols( vol_even, vol_odd )
        class(image), intent(out) :: vol_even, vol_odd
    end subroutine nu_filter_vols

    module subroutine nu_filter_vol( vol_in, vol_out )
        class(image), intent(in)  :: vol_in
        class(image), intent(out) :: vol_out
    end subroutine nu_filter_vol

    module subroutine nu_postprocess_vol( vol_in, vol_lp, vol_pproc, global_lp, global_bfac, fsc_filter )
        class(image), intent(in)  :: vol_in
        class(image), intent(out) :: vol_lp, vol_pproc
        real,         intent(in)  :: global_lp, global_bfac
        real, optional, intent(in) :: fsc_filter(:)
    end subroutine nu_postprocess_vol

    ! In submodule: simple_nu_filter_stats.f90
    module subroutine pack_filtmap_lowpass_limits( lowpass_vals, mask )
        real, allocatable, intent(inout) :: lowpass_vals(:)
        logical, intent(in)              :: mask(:,:,:)
    end subroutine pack_filtmap_lowpass_limits

    module subroutine calc_filtmap_lowpass_stats( statvars, mask )
        type(stats_struct), intent(out) :: statvars
        logical, intent(in)             :: mask(:,:,:)
    end subroutine calc_filtmap_lowpass_stats

    module subroutine calc_filtmap_lowpass_histogram( counts, percentages, mask )
        integer, intent(out) :: counts(:)
        real,    intent(out) :: percentages(:)
        logical, intent(in)  :: mask(:,:,:)
    end subroutine calc_filtmap_lowpass_histogram

    module subroutine print_filtmap_lowpass_histogram( mask, aux_resolutions )
        logical,        intent(in) :: mask(:,:,:)
        real, optional, intent(in) :: aux_resolutions(:)
    end subroutine print_filtmap_lowpass_histogram

    module subroutine print_nu_filtmap_lowpass_stats( mask, aux_resolutions )
        logical,        intent(in) :: mask(:,:,:)
        real, optional, intent(in) :: aux_resolutions(:)
    end subroutine print_nu_filtmap_lowpass_stats

    module subroutine analyze_filtmap_neighbor_continuity( mask )
        logical, intent(in) :: mask(:,:,:)
    end subroutine analyze_filtmap_neighbor_continuity

end interface

end module simple_nu_filter
