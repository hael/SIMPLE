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
! accepted. The challenger starts from the finest populated base-bank label, so
! empty finer discrete labels do not block shell refinement. Refinement-style
! callers accept the next shell only when the challenger wins enough absolute
! support and at least 5% of the tested frontier. Diagnostic callers can pass
! accept_pct=0. for an explicitly permissive shell walk. Accepted extension
! shells are challenged at the full Fourier sampling rate, but the retained
! high-resolution bank is thinned to every second shell plus the current
! terminal shell. After the shell walk stops, the final accepted label map is
! cleaned with the ordered-label Potts prior over the expanded retained bank.
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

public :: setup_nu_dmats, optimize_nu_cutoff_finds, nu_filter_vols, nu_filter_vol, &
          cleanup_nu_filter, pack_filtmap_lowpass_limits,&
          calc_filtmap_lowpass_stats, print_nu_filtmap_lowpass_stats, calc_filtmap_lowpass_histogram,&
          print_filtmap_lowpass_histogram, extend_nu_filter_highres_shell_next, extend_nu_filter_highres_shells,&
          refine_nu_extension_filtmap_ordered_labels, analyze_filtmap_neighbor_continuity,&
          nu_highres_extension_stats, get_nu_filter_bank_finest_lp, get_nu_filtmap_finest_selected_lp,&
          get_nu_filtmap_highres_shell_depth
private
#include "simple_local_flags.inc"

real,             parameter   :: lowpass_limits(8) = [20.,15.,12.,10.,8.,6.,5.,4.]
! Minimum finest-frontier fraction of the NU mask required before testing a
! finer shell. Zero means challenge whenever at least one frontier voxel exists.
real,             parameter   :: NU_HIGHRES_EXTENSION_THRESHOLD_PCT  = 0.
! Percentage of the tested frontier that must select a challenger before it is
! accepted into refinement-style NU banks. Diagnostic callers can pass
! accept_pct=0. to request a permissive one-voxel shell walk.
real,             parameter   :: NU_HIGHRES_EXTENSION_ACCEPT_PCT     = 5.0
! Refinement-style shell acceptance also requires enough absolute support so a
! tiny frontier cannot march indefinitely. Diagnostic callers can pass
! accept_pct=0. to bypass both the frontier fraction and this seed floor.
integer,          parameter   :: NU_HIGHRES_EXTENSION_MIN_SEED_VOXELS = 32
! Hard cap on mask-packed distance-matrix columns retained for NU optimization.
! When this fills, unselected high-resolution labels are compacted away before
! another shell is accepted.
integer,          parameter   :: NU_DMAT_CANDIDATE_CAP                = 24
! High-resolution extension still challenges every shell, but retained
! extension candidates are thinned to this stride. With the default stride of
! two, odd shell steps are kept only as the temporary frontier needed to test
! the next even shell step.
integer,          parameter   :: NU_HIGHRES_EXTENSION_RETAIN_STRIDE    = 2
integer,          parameter   :: NU_DMAT_CANDIDATE_HEADROOM           = 2
! Candidate-scale objective smoothing. The normalized unary objective for a
! candidate with low-pass L is averaged over an AWF-like local support:
! radius_A = 0.5 * NU_OBJECTIVE_SMOOTH_AWF * L, capped below. Increasing AWF
! makes local evidence more collective; lowering it makes assignments more
! voxel-local. The cap prevents very coarse candidates from washing out the
! objective over an unrealistically large support.
real,             parameter   :: NU_OBJECTIVE_SMOOTH_AWF          = 3.0
real,             parameter   :: NU_OBJECTIVE_SMOOTH_RADIUS_FRAC  = 0.5
real,             parameter   :: NU_OBJECTIVE_SMOOTH_MAX_RADIUS_A = 30.0
! Report every nonzero retained-bank coordinate jump as a continuity boundary.
! Larger jumps are still separated in the diagnostics.
integer,          parameter   :: DISCONT_STEP_THRESH         = 0
integer,          parameter   :: NU_CONTINUITY_LARGE_STEP_THRESH = 1
integer,          parameter   :: NU_LABEL_SMOOTH_MAXITS      = 6
! Adjacent retained-bank coordinate jumps are tolerated by the ordered-label
! Potts prior; the quadratic hinge makes larger jumps increasingly expensive.
integer,          parameter   :: NU_LABEL_SMOOTH_STEP_TOL    = 1
integer,          parameter   :: NU_LABEL_SMOOTH_NNEIGH      = 26
integer,          parameter   :: NU_LABEL_SMOOTH_NCOLORS     = 8
real,             parameter   :: NU_LABEL_SMOOTH_BETA_FRAC   = 2.0
real,             parameter   :: NU_LABEL_SMOOTH_QUAD_FRAC   = 1.0
real,             parameter   :: NU_LABEL_SMOOTH_TIE_EPS     = 1.e-6
! Optional unordered auxiliary sources are source alternatives rather than
! ordered resolution-bin labels. Keep this at zero to let the unary objective
! decide aux-vs-base assignments without frequency-ladder Potts bias; raise
! slightly if a future dataset needs source-boundary cleanup.
real,             parameter   :: NU_AUX_SOURCE_BOUNDARY_COST = 0.0
integer,          parameter   :: NU_LABEL_KIND = selected_int_kind(4)
character(len=*), parameter   :: NU_FILTER_CACHE_EVEN        = 'nu_filter_cache_even'
character(len=*), parameter   :: NU_FILTER_CACHE_ODD         = 'nu_filter_cache_odd'
real,             allocatable :: dmats_mask(:,:)
real,             allocatable :: dmats_aux_mask(:,:)
real,             allocatable :: bwfilters(:,:)
real,             allocatable :: candidate_coords(:)
real,             allocatable :: aux_candidate_resolutions(:)
integer(kind=NU_LABEL_KIND), allocatable :: filtmap(:,:,:)
integer(kind=NU_LABEL_KIND), allocatable :: srcmap(:,:,:)
integer,          allocatable :: cutoff_finds(:)
real,             allocatable :: dmat_finest_cached(:)
logical,          allocatable :: nu_lmask(:,:,:)
integer,          allocatable :: nu_mask_vox(:,:)
real,             allocatable :: nu_smooth_norm(:,:,:)
type(image),      allocatable :: aux_even_bank(:), aux_odd_bank(:)
integer :: ldim(3), box
integer :: n_nu_mask = 0
integer :: nu_smooth_norm_radius = -1
real    :: smpd
! Cache of the raw E/O noise scale used to normalize the Huber unary objective.
! Computed once per setup_nu_dmats and reused across all shell-extension
! challenges so the per-shell median/MAD pass is not paid repeatedly.
real    :: nu_noise_sigma_cached = 0.
logical :: l_aux_source_unordered_potts = .false.

type :: nu_highres_extension_stats
    logical :: attempted    = .false.
    logical :: applied      = .false.
    logical :: promote_next = .false.
    real    :: old_limit = 0.
    real    :: new_limit = 0.
    integer :: old_find  = 0
    integer :: new_find  = 0
    integer :: n_mask     = 0
    integer :: n_tested   = 0
    integer :: n_unary_wins = 0
    integer :: n_extended = 0
    integer :: n_seed_min = 0
    real    :: pct_tested_mask     = 0.
    real    :: pct_unary_wins_tested = 0.
    real    :: pct_unary_wins_mask = 0.
    real    :: pct_extended_tested = 0.
    logical :: accepted_by_frontier = .false.
    logical :: memory_limited       = .false.
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

    module logical function keep_nu_highres_extension_step( istep, finest_step )
        integer, intent(in) :: istep, finest_step
    end function keep_nu_highres_extension_step

    module integer function count_nu_highres_extension_retained_steps( nsteps )
        integer, intent(in) :: nsteps
    end function count_nu_highres_extension_retained_steps

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

    module subroutine release_nu_filter_unary_storage
    end subroutine release_nu_filter_unary_storage

    module subroutine release_nu_smooth_norm
    end subroutine release_nu_smooth_norm

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

    module subroutine setup_nu_mask_voxels
    end subroutine setup_nu_mask_voxels

    module real function nu_objective_smooth_radius_angstrom( lp_angstrom )
        real, intent(in) :: lp_angstrom
    end function nu_objective_smooth_radius_angstrom

    module integer function nu_objective_smooth_radius_pixels( lp_angstrom )
        real, intent(in) :: lp_angstrom
    end function nu_objective_smooth_radius_pixels

    module subroutine smooth_nu_objective( dmat, tmp, lp_angstrom )
        real, intent(inout) :: dmat(:,:,:)
        real, intent(inout) :: tmp(:,:,:)
        real, intent(in)    :: lp_angstrom
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
    module subroutine setup_nu_dmats( vol_even, vol_odd, l_mask, aux_resolutions, aux_even, aux_odd, &
            &n_highres_steps, l_aux_source_unordered )
        class(image),          intent(in) :: vol_even, vol_odd
        logical,               intent(in) :: l_mask(:,:,:)
        real,                  intent(in) :: aux_resolutions(:)
        type(image), optional, intent(in) :: aux_even(:), aux_odd(:)
        integer,     optional, intent(in) :: n_highres_steps
        logical,     optional, intent(in) :: l_aux_source_unordered
    end subroutine setup_nu_dmats

    module subroutine setup_nu_candidate_coords( n_candidates, aux_resolutions )
        integer, intent(in) :: n_candidates
        real,    intent(in) :: aux_resolutions(:)
    end subroutine setup_nu_candidate_coords

    module subroutine refresh_nu_aux_candidate_coords()
    end subroutine refresh_nu_aux_candidate_coords

    module real function lowpass_limit_to_candidate_coord( lp_angstrom )
        real, intent(in) :: lp_angstrom
    end function lowpass_limit_to_candidate_coord

    module real function get_nu_filter_bank_finest_lp()
    end function get_nu_filter_bank_finest_lp

    module integer function get_nu_filtmap_highres_shell_depth()
    end function get_nu_filtmap_highres_shell_depth

    module subroutine optimize_nu_cutoff_finds()
    end subroutine optimize_nu_cutoff_finds

    module subroutine candidate_map_to_filt_and_src( candmap, n_base )
        integer(kind=NU_LABEL_KIND), intent(in) :: candmap(:,:,:)
        integer, intent(in) :: n_base
    end subroutine candidate_map_to_filt_and_src

    module subroutine log_nu_candidate_selection_counts( candmap, n_base, stage )
        integer(kind=NU_LABEL_KIND), intent(in) :: candmap(:,:,:)
        integer,          intent(in) :: n_base
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
        integer(kind=NU_LABEL_KIND), intent(inout) :: candmap(:,:,:)
        integer, intent(in)    :: n_candidates
    end subroutine refine_nu_candidate_map_ordered_labels

    module real function estimate_nu_label_smooth_beta( n_candidates )
        integer, intent(in) :: n_candidates
    end function estimate_nu_label_smooth_beta

    module real function nu_label_smooth_neighborhood_cost( icand, candmap, neigh, nsz )
        integer, intent(in) :: icand, neigh(3,NU_LABEL_SMOOTH_NNEIGH), nsz
        integer(kind=NU_LABEL_KIND), intent(in) :: candmap(:,:,:)
    end function nu_label_smooth_neighborhood_cost

    module real function nu_label_smooth_pair_cost( icand, jcand )
        integer, intent(in) :: icand, jcand
    end function nu_label_smooth_pair_cost

    module real function nu_label_smooth_coord_pair_cost( icoord, jcoord )
        real, intent(in) :: icoord, jcoord
    end function nu_label_smooth_coord_pair_cost

    module real function nu_label_smooth_source_pair_cost( icoord, jcoord, l_aux_i, l_aux_j )
        real,    intent(in) :: icoord, jcoord
        logical, intent(in) :: l_aux_i, l_aux_j
    end function nu_label_smooth_source_pair_cost

    module integer function nu_label_smooth_color( i, j, k )
        integer, intent(in) :: i, j, k
    end function nu_label_smooth_color

    module logical function nu_label_smooth_is_better( e, best_e )
        real, intent(in) :: e, best_e
    end function nu_label_smooth_is_better

    module real function calc_nu_label_smooth_site_energy( candmap, beta )
        integer(kind=NU_LABEL_KIND), intent(in) :: candmap(:,:,:)
        real,    intent(in) :: beta
    end function calc_nu_label_smooth_site_energy

    ! In submodule: simple_nu_filter_extend.f90
    module subroutine extend_nu_filter_highres( vol_even, vol_odd, threshold_pct, new_limit, stats, accept_pct )
        class(image), intent(in) :: vol_even, vol_odd
        real,         intent(in) :: threshold_pct   ! e.g. 10.0
        real,         intent(in) :: new_limit        ! Angstrom limit for the proposed shell
        type(nu_highres_extension_stats), optional, intent(out) :: stats
        real, optional, intent(in) :: accept_pct
    end subroutine extend_nu_filter_highres

    module subroutine extend_nu_filter_highres_shell_next( vol_even, vol_odd, stats, accept_pct )
        class(image),                               intent(in)  :: vol_even, vol_odd
        type(nu_highres_extension_stats), optional, intent(out) :: stats
        real, optional, intent(in) :: accept_pct
    end subroutine extend_nu_filter_highres_shell_next

    module subroutine extend_nu_filter_highres_shells( vol_even, vol_odd, nsteps, accept_pct )
        class(image), intent(in) :: vol_even, vol_odd
        integer, optional, intent(out) :: nsteps
        real, optional, intent(in) :: accept_pct
    end subroutine extend_nu_filter_highres_shells

    module subroutine refine_nu_extension_filtmap_ordered_labels
    end subroutine refine_nu_extension_filtmap_ordered_labels

    module subroutine init_nu_highres_extension_selection( frontier_vox, dmat_old, dmat_new, &
            &extend_choice, n_extended )
        integer, intent(in)    :: frontier_vox(:)
        real,    intent(in)    :: dmat_old(:), dmat_new(:,:,:)
        integer(kind=NU_LABEL_KIND), intent(inout) :: extend_choice(:)
        integer, intent(out)   :: n_extended
    end subroutine init_nu_highres_extension_selection

    module subroutine init_nu_highres_extension_selection_aux( frontier_vox, dmat_old, dmat_new, &
            &extend_choice, n_extended )
        integer, intent(in)    :: frontier_vox(:)
        real,    intent(in)    :: dmat_old(:), dmat_new(:,:,:)
        integer(kind=NU_LABEL_KIND), intent(inout) :: extend_choice(:)
        integer, intent(out)   :: n_extended
    end subroutine init_nu_highres_extension_selection_aux

    module subroutine apply_nu_highres_extension_selection( frontier_vox, extend_choice, old_label, new_label )
        integer, intent(in) :: frontier_vox(:)
        integer(kind=NU_LABEL_KIND), intent(in) :: extend_choice(:)
        integer, intent(in) :: old_label, new_label
    end subroutine apply_nu_highres_extension_selection

    ! In submodule: simple_nu_filter_apply.f90
    module subroutine nu_filter_vols( vol_even, vol_odd )
        class(image), intent(out) :: vol_even, vol_odd
    end subroutine nu_filter_vols

    module subroutine nu_filter_vol( vol_in, vol_out )
        class(image), intent(in)  :: vol_in
        class(image), intent(out) :: vol_out
    end subroutine nu_filter_vol

    ! In submodule: simple_nu_filter_stats.f90
    module subroutine pack_filtmap_lowpass_limits( lowpass_vals, mask )
        real, allocatable, intent(inout) :: lowpass_vals(:)
        logical, optional, intent(in)    :: mask(:,:,:)
    end subroutine pack_filtmap_lowpass_limits

    module subroutine calc_filtmap_lowpass_stats( statvars, mask )
        type(stats_struct), intent(out) :: statvars
        logical, optional, intent(in)   :: mask(:,:,:)
    end subroutine calc_filtmap_lowpass_stats

    module subroutine calc_filtmap_lowpass_histogram( counts, percentages, mask )
        integer, intent(out) :: counts(:)
        real,    intent(out) :: percentages(:)
        logical, optional, intent(in) :: mask(:,:,:)
    end subroutine calc_filtmap_lowpass_histogram

    module real function get_nu_filtmap_finest_selected_lp( mask, aux_resolutions )
        logical, optional, intent(in) :: mask(:,:,:)
        real, optional, intent(in) :: aux_resolutions(:)
    end function get_nu_filtmap_finest_selected_lp

    module subroutine print_filtmap_lowpass_histogram( mask, aux_resolutions )
        logical, optional, intent(in) :: mask(:,:,:)
        real, optional, intent(in) :: aux_resolutions(:)
    end subroutine print_filtmap_lowpass_histogram

    module subroutine print_nu_filtmap_lowpass_stats( mask, aux_resolutions )
        logical, optional, intent(in) :: mask(:,:,:)
        real, optional, intent(in) :: aux_resolutions(:)
    end subroutine print_nu_filtmap_lowpass_stats

    module subroutine analyze_filtmap_neighbor_continuity( mask )
        logical, optional, intent(in) :: mask(:,:,:)
    end subroutine analyze_filtmap_neighbor_continuity

end interface

end module simple_nu_filter
