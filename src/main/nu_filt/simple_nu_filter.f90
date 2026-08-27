!@descr: volume-domain nonuniform filtering of even/odd volumes
!
! Static low-pass bank sequence:
!    call setup_nu_dmats(vol_even, vol_odd, mskdiam, [real ::])
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
!    call setup_nu_dmats(vol_even, vol_odd, mskdiam, [real ::])
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
! Auxiliary pairs supplied through setup_nu_dmats are only used when
! they extend beyond the finest discrete member. In that case the auxiliary pair
! replaces the finest discrete label instead of appending a sidecar candidate.
!
module simple_nu_filter
use simple_core_module_api
use simple_image, only: image
use simple_butterworth
use simple_tent_smooth, only: tent_smooth_3d
use simple_neighs,      only: neigh_8_3D
use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
implicit none

public :: setup_nu_dmats, optimize_nu_cutoff_finds, nu_filter_vols, nu_filter_vol, &
          cleanup_nu_filter, pack_filtmap_lowpass_limits,&
          calc_filtmap_lowpass_stats, print_nu_filtmap_lowpass_stats, calc_filtmap_lowpass_histogram,&
          print_filtmap_lowpass_histogram, extend_nu_filter_highres_shell_next, extend_nu_filter_highres_shells,&
          refine_nu_extension_filtmap_ordered_labels, analyze_filtmap_neighbor_continuity,&
          nu_highres_extension_stats, get_nu_filter_bank_finest_lp, get_nu_filtmap_finest_selected_lp,&
          get_nu_filtmap_highres_shell_depth, write_nu_local_resolution_map, set_nu_filter_report, NU_DEV_OUTPUT,&
          nu_envmask_params, nu_envmask_stats, nu_evidence_envelope, calc_nu_evidence_margin,&
          write_nu_evidence_map, print_nu_envmask_stats, NU_ENVMASK_BETA, NU_ENVMASK_DENS_WEIGHT,&
          NU_ENVMASK_RELATIVE, NU_ENVMASK_MINVOL_FRAC, NU_ENVMASK_GROW_A, NU_ENVMASK_EDGE_A,&
          nu_evidence_state, nu_evidence_summary, build_nu_evidence_state, unpack_nu_evidence_state,&
          get_nu_evidence_summary, nu_evidence_state_is_valid, print_nu_evidence_summary,&
          expand_nu_evidence_band_weights, assert_nu_evidence_replay_ready,&
          NU_EVIDENCE_NBANDS, NU_EVIDENCE_BAND_LIMITS, NU_EVIDENCE_MIN_NULL_FRAC
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
integer,          parameter   :: NU_DMAT_CANDIDATE_HEADROOM            = 2
! Candidate-scale objective smoothing. The normalized unary objective for a
! candidate with low-pass L is averaged over an AWF-like local support:
! radius_A = 0.5 * NU_OBJECTIVE_SMOOTH_AWF * L, capped below. Increasing AWF
! makes local evidence more collective; lowering it makes assignments more
! voxel-local. The cap prevents very coarse candidates from washing out the
! objective over an unrealistically large support.
real,             parameter   :: NU_OBJECTIVE_SMOOTH_AWF          = 3.0
real,             parameter   :: NU_OBJECTIVE_SMOOTH_RADIUS_FRAC  = 0.5
real,             parameter   :: NU_OBJECTIVE_SMOOTH_MAX_RADIUS_A = 30.0
! Report continuity health using the same hinge as the ordered-label prior:
! one-step retained-bank transitions are tolerated, larger jumps are penalized.
integer,          parameter   :: DISCONT_STEP_THRESH         = 1
integer,          parameter   :: NU_LABEL_SMOOTH_MAXITS      = 6
! Adjacent retained-bank coordinate jumps are tolerated by the ordered-label
! Potts prior; the quadratic hinge makes larger jumps increasingly expensive.
integer,          parameter   :: NU_LABEL_SMOOTH_STEP_TOL    = 1
integer,          parameter   :: NU_LABEL_SMOOTH_NNEIGH      = 26
integer,          parameter   :: NU_LABEL_SMOOTH_NCOLORS     = 8
real,             parameter   :: NU_LABEL_SMOOTH_BETA_FRAC   = 2.0
real,             parameter   :: NU_LABEL_SMOOTH_QUAD_FRAC   = 1.0
real,             parameter   :: NU_LABEL_SMOOTH_TIE_EPS     = 1.e-6
integer,          parameter   :: NU_LABEL_KIND               = selected_int_kind(4)
! Fixed standalone NU-envelope policy. Only the evidence threshold and scale
! remain user-tunable; morphology lengths are physical and converted to pixels
! by nu_filt3D at the input-map sampling distance.
! Binary MRF boundary smoothness; higher values give smoother boundaries.
real,             parameter   :: NU_ENVMASK_BETA              = 1.0
! Weight of local density evidence, which can retain strong but poorly ordered density.
real,             parameter   :: NU_ENVMASK_DENS_WEIGHT       = 0.0
! A scale-free ratio can stop a high-contrast core outvoting weak, ordered density.
logical,          parameter   :: NU_ENVMASK_RELATIVE          = .false.
! Smallest connected component kept, expressed as a fraction of the largest.
real,             parameter   :: NU_ENVMASK_MINVOL_FRAC       = 0.1
! Physical binary growth applied before the soft edge. Callers floor the pixel
! conversion at one layer, so at any sampling coarser than ~1 A this is a
! single-voxel dilation rather than a true 1 A margin. Raise it if envelopes are
! observed to clip flexible periphery: reference masking suppresses whatever it
! excludes, so a slightly generous margin is much cheaper than a clipped domain.
real,             parameter   :: NU_ENVMASK_GROW_A            = 1.0
! Physical cosine-edge width used to soften the molecular envelope.
real,             parameter   :: NU_ENVMASK_EDGE_A            = 6.0
! Sentinel floor above which a mask-packed unary entry is treated as unpopulated.
! Columns are allocated with huge() and compaction can leave stale members behind,
! so evidence comparisons must ignore anything at that magnitude.
real,             parameter   :: NU_EVIDENCE_INVALID         = 0.5 * huge(1.)
! Four broad support bands are the compact interface between the full NU unary
! bank and the future replay precision.  Confidence in successive entries means
! reproducible detail is supported through 20, 12, 8, and 5 A respectively.
! The sets are nested, so the packed confidences must be monotonically
! non-increasing from coarse to fine.
integer,          parameter   :: NU_EVIDENCE_NBANDS = 4
real,             parameter   :: NU_EVIDENCE_BAND_LIMITS(NU_EVIDENCE_NBANDS) = [20., 12., 8., 5.]
real,             parameter   :: NU_EVIDENCE_UNCERTAIN_ENTROPY = 0.5
real,             parameter   :: NU_EVIDENCE_NULL_NSIGMA = 3.0
! Replay-readiness contract: the generous spherical support always contains
! substantial solvent, so a compact state whose explicit null wins less than
! this fraction of the support is evidence of a failed null calibration (the
! observed zero-null failure mode), not of a solvent-free box. The replay must
! hard-error on it rather than attach an uncalibrated precision. Provisional
! (R9), anchored to the simple_test_nu_envmask fixture; recalibrate against
! real-data operating points before relaxing.
real,             parameter   :: NU_EVIDENCE_MIN_NULL_FRAC = 0.01
character(len=*), parameter   :: NU_EVIDENCE_SOURCE_BASE = 'base_unfil'
character(len=*), parameter   :: NU_EVIDENCE_ALGORITHM = 'nu_evidence_v1'
character(len=*), parameter   :: NU_FILTER_CACHE_EVEN        = 'nu_filter_cache_even'
character(len=*), parameter   :: NU_FILTER_CACHE_ODD         = 'nu_filter_cache_odd'
real,             allocatable :: dmats_mask(:,:)
real,             allocatable :: bwfilters(:,:)
real,             allocatable :: candidate_coords(:)
integer(kind=NU_LABEL_KIND), allocatable :: filtmap(:,:,:)
integer,          allocatable :: cutoff_finds(:)
real,             allocatable :: dmat_finest_cached(:)
! Raw, unsmoothed unary costs kept for envelope masking. dmats_mask is smoothed at
! candidate-dependent radii (30 A for the coarsest member, 6 A for the finest),
! which is right for label selection but wrong for locating a boundary: it blurs
! the coarse baseline five times harder than its competitors and erodes the
! envelope inward on a 30 A scale. The envelope therefore gets its own evidence
! from these, smoothed once, symmetrically, at a scale the caller chooses.
real,             allocatable :: nu_ev_base(:)
real,             allocatable :: nu_ev_best(:)
logical,          allocatable :: nu_lmask(:,:,:)
integer,          allocatable :: nu_mask_vox(:,:)
real,             allocatable :: nu_smooth_norm(:,:,:)
type(image),      allocatable :: aux_even_bank(:), aux_odd_bank(:)
integer :: ldim(3), box
integer :: n_nu_mask = 0
integer :: nu_smooth_norm_radius = -1
real    :: smpd, nu_support_mskdiam = 0.
logical :: nu_l_report = .true.
! Opt-in diagnostics for NU-filter development. Keep normal runs concise; this
! flag restores the detailed candidate, shell-extension, and continuity logs.
logical :: NU_DEV_OUTPUT = .false.
! Cache of the radial raw E/O noise profile that whitens the Huber unary
! objective (see image::nu_objective_noise_profile). Computed once per
! setup_nu_dmats and reused across all shell-extension challenges so the
! per-shell median/MAD pass is not paid repeatedly.
real, allocatable :: nu_noise_profile_cached(:)
real    :: nu_noise_rmax_cached = 0.
integer :: nu_aux_replacement_label = 0
real    :: nu_aux_replacement_resolution = 0.
logical :: nu_evidence_requested = .false.
character(len=32) :: nu_evidence_source = ''
real(kind=8) :: nu_evidence_source_fingerprint(6) = 0.d0

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

! Controls for NU-evidence-driven envelope masking. beta regularizes boundary
! area only; connectivity and hole filling are the caller's responsibility.
type :: nu_envmask_params
    real    :: nsigma      = 3.0   ! threshold, in null MADs above the null median
    ! Binary MRF boundary smoothness; higher values give smoother boundaries.
    real    :: beta        = NU_ENVMASK_BETA
    ! Local density weight; positive values retain strong but poorly ordered density.
    real    :: dens_weight = NU_ENVMASK_DENS_WEIGHT
    real    :: lp_smooth   = 8.0   ! scale, in Angstrom, at which the envelope is defined
    ! Use a scale-free baseline-to-best ratio so weak ordered density is not outvoted.
    logical :: l_relative  = NU_ENVMASK_RELATIVE
    integer :: maxits      = 6     ! ICM sweeps
end type nu_envmask_params

type :: nu_envmask_stats
    real    :: null_med    = 0.
    real    :: null_mad    = 0.
    real    :: thres       = 0.
    real    :: nsigma      = 0.
    real    :: beta_used   = 0.
    real    :: dens_med    = 0.
    real    :: dens_mad    = 0.
    real    :: dens_weight = 0.
    integer :: n_support   = 0
    integer :: n_seed      = 0
    integer :: n_signal    = 0
    integer :: nits        = 0
    real    :: pct_seed    = 0.
    real    :: pct_signal  = 0.
    real    :: lp_smooth   = 0.
    logical :: l_relative  = .false.
end type nu_envmask_stats

! Public scalar metadata for a frozen NU evidence state.  The large packed
! arrays remain private in nu_evidence_state and can only be copied out through
! unpack_nu_evidence_state, preventing accidental mutation between half replays.
type :: nu_evidence_summary
    logical :: valid = .false.
    integer :: ldim(3) = 0
    integer :: n_support = 0
    integer :: n_candidates = 0
    integer :: n_bands = NU_EVIDENCE_NBANDS
    real    :: smpd = 0.
    real    :: mskdiam = 0.
    real    :: null_fraction = 0.
    real    :: uncertain_fraction = 0.
    real    :: calibration_temperature = 0.
    real    :: spatial_beta = 0.
    real    :: null_cost_mean = 0.
    real    :: null_bias_median = 0.
    real    :: null_bias_mad = 0.
    real    :: null_bias_threshold = 0.
    real    :: supported_fraction(NU_EVIDENCE_NBANDS) = 0.
    real    :: band_limits(NU_EVIDENCE_NBANDS) = NU_EVIDENCE_BAND_LIMITS
    character(len=32) :: source = ''
    character(len=16) :: identity = ''
    character(len=LONGSTRLEN) :: provenance = ''
end type nu_evidence_summary

type :: nu_evidence_state
    private
    type(nu_evidence_summary) :: summary
    integer(kind=NU_LABEL_KIND), allocatable :: selected_label(:)
    real, allocatable :: selected_cutoff(:)
    real, allocatable :: uncertainty(:)
    real, allocatable :: band_support(:,:)
end type nu_evidence_state

interface

    ! In submodule: simple_nu_filter_state.f90
    module real function cutoff_find_to_lowpass_limit( icut )
        integer, intent(in) :: icut
    end function cutoff_find_to_lowpass_limit

    module real function nu_label_lowpass_limit( ilabel )
        integer, intent(in) :: ilabel
    end function nu_label_lowpass_limit

    module logical function nu_label_is_aux_replacement( ilabel )
        integer, intent(in) :: ilabel
    end function nu_label_is_aux_replacement

    module subroutine init_nu_filter( vol_even, vol_odd, n_highres_steps )
        class(image), intent(in) :: vol_even, vol_odd
        integer, optional, intent(in) :: n_highres_steps
    end subroutine init_nu_filter

    module subroutine set_nu_filter_report( l_report )
        logical, intent(in) :: l_report
    end subroutine set_nu_filter_report

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
    module subroutine setup_nu_dmats( vol_even, vol_odd, mskdiam, aux_resolutions, aux_even, aux_odd, &
            &n_highres_steps, evidence_source )
        class(image),          intent(in) :: vol_even, vol_odd
        real,                  intent(in) :: mskdiam
        real,                  intent(in) :: aux_resolutions(:)
        type(image), optional, intent(in) :: aux_even(:), aux_odd(:)
        integer,     optional, intent(in) :: n_highres_steps
        character(len=*), optional, intent(in) :: evidence_source
    end subroutine setup_nu_dmats

    module subroutine setup_nu_candidate_coords( n_candidates )
        integer, intent(in) :: n_candidates
    end subroutine setup_nu_candidate_coords

    module real function get_nu_filter_bank_finest_lp()
    end function get_nu_filter_bank_finest_lp

    module integer function get_nu_filtmap_highres_shell_depth()
    end function get_nu_filtmap_highres_shell_depth

    module subroutine optimize_nu_cutoff_finds()
    end subroutine optimize_nu_cutoff_finds

    module subroutine clamp_nu_filtmap_labels( n_base )
        integer, intent(in) :: n_base
    end subroutine clamp_nu_filtmap_labels

    module subroutine log_nu_candidate_selection_counts( candmap, n_base, stage )
        integer(kind=NU_LABEL_KIND), intent(in) :: candmap(:,:,:)
        integer,          intent(in) :: n_base
        character(len=*), intent(in) :: stage
    end subroutine log_nu_candidate_selection_counts

    module subroutine log_nu_candidate_coords
    end subroutine log_nu_candidate_coords

    module real function nu_candidate_coord_for_label( ilabel )
        integer, intent(in) :: ilabel
    end function nu_candidate_coord_for_label

    module integer function nu_effective_base_label_for_candidate( icand, n_base )
        integer, intent(in) :: icand, n_base
    end function nu_effective_base_label_for_candidate

    ! In submodule: simple_nu_filter_evidence.f90
    module subroutine build_nu_evidence_state( vol_even, vol_odd, state )
        class(image),            intent(in)  :: vol_even, vol_odd
        type(nu_evidence_state), intent(out) :: state
    end subroutine build_nu_evidence_state

    module subroutine calculate_nu_source_fingerprint( vol_even, vol_odd, fingerprint )
        class(image), intent(in) :: vol_even, vol_odd
        real(kind=8), intent(out) :: fingerprint(6)
    end subroutine calculate_nu_source_fingerprint

    module logical function nu_evidence_state_is_valid( state )
        type(nu_evidence_state), intent(in) :: state
    end function nu_evidence_state_is_valid

    module subroutine get_nu_evidence_summary( state, summary )
        type(nu_evidence_state),   intent(in)  :: state
        type(nu_evidence_summary), intent(out) :: summary
    end subroutine get_nu_evidence_summary

    module subroutine unpack_nu_evidence_state( state, selected_label, selected_cutoff, uncertainty, band_support )
        type(nu_evidence_state), intent(in) :: state
        integer, allocatable, optional, intent(out) :: selected_label(:)
        real,    allocatable, optional, intent(out) :: selected_cutoff(:), uncertainty(:), band_support(:,:)
    end subroutine unpack_nu_evidence_state

    module subroutine print_nu_evidence_summary( state )
        type(nu_evidence_state), intent(in) :: state
    end subroutine print_nu_evidence_summary

    module subroutine expand_nu_evidence_band_weights( state, band_w )
        type(nu_evidence_state), intent(in)  :: state
        real, allocatable,       intent(out) :: band_w(:,:,:,:)
    end subroutine expand_nu_evidence_band_weights

    module subroutine assert_nu_evidence_replay_ready( state )
        type(nu_evidence_state), intent(in) :: state
    end subroutine assert_nu_evidence_replay_ready

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

    module real function get_nu_filtmap_finest_selected_lp( mask )
        logical, optional, intent(in) :: mask(:,:,:)
    end function get_nu_filtmap_finest_selected_lp

    module subroutine print_filtmap_lowpass_histogram( mask )
        logical, optional, intent(in) :: mask(:,:,:)
    end subroutine print_filtmap_lowpass_histogram

    module subroutine print_nu_filtmap_lowpass_stats( mask )
        logical, optional, intent(in) :: mask(:,:,:)
    end subroutine print_nu_filtmap_lowpass_stats

    module subroutine analyze_filtmap_neighbor_continuity( mask )
        logical, optional, intent(in) :: mask(:,:,:)
    end subroutine analyze_filtmap_neighbor_continuity

    module subroutine write_nu_local_resolution_map( fname, mask )
        class(string), intent(in) :: fname
        logical, optional, intent(in) :: mask(:,:,:)
    end subroutine write_nu_local_resolution_map

    module subroutine accumulate_nu_evidence_raw( dmat_full, icand )
        real,    intent(in) :: dmat_full(:,:,:)
        integer, intent(in) :: icand
    end subroutine accumulate_nu_evidence_raw

    ! In submodule: simple_nu_filter_envmask.f90
    module subroutine calc_nu_evidence_margin( margin, lp_smooth, l_relative )
        real, allocatable, intent(inout) :: margin(:)
        real,    optional, intent(in)    :: lp_smooth
        logical, optional, intent(in)    :: l_relative
    end subroutine calc_nu_evidence_margin

    module real function nu_evidence_baseline_floor( base_full ) result( floor_val )
        real, intent(in) :: base_full(:,:,:)
    end function nu_evidence_baseline_floor

    module subroutine calc_nu_evidence_score( margin, nsigma, score, stats )
        real,                   intent(in)    :: margin(:)
        real,                   intent(in)    :: nsigma
        real, allocatable,      intent(inout) :: score(:)
        type(nu_envmask_stats), intent(inout) :: stats
    end subroutine calc_nu_evidence_score

    module subroutine add_nu_evidence_density( vol_dens, weight, score, stats )
        class(image), target,   intent(in)    :: vol_dens
        real,                   intent(in)    :: weight
        real, allocatable,      intent(inout) :: score(:)
        type(nu_envmask_stats), intent(inout) :: stats
    end subroutine add_nu_evidence_density

    module subroutine segment_nu_evidence( score, p, lmask, stats )
        real,                    intent(in)    :: score(:)
        type(nu_envmask_params), intent(in)    :: p
        logical, allocatable,    intent(inout) :: lmask(:,:,:)
        type(nu_envmask_stats),  intent(inout) :: stats
    end subroutine segment_nu_evidence

    module subroutine nu_evidence_envelope( p, lmask, stats, vol_dens )
        type(nu_envmask_params), intent(in)    :: p
        logical, allocatable,    intent(inout) :: lmask(:,:,:)
        type(nu_envmask_stats),  intent(inout) :: stats
        class(image), optional, target, intent(in) :: vol_dens
    end subroutine nu_evidence_envelope

    module subroutine write_nu_evidence_map( fname, lp_smooth, l_relative )
        class(string),     intent(in) :: fname
        real,    optional, intent(in) :: lp_smooth
        logical, optional, intent(in) :: l_relative
    end subroutine write_nu_evidence_map

    module subroutine print_nu_envmask_stats( stats )
        type(nu_envmask_stats), intent(in) :: stats
    end subroutine print_nu_envmask_stats

end interface

end module simple_nu_filter
