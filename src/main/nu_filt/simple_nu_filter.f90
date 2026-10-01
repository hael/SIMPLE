!@descr: volume-domain nonuniform filtering of even/odd volumes
! Sequence: setup_nu_dmats -> optimize_nu_cutoff_finds -> nu_filter_vols -> cleanup_nu_filter.
! Bank: static ladder lowpass_limits, cut at fsc_res/NU_BANK_FSC_HEADROOM when given (>= 2 rungs);
! an auxiliary (ML) pair at or beyond the finest retained rung is appended at that rung's Potts coordinate.
! The last label is the finest member and the matching low-pass handoff.
module simple_nu_filter
use simple_core_module_api
use simple_image, only: image
use simple_butterworth
use simple_tent_smooth, only: tent_smooth_3d
use simple_neighs,      only: neigh_8_3D
use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
implicit none

public :: setup_nu_dmats, optimize_nu_cutoff_finds, nu_filter_vols, nu_filter_vol, &
          set_nu_solvent_envelope, clear_nu_solvent_envelope, &
          set_nu_evidence_null_shell, clear_nu_evidence_null_shell, &
          cleanup_nu_filter, pack_filtmap_lowpass_limits,&
          calc_filtmap_lowpass_stats, print_nu_filtmap_lowpass_stats, calc_filtmap_lowpass_histogram,&
          print_filtmap_lowpass_histogram, analyze_filtmap_neighbor_continuity,&
          get_nu_filter_bank_finest_lp, get_nu_filtmap_finest_selected_lp,&
          write_nu_local_resolution_map, set_nu_filter_report, NU_DEV_OUTPUT,&
          nu_envmask_params, nu_envmask_stats, nu_evidence_envelope, calc_nu_evidence_margin,&
          write_nu_evidence_map, write_nu_evidence_envmask, print_nu_envmask_stats, NU_ENVMASK_BETA, NU_ENVMASK_DENS_WEIGHT,&
          NU_ENVMASK_RELATIVE, NU_ENVMASK_MINVOL_FRAC, NU_ENVMASK_GROW_A, NU_ENVMASK_EDGE_A,&
          nu_evidence_state, nu_evidence_summary, build_nu_evidence_state, unpack_nu_evidence_state,&
          NU_ALIGN_LP_MIN_ASSIGNED_PCT, NU_ALIGN_LP_MIN_SIGNAL_PCT,&
          get_nu_evidence_summary, nu_evidence_state_is_valid, print_nu_evidence_summary,&
          assert_nu_evidence_replay_ready,&
          nu_evidence_sharpen_vol,&
          NU_EVIDENCE_NBANDS, NU_EVIDENCE_BAND_LIMITS, NU_EVIDENCE_MIN_NULL_FRAC,&
          NU_EVIDENCE_MAX_NULL_FRAC, NU_EVIDENCE_SOURCE_BASE, NU_EVIDENCE_SOURCE_PREV
private
#include "simple_local_flags.inc"

real,             parameter   :: lowpass_limits(8) = [20.,15.,12.,10.,8.,6.,5.,4.]
! Static-bank cap: with fsc_res, keep rungs at or coarser than fsc_res/NU_BANK_FSC_HEADROOM (>= 2);
! without it (nu_filt3D, flex_pca) the bank is uncapped. Above the two-rung floor the cap bounds the
! handoff's lead over the FSC; rationale in doc/policies/NU/nonuniform_filtering_policy.md sections 8, 12.
real,             parameter   :: NU_BANK_FSC_HEADROOM = 1.5
! Hard cap on mask-packed distance-matrix columns retained for NU optimization.
integer,          parameter   :: NU_DMAT_CANDIDATE_CAP                = 24
! Candidate-scale objective smoothing. The normalized unary objective for a
! candidate with low-pass L is averaged over an AWF-like local support:
! radius_A = 0.5 * NU_OBJECTIVE_SMOOTH_AWF * L, capped below. Increasing AWF
! makes local evidence more collective; lowering it makes assignments more
! voxel-local. The cap prevents very coarse candidates from washing out the
! objective over an unrealistically large support.
real,             parameter   :: NU_OBJECTIVE_SMOOTH_AWF          = 3.0
real,             parameter   :: NU_OBJECTIVE_SMOOTH_RADIUS_FRAC  = 0.5
real,             parameter   :: NU_OBJECTIVE_SMOOTH_MAX_RADIUS_A = 30.0
! Largest representable NU unary with useful relative precision. The image
! objective uses the same bound before its single-precision smoothing volume
! is formed. This is still overwhelmingly unfavorable relative to ordinary
! noise-normalized Huber costs, without behaving like a numeric sentinel.
real,             parameter   :: NU_OBJECTIVE_UNARY_CAP           = 1. / epsilon(1.)
! Report continuity health using the same hinge as the ordered-label prior:
! one-step retained-bank transitions are tolerated, larger jumps are penalized.
integer,          parameter   :: DISCONT_STEP_THRESH         = 1
integer,          parameter   :: NU_LABEL_SMOOTH_MAXITS      = 6
! Adjacent ladder-coordinate jumps are tolerated by the ordered-label Potts prior; the linear-quadratic
! hinge prices larger ones. The auxiliary label shares the finest rung's coordinate (setup_nu_candidate_coords).
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
!> score, in null MADs, assigned to voxels outside the evidence calibration
!! domain: a fixed solvent label the ICM boundary term (beta ~ 1) cannot flip
real,             parameter   :: NU_ENVMASK_EXCLUDED_SCORE    = -1.0e3
!> a Euclidean null shell must carry at least this many voxels, and this
!! fraction of the label domain, for its median/MAD to be trusted
integer,          parameter   :: NU_ENVMASK_MIN_NULL_VOX      = 1000
real,             parameter   :: NU_ENVMASK_MIN_NULL_FRAC     = 0.02
! Sentinel floor above which a mask-packed unary entry is treated as unpopulated.
! Columns are allocated with huge() and compaction can leave stale members behind,
! so evidence comparisons must ignore anything at that magnitude.
real,             parameter   :: NU_EVIDENCE_INVALID         = 0.5 * huge(1.)
! Static nested support bands of the compact evidence state: entry b is the confidence that reproducible
! detail is supported through 20, 12, 8, 5 A, so entries are non-increasing from coarse to fine.
integer,          parameter   :: NU_EVIDENCE_NBANDS = 4
real,             parameter   :: NU_EVIDENCE_BAND_LIMITS(NU_EVIDENCE_NBANDS) = [20., 12., 8., 5.]
! Evidence bands beyond the static four are appended at NU_EVIDENCE_BAND_RATIO steps while the finest
! signal candidate reaches the next boundary, up to NU_EVIDENCE_MAX_NBANDS; the 4 A ladder floor lies
! above the first one (3.2 A), so in practice the four static bands apply.
integer,          parameter   :: NU_EVIDENCE_MAX_NBANDS      = 8    !< band-count cap
real,             parameter   :: NU_EVIDENCE_BAND_RATIO      = 0.64 !< geometric step, matching the 20->12->8->5 spacing
!> an appended band is kept only if its mean support reaches this fraction (pruned finest-first)
real,             parameter   :: NU_EVIDENCE_MIN_BAND_SUPPORT = 0.01
!> Default gate of get_nu_filtmap_finest_selected_lp: the finest cutoff selected (it or finer) by at
!! least this percentage of the assigned voxels. Diagnostic; the flex_pca nufilt report uses this default.
real,             parameter   :: NU_ALIGN_LP_MIN_ASSIGNED_PCT = 5.0
! Diagnostic only: finest label whose cumulative population reaches this % of the signal voxels
! (mask minus the background clamp), quoted on the NU MATCHING LOW-PASS HANDOFF line.
! The handoff itself is the finest bank member (get_nu_filter_bank_finest_lp).
real,             parameter   :: NU_ALIGN_LP_MIN_SIGNAL_PCT   = 1.0
real,             parameter   :: NU_EVIDENCE_UNCERTAIN_ENTROPY = 0.5
! NU-evidence nonuniform postprocessing v2 (nu_evidence_local_sharpening.md,
! postprocess_nu commander): classical shrink-then-sharpen, localized by the
! evidence. One Guinier B-factor inside the evidenced local passband; no user
! gain knob.
real,             parameter   :: NU_SHARP_BFAC_FINEST_A = 5.0 !< Guinier sharpening only when the finest evidenced cutoff is finer (mirrors the standard postprocess lp<5A gate)
! Readiness bounds of the compact evidence state (assert_nu_evidence_replay_ready): the generous
! sphere always holds solvent and molecule, so an explicit-null fraction outside [MIN, MAX] is a
! failed null calibration. Provisional; recalibrate against real-data operating points before relaxing.
real,             parameter   :: NU_EVIDENCE_MIN_NULL_FRAC = 0.01
real,             parameter   :: NU_EVIDENCE_MAX_NULL_FRAC = 0.90
character(len=*), parameter   :: NU_EVIDENCE_SOURCE_BASE = 'base_unfil'
! lag-one evidence source for the PCG trailing bootstrap: the previous
! iteration's shipped half pair, the same pair that supplies the bootstrap FSC
character(len=*), parameter   :: NU_EVIDENCE_SOURCE_PREV = 'previous_shipped'
character(len=*), parameter   :: NU_EVIDENCE_ALGORITHM = 'nu_evidence_v1'
character(len=*), parameter   :: NU_FILTER_CACHE_EVEN        = 'nu_filter_cache_even'
character(len=*), parameter   :: NU_FILTER_CACHE_ODD         = 'nu_filter_cache_odd'
real,             allocatable :: dmats_mask(:,:)
! Raw (unsmoothed) mask-packed unaries of every candidate, for the like-for-like selection in
! optimize_nu_cutoff_finds, which dmats_mask's per-candidate radii cannot give. Released with dmats_mask.
real,             allocatable :: raw_dmats_mask(:,:)
real,             allocatable :: bwfilters(:,:)
real,             allocatable :: candidate_coords(:)
integer(kind=NU_LABEL_KIND), allocatable :: filtmap(:,:,:)
integer,          allocatable :: cutoff_finds(:)
! Raw coarsest-candidate (nu_ev_base) and best-candidate (nu_ev_best) unaries for envelope masking:
! dmats_mask's per-candidate radii (30 A for the coarsest) would erode the boundary, so
! calc_nu_evidence_margin smooths their difference once at the caller's scale.
real,             allocatable :: nu_ev_base(:)
real,             allocatable :: nu_ev_best(:)
logical,          allocatable :: nu_lmask(:,:,:)
integer,          allocatable :: nu_mask_vox(:,:)
! Packed (nu_mask_vox order) observation mask of the evidence pair: .false.
! where both half maps are exactly zero. A density-constrained PCG solve
! leaves such voxels inside the broader spherical NU domain; they are
! boundary conditions, not measurements, so the evidence compaction keeps
! them out of every calibration statistic and freezes them at the null label.
logical,          allocatable :: nu_observed_mask(:)
integer :: n_nu_observed = 0
real,             allocatable :: nu_smooth_norm(:,:,:)
! Solvent-constraint clamp (pcg_priors_history.md dev item 4): the objective and null
! competition retain the broad spherical domain. The whitening fit omits exact
! zero/zero samples introduced by a narrower PCG solve support because those
! are boundary conditions, not noise observations. Solvent voxels are fixed to
! the coarsest signal candidate in both the filter and evidence Potts fields,
! preserving coarse support while suppressing unsupported finer bands. Set by
! the caller after setup_nu_dmats; cleared by cleanup_nu_filter; absent = no clamp.
logical,          allocatable :: nu_solvent_lmask(:,:,:) !< .true. = outside the envelope (solvent)
logical :: nu_l_solvent_clamp = .false.
! Which envelope defines the filter-field background: the conservative
! density automask (automsk=no solvent constraint) or the NU evidence
! envelope (automsk=yes background policy, derived in the same evidence
! pass). Part of the frozen-evidence identity via the provenance string.
character(len=32) :: nu_solvent_clamp_source = 'density_envelope'
! Evidence-envelope null. Spherical base pair: margin median/MAD over the observed support (no shell).
! Envelope-constrained pair (set_nu_evidence_null_shell): null on the density envelope's dilation
! ring at full base-support weight; labels free only on the observed density envelope.
logical, allocatable :: nu_calib_lmask(:) !< labels free here; fixed solvent elsewhere
logical, allocatable :: nu_null_lmask(:)  !< null statistics estimated here (the Euclidean shell)
integer :: n_nu_calib = 0
integer :: n_nu_null  = 0
integer :: nu_bank_cap_find = 0 !< FSC-anchored static-bank cap in Fourier shells; 0 = uncapped (no fsc_res)
type(image),      allocatable :: aux_even_bank(:), aux_odd_bank(:)
integer :: ldim(3), box
integer :: n_nu_mask = 0
integer :: nu_smooth_norm_radius = -1
real    :: smpd, nu_support_mskdiam = 0.
logical :: nu_l_report = .true.
! Opt-in diagnostics for NU-filter development. Keep normal runs concise; this
! flag restores the detailed candidate, smoothing, and continuity logs.
logical :: NU_DEV_OUTPUT = .false.
! Radial raw E/O noise profile that whitens the Huber unary (image::nu_objective_noise_profile),
! computed once per setup_nu_dmats and reused by build_nu_evidence_state for its null candidate.
real, allocatable :: nu_noise_profile_cached(:)
real    :: nu_noise_rmax_cached = 0.
integer :: nu_aux_replacement_label = 0
real    :: nu_aux_replacement_resolution = 0.
logical :: nu_evidence_requested = .false.
character(len=32) :: nu_evidence_source = ''
real(kind=8) :: nu_evidence_source_fingerprint(6) = 0.d0


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
    ! label domain (observed density envelope, or the observed support) and null set
    integer :: n_calib          = 0
    real    :: pct_calib        = 0.  !< domain as a percentage of the support
    real    :: pct_signal_calib = 0.  !< signal as a percentage of the domain
    logical :: l_null_majority  = .true. !< solvent held the majority of the domain (mixture regime validity)
    logical :: l_null_shell     = .false. !< null designated by the Euclidean shell rather than estimated from the mixture
    integer :: n_null           = 0   !< voxels the null statistics were estimated on
    real    :: pct_null         = 0.  !< null set as a percentage of the domain
    real    :: core_med         = 0.  !< median margin inside the core (domain minus shell), separation diagnostic
    real    :: pct_signal_null  = 0.  !< signal as a percentage of the null shell: the dilation ring carrying evidence
    logical :: l_null_valid     = .true. !< the null is trusted: majority (mixture) or sufficient shell (Euclidean)
end type nu_envmask_stats

! Public scalar metadata for a frozen NU evidence state. The packed arrays stay private in
! nu_evidence_state; unpack_nu_evidence_state copies them out, so consumers cannot mutate the state.
type :: nu_evidence_summary
    logical :: valid = .false.
    integer :: ldim(3) = 0
    integer :: n_support = 0
    integer :: n_candidates = 0
    integer :: n_bands = 0
    real    :: smpd = 0.
    real    :: mskdiam = 0.
    real    :: null_fraction = 0.
    real    :: uncertain_fraction = 0.
    real    :: observed_fraction = 1. !< fraction of the packed support carrying a non-zero half-map pair
    real    :: calibration_temperature = 0.
    real    :: spatial_beta = 0.
    real    :: null_cost_mean = 0.
    real    :: null_bias_median = 0.
    real    :: null_bias_mad = 0.
    real    :: null_bias_threshold = 0.
    ! NU_EVIDENCE_NBANDS static entries plus appended bands, up to NU_EVIDENCE_MAX_NBANDS
    real, allocatable :: supported_fraction(:)
    real, allocatable :: band_limits(:)
    character(len=32) :: source = ''
    character(len=16) :: identity = ''
    character(len=XLONGSTRLEN) :: provenance = ''
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

    module subroutine init_nu_filter( vol_even, vol_odd, fsc_res )
        class(image), intent(in) :: vol_even, vol_odd
        real,    optional, intent(in) :: fsc_res
    end subroutine init_nu_filter

    module subroutine set_nu_filter_report( l_report )
        logical, intent(in) :: l_report
    end subroutine set_nu_filter_report

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

    module subroutine set_nu_solvent_envelope( envmask, source )
        class(image),               intent(in) :: envmask
        character(len=*), optional, intent(in) :: source
    end subroutine set_nu_solvent_envelope

    module subroutine clear_nu_solvent_envelope
    end subroutine clear_nu_solvent_envelope

    module subroutine set_nu_evidence_null_shell( envmask, core, dilated, base_support )
        class(image), target, intent(in) :: envmask, core, dilated, base_support
    end subroutine set_nu_evidence_null_shell

    module subroutine clear_nu_evidence_null_shell
    end subroutine clear_nu_evidence_null_shell

    module subroutine apply_nu_solvent_clamp( n_clamped )
        integer, optional, intent(out) :: n_clamped
    end subroutine apply_nu_solvent_clamp

    module subroutine cleanup_aux_bank
    end subroutine cleanup_aux_bank

    module subroutine validate_aux_volumes( aux_even, aux_odd )
        type(image), intent(in) :: aux_even(:), aux_odd(:)
    end subroutine validate_aux_volumes

    module subroutine stash_aux_volumes( aux_even, aux_odd )
        type(image), intent(in) :: aux_even(:), aux_odd(:)
    end subroutine stash_aux_volumes

    module subroutine setup_nu_observed_mask( vol_even, vol_odd )
        class(image), intent(in) :: vol_even, vol_odd
    end subroutine setup_nu_observed_mask

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

    module subroutine pack_nu_raw_candidate( dmat_full, icand )
        real,    intent(in) :: dmat_full(:,:,:)
        integer, intent(in) :: icand
    end subroutine pack_nu_raw_candidate

    module subroutine unpack_nu_raw_candidate( icand, dmat_full )
        integer, intent(in)    :: icand
        real,    intent(inout) :: dmat_full(:,:,:)
    end subroutine unpack_nu_raw_candidate

    module subroutine pack_nu_full_to_mask( dmat_full, packed )
        real, intent(in)    :: dmat_full(:,:,:)
        real, intent(inout) :: packed(:)
    end subroutine pack_nu_full_to_mask

    module subroutine cache_filtered_vols( vol_even, vol_odd )
        class(image), intent(in) :: vol_even, vol_odd
    end subroutine cache_filtered_vols

    ! In submodule: simple_nu_filter_bank.f90
    module subroutine setup_nu_dmats( vol_even, vol_odd, mskdiam, aux_resolutions, aux_even, aux_odd, &
            &evidence_source, fsc_res )
        class(image),          intent(in) :: vol_even, vol_odd
        real,                  intent(in) :: mskdiam
        real,                  intent(in) :: aux_resolutions(:)
        type(image), optional, intent(in) :: aux_even(:), aux_odd(:)
        character(len=*), optional, intent(in) :: evidence_source
        real,        optional, intent(in) :: fsc_res !< pair FSC=0.143 resolution in A; caps the bank
    end subroutine setup_nu_dmats

    module subroutine setup_nu_candidate_coords( n_candidates )
        integer, intent(in) :: n_candidates
    end subroutine setup_nu_candidate_coords

    module integer function nu_static_ladder_count( n_base )
        integer, intent(in) :: n_base
    end function nu_static_ladder_count

    module real function nu_potts_coord_for_label( ilabel, n_base )
        integer, intent(in) :: ilabel, n_base
    end function nu_potts_coord_for_label

    module real function get_nu_filter_bank_finest_lp()
    end function get_nu_filter_bank_finest_lp

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

    module subroutine assert_nu_evidence_replay_ready( state )
        type(nu_evidence_state), intent(in) :: state
    end subroutine assert_nu_evidence_replay_ready

    ! In submodule: simple_nu_filter_sharpen.f90
    module subroutine nu_evidence_sharpen_vol( state, vol_even, vol_odd, fsc, vol_sharp, apply_even, apply_odd )
        type(nu_evidence_state), intent(in)    :: state
        type(image),             intent(in)    :: vol_even, vol_odd
        real,                    intent(in)    :: fsc(:)
        type(image),             intent(inout) :: vol_sharp
        type(image), optional,   intent(in)    :: apply_even, apply_odd
    end subroutine nu_evidence_sharpen_vol

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


    ! In submodule: simple_nu_filter_apply.f90
    !> Filtered even/odd references from the current label field: composed from the cached candidates of
    !! the competition pair or, given vol_apply_even/odd (pcg_solvent=yes), from that pair filtered per label.
    !! The auxiliary label is filled from the auxiliary pair in both cases.
    module subroutine nu_filter_vols( vol_even, vol_odd, vol_apply_even, vol_apply_odd )
        class(image),           intent(out) :: vol_even, vol_odd
        class(image), optional, intent(in)  :: vol_apply_even, vol_apply_odd
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

    module real function get_nu_filtmap_finest_selected_lp( mask, min_assigned_pct, min_signal_pct, n_signal )
        logical, optional, intent(in)  :: mask(:,:,:)
        real,    optional, intent(in)  :: min_assigned_pct
        real,    optional, intent(in)  :: min_signal_pct
        integer, optional, intent(out) :: n_signal
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

    module subroutine write_nu_evidence_envmask( nsigma, lp_smooth, smpd, state, fname, mask_out, l_valid )
        real,              intent(in)  :: nsigma, lp_smooth, smpd
        integer,           intent(in)  :: state
        class(string),     intent(in)  :: fname
        class(image), optional, intent(inout) :: mask_out
        logical,      optional, intent(out)   :: l_valid
    end subroutine write_nu_evidence_envmask

    module subroutine print_nu_envmask_stats( stats )
        type(nu_envmask_stats), intent(in) :: stats
    end subroutine print_nu_envmask_stats

end interface

end module simple_nu_filter
