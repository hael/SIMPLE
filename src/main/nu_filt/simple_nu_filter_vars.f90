!@descr: working state of simple_nu_filter, shared by its submodules (the filter bank, label field and evidence)
module simple_nu_filter_vars
use simple_image, only: image
implicit none
public

integer,          parameter   :: NU_LABEL_KIND     = selected_int_kind(4)
real,             parameter   :: lowpass_limits(8) = [20.,15.,12.,10.,8.,6.,5.,4.]
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
! Radial raw E/O noise profile that whitens the Huber unary (image::nu_objective_noise_profile),
! computed once per setup_nu_dmats and reused by build_nu_evidence_state for its null candidate.
real, allocatable :: nu_noise_profile_cached(:)
real    :: nu_noise_rmax_cached = 0.
integer :: nu_aux_replacement_label = 0
real    :: nu_aux_replacement_resolution = 0.
logical :: nu_evidence_requested = .false.
character(len=32) :: nu_evidence_source = ''
real(kind=8) :: nu_evidence_source_fingerprint(6) = 0.d0

end module simple_nu_filter_vars
