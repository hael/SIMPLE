!@descr: EM estimation of the low-rank flex_pca covariance model: basis fit, latent embedding and mean estimation
module simple_flex_pca_em
use simple_core_module_api
use simple_builder,         only: builder
use simple_image,           only: image
use simple_parameters,      only: parameters
use simple_reconstructor,   only: reconstructor
use simple_kbinterpol,      only: kbinterpol
use simple_gridding,        only: prep3D_inv_kbenvelope4mul
use simple_linalg,          only: jacobi, eigsrt
use simple_math,            only: ceil_div, floor_div
use simple_srch_sort_loc,   only: hpsort
use simple_math_ft,         only: resample_sigma2
use simple_ftiter,          only: ftiter
use simple_flex_pca_polar,  only: polar_grid_t
use simple_flex_pca_crossfsc, only: crossfsc_file, crossfsc_record
use simple_flex_reconstructor_latent_ops, only: pair_index
use simple_flex_pca_pcg,    only: flex_pcg_t, flex_pcg_outcome_t, flex_window_apply, flex_window_apply_rec, flex_pcg_install_window, flex_env_init, flex_env_active
use simple_ori,             only: ori
use simple_flex_pca_rounds, only: flex_pca_rounds, flex_pca_rounds_shmem, flex_pca_half_of, flex_pca_part_fname, &
    &FLEX_PCA_PART_MAGIC, PROBE_PART_VERSION, PROBE_PART_VERSION5, EMBED_STATS_VERSION, &
    &PCA_STAGE_PROBE, PCA_STAGE_POLISH, PCA_STAGE_EMBED, PCA_STAGE_STATES, FLEX_FIT_ALL, FLEX_FIT_A, FLEX_FIT_B
implicit none
private
#include "simple_local_flags.inc"

public :: build_covariance_eigenbasis, embed_latents_with_contrast, estimate_covariance_mean
public :: probe_subspace_iteration, align_basis_to_reference
public :: init_basis_reconstructor
public :: bag_basis_pool, basis_recs_from_images
public :: cov_env_int_pub, save_probe_state
public :: compose_basis_from_runs, compose_cut_reembed
public :: probe_worker_pass, embed_worker_pass
public :: probe_fit_t, new_probe_fit, kill_probe_fit, xfsc_ctx_t
public :: run_flex_pca_paired, run_flex_pca_paired_worker

! Density observability floor, matching simple_image_arith::div_cmat_at_1 and the projected-latent coupled
! solve.
real(dp), parameter :: COV_DENSITY_FLOOR = 1.0d-6
! Resident-volume capacity of the fused device E-step. MUST match the `ncomp1 > 64` guard and
! the u_re/u_im extents in cuda/simple_flex_gpu_kernels.cu; exceeding it returns a hard error.
! Raised 24 -> 64 (2026-08-18): the per-thread arrays are dynamically indexed, hence in local
! memory either way -- 24 was a stack-frame guess, not a register wall. neigs=40 now rides the
! device path; parity/speed vs the CPU E-step is validated per-rank by the phase-12 A/B.
integer,  parameter :: COV_GPU_ESTEP_MAXVOLS = 64
! Relative ridge used ONLY for the covariance diagonal in the S.B SNR proxy, which runs before the S.C
! weights exist (Algorithm 1 precedes Algorithm 2).
real,     parameter :: COV_RIDGE_REL     = 5.0e-2
! Relative eigenvalue floor for retaining direct-column PCA components.
real(dp), parameter :: COV_EIG_REL_FLOOR = 1.0d-6
! Cap on the column-subspace dimension. The accumulation is a batched dsyrk on the Van Loan-Pitsianis
! rearrangement (see unrearrange_kron_selfsum), which needs ONE shared d^4 array regardless of thread count.
integer,  parameter :: COV_MAX_DTILDE    = 320
! Default column-subspace dimension, applied as a min against the memory budget so the rank follows
! the data rather than free RAM. Override with SIMPLE_COV_DTILDE.
integer,  parameter :: COV_DEFAULT_DTILDE = 128
! Total particles the probe / SNR / column-accumulation initialiser fit on, summed across processes.
! 0 = OFF: capping traded a recovered conformation for speed, which is not a trade worth taking.
! See doc/policies/flex_pca_policy.md. Enable per-run with SIMPLE_COV_PROBE_MAX / SIMPLE_COV_BASIS_MAX.
integer,  parameter :: COV_PROBE_MAX_PTCLS = 0
integer,  parameter :: COV_BASIS_MAX_PTCLS = 0
! Particle budget for em_calibrate_noise_prior. That pass estimates exactly TWO global
! scalars (the whitened-noise constant sig2 and the initial prior variance Gamma^0), whose
! precision improves only as 1/sqrt(N) -- 20k particles already give ~0.7% on sig2, far
! tighter than the EM needs. It previously shared SIMPLE_COV_BASIS_MAX, whose default of 0
! means UNCAPPED, so on a full dataset it read every particle on the MASTER ALONE (workers
! idle, no stage guard covers it): measured ~5 min of 1/10-capacity time at 105k/box64
! before the first distributed round. Capping removes the work rather than parallelising it.
integer,  parameter :: COV_CALIB_MAX_PTCLS = 20000
! How far above the spectrum's noise bulk a direction must stand to count as signal. Loose by design:
! keeping a noise direction costs one rank, dropping a real one costs a conformation.
real(dp), parameter :: COV_SIGNAL_FACTOR = 4.0d0
! Samples per free parameter for the rank bound d ~ sqrt(2N/R). REPORT ONLY.
real(dp), parameter :: COV_SAMPLES_PER_PARAM = 10.0d0
!> probe stops when successive bases agree to this mean principal-angle cosine: tight enough that the
!! remaining rotation cannot move a state target, and it fires at the measured knee rather than
!! running the tuned count out.
!> The fit runs its full n_probe_iters budget: an in-loop convergence stop was measured to fire on
!> semi-convergence, and the paired merge needs the last two iterations' frames. Kept as a ceiling
!> that a rank-1 fit can still trip (a lone component that stops moving has nothing left to learn).
real(dp), parameter :: COV_PROBE_CONV    = 0.999999d0
!> mixture width of the MCFA E-step: the deconvolution picks the macro-clusters downstream, so this
!> is only the E-step's flexibility budget (measured stable 8-32 on 10028/10076)
integer,  parameter :: COV_EM_MIX        = 16
!> nuisance shells deflated out of every basis volume each M-step (envelope + background + the
!> consensus dilation): without it ~86% of the 10028 basis was solvent and envelope
integer,  parameter :: COV_EM_DEFLATE    = 4
!> Convergence is declared when the REPRODUCIBLE dimension -- the sum of principal-angle cosines
!! between the even and odd half-bases, which the M-step already produces every iteration -- has
!! not improved by more than COV_EO_TOL for COV_EO_PATIENCE consecutive iterations. This is the
!! only criterion here with no dataset-specific constant in it: it asks the data how much of the
!! basis survives a change of particles, and stops when that stops growing.
real(dp), parameter :: COV_EO_TOL        = 0.02d0
integer,  parameter :: COV_EO_PATIENCE   = 3
! Memory budget for the shared A accumulator, in bytes.
real(dp), parameter :: COV_ATHR_BUDGET   = 8.0d9
! Accumulate the columns against the unscaled mean (a==1) rather than the per-particle ML contrast a_i.
! Subtracting a_i*T*mu also deletes the component of the conformational signal parallel to T*mu.
logical,  parameter :: COV_UNIT_CONTRAST  = .true.
! Grid-search the per-particle contrast in the embedding instead of using the closed-form estimate.
logical,  parameter :: COV_EMBED_CONTRAST_GRID = .false.
integer,  parameter :: COV_CG_MAXIT = 2000     ! CG iteration cap; convergence is reported, not assumed
real(dp), parameter :: COV_CG_TOL   = 1.d-10   ! relative residual target
integer,  parameter :: GRAM_DIAG_STRIDE = 200   ! subsample for the projected-Gram spectrum
integer,  parameter :: NCONTRAST_GRID = 50
real(dp), parameter :: A_GRID_HI = 2.0d0
! bracket for the fitted per-particle contrast; wide enough not to bind on real amplitude spread,
! tight enough that a particle the mask or the band has emptied cannot drag its latent to infinity
real(dp), parameter :: A_CONTRAST_LO = 0.2d0
real(dp), parameter :: A_CONTRAST_HI = 3.0d0
real(dp), parameter :: COV_PINV_RCOND = 1.0d-6

! Source of the covariance mean mu.
logical, parameter :: COV_MEAN_FROM_DATA = .false.

! Soft-mask each particle image to the projected molecular envelope before the column accumulation, so
! solvent (pure noise) does not enter the inner products.
logical, parameter :: COV_MASK_IMAGES = .false.
! Radial margin on that disc, as a multiple of the model radius. It must cover CTF delocalisation,
! lambda*defocus/d, which reaches ~70 A at 5.5 um defocus and 15 A resolution.
real, parameter :: COV_MASK_MARGIN = 1.4

! Subtract the analytic per-sample noise bias K_R(.,q_s)|T|^2 from the column numerator. Without it the
! bias survives into the half-set column FSC and the Wiener shrinkage deletes the low-frequency band.
logical, parameter :: COV_COLUMN_NOISE_DEBIAS = .true.

character(len=*), parameter :: COV_UTILDE_FBODY = 'flex_pca_utilde'
character(len=*), parameter :: COV_UTILDE_META  = 'flex_pca_utilde.txt'
!> master -> probe-worker handoff: the basis dimension, its prior variances and the whitened-noise
!! level. The basis volumes themselves are already on disk as flex_pca_pc*.mrc.
character(len=*), parameter :: COV_PROBE_META   = 'flex_pca_probe.txt'
! Half-width of the KB backprojection stencil in grid units, as cov_kb_weights derives it.
integer, parameter :: COV_KB_IWINSZ = ceiling(KBWINSZ - 0.5)

!> Per-fit state of one probe EM fit (the paired-engine state hoist, proposal §4 step 0/1).
!!
!! probe_subspace_iteration used to keep everything that survives from EM iteration N to N+1 as
!! subroutine-local variables of the one long call -- which is exactly why a distributed probe
!! worker (relaunched with niters=1 every round) destroys its MCFA state, prev_real
!! and polar caches each round, and why two alternating single-fit calls could never
!! implement the paired engine. This type is that state, hoisted, so one master-loop iteration
!! can advance two resident fits.
!!
!! LIFECYCLE CONTRACT (the MCFA free-on-iteration crash, em_iter teardown note, verbatim
!! archetype): the mixture state (mix_*) and its work arrays (rhs0th etc) must survive from one
!! iteration's M-step to the next iteration's E-step. They are freed ONLY by kill_probe_fit or
!! by the basis-rank-change resize at iteration start -- never by per-iteration cleanup. The
!! historical bug: iteration 3 initialised the mixture, per-iteration cleanup freed it,
!! l_mix_active stayed true, and iteration 4's E-step walked into an unallocated mix_Ominv (a
!! null descriptor, invisible to -fcheck=bounds, every thread segfaulting on the same line).
type :: probe_fit_t
    ! ---- identity, selection, files ----
    integer      :: id = 0                  !< 1 = fit A, 2 = fit B, 0 = single-fit (legacy)
    type(string) :: fprefix                 !< eigenvolume namespace (default 'flex_pca_pc')
    type(string) :: meta_fname              !< probe-state file (default COV_PROBE_META)
    integer, allocatable :: pinds(:)        !< this fit's half of the master's selection
    integer      :: nptcls = 0
    integer, allocatable :: ppinds(:)       !< probe-stage subsample; owned per fit (the banked-
                                            !! adjoint path permutes it in place -- see hazards)
    integer      :: npp = 0
    ! ---- model handles ----
    type(reconstructor), allocatable :: basis_recs(:)
    real(dp),            allocatable :: eigvals(:)
    real(dp)     :: sig2_eff = 0.d0         !< calibrated whitened-noise level (raw)
    real(dp)     :: sig2     = 0.d0         !< floored copy the solves consume
    integer      :: ncomp    = 0            !< can shrink at the M-step swap; per-fit
    type(reconstructor)   :: mean_rec       !< per-fit mean copy carrying the per-fit mean scale
                                            !! (paired driver only; the single-fit path keeps the
                                            !! caller's mean_rec argument)
    real(dp), allocatable :: z(:,:)         !< (npp,ncomp at stage entry); never reallocated
                                            !! mid-stage even when ncomp shrinks (rank never grows)
    ! ---- convergence state (cross-iteration) ----
    real(dp)     :: nll_prev = 0.d0
    real(dp)     :: eo_best  = -1.d0
    integer      :: eo_stall = 0, eo_patience = 0
    logical      :: l_converged = .false.
    type(image), allocatable :: prev_real(:)   !< previous iteration's orthonormal basis
    real(dp)     :: conv_thresh = 0.d0
    ! ---- MCFA mixture state (cross-iteration; see the lifecycle contract above) ----
    logical      :: l_mix_req = .false., l_mix_active = .false., l_mix_used = .false.
    integer      :: kmix = 0, n_mix_warm = 3
    real(dp)     :: ldOm_mix = 0.d0, ldOm_used = 0.d0
    real(dp), allocatable :: mix_xi(:,:), mix_Om(:,:), mix_Ominv(:,:), mix_pi(:)
    real(dp), allocatable :: mix_Omxi(:,:), mix_xiOx(:), mix_lpi(:)
    real(dp), allocatable :: rhs0th(:,:), mkth(:,:,:), lwth(:,:), rkth(:,:)
    real(dp), allocatable :: mxa_sr(:,:), mxa_sm(:,:,:), mxa_smm(:,:,:,:), mxa_sainv(:,:,:)
    !> distributed-reduce buffers (phase 2 dispatch; fields live here from the start)
    real(dp), allocatable :: dm_sr(:), dm_sm(:,:), dm_smm(:,:,:), dm_sai(:,:), dm_z(:,:)
    integer      :: dm_nz = 0
    ! ---- polar E-step bank state (per stage, pose-keyed to this fit's ppinds) ----
    logical      :: l_pol_es = .false., l_pol_grid = .false.
    logical      :: l_pol_bank_it = .false., l_pol_hyb = .false., l_rhyb_off = .false.
    integer      :: rhyb_req = 0, osamp_pol = 1
    integer      :: ndir_es = 0, nsamp_es = 0, nsamp2_es = 0, nk_es = 0
    integer      :: ph0_es = 0, pk0_es = 0, hlo_es = 0, hhi_es = 0, klo_es = 0
    integer      :: nyqr_es = 0, nyqb_es = 0, rhyb_es = 0, npos_es = 0
    integer,  allocatable :: hex_es(:), kex_es(:)
    type(polar_grid_t)    :: pg_es
    real,     allocatable :: rmatb_es(:,:,:), nrmb_es(:,:)
    real,     allocatable :: cae(:), sae(:)
    integer,  allocatable :: dir_es(:)
    logical,  allocatable :: dused_es(:)
    ! per-iteration bank (rebuilt from the refreshed basis) + thread scratch (ncomp-shaped)
    real,     allocatable :: UsallE(:,:,:)
    real(dp), allocatable :: CfE(:,:,:), Cm0E(:,:,:), c00E(:,:)
    complex,  allocatable :: UbankE(:,:,:)
    real,     allocatable :: CspE(:,:,:)
    real,     allocatable :: xws_es(:,:), wr_es(:,:), Reb_es(:,:)
    real(dp), allocatable :: wrd_es(:,:)
    real         :: sec_bank = 0.
    ! ---- M-step accumulators and per-batch scratch (per-iteration, per-fit within it) ----
    type(reconstructor), allocatable :: Yeven(:), Yodd(:)
    real,     allocatable :: rho_e(:,:,:,:), rho_o(:,:,:,:)
    real(dp), allocatable :: prior(:)
    real(dp), allocatable :: Gth(:,:,:), Ath(:,:,:), bth(:,:), cth(:,:), zth(:,:)
    real(dp), allocatable :: Ainvth(:,:,:), Acpth(:,:,:)
    real(dp), allocatable :: hth(:,:), nll_thr(:), gam_dbg(:,:)
    real(dp), allocatable :: zbatch(:,:), dens(:,:,:)
    logical,  allocatable :: valid(:), valid_e(:), valid_o(:)
    real(dp), allocatable :: gam_thr(:,:), gam_acc(:), gam_sum(:)
    integer,  allocatable :: nval_thr(:)
    real(dp)     :: nll_tot = 0.d0
    integer      :: nval = 0, npairs = 0, es(3) = 0
    type(fplane_type), allocatable :: basis_fpls(:,:), mean_fpl(:)
    real(dp), allocatable :: sec_proj_thr(:), sec_gram_thr(:)
    integer      :: khi_fit = 0             !< per-fit band, written by the BAND/RANK diagnostic;
                                            !! the future per-fit band schedule's home (unused v1)
    ! ---- band-anneal schedule ----
    integer      :: khi_full = 0, kfr_ann(2) = 0
    real         :: dstep_ann = 0., lp_it = 0.
    ! ---- per-fit copies of process-wide config decisions (self-description) ----
    integer      :: n_probe_cm = 0, nml_plain = 0
    logical      :: l_probe_mls = .false.
    logical      :: l_deflate_mean = .false.
    integer      :: vdfl = 0
    ! ---- cross-fit-FSC (crossfsc) per-fit hooks ----
    ! The driver owns the artifact and all consumer decisions (see xfsc_ctx_t in em_iter);
    ! fit_iter_finish only (a) harvests the writer payloads below when l_xf_harvest is set --
    ! the H profiles from the rho_e/rho_o pair diagonals BEFORE any invtau2 mutates them and
    ! before their deallocation, the internal e/o FSC curves and Gamma at the BAND/RANK site --
    ! and (b) applies xf_invtau2 (built by the driver from the previous iteration's record) to
    ! the diagonal rows of this fit's rho_e AND rho_o immediately before the two coupled solves,
    ! deallocating it after (one-shot per iteration). With no XFSC gate live, l_xf_harvest stays
    ! false and xf_invtau2 unallocated: both hooks are inert.
    logical  :: l_xf_harvest = .false.
    real,     allocatable :: xf_h_e(:,:), xf_h_o(:,:)  !< (filtsz,ncomp) per-shell H, pre-ridge
    integer,  allocatable :: xf_cnt(:)                 !< (filtsz) shared per-shell voxel counts
    real,     allocatable :: xf_fscq(:,:)              !< (filtsz,ncomp) internal e/o FSC curves
    real(dp), allocatable :: xf_gam(:)                 !< (ncomp) Gamma at the writer site
    real,     allocatable :: xf_invtau2(:,:)           !< (ncomp,filtsz) ridge from record t-1
    ! ---- paired-merge stash (SIMPLE_COV_PAIRED_MERGE=1; proposal par.7 step 2) ----
    ! The final stage merges the two fits' LAST-iteration M-step sufficient statistics; those
    ! live only inside one iteration (fit_iter_finish mutates the numerators in the coupled
    ! solve and frees everything), so under the merge gate the paired driver snapshots them
    ! RAW -- after the batch loop + reductions, before fit_iter_finish (pre-ridge, pre-solve)
    ! -- every iteration, overwriting: any iteration can turn out to be the last (convergence
    ! or the crossfsc stop are evaluated after the tails). mg_prev is the fit's ENTRY-frame
    ! orthonormal basis (prev_real at stash time = the previous iteration's delivered basis),
    ! which is the latent frame the stashed statistics are expressed in -- the frame map R is
    ! computed between the two fits' mg_prev sets, NOT between the post-solve delivered bases.
    logical  :: l_mg_stash = .false.               !< stash present (driver-gated)
    integer  :: mg_ncomp = 0, mg_npairs = 0        !< entry rank of the stashed iteration
    complex,  allocatable :: mg_ye(:,:,:,:), mg_yo(:,:,:,:)   !< (es) x ncomp numerators
    real,     allocatable :: mg_rhe(:,:,:,:), mg_rho(:,:,:,:) !< (npairs, es) packed densities
    real,     allocatable :: mg_kpe(:,:), mg_kpo(:,:)     !< (npairs, npk) packed PCG pair kernels on the band list (rec_backend=pcg)
    complex,  allocatable :: mg_rpe(:,:), mg_rpo(:,:)     !< (ncomp, npk) packed PCG right-hand sides on the band list (rec_backend=pcg)
    ! ---- rec_backend=pcg: the coupled M-step on the PCG operator (simple_flex_pca_pcg) ----
    logical  :: l_pcg = .false.                    !< this fit solves its M-step by PCG
    type(flex_pcg_t) :: pcg                        !< lattice, envelopes, support, finalized kernels
    real,     allocatable :: kacc_e(:,:), kacc_o(:,:)     !< expanded-lattice accumulators on the band list (inserting processes)
    real,     allocatable :: kpk_e(:,:), kpk_o(:,:)       !< packed pair kernels on the physical band list: transport, reduction and merge form
    complex,  allocatable :: racc_e(:,:), racc_o(:,:)     !< right-hand-side accumulators on the band list
    complex,  allocatable :: rpk_e(:,:), rpk_o(:,:)       !< packed right-hand sides on the physical band list
    type(image), allocatable :: mg_prev(:)         !< entry-frame orthonormal basis images
end type probe_fit_t


!> ---- CROSS-FIT-FSC driver context (artifact writer + SSNR ridge) ----
!! Spec: doc/for_developers/ideas/flex_pca_crossfsc_shrinkage_marching_spec.md (the marching and
!! stopping consumers of that spec were removed 2026-09-07; only the artifact and the ridge remain).
!! The ridge defaults OFF (SIMPLE_COV_XFSC_REG=0). One context is owned by each driver loop
!! (single-fit probe_subspace_iteration, paired probe_subspace_paired); the per-fit payloads the
!! writer needs (H profiles harvested from the rho pair diagonals, internal FSC curves, Gamma) are
!! stashed in probe_fit_t by fit_iter_finish, so the phase procedure stays ignorant of the artifact
!! and the driver owns every consumer decision. Only the paired engine writes records (paired=1,
!! one per iteration: the artifact is its lasting product); the single-fit engine can consume a
!! paired artifact left in the directory (pcafit names its side).
type :: xfsc_ctx_t
    logical  :: l_writer   = .false.   !< this driver writes records (paired master)
    logical  :: l_any      = .false.   !< ridge or writer live -> artifact state resident
    logical  :: l_loaded   = .false.   !< artifact reloaded from disk (restart-complete series)
    logical  :: l_paired   = .false.   !< this driver is the paired engine
    integer  :: v_reg      = 0         !< par.4 arms: 0 internal / 1 cross / 2 blend
    integer  :: reg_active = 0         !< arm ACTIVE this iteration (0 when degraded)
    integer  :: klo        = 6         !< low-resolution exemption index (reslim_ind analog)
    integer  :: filtsz     = 0
    integer  :: pairing_id = 0         !< balanced mod-4 pairing id (paired engine; header field)
    type(crossfsc_file) :: xf
end type xfsc_ctx_t


interface

    module subroutine write_probe_part( fname, cmat_e, rho_ex, cmat_o, rho_ox, rho_e, rho_o, gam_sum, nll_sum, nval, ncomp, kpk_e, kpk_o, rpk_e, rpk_o, mix_sr, mix_sm, mix_smm, mix_sainv, z_sub )
        class(string), intent(in) :: fname
        complex,       intent(in) :: cmat_e(:,:,:,:), cmat_o(:,:,:,:)
        real,          intent(in) :: rho_ex(:,:,:,:), rho_ox(:,:,:,:)
        real,          intent(in) :: rho_e(:,:,:,:),  rho_o(:,:,:,:)
        real(dp),      intent(in) :: gam_sum(:)
        real(dp),      intent(in) :: nll_sum
        integer,       intent(in) :: nval, ncomp
        real, allocatable, intent(in) :: kpk_e(:,:), kpk_o(:,:)
        complex, allocatable, intent(in) :: rpk_e(:,:), rpk_o(:,:)
        real(dp), optional, intent(in) :: mix_sr(:), mix_sm(:,:), mix_smm(:,:,:), mix_sainv(:,:)
        real(dp), optional, intent(in) :: z_sub(:,:)
    end subroutine write_probe_part

    module subroutine reduce_probe_parts( params, nparts, cmat_e, rho_ex, cmat_o, rho_ox, rho_e, rho_o, gam_sum, nll_sum, nval, ncomp, kpk_e, kpk_o, rpk_e, rpk_o, mix_sr, mix_sm, mix_smm, mix_sainv, z_pool, nz_pool )
        class(parameters), intent(in)    :: params
        integer,           intent(in)    :: nparts, ncomp
        complex,           intent(inout) :: cmat_e(:,:,:,:), cmat_o(:,:,:,:)
        real,              intent(inout) :: rho_ex(:,:,:,:), rho_ox(:,:,:,:)
        real,              intent(inout) :: rho_e(:,:,:,:),  rho_o(:,:,:,:)
        real(dp),          intent(inout) :: gam_sum(:)
        real(dp),          intent(inout) :: nll_sum
        integer,           intent(inout) :: nval
        real, allocatable, intent(inout) :: kpk_e(:,:), kpk_o(:,:)
        complex, allocatable, intent(inout) :: rpk_e(:,:), rpk_o(:,:)
        real(dp), optional, intent(inout) :: mix_sr(:), mix_sm(:,:), mix_smm(:,:,:), mix_sainv(:,:)
        real(dp), optional, intent(inout) :: z_pool(:,:)
        integer,  optional, intent(inout) :: nz_pool
    end subroutine reduce_probe_parts

    module subroutine open_probe_part_v5_write( fname, nfits, funit, tmp_fname )
        class(string), intent(in)  :: fname
        integer,       intent(in)  :: nfits
        integer,       intent(out) :: funit
        type(string),  intent(out) :: tmp_fname
    end subroutine open_probe_part_v5_write

    module subroutine write_probe_part_v5_fit( funit, cmat_e, rho_ex, cmat_o, rho_ox, rho_e, rho_o, gam_sum, nll_sum, nval, ncomp, kpk_e, kpk_o, rpk_e, rpk_o, mix_sr, mix_sm, mix_smm, mix_sainv, z_sub )
        integer,       intent(in) :: funit
        complex,       intent(in) :: cmat_e(:,:,:,:), cmat_o(:,:,:,:)
        real,          intent(in) :: rho_ex(:,:,:,:), rho_ox(:,:,:,:)
        real,          intent(in) :: rho_e(:,:,:,:),  rho_o(:,:,:,:)
        real(dp),      intent(in) :: gam_sum(:)
        real(dp),      intent(in) :: nll_sum
        integer,       intent(in) :: nval, ncomp
        real, allocatable, intent(in) :: kpk_e(:,:), kpk_o(:,:)
        complex, allocatable, intent(in) :: rpk_e(:,:), rpk_o(:,:)
        real(dp), optional, intent(in) :: mix_sr(:), mix_sm(:,:), mix_smm(:,:,:), mix_sainv(:,:)
        real(dp), optional, intent(in) :: z_sub(:,:)
    end subroutine write_probe_part_v5_fit

    module subroutine close_probe_part_v5_write( funit, tmp_fname, fname )
        integer,       intent(in)    :: funit
        type(string),  intent(inout) :: tmp_fname
        class(string), intent(in)    :: fname
    end subroutine close_probe_part_v5_write

    module subroutine open_probe_part_v5_read( fname, nfits, funit )
        class(string), intent(in)  :: fname
        integer,       intent(in)  :: nfits
        integer,       intent(out) :: funit
    end subroutine open_probe_part_v5_read

    module subroutine fold_probe_part_v5_fit( funit, cmat_e, rho_ex, cmat_o, rho_ox, rho_e, rho_o, gam_sum, nll_sum, nval, ncomp, kpk_e, kpk_o, rpk_e, rpk_o, mix_sr, mix_sm, mix_smm, mix_sainv, z_pool, nz_pool )
        integer,           intent(in)    :: funit
        complex,           intent(inout) :: cmat_e(:,:,:,:), cmat_o(:,:,:,:)
        real,              intent(inout) :: rho_ex(:,:,:,:), rho_ox(:,:,:,:)
        real,              intent(inout) :: rho_e(:,:,:,:),  rho_o(:,:,:,:)
        real(dp),          intent(inout) :: gam_sum(:)
        real(dp),          intent(inout) :: nll_sum
        integer,           intent(inout) :: nval
        integer,           intent(in)    :: ncomp
        real, allocatable, intent(inout) :: kpk_e(:,:), kpk_o(:,:)
        complex, allocatable, intent(inout) :: rpk_e(:,:), rpk_o(:,:)
        real(dp), optional, intent(inout) :: mix_sr(:), mix_sm(:,:), mix_smm(:,:,:), mix_sainv(:,:)
        real(dp), optional, intent(inout) :: z_pool(:,:)
        integer,  optional, intent(inout) :: nz_pool
    end subroutine fold_probe_part_v5_fit

    module subroutine close_probe_part_v5_read( funit, fname )
        integer,       intent(in) :: funit
        class(string), intent(in) :: fname
    end subroutine close_probe_part_v5_read

    module subroutine fit_iter_begin( params, build, fit, mean_rec, it_eff, niters_eff, nthr )
        class(parameters),   intent(inout) :: params
        type(builder),       intent(inout) :: build
        type(probe_fit_t),   intent(inout) :: fit
        type(reconstructor), intent(inout) :: mean_rec
        integer,             intent(in)    :: it_eff, niters_eff, nthr
    end subroutine fit_iter_begin

    module subroutine fit_polar_bank_build( params, build, fit, mean_rec, fpl1, nthr, l_dev_geom )
        class(parameters),   intent(inout) :: params
        type(builder),       intent(inout) :: build
        type(probe_fit_t),   intent(inout) :: fit
        type(reconstructor), intent(inout) :: mean_rec
        type(fplane_type),   intent(in)    :: fpl1
        integer,             intent(in)    :: nthr
        logical,             intent(in)    :: l_dev_geom
    end subroutine fit_polar_bank_build

    module subroutine fit_estep_former_polar( fit, mean_rec, o, fpl, row, ithr, a, aa, e_mm, myv )
        type(probe_fit_t),   intent(inout) :: fit
        type(reconstructor), intent(inout) :: mean_rec
        class(ori),          intent(inout) :: o
        type(fplane_type),   intent(inout) :: fpl
        integer,             intent(in)    :: row, ithr
        real(dp),            intent(out)   :: a, aa, e_mm, myv
    end subroutine fit_estep_former_polar

    module subroutine fit_estep_former_cart( fit, mean_rec, o, fpl, ithr, a, aa, e_mm, myv )
        type(probe_fit_t),   intent(inout) :: fit
        type(reconstructor), intent(inout) :: mean_rec
        class(ori),          intent(inout) :: o
        type(fplane_type),   intent(inout) :: fpl
        integer,             intent(in)    :: ithr
        real(dp),            intent(out)   :: a, aa, e_mm, myv
    end subroutine fit_estep_former_cart

    ! ===== implemented in simple_flex_pca_em_compose =====
    module subroutine compose_basis_from_runs( params, build, basis_recs, eigvals, ncomp, sig2_eff )
        class(parameters),   intent(inout) :: params
        type(builder),       intent(inout) :: build
        type(reconstructor), allocatable, intent(out) :: basis_recs(:)
        real(dp),            allocatable, intent(out) :: eigvals(:)
        integer,             intent(out)   :: ncomp
        real(dp),            intent(out)   :: sig2_eff
    end subroutine compose_basis_from_runs

    module subroutine compose_cut_reembed( params, build, mean_rec, basis_recs, eigvals, ncomp, sig2_eff, &
        &pinds, nptcls, rounds )
        class(flex_pca_rounds), intent(inout) :: rounds
        class(parameters),   intent(inout) :: params
        type(builder),       intent(inout) :: build
        type(reconstructor), intent(inout) :: mean_rec
        type(reconstructor), allocatable, intent(inout) :: basis_recs(:)
        real(dp),            allocatable, intent(inout) :: eigvals(:)
        integer,             intent(inout) :: ncomp
        real(dp),            intent(in)    :: sig2_eff
        integer,             intent(in)    :: pinds(:), nptcls
    end subroutine compose_cut_reembed

    module subroutine fit_estep_solve_stats( fit, fpl, i, row, ithr, a, aa, e_mm, myv )
        type(probe_fit_t), intent(inout) :: fit
        type(fplane_type), intent(inout) :: fpl
        integer,           intent(in)    :: i, row, ithr
        real(dp),          intent(inout) :: a, aa, e_mm, myv
    end subroutine fit_estep_solve_stats

    module subroutine fit_batch_insert( build, fit, orientations, fpls, eo, batchsz )
        type(builder),     intent(inout) :: build
        type(probe_fit_t), intent(inout) :: fit
        type(ori),         intent(inout) :: orientations(:)
        type(fplane_type), intent(inout) :: fpls(:)
        integer,           intent(in)    :: eo(:), batchsz
    end subroutine fit_batch_insert

    module subroutine fit_iter_reduce( fit, it_eff, nthr , rounds)
        class(flex_pca_rounds), intent(inout) :: rounds
        type(probe_fit_t), intent(inout) :: fit
        integer,           intent(in)    :: it_eff, nthr
    end subroutine fit_iter_reduce

    module subroutine fit_iter_finish( params, build, fit, it_eff, nthr )
        class(parameters),  intent(inout) :: params
        type(builder),      intent(inout) :: build
        type(probe_fit_t),  intent(inout) :: fit
        integer,            intent(in)    :: it_eff, nthr
    end subroutine fit_iter_finish

    module subroutine xfsc_setup( ctx, params, kfr_ann, l_paired, l_master )
        type(xfsc_ctx_t),  intent(inout) :: ctx
        class(parameters), intent(in)    :: params
        integer,           intent(in)    :: kfr_ann(2)
        logical,           intent(in)    :: l_paired, l_master
    end subroutine xfsc_setup

    module subroutine xfsc_prep_iter( ctx, params, fit, it_eff, tag )
        type(xfsc_ctx_t),  intent(inout) :: ctx
        class(parameters), intent(in)    :: params
        type(probe_fit_t), intent(inout) :: fit
        integer,           intent(in)    :: it_eff
        character(len=*),  intent(in)    :: tag
    end subroutine xfsc_prep_iter

    module subroutine xfsc_teardown( ctx )
        type(xfsc_ctx_t), intent(inout) :: ctx
    end subroutine xfsc_teardown

    module subroutine xfsc_paired_record( ctx, params, fits, it_eff )
        type(xfsc_ctx_t),  intent(inout) :: ctx
        class(parameters), intent(in)    :: params
        type(probe_fit_t), intent(inout) :: fits(2)
        integer,           intent(in)    :: it_eff
    end subroutine xfsc_paired_record

    module subroutine paired_estep_pass( params, build, fits, it_eff, nthr )
        class(parameters), intent(inout) :: params
        type(builder),     intent(inout) :: build
        type(probe_fit_t), intent(inout) :: fits(2)
        integer,           intent(in)    :: it_eff, nthr
    end subroutine paired_estep_pass

    module subroutine paired_reduce_parts_v5( params, fits , rounds)
        class(flex_pca_rounds), intent(inout) :: rounds
        class(parameters), intent(in)    :: params
        type(probe_fit_t), intent(inout) :: fits(2)
    end subroutine paired_reduce_parts_v5

    ! ===== implemented in simple_flex_pca_em_env =====

    real(dp) module function corr_dp( a, b, n ) result( r )
        integer,  intent(in) :: n
        real(dp), intent(in) :: a(n), b(n)
    end function corr_dp

    logical module function cov_env_int_off( name ) result(off)
        character(len=*), intent(in) :: name
    end function cov_env_int_off

    logical module function cov_env_int_on( name ) result(on)
        character(len=*), intent(in) :: name
    end function cov_env_int_on

    pure integer module function cov_signal_rank( eval, n ) result( d )
        integer,  intent(in) :: n
        real(dp), intent(in) :: eval(n)          !< DESCENDING eigenvalues
    end function cov_signal_rank

    module subroutine cov_stage_subsample( build, pinds, nptcls, nparts, maxtot, env_max, &
        &label, spinds, nsel )
        type(builder),        intent(inout) :: build
        integer,              intent(in)    :: pinds(:), nptcls, nparts, maxtot
        character(len=*),     intent(in)    :: env_max, label
        integer, allocatable, intent(out)   :: spinds(:)
        integer,              intent(out)   :: nsel
    end subroutine cov_stage_subsample

    module subroutine cov_env_int_pub( name, val )
        character(len=*), intent(in)    :: name
        integer,          intent(inout) :: val
    end subroutine cov_env_int_pub

    logical module function cov_env_is_set( name )
        character(len=*), intent(in) :: name
    end function cov_env_is_set

    module subroutine cov_env_int( name, val )
        character(len=*), intent(in)    :: name
        integer,          intent(inout) :: val
    end subroutine cov_env_int

    pure real(dp) module function cov_accum_bytes( d ) result( nbytes )
        integer, intent(in) :: d
    end function cov_accum_bytes

    pure integer module function cov_dim_budget() result( d )
    end function cov_dim_budget

    module subroutine cov_dev_prep_start( params, build, l_on )
        class(parameters), intent(in)    :: params
        class(builder),    intent(inout) :: build
        logical,           intent(out)   :: l_on
    end subroutine cov_dev_prep_start

    module subroutine cov_dev_prep_stop( l_on )
        logical, intent(in) :: l_on
    end subroutine cov_dev_prep_stop

    module subroutine map_sampling_precision( Gtil, prior, n, Qout )
        integer,  intent(in)  :: n
        real(dp), intent(in)  :: Gtil(n,n), prior(n)
        real(dp), intent(out) :: Qout(n,n)
    end subroutine map_sampling_precision

    ! ===== implemented in simple_flex_pca_em_state =====

    module subroutine new_probe_fit( fit, id, fprefix, meta_fname, pinds, nptcls )
        type(probe_fit_t), intent(inout) :: fit
        integer,           intent(in)    :: id
        character(len=*),  intent(in)    :: fprefix, meta_fname
        integer,           intent(in)    :: pinds(:)
        integer,           intent(in)    :: nptcls
    end subroutine new_probe_fit

    module subroutine kill_probe_fit( fit )
        type(probe_fit_t), intent(inout) :: fit
    end subroutine kill_probe_fit

    ! ===== implemented in simple_flex_pca_em_pairmerge =====

    module subroutine probe_fit_merge_stash( fit )
        type(probe_fit_t), intent(inout) :: fit
    end subroutine probe_fit_merge_stash

    module subroutine probe_paired_merge( params, build, fits, merge_mode, m_basis_recs, &
        &m_eigvals, m_ncomp, m_sig2, m_matchcos )
        class(parameters),   intent(inout) :: params
        type(builder),       intent(inout) :: build
        type(probe_fit_t),   intent(inout) :: fits(2)
        integer,             intent(in)    :: merge_mode
        type(reconstructor), allocatable, intent(out) :: m_basis_recs(:)
        real(dp),            allocatable, intent(out) :: m_eigvals(:)
        integer,             intent(out)              :: m_ncomp
        real(dp),            intent(out)              :: m_sig2
        real(dp), allocatable, optional, intent(out)  :: m_matchcos(:)
    end subroutine probe_paired_merge

    module subroutine run_flex_pca_paired_worker( params, build, pinds, nptcls, it_stamp, &
        &niters_stamp, vpair , rounds)
        class(flex_pca_rounds), intent(inout) :: rounds
        class(parameters), intent(inout) :: params
        type(builder),     intent(inout) :: build
        integer,           intent(in)    :: pinds(:), nptcls, it_stamp, niters_stamp, vpair
    end subroutine run_flex_pca_paired_worker

    ! ===== implemented in simple_flex_pca_em_fit =====

    module subroutine build_covariance_eigenbasis( params, build, mean_rec, pinds, nptcls, &
        &col_sep, neigs_req, basis_recs, eigvals, ncomp_out, sig2_out, fprefix, rounds)
        class(flex_pca_rounds), intent(inout) :: rounds
        class(parameters),   intent(inout) :: params
        type(builder),       intent(inout) :: build
        type(reconstructor), intent(inout) :: mean_rec
        integer,             intent(in)    :: pinds(:), nptcls, col_sep, neigs_req
        type(reconstructor), allocatable, intent(out) :: basis_recs(:)
        real(dp),            allocatable, intent(out) :: eigvals(:)
        integer,             intent(out)   :: ncomp_out
        real(dp),            intent(out)   :: sig2_out
        character(len=*),        optional, intent(in)  :: fprefix
    end subroutine build_covariance_eigenbasis

    module subroutine init_basis_datafree( params, build, mean_rec, pinds, nptcls, col_sep, neigs_req, &
        &basis_recs, eigvals, ncomp_out, sig2_out, fprefix , rounds)
        class(flex_pca_rounds), intent(inout) :: rounds
        class(parameters),   intent(inout) :: params
        type(builder),       intent(inout) :: build
        type(reconstructor), intent(inout) :: mean_rec
        integer,             intent(in)    :: pinds(:), nptcls, col_sep, neigs_req
        type(reconstructor), allocatable, intent(out) :: basis_recs(:)
        real(dp),            allocatable, intent(out) :: eigvals(:)
        integer,             intent(out)   :: ncomp_out
        real(dp),            intent(out)   :: sig2_out
        !> per-fit eigenvolume namespace (paired engine); default 'flex_pca_pc'
        character(len=*), optional, intent(in) :: fprefix
    end subroutine init_basis_datafree

    module subroutine em_calibrate_noise_prior( params, build, mean_rec, basis_recs, ncomp, pinds, nptcls, &
        &sig2_out, gam0 )
        class(parameters),   intent(inout) :: params
        type(builder),       intent(inout) :: build
        type(reconstructor), intent(inout) :: mean_rec
        type(reconstructor), intent(inout) :: basis_recs(:)
        integer,             intent(in)    :: ncomp, pinds(:), nptcls
        real(dp),            intent(out)   :: sig2_out, gam0
    end subroutine em_calibrate_noise_prior

    module subroutine save_probe_state( ncomp, eigvals, sig2_eff, fname )
        integer,  intent(in) :: ncomp
        real(dp), intent(in) :: eigvals(:), sig2_eff
        character(len=*), optional, intent(in) :: fname
    end subroutine save_probe_state

    module subroutine load_probe_state( ncomp, eigvals, sig2_eff, fname )
        integer,               intent(out) :: ncomp
        real(dp), allocatable, intent(out) :: eigvals(:)
        real(dp),              intent(out) :: sig2_eff
        character(len=*), optional, intent(in) :: fname
    end subroutine load_probe_state

    !> Distributed worker stage bodies. One E-step pass over the worker's particle list against
    !! the basis the master left on disk (PROBE/POLISH: single or paired, chosen by params%nfits,
    !! writes this part's v5 file), or one embedding-statistics pass (EMBED). The round state
    !! (stage, which_iter, maxits, nfits) arrives in params from job_descr.
    module subroutine probe_worker_pass( params, build, mean_rec, pinds, nptcls , rounds)
        class(flex_pca_rounds), intent(inout) :: rounds
        class(parameters),   intent(inout) :: params
        type(builder),       intent(inout) :: build
        type(reconstructor), intent(inout) :: mean_rec
        integer,             intent(in)    :: pinds(:), nptcls
    end subroutine probe_worker_pass

    module subroutine embed_worker_pass( params, build, mean_rec, pinds, nptcls , rounds)
        class(flex_pca_rounds), intent(inout) :: rounds
        class(parameters),   intent(inout) :: params
        type(builder),       intent(inout) :: build
        type(reconstructor), intent(inout) :: mean_rec
        integer,             intent(in)    :: pinds(:), nptcls
    end subroutine embed_worker_pass

    module subroutine load_probe_basis( params, build, ncomp, basis_recs, fprefix )
        class(parameters),   intent(inout) :: params
        type(builder),       intent(inout) :: build
        integer,             intent(in)    :: ncomp
        type(reconstructor), allocatable, intent(out) :: basis_recs(:)
        character(len=*), optional, intent(in) :: fprefix
    end subroutine load_probe_basis

    module subroutine select_frequencies_lowfreq( params, ncols_req, col_sep, col_hkl, ncol )
        class(parameters),    intent(in)  :: params
        integer,              intent(in)  :: ncols_req, col_sep
        integer, allocatable, intent(out) :: col_hkl(:,:)
        integer,              intent(out) :: ncol
    end subroutine select_frequencies_lowfreq

    module subroutine pick_next_lowfreq( cand, ncand, chosen, nchosen, sep, best )
        integer, intent(in)  :: cand(:,:), ncand, chosen(:,:), nchosen, sep
        integer, intent(out) :: best
    end subroutine pick_next_lowfreq

    module function covariance_kfromto( params ) result( kfromto )
        class(parameters), intent(in) :: params
        integer :: kfromto(2)
    end function covariance_kfromto

    module subroutine init_basis_reconstructor( params, build, rec )
        class(parameters),   intent(inout) :: params
        type(builder),       intent(inout) :: build
        type(reconstructor), intent(inout) :: rec
    end subroutine init_basis_reconstructor

    module function cov_image_mask_radius( params ) result( r )
        class(parameters), intent(in) :: params
        real :: r
    end function cov_image_mask_radius

    ! ===== implemented in simple_flex_pca_em_mean =====

    module subroutine estimate_covariance_mean( params, build, mean_rec, pinds, nptcls , rounds)
        class(flex_pca_rounds), intent(inout) :: rounds
        class(parameters),   intent(inout) :: params
        type(builder),       intent(inout) :: build
        type(reconstructor), intent(inout) :: mean_rec
        integer,             intent(in)    :: pinds(:), nptcls
    end subroutine estimate_covariance_mean

    module subroutine estimate_mean_from_data( params, build, mean_rec, pinds, nptcls )
        class(parameters),   intent(inout) :: params
        type(builder),       intent(inout) :: build
        type(reconstructor), intent(inout) :: mean_rec
        integer,             intent(in)    :: pinds(:), nptcls
    end subroutine estimate_mean_from_data

    module subroutine apply_cached_mean_scale( params, mean_rec, cache_fname )
        class(parameters),   intent(inout) :: params
        type(reconstructor), intent(inout) :: mean_rec
        !> per-fit namespace (paired engine); default flex_pca_mean_scale.bin
        character(len=*), optional, intent(in) :: cache_fname
    end subroutine apply_cached_mean_scale

    module subroutine estimate_mean_scale( params, build, mean_rec, pinds, nptcls, cache_fname , rounds)
        class(flex_pca_rounds), intent(inout) :: rounds
        class(parameters),   intent(inout) :: params
        type(builder),       intent(inout) :: build
        type(reconstructor), intent(inout) :: mean_rec
        integer,             intent(in)    :: pinds(:), nptcls
        !> per-fit cache namespace (paired distributed master writes one per fit)
        character(len=*), optional, intent(in) :: cache_fname
    end subroutine estimate_mean_scale

    module subroutine plane_shell_cross_accum( mean_fpl, fpl, nyq, my_sh, mm_sh )
        type(fplane_type), intent(in)    :: mean_fpl, fpl
        integer,           intent(in)    :: nyq
        real(dp),          intent(inout) :: my_sh(0:), mm_sh(0:)
    end subroutine plane_shell_cross_accum

    module subroutine init_mean_reconstructor( params, build, mean_rec )
        class(parameters),  intent(inout) :: params
        type(builder),      intent(inout) :: build
        type(reconstructor),intent(inout) :: mean_rec
    end subroutine init_mean_reconstructor

    module subroutine form_reconstruction_plane( fpl, num )
        type(fplane_type), intent(in)    :: fpl
        type(fplane_type), intent(inout) :: num
    end subroutine form_reconstruction_plane

    module subroutine cov_herm_sample_list( ref, slist, nsamp )
        type(fplane_type), intent(in)  :: ref
        integer,           intent(out) :: slist(:,:)
        integer,           intent(out) :: nsamp
    end subroutine cov_herm_sample_list

    module function cov_herm_inner( lhs, rhs, half ) result( val )
        type(fplane_type), intent(in) :: lhs, rhs
        integer, optional, intent(in) :: half
        complex(dp) :: val
    end function cov_herm_inner

    real module function particle_contrast( mean_fpl, fpl )
        type(fplane_type), intent(in) :: mean_fpl, fpl
    end function particle_contrast

    module subroutine cov_herm_selfpower( fpl, pw, cnt )
        type(fplane_type), intent(in)  :: fpl
        real(dp),          intent(out) :: pw, cnt
    end subroutine cov_herm_selfpower

    module subroutine plane_hf_power( fpl, nyq, frac, pw, cnt )
        type(fplane_type), intent(in)  :: fpl
        integer,           intent(in)  :: nyq
        real,              intent(in)  :: frac
        real(dp),          intent(out) :: pw, cnt
    end subroutine plane_hf_power

    module subroutine plane_shell_power_accum( fpl, nyq, pw_sh, cnt_sh )
        type(fplane_type), intent(in)    :: fpl
        integer,           intent(in)    :: nyq
        real(dp),          intent(inout) :: pw_sh(0:), cnt_sh(0:)
    end subroutine plane_shell_power_accum

    module subroutine form_residual_plane( fpl, mean_fpl, resid, contrast )
        type(fplane_type), intent(in)    :: fpl, mean_fpl
        type(fplane_type), intent(inout) :: resid
        real, optional,    intent(in)    :: contrast
    end subroutine form_residual_plane

    module subroutine cleanup_plane( fpl )
        type(fplane_type), intent(inout) :: fpl
    end subroutine cleanup_plane

    ! ===== implemented in simple_flex_pca_em_basis =====

    module subroutine basis_to_real_representatives( params, work, colvol, ncol, lb, ub, realvols, nreal )
        class(parameters),   intent(inout) :: params
        type(reconstructor), intent(inout) :: work
        complex,             intent(in)    :: colvol(:,:,:,:)
        integer,             intent(in)    :: ncol, lb(3), ub(3)
        type(image), allocatable, intent(out) :: realvols(:)
        integer,                  intent(out) :: nreal
    end subroutine basis_to_real_representatives

    module subroutine realize_hermitian_volume( params, work, vherm, gridcorr_img, energy )
        class(parameters),   intent(in)    :: params
        type(reconstructor), intent(inout) :: work
        complex,             intent(in)    :: vherm(:,:,:)
        type(image),         intent(inout) :: gridcorr_img
        real,                intent(out)   :: energy
    end subroutine realize_hermitian_volume

    module subroutine orthonormalize_representatives( params, build, realvols, nreal, utilde, utilde_real, d_tilde, svals, &
        &nptcls_basis )
        class(parameters),   intent(inout) :: params
        type(builder),       intent(inout) :: build
        type(image),         intent(inout) :: realvols(:)
        integer,             intent(in)    :: nreal
        type(reconstructor), allocatable, intent(out) :: utilde(:)
        type(image),         allocatable, intent(out) :: utilde_real(:)
        integer,             intent(out)   :: d_tilde
        real(dp), allocatable, optional, intent(out) :: svals(:)
        integer, optional, intent(in) :: nptcls_basis
    end subroutine orthonormalize_representatives

    module subroutine align_basis_to_reference( ref_imgs, nref_c, tgt_imgs, ntgt_c, M, svals )
        integer,     intent(in)    :: nref_c, ntgt_c
        type(image), intent(inout) :: ref_imgs(nref_c), tgt_imgs(ntgt_c)
        real(dp), allocatable, intent(out) :: M(:,:), svals(:)
    end subroutine align_basis_to_reference

    module subroutine bag_basis_pool( imgs_a, na_c, eig_a, imgs_b, nb_c, eig_b, ncomp_out, pooled, eig_pooled )
        integer,     intent(in)    :: na_c, nb_c, ncomp_out
        type(image), intent(inout) :: imgs_a(na_c), imgs_b(nb_c)
        real(dp),    intent(in)    :: eig_a(na_c), eig_b(nb_c)
        type(image), allocatable, intent(out) :: pooled(:)
        real(dp),    allocatable, intent(out) :: eig_pooled(:)
    end subroutine bag_basis_pool

    module subroutine basis_recs_from_images( params, build, imgs, ncomp, basis_recs )
        class(parameters),   intent(inout) :: params
        type(builder),       intent(inout) :: build
        integer,             intent(in)    :: ncomp
        type(image),         intent(inout) :: imgs(ncomp)
        type(reconstructor), allocatable, intent(out) :: basis_recs(:)
    end subroutine basis_recs_from_images

    ! ===== implemented in simple_flex_pca_em_pose =====

    logical module function cov_polar_enabled()
    end function cov_polar_enabled


    integer module function cov_polar_ndir( nptcls )
        integer, intent(in) :: nptcls
    end function cov_polar_ndir










    module subroutine polar_ring_gram( Us, ldu, ncomp, row0, nrow, Csp, Cout, Mout )
        integer,  intent(in)    :: ldu, ncomp, row0, nrow
        real,     intent(in)    :: Us(ldu,0:ncomp)
        real,     intent(inout) :: Csp(0:ncomp,0:ncomp)      !< caller-owned scratch
        real(dp), intent(out)   :: Cout(ncomp*ncomp), Mout(ncomp)
    end subroutine polar_ring_gram

    real(dp) module function polar_ring_selfpower( Us, ldu, row0, nrow )
        integer, intent(in) :: ldu, row0, nrow
        real,    intent(in) :: Us(ldu,0:*)
    end function polar_ring_selfpower

    real(dp) module function polar_self_energy( xws, wr, pg ) result( e )
        real,               intent(in) :: xws(:), wr(:)
        type(polar_grid_t), intent(in) :: pg
    end function polar_self_energy

    pure real(dp) module function sum_dp_safe( acc, n ) result( v )
        real(dp), intent(in) :: acc
        integer,  intent(in) :: n
    end function sum_dp_safe

    module subroutine align_halfplane_inplane( frlims, nyq_eff, src, ca, sa, dst )
        integer, intent(in)  :: frlims(3,2), nyq_eff
        complex, intent(in)  :: src(frlims(1,1):frlims(1,2), frlims(2,1):0)
        real,    intent(in)  :: ca, sa
        complex, intent(out) :: dst(frlims(1,1):frlims(1,2), frlims(2,1):0)
    end subroutine align_halfplane_inplane

    module subroutine polar_sample_particle_packed( fpl, pg, ca, sa, xws, wr, hfpw, hfcnt, tazim, xws1, xws2 )
        type(fplane_type),  intent(in)    :: fpl
        type(polar_grid_t), intent(in)    :: pg
        real,               intent(in)    :: ca, sa
        real,               intent(out)   :: xws(:)
        real,               intent(out)   :: wr(:)
        real(dp),           intent(inout) :: hfpw, hfcnt
        real,               intent(out)   :: tazim
        real, optional,     intent(out)   :: xws1(:), xws2(:)
    end subroutine polar_sample_particle_packed

    module subroutine project_fplane_mean_banded( rec, o, fpl_ref, fpl_out )
        type(reconstructor), intent(in)    :: rec
        class(ori),          intent(inout) :: o
        type(fplane_type),   intent(in)    :: fpl_ref
        type(fplane_type),   intent(inout) :: fpl_out
    end subroutine project_fplane_mean_banded

    module subroutine polar_hybrid_exact_accum( rec0, recs, ncomp, o, fpl, hex, kex, npos, &
            &Gd, bd, cd, e_mm, myv )
        type(reconstructor), intent(in)    :: rec0
        type(reconstructor), intent(in)    :: recs(ncomp)
        integer,             intent(in)    :: ncomp, npos
        class(ori),          intent(inout) :: o
        type(fplane_type),   intent(in)    :: fpl
        integer,             intent(in)    :: hex(npos), kex(npos)
        real(dp),            intent(inout) :: Gd(ncomp,ncomp), bd(ncomp), cd(ncomp)
        real(dp),            intent(inout) :: e_mm, myv
    end subroutine polar_hybrid_exact_accum

    module subroutine subtract_mean_banded( fpl, mean_fpl, a, rec_nyq )
        type(fplane_type), intent(inout) :: fpl
        type(fplane_type), intent(in)    :: mean_fpl
        real,              intent(in)    :: a
        integer,           intent(in)    :: rec_nyq
    end subroutine subtract_mean_banded

    ! ===== implemented in simple_flex_pca_em_embed =====

    module subroutine embed_latents_with_contrast( params, build, mean_rec, basis_recs, ncomp, eigvals, sig2_eff, &
        &pinds, nptcls, z, contrast, precision, resid_energy, resid_mean_energy, rho_out, stats_only, &
        &from_parts, rounds, zhalf_out)
        class(flex_pca_rounds), intent(inout) :: rounds
        class(parameters),   intent(inout) :: params
        type(builder),       intent(inout) :: build
        type(reconstructor), intent(inout) :: mean_rec
        type(reconstructor), intent(inout) :: basis_recs(ncomp)
        integer,             intent(in)    :: ncomp, pinds(:), nptcls
        real(dp),            intent(in)    :: eigvals(ncomp)
        real(dp),            intent(in)    :: sig2_eff
        real(dp),            intent(out)   :: z(nptcls,ncomp), contrast(nptcls)
        real(dp),            intent(out)   :: precision(ncomp,ncomp,nptcls)
        real(dp),            intent(out)   :: resid_energy(nptcls), resid_mean_energy(nptcls)
        real(dp), optional,  intent(out)   :: rho_out(ncomp)
        logical,  optional,  intent(in)    :: stats_only
        logical,  optional,  intent(in)    :: from_parts
        real(dp), optional,  intent(out)   :: zhalf_out(:,:,:)
    end subroutine embed_latents_with_contrast

    ! ===== implemented in simple_flex_pca_em_iter =====

    module subroutine run_flex_pca_paired( params, build, pinds, nptcls, col_sep, neigs_req, &
        &m_basis_recs, m_eigvals, m_ncomp, m_sig2, m_matchcos , rounds)
        class(flex_pca_rounds), intent(inout) :: rounds
        class(parameters),   intent(inout) :: params
        type(builder),       intent(inout) :: build
        integer,             intent(in)    :: pinds(:), nptcls, col_sep, neigs_req
        !> merged final-stage product (SIMPLE_COV_PAIRED_MERGE>=1 only; stay unallocated/0
        !! otherwise) -- the caller embeds ALL particles against this basis
        type(reconstructor), allocatable, optional, intent(out) :: m_basis_recs(:)
        real(dp),            allocatable, optional, intent(out) :: m_eigvals(:)
        integer,             optional,    intent(out)           :: m_ncomp
        real(dp),            optional,    intent(out)           :: m_sig2
        !> per merged component: cross-half match |cos| (the CV rank-cut signal)
        real(dp), allocatable, optional,  intent(out)           :: m_matchcos(:)
    end subroutine run_flex_pca_paired

    module subroutine probe_subspace_iteration( params, build, mean_rec, basis_recs, eigvals, sig2_eff, &
        &pinds, nptcls, ncomp, niters, it_glob, niters_glob, fprefix, meta_fname, rounds)
        class(flex_pca_rounds), intent(inout) :: rounds
        class(parameters),   intent(inout) :: params
        type(builder),       intent(inout) :: build
        type(reconstructor), intent(inout) :: mean_rec
        type(reconstructor), allocatable, intent(inout) :: basis_recs(:)
        real(dp),            allocatable, intent(inout) :: eigvals(:)
        real(dp),            intent(in)    :: sig2_eff
        integer,             intent(in)    :: pinds(:), nptcls, niters
        integer,             intent(inout) :: ncomp
        integer, optional,   intent(in)    :: it_glob, niters_glob
        character(len=*), optional, intent(in) :: fprefix, meta_fname
    end subroutine probe_subspace_iteration

    ! ===== implemented in simple_flex_pca_em_solve =====

    pure module subroutine spd_logdet_dp( A, n, logdet, ok )
        integer,  intent(in)  :: n
        real(dp), intent(in)  :: A(n,n)
        real(dp), intent(out) :: logdet
        logical,  intent(out) :: ok
    end subroutine spd_logdet_dp

    module subroutine deflate_against_basis( imgs, n, basis, nb )
        integer,     intent(in)    :: n, nb
        type(image), intent(inout) :: imgs(n), basis(nb)
    end subroutine deflate_against_basis

    module subroutine cross_half_subspace_angles( eimgs, oimgs, n, svals )
        integer,     intent(in)    :: n
        type(image), intent(inout) :: eimgs(n), oimgs(n)
        real(dp), allocatable, intent(out) :: svals(:)
    end subroutine cross_half_subspace_angles

    pure module function quad_form( M, z, n ) result( val )
        integer,  intent(in) :: n
        real(dp), intent(in) :: M(n,n), z(n)
        real(dp) :: val
    end function quad_form

    module subroutine spd_solve_dp( A, b, n )
        integer,  intent(in)    :: n
        real(dp), intent(inout) :: A(n,n), b(n)
    end subroutine spd_solve_dp

    module subroutine probe_solve_ecm( n, G, b, c, myv, e_mm, prior_, sig2, nml, a, z_, Ainv_, ldA, lok, quad )
        integer,  intent(in)    :: n, nml
        real(dp), intent(in)    :: G(n,n), b(n), c(n), myv, e_mm, prior_(n), sig2
        real(dp), intent(inout) :: a
        real(dp), intent(out)   :: z_(n), Ainv_(n,n), ldA
        logical,  intent(out)   :: lok
        real(dp), intent(out)   :: quad
    end subroutine probe_solve_ecm

    module subroutine probe_solve_mix( n, kmix, G, b, c, myv, e_mm, nml, a, sig2, Ominv, Omxi, lpi, xiOx, &
            &zbar, Ainv_, dens_, ldA, lok, nll_add, sr_acc, sm_acc, smm_acc, sainv_acc )
        integer,  intent(in)    :: n, kmix, nml
        real(dp), intent(in)    :: G(n,n), b(n), c(n), myv, e_mm, sig2
        real(dp), intent(inout) :: a
        real(dp), intent(in)    :: Ominv(n,n), Omxi(n,kmix), lpi(kmix), xiOx(kmix)
        real(dp), intent(out)   :: zbar(n), Ainv_(n,n), dens_(n,n), ldA, nll_add
        logical,  intent(out)   :: lok
        real(dp), intent(inout) :: sr_acc(kmix), sm_acc(n,kmix), smm_acc(n,n,kmix), sainv_acc(n,n)
    end subroutine probe_solve_mix

    module subroutine mcfa_init( z, nptcls, ncomp, kmix, gam_sum, nval, xi, ppi, Om )
        integer,  intent(in)  :: nptcls, ncomp, kmix, nval
        real(dp), intent(in)  :: z(nptcls,ncomp), gam_sum(ncomp)
        real(dp), intent(out) :: xi(ncomp,kmix), ppi(kmix), Om(ncomp,ncomp)
    end subroutine mcfa_init

    module subroutine mcfa_condition( ncomp, diag_only, Om, Ominv, ldOm )
        integer,  intent(in)    :: ncomp
        logical,  intent(in)    :: diag_only
        real(dp), intent(inout) :: Om(ncomp,ncomp)
        real(dp), intent(out)   :: Ominv(ncomp,ncomp), ldOm
    end subroutine mcfa_condition

    module subroutine mcfa_mstep( ncomp, kmix, nval, sr, sm, smm, sai, &
        &pin_origin, ppi, xi, Om, Ominv, ldOm )
        integer,  intent(in)    :: ncomp, kmix, nval
        real(dp), intent(in)    :: sr(kmix), sm(ncomp,kmix)
        real(dp), intent(in)    :: smm(ncomp,ncomp,kmix), sai(ncomp,ncomp)
        logical,  intent(in)    :: pin_origin
        real(dp), intent(inout) :: ppi(kmix), xi(ncomp,kmix), Om(ncomp,ncomp)
        real(dp), intent(out)   :: Ominv(ncomp,ncomp), ldOm
    end subroutine mcfa_mstep

    module subroutine spd_inv_dp( A, Ainv, n )
        integer,  intent(in)    :: n
        real(dp), intent(inout) :: A(n,n)
        real(dp), intent(out)   :: Ainv(n,n)
    end subroutine spd_inv_dp

end interface

end module simple_flex_pca_em
