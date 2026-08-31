!@descr: flex_pca EM: the subspace EM iteration (E-step accumulation, M-step update)
submodule (simple_flex_pca_em) simple_flex_pca_em_iter
use simple_flex_pca_distr,  only: flex_pca_is_master, flex_pca_is_worker, flex_pca_nparts,&
    &flex_pca_run_stage, PCA_STAGE_PROBE
use simple_flex_pca_parts,  only: flex_pca_part_fname, write_probe_part, reduce_probe_parts
use simple_matcher_3Drec,   only: init_rec, cleanup_rec_buffers
use simple_matcher_ptcl_io, only: discrete_read_imgbatch, prepimgbatch
use simple_flex_reconstructor_latent_ops, only: project_fplane_mean, project_fplanes_mean_basis,&
    &insert_planes_oversamp_coupled_batch_scaled
use simple_flex_projected_latent_model, only: prep_imgs4projected_model, solve_coupled_basis_exp,&
    &projected_model_kfromto
use simple_flex_gpu,        only: flex_gpu_available, flex_gpu_coupled_begin_f,&
    &flex_gpu_coupled_batch_raw_f, flex_gpu_coupled_end_f, flex_gpu_coupled_bank_f,&
    &flex_gpu_coupled_batch_banked_f, flex_gpu_coupled_bank_free_f, flex_gpu_estep_vols_f,&
    &flex_gpu_estep_batch_f, flex_gpu_estep_resid_f, flex_gpu_estep_free_f,&
    &flex_gpu_coupled_batch_banked_res_f, flex_gpu_prep_begin_f, flex_gpu_prep_batch_f,&
    &flex_gpu_prep_free_f, flex_gpu_estep_batch_res_f, flex_gpu_prep_check_f,&
    &flex_gpu_poles_begin_f, flex_gpu_poles_bank_f, flex_gpu_poles_batch_f, flex_gpu_poles_free_f
use simple_flex_pca_polar,  only: polar_grid_build, polar_grid_kill, polar_project_recs,&
    &polar_relative_inplane, polar_assign_directions, polar_sample_particle_fused
implicit none
#include "simple_local_flags.inc"

contains

    !> Probe-based subspace iteration: alternate a Wiener E-step (per-particle latents in the current
    !! basis) with a weighted-backprojection M-step (Y_q += sum_i z_iq * backproject(r_i)), then
    !! orthonormalize the refined probe volumes into the next basis.
    module subroutine probe_subspace_iteration( params, build, mean_rec, basis_recs, eigvals, sig2_eff, &
        &pinds, nptcls, ncomp, niters )
        use simple_ptcl_cache, only: ptcl_cache_in_use, ptcl_cache_read_batch
        class(parameters),   intent(inout) :: params
        type(builder),       intent(inout) :: build
        type(reconstructor), intent(inout) :: mean_rec
        type(reconstructor), allocatable, intent(inout) :: basis_recs(:)
        real(dp),            allocatable, intent(inout) :: eigvals(:)
        real(dp),            intent(in)    :: sig2_eff
        integer,             intent(in)    :: pinds(:), nptcls, niters
        integer,             intent(inout) :: ncomp
        type(fplane_type), allocatable :: fpls(:), basis_fpls(:,:), mean_fpl(:)
        type(ori),           allocatable :: orientations(:)
        type(reconstructor), allocatable :: Yeven(:), Yodd(:), utilde(:)
        type(image),         allocatable :: realvols(:), utilde_real(:)
        type(image) :: img_o, mstep_gridcorr
        real,                allocatable :: rho_e(:,:,:,:), rho_o(:,:,:,:), filt(:), corrs(:)
        ! GPU M-step path (SIMPLE_COV_GPU=1): one coupled device accumulation over the combined
        ! [Yeven, Yodd] x [rho_e, rho_o] layout, halfsets routed through the scale slots
        logical  :: l_gpu_probe
        integer  :: vgpu, gq, gr
        real(dp), allocatable :: gdsc(:,:), grsc(:,:)
        ! banked coupled adjoint (SIMPLE_COV_GPU_BANKADJ, default on with the GPU M-step):
        ! probe records sorted by bank direction, in-plane aligned on device, GEMM-combined per
        ! direction, one splat per (direction, slot) instead of per particle
        logical  :: l_bank_adj
        integer  :: vbank, ndir_pb, jdir
        type(oris) :: dirs_pb
        type(ori)  :: o_pb
        real     :: ca1, sa1
        real,    allocatable :: rmat_pb(:,:,:), nrm_pb(:,:), rmat_pp(:,:,:), nrm_pp(:,:)
        real,    allocatable :: cap(:), sap(:)
        integer, allocatable :: dirof_pb(:), cnt_pb(:), perm_pb(:), pptmp(:)
        real(dp),            allocatable :: zbatch(:,:), dens(:,:,:), z(:,:), prior(:)
        real(dp),            allocatable :: Gth(:,:,:), Ath(:,:,:), bth(:,:), cth(:,:), zth(:,:)
        !> per-thread posterior covariance A^-1 and the untouched copy of A it is taken from
        real(dp),            allocatable :: Ainvth(:,:,:), Acpth(:,:,:)
        !> marginal-likelihood scratch: h_i before the solve, and the per-thread accumulator of
        !! log det A_i - h_i' A_i^-1 h_i. The remaining terms of -2 log p(y_i) are the per-particle
        !! data energy (fixed across iterations, since a_i and mu are) and N*log det Gamma, added at
        !! the reduction. See spd_logdet_dp for why resid_energy could never serve this purpose.
        real(dp),            allocatable :: hth(:,:), nll_thr(:), gam_dbg(:,:)
        real(dp) :: ldA, nll_tot, nll_prev
        !> principal-angle convergence threshold. Overridable because the default 0.97 was tuned for a
        !! probe that STARTED from the covariance eigenbasis and only had to clean it: from a data-free
        !! start the subspace stops rotating several iterations before Gamma stops moving (measured on
        !! 10076: cos crosses 0.97 at iteration 5 while max var is still falling ~20 % per iteration),
        !! so a pure-EM fit needs a tighter bar. SIMPLE_COV_PROBE_CONV is in PER MILLE (995 = 0.995).
        real(dp) :: conv_thresh
        integer  :: vconv
        !> even/odd half-basis agreement: the dataset-agnostic convergence signal
        type(image), allocatable :: eimgs(:), oimgs(:)
        real(dp),    allocatable :: sv_eo(:)
        real(dp) :: eo_dim, eo_best
        integer  :: eo_stall, vpat, eo_patience
        !> BAND ANNEALING, SIMPLE_COV_EM_ANNEAL = number of iterations to reach the full band.
        !! Expressed as a k_hi schedule and converted to lp, NOT set through lp directly: on this
        !! geometry covariance_kfromto only narrows the band when lp > 2*smpd_crop, so a schedule
        !! written in lp is silently a no-op at common settings. Rationale is that a from-scratch
        !! basis has nothing to constrain the high shells in the first iterations, so letting them
        !! in early spends rank on noise. NOTE: this schedule is the proposal's own invention -- it
        !! is not from 3DVA, which uses a single fixed filter -- so it is OFF by default and is a
        !! thing to A/B, not a thing to assume.
        integer  :: vanneal, khi_full, khi_it, kfr_ann(2)
        integer, parameter :: MIX_ZSUB_MAX = 2000   ! per-part latent subsample shipped for mcfa_init
        real(dp), allocatable :: dm_sr(:), dm_sm(:,:), dm_smm(:,:,:), dm_sai(:,:), dm_z(:,:)
        integer :: dm_nz = 0
        real     :: lp_it, dstep_ann
        !> POLAR EM E-STEP, SIMPLE_COV_EM_POLAR=1.
        !! probe_subspace_iteration has never had a polar former: its three E-step branches are all
        !! Cartesian (GPU-fused exact, banked-forward, CPU exact), and cov_polar_enabled() is only
        !! consulted by reduced_covariance_solve and embed_latents_with_contrast. This routes the
        !! per-particle sufficient statistics through embed_accumulate_polar -- the SAME polar
        !! former the embedding uses, which fills exactly Gcache/bcache/ccache/contrast -- and then
        !! keeps the Cartesian residual pass for the M-step, because there is no polar->volume
        !! adjoint anywhere in the tree and there is no prospect of one.
        !! The residual pass projects the MEAN ONLY, not the r+1 volumes, so at large r this is
        !! cheaper than the Cartesian E-step it replaces as well as being the right representation.
        logical  :: l_em_polar
        integer  :: vempol, nvalid_pol
        real(dp), allocatable :: pGc(:,:,:), pbc(:,:), pcc(:,:), pzh(:,:,:)
        real(dp), allocatable :: pcon(:), pre(:), prme(:)
        !> POLAR SHARED-DIRECTION E-STEP, SIMPLE_COV_POLAR_ESTEP=1 (stage 1).
        !! IN the batch loop, unlike SIMPLE_COV_EM_POLAR (which routes the whole statistics pass
        !! through embed_accumulate_polar up front and is windows/mixture-incompatible): once per
        !! EM iteration the mean + ncomp basis volumes are projected onto the polar rings of a
        !! shared direction table (the bank); per particle, the bank rings at (direction, psi)
        !! replace the ncomp+1 per-particle Cartesian projections in forming G/b/c/e_mm/myv, and
        !! everything downstream -- the ECM/MCFA solve, mixture accumulators, ONLINE windows,
        !! e/o halfsets, the Cartesian M-step insertion -- is untouched. CTF amplitude, shift
        !! phase and per-shell whitening ride in exactly as in the Cartesian former, because the
        !! SAME prepped cmplx/transfer planes are polar-sampled (the validated embed formulation).
        !! The mean-shaped deflation needs no special handling here: dfl_basis deflates the
        !! refined basis VOLUMES at the end of each M-step, and the bank is rebuilt from those
        !! same basis_recs at the next iteration, so the deflation enters the bank exactly as it
        !! enters the per-particle Cartesian projections.
        logical  :: l_pol_es, l_pol_grid, l_pol_bank_it
        integer  :: vpoles, n_pol_check, ndir_es, nsamp_es, nsamp2_es, nk_es, id_es, idp_es, ir
        integer  :: ph0_es, pk0_es, hlo_es, hhi_es, klo_es, nyqr_es, nyqb_es
        !> HYBRID exact/ring quadrature (accuracy): the shells 0..rhyb_es are accumulated as
        !! EXACT Cartesian lattice statistics per particle (data, CTF and model read exactly as
        !! the Cartesian former reads them -- including the DC sample the rings never had) and
        !! the shared-direction rings cover only the annulus above. Measured on the calibrated
        !! synthetic dissector: the pure ring quadrature carries a MULTIPLICATIVE low bias of
        !! the posterior latents (z-variance ratio ~0.85 vs Cartesian, the eigen-spectrum
        !! collapse mechanism of the first science A/B), which angular/radial oversampling and
        !! the DC sample alone do NOT fix; the hybrid at rhyb ~ 0.55*band restores ratio ~1.00
        !! with G err 0.4% and b err ~3%. SIMPLE_COV_POLAR_RHYB overrides (0=off: pure rings,
        !! tonight's baseline); SIMPLE_COV_POLAR_OSAMP multiplies ring angular sampling.
        logical  :: l_pol_hyb, l_rhyb_off
        integer  :: vpol_rhyb, vpol_osamp, rhyb_es, npos_es, jx_es, kx_es
        integer,  allocatable :: hex_es(:), kex_es(:)
        type(polar_grid_t)    :: pg_es
        type(oris)            :: dirs_es
        type(ori)             :: o_es
        real,     allocatable :: rmatp_es(:,:,:), nrmp_es(:,:), rmatb_es(:,:,:), nrmb_es(:,:)
        real,     allocatable :: cae(:), sae(:)
        integer,  allocatable :: dir_es(:)
        logical,  allocatable :: dused_es(:)
        real,     allocatable :: UsallE(:,:,:)                       !< (nsamp2, 0:ncomp, ndir) the bank
        real(dp), allocatable :: CfE(:,:,:), Cm0E(:,:,:), c00E(:,:)  !< ring Gram tables per direction
        complex,  allocatable :: UbankE(:,:,:)                       !< per-thread projection scratch
        real,     allocatable :: CspE(:,:,:)                         !< per-thread ring-Gram scratch
        real,     allocatable :: xws_es(:,:), wr_es(:,:), Reb_es(:,:)
        real(dp), allocatable :: wrd_es(:,:)
        real(dp) :: pw_es, cnt_es
        real     :: taz_es, sec_bank
        integer(timer_int_kind) :: t_bank
        !> DEVICE polar E-step (stage 2), SIMPLE_COV_POLAR_GPU=1 (requires POLAR_ESTEP=1): the
        !! host-built bank + ring Gram tables upload once per iteration; per batch the RAW
        !! (y, T) plane sections go up and per-particle G/b/c/e_mm/myv come back in the fused
        !! E-step's stat layout -- the ring part via a device sampler + bank/table contractions,
        !! the hybrid low-k exact part via the unmodified estep_stats_kernel with the disc
        !! capped at rhyb*(rhyb+1) against the resident mean+basis volumes. The solve (plain or
        !! mixture), Gamma/nll accumulation, ONLINE windows, e/o bookkeeping, the Cartesian
        !! mean projection/residual and the M-step insertion are the CPU polar path's,
        !! untouched. All device state is per-stage (ring geometry), per-iteration (bank,
        !! volumes) or per-batch (planes), so the path is window-safe under SIMPLE_COV_ONLINE
        !! by construction (the ONLINE x fused-device breakage was in the prep/banked caches,
        !! which this path does not use).
        logical  :: l_pol_gpu, l_pol_gpu_any
        integer  :: vpolg, nst_pg
        real(dp), allocatable :: stats_pg(:,:)
        real,     allocatable :: ca_pg(:), sa_pg(:)
        integer,  allocatable :: dir_pg(:)
        logical,  allocatable :: vld_pg(:)
        !> gate-1 GPU parity arrays: the CPU polar former re-run on the check particles
        real(dp), allocatable :: zchk_g(:,:), gerrg_chk(:), berrg_chk(:)
        real(dp), allocatable :: Gpc_chk(:,:,:), bpc_chk(:,:), cpc_chk(:,:), zpc_chk(:,:)
        real(dp) :: emmg_chk, myvg_chk, ag_chk
        !> gate-1 parity diagnostic, SIMPLE_COV_POLAR_CHECK=N: the first N particles of EM
        !! iteration 1 run BOTH formers; per-component z correlation and the median relative
        !! G/b errors are printed, then the run continues on the polar path.
        logical  :: l_chk_act, l_chk_p, lok_chk
        integer  :: nchk_max, ichk
        real(dp), allocatable :: zchk_c(:,:), zchk_p(:,:), gerr_chk(:), berr_chk(:)
        logical,  allocatable :: lchk(:)
        real(dp), allocatable :: Gc_chk(:,:,:), bc_chk(:,:), cc_chk(:,:), zc_chk(:,:), Ai_chk(:,:,:)
        real(dp) :: emmc_chk, myvc_chk, ac_chk, lda_chk, qml_chk
        real(dp) :: gnum_chk, gden_chk, bnum_chk, bden_chk
        !> exact-direction polar arm of the check: the SAME quadrature at the particle's own
        !! direction (no snap), which attributes any disagreement between direction quantization
        !! (banked vs exact) and the polar formulation itself (exact vs Cartesian). tazim -- the
        !! mean within-ring relative spread of |T|^2 -- bounds the radial-factorisation error in G.
        real(dp), allocatable :: zchk_x(:,:), gerrx_chk(:), berrx_chk(:)
        real,     allocatable :: UsX_chk(:,:,:), xwsX_chk(:,:), wrX_chk(:,:), RebX_chk(:,:)
        real(dp), allocatable :: CfX_chk(:,:,:), Cm0X_chk(:,:,:), c00X_chk(:,:), wrdX_chk(:,:)
        real(dp), allocatable :: GX_chk(:,:,:), bX_chk(:,:), cX_chk(:,:), zx_chk(:,:)
        real(dp), allocatable :: tazsum_thr(:)
        integer,  allocatable :: tazcnt_thr(:)
        real     :: rmx_chk(3,3), tazx_chk
        real(dp) :: pwx_chk, cntx_chk, emmx_chk, myvx_chk, ax_chk
        !> DC-less Cartesian arm of the check: cov_herm_inner includes the (0,0) sample, the polar
        !! quadrature starts at ring 1 and has none. If the exact-direction polar arm agrees with
        !! THIS reference, the polar/Cartesian gap is the DC term, not the ring quadrature.
        real(dp), allocatable :: zchk_d(:,:), gerrd_chk(:), berrd_chk(:)
        real(dp), allocatable :: Gd_chk(:,:,:), bd_chk(:,:), cd_chk(:,:), zd_chk(:,:)
        complex,  allocatable :: dcb_chk(:,:)
        complex  :: dcy_chk, dcm_chk
        real(dp) :: emmd_chk, myvd_chk, ad_chk
        logical  :: l_probe_distr_pre
        !> mean-shaped (contrast) deflation of the refined basis, SIMPLE_COV_EM_DEFLATE
        logical       :: l_deflate_mean
        integer       :: vdfl, ndfl, ndfl_sh, idfl, jdfl, nkeep_dfl, kfr_dfl(2)
        logical       :: l_dfl_bg
        logical       :: l_dfl_pose
        integer       :: ipdfl, ixp, iyp, izp
        real(dp)      :: gpx, gpy, gpz, cp_dfl
        type(image)   :: mvol_dfl
        type(image), allocatable :: dfl_basis(:)
        real, pointer :: rm_dfl(:,:,:), rv_dfl(:,:,:)
        real          :: res_lo, res_hi
        real(dp)      :: mm_dfl, mv_dfl, rem_dfl, tot_dfl, mnorm_dfl
        logical  :: lok
        !> per-thread accumulator for the EM Gamma update, and its reduction
        real(dp),            allocatable :: gam_thr(:,:), gam_acc(:)
        integer,             allocatable :: nval_thr(:)
        logical,             allocatable :: valid(:), valid_e(:), valid_o(:)
        integer,             allocatable :: eo(:)
        real, pointer :: rmatp(:,:,:)
        real     :: fc
        real(dp) :: sig2, a, aa, e_mm, myv, mu_q, sd_q
        integer  :: it, q, r, i, ithr, nthr, batchlims(2), batchsz, ibatch, row, d_new, es(3), filtsz, sh
        integer  :: npairs, nval, npp, nparts_sub
        logical  :: l_probe_distr
        real(dp),            allocatable :: gam_sum(:)
        complex,             allocatable :: cme(:,:,:,:), cmo(:,:,:,:)
        real,                allocatable :: rhe(:,:,:,:), rhoo(:,:,:,:)
        integer,             allocatable :: ppinds(:)
        type(image),         allocatable :: prev_real(:)      ! previous iteration's orthonormal basis
        real(dp),            allocatable :: Mconv(:,:), sconv(:)
        real(dp) :: cos_mean
        logical  :: l_converged
        integer  :: n_probe_cm, nml_plain
        real,     allocatable :: fscq_dg(:,:)
        real(dp) :: fmean_dg(512), fbest_dg
        real     :: res_dg
        integer  :: khi_dg, nsig_dg, ntop_dg, sel_dg(4), tq_dg, bq_dg
        ! ---- ONLINE (minibatch) EM: SIMPLE_COV_ONLINE=1 + SIMPLE_COV_BATCH=B ----
        ! Each iteration E-steps ONE rotating contiguous window of the (full) probe set and the
        ! M-step consumes RUNNING-AVERAGED sufficient statistics S_t = (1-g_t) R(S_{t-1}) + g_t s_t
        ! (Cappe-Moulines), with R the Procrustes rotation between successive basis frames. The
        ! Wiener merge therefore sees an effective sample that GROWS with iterations while each
        ! iteration costs one batch -- "the whole dataset informs the basis without increasing
        ! per-step compute". Naive per-iteration resampling without averaging is measured negative;
        ! full-N at fixed regularisation is measured negative (subsampling is regularisation) --
        ! this is the formulation that addresses both. v1: plain CPU body, nparts=1, contiguous
        ! project-order windows (sequential IO); rho/gamma rotate with R^2 weights (exact only for
        ! the map stats; second-order error near convergence, diluted by g early).
        logical  :: l_online, l_focus_annulus
        integer  :: vonl_pr, nb_online, ob_lo, ob_hi, b_online, b_win, it_win, vfoc_r1, vfoc_ed
        real     :: foc_r1, foc_edge
        real(dp) :: gam_onl
        complex, allocatable :: oh_cme(:,:,:,:), oh_cmo(:,:,:,:)
        real,    allocatable :: oh_rhe(:,:,:,:), oh_rho(:,:,:,:)
        real(dp), allocatable :: oh_gam(:), ut_hist(:,:), r_onl(:,:)
        logical  :: l_oh_live
        logical  :: l_probe_mls
        real(dp) :: qml, nll_mix_add
        ! ---- MCFA state: tied-covariance mixture prior over the latents ----
        integer  :: kmix_pr, n_mix_warm, kk2
        logical  :: l_mix_req, l_mix_active, l_mix_used
        real(dp) :: ldOm_mix, ldOm_used, lwm, wsm
        real(dp), allocatable :: mix_xi(:,:), mix_Om(:,:), mix_Ominv(:,:), mix_pi(:)
        real(dp), allocatable :: mix_Omxi(:,:), mix_xiOx(:), mix_lpi(:)
        real(dp), allocatable :: rhs0th(:,:), mkth(:,:,:), lwth(:,:), rkth(:,:)
        real(dp), allocatable :: mxa_sr(:,:), mxa_sm(:,:,:), mxa_smm(:,:,:,:), mxa_sainv(:,:,:)
        ! unified-mixture state: frame rotation (Gap B), full-N running average of the
        ! reduced mixture statistics (Gap A), starved-component reseeding, final full pass
        logical  :: l_onl_final
        type(string) :: fname
        integer(timer_int_kind) :: t_it, t_sec
        real    :: sec_read, sec_prep, sec_estep, sec_ins
        real(dp) :: twp0, twp1, twp2
        real(dp), allocatable :: sec_proj_thr(:), sec_gram_thr(:)
        ! banked FORWARD: shared raw projections per direction segment; data aligned into the
        ! bank frame once (which also hands the M-step a pre-aligned residual).
        ! MEASURED REFUTED 2026-08-14 (gpu15): science ARI 0.836 vs the 0.95 band -- the E-step
        ! latents cannot take direction-quantized projections (the M-step adjoint can) -- and
        ! the full-array CTF multiplies are memory-traffic-bound (gram bucket 30->200 thread-s).
        ! Kept opt-in behind SIMPLE_COV_GPU_BANKFWD for reference; default is the exact
        ! per-particle forward with the banked adjoint.
        type(fplane_type), allocatable :: mean_fplC(:), basis_fplsC(:,:)
        type(ori),         allocatable :: o_pb_thr(:)
        complex,           allocatable :: tmpc(:,:,:)
        real,              allocatable :: cap1(:), sap0(:)
        integer,           allocatable :: seg_dir_b(:), seg_beg_b(:), seg_cnt_b(:)
        integer  :: iseg, ii, nseg_b, nyqb
        logical  :: l_bank_fwd
        ! fused device E-step (exact per-particle directions, volumes resident)
        logical  :: l_gpu_estep, l_deflate
        integer  :: ves, nst_es, off_es
        real(dp), allocatable :: stats_es(:,:)
        real,     allocatable :: avec_es(:)
        logical,  allocatable :: vld_es(:)
        ! device prep (raw images up; packed planes born on device)
        logical  :: l_gpu_prep, l_pcache
        integer  :: vprep, frlims_pb(3,2), nyqpd_pb, kf_pb(2)
        type(ftiter) :: fit_pb
        type(ctfparams), allocatable :: ctfp_arr(:)
        real,            allocatable :: shf_pb(:,:), sig2_pb(:,:)
        integer :: signyq_pb, nyqfull_pb
        nthr = nthr_glob
        sig2 = max(sig2_eff, DTINY)
        ! ---- optional STRIDE subsample, for the basis refinement only ----
        ! The probe refines ncomp band-limited, FSC-regularised volumes, and every iteration costs a
        ! full pass over the data -- far more particles than that many parameters need. A stride keeps
        ! both halfsets and every state proportionally represented; fromp/top would NOT, because
        ! particles are commonly ordered by state (on Ribosembly a contiguous window selects whole
        ! states). The embedding stage that follows still uses every particle: only the basis
        ! refinement is subsampled.
        ! The stride MUST be applied within each halfset, not across the particle list. `eo` alternates
        ! strictly by particle index (0,1,0,1,...), so a plain stride of 2 selects one halfset entirely
        ! and leaves the other empty -- and every probe M-step is regularised by an even/odd FSC, which
        ! is then computed against nothing. Measured: the Wiener filter kills the basis and the run dies
        ! at the "embedding collapsed" guard. Striding per halfset keeps both populated at any stride.
        ! ---- absolute cap, not a fixed ratio ----
        ! The probe refines ncomp band-limited, FSC-regularised volumes, and that parameter count
        ! does not grow with the dataset, so the particles needed to determine it do not either. A
        ! constant stride would leave the probe scaling linearly and dominating the run; capping the
        ! count makes it O(1) in dataset size.
        !
        ! COV_PROBE_MAX_PTCLS is the total across all processes, so each takes its share -- a worker
        ! sees only its own partition and would otherwise take the whole budget nparts times over.
        ! SIMPLE_COV_PROBE_STRIDE still wins if set.
        ! Only a WORKER divides the total by nparts: it holds one fromp/top partition. The master
        ! holds every particle, so passing nparts there divides twice and inflates the stride by
        ! exactly nparts. See cov_stage_subsample, which the initialiser shares.
        nparts_sub = 1
        if( flex_pca_is_worker() ) nparts_sub = params%nparts
        call cov_stage_subsample(build, pinds, nptcls, nparts_sub, COV_PROBE_MAX_PTCLS, &
            &'SIMPLE_COV_PROBE_STRIDE', 'SIMPLE_COV_PROBE_MAX', 'PROBE', ppinds, npp)
        allocate(z(npp,ncomp))
        ! ---- banked-adjoint setup: pose geometry is fixed across EM iterations, so the bank,
        ! the direction assignment and the direction sort are built ONCE per stage. Sorting
        ! The polar and mixture E-step selectors are read BEFORE the GPU E-step gating below,
        ! which has to stand those paths down. The EM_POLAR read used to sit a hundred lines
        ! AFTER the gate that tests it, so the gate read an UNINITIALISED logical -- usually
        ! false on this stack, which is why it never bit in practice, but it was undefined
        ! behaviour standing in for "off".
        vempol     = 0
        call cov_env_int('SIMPLE_COV_EM_POLAR', vempol)
        l_em_polar = vempol > 0
        ! ---- POLAR E-STEP (SIMPLE_COV_POLAR_ESTEP=1): in-batch-loop shared-direction bank ----
        vpoles   = 0
        call cov_env_int('SIMPLE_COV_POLAR_ESTEP', vpoles)
        l_pol_es = vpoles > 0
        if( l_pol_es .and. l_em_polar ) THROW_HARD('SIMPLE_COV_POLAR_ESTEP + SIMPLE_COV_EM_POLAR &
            &are mutually exclusive: both replace the E-step former')
        n_pol_check = 0
        call cov_env_int('SIMPLE_COV_POLAR_CHECK', n_pol_check)
        if( n_pol_check > 0 .and. .not. l_pol_es )then
            write(logfhandle,'(A)') '>>> FLEX_PCA POLAR CHECK ignored: needs SIMPLE_COV_POLAR_ESTEP=1'
            n_pol_check = 0
        endif
        ! hybrid/oversampling accuracy knobs (see the declaration block for the measurements)
        vpol_rhyb  = 0
        call cov_env_int('SIMPLE_COV_POLAR_RHYB', vpol_rhyb)
        l_rhyb_off = cov_env_int_off('SIMPLE_COV_POLAR_RHYB')
        vpol_osamp = 1
        call cov_env_int('SIMPLE_COV_POLAR_OSAMP', vpol_osamp)
        vpol_osamp = max(1, min(8, vpol_osamp))
        l_pol_hyb  = .false.
        rhyb_es    = 0
        npos_es    = 0
        l_pol_grid    = .false.
        l_pol_bank_it = .false.
        l_chk_act     = .false.
        nchk_max      = 0
        sec_bank      = 0.
        if( l_pol_es )then
            write(logfhandle,'(A)') '>>> FLEX_PCA POLAR E-STEP ON (stage 1): shared-direction bank &
                &supplies G/b/c; data-plane prep and M-step insertion stay Cartesian'
            if( n_pol_check > 0 ) write(logfhandle,'(A,I0,A)') '>>> FLEX_PCA POLAR CHECK: dual &
                &former on the first ',n_pol_check,' particles of iteration 1'
            call flush(logfhandle)
        endif
        ! ---- DEVICE polar E-step (stage 2), SIMPLE_COV_POLAR_GPU=1: opt-in, composes with the
        ! production env set; falls back to the CPU polar path (same answer) when the card or
        ! the CUDA build is absent, or when the basis exceeds the resident-volume cap below.
        l_pol_gpu = .false.
        if( l_pol_es )then
            vpolg = 0
            call cov_env_int('SIMPLE_COV_POLAR_GPU', vpolg)
            if( vpolg > 0 )then
                if( flex_gpu_available() )then
                    l_pol_gpu = .true.
                    write(logfhandle,'(A)') '>>> FLEX_PCA POLAR E-STEP DEVICE ON (stage 2): bank + &
                        &ring tables resident per iteration, per-batch ring/hybrid statistics on &
                        &device; solve, windows, mixture and M-step unchanged'
                else
                    write(logfhandle,'(A)') '>>> FLEX_PCA POLAR GPU requested but no device/CUDA &
                        &build available: using the CPU polar E-step'
                endif
                call flush(logfhandle)
            endif
        endif
        l_pol_gpu_any = l_pol_gpu
        ! ---- MCFA (SIMPLE_COV_EM_MIX=K): tied-covariance K-component mixture latent prior ----
        ! Mixtures of common factor analysers (Baek, McLachlan & Flack, IEEE TPAMI 32:1298,
        ! 2010): ONE shared basis, K Gaussian components in the latent space, fitted in the same
        ! EM as the basis itself. The structural point: a direction earns rank in U only if it
        ! improves the fit under a MULTI-MODAL prior, so purely continuous nuisance variation --
        ! fit equally well by any single component -- gains nothing. This is the one thing a
        ! moment estimator cannot reproduce, because responsibilities are posterior objects.
        ! K=1 pins the component at the origin with a diagonal covariance, which reduces EXACTLY
        ! to the plain PPCA EM -- that is the mandated regression test, not a feature.
        kmix_pr = 0
        call cov_env_int('SIMPLE_COV_EM_MIX', kmix_pr)
        kmix_pr   = max(0, min(64, kmix_pr))
        l_mix_req = kmix_pr >= 1
        ! v2 (2026-08-21): the four MCFA accumulators are additive sufficient statistics and now
        ! ride in the probe part files (PROBE_PART_VERSION 3), together with a bounded latent
        ! subsample so the master can seed mcfa_init. nparts>1 is therefore supported.
        if( l_mix_req .and. l_em_polar ) THROW_HARD('SIMPLE_COV_EM_MIX + SIMPLE_COV_EM_POLAR &
            &are mutually exclusive')
        vonl_pr = 0
        call cov_env_int('SIMPLE_COV_ONLINE', vonl_pr)
        l_online = vonl_pr > 0
        vfoc_r1 = 0
        call cov_env_int('SIMPLE_COV_FOCUS_R1', vfoc_r1)
        l_focus_annulus = vfoc_r1 > 0
        foc_r1 = real(vfoc_r1)
        vfoc_ed = 10
        call cov_env_int('SIMPLE_COV_FOCUS_EDGE', vfoc_ed)
        foc_edge = real(max(1, vfoc_ed))
        if( l_focus_annulus )then
            write(logfhandle,'(A,F6.1,A,F5.1,A)') '>>> FLEX_PCA FOCUSED HETEROGENEITY: basis &
                &restricted to the shell beyond r=', foc_r1, ' A (soft edge ', foc_edge, ' A)'
            call flush(logfhandle)
        endif
        l_oh_live = .false.
        b_online = 40000
        call cov_env_int('SIMPLE_COV_BATCH', b_online)
        ! Under windows-only ONLINE the delivered basis is the LAST window's M-step --
        ! and which window that is falls out of mod(iters, nwindows). ONLINE_FINAL runs
        ! the final iteration over the FULL probe set, so the endpoint estimator is
        ! built from every particle with no history blending involved (the
        ! full-N-under-mixture-prior configuration, measured good at the ceiling).
        l_onl_final = l_online .and. cov_env_int_on('SIMPLE_COV_ONLINE_FINAL')
        ! assigned BEFORE its first read: this guard used to read l_probe_distr_pre ~200 lines
        ! before the assignment below -- undefined behaviour standing in for "off", the same trap
        ! the EM_POLAR selector note above documents
        l_probe_distr_pre = flex_pca_is_master() .and. flex_pca_nparts() > 1
        if( l_online )then
            ! v2 (2026-08-21): ONLINE distributes. The recurrence itself -- the running average
            ! S_t = (1-g_t) R(S_{t-1}) + g_t s_t and its Procrustes frame rotation -- is MASTER-side
            ! and operates on the REDUCED statistics, so workers stay stateless: each contributes
            ! s_t for its slice of the current window. The only thing that has to be made
            ! distribution-aware is the window walk itself (below): each rank walks its own shard
            ! with a per-rank window of b_online/nparts, so the UNION across ranks is one global
            ! window of the intended size b_online, evenly spread and load-balanced. At nparts=1
            ! this reduces identically to the shared-memory walk.
        endif
        n_mix_warm = 3
        call cov_env_int('SIMPLE_COV_EM_MIX_WARM', n_mix_warm)
        n_mix_warm   = max(1, min(20, n_mix_warm))
        l_mix_active = .false.
        l_mix_used   = .false.
        ldOm_mix     = 0.d0
        ldOm_used    = 0.d0
        if( l_mix_req ) write(logfhandle,'(A,I0,A,I0,A)') '>>> FLEX_PCA MIX requested: K=', &
            &kmix_pr,'  (plain warm-up: ',n_mix_warm,' iterations)'
        ! ppinds is safe: z is only ever consumed through order-insensitive sums, and z/ppinds
        ! stay aligned. Gated to nsym==1 (the bank factorization R = R_dir * Rz(psi) has no
        ! symmetry replication) and to the in-process E-step (the distributed master never runs
        ! the batch loop).
        l_bank_adj = .false.
        l_bank_fwd = .false.
        l_gpu_estep = .false.
        ! initialized HERE, not inside the SIMPLE_COV_GPU block below: the batch loop reads it
        ! on the CPU path too (skips the CPU prep and touches unallocated ctfp_arr if garbage-true)
        l_gpu_prep = .false.
        vgpu = 0;  call cov_env_int('SIMPLE_COV_GPU', vgpu)
        vbank = 1; call cov_env_int('SIMPLE_COV_GPU_BANKADJ', vbank)
        if( cov_env_int_off('SIMPLE_COV_GPU_BANKADJ') ) vbank = 0
        if( vgpu > 0 .and. vbank > 0 .and. flex_gpu_available() .and. &
            &build%pgrpsyms%get_nsym() == 1 .and. &
            &.not. (flex_pca_is_master() .and. flex_pca_nparts() > 1) )then
            l_bank_adj = .true.
            ndir_pb = cov_polar_ndir(npp)
            call dirs_pb%new(ndir_pb, is_ptcl=.false.)
            call build%pgrpsyms%build_refspiral(dirs_pb)
            allocate(rmat_pb(3,3,ndir_pb), nrm_pb(3,ndir_pb))
            do jdir = 1, ndir_pb
                rmat_pb(:,:,jdir) = dirs_pb%get_mat(jdir)
                nrm_pb(:,jdir)    = rmat_pb(3,:,jdir)
            end do
            ! dirs_pb stays alive: the banked forward projects the shared basis at bank
            ! orientations pulled from it per segment
            vbank = 0
            call cov_env_int('SIMPLE_COV_GPU_BANKFWD', vbank)
            l_bank_fwd = vbank > 0
            ves = 1
            call cov_env_int('SIMPLE_COV_GPU_ESTEP', ves)
            ! cov_env_int only assigns for values > 0, so it cannot express "off" and =0 was
            ! silently ignored -- which voided the first MCFA reduction test: the plain arm ran
            ! the FUSED body while the mixture arm ran the CPU body, and the pair compared
            ! former numerics instead of the mixture algebra. Same trap, same fix, as
            ! SIMPLE_COV_EM_DEFLATE=0.
            l_gpu_estep = ves > 0 .and. .not. cov_env_int_off('SIMPLE_COV_GPU_ESTEP')
            ! the polar former supplies G, b, c and the contrast itself, so the device E-step and
            ! the banked forward must both stand down or their statistics would be the ones used
            if( l_em_polar .or. l_pol_es )then
                l_gpu_estep = .false.
                l_bank_fwd  = .false.
            endif
            ! the banked forward still lacks a mixture solve; the fused device E-step COMPOSES
            ! with the mixture since 2026-08-19 -- the device fetches the same G/b/c sufficient
            ! statistics and the host branches to probe_solve_mix (one shared solver, both bodies)
            if( l_mix_req )then
                l_bank_fwd  = .false.
            endif
            ! CAPACITY GATE. The device holds mean+basis as ONE resident allocation and the CUDA
            ! entry point refuses ncomp+1 > COV_GPU_ESTEP_MAXVOLS. Without this check the worker
            ! dies inside flex_gpu_estep_vols_f (THROW_HARD) and, because a dead worker never writes
            ! JOB_FINISHED, the distributed MASTER then waits forever -- an 8-hour silent hang
            ! observed at neigs=40. Fall back to the CPU E-step instead: slower, same answer.
            if( l_gpu_estep .and. (ncomp + 1) > COV_GPU_ESTEP_MAXVOLS )then
                l_gpu_estep = .false.
                write(logfhandle,'(A,I0,A,I0,A)') '>>> FLEX_PCA PROBE fused device E-STEP OFF: needs ', &
                    &ncomp+1,' resident volumes, device entry point caps at ',COV_GPU_ESTEP_MAXVOLS, &
                    &' -- using the CPU E-step for this basis size'
                call flush(logfhandle)
            endif
            ! ONLINE x device is MEASURED BROKEN (2026-08-19, onlk8_dev: ceiling 0.2087 vs
            ! 0.2742 CPU at the identical config -- the representation itself collapses).
            ! Some device-side state IS keyed to a fixed particle set across iterations
            ! (candidate: the banked direction segments / prep caches; the streaming argument
            ! that justified removing this guard was wrong). CPU body under windows until the
            ! device path is made window-safe and re-validated.
            if( l_online .and. (l_gpu_estep .or. l_bank_adj) )then
                l_gpu_estep = .false.
                l_bank_adj  = .false.
                write(logfhandle,'(A)') '>>> FLEX_PCA ONLINE: plain CPU E-step body (device path &
                    &measured broken under rotating windows; see 2026-08-19 note)'
                call flush(logfhandle)
            endif
            if( l_gpu_estep )then
                allocate(avec_es(MAXIMGBATCHSZ), source=0.0)
                allocate(vld_es(MAXIMGBATCHSZ), source=.false.)
            endif
            ! mean deflation of the M-step residual (SIMPLE_COV_PROBE_DEFLATE=0 to skip):
            ! a geometry defect used to skip it silently on the fused path, and the
            ! no-deflation variant measured BETTER twice (ARI 0.970/0.981 incl. 16/16 GT
            ! states, +2 probe rounds) -- kept as an explicit arm until the A/B settles it.
            ! presence-and-zero test: cov_env_int silently ignores <= 0 (the CGSOLVE trap)
            l_deflate = .not. cov_env_int_off('SIMPLE_COV_PROBE_DEFLATE')
            if( l_gpu_estep .and. .not. l_deflate )then
                write(logfhandle,'(A)') '>>> FLEX_PCA PROBE M-STEP DEFLATION OFF (raw-data subspace iteration)'
                call flush(logfhandle)
            endif
            ! device prep: mask path only, no sigma2 whitening; planes born on device
            l_gpu_prep = .false.
            if( l_gpu_estep )then
                vprep = 1
                call cov_env_int('SIMPLE_COV_GPU_PREP', vprep)
                l_gpu_prep = vprep > 0
                if( l_gpu_prep .and. params%l_ml_reg )then
                    ! whitening path: needs the loaded sigma2 spectra (the CPU prep THROWs
                    ! without them, so this mirrors its requirement)
                    l_gpu_prep = allocated(build%esig%sigma2_noise)
                endif
            endif
            if( l_gpu_prep )then
                kf_pb = projected_model_kfromto(params)
                call fit_pb%new([params%boxpd, params%boxpd, 1], params%smpd_crop)
                frlims_pb = fit_pb%loop_lims(3)
                ! fill to the FULL padded band like the CPU generator (the working-band cap
                ! acts downstream through fpl%nyq, not through the stored values)
                nyqpd_pb  = fit_pb%get_lfny(1)
                call flex_gpu_prep_begin_f(build%lmsk, params%box, params%boxpd, &
                    &MAXIMGBATCHSZ, cov_image_mask_radius(params), .true.)
                allocate(ctfp_arr(MAXIMGBATCHSZ), shf_pb(2,MAXIMGBATCHSZ))
                if( params%l_ml_reg )then
                    block
                        type(ftiter) :: fit_cr
                        call fit_cr%new([params%box_crop, params%box_crop, 1], params%smpd_crop)
                        signyq_pb = fit_cr%get_lfny(1)
                    end block
                    nyqfull_pb = fit_pb%get_lfny(1)
                    allocate(sig2_pb(0:nyqfull_pb, MAXIMGBATCHSZ), source=1.0)
                endif
                write(logfhandle,'(A)') '>>> FLEX_PCA PROBE DEVICE PREP ON (planes born on device)'
                call flush(logfhandle)
            endif
            allocate(o_pb_thr(nthr))
            allocate(cap1(MAXIMGBATCHSZ), source=1.0)
            allocate(sap0(MAXIMGBATCHSZ), source=0.0)
            allocate(seg_dir_b(MAXIMGBATCHSZ), seg_beg_b(MAXIMGBATCHSZ), seg_cnt_b(MAXIMGBATCHSZ))
            allocate(rmat_pp(3,3,npp), nrm_pp(3,npp), dirof_pb(npp), cap(npp), sap(npp))
            do i = 1, npp
                call build%spproj_field%get_ori(ppinds(i), o_pb)
                rmat_pp(:,:,i) = o_pb%get_mat()
                nrm_pp(:,i)    = rmat_pp(3,:,i)
            end do
            call o_pb%kill
            call polar_assign_directions(nrm_pp, npp, nrm_pb, ndir_pb, dirof_pb)
            do i = 1, npp
                call polar_relative_inplane(rmat_pp(:,:,i), rmat_pb(:,:,dirof_pb(i)), ca1, sa1)
                cap(i) = ca1; sap(i) = sa1
            end do
            deallocate(rmat_pp, nrm_pp)
            ! stable counting sort by bank direction: contiguous same-direction runs inside each
            ! batch are what turn per-particle splats into per-direction splats
            allocate(cnt_pb(ndir_pb), source=0)
            do i = 1, npp
                cnt_pb(dirof_pb(i)) = cnt_pb(dirof_pb(i)) + 1
            end do
            allocate(perm_pb(npp))
            block
                integer, allocatable :: pos(:)
                allocate(pos(ndir_pb)); pos(1) = 1
                do jdir = 2, ndir_pb
                    pos(jdir) = pos(jdir-1) + cnt_pb(jdir-1)
                end do
                do i = 1, npp
                    perm_pb(pos(dirof_pb(i))) = i
                    pos(dirof_pb(i)) = pos(dirof_pb(i)) + 1
                end do
                deallocate(pos)
            end block
            allocate(pptmp(npp))
            pptmp = ppinds;   ppinds   = pptmp(perm_pb)
            pptmp = dirof_pb; dirof_pb = pptmp(perm_pb)
            deallocate(pptmp)
            block
                real, allocatable :: rtmp(:)
                allocate(rtmp(npp))
                rtmp = cap; cap = rtmp(perm_pb)
                rtmp = sap; sap = rtmp(perm_pb)
                deallocate(rtmp)
            end block
            deallocate(perm_pb, cnt_pb)
            write(logfhandle,'(A,I0,A,I0,A)') '>>> FLEX_PCA PROBE BANKED ADJOINT ON (', &
                &ndir_pb,' bank directions, ',npp,' records)'
            call flush(logfhandle)
        endif
        if( .not. allocated(sec_proj_thr) )then
            allocate(sec_proj_thr(nthr), sec_gram_thr(nthr), source=0.d0)
        endif
        nll_prev    = 0.d0
        eo_best     = -1.d0
        eo_stall    = 0
        eo_patience = COV_EO_PATIENCE
        vpat        = 0
        call cov_env_int('SIMPLE_COV_EO_PATIENCE', vpat)
        if( vpat > 0 ) eo_patience = vpat
        ! ON by default. The consensus-shaped term is one the model provably cannot represent:
        ! the contrast coefficient is pinned at 1, so a particle of true amplitude a_i leaves
        ! (a_i - 1) * T_i P(R_i) mu behind, rank one along the consensus and identical in every
        ! particle. Deflating it costs a projection per basis volume per iteration and was
        ! positive in all four cells measured at K=20 on EMPIAR-10076 -- moment 0.1756 -> 0.1906,
        ! EM 0.1193 -> 0.1764, and positive again with latent whitening layered on both arms.
        ! It removes 34.6% of the basis energy on the first EM iteration and 16.2% at convergence.
        ! Validated on 10076 only so far; SIMPLE_COV_EM_DEFLATE=0 restores the old behaviour.
        vdfl        = 1
        call cov_env_int('SIMPLE_COV_EM_DEFLATE', vdfl)
        ! cov_env_int only accepts values > 0, so it cannot express "off" and =0 would silently
        ! leave the default standing. The explicit off-check is what makes the escape hatch real.
        l_deflate_mean = vdfl > 0 .and. .not. cov_env_int_off('SIMPLE_COV_EM_DEFLATE')
        ! ---- per-particle contrast INSIDE the EM (a_i in the basis loop) ----
        ! The probe has always fit a scale a_i = <m,y>/||m||^2 per particle per iteration and
        ! subtracted a_i*(T mu) from the residual -- the claim that the EM path has no scale
        ! parameter was about a different routine. What the historical path does NOT do:
        !   1. refine a_i against the basis: the projection fit ignores B z, which biases a_i
        !      wherever the basis is not orthogonal to the consensus (deflation removes most of
        !      that, which is one reason the two compose), and
        !   2. carry a_i into the M-step: the residual y - a m ~ a B z + n is inserted with
        !      weight z and density E[zz'], where the ML weighting is a z and a^2 E[zz'].
        ! 3DVA fits its alpha_i jointly in the M-step, which is both of these at once.
        ! SIMPLE_COV_PROBE_CONTRAST=n enables n ECM alternations per particle (fixes 1);
        ! SIMPLE_COV_PROBE_MLSCALE=1 enables the consistent M-step weighting (fixes 2).
        ! Both off by default; the polar statistics path keeps the fixed polar-fit contrast.
        n_probe_cm = 0
        call cov_env_int('SIMPLE_COV_PROBE_CONTRAST', n_probe_cm)
        n_probe_cm = max(0, min(8, n_probe_cm))
        nml_plain = n_probe_cm
        if( l_em_polar ) nml_plain = 0
        l_probe_mls = cov_env_int_on('SIMPLE_COV_PROBE_MLSCALE')
        if( n_probe_cm > 0 ) write(logfhandle,'(A,I0,A)') &
            &'>>> FLEX_PCA PROBE per-particle contrast ECM ON (', n_probe_cm, ' alternations)'
        if( l_probe_mls ) write(logfhandle,'(A)') &
            &'>>> FLEX_PCA PROBE ML-scaled M-step ON (a z insertion, a^2 density)'
        vanneal   = 0
        call cov_env_int('SIMPLE_COV_EM_ANNEAL', vanneal)
        kfr_ann   = covariance_kfromto(params)
        khi_full  = max(1, kfr_ann(2))
        dstep_ann = real(max(1, params%box_crop - 1)) * params%smpd_crop
        conv_thresh = COV_PROBE_CONV
        vconv       = 0
        call cov_env_int('SIMPLE_COV_PROBE_CONV', vconv)
        if( vconv > 0 ) conv_thresh = min(0.999999d0, real(vconv,dp)/1000.d0)
        ! ---- DOWNSCALED-PARTICLE CACHE ----
        ! Every EM iteration re-reads and re-preps the SAME particles, and one qsys round is
        ! launched per iteration, so the workers are fresh processes each time and nothing held in
        ! memory survives them. The on-disk cache does survive: it stores the iteration-independent
        ! prefix (noise normalisation, FFT, crop to box_crop), which at box=360/box_crop=64 is where
        ! the measured 48% of probe time goes. Masking and the device prep both keep the full box.
        l_pcache = ptcl_cache_in_use(params, build) .and. .not. l_gpu_prep .and. &
            &cov_image_mask_radius(params) <= 0.0
        if( l_pcache )then
            write(logfhandle,'(A)') '>>> FLEX_PCA PROBE: particles served from the downscaled cache'
            call flush(logfhandle)
        endif
        do it = 1, niters
            lp_it = params%lp
            if( vanneal > 0 )then
                khi_it = max(1, min(khi_full, 1 + ceiling(real(khi_full - 1) * &
                    &min(1.0, real(it) / real(vanneal)))))
                lp_it  = max(params%lp, dstep_ann / real(khi_it))
                write(logfhandle,'(A,I0,A,I0,A,I0,A,F7.2)') '>>> FLEX_PCA PROBE ITER ',it, &
                    &'  band anneal k_hi=',khi_it,' / ',khi_full,'  lp=',lp_it
            endif
            t_it = tic()
            sec_read = 0.; sec_prep = 0.; sec_estep = 0.; sec_ins = 0.
            sec_proj_thr = 0.d0; sec_gram_thr = 0.d0
            sec_bank = 0.; l_pol_bank_it = .false.   ! the polar E-step bank is rebuilt every iteration
            write(logfhandle,'(A,I0,A,I0,A,I0)') '>>> FLEX_PCA PROBE SUBSPACE ITERATION ',it,' / ',niters, &
                &'  basis dim=',ncomp
            call flush(logfhandle)
            allocate(prior(ncomp))
            do q = 1, ncomp
                prior(q) = 1.d0 / max(eigvals(q), DTINY)
            end do
            ! MCFA bookkeeping. l_mix_used freezes which prior THIS iteration's E-step runs
            ! under (the M-step hook below flips l_mix_active mid-iteration); a basis dimension
            ! change invalidates xi/Omega, so the mixture re-warms and re-initialises.
            l_mix_used = l_mix_active
            ldOm_used  = ldOm_mix
            if( l_mix_req )then
                if( allocated(rhs0th) )then
                    if( size(rhs0th,1) /= ncomp )then
                        deallocate(rhs0th, mkth, lwth, rkth, mxa_sr, mxa_sm, mxa_smm, mxa_sainv)
                        if( allocated(mix_xi) ) deallocate(mix_xi, mix_Om, mix_Ominv, mix_pi, &
                            &mix_Omxi, mix_xiOx, mix_lpi)
                        if( l_mix_active ) write(logfhandle,'(A)') &
                            &'>>> FLEX_PCA MIX re-initialising: basis dimension changed'
                        l_mix_active = .false.
                        l_mix_used   = .false.
                    endif
                endif
                if( .not. allocated(rhs0th) )then
                    allocate(rhs0th(ncomp,nthr), mkth(ncomp,kmix_pr,nthr), lwth(kmix_pr,nthr), &
                        &rkth(kmix_pr,nthr), mxa_sr(kmix_pr,nthr), mxa_sm(ncomp,kmix_pr,nthr), &
                        &mxa_smm(ncomp,ncomp,kmix_pr,nthr), mxa_sainv(ncomp,ncomp,nthr))
                endif
                mxa_sr = 0.d0; mxa_sm = 0.d0; mxa_smm = 0.d0; mxa_sainv = 0.d0
            endif
            ! even/odd Y_q accumulators (half-set FSC regularization) + the COUPLED latent normal matrix.
            ! rho carries one entry per (q,r) pair, not one shared density: the M-step below solves the
            ! components together at every grid point.
            allocate(Yeven(ncomp), Yodd(ncomp))
            do q = 1, ncomp
                call init_basis_reconstructor(params, build, Yeven(q)); call Yeven(q)%reset; call Yeven(q)%reset_exp
                call init_basis_reconstructor(params, build, Yodd(q));  call Yodd(q)%reset;  call Yodd(q)%reset_exp
            end do
            es     = shape(Yeven(1)%cmat_exp)
            ! The per-voxel coupled normal matrix is always the FULL ncomp x ncomp SPD system.
            ! A diagonal approximation used to be selectable here; it dropped the cross-component
            ! terms, and on EMPIAR-10028 that cost the 40S rotation outright -- the mode only
            ! appeared once the full solve was restored. It is not a speed/accuracy trade worth
            ! offering, so the option is gone rather than merely defaulted off.
            npairs = (ncomp*(ncomp+1))/2
            allocate(rho_e(npairs,es(1),es(2),es(3)), rho_o(npairs,es(1),es(2),es(3)), source=0.)
            write(logfhandle,'(A,I0,A,F8.2,A)') '>>> FLEX_PCA PROBE coupled normal matrix rows=',npairs, &
                &' (full)  rho even+odd ', 8.d0*real(npairs,dp)*real(es(1),dp)*real(es(2),dp)*real(es(3),dp)/1.d9,' GB'
            call flush(logfhandle)
            allocate(Gth(ncomp,ncomp,nthr), Ath(ncomp,ncomp,nthr), bth(ncomp,nthr), cth(ncomp,nthr), zth(ncomp,nthr))
            allocate(Ainvth(ncomp,ncomp,nthr), Acpth(ncomp,ncomp,nthr))
            allocate(basis_fpls(ncomp,nthr), mean_fpl(nthr), orientations(MAXIMGBATCHSZ))
            if( l_bank_fwd ) allocate(basis_fplsC(ncomp,nthr), mean_fplC(nthr))
            allocate(zbatch(ncomp,MAXIMGBATCHSZ), dens(ncomp,ncomp,MAXIMGBATCHSZ))
            allocate(valid(MAXIMGBATCHSZ), valid_e(MAXIMGBATCHSZ), valid_o(MAXIMGBATCHSZ), eo(MAXIMGBATCHSZ))
            allocate(gam_thr(ncomp,nthr), source=0.d0)
            allocate(gam_acc(ncomp), source=0.d0)
            allocate(nval_thr(nthr), source=0)
            allocate(hth(ncomp,nthr), source=0.d0)
            allocate(nll_thr(nthr), source=0.d0)
            allocate(gam_dbg(4,nthr), source=0.d0)
            nll_tot = 0.d0
            dens = 0.d0
            call mean_rec%expand_exp
            do q = 1, ncomp
                call basis_recs(q)%expand_exp
            end do
            ! ---- POLAR E-STEP (SIMPLE_COV_EM_POLAR=1) ----
            ! Must run BEFORE init_rec: embed_accumulate_polar ends with cleanup_rec_buffers, which
            ! calls killimgbatch / dealloc_imgarr(img_pad_heap) / forget_ft_maps -- shared BUILD
            ! state, not just its own planes -- so running it after init_rec pulls the batch loop's
            ! buffers out from under it and segfaults on the first insertion.
            ! One pass through the shared-direction bank fills every particle's G, b, c and
            ! contrast; the batch loop below then projects the MEAN ONLY, to form the Cartesian
            ! residual the M-step adjoint needs. There is no polar->volume adjoint in the tree, so
            ! the M-step stays Cartesian by necessity.
            if( l_em_polar )then
                if( l_probe_distr_pre ) THROW_HARD('SIMPLE_COV_EM_POLAR is shared-memory only for &
                    &now; the polar former would have to run inside the worker. Use nparts=1.')
                allocate(pGc(ncomp,ncomp,npp), source=0.d0)
                allocate(pbc(ncomp,npp), pcc(ncomp,npp), source=0.d0)
                allocate(pzh(npp,ncomp,2), source=0.d0)
                allocate(pcon(npp), pre(npp), prme(npp), source=0.d0)
                call embed_accumulate_polar(params, build, mean_rec, basis_recs, ncomp, eigvals, &
                    &sig2, ppinds, npp, pGc, pbc, pcc, pzh, pcon, pre, prme, prior, nvalid_pol)
                write(logfhandle,'(A,I0,A,I0)') '>>> FLEX_PCA PROBE POLAR E-STEP: statistics from &
                    &the shared-direction bank for ',nvalid_pol,' of ',npp
                call flush(logfhandle)
            endif
            call init_rec(params, build, MAXIMGBATCHSZ, fpls, cached=l_pcache)
            ! ---- ACCUMULATE: distributed when the master has parts, in-process otherwise ----
            ! One qsys round per EM iteration, because the basis the E-step projects changes every
            ! iteration: workers are relaunched against the master's refreshed flex_pca_pc*.mrc
            ! rather than looping locally. Everything below the reduction -- the coupled M-step
            ! solve, the FSC merge, the re-orthonormalisation -- stays master-only and unchanged, so
            ! the shared-memory result remains the reference.
            l_probe_distr = flex_pca_is_master() .and. flex_pca_nparts() > 1
            allocate(gam_sum(ncomp), source=0.d0)
            nval = 0
            if( l_probe_distr )then
                call save_probe_state(ncomp, eigvals, sig2_eff)
                call flex_pca_run_stage(PCA_STAGE_PROBE, 'probe iteration')
                allocate(cme(es(1),es(2),es(3),ncomp), cmo(es(1),es(2),es(3),ncomp), source=(0.,0.))
                allocate(rhe(es(1),es(2),es(3),ncomp), rhoo(es(1),es(2),es(3),ncomp), source=0.)
                if( l_mix_req )then
                    if( .not. allocated(dm_sr) )then
                        allocate(dm_sr(kmix_pr), dm_sm(ncomp,kmix_pr), &
                            &dm_smm(ncomp,ncomp,kmix_pr), dm_sai(ncomp,ncomp))
                        allocate(dm_z(MIX_ZSUB_MAX*flex_pca_nparts(), ncomp))
                    endif
                    dm_sr = 0.d0; dm_sm = 0.d0; dm_smm = 0.d0; dm_sai = 0.d0; dm_nz = 0
                    call reduce_probe_parts(params, flex_pca_nparts(), cme, rhe, cmo, rhoo, &
                        &rho_e, rho_o, gam_sum, nll_tot, nval, ncomp, &
                        &dm_sr, dm_sm, dm_smm, dm_sai, dm_z, dm_nz)
                    write(logfhandle,'(A,I0,A,I0)') '>>> FLEX_PCA MIX distributed reduce: kmix=', &
                        &kmix_pr,'  pooled latent subsample rows=',dm_nz
                    call flush(logfhandle)
                else
                    call reduce_probe_parts(params, flex_pca_nparts(), cme, rhe, cmo, rhoo, &
                        &rho_e, rho_o, gam_sum, nll_tot, nval, ncomp)
                endif
                do q = 1, ncomp
                    Yeven(q)%cmat_exp = cme(:,:,:,q); Yeven(q)%rho_exp = rhe(:,:,:,q)
                    Yodd(q)%cmat_exp  = cmo(:,:,:,q); Yodd(q)%rho_exp  = rhoo(:,:,:,q)
                end do
                deallocate(cme, cmo, rhe, rhoo)
            else
            if( l_pcache )then
                call prepimgbatch(params, build, MAXIMGBATCHSZ, box=params%box_crop, smpd=params%smpd_crop)
            else
                call prepimgbatch(params, build, MAXIMGBATCHSZ)
            endif
            vgpu = 0
            call cov_env_int('SIMPLE_COV_GPU', vgpu)
            l_gpu_probe  = vgpu > 0 .and. flex_gpu_available()
            if( l_gpu_probe )then
                write(logfhandle,'(A,I0,A,I0,A)') '>>> FLEX_PCA PROBE M-step GPU insertion ON (', &
                    &2*ncomp, ' components, ', 2*npairs, ' density rows)'
                call flush(logfhandle)
                call flex_gpu_coupled_begin_f(Yeven, 2*ncomp, 2*npairs)
                allocate(gdsc(2*ncomp,MAXIMGBATCHSZ), grsc(2*npairs,MAXIMGBATCHSZ), source=0.d0)
                if( l_bank_adj ) call flex_gpu_coupled_bank_f(rmat_pb, ndir_pb)
                if( l_gpu_estep .and. l_bank_adj )then
                    call flex_gpu_estep_vols_f(mean_rec, basis_recs, ncomp)
                    nst_es = 2 + 2*ncomp + (ncomp*(ncomp+1))/2
                    allocate(stats_es(nst_es, MAXIMGBATCHSZ), source=0.d0)
                    write(logfhandle,'(A,I0,A)') '>>> FLEX_PCA PROBE FUSED E-STEP ON (', &
                        &ncomp+1, ' volumes resident, exact directions)'
                    call flush(logfhandle)
                endif
            endif
            z = 0.d0
            ! parity-check bookkeeping: only iteration 1, only when the polar former runs
            l_chk_act = l_pol_es .and. n_pol_check > 0 .and. it == 1
            if( l_chk_act .and. .not. allocated(zchk_c) )then
                nchk_max = min(n_pol_check, npp)
                allocate(zchk_c(ncomp,nchk_max), zchk_p(ncomp,nchk_max), source=0.d0)
                allocate(gerr_chk(nchk_max), berr_chk(nchk_max), source=0.d0)
                allocate(lchk(nchk_max), source=.false.)
                allocate(Gc_chk(ncomp,ncomp,nthr), bc_chk(ncomp,nthr), cc_chk(ncomp,nthr))
                allocate(zc_chk(ncomp,nthr), Ai_chk(ncomp,ncomp,nthr))
                allocate(zchk_d(ncomp,nchk_max), source=0.d0)
                allocate(gerrd_chk(nchk_max), berrd_chk(nchk_max), source=0.d0)
                allocate(Gd_chk(ncomp,ncomp,nthr), bd_chk(ncomp,nthr), cd_chk(ncomp,nthr))
                allocate(zd_chk(ncomp,nthr), dcb_chk(ncomp,nthr))
            endif
            ob_lo = 1
            ob_hi = npp
            if( l_online )then
                b_win     = max(1, b_online / max(1, params%nparts))   ! params%nparts is worker-visible; flex_pca_nparts() is master-only
                nb_online = max(1, (npp + b_win - 1) / b_win)
                it_win    = mod(it - 1, nb_online)
                ob_lo     = 1 + it_win * b_win
                ob_hi     = min(npp, ob_lo + b_win - 1)
                if( l_onl_final .and. it == niters )then
                    ob_lo = 1
                    ob_hi = npp
                    write(logfhandle,'(A,I0,A,I0,A)') '>>> FLEX_PCA ONLINE FINAL it=', it, &
                        &'  full-set pass over ', npp, ' particles'
                else
                    write(logfhandle,'(A,I0,A,I0,A,I0,A,I0)') '>>> FLEX_PCA ONLINE it=', it, &
                        &'  window ', ob_lo, '..', ob_hi, '  of ', npp
                endif
                call flush(logfhandle)
            endif
            do ibatch = ob_lo, ob_hi, MAXIMGBATCHSZ
                batchlims = [ibatch, min(ob_hi, ibatch + MAXIMGBATCHSZ - 1)]
                batchsz   = batchlims(2) - batchlims(1) + 1
                t_sec = tic()
                if( l_pcache )then
                    call ptcl_cache_read_batch(params, build, npp, ppinds, batchlims)
                else
                    call discrete_read_imgbatch(params, build, npp, ppinds, batchlims)
                endif
                sec_read = sec_read + toc(t_sec)
                t_sec = tic()
                if( .not. l_gpu_prep )then
                    call prep_imgs4projected_model(params, build, batchsz, build%imgbatch(:batchsz), &
                        &ppinds(batchlims(1):batchlims(2)), fpls(:batchsz), &
                        &mskrad=cov_image_mask_radius(params), cached=l_pcache)
                endif
                sec_prep = sec_prep + toc(t_sec)
                do i = 1, batchsz
                    call build%spproj_field%get_ori(ppinds(batchlims(1)+i-1), orientations(i))
                    eo(i) = build%spproj_field%get_eo(ppinds(batchlims(1)+i-1))
                    if( l_gpu_prep )then
                        ctfp_arr(i)  = build%spproj%get_ctfparams(params%oritype, ppinds(batchlims(1)+i-1))
                        shf_pb(:,i)  = build%spproj_field%get_2Dshift(ppinds(batchlims(1)+i-1))
                        if( params%l_ml_reg ) call resample_sigma2(kf_pb(1), signyq_pb, &
                            &build%esig%sigma2_noise(kf_pb(1):kf_pb(2), ppinds(batchlims(1)+i-1)), &
                            &nyqfull_pb, real(signyq_pb)/real(nyqfull_pb), sig2_pb(:,i))
                    endif
                end do
                ! ---- POLAR E-STEP BANK, built once per EM iteration at the first prepped batch
                ! (basis_recs change every M-step, so the bank cannot be cached across iterations;
                ! the grid and the pose-fixed direction assignment are built once per stage) ----
                if( l_pol_es .and. .not. l_pol_bank_it )then
                    t_bank = tic()
                    if( .not. l_pol_grid )then
                        ! grid geometry from the first prepped plane: the identical derivation
                        ! embed_accumulate_polar uses, so polar embed and polar E-step share
                        ! band/quadrature conventions. No noise rings: sig2 arrives as sig2_eff.
                        ph0_es  = lbound(fpls(1)%cmplx_plane,1)
                        pk0_es  = lbound(fpls(1)%cmplx_plane,2)
                        hlo_es  = ceil_div (lbound(fpls(1)%cmplx_plane,1), OSMPL_PAD_FAC)
                        hhi_es  = floor_div(ubound(fpls(1)%cmplx_plane,1), OSMPL_PAD_FAC)
                        klo_es  = ceil_div (lbound(fpls(1)%cmplx_plane,2), OSMPL_PAD_FAC)
                        nyqr_es = mean_rec%get_lfny(1)
                        nyqb_es = nyqr_es
                        if( fpls(1)%nyq > 0 ) nyqb_es = min(nyqb_es, max(1, fpls(1)%nyq / OSMPL_PAD_FAC))
                        ! hybrid split point: auto at 0.72*band -- the measured knee of the
                        ! real-data ladder (10049 gate, band 11: rhyb 6 -> b err 8.9%, min z
                        ! corr 0.949, lambda head 153/282/326; rhyb 8 -> b err 3.1%, min z
                        ! corr 0.990, lambda 166/273/394 vs Cartesian 166/264/382). Env
                        ! override; explicit 0 = pure rings (the pre-hybrid baseline).
                        rhyb_es = nint(0.72*real(nyqb_es))
                        if( vpol_rhyb > 0 ) rhyb_es = vpol_rhyb
                        if( l_rhyb_off )    rhyb_es = 0
                        rhyb_es   = max(0, min(rhyb_es, nyqb_es-1))
                        l_pol_hyb = rhyb_es > 0
                        if( l_pol_hyb )then
                            call polar_grid_build(pg_es, rhyb_es+1, nyqb_es, nyqb_es+1, nyqb_es, &
                                &hlo_es, hhi_es, klo_es, ph0_es, pk0_es, ang_osamp=vpol_osamp, &
                                &gate_lo=rhyb_es*(rhyb_es+1))
                            ! exact-part lattice positions, cov_herm_inner's half-plane rule
                            ! (k<=0; on the k=0 line only h<=0), shells 0..rhyb by the nint
                            ! convention: h^2+k^2 <= rhyb*(rhyb+1). Raster order k-outer.
                            npos_es = 0
                            do kx_es = -rhyb_es, 0
                                do jx_es = -rhyb_es, merge(0, rhyb_es, kx_es == 0)
                                    if( jx_es*jx_es + kx_es*kx_es > rhyb_es*(rhyb_es+1) ) cycle
                                    npos_es = npos_es + 1
                                end do
                            end do
                            allocate(hex_es(npos_es), kex_es(npos_es))
                            npos_es = 0
                            do kx_es = -rhyb_es, 0
                                do jx_es = -rhyb_es, merge(0, rhyb_es, kx_es == 0)
                                    if( jx_es*jx_es + kx_es*kx_es > rhyb_es*(rhyb_es+1) ) cycle
                                    npos_es = npos_es + 1
                                    hex_es(npos_es) = jx_es
                                    kex_es(npos_es) = kx_es
                                end do
                            end do
                            write(logfhandle,'(A,I0,A,I0,A,I0,A,I0)') '>>> FLEX_PCA POLAR HYBRID: exact &
                                &Cartesian statistics for shells 0..',rhyb_es,' (',npos_es,' lattice &
                                &points/particle incl. DC), rings ',rhyb_es+1,'..',nyqb_es
                            if( vpol_osamp > 1 ) write(logfhandle,'(A,I0)') &
                                &'>>> FLEX_PCA POLAR OSAMP: ring angular oversampling x', vpol_osamp
                            call flush(logfhandle)
                        else
                            call polar_grid_build(pg_es, 1, nyqb_es, nyqb_es+1, nyqb_es, hlo_es, hhi_es, &
                                &klo_es, ph0_es, pk0_es, ang_osamp=vpol_osamp)
                            if( vpol_osamp > 1 )then
                                write(logfhandle,'(A,I0)') &
                                    &'>>> FLEX_PCA POLAR OSAMP: ring angular oversampling x', vpol_osamp
                                call flush(logfhandle)
                            endif
                        endif
                        nsamp_es = pg_es%nsamp; nsamp2_es = 2*nsamp_es; nk_es = pg_es%nk
                        ! shared direction table: same refspiral + count derivation as the polar embed
                        ndir_es = cov_polar_ndir(npp)
                        call dirs_es%new(ndir_es, is_ptcl=.false.)
                        call build%pgrpsyms%build_refspiral(dirs_es)
                        allocate(rmatb_es(3,3,ndir_es), nrmb_es(3,ndir_es))
                        do id_es = 1, ndir_es
                            rmatb_es(:,:,id_es) = dirs_es%get_mat(id_es)
                            nrmb_es(:,id_es)    = rmatb_es(3,:,id_es)
                        end do
                        call dirs_es%kill
                        ! pose-fixed per-particle (direction, in-plane) assignment, once per stage
                        allocate(rmatp_es(3,3,npp), nrmp_es(3,npp), dir_es(npp), cae(npp), sae(npp))
                        do i = 1, npp
                            call build%spproj_field%get_ori(ppinds(i), o_es)
                            rmatp_es(:,:,i) = o_es%get_mat()
                            nrmp_es(:,i)    = rmatp_es(3,:,i)
                        end do
                        call o_es%kill
                        call polar_assign_directions(nrmp_es, npp, nrmb_es, ndir_es, dir_es)
                        do i = 1, npp
                            call polar_relative_inplane(rmatp_es(:,:,i), rmatb_es(:,:,dir_es(i)), ca1, sa1)
                            cae(i) = ca1; sae(i) = sa1
                        end do
                        deallocate(rmatp_es, nrmp_es)
                        allocate(dused_es(ndir_es))
                        ! device ring geometry, once per stage (the grid is pose- and band-fixed)
                        if( l_pol_gpu ) call flex_gpu_poles_begin_f(pg_es%rad, pg_es%cs, &
                            &pg_es%sn, pg_es%sqwq, pg_es%rbeg, pg_es%rend, nsamp_es, nk_es)
                        l_pol_grid = .true.
                        write(logfhandle,'(A,I0,A,I0,A,I0,A,F8.1,A,F8.1,A)') &
                            &'>>> FLEX_PCA POLAR ESTEP BANK: ',ncomp+1,' volumes x ',ndir_es, &
                            &' directions x ',nsamp_es,' ring samples = ', &
                            &4.d0*real(nsamp2_es,dp)*real(ncomp+1,dp)*real(ndir_es,dp)/1.d6, &
                            &' MB (+ ring tables ', &
                            &8.d0*real(ncomp*ncomp+ncomp+1,dp)*real(nk_es,dp)*real(ndir_es,dp)/1.d6,' MB)'
                        call flush(logfhandle)
                    endif
                    ! (re)allocate at this iteration's rank (ncomp can change between iterations)
                    if( allocated(UsallE) )then
                        if( size(UsallE,2) /= ncomp+1 ) deallocate(UsallE, CfE, Cm0E, c00E, &
                            &UbankE, CspE, xws_es, wr_es, wrd_es, Reb_es)
                    endif
                    if( .not. allocated(UsallE) )then
                        allocate(UsallE(nsamp2_es,0:ncomp,ndir_es), CfE(ncomp*ncomp,nk_es,ndir_es), &
                            &Cm0E(ncomp,nk_es,ndir_es), c00E(nk_es,ndir_es))
                        allocate(UbankE(nsamp_es,0:ncomp,nthr), CspE(0:ncomp,0:ncomp,nthr))
                        allocate(xws_es(nsamp2_es,nthr), wr_es(nk_es,nthr), wrd_es(nk_es,nthr), &
                            &Reb_es(0:ncomp,nthr))
                    endif
                    ! GPU parity arrays (gate 1): the CPU polar former re-run on the check
                    ! particles; the exact-direction/DC-less attribution arms stand down under
                    ! the device path (they answer the polar-vs-Cartesian question, already
                    ! measured on this recipe; the device gate is GPU-vs-CPU-polar)
                    if( l_chk_act .and. l_pol_gpu .and. .not. allocated(zchk_g) )then
                        allocate(zchk_g(ncomp,nchk_max), source=0.d0)
                        allocate(gerrg_chk(nchk_max), berrg_chk(nchk_max), source=0.d0)
                        allocate(Gpc_chk(ncomp,ncomp,nthr), bpc_chk(ncomp,nthr), &
                            &cpc_chk(ncomp,nthr), zpc_chk(ncomp,nthr))
                    endif
                    ! exact-direction check arm scratch (grid-shaped, so allocated here)
                    if( l_chk_act .and. .not. l_pol_gpu .and. .not. allocated(UsX_chk) )then
                        allocate(zchk_x(ncomp,nchk_max), source=0.d0)
                        allocate(gerrx_chk(nchk_max), berrx_chk(nchk_max), source=0.d0)
                        allocate(UsX_chk(nsamp2_es,0:ncomp,nthr), xwsX_chk(nsamp2_es,nthr), &
                            &wrX_chk(nk_es,nthr), RebX_chk(0:ncomp,nthr))
                        allocate(CfX_chk(ncomp*ncomp,nk_es,nthr), Cm0X_chk(ncomp,nk_es,nthr), &
                            &c00X_chk(nk_es,nthr), wrdX_chk(nk_es,nthr))
                        allocate(GX_chk(ncomp,ncomp,nthr), bX_chk(ncomp,nthr), cX_chk(ncomp,nthr), &
                            &zx_chk(ncomp,nthr))
                        allocate(tazsum_thr(nthr), source=0.d0)
                        allocate(tazcnt_thr(nthr), source=0)
                    endif
                    ! only the directions this iteration's window touches
                    dused_es = .false.
                    do i = ob_lo, ob_hi
                        if( dir_es(i) > 0 ) dused_es(dir_es(i)) = .true.
                    end do
                    !$omp parallel do default(shared) schedule(dynamic) proc_bind(close) &
                    !$omp& private(id_es,ithr,q,r,ir)
                    do id_es = 1, ndir_es
                        if( .not. dused_es(id_es) ) cycle
                        ithr = omp_get_thread_num() + 1
                        call polar_project_recs(mean_rec, basis_recs, ncomp, rmatb_es(:,:,id_es), &
                            &pg_es, UbankE(:,:,ithr))
                        do q = 0, ncomp
                            do r = 1, nsamp_es
                                UsallE(2*r-1,q,id_es) = pg_es%sqwq(r)*real (UbankE(r,q,ithr))
                                UsallE(2*r,  q,id_es) = pg_es%sqwq(r)*aimag(UbankE(r,q,ithr))
                            end do
                        end do
                        do ir = 1, nk_es
                            call polar_ring_gram(UsallE(1,0,id_es), nsamp2_es, ncomp, pg_es%rbeg(ir), &
                                &pg_es%rend(ir)-pg_es%rbeg(ir)+1, CspE(0,0,ithr), CfE(1,ir,id_es), &
                                &Cm0E(1,ir,id_es))
                            c00E(ir,id_es) = polar_ring_selfpower(UsallE(1,0,id_es), nsamp2_es, &
                                &pg_es%rbeg(ir), pg_es%rend(ir)-pg_es%rbeg(ir)+1)
                        end do
                    end do
                    !$omp end parallel do
                    ! ---- device upload (stage 2): the bank + ring tables once per iteration, and
                    ! the mean+basis volumes for the hybrid-exact kernel. The rank can change
                    ! between iterations, so the capacity gate is re-checked here; falling back
                    ! mid-run is seamless because every consumer below branches per iteration.
                    if( l_pol_gpu .and. l_pol_hyb .and. (ncomp + 1) > COV_GPU_ESTEP_MAXVOLS )then
                        l_pol_gpu = .false.
                        write(logfhandle,'(A,I0,A,I0,A)') '>>> FLEX_PCA POLAR E-STEP DEVICE OFF: needs ', &
                            &ncomp+1,' resident volumes, device entry point caps at ', &
                            &COV_GPU_ESTEP_MAXVOLS,' -- using the CPU polar E-step'
                        call flush(logfhandle)
                    endif
                    if( l_pol_gpu )then
                        call flex_gpu_poles_bank_f(UsallE, CfE, Cm0E, c00E, ncomp, ndir_es)
                        if( l_pol_hyb ) call flex_gpu_estep_vols_f(mean_rec, basis_recs, ncomp)
                        nst_pg = 2 + 2*ncomp + (ncomp*(ncomp+1))/2
                        if( allocated(stats_pg) )then
                            if( size(stats_pg,1) /= nst_pg ) deallocate(stats_pg)
                        endif
                        if( .not. allocated(stats_pg) ) allocate(stats_pg(nst_pg, MAXIMGBATCHSZ), &
                            &source=0.d0)
                        if( .not. allocated(vld_pg) ) allocate(vld_pg(MAXIMGBATCHSZ), &
                            &dir_pg(MAXIMGBATCHSZ), ca_pg(MAXIMGBATCHSZ), sa_pg(MAXIMGBATCHSZ))
                    endif
                    sec_bank      = sec_bank + toc(t_bank)
                    l_pol_bank_it = .true.
                    write(logfhandle,'(A,I0,A,I0,A,I0,A,F7.1)') '>>> FLEX_PCA POLAR ESTEP BANK it=', &
                        &it,'  directions built=',count(dused_es),' of ',ndir_es, &
                        &'  build seconds=',sec_bank
                    call flush(logfhandle)
                endif
                valid(:batchsz) = .false.
                zbatch(:,:batchsz) = 0.d0
                dens(:,:,:batchsz) = 0.d0
                t_sec = tic()
                if( l_gpu_estep .and. l_gpu_probe .and. l_bank_adj )then
                    ! ---- fused device E-step: stats on device at EXACT directions, tiny fetch,
                    ! host posterior solves, residual formed on device and fetched packed
                    twp0 = omp_get_wtime()
                    do i = 1, batchsz
                        vld_es(i) = .not. orientations(i)%isstatezero()
                    end do
                    if( l_gpu_prep )then
                        if( params%l_ml_reg )then
                            call flex_gpu_prep_batch_f(build%imgbatch(:batchsz), ctfp_arr(:batchsz), &
                                &shf_pb(:,:batchsz), vld_es(:batchsz), batchsz, params%box, &
                                &frlims_pb, nyqpd_pb, sig2_ups=sig2_pb(:,:batchsz))
                        else
                            call flex_gpu_prep_batch_f(build%imgbatch(:batchsz), ctfp_arr(:batchsz), &
                                &shf_pb(:,:batchsz), vld_es(:batchsz), batchsz, params%box, &
                                &frlims_pb, nyqpd_pb)
                        endif
                        ! one-time in-situ cross-check against the CPU prep on real particles
                        if( it == 1 .and. ibatch == 1 )then
                            vprep = 0
                            call cov_env_int('SIMPLE_COV_PREP_CHECK', vprep)
                            if( vprep > 0 )then
                                call prep_imgs4projected_model(params, build, batchsz, &
                                    &build%imgbatch(:batchsz), ppinds(batchlims(1):batchlims(2)), &
                                    &fpls(:batchsz), mskrad=cov_image_mask_radius(params), &
                                    &force_cpu=.true.)
                                call flex_gpu_prep_check_f(fpls(:batchsz), vld_es(:batchsz), batchsz)
                            endif
                        endif
                        call flex_gpu_estep_batch_res_f(orientations(:batchsz), batchsz, &
                            &frlims_pb, kf_pb(2), stats_es(:,:batchsz))
                    else
                        call flex_gpu_estep_batch_f(fpls(:batchsz), orientations(:batchsz), &
                            &vld_es(:batchsz), batchsz, stats_es(:,:batchsz))
                    endif
                    sec_proj_thr(1) = sec_proj_thr(1) + (omp_get_wtime() - twp0)
                    twp0 = omp_get_wtime()
                    !$omp parallel do default(shared) schedule(static) proc_bind(close) &
                    !$omp& private(i,ithr,row,q,r,off_es,a,aa,e_mm,myv,ldA,lok,qml)
                    do i = 1, batchsz
                        if( .not. vld_es(i) ) cycle
                        ithr = omp_get_thread_num() + 1
                        row  = batchlims(1) + i - 1
                        e_mm = stats_es(1,i)
                        myv  = stats_es(2,i)
                        a    = max(0.1d0, min(5.0d0, myv / max(e_mm, DTINY)))
                        aa   = a*a
                        do q = 1, ncomp
                            bth(q,ithr) = stats_es(2+q,i)
                            cth(q,ithr) = stats_es(2+ncomp+q,i)
                        end do
                        off_es = 2 + 2*ncomp
                        do q = 1, ncomp
                            do r = q, ncomp
                                off_es = off_es + 1
                                Gth(q,r,ithr) = stats_es(off_es,i)
                                Gth(r,q,ithr) = Gth(q,r,ithr)
                            end do
                        end do
                        if( l_mix_used )then
                            ! same shared MCFA solver as the plain body: the device fetched the
                            ! identical G/b/c sufficient statistics, only the host solve differs
                            call probe_solve_mix(ncomp, kmix_pr, Gth(:,:,ithr), bth(:,ithr), &
                                &cth(:,ithr), myv, e_mm, n_probe_cm, a, sig2, mix_Ominv, mix_Omxi, &
                                &mix_lpi, mix_xiOx, zth(:,ithr), Ainvth(:,:,ithr), dens(:,:,i), &
                                &ldA, lok, nll_mix_add, mxa_sr(:,ithr), mxa_sm(:,:,ithr), &
                                &mxa_smm(:,:,:,ithr), mxa_sainv(:,:,ithr))
                            if( lok ) nll_thr(ithr) = nll_thr(ithr) + nll_mix_add
                        else
                            call probe_solve_ecm(ncomp, Gth(:,:,ithr), bth(:,ithr), cth(:,ithr), &
                                &myv, e_mm, prior, sig2, n_probe_cm, a, zth(:,ithr), &
                                &Ainvth(:,:,ithr), ldA, lok, qml)
                            if( lok ) nll_thr(ithr) = nll_thr(ithr) + ldA - qml
                        endif
                        aa = a*a
                        z(row,:)    = zth(:,ithr)
                        zbatch(:,i) = zth(:,ithr)
                        if( .not. l_mix_used )then
                            do r = 1, ncomp
                                do q = 1, ncomp
                                    dens(q,r,i) = zth(q,ithr)*zth(r,ithr) + Ainvth(q,r,ithr)
                                end do
                            end do
                        endif
                        do q = 1, ncomp
                            gam_thr(q,ithr) = gam_thr(q,ithr) + dens(q,q,i)
                        end do
                        ! Gamma above reads the UNSCALED E[zz'] (prior on z, not on a z);
                        ! everything below feeds the M-step, which under the fitted scale
                        ! solves  r ~ a B z + n
                        if( l_probe_mls )then
                            zbatch(:,i) = a*zbatch(:,i)
                            dens(:,:,i) = aa*dens(:,:,i)
                        endif
                        nval_thr(ithr) = nval_thr(ithr) + 1
                        valid(i)       = .true.
                        avec_es(i)     = real(a)
                    end do
                    !$omp end parallel do
                    ! residual stays on device; the M-step consumes it via the resident entry
                    if( l_deflate )then
                        call flex_gpu_estep_resid_f(fpls(:batchsz), avec_es(:batchsz), &
                            &vld_es(:batchsz), batchsz, fetch=.false.)
                    endif
                    sec_gram_thr(1) = sec_gram_thr(1) + (omp_get_wtime() - twp0)
                elseif( l_pol_gpu )then
                    ! ---- POLAR E-STEP ON DEVICE (stage 2): per-batch ring + hybrid-exact
                    ! statistics in the fused stat layout; everything downstream of the
                    ! statistics -- solve, mixture, Gamma, windows, e/o, mean projection,
                    ! residual, insertion -- is the CPU polar path's, verbatim.
                    twp0 = omp_get_wtime()
                    do i = 1, batchsz
                        row = batchlims(1) + i - 1
                        vld_pg(i) = .not. orientations(i)%isstatezero()
                        dir_pg(i) = dir_es(row)
                        ca_pg(i)  = cae(row)
                        sa_pg(i)  = sae(row)
                    end do
                    call flex_gpu_poles_batch_f(fpls(:batchsz), orientations(:batchsz), &
                        &ca_pg(:batchsz), sa_pg(:batchsz), dir_pg(:batchsz), vld_pg(:batchsz), &
                        &batchsz, merge(rhyb_es*(rhyb_es+1), 0, l_pol_hyb), stats_pg(:,:batchsz))
                    sec_proj_thr(1) = sec_proj_thr(1) + (omp_get_wtime() - twp0)
                    twp0 = omp_get_wtime()
                    !$omp parallel do default(shared) schedule(dynamic) proc_bind(close) &
                    !$omp& private(i,ithr,row,q,r,off_es,a,aa,e_mm,myv,ldA,lok,qml,nll_mix_add) &
                    !$omp& private(idp_es,taz_es,l_chk_p,ichk,emmc_chk,myvc_chk,ac_chk,lda_chk) &
                    !$omp& private(qml_chk,lok_chk,gnum_chk,gden_chk,bnum_chk,bden_chk) &
                    !$omp& private(emmg_chk,myvg_chk,ag_chk)
                    do i = 1, batchsz
                        if( .not. vld_pg(i) ) cycle
                        ithr = omp_get_thread_num() + 1
                        row  = batchlims(1) + i - 1
                        e_mm = stats_pg(1,i)
                        myv  = stats_pg(2,i)
                        a    = max(0.1d0, min(5.0d0, myv / max(e_mm, DTINY)))
                        aa   = a*a
                        do q = 1, ncomp
                            bth(q,ithr) = stats_pg(2+q,i)
                            cth(q,ithr) = stats_pg(2+ncomp+q,i)
                        end do
                        off_es = 2 + 2*ncomp
                        do q = 1, ncomp
                            do r = q, ncomp
                                off_es = off_es + 1
                                Gth(q,r,ithr) = stats_pg(off_es,i)
                                Gth(r,q,ithr) = Gth(q,r,ithr)
                            end do
                        end do
                        l_chk_p = l_chk_act .and. row <= nchk_max
                        if( l_chk_p )then
                            ! Cartesian reference former (its mean plane doubles as the
                            ! residual's mean projection below)
                            call project_fplanes_mean_basis(mean_rec, basis_recs, orientations(i), &
                                &fpls(i), mean_fpl(ithr), basis_fpls(:,ithr), apply_ctf_amp=.true.)
                            emmc_chk = real(cov_herm_inner(mean_fpl(ithr), mean_fpl(ithr)), dp)
                            myvc_chk = real(cov_herm_inner(mean_fpl(ithr), fpls(i)), dp)
                            do q = 1, ncomp
                                bc_chk(q,ithr) = real(cov_herm_inner(basis_fpls(q,ithr), fpls(i)), dp)
                                cc_chk(q,ithr) = real(cov_herm_inner(basis_fpls(q,ithr), mean_fpl(ithr)), dp)
                                do r = q, ncomp
                                    Gc_chk(q,r,ithr) = real(cov_herm_inner(basis_fpls(q,ithr), &
                                        &basis_fpls(r,ithr)), dp)
                                    Gc_chk(r,q,ithr) = Gc_chk(q,r,ithr)
                                end do
                            end do
                            ac_chk = max(0.1d0, min(5.0d0, myvc_chk / max(emmc_chk, DTINY)))
                            call probe_solve_ecm(ncomp, Gc_chk(:,:,ithr), bc_chk(:,ithr), &
                                &cc_chk(:,ithr), myvc_chk, emmc_chk, prior, sig2, 0, ac_chk, &
                                &zc_chk(:,ithr), Ai_chk(:,:,ithr), lda_chk, lok_chk, qml_chk)
                            ichk = row
                            zchk_c(:,ichk) = zc_chk(:,ithr)
                            ! GPU-vs-Cartesian statistic errors (the polar-vs-Cartesian bar)
                            gnum_chk = 0.d0; gden_chk = 0.d0; bnum_chk = 0.d0; bden_chk = 0.d0
                            do r = 1, ncomp
                                bnum_chk = bnum_chk + (bth(r,ithr) - bc_chk(r,ithr))**2
                                bden_chk = bden_chk + bc_chk(r,ithr)**2
                                do q = 1, ncomp
                                    gnum_chk = gnum_chk + (Gth(q,r,ithr) - Gc_chk(q,r,ithr))**2
                                    gden_chk = gden_chk + Gc_chk(q,r,ithr)**2
                                end do
                            end do
                            gerr_chk(ichk) = sqrt(gnum_chk / max(gden_chk, DTINY))
                            berr_chk(ichk) = sqrt(bnum_chk / max(bden_chk, DTINY))
                            ! the CPU polar former, the production stage-1 path verbatim,
                            ! for the GPU parity gate (same bank, same tables, same hybrid)
                            idp_es = dir_es(row)
                            call polar_sample_particle_fused(fpls(i)%cmplx_plane, &
                                &fpls(i)%transfer_plane, pg_es, cae(row), sae(row), &
                                &xws_es(:,ithr), wr_es(:,ithr), taz_es)
                            wrd_es(:,ithr) = real(wr_es(:,ithr), dp)
                            call sgemv('T', nsamp2_es, ncomp+1, 1.0, UsallE(1,0,idp_es), &
                                &nsamp2_es, xws_es(1,ithr), 1, 0.0, Reb_es(0,ithr), 1)
                            call dgemv('N', ncomp*ncomp, nk_es, 1.d0, CfE(1,1,idp_es), &
                                &ncomp*ncomp, wrd_es(1,ithr), 1, 0.d0, Gpc_chk(1,1,ithr), 1)
                            call dgemv('N', ncomp, nk_es, 1.d0, Cm0E(1,1,idp_es), ncomp, &
                                &wrd_es(1,ithr), 1, 0.d0, cpc_chk(1,ithr), 1)
                            emmg_chk = dot_product(c00E(:,idp_es), wrd_es(:,ithr))
                            myvg_chk = real(Reb_es(0,ithr), dp)
                            do q = 1, ncomp
                                bpc_chk(q,ithr) = real(Reb_es(q,ithr), dp)
                            end do
                            if( l_pol_hyb ) call polar_hybrid_exact_accum(mean_rec, basis_recs, &
                                &ncomp, orientations(i), fpls(i), hex_es, kex_es, npos_es, &
                                &Gpc_chk(:,:,ithr), bpc_chk(:,ithr), cpc_chk(:,ithr), emmg_chk, &
                                &myvg_chk)
                            ! GPU-vs-CPU-polar statistic errors (the gate-1 numbers)
                            gnum_chk = 0.d0; gden_chk = 0.d0; bnum_chk = 0.d0; bden_chk = 0.d0
                            do r = 1, ncomp
                                bnum_chk = bnum_chk + (bth(r,ithr) - bpc_chk(r,ithr))**2
                                bden_chk = bden_chk + bpc_chk(r,ithr)**2
                                do q = 1, ncomp
                                    gnum_chk = gnum_chk + (Gth(q,r,ithr) - Gpc_chk(q,r,ithr))**2
                                    gden_chk = gden_chk + Gpc_chk(q,r,ithr)**2
                                end do
                            end do
                            gerrg_chk(ichk) = sqrt(gnum_chk / max(gden_chk, DTINY))
                            berrg_chk(ichk) = sqrt(bnum_chk / max(bden_chk, DTINY))
                            ! CPU-polar z through the SAME solver settings as the production
                            ! solve below, so the correlation isolates the statistics
                            ag_chk = max(0.1d0, min(5.0d0, myvg_chk / max(emmg_chk, DTINY)))
                            call probe_solve_ecm(ncomp, Gpc_chk(:,:,ithr), bpc_chk(:,ithr), &
                                &cpc_chk(:,ithr), myvg_chk, emmg_chk, prior, sig2, n_probe_cm, &
                                &ag_chk, zpc_chk(:,ithr), Ai_chk(:,:,ithr), lda_chk, lok_chk, &
                                &qml_chk)
                            zchk_g(:,ichk) = zpc_chk(:,ithr)
                            lchk(ichk)     = .true.
                        endif
                        if( l_mix_used )then
                            ! same shared MCFA solver as the CPU bodies: the device fetched the
                            ! identical G/b/c sufficient statistics, only the former differs
                            call probe_solve_mix(ncomp, kmix_pr, Gth(:,:,ithr), bth(:,ithr), &
                                &cth(:,ithr), myv, e_mm, n_probe_cm, a, sig2, mix_Ominv, mix_Omxi, &
                                &mix_lpi, mix_xiOx, zth(:,ithr), Ainvth(:,:,ithr), dens(:,:,i), &
                                &ldA, lok, nll_mix_add, mxa_sr(:,ithr), mxa_sm(:,:,ithr), &
                                &mxa_smm(:,:,:,ithr), mxa_sainv(:,:,ithr))
                            if( lok ) nll_thr(ithr) = nll_thr(ithr) + nll_mix_add
                        else
                            call probe_solve_ecm(ncomp, Gth(:,:,ithr), bth(:,ithr), cth(:,ithr), &
                                &myv, e_mm, prior, sig2, n_probe_cm, a, zth(:,ithr), &
                                &Ainvth(:,:,ithr), ldA, lok, qml)
                            if( lok ) nll_thr(ithr) = nll_thr(ithr) + ldA - qml
                        endif
                        aa = a*a
                        z(row,:)    = zth(:,ithr)
                        zbatch(:,i) = zth(:,ithr)
                        ! parity check: the device-path z actually used, after the shared solve
                        if( l_chk_act )then
                            if( row <= nchk_max ) zchk_p(:,row) = zth(:,ithr)
                        endif
                        if( .not. l_mix_used )then
                            do r = 1, ncomp
                                do q = 1, ncomp
                                    dens(q,r,i) = zth(q,ithr)*zth(r,ithr) + Ainvth(q,r,ithr)
                                end do
                            end do
                        endif
                        do q = 1, ncomp
                            gam_thr(q,ithr) = gam_thr(q,ithr) + dens(q,q,i)
                        end do
                        if( l_probe_mls )then
                            zbatch(:,i) = a*zbatch(:,i)
                            dens(:,:,i) = aa*dens(:,:,i)
                        endif
                        nval_thr(ithr) = nval_thr(ithr) + 1
                        valid(i)       = .true.
                        gam_dbg(1,ithr) = gam_dbg(1,ithr) + sum([(Gth(q,q,ithr), q=1,ncomp)])
                        gam_dbg(2,ithr) = gam_dbg(2,ithr) + dot_product(bth(:,ithr), bth(:,ithr))
                        gam_dbg(3,ithr) = gam_dbg(3,ithr) + dot_product(cth(:,ithr), cth(:,ithr))
                        gam_dbg(4,ithr) = gam_dbg(4,ithr) + a
                        ! residual r_i = y - a*(T mu) in place, banded, exactly the CPU polar
                        ! path's (the check arm's Cartesian mean plane serves where present)
                        if( .not. l_chk_p )then
                            call project_fplane_mean_banded(mean_rec, orientations(i), fpls(i), &
                                &mean_fpl(ithr))
                        endif
                        call subtract_mean_banded(fpls(i), mean_fpl(ithr), real(a), nyqr_es)
                    end do
                    !$omp end parallel do
                    sec_gram_thr(1) = sec_gram_thr(1) + (omp_get_wtime() - twp0)
                elseif( l_bank_fwd )then
                ! ---- banked forward: one raw (no-CTF) projection set per direction segment;
                ! per particle, y/T are RESAMPLED into the bank frame (the same unit-tap scheme
                ! the polar former and the device align use), the CTF is applied to the shared
                ! projections as an elementwise multiply, and the residual is formed already
                ! aligned -- the banked M-step then runs at unit in-plane rotation.
                if( .not. allocated(tmpc) ) allocate(tmpc(fpls(1)%frlims(1,1):fpls(1)%frlims(1,2), &
                    &fpls(1)%frlims(2,1):0, nthr))
                nseg_b = 0
                do i = 1, batchsz
                    row = batchlims(1) + i - 1
                    if( nseg_b > 0 )then
                        if( dirof_pb(row) == seg_dir_b(nseg_b) )then
                            seg_cnt_b(nseg_b) = seg_cnt_b(nseg_b) + 1
                            cycle
                        endif
                    endif
                    nseg_b = nseg_b + 1
                    seg_dir_b(nseg_b) = dirof_pb(row)
                    seg_beg_b(nseg_b) = i
                    seg_cnt_b(nseg_b) = 1
                end do
                !$omp parallel do default(shared) schedule(dynamic) proc_bind(close) &
                !$omp& private(iseg,ii,i,ithr,q,r,a,aa,e_mm,myv,row,twp0,twp1,twp2,nyqb,ldA,lok,qml)
                do iseg = 1, nseg_b
                    ithr = omp_get_thread_num() + 1
                    call dirs_pb%get_ori(seg_dir_b(iseg), o_pb_thr(ithr))
                    twp0 = omp_get_wtime()
                    call project_fplanes_mean_basis(mean_rec, basis_recs, o_pb_thr(ithr), &
                        &fpls(seg_beg_b(iseg)), mean_fpl(ithr), basis_fpls(:,ithr), &
                        &apply_ctf_amp=.false.)
                    twp1 = omp_get_wtime()
                    sec_proj_thr(ithr) = sec_proj_thr(ithr) + (twp1 - twp0)
                    do ii = 1, seg_cnt_b(iseg)
                        i = seg_beg_b(iseg) + ii - 1
                        if( orientations(i)%isstatezero() ) cycle
                        row  = batchlims(1) + i - 1
                        twp1 = omp_get_wtime()
                        nyqb = max(1, fpls(i)%nyq / OSMPL_PAD_FAC)
                        call align_halfplane_inplane(fpls(i)%frlims, nyqb, fpls(i)%cmplx_plane, &
                            &cap(row), sap(row), tmpc(:,:,ithr))
                        fpls(i)%cmplx_plane = tmpc(:,:,ithr)
                        call align_halfplane_inplane(fpls(i)%frlims, nyqb, fpls(i)%transfer_plane, &
                            &cap(row), sap(row), tmpc(:,:,ithr))
                        fpls(i)%transfer_plane = tmpc(:,:,ithr)
                        fpls(i)%ctfsq_plane = real(fpls(i)%transfer_plane*conjg(fpls(i)%transfer_plane))
                        call ensure_bank_ctf_planes(fpls(i), mean_fplC(ithr), basis_fplsC(:,ithr), ncomp)
                        mean_fplC(ithr)%cmplx_plane = fpls(i)%transfer_plane * mean_fpl(ithr)%cmplx_plane
                        do q = 1, ncomp
                            basis_fplsC(q,ithr)%cmplx_plane = fpls(i)%transfer_plane * &
                                &basis_fpls(q,ithr)%cmplx_plane
                        end do
                        e_mm = real(cov_herm_inner(mean_fplC(ithr), mean_fplC(ithr)), dp)
                        myv  = real(cov_herm_inner(mean_fplC(ithr), fpls(i)), dp)
                        a    = max(0.1d0, min(5.0d0, myv / max(e_mm, DTINY)))
                        aa   = a*a
                        do q = 1, ncomp
                            bth(q,ithr) = real(cov_herm_inner(basis_fplsC(q,ithr), fpls(i)), dp)
                            cth(q,ithr) = real(cov_herm_inner(basis_fplsC(q,ithr), mean_fplC(ithr)), dp)
                            do r = q, ncomp
                                Gth(q,r,ithr) = real(cov_herm_inner(basis_fplsC(q,ithr), &
                                    &basis_fplsC(r,ithr)), dp)
                                Gth(r,q,ithr) = Gth(q,r,ithr)
                            end do
                        end do
                        call probe_solve_ecm(ncomp, Gth(:,:,ithr), bth(:,ithr), cth(:,ithr), &
                            &myv, e_mm, prior, sig2, n_probe_cm, a, zth(:,ithr), &
                            &Ainvth(:,:,ithr), ldA, lok, qml)
                        aa = a*a
                        if( lok ) nll_thr(ithr) = nll_thr(ithr) + ldA - qml
                        z(row,:)          = zth(:,ithr)
                        zbatch(:,i)       = zth(:,ithr)
                        do r = 1, ncomp
                            do q = 1, ncomp
                                dens(q,r,i) = zth(q,ithr)*zth(r,ithr) + Ainvth(q,r,ithr)
                            end do
                        end do
                        do q = 1, ncomp
                            gam_thr(q,ithr) = gam_thr(q,ithr) + dens(q,q,i)
                        end do
                        if( l_probe_mls )then
                            zbatch(:,i) = a*zbatch(:,i)
                            dens(:,:,i) = aa*dens(:,:,i)
                        endif
                        nval_thr(ithr)    = nval_thr(ithr) + 1
                        valid(i)          = .true.
                        ! residual in the BANK frame (transfer/ctfsq already aligned above)
                        fpls(i)%cmplx_plane = fpls(i)%cmplx_plane - real(a)*mean_fplC(ithr)%cmplx_plane
                        twp2 = omp_get_wtime()
                        sec_gram_thr(ithr) = sec_gram_thr(ithr) + (twp2 - twp1)
                    end do
                end do
                !$omp end parallel do
                else
                !$omp parallel do default(shared) schedule(dynamic) proc_bind(close) &
                !$omp& private(i,ithr,q,r,a,aa,e_mm,myv,row,twp0,twp1,twp2,ldA,lok,qml,kk2,lwm,wsm) &
                !$omp& private(idp_es,pw_es,cnt_es,taz_es,l_chk_p,ichk,emmc_chk,myvc_chk,ac_chk) &
                !$omp& private(lda_chk,qml_chk,lok_chk,gnum_chk,gden_chk,bnum_chk,bden_chk) &
                !$omp& private(ir,rmx_chk,tazx_chk,pwx_chk,cntx_chk,emmx_chk,myvx_chk,ax_chk) &
                !$omp& private(dcy_chk,dcm_chk,emmd_chk,myvd_chk,ad_chk)
                do i = 1, batchsz
                    if( orientations(i)%isstatezero() ) cycle
                    ithr = omp_get_thread_num() + 1
                    row  = batchlims(1) + i - 1
                    twp0 = omp_get_wtime()
                    if( l_em_polar )then
                        ! statistics already formed in the polar quadrature; the mean plane is still
                        ! needed, in Cartesian, because the M-step backprojects y - a*T*mu
                        call project_fplane_mean(mean_rec, orientations(i), fpls(i), &
                            &mean_fpl(ithr), apply_ctf_amp=.true.)
                        a  = pcon(row)
                        aa = a*a
                        e_mm = 0.d0; myv = 0.d0   ! unused at nml=0; never pass uninitialised
                        Gth(:,:,ithr) = pGc(:,:,row)
                        bth(:,ithr)   = pbc(:,row)
                        cth(:,ithr)   = pcc(:,row)
                        twp1 = omp_get_wtime()
                        sec_proj_thr(ithr) = sec_proj_thr(ithr) + (twp1 - twp0)
                    elseif( l_pol_es )then
                        ! ---- POLAR SHARED-DIRECTION FORMER (stage 1) ----
                        idp_es  = dir_es(row)
                        l_chk_p = l_chk_act .and. row <= nchk_max
                        if( l_chk_p )then
                            ! parity check: the production Cartesian former on this particle; its
                            ! mean plane doubles as the residual's mean projection below
                            call project_fplanes_mean_basis(mean_rec, basis_recs, orientations(i), &
                                &fpls(i), mean_fpl(ithr), basis_fpls(:,ithr), apply_ctf_amp=.true.)
                            emmc_chk = real(cov_herm_inner(mean_fpl(ithr), mean_fpl(ithr)), dp)
                            myvc_chk = real(cov_herm_inner(mean_fpl(ithr), fpls(i)), dp)
                            do q = 1, ncomp
                                bc_chk(q,ithr) = real(cov_herm_inner(basis_fpls(q,ithr), fpls(i)), dp)
                                cc_chk(q,ithr) = real(cov_herm_inner(basis_fpls(q,ithr), mean_fpl(ithr)), dp)
                                do r = q, ncomp
                                    Gc_chk(q,r,ithr) = real(cov_herm_inner(basis_fpls(q,ithr), &
                                        &basis_fpls(r,ithr)), dp)
                                    Gc_chk(r,q,ithr) = Gc_chk(q,r,ithr)
                                end do
                            end do
                        else
                            ! mean only, in Cartesian: the M-step backprojects y - a*(T mu) and
                            ! there is no polar->volume adjoint. Banded variant: identical
                            ! interpolation, none of project_fplane's per-call full-plane
                            ! zero-fill + ctfsq/transfer copies (measured as the bulk of the
                            ! polar project bucket at the native padded plane).
                            call project_fplane_mean_banded(mean_rec, orientations(i), fpls(i), &
                                &mean_fpl(ithr))
                        endif
                        ! polar-sample the prepped data plane ONCE at (bank direction, relative
                        ! in-plane angle): CTF amplitude, shift phase and per-shell whitening all
                        ! ride in from the same cmplx/transfer planes the Cartesian former reads.
                        ! Fused sampler: one KB geometry per ring sample shared by both plane
                        ! gathers, packed output written in place, no per-call allocations --
                        ! bit-identical statistics (see polar_sample_particle_fused).
                        call polar_sample_particle_fused(fpls(i)%cmplx_plane, fpls(i)%transfer_plane, &
                            &pg_es, cae(row), sae(row), xws_es(:,ithr), wr_es(:,ithr), taz_es)
                        wrd_es(:,ithr) = real(wr_es(:,ithr), dp)
                        ! b and the mean row, exact per-sample CTF: one GEMV against the bank
                        call sgemv('T', nsamp2_es, ncomp+1, 1.0, UsallE(1,0,idp_es), nsamp2_es, &
                            &xws_es(1,ithr), 1, 0.0, Reb_es(0,ithr), 1)
                        ! G, c, e_mm by the radial factorisation: ring Grams x per-ring mean |T|^2
                        call dgemv('N', ncomp*ncomp, nk_es, 1.d0, CfE(1,1,idp_es), ncomp*ncomp, &
                            &wrd_es(1,ithr), 1, 0.d0, Gth(1,1,ithr), 1)
                        call dgemv('N', ncomp, nk_es, 1.d0, Cm0E(1,1,idp_es), ncomp, &
                            &wrd_es(1,ithr), 1, 0.d0, cth(1,ithr), 1)
                        e_mm = dot_product(c00E(:,idp_es), wrd_es(:,ithr))
                        myv  = real(Reb_es(0,ithr), dp)
                        do q = 1, ncomp
                            bth(q,ithr) = real(Reb_es(q,ithr), dp)
                        end do
                        ! hybrid: the low-k shells enter as exact Cartesian statistics
                        if( l_pol_hyb ) call polar_hybrid_exact_accum(mean_rec, basis_recs, &
                            &ncomp, orientations(i), fpls(i), hex_es, kex_es, npos_es, &
                            &Gth(:,:,ithr), bth(:,ithr), cth(:,ithr), e_mm, myv)
                        a    = max(0.1d0, min(5.0d0, myv / max(e_mm, DTINY)))
                        aa   = a*a
                        twp1 = omp_get_wtime()
                        sec_proj_thr(ithr) = sec_proj_thr(ithr) + (twp1 - twp0)
                        if( l_chk_p )then
                            ! Cartesian z at this particle (plain fixed-a solve; iteration 1 only,
                            ! where the mixture is never active, so both paths solve the same system)
                            ac_chk = max(0.1d0, min(5.0d0, myvc_chk / max(emmc_chk, DTINY)))
                            call probe_solve_ecm(ncomp, Gc_chk(:,:,ithr), bc_chk(:,ithr), &
                                &cc_chk(:,ithr), myvc_chk, emmc_chk, prior, sig2, 0, ac_chk, &
                                &zc_chk(:,ithr), Ai_chk(:,:,ithr), lda_chk, lok_chk, qml_chk)
                            ichk = row
                            zchk_c(:,ichk) = zc_chk(:,ithr)
                            gnum_chk = 0.d0; gden_chk = 0.d0; bnum_chk = 0.d0; bden_chk = 0.d0
                            do r = 1, ncomp
                                bnum_chk = bnum_chk + (bth(r,ithr) - bc_chk(r,ithr))**2
                                bden_chk = bden_chk + bc_chk(r,ithr)**2
                                do q = 1, ncomp
                                    gnum_chk = gnum_chk + (Gth(q,r,ithr) - Gc_chk(q,r,ithr))**2
                                    gden_chk = gden_chk + Gc_chk(q,r,ithr)**2
                                end do
                            end do
                            gerr_chk(ichk) = sqrt(gnum_chk / max(gden_chk, DTINY))
                            berr_chk(ichk) = sqrt(bnum_chk / max(bden_chk, DTINY))
                            lchk(ichk)     = .true.
                            ! exact-direction polar arm: the same quadrature at the particle's OWN
                            ! direction (identity in-plane), separating the direction snap from the
                            ! polar formulation itself
                            rmx_chk = orientations(i)%get_mat()
                            call polar_project_recs(mean_rec, basis_recs, ncomp, rmx_chk, pg_es, &
                                &UbankE(:,:,ithr))
                            do q = 0, ncomp
                                do r = 1, nsamp_es
                                    UsX_chk(2*r-1,q,ithr) = pg_es%sqwq(r)*real (UbankE(r,q,ithr))
                                    UsX_chk(2*r,  q,ithr) = pg_es%sqwq(r)*aimag(UbankE(r,q,ithr))
                                end do
                            end do
                            do ir = 1, nk_es
                                call polar_ring_gram(UsX_chk(1,0,ithr), nsamp2_es, ncomp, &
                                    &pg_es%rbeg(ir), pg_es%rend(ir)-pg_es%rbeg(ir)+1, &
                                    &CspE(0,0,ithr), CfX_chk(1,ir,ithr), Cm0X_chk(1,ir,ithr))
                                c00X_chk(ir,ithr) = polar_ring_selfpower(UsX_chk(1,0,ithr), &
                                    &nsamp2_es, pg_es%rbeg(ir), pg_es%rend(ir)-pg_es%rbeg(ir)+1)
                            end do
                            call polar_sample_particle_fused(fpls(i)%cmplx_plane, &
                                &fpls(i)%transfer_plane, pg_es, 1.0, 0.0, xwsX_chk(:,ithr), &
                                &wrX_chk(:,ithr), tazx_chk)
                            wrdX_chk(:,ithr) = real(wrX_chk(:,ithr), dp)
                            call sgemv('T', nsamp2_es, ncomp+1, 1.0, UsX_chk(1,0,ithr), nsamp2_es, &
                                &xwsX_chk(1,ithr), 1, 0.0, RebX_chk(0,ithr), 1)
                            call dgemv('N', ncomp*ncomp, nk_es, 1.d0, CfX_chk(1,1,ithr), ncomp*ncomp, &
                                &wrdX_chk(1,ithr), 1, 0.d0, GX_chk(1,1,ithr), 1)
                            call dgemv('N', ncomp, nk_es, 1.d0, Cm0X_chk(1,1,ithr), ncomp, &
                                &wrdX_chk(1,ithr), 1, 0.d0, cX_chk(1,ithr), 1)
                            emmx_chk = dot_product(c00X_chk(:,ithr), wrdX_chk(:,ithr))
                            myvx_chk = real(RebX_chk(0,ithr), dp)
                            do q = 1, ncomp
                                bX_chk(q,ithr) = real(RebX_chk(q,ithr), dp)
                            end do
                            if( l_pol_hyb ) call polar_hybrid_exact_accum(mean_rec, basis_recs, &
                                &ncomp, orientations(i), fpls(i), hex_es, kex_es, npos_es, &
                                &GX_chk(:,:,ithr), bX_chk(:,ithr), cX_chk(:,ithr), emmx_chk, myvx_chk)
                            ax_chk = max(0.1d0, min(5.0d0, myvx_chk / max(emmx_chk, DTINY)))
                            call probe_solve_ecm(ncomp, GX_chk(:,:,ithr), bX_chk(:,ithr), &
                                &cX_chk(:,ithr), myvx_chk, emmx_chk, prior, sig2, 0, ax_chk, &
                                &zx_chk(:,ithr), Ai_chk(:,:,ithr), lda_chk, lok_chk, qml_chk)
                            zchk_x(:,ichk) = zx_chk(:,ithr)
                            gnum_chk = 0.d0; gden_chk = 0.d0; bnum_chk = 0.d0; bden_chk = 0.d0
                            do r = 1, ncomp
                                bnum_chk = bnum_chk + (bX_chk(r,ithr) - bc_chk(r,ithr))**2
                                bden_chk = bden_chk + bc_chk(r,ithr)**2
                                do q = 1, ncomp
                                    gnum_chk = gnum_chk + (GX_chk(q,r,ithr) - Gc_chk(q,r,ithr))**2
                                    gden_chk = gden_chk + Gc_chk(q,r,ithr)**2
                                end do
                            end do
                            gerrx_chk(ichk)  = sqrt(gnum_chk / max(gden_chk, DTINY))
                            berrx_chk(ichk)  = sqrt(bnum_chk / max(bden_chk, DTINY))
                            tazsum_thr(ithr) = tazsum_thr(ithr) + real(tazx_chk, dp)
                            tazcnt_thr(ithr) = tazcnt_thr(ithr) + 1
                            ! DC-less Cartesian arm: subtract the single (0,0) sample from every
                            ! Cartesian statistic (the polar quadrature has no r=0 sample). Exact-
                            ! polar agreement with THIS reference convicts or clears the DC term.
                            dcy_chk = fpls(i)%cmplx_plane(0,0)
                            dcm_chk = mean_fpl(ithr)%cmplx_plane(0,0)
                            do q = 1, ncomp
                                dcb_chk(q,ithr) = basis_fpls(q,ithr)%cmplx_plane(0,0)
                            end do
                            emmd_chk = emmc_chk - real(conjg(dcm_chk)*dcm_chk, dp)
                            myvd_chk = myvc_chk - real(conjg(dcm_chk)*dcy_chk, dp)
                            do q = 1, ncomp
                                bd_chk(q,ithr) = bc_chk(q,ithr) - real(conjg(dcb_chk(q,ithr))*dcy_chk, dp)
                                cd_chk(q,ithr) = cc_chk(q,ithr) - real(conjg(dcb_chk(q,ithr))*dcm_chk, dp)
                                do r = q, ncomp
                                    Gd_chk(q,r,ithr) = Gc_chk(q,r,ithr) &
                                        &- real(conjg(dcb_chk(q,ithr))*dcb_chk(r,ithr), dp)
                                    Gd_chk(r,q,ithr) = Gd_chk(q,r,ithr)
                                end do
                            end do
                            ad_chk = max(0.1d0, min(5.0d0, myvd_chk / max(emmd_chk, DTINY)))
                            call probe_solve_ecm(ncomp, Gd_chk(:,:,ithr), bd_chk(:,ithr), &
                                &cd_chk(:,ithr), myvd_chk, emmd_chk, prior, sig2, 0, ad_chk, &
                                &zd_chk(:,ithr), Ai_chk(:,:,ithr), lda_chk, lok_chk, qml_chk)
                            zchk_d(:,ichk) = zd_chk(:,ithr)
                            ! exact-direction polar vs DC-less Cartesian statistic errors
                            gnum_chk = 0.d0; gden_chk = 0.d0; bnum_chk = 0.d0; bden_chk = 0.d0
                            do r = 1, ncomp
                                bnum_chk = bnum_chk + (bX_chk(r,ithr) - bd_chk(r,ithr))**2
                                bden_chk = bden_chk + bd_chk(r,ithr)**2
                                do q = 1, ncomp
                                    gnum_chk = gnum_chk + (GX_chk(q,r,ithr) - Gd_chk(q,r,ithr))**2
                                    gden_chk = gden_chk + Gd_chk(q,r,ithr)**2
                                end do
                            end do
                            gerrd_chk(ichk) = sqrt(gnum_chk / max(gden_chk, DTINY))
                            berrd_chk(ichk) = sqrt(bnum_chk / max(bden_chk, DTINY))
                        endif
                    else
                    call project_fplanes_mean_basis(mean_rec, basis_recs, orientations(i), fpls(i), &
                        &mean_fpl(ithr), basis_fpls(:,ithr), apply_ctf_amp=.true.)
                    twp1 = omp_get_wtime()
                    sec_proj_thr(ithr) = sec_proj_thr(ithr) + (twp1 - twp0)
                    e_mm = real(cov_herm_inner(mean_fpl(ithr), mean_fpl(ithr)), dp)
                    myv  = real(cov_herm_inner(mean_fpl(ithr), fpls(i)), dp)
                    a    = max(0.1d0, min(5.0d0, myv / max(e_mm, DTINY)))
                    aa   = a*a
                    do q = 1, ncomp
                        bth(q,ithr) = real(cov_herm_inner(basis_fpls(q,ithr), fpls(i)), dp)
                        cth(q,ithr) = real(cov_herm_inner(basis_fpls(q,ithr), mean_fpl(ithr)), dp)
                        do r = q, ncomp
                            Gth(q,r,ithr) = real(cov_herm_inner(basis_fpls(q,ithr), basis_fpls(r,ithr)), dp)
                            Gth(r,q,ithr) = Gth(q,r,ithr)
                        end do
                    end do
                    endif
                    ! Posterior precision A = (a^2/sig2) G + Gamma^-1. The whole normal system is already
                    ! scaled by 1/sig2, so Cov[z|y] = A^-1 exactly -- no further sig2 factor.
                    if( l_mix_used )then
                        ! ---- MCFA E-step via the ONE shared solver (probe_solve_mix) ----
                        call probe_solve_mix(ncomp, kmix_pr, Gth(:,:,ithr), bth(:,ithr), &
                            &cth(:,ithr), myv, e_mm, nml_plain, a, sig2, mix_Ominv, mix_Omxi, &
                            &mix_lpi, mix_xiOx, zth(:,ithr), Ainvth(:,:,ithr), dens(:,:,i), &
                            &ldA, lok, nll_mix_add, mxa_sr(:,ithr), mxa_sm(:,:,ithr), &
                            &mxa_smm(:,:,:,ithr), mxa_sainv(:,:,ithr))
                        if( lok ) nll_thr(ithr) = nll_thr(ithr) + nll_mix_add
                    else
                        call probe_solve_ecm(ncomp, Gth(:,:,ithr), bth(:,ithr), cth(:,ithr), &
                            &myv, e_mm, prior, sig2, nml_plain, a, zth(:,ithr), &
                            &Ainvth(:,:,ithr), ldA, lok, qml)
                        if( lok ) nll_thr(ithr) = nll_thr(ithr) + ldA - qml
                    endif
                    aa = a*a
                    z(row,:)          = zth(:,ithr)
                    zbatch(:,i)       = zth(:,ithr)
                    ! parity check: the polar-path z actually used, captured after the shared solve
                    if( l_chk_act )then
                        if( row <= nchk_max ) zchk_p(:,row) = zth(:,ithr)
                    endif
                    ! EM sufficient statistic E[z z'|y] = z z' + Cov[z|y]. BOTH the coupled M-step normal
                    ! matrix and the Gamma update below need it. Dropping Cov underestimates Gamma, which
                    ! tightens the prior, which shrinks z further: the bias compounds across iterations.
                    ! Under the mixture E[zz'|y] = A^-1 + sum_k r_k m_k m_k', NOT A^-1 + E[z]E[z]'
                    ! -- the between-component spread is real posterior variance; probe_solve_mix
                    ! has already written dens(:,:,i) in that case.
                    if( .not. l_mix_used )then
                        do r = 1, ncomp
                            do q = 1, ncomp
                                dens(q,r,i) = zth(q,ithr)*zth(r,ithr) + Ainvth(q,r,ithr)
                            end do
                        end do
                    endif
                    do q = 1, ncomp
                        gam_thr(q,ithr) = gam_thr(q,ithr) + dens(q,q,i)
                    end do
                    if( l_probe_mls )then
                        zbatch(:,i) = a*zbatch(:,i)
                        dens(:,:,i) = aa*dens(:,:,i)
                    endif
                    nval_thr(ithr)    = nval_thr(ithr) + 1
                    valid(i)          = .true.
                    gam_dbg(1,ithr) = gam_dbg(1,ithr) + sum([(Gth(q,q,ithr), q=1,ncomp)])
                    gam_dbg(2,ithr) = gam_dbg(2,ithr) + dot_product(bth(:,ithr), bth(:,ithr))
                    gam_dbg(3,ithr) = gam_dbg(3,ithr) + dot_product(cth(:,ithr), cth(:,ithr))
                    gam_dbg(4,ithr) = gam_dbg(4,ithr) + a
                    ! residual observation r_i = y - a*(T mu) in place (transfer/ctfsq intact for backprojection)
                    if( l_pol_es )then
                        ! banded: the mean plane is zero outside the working disc (both formers
                        ! write the same disc), so the full-array statement only rewrote
                        ! unchanged values there -- at the native padded lattice that traffic
                        ! was most of the polar path's gram+solve bucket
                        call subtract_mean_banded(fpls(i), mean_fpl(ithr), real(a), nyqr_es)
                    else
                        fpls(i)%cmplx_plane = fpls(i)%cmplx_plane - real(a)*mean_fpl(ithr)%cmplx_plane
                    endif
                    twp2 = omp_get_wtime()
                    sec_gram_thr(ithr) = sec_gram_thr(ithr) + (twp2 - twp1)
                end do
                !$omp end parallel do
                endif
                sec_estep = sec_estep + toc(t_sec)
                t_sec = tic()
                ! M-step by halfset: Y_q += sum_i z_iq * backproject(r_i), and the coupled normal matrix
                ! rho(q,r) += sum_i |CTF|^2 E[z_iq z_ir]   (batched KB)
                do i = 1, batchsz
                    valid_e(i) = valid(i) .and. eo(i)==0
                    valid_o(i) = valid(i) .and. eo(i)==1
                end do
                if( l_gpu_probe )then
                    ! halfset routing through the slots: even records fill 1..ncomp / 1..npairs,
                    ! odd records the upper halves; density slots honor the diag/full mode
                    gdsc(:,:batchsz) = 0.d0
                    grsc(:,:batchsz) = 0.d0
                    do i = 1, batchsz
                        if( .not. (valid_e(i) .or. valid_o(i)) ) cycle
                        gq = merge(0, ncomp,  valid_e(i))
                        gr = merge(0, npairs, valid_e(i))
                        gdsc(gq+1:gq+ncomp,i) = zbatch(:,i)
                        do r = 1, ncomp
                            do q = 1, r
                                grsc(gr+(r*(r-1))/2+q,i) = dens(q,r,i)
                            end do
                        end do
                    end do
                    if( l_bank_adj )then
                        if( l_bank_fwd )then
                            ! residuals were built in the bank frame by the forward -- unit rotation
                            call flex_gpu_coupled_batch_banked_f(fpls(:batchsz), gdsc(:,:batchsz), &
                                &grsc(:,:batchsz), valid_e(:batchsz) .or. valid_o(:batchsz), batchsz, &
                                &dirof_pb(batchlims(1):batchlims(2)), cap1(:batchsz), sap0(:batchsz))
                        else if( l_gpu_estep )then
                            ! fused path: the residual is already resident on device
                            call flex_gpu_coupled_batch_banked_res_f(gdsc(:,:batchsz), &
                                &grsc(:,:batchsz), valid_e(:batchsz) .or. valid_o(:batchsz), batchsz, &
                                &dirof_pb(batchlims(1):batchlims(2)), cap(batchlims(1):batchlims(2)), &
                                &sap(batchlims(1):batchlims(2)))
                        else
                            call flex_gpu_coupled_batch_banked_f(fpls(:batchsz), gdsc(:,:batchsz), &
                                &grsc(:,:batchsz), valid_e(:batchsz) .or. valid_o(:batchsz), batchsz, &
                                &dirof_pb(batchlims(1):batchlims(2)), cap(batchlims(1):batchlims(2)), &
                                &sap(batchlims(1):batchlims(2)))
                        endif
                    else
                        call flex_gpu_coupled_batch_raw_f(build%pgrpsyms, orientations(:batchsz), &
                            &fpls(:batchsz), gdsc(:,:batchsz), grsc(:,:batchsz), &
                            &valid_e(:batchsz) .or. valid_o(:batchsz), batchsz)
                    endif
                else
                    call insert_planes_oversamp_coupled_batch_scaled(Yeven, rho_e, build%pgrpsyms, &
                        &orientations(:batchsz), fpls(:batchsz), zbatch(:,:batchsz), dens(:,:,:batchsz), &
                        &valid_e(:batchsz), batchsz)
                    call insert_planes_oversamp_coupled_batch_scaled(Yodd, rho_o, build%pgrpsyms, &
                        &orientations(:batchsz), fpls(:batchsz), zbatch(:,:batchsz), dens(:,:,:batchsz), &
                        &valid_o(:batchsz), batchsz)
                endif
                sec_ins = sec_ins + toc(t_sec)
                if( batchlims(2)==nptcls .or. mod(batchlims(2), 5*MAXIMGBATCHSZ)==0 )then
                    write(logfhandle,'(A,I0,A,I0)') '>>> FLEX_PCA PROBE PASS PARTICLES: ',batchlims(2),' / ',npp
                    call flush(logfhandle)
                endif
            end do
            if( l_gpu_probe )then
                call flex_gpu_coupled_end_f(Yeven, rho_e, Yodd, rho_o)
                deallocate(gdsc, grsc)
            endif
            write(logfhandle,'(A,F7.1,A,F7.1,A,F7.1,A,F7.1)') '>>> FLEX_PCA PROBE E-STEP SPLIT (seconds): read=', &
                &sec_read,'  prep=',sec_prep,'  project+solve=',sec_estep,'  insert=',sec_ins
            write(logfhandle,'(A,F7.1,A,F7.1)') '>>> FLEX_PCA PROBE E-STEP INNER (thread-seconds): project=', &
                &sum(sec_proj_thr),'  gram+solve=',sum(sec_gram_thr)
            call flush(logfhandle)
            if( l_pol_es )then
                ! the two numbers the polar A/B reads from run.log
                write(logfhandle,'(A,I0,A,F7.1,A,F7.1)') '>>> FLEX_PCA POLAR ESTEP it=',it, &
                    &'  bank build seconds=',sec_bank,'  estep seconds=',sec_estep
                call flush(logfhandle)
            endif
            ! ---- POLAR CHECK SUMMARY (iteration 1): three formers ran on the first N particles:
            ! Cartesian (reference), polar at the bank direction (the production stage-1 path) and
            ! polar at the exact direction. banked-vs-exact isolates the direction quantization;
            ! exact-vs-Cartesian isolates the polar formulation (quadrature + radial |T|^2). ----
            if( l_chk_act )then
                block
                    integer  :: nck, qq2, j2
                    real(dp), allocatable :: va_ck(:), vb_ck(:), vx_ck(:)
                    real,     allocatable :: se_ck(:)
                    real(dp) :: cq_ck(ncomp), cqx_ck(ncomp), cqbx_ck(ncomp)
                    real(dp) :: gmed_ck, bmed_ck, gmedx_ck, bmedx_ck, taz_ck
                    nck = count(lchk)
                    if( nck >= 3 )then
                        allocate(va_ck(nck), vb_ck(nck), vx_ck(nck), se_ck(nck))
                        do qq2 = 1, ncomp
                            j2 = 0
                            do ichk = 1, nchk_max
                                if( .not. lchk(ichk) ) cycle
                                j2 = j2 + 1
                                va_ck(j2) = zchk_c(qq2,ichk)
                                vb_ck(j2) = zchk_p(qq2,ichk)
                                vx_ck(j2) = 0.d0
                                if( allocated(zchk_x) ) vx_ck(j2) = zchk_x(qq2,ichk)
                            end do
                            cq_ck(qq2)   = corr_dp(va_ck, vb_ck, nck)
                            cqx_ck(qq2)  = 0.d0
                            cqbx_ck(qq2) = 0.d0
                            if( allocated(zchk_x) )then
                                cqx_ck(qq2)  = corr_dp(va_ck, vx_ck, nck)
                                cqbx_ck(qq2) = corr_dp(vb_ck, vx_ck, nck)
                            endif
                        end do
                        j2 = 0
                        do ichk = 1, nchk_max
                            if( .not. lchk(ichk) ) cycle
                            j2 = j2 + 1
                            se_ck(j2) = real(gerr_chk(ichk))
                        end do
                        call hpsort(se_ck)
                        gmed_ck = real(se_ck((nck+1)/2), dp)
                        j2 = 0
                        do ichk = 1, nchk_max
                            if( .not. lchk(ichk) ) cycle
                            j2 = j2 + 1
                            se_ck(j2) = real(berr_chk(ichk))
                        end do
                        call hpsort(se_ck)
                        bmed_ck = real(se_ck((nck+1)/2), dp)
                        gmedx_ck = 0.d0
                        bmedx_ck = 0.d0
                        if( allocated(gerrx_chk) )then
                            j2 = 0
                            do ichk = 1, nchk_max
                                if( .not. lchk(ichk) ) cycle
                                j2 = j2 + 1
                                se_ck(j2) = real(gerrx_chk(ichk))
                            end do
                            call hpsort(se_ck)
                            gmedx_ck = real(se_ck((nck+1)/2), dp)
                            j2 = 0
                            do ichk = 1, nchk_max
                                if( .not. lchk(ichk) ) cycle
                                j2 = j2 + 1
                                se_ck(j2) = real(berrx_chk(ichk))
                            end do
                            call hpsort(se_ck)
                            bmedx_ck = real(se_ck((nck+1)/2), dp)
                        endif
                        write(logfhandle,'(A,I0,A)') '>>> FLEX_PCA POLAR CHECK: Cartesian vs polar &
                            &E-step statistics on ',nck,' particles (iteration 1)'
                        if( allocated(zchk_g) )then
                            write(logfhandle,'(A)') '>>>   z correlation per component (DEVICE polar vs Cartesian):'
                        else
                            write(logfhandle,'(A)') '>>>   z correlation per component (BANKED polar vs Cartesian):'
                        endif
                        do qq2 = 1, ncomp, 8
                            write(logfhandle,'(A,I3,A,8F8.4)') '>>>     z',qq2,'..', &
                                &cq_ck(qq2:min(ncomp,qq2+7))
                        end do
                        write(logfhandle,'(A,F8.4,A,ES10.2,A,ES10.2)') '>>>   BANKED vs Cartesian:   min z corr=', &
                            &minval(cq_ck),'  median rel |G-Gc|_F/|Gc|_F=',gmed_ck, &
                            &'  median rel |b-bc|/|bc|=',bmed_ck
                        ! ---- GPU parity summary (gate 1): device former vs the CPU polar former,
                        ! identical bank/tables/hybrid and identical solver settings ----
                        if( allocated(zchk_g) )then
                            block
                                real(dp) :: cqg_ck(ncomp), gmedg_ck, bmedg_ck
                                do qq2 = 1, ncomp
                                    j2 = 0
                                    do ichk = 1, nchk_max
                                        if( .not. lchk(ichk) ) cycle
                                        j2 = j2 + 1
                                        va_ck(j2) = zchk_g(qq2,ichk)
                                        vb_ck(j2) = zchk_p(qq2,ichk)
                                    end do
                                    cqg_ck(qq2) = corr_dp(va_ck, vb_ck, nck)
                                end do
                                j2 = 0
                                do ichk = 1, nchk_max
                                    if( .not. lchk(ichk) ) cycle
                                    j2 = j2 + 1
                                    se_ck(j2) = real(gerrg_chk(ichk))
                                end do
                                call hpsort(se_ck)
                                gmedg_ck = real(se_ck((nck+1)/2), dp)
                                j2 = 0
                                do ichk = 1, nchk_max
                                    if( .not. lchk(ichk) ) cycle
                                    j2 = j2 + 1
                                    se_ck(j2) = real(berrg_chk(ichk))
                                end do
                                call hpsort(se_ck)
                                bmedg_ck = real(se_ck((nck+1)/2), dp)
                                write(logfhandle,'(A)') '>>>   z correlation per component (DEVICE &
                                    &polar vs CPU polar):'
                                do qq2 = 1, ncomp, 8
                                    write(logfhandle,'(A,I3,A,8F8.4)') '>>>     z',qq2,'..', &
                                        &cqg_ck(qq2:min(ncomp,qq2+7))
                                end do
                                write(logfhandle,'(A,F8.4,A,F8.4,A,ES10.2,A,ES10.2)') &
                                    &'>>>   DEVICE vs CPU-POLAR:   min z corr=',minval(cqg_ck), &
                                    &'  mean=',sum(cqg_ck)/real(ncomp,dp), &
                                    &'  median rel |G-Gp|_F/|Gp|_F=',gmedg_ck, &
                                    &'  median rel |b-bp|/|bp|=',bmedg_ck
                            end block
                        endif
                        if( allocated(zchk_x) )then
                        write(logfhandle,'(A)') '>>>   z correlation per component (EXACT-DIR polar vs Cartesian):'
                        do qq2 = 1, ncomp, 8
                            write(logfhandle,'(A,I3,A,8F8.4)') '>>>     z',qq2,'..', &
                                &cqx_ck(qq2:min(ncomp,qq2+7))
                        end do
                        write(logfhandle,'(A,F8.4,A,ES10.2,A,ES10.2)') '>>>   EXACT-DIR vs Cartesian: min z corr=', &
                            &minval(cqx_ck),'  median rel |G-Gc|_F/|Gc|_F=',gmedx_ck, &
                            &'  median rel |b-bc|/|bc|=',bmedx_ck
                        write(logfhandle,'(A,F8.4,A,F8.4)') '>>>   BANKED vs EXACT-DIR (direction &
                            &quantization alone): min z corr=',minval(cqbx_ck),'  mean=', &
                            &sum(cqbx_ck)/real(ncomp,dp)
                        taz_ck = sum(tazsum_thr) / real(max(1, sum(tazcnt_thr)), dp)
                        write(logfhandle,'(A,F8.4,A)') '>>>   mean within-ring |T|^2 spread (tazim)=', &
                            &taz_ck,'  (bounds the radial-factorisation error in G)'
                        endif
                        ! DC attribution: exact-direction polar against the DC-less Cartesian
                        ! (CPU attribution arm only; stands down under the device path)
                        if( allocated(zchk_x) )then
                        block
                            real(dp) :: cqd_ck(ncomp), cqcd_ck(ncomp), gdm_ck, bdm_ck
                            do qq2 = 1, ncomp
                                j2 = 0
                                do ichk = 1, nchk_max
                                    if( .not. lchk(ichk) ) cycle
                                    j2 = j2 + 1
                                    va_ck(j2) = zchk_d(qq2,ichk)
                                    vx_ck(j2) = zchk_x(qq2,ichk)
                                    vb_ck(j2) = zchk_c(qq2,ichk)
                                end do
                                cqd_ck(qq2)  = corr_dp(va_ck, vx_ck, nck)   ! exact polar vs DC-less cart
                                cqcd_ck(qq2) = corr_dp(va_ck, vb_ck, nck)   ! DC-less cart vs full cart
                            end do
                            j2 = 0
                            do ichk = 1, nchk_max
                                if( .not. lchk(ichk) ) cycle
                                j2 = j2 + 1
                                se_ck(j2) = real(gerrd_chk(ichk))
                            end do
                            call hpsort(se_ck)
                            gdm_ck = real(se_ck((nck+1)/2), dp)
                            j2 = 0
                            do ichk = 1, nchk_max
                                if( .not. lchk(ichk) ) cycle
                                j2 = j2 + 1
                                se_ck(j2) = real(berrd_chk(ichk))
                            end do
                            call hpsort(se_ck)
                            bdm_ck = real(se_ck((nck+1)/2), dp)
                            write(logfhandle,'(A)') '>>>   z correlation per component (EXACT-DIR polar &
                                &vs DC-LESS Cartesian):'
                            do qq2 = 1, ncomp, 8
                                write(logfhandle,'(A,I3,A,8F8.4)') '>>>     z',qq2,'..', &
                                    &cqd_ck(qq2:min(ncomp,qq2+7))
                            end do
                            write(logfhandle,'(A,F8.4,A,ES10.2,A,ES10.2)') &
                                &'>>>   EXACT-DIR vs DC-LESS Cartesian: min z corr=',minval(cqd_ck), &
                                &'  median rel G err=',gdm_ck,'  median rel b err=',bdm_ck
                            write(logfhandle,'(A,F8.4,A,F8.4)') '>>>   DC-LESS vs FULL Cartesian &
                                &(the DC term alone): min z corr=',minval(cqcd_ck),'  mean=', &
                                &sum(cqcd_ck)/real(ncomp,dp)
                        end block
                        endif
                        call flush(logfhandle)
                        deallocate(va_ck, vb_ck, vx_ck, se_ck)
                    else
                        write(logfhandle,'(A,I0,A)') '>>> FLEX_PCA POLAR CHECK: only ',nck, &
                            &' particles checked -- too few for a summary'
                        call flush(logfhandle)
                    endif
                end block
                deallocate(zchk_c, zchk_p, gerr_chk, berr_chk, lchk, Gc_chk, bc_chk, cc_chk, &
                    &zc_chk, Ai_chk)
                deallocate(zchk_d, gerrd_chk, berrd_chk, Gd_chk, bd_chk, cd_chk, zd_chk, dcb_chk)
                if( allocated(zchk_x) ) deallocate(zchk_x, gerrx_chk, berrx_chk, UsX_chk, &
                    &xwsX_chk, wrX_chk, RebX_chk, CfX_chk, Cm0X_chk, c00X_chk, wrdX_chk, &
                    &GX_chk, bX_chk, cX_chk, zx_chk, tazsum_thr, tazcnt_thr)
                if( allocated(zchk_g) ) deallocate(zchk_g, gerrg_chk, berrg_chk, Gpc_chk, &
                    &bpc_chk, cpc_chk, zpc_chk)
                l_chk_act   = .false.
                n_pol_check = 0
            endif
            ! reduce the EM Gamma accumulator BEFORE ncomp is replaced by d_new below.
            ! Gamma travels between parts as a SUM and is divided by the REDUCED nval: dividing per
            ! part would weight a small part equally with a large one.
            nval = sum(nval_thr)
            ! E-STEP STATISTICS SUMMARY. The polar and Cartesian formers are supposed to produce
            ! the SAME G, b, c on the same basis; at iteration 1 the basis is the deterministic
            ! data-free init, so any difference here is the former, not the fit. Printing the
            ! actual moments is the only way to tell a scale mismatch from a genuine difference.
            if( nval > 0 )then
                write(logfhandle,'(A,I0,A,ES12.4,A,ES12.4,A,ES12.4,A,F7.4)') &
                    &'>>> FLEX_PCA PROBE ESTAT it=',it,'  <trG>=',sum(gam_dbg(1,:))/real(nval,dp), &
                    &'  <b.b>=',sum(gam_dbg(2,:))/real(nval,dp), &
                    &'  <c.c>=',sum(gam_dbg(3,:))/real(nval,dp), &
                    &'  <a>=',real(sum(gam_dbg(4,:))/real(nval,dp))
                call flush(logfhandle)
            endif
            ! In-process reduction of the likelihood accumulator. This was MISSING: nll_tot was
            ! summed only in the worker branch, so a shared-memory run (nparts=1) reported just the
            ! N*log det Gamma term and none of the per-particle part -- which silently invalidated
            ! every nparts=1 likelihood comparison. Distributed runs were unaffected, since
            ! reduce_probe_parts sums the workers' contributions.
            nll_tot = sum(nll_thr)
            do q = 1, ncomp
                gam_sum(q) = sum(gam_thr(q,:))
            end do
            ! ---- MCFA: mixture M-step, or its (re)initialisation from the current latents ----
            ! nparts>1 is excluded at the flag read, so this always runs on the in-process
            ! master with this iteration's accumulators complete.
            if( l_mix_req .and. .not. flex_pca_is_worker() )then
                if( l_mix_active )then
                    block
                        real(dp), allocatable :: rr_sr(:), rr_sm(:,:), rr_smm(:,:,:), rr_sai(:,:)
                        real(dp) :: g_mx, pi_floor, best_d, dmin, dsq
                        integer  :: tt, kk3, kkp, ii2, best_i
                        allocate(rr_sr(kmix_pr), rr_sm(ncomp,kmix_pr), &
                            &rr_smm(ncomp,ncomp,kmix_pr), rr_sai(ncomp,ncomp))
                        rr_sr = 0.d0; rr_sm = 0.d0; rr_smm = 0.d0; rr_sai = 0.d0
                        do tt = 1, nthr
                            rr_sr  = rr_sr  + mxa_sr(:,tt)
                            rr_sm  = rr_sm  + mxa_sm(:,:,tt)
                            rr_sai = rr_sai + mxa_sainv(:,:,tt)
                            do kk3 = 1, kmix_pr
                                rr_smm(:,:,kk3) = rr_smm(:,:,kk3) + mxa_smm(:,:,kk3,tt)
                            end do
                        end do
                        if( allocated(dm_sr) .and. flex_pca_nparts() > 1 )then
                            ! distributed: the reduce already summed every part's thread-sums
                            call mcfa_mstep(ncomp, kmix_pr, nval, dm_sr, dm_sm, dm_smm, dm_sai, &
                                &kmix_pr == 1, mix_pi, mix_xi, mix_Om, mix_Ominv, ldOm_mix)
                        else
                            call mcfa_mstep(ncomp, kmix_pr, nval, rr_sr, rr_sm, rr_smm, rr_sai, &
                                &kmix_pr == 1, mix_pi, mix_xi, mix_Om, mix_Ominv, ldOm_mix)
                        endif
                        deallocate(rr_sr, rr_sm, rr_smm, rr_sai)
                    end block
                else if( it >= n_mix_warm )then
                    allocate(mix_xi(ncomp,kmix_pr), mix_Om(ncomp,ncomp), mix_Ominv(ncomp,ncomp), &
                        &mix_pi(kmix_pr), mix_Omxi(ncomp,kmix_pr), mix_xiOx(kmix_pr), mix_lpi(kmix_pr))
                    if( flex_pca_nparts() > 1 .and. dm_nz > 0 )then
                        ! distributed: seed from the pooled per-part latent subsample
                        call mcfa_init(dm_z(:dm_nz,:), dm_nz, ncomp, kmix_pr, gam_sum, nval, &
                            &mix_xi, mix_pi, mix_Om)
                    else
                        call mcfa_init(z, size(z,1), ncomp, kmix_pr, gam_sum, nval, mix_xi, mix_pi, mix_Om)
                    endif
                    call mcfa_condition(ncomp, kmix_pr == 1, mix_Om, mix_Ominv, ldOm_mix)
                    l_mix_active = .true.
                    write(logfhandle,'(A,I0,A,I0)') '>>> FLEX_PCA MIX initialised: K=',kmix_pr, &
                        &'  at iteration ',it
                endif
                if( l_mix_active )then
                    do kk2 = 1, kmix_pr
                        mix_Omxi(:,kk2) = matmul(mix_Ominv, mix_xi(:,kk2))
                        mix_xiOx(kk2)   = dot_product(mix_xi(:,kk2), mix_Omxi(:,kk2))
                        mix_lpi(kk2)    = log(max(mix_pi(kk2), 1.d-12))
                    end do
                    write(logfhandle,'(A,I0,A,ES11.3,A,F7.4,A,F8.3)') '>>> FLEX_PCA MIX it=',it, &
                        &'  logdetOm=',ldOm_mix,'  pi_max=',maxval(mix_pi), &
                        &'  |xi|_max=',sqrt(maxval(sum(mix_xi**2, dim=1)))
                    call flush(logfhandle)
                endif
            endif
            if( flex_pca_is_worker() )then
                allocate(cme(es(1),es(2),es(3),ncomp), cmo(es(1),es(2),es(3),ncomp))
                allocate(rhe(es(1),es(2),es(3),ncomp), rhoo(es(1),es(2),es(3),ncomp))
                do q = 1, ncomp
                    cme(:,:,:,q) = Yeven(q)%cmat_exp; rhe(:,:,:,q)  = Yeven(q)%rho_exp
                    cmo(:,:,:,q) = Yodd(q)%cmat_exp;  rhoo(:,:,:,q) = Yodd(q)%rho_exp
                end do
                nll_tot = sum(nll_thr)
                if( l_mix_req )then
                    block
                        real(dp), allocatable :: w_sr(:), w_sm(:,:), w_smm(:,:,:), w_sai(:,:), w_z(:,:)
                        integer :: tt2, kk4, nzs, izs, istep
                        allocate(w_sr(kmix_pr), w_sm(ncomp,kmix_pr), &
                            &w_smm(ncomp,ncomp,kmix_pr), w_sai(ncomp,ncomp))
                        w_sr = 0.d0; w_sm = 0.d0; w_smm = 0.d0; w_sai = 0.d0
                        do tt2 = 1, nthr
                            w_sr  = w_sr  + mxa_sr(:,tt2)
                            w_sm  = w_sm  + mxa_sm(:,:,tt2)
                            w_sai = w_sai + mxa_sainv(:,:,tt2)
                            do kk4 = 1, kmix_pr
                                w_smm(:,:,kk4) = w_smm(:,:,kk4) + mxa_smm(:,:,kk4,tt2)
                            end do
                        end do
                        ! deterministic stride subsample of this part's latents (bounded)
                        nzs   = min(MIX_ZSUB_MAX, npp)
                        istep = max(1, npp / max(1,nzs))
                        nzs   = min(nzs, (npp + istep - 1)/istep)
                        allocate(w_z(nzs,ncomp))
                        do izs = 1, nzs
                            w_z(izs,:) = z(min(npp, 1 + (izs-1)*istep), :)
                        end do
                        call write_probe_part(flex_pca_part_fname('probe', params%part, params%numlen), &
                            &cme, rhe, cmo, rhoo, rho_e, rho_o, gam_sum, nll_tot, nval, ncomp, &
                            &w_sr, w_sm, w_smm, w_sai, w_z)
                        deallocate(w_sr, w_sm, w_smm, w_sai, w_z)
                    end block
                else
                    call write_probe_part(flex_pca_part_fname('probe', params%part, params%numlen), &
                        &cme, rhe, cmo, rhoo, rho_e, rho_o, gam_sum, nll_tot, nval, ncomp)
                endif
                deallocate(cme, cmo, rhe, rhoo, gam_sum)
                call cleanup_rec_buffers(build, fpls)
                do q = 1, size(Yeven)
                    call Yeven(q)%dealloc_rho; call Yeven(q)%kill
                    call Yodd(q)%dealloc_rho;  call Yodd(q)%kill
                end do
                deallocate(Yeven, Yodd, rho_e, rho_o, prior)
                deallocate(Gth, Ath, bth, cth, zth, basis_fpls, mean_fpl, orientations, zbatch, dens)
                deallocate(valid, valid_e, valid_o, eo, Ainvth, Acpth, gam_thr, gam_acc, nval_thr)
                deallocate(hth, nll_thr)
                if( l_bank_fwd )then
                    do ithr = 1, nthr
                        call cleanup_plane(mean_fplC(ithr))
                        do q = 1, size(basis_fplsC,1); call cleanup_plane(basis_fplsC(q,ithr)); end do
                    end do
                    deallocate(mean_fplC, basis_fplsC)
                endif
                if( l_bank_adj )then
                    call flex_gpu_coupled_bank_free_f
                    call flex_gpu_estep_free_f
                    if( l_gpu_prep )then
                        call flex_gpu_prep_free_f
                        deallocate(ctfp_arr, shf_pb)
                        if( allocated(sig2_pb) ) deallocate(sig2_pb)
                    endif
                    call dirs_pb%kill
                    deallocate(rmat_pb, nrm_pb, dirof_pb, cap, sap, o_pb_thr, cap1, sap0)
                    deallocate(seg_dir_b, seg_beg_b, seg_cnt_b)
                    if( allocated(tmpc) ) deallocate(tmpc)
                    if( allocated(stats_es) ) deallocate(stats_es)
                    if( allocated(avec_es) ) deallocate(avec_es, vld_es)
                endif
                deallocate(z, ppinds)
                return
            endif
            endif
            do q = 1, ncomp
                gam_acc(q) = gam_sum(q) / real(max(1,nval),dp)
            end do
            deallocate(gam_sum)
            ! ---- ONLINE EM: blend this batch's sufficient statistics into the running average,
            ! rotating the history into the current basis frame first (rationale at the flag).
            ! MEASURED 2026-08-19 on 10076: the averaging HURTS (0.1936/0.3278/ceil 0.2537 vs
            ! windows-only 0.2361/0.3509/0.2648) -- stale M-step inputs + double regularization
            ! on top of the Wiener. Windowed rotation without averaging is the winner, so
            ! SIMPLE_COV_ONLINE=1 is windows-only and the blend is opt-in at =2 (research).
            if( l_online .and. vonl_pr >= 2 )then
                block
                    real,     pointer     :: onl_pr(:,:,:)
                    real(dp), allocatable :: Mo(:,:), MtMo(:,:), Vo(:,:), evo(:)
                    real,     allocatable :: ut_now(:,:,:,:)
                    complex,  allocatable :: tcme(:,:,:,:), tcmo(:,:,:,:)
                    real,     allocatable :: trhe(:,:,:,:), trho(:,:,:,:)
                    real(dp), allocatable :: tgam(:)
                    real(dp) :: g_onl, wpq
                    integer  :: p2, q2, nrot_o
                    logical  :: l_rot_ok
                    g_onl = (2.d0 + real(it,dp))**(-0.7d0)
                    if( it == 1 .or. .not. l_oh_live ) g_onl = 1.d0
                    ! the CURRENT E-step frame lives in prev_real (utilde_real is per-iteration
                    ! scratch, killed after propagation -- it is NEVER alive at blend time)
                    l_rot_ok = .false.
                    if( allocated(prev_real) .and. l_oh_live .and. allocated(ut_hist) )then
                        if( size(prev_real) == ncomp )then
                            call prev_real(1)%get_rmat_ptr(onl_pr)
                            if( size(ut_hist,1) == size(onl_pr) ) l_rot_ok = .true.
                        endif
                    endif
                    if( g_onl < 1.d0 .and. .not. l_rot_ok )then
                        ! rank changed or first usable frame: restart the average rather than mix frames
                        g_onl = 1.d0
                        write(logfhandle,'(A)') '>>> FLEX_PCA ONLINE: history reset (frame not alignable)'
                    endif
                    if( g_onl < 1.d0 )then
                        ! Procrustes rotation hist-frame -> current frame, from the basis overlap
                        allocate(Mo(ncomp,ncomp), MtMo(ncomp,ncomp), Vo(ncomp,ncomp), evo(ncomp))
                        block
                            real(dp), allocatable :: vnow(:)
                            allocate(vnow(size(ut_hist,1)))
                            do q2 = 1, ncomp
                                call prev_real(q2)%get_rmat_ptr(onl_pr)
                                vnow = real(reshape(onl_pr, [size(onl_pr)]), dp)
                                do p2 = 1, ncomp
                                    Mo(p2,q2) = sum(ut_hist(:,p2) * vnow)
                                end do
                            end do
                        end block
                        MtMo = matmul(transpose(Mo), Mo)
                        call jacobi(MtMo, ncomp, ncomp, evo, Vo, nrot_o)
                        do q2 = 1, ncomp
                            evo(q2) = 1.d0/sqrt(max(evo(q2), 1.d-12))
                        end do
                        do q2 = 1, ncomp
                            MtMo(:,q2) = Vo(:,q2)*evo(q2)     ! scratch: V * diag(1/sqrt(ev))
                        end do
                        Mo = matmul(Mo, matmul(MtMo, transpose(Vo)))
                        ! rotate history: linear for the map sums, R^2 weights for rho/gamma
                        allocate(tcme, mold=oh_cme); allocate(tcmo, mold=oh_cmo)
                        allocate(trhe, mold=oh_rhe); allocate(trho, mold=oh_rho)
                        allocate(tgam(ncomp))
                        tcme = cmplx(0.,0.); tcmo = cmplx(0.,0.); trhe = 0.; trho = 0.; tgam = 0.d0
                        do q2 = 1, ncomp
                            do p2 = 1, ncomp
                                wpq = Mo(p2,q2)
                                if( abs(wpq) < 1.d-6 ) cycle
                                tcme(:,:,:,q2) = tcme(:,:,:,q2) + real(wpq)*oh_cme(:,:,:,p2)
                                tcmo(:,:,:,q2) = tcmo(:,:,:,q2) + real(wpq)*oh_cmo(:,:,:,p2)
                                trhe(:,:,:,q2) = trhe(:,:,:,q2) + real(wpq*wpq)*oh_rhe(:,:,:,p2)
                                trho(:,:,:,q2) = trho(:,:,:,q2) + real(wpq*wpq)*oh_rho(:,:,:,p2)
                                tgam(q2)       = tgam(q2)       + wpq*wpq*oh_gam(p2)
                            end do
                        end do
                        call move_alloc(tcme, oh_cme); call move_alloc(tcmo, oh_cmo)
                        call move_alloc(trhe, oh_rhe); call move_alloc(trho, oh_rho)
                        oh_gam = tgam
                        deallocate(tgam, Mo, MtMo, Vo, evo)
                    endif
                    if( .not. allocated(oh_cme) )then
                        associate( ce => Yeven(1)%cmat_exp )
                            allocate(oh_cme(size(ce,1),size(ce,2),size(ce,3),ncomp))
                            allocate(oh_cmo(size(ce,1),size(ce,2),size(ce,3),ncomp))
                        end associate
                        associate( re => Yeven(1)%rho_exp )
                            allocate(oh_rhe(size(re,1),size(re,2),size(re,3),ncomp))
                            allocate(oh_rho(size(re,1),size(re,2),size(re,3),ncomp))
                        end associate
                        allocate(oh_gam(ncomp))
                    endif
                    ! blend and write back so the M-step consumes the running average
                    do q2 = 1, ncomp
                        if( g_onl >= 1.d0 )then
                            oh_cme(:,:,:,q2) = Yeven(q2)%cmat_exp
                            oh_cmo(:,:,:,q2) = Yodd(q2)%cmat_exp
                            oh_rhe(:,:,:,q2) = Yeven(q2)%rho_exp
                            oh_rho(:,:,:,q2) = Yodd(q2)%rho_exp
                            oh_gam(q2)       = gam_acc(q2)
                        else
                            oh_cme(:,:,:,q2) = real(1.d0-g_onl)*oh_cme(:,:,:,q2) + real(g_onl)*Yeven(q2)%cmat_exp
                            oh_cmo(:,:,:,q2) = real(1.d0-g_onl)*oh_cmo(:,:,:,q2) + real(g_onl)*Yodd(q2)%cmat_exp
                            oh_rhe(:,:,:,q2) = real(1.d0-g_onl)*oh_rhe(:,:,:,q2) + real(g_onl)*Yeven(q2)%rho_exp
                            oh_rho(:,:,:,q2) = real(1.d0-g_onl)*oh_rho(:,:,:,q2) + real(g_onl)*Yodd(q2)%rho_exp
                            oh_gam(q2)       = (1.d0-g_onl)*oh_gam(q2) + g_onl*gam_acc(q2)
                        endif
                        Yeven(q2)%cmat_exp = oh_cme(:,:,:,q2)
                        Yodd(q2)%cmat_exp  = oh_cmo(:,:,:,q2)
                        Yeven(q2)%rho_exp  = oh_rhe(:,:,:,q2)
                        Yodd(q2)%rho_exp   = oh_rho(:,:,:,q2)
                        gam_acc(q2)        = oh_gam(q2)
                    end do
                    l_oh_live = .true.
                    ! snapshot the CURRENT frame as the one the history now lives in
                    if( allocated(prev_real) )then
                        if( size(prev_real) == ncomp )then
                            call prev_real(1)%get_rmat_ptr(onl_pr)
                            if( .not. allocated(ut_hist) ) allocate(ut_hist(size(onl_pr), ncomp))
                            if( size(ut_hist,1) /= size(onl_pr) )then
                                deallocate(ut_hist); allocate(ut_hist(size(onl_pr), ncomp))
                            endif
                            do q2 = 1, ncomp
                                call prev_real(q2)%get_rmat_ptr(onl_pr)
                                ut_hist(:,q2) = real(reshape(onl_pr, [size(onl_pr)]), dp)
                            end do
                        endif
                    endif
                    write(logfhandle,'(A,I0,A,F6.3,A,F9.0)') '>>> FLEX_PCA ONLINE it=', it, &
                        &'  gamma=', g_onl, '  effective N~', real(nval,dp)/max(g_onl,1.d-3)
                    call flush(logfhandle)
                end block
            endif
            ! ---- MARGINAL LIKELIHOOD (the EM objective; NOT resid_energy) ----
            !   -2 log p(y_i) = ||y_i - a_i T_i u_0||^2/sig2  -  h_i' A_i^-1 h_i
            !                   + log det A_i + log det Gamma + const
            ! The first term is FIXED across iterations -- a_i and mu are both held fixed -- so only
            ! the rest is accumulated, and differences between iterations are exact. log det Gamma is
            ! taken from `prior`, which is the Gamma the E-step just USED, not the one it just
            ! produced. Reported per particle so the number is comparable across N.
            ! The filtered M-step (FSC-Wiener + band-limit + mask + re-orthonormalisation) sits
            ! OUTSIDE the likelihood, so this is a DIAGNOSTIC, not a monotonicity guarantee: a rise
            ! means the regularisation is costing more than it buys, which is worth seeing.
            if( l_mix_used )then
                ! the mixture marginal's Omega log det, with the Omega this iteration's E-step
                ! actually USED (the M-step above has already moved mix_Om on)
                nll_tot = nll_tot + real(nval,dp)*ldOm_used
            else
                nll_tot = nll_tot - real(nval,dp)*sum(log(max(prior(1:ncomp), DTINY)))
            endif
            write(logfhandle,'(A,I0,A,ES14.6,A,ES12.4)') '>>> FLEX_PCA PROBE ITER ',it, &
                &'  -2logL/N (varying part)=',nll_tot/real(max(1,nval),dp), &
                &'  delta=',(nll_tot/real(max(1,nval),dp)) - nll_prev
            nll_prev = nll_tot/real(max(1,nval),dp)
            call flush(logfhandle)
            ! COUPLED M-step solve. At every grid point the components share ONE k x k normal matrix
            ! sum_i |CTF|^2 E[z_i z_i'], so the basis volumes have to be solved together. Dividing each
            ! Y_q by the scalar density sum_i |CTF|^2 instead drops the cross-component coupling and
            ! degrades the update to plain subspace iteration. Solve per halfset, then hand a unit
            ! density to compress_exp and skip sampl_dens_correct -- the divide has already happened.
            call solve_coupled_basis_exp(Yeven, rho_e, ncomp)
            call solve_coupled_basis_exp(Yodd,  rho_o, ncomp)
            ! finalize even/odd Y_q, half-set FSC Wiener-merge -> clean band-limited masked basis volume.
            allocate(realvols(ncomp))
            ! the two half-bases are kept, not discarded: comparing them is the only stopping
            ! criterion here that does not embed a dataset-specific constant (see below)
            allocate(eimgs(ncomp), oimgs(ncomp))
            filtsz = max(1, fdim(params%box_crop) - 1)
            allocate(filt(filtsz), corrs(filtsz))
            ! native-lattice inverse KB envelope, applied to each half before FSC/merge/mask
            mstep_gridcorr = prep3D_inv_kbenvelope4mul([params%box_crop,params%box_crop,params%box_crop], params%smpd_crop)
            allocate(fscq_dg(filtsz, ncomp), source=0.)
            do q = 1, ncomp
                ! even half -> band-limited real image (UNmasked, for an unbiased FSC)
                Yeven(q)%rho_exp = 1.0
                call Yeven(q)%compress_exp; call Yeven(q)%ifft
                call realvols(q)%new([params%box_crop,params%box_crop,params%box_crop], params%smpd_crop)
                call Yeven(q)%get_rmat_ptr(rmatp); call realvols(q)%set_rmat(rmatp, .false.)
                call realvols(q)%mul(mstep_gridcorr)
                ! odd half
                Yodd(q)%rho_exp = 1.0
                call Yodd(q)%compress_exp; call Yodd(q)%ifft
                call img_o%new([params%box_crop,params%box_crop,params%box_crop], params%smpd_crop)
                call Yodd(q)%get_rmat_ptr(rmatp); call img_o%set_rmat(rmatp, .false.)
                call img_o%mul(mstep_gridcorr)
                ! half-set FSC -> Wiener filter 2F/(1+F)
                call realvols(q)%fft; call img_o%fft
                call realvols(q)%fsc(img_o, corrs)
                fscq_dg(:,q) = corrs
                do sh = 1, filtsz
                    fc = max(0., min(0.999, corrs(sh)))
                    filt(sh) = 2.*fc/(1.+fc)
                end do
                ! merged, FSC-Wiener + band-limit filtered, back to real, masked
                call realvols(q)%add(img_o); call realvols(q)%mul(0.5)
                call realvols(q)%apply_filter(filt)
                if( lp_it > 2.0*params%smpd_crop + TINY ) call realvols(q)%bp(0., lp_it)
                call realvols(q)%ifft
                if( params%msk_crop > TINY ) call realvols(q)%mask3D_soft(params%msk_crop, backgr=0.)
                ! ---- FOCUSED HETEROGENEITY (SIMPLE_COV_FOCUS_R1, Angstrom): soft inner
                ! exclusion so the basis explains ONLY the peripheral shell (arms/NBD on
                ! 10049 — 10-24% of motion variance lives beyond r=60 A but never owns an
                ! axis against core+background variance). The mean stays global; only the
                ! heterogeneity subspace is focused — focused-classification semantics.
                if( l_focus_annulus )then
                    block
                        real, pointer :: fpr(:,:,:)
                        integer :: fi, fj, fk
                        real    :: frad, fw
                        call realvols(q)%get_rmat_ptr(fpr)
                        do fk = 1, params%box_crop
                            do fj = 1, params%box_crop
                                do fi = 1, params%box_crop
                                    frad = sqrt(real(fi-params%box_crop/2-1)**2 + &
                                        &real(fj-params%box_crop/2-1)**2 + &
                                        &real(fk-params%box_crop/2-1)**2)*params%smpd_crop
                                    fw = min(1.0, max(0.0, (frad - foc_r1)/max(foc_edge, 1.0)))
                                    fpr(fi,fj,fk) = fpr(fi,fj,fk)*fw
                                end do
                            end do
                        end do
                    end block
                endif
                ! stash both halves, filtered exactly as the merged basis is, so the comparison
                ! below is between two bases that have had identical treatment
                call eimgs(q)%copy(realvols(q))
                call img_o%ifft
                if( lp_it > 2.0*params%smpd_crop + TINY )then
                    call img_o%fft; call img_o%bp(0., lp_it); call img_o%ifft
                endif
                if( params%msk_crop > TINY ) call img_o%mask3D_soft(params%msk_crop, backgr=0.)
                if( l_focus_annulus )then
                    block
                        real, pointer :: fpr(:,:,:)
                        integer :: fi, fj, fk
                        real    :: frad, fw
                        call img_o%get_rmat_ptr(fpr)
                        do fk = 1, params%box_crop
                            do fj = 1, params%box_crop
                                do fi = 1, params%box_crop
                                    frad = sqrt(real(fi-params%box_crop/2-1)**2 + &
                                        &real(fj-params%box_crop/2-1)**2 + &
                                        &real(fk-params%box_crop/2-1)**2)*params%smpd_crop
                                    fw = min(1.0, max(0.0, (frad - foc_r1)/max(foc_edge, 1.0)))
                                    fpr(fi,fj,fk) = fpr(fi,fj,fk)*fw
                                end do
                            end do
                        end do
                    end block
                endif
                call oimgs(q)%copy(img_o)
                call img_o%kill
            end do
            call mstep_gridcorr%kill
            ! ---- AUTO-lp / AUTO-neigs CANDIDATE DIAGNOSTIC (report only) ----
            ! Both knobs are derivable from this one measurement, which the Wiener merge
            ! computes anyway. The HETEROGENEITY resolution -- where the leading components'
            ! split-half FSC crosses 0.143 -- is what lp should encode (NOT the consensus
            ! resolution: compositional signal dies far coarser than the consensus map).
            ! The count of components with significant in-band FSC is what neigs should
            ! encode. Validation bar for any rule built on these: reproduce the measured
            ! ladders (band: lp16 beats lp20 delivered; rank: delivered peaks near 40).
            khi_dg = filtsz
            if( params%lp > 2.0*params%smpd_crop + TINY ) &
                &khi_dg = max(2, min(filtsz, int(dstep_ann/params%lp)))
            nsig_dg = 0
            do q = 1, ncomp
                fmean_dg(q) = sum(fscq_dg(1:khi_dg,q))/real(khi_dg)
                if( fmean_dg(q) > 0.143 ) nsig_dg = nsig_dg + 1
            end do
            ! heterogeneity resolution off the STRONGEST components. Ranked by in-band mean
            ! FSC, NOT by component index: the EM basis is orthonormalised in fixed order and
            ! never variance-sorted (measured on 10076: the 39th/40th components carried the
            ! 2nd/3rd largest latent sd), so "first 4" would read arbitrary components.
            ntop_dg = min(4, ncomp)
            sel_dg  = 0
            do tq_dg = 1, ntop_dg
                fbest_dg = -huge(1.d0); bq_dg = 1
                do q = 1, ncomp
                    if( any(sel_dg(1:tq_dg-1) == q) ) cycle
                    if( fmean_dg(q) > fbest_dg )then
                        fbest_dg = fmean_dg(q); bq_dg = q
                    endif
                end do
                sel_dg(tq_dg) = bq_dg
            end do
            res_dg  = dstep_ann/real(khi_dg)
            do sh = 2, filtsz
                fc = 0.
                do tq_dg = 1, ntop_dg
                    fc = fc + fscq_dg(sh,sel_dg(tq_dg))
                end do
                fc = fc/real(ntop_dg)
                if( fc < 0.143 )then
                    res_dg = dstep_ann/real(sh)
                    exit
                endif
            end do
            write(logfhandle,'(A,I0,A,F7.2,A,I0,A,I0,A)') '>>> FLEX_PCA BAND/RANK it=',it, &
                &'  het-resolution(FSC 0.143)=',res_dg,' A   components with in-band FSC>0.143: ', &
                &nsig_dg,' of ',ncomp,'   (auto-lp / auto-neigs candidates)'
            call flush(logfhandle)
            deallocate(fscq_dg)
            deallocate(filt, corrs)
            ! ---- DATASET-AGNOSTIC STOPPING RULE ----
            ! Every fixed iteration count is a tuned constant standing in for a convergence test,
            ! and the number that is right for one dataset is wrong for the next. The successive-
            ! basis principal angle is not the fix either: it averages over ALL components, so it
            ! is dominated by the tail directions that never reproduce, and it fires at iteration
            ! 3-5 while the leading directions are still moving.
            !
            ! What the algorithm already has, and throws away, is TWO INDEPENDENT HALF-FITS: Yeven
            ! and Yodd are solved separately every iteration and only then merged. The principal
            ! angles between those two spans measure the REPRODUCIBLE dimension of the current
            ! basis directly, and their sum is a soft count of it. Iterating past the point where
            ! that stops rising is refining structure that does not survive a change of particles
            ! -- measured on 10076: the leading eight directions are converged by iteration ~20
            ! while the likelihood keeps falling to iteration 80, and only ~8 reproduce.
            !
            ! So: stop when the reproducible dimension has not improved for COV_EO_PATIENCE
            ! iterations. No particle count, no band, no dataset in the rule.
            ! Comparing the half-bases DIRECTLY does not work, and the reason is worth recording:
            ! both halves are updates to the SAME current basis, so both inherit it and the
            ! agreement saturates immediately -- measured 15.3 of 16 at iteration 1, which is the
            ! shared basis showing through, not reproducible signal.
            !
            ! The quantity that is actually informative is whether the two halves agree about what
            ! to ADD. Projecting each half-basis onto the orthogonal complement of the previous
            ! basis leaves exactly the new directions each half wants, and the principal angles
            ! between THOSE say whether the current update carries signal or noise. When the
            ! updates stop agreeing, further iterations are fitting the half they happen to see.
            if( allocated(prev_real) )then
                call deflate_against_basis(eimgs, ncomp, prev_real, size(prev_real))
                call deflate_against_basis(oimgs, ncomp, prev_real, size(prev_real))
                call cross_half_subspace_angles(eimgs, oimgs, ncomp, sv_eo)
                eo_dim = sum(sv_eo)
            else
                allocate(sv_eo(ncomp), source=0.d0)
                eo_dim = -1.d0            ! first iteration: no previous basis to deflate against
            endif
            write(logfhandle,'(A,I0,A,F8.3,A,I0,A,I0)') '>>> FLEX_PCA PROBE ITER ',it, &
                &'  update agreement (even|odd, prev basis deflated)=',eo_dim, &
                &'   n>=0.9: ',count(sv_eo >= 0.9d0),' of ',ncomp
            call flush(logfhandle)
            if( eo_dim < 0.d0 )then
                continue                  ! nothing to judge yet
            else if( eo_dim > eo_best + COV_EO_TOL )then
                eo_best = eo_dim; eo_stall = 0
            else
                eo_stall = eo_stall + 1
            endif
            deallocate(sv_eo)
            do q = 1, ncomp
                call eimgs(q)%kill; call oimgs(q)%kill
            end do
            deallocate(eimgs, oimgs)
            ! ---- MEAN-SHAPED (CONTRAST) DEFLATION, SIMPLE_COV_EM_DEFLATE=1 ----
            ! a_i is a scalar ML fit of the contrast against the mean, so any FREQUENCY-dependent
            ! part of the per-image scale (envelope/B-factor spread) survives it and lands in the
            ! residual as a term proportional to the consensus shape, coherent across every particle.
            ! It is therefore the strongest single direction available and the EM spends a whole
            ! component on it: measured on 10076, pc001 is a featureless blob in the EM arm AND in
            ! the moment arm. The column path already deflates it in volume space (SIMPLE_COV_RANK1);
            ! this is the same operation applied to the refined basis each iteration.
            ! Only the GLOBAL mean-shaped component is removed -- a domain-local mass change keeps
            ! whatever part of it is orthogonal to the consensus, which is most of it.
            ! SIMPLE_COV_EM_DEFLATE=n with n>1 deflates an n-dimensional ENVELOPE subspace instead
            ! of the single mean direction, which is the gap the paragraph above names. Splitting
            ! the fitted band into n equal shells of the consensus and removing all of them takes
            ! out any per-particle scale that varies smoothly with frequency -- a B-factor spread,
            ! an envelope difference, an ice-thickness gradient -- not just the one scalar a_i
            ! already absorbs. Disjoint Fourier bands are orthogonal by Parseval, but soft-masking
            ! in real space breaks that, so the shells are explicitly orthonormalised rather than
            ! assumed independent.
            if( l_deflate_mean )then
                ndfl = max(1, vdfl)
                ! SIMPLE_COV_DEFLATE_BG=1 appends a FLAT-IN-MASK background template to the
                ! deflation set: the ice/background mode is NOT consensus-shaped (measured
                ! overlap 0.034 on 10049) so consensus shells never remove it, and it owns
                ! the leading axis (27x) on real data whenever the mask is wide enough to
                ! admit solvent. Deflating an explicit background direction lets the mask
                ! OPEN (arms/NBD variance survives) without the nuisance eating the basis.
                if( cov_env_int_on('SIMPLE_COV_DEFLATE_BG') )then
                    ndfl = ndfl + 1
                    l_dfl_bg = .true.
                else
                    l_dfl_bg = .false.
                endif
                ! SIMPLE_COV_DEFLATE_POSE=1: six global pose-derivative templates (below)
                if( cov_env_int_on('SIMPLE_COV_DEFLATE_POSE') )then
                    ndfl = ndfl + 6
                    l_dfl_pose = .true.
                else
                    l_dfl_pose = .false.
                endif
                ndfl_sh = ndfl
                if( l_dfl_bg )   ndfl_sh = ndfl_sh - 1
                if( l_dfl_pose ) ndfl_sh = ndfl_sh - 6
                allocate(dfl_basis(ndfl))
                call mvol_dfl%read_and_crop(params%vols(1), params%smpd, params%box_crop, params%smpd_crop)
                kfr_dfl = covariance_kfromto(params)
                if( l_dfl_bg )then
                    ! last slot: uniform density under the soft spherical mask
                    call dfl_basis(ndfl)%copy(mvol_dfl)
                    call dfl_basis(ndfl)%get_rmat_ptr(rm_dfl)
                    rm_dfl = 0.
                    rm_dfl(1:params%box_crop, 1:params%box_crop, 1:params%box_crop) = 1.
                    if( params%msk_crop > TINY ) call dfl_basis(ndfl)%mask3D_soft(params%msk_crop, backgr=0.)
                    write(logfhandle,'(A)') '>>> FLEX_PCA DEFLATE_BG: flat-in-mask background template &
                        &appended to the deflation set'
                    call flush(logfhandle)
                endif
                do idfl = 1, ndfl_sh
                    call dfl_basis(idfl)%copy(mvol_dfl)
                    if( ndfl_sh > 1 )then
                        ! shell idfl of ndfl, in resolution not shell index, so the bands carry
                        ! comparable signal rather than the outermost one carrying nearly all of it
                        ! dstep_crop / shell, matching covariance_kfromto's own convention.
                        ! The FIRST shell is a pure low-pass: a high-pass edge at k~1 would shave
                        ! the shells where nearly all the consensus power lives.
                        if( idfl == 1 )then
                            res_lo = 0.
                        else
                            res_lo = dstep_ann / &
                                &max(1., real(kfr_dfl(1)) + real(idfl-1)*real(max(1,kfr_dfl(2)-kfr_dfl(1)))/real(ndfl_sh))
                        endif
                        res_hi = dstep_ann / &
                            &max(1., real(kfr_dfl(1)) + real(idfl)  *real(max(1,kfr_dfl(2)-kfr_dfl(1)))/real(ndfl_sh))
                        call dfl_basis(idfl)%fft
                        ! width=1: bp's DEFAULT cosine edge is 10 shells wide PER SIDE, and at
                        ! ndfl=2 each shell of the fitted band is only ~8 shells, so the default
                        ! obliterated the shells and unit-normalisation amplified the residue into
                        ! directions near-orthogonal to everything: rank-2 deflation removed 0.12%
                        ! of basis energy where rank-1 removed 34.6% (measured, hd_defl2)
                        call dfl_basis(idfl)%bp(res_lo, res_hi, width=1.0)
                        call dfl_basis(idfl)%ifft
                        ! real-space padding columns are not guaranteed zero after an ifft, and
                        ! every inner product below runs over the whole padded array
                        call dfl_basis(idfl)%get_rmat_ptr(rm_dfl)
                        rm_dfl(params%box_crop+1:,:,:) = 0.
                        ! band-passing spreads density outside the particle, and the volumes being
                        ! deflated are soft-masked, so the shells must be too or the projection is
                        ! taken against solvent. Deliberately NOT applied at ndfl=1: that is the
                        ! arm measured at ARI 0.176 and it stays bit-identical.
                        if( params%msk_crop > TINY ) call dfl_basis(idfl)%mask3D_soft(params%msk_crop, backgr=0.)
                    endif
                end do
                ! SIMPLE_COV_DEFLATE_POSE=1 appends the six GLOBAL pose-derivative
                ! templates: translation gradients d(mu)/dx,dy,dz and whole-map rotation
                ! derivatives -(a x r).grad(mu) about x/y/z. Per-particle alignment error
                ! IS a small global rotation/shift, so its variance is nuisance in every
                ! dataset regardless of heterogeneity type; DIFFERENTIAL (subunit-local)
                ! motion is orthogonal to whole-map rotation and survives the deflation.
                ! Measured need: on EMPIAR-10028 the delivered "ratcheting" states were
                ! globally rotated consensus copies (pose jitter), and on 3.42 A cryoSPARC
                ! poses that mode collapses -- so deflate it at source. Slots ndfl_sh+1..+6;
                ! Gram-Schmidt below orthonormalises them against the consensus shells.
                if( l_dfl_pose )then
                    call mvol_dfl%get_rmat_ptr(rv_dfl)
                    cp_dfl = real(params%box_crop/2 + 1, dp)
                    do ipdfl = 1, 6
                        call dfl_basis(ndfl_sh+ipdfl)%copy(mvol_dfl)
                        call dfl_basis(ndfl_sh+ipdfl)%get_rmat_ptr(rm_dfl)
                        rm_dfl = 0.
                        do izp = 2, params%box_crop-1
                            do iyp = 2, params%box_crop-1
                                do ixp = 2, params%box_crop-1
                                    gpx = 0.5d0*(real(rv_dfl(ixp+1,iyp,izp),dp) - real(rv_dfl(ixp-1,iyp,izp),dp))
                                    gpy = 0.5d0*(real(rv_dfl(ixp,iyp+1,izp),dp) - real(rv_dfl(ixp,iyp-1,izp),dp))
                                    gpz = 0.5d0*(real(rv_dfl(ixp,iyp,izp+1),dp) - real(rv_dfl(ixp,iyp,izp-1),dp))
                                    select case(ipdfl)
                                        case(1); rm_dfl(ixp,iyp,izp) = real(gpx)
                                        case(2); rm_dfl(ixp,iyp,izp) = real(gpy)
                                        case(3); rm_dfl(ixp,iyp,izp) = real(gpz)
                                        case(4); rm_dfl(ixp,iyp,izp) = &
                                            &real((real(izp,dp)-cp_dfl)*gpy - (real(iyp,dp)-cp_dfl)*gpz)
                                        case(5); rm_dfl(ixp,iyp,izp) = &
                                            &real((real(ixp,dp)-cp_dfl)*gpz - (real(izp,dp)-cp_dfl)*gpx)
                                        case(6); rm_dfl(ixp,iyp,izp) = &
                                            &real((real(iyp,dp)-cp_dfl)*gpx - (real(ixp,dp)-cp_dfl)*gpy)
                                    end select
                                end do
                            end do
                        end do
                        if( params%msk_crop > TINY ) call dfl_basis(ndfl_sh+ipdfl)%mask3D_soft(params%msk_crop, backgr=0.)
                    end do
                    write(logfhandle,'(A)') '>>> FLEX_PCA DEFLATE_POSE: 6 global pose-derivative templates &
                        &appended to the deflation set'
                    call flush(logfhandle)
                endif
                ! shell health: a shell whose norm is tiny relative to the consensus was emptied
                ! by the band edges or the mask, and unit-normalising it manufactures a spurious
                ! deflation direction out of numerical residue -- which is exactly how rank-2
                ! went inert. Print the norms, and keep only shells above a RELATIVE floor.
                call mvol_dfl%get_rmat_ptr(rv_dfl)
                mnorm_dfl = sqrt(sum(real(rv_dfl,dp)**2))
                if( ndfl > 1 )then
                    write(logfhandle,'(A)',advance='no') '>>> FLEX_PCA deflation shell norms / |consensus|:'
                    do idfl = 1, ndfl
                        call dfl_basis(idfl)%get_rmat_ptr(rm_dfl)
                        write(logfhandle,'(A,ES9.2)',advance='no') ' ', &
                            &sqrt(sum(real(rm_dfl,dp)**2))/max(mnorm_dfl,DTINY)
                    end do
                    write(logfhandle,*)
                endif
                ! modified Gram-Schmidt; a shell that the band edges or the mask emptied drops out
                nkeep_dfl = 0
                do idfl = 1, ndfl
                    call dfl_basis(idfl)%get_rmat_ptr(rm_dfl)
                    do jdfl = 1, nkeep_dfl
                        call dfl_basis(jdfl)%get_rmat_ptr(rv_dfl)
                        mv_dfl = sum(real(rm_dfl,dp)*real(rv_dfl,dp))
                        rm_dfl = rm_dfl - real(mv_dfl)*rv_dfl
                    end do
                    mm_dfl = sqrt(sum(real(rm_dfl,dp)*real(rm_dfl,dp)))
                    if( mm_dfl <= 1.d-3*mnorm_dfl ) cycle
                    rm_dfl    = rm_dfl / real(mm_dfl)
                    nkeep_dfl = nkeep_dfl + 1
                    if( nkeep_dfl /= idfl ) call dfl_basis(nkeep_dfl)%copy(dfl_basis(idfl))
                end do
                if( nkeep_dfl > 0 )then
                    rem_dfl = 0.d0; tot_dfl = 0.d0
                    do q = 1, ncomp
                        call realvols(q)%get_rmat_ptr(rv_dfl)
                        tot_dfl = tot_dfl + sum(real(rv_dfl,dp)**2)
                        do idfl = 1, nkeep_dfl
                            call dfl_basis(idfl)%get_rmat_ptr(rm_dfl)
                            mv_dfl  = sum(real(rm_dfl,dp)*real(rv_dfl,dp))   ! unit-norm already
                            rem_dfl = rem_dfl + mv_dfl*mv_dfl
                            rv_dfl  = rv_dfl - real(mv_dfl)*rm_dfl
                        end do
                    end do
                    write(logfhandle,'(A,I0,A,I0,A,F6.2,A)') '>>> FLEX_PCA PROBE ITER ',it, &
                        &'  mean-shaped deflation (rank ',nkeep_dfl,') removed ', &
                        &100.d0*rem_dfl/max(tot_dfl,DTINY),' % of basis energy'
                    call flush(logfhandle)
                endif
                do idfl = 1, ndfl
                    call dfl_basis(idfl)%kill
                end do
                deallocate(dfl_basis)
                call mvol_dfl%kill
            endif
            ! orthonormalize the probe volumes -> refined basis
            call orthonormalize_representatives(params, build, realvols, ncomp, utilde, utilde_real, d_new)
            ! ---- CONVERGENCE: principal angles between successive bases ----
            ! This is what n_probe_iters should have been. The measured ladder on Ribosembly saturates
            ! at 2 (25-NN 0.7712 -> 0.7715 from 2 to 3, ceiling pinned at 14/16) while the max-variance
            ! that gets logged keeps oscillating 1.96e4 / 6.71e3 / 1.40e4 -- so the basis had settled
            ! while the quantity being watched had not, and the iteration count was a tuned constant
            ! standing in for a convergence test. Both bases come out of orthonormalize_representatives
            ! ORTHONORMAL, so the singular values of their cross-Gram really are principal-angle
            ! cosines here; align_basis_to_reference only per-vector normalises, which is a no-op on an
            ! already-orthonormal set (it is NOT safe on arbitrary bases -- see its own caveat).
            l_converged = .false.
            if( allocated(prev_real) )then
                call align_basis_to_reference(prev_real, size(prev_real), utilde_real, d_new, Mconv, sconv)
                cos_mean = sum(sconv) / real(max(1,size(sconv)),dp)
                write(logfhandle,'(A,I0,A,F9.6)') '>>> FLEX_PCA PROBE ITER ',it, &
                    &'  mean principal-angle cosine vs previous basis=',cos_mean
                call flush(logfhandle)
                ! The successive-basis cosine is kept as a REPORTED diagnostic only. It averages
                ! over all components, so the never-reproducing tail dominates it and it fires far
                ! too early -- on 10076 at iteration 3-5, with the leading directions still moving
                ! and the likelihood an order of magnitude from its plateau. It becomes the
                ! criterion again only if the even/odd signal is unavailable.
                if( ncomp < 2 .and. cos_mean >= conv_thresh ) l_converged = .true.
                deallocate(Mconv, sconv)
                do q = 1, size(prev_real)
                    call prev_real(q)%kill
                end do
                deallocate(prev_real)
            endif
            allocate(prev_real(d_new))
            do q = 1, d_new
                call prev_real(q)%copy(utilde_real(q))
            end do
            ! replace basis_recs with the refined (projection-ready) basis; eigvals = latent variances
            do q = 1, size(basis_recs)
                call basis_recs(q)%dealloc_rho; call basis_recs(q)%kill
            end do
            deallocate(basis_recs); allocate(basis_recs(d_new))
            if( allocated(eigvals) ) deallocate(eigvals); allocate(eigvals(d_new))
            do q = 1, d_new
                ! projection-ready basis reconstructor from the clean real basis image (mean_rec idiom)
                call init_basis_reconstructor(params, build, basis_recs(q))
                call basis_recs(q)%set_rmat(utilde_real(q)%get_rmat(), .false.)   ! NOT add(): see above
                call basis_recs(q)%fft
                call basis_recs(q)%expand_exp
                ! EM Gamma update: the POSTERIOR second moment (1/n) sum_i (z_iq^2 + [A_i^-1]_qq), not the
                ! sample variance of the MAP point estimates. MAP latents are shrunk toward zero by the
                ! prior, so the point-estimate variance underestimates Gamma, which sets a tighter prior
                ! next iteration, which shrinks them further -- a self-reinforcing collapse.
                eigvals(q) = max(gam_acc(min(q,ncomp)), DTINY)
                ! overwrite the eigenvolume MRC with the refined basis vector
                fname = 'flex_pca_pc'//int2str_pad(q,3)//MRC_EXT
                call utilde_real(q)%write(fname, del_if_exists=.true.); call fname%kill
            end do
            ! posterior_frac is the share of Gamma that the old point-estimate update was throwing away;
            ! mean(z) should sit near zero -- drift there means the consensus mean mu is off.
            write(logfhandle,'(A,I0,A)') '>>> FLEX_PCA PROBE EM Gamma update (n=',nval,' valid particles):'
            do q = 1, min(ncomp,10)
                mu_q = sum(z(:,q)) / real(max(1,npp),dp)
                sd_q = sum(z(:,q)**2) / real(max(1,nval),dp)
                write(logfhandle,'(A,I3,A,ES11.3,A,ES11.3,A,F7.3,A,ES10.2)') '>>>   z',q, &
                    &'  <z^2>=',sd_q,'  Gamma=',gam_acc(q), &
                    &'  posterior_frac=',real(1.d0 - sd_q/max(gam_acc(q),DTINY)), &
                    &'  mean(z)=',mu_q
            end do
            call flush(logfhandle)
            ncomp = d_new
            write(logfhandle,'(A,I0,A,ES12.4,A,ES12.4,A,F8.1)') '>>> FLEX_PCA PROBE ITER ',it, &
                &' refined dim=',real(ncomp),' max var=',maxval(eigvals),' seconds=',toc(t_it)
            call flush(logfhandle)
            ! cleanup per-iteration scratch
            do i = 1, size(orientations); call orientations(i)%kill; end do
            do ithr = 1, nthr
                call cleanup_plane(mean_fpl(ithr))
                do q = 1, size(basis_fpls,1); call cleanup_plane(basis_fpls(q,ithr)); end do
                if( l_bank_fwd )then
                    call cleanup_plane(mean_fplC(ithr))
                    do q = 1, size(basis_fplsC,1); call cleanup_plane(basis_fplsC(q,ithr)); end do
                endif
            end do
            if( l_bank_fwd ) deallocate(mean_fplC, basis_fplsC)
            if( allocated(stats_es) ) deallocate(stats_es)
            call cleanup_rec_buffers(build, fpls)
            do q = 1, size(Yeven); call Yeven(q)%dealloc_rho; call Yeven(q)%kill; call Yodd(q)%dealloc_rho; call Yodd(q)%kill; end do
            do q = 1, size(utilde); call utilde(q)%dealloc_rho; call utilde(q)%kill; call utilde_real(q)%kill; end do
            do q = 1, size(realvols); call realvols(q)%kill; end do
            deallocate(Yeven, Yodd, utilde, utilde_real, realvols, rho_e, rho_o, prior)
            deallocate(Gth, Ath, bth, cth, zth, basis_fpls, mean_fpl, orientations, zbatch, dens, valid, valid_e, valid_o, eo)
            deallocate(Ainvth, Acpth, gam_thr, gam_acc, nval_thr, hth, nll_thr)
            ! The MCFA state (mix_*) and its work arrays (rhs0th etc) are deliberately NOT
            ! deallocated here: this block runs at the end of EVERY iteration, and the mixture
            ! must survive from one iteration's M-step to the next iteration's E-step. Putting
            ! their cleanup here was the crash: iteration 3 initialised the mixture, this block
            ! freed it, l_mix_active stayed true, and iteration 4's E-step walked into an
            ! unallocated mix_Ominv -- a null descriptor, which is why -fcheck=bounds stayed
            ! silent while every thread segfaulted on the same line. They are subroutine-local
            ! allocatables, so Fortran frees them at exit; the resize block at iteration start
            ! handles dimension changes.
            if( allocated(gam_dbg) ) deallocate(gam_dbg)
            if( allocated(pGc) ) deallocate(pGc, pbc, pcc, pzh, pcon, pre, prme)
            ! NOT used to stop. The update agreement decays monotonically from the first
            ! iteration (measured 15.0 -> 13.7 -> 10.6 -> 8.2, strongly-reproducible update
            ! directions 15 -> 2 -> 0), so any "peaked and stalled" rule fires at iteration 4-5.
            ! But the basis demonstrably keeps improving to ~20 on the same data, because an
            ! update whose DIRECTION is not reproducible can still carry a reproducible bias that
            ! accumulates. So this is a genuine diagnostic of update quality and NOT a stopping
            ! rule, and pretending otherwise would just be a new tuned constant wearing a
            ! principled hat. The defensible in-run rule is cross-fit agreement between two
            ! INDEPENDENT half-fits, which costs 2x; see the acceleration note.
            if( .false. .and. eo_stall >= eo_patience ) l_converged = .true.
            if( l_converged )then
                write(logfhandle,'(A,I0,A,F8.3,A,I0,A)') '>>> FLEX_PCA PROBE converged after ',it, &
                    &' iterations: reproducible dim (even|odd) peaked at ',eo_best, &
                    &' and did not improve for ',eo_patience,' iterations'
                call flush(logfhandle)
                exit
            endif
        end do
        if( allocated(prev_real) )then
            do q = 1, size(prev_real)
                call prev_real(q)%kill
            end do
            deallocate(prev_real)
        endif
        if( l_pol_es )then
            call polar_grid_kill(pg_es)
            if( allocated(UsallE) ) deallocate(UsallE, CfE, Cm0E, c00E, UbankE, CspE, &
                &xws_es, wr_es, wrd_es, Reb_es)
            if( allocated(dir_es) ) deallocate(dir_es, cae, sae, dused_es, rmatb_es, nrmb_es)
            if( allocated(hex_es) ) deallocate(hex_es, kex_es)
            if( l_pol_gpu_any )then
                call flex_gpu_poles_free_f
                ! resident hybrid volumes (idempotent; also freed by the banked-adjoint teardown)
                call flex_gpu_estep_free_f
            endif
            if( allocated(stats_pg) ) deallocate(stats_pg)
            if( allocated(vld_pg)   ) deallocate(vld_pg, dir_pg, ca_pg, sa_pg)
        endif
        if( l_bank_adj )then
            call flex_gpu_coupled_bank_free_f
            call flex_gpu_estep_free_f
            if( l_gpu_prep )then
                call flex_gpu_prep_free_f
                deallocate(ctfp_arr, shf_pb)
                if( allocated(sig2_pb) ) deallocate(sig2_pb)
            endif
            call dirs_pb%kill
            deallocate(rmat_pb, nrm_pb, dirof_pb, cap, sap, o_pb_thr, cap1, sap0)
            deallocate(seg_dir_b, seg_beg_b, seg_cnt_b)
            if( allocated(tmpc) ) deallocate(tmpc)
            if( allocated(avec_es) ) deallocate(avec_es, vld_es)
        endif
        deallocate(z, ppinds)

    contains

        subroutine ensure_bank_ctf_planes( ref, mfc, bfc, nc )
            type(fplane_type), intent(in)    :: ref
            type(fplane_type), intent(inout) :: mfc
            integer,           intent(in)    :: nc
            type(fplane_type), intent(inout) :: bfc(nc)
            integer :: qq
            if( .not. allocated(mfc%cmplx_plane) )then
                allocate(mfc%cmplx_plane, mold=ref%cmplx_plane)
                mfc%cmplx_plane = CMPLX_ZERO
            endif
            mfc%frlims = ref%frlims; mfc%nyq = ref%nyq
            do qq = 1, nc
                if( .not. allocated(bfc(qq)%cmplx_plane) )then
                    allocate(bfc(qq)%cmplx_plane, mold=ref%cmplx_plane)
                    bfc(qq)%cmplx_plane = CMPLX_ZERO
                endif
                bfc(qq)%frlims = ref%frlims; bfc(qq)%nyq = ref%nyq
            end do
        end subroutine ensure_bank_ctf_planes

    end subroutine probe_subspace_iteration

end submodule simple_flex_pca_em_iter
