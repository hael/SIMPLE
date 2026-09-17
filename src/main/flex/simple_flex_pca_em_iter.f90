!@descr: flex_pca EM: the subspace EM iteration (E-step accumulation, M-step update)
submodule (simple_flex_pca_em) simple_flex_pca_em_iter
use simple_matcher_3Drec,   only: init_rec, cleanup_rec_buffers
use simple_matcher_ptcl_io, only: discrete_read_imgbatch, prepimgbatch
use simple_flex_reconstructor_latent_ops, only: project_fplane_mean, project_fplanes_mean_basis,&
    &insert_planes_oversamp_coupled_batch_scaled
use simple_flex_reconstructor_latent_ops, only: prep_imgs4projected_model, solve_coupled_basis_exp,&
    &projected_model_kfromto, add_invtausq2rho_coupled
use simple_flex_pca_crossfsc, only: crossfsc_file, crossfsc_record, crossfsc_load, crossfsc_write,&
    &crossfsc_append, crossfsc_latest_upto, crossfsc_kill, crossfsc_kill_record, crossfsc_to_invtau2,&
    &crossfsc_harvest_h, crossfsc_stop_stat, crossfsc_inband_mean, crossfsc_khi_deepest,&
    &crossfsc_assert_paired, COV_XFSC_FNAME
use simple_flex_gpu,        only: flex_gpu_available, flex_gpu_coupled_begin_f,&
    &flex_gpu_coupled_batch_raw_f, flex_gpu_coupled_end_f, flex_gpu_coupled_bank_f,&
    &flex_gpu_coupled_batch_banked_f, flex_gpu_coupled_bank_free_f, flex_gpu_estep_vols_f,&
    &flex_gpu_estep_batch_f, flex_gpu_estep_resid_f, flex_gpu_estep_free_f,&
    &flex_gpu_coupled_batch_banked_res_f, flex_gpu_prep_begin_f, flex_gpu_prep_batch_f,&
    &flex_gpu_prep_free_f, flex_gpu_estep_batch_res_f, flex_gpu_prep_check_f,&
    &flex_gpu_poles_begin_f, flex_gpu_poles_bank_f, flex_gpu_poles_batch_f, flex_gpu_poles_free_f
use simple_flex_pca_util,   only: cov_env_flag_on, cov_env_flag_off
use simple_flex_pca_polar,  only: polar_grid_build, polar_grid_kill, polar_project_recs,&
    &polar_relative_inplane, polar_assign_directions, polar_sample_particle_fused
implicit none
#include "simple_local_flags.inc"

contains

    !> Probe-based subspace iteration: alternate a Wiener E-step (per-particle latents in the current
    !! basis) with a weighted-backprojection M-step (Y_q += sum_i z_iq * backproject(r_i)), then
    !! orthonormalize the refined probe volumes into the next basis.
    module subroutine probe_subspace_iteration( params, build, mean_rec, basis_recs, eigvals, sig2_eff, &
        &pinds, nptcls, ncomp, niters, it_glob, niters_glob, fprefix, meta_fname, rounds)
        use simple_flex_pca_plane_cache, only: plane_cache_in_use, plane_cache_read_batch
        use simple_flex_pca_planes, only: planes_batch_load
        class(flex_pca_rounds), intent(inout) :: rounds
        class(parameters),   intent(inout) :: params
        type(builder),       intent(inout) :: build
        type(reconstructor), intent(inout) :: mean_rec
        type(reconstructor), allocatable, intent(inout) :: basis_recs(:)
        real(dp),            allocatable, intent(inout) :: eigvals(:)
        real(dp),            intent(in)    :: sig2_eff
        integer,             intent(in)    :: pinds(:), nptcls, niters
        integer,             intent(inout) :: ncomp
        !> master's global-iteration stamp, passed by a distributed probe worker whose own loop
        !! runs once per relaunch; absent (or 0) in shared-memory runs and on the master itself
        integer, optional,   intent(in)    :: it_glob, niters_glob
        !> file namespace override (merge-polish pass: never clobber fit A's delivery)
        character(len=*), optional, intent(in) :: fprefix, meta_fname
        type(fplane_type), allocatable :: fpls(:)
        type(ori),           allocatable :: orientations(:)
        type(reconstructor), allocatable :: utilde(:)
        type(image),         allocatable :: realvols(:), utilde_real(:)
        type(image) :: img_o, mstep_gridcorr
        real,                allocatable :: filt(:), corrs(:)
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
        !> marginal-likelihood scratch: h_i before the solve, and the per-thread accumulator of
        !! log det A_i - h_i' A_i^-1 h_i. The remaining terms of -2 log p(y_i) are the per-particle
        !! data energy (fixed across iterations, since a_i and mu are) and N*log det Gamma, added at
        !! the reduction. See spd_logdet_dp for why resid_energy could never serve this purpose.
        real(dp) :: ldA
        !> principal-angle convergence threshold. Overridable because the default 0.97 was tuned for a
        !! probe that STARTED from the covariance eigenbasis and only had to clean it: from a data-free
        !! start the subspace stops rotating several iterations before Gamma stops moving (measured on
        !! 10076: cos crosses 0.97 at iteration 5 while max var is still falling ~20 % per iteration),
        !! so a pure-EM fit needs a tighter bar. SIMPLE_COV_PROBE_CONV is in PER MILLE (995 = 0.995).
        !> even/odd half-basis agreement: the dataset-agnostic convergence signal
        type(image), allocatable :: eimgs(:), oimgs(:)
        real(dp),    allocatable :: sv_eo(:)
        real(dp) :: eo_dim
        integer, parameter :: MIX_ZSUB_MAX = 2000   ! per-part latent subsample shipped for mcfa_init
        !> POLAR SHARED-DIRECTION E-STEP, SIMPLE_COV_POLAR_ESTEP=1 (stage 1).
        !! IN the batch loop: once per
        !! EM iteration the mean + ncomp basis volumes are projected onto the polar rings of a
        !! shared direction table (the bank); per particle, the bank rings at (direction, psi)
        !! replace the ncomp+1 per-particle Cartesian projections in forming G/b/c/e_mm/myv, and
        !! everything downstream -- the ECM/MCFA solve, mixture accumulators,
        !! e/o halfsets, the Cartesian M-step insertion -- is untouched. CTF amplitude, shift
        !! phase and per-shell whitening ride in exactly as in the Cartesian former, because the
        !! SAME prepped cmplx/transfer planes are polar-sampled (the validated embed formulation).
        !! The mean-shaped deflation needs no special handling here: dfl_basis deflates the
        !! refined basis VOLUMES at the end of each M-step, and the bank is rebuilt from those
        !! same basis_recs at the next iteration, so the deflation enters the bank exactly as it
        !! enters the per-particle Cartesian projections.
        integer  :: id_es, idp_es, ir
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
        integer  :: jx_es, kx_es
        type(oris)            :: dirs_es
        type(ori)             :: o_es
        real,     allocatable :: rmatp_es(:,:,:), nrmp_es(:,:)
        real(dp) :: pw_es, cnt_es
        real     :: taz_es
        integer(timer_int_kind) :: t_bank
        !> DEVICE polar E-step (stage 2), SIMPLE_COV_POLAR_GPU=1 (requires POLAR_ESTEP=1): the
        !! host-built bank + ring Gram tables upload once per iteration; per batch the RAW
        !! (y, T) plane sections go up and per-particle G/b/c/e_mm/myv come back in the fused
        !! E-step's stat layout -- the ring part via a device sampler + bank/table contractions,
        !! the hybrid low-k exact part via the unmodified estep_stats_kernel with the disc
        !! capped at rhyb*(rhyb+1) against the resident mean+basis volumes. The solve (plain or
        !! mixture), Gamma/nll accumulation, e/o bookkeeping, the Cartesian
        !! mean projection/residual and the M-step insertion are the CPU polar path's,
        !! untouched. All device state is per-stage (ring geometry), per-iteration (bank,
        !! volumes) or per-batch (planes).
        logical  :: l_pol_gpu, l_pol_gpu_any
        integer  :: vpolg, nst_pg
        real(dp), allocatable :: stats_pg(:,:)
        real,     allocatable :: ca_pg(:), sa_pg(:)
        integer,  allocatable :: dir_pg(:)
        logical,  allocatable :: vld_pg(:)
        !> gate-1 GPU parity arrays: the CPU polar former re-run on the check particles
        !! iteration 1 run BOTH formers; per-component z correlation and the median relative
        !! G/b errors are printed, then the run continues on the polar path.
        !> exact-direction polar arm of the check: the SAME quadrature at the particle's own
        !! direction (no snap), which attributes any disagreement between direction quantization
        !! (banked vs exact) and the polar formulation itself (exact vs Cartesian). tazim -- the
        !! mean within-ring relative spread of |T|^2 -- bounds the radial-factorisation error in G.
        !> DC-less Cartesian arm of the check: cov_herm_inner includes the (0,0) sample, the polar
        !! quadrature starts at ring 1 and has none. If the exact-direction polar arm agrees with
        !! THIS reference, the polar/Cartesian gap is the DC term, not the ring quadrature.
        logical  :: l_probe_distr_pre
        !> mean-shaped (contrast) deflation of the refined basis, SIMPLE_COV_EM_DEFLATE
        integer       :: ndfl, ndfl_sh, idfl, jdfl, nkeep_dfl, kfr_dfl(2)
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
        integer,             allocatable :: eo(:)
        real, pointer :: rmatp(:,:,:)
        real     :: fc
        real(dp) :: a, aa, e_mm, myv, mu_q, sd_q
        integer  :: it, q, r, i, ithr, nthr, batchlims(2), batchsz, ibatch, row, d_new, filtsz, sh
        !> effective (global) iteration numbering -- the ONLY counters iteration-keyed schedules
        !! and iteration logs may use; equal to it/niters except on a distributed probe worker
        integer  :: it_eff, niters_eff
        integer  :: nparts_sub
        logical  :: l_probe_distr
        complex,             allocatable :: cme(:,:,:,:), cmo(:,:,:,:)
        real,                allocatable :: rhe(:,:,:,:), rhoo(:,:,:,:)
        real(dp),            allocatable :: Mconv(:,:), sconv(:)
        real(dp) :: cos_mean
        real,     allocatable :: fscq_dg(:,:)
        real(dp) :: fmean_dg(512), fbest_dg
        real     :: res_dg
        integer  :: khi_dg, nsig_dg, ntop_dg, sel_dg(4), tq_dg, bq_dg
        real(dp) :: qml, nll_mix_add
        ! ---- MCFA state: tied-covariance mixture prior over the latents ----
        integer  :: kk2
        real(dp) :: lwm, wsm
        ! unified-mixture state: frame rotation (Gap B), full-N running average of the
        ! reduced mixture statistics (Gap A), starved-component reseeding, final full pass
        type(string) :: fname
        integer(timer_int_kind) :: t_it, t_sec
        real(timer_int_kind) :: sec_read, sec_prep, sec_estep, sec_ins
        real(dp) :: twp0, twp1, twp2
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
        !> the hoisted per-fit state (M2): every cross-iteration and per-iteration-per-fit
        !! variable of this routine now lives in one probe_fit_t owned by the driver
        type(probe_fit_t) :: fit
        !> cross-fit-FSC driver context (inert unless a SIMPLE_COV_XFSC_* gate is set)
        type(xfsc_ctx_t)  :: xfctx
        nthr = omp_get_max_threads()
        ! ---- M2 state hoist: populate the fit object from the arguments ----
        block
            character(len=:), allocatable :: pfx_l, meta_l
            pfx_l  = 'flex_pca_pc'
            meta_l = COV_PROBE_META
            if( present(fprefix) )    pfx_l  = trim(fprefix)
            if( present(meta_fname) ) meta_l = trim(meta_fname)
            call new_probe_fit(fit, 0, pfx_l, meta_l, pinds, nptcls)
        end block
        call move_alloc(basis_recs, fit%basis_recs)
        call move_alloc(eigvals,    fit%eigvals)
        fit%ncomp    = ncomp
        fit%sig2_eff = sig2_eff
        fit%sig2     = max(sig2_eff, DTINY)
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
        ! Only a WORKER divides the total by nparts: it holds one fromp/top partition. The master
        ! holds every particle, so passing nparts there divides twice and inflates the stride by
        ! exactly nparts. See cov_stage_subsample, which the initialiser shares.
        nparts_sub = 1
        if( rounds%is_worker() ) nparts_sub = params%nparts
        call cov_stage_subsample(build, fit%pinds, fit%nptcls, nparts_sub, COV_PROBE_MAX_PTCLS, &
            &'SIMPLE_COV_PROBE_MAX', 'PROBE', fit%ppinds, fit%npp)
        allocate(fit%z(fit%npp,fit%ncomp))
        ! ---- banked-adjoint setup: pose geometry is fixed across EM iterations, so the bank,
        ! the direction assignment and the direction sort are built ONCE per stage. Sorting
        ! per-fit stage configuration: every SIMPLE_COV_* selector that becomes fit state
        call fit_stage_config(params, fit, nthr)
        ! ---- CROSS-FIT-FSC setup: the ridge defaults OFF; the single-fit engine writes no records ----
        call xfsc_setup(xfctx, params, fit%kfr_ann, .false., .not. rounds%is_worker())
        l_probe_distr_pre = rounds%distributed()
        ! ---- DEVICE polar E-step (stage 2), SIMPLE_COV_POLAR_GPU=1: opt-in, composes with the
        ! production env set; falls back to the CPU polar path (same answer) when the card or
        ! the CUDA build is absent, or when the basis exceeds the resident-volume cap below.
        l_pol_gpu = .false.
        if( fit%l_pol_es )then
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
            &.not. (rounds%distributed()) )then
            l_bank_adj = .true.
            ndir_pb = cov_polar_ndir(fit%npp)
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
            ! rec_backend=pcg accumulates its pair kernels on the CPU insert path only
            if( trim(params%rec_backend) == 'pcg' )then
                l_gpu_estep = .false.
                l_bank_fwd  = .false.
            endif
            ! the polar former supplies G, b, c and the contrast itself, so the device E-step and
            ! the banked forward must both stand down or their statistics would be the ones used
            if( fit%l_pol_es )then
                l_gpu_estep = .false.
                l_bank_fwd  = .false.
            endif
            ! the banked forward still lacks a mixture solve; the fused device E-step COMPOSES
            ! with the mixture since 2026-08-19 -- the device fetches the same G/b/c sufficient
            ! statistics and the host branches to probe_solve_mix (one shared solver, both bodies)
            if( fit%l_mix_req )then
                l_bank_fwd  = .false.
            endif
            ! CAPACITY GATE. The device holds mean+basis as ONE resident allocation and the CUDA
            ! entry point refuses ncomp+1 > COV_GPU_ESTEP_MAXVOLS. Without this check the worker
            ! dies inside flex_gpu_estep_vols_f (THROW_HARD) and, because a dead worker never writes
            ! JOB_FINISHED, the distributed MASTER then waits forever -- an 8-hour silent hang
            ! observed at neigs=40. Fall back to the CPU E-step instead: slower, same answer.
            if( l_gpu_estep .and. (fit%ncomp + 1) > COV_GPU_ESTEP_MAXVOLS )then
                l_gpu_estep = .false.
                write(logfhandle,'(A,I0,A,I0,A)') '>>> FLEX_PCA PROBE fused device E-STEP OFF: needs ', &
                    &fit%ncomp+1,' resident volumes, device entry point caps at ',COV_GPU_ESTEP_MAXVOLS, &
                    &' -- using the CPU E-step for this basis size'
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
            l_deflate = .true.
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
            allocate(rmat_pp(3,3,fit%npp), nrm_pp(3,fit%npp), dirof_pb(fit%npp), cap(fit%npp), sap(fit%npp))
            do i = 1, fit%npp
                call build%spproj_field%get_ori(fit%ppinds(i), o_pb)
                rmat_pp(:,:,i) = o_pb%get_mat()
                nrm_pp(:,i)    = rmat_pp(3,:,i)
            end do
            call o_pb%kill
            call polar_assign_directions(nrm_pp, fit%npp, nrm_pb, ndir_pb, dirof_pb)
            do i = 1, fit%npp
                call polar_relative_inplane(rmat_pp(:,:,i), rmat_pb(:,:,dirof_pb(i)), ca1, sa1)
                cap(i) = ca1; sap(i) = sa1
            end do
            deallocate(rmat_pp, nrm_pp)
            ! stable counting sort by bank direction: contiguous same-direction runs inside each
            ! batch are what turn per-particle splats into per-direction splats
            allocate(cnt_pb(ndir_pb), source=0)
            do i = 1, fit%npp
                cnt_pb(dirof_pb(i)) = cnt_pb(dirof_pb(i)) + 1
            end do
            allocate(perm_pb(fit%npp))
            block
                integer, allocatable :: pos(:)
                allocate(pos(ndir_pb)); pos(1) = 1
                do jdir = 2, ndir_pb
                    pos(jdir) = pos(jdir-1) + cnt_pb(jdir-1)
                end do
                do i = 1, fit%npp
                    perm_pb(pos(dirof_pb(i))) = i
                    pos(dirof_pb(i)) = pos(dirof_pb(i)) + 1
                end do
                deallocate(pos)
            end block
            allocate(pptmp(fit%npp))
            pptmp = fit%ppinds;   fit%ppinds   = pptmp(perm_pb)
            pptmp = dirof_pb; dirof_pb = pptmp(perm_pb)
            deallocate(pptmp)
            block
                real, allocatable :: rtmp(:)
                allocate(rtmp(fit%npp))
                rtmp = cap; cap = rtmp(perm_pb)
                rtmp = sap; sap = rtmp(perm_pb)
                deallocate(rtmp)
            end block
            deallocate(perm_pb, cnt_pb)
            write(logfhandle,'(A,I0,A,I0,A)') '>>> FLEX_PCA PROBE BANKED ADJOINT ON (', &
                &ndir_pb,' bank directions, ',fit%npp,' records)'
            call flush(logfhandle)
        endif
        ! ---- DOWNSCALED-PARTICLE CACHE ----
        ! Every EM iteration re-reads and re-preps the SAME particles, and one qsys round is
        ! launched per iteration, so the workers are fresh processes each time and nothing held in
        ! memory survives them. The on-disk cache does survive: it stores the iteration-independent
        ! prefix (noise normalisation, FFT, crop to box_crop), which at box=360/box_crop=64 is where
        ! the measured 48% of probe time goes. Masking and the device prep both keep the full box.
        ! A cache-served batch is read at box_crop into cropped planes (init_rec cropped=,
        ! prepimgbatch at box_crop) and prepped with cached=.true. (no second noise renorm).
        l_pcache = plane_cache_in_use(params, build)
        if( l_pcache )then
            write(logfhandle,'(A)') '>>> FLEX_PCA PROBE: particles served from the downscaled cache'
            call flush(logfhandle)
        endif
        ! ---- EFFECTIVE (GLOBAL) ITERATION NUMBERING ----
        ! A distributed probe worker is relaunched once per master EM iteration with niters=1,
        ! so the local counter it is pinned at 1 in every round. Every iteration-keyed schedule
        ! (mixture warm-up, one-time checks) and every iteration log must key off
        ! it_eff/niters_eff, never it/niters (a schedule keyed off the local counter once left
        ! ~24% of every shard out of the M-step under nparts>1). The master stamps
        ! the true iteration through the probe-state file; shared memory reduces to it_eff == it.
        niters_eff = niters
        if( present(niters_glob) )then
            if( niters_glob > 0 ) niters_eff = niters_glob
        endif
        do it = 1, niters
            it_eff = it
            if( present(it_glob) )then
                if( it_glob > 0 ) it_eff = it_glob
            endif
            t_it = tic()
            sec_read = 0.; sec_prep = 0.; sec_estep = 0.; sec_ins = 0.
            call fit_iter_begin(params, build, fit, mean_rec, it_eff, niters_eff, nthr)
            allocate(orientations(MAXIMGBATCHSZ), eo(MAXIMGBATCHSZ))
            if( l_bank_fwd ) allocate(basis_fplsC(fit%ncomp,nthr), mean_fplC(nthr))
            call init_rec(params, build, MAXIMGBATCHSZ, fpls, cropped=l_pcache)
            ! ---- ACCUMULATE: distributed when the master has parts, in-process otherwise ----
            ! One qsys round per EM iteration, because the basis the E-step projects changes every
            ! iteration: workers are relaunched against the master's refreshed flex_pca_pc*.mrc
            ! rather than looping locally. Everything below the reduction -- the coupled M-step
            ! solve, the FSC merge, the re-orthonormalisation -- stays master-only and unchanged, so
            ! the shared-memory result remains the reference.
            l_probe_distr = rounds%distributed()
            if( l_probe_distr )then
                ! the round's control state (iteration, budget, one fit) travels in job_descr; the
                ! stage selects the basis namespace a relaunched worker loads -- the joint fit
                ! after the merge refines the polished basis, not fit A's
                call save_probe_state(fit%ncomp, fit%eigvals, fit%sig2_eff)
                if( fit%fprefix%to_char() == 'flex_pca_polished_pc' )then
                    call rounds%run_stage(params, PCA_STAGE_POLISH, 'joint-fit iteration', &
                        &which_iter=it_eff, maxits=niters_eff, nfits=1)
                else
                    call rounds%run_stage(params, PCA_STAGE_PROBE, 'probe iteration', &
                        &which_iter=it_eff, maxits=niters_eff, nfits=1)
                endif
                allocate(cme(fit%es(1),fit%es(2),fit%es(3),fit%ncomp), cmo(fit%es(1),fit%es(2),fit%es(3),fit%ncomp), source=(0.,0.))
                allocate(rhe(fit%es(1),fit%es(2),fit%es(3),fit%ncomp), rhoo(fit%es(1),fit%es(2),fit%es(3),fit%ncomp), source=0.)
                if( fit%l_mix_req )then
                    if( allocated(fit%dm_sm) )then
                        if( size(fit%dm_sm,1) /= fit%ncomp .or. size(fit%dm_sr) /= fit%kmix ) &
                            &deallocate(fit%dm_sr, fit%dm_sm, fit%dm_smm, fit%dm_sai, fit%dm_z)
                    endif
                    if( .not. allocated(fit%dm_sr) )then
                        allocate(fit%dm_sr(fit%kmix), fit%dm_sm(fit%ncomp,fit%kmix), &
                            &fit%dm_smm(fit%ncomp,fit%ncomp,fit%kmix), fit%dm_sai(fit%ncomp,fit%ncomp))
                        allocate(fit%dm_z(MIX_ZSUB_MAX*rounds%nparts(), fit%ncomp))
                    endif
                    fit%dm_sr = 0.d0; fit%dm_sm = 0.d0; fit%dm_smm = 0.d0; fit%dm_sai = 0.d0; fit%dm_nz = 0
                    call reduce_probe_parts(params, rounds%nparts(), cme, rhe, cmo, rhoo, &
                        &fit%rho_e, fit%rho_o, fit%gam_sum, fit%nll_tot, fit%nval, fit%ncomp, fit%kpk_e, fit%kpk_o, fit%rpk_e, fit%rpk_o, &
                        &fit%dm_sr, fit%dm_sm, fit%dm_smm, fit%dm_sai, fit%dm_z, fit%dm_nz)
                    write(logfhandle,'(A,I0,A,I0)') '>>> FLEX_PCA MIX distributed reduce: kmix=', &
                        &fit%kmix,'  pooled latent subsample rows=',fit%dm_nz
                    call flush(logfhandle)
                else
                    call reduce_probe_parts(params, rounds%nparts(), cme, rhe, cmo, rhoo, &
                        &fit%rho_e, fit%rho_o, fit%gam_sum, fit%nll_tot, fit%nval, fit%ncomp, fit%kpk_e, fit%kpk_o, fit%rpk_e, fit%rpk_o)
                endif
                do q = 1, fit%ncomp
                    fit%Yeven(q)%cmat_exp = cme(:,:,:,q); fit%Yeven(q)%rho_exp = rhe(:,:,:,q)
                    fit%Yodd(q)%cmat_exp  = cmo(:,:,:,q); fit%Yodd(q)%rho_exp  = rhoo(:,:,:,q)
                end do
                deallocate(cme, cmo, rhe, rhoo)
            else
            if( l_pcache )then
                call prepimgbatch(params, build, MAXIMGBATCHSZ, box=params%box_crop, smpd=params%smpd_crop)
            else
                call prepimgbatch(params, build, MAXIMGBATCHSZ)
            endif
            ! default-ON whenever the CUDA build sees a device (measured 1.2 s/iter vs 11-12 s
            ! CPU insert); SIMPLE_COV_GPU=0 opts out. cov_env_int cannot express "off" (=0 is
            ! indistinguishable from unset), so the opt-out goes through the presence-and-zero
            ! reader. The device E-step and the bank adjoint stay strictly opt-in above.
            l_gpu_probe  = flex_gpu_available() .and. .not. cov_env_int_off('SIMPLE_COV_GPU')
            if( l_gpu_probe )then
                write(logfhandle,'(A,I0,A,I0,A)') '>>> FLEX_PCA PROBE M-step GPU insertion ON (', &
                    &2*fit%ncomp, ' components, ', 2*fit%npairs, ' density rows)'
                call flush(logfhandle)
                call flex_gpu_coupled_begin_f(fit%Yeven, 2*fit%ncomp, 2*fit%npairs)
                allocate(gdsc(2*fit%ncomp,MAXIMGBATCHSZ), grsc(2*fit%npairs,MAXIMGBATCHSZ), source=0.d0)
                if( l_bank_adj ) call flex_gpu_coupled_bank_f(rmat_pb, ndir_pb)
                if( l_gpu_estep .and. l_bank_adj )then
                    call flex_gpu_estep_vols_f(mean_rec, fit%basis_recs, fit%ncomp)
                    nst_es = 2 + 2*fit%ncomp + (fit%ncomp*(fit%ncomp+1))/2
                    allocate(stats_es(nst_es, MAXIMGBATCHSZ), source=0.d0)
                    write(logfhandle,'(A,I0,A)') '>>> FLEX_PCA PROBE FUSED E-STEP ON (', &
                        &fit%ncomp+1, ' volumes resident, exact directions)'
                    call flush(logfhandle)
                endif
            endif
            fit%z = 0.d0
            do ibatch = 1, fit%npp, MAXIMGBATCHSZ
                batchlims = [ibatch, min(fit%npp, ibatch + MAXIMGBATCHSZ - 1)]
                batchsz   = batchlims(2) - batchlims(1) + 1
                if( l_gpu_prep )then
                    t_sec = tic()
                    if( l_pcache )then
                        call plane_cache_read_batch(params, fit%npp, fit%ppinds, batchlims)
                    else
                        call discrete_read_imgbatch(params, build, fit%npp, fit%ppinds, batchlims)
                    endif
                    sec_read = sec_read + toc(t_sec)
                else
                    call planes_batch_load(params, build, fit%npp, fit%ppinds, batchlims, fpls, &
                        &cov_image_mask_radius(params), l_pcache, sec_read, sec_prep)
                endif
                do i = 1, batchsz
                    call build%spproj_field%get_ori(fit%ppinds(batchlims(1)+i-1), orientations(i))
                    eo(i) = build%spproj_field%get_eo(fit%ppinds(batchlims(1)+i-1))
                    if( l_gpu_prep )then
                        ctfp_arr(i)  = build%spproj%get_ctfparams(params%oritype, fit%ppinds(batchlims(1)+i-1))
                        shf_pb(:,i)  = build%spproj_field%get_2Dshift(fit%ppinds(batchlims(1)+i-1))
                        if( params%l_ml_reg ) call resample_sigma2(kf_pb(1), signyq_pb, &
                            &build%esig%sigma2_noise(kf_pb(1):kf_pb(2), fit%ppinds(batchlims(1)+i-1)), &
                            &nyqfull_pb, real(signyq_pb)/real(nyqfull_pb), sig2_pb(:,i))
                    endif
                end do
                ! ---- POLAR E-STEP BANK, built once per EM iteration at the first prepped batch
                ! (basis_recs change every M-step, so the bank cannot be cached across iterations;
                ! the grid and the pose-fixed direction assignment are built once per stage) ----
                if( fit%l_pol_es .and. .not. fit%l_pol_bank_it )then
                    t_bank = tic()
                    call fit_polar_bank_build(params, build, fit, mean_rec, fpls(1), nthr, l_pol_gpu)
                    ! GPU parity arrays (gate 1): the CPU polar former re-run on the check
                    ! particles; the exact-direction/DC-less attribution arms stand down under
                    ! the device path (they answer the polar-vs-Cartesian question, already
                    ! measured on this recipe; the device gate is GPU-vs-CPU-polar)
                    ! exact-direction check arm scratch (grid-shaped, so allocated here)
                    ! ---- device upload (stage 2): the bank + ring tables once per iteration, and
                    ! the mean+basis volumes for the hybrid-exact kernel. The rank can change
                    ! between iterations, so the capacity gate is re-checked here; falling back
                    ! mid-run is seamless because every consumer below branches per iteration.
                    if( l_pol_gpu .and. fit%l_pol_hyb .and. (fit%ncomp + 1) > COV_GPU_ESTEP_MAXVOLS )then
                        l_pol_gpu = .false.
                        write(logfhandle,'(A,I0,A,I0,A)') '>>> FLEX_PCA POLAR E-STEP DEVICE OFF: needs ', &
                            &fit%ncomp+1,' resident volumes, device entry point caps at ', &
                            &COV_GPU_ESTEP_MAXVOLS,' -- using the CPU polar E-step'
                        call flush(logfhandle)
                    endif
                    if( l_pol_gpu )then
                        call flex_gpu_poles_bank_f(fit%UsallE, fit%CfE, fit%Cm0E, fit%c00E, fit%ncomp, fit%ndir_es)
                        if( fit%l_pol_hyb ) call flex_gpu_estep_vols_f(mean_rec, fit%basis_recs, fit%ncomp)
                        nst_pg = 2 + 2*fit%ncomp + (fit%ncomp*(fit%ncomp+1))/2
                        if( allocated(stats_pg) )then
                            if( size(stats_pg,1) /= nst_pg ) deallocate(stats_pg)
                        endif
                        if( .not. allocated(stats_pg) ) allocate(stats_pg(nst_pg, MAXIMGBATCHSZ), &
                            &source=0.d0)
                        if( .not. allocated(vld_pg) ) allocate(vld_pg(MAXIMGBATCHSZ), &
                            &dir_pg(MAXIMGBATCHSZ), ca_pg(MAXIMGBATCHSZ), sa_pg(MAXIMGBATCHSZ))
                    endif
                    fit%sec_bank      = fit%sec_bank + toc(t_bank)
                    fit%l_pol_bank_it = .true.
                    write(logfhandle,'(A,I0,A,I0,A,I0,A,F7.1)') '>>> FLEX_PCA POLAR ESTEP BANK it=', &
                        &it_eff,'  directions built=',count(fit%dused_es),' of ',fit%ndir_es, &
                        &'  build seconds=',fit%sec_bank
                    call flush(logfhandle)
                endif
                fit%valid(:batchsz) = .false.
                fit%zbatch(:,:batchsz) = 0.d0
                fit%dens(:,:,:batchsz) = 0.d0
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
                        if( it_eff == 1 .and. ibatch == 1 )then
                            vprep = 0
                            call cov_env_int('SIMPLE_COV_PREP_CHECK', vprep)
                            if( vprep > 0 )then
                                call prep_imgs4projected_model(params, build, batchsz, &
                                    &build%imgbatch(:batchsz), fit%ppinds(batchlims(1):batchlims(2)), &
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
                    fit%sec_proj_thr(1) = fit%sec_proj_thr(1) + (omp_get_wtime() - twp0)
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
                        do q = 1, fit%ncomp
                            fit%bth(q,ithr) = stats_es(2+q,i)
                            fit%cth(q,ithr) = stats_es(2+fit%ncomp+q,i)
                        end do
                        off_es = 2 + 2*fit%ncomp
                        do q = 1, fit%ncomp
                            do r = q, fit%ncomp
                                off_es = off_es + 1
                                fit%Gth(q,r,ithr) = stats_es(off_es,i)
                                fit%Gth(r,q,ithr) = fit%Gth(q,r,ithr)
                            end do
                        end do
                        if( fit%l_mix_used )then
                            ! same shared MCFA solver as the plain body: the device fetched the
                            ! identical G/b/c sufficient statistics, only the host solve differs
                            call probe_solve_mix(fit%ncomp, fit%kmix, fit%Gth(:,:,ithr), fit%bth(:,ithr), &
                                &fit%cth(:,ithr), myv, e_mm, fit%n_probe_cm, a, fit%sig2, fit%mix_Ominv, fit%mix_Omxi, &
                                &fit%mix_lpi, fit%mix_xiOx, fit%zth(:,ithr), fit%Ainvth(:,:,ithr), fit%dens(:,:,i), &
                                &ldA, lok, nll_mix_add, fit%mxa_sr(:,ithr), fit%mxa_sm(:,:,ithr), &
                                &fit%mxa_smm(:,:,:,ithr), fit%mxa_sainv(:,:,ithr))
                            if( lok ) fit%nll_thr(ithr) = fit%nll_thr(ithr) + nll_mix_add
                        else
                            call probe_solve_ecm(fit%ncomp, fit%Gth(:,:,ithr), fit%bth(:,ithr), fit%cth(:,ithr), &
                                &myv, e_mm, fit%prior, fit%sig2, fit%n_probe_cm, a, fit%zth(:,ithr), &
                                &fit%Ainvth(:,:,ithr), ldA, lok, qml)
                            if( lok ) fit%nll_thr(ithr) = fit%nll_thr(ithr) + ldA - qml
                        endif
                        aa = a*a
                        fit%z(row,:)    = fit%zth(:,ithr)
                        fit%zbatch(:,i) = fit%zth(:,ithr)
                        if( .not. fit%l_mix_used )then
                            do r = 1, fit%ncomp
                                do q = 1, fit%ncomp
                                    fit%dens(q,r,i) = fit%zth(q,ithr)*fit%zth(r,ithr) + fit%Ainvth(q,r,ithr)
                                end do
                            end do
                        endif
                        do q = 1, fit%ncomp
                            fit%gam_thr(q,ithr) = fit%gam_thr(q,ithr) + fit%dens(q,q,i)
                        end do
                        ! Gamma above reads the UNSCALED E[zz'] (prior on z, not on a z);
                        ! everything below feeds the M-step, which under the fitted scale
                        ! solves  r ~ a B z + n
                        if( fit%l_probe_mls )then
                            fit%zbatch(:,i) = a*fit%zbatch(:,i)
                            fit%dens(:,:,i) = aa*fit%dens(:,:,i)
                        endif
                        fit%nval_thr(ithr) = fit%nval_thr(ithr) + 1
                        fit%valid(i)       = .true.
                        avec_es(i)     = real(a)
                    end do
                    !$omp end parallel do
                    ! residual stays on device; the M-step consumes it via the resident entry
                    if( l_deflate )then
                        call flex_gpu_estep_resid_f(fpls(:batchsz), avec_es(:batchsz), &
                            &vld_es(:batchsz), batchsz, fetch=.false.)
                    endif
                    fit%sec_gram_thr(1) = fit%sec_gram_thr(1) + (omp_get_wtime() - twp0)
                elseif( l_pol_gpu )then
                    ! ---- POLAR E-STEP ON DEVICE (stage 2): per-batch ring + hybrid-exact
                    ! statistics in the fused stat layout; everything downstream of the
                    ! statistics -- solve, mixture, Gamma, windows, e/o, mean projection,
                    ! residual, insertion -- is the CPU polar path's, verbatim.
                    twp0 = omp_get_wtime()
                    do i = 1, batchsz
                        row = batchlims(1) + i - 1
                        vld_pg(i) = .not. orientations(i)%isstatezero()
                        dir_pg(i) = fit%dir_es(row)
                        ca_pg(i)  = fit%cae(row)
                        sa_pg(i)  = fit%sae(row)
                    end do
                    call flex_gpu_poles_batch_f(fpls(:batchsz), orientations(:batchsz), &
                        &ca_pg(:batchsz), sa_pg(:batchsz), dir_pg(:batchsz), vld_pg(:batchsz), &
                        &batchsz, merge(fit%rhyb_es*(fit%rhyb_es+1), 0, fit%l_pol_hyb), stats_pg(:,:batchsz))
                    fit%sec_proj_thr(1) = fit%sec_proj_thr(1) + (omp_get_wtime() - twp0)
                    twp0 = omp_get_wtime()
                    !$omp parallel do default(shared) schedule(dynamic) proc_bind(close) &
                    !$omp& private(i,ithr,row,q,r,off_es,a,aa,e_mm,myv,ldA,lok,qml,nll_mix_add) &
                    !$omp& private(idp_es,taz_es)
                    do i = 1, batchsz
                        if( .not. vld_pg(i) ) cycle
                        ithr = omp_get_thread_num() + 1
                        row  = batchlims(1) + i - 1
                        e_mm = stats_pg(1,i)
                        myv  = stats_pg(2,i)
                        a    = max(0.1d0, min(5.0d0, myv / max(e_mm, DTINY)))
                        aa   = a*a
                        do q = 1, fit%ncomp
                            fit%bth(q,ithr) = stats_pg(2+q,i)
                            fit%cth(q,ithr) = stats_pg(2+fit%ncomp+q,i)
                        end do
                        off_es = 2 + 2*fit%ncomp
                        do q = 1, fit%ncomp
                            do r = q, fit%ncomp
                                off_es = off_es + 1
                                fit%Gth(q,r,ithr) = stats_pg(off_es,i)
                                fit%Gth(r,q,ithr) = fit%Gth(q,r,ithr)
                            end do
                        end do
                        if( fit%l_mix_used )then
                            ! same shared MCFA solver as the CPU bodies: the device fetched the
                            ! identical G/b/c sufficient statistics, only the former differs
                            call probe_solve_mix(fit%ncomp, fit%kmix, fit%Gth(:,:,ithr), fit%bth(:,ithr), &
                                &fit%cth(:,ithr), myv, e_mm, fit%n_probe_cm, a, fit%sig2, fit%mix_Ominv, fit%mix_Omxi, &
                                &fit%mix_lpi, fit%mix_xiOx, fit%zth(:,ithr), fit%Ainvth(:,:,ithr), fit%dens(:,:,i), &
                                &ldA, lok, nll_mix_add, fit%mxa_sr(:,ithr), fit%mxa_sm(:,:,ithr), &
                                &fit%mxa_smm(:,:,:,ithr), fit%mxa_sainv(:,:,ithr))
                            if( lok ) fit%nll_thr(ithr) = fit%nll_thr(ithr) + nll_mix_add
                        else
                            call probe_solve_ecm(fit%ncomp, fit%Gth(:,:,ithr), fit%bth(:,ithr), fit%cth(:,ithr), &
                                &myv, e_mm, fit%prior, fit%sig2, fit%n_probe_cm, a, fit%zth(:,ithr), &
                                &fit%Ainvth(:,:,ithr), ldA, lok, qml)
                            if( lok ) fit%nll_thr(ithr) = fit%nll_thr(ithr) + ldA - qml
                        endif
                        aa = a*a
                        fit%z(row,:)    = fit%zth(:,ithr)
                        fit%zbatch(:,i) = fit%zth(:,ithr)
                        ! parity check: the device-path z actually used, after the shared solve
                        if( .not. fit%l_mix_used )then
                            do r = 1, fit%ncomp
                                do q = 1, fit%ncomp
                                    fit%dens(q,r,i) = fit%zth(q,ithr)*fit%zth(r,ithr) + fit%Ainvth(q,r,ithr)
                                end do
                            end do
                        endif
                        do q = 1, fit%ncomp
                            fit%gam_thr(q,ithr) = fit%gam_thr(q,ithr) + fit%dens(q,q,i)
                        end do
                        if( fit%l_probe_mls )then
                            fit%zbatch(:,i) = a*fit%zbatch(:,i)
                            fit%dens(:,:,i) = aa*fit%dens(:,:,i)
                        endif
                        fit%nval_thr(ithr) = fit%nval_thr(ithr) + 1
                        fit%valid(i)       = .true.
                        fit%gam_dbg(1,ithr) = fit%gam_dbg(1,ithr) + sum([(fit%Gth(q,q,ithr), q=1,fit%ncomp)])
                        fit%gam_dbg(2,ithr) = fit%gam_dbg(2,ithr) + dot_product(fit%bth(:,ithr), fit%bth(:,ithr))
                        fit%gam_dbg(3,ithr) = fit%gam_dbg(3,ithr) + dot_product(fit%cth(:,ithr), fit%cth(:,ithr))
                        fit%gam_dbg(4,ithr) = fit%gam_dbg(4,ithr) + a
                        ! residual r_i = y - a*(T mu) in place, banded, exactly the CPU polar
                        ! path's (the check arm's Cartesian mean plane serves where present)
                        call project_fplane_mean_banded(mean_rec, orientations(i), fpls(i), &
                            &fit%mean_fpl(ithr))
                        call subtract_mean_banded(fpls(i), fit%mean_fpl(ithr), real(a), fit%nyqr_es)
                    end do
                    !$omp end parallel do
                    fit%sec_gram_thr(1) = fit%sec_gram_thr(1) + (omp_get_wtime() - twp0)
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
                    call project_fplanes_mean_basis(mean_rec, fit%basis_recs, o_pb_thr(ithr), &
                        &fpls(seg_beg_b(iseg)), fit%mean_fpl(ithr), fit%basis_fpls(:,ithr), &
                        &apply_ctf_amp=.false.)
                    twp1 = omp_get_wtime()
                    fit%sec_proj_thr(ithr) = fit%sec_proj_thr(ithr) + (twp1 - twp0)
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
                        call ensure_bank_ctf_planes(fpls(i), mean_fplC(ithr), basis_fplsC(:,ithr), fit%ncomp)
                        mean_fplC(ithr)%cmplx_plane = fpls(i)%transfer_plane * fit%mean_fpl(ithr)%cmplx_plane
                        do q = 1, fit%ncomp
                            basis_fplsC(q,ithr)%cmplx_plane = fpls(i)%transfer_plane * &
                                &fit%basis_fpls(q,ithr)%cmplx_plane
                        end do
                        e_mm = real(cov_herm_inner(mean_fplC(ithr), mean_fplC(ithr)), dp)
                        myv  = real(cov_herm_inner(mean_fplC(ithr), fpls(i)), dp)
                        a    = max(0.1d0, min(5.0d0, myv / max(e_mm, DTINY)))
                        aa   = a*a
                        do q = 1, fit%ncomp
                            fit%bth(q,ithr) = real(cov_herm_inner(basis_fplsC(q,ithr), fpls(i)), dp)
                            fit%cth(q,ithr) = real(cov_herm_inner(basis_fplsC(q,ithr), mean_fplC(ithr)), dp)
                            do r = q, fit%ncomp
                                fit%Gth(q,r,ithr) = real(cov_herm_inner(basis_fplsC(q,ithr), &
                                    &basis_fplsC(r,ithr)), dp)
                                fit%Gth(r,q,ithr) = fit%Gth(q,r,ithr)
                            end do
                        end do
                        call probe_solve_ecm(fit%ncomp, fit%Gth(:,:,ithr), fit%bth(:,ithr), fit%cth(:,ithr), &
                            &myv, e_mm, fit%prior, fit%sig2, fit%n_probe_cm, a, fit%zth(:,ithr), &
                            &fit%Ainvth(:,:,ithr), ldA, lok, qml)
                        aa = a*a
                        if( lok ) fit%nll_thr(ithr) = fit%nll_thr(ithr) + ldA - qml
                        fit%z(row,:)          = fit%zth(:,ithr)
                        fit%zbatch(:,i)       = fit%zth(:,ithr)
                        do r = 1, fit%ncomp
                            do q = 1, fit%ncomp
                                fit%dens(q,r,i) = fit%zth(q,ithr)*fit%zth(r,ithr) + fit%Ainvth(q,r,ithr)
                            end do
                        end do
                        do q = 1, fit%ncomp
                            fit%gam_thr(q,ithr) = fit%gam_thr(q,ithr) + fit%dens(q,q,i)
                        end do
                        if( fit%l_probe_mls )then
                            fit%zbatch(:,i) = a*fit%zbatch(:,i)
                            fit%dens(:,:,i) = aa*fit%dens(:,:,i)
                        endif
                        fit%nval_thr(ithr)    = fit%nval_thr(ithr) + 1
                        fit%valid(i)          = .true.
                        ! residual in the BANK frame (transfer/ctfsq already aligned above)
                        fpls(i)%cmplx_plane = fpls(i)%cmplx_plane - real(a)*mean_fplC(ithr)%cmplx_plane
                        twp2 = omp_get_wtime()
                        fit%sec_gram_thr(ithr) = fit%sec_gram_thr(ithr) + (twp2 - twp1)
                    end do
                end do
                !$omp end parallel do
                else
                !$omp parallel do default(shared) schedule(dynamic) proc_bind(close) &
                !$omp& private(i,ithr,q,r,a,aa,e_mm,myv,row,twp0,twp1,twp2,ldA,lok,qml,kk2,lwm,wsm) &
                !$omp& private(idp_es,pw_es,cnt_es,taz_es,ir)
                do i = 1, batchsz
                    if( orientations(i)%isstatezero() ) cycle
                    ithr = omp_get_thread_num() + 1
                    row  = batchlims(1) + i - 1
                    twp0 = omp_get_wtime()
                    if( fit%l_pol_es )then
                        call fit_estep_former_polar(fit, mean_rec, orientations(i), fpls(i), row, ithr, a, aa, e_mm, myv)
                    else
                        call fit_estep_former_cart(fit, mean_rec, orientations(i), fpls(i), ithr, a, aa, e_mm, myv)
                    endif
                    call fit_estep_solve_stats(fit, fpls(i), i, row, ithr, a, aa, e_mm, myv)
                end do
                !$omp end parallel do
                endif
                sec_estep = sec_estep + toc(t_sec)
                t_sec = tic()
                ! M-step by halfset: Y_q += sum_i z_iq * backproject(r_i), and the coupled normal matrix
                ! rho(q,r) += sum_i |CTF|^2 E[z_iq z_ir]   (batched KB)
                if( l_gpu_probe )then
                    do i = 1, batchsz
                        fit%valid_e(i) = fit%valid(i) .and. eo(i)==0
                        fit%valid_o(i) = fit%valid(i) .and. eo(i)==1
                    end do
                    ! halfset routing through the slots: even records fill 1..ncomp / 1..npairs,
                    ! odd records the upper halves; density slots honor the diag/full mode
                    gdsc(:,:batchsz) = 0.d0
                    grsc(:,:batchsz) = 0.d0
                    do i = 1, batchsz
                        if( .not. (fit%valid_e(i) .or. fit%valid_o(i)) ) cycle
                        gq = merge(0, fit%ncomp,  fit%valid_e(i))
                        gr = merge(0, fit%npairs, fit%valid_e(i))
                        gdsc(gq+1:gq+fit%ncomp,i) = fit%zbatch(:,i)
                        do r = 1, fit%ncomp
                            do q = 1, r
                                grsc(gr+(r*(r-1))/2+q,i) = fit%dens(q,r,i)
                            end do
                        end do
                    end do
                    if( l_bank_adj )then
                        if( l_bank_fwd )then
                            ! residuals were built in the bank frame by the forward -- unit rotation
                            call flex_gpu_coupled_batch_banked_f(fpls(:batchsz), gdsc(:,:batchsz), &
                                &grsc(:,:batchsz), fit%valid_e(:batchsz) .or. fit%valid_o(:batchsz), batchsz, &
                                &dirof_pb(batchlims(1):batchlims(2)), cap1(:batchsz), sap0(:batchsz))
                        else if( l_gpu_estep )then
                            ! fused path: the residual is already resident on device
                            call flex_gpu_coupled_batch_banked_res_f(gdsc(:,:batchsz), &
                                &grsc(:,:batchsz), fit%valid_e(:batchsz) .or. fit%valid_o(:batchsz), batchsz, &
                                &dirof_pb(batchlims(1):batchlims(2)), cap(batchlims(1):batchlims(2)), &
                                &sap(batchlims(1):batchlims(2)))
                        else
                            call flex_gpu_coupled_batch_banked_f(fpls(:batchsz), gdsc(:,:batchsz), &
                                &grsc(:,:batchsz), fit%valid_e(:batchsz) .or. fit%valid_o(:batchsz), batchsz, &
                                &dirof_pb(batchlims(1):batchlims(2)), cap(batchlims(1):batchlims(2)), &
                                &sap(batchlims(1):batchlims(2)))
                        endif
                    else
                        call flex_gpu_coupled_batch_raw_f(build%pgrpsyms, orientations(:batchsz), &
                            &fpls(:batchsz), gdsc(:,:batchsz), grsc(:,:batchsz), &
                            &fit%valid_e(:batchsz) .or. fit%valid_o(:batchsz), batchsz)
                    endif
                else
                    call fit_batch_insert(build, fit, orientations, fpls, eo, batchsz)
                endif
                sec_ins = sec_ins + toc(t_sec)
                if( batchlims(2)==fit%nptcls .or. mod(batchlims(2), 5*MAXIMGBATCHSZ)==0 )then
                    write(logfhandle,'(A,I0,A,I0)') '>>> FLEX_PCA PROBE PASS PARTICLES: ',batchlims(2),' / ',fit%npp
                    call flush(logfhandle)
                endif
            end do
            if( l_gpu_probe )then
                call flex_gpu_coupled_end_f(fit%Yeven, fit%rho_e, fit%Yodd, fit%rho_o)
                deallocate(gdsc, grsc)
            endif
            write(logfhandle,'(A,F7.1,A,F7.1,A,F7.1,A,F7.1)') '>>> FLEX_PCA PROBE E-STEP SPLIT (seconds): read=', &
                &sec_read,'  prep=',sec_prep,'  project+solve=',sec_estep,'  insert=',sec_ins
            write(logfhandle,'(A,F7.1,A,F7.1)') '>>> FLEX_PCA PROBE E-STEP INNER (thread-seconds): project=', &
                &sum(fit%sec_proj_thr),'  gram+solve=',sum(fit%sec_gram_thr)
            call flush(logfhandle)
            if( fit%l_pol_es )then
                ! the two numbers the polar A/B reads from run.log
                write(logfhandle,'(A,I0,A,F7.1,A,F7.1)') '>>> FLEX_PCA POLAR ESTEP it=',it_eff, &
                    &'  bank build seconds=',fit%sec_bank,'  estep seconds=',sec_estep
                call flush(logfhandle)
            endif
            ! ---- POLAR CHECK SUMMARY (iteration 1): three formers ran on the first N particles:
            ! Cartesian (reference), polar at the bank direction (the production stage-1 path) and
            ! polar at the exact direction. banked-vs-exact isolates the direction quantization;
            ! exact-vs-Cartesian isolates the polar formulation (quadrature + radial |T|^2). ----
            call fit_iter_reduce(fit, it_eff, nthr, rounds=rounds)
            if( rounds%is_worker() )then
                allocate(cme(fit%es(1),fit%es(2),fit%es(3),fit%ncomp), cmo(fit%es(1),fit%es(2),fit%es(3),fit%ncomp))
                allocate(rhe(fit%es(1),fit%es(2),fit%es(3),fit%ncomp), rhoo(fit%es(1),fit%es(2),fit%es(3),fit%ncomp))
                do q = 1, fit%ncomp
                    cme(:,:,:,q) = fit%Yeven(q)%cmat_exp; rhe(:,:,:,q)  = fit%Yeven(q)%rho_exp
                    cmo(:,:,:,q) = fit%Yodd(q)%cmat_exp;  rhoo(:,:,:,q) = fit%Yodd(q)%rho_exp
                end do
                fit%nll_tot = sum(fit%nll_thr)
                if( fit%l_mix_req )then
                    block
                        real(dp), allocatable :: w_sr(:), w_sm(:,:), w_smm(:,:,:), w_sai(:,:), w_z(:,:)
                        integer :: tt2, kk4, nzs, izs, istep
                        allocate(w_sr(fit%kmix), w_sm(fit%ncomp,fit%kmix), &
                            &w_smm(fit%ncomp,fit%ncomp,fit%kmix), w_sai(fit%ncomp,fit%ncomp))
                        w_sr = 0.d0; w_sm = 0.d0; w_smm = 0.d0; w_sai = 0.d0
                        do tt2 = 1, nthr
                            w_sr  = w_sr  + fit%mxa_sr(:,tt2)
                            w_sm  = w_sm  + fit%mxa_sm(:,:,tt2)
                            w_sai = w_sai + fit%mxa_sainv(:,:,tt2)
                            do kk4 = 1, fit%kmix
                                w_smm(:,:,kk4) = w_smm(:,:,kk4) + fit%mxa_smm(:,:,kk4,tt2)
                            end do
                        end do
                        ! deterministic stride subsample of this part's latents (bounded)
                        nzs   = min(MIX_ZSUB_MAX, fit%npp)
                        istep = max(1, fit%npp / max(1,nzs))
                        nzs   = min(nzs, (fit%npp + istep - 1)/istep)
                        allocate(w_z(nzs,fit%ncomp))
                        do izs = 1, nzs
                            w_z(izs,:) = fit%z(min(fit%npp, 1 + (izs-1)*istep), :)
                        end do
                        call write_probe_part(flex_pca_part_fname('probe', params%part, params%numlen), &
                            &cme, rhe, cmo, rhoo, fit%rho_e, fit%rho_o, fit%gam_sum, fit%nll_tot, fit%nval, fit%ncomp, &
                            &fit%kpk_e, fit%kpk_o, fit%rpk_e, fit%rpk_o, w_sr, w_sm, w_smm, w_sai, w_z)
                        deallocate(w_sr, w_sm, w_smm, w_sai, w_z)
                    end block
                else
                    call write_probe_part(flex_pca_part_fname('probe', params%part, params%numlen), &
                        &cme, rhe, cmo, rhoo, fit%rho_e, fit%rho_o, fit%gam_sum, fit%nll_tot, fit%nval, fit%ncomp, &
                        &fit%kpk_e, fit%kpk_o, fit%rpk_e, fit%rpk_o)
                endif
                deallocate(cme, cmo, rhe, rhoo, fit%gam_sum)
                call cleanup_rec_buffers(build, fpls)
                do q = 1, size(fit%Yeven)
                    call fit%Yeven(q)%dealloc_rho; call fit%Yeven(q)%kill
                    call fit%Yodd(q)%dealloc_rho;  call fit%Yodd(q)%kill
                end do
                deallocate(fit%Yeven, fit%Yodd, fit%rho_e, fit%rho_o, fit%prior)
                if( allocated(fit%kacc_e) ) deallocate(fit%kacc_e)
                if( allocated(fit%kacc_o) ) deallocate(fit%kacc_o)
                if( allocated(fit%kpk_e) ) deallocate(fit%kpk_e)
                if( allocated(fit%kpk_o) ) deallocate(fit%kpk_o)
                if( allocated(fit%racc_e) ) deallocate(fit%racc_e)
                if( allocated(fit%racc_o) ) deallocate(fit%racc_o)
                if( allocated(fit%rpk_e) ) deallocate(fit%rpk_e)
                if( allocated(fit%rpk_o) ) deallocate(fit%rpk_o)
                deallocate(fit%Gth, fit%Ath, fit%bth, fit%cth, fit%zth, fit%basis_fpls, fit%mean_fpl, orientations, fit%zbatch, fit%dens)
                deallocate(fit%valid, fit%valid_e, fit%valid_o, eo, fit%Ainvth, fit%Acpth, fit%gam_thr, fit%gam_acc, fit%nval_thr)
                deallocate(fit%hth, fit%nll_thr)
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
                deallocate(fit%z, fit%ppinds)
                ! hoisted model handles back to the caller's arguments
                call move_alloc(fit%basis_recs, basis_recs)
                call move_alloc(fit%eigvals,    eigvals)
                ncomp = fit%ncomp
                return
            endif
            endif
            ! crossfsc per-iteration prep (master): set the harvest flag for the writer payloads
            ! and build the SSNR ridge from record t-1 (fit_iter_finish applies it just before
            ! the coupled solves; on scaffolding records the arms degrade to arm 0 and log it)
            if( xfctx%l_any ) call xfsc_prep_iter(xfctx, params, fit, it_eff, '')
            call fit_iter_finish(params, build, fit, it_eff, nthr)
            write(logfhandle,'(A,I0,A,ES12.4,A,ES12.4,A,F8.1)') '>>> FLEX_PCA PROBE ITER ',it_eff, &
                &' refined dim=',real(fit%ncomp),' max var=',maxval(fit%eigvals),' seconds=',toc(t_it)
            call flush(logfhandle)
            ! cleanup driver-owned per-iteration scratch (fit-owned scratch was freed by fit_iter_finish)
            do i = 1, size(orientations); call orientations(i)%kill; end do
            if( l_bank_fwd )then
                do ithr = 1, nthr
                    call cleanup_plane(mean_fplC(ithr))
                    do q = 1, size(basis_fplsC,1); call cleanup_plane(basis_fplsC(q,ithr)); end do
                end do
                deallocate(mean_fplC, basis_fplsC)
            endif
            if( allocated(stats_es) ) deallocate(stats_es)
            call cleanup_rec_buffers(build, fpls)
            deallocate(orientations, eo)
            if( fit%l_converged )then
                write(logfhandle,'(A,I0,A,F8.3,A,I0,A)') '>>> FLEX_PCA PROBE converged after ',it_eff, &
                    &' iterations: reproducible dim (even|odd) peaked at ',fit%eo_best, &
                    &' and did not improve for ',fit%eo_patience,' iterations'
                call flush(logfhandle)
                exit
            endif
        end do
        call xfsc_teardown(xfctx)
        if( allocated(fit%prev_real) )then
            do q = 1, size(fit%prev_real)
                call fit%prev_real(q)%kill
            end do
            deallocate(fit%prev_real)
        endif
        if( fit%l_pol_es )then
            call polar_grid_kill(fit%pg_es)
            if( allocated(fit%UsallE) ) deallocate(fit%UsallE, fit%CfE, fit%Cm0E, fit%c00E, fit%UbankE, fit%CspE, &
                &fit%xws_es, fit%wr_es, fit%wrd_es, fit%Reb_es)
            if( allocated(fit%dir_es) ) deallocate(fit%dir_es, fit%cae, fit%sae, fit%dused_es, fit%rmatb_es, fit%nrmb_es)
            if( allocated(fit%hex_es) ) deallocate(fit%hex_es, fit%kex_es)
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
        deallocate(fit%z, fit%ppinds)
        ! hoisted model handles back to the caller's arguments
        call move_alloc(fit%basis_recs, basis_recs)
        call move_alloc(fit%eigvals,    eigvals)
        ncomp = fit%ncomp

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

    !> Per-fit stage configuration: every SIMPLE_COV_* selector that becomes fit state, read in
    !! the single-fit entry order. Driver-only selectors (GPU paths, parity check, particle
    !! cache) stay with the callers. Under the paired driver this runs once per fit; the env
    !! values are process-wide constants, so the two reads agree by construction (hazard note:
    !! the moment any of these becomes per-fit, hoist the read to the driver and stamp copies).
    subroutine fit_stage_config( params, fit, nthr )
        class(parameters),  intent(inout) :: params
        type(probe_fit_t),  intent(inout) :: fit
        integer,            intent(in)    :: nthr
        ! ---- POLAR E-STEP: in-batch-loop shared-direction bank (the E-step) ----
        fit%l_pol_es = .true.
        ! hybrid/oversampling accuracy: the derived hybrid radius, no angular oversampling
        fit%rhyb_req   = 0
        fit%l_rhyb_off = .false.
        fit%osamp_pol  = 1
        fit%l_pol_hyb  = .false.
        fit%rhyb_es    = 0
        fit%npos_es    = 0
        fit%l_pol_grid    = .false.
        fit%l_pol_bank_it = .false.
        fit%sec_bank      = 0.
        if( fit%l_pol_es )then
            write(logfhandle,'(A)') '>>> FLEX_PCA POLAR E-STEP ON (stage 1): shared-direction bank &
                &supplies G/b/c; data-plane prep and M-step insertion stay Cartesian'
            call flush(logfhandle)
        endif
        ! ---- MCFA (SIMPLE_COV_EM_MIX=K): tied-covariance K-component mixture latent prior ----
        ! Mixtures of common factor analysers (Baek, McLachlan & Flack, IEEE TPAMI 32:1298,
        ! 2010): ONE shared basis, K Gaussian components in the latent space, fitted in the same
        ! EM as the basis itself. The structural point: a direction earns rank in U only if it
        ! improves the fit under a MULTI-MODAL prior, so purely continuous nuisance variation --
        ! fit equally well by any single component -- gains nothing. This is the one thing a
        ! moment estimator cannot reproduce, because responsibilities are posterior objects.
        ! K=1 pins the component at the origin with a diagonal covariance, which reduces EXACTLY
        ! to the plain PPCA EM -- that is the mandated regression test, not a feature.
        fit%kmix   = COV_EM_MIX
        fit%l_mix_req = fit%kmix >= 1
        ! v2 (2026-08-21): the four MCFA accumulators are additive sufficient statistics and now
        ! ride in the probe part files (PROBE_PART_VERSION 3), together with a bounded latent
        ! subsample so the master can seed mcfa_init. nparts>1 is therefore supported.
        fit%n_mix_warm = 1
        fit%l_mix_active = .false.
        fit%l_mix_used   = .false.
        fit%ldOm_mix     = 0.d0
        fit%ldOm_used    = 0.d0
        if( fit%l_mix_req ) write(logfhandle,'(A,I0,A,I0,A)') '>>> FLEX_PCA MIX requested: K=', &
            &fit%kmix,'  (plain warm-up: ',fit%n_mix_warm,' iterations)'
        if( .not. allocated(fit%sec_proj_thr) )then
            allocate(fit%sec_proj_thr(nthr), fit%sec_gram_thr(nthr), source=0.d0)
        endif
        fit%nll_prev    = 0.d0
        fit%eo_best     = -1.d0
        fit%eo_stall    = 0
        fit%eo_patience = COV_EO_PATIENCE
        ! ON by default. The consensus-shaped term is one the model provably cannot represent:
        ! the contrast coefficient is pinned at 1, so a particle of true amplitude a_i leaves
        ! (a_i - 1) * T_i P(R_i) mu behind, rank one along the consensus and identical in every
        ! particle. Deflating it costs a projection per basis volume per iteration and was
        ! positive in all four cells measured at K=20 on EMPIAR-10076 -- moment 0.1756 -> 0.1906,
        ! EM 0.1193 -> 0.1764, and positive again with latent whitening layered on both arms.
        ! It removes 34.6% of the basis energy on the first EM iteration and 16.2% at convergence.
        ! Validated on 10076 only so far; SIMPLE_COV_EM_DEFLATE=0 restores the old behaviour.
        fit%vdfl           = COV_EM_DEFLATE
        fit%l_deflate_mean = .true.
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
        fit%n_probe_cm  = 0
        fit%nml_plain   = 0
        fit%l_probe_mls = .false.
        if( fit%n_probe_cm > 0 ) write(logfhandle,'(A,I0,A)') &
            &'>>> FLEX_PCA PROBE per-particle contrast ECM ON (', fit%n_probe_cm, ' alternations)'
        if( fit%l_probe_mls ) write(logfhandle,'(A)') &
            &'>>> FLEX_PCA PROBE ML-scaled M-step ON (a z insertion, a^2 density)'
        fit%kfr_ann   = covariance_kfromto(params)
        fit%khi_full  = max(1, fit%kfr_ann(2))
        fit%dstep_ann = real(max(1, params%box_crop - 1)) * params%smpd_crop
        fit%conv_thresh = COV_PROBE_CONV
    end subroutine fit_stage_config

    !> Per-iteration begin: band/lp schedule, prior from the current Gamma, MCFA freeze +
    !! rank-change resize, the per-iteration accumulators and thread scratch, and the
    !! projection-ready expansion of the fit's mean + basis.
    module subroutine fit_iter_begin( params, build, fit, mean_rec, it_eff, niters_eff, nthr )
        class(parameters),   intent(inout) :: params
        type(builder),       intent(inout) :: build
        type(probe_fit_t),   intent(inout) :: fit
        type(reconstructor), intent(inout) :: mean_rec
        integer,             intent(in)    :: it_eff, niters_eff, nthr
        integer :: q
            fit%lp_it  = params%lp
            fit%sec_proj_thr = 0.d0; fit%sec_gram_thr = 0.d0
            fit%sec_bank = 0.; fit%l_pol_bank_it = .false.   ! the polar E-step bank is rebuilt every iteration
            write(logfhandle,'(A,I0,A,I0,A,I0)') '>>> FLEX_PCA PROBE SUBSPACE ITERATION ',it_eff,' / ',niters_eff, &
                &'  basis dim=',fit%ncomp
            call flush(logfhandle)
            allocate(fit%prior(fit%ncomp))
            do q = 1, fit%ncomp
                fit%prior(q) = 1.d0 / max(fit%eigvals(q), DTINY)
            end do
            ! MCFA bookkeeping. l_mix_used freezes which prior THIS iteration's E-step runs
            ! under (the M-step hook below flips l_mix_active mid-iteration); a basis dimension
            ! change invalidates xi/Omega, so the mixture re-warms and re-initialises.
            fit%l_mix_used = fit%l_mix_active
            fit%ldOm_used  = fit%ldOm_mix
            if( fit%l_mix_req )then
                if( allocated(fit%rhs0th) )then
                    if( size(fit%rhs0th,1) /= fit%ncomp )then
                        deallocate(fit%rhs0th, fit%mkth, fit%lwth, fit%rkth, fit%mxa_sr, fit%mxa_sm, fit%mxa_smm, fit%mxa_sainv)
                        if( allocated(fit%mix_xi) ) deallocate(fit%mix_xi, fit%mix_Om, fit%mix_Ominv, fit%mix_pi, &
                            &fit%mix_Omxi, fit%mix_xiOx, fit%mix_lpi)
                        if( fit%l_mix_active ) write(logfhandle,'(A)') &
                            &'>>> FLEX_PCA MIX re-initialising: basis dimension changed'
                        fit%l_mix_active = .false.
                        fit%l_mix_used   = .false.
                    endif
                endif
                if( .not. allocated(fit%rhs0th) )then
                    allocate(fit%rhs0th(fit%ncomp,nthr), fit%mkth(fit%ncomp,fit%kmix,nthr), fit%lwth(fit%kmix,nthr), &
                        &fit%rkth(fit%kmix,nthr), fit%mxa_sr(fit%kmix,nthr), fit%mxa_sm(fit%ncomp,fit%kmix,nthr), &
                        &fit%mxa_smm(fit%ncomp,fit%ncomp,fit%kmix,nthr), fit%mxa_sainv(fit%ncomp,fit%ncomp,nthr))
                endif
                fit%mxa_sr = 0.d0; fit%mxa_sm = 0.d0; fit%mxa_smm = 0.d0; fit%mxa_sainv = 0.d0
            endif
            ! even/odd Y_q accumulators (half-set FSC regularization) + the COUPLED latent normal matrix.
            ! rho carries one entry per (q,r) pair, not one shared density: the M-step below solves the
            ! components together at every grid point.
            allocate(fit%Yeven(fit%ncomp), fit%Yodd(fit%ncomp))
            do q = 1, fit%ncomp
                call init_basis_reconstructor(params, build, fit%Yeven(q)); call fit%Yeven(q)%reset; call fit%Yeven(q)%reset_exp
                call init_basis_reconstructor(params, build, fit%Yodd(q));  call fit%Yodd(q)%reset;  call fit%Yodd(q)%reset_exp
            end do
            fit%es     = shape(fit%Yeven(1)%cmat_exp)
            ! The per-voxel coupled normal matrix is always the FULL ncomp x ncomp SPD system.
            ! A diagonal approximation used to be selectable here; it dropped the cross-component
            ! terms, and on EMPIAR-10028 that cost the 40S rotation outright -- the mode only
            ! appeared once the full solve was restored. It is not a speed/accuracy trade worth
            ! offering, so the option is gone rather than merely defaulted off.
            fit%npairs = (fit%ncomp*(fit%ncomp+1))/2
            allocate(fit%rho_e(fit%npairs,fit%es(1),fit%es(2),fit%es(3)), fit%rho_o(fit%npairs,fit%es(1),fit%es(2),fit%es(3)), source=0.)
            write(logfhandle,'(A,I0,A,F8.2,A)') '>>> FLEX_PCA PROBE coupled normal matrix rows=',fit%npairs, &
                &' (full)  rho even+odd ', 8.d0*real(fit%npairs,dp)*real(fit%es(1),dp)*real(fit%es(2),dp)*real(fit%es(3),dp)/1.d9,' GB'
            call flush(logfhandle)
            ! rec_backend=pcg: the solver's lattice and support at this rank, the packed kernel sums zeroed
            ! (the full-range accumulators are allocated by the first batch insert of an inserting process)
            fit%l_pcg = trim(params%rec_backend) == 'pcg'
            if( allocated(fit%kacc_e) ) deallocate(fit%kacc_e)
            if( allocated(fit%kacc_o) ) deallocate(fit%kacc_o)
            if( allocated(fit%kpk_e) ) deallocate(fit%kpk_e)
            if( allocated(fit%kpk_o) ) deallocate(fit%kpk_o)
            if( allocated(fit%racc_e) ) deallocate(fit%racc_e)
            if( allocated(fit%racc_o) ) deallocate(fit%racc_o)
            if( allocated(fit%rpk_e) ) deallocate(fit%rpk_e)
            if( allocated(fit%rpk_o) ) deallocate(fit%rpk_o)
            if( fit%l_pcg )then
                call fit%pcg%new(params%box_crop, params%smpd_crop, fit%ncomp)
                call flex_pcg_install_window(fit%pcg, params)
                call fit%pcg%alloc_packed(fit%kpk_e)
                call fit%pcg%alloc_packed(fit%kpk_o)
                call fit%pcg%alloc_rhs_packed(fit%rpk_e)
                call fit%pcg%alloc_rhs_packed(fit%rpk_o)
                write(logfhandle,'(A,F8.2,A,F8.2,A,F8.2,A,F8.2,A)') '>>> FLEX_PCA PROBE PCG pair kernels: packed even+odd ', &
                    &2.d0*fit%pcg%bytes_packed()/1.d9, ' GB, accumulators even+odd ', 2.d0*fit%pcg%bytes_accum()/1.d9, &
                    &' GB; right-hand sides: packed ', 2.d0*fit%pcg%bytes_rhs_packed()/1.d9, ' GB, accumulators ', &
                    &2.d0*fit%pcg%bytes_rhs_accum()/1.d9, ' GB'
                call flush(logfhandle)
            endif
            allocate(fit%Gth(fit%ncomp,fit%ncomp,nthr), fit%Ath(fit%ncomp,fit%ncomp,nthr), fit%bth(fit%ncomp,nthr), fit%cth(fit%ncomp,nthr), fit%zth(fit%ncomp,nthr))
            allocate(fit%Ainvth(fit%ncomp,fit%ncomp,nthr), fit%Acpth(fit%ncomp,fit%ncomp,nthr))
            allocate(fit%basis_fpls(fit%ncomp,nthr), fit%mean_fpl(nthr))
            allocate(fit%zbatch(fit%ncomp,MAXIMGBATCHSZ), fit%dens(fit%ncomp,fit%ncomp,MAXIMGBATCHSZ))
            allocate(fit%valid(MAXIMGBATCHSZ), fit%valid_e(MAXIMGBATCHSZ), fit%valid_o(MAXIMGBATCHSZ))
            allocate(fit%gam_thr(fit%ncomp,nthr), source=0.d0)
            allocate(fit%gam_acc(fit%ncomp), source=0.d0)
            allocate(fit%nval_thr(nthr), source=0)
            allocate(fit%hth(fit%ncomp,nthr), source=0.d0)
            allocate(fit%nll_thr(nthr), source=0.d0)
            allocate(fit%gam_dbg(4,nthr), source=0.d0)
            fit%nll_tot = 0.d0
            fit%dens = 0.d0
            ! the Gamma reduction target (filled by fit_iter_reduce or the distributed
            ! master's part reduce; freed by fit_iter_finish at the gam_acc step)
            allocate(fit%gam_sum(fit%ncomp), source=0.d0)
            fit%nval = 0
            call mean_rec%expand_exp
            do q = 1, fit%ncomp
                call fit%basis_recs(q)%expand_exp
            end do
    end subroutine fit_iter_begin

    !> PAIRED-ENGINE ENTRY (SIMPLE_COV_PAIRED=1; plan step 1-2, shared-memory probe-only v1):
    !! two resident probe fits over the mod-4 disjoint halves of the master's selection,
    !! advanced by ONE shared master loop. Per-fit initialisation mirrors the single-fit order:
    !! mean copy + mean scale on the fit's own half, deterministic data-free basis (identical
    !! geometry both fits) with per-fit noise/prior calibration, probe-stage subsample.
    !! Delivery is probe-only: fit A keeps the legacy namespaces (flex_pca_pc*/flex_pca_probe.txt)
    !! so every downstream consumer keeps working; fit B writes flex_pca_fitB_pc* +
    !! flex_pca_probe_fitB.txt (naming precedent: the bagA/bagB pools).
    module subroutine run_flex_pca_paired( params, build, pinds, nptcls, col_sep, neigs_req, &
        &m_basis_recs, m_eigvals, m_ncomp, m_sig2, m_matchcos , rounds)
        class(flex_pca_rounds), intent(inout) :: rounds
        class(parameters),   intent(inout) :: params
        type(builder),       intent(inout) :: build
        integer,             intent(in)    :: pinds(:), nptcls, col_sep, neigs_req
        type(reconstructor), allocatable, optional, intent(out) :: m_basis_recs(:)
        real(dp),            allocatable, optional, intent(out) :: m_eigvals(:)
        integer,             optional,    intent(out)           :: m_ncomp
        real(dp),            optional,    intent(out)           :: m_sig2
        real(dp), allocatable, optional,  intent(out)           :: m_matchcos(:)
        type(probe_fit_t) :: fits(2)
        integer, allocatable :: half_pinds(:)
        integer :: vpair, f, i, cnt, nhalf(2), u, vmerge
        ! ---- FINAL STAGE (par.7 "merge, don't refit"): the accumulator-level merge is the delivery.
        ! vmerge is kept as the mode selector consumed by the stash hook and the merge itself
        ! (1 = accumulators, 2 = the pooled-Gram fallback taken when the frame map degenerates).
        vmerge = 1
        if( .not. (present(m_basis_recs) .and. present(m_eigvals) .and. present(m_ncomp) &
            &.and. present(m_sig2)) ) THROW_HARD('the paired merge requires the merged-product arguments')
        if( params%n_probe_iters < 2 ) THROW_HARD('the paired merge needs n_probe_iters >= 2: the final iteration''s statistics are expressed in the previous iteration''s delivered frame')
        if( present(m_ncomp) ) m_ncomp = 0
        if( present(m_sig2)  ) m_sig2  = 0.d0
        ! ---- the split: the ONE rule, shared with the two-job pcafit harness ----
        vpair = 1
        call cov_env_int_pub('SIMPLE_COV_MOD4_PAIRING', vpair)
        if( vpair == 2 ) THROW_HARD('SIMPLE_COV_MOD4_PAIRING=2 groups same-parity rows: both &
            &halves lose one eo class under the row-alternating project eo split. Use pairing 1 or 3.')
        if( vpair /= 1 .and. vpair /= 3 ) THROW_HARD('SIMPLE_COV_MOD4_PAIRING must be 1, 2 or 3')
        do f = 1, 2
            cnt = 0
            do i = 1, nptcls
                if( flex_pca_half_of(pinds(i), vpair) == f ) cnt = cnt + 1
            end do
            nhalf(f) = cnt
            if( cnt < 100 ) THROW_HARD('paired engine: a mod-4 half has fewer than 100 particles')
        end do
        write(logfhandle,'(A,I0,A,I0,A,I0,A)') '>>> FLEX_PCA PAIRED ENGINE: two resident fits, &
            &one shared loop; half A=',nhalf(1),'  half B=',nhalf(2),'  (mod-4 pairing ',vpair,')'
        call flush(logfhandle)
        do f = 1, 2
            allocate(half_pinds(nhalf(f)))
            cnt = 0
            do i = 1, nptcls
                if( flex_pca_half_of(pinds(i), vpair) == f )then
                    cnt = cnt + 1; half_pinds(cnt) = pinds(i)
                endif
            end do
            if( f == FLEX_FIT_A )then
                call new_probe_fit(fits(f), f, 'flex_pca_pc', COV_PROBE_META, half_pinds, nhalf(f))
            else
                call new_probe_fit(fits(f), f, 'flex_pca_fitB_pc', 'flex_pca_probe_fitB.txt', &
                    &half_pinds, nhalf(f))
            endif
            deallocate(half_pinds)
        end do
        ! ---- per-fit initialisation, mirroring the single-fit order ----
        do f = 1, 2
            ! per-fit mean copy carrying the per-fit mean scale, fitted on the fit's own half
            ! (the two-job harness scales the mean on its own half, so the paired engine must
            ! too or T0 compares different means). Shared memory: estimate_mean_scale writes
            ! no cache file here (rounds%is_master() is false), so no fname clash.
            call init_mean_reconstructor(params, build, fits(f)%mean_rec)
            ! per-fit cache namespace: on the DISTRIBUTED master estimate_mean_scale writes the
            ! fitted radial scale for the workers to apply (per-fit files, or the two fits clash)
            if( f == FLEX_FIT_A )then
                call estimate_mean_scale(params, build, fits(f)%mean_rec, fits(f)%pinds, &
                    &fits(f)%nptcls, cache_fname='flex_pca_mean_scale.bin', rounds=rounds)
            else
                call estimate_mean_scale(params, build, fits(f)%mean_rec, fits(f)%pinds, &
                    &fits(f)%nptcls, cache_fname='flex_pca_mean_scale_fitB.bin', rounds=rounds)
            endif
            ! deterministic data-free basis; its calibration pass runs on the fit's half ->
            ! per-fit sig2/gamma0; the init eigenvolume write uses the fit's namespace
            call init_basis_datafree(params, build, fits(f)%mean_rec, fits(f)%pinds, fits(f)%nptcls, &
                &col_sep, neigs_req, fits(f)%basis_recs, fits(f)%eigvals, fits(f)%ncomp, &
                &fits(f)%sig2_eff, fprefix=fits(f)%fprefix%to_char(), rounds=rounds)
            fits(f)%sig2 = max(fits(f)%sig2_eff, DTINY)
            ! probe-stage subsample of the fit's half (master process: no nparts division)
            call cov_stage_subsample(build, fits(f)%pinds, fits(f)%nptcls, 1, COV_PROBE_MAX_PTCLS, &
                &'SIMPLE_COV_PROBE_MAX', 'PROBE', &
                &fits(f)%ppinds, fits(f)%npp)
            allocate(fits(f)%z(fits(f)%npp, fits(f)%ncomp))
        end do
        call probe_subspace_paired(params, build, fits, params%n_probe_iters, vmerge == 1, vpair, rounds=rounds)
        ! ---- delivery (v1, probe-only): both fits' probe metas + eigenvalue tables; the
        ! eigenvolumes are already on disk from the last iteration under each fit's prefix ----
        do f = 1, 2
            call save_probe_state(fits(f)%ncomp, fits(f)%eigvals, fits(f)%sig2_eff, &
                &fname=fits(f)%meta_fname%to_char())
            call write_paired_eigen_table(fits(f))
        end do
        ! the paired record (probe-only mode never reaches the covariance manifest)
        call del_file('flex_pca_paired.txt')
        open(newunit=u, file='flex_pca_paired.txt', status='replace', action='write')
        write(u,'(A)') '# paired-engine record: paired  pairing  fit  nptcls  npp  ncomp  sig2'
        do f = 1, 2
            write(u,'(I2,1X,I2,1X,A1,1X,I10,1X,I10,1X,I5,1X,ES16.8)') 1, vpair, &
                &merge('A','B',f==FLEX_FIT_A), fits(f)%nptcls, fits(f)%npp, fits(f)%ncomp, &
                &fits(f)%sig2_eff
        end do
        close(u)
        write(logfhandle,'(A,I0,A,I0,A,ES11.4,A,ES11.4)') '>>> FLEX_PCA PAIRED delivered: &
            &ncomp A=',fits(1)%ncomp,' B=',fits(2)%ncomp,'  sig2 A=',fits(1)%sig2_eff, &
            &' B=',fits(2)%sig2_eff
        call flush(logfhandle)
        ! ---- FINAL STAGE (par.7): frame-align + accumulator merge + ONE joint solve (mode 1)
        ! or the pooled-Gram fallback (mode 2); merged eigenvolumes/meta/manifest written inside.
        ! The fits' own delivery above is untouched -- the merge is an ADDITIONAL product.
        if( vmerge > 0 ) call probe_paired_merge(params, build, fits, vmerge, m_basis_recs, &
            &m_eigvals, m_ncomp, m_sig2, m_matchcos=m_matchcos)
        do f = 1, 2
            call kill_probe_fit(fits(f))
        end do
    end subroutine run_flex_pca_paired

    !> Per-fit eigenvalue table (the paired-local equivalent of write_covariance_eigenvolumes'
    !! table): fit A keeps the legacy name, fit B gets the fitB namespace.
    subroutine write_paired_eigen_table( fit )
        type(probe_fit_t), intent(in) :: fit
        character(len=:), allocatable :: fn
        integer :: q, u
        if( fit%id == FLEX_FIT_A )then
            fn = 'flex_pca_eigenvalues.txt'
        else
            fn = 'flex_pca_eigenvalues_fitB.txt'
        endif
        call del_file(fn)
        open(newunit=u, file=fn, status='replace', action='write')
        write(u,'(A)') '# component eigenvalue'
        do q = 1, fit%ncomp
            write(u,'(I6,1X,ES20.10)') q, fit%eigvals(q)
        end do
        close(u)
    end subroutine write_paired_eigen_table

    !> THE SHARED-MEMORY PAIRED MASTER LOOP (plan §3.3): one it_eff advances BOTH fits. Each
    !! iteration merges the two fits' per-fit probe subsets into one sorted-by-project-row read
    !! list, runs ONE read/prep pass, E-steps each particle against exactly the basis of the fit
    !! that owns it (per-particle cost unchanged; what doubles is basis memory, per-fit
    !! accumulators, the master's per-voxel solves and the per-iteration bank build), then runs
    !! the per-fit CPU inserts, reductions and master tails. Supported E-step bodies are the CPU
    !! polar former and the plain Cartesian former (plan §3.5); every GPU path stands down --
    !! the banked-adjoint ppinds permutation therefore never runs, which the merged read list's
    !! ascending-order assumption depends on (hazard 1).
    subroutine probe_subspace_paired( params, build, fits, niters, l_merge_stash, vpair , rounds)
        class(flex_pca_rounds), intent(inout) :: rounds
        class(parameters), intent(inout) :: params
        type(builder),     intent(inout) :: build
        type(probe_fit_t), intent(inout) :: fits(2)
        integer,           intent(in)    :: niters
        !> par.7 final stage (SIMPLE_COV_PAIRED_MERGE=1): snapshot each fit's raw M-step
        !! sufficient statistics + entry frame every iteration, pre-ridge pre-solve
        logical,           intent(in)    :: l_merge_stash
        !> the mod-4 pairing, stamped into both probe-state files for distributed dispatch
        integer,           intent(in)    :: vpair
        integer  :: it_eff, niters_eff, f, i, nthr
        integer(timer_int_kind) :: t_it
        !> cross-fit-FSC driver context: the paired master writes honest paired=1 records every
        !! iteration; the ridge/marching/stopping consumers act per their own gates
        type(xfsc_ctx_t) :: xfctx
        !> phase-2 distributed dispatch: the E-step fans out over parts, one qsys round per
        !! iteration, one v5 part per worker carrying BOTH fits' blocks
        logical :: l_pdistr
        nthr = omp_get_max_threads()
        l_pdistr = rounds%distributed()
        if( l_pdistr )then
            write(logfhandle,'(A,I0,A)') '>>> FLEX_PCA PAIRED DISTRIBUTED: E-step over ', &
                &rounds%nparts(), ' parts, one v5 part per worker, per-fit reduce'
            call flush(logfhandle)
        endif
        ! ---- per-fit stage config: process-wide env values, read once per fit ----
        do f = 1, 2
            call fit_stage_config(params, fits(f), nthr)
        end do
        ! ---- CROSS-FIT-FSC setup (spec par.7): paired master is always a writer ----
        call xfsc_setup(xfctx, params, fits(1)%kfr_ann, .true., .true.)
        ! ---- refusals and stand-downs (plan §3.2/§3.5), printed once at driver start ----
        if( flex_gpu_available() .and. .not. cov_env_int_off('SIMPLE_COV_GPU') ) &
            &write(logfhandle,'(A)') '>>> FLEX_PCA PAIRED: CPU M-step insertion (single-fit &
            &device accumulator)'
        write(logfhandle,'(A,I0,A,I0,A)') '>>> FLEX_PCA PAIRED: basis dim A=',fits(1)%ncomp, &
            &' B=',fits(2)%ncomp,'  CPU+POLAR E-step, device E-step/bank/banked-adjoint down'
        call flush(logfhandle)
        niters_eff = niters
        do it_eff = 1, niters
            t_it = tic()
            do f = 1, 2
                call fit_iter_begin(params, build, fits(f), fits(f)%mean_rec, it_eff, niters_eff, nthr)
                fits(f)%z = 0.d0
            end do
            if( l_pdistr )then
                ! ---- phase-2 dispatch: stamp BOTH fits' probe-state files with the global
                ! iteration AND the pairing (the worker's dispatch flag -- plan par.7.1: not via
                ! env), refresh nothing else (each fit's eigenvolumes are already on disk under
                ! its own prefix from the previous fit_iter_finish / the data-free init), run ONE
                ! qsys round, then fold every worker's v5 part into the per-fit accumulators ----
                do f = 1, 2
                    call save_probe_state(fits(f)%ncomp, fits(f)%eigvals, fits(f)%sig2_eff, &
                        &fname=fits(f)%meta_fname%to_char())
                end do
                call rounds%run_stage(params, PCA_STAGE_PROBE, 'paired probe iteration', &
                    &which_iter=it_eff, maxits=niters_eff, nfits=2)
                call paired_reduce_parts_v5(params, fits, rounds=rounds)
            else
                call paired_estep_pass(params, build, fits, it_eff, nthr)
            endif
            do f = 1, 2
                write(logfhandle,'(A,A,A,I0)') '>>> FLEX_PCA PAIRED FIT ', merge('A','B',f==1), &
                    &' master tail, it=', it_eff
                call flush(logfhandle)
                ! in-process: fold the thread accumulators; distributed: the v5 part reduce
                ! above already summed every worker's contribution into the fit accumulators
                if( .not. l_pdistr ) call fit_iter_reduce(fits(f), it_eff, nthr, rounds=rounds)
                ! crossfsc per-iteration prep: harvest flag for the writer payloads (H from the
                ! rho pair diagonals pre-ridge, internal FSC + Gamma) and the per-fit SSNR ridge
                ! from record t-1 (applied by fit_iter_finish just before the coupled solves;
                ! fit X's ridge maps the record through fit X's side of the signed permutation
                ! with fit X's own H -- the spec par.2.3 invariant, keyed off fit%id here)
                call xfsc_prep_iter(xfctx, params, fits(f), it_eff, &
                    &merge('  fit=A','  fit=B', f == 1))
                ! par.7 merge stash: raw last-iteration statistics, BEFORE fit_iter_finish
                ! ridges rho / solves the numerators in place / frees the accumulators. Every
                ! iteration overwrites -- any iteration can turn out to be the last.
                if( l_merge_stash ) call probe_fit_merge_stash(fits(f))
                call fit_iter_finish(params, build, fits(f), it_eff, nthr)
                write(logfhandle,'(A,A,A,I0,A,I0,A,ES12.4)') '>>> FLEX_PCA PAIRED FIT ', &
                    &merge('A','B',f==1), ' it=', it_eff, '  refined dim=', fits(f)%ncomp, &
                    &'  max var=', maxval(fits(f)%eigvals)
                call flush(logfhandle)
            end do
            ! impl-map step 3 hook: the cross-fit comparison + artifact record + stopping series.
            ! One honest paired=1 record per iteration from the two fits' realized bases and
            ! stashed payloads (v1-minimal matched blocks; see xfsc_paired_record).
            if( xfctx%l_writer ) call xfsc_paired_record(xfctx, params, fits, it_eff)
            write(logfhandle,'(A,I0,A,I0,A,F8.1)') '>>> FLEX_PCA PAIRED ITER ', it_eff, ' / ', &
                &niters_eff, '  seconds=', toc(t_it)
            call flush(logfhandle)
            ! v1: exit when BOTH fits have converged (a converged fit keeps updating until then;
            ! freezing it while resident is a step-3 concern together with the cross-fit series)
            if( fits(1)%l_converged .and. fits(2)%l_converged )then
                write(logfhandle,'(A,I0,A)') '>>> FLEX_PCA PAIRED converged after ', it_eff, &
                    &' iterations (both fits)'
                call flush(logfhandle)
                exit
            endif
        end do
        call xfsc_teardown(xfctx)
    end subroutine probe_subspace_paired

    !> DISTRIBUTED PAIRED WORKER (plan par.7.2): one relaunch per master iteration. Loads BOTH
    !! fits' metas and bases (per-fit namespaces), splits this worker's fromp/top shard by the
    !! ONE mod-4 rule, runs the same shared E-step pass the shared-memory master runs, and
    !! writes ONE v5 part carrying both fits' accumulator blocks. No master tail runs here.
    module subroutine run_flex_pca_paired_worker( params, build, pinds, nptcls, it_stamp, &
        &niters_stamp, vpair , rounds)
        class(flex_pca_rounds), intent(inout) :: rounds
        class(parameters), intent(inout) :: params
        type(builder),     intent(inout) :: build
        integer,           intent(in)    :: pinds(:), nptcls, it_stamp, niters_stamp, vpair
        integer, parameter :: MIX_ZSUB_MAX5 = 2000  ! keep equal to probe_subspace_iteration's
        type(probe_fit_t) :: fits(2)
        integer,  allocatable :: half_pinds(:)
        real(dp), allocatable :: ev_load(:)
        complex,  allocatable :: cme(:,:,:,:), cmo(:,:,:,:)
        real,     allocatable :: rhe(:,:,:,:), rhoo(:,:,:,:)
        type(string) :: pfname, tmp_fname
        integer  :: nhalf(2), f, i, cnt, nthr, q, funit, nc_load, it_eff, niters_eff
        real(dp) :: s2_load
        nthr = omp_get_max_threads()
        it_eff     = max(1, it_stamp)
        niters_eff = max(1, niters_stamp)
        if( vpair /= 1 .and. vpair /= 3 ) THROW_HARD('paired worker: invalid mod-4 pairing stamp')
        do f = 1, 2
            cnt = 0
            do i = 1, nptcls
                if( flex_pca_half_of(pinds(i), vpair) == f ) cnt = cnt + 1
            end do
            nhalf(f) = cnt
            if( cnt < 2 ) THROW_HARD('paired worker: a mod-4 half of this shard is (near) empty')
        end do
        write(logfhandle,'(A,I0,A,I0,A,I0,A,I0)') '>>> FLEX_PCA PAIRED WORKER part=', params%part, &
            &'  it_eff=', it_eff, '  shard half A=', nhalf(1), '  B=', nhalf(2)
        call flush(logfhandle)
        do f = 1, 2
            allocate(half_pinds(nhalf(f)))
            cnt = 0
            do i = 1, nptcls
                if( flex_pca_half_of(pinds(i), vpair) == f )then
                    cnt = cnt + 1; half_pinds(cnt) = pinds(i)
                endif
            end do
            if( f == FLEX_FIT_A )then
                call new_probe_fit(fits(f), f, 'flex_pca_pc', COV_PROBE_META, half_pinds, nhalf(f))
            else
                call new_probe_fit(fits(f), f, 'flex_pca_fitB_pc', 'flex_pca_probe_fitB.txt', &
                    &half_pinds, nhalf(f))
            endif
            deallocate(half_pinds)
            ! per-fit mean: deterministic rebuild from vol1 + the MASTER-fitted per-fit radial
            ! scale (a worker must never re-fit the scale on its own shard)
            call init_mean_reconstructor(params, build, fits(f)%mean_rec)
            if( f == FLEX_FIT_A )then
                call apply_cached_mean_scale(params, fits(f)%mean_rec, &
                    &cache_fname='flex_pca_mean_scale.bin')
            else
                call apply_cached_mean_scale(params, fits(f)%mean_rec, &
                    &cache_fname='flex_pca_mean_scale_fitB.bin')
            endif
            ! per-fit basis + meta, refreshed by the master before this round was scheduled
            call load_probe_state(nc_load, ev_load, s2_load, fname=fits(f)%meta_fname%to_char())
            fits(f)%ncomp    = nc_load
            call move_alloc(ev_load, fits(f)%eigvals)
            fits(f)%sig2_eff = s2_load
            fits(f)%sig2     = max(s2_load, DTINY)
            call load_probe_basis(params, build, fits(f)%ncomp, fits(f)%basis_recs, &
                &fprefix=fits(f)%fprefix%to_char())
            ! probe-stage subsample of this shard's half; the budget is a TOTAL across processes,
            ! so a worker divides by nparts (cov_stage_subsample contract)
            call cov_stage_subsample(build, fits(f)%pinds, fits(f)%nptcls, params%nparts, &
                &COV_PROBE_MAX_PTCLS, 'SIMPLE_COV_PROBE_MAX', 'PROBE', &
                &fits(f)%ppinds, fits(f)%npp)
            allocate(fits(f)%z(fits(f)%npp, fits(f)%ncomp))
            call fit_stage_config(params, fits(f), nthr)
        end do
        ! ---- one accumulate iteration keyed to the master's global stamp ----
        do f = 1, 2
            call fit_iter_begin(params, build, fits(f), fits(f)%mean_rec, it_eff, niters_eff, nthr)
            fits(f)%z = 0.d0
        end do
        call paired_estep_pass(params, build, fits, it_eff, nthr)
        do f = 1, 2
            call fit_iter_reduce(fits(f), it_eff, nthr, rounds=rounds)
            fits(f)%nll_tot = sum(fits(f)%nll_thr)
        end do
        ! ---- ONE v5 part: both fits' blocks, in fit order ----
        pfname = flex_pca_part_fname('probe', params%part, params%numlen)
        call open_probe_part_v5_write(pfname, 2, funit, tmp_fname)
        do f = 1, 2
            allocate(cme(fits(f)%es(1),fits(f)%es(2),fits(f)%es(3),fits(f)%ncomp), &
                &cmo(fits(f)%es(1),fits(f)%es(2),fits(f)%es(3),fits(f)%ncomp))
            allocate(rhe(fits(f)%es(1),fits(f)%es(2),fits(f)%es(3),fits(f)%ncomp), &
                &rhoo(fits(f)%es(1),fits(f)%es(2),fits(f)%es(3),fits(f)%ncomp))
            do q = 1, fits(f)%ncomp
                cme(:,:,:,q) = fits(f)%Yeven(q)%cmat_exp; rhe(:,:,:,q)  = fits(f)%Yeven(q)%rho_exp
                cmo(:,:,:,q) = fits(f)%Yodd(q)%cmat_exp;  rhoo(:,:,:,q) = fits(f)%Yodd(q)%rho_exp
            end do
            if( fits(f)%l_mix_req )then
                block
                    real(dp), allocatable :: w_sr(:), w_sm(:,:), w_smm(:,:,:), w_sai(:,:), w_z(:,:)
                    integer :: tt2, kk4, nzs, izs, istep
                    allocate(w_sr(fits(f)%kmix), w_sm(fits(f)%ncomp,fits(f)%kmix), &
                        &w_smm(fits(f)%ncomp,fits(f)%ncomp,fits(f)%kmix), &
                        &w_sai(fits(f)%ncomp,fits(f)%ncomp))
                    w_sr = 0.d0; w_sm = 0.d0; w_smm = 0.d0; w_sai = 0.d0
                    do tt2 = 1, nthr
                        w_sr  = w_sr  + fits(f)%mxa_sr(:,tt2)
                        w_sm  = w_sm  + fits(f)%mxa_sm(:,:,tt2)
                        w_sai = w_sai + fits(f)%mxa_sainv(:,:,tt2)
                        do kk4 = 1, fits(f)%kmix
                            w_smm(:,:,kk4) = w_smm(:,:,kk4) + fits(f)%mxa_smm(:,:,kk4,tt2)
                        end do
                    end do
                    ! deterministic stride subsample of this part's latents (bounded)
                    nzs   = min(MIX_ZSUB_MAX5, fits(f)%npp)
                    istep = max(1, fits(f)%npp / max(1,nzs))
                    nzs   = min(nzs, (fits(f)%npp + istep - 1)/istep)
                    allocate(w_z(nzs,fits(f)%ncomp))
                    do izs = 1, nzs
                        w_z(izs,:) = fits(f)%z(min(fits(f)%npp, 1 + (izs-1)*istep), :)
                    end do
                    call write_probe_part_v5_fit(funit, cme, rhe, cmo, rhoo, fits(f)%rho_e, &
                        &fits(f)%rho_o, fits(f)%gam_sum, fits(f)%nll_tot, fits(f)%nval, &
                        &fits(f)%ncomp, fits(f)%kpk_e, fits(f)%kpk_o, fits(f)%rpk_e, fits(f)%rpk_o, w_sr, w_sm, w_smm, w_sai, w_z)
                    deallocate(w_sr, w_sm, w_smm, w_sai, w_z)
                end block
            else
                call write_probe_part_v5_fit(funit, cme, rhe, cmo, rhoo, fits(f)%rho_e, &
                    &fits(f)%rho_o, fits(f)%gam_sum, fits(f)%nll_tot, fits(f)%nval, fits(f)%ncomp, &
                    &fits(f)%kpk_e, fits(f)%kpk_o, fits(f)%rpk_e, fits(f)%rpk_o)
            endif
            deallocate(cme, cmo, rhe, rhoo)
        end do
        call close_probe_part_v5_write(funit, tmp_fname, pfname)
        call pfname%kill
        write(logfhandle,'(A,I0,A,I0,A,I0)') '>>> FLEX_PCA PAIRED WORKER part=', params%part, &
            &'  wrote v5 part: valid A=', fits(1)%nval, '  B=', fits(2)%nval
        call flush(logfhandle)
        do f = 1, 2
            call kill_probe_fit(fits(f))
        end do
    end subroutine run_flex_pca_paired_worker

end submodule simple_flex_pca_em_iter
