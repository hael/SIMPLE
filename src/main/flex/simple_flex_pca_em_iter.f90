!@descr: flex_pca EM: the subspace EM iteration (E-step accumulation, M-step update)
submodule (simple_flex_pca_em) simple_flex_pca_em_iter
use simple_matcher_3Drec,   only: init_rec, cleanup_rec_buffers
use simple_matcher_ptcl_io, only: prepimgbatch
use simple_flex_reconstructor_latent_ops, only: project_fplane_mean,&
    &insert_planes_oversamp_coupled_batch_scaled
use simple_flex_reconstructor_latent_ops, only: solve_coupled_basis_exp, add_invtausq2rho_coupled
use simple_flex_pca_crossfsc, only: crossfsc_file, crossfsc_record, crossfsc_load, crossfsc_write,&
    &crossfsc_append, crossfsc_latest_upto, crossfsc_kill, crossfsc_kill_record, crossfsc_to_invtau2,&
    &crossfsc_harvest_h, crossfsc_stop_stat, crossfsc_inband_mean, crossfsc_khi_deepest,&
    &COV_XFSC_FNAME
use simple_flex_pca_util,   only: cov_env_flag_on, cov_env_flag_off
use simple_flex_pca_polar,  only: polar_grid_build, polar_grid_kill, polar_project_recs,&
    &polar_sample_particle_fused
implicit none
#include "simple_local_flags.inc"

contains

    !> Probe-based subspace iteration: alternate a Wiener E-step (per-particle latents in the current
    !! basis) with a weighted-backprojection M-step (Y_q += sum_i z_iq * backproject(r_i)), then
    !! orthonormalize the refined probe volumes into the next basis.
    module subroutine probe_subspace_iteration( params, build, mean_rec, basis_recs, eigvals, sig2_eff, &
        &pinds, nptcls, ncomp, niters, it_glob, niters_glob, fprefix, meta_fname, rounds)
        use simple_flex_pca_plane_cache, only: plane_cache_in_use
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
        integer, parameter :: MIX_ZSUB_MAX = 2000   ! per-part latent subsample shipped for mcfa_init
        integer(timer_int_kind) :: t_bank
        logical  :: l_probe_distr_pre
        integer,             allocatable :: eo(:)
        real(dp) :: a, aa, e_mm, myv
        integer  :: it, q, i, ithr, nthr, batchlims(2), batchsz, ibatch, row
        !> effective (global) iteration numbering -- the ONLY counters iteration-keyed schedules
        !! and iteration logs may use; equal to it/niters except on a distributed probe worker
        integer  :: it_eff, niters_eff
        integer  :: nparts_sub
        logical  :: l_probe_distr
        complex,             allocatable :: cme(:,:,:,:), cmo(:,:,:,:)
        real,                allocatable :: rhe(:,:,:,:), rhoo(:,:,:,:)
        integer(timer_int_kind) :: t_it, t_sec
        real(timer_int_kind) :: sec_read, sec_prep, sec_estep, sec_ins
        real(dp) :: twp0
        logical  :: l_pcache
        !> the hoisted per-fit state (M2): every cross-iteration and per-iteration-per-fit
        !! variable of this routine now lives in one probe_fit_t owned by the driver
        type(probe_fit_t) :: fit
        !> cross-fit-FSC driver context (ridge only: the single-fit engine writes no records)
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
        ! Optional probe-stage cap (SIMPLE_COV_PROBE_MAX, default off), applied per halfset; embedding
        ! still uses every particle. The cap is a cross-process total: only a worker divides it by nparts.
        nparts_sub = 1
        if( rounds%is_worker() ) nparts_sub = params%nparts
        call cov_stage_subsample(build, fit%pinds, fit%nptcls, nparts_sub, COV_PROBE_MAX_PTCLS, &
            &'SIMPLE_COV_PROBE_MAX', 'PROBE', fit%ppinds, fit%npp)
        allocate(fit%z(fit%npp,fit%ncomp))
        fit%z = 0.d0
        ! per-fit stage configuration (hard-wired defaults)
        call fit_stage_config(params, fit, nthr)
        ! ---- CROSS-FIT-FSC setup: the single-fit engine writes no records; its ridge reads any paired artifact ----
        call xfsc_setup(xfctx, params, fit%kfr_ann, .false., .not. rounds%is_worker())
        l_probe_distr_pre = rounds%distributed()
        ! ---- DOWNSCALED-PARTICLE CACHE ----
        ! Every EM iteration re-reads and re-preps the SAME particles, and one qsys round is
        ! launched per iteration, so the workers are fresh processes each time and nothing held in
        ! memory survives them. The on-disk cache does survive: it stores the iteration-independent
        ! prefix (noise normalisation, FFT, crop to box_crop), which at box=360/box_crop=64 is where
        ! the measured 48% of probe time goes. Masking keeps the full box.
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
            fit%z = 0.d0
            do ibatch = 1, fit%npp, MAXIMGBATCHSZ
                batchlims = [ibatch, min(fit%npp, ibatch + MAXIMGBATCHSZ - 1)]
                batchsz   = batchlims(2) - batchlims(1) + 1
                call planes_batch_load(params, build, fit%npp, fit%ppinds, batchlims, fpls, &
                    &cov_image_mask_radius(params), l_pcache, sec_read, sec_prep)
                do i = 1, batchsz
                    call build%spproj_field%get_ori(fit%ppinds(batchlims(1)+i-1), orientations(i))
                    eo(i) = build%spproj_field%get_eo(fit%ppinds(batchlims(1)+i-1))
                end do
                ! ---- POLAR E-STEP BANK, built once per EM iteration at the first prepped batch
                ! (basis_recs change every M-step, so the bank cannot be cached across iterations;
                ! the grid and the pose-fixed direction assignment are built once per stage) ----
                if( fit%l_pol_es .and. .not. fit%l_pol_bank_it )then
                    t_bank = tic()
                    call fit_polar_bank_build(params, build, fit, mean_rec, fpls(1), nthr)
                    fit%sec_bank      = fit%sec_bank + real(toc(t_bank))
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
                !$omp parallel do default(shared) schedule(dynamic) proc_bind(close) &
                !$omp& private(i,ithr,q,a,aa,e_mm,myv,row,twp0)
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
                sec_estep = sec_estep + toc(t_sec)
                t_sec = tic()
                ! M-step by halfset: Y_q += sum_i z_iq * backproject(r_i), and the coupled normal matrix
                ! rho(q,r) += sum_i |CTF|^2 E[z_iq z_ir]   (batched KB)
                call fit_batch_insert(build, fit, orientations, fpls, eo, batchsz)
                sec_ins = sec_ins + toc(t_sec)
                if( batchlims(2)==fit%nptcls .or. mod(batchlims(2), 5*MAXIMGBATCHSZ)==0 )then
                    write(logfhandle,'(A,I0,A,I0)') '>>> FLEX_PCA PROBE PASS PARTICLES: ',batchlims(2),' / ',fit%npp
                    call flush(logfhandle)
                endif
            end do
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
            call cleanup_rec_buffers(build, fpls)
            deallocate(orientations, eo)
            if( fit%l_converged )then
                write(logfhandle,'(A,I0,A,F9.6)') '>>> FLEX_PCA PROBE converged after ',it_eff, &
                    &' iterations: rank-1 basis cosine vs previous >= ',fit%conv_thresh
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
        endif
        deallocate(fit%z, fit%ppinds)
        ! hoisted model handles back to the caller's arguments
        call move_alloc(fit%basis_recs, basis_recs)
        call move_alloc(fit%eigvals,    eigvals)
        ncomp = fit%ncomp
    end subroutine probe_subspace_iteration

    !> Per-fit stage defaults, all hard-wired: polar E-step, MCFA prior (K=COV_EM_MIX), mean-shaped
    !! deflation, fixed per-particle contrast (no ECM, no a-scaled M-step), band and convergence settings.
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
        ! MCFA: one basis, K latent Gaussians (Baek et al., IEEE TPAMI 2010); K=1 is plain PPCA EM.
        fit%kmix   = COV_EM_MIX
        fit%l_mix_req = fit%kmix >= 1
        ! the MCFA accumulators and a bounded latent subsample (mcfa_init's seed) ride in the probe part files
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
        ! mean-shaped deflation: the per-particle scalar contrast cannot absorb a frequency-dependent
        ! scale, which would otherwise take a whole consensus-shaped component
        fit%vdfl           = COV_EM_DEFLATE
        fit%l_deflate_mean = .true.
        ! per-particle contrast a_i = <m,y>/||m||^2 (clamped to [0.1,5]) is held fixed: no ECM
        ! refinement against the basis, no a-scaled M-step insertion
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

    !> Paired engine: two resident probe fits over the mod-4 halves of the selection, advanced by one
    !! shared loop from identical data-free bases, then merged by probe_paired_merge. Fit A writes
    !! flex_pca_pc*/flex_pca_probe.txt, fit B flex_pca_fitB_pc*/flex_pca_probe_fitB.txt.
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
        integer :: vpair, f, i, cnt, nhalf(2), u
        if( .not. (present(m_basis_recs) .and. present(m_eigvals) .and. present(m_ncomp) &
            &.and. present(m_sig2)) ) THROW_HARD('the paired merge requires the merged-product arguments')
        if( params%n_probe_iters < 2 ) THROW_HARD('the paired merge needs n_probe_iters >= 2: the final iteration''s statistics are expressed in the previous iteration''s delivered frame')
        if( present(m_ncomp) ) m_ncomp = 0
        if( present(m_sig2)  ) m_sig2  = 0.d0
        ! ---- the split: the ONE rule (flex_pca_half_of), shared with the paired workers ----
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
        call probe_subspace_paired(params, build, fits, params%n_probe_iters, vpair, rounds=rounds)
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
        ! ---- final stage: frame-align, accumulator merge and one joint solve; the merged eigenvolumes,
        ! meta and manifest lines are an additional product next to the fits' own delivery ----
        call probe_paired_merge(params, build, fits, m_basis_recs, &
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
    !! polar former and the plain Cartesian former (plan §3.5); the merged read list assumes
    !! ascending ppinds order (hazard 1).
    subroutine probe_subspace_paired( params, build, fits, niters, vpair , rounds)
        class(flex_pca_rounds), intent(inout) :: rounds
        class(parameters), intent(inout) :: params
        type(builder),     intent(inout) :: build
        type(probe_fit_t), intent(inout) :: fits(2)
        integer,           intent(in)    :: niters
        !> the mod-4 pairing, stamped into both probe-state files for distributed dispatch
        integer,           intent(in)    :: vpair
        integer  :: it_eff, niters_eff, f, nthr
        integer(timer_int_kind) :: t_it
        !> cross-fit-FSC driver context: the paired master writes paired=1 records every
        !! iteration; the ridge is their only consumer
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
        ! ---- per-fit stage config (hard-wired defaults) ----
        do f = 1, 2
            call fit_stage_config(params, fits(f), nthr)
        end do
        ! ---- CROSS-FIT-FSC setup (spec par.7): paired master is always a writer ----
        call xfsc_setup(xfctx, params, fits(1)%kfr_ann, .true., .true.)
        ! ---- refusals and stand-downs (plan §3.2/§3.5), printed once at driver start ----
        write(logfhandle,'(A,I0,A,I0,A)') '>>> FLEX_PCA PAIRED: basis dim A=',fits(1)%ncomp, &
            &' B=',fits(2)%ncomp,'  CPU+POLAR E-step'
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
                call probe_fit_merge_stash(fits(f))
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
