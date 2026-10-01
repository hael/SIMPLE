!@descr: flex_pca: the fit drivers: the paired master run (two resident fits over the mod-4 halves, merge) and the worker stage bodies
module simple_flex_pca_fit_driver
use simple_core_module_api
use simple_flex_pca_records, only: flex_fit_model, flex_selection, flex_latent
use simple_builder, only: builder
use simple_parameters, only: parameters
use simple_reconstructor, only: reconstructor
use simple_image, only: image
use simple_flex_pca_rounds, only: flex_pca_rounds
use simple_flex_pca_stages, only: flex_pca_half_of, PCA_STAGE_POLISH, FLEX_FIT_A
use simple_flex_pca_run_types, only: flex_run_settings
use simple_flex_pca_basis, only: apply_cached_mean_scale, estimate_mean_scale, init_mean_reconstructor,&
    &init_basis_datafree, save_probe_state, load_probe_state, load_probe_basis, COV_PROBE_META
use simple_flex_pca_embed, only: embed_latents_with_contrast
use simple_flex_pca_util, only: cov_stage_subsample
use simple_flex_probe_fit, only: flex_probe_fit, flex_mean_ref, probe_subspace_iteration, fit_engine_iterate
use simple_flex_pca_pairmerge, only: probe_paired_merge
implicit none
private
#include "simple_local_flags.inc"

public :: run_flex_pca_paired, probe_worker_pass, embed_worker_pass

contains

    !> PAIRED-ENGINE ENTRY (SIMPLE_COV_PAIRED=1; plan step 1-2, shared-memory probe-only v1):
    !! two resident probe fits over the mod-4 disjoint halves of the master's selection,
    !! advanced by ONE shared master loop (fit_engine_iterate with nfits=2). Per-fit initialisation mirrors the single-fit order:
    !! mean copy + mean scale on the fit's own half, deterministic data-free basis (identical
    !! geometry both fits) with per-fit noise/prior calibration, probe-stage subsample.
    !! Delivery is probe-only: fit A keeps the legacy namespaces (flex_pca_pc*/flex_pca_probe.txt)
    !! so every downstream consumer keeps working; fit B writes flex_pca_fitB_pc* +
    !! flex_pca_probe_fitB.txt (naming precedent: the bagA/bagB pools).
    subroutine run_flex_pca_paired( params, cfg, build, sel, col_sep, neigs_req, model , rounds)
        type(flex_selection), intent(in)    :: sel
        type(flex_fit_model), intent(inout) :: model   !< the merged product: basis, prior variances, rank, noise level
        class(flex_pca_rounds), intent(inout) :: rounds
        class(parameters),   intent(inout) :: params
        type(flex_run_settings), intent(in)    :: cfg
        type(builder),       intent(inout) :: build
        integer,             intent(in) :: col_sep, neigs_req
        !> per merged component: the cross-half match |cos| (the axis weight below)
        real(dp), allocatable :: mergecos(:)
        type(flex_probe_fit), target :: fits(2)
        type(flex_mean_ref) :: means(2)
        type(flex_selection) :: half
        integer :: vpair, f, i, cnt, nhalf(2), u
        ! ---- FINAL STAGE (par.7 "merge, don't refit"): the accumulator-level merge is the delivery ----
        if( params%n_probe_iters < 2 ) THROW_HARD('the paired merge needs n_probe_iters >= 2: the final iteration''s statistics are expressed in the previous iteration''s delivered frame')
        model%ncomp = 0
        model%sig2_eff  = 0.d0
        ! ---- the split: the ONE rule, shared with the two-job pcafit harness ----
        vpair = cfg%mod4_pairing
        if( vpair == 2 ) THROW_HARD('SIMPLE_COV_MOD4_PAIRING=2 groups same-parity rows: both halves lose one eo class under the row-alternating project eo split. Use pairing 1 or 3.')
        if( vpair /= 1 .and. vpair /= 3 ) THROW_HARD('SIMPLE_COV_MOD4_PAIRING must be 1, 2 or 3')
        do f = 1, 2
            cnt = 0
            do i = 1, sel%nptcls
                if( flex_pca_half_of(sel%pinds(i), vpair) == f ) cnt = cnt + 1
            end do
            nhalf(f) = cnt
            if( cnt < 100 ) THROW_HARD('paired engine: a mod-4 half has fewer than 100 particles')
        end do
        write(logfhandle,'(A,I0,A,I0,A,I0,A)') '>>> FLEX_PCA PAIRED ENGINE: two resident fits, &
            &one shared loop; half A=',nhalf(1),'  half B=',nhalf(2),'  (mod-4 pairing ',vpair,')'
        call flush(logfhandle)
        do f = 1, 2
            allocate(half%pinds(nhalf(f)))
            half%nptcls = nhalf(f)
            cnt = 0
            do i = 1, sel%nptcls
                if( flex_pca_half_of(sel%pinds(i), vpair) == f )then
                    cnt = cnt + 1; half%pinds(cnt) = sel%pinds(i)
                endif
            end do
            if( f == FLEX_FIT_A )then
                call fits(f)%new(cfg, f, 'flex_pca_pc', COV_PROBE_META, half)
            else
                call fits(f)%new(cfg, f, 'flex_pca_fitB_pc', 'flex_pca_probe_fitB.txt', half)
            endif
            call half%kill
        end do
        ! ---- per-fit initialisation, mirroring the single-fit order ----
        do f = 1, 2
            ! per-fit mean copy carrying the per-fit mean scale, fitted on the fit's own half
            ! (the two-job harness scales the mean on its own half, so the paired engine must
            ! too or T0 compares different means). Shared memory: estimate_mean_scale writes
            ! no cache file here (rounds%is_master() is false), so no fname clash.
            call init_mean_reconstructor(params, build, fits(f)%model%mean_rec)
            ! per-fit cache namespace: on the DISTRIBUTED master estimate_mean_scale writes the
            ! fitted radial scale for the workers to apply (per-fit files, or the two fits clash)
            if( f == FLEX_FIT_A )then
                call estimate_mean_scale(params, build, fits(f)%model%mean_rec, fits(f)%spec%sel%pinds, &
                    &fits(f)%spec%sel%nptcls, cache_fname='flex_pca_mean_scale.bin', rounds=rounds)
            else
                call estimate_mean_scale(params, build, fits(f)%model%mean_rec, fits(f)%spec%sel%pinds, &
                    &fits(f)%spec%sel%nptcls, cache_fname='flex_pca_mean_scale_fitB.bin', rounds=rounds)
            endif
            ! deterministic data-free basis; its calibration pass runs on the fit's half ->
            ! per-fit sig2/gamma0; the init eigenvolume write uses the fit's namespace
            call init_basis_datafree(params, cfg, build, fits(f)%model, fits(f)%spec%sel, col_sep, neigs_req, &
                &fprefix=fits(f)%spec%fprefix%to_char(), rounds=rounds)
            fits(f)%model%sig2 = max(fits(f)%model%sig2_eff, DTINY)
            ! probe-stage subsample of the fit's half (master process: no nparts division)
            call cov_stage_subsample(build, fits(f)%spec%sel%pinds, fits(f)%spec%sel%nptcls, 1, cfg%probe_max, 'PROBE', &
                &fits(f)%spec%ppinds, fits(f)%spec%npp)
            allocate(fits(f)%model%z(fits(f)%spec%npp, fits(f)%model%ncomp))
        end do
        ! ---- the shared engine loop: one it_eff advances BOTH fits ----
        do f = 1, 2
            means(f)%p => fits(f)%model%mean_rec
        end do
        call fit_engine_iterate(params, build, fits, means, 2, params%n_probe_iters, 0, 0, .true., rounds)
        ! ---- delivery (v1, probe-only): both fits' probe metas + eigenvalue tables; the
        ! eigenvolumes are already on disk from the last iteration under each fit's prefix ----
        do f = 1, 2
            call save_probe_state(fits(f)%model%ncomp, fits(f)%model%eigvals, fits(f)%model%sig2_eff, &
                &fname=fits(f)%spec%meta_fname%to_char())
            call write_paired_eigen_table(fits(f))
        end do
        ! the paired record (probe-only mode never reaches the covariance manifest)
        call del_file('flex_pca_paired.txt')
        open(newunit=u, file='flex_pca_paired.txt', status='replace', action='write')
        write(u,'(A)') '# paired-engine record: paired  pairing  fit  nptcls  npp  ncomp  sig2'
        do f = 1, 2
            write(u,'(I2,1X,I2,1X,A1,1X,I10,1X,I10,1X,I5,1X,ES16.8)') 1, vpair, &
                &merge('A','B',f==FLEX_FIT_A), fits(f)%spec%sel%nptcls, fits(f)%spec%npp, fits(f)%model%ncomp, &
                &fits(f)%model%sig2_eff
        end do
        close(u)
        write(logfhandle,'(A,I0,A,I0,A,ES11.4,A,ES11.4)') '>>> FLEX_PCA PAIRED delivered: &
            &ncomp A=',fits(1)%model%ncomp,' B=',fits(2)%model%ncomp,'  sig2 A=',fits(1)%model%sig2_eff, &
            &' B=',fits(2)%model%sig2_eff
        call flush(logfhandle)
        ! ---- FINAL STAGE (par.7): frame-align + accumulator merge + ONE joint solve (mode 1)
        ! merged eigenvolumes/meta/manifest written inside.
        ! The fits' own delivery above is untouched -- the merge is an ADDITIONAL product.
        call probe_paired_merge(params, build, fits, model, mergecos)
        if( .not. allocated(model%basis_recs) .or. model%ncomp < 1 ) THROW_HARD('paired merge returned no merged basis')
        ! ---- FSC-DOCTRINE AXIS WEIGHTING: the per-axis cross-half match cosine is an FSC-per-component,
        ! and SIMPLE's house treatment of reproducibility is per-shell WEIGHTING,
        ! never hard truncation -- fsc2optlp doctrine, merged-estimate correction
        ! w = 2c/(1+c) (Rosenthal-Henderson). Axes below the 0.143 information
        ! floor are dropped (nothing is there); everything else is down-weighted
        ! through its prior variance, which handles a smooth no-cliff spectrum
        ! (measured on 10028) the way hard thresholds cannot. 0.5 remains the
        ! REPORTING bar for interpretable axes, logged per component below.
        block
            integer  :: qq
            real(dp) :: cq, wq
            if( allocated(mergecos) )then
                write(logfhandle,'(A)') '>>> FLEX_PCA AXIS RELIABILITY (FSC doctrine): &
                    &component | match c | weight 2c/(1+c) | class'
                do qq = 1, model%ncomp
                    cq = max(0.d0, min(1.d0, mergecos(qq)))
                    if( cq < 0.143d0 )then
                        wq = 0.d0
                    else
                        wq = 2.d0*cq/(1.d0 + cq)
                    endif
                    if( cq >= 0.5d0 )then
                        write(logfhandle,'(A,I3,A,F7.3,A,F7.3,A)') '>>>   ', qq, &
                            &'  c=', real(cq), '  w=', real(wq), '  interpretable (c>=0.5)'
                    else if( cq >= 0.143d0 )then
                        write(logfhandle,'(A,I3,A,F7.3,A,F7.3,A)') '>>>   ', qq, &
                            &'  c=', real(cq), '  w=', real(wq), '  weighted (0.143<=c<0.5)'
                    else
                        write(logfhandle,'(A,I3,A,F7.3,A,F7.3,A)') '>>>   ', qq, &
                            &'  c=', real(cq), '  w=', real(wq), '  DROPPED (c<0.143)'
                    endif
                    model%eigvals(qq) = model%eigvals(qq)*max(wq, 1.d-6)
                end do
                call flush(logfhandle)
            endif
        end block
        ! ---- FULL-SET POLISH: after the merge froze rank, axes and convergence on the
        ! half pair, ONE unmonitored full-selection EM iteration from the merged basis
        ! (the eo-combine-then-polish convention). Depth stays at 1: deeper polish has no
        ! held-out signal and was measured to erode the island (2026-09-01 harness).
        ! Writes under the polished namespace so the per-fit deliveries stay intact.
        block
            integer, parameter :: NPOLISH = 1
            write(logfhandle,'(A,I0,A,I0,A)') '>>> FLEX_PCA MERGE POLISH: ', &
                &NPOLISH,' full-selection EM iteration(s) over ',sel%nptcls, &
                &' particles from the merged basis'
            call flush(logfhandle)
            ! materialize the ENTRY basis under the polished namespace BEFORE the
            ! loop: a distributed round-1 worker loads the current basis from disk
            ! (the stamped prefix), and the merged basis lives only in memory here
            ! -- the merged eigenvolumes on disk carry the same content, so copy
            ! them into the polished names (they are overwritten every iteration
            ! by fit_iter_finish thereafter)
            block
                type(image)  :: vcp
                type(string) :: fsrc, fdst
                integer :: qcp
                call vcp%new([params%box_crop,params%box_crop,params%box_crop], &
                    &params%smpd_crop)
                do qcp = 1, model%ncomp
                    fsrc = 'flex_pca_merged_pc'//int2str_pad(qcp,3)//MRC_EXT
                    fdst = 'flex_pca_polished_pc'//int2str_pad(qcp,3)//MRC_EXT
                    if( file_exists(fsrc) )then
                        call vcp%read(fsrc)
                        call vcp%write(fdst, del_if_exists=.true.)
                    endif
                    call fsrc%kill; call fdst%kill
                end do
                call vcp%kill
            end block
            call probe_subspace_iteration(params, cfg, build, model, sel, NPOLISH, &
                &fprefix='flex_pca_polished_pc', &
                &meta_fname='flex_pca_probe_polished.txt', rounds=rounds)
            call save_probe_state(model%ncomp, model%eigvals, model%sig2_eff, &
                &fname='flex_pca_probe_polished.txt')
        end block
        if( allocated(mergecos) ) deallocate(mergecos)
        do f = 1, 2
            call fits(f)%kill
        end do
    end subroutine run_flex_pca_paired

    !> Per-fit eigenvalue table (the paired-local equivalent of write_covariance_eigenvolumes'
    !! table): fit A keeps the legacy name, fit B gets the fitB namespace.
    subroutine write_paired_eigen_table( fit )
        type(flex_probe_fit), intent(in) :: fit
        character(len=:), allocatable :: fn
        integer :: q, u
        if( fit%spec%id == FLEX_FIT_A )then
            fn = 'flex_pca_eigenvalues.txt'
        else
            fn = 'flex_pca_eigenvalues_fitB.txt'
        endif
        call del_file(fn)
        open(newunit=u, file=fn, status='replace', action='write')
        write(u,'(A)') '# component eigenvalue'
        do q = 1, fit%model%ncomp
            write(u,'(I6,1X,ES20.10)') q, fit%model%eigvals(q)
        end do
        close(u)
    end subroutine write_paired_eigen_table

    !> DISTRIBUTED PAIRED WORKER (plan par.7.2): one relaunch per master iteration. Loads BOTH
    !! fits' metas and bases (per-fit namespaces), splits this worker's fromp/top shard by the
    !! ONE mod-4 rule and runs the engine for one iteration keyed to the master's global stamp:
    !! the same shared E-step pass the shared-memory master runs, then ONE v5 part carrying both
    !! fits' accumulator blocks. No master tail runs here (the engine returns on a worker).
    subroutine run_flex_pca_paired_worker( params, cfg, build, sel, it_stamp, &
        &niters_stamp, vpair , rounds)
        type(flex_selection), intent(in) :: sel
        class(flex_pca_rounds), intent(inout) :: rounds
        class(parameters), intent(inout) :: params
        type(flex_run_settings), intent(in)    :: cfg
        type(builder),     intent(inout) :: build
        integer,           intent(in) :: it_stamp, niters_stamp, vpair
        type(flex_probe_fit), target :: fits(2)
        type(flex_mean_ref) :: means(2)
        type(flex_selection) :: half
        real(dp), allocatable :: ev_load(:)
        integer  :: nhalf(2), f, i, cnt, nc_load, it_eff
        real(dp) :: s2_load
        it_eff = max(1, it_stamp)
        if( vpair /= 1 .and. vpair /= 3 ) THROW_HARD('paired worker: invalid mod-4 pairing stamp')
        do f = 1, 2
            cnt = 0
            do i = 1, sel%nptcls
                if( flex_pca_half_of(sel%pinds(i), vpair) == f ) cnt = cnt + 1
            end do
            nhalf(f) = cnt
            if( cnt < 2 ) THROW_HARD('paired worker: a mod-4 half of this shard is (near) empty')
        end do
        write(logfhandle,'(A,I0,A,I0,A,I0,A,I0)') '>>> FLEX_PCA PAIRED WORKER part=', params%part, &
            &'  it_eff=', it_eff, '  shard half A=', nhalf(1), '  B=', nhalf(2)
        call flush(logfhandle)
        do f = 1, 2
            allocate(half%pinds(nhalf(f)))
            half%nptcls = nhalf(f)
            cnt = 0
            do i = 1, sel%nptcls
                if( flex_pca_half_of(sel%pinds(i), vpair) == f )then
                    cnt = cnt + 1; half%pinds(cnt) = sel%pinds(i)
                endif
            end do
            if( f == FLEX_FIT_A )then
                call fits(f)%new(cfg, f, 'flex_pca_pc', COV_PROBE_META, half)
            else
                call fits(f)%new(cfg, f, 'flex_pca_fitB_pc', 'flex_pca_probe_fitB.txt', half)
            endif
            call half%kill
            ! per-fit mean: deterministic rebuild from vol1 + the MASTER-fitted per-fit radial
            ! scale (a worker must never re-fit the scale on its own shard)
            call init_mean_reconstructor(params, build, fits(f)%model%mean_rec)
            if( f == FLEX_FIT_A )then
                call apply_cached_mean_scale(params, fits(f)%model%mean_rec, &
                    &cache_fname='flex_pca_mean_scale.bin')
            else
                call apply_cached_mean_scale(params, fits(f)%model%mean_rec, &
                    &cache_fname='flex_pca_mean_scale_fitB.bin')
            endif
            ! per-fit basis + meta, refreshed by the master before this round was scheduled
            call load_probe_state(nc_load, ev_load, s2_load, fname=fits(f)%spec%meta_fname%to_char())
            fits(f)%model%ncomp    = nc_load
            call move_alloc(ev_load, fits(f)%model%eigvals)
            fits(f)%model%sig2_eff = s2_load
            fits(f)%model%sig2     = max(s2_load, DTINY)
            call load_probe_basis(params, build, fits(f)%model%ncomp, fits(f)%model%basis_recs, &
                &fprefix=fits(f)%spec%fprefix%to_char())
            ! probe-stage subsample of this shard's half; the budget is a TOTAL across processes,
            ! so a worker divides by nparts (cov_stage_subsample contract)
            call cov_stage_subsample(build, fits(f)%spec%sel%pinds, fits(f)%spec%sel%nptcls, params%nparts, &
                &cfg%probe_max, 'PROBE', fits(f)%spec%ppinds, fits(f)%spec%npp)
            allocate(fits(f)%model%z(fits(f)%spec%npp, fits(f)%model%ncomp))
            means(f)%p => fits(f)%model%mean_rec
        end do
        ! ---- one accumulate iteration keyed to the master's global stamp; the engine writes
        ! the v5 part and returns on a worker ----
        call fit_engine_iterate(params, build, fits, means, 2, 1, it_stamp, niters_stamp, .false., rounds)
        do f = 1, 2
            call fits(f)%kill
        end do
    end subroutine run_flex_pca_paired_worker

    subroutine probe_worker_pass( params, cfg, build, model, sel , rounds)
        type(flex_fit_model), intent(inout) :: model   !< the workers mean in; the basis of this round loaded here
        type(flex_selection), intent(in)    :: sel
        class(flex_pca_rounds), intent(inout) :: rounds
        class(parameters),   intent(inout) :: params
        type(flex_run_settings), intent(in)    :: cfg
        type(builder),       intent(inout) :: build
        integer :: vpair, s
        ! the master refreshed the basis volumes and the probe-state file before scheduling this
        ! round; which_iter keys every iteration schedule inside probe_subspace_iteration off the
        ! master's true iteration (this worker's own loop runs exactly once per relaunch), maxits
        ! is the budget, nfits selects the paired pass, the POLISH stage the polished namespace
        call load_probe_state(model%ncomp, model%eigvals, model%sig2_eff)
        if( params%nfits == 2 )then
            ! paired pass: BOTH fits over this worker's particle list, split by the mod-4 rule
            ! (the pairing is a pinned environment constant, SIMPLE_COV_MOD4_PAIRING), one v5 part
            vpair = cfg%mod4_pairing
            if( allocated(model%eigvals) ) deallocate(model%eigvals)
            call run_flex_pca_paired_worker(params, cfg, build, sel, params%which_iter, &
                &params%maxits, vpair, rounds=rounds)
            return
        endif
        if( params%stage == PCA_STAGE_POLISH )then
            call load_probe_basis(params, build, model%ncomp, model%basis_recs, fprefix='flex_pca_polished_pc')
        else
            call load_probe_basis(params, build, model%ncomp, model%basis_recs)
        endif
        call probe_subspace_iteration(params, cfg, build, model, sel, 1, it_glob=params%which_iter, niters_glob=params%maxits, &
            &rounds=rounds)
        do s = 1, size(model%basis_recs)
            call model%basis_recs(s)%dealloc_rho; call model%basis_recs(s)%kill
        end do
        deallocate(model%basis_recs)
        if( allocated(model%eigvals) ) deallocate(model%eigvals)
    end subroutine probe_worker_pass

    !> Same handoff as the probe worker (basis on disk as flex_pca_pc*.mrc, dimension/prior
    !! variances/noise level in flex_pca_probe.txt) but one round: the basis is final by now.
    subroutine embed_worker_pass( params, build, model, sel , rounds)
        type(flex_fit_model), intent(inout) :: model
        type(flex_selection), intent(in)    :: sel
        type(flex_latent) :: lat
        class(flex_pca_rounds), intent(inout) :: rounds
        class(parameters),   intent(inout) :: params
        type(builder),       intent(inout) :: build
        integer :: s
        call load_probe_state(model%ncomp, model%eigvals, model%sig2_eff)
        call load_probe_basis(params, build, model%ncomp, model%basis_recs)
        call embed_latents_with_contrast(params, build, model, sel, lat, rounds, stats_only=.true., l_zhalf=.false.)
        call lat%kill
        do s = 1, size(model%basis_recs)
            call model%basis_recs(s)%dealloc_rho; call model%basis_recs(s)%kill
        end do
        deallocate(model%basis_recs)
        if( allocated(model%eigvals) ) deallocate(model%eigvals)
    end subroutine embed_worker_pass

end module simple_flex_pca_fit_driver
