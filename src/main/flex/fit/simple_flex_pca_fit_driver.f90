!@descr: flex_pca: the fit drivers: the paired master run (two resident fits over the mod-4 halves, merge) and the worker stage bodies
module simple_flex_pca_fit_driver
use simple_core_module_api, only: del_file, dp, dtiny, file_exists, int2str_pad, logfhandle, mrc_ext, &
    &simple_exception, string
use simple_defs_flex,          only: FLEX_FSC_SIGNAL_THRESHOLD
use simple_flex_pca_records,   only: flex_fit_model, flex_selection, flex_latent
use simple_builder,            only: builder
use simple_parameters,         only: parameters
use simple_reconstructor,      only: reconstructor
use simple_image,              only: image
use simple_flex_pca_rounds,    only: flex_pca_rounds
use simple_flex_pca_stages,    only: flex_pca_half_of, PCA_STAGE_POLISH, FLEX_FIT_A, FLEX_MOD4_PAIRING
use simple_flex_pca_pcg,       only: flex_pcg_environment
use simple_flex_pca_basis,     only: apply_cached_mean_scale, estimate_mean_scale, init_mean_reconstructor, &
    &init_basis_datafree, save_probe_state, load_probe_state, load_probe_basis, COV_PROBE_META
use simple_flex_pca_embed,     only: embed_latents_with_contrast
use simple_flex_pca_util,      only: cov_stage_subsample
use simple_flex_probe_fit,     only: flex_probe_fit, flex_mean_ref, probe_subspace_iteration, fit_engine_iterate
use simple_flex_pca_pairmerge, only: probe_paired_merge
use simple_flex_pca_planes,    only: flex_plane_store
implicit none
private
#include "simple_local_flags.inc"

public :: run_flex_pca_paired, probe_worker_pass, embed_worker_pass

contains

    !> Paired engine: two resident probe fits over the mod-4 halves of the selection, advanced by one
    !! shared loop (fit_engine_iterate with nfits=2) from identical data-free bases, then merged by
    !! probe_paired_merge. Fit A writes flex_pca_pc*/flex_pca_probe.txt, fit B flex_pca_fitB_pc*/flex_pca_probe_fitB.txt.
    subroutine run_flex_pca_paired( params, build, plane_store, pcg_env, sel, col_sep, neigs_req, model, rounds )
        type(flex_selection),   intent(in)    :: sel
        type(flex_fit_model),   intent(inout) :: model   !< the merged product: basis, prior variances, rank, noise level
        class(flex_pca_rounds),  intent(inout) :: rounds
        class(parameters),       intent(inout) :: params
        type(builder),           intent(inout) :: build
        class(flex_plane_store), intent(inout) :: plane_store
        class(flex_pcg_environment), intent(in) :: pcg_env
        integer,                 intent(in)    :: col_sep, neigs_req
        !> per merged component: the cross-half match |cos| (the axis weight below)
        real(dp), allocatable :: mergecos(:)
        type(flex_probe_fit), target :: fits(2)
        type(flex_mean_ref) :: means(2)
        type(flex_selection) :: half
        integer :: f, i, cnt, nhalf(2), u
        ! ---- FINAL STAGE (par.7 "merge, don't refit"): the accumulator-level merge is the delivery ----
        if( params%n_probe_iters < 2 ) THROW_HARD('the paired merge needs n_probe_iters >= 2: the final iteration''s statistics are expressed in the previous iteration''s delivered frame')
        model%ncomp = 0
        model%sig2_eff  = 0.d0
        ! ---- the split: the ONE rule (flex_pca_half_of), shared with the paired workers ----
        do f = 1, 2
            cnt = 0
            do i = 1, sel%nptcls
                if( flex_pca_half_of(sel%pinds(i)) == f ) cnt = cnt + 1
            end do
            nhalf(f) = cnt
            if( cnt < 100 ) THROW_HARD('paired engine: a mod-4 half has fewer than 100 particles')
        end do
        write(logfhandle,'(A,I0,A,I0,A,I0,A)') '>>> FLEX_PCA PAIRED ENGINE: two resident fits, &
            &one shared loop; half A=',nhalf(1),'  half B=',nhalf(2),'  (mod-4 pairing ',FLEX_MOD4_PAIRING,')'
        call flush(logfhandle)
        do f = 1, 2
            allocate(half%pinds(nhalf(f)))
            half%nptcls = nhalf(f)
            cnt = 0
            do i = 1, sel%nptcls
                if( flex_pca_half_of(sel%pinds(i)) == f )then
                    cnt = cnt + 1; half%pinds(cnt) = sel%pinds(i)
                endif
            end do
            if( f == FLEX_FIT_A )then
                call fits(f)%new(f, 'flex_pca_pc', COV_PROBE_META, half)
            else
                call fits(f)%new(f, 'flex_pca_fitB_pc', 'flex_pca_probe_fitB.txt', half)
            endif
            call half%kill
        end do
        call prepare_fit_environments(pcg_env, fits)
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
                call estimate_mean_scale(params, build, plane_store, fits(f)%model%mean_rec, fits(f)%spec%sel%pinds, &
                    &fits(f)%spec%sel%nptcls, cache_fname='flex_pca_mean_scale.bin', rounds=rounds)
            else
                call estimate_mean_scale(params, build, plane_store, fits(f)%model%mean_rec, fits(f)%spec%sel%pinds, &
                    &fits(f)%spec%sel%nptcls, cache_fname='flex_pca_mean_scale_fitB.bin', rounds=rounds)
            endif
            ! deterministic data-free basis; its calibration pass runs on the fit's half ->
            ! per-fit sig2/gamma0; the init eigenvolume write uses the fit's namespace
            call init_basis_datafree(params, build, fits(f)%mstep%env, fits(f)%model, fits(f)%spec%sel, col_sep, neigs_req, &
                &fprefix=fits(f)%spec%fprefix%to_char(), rounds=rounds)
            fits(f)%model%sig2 = max(fits(f)%model%sig2_eff, DTINY)
            ! probe-stage subsample of the fit's half (master process: no nparts division)
            call cov_stage_subsample(build%spproj_field, fits(f)%spec%sel%pinds, fits(f)%spec%sel%nptcls, 1, 0, 'PROBE', &
                &fits(f)%spec%ppinds, fits(f)%spec%npp)
            allocate(fits(f)%model%z(fits(f)%spec%npp, fits(f)%model%ncomp))
        end do
        ! ---- the shared engine loop: one it_eff advances BOTH fits ----
        do f = 1, 2
            means(f)%p => fits(f)%model%mean_rec
        end do
        call fit_engine_iterate(params, build, plane_store, fits, means, 2, params%n_probe_iters, 0, 0, .true., rounds)
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
            write(u,'(I2,1X,I2,1X,A1,1X,I10,1X,I10,1X,I5,1X,ES16.8)') 1, FLEX_MOD4_PAIRING, &
                &merge('A','B',f==FLEX_FIT_A), fits(f)%spec%sel%nptcls, fits(f)%spec%npp, fits(f)%model%ncomp, &
                &fits(f)%model%sig2_eff
        end do
        close(u)
        write(logfhandle,'(A,I0,A,I0,A,ES11.4,A,ES11.4)') '>>> FLEX_PCA PAIRED delivered: &
            &ncomp A=',fits(1)%model%ncomp,' B=',fits(2)%model%ncomp,'  sig2 A=',fits(1)%model%sig2_eff, &
            &' B=',fits(2)%model%sig2_eff
        call flush(logfhandle)
        ! ---- final stage: frame-align, accumulator merge and one joint solve; the merged eigenvolumes,
        ! meta and manifest lines are an additional product next to the fits' own delivery ----
        call probe_paired_merge(params, build, fits, model, mergecos)
        if( .not. allocated(model%basis_recs) .or. model%ncomp < 1 ) THROW_HARD('paired merge returned no merged basis')
        ! Axis weight from cross-half match cosine c: prior variance x 2c/(1+c), ~0 below 0.143.
        block
            integer  :: qq
            real(dp) :: cq, wq
            if( allocated(mergecos) )then
                write(logfhandle,'(A)') '>>> FLEX_PCA AXIS RELIABILITY (FSC doctrine): &
                    &component | match c | weight 2c/(1+c) | class'
                do qq = 1, model%ncomp
                    cq = max(0.d0, min(1.d0, mergecos(qq)))
                    if( cq < FLEX_FSC_SIGNAL_THRESHOLD )then
                        wq = 0.d0
                    else
                        wq = 2.d0*cq/(1.d0 + cq)
                    endif
                    if( cq >= 0.5d0 )then
                        write(logfhandle,'(A,I3,A,F7.3,A,F7.3,A)') '>>>   ', qq, &
                            &'  c=', real(cq), '  w=', real(wq), '  interpretable (c>=0.5)'
                    else if( cq >= FLEX_FSC_SIGNAL_THRESHOLD )then
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
        ! Full-set polish: one unmonitored full-selection EM iteration from the merged basis, written
        ! under the polished namespace so the per-fit deliveries stay intact.
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
            call probe_subspace_iteration(params, build, plane_store, pcg_env, model, sel, NPOLISH, &
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
    !! the same shared E-step pass the shared-memory master runs, then one shared-codec part carrying both
    !! fits' accumulator blocks. No master tail runs here (the engine returns on a worker).
    subroutine run_flex_pca_paired_worker( params, build, plane_store, pcg_env, sel, it_stamp, niters_stamp, rounds )
        type(flex_selection),   intent(in)    :: sel
        class(flex_pca_rounds),  intent(inout) :: rounds
        class(parameters),       intent(inout) :: params
        type(builder),           intent(inout) :: build
        class(flex_plane_store), intent(inout) :: plane_store
        class(flex_pcg_environment), intent(in) :: pcg_env
        integer,                 intent(in)    :: it_stamp, niters_stamp
        type(flex_probe_fit), target :: fits(2)
        type(flex_mean_ref) :: means(2)
        type(flex_selection) :: half
        real(dp), allocatable :: ev_load(:)
        integer  :: nhalf(2), f, i, cnt, nc_load, it_eff
        real(dp) :: s2_load
        it_eff = max(1, it_stamp)
        do f = 1, 2
            cnt = 0
            do i = 1, sel%nptcls
                if( flex_pca_half_of(sel%pinds(i)) == f ) cnt = cnt + 1
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
                if( flex_pca_half_of(sel%pinds(i)) == f )then
                    cnt = cnt + 1; half%pinds(cnt) = sel%pinds(i)
                endif
            end do
            if( f == FLEX_FIT_A )then
                call fits(f)%new(f, 'flex_pca_pc', COV_PROBE_META, half)
            else
                call fits(f)%new(f, 'flex_pca_fitB_pc', 'flex_pca_probe_fitB.txt', half)
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
            ! the worker probe pass uses every particle in this shard's half
            call cov_stage_subsample(build%spproj_field, fits(f)%spec%sel%pinds, fits(f)%spec%sel%nptcls, 1, &
                &0, 'PROBE', fits(f)%spec%ppinds, fits(f)%spec%npp)
            allocate(fits(f)%model%z(fits(f)%spec%npp, fits(f)%model%ncomp))
            means(f)%p => fits(f)%model%mean_rec
        end do
        call prepare_fit_environments(pcg_env, fits)
        ! ---- one accumulate iteration keyed to the master's global stamp; the engine writes
        ! the paired part and returns on a worker ----
        call fit_engine_iterate(params, build, plane_store, fits, means, 2, 1, it_stamp, niters_stamp, .false., rounds)
        do f = 1, 2
            call fits(f)%kill
        end do
    end subroutine run_flex_pca_paired_worker

    subroutine probe_worker_pass( params, build, plane_store, pcg_env, model, sel, rounds )
        type(flex_fit_model),   intent(inout) :: model   !< the workers mean in; the basis of this round loaded here
        type(flex_selection),   intent(in)    :: sel
        class(flex_pca_rounds),  intent(inout) :: rounds
        class(parameters),       intent(inout) :: params
        type(builder),           intent(inout) :: build
        class(flex_plane_store), intent(inout) :: plane_store
        class(flex_pcg_environment), intent(in) :: pcg_env
        integer :: s
        ! the master refreshed the basis volumes and the probe-state file before scheduling this
        ! round; which_iter keys every iteration schedule inside probe_subspace_iteration off the
        ! master's true iteration (this worker's own loop runs exactly once per relaunch), maxits
        ! is the budget, nfits selects the paired pass, the POLISH stage the polished namespace
        call load_probe_state(model%ncomp, model%eigvals, model%sig2_eff)
        if( params%nfits == 2 )then
            ! paired pass: BOTH fits over this worker's particle list, split by the mod-4 rule
            ! using the same pinned balanced split as the master, one paired part
            if( allocated(model%eigvals) ) deallocate(model%eigvals)
            call run_flex_pca_paired_worker(params, build, plane_store, pcg_env, sel, params%which_iter, &
                &params%maxits, rounds=rounds)
            return
        endif
        if( params%stage == PCA_STAGE_POLISH )then
            call load_probe_basis(params, build, model%ncomp, model%basis_recs, fprefix='flex_pca_polished_pc')
        else
            call load_probe_basis(params, build, model%ncomp, model%basis_recs)
        endif
        call probe_subspace_iteration(params, build, plane_store, pcg_env, model, sel, 1, &
            &it_glob=params%which_iter, niters_glob=params%maxits, rounds=rounds)
        do s = 1, size(model%basis_recs)
            call model%basis_recs(s)%dealloc_rho; call model%basis_recs(s)%kill
        end do
        deallocate(model%basis_recs)
        if( allocated(model%eigvals) ) deallocate(model%eigvals)
    end subroutine probe_worker_pass

    !> Give every resident fit an exact copy of the run's already-resampled support.
    subroutine prepare_fit_environments( pcg_env, fits )
        class(flex_pcg_environment), intent(in)    :: pcg_env
        type(flex_probe_fit),        intent(inout) :: fits(:)
        integer :: f
        do f = 1, size(fits)
            if( .not. allocated(fits(f)%mstep%env) ) allocate(fits(f)%mstep%env)
            call fits(f)%mstep%env%copy_from(pcg_env)
        end do
    end subroutine prepare_fit_environments

    !> Same handoff as the probe worker (basis on disk as flex_pca_pc*.mrc, dimension/prior
    !! variances/noise level in flex_pca_probe.txt) but one round: the basis is final by now.
    subroutine embed_worker_pass( params, build, plane_store, model, sel , rounds)
        type(flex_fit_model),   intent(inout) :: model
        type(flex_selection),   intent(in)    :: sel
        type(flex_latent) :: lat
        class(flex_pca_rounds),  intent(inout) :: rounds
        class(parameters),       intent(inout) :: params
        type(builder),           intent(inout) :: build
        class(flex_plane_store), intent(inout) :: plane_store
        integer :: s
        call load_probe_state(model%ncomp, model%eigvals, model%sig2_eff)
        call load_probe_basis(params, build, model%ncomp, model%basis_recs)
        call embed_latents_with_contrast(params, build, plane_store, model, sel, lat, rounds, stats_only=.true., l_zhalf=.false.)
        call lat%kill
        do s = 1, size(model%basis_recs)
            call model%basis_recs(s)%dealloc_rho; call model%basis_recs(s)%kill
        end do
        deallocate(model%basis_recs)
        if( allocated(model%eigvals) ) deallocate(model%eigvals)
    end subroutine embed_worker_pass

end module simple_flex_pca_fit_driver
