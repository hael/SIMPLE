!@descr: Standalone projection-aware low-rank covariance workflow for heterogeneous SPA data
module simple_flex_pca_model
use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
use simple_core_module_api
use simple_builder,                    only: builder
use simple_cmdline,                    only: cmdline
use simple_flex_pca_rec3D,         only: reconstruct_flex_weighted_states, flex_rec_box, flex_rec_smpd, &
    &read_state_weights_round
use simple_flex_pca_planes,        only: planes_enable, planes_kill
use simple_flex_pca_plane_cache,   only: plane_cache_ensure, plane_cache_in_use, plane_cache_close
use simple_flex_pca_deconv,        only: calibrate_noise_scale, deconvolve_latent
use simple_flex_pca_rounds,        only: flex_pca_rounds, flex_pca_rounds_shmem, flex_pca_half_of, &
    &PCA_STAGE_STATES, PCA_STAGE_EMBED, PCA_STAGE_PROBE, PCA_STAGE_POLISH, FLEX_FIT_ALL, FLEX_FIT_A, FLEX_FIT_B
use simple_flex_pca_em,         only: cov_env_int_pub, compose_basis_from_runs, compose_cut_reembed, &
    &build_covariance_eigenbasis, embed_latents_with_contrast, probe_worker_pass, embed_worker_pass, &
    &estimate_covariance_mean, probe_subspace_iteration, align_basis_to_reference, &
    &save_probe_state, run_flex_pca_paired, &
    &init_basis_reconstructor
use simple_image,                      only: image
use simple_parameters,                 only: parameters
use simple_reconstructor,              only: reconstructor
use simple_sigma2_files,               only: load_sigma2_groups
use simple_sp_project,                 only: sp_project
use simple_srch_sort_loc,              only: hpsort
use simple_umap,                       only: umap_embed, umap_subsample
use simple_flex_pca_plot,              only: flex_plot_latent_jpg
use simple_finch,                      only: finch_hierarchy, fit_finch, finch_representatives, &
    &                                          select_finch_level, refine_finch_level
use simple_kd_tree,                    only: kd_tree, knn_table
use simple_linalg,                     only: jacobi, eigsrt, matinv
use simple_flex_pca_merge,             only: flex_pca_merge_enabled, two_gate_state_merge, &
    &flex_pca_merge_force_on
use simple_flex_pca_util,  only: cov_env_flag_on, cov_env_flag_off, cov_env_dp
use simple_flex_pca_weights, only: build_covariance_state_weights, cv_select_bandwidths
use simple_flex_pca_targets, only: component_reliability_proxy
use simple_flex_weights_state, only: flex_weights_deliver, flex_weights_state_fname, FLEX_WEIGHTS_STALE_SCAN
use simple_flex_weights_file,  only: FLEX_WEIGHTS_PROV_FLEX_PCA, FLEX_WEIGHTS_PROV_MERGED
implicit none
private
character(len=*), parameter :: SIGMA_STATE_FNAME= 'flex_pca_sigma_state.txt'
#include "simple_local_flags.inc"

public :: run_flex_pca, run_flex_pca_worker
public :: test_flex_pca_embedding_cache_io
public :: test_flex_pca_auto_settings
public :: test_flex_pca_population_floor
public :: auto_box_crop, auto_min_neff, auto_state_count

!> Over-provisioning level for npreimages=0. Bounded by cost, not accuracy: 24 and 32 converge to
!! the same answer, while gate 2 compares K(K-1)/2 map pairs.
integer, parameter :: FLEX_AUTO_K_START = 32
integer, parameter :: FLEX_AUTO_K_MIN   = 8
!> Nyquist margin for a derived box_crop; columns are selected inside that band.
real,    parameter :: FLEX_AUTO_BOX_SAFETY = 1.25
!> Occupancy share of an equal split. Paired with an SNR term because neither alone reproduces both
!! validation datasets.
real,    parameter :: FLEX_AUTO_NEFF_OCCUPANCY = 0.10

character(len=8), parameter :: COV_CACHE_MAGIC   = 'SIMPLFXC'
! Bump whenever the cache layout changes; read_embedding_cache rejects any other version.
integer,          parameter :: COV_CACHE_VERSION = 4   ! 4: box_crop/smpd_crop provenance after (nptcls, ncomp)
!> optional trailing block of the SAME file: the deconvolved coordinates and precisions, the mixture
!! labels and the noise scale; a resume adopts it from infile, so no second cache file and no K ladder
character(len=8), parameter :: COV_DECONV_MAGIC  = 'SIMPLFXD'
!> the run's working lattice, set at the top of run_flex_pca; the cache I/O stamps and checks it
integer :: cov_box_crop_glob  = 0
!> UMAP coordinates of the last delivery readout (and their particle indices), kept so the state
!! stage can draw the same embedding coloured by the delivered states (flex_pca_umap_states.jpg)
real,    allocatable :: umap_plot_xy(:,:)
integer, allocatable :: umap_plot_pind(:)
real    :: cov_smpd_crop_glob = 0.
! Safety cap on the bandwidth widening loop.
integer,          parameter :: COV_MAX_BW_GROW   = 4
integer,          parameter :: MIN_NSTATES       = 3
! npreimages is a PROVISION CEILING, not a target: state placement lays down that many kernels and
! the two-gate merge collapses the indistinct ones, so the recovered K is only ever <= it.
! preimage_auto=yes raises that ceiling to AUTO_NSTATES and turns the merge on, since over-provisioning
! is the only regime in which the merge can recover K at all. Measured on Ribosembly: ceiling 32 with
! the complete-linkage merge recovered 14 states at ARI 0.947 (14/16 GT states covered).
integer,          parameter :: AUTO_NSTATES      = 32

contains

    subroutine run_flex_pca( params, build, cline , rounds)
        use simple_qsys_funs, only: qsys_job_finished
        class(flex_pca_rounds), intent(inout) :: rounds
        type(parameters), intent(inout) :: params
        type(builder),    intent(inout) :: build
        class(cmdline),   intent(inout) :: cline
        type(reconstructor) :: mean_rec
        type(reconstructor), allocatable :: basis_recs(:)
        integer, allocatable :: pinds(:), labels(:)
        real(dp), allocatable :: z(:,:), eigvals(:), prior_precision(:), contrast(:)
        real(dp), allocatable :: latent_second(:,:,:)
        real(dp), allocatable :: zhalf(:,:,:)
        integer,  allocatable :: deconv_labels(:)
        logical  :: l_deconv_applied, l_deconv_adopted, l_state_rec
        logical  :: l_pop_floor
        logical  :: l_merged          ! the two-gate merge collapsed states: the delivered table is a merged one
        logical  :: l_paired_states   ! paired merge delivered the basis+embedding; fall through to the state stage
        real(dp), allocatable :: resid_energy(:), resid_mean_energy(:)
        real, allocatable :: state_weights(:,:), half_weights(:,:), targets(:,:), bandwidths(:), neff(:)
        integer :: nptcls, ncomp, nstates, min_neff, state_axis, col_sep, neigs_req, nkern
        integer :: q, i, r, s
        integer :: nstates_merged
        integer, allocatable :: merge_label(:)
        real,    allocatable :: merged_weights(:,:), merged_targets(:,:), merged_bw(:)
        real(dp), allocatable :: state_mass(:), merged_mass(:)
        real(dp) :: sumw_s, sumw2_s
        real(dp), allocatable :: kdist(:,:), kfloor(:)
        real(dp), allocatable :: comp_rho(:)     ! per-component reliability, drives state-target ordering
        ! ---- B1: responsibility-delivered states from the probe-fitted mixture ----
        real(dp), allocatable :: pviews(:,:)     ! per-particle viewing AXIS, for the GMM orientation term
        character(len=:), allocatable :: cachedir, cachestr
        real(dp) :: sig2_eff
        logical :: sigma_loaded, l_resume
        logical :: l_compose
        character(len=XLONGSTRLEN) :: envc
        integer :: envlen_c, envstat_c
        integer(timer_int_kind) :: t_blk

        call validate_covariance_inputs(params, build, cline, pinds, nptcls, rounds=rounds)
        l_paired_states = .false.
        l_pop_floor     = .false.
        l_merged        = .false.
        cov_box_crop_glob  = params%box_crop
        cov_smpd_crop_glob = params%smpd_crop
        ! distributed master: one particle-index list per part, shipped as pindfile= (no-op otherwise)
        call rounds%plan_partitions(params, pinds)
        call load_and_validate_sigma(params, build, cline, pinds, sigma_loaded, rounds=rounds)
        ! a shared-memory run keeps its prepped planes across every pass; the distributed master
        ! never runs a pass itself and a worker is a fresh process per round
        ! cache=yes: the plane cache is built once by the process that owns the run (workers adopt it)
        if( .not. rounds%is_worker() ) call plane_cache_ensure(params, build, pinds, nptcls)
        if( .not. rounds%distributed() .and. .not. rounds%is_worker() )then
            call planes_enable(build%spproj_field%get_noris(), plane_cache_in_use(params, build))
        endif

        neigs_req  = max(1, min(48, params%neigs))
        ! No ceiling on the state count: over-provisioning is the only regime in which the merge recovers K.
        ! npreimages=0 selects the automatic ceiling; anything else is taken as the requested ceiling.
        ! npreimages is the ceiling; preimage_auto raises it and enables the collapse that makes a
        ! ceiling meaningful. An explicit npreimages alongside preimage_auto=yes is honoured as the
        ! ceiling -- auto then contributes only the merge.
        nstates = max(MIN_NSTATES, params%npreimages)
        ! population floor (min_state_frac > 0): exactly npreimages hard-labelled states, each above
        ! the floor; it cannot be combined with the automatic ceiling or the two-gate merge
        if( params%min_state_frac > 0. )then
            if( params%l_preimage_auto .or. flex_pca_merge_enabled() )then
                THROW_HARD('min_state_frac delivers exactly npreimages states; it cannot be combined with preimage_auto or the merge')
            endif
        endif
        if( params%l_preimage_auto )then
            if( .not. cline%defined('npreimages') ) nstates = AUTO_NSTATES
            ! a ceiling without the collapse is just a large state count, so auto implies the merge
            call flex_pca_merge_force_on
            write(logfhandle,'(A,I0,A)') '>>> FLEX_PCA preimage_auto=yes: state count is a CEILING of ', &
                &nstates, '; the two-gate merge collapses these to the recovered count'
            if( .not. flex_pca_merge_enabled() ) write(logfhandle,'(A)') &
                &'>>> FLEX_PCA WARNING: SIMPLE_COV_MERGE=0 disables the merge, so the ceiling will be &
                &delivered as-is; drop preimage_auto and set npreimages explicitly instead'
        endif
        call report_state_memory(params, nstates)
        ! the minimum effective sample size of a delivered state: the kernel bandwidth floor AND
        ! the occupancy floor that decides whether a state is reconstructed at all
        min_neff = max(20, min(nptcls, params%min_neff))
        state_axis = params%state_axis      ! <0 path, 0 k-means, >=1 legacy single axis
        ! nkern decouples the number of components the STATE STAGE uses from neigs, the number estimated.
        nkern      = params%nkern
        if( nkern <= 0 ) nkern = ishft(huge(1), -1) ! clamped against ncomp once the fit is known
        col_sep    = max(1, params%column_separation)

        write(logfhandle,'(A,I0,A,I0)') '>>> FLEX_PCA particles=',nptcls, &
            &' requested_components=',neigs_req
        write(logfhandle,'(A,L1,A,I0,A,I0)') '>>> FLEX_PCA sigma_whitened=',sigma_loaded, &
            &' state_axis=',state_axis,' minimum_state_neff=',min_neff
        ! A low-pass finer than the working Nyquist (2*smpd_crop) cannot be honoured
        if( params%lp <= 2.0*params%smpd_crop + TINY )then
            write(logfhandle,'(A,F8.3,A,F8.3,A)') '>>> FLEX_PCA WARNING: lp=',params%lp, &
                &' A is at or beyond the working Nyquist (2*smpd_crop=',2.0*params%smpd_crop, &
                &' A); NO low-pass will be applied and the covariance is estimated to Nyquist.'
            write(logfhandle,'(A)') '>>> FLEX_PCA WARNING: increase lp or decrease box_crop/smpd_crop.'
        else
            write(logfhandle,'(A,F8.3,A,I0)') '>>> FLEX_PCA covariance band: lp=',params%lp, &
                &' A, kstop=',max(1,min(max(1,fdim(params%box_crop)-1), &
                &int(real(max(1,params%box_crop-1))*params%smpd_crop/params%lp)))
        endif
        call flush(logfhandle)

        ! RESUME MODE. The basis and embedding dominate runtime and do not depend on the state stage,
        ! so infile= re-runs only that stage.
        l_resume = cline%defined('infile')
        cachedir = ''
        if( l_resume )then
            call read_embedding_cache(params%infile%to_char(), pinds, nptcls, ncomp, z, eigvals, &
                &contrast, resid_energy, resid_mean_energy, latent_second, sig2_eff)
            allocate(prior_precision(ncomp))
            do q = 1, ncomp
                prior_precision(q) = 1.d0 / max(eigvals(q), DTINY)
            end do
            if( state_axis > 0 ) state_axis = min(state_axis, min(ncomp, nkern))
            cachestr = params%infile%to_char()
            cachedir = ''
            do i = len_trim(cachestr), 1, -1
                if( cachestr(i:i) == '/' )then
                    cachedir = cachestr(1:i)
                    exit
                endif
            end do
            write(logfhandle,'(A,A)') '>>> FLEX_PCA RESUMED from embedding cache ',trim(cachestr)
            write(logfhandle,'(A,I0,A,I0)') '>>> FLEX_PCA cached particles=',nptcls,' components=',ncomp
            call flush(logfhandle)
        else

        call estimate_covariance_mean(params, build, mean_rec, pinds, nptcls, rounds=rounds)

        ! ---- BASIS COMPOSITION (SIMPLE_COV_COMPOSE=<dir>,...): the basis is the coarse-first
        ! orthonormalised union of finished runs' eigenvolumes (see simple_flex_pca_em_compose);
        ! no covariance moment, no EM here -- straight to the embedding
        envc = ''
        call get_environment_variable('SIMPLE_COV_COMPOSE', envc, envlen_c, envstat_c)
        l_compose = envstat_c == 0 .and. envlen_c > 0

        ! ---- THE FIT: two disjoint mod-4-half fits, resident together and advanced by one shared
        ! master loop, then the merge, then ONE all-N embedding against the merged basis (the
        ! combine-then-polish convention: no restart, frozen rank and axes). Placed after
        ! estimate_covariance_mean so the read-once env caches are set exactly as on the compose
        ! path; the full-selection mean_rec built above is not consumed -- each fit owns a mean
        ! scaled on its own half. A distributed WORKER never enters the driver: it falls through to
        ! the worker control flow, routed by the probe-state file's paired dispatch stamp.
        block
            real(dp), allocatable :: mergecos(:)
            if( .not. l_compose .and. .not. rounds%is_worker() )then
                if( rounds%nparts() > 1 )then
                    write(logfhandle,'(A,I0,A)') '>>> FLEX_PCA PAIRED DISTRIBUTED: ', &
                        &rounds%nparts(), ' parts, one part per worker per iteration'
                    call flush(logfhandle)
                endif
                    call run_flex_pca_paired(params, build, pinds, nptcls, col_sep, neigs_req, &
                        &m_basis_recs=basis_recs, m_eigvals=eigvals, m_ncomp=ncomp, m_sig2=sig2_eff, &
                        &m_matchcos=mergecos, rounds=rounds)
                    if( .not. allocated(basis_recs) .or. ncomp < 1 ) &
                        &THROW_HARD('paired merge returned no merged basis')
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
                            do qq = 1, ncomp
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
                                eigvals(qq) = eigvals(qq)*max(wq, 1.d-6)
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
                        if( .true. )then
                            write(logfhandle,'(A,I0,A,I0,A)') '>>> FLEX_PCA MERGE POLISH: ', &
                                &NPOLISH,' full-selection EM iteration(s) over ',nptcls, &
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
                                do qcp = 1, ncomp
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
                            call probe_subspace_iteration(params, build, mean_rec, basis_recs, &
                                &eigvals, sig2_eff, pinds, nptcls, ncomp, NPOLISH, &
                                &fprefix='flex_pca_polished_pc', &
                                &meta_fname='flex_pca_probe_polished.txt', rounds=rounds)
                            call save_probe_state(ncomp, eigvals, sig2_eff, &
                                &fname='flex_pca_probe_polished.txt')
                        endif
                    end block
                    write(logfhandle,'(A,I0,A,I0,A)') '>>> FLEX_PCA MERGED EMBED: all-N pass over ', &
                        &nptcls,' selected particles (both halves) against the ',ncomp, &
                        &'-component merged basis'
                    call flush(logfhandle)
                    allocate(z(nptcls,ncomp), latent_second(ncomp,ncomp,nptcls))
                    allocate(resid_energy(nptcls), resid_mean_energy(nptcls), contrast(nptcls))
                    allocate(zhalf(nptcls,ncomp,2), source=0.d0)
                    write(logfhandle,'(A)') '>>> FLEX_PCA MERGED EMBED: plain prior; axis reliability is &
                        &the merge match cosine already in eigvals (split-half rho not applied)'
                    call flush(logfhandle)
                    call embed_latents_with_contrast(params, build, mean_rec, basis_recs, ncomp, &
                        &eigvals, sig2_eff, pinds, nptcls, z, contrast, latent_second, &
                        &resid_energy, resid_mean_energy, rounds=rounds, &
                        &zhalf_out=zhalf)
                    write(logfhandle,'(A,F7.3,A,F7.3)') '>>> FLEX_PCA per-particle contrast: mean=', &
                        &real(sum(contrast)/real(nptcls,dp)),' sd=', &
                        &real(sqrt(max(sum((contrast-sum(contrast)/nptcls)**2)/real(nptcls,dp),DTINY)))
                    call flush(logfhandle)
                    call write_embedding_cache('flex_pca_embedding.bin', pinds, nptcls, ncomp, z, &
                        &eigvals, contrast, resid_energy, resid_mean_energy, latent_second, sig2_eff)
                    ! raw readouts first (the cache holds raw z), then the deconvolved coordinates are
                    ! what the run delivers and what the state stage re-derives on resume
                    call deliver_latent_readouts(pinds, nptcls, ncomp, z, params%l_umap)
                    ! the merged basis, the raw embedding and its posterior moments stay allocated: the
                    ! run falls through to the common deconvolution + state stage below, exactly what a
                    ! resume from this cache would do (the cache holds raw z)
                    allocate(prior_precision(ncomp))
                    do q = 1, ncomp
                        prior_precision(q) = 1.d0 / max(eigvals(q), DTINY)
                    end do
                    if( state_axis > 0 ) state_axis = min(state_axis, min(ncomp, nkern))
                    l_paired_states = .true.
                    write(logfhandle,'(A,I0,A)') '>>> FLEX_PCA PAIRED MERGE delivered; continuing into &
                        &the state stage with the ', ncomp, '-component merged basis'
                    call flush(logfhandle)
            endif
        end block

        if( .not. l_paired_states )then
        block_eigenbasis: block

            if( l_compose )then
                call compose_basis_from_runs(params, build, basis_recs, eigvals, ncomp, sig2_eff)
                if( .not. cov_env_flag_off('SIMPLE_COV_COMPOSE_CUT') ) call compose_cut_reembed(params, build, &
                    &mean_rec, basis_recs, eigvals, ncomp, sig2_eff, pinds, nptcls, rounds)
                if( state_axis > 0 ) state_axis = min(state_axis, min(ncomp, nkern))
                write(logfhandle,'(A,I0,A)') '>>> FLEX_PCA composed basis: ',ncomp, &
                    &' components; no covariance moment, no EM (n_probe_iters ignored)'
                call flush(logfhandle)
            else
            call build_covariance_eigenbasis(params, build, mean_rec, pinds, nptcls, &
                &col_sep, neigs_req, basis_recs, eigvals, ncomp, sig2_eff, rounds=rounds)
            if( state_axis > 0 ) state_axis = min(state_axis, min(ncomp, nkern))
            write(logfhandle,'(A,I0)') '>>> FLEX_PCA retained covariance components=',ncomp
            call flush(logfhandle)

            ! Optional Wiener E-step / weighted-backprojection M-step, to clean the noisy column directions.
            if( params%n_probe_iters > 0 )then
                call probe_subspace_iteration(params, build, mean_rec, basis_recs, eigvals, sig2_eff, &
                    &pinds, nptcls, ncomp, params%n_probe_iters, rounds=rounds)
                if( state_axis > 0 ) state_axis = min(state_axis, min(ncomp, nkern))
                write(logfhandle,'(A,I0)') '>>> FLEX_PCA probe-refined components=',ncomp
                call flush(logfhandle)
            endif
            endif   ! l_compose
            allocate(z(nptcls,ncomp), prior_precision(ncomp))
            allocate(latent_second(ncomp,ncomp,nptcls))
            allocate(resid_energy(nptcls), resid_mean_energy(nptcls))

            do q = 1, ncomp
                prior_precision(q) = 1.d0 / max(eigvals(q), DTINY)
            end do
            write(logfhandle,'(A,ES12.4,A,ES12.4)') '>>> FLEX_PCA covariance eigenvalues: max=', &
                &maxval(eigvals),' min=',minval(eigvals)
            call flush(logfhandle)
            allocate(contrast(nptcls))
            allocate(comp_rho(ncomp), source=1.d0)
            ! One qsys round: the basis is fixed. Workers cannot finish -- the reliability prior couples every
            ! particle -- so they ship sufficient statistics and the master owns rho and the re-solve.
            allocate(zhalf(nptcls,ncomp,2), source=0.d0)
            if( rounds%distributed() )then
                call save_probe_state(ncomp, eigvals, sig2_eff)
                call rounds%run_stage(params, PCA_STAGE_EMBED, 'embedding')
                call embed_latents_with_contrast(params, build, mean_rec, basis_recs, ncomp, eigvals, &
                    &sig2_eff, pinds, nptcls, z, contrast, latent_second, resid_energy, &
                    &resid_mean_energy, rho_out=comp_rho, from_parts=.true., rounds=rounds, zhalf_out=zhalf)
            else
                call embed_latents_with_contrast(params, build, mean_rec, basis_recs, ncomp, eigvals, &
                    &sig2_eff, pinds, nptcls, z, contrast, latent_second, resid_energy, &
                    &resid_mean_energy, rho_out=comp_rho, rounds=rounds, zhalf_out=zhalf)
            endif
        end block block_eigenbasis
        write(logfhandle,'(A,F7.3,A,F7.3)') '>>> FLEX_PCA per-particle contrast: mean=', &
            &real(sum(contrast)/real(nptcls,dp)),' sd=', &
            &real(sqrt(max(sum((contrast-sum(contrast)/nptcls)**2)/real(nptcls,dp),DTINY)))
        call flush(logfhandle)
        ! z is left in the physical units the MAP solve returns; the kernel metric does the weighting
        call write_covariance_eigenvolumes(basis_recs, eigvals, ncomp)
        call write_embedding_cache('flex_pca_embedding.bin', pinds, nptcls, ncomp, z, eigvals, &
            &contrast, resid_energy, resid_mean_energy, latent_second, sig2_eff)
        call deliver_latent_readouts(pinds, nptcls, ncomp, z, params%l_umap)
        endif   ! .not. l_paired_states

        endif   ! .not. l_resume
        ! Resuming skips the split-half solve, so fall back to the spread-over-posterior-variance proxy.
        if( .not. allocated(comp_rho) )then
            allocate(comp_rho(ncomp))
            call component_reliability_proxy(z, latent_second, nptcls, ncomp, comp_rho)
            write(logfhandle,'(A)') '>>> FLEX_PCA resumed embedding: component reliability from the &
                &spread/posterior-variance proxy (no cached split-half rho)'
        endif
        ! Latent deconvolution. Must stay after the cache write: the cache stores raw z, so a
        ! resume never deconvolves twice.
        call apply_latent_deconvolution(z, latent_second, eigvals, zhalf, pinds, nptcls, ncomp, l_deconv_applied, &
            &labels=deconv_labels, contrast=contrast, resid_energy=resid_energy, &
            &resid_mean_energy=resid_mean_energy, sig2_eff=sig2_eff, resume=l_resume, adopted=l_deconv_adopted, &
            &srcdir=infile_dir(params, l_resume), srcfile=infile_path(params, l_resume))
        if( l_deconv_applied )then
            ! a resume that adopted the cache already has the readouts from the original run
            if( .not. l_deconv_adopted )then
                if( allocated(deconv_labels) )then
                    call deliver_latent_readouts(pinds, nptcls, ncomp, z, params%l_umap, tag='_deconv', labels=deconv_labels)
                else
                    call deliver_latent_readouts(pinds, nptcls, ncomp, z, params%l_umap, tag='_deconv')
                endif
            endif
        endif

        ! Per-particle viewing AXIS, folded antipodally: +n and -n are mirror projections, so orientation
        ! bias lives in the axis and never in the mean resultant. Only the GMM's coverage term reads it.
        allocate(pviews(3,nptcls))
        do i = 1, nptcls
            pviews(:,i) = real(build%spproj_field%get_normal(pinds(i)), dp)
            if( pviews(3,i) < 0.d0 ) pviews(:,i) = -pviews(:,i)
        end do
        ! the deconvolution's mixture components are the macro-clusters (full covariances,
        ! per-particle noise, held-out K); GMM AUTO's tied-covariance discovery only runs
        ! when no deconvolution happened
        if( params%min_state_frac > 0. )then
            ! population floor: exactly nstates hard-labelled states, every one at or above the floor
            ! (refine3D_states flex=yes relies on it); the bandwidth CV below does not apply
            l_pop_floor = .true.
            call place_states_with_population_floor(z, nptcls, ncomp, nkern, nstates, state_axis, min_neff, &
                &params%min_state_frac, eigvals, latent_second, state_weights, targets, bandwidths, neff, &
                &labels, comp_rho=comp_rho)
        else
            call build_covariance_state_weights(z, nptcls, ncomp, nkern, nstates, state_axis, min_neff, &
                &eigvals, latent_second, state_weights, targets, bandwidths, neff, labels, &
                &dist_out=kdist, bfloor_out=kfloor, comp_rho=comp_rho, macro_in=deconv_labels)
        endif

        ! ---- OCCUPANCY FLOOR ----
        ! min_neff is the minimum EFFECTIVE sample size of a state; until 2026-09-12 it only set the
        ! kernel bandwidth, so a seat that ended up with a handful of particles was still reconstructed.
        ! Measured on 10028: two of eight delivered states held 232 and 284 particles and their box-360
        ! maps were reconstruction artifacts (directional striping, a uniform amplitude offset), while
        ! every state above the floor was a clean map. A state below the floor carries no information
        ! its macro-cluster does not already carry, so it is dropped before any map is reconstructed.
        call prune_underpopulated_states(nptcls, nstates, min_neff, state_weights, targets, &
            &bandwidths, neff, labels, kdist, kfloor)

        ! MUST precede cv_select_bandwidths, which reconstructs trial half maps through the same backend.
        params%l_ml_reg = .false.
        params%ml_reg   = 'no'
        call cline%set('ml_reg','no')
        ! The reconstruction backend (prep_imgs4rec) reads its Fourier band from build%esig%get_kfromto().
        call build%esig%set_kfromto([1, max(1, fdim(flex_rec_box(params)) - 1)])
        if( flex_rec_box(params) /= params%box_crop )then
            write(logfhandle,'(A,I0,A,F6.3,A,I0,A,F6.3,A)') '>>> FLEX_PCA state maps decoupled from covariance box: &
                &box_rec=',flex_rec_box(params),' smpd_rec=',flex_rec_smpd(params),' A vs box_crop=',params%box_crop, &
                &' smpd_crop=',params%smpd_crop,' A'
        endif

        if( params%nbins > 1 .and. .not. l_pop_floor )then
            call cv_select_bandwidths(params, build, pinds, nptcls, nstates, params%nbins, min_neff, &
                &kdist, kfloor, state_weights, bandwidths, neff, rounds=rounds)
        endif
        ! rec_states=no: deliver the placement, labels and readouts without reconstructing any
        ! map -- the states stage is then seconds instead of minutes
        l_state_rec = trim(params%rec_states) /= 'no'
        if( .not. l_state_rec )then
            write(logfhandle,'(A)') '>>> FLEX_PCA STATE RECONSTRUCTION OFF (rec_states=no): labels and &
                &tables only, no maps, no merge gate'
            call flush(logfhandle)
        endif
        if( l_state_rec )then
        ! combined states and both halfsets in ONE pass; combined == even + odd exactly
        params%outvol = 'flex_pca_state_001.mrc'
        t_blk = tic()
        ! states of one macro-cluster share a delivery pool (per-shell adaptive kernel in rec3D)
        call reconstruct_flex_weighted_states(params, build, pinds, state_weights, nstates, &
            &floor_rho=.true., outvol_even=string('flex_pca_even_state_001.mrc'), &
            &outvol_odd=string('flex_pca_odd_state_001.mrc'), rounds=rounds)
        write(logfhandle,'(A,F9.1)') '>>> FLEX_PCA STAGE states_combined_eo seconds=', toc(t_blk)
        ! Collapse indistinct states and reconstruct once at the surviving count. Needs the half maps above.
        if( flex_pca_merge_enabled() .and. nstates > 1 )then
            t_blk = tic()
            allocate(merge_label(nstates))
            call two_gate_state_merge(params, pviews, state_weights, nptcls, nstates, &
                &merge_label, nstates_merged)
            if( nstates_merged < nstates )then
                allocate(state_mass(nstates))
                do s = 1, nstates
                    state_mass(s) = sum(real(state_weights(:,s), dp))
                end do
                allocate(merged_weights(nptcls,nstates_merged), source=0.)
                do s = 1, nstates
                    merged_weights(:,merge_label(s)) = merged_weights(:,merge_label(s)) + state_weights(:,s)
                end do
                call move_alloc(merged_weights, state_weights)
                ! RECOMPUTE the label, do not remap it: the merge SUMS columns, and the pre-merge argmax
                ! can lose to a combined rival it beat individually (0.40 vs 0.35+0.25). Remapping would
                ! label the particle to a map it is no longer the largest contributor to, so the hard
                ! assignment and the delivered maps would describe different partitions. 0 is preserved:
                ! summing cannot lift a particle that was outside every kernel support.
                do i = 1, nptcls
                    if( labels(i) >= 1 ) labels(i) = maxloc(state_weights(i,:), dim=1)
                end do
                ! collapse the per-state tables by the same mass, else they still describe the pre-merge nstates
                allocate(merged_targets(size(targets,1),nstates_merged), source=0.)
                allocate(merged_bw(nstates_merged), source=0.)
                allocate(merged_mass(nstates_merged), source=0.d0)
                do s = 1, nstates
                    r = merge_label(s)
                    merged_targets(:,r) = merged_targets(:,r) + real(state_mass(s))*targets(:,s)
                    merged_bw(r)        = merged_bw(r)        + real(state_mass(s))*bandwidths(s)
                    merged_mass(r)      = merged_mass(r)      + state_mass(s)
                end do
                do r = 1, nstates_merged
                    if( merged_mass(r) > DTINY )then
                        merged_targets(:,r) = merged_targets(:,r) / real(merged_mass(r))
                        merged_bw(r)        = merged_bw(r)        / real(merged_mass(r))
                    endif
                end do
                call move_alloc(merged_targets, targets)
                call move_alloc(merged_bw,      bandwidths)
                deallocate(neff); allocate(neff(nstates_merged))
                do r = 1, nstates_merged
                    sumw_s  = sum(real(state_weights(:,r), dp))
                    sumw2_s = sum(real(state_weights(:,r), dp)**2)
                    neff(r) = real(sumw_s*sumw_s / max(sumw2_s, DTINY))
                end do
                deallocate(state_mass, merged_mass)
                ! downstream addresses the maps as a contiguous run, so the stale tail would deliver merged-away states
                do s = nstates_merged + 1, nstates
                    call del_file('flex_pca_state_'     //int2str_pad(s,3)//MRC_EXT)
                    call del_file('flex_pca_even_state_'//int2str_pad(s,3)//MRC_EXT)
                    call del_file('flex_pca_odd_state_' //int2str_pad(s,3)//MRC_EXT)
                end do
                nstates  = nstates_merged
                l_merged = .true.
                call reconstruct_flex_weighted_states(params, build, pinds, state_weights, nstates, &
                    &floor_rho=.true., outvol_even=string('flex_pca_even_state_001.mrc'), &
                    &outvol_odd=string('flex_pca_odd_state_001.mrc'), rounds=rounds)
            endif
            deallocate(merge_label)
            write(logfhandle,'(A,F9.1)') '>>> FLEX_PCA STAGE two_gate_merge seconds=', toc(t_blk)
        endif
        endif   ! l_state_rec
        call write_covariance_tables(build, pinds, z, eigvals, prior_precision, state_weights, labels, &
            &targets, bandwidths, neff, resid_energy, resid_mean_energy, contrast)
        ! The delivered weight table into the project-registered store, then the hard labels into
        ! the project itself, so the assignment can be judged by an INDEPENDENT reconstructor. Every
        ! non-worker delivers (shared memory, nparts=1 and the distributed master alike); a worker
        ! shares the master's projfile and must not write it. The store goes first: it validates
        ! against the field's activity as the run saw it, before the labels overwrite `state`.
        if( .not. rounds%is_worker() )then
            call write_flex_weights_store(params, build, pinds, state_weights, labels, targets, bandwidths, &
                &l_merged)
            call write_discrete_state_project(build%spproj, pinds, labels, nstates, params%projfile)
        endif
        allocate(half_weights(nptcls,nstates), source=state_weights)
        ! nonuniform filtering LAST, so every delivered map is filtered the same way
        if( l_state_rec .and. trim(params%nufilt) == 'yes' ) call apply_consensus_nu_filter(params, nstates)
        call write_covariance_manifest(params, nptcls, ncomp, nstates, state_axis, min_neff, sigma_loaded)

        ! in resume mode no reconstructor was ever built -- only the cached embedding was read
        if( .not. l_resume )then
            call mean_rec%dealloc_rho
            call mean_rec%kill
            do q = 1, ncomp
                call basis_recs(q)%dealloc_rho
                call basis_recs(q)%kill
            end do
            deallocate(basis_recs)
        endif
        deallocate(pinds, labels, z, eigvals, prior_precision, latent_second, contrast)
        deallocate(resid_energy, resid_mean_energy)
        deallocate(state_weights, half_weights, targets, bandwidths, neff)
        if( allocated(pviews) ) deallocate(pviews)
        if( allocated(kdist) ) deallocate(kdist, kfloor)
        if( allocated(deconv_labels) ) deallocate(deconv_labels)
        call planes_kill()
        call plane_cache_close()
    end subroutine run_flex_pca

    !> The distributed worker's one entry: the shared preparation (its particle list, poses, the
    !! pinned sigma decision), then exactly one stage body. Round state (stage, which_iter,
    !! maxits, nfits, pcafit) arrives in params from the master's job_descr. The caller (the
    !! worker strategy) signals completion with qsys_job_finished.
    subroutine run_flex_pca_worker( params, build, cline , rounds)
        class(flex_pca_rounds), intent(inout) :: rounds
        type(parameters), intent(inout) :: params
        type(builder),    intent(inout) :: build
        class(cmdline),   intent(inout) :: cline
        type(reconstructor)   :: mean_rec
        integer,  allocatable :: pinds(:)
        real,     allocatable :: state_weights(:,:)
        integer :: nptcls, nstates
        logical :: sigma_loaded, l_split_eo
        call validate_covariance_inputs(params, build, cline, pinds, nptcls, rounds=rounds)
        call load_and_validate_sigma(params, build, cline, pinds, sigma_loaded, rounds=rounds)
        select case(params%stage)
            case(PCA_STAGE_STATES)
                ! the master produced the weight table; replicate its operator setup (ml_reg off,
                ! reconstruction Fourier band) and accumulate this part
                params%l_ml_reg = .false.
                params%ml_reg   = 'no'
                call cline%set('ml_reg','no')
                call build%esig%set_kfromto([1, max(1, fdim(flex_rec_box(params)) - 1)])
                call read_state_weights_round(pinds, nptcls, state_weights, nstates, l_split_eo)
                call reconstruct_flex_weighted_states(params, build, pinds, state_weights, &
                    &nstates, floor_rho=.true., split_eo=l_split_eo, rounds=rounds)
                deallocate(state_weights)
            case(PCA_STAGE_PROBE, PCA_STAGE_POLISH, PCA_STAGE_EMBED)
                ! the mean's scale comes from the master's cache; the read-once env caches are set
                ! exactly as on the master path
                call estimate_covariance_mean(params, build, mean_rec, pinds, nptcls, rounds=rounds)
                if( params%stage == PCA_STAGE_EMBED )then
                    call embed_worker_pass(params, build, mean_rec, pinds, nptcls, rounds=rounds)
                else
                    call probe_worker_pass(params, build, mean_rec, pinds, nptcls, rounds=rounds)
                endif
                call mean_rec%dealloc_rho; call mean_rec%kill
            case default
                THROW_HARD('unknown flex_pca worker stage')
        end select
        deallocate(pinds)
    end subroutine run_flex_pca_worker

    !> Read a part's particle-index list (one integer per line, as written by the master's
    !! flex_pca_plan_partitions through arr2txtfile).
    subroutine read_pind_list( fname, pinds, nptcls )
        character(len=*),     intent(in)  :: fname
        integer, allocatable, intent(out) :: pinds(:)
        integer,              intent(out) :: nptcls
        type(string) :: fn
        integer :: funit, io_stat, i, n
        fn = trim(fname)
        if( .not. file_exists(fn) ) THROW_HARD('flex_pca worker: particle list not found: '//trim(fname))
        n = nlines(fn)
        if( n < 1 ) THROW_HARD('flex_pca worker: empty particle list: '//trim(fname))
        allocate(pinds(n))
        call fopen(funit, file=fn, action='READ', status='OLD', iostat=io_stat)
        call fileiochk('read_pind_list', io_stat)
        do i = 1, n
            read(funit,*) pinds(i)
        end do
        call fclose(funit)
        call fn%kill
        nptcls = n
    end subroutine read_pind_list

    subroutine validate_covariance_inputs( params, build, cline, pinds, nptcls , rounds)
        class(flex_pca_rounds), intent(inout) :: rounds
        type(parameters), intent(inout) :: params
        type(builder),    intent(inout) :: build
        class(cmdline),   intent(in)    :: cline
        integer, allocatable, intent(out) :: pinds(:)
        integer, intent(out) :: nptcls
        integer :: q, i, cnt
        integer, allocatable :: sel(:)
        if( trim(params%oritype) /= 'ptcl3D' ) THROW_HARD('flex_pca requires oritype=ptcl3D')
        if( .not. cline%defined('vol1') )then
            THROW_HARD('flex_pca requires a consensus mean map: pass vol1 or register one in the project out segment')
        endif
        if( trim(params%ptcl_src) /= 'raw' ) THROW_HARD('flex_pca currently requires ptcl_src=raw')
        ! a run writes its state labels into its project copy, so a rerun must start from the original
        ! project (a subset is selected with pindfile=, never with the labels)
        if( build%spproj_field%get_n('state') /= 1 ) THROW_HARD('flex_pca works on one population but the project &
            &carries several state labels (a delivered copy): rerun from the original project, selecting particles &
            &with pindfile= if needed')
        if( cline%defined('pindfile') )then
            ! a distributed worker takes the master's partition as it was planned: its own
            ! particle-index list, never a re-derivation of the selection from fromp/top
            call read_pind_list(params%pindfile%to_char(), pinds, nptcls)
            write(logfhandle,'(A,I0,A,A)') '>>> FLEX_PCA (WORKER) ',nptcls,' particles from ', &
                &params%pindfile%to_char()
            call flush(logfhandle)
        else
            ! not sample4rec: its updatecnt > 0 condition belongs to trailing reconstruction
            allocate(sel(max(0, params%top - params%fromp + 1)))
            cnt = 0
            do i = params%fromp, params%top
                if( build%spproj_field%get_state(i) > 0 )then
                    cnt = cnt + 1; sel(cnt) = i
                endif
            end do
            nptcls = cnt
            if( allocated(pinds) ) deallocate(pinds)
            allocate(pinds(nptcls), source=sel(1:nptcls))
            deallocate(sel)
        endif
        if( nptcls < 100 ) THROW_HARD('flex_pca requires at least 100 active particles')
        if( build%spproj%os_ptcl2D%get_noris() > 0 .and. &
            &build%spproj%os_ptcl2D%get_noris() /= build%spproj%os_ptcl3D%get_noris() ) &
            &THROW_HARD('flex_pca requires matching ptcl2D and ptcl3D rows')
        ! eo split by index parity when the project's is degenerate. MASTER-ONLY and must be persisted: a
        ! worker counts over its own range, so it would silently build a different -- equally valid -- split.
        if( rounds%is_worker() )then
            if( count([(build%spproj_field%get_eo(pinds(q))==0,q=1,nptcls)]) < 1 .or. &
                &count([(build%spproj_field%get_eo(pinds(q))==1,q=1,nptcls)]) < 1 )then
                THROW_HARD('flex_pca worker sees a degenerate eo split; the master did not persist its repair')
            endif
        else if( count([(build%spproj_field%get_eo(pinds(q))==0,q=1,nptcls)]) < 20 .or. &
            &count([(build%spproj_field%get_eo(pinds(q))==1,q=1,nptcls)]) < 20 )then
            write(logfhandle,'(A)') '>>> FLEX_PCA assigning alternating even/odd halfsets (project eo was degenerate)'
            call build%spproj_field%partition_eo
            if( rounds%is_master() )then
                call build%spproj%write_segment_inside(params%oritype, params%projfile)
                write(logfhandle,'(A)') '>>> FLEX_PCA persisted the repaired eo split for the workers'
            endif
        endif
        if( count([(build%spproj_field%get_eo(pinds(q))==0,q=1,nptcls)]) < 20 .or. &
            &count([(build%spproj_field%get_eo(pinds(q))==1,q=1,nptcls)]) < 20 ) &
            &THROW_HARD('flex_pca requires populated even and odd halfsets')
    end subroutine validate_covariance_inputs

    !> Cost of the requested state count in resident reconstructors. REPORT ONLY: a former hard 64 GB
    !! refusal blocked runs that fit and missed ones that did not.
    subroutine report_state_memory( params, nstates )
        class(parameters), intent(in) :: params
        integer,           intent(in) :: nstates
        real(dp) :: gb, gb_proc, nexp
        integer  :: box_rec, dim_exp, nproc
        box_rec = flex_rec_box(params)
        ! expanded grid half-width, as reconstructor::alloc_rho derives it: |lims| + ceiling(KBWINSZ)
        dim_exp = box_rec/2 + ceiling(KBWINSZ) + 1
        nexp    = (2.d0*real(dim_exp,dp) + 1.d0)**3
        ! complex cmat_exp (8 B) + real rho_exp (4 B) per grid point, x2 when even and odd coexist
        gb_proc = 2.d0 * real(nstates,dp) * nexp * 12.d0 / 1.d9
        ! every process pays this, so the machine-wide peak is nparts+1 times the per-process figure
        nproc   = max(1, params%nparts) + 1
        gb      = gb_proc * real(nproc,dp)
        write(logfhandle,'(A,I0,A,I0,A,F8.2,A,I0,A,F8.2,A)') '>>> FLEX_PCA states=',nstates, &
            &' at reconstruction box ',box_rec,' -> approx ',gb_proc,' GB of reconstructors per &
            &process x ',nproc,' processes = ',gb,' GB machine-wide'
        write(logfhandle,'(A)') '>>> FLEX_PCA the knobs that move it are npreimages (linear) and &
            &box_crop (cubic)'
        ! the reconstructors rarely hurt: the reduced solve's accumulator is sized against COV_ATHR_BUDGET,
        ! likewise per process, so a distributed run multiplies it by nproc. SIMPLE_COV_DTILDE moves it.
        if( params%nparts > 1 ) write(logfhandle,'(A,I0,A)') '>>> FLEX_PCA NOTE: the reduced-solve &
            &accumulator is also per process; at nparts=',params%nparts,' it is paid that many times &
            &over. Cap it with SIMPLE_COV_DTILDE if the machine is tight.'
        call flush(logfhandle)
    end subroutine report_state_memory

    subroutine load_and_validate_sigma( params, build, cline, pinds, loaded , rounds)
        class(flex_pca_rounds), intent(inout) :: rounds
        type(parameters), intent(inout) :: params
        type(builder),    intent(inout) :: build
        class(cmdline),   intent(inout) :: cline
        integer,          intent(in)    :: pinds(:)
        logical,          intent(out)   :: loaded
        integer :: i, k, iptcl, noris, fromp_save, top_save
        ! The sigma table is allocated on params%fromp:top. A distributed worker's particle list is the
        ! master's partition (pindfile), not the fromp/top split of the project, so widen the range to
        ! every particle for the load (the canonical reader fills rows for state>0 within the range) and
        ! restore after. The unit fallback below has always done the same.
        noris      = build%spproj_field%get_noris()
        fromp_save = params%fromp
        top_save   = params%top
        params%fromp = 1
        params%top   = noris
        call load_sigma2_groups(params, build%pftc, build%esig, build%spproj, build%spproj_field, loaded)
        params%fromp = fromp_save
        params%top   = top_save
        if( .not. loaded )then
            ! gen_fplane4rec follows norm_noise_taper_edge_pad_fft, so a unit spectrum is the correct fallback
            ! Construct the sigma object through its own API (the legacy loader's idiom: widen
            ! fromp/top to every particle so the table covers the whole project, restore after):
            ! esig%new also registers the table with the polar calculator; a bare allocate of the
            ! component inside a never-constructed object leaves the consumers with an
            ! inconsistent object (gate crash in the mean-stage prep, 2026-09-08)
            params%fromp = 1
            params%top   = noris
            call build%esig%new(params, build%pftc, string('flex_pca_unit_sigma2.dat'), params%box)
            params%fromp = fromp_save
            params%top   = top_save
            build%esig%sigma2_noise = 1.0
            loaded = .true.
            write(logfhandle,'(A)') '>>> FLEX_PCA SIGMA: using unit white-noise spectrum after background normalization'
            write(logfhandle,'(A)') '>>> FLEX_PCA SIGMA: provide sigma2_group*.star for coloured experimental noise'
        endif
        do i = 1, size(pinds)
            iptcl = pinds(i)
            do k = lbound(build%esig%sigma2_noise,1), ubound(build%esig%sigma2_noise,1)
                if( .not. ieee_is_finite(build%esig%sigma2_noise(k,iptcl)) .or. &
                    &build%esig%sigma2_noise(k,iptcl) <= TINY )then
                    THROW_HARD('flex_pca found a nonpositive or nonfinite sigma2 value')
                endif
            end do
        end do
        params%ml_reg   = 'yes'
        params%l_ml_reg = .true.
        call cline%set('ml_reg','yes')
        ! Pin the sigma2 decision across the master/worker boundary (see write_sigma_state).
        if( rounds%is_master() )then
            call write_sigma_state(loaded)
        else if( rounds%is_worker() )then
            call check_sigma_state(loaded)
        endif
    end subroutine load_and_validate_sigma


    !> Calibrate the per-particle noise (from the even/odd half solutions when the run has them,
    !! else from the scale file the original run wrote) and replace z / precision by the posterior
    !! means / precisions under the deconvolved mixture prior. applied=.false. when SIMPLE_COV_DECONV=0.
    subroutine apply_latent_deconvolution( z, precision, eigvals, zhalf, pinds, nptcls, ncomp, applied, labels, &
        &contrast, resid_energy, resid_mean_energy, sig2_eff, resume, adopted, srcdir, srcfile )
        integer,  intent(in)    :: pinds(:), nptcls, ncomp
        integer, allocatable, optional, intent(inout) :: labels(:)   !< mixture component per particle
        !> with these present the deconvolved coordinates are cached to flex_pca_embedding_deconv.bin
        !! (same layout as the raw cache) so a resume adopts them instead of deconvolving again
        real(dp), optional, intent(in)  :: contrast(:), resid_energy(:), resid_mean_energy(:), sig2_eff
        logical,  optional, intent(in)  :: resume
        logical,  optional, intent(out) :: adopted
        !> directory of the embedding a resume was given (infile): the deconvolved cache and the labels are
        !! looked up there when the run directory has none, so a states-only resume in a fresh directory
        !! adopts the original run's deconvolution instead of re-running the K ladder
        character(len=*), optional, intent(in) :: srcdir
        character(len=*), optional, intent(in) :: srcfile   !< the infile itself: its trailing deconvolved block is adopted first
        real(dp), intent(inout) :: z(:,:), precision(:,:,:)
        real(dp), intent(in)    :: eigvals(:)
        real(dp), allocatable, intent(inout) :: zhalf(:,:,:)
        logical,  intent(out)   :: applied
        real(dp) :: prior(ncomp), a_comp(ncomp), noise_scale
        integer  :: q, k_deconv, u_ns, io_ns, ncomp_c, i, pind_c, comp_c
        real(dp) :: resp_c
        real(dp), allocatable :: z_c(:,:), eig_c(:), con_c(:), re_c(:), rme_c(:), prec_c(:,:,:)
        real(dp) :: sig2_c
        logical  :: l_resume
        character(len=:), allocatable :: cache_fname, labels_fname, ns_fname
        applied = .false.
        if( present(adopted) ) adopted = .false.
        l_resume = .false.
        if( present(resume) ) l_resume = resume
        ! ---- resume: adopt the deconvolved block of the infile cache (one file carries raw + deconvolved) ----
        if( l_resume .and. present(srcfile) )then
            if( len_trim(srcfile) > 0 )then
                if( file_exists(srcfile) )then
                    block
                        real(dp), allocatable :: zb(:,:), pb(:,:,:)
                        integer,  allocatable :: lb(:)
                        real(dp) :: nsb
                        logical  :: l_found
                        call read_deconv_block(srcfile, nptcls, ncomp, zb, pb, lb, nsb, l_found)
                        if( l_found )then
                            z(1:nptcls,1:ncomp) = zb
                            precision(1:ncomp,1:ncomp,1:nptcls) = pb
                            if( present(labels) )then
                                if( allocated(labels) ) deallocate(labels)
                                allocate(labels(nptcls), source=lb)
                                if( any(labels < 1) ) deallocate(labels)
                            endif
                            write(logfhandle,'(A)') '>>> FLEX_PCA resumed embedding: deconvolved block adopted from '//&
                                &trim(srcfile)//' (no re-deconvolution)'
                            call flush(logfhandle)
                            applied = .true.
                            if( present(adopted) ) adopted = .true.
                            deallocate(zb, pb, lb)
                            return
                        endif
                    end block
                endif
            endif
        endif
        ! ---- legacy: the separate deconvolved cache older runs wrote (here, or next to infile) ----
        cache_fname  = 'flex_pca_embedding_deconv.bin'
        labels_fname = 'flex_pca_deconv_labels.txt'
        if( l_resume .and. .not. file_exists(cache_fname) .and. present(srcdir) )then
            if( len_trim(srcdir) > 0 )then
                if( file_exists(trim(srcdir)//'/flex_pca_embedding_deconv.bin') )then
                    cache_fname  = trim(srcdir)//'/flex_pca_embedding_deconv.bin'
                    labels_fname = trim(srcdir)//'/flex_pca_deconv_labels.txt'
                    write(logfhandle,'(A)') '>>> FLEX_PCA resumed embedding: deconvolved cache found next to infile: '//cache_fname
                endif
            endif
        endif
        if( l_resume .and. file_exists(cache_fname) )then
            call read_embedding_cache(cache_fname, pinds, nptcls, ncomp_c, z_c, eig_c, &
                &con_c, re_c, rme_c, prec_c, sig2_c)
            if( ncomp_c == ncomp .and. size(z_c,1) == nptcls )then
                z(1:nptcls,1:ncomp) = z_c
                precision(1:ncomp,1:ncomp,1:nptcls) = prec_c
                if( present(labels) )then
                    if( allocated(labels) ) deallocate(labels)
                    allocate(labels(nptcls), source=0)
                    open(newunit=u_ns, file=labels_fname, status='old', action='read', iostat=io_ns)
                    if( io_ns == 0 )then
                        read(u_ns,*,iostat=io_ns)      ! header
                        do i = 1, nptcls
                            read(u_ns,*,iostat=io_ns) pind_c, comp_c, resp_c
                            if( io_ns /= 0 .or. pind_c /= pinds(i) ) exit
                            labels(i) = comp_c
                        end do
                        close(u_ns)
                    endif
                    if( io_ns /= 0 .or. any(labels < 1) )then
                        write(logfhandle,'(A)') '>>> FLEX_PCA resumed deconvolution: labels file missing or &
                            &inconsistent; the state stage will rediscover its clusters'
                        deallocate(labels)
                    endif
                endif
                write(logfhandle,'(A)') '>>> FLEX_PCA resumed embedding: deconvolved coordinates adopted from '//&
                    &cache_fname//' (no re-deconvolution)'
                call flush(logfhandle)
                applied = .true.
                if( present(adopted) ) adopted = .true.
                deallocate(z_c, eig_c, con_c, re_c, rme_c, prec_c)
                return
            endif
            deallocate(z_c, eig_c, con_c, re_c, rme_c, prec_c)
        endif
        do q = 1, ncomp
            prior(q) = 1.d0 / max(eigvals(q), DTINY)
        end do
        if( allocated(zhalf) )then
            call calibrate_noise_scale(z, zhalf, precision, prior, nptcls, ncomp, noise_scale, a_comp)
            open(newunit=u_ns, file='flex_pca_noise_scale.txt', status='replace', action='write')
            write(u_ns,'(ES16.8)') noise_scale
            close(u_ns)
            deallocate(zhalf)
        else
            noise_scale = 1.d0
            ! the calibrated noise scale of the fit: here, or next to infile (a resume in a fresh directory
            ! used to silently fall back to 1.0 and deconvolve differently from the fit)
            ns_fname = 'flex_pca_noise_scale.txt'
            if( .not. file_exists(ns_fname) .and. present(srcdir) )then
                if( len_trim(srcdir) > 0 )then
                    if( file_exists(trim(srcdir)//'/flex_pca_noise_scale.txt') ) ns_fname = trim(srcdir)//'/flex_pca_noise_scale.txt'
                endif
            endif
            open(newunit=u_ns, file=ns_fname, status='old', action='read', iostat=io_ns)
            if( io_ns == 0 )then
                read(u_ns,*,iostat=io_ns) noise_scale
                close(u_ns)
            endif
            if( io_ns /= 0 ) noise_scale = 1.d0
            write(logfhandle,'(A,A,A,F8.3)') '>>> FLEX_PCA resumed embedding: noise scale from ', ns_fname, &
                &' (1.0 when absent) =', noise_scale
        endif
        call deconvolve_latent(z, precision, prior, nptcls, ncomp, noise_scale, 16, k_deconv, &
            &prior_fname='flex_pca_deconv_prior.txt', labels_fname='flex_pca_deconv_labels.txt', pinds=pinds, &
            &labels_out=labels)
        applied = .true.
        ! the deconvolved coordinates join the raw cache of this run as a trailing block (one file)
        if( file_exists('flex_pca_embedding.bin') )then
            if( present(labels) )then
                if( allocated(labels) )then
                    call append_deconv_block('flex_pca_embedding.bin', nptcls, ncomp, z, precision, labels, noise_scale)
                else
                    call append_deconv_block('flex_pca_embedding.bin', nptcls, ncomp, z, precision, noise_scale=noise_scale)
                endif
            else
                call append_deconv_block('flex_pca_embedding.bin', nptcls, ncomp, z, precision, noise_scale=noise_scale)
            endif
        else
            write(logfhandle,'(A)') '>>> FLEX_PCA deconvolved coordinates not cached: no flex_pca_embedding.bin in this directory'
        endif
    end subroutine apply_latent_deconvolution

    !> Byte offset of the first byte after the raw block of a cache written by write_embedding_cache
    !! (stream access: magic(8) ver(4) nptcls(4) ncomp(4) box(4) smpd(4) pinds(4n) z(8nk) eigvals(8k)
    !! contrast(8n) resid(8n) resid_mean(8n) precision(8kkn) sig2(8)).
    pure function raw_block_end( nptcls, ncomp ) result( pos )
        integer, intent(in) :: nptcls, ncomp
        integer(kind=8) :: pos, n, k
        n = int(nptcls, 8); k = int(ncomp, 8)
        pos = 8_8 + 4_8 + 4_8 + 4_8 + 4_8 + 4_8 + 4_8*n + 8_8*n*k + 8_8*k + 3_8*8_8*n + 8_8*k*k*n + 8_8 + 1_8
    end function raw_block_end

    !> Append the deconvolved block to an existing raw cache. Anything already after the raw block
    !! (an older deconvolved block) is overwritten: the file is truncated at the raw block's end.
    subroutine append_deconv_block( fname, nptcls, ncomp, z, precision, labels, noise_scale )
        character(len=*),  intent(in) :: fname
        integer,           intent(in) :: nptcls, ncomp
        real(dp),          intent(in) :: z(:,:), precision(:,:,:)
        integer, optional, intent(in) :: labels(:)
        real(dp), optional, intent(in) :: noise_scale
        integer, allocatable :: lab(:)
        real(dp) :: ns
        integer  :: u, io
        integer(kind=8) :: pos
        allocate(lab(nptcls), source=0)
        if( present(labels) ) lab(1:nptcls) = labels(1:nptcls)
        ns = 1.d0; if( present(noise_scale) ) ns = noise_scale
        pos = raw_block_end(nptcls, ncomp)
        open(newunit=u, file=fname, status='old', action='readwrite', access='stream', form='unformatted', iostat=io)
        if( io /= 0 ) THROW_HARD('append_deconv_block: cannot open '//trim(fname))
        write(u, pos=pos) COV_DECONV_MAGIC
        write(u) nptcls, ncomp
        write(u) z(1:nptcls,1:ncomp)
        write(u) precision(1:ncomp,1:ncomp,1:nptcls)
        write(u) lab
        write(u) ns
        endfile(u)
        close(u)
        write(logfhandle,'(A,A)') '>>> FLEX_PCA deconvolved coordinates, precisions and labels appended to ', trim(fname)
        call flush(logfhandle)
        deallocate(lab)
    end subroutine append_deconv_block

    !> Read the deconvolved block of a cache if present and consistent with (nptcls, ncomp); found=.false.
    !! otherwise (older caches end after the raw block).
    subroutine read_deconv_block( fname, nptcls, ncomp, z, precision, labels, noise_scale, found )
        character(len=*),      intent(in)  :: fname
        integer,               intent(in)  :: nptcls, ncomp
        real(dp), allocatable, intent(out) :: z(:,:), precision(:,:,:)
        integer,  allocatable, intent(out) :: labels(:)
        real(dp),              intent(out) :: noise_scale
        logical,               intent(out) :: found
        character(len=len(COV_DECONV_MAGIC)) :: magic
        character(len=len(COV_CACHE_MAGIC))  :: magic0
        integer :: u, io, ver, n_c, k_c, n_b, k_b
        integer(kind=8) :: pos
        found = .false.; noise_scale = 1.d0
        if( .not. file_exists(fname) ) return
        open(newunit=u, file=fname, status='old', action='read', access='stream', form='unformatted', iostat=io)
        if( io /= 0 ) return
        read(u, iostat=io) magic0, ver
        if( io /= 0 .or. magic0 /= COV_CACHE_MAGIC .or. ver /= COV_CACHE_VERSION )then; close(u); return; endif
        read(u, iostat=io) n_c, k_c
        if( io /= 0 .or. n_c /= nptcls .or. k_c /= ncomp )then; close(u); return; endif
        pos = raw_block_end(nptcls, ncomp)
        read(u, pos=pos, iostat=io) magic
        if( io /= 0 .or. magic /= COV_DECONV_MAGIC )then; close(u); return; endif
        read(u, iostat=io) n_b, k_b
        if( io /= 0 .or. n_b /= nptcls .or. k_b /= ncomp )then; close(u); return; endif
        allocate(z(nptcls,ncomp), precision(ncomp,ncomp,nptcls), labels(nptcls))
        read(u, iostat=io) z
        if( io == 0 ) read(u, iostat=io) precision
        if( io == 0 ) read(u, iostat=io) labels
        if( io == 0 ) read(u, iostat=io) noise_scale
        close(u)
        if( io /= 0 )then
            deallocate(z, precision, labels); return
        endif
        found = .true.
    end subroutine read_deconv_block



    !> Binary cache of the state stages' inputs, so a different state count / bandwidth / placement can
    !! be tried without re-fitting the basis and re-embedding every particle.
    subroutine write_embedding_cache( fname, pinds, nptcls, ncomp, z, eigvals, contrast, &
        &resid_energy, resid_mean_energy, precision, sig2_eff )
        character(len=*), intent(in) :: fname
        integer,          intent(in) :: pinds(:), nptcls, ncomp
        real(dp),         intent(in) :: z(nptcls,ncomp), eigvals(ncomp), contrast(nptcls)
        real(dp),         intent(in) :: resid_energy(nptcls), resid_mean_energy(nptcls)
        real(dp),         intent(in) :: precision(ncomp,ncomp,nptcls), sig2_eff
        integer :: u
        call del_file(fname)
        open(newunit=u, file=fname, status='replace', action='write', access='stream', form='unformatted')
        write(u) COV_CACHE_MAGIC, COV_CACHE_VERSION
        write(u) nptcls, ncomp
        write(u) cov_box_crop_glob, cov_smpd_crop_glob
        write(u) pinds(1:nptcls)
        write(u) z
        write(u) eigvals
        write(u) contrast
        write(u) resid_energy
        write(u) resid_mean_energy
        write(u) precision
        write(u) sig2_eff
        close(u)
        write(logfhandle,'(A,A,A,F8.1,A)') '>>> FLEX_PCA embedding cached to ',trim(fname), &
            &' (',real(8*(nptcls*ncomp + ncomp*ncomp*nptcls))/1048576.0,' MB); reuse with infile=<path>'
        call flush(logfhandle)
    end subroutine write_embedding_cache

    subroutine read_embedding_cache( fname, pinds, nptcls, ncomp, z, eigvals, contrast, &
        &resid_energy, resid_mean_energy, precision, sig2_eff )
        character(len=*), intent(in)    :: fname
        integer,          intent(in)    :: pinds(:)
        integer,          intent(in)    :: nptcls
        integer,          intent(out)   :: ncomp
        real(dp), allocatable, intent(out) :: z(:,:), eigvals(:), contrast(:)
        real(dp), allocatable, intent(out) :: resid_energy(:), resid_mean_energy(:)
        real(dp), allocatable, intent(out) :: precision(:,:,:)
        real(dp),         intent(out)   :: sig2_eff
        integer, allocatable :: pinds_cached(:)
        character(len=len(COV_CACHE_MAGIC)) :: magic
        integer :: u, ver, nptcls_c, i, box_c
        real    :: smpd_c
        if( .not. file_exists(fname) ) THROW_HARD('flex_pca embedding cache not found: '//trim(fname))
        open(newunit=u, file=fname, status='old', action='read', access='stream', form='unformatted')
        read(u) magic, ver
        if( magic /= COV_CACHE_MAGIC ) THROW_HARD('not a flex_pca embedding cache: '//trim(fname))
        if( ver /= COV_CACHE_VERSION ) THROW_HARD('flex_pca embedding cache version mismatch; re-run the fit')
        read(u) nptcls_c, ncomp
        if( nptcls_c /= nptcls ) THROW_HARD('flex_pca embedding cache particle count does not match the project')
        read(u) box_c, smpd_c
        if( cov_box_crop_glob > 0 .and. (box_c /= cov_box_crop_glob .or. &
            &abs(smpd_c - cov_smpd_crop_glob) > 1.e-3*cov_smpd_crop_glob) )then
            write(logfhandle,'(A,I0,A,F7.3,A,I0,A,F7.3,A)') '>>> FLEX_PCA embedding cache lattice box_crop=', &
                &box_c, ' @ ', smpd_c, ' A; this run: box_crop=', cov_box_crop_glob, ' @ ', cov_smpd_crop_glob, ' A'
            THROW_HARD('flex_pca embedding cache was built at a different working sampling (box_crop/lp); &
                &resume at the cached lattice or re-run the fit')
        endif
        allocate(pinds_cached(nptcls_c))
        read(u) pinds_cached
        do i = 1, nptcls
            if( pinds_cached(i) /= pinds(i) ) &
                &THROW_HARD('flex_pca embedding cache was built from a different particle selection')
        end do
        allocate(z(nptcls,ncomp), eigvals(ncomp), contrast(nptcls), resid_energy(nptcls), &
            &resid_mean_energy(nptcls), precision(ncomp,ncomp,nptcls))
        read(u) z
        read(u) eigvals
        read(u) contrast
        read(u) resid_energy
        read(u) resid_mean_energy
        read(u) precision
        read(u) sig2_eff
        close(u)
        deallocate(pinds_cached)
    end subroutine read_embedding_cache

    !> Drop every state whose effective sample size is below min_neff and compact the per-state
    !! arrays, so the reconstruction, the bandwidth CV and the merge all see only states that can
    !! support a map. Particles whose argmax state is dropped become unassigned (label 0) and feed no
    !! map -- they are, by construction, the particles the placement could not commit anywhere. At
    !! least two states always survive: with fewer the run has no heterogeneity to deliver.
    subroutine prune_underpopulated_states( nptcls, nstates, min_neff, weights, targets, bandwidths, &
        &neff, labels, dist, bfloor )
        integer,              intent(in)    :: nptcls, min_neff
        integer,              intent(inout) :: nstates
        real,    allocatable, intent(inout) :: weights(:,:), targets(:,:), bandwidths(:), neff(:)
        integer, allocatable, intent(inout) :: labels(:)
        real(dp), allocatable, intent(inout) :: dist(:,:), bfloor(:)
        logical,  allocatable :: keep(:)
        integer,  allocatable :: map(:), occ(:), ord(:)
        real,     allocatable :: w2(:,:), t2(:,:), b2(:), n2(:)
        real(dp), allocatable :: d2(:,:), f2(:)
        real,     allocatable :: key(:)
        integer :: s, i, nkeep, ndrop, nlost, snew
        if( nstates < 2 ) return
        allocate(keep(nstates), source=.true.)
        allocate(map(nstates), occ(nstates), source=0)
        do i = 1, nptcls
            if( labels(i) >= 1 .and. labels(i) <= nstates ) occ(labels(i)) = occ(labels(i)) + 1
        end do
        do s = 1, nstates
            keep(s) = neff(s) >= real(min_neff)
        end do
        nkeep = count(keep)
        if( nkeep >= nstates ) return
        if( nkeep < 2 )then
            ! nothing clears the floor: keep the two best-supported seats rather than delivering none
            allocate(key(nstates), ord(nstates))
            key = neff
            do s = 1, nstates
                ord(s) = s
            end do
            call hpsort(key, ord)
            keep = .false.
            keep(ord(nstates))   = .true.
            keep(ord(nstates-1)) = .true.
            nkeep = 2
            deallocate(key, ord)
        endif
        ndrop = nstates - nkeep
        write(logfhandle,'(A,I0,A,I0,A,I0,A)') '>>> FLEX_PCA OCCUPANCY FLOOR: ', ndrop, ' of ', nstates, &
            &' states have fewer than ', min_neff, ' effective particles and are dropped before reconstruction'
        do s = 1, nstates
            if( keep(s) ) cycle
            write(logfhandle,'(A,I3,A,F10.1,A,I0)') '>>>   dropped state=', s, '  neff=', neff(s), &
                &'  particles=', occ(s)
        end do
        snew = 0
        do s = 1, nstates
            if( .not. keep(s) ) cycle
            snew   = snew + 1
            map(s) = snew
        end do
        allocate(w2(nptcls,nkeep), t2(size(targets,1),nkeep), b2(nkeep), n2(nkeep))
        snew = 0
        do s = 1, nstates
            if( .not. keep(s) ) cycle
            snew = snew + 1
            w2(:,snew) = weights(:,s)
            t2(:,snew) = targets(:,s)
            b2(snew)   = bandwidths(s)
            n2(snew)   = neff(s)
        end do
        call move_alloc(w2, weights)
        call move_alloc(t2, targets)
        call move_alloc(b2, bandwidths)
        call move_alloc(n2, neff)
        if( allocated(dist) )then
            allocate(d2(nptcls,nkeep))
            snew = 0
            do s = 1, nstates
                if( .not. keep(s) ) cycle
                snew = snew + 1
                d2(:,snew) = dist(:,s)
            end do
            call move_alloc(d2, dist)
        endif
        if( allocated(bfloor) )then
            allocate(f2(nkeep))
            snew = 0
            do s = 1, nstates
                if( .not. keep(s) ) cycle
                snew = snew + 1
                f2(snew) = bfloor(s)
            end do
            call move_alloc(f2, bfloor)
        endif
        nlost = 0
        do i = 1, nptcls
            if( labels(i) < 1 .or. labels(i) > nstates ) cycle
            if( keep(labels(i)) )then
                labels(i) = map(labels(i))
            else
                labels(i) = 0
                nlost     = nlost + 1
            endif
        end do
        write(logfhandle,'(A,I0,A,F6.2,A,I0,A)') '>>> FLEX_PCA OCCUPANCY FLOOR: ', nlost, ' particles (', &
            &100.0*real(nlost)/real(max(nptcls,1)), '%) lost their state and feed no map; ', nkeep, ' states remain'
        call flush(logfhandle)
        nstates = nkeep
        deallocate(keep, map, occ)
    end subroutine prune_underpopulated_states

    subroutine write_covariance_eigenvolumes( basis_recs, eigvals, ncomp, fname )
        type(reconstructor), intent(inout) :: basis_recs(ncomp)
        real(dp),            intent(in)    :: eigvals(ncomp)
        integer,             intent(in)    :: ncomp
        !> per-fit eigenvalue-table namespace (paired engine); default flex_pca_eigenvalues.txt
        character(len=*), optional, intent(in) :: fname
        character(len=:), allocatable :: fn
        integer :: q, u
        fn = 'flex_pca_eigenvalues.txt'
        if( present(fname) ) fn = trim(fname)
        ! the eigenvolume MRCs are written by form_eigenbasis_from_reduced; only the table is written here
        call del_file(fn)
        open(newunit=u,file=fn,status='replace',action='write')
        write(u,'(A)') '# component eigenvalue'
        do q = 1, ncomp
            write(u,'(I6,1X,ES20.10)') q,eigvals(q)
        end do
        close(u)
    end subroutine write_covariance_eigenvolumes


    !> Nonuniform filtering of the delivered state maps, using ONE filter from the CONSENSUS half maps.
    subroutine apply_consensus_nu_filter( params, nstates )
        use simple_nu_filter, only: setup_nu_dmats, optimize_nu_cutoff_finds, nu_filter_vol, &
            &cleanup_nu_filter, write_nu_local_resolution_map, get_nu_filtmap_finest_selected_lp
        class(parameters), intent(in) :: params
        integer,           intent(in) :: nstates
        type(image)  :: vol_e, vol_o, vin, vout
        type(string) :: fn, fn_out, spe, spo
        character(len=:), allocatable :: pe, po, base
        integer :: s, ldim(3)
        ! consensus half maps: explicit vol_even/vol_odd if given, else derived from the vol1 stem
        base = params%vols(1)%to_char()
        if( len_trim(params%vol_even%to_char()) > 0 .and. len_trim(params%vol_odd%to_char()) > 0 )then
            pe = params%vol_even%to_char(); po = params%vol_odd%to_char()
        else
            if( index(base, MRC_EXT) < 1 )then
                write(logfhandle,'(A)') '>>> FLEX_PCA nufilt: cannot derive half-map names from vol1; &
                    &pass vol_even and vol_odd explicitly'
                return
            endif
            pe = base(:index(base, MRC_EXT, back=.true.)-1)//'_even'//MRC_EXT
            po = base(:index(base, MRC_EXT, back=.true.)-1)//'_odd'//MRC_EXT
        endif
        if( .not. file_exists(pe) .or. .not. file_exists(po) )then
            write(logfhandle,'(A)') '>>> FLEX_PCA nufilt: consensus half maps not found ('//pe//', '//po// &
                &'); skipping nonuniform filtering'
            return
        endif
        write(logfhandle,'(A)') '>>> FLEX_PCA nonuniform filter from the consensus half maps:'
        write(logfhandle,'(A)') '>>>   '//pe
        write(logfhandle,'(A)') '>>>   '//po
        call vol_e%new([params%box_crop,params%box_crop,params%box_crop], params%smpd_crop)
        call vol_o%new([params%box_crop,params%box_crop,params%box_crop], params%smpd_crop)
        spe = pe; spo = po
        call vol_e%read_and_crop(spe, params%smpd, params%box_crop, params%smpd_crop)
        call vol_o%read_and_crop(spo, params%smpd, params%box_crop, params%smpd_crop)
        call spe%kill; call spo%kill
        ldim = [params%box_crop,params%box_crop,params%box_crop]
        call setup_nu_dmats(vol_e, vol_o, params%mskdiam, [real ::])
        call optimize_nu_cutoff_finds()
        write(logfhandle,'(A,F8.2,A)') '>>> FLEX_PCA nufilt: finest selected local resolution ', &
            &get_nu_filtmap_finest_selected_lp(),' A'
        fn = 'flex_pca_nu_locres'//MRC_EXT
        call write_nu_local_resolution_map(fn)
        call fn%kill
        ! apply the SAME filter map to every delivered volume
        do s = 1, nstates
            fn     = 'flex_pca_state_'//int2str_pad(s,3)//MRC_EXT
            fn_out = 'flex_pca_state_'//int2str_pad(s,3)//'_nu'//MRC_EXT
            if( .not. file_exists(fn%to_char()) )then
                call fn%kill; call fn_out%kill; cycle
            endif
            call vin%new(ldim, params%smpd_crop)
            call vin%read(fn)
            call nu_filter_vol(vin, vout)
            call vout%write(fn_out, del_if_exists=.true.)
            call vin%kill; call vout%kill
            call fn%kill; call fn_out%kill
        end do
        call cleanup_nu_filter()
        call vol_e%kill; call vol_o%kill
        write(logfhandle,'(A)') '>>> FLEX_PCA nonuniform-filtered maps written as *_nu.mrc &
            &(originals retained); local resolution map: flex_pca_nu_locres.mrc'
        call flush(logfhandle)
    end subroutine apply_consensus_nu_filter


    !> Coordinates table for the paired-merge delivery (par.7 step 3): the standard
    !! flex_pca_coordinates.txt column layout -- consumers slice the latents as "everything
    !! from column 6 on" -- with label 0 everywhere, because the state stage never ran.
    subroutine write_merged_coordinates( build, pinds, z, resid_energy, resid_mean_energy )
        type(builder), intent(inout) :: build
        integer,       intent(in)    :: pinds(:)
        real(dp),      intent(in)    :: z(:,:)
        real(dp),      intent(in)    :: resid_energy(:), resid_mean_energy(:)
        integer :: u, i, q
        call del_file('flex_pca_coordinates.txt')
        open(newunit=u,file='flex_pca_coordinates.txt',status='replace',action='write')
        write(u,'(A)',advance='no') '# particle eo label residual mean_residual'
        do q=1,size(z,2); write(u,'(A,I0)',advance='no') ' z',q; end do
        write(u,*)
        do i=1,size(pinds)
            write(u,'(I10,1X,I1,1X,I4,2(1X,ES16.8))',advance='no') pinds(i), &
                &build%spproj_field%get_eo(pinds(i)), 0, resid_energy(i), resid_mean_energy(i)
            do q=1,size(z,2); write(u,'(1X,ES16.8)',advance='no') z(i,q); end do
            write(u,*)
        end do
        close(u)
        write(logfhandle,'(A,I0,A)') '>>> FLEX_PCA MERGED coordinates written for ',size(pinds), &
            &' particles (label=0: no state stage in the merge delivery)'
        call flush(logfhandle)
    end subroutine write_merged_coordinates

    subroutine write_covariance_tables( build, pinds, z, eigvals, prior_precision, weights, labels, &
        &targets, bandwidths, neff, resid_energy, resid_mean_energy, contrast )
        type(builder), intent(inout) :: build
        integer,       intent(in) :: pinds(:), labels(:)
        real(dp),      intent(in) :: z(:,:), eigvals(:), prior_precision(:)
        real,          intent(in) :: weights(:,:), targets(:,:), bandwidths(:), neff(:)
        real(dp),      intent(in) :: resid_energy(:), resid_mean_energy(:)
        !> fitted per-particle contrast. Written to its own file rather than appended as a column
        !! of flex_pca_coordinates.txt, because every consumer of that file slices the latents as
        !! "everything from column 6 on" and an extra column would silently become a latent.
        real(dp), optional, intent(in) :: contrast(:)
        integer :: u, i, q, state
        call del_file('flex_pca_coordinates.txt')
        open(newunit=u,file='flex_pca_coordinates.txt',status='replace',action='write')
        write(u,'(A)',advance='no') '# particle eo label residual mean_residual'
        do q=1,size(z,2); write(u,'(A,I0)',advance='no') ' z',q; end do
        write(u,*)
        do i=1,size(pinds)
            write(u,'(I10,1X,I1,1X,I4,2(1X,ES16.8))',advance='no') pinds(i), &
                &build%spproj_field%get_eo(pinds(i)),labels(i),resid_energy(i),resid_mean_energy(i)
            do q=1,size(z,2); write(u,'(1X,ES16.8)',advance='no') z(i,q); end do
            write(u,*)
        end do
        close(u)
        if( present(contrast) )then
            call del_file('flex_pca_contrast.txt')
            open(newunit=u,file='flex_pca_contrast.txt',status='replace',action='write')
            write(u,'(A)') '# particle contrast'
            do i=1,min(size(pinds),size(contrast))
                write(u,'(I10,1X,ES16.8)') pinds(i), contrast(i)
            end do
            close(u)
        endif
        call del_file('flex_pca_state_weights.txt')
        open(newunit=u,file='flex_pca_state_weights.txt',status='replace',action='write')
        write(u,'(A)',advance='no') '# particle'
        do state=1,size(weights,2); write(u,'(A,I3.3)',advance='no') ' w',state; end do
        write(u,*)
        do i=1,size(pinds)
            write(u,'(I10)',advance='no') pinds(i)
            do state=1,size(weights,2); write(u,'(1X,ES16.8)',advance='no') weights(i,state); end do
            write(u,*)
        end do
        close(u)
        call write_states_umap_figure(pinds, z, weights)
        call del_file('flex_pca_state_targets.txt')
        open(newunit=u,file='flex_pca_state_targets.txt',status='replace',action='write')
        ! targets are full latent-space points, so every coordinate is written
        write(u,'(A)',advance='no') '# state bandwidth effective_particles'
        do q=1,size(targets,1); write(u,'(A,I0)',advance='no') ' t',q; end do
        write(u,*)
        do state=1,size(targets,2)
            write(u,'(I5,2(1X,ES16.8))',advance='no') state,bandwidths(state),neff(state)
            do q=1,size(targets,1); write(u,'(1X,ES16.8)',advance='no') targets(q,state); end do
            write(u,*)
        end do
        close(u)
        call del_file('flex_pca_map_prior.txt')
        open(newunit=u,file='flex_pca_map_prior.txt',status='replace',action='write')
        write(u,'(A)') '# component covariance_eigenvalue prior_precision'
        do q=1,size(eigvals)
            write(u,'(I5,2(1X,ES20.10))') q,eigvals(q),prior_precision(q)
        end do
        close(u)
    end subroutine write_covariance_tables

    !>  The delivered weight table into one file per state (flex_weights_state_NNN.bin: every physical
    !!  row, rows outside the selection zero, that state's scalars alongside), each registered in the
    !!  out segment of the run's own project copy as imgkind flex_weights, state NNN, beside vol_flex
    !!  state NNN. The files validate their rows against the field's `state` as the run saw it, so
    !!  this runs BEFORE write_discrete_state_project overwrites those labels.
    subroutine write_flex_weights_store( params, build, pinds, weights, labels, targets, bandwidths, l_merged )
        type(parameters), intent(in)    :: params
        type(builder),    intent(inout) :: build
        integer,          intent(in)    :: pinds(:), labels(:)
        real,             intent(in)    :: weights(:,:), targets(:,:), bandwidths(:)
        logical,          intent(in)    :: l_merged
        character(len=STDLEN) :: message
        integer :: status, s, nstates
        nstates = size(weights,2)
        call flex_weights_deliver(build%spproj, build%spproj_field, params%box, params%smpd, &
            &params%box_crop, params%smpd_crop, pinds, weights, labels, targets, bandwidths, &
            &merge(FLEX_WEIGHTS_PROV_MERGED, FLEX_WEIGHTS_PROV_FLEX_PCA, l_merged), status, message)
        if( status /= 0 ) THROW_HARD('flex_pca could not deliver the state weights: '//trim(message))
        do s = 1, nstates
            call build%spproj%add_flex_weights2os_out(flex_weights_state_fname(s), s, params%box, params%smpd)
        end do
        ! entries of a previous delivery with more states would point at files the delivery removed
        do s = nstates+1, nstates+FLEX_WEIGHTS_STALE_SCAN
            if( build%spproj%isthere_in_osout('flex_weights', s) ) call build%spproj%remove_entry_from_osout('flex_weights', s)
        end do
        call build%spproj%write_segment_inside('out', params%projfile)
        write(logfhandle,'(A,I0,A,I0,A)') '>>> FLEX_PCA STATE WEIGHTS WRITTEN: flex_weights_state_001..'// &
            &int2str_pad(nstates,3)//'.bin (', size(weights,1), ' particles x ', nstates, &
            &' states; registered in the out segment as flex_weights per state)'
        call flush(logfhandle)
    end subroutine write_flex_weights_store

    !>  Write the hard state assignment INTO the run's own project: ptcl3D/state carries each embedded
    !!  particle's label, 0 elsewhere. Judge the clusters independently of the kernel-weighted backend with
    !!      simple_exec prg=reconstruct3D projfile=<projfile> nstates=<nstates>
    !!  mkdir=yes already gave the master a private copy of the project, so this rewrites that copy and
    !!  never the project the user pointed at. MUTATES the live field, so it must run after every stage
    !!  that reads the input particle selection, and on the master only.
    subroutine write_discrete_state_project( spproj, pinds, labels, nstates, projfile )
        type(sp_project), intent(inout) :: spproj
        integer,          intent(in)    :: pinds(:), labels(:), nstates
        type(string),     intent(in)    :: projfile
        logical, allocatable :: assigned(:)
        integer :: i, iptcl, state, nptcls, nexcluded
        if( size(pinds) < 1 .or. size(labels) /= size(pinds) .or. nstates < 2 ) &
            &THROW_HARD('invalid flex_pca discrete-state assignment')
        if( len_trim(projfile%to_char()) == 0 ) THROW_HARD('flex_pca discrete-state project file is empty')
        nptcls = spproj%os_ptcl3D%get_noris()
        ! validate BEFORE mutating: this overwrites the live project field rather than a private copy,
        ! so a mid-loop abort would leave the input selection half-replaced
        allocate(assigned(nptcls), source=.false.)
        nexcluded = 0
        do i = 1,size(pinds)
            iptcl = pinds(i)
            state = labels(i)
            if( iptcl < 1 .or. iptcl > nptcls ) THROW_HARD('flex_pca discrete-state particle index outside project')
            if( assigned(iptcl) ) THROW_HARD('duplicate particle in flex_pca discrete-state assignment')
            if( state > nstates ) THROW_HARD('flex_pca discrete-state label outside state range')
            if( state < 1 ) nexcluded = nexcluded + 1
            assigned(iptcl) = .true.
        end do
        do iptcl = 1,nptcls
            call spproj%os_ptcl3D%set_state(iptcl,0)
        end do
        do i = 1,size(pinds)
            if( labels(i) >= 1 ) call spproj%os_ptcl3D%set_state(pinds(i),labels(i))
        end do
        ! ptcl3D only, as refine3D does for nstates>1: a state INDEX carries no ptcl2D meaning, and only
        ! selection (0/1) is mirrored across the two segments
        call spproj%write_segment_inside('ptcl3D', projfile)
        write(logfhandle,'(A,A)') '>>> FLEX_PCA HARD STATES WRITTEN TO: ',projfile%to_char()
        do state = 1,nstates
            write(logfhandle,'(A,I0,A,I0)') '>>> FLEX_PCA DISCRETE-STATE state=',state, &
                &' population=',count(labels==state)
        end do
        if( nexcluded > 0 ) write(logfhandle,'(A,I0)') &
            &'>>> FLEX_PCA DISCRETE-STATE unassigned particles left at state=0: ',nexcluded
        write(logfhandle,'(A,A,A,I0)') '>>> RECONSTRUCT WITH: simple_exec prg=reconstruct3D projfile=', &
            &projfile%to_char(),' nstates=',nstates
        call flush(logfhandle)
        deallocate(assigned)
    end subroutine write_discrete_state_project

    subroutine write_covariance_manifest( params, nptcls, ncomp, nstates, axis, min_neff, sigma_loaded )
        type(parameters), intent(in) :: params
        integer, intent(in) :: nptcls, ncomp, nstates, axis, min_neff
        logical, intent(in) :: sigma_loaded
        integer :: u
        call del_file('flex_pca_manifest.txt')
        open(newunit=u,file='flex_pca_manifest.txt',status='replace',action='write')
        write(u,'(A)') 'method=matched_kb_selected_column_covariance'
        write(u,'(A)') 'diffusion_map_dependency=no'
        write(u,'(A)') 'interpolation=matched_simple_kaiser_bessel'
        write(u,'(A,L1)') 'sigma_whitened=',sigma_loaded
        write(u,'(A,I0)') 'particles=',nptcls
        write(u,'(A,I0)') 'box_crop=',params%box_crop
        write(u,'(A,F10.4)') 'smpd_crop=',params%smpd_crop
        write(u,'(A,F10.4)') 'lowpass_angstrom=',params%lp
        write(u,'(A,I0)') 'column_separation=',params%column_separation
        write(u,'(A,I0)') 'probe_iters=',params%n_probe_iters
        ! the covariance band is capped at fdim(box_crop)-1
        write(u,'(A,L1)') 'lowpass_active=',(params%lp > 2.0*params%smpd_crop)
        write(u,'(A,I0)') 'components=',ncomp
        write(u,'(A,I0)') 'states=',nstates
        if( axis <= 0 )then
            write(u,'(A)')    'state_placement=kmeans_full_latent_space'
        else
            write(u,'(A)')    'state_placement=single_axis_quantiles'
        endif
        write(u,'(A,I0)') 'state_axis=',axis
        write(u,'(A,I0)') 'minimum_state_neff=',min_neff
        write(u,'(A)') 'half_maps=combined_even_odd'
        write(u,'(A)') 'validation=inspect_half_map_agreement_nuisance_correlations_and_heldout_residuals'
        ! ---- cross-fit-FSC provenance (spec par.4.2): a run whose fits consumed the crossfsc
        ! ridge is marked COUPLED, because the coupling deflates the cross-fit statistic's
        ! independence the same bounded way gold-standard FSC regularization does. Lines appear
        ! only when SIMPLE_COV_XFSC_REG is set, so the default manifest is unchanged.
        write(u,'(A,I0)') 'crossfsc_reg_mode=',1
        write(u,'(A,L1)') 'crossfsc_coupled=',.true.
        close(u)
    end subroutine write_covariance_manifest

    ! ============ SELF-CONTAINED TESTS ============ No project, no images, no data files.

    subroutine test_flex_pca_embedding_cache_io()
        integer,  parameter :: NP = 37, NC = 4
        character(len=*), parameter :: FN = 'test_flex_pca_cache.bin'
        integer  :: pinds(NP), i, q, r, ncomp_rd
        real(dp) :: z(NP,NC), eigvals(NC), contrast(NP), re(NP), rme(NP)
        real(dp) :: prec(NC,NC,NP), sig2, sig2_rd
        real(dp), allocatable :: z_rd(:,:), eig_rd(:), con_rd(:), re_rd(:), rme_rd(:), prec_rd(:,:,:)
        write(logfhandle,'(A)') '>>> TEST flex_pca embedding cache I/O'
        do i = 1, NP
            pinds(i) = 3*i + 1                      ! non-contiguous, as a real selection is
            contrast(i) = 0.5d0 + 0.01d0*real(i,dp)
            re(i)       = real(i,dp)
            rme(i)      = 2.d0*real(i,dp)
            do q = 1, NC
                z(i,q) = sin(real(i*q,dp))          ! deterministic, no RNG
                do r = 1, NC
                    prec(q,r,i) = 1.d0/real(q+r+i,dp)
                end do
            end do
        end do
        do q = 1, NC
            eigvals(q) = 10.d0/real(q,dp)
        end do
        sig2 = 0.137d0
        call write_embedding_cache(FN, pinds, NP, NC, z, eigvals, contrast, re, rme, prec, sig2)
        call read_embedding_cache(FN, pinds, NP, ncomp_rd, z_rd, eig_rd, con_rd, &
            &re_rd, rme_rd, prec_rd, sig2_rd)
        if( ncomp_rd /= NC ) THROW_HARD('cache round trip: component count changed')
        if( maxval(abs(z    - z_rd   )) > 0.d0 ) THROW_HARD('cache round trip: latents differ')
        if( maxval(abs(eigvals - eig_rd)) > 0.d0 ) THROW_HARD('cache round trip: eigenvalues differ')
        if( maxval(abs(contrast - con_rd)) > 0.d0 ) THROW_HARD('cache round trip: contrast differs')
        if( maxval(abs(re   - re_rd  )) > 0.d0 ) THROW_HARD('cache round trip: residual energy differs')
        if( maxval(abs(rme  - rme_rd )) > 0.d0 ) THROW_HARD('cache round trip: mean residual energy differs')
        if( maxval(abs(prec - prec_rd)) > 0.d0 ) THROW_HARD('cache round trip: precision differs')
        if( abs(sig2 - sig2_rd)         > 0.d0 ) THROW_HARD('cache round trip: sig2_eff differs')
        call del_file(FN)
        write(logfhandle,'(A)') '>>>   PASSED (bit-exact for all seven payloads)'
    end subroutine test_flex_pca_embedding_cache_io


    !> Smallest even crop that still resolves lp with margin: smpd_crop = smpd*box/box_crop and the
    !! crop's Nyquist is 2*smpd_crop, so lp needs box_crop > 2*box*smpd/lp.
    pure integer function auto_box_crop( box, smpd, lp ) result( bc )
        integer, intent(in) :: box
        real,    intent(in) :: smpd, lp
        if( lp <= 0. .or. smpd <= 0. .or. box <= 0 )then
            bc = box
            return
        endif
        bc = 2*nint(0.5*FLEX_AUTO_BOX_SAFETY*2.0*real(box)*smpd/lp)   ! nearest even
        bc = max(32, min(box, bc))
    end function auto_box_crop

    !> Minimum effective particles per state: the larger of the SNR requirement (~1/s particles for
    !! unit conformational SNR) and an occupancy floor. IgG is limited by the first, Ribosembly the
    !! second, so neither term alone suffices.
    pure integer function auto_min_neff( nptcls, nstates, snr_best ) result( mn )
        integer,  intent(in) :: nptcls, nstates
        real(dp), intent(in) :: snr_best          !< best per-component conformational SNR, 0 if unknown
        integer :: n_snr, n_occ
        n_snr = 20
        if( snr_best > 0.d0 ) n_snr = max(20, nint(1.d0/snr_best))
        n_occ = 20
        if( nstates > 0 ) n_occ = nint(FLEX_AUTO_NEFF_OCCUPANCY*real(nptcls)/real(nstates))
        mn = max(20, min(nptcls, max(n_snr, n_occ)))
    end function auto_min_neff

    !> Over-provisioned starting count for npreimages=0; see the call site for why 32 and not more.
    pure integer function auto_state_count( nptcls, min_neff ) result( k )
        integer, intent(in) :: nptcls, min_neff
        k = FLEX_AUTO_K_START
        if( min_neff > 0 ) k = min(k, nptcls/(4*min_neff))
        k = max(FLEX_AUTO_K_MIN, k)
    end function auto_state_count

    !> The derived settings must reproduce what the validation datasets were actually run at.
    subroutine test_flex_pca_auto_settings()
        integer :: bc, mn, k
        write(logfhandle,'(A)') '>>> TEST flex_pca derived settings'
        ! box_crop: both IgG-RL and Ribosembly are box 128 at 3.0 A/px run at lp=15, box_crop=64
        bc = auto_box_crop(128, 3.0, 15.0)
        if( bc /= 64 ) THROW_HARD('auto box_crop did not reproduce the validated 64')
        ! finer lp must not silently keep a crop that cannot resolve it
        if( auto_box_crop(128, 3.0, 8.0) <= 64 ) THROW_HARD('auto box_crop did not grow with finer lp')
        if( auto_box_crop(128, 3.0, 30.0) >= 64 ) THROW_HARD('auto box_crop did not shrink with coarser lp')
        if( auto_box_crop(128, 3.0, 15.0) > 128 ) THROW_HARD('auto box_crop exceeded the native box')
        ! min_neff: Ribosembly is occupancy-limited, IgG is SNR-limited
        mn = auto_min_neff(335240, 16, 0.d0)
        if( abs(mn - 2095) > 50 ) THROW_HARD('auto min_neff did not reproduce the Ribosembly scale')
        mn = auto_min_neff(100000, 20, 0.0178d0)
        if( mn < 56 ) THROW_HARD('auto min_neff fell below the IgG SNR requirement')
        ! the SNR term must be able to dominate when the signal is weak
        if( auto_min_neff(10000, 20, 1.d-3) <= auto_min_neff(10000, 20, 1.d-1) ) &
            &THROW_HARD('auto min_neff did not grow as conformational SNR fell')
        ! state count: over-provision, but never below the floor or above the validated level
        k = auto_state_count(335240, 2000)
        if( k /= 32 ) THROW_HARD('auto state count did not over-provision to the validated level')
        if( auto_state_count(1000, 2000) /= FLEX_AUTO_K_MIN ) &
            &THROW_HARD('auto state count ignored the small-dataset floor')
        if( auto_state_count(100000000, 100) /= FLEX_AUTO_K_START ) &
            &THROW_HARD('auto state count exceeded the validated over-provision cap')
        write(logfhandle,'(A,I0,A,I0,A)') '>>>   PASSED (box_crop=',auto_box_crop(128,3.0,15.0), &
            &', state count=',auto_state_count(335240,2000),', min_neff both regimes)'
    end subroutine test_flex_pca_auto_settings

    subroutine place_states_with_population_floor( z, nptcls, ncomp, nkern, nstates_req, axis, min_neff, &
        &min_state_frac, eigvals, precision, weights, targets, bandwidths, neff, labels, comp_rho )
        use simple_rnd, only: irnd_uni
        integer,  intent(in) :: nptcls, ncomp, nkern, nstates_req, axis, min_neff
        real,     intent(in) :: min_state_frac
        real(dp), intent(in) :: z(nptcls,ncomp), eigvals(ncomp), precision(ncomp,ncomp,nptcls)
        real,    allocatable, intent(out) :: weights(:,:), targets(:,:), bandwidths(:), neff(:)
        integer, allocatable, intent(out) :: labels(:)
        real(dp), optional,   intent(in)  :: comp_rho(ncomp)
        integer,  parameter :: ROUND_CAP = 8
        real,     allocatable :: w_r(:,:), t_r(:,:), bw_r(:), nf_r(:)
        real(dp), allocatable :: z_r(:,:), p_r(:,:,:), sdv(:)
        integer,  allocatable :: idx(:), lab_r(:), occ(:), order(:), kept(:), deliver(:)
        logical,  allocatable :: retained(:), qualifies(:)
        integer  :: nmin, K, round, nret, nk, i, q, s, t, nqual, nrand, nsurplus, kbest, itmp
        real(dp) :: d2, dbest, zbar
        logical  :: l_success
        nk   = max(1, min(ncomp, nkern))
        nmin = max(1, nint(min_state_frac * real(nptcls)))
        if( nstates_req * nmin > nptcls )then
            THROW_HARD('min_state_frac is too large for the requested state count: the floors exceed the particle count')
        endif
        allocate(retained(nptcls), source=.true.)
        K         = nstates_req
        nret      = nptcls
        nqual     = 0
        l_success = .false.
        do round = 1, ROUND_CAP
            nret = count(retained)
            if( allocated(idx) ) deallocate(idx)
            allocate(idx(nret))
            t = 0
            do i = 1, nptcls
                if( retained(i) )then
                    t = t + 1
                    idx(t) = i
                endif
            end do
            allocate(z_r(nret,ncomp), p_r(ncomp,ncomp,nret))
            do i = 1, nret
                z_r(i,:)   = z(idx(i),:)
                p_r(:,:,i) = precision(:,:,idx(i))
            end do
                call build_covariance_state_weights(z_r, nret, ncomp, nkern, K, axis, &
                    &max(20, min(min_neff, nret/2)), eigvals, p_r, w_r, t_r, bw_r, nf_r, lab_r, &
                    &comp_rho=comp_rho)
            deallocate(z_r, p_r)
            if( allocated(occ) ) deallocate(occ, qualifies)
            allocate(occ(K), source=0)
            allocate(qualifies(K), source=.false.)
            do i = 1, nret
                if( lab_r(i) >= 1 ) occ(lab_r(i)) = occ(lab_r(i)) + 1
            end do
            qualifies = occ >= nmin
            nqual     = count(qualifies)
            write(logfhandle,'(A,I0,A,I0,A,I0,A,I0,A,I0)') '>>> FLEX_PCA POPULATION FLOOR round=', round, &
                &' provisioned=', K, ' retained=', nret, ' qualifying=', nqual, ' floor=', nmin
            if( nqual >= nstates_req )then
                l_success = .true.
                exit
            endif
            ! peel: members of under-populated clusters and particles outside every kernel support
            ! leave the placement mass, so the next placement spends its centers on the retained mass
            do i = 1, nret
                if( lab_r(i) < 1 )then
                    retained(idx(i)) = .false.
                else if( .not. qualifies(lab_r(i)) )then
                    retained(idx(i)) = .false.
                endif
            end do
            if( count(retained) < nstates_req * nmin )then
                write(logfhandle,'(A)') '>>> FLEX_PCA POPULATION FLOOR: the retained mass can no longer hold the floors'
                exit
            endif
            if( K >= AUTO_NSTATES )then
                write(logfhandle,'(A)') '>>> FLEX_PCA POPULATION FLOOR: provision cap reached'
                exit
            endif
            K = min(AUTO_NSTATES, K + (nstates_req - nqual))
        end do
        if( .not. l_success )then
            THROW_WARN('flex_pca population floor not reached for every requested state; keeping the most populated clusters')
        endif
        ! clusters ordered by population, descending (K <= AUTO_NSTATES, selection sort)
        allocate(order(K))
        order = [(s, s=1,K)]
        do s = 1, K-1
            kbest = s
            do t = s+1, K
                if( occ(order(t)) > occ(order(kbest)) ) kbest = t
            end do
            if( kbest /= s )then
                itmp         = order(s)
                order(s)     = order(kbest)
                order(kbest) = itmp
            endif
        end do
        allocate(kept(nstates_req), deliver(K))
        kept    = order(1:nstates_req)
        deliver = 0
        do s = 1, nstates_req
            deliver(kept(s)) = s
        end do
        ! standardized latent metric on the placement components, for attaching surplus clusters
        allocate(sdv(nk))
        do q = 1, nk
            zbar   = sum(z(:,q)) / real(nptcls,dp)
            sdv(q) = max(sqrt(sum((z(:,q) - zbar)**2) / real(nptcls,dp)), 1.d-12)
        end do
        allocate(labels(nptcls), source=0)
        nsurplus = 0
        do i = 1, nret
            s = lab_r(i)
            if( s < 1 ) cycle
            if( deliver(s) > 0 )then
                labels(idx(i)) = deliver(s)
            else if( qualifies(s) )then
                ! surplus qualifying cluster: real mass, attached to the nearest delivered target
                dbest = huge(1.d0)
                kbest = 1
                do t = 1, nstates_req
                    d2 = 0.d0
                    do q = 1, nk
                        d2 = d2 + ((z(idx(i),q) - real(t_r(q,kept(t)),dp)) / sdv(q))**2
                    end do
                    if( d2 < dbest )then
                        dbest = d2
                        kbest = t
                    endif
                end do
                labels(idx(i)) = kbest
                nsurplus       = nsurplus + 1
            endif
        end do
        ! members of dropped clusters, particles outside every kernel support and particles peeled in
        ! earlier rounds receive a uniformly random delivered label
        nrand = 0
        do i = 1, nptcls
            if( labels(i) < 1 )then
                labels(i) = irnd_uni(nstates_req)
                nrand     = nrand + 1
            endif
        end do
        ! delivered tables: hard-label indicator weights, so the state maps are ordinary
        ! reconstructions of the labelled particles
        allocate(weights(nptcls,nstates_req), source=0.)
        do i = 1, nptcls
            weights(i,labels(i)) = 1.
        end do
        allocate(targets(ncomp,nstates_req), bandwidths(nstates_req), neff(nstates_req))
        do s = 1, nstates_req
            targets(:,s)  = t_r(:,kept(s))
            bandwidths(s) = bw_r(kept(s))
            neff(s)       = real(count(labels == s))
        end do
        write(logfhandle,'(A,I0,A,I0,A,I0)') '>>> FLEX_PCA POPULATION FLOOR delivered states=', nstates_req, &
            &' surplus-attached=', nsurplus, ' randomized=', nrand
        do s = 1, nstates_req
            write(logfhandle,'(A,I3,A,I9,A,I9)') '>>>   state=', s, '  particles=', nint(neff(s)), '  floor=', nmin
            if( nint(neff(s)) < nmin ) THROW_WARN('flex_pca delivered a state below the population floor')
        end do
        call flush(logfhandle)
        deallocate(retained, idx, w_r, t_r, bw_r, nf_r, lab_r, occ, qualifies, order, kept, deliver, sdv)
    end subroutine place_states_with_population_floor

    subroutine test_flex_pca_population_floor()
        use simple_rnd, only: seed_rnd
        integer,  parameter :: NA = 300, NB = 250, NO = 20, NP = NA+NB+NO, NC = 2, NST = 2
        real,     parameter :: FRAC = 0.2
        integer  :: i, q, state, nmin, nlab(NST)
        real(dp) :: z(NP,NC), eigvals(NC), prec(NC,NC,NP)
        real,     allocatable :: weights(:,:), targets(:,:), bandwidths(:), neff(:)
        integer,  allocatable :: labels(:)
        write(logfhandle,'(A)') '>>> TEST flex_pca population floor on two clusters plus outliers'
        call seed_rnd
        do i = 1, NA
            z(i,1) = -3.d0 + 0.01d0*real(mod(i,7),dp)
            z(i,2) =  0.02d0*real(mod(i,5),dp)
        end do
        do i = 1, NB
            z(NA+i,1) = 3.d0 + 0.01d0*real(mod(i,7),dp)
            z(NA+i,2) = 0.02d0*real(mod(i,5),dp)
        end do
        do i = 1, NO
            z(NA+NB+i,1) = 60.d0 + 0.05d0*real(mod(i,3),dp)
            z(NA+NB+i,2) = 60.d0 + 0.05d0*real(mod(i,4),dp)
        end do
        eigvals = 1.d0
        prec    = 0.d0
        do i = 1, NP
            do q = 1, NC
                prec(q,q,i) = 1.d0
            end do
        end do
        call place_states_with_population_floor(z, NP, NC, NC, NST, 0, 10, FRAC, eigvals, prec, &
            &weights, targets, bandwidths, neff, labels)
        nmin = max(1, nint(FRAC*real(NP)))
        if( size(labels) /= NP ) THROW_HARD('labels shape wrong')
        if( any(labels < 1) .or. any(labels > NST) ) THROW_HARD('a particle was left without a delivered state')
        if( size(weights,1) /= NP .or. size(weights,2) /= NST ) THROW_HARD('weights shape wrong')
        do i = 1, NP
            if( abs(sum(weights(i,:)) - 1.) > 1.e-6 ) THROW_HARD('delivered weights are not hard-label indicators')
        end do
        nlab = 0
        do i = 1, NP
            nlab(labels(i)) = nlab(labels(i)) + 1
        end do
        do state = 1, NST
            if( nlab(state) < nmin ) THROW_HARD('a delivered state is below the population floor')
            if( nint(neff(state)) /= nlab(state) ) THROW_HARD('neff does not report the delivered population')
        end do
        deallocate(weights, targets, bandwidths, neff, labels)
        write(logfhandle,'(A)') '>>>   PASSED (every particle labelled, every state at or above the floor)'
    end subroutine test_flex_pca_population_floor


    !> The resume embedding path itself ('' when not resuming)
    function infile_path( params, l_resume ) result( f )
        type(parameters), intent(in) :: params
        logical,          intent(in) :: l_resume
        character(len=:), allocatable :: f
        f = ''
        if( l_resume ) f = params%infile%to_char()
    end function infile_path

    !> Directory of the resume embedding (params%infile), '' when not resuming or when the path has no
    !! directory part; used to adopt the original run's deconvolved cache from a fresh run directory.
    function infile_dir( params, l_resume ) result( d )
        type(parameters), intent(in) :: params
        logical,          intent(in) :: l_resume
        character(len=:), allocatable :: d, f
        integer :: k
        d = ''
        if( .not. l_resume ) return
        f = params%infile%to_char()
        k = index(f, '/', back=.true.)
        if( k > 1 ) d = f(1:k-1)
    end function infile_dir

    !> The delivery UMAP coloured by the delivered states (argmax state weight; 0 = no weight), written
    !! as flex_pca_umap_states.jpg next to the state tables. Needs the coordinates kept by the last
    !! deliver_latent_readouts of this process; a states-only worker that never ran it writes nothing.
    subroutine write_states_umap_figure( pinds, z, weights )
        integer,  intent(in) :: pinds(:)
        real(dp), intent(in) :: z(:,:)
        real,     intent(in) :: weights(:,:)
        integer, allocatable :: lut(:), lab(:)
        real,    allocatable :: pz1(:), pz2(:)
        integer :: i, row, nsub, pmax, s, sbest
        real    :: wbest
        if( .not. allocated(umap_plot_pind) ) return
        if( size(z,2) < 2 ) return
        nsub = size(umap_plot_pind); pmax = max(maxval(pinds), maxval(umap_plot_pind))
        allocate(lut(pmax), source=0)
        do i = 1, size(pinds); lut(pinds(i)) = i; end do
        allocate(lab(nsub), pz1(nsub), pz2(nsub))
        do i = 1, nsub
            row = lut(umap_plot_pind(i))
            if( row < 1 )then
                lab(i) = 0; pz1(i) = 0.0; pz2(i) = 0.0; cycle
            endif
            pz1(i) = real(z(row,1)); pz2(i) = real(z(row,2))
            sbest = 0; wbest = 0.0
            do s = 1, size(weights,2)
                if( weights(row,s) > wbest )then; wbest = weights(row,s); sbest = s; endif
            end do
            lab(i) = sbest
        end do
        call flex_plot_latent_jpg('flex_pca_umap_states.jpg', umap_plot_xy(1,:), umap_plot_xy(2,:), &
            &pz1, pz2, nsub, labels=lab, nlab=size(weights,2))
        write(logfhandle,'(A)') '>>> FLEX_PCA UMAP figure by delivered state: flex_pca_umap_states.jpg'
        deallocate(lut, lab, pz1, pz2)
    end subroutine write_states_umap_figure

    !> Delivery readout on the all-N latents: UMAP plot coordinates of a bounded subsample
    !! (flex_pca_umap.txt), opt-in through umap=yes. Dead axes (near-zero latent variance) are
    !! excluded, since whitening explodes them into pure noise.
    subroutine deliver_latent_readouts( pinds, nptcls, ncomp, z, l_umap, tag, labels )
        integer,  intent(in) :: pinds(:), nptcls, ncomp
        real(dp), intent(in) :: z(:,:)
        logical,  intent(in) :: l_umap
        character(len=*), optional, intent(in) :: tag   !< suffix on the output names (e.g. '_deconv')
        integer,          optional, intent(in) :: labels(:) !< per-row labels (e.g. the deconvolution's mixture component) for the figure
        real,    allocatable :: pz1(:), pz2(:)
        integer, allocatable :: plab(:)
        character(len=:), allocatable :: sfx
        integer,  allocatable :: hcols(:), uinds(:)
        real,     allocatable :: Xu(:,:), Yu(:,:)
        real(dp) :: colvar(ncomp), vmax, mu_c
        integer  :: vumap, nhc, q, i, u, nsub
        sfx = ''
        if( present(tag) ) sfx = tag
        vumap = merge(1, 0, l_umap)
        if( vumap < 1 ) return
        ! healthy columns: latent variance above 1e-4 of the max
        vmax = 0.d0
        do q = 1, ncomp
            mu_c      = sum(z(1:nptcls,q))/real(nptcls,dp)
            colvar(q) = sum((z(1:nptcls,q) - mu_c)**2)/real(max(nptcls-1,1),dp)
            vmax      = max(vmax, colvar(q))
        end do
        allocate(hcols(ncomp))
        nhc = 0
        do q = 1, ncomp
            if( colvar(q) > 1.d-4*vmax )then
                nhc = nhc + 1
                hcols(nhc) = q
            endif
        end do
        if( nhc < ncomp ) write(logfhandle,'(A,I0,A,I0,A)') '>>> FLEX_PCA delivery readouts: ', &
            &ncomp - nhc, ' dead/retired axis(es) excluded (', nhc, ' healthy)'
        if( vumap > 0 .and. nhc >= 3 .and. nptcls > 100 )then
            nsub = min(1000000, nptcls)   ! kd-tree kNN + 200 SGD epochs: every particle, not a subsample
            call umap_subsample(nptcls, nsub, 1234, uinds)
            allocate(Xu(nsub,nhc))
            do q = 1, nhc
                do i = 1, nsub
                    Xu(i,q) = real(z(uinds(i),hcols(q)))
                end do
            end do
            call umap_embed(Xu, 1234, Yu)
            write(logfhandle,'(A,I0,A)') '>>> FLEX_PCA UMAP delivery plot: ', nsub, &
                &' latents -> flex_pca_umap.txt'
            call flush(logfhandle)
            call del_file('flex_pca_umap'//sfx//'.txt')
            open(newunit=u, file='flex_pca_umap'//sfx//'.txt', status='replace', action='write')
            write(u,'(A)') '# particle  umap1  umap2'
            do i = 1, nsub
                write(u,'(I10,2(1X,ES14.6))') pinds(uinds(i)), Yu(i,1), Yu(i,2)
            end do
            close(u)
            ! the figure: log-density UMAP core, UMAP by label, z1 vs z2 by label (flex_pca_umap<sfx>.jpg);
            ! the coordinates are kept for the state stage's own colouring
            if( allocated(umap_plot_xy) )   deallocate(umap_plot_xy)
            if( allocated(umap_plot_pind) ) deallocate(umap_plot_pind)
            allocate(umap_plot_xy(2,nsub), umap_plot_pind(nsub), pz1(nsub), pz2(nsub), plab(nsub))
            do i = 1, nsub
                umap_plot_xy(1,i) = Yu(i,1); umap_plot_xy(2,i) = Yu(i,2)
                umap_plot_pind(i) = pinds(uinds(i))
                pz1(i) = real(z(uinds(i),hcols(1))); pz2(i) = real(z(uinds(i),hcols(2)))
                plab(i) = 0
                if( present(labels) ) plab(i) = max(0, labels(uinds(i)))
            end do
            call flex_plot_latent_jpg('flex_pca_umap'//sfx//'.jpg', umap_plot_xy(1,:), umap_plot_xy(2,:), &
                &pz1, pz2, nsub, labels=plab)
            write(logfhandle,'(A)') '>>> FLEX_PCA UMAP figure: flex_pca_umap'//sfx//'.jpg'
            deallocate(Xu, Yu, uinds, pz1, pz2, plab)
        endif
        deallocate(hcols)
    end subroutine deliver_latent_readouts

    !> Whether the master resolved real sigma2 spectra or fell back to the unit white-noise spectrum.
    !! Discovery is re-run independently in every worker, and a worker that resolves it differently
    !! whitens against a different noise model -- which changes every inner product without any error
    !! being raised, because both outcomes are individually legitimate. Pin it and make workers assert.
    subroutine write_sigma_state( loaded )
        logical, intent(in) :: loaded
        integer :: funit, io_stat
        call fopen(funit, file=string(SIGMA_STATE_FNAME), action='WRITE', status='REPLACE', iostat=io_stat)
        call fileiochk('write_sigma_state; open', io_stat)
        write(funit,'(L1)') loaded
        call fclose(funit)
    end subroutine write_sigma_state

    subroutine check_sigma_state( loaded )
        logical, intent(in) :: loaded
        integer :: funit, io_stat
        logical :: loaded_master
        if( .not. file_exists(string(SIGMA_STATE_FNAME)) )then
            THROW_HARD('flex_pca worker found no flex_pca_sigma_state.txt from the master')
        endif
        call fopen(funit, file=string(SIGMA_STATE_FNAME), action='READ', status='OLD', iostat=io_stat)
        call fileiochk('check_sigma_state; open', io_stat)
        read(funit,'(L1)') loaded_master
        call fclose(funit)
        if( loaded_master .neqv. loaded )then
            THROW_HARD('flex_pca worker resolved sigma2 differently from the master; the parts would be &
                &whitened against different noise models')
        endif
    end subroutine check_sigma_state

end module simple_flex_pca_model
