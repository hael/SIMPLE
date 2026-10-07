!@descr: flex_pca application: the phase order of a run and the worker's stage dispatch
!!
!! One object per process. The master/shared-memory run is prepare -> obtain_embedding ->
!! infer_states -> reconstruct_states -> publish over one flex_run_session, which owns every
!! value that outlives a phase (selection, mean and basis, embedding and its moments, the state
!! table) and frees them once. The worker runs exactly one stage body. Distribution goes through
!! the rounds object the strategy hands in; every switch comes from the run settings.
module simple_flex_pca_application
use simple_core_module_api
use simple_flex_pca_records, only: flex_fit_model, flex_selection
use simple_builder,                    only: builder
use simple_cmdline,                    only: cmdline
use simple_flex_pca_rec3D,         only: reconstruct_flex_weighted_states, flex_rec_box, flex_rec_smpd, &
    &read_state_weights_round
use simple_flex_pca_planes,        only: planes_enable, planes_kill
use simple_flex_pca_plane_cache,   only: plane_cache_ensure, plane_cache_in_use, plane_cache_close
use simple_flex_pca_rounds,        only: flex_pca_rounds
use simple_flex_pca_embedding_io,  only: write_embedding_cache, read_embedding_cache
use simple_flex_pca_state_service, only: apply_latent_deconvolution, prune_underpopulated_states, &
    &place_states_with_population_floor, infile_path, infile_dir, AUTO_NSTATES
use simple_flex_pca_delivery_3d,   only: deliver_latent_readouts, write_covariance_tables, &
    &write_covariance_eigenvolumes, write_covariance_manifest, apply_consensus_nu_filter
use simple_flex_pca_project_gateway, only: validate_covariance_inputs, load_and_validate_sigma, &
    &write_flex_weights_store, write_discrete_state_project
use simple_flex_pca_stages,        only: flex_stage_request, PCA_STAGE_STATES, PCA_STAGE_EMBED, PCA_STAGE_PROBE, &
    &PCA_STAGE_POLISH
use simple_flex_pca_run_types,     only: flex_run_settings, flex_run_session
use simple_flex_pca_mstep,         only: init_basis_reconstructor
use simple_flex_pca_basis,         only: estimate_covariance_mean, &
    &align_basis_to_reference, save_probe_state
use simple_flex_pca_embed,         only: embed_latents_with_contrast, compose_basis_from_runs, compose_cut_reembed
use simple_flex_probe_fit,         only: probe_subspace_iteration
use simple_flex_pca_fit_driver,    only: run_flex_pca_paired, probe_worker_pass, embed_worker_pass
use simple_image,                      only: image
use simple_parameters,                 only: parameters
use simple_reconstructor,              only: reconstructor
use simple_flex_pca_merge,             only: two_gate_state_merge
use simple_flex_pca_weights, only: build_covariance_state_weights, cv_select_bandwidths
use simple_flex_pca_targets, only: component_reliability_proxy
implicit none
private
#include "simple_local_flags.inc"

public :: flex_pca_application

integer, parameter :: MIN_NSTATES = 3

type :: flex_pca_application
    type(flex_run_session) :: sess
  contains
    procedure :: run        => app_run
    procedure :: run_worker => execute_worker_stage
    procedure :: kill       => app_kill
    procedure, private :: prepare, obtain_embedding, infer_states, reconstruct_states, publish
    procedure, private :: resume_from_cache, embed_all
end type flex_pca_application

contains

    !> The master / shared-memory run: the five phases in order, then the session is released
    subroutine app_run( self, params, build, cline, rounds )
        class(flex_pca_application), intent(inout) :: self
        class(flex_pca_rounds), intent(inout) :: rounds
        type(parameters), intent(inout) :: params
        type(builder),    intent(inout) :: build
        class(cmdline),   intent(inout) :: cline
        call self%prepare(params, build, cline, rounds)
        call self%obtain_embedding(params, build, rounds)
        call self%infer_states(params, build, cline, rounds)
        call self%reconstruct_states(params, build, rounds)
        call self%publish(params, build, rounds)
        call self%kill
    end subroutine app_run

    subroutine app_kill( self )
        class(flex_pca_application), intent(inout) :: self
        call planes_kill()
        call plane_cache_close()
        call self%sess%kill
    end subroutine app_kill

    !> Selection, sigma, planes and the derived state settings; the run settings are resolved here
    subroutine prepare( self, params, build, cline, rounds )
        class(flex_pca_application), intent(inout) :: self
        class(flex_pca_rounds), intent(inout) :: rounds
        type(parameters), intent(inout) :: params
        type(builder),    intent(inout) :: build
        class(cmdline),   intent(inout) :: cline
        associate( sess => self%sess )
        ! every switch of the run, resolved once: the command-line facts params drops and the
        ! SIMPLE_COV_* environment overrides
        call sess%cfg%new(params, cline)
        call validate_covariance_inputs(params, sess%cfg, build, sess%sel%pinds, sess%sel%nptcls, rounds=rounds)
        sess%l_paired_states = .false.
        sess%l_pop_floor     = .false.
        sess%l_merged        = .false.
        ! distributed master: one particle-index list per part, shipped as pindfile= (no-op otherwise)
        call rounds%plan_partitions(params, sess%sel%pinds)
        call load_and_validate_sigma(params, build, cline, sess%sel%pinds, sess%sigma_loaded, rounds=rounds)
        ! a shared-memory run keeps its prepped planes across every pass; the distributed master
        ! never runs a pass itself and a worker is a fresh process per round
        ! cache=yes: the plane cache is built once by the process that owns the run (workers adopt it)
        if( .not. rounds%is_worker() ) call plane_cache_ensure(params, build, sess%sel%pinds, sess%sel%nptcls)
        if( .not. rounds%distributed() .and. .not. rounds%is_worker() )then
            call planes_enable(build%spproj_field%get_noris(), plane_cache_in_use(params, build))
        endif

        sess%neigs_req  = max(1, min(48, params%neigs))
        ! npreimages = state ceiling (>= MIN_NSTATES); preimage_auto=yes raises it to AUTO_NSTATES unless
        ! npreimages is given, and turns the merge on.
        sess%states%nstates = max(MIN_NSTATES, params%npreimages)
        ! population floor (min_state_frac > 0): exactly npreimages hard-labelled states, each above
        ! the floor; it cannot be combined with the automatic ceiling or the two-gate merge
        if( params%min_state_frac > 0. )then
            if( params%l_preimage_auto .or. sess%cfg%l_merge )then
                THROW_HARD('min_state_frac delivers exactly npreimages states; it cannot be combined with preimage_auto or the merge')
            endif
        endif
        if( params%l_preimage_auto )then
            if( .not. sess%cfg%l_npreimages_explicit ) sess%states%nstates = AUTO_NSTATES
            ! a ceiling without the collapse is just a large state count, so auto implies the merge
            ! (resolved in the run settings; an explicit SIMPLE_COV_MERGE=0 still wins)
            write(logfhandle,'(A,I0,A)') '>>> FLEX_PCA preimage_auto=yes: state count is a CEILING of ', &
                &sess%states%nstates, '; the two-gate merge collapses these to the recovered count'
            if( .not. sess%cfg%l_merge ) write(logfhandle,'(A)') &
                &'>>> FLEX_PCA WARNING: SIMPLE_COV_MERGE=0 disables the merge, so the ceiling will be &
                &delivered as-is; drop preimage_auto and set npreimages explicitly instead'
        endif
        call report_state_memory(params, sess%states%nstates)
        ! the minimum effective sample size of a delivered state: the kernel bandwidth floor AND
        ! the occupancy floor that decides whether a state is reconstructed at all
        sess%min_neff = max(20, min(sess%sel%nptcls, params%min_neff))
        sess%state_axis = params%state_axis      ! <0 path, 0 state_placement (diffusion k-center | equal_occ), >=1 single axis
        ! nkern decouples the number of components the STATE STAGE uses from neigs, the number estimated.
        sess%nkern      = params%nkern
        if( sess%nkern <= 0 ) sess%nkern = ishft(huge(1), -1) ! clamped against ncomp once the fit is known
        sess%col_sep    = max(1, params%column_separation)

        write(logfhandle,'(A,I0,A,I0)') '>>> FLEX_PCA particles=',sess%sel%nptcls, &
            &' requested_components=',sess%neigs_req
        write(logfhandle,'(A,L1,A,I0,A,A,A,I0)') '>>> FLEX_PCA sigma_whitened=',sess%sigma_loaded, &
            &' state_axis=',sess%state_axis,' state_placement=',trim(params%state_placement),' minimum_state_neff=',sess%min_neff
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
        end associate
    end subroutine prepare

    !> The basis and the raw embedding: from the cache (resume), from finished runs (composition)
    !! or from the paired fit with its merge and polish; then one all-N embedding, cached and read
    !! out; then the component reliability and the latent deconvolution.
    subroutine obtain_embedding( self, params, build, rounds )
        class(flex_pca_application), intent(inout) :: self
        class(flex_pca_rounds), intent(inout) :: rounds
        type(parameters), intent(inout) :: params
        type(builder),    intent(inout) :: build
        associate( sess => self%sess )
        sess%l_resume  = sess%cfg%l_resume
        sess%l_compose = sess%cfg%l_compose
        if( sess%l_resume )then
            ! the basis and embedding dominate runtime and do not depend on the state stage, so
            ! infile= re-runs only that stage
            call self%resume_from_cache(params)
        else
            call estimate_covariance_mean(params, build, sess%model%mean_rec, sess%sel%pinds, sess%sel%nptcls, rounds=rounds)
            if( sess%l_compose )then
                ! the coarse-first orthonormalised union of finished runs' eigenvolumes: no
                ! covariance moment, no EM, straight to the embedding
                call compose_basis_from_runs(params, sess%cfg, build, sess%model)
                if( sess%cfg%l_compose_cut ) call compose_cut_reembed(params, sess%cfg, build, sess%model, sess%sel, rounds)
                write(logfhandle,'(A,I0,A)') '>>> FLEX_PCA composed basis: ',sess%model%ncomp, &
                    &' components; no covariance moment, no EM (n_probe_iters ignored)'
                call flush(logfhandle)
            else
                ! two disjoint mod-4-half fits advanced by one shared loop, merged, axis-weighted
                ! and polished (combine-then-polish: no restart, frozen rank and axes). The
                ! full-selection mean built above is not consumed: each fit scales its own.
                if( rounds%nparts() > 1 )then
                    write(logfhandle,'(A,I0,A)') '>>> FLEX_PCA PAIRED DISTRIBUTED: ', &
                        &rounds%nparts(), ' parts, one part per worker per iteration'
                    call flush(logfhandle)
                endif
                ! Paired fit (every non-compose run): two mod-4-half fits advanced together, merged, polished once,
                ! embedded on all N.
                call run_flex_pca_paired(params, sess%cfg, build, sess%sel, sess%col_sep, sess%neigs_req, sess%model, rounds=rounds)
                sess%l_paired_states = .true.
            endif
            call sess%clamp_state_axis
            call self%embed_all(params, build, rounds)
        endif
        ! resuming skips the split-half solve: fall back to the spread-over-posterior-variance proxy
        if( .not. allocated(sess%latent%comp_rho) )then
            allocate(sess%latent%comp_rho(sess%model%ncomp))
            call component_reliability_proxy(sess%latent%z, sess%latent%precision, sess%sel%nptcls, sess%model%ncomp, sess%latent%comp_rho)
            write(logfhandle,'(A)') '>>> FLEX_PCA resumed embedding: component reliability from the &
                &spread/posterior-variance proxy (no cached split-half rho)'
        endif
        ! Latent deconvolution. Must stay after the cache write: the cache stores raw z, so a
        ! resume never deconvolves twice.
        call apply_latent_deconvolution(sess%latent, sess%model, sess%sel, &
            &sess%l_deconv_applied, sess%deconv_labels, sess%l_resume, sess%l_deconv_adopted, &
            &infile_dir(params, sess%l_resume), infile_path(params, sess%l_resume))
        if( sess%l_deconv_applied )then
            ! a resume that adopted the cache already has the readouts from the original run
            if( .not. sess%l_deconv_adopted )then
                if( allocated(sess%deconv_labels) )then
                    call deliver_latent_readouts(sess%readout, sess%sel, sess%model, sess%latent, params%l_umap, tag='_deconv', labels=sess%deconv_labels)
                else
                    call deliver_latent_readouts(sess%readout, sess%sel, sess%model, sess%latent, params%l_umap, tag='_deconv')
                endif
            endif
        endif
        end associate
    end subroutine obtain_embedding

    !> Resume: the raw embedding, its basis and statistics from the cache a previous run wrote.
    subroutine resume_from_cache( self, params )
        class(flex_pca_application), intent(inout) :: self
        type(parameters), intent(inout) :: params
        integer :: q
        associate( sess => self%sess )
        call read_embedding_cache(params%infile%to_char(), params%box_crop, params%smpd_crop, sess%sel, sess%model, sess%latent)
        allocate(sess%latent%prior_precision(sess%model%ncomp))
        do q = 1, sess%model%ncomp
            sess%latent%prior_precision(q) = 1.d0 / max(sess%model%eigvals(q), DTINY)
        end do
        call sess%clamp_state_axis
        write(logfhandle,'(A,A)') '>>> FLEX_PCA RESUMED from embedding cache ',params%infile%to_char()
        write(logfhandle,'(A,I0,A,I0)') '>>> FLEX_PCA cached particles=',sess%sel%nptcls,' components=',sess%model%ncomp
        call flush(logfhandle)
        end associate
    end subroutine resume_from_cache

    !> One all-N embedding against the session's basis: the prior precision, the MAP solve with
    !! per-particle contrast (a qsys round when distributed and the basis is on disk under the plain
    !! namespace; in-process after the paired merge, whose polished basis the embed workers do not
    !! load), the cache and the raw readouts.
    subroutine embed_all( self, params, build, rounds )
        class(flex_pca_application), intent(inout) :: self
        class(flex_pca_rounds), intent(inout) :: rounds
        type(parameters), intent(inout) :: params
        type(builder),    intent(inout) :: build
        integer :: q
        associate( sess => self%sess )
        if( allocated(sess%latent%prior_precision) ) deallocate(sess%latent%prior_precision)
        allocate(sess%latent%prior_precision(sess%model%ncomp))
        do q = 1, sess%model%ncomp
            sess%latent%prior_precision(q) = 1.d0 / max(sess%model%eigvals(q), DTINY)
        end do
        if( sess%l_paired_states )then
            write(logfhandle,'(A,I0,A,I0,A)') '>>> FLEX_PCA MERGED EMBED: all-N pass over ', &
                &sess%sel%nptcls,' selected particles (both halves) against the ',sess%model%ncomp, &
                &'-component merged basis'
            write(logfhandle,'(A)') '>>> FLEX_PCA MERGED EMBED: plain prior; axis reliability is &
                &the merge match cosine already in eigvals (split-half rho not applied)'
            call flush(logfhandle)
            call embed_latents_with_contrast(params, build, sess%model, sess%sel, sess%latent, rounds, l_zhalf=.true.)
        else
            write(logfhandle,'(A,ES12.4,A,ES12.4)') '>>> FLEX_PCA covariance eigenvalues: max=', &
                &maxval(sess%model%eigvals),' min=',minval(sess%model%eigvals)
            call flush(logfhandle)
            ! one qsys round: the basis is fixed. Workers cannot finish -- the reliability prior
            ! couples every particle -- so they ship sufficient statistics and the master owns rho
            ! and the re-solve.
            if( rounds%distributed() )then
                call save_probe_state(sess%model%ncomp, sess%model%eigvals, sess%model%sig2_eff)
                call rounds%run_stage(params, flex_stage_request(stage=PCA_STAGE_EMBED, label='embedding'))
                call embed_latents_with_contrast(params, build, sess%model, sess%sel, sess%latent, rounds, from_parts=.true., l_zhalf=.true.)
            else
                call embed_latents_with_contrast(params, build, sess%model, sess%sel, sess%latent, rounds, l_zhalf=.true.)
            endif
        endif
        write(logfhandle,'(A,F7.3,A,F7.3)') '>>> FLEX_PCA per-particle contrast: mean=', &
            &real(sum(sess%latent%contrast)/real(sess%sel%nptcls,dp)),' sd=', &
            &real(sqrt(max(sum((sess%latent%contrast-sum(sess%latent%contrast)/sess%sel%nptcls)**2)/real(sess%sel%nptcls,dp),DTINY)))
        call flush(logfhandle)
        ! z is left in the physical units the MAP solve returns; the kernel metric does the weighting
        if( .not. sess%l_paired_states ) call write_covariance_eigenvolumes(sess%model%eigvals, sess%model%ncomp)
        call write_embedding_cache('flex_pca_embedding.bin', params%box_crop, params%smpd_crop, sess%sel, sess%model, sess%latent)
        ! raw readouts (the cache holds raw z); the deconvolved coordinates are what the run
        ! delivers and what the state stage re-derives on resume
        call deliver_latent_readouts(sess%readout, sess%sel, sess%model, sess%latent, params%l_umap)
        if( sess%l_paired_states )then
            write(logfhandle,'(A,I0,A)') '>>> FLEX_PCA PAIRED MERGE delivered; continuing into &
                &the state stage with the ', sess%model%ncomp, '-component merged basis'
            call flush(logfhandle)
        endif
        end associate
    end subroutine embed_all

    !> State placement on the (deconvolved) embedding, the occupancy floor, the reconstruction band and the bandwidth CV
    subroutine infer_states( self, params, build, cline, rounds )
        class(flex_pca_application), intent(inout) :: self
        class(flex_pca_rounds), intent(inout) :: rounds
        type(parameters), intent(inout) :: params
        type(builder),    intent(inout) :: build
        class(cmdline),   intent(inout) :: cline
        integer :: i
        associate( sess => self%sess )
        ! Per-particle viewing AXIS, folded antipodally: +n and -n are mirror projections, so orientation
        ! bias lives in the axis and never in the mean resultant. Only the GMM's coverage term reads it.
        allocate(sess%pviews(3,sess%sel%nptcls))
        do i = 1, sess%sel%nptcls
            sess%pviews(:,i) = real(build%spproj_field%get_normal(sess%sel%pinds(i)), dp)
            if( sess%pviews(3,i) < 0.d0 ) sess%pviews(:,i) = -sess%pviews(:,i)
        end do
        ! the deconvolution's mixture components are the macro-clusters (full covariances,
        ! per-particle noise, held-out K); GMM AUTO's tied-covariance discovery only runs
        ! when no deconvolution happened
        if( params%min_state_frac > 0. )then
            ! population floor: exactly nstates hard-labelled states, every one at or above the floor
            ! (refine3D_states flex initialization relies on it); the bandwidth CV below does not apply
            sess%l_pop_floor = .true.
            call place_states_with_population_floor(sess%latent, sess%model, sess%nkern, sess%state_axis, sess%min_neff, &
                &params%min_state_frac, sess%states, equal_occ=trim(params%state_placement) == 'equal_occ')
        else
            call build_covariance_state_weights(sess%latent, sess%nkern, sess%state_axis, sess%min_neff, sess%states, &
                &macro_in=sess%deconv_labels, equal_occ=trim(params%state_placement) == 'equal_occ')
        endif

        ! Drop states below min_neff before reconstruction (they gave artefact maps).
        call prune_underpopulated_states(sess%min_neff, sess%states)

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

        if( params%nbins > 1 .and. .not. sess%l_pop_floor )then
            call cv_select_bandwidths(params, sess%cfg, build, sess%sel, params%nbins, sess%min_neff, sess%states, rounds=rounds)
        endif
        ! rec_states=no: deliver the placement, labels and readouts without reconstructing any
        ! map -- the states stage is then seconds instead of minutes
        sess%l_state_rec = trim(params%rec_states) /= 'no'
        if( .not. sess%l_state_rec )then
            write(logfhandle,'(A)') '>>> FLEX_PCA STATE RECONSTRUCTION OFF (rec_states=no): labels and &
                &tables only, no maps, no merge gate'
            call flush(logfhandle)
        endif
        end associate
    end subroutine infer_states

    !> The state maps (combined and per halfset), the two-gate merge and the re-reconstruction at the merged count
    subroutine reconstruct_states( self, params, build, rounds )
        class(flex_pca_application), intent(inout) :: self
        class(flex_pca_rounds), intent(inout) :: rounds
        type(parameters), intent(inout) :: params
        type(builder),    intent(inout) :: build
        integer :: i, r, s, nstates_merged
        integer, allocatable :: merge_label(:)
        real,    allocatable :: merged_weights(:,:), merged_targets(:,:), merged_bw(:)
        real(dp), allocatable :: state_mass(:), merged_mass(:)
        real(dp) :: sumw_s, sumw2_s
        integer(timer_int_kind) :: t_blk
        associate( sess => self%sess )
        if( sess%l_state_rec )then
        ! combined states and both halfsets in ONE pass; combined == even + odd exactly
        params%outvol = 'flex_pca_state_001.mrc'
        t_blk = tic()
        ! states of one macro-cluster share a delivery pool (per-shell adaptive kernel in rec3D)
        call reconstruct_flex_weighted_states(params, sess%cfg, build, sess%sel%pinds, sess%states%weights, sess%states%nstates, &
            &floor_rho=.true., outvol_even=string('flex_pca_even_state_001.mrc'), &
            &outvol_odd=string('flex_pca_odd_state_001.mrc'), rounds=rounds)
        write(logfhandle,'(A,F9.1)') '>>> FLEX_PCA STAGE states_combined_eo seconds=', toc(t_blk)
        ! Collapse indistinct states and reconstruct once at the surviving count. Needs the half maps above.
        if( sess%cfg%l_merge .and. sess%states%nstates > 1 )then
            t_blk = tic()
            allocate(merge_label(sess%states%nstates))
            call two_gate_state_merge(params, sess%cfg, sess%pviews, sess%states, merge_label, nstates_merged)
            if( nstates_merged < sess%states%nstates )then
                allocate(state_mass(sess%states%nstates))
                do s = 1, sess%states%nstates
                    state_mass(s) = sum(real(sess%states%weights(:,s), dp))
                end do
                allocate(merged_weights(sess%sel%nptcls,nstates_merged), source=0.)
                do s = 1, sess%states%nstates
                    merged_weights(:,merge_label(s)) = merged_weights(:,merge_label(s)) + sess%states%weights(:,s)
                end do
                call move_alloc(merged_weights, sess%states%weights)
                ! RECOMPUTE the label, do not remap it: the merge SUMS columns, and the pre-merge argmax
                ! can lose to a combined rival it beat individually (0.40 vs 0.35+0.25). Remapping would
                ! label the particle to a map it is no longer the largest contributor to, so the hard
                ! assignment and the delivered maps would describe different partitions. 0 is preserved:
                ! summing cannot lift a particle that was outside every kernel support.
                do i = 1, sess%sel%nptcls
                    if( sess%states%labels(i) >= 1 ) sess%states%labels(i) = maxloc(sess%states%weights(i,:), dim=1)
                end do
                ! collapse the per-state tables by the same mass, else they still describe the pre-merge nstates
                allocate(merged_targets(size(sess%states%targets,1),nstates_merged), source=0.)
                allocate(merged_bw(nstates_merged), source=0.)
                allocate(merged_mass(nstates_merged), source=0.d0)
                do s = 1, sess%states%nstates
                    r = merge_label(s)
                    merged_targets(:,r) = merged_targets(:,r) + real(state_mass(s))*sess%states%targets(:,s)
                    merged_bw(r)        = merged_bw(r)        + real(state_mass(s))*sess%states%bandwidths(s)
                    merged_mass(r)      = merged_mass(r)      + state_mass(s)
                end do
                do r = 1, nstates_merged
                    if( merged_mass(r) > DTINY )then
                        merged_targets(:,r) = merged_targets(:,r) / real(merged_mass(r))
                        merged_bw(r)        = merged_bw(r)        / real(merged_mass(r))
                    endif
                end do
                call move_alloc(merged_targets, sess%states%targets)
                call move_alloc(merged_bw,      sess%states%bandwidths)
                deallocate(sess%states%neff); allocate(sess%states%neff(nstates_merged))
                do r = 1, nstates_merged
                    sumw_s  = sum(real(sess%states%weights(:,r), dp))
                    sumw2_s = sum(real(sess%states%weights(:,r), dp)**2)
                    sess%states%neff(r) = real(sumw_s*sumw_s / max(sumw2_s, DTINY))
                end do
                deallocate(state_mass, merged_mass)
                ! downstream addresses the maps as a contiguous run, so the stale tail would deliver merged-away states
                do s = nstates_merged + 1, sess%states%nstates
                    call del_file('flex_pca_state_'     //int2str_pad(s,3)//MRC_EXT)
                    call del_file('flex_pca_even_state_'//int2str_pad(s,3)//MRC_EXT)
                    call del_file('flex_pca_odd_state_' //int2str_pad(s,3)//MRC_EXT)
                end do
                sess%states%nstates  = nstates_merged
                sess%l_merged = .true.
                call reconstruct_flex_weighted_states(params, sess%cfg, build, sess%sel%pinds, sess%states%weights, sess%states%nstates, &
                    &floor_rho=.true., outvol_even=string('flex_pca_even_state_001.mrc'), &
                    &outvol_odd=string('flex_pca_odd_state_001.mrc'), rounds=rounds)
            endif
            deallocate(merge_label)
            write(logfhandle,'(A,F9.1)') '>>> FLEX_PCA STAGE two_gate_merge seconds=', toc(t_blk)
        endif
        endif   ! l_state_rec
        end associate
    end subroutine reconstruct_states

    !> Tables, the weight store and hard labels into the project, the consensus nonuniform filter and the manifest
    subroutine publish( self, params, build, rounds )
        class(flex_pca_application), intent(inout) :: self
        class(flex_pca_rounds), intent(inout) :: rounds
        type(parameters), intent(inout) :: params
        type(builder),    intent(inout) :: build
        associate( sess => self%sess )
        call write_covariance_tables(sess%readout, build, sess%sel, sess%model, sess%latent, sess%states)
        ! The delivered weight table into the project-registered store, then the hard labels into
        ! the project itself, so the assignment can be judged by an INDEPENDENT reconstructor. Every
        ! non-worker delivers (shared memory, nparts=1 and the distributed master alike); a worker
        ! shares the master's projfile and must not write it. The store goes first: it validates
        ! against the field's activity as the run saw it, before the labels overwrite `state`.
        if( .not. rounds%is_worker() )then
            call write_flex_weights_store(params, build, sess%sel%pinds, sess%states%weights, sess%states%labels, sess%states%targets, sess%states%bandwidths, &
                &sess%l_merged)
            call write_discrete_state_project(build%spproj, sess%sel%pinds, sess%states%labels, sess%states%nstates, params%projfile)
        endif
        allocate(sess%states%half_weights(sess%sel%nptcls,sess%states%nstates), source=sess%states%weights)
        ! nonuniform filtering LAST, so every delivered map is filtered the same way
        if( sess%l_state_rec .and. trim(params%nufilt) == 'yes' ) call apply_consensus_nu_filter(params, sess%states%nstates)
        call write_covariance_manifest(params, sess%sel%nptcls, sess%model%ncomp, sess%states%nstates, sess%state_axis, sess%min_neff, sess%sigma_loaded)
        end associate
    end subroutine publish

    !> The distributed worker's one entry: the shared preparation (its particle list, poses, the
    !! pinned sigma decision), then exactly one stage body. Round state (stage, which_iter,
    !! maxits, nfits, pcafit) arrives in params from the master's job_descr. The caller (the
    !! worker strategy) signals completion with qsys_declare_part_finished.
    subroutine execute_worker_stage( self, params, build, cline, rounds )
        type(flex_fit_model)  :: model
        type(flex_selection)  :: sel
        class(flex_pca_application), intent(inout) :: self
        class(flex_pca_rounds), intent(inout) :: rounds
        type(parameters), intent(inout) :: params
        type(builder),    intent(inout) :: build
        class(cmdline),   intent(inout) :: cline
        real,     allocatable :: state_weights(:,:)
        integer :: nstates
        logical :: sigma_loaded, l_split_eo
        type(flex_run_settings) :: cfg
        call cfg%new(params, cline)
        call validate_covariance_inputs(params, cfg, build, sel%pinds, sel%nptcls, rounds=rounds)
        call load_and_validate_sigma(params, build, cline, sel%pinds, sigma_loaded, rounds=rounds)
        select case(params%stage)
            case(PCA_STAGE_STATES)
                ! the master produced the weight table; replicate its operator setup (ml_reg off,
                ! reconstruction Fourier band) and accumulate this part
                params%l_ml_reg = .false.
                params%ml_reg   = 'no'
                call cline%set('ml_reg','no')
                call build%esig%set_kfromto([1, max(1, fdim(flex_rec_box(params)) - 1)])
                call read_state_weights_round(sel%pinds, sel%nptcls, state_weights, nstates, l_split_eo)
                call reconstruct_flex_weighted_states(params, cfg, build, sel%pinds, state_weights, &
                    &nstates, floor_rho=.true., split_eo=l_split_eo, rounds=rounds)
                deallocate(state_weights)
            case(PCA_STAGE_PROBE, PCA_STAGE_POLISH, PCA_STAGE_EMBED)
                ! the mean's scale comes from the master's cache; the read-once env caches are set
                ! exactly as on the master path
                call estimate_covariance_mean(params, build, model%mean_rec, sel%pinds, sel%nptcls, rounds=rounds)
                if( params%stage == PCA_STAGE_EMBED )then
                    call embed_worker_pass(params, build, model, sel, rounds=rounds)
                else
                    call probe_worker_pass(params, cfg, build, model, sel, rounds=rounds)
                endif
                call model%mean_rec%dealloc_rho; call model%mean_rec%kill
            case default
                THROW_HARD('unknown flex_pca worker stage')
        end select
        deallocate(sel%pinds)
        call cfg%kill
    end subroutine execute_worker_stage

    !> Cost of the requested state count in resident reconstructors. REPORT ONLY.
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
        call flush(logfhandle)
    end subroutine report_state_memory

end module simple_flex_pca_application
