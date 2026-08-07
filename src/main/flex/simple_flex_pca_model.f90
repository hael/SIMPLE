!@descr: Standalone projection-aware low-rank covariance workflow for heterogeneous SPA data
module simple_flex_pca_model
use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
use simple_core_module_api
use simple_builder,                    only: builder
use simple_cmdline,                    only: cmdline
use simple_flex_pca_rec3D,         only: reconstruct_flex_weighted_states, flex_rec_box, flex_rec_smpd
use simple_flex_pca_columns,    only: cov_env_int_pub, build_covariance_eigenbasis, embed_latents_with_contrast, &
    &estimate_covariance_mean, probe_subspace_iteration, align_basis_to_reference, probe_external_basis, &
    &save_probe_state
use simple_image,                      only: image
use simple_parameters,                 only: parameters
use simple_reconstructor,              only: reconstructor
use simple_sigma2_files,               only: load_sigma2_groups
use simple_sp_project,                 only: sp_project
use simple_srch_sort_loc,              only: hpsort
use simple_flex_pca_distr,             only: flex_pca_is_worker, flex_pca_is_master, flex_pca_nparts, flex_pca_run_stage
use simple_flex_pca_parts,             only: write_sigma_state, check_sigma_state, &
    &read_state_weights_round, write_embed_part
use simple_flex_pca_distr,             only: PCA_STAGE_STATES, PCA_STAGE_EMBED
use simple_finch,                      only: finch_hierarchy, fit_finch, finch_representatives, &
    &                                          select_finch_level, refine_finch_level
use simple_kd_tree,                    only: kd_tree, knn_table
use simple_linalg,                     only: jacobi, eigsrt, matinv
use simple_flex_pca_merge,             only: flex_pca_merge_enabled, two_gate_state_merge
implicit none
private
#include "simple_local_flags.inc"

public :: run_flex_pca
public :: test_flex_pca_embedding_cache_io, test_flex_pca_kernel_bandwidth
public :: test_flex_pca_state_weights

character(len=8), parameter :: COV_CACHE_MAGIC   = 'SIMPLFXC'
! Bump whenever the cache layout changes; read_embedding_cache rejects any other version.
integer,          parameter :: COV_CACHE_VERSION = 3
! Safety cap on the bandwidth widening loop.
integer,          parameter :: COV_MAX_BW_GROW   = 4

contains

    subroutine run_flex_pca( params, build, cline )
        use simple_qsys_funs, only: qsys_job_finished
        type(parameters), intent(inout) :: params
        type(builder),    intent(inout) :: build
        class(cmdline),   intent(inout) :: cline
        type(reconstructor) :: mean_rec
        type(reconstructor), allocatable :: basis_recs(:)
        integer, allocatable :: pinds(:), labels(:)
        real(dp), allocatable :: z(:,:), eigvals(:), prior_precision(:), contrast(:)
        real(dp), allocatable :: latent_second(:,:,:)
        real(dp), allocatable :: resid_energy(:), resid_mean_energy(:)
        real, allocatable :: state_weights(:,:), half_weights(:,:), targets(:,:), bandwidths(:), neff(:)
        real, allocatable :: traj_weights(:,:), traj_targets(:,:), traj_bandwidths(:), traj_neff(:)
        real(dp), allocatable :: geo_targets(:,:)
        integer, allocatable :: traj_labels(:), comp_rank(:)
        integer :: nptcls, ncomp, nstates, min_neff, state_axis, ncols_req, col_sep, neigs_req, nkern
        integer :: q, i, r, s, metric_valid_count, axis_sel, nfinch
        integer :: nstates_merged
        integer, allocatable :: merge_label(:)
        real,    allocatable :: merged_weights(:,:)
        real(dp), allocatable :: finch_targets(:,:)
        logical :: l_finch_states
        logical :: l_geo, l_rot
        real(dp), allocatable :: kdist(:,:), kfloor(:), tdist(:,:), tfloor(:)
        real(dp), allocatable :: comp_rho(:)     ! per-component reliability, drives state-target ordering
        real(dp), allocatable :: pviews(:,:)     ! per-particle viewing AXIS, for the GMM orientation term
        character(len=:), allocatable :: cachedir, cachestr
        real(dp) :: sig2_eff
        logical :: sigma_loaded, l_resume, l_split_eo
        integer(timer_int_kind) :: t_blk

        call validate_covariance_inputs(params, build, cline, pinds, nptcls)
        call load_and_validate_sigma(params, build, cline, pinds, sigma_loaded)

        ! STATE-RECONSTRUCTION WORKER. Nothing upstream is needed: the master has already produced the
        ! weight table for this round. Replicate the operator setup the master applies before its own
        ! state calls (ml_reg off, the reconstruction Fourier band) and go straight to the accumulation.
        if( flex_pca_is_worker() .and. params%stage == PCA_STAGE_STATES )then
            params%l_ml_reg = .false.
            params%ml_reg   = 'no'
            call cline%set('ml_reg','no')
            call build%esig%set_kfromto([1, max(1, fdim(flex_rec_box(params)) - 1)])
            call read_state_weights_round(pinds, nptcls, state_weights, nstates, l_split_eo)
            call reconstruct_flex_weighted_states(params, build, pinds, state_weights, &
                &nstates, floor_rho=.true., split_eo=l_split_eo)
            call qsys_job_finished(params, string('simple_flex_pca_model :: run_flex_pca states'))
            deallocate(pinds, state_weights)
            return
        endif

        neigs_req  = max(1, min(48, params%neigs))
        ! No fixed ceiling on the state count: a compositionally rich specimen can carry far more
        ! states than any constant would allow, and an over-provisioned run is the only regime in
        ! which the two-gate merge can recover K. What actually bounds it is memory, which is
        ! REPORTED below rather than enforced.
        nstates    = max(3, params%npreimages)
        call report_state_memory(params, nstates)
        ! env-only: inert on the default path (the GMM replaces the kernel weights and bandwidth),
        ! live on the nbins>1 and SIMPLE_COV_GMM=0 opt-outs
        min_neff = params%min_neff
        call cov_env_int_pub('SIMPLE_COV_MIN_NEFF', min_neff)
        min_neff = max(20, min(nptcls, min_neff))
        state_axis = params%state_axis      ! <0 path, 0 k-means, >=1 legacy single axis
        ! nkern decouples the number of components the STATE STAGE uses from neigs, the number estimated.
        nkern      = params%nkern
        if( nkern <= 0 ) nkern = huge(1)/2   ! clamped against ncomp once the fit is known
        ncols_req  = max(4, params%ncols)
        col_sep    = max(1, params%column_separation)

        write(logfhandle,'(A,I0,A,I0,A,I0)') '>>> FLEX_PCA particles=',nptcls, &
            &' requested_columns=',ncols_req,' requested_components=',neigs_req
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

        ! RESUME MODE. The covariance basis and the per-particle embedding dominate the runtime and are
        ! COMPLETELY INDEPENDENT of the state count, the kernel bandwidth and the target placement, so
        ! infile= resumes from a cached embedding and re-runs only the state stages.
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

        ! consensus mean mu (eq. S.1)
        call estimate_covariance_mean(params, build, mean_rec, pinds, nptcls)

        if( params%l_heldout )then
            ! each halfset is embedded in the OTHER halfset's basis
            call heldout_embedding(params, build, mean_rec, pinds, nptcls, ncols_req, col_sep, neigs_req, &
                &basis_recs, eigvals, ncomp, sig2_eff, z, contrast, latent_second, resid_energy, resid_mean_energy)
            if( state_axis > 0 ) state_axis = min(state_axis, min(ncomp, nkern))
            allocate(prior_precision(ncomp))
            do q = 1, ncomp
                prior_precision(q) = 1.d0 / max(eigvals(q), DTINY)
            end do
            write(logfhandle,'(A,I0)') '>>> FLEX_PCA held-out retained covariance components=',ncomp
            call flush(logfhandle)
        else

            call build_covariance_eigenbasis(params, build, mean_rec, pinds, nptcls, &
                &ncols_req, col_sep, neigs_req, basis_recs, eigvals, ncomp, sig2_eff)
            if( flex_pca_is_worker() )then
                ! this worker's part file for the requested stage is on disk; the master reduces it
                call qsys_job_finished(params, string('simple_flex_pca_model :: run_flex_pca'))
                call mean_rec%dealloc_rho; call mean_rec%kill
                deallocate(pinds)
                return
            endif
            if( state_axis > 0 ) state_axis = min(state_axis, min(ncomp, nkern))
            write(logfhandle,'(A,I0)') '>>> FLEX_PCA retained covariance components=',ncomp
            call flush(logfhandle)

            ! Optional refinement: alternating the Wiener E-step (per-particle latents, with the anisotropic
            ! posterior second moment) with the weighted-backprojection M-step cleans the per-particle
            ! projection directions the covariance columns give, which carry a large noise fraction.
            if( params%n_probe_iters > 0 )then
                call probe_subspace_iteration(params, build, mean_rec, basis_recs, eigvals, sig2_eff, &
                    &pinds, nptcls, ncomp, params%n_probe_iters)
                if( state_axis > 0 ) state_axis = min(state_axis, min(ncomp, nkern))
                write(logfhandle,'(A,I0)') '>>> FLEX_PCA probe-refined components=',ncomp
                call flush(logfhandle)
            endif

            allocate(z(nptcls,ncomp), prior_precision(ncomp))
            allocate(latent_second(ncomp,ncomp,nptcls))
            allocate(resid_energy(nptcls), resid_mean_energy(nptcls))

            ! MAP prior precision = 1/Gamma_q
            do q = 1, ncomp
                prior_precision(q) = 1.d0 / max(eigvals(q), DTINY)
            end do
            write(logfhandle,'(A,ES12.4,A,ES12.4)') '>>> FLEX_PCA covariance eigenvalues: max=', &
                &maxval(eigvals),' min=',minval(eigvals)
            call flush(logfhandle)
            ! contrast-aware MAP embedding (S.D/S.E)
            allocate(contrast(nptcls))
            allocate(comp_rho(ncomp), source=1.d0)
            ! ---- embedding: distributed when the master has parts ----
            ! One qsys round, unlike the probe, because the basis is fixed here and the workers need
            ! no per-iteration refresh. What they cannot do is finish, since the reliability prior
            ! couples every particle, so they ship sufficient statistics and the master owns rho and
            ! the re-solve. See write_embed_stats_part.
            if( flex_pca_is_master() .and. flex_pca_nparts() > 1 )then
                call save_probe_state(ncomp, eigvals, sig2_eff)
                call flex_pca_run_stage(PCA_STAGE_EMBED, 'embedding')
                call embed_latents_with_contrast(params, build, mean_rec, basis_recs, ncomp, eigvals, &
                    &sig2_eff, pinds, nptcls, z, contrast, latent_second, resid_energy, &
                    &resid_mean_energy, rho_out=comp_rho, from_parts=.true.)
            else
                call embed_latents_with_contrast(params, build, mean_rec, basis_recs, ncomp, eigvals, &
                    &sig2_eff, pinds, nptcls, z, contrast, latent_second, resid_energy, &
                    &resid_mean_energy, rho_out=comp_rho)
            endif
        endif
        write(logfhandle,'(A,F7.3,A,F7.3)') '>>> FLEX_PCA per-particle contrast: mean=', &
            &real(sum(contrast)/real(nptcls,dp)),' sd=', &
            &real(sqrt(max(sum((contrast-sum(contrast)/nptcls)**2)/real(nptcls,dp),DTINY)))
        call flush(logfhandle)
        ! z is left in the physical units the MAP solve returns; the kernel metric does the weighting
        call write_covariance_eigenvolumes(basis_recs, eigvals, ncomp)
        call write_embedding_cache('flex_pca_embedding.bin', pinds, nptcls, ncomp, z, eigvals, &
            &contrast, resid_energy, resid_mean_energy, latent_second, sig2_eff)

        endif   ! .not. l_resume
        ! Resuming from a cached embedding skips the split-half solve, so no measured rho exists.
        ! Fall back to the spread-over-posterior-variance proxy, which is computable from the cache
        ! alone and is the same ranking the nkern documentation prescribes.
        if( .not. allocated(comp_rho) )then
            allocate(comp_rho(ncomp))
            call component_reliability_proxy(z, latent_second, nptcls, ncomp, comp_rho)
            write(logfhandle,'(A)') '>>> FLEX_PCA resumed embedding: component reliability from the &
                &spread/posterior-variance proxy (no cached split-half rho)'
        endif
        ! optional external-basis probe (SIMPLE_COV_PROBEBASIS / SIMPLE_COV_PROBEN / SIMPLE_COV_PROBEEIG)
        call run_external_basis_probe(params, build, mean_rec, pinds, nptcls, ncomp, eigvals, &
            &sig2_eff, l_resume)

        ! Optional rotation of the eigenbasis toward spatially coherent components, BEFORE any state
        ! placement so every downstream stage sees one consistent frame.
        if( params%pcrot > 0. )then
            if( l_resume .and. .not. file_exists('flex_pca_pc001.mrc') )then
                call rotate_basis_by_smoothness(params, ncomp, params%pcrot, z, nptcls, latent_second, l_rot, &
                    &eigdir=cachedir)
            else
                call rotate_basis_by_smoothness(params, ncomp, params%pcrot, z, nptcls, latent_second, l_rot)
            endif
            if( .not. l_rot ) write(logfhandle,'(A)') &
                &'>>> FLEX_PCA smoothness rotation requested but not applied; continuing in the PCA frame'
        endif
        ! external targets take precedence, then FINCH (SIMPLE_COV_FINCH_STATES=1), then k-means
        l_finch_states = .false.
        call read_external_targets(ncomp, nstates, finch_targets, l_finch_states)
        if( l_finch_states )then
            write(logfhandle,'(A,I0,A,I0,A)') '>>> FLEX_PCA using EXTERNAL latent targets: ', &
                &nstates,' states over ',ncomp,' components'
        else if( cov_env_flag_on('SIMPLE_COV_FINCH_STATES') )then
            call finch_state_targets(z, nptcls, ncomp, nstates, finch_targets, nfinch, l_finch_states)
            if( l_finch_states )then
                write(logfhandle,'(A,I0,A,I0)') '>>> FLEX_PCA state count set by FINCH: ',nfinch, &
                    &'  (npreimages was an upper bound of ',nstates
                nstates = nfinch
            else
                write(logfhandle,'(A)') '>>> FLEX_PCA FINCH state placement failed; falling back to k-means'
            endif
        endif
        ! Per-particle viewing AXIS, folded antipodally: +n and -n are mirror projections carrying the
        ! same information, so orientation bias lives in the axis and never in the mean resultant.
        ! Only the GMM's optional orientation-coverage term consumes this; it is cheap to always build.
        allocate(pviews(3,nptcls))
        do i = 1, nptcls
            pviews(:,i) = real(build%spproj_field%get_normal(pinds(i)), dp)
            if( pviews(3,i) < 0.d0 ) pviews(:,i) = -pviews(:,i)
        end do
        if( l_finch_states )then
            call build_covariance_state_weights(z, nptcls, ncomp, nkern, nstates, state_axis, min_neff, &
                &eigvals, latent_second, state_weights, targets, bandwidths, neff, labels, &
                &dist_out=kdist, bfloor_out=kfloor, targets_in=finch_targets, views=pviews)
            deallocate(finch_targets)
        else
            call build_covariance_state_weights(z, nptcls, ncomp, nkern, nstates, state_axis, min_neff, &
                &eigvals, latent_second, state_weights, targets, bandwidths, neff, labels, &
                &dist_out=kdist, bfloor_out=kfloor, comp_rho=comp_rho, views=pviews)
        endif
        ! ---- STATE-RECONSTRUCTION OPERATOR SETUP ---- MUST precede cv_select_bandwidths, which
        ! reconstructs trial half maps through the same backend.
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

        if( params%nbins > 1 )then
            call cv_select_bandwidths(params, build, pinds, nptcls, nstates, params%nbins, min_neff, &
                &kdist, kfloor, state_weights, bandwidths, neff)
        endif
        call write_covariance_tables(build, pinds, z, eigvals, prior_precision, state_weights, labels, &
            &targets, bandwidths, neff, resid_energy, resid_mean_energy)
        ! Hard state labels as a runnable project, so the embedding and its state assignment can be
        ! judged by an INDEPENDENT reconstructor (plain reconstruct3D) rather than only through the
        ! kernel-weighted backend below, which shares every upstream assumption with the embedding.
        call write_discrete_state_project(build%spproj, pinds, labels, nstates, params%outfile)

        ! combined states and both halfsets in ONE pass through the gridding reconstructor
        ! (see reconstruct_flex_weighted_states: combined == even + odd exactly)
        params%outvol = 'flex_pca_state_001.mrc'
        t_blk = tic()
        call reconstruct_flex_weighted_states(params, build, pinds, state_weights, nstates, &
            &floor_rho=.true., outvol_even=string('flex_pca_even_state_001.mrc'), &
            &outvol_odd=string('flex_pca_odd_state_001.mrc'))
        write(logfhandle,'(A,F9.1)') '>>> FLEX_PCA STAGE states_combined_eo seconds=', toc(t_blk)
        ! ---- two-gate merge ----
        ! Collapse states the orientations or the maps say are not distinct, then reconstruct once
        ! at the surviving count. Runs here because it needs the half maps the pass above just wrote,
        ! and before the trajectory ordering, which should order the states actually delivered.
        if( flex_pca_merge_enabled() .and. nstates > 1 )then
            t_blk = tic()
            allocate(merge_label(nstates))
            call two_gate_state_merge(params, pviews, state_weights, nptcls, nstates, &
                &merge_label, nstates_merged)
            if( nstates_merged < nstates )then
                allocate(merged_weights(nptcls,nstates_merged), source=0.)
                do s = 1, nstates
                    merged_weights(:,merge_label(s)) = merged_weights(:,merge_label(s)) + state_weights(:,s)
                end do
                call move_alloc(merged_weights, state_weights)
                do i = 1, nptcls
                    if( labels(i) >= 1 .and. labels(i) <= nstates ) labels(i) = merge_label(labels(i))
                end do
                ! the re-reconstruction below writes 001..nstates_merged, so the maps above that
                ! index are stale. Everything downstream addresses the state maps as the contiguous
                ! run flex_pca_state_001..NNN, so leaving the tail behind delivers states the merge
                ! just decided do not exist.
                do s = nstates_merged + 1, nstates
                    call del_file('flex_pca_state_'     //int2str_pad(s,3)//MRC_EXT)
                    call del_file('flex_pca_even_state_'//int2str_pad(s,3)//MRC_EXT)
                    call del_file('flex_pca_odd_state_' //int2str_pad(s,3)//MRC_EXT)
                end do
                nstates = nstates_merged
                call reconstruct_flex_weighted_states(params, build, pinds, state_weights, nstates, &
                    &floor_rho=.true., outvol_even=string('flex_pca_even_state_001.mrc'), &
                    &outvol_odd=string('flex_pca_odd_state_001.mrc'))
                call write_discrete_state_project(build%spproj, pinds, labels, nstates, params%outfile)
            endif
            deallocate(merge_label)
            write(logfhandle,'(A,F9.1)') '>>> FLEX_PCA STAGE two_gate_merge seconds=', toc(t_blk)
        endif
        allocate(half_weights(nptcls,nstates), source=state_weights)
        ! Order the reconstructed states along the dominant conformational motion in VOLUME space so the
        ! trajectory endpoints span the actual change; the covariance state axis need not order the
        ! states along the motion.
        axis_sel = 0
        allocate(comp_rank(ncomp))
        do q = 1, ncomp
            comp_rank(q) = q
        end do
        if( file_exists('flex_pca_pc001.mrc') )then
            t_blk = tic()
            call order_states_by_volume_trajectory(params, nstates, ncomp, axis_sel=axis_sel, comp_rank=comp_rank)
            write(logfhandle,'(A,F9.1)') '>>> FLEX_PCA STAGE traj_ordering seconds=', toc(t_blk)
        else if( l_resume .and. file_exists(cachedir//'flex_pca_pc001.mrc') )then
            ! resumed run: the eigenvolumes sit beside the cache we resumed from
            call order_states_by_volume_trajectory(params, nstates, ncomp, eigdir=cachedir, axis_sel=axis_sel, &
                &comp_rank=comp_rank)
        else
            write(logfhandle,'(A)') '>>> FLEX_PCA skipping volume-space trajectory ordering: no &
                &flex_pca_pc###.mrc in the working directory or beside the embedding cache'
        endif
        ! ---- AUTOMATIC CONFORMATIONAL TRAJECTORY ---- Re-place the targets ALONG the component just
        ! identified as carrying the reproducible motion, spanning its occupied range, and reconstruct.
        if( axis_sel > min(ncomp, nkern) )then
            write(logfhandle,'(A,I0,A,I0,A)') '>>> FLEX_PCA most reproducible component (',axis_sel, &
                &') lies outside the ',min(ncomp,nkern),' components retained by nkern; keeping the &
                &re-ordered trajectory. Raise nkern to place targets along it'
            axis_sel = 0
        endif
        if( axis_sel > 0 )then
            ! preferred placement is the FINCH manifold geodesic; fall back to the single axis
            l_geo = .false.
            if( params%ngeo /= 0 )then
                allocate(geo_targets(ncomp,nstates))
                call finch_geodesic_targets(z, nptcls, ncomp, nstates, comp_rank, &
                    &merge(params%ngeo, min(ncomp,6), params%ngeo > 0), axis_sel, geo_targets, l_geo)
                if( .not. l_geo ) deallocate(geo_targets)
            endif
            if( l_geo )then
                write(logfhandle,'(A)') '>>> FLEX_PCA placing trajectory targets on the FINCH manifold &
                    &geodesic -> flex_pca_traj_###.mrc'
                call build_covariance_state_weights(z, nptcls, ncomp, nkern, nstates, axis_sel, min_neff, &
                    &eigvals, latent_second, traj_weights, traj_targets, traj_bandwidths, traj_neff, traj_labels, &
                    &dist_out=tdist, bfloor_out=tfloor, targets_in=geo_targets)
                deallocate(geo_targets)
            else
                write(logfhandle,'(A,I0,A)') '>>> FLEX_PCA placing trajectory targets along component ', &
                    &axis_sel,' (selected by half-set reproducibility) -> flex_pca_traj_###.mrc'
                call build_covariance_state_weights(z, nptcls, ncomp, nkern, nstates, axis_sel, min_neff, &
                    &eigvals, latent_second, traj_weights, traj_targets, traj_bandwidths, traj_neff, traj_labels, &
                    &dist_out=tdist, bfloor_out=tfloor, comp_rho=comp_rho)
            endif
            ! The frames are a SEQUENCE -- flex_pca_traj_001..NNN is read as one -- but neither
            ! placement guarantees the emitted order runs along the path. Sort before anything
            ! downstream indexes them.
            call order_trajectory_targets(traj_weights, traj_targets, traj_bandwidths, traj_neff, &
                &traj_labels, min(ncomp,nkern), nstates)
            ! Cross-validate the bandwidth PER TRAJECTORY FRAME. Placement fixes each frame's target but
            ! not how much of the manifold it averages, and a density-blind floor gives dense middle frames
            ! several times the effective particle count of the sparse ends -- the smeared-middle failure.
            if( params%nbins > 1 )then
                call cv_select_bandwidths(params, build, pinds, nptcls, nstates, params%nbins, min_neff, &
                    &tdist, tfloor, traj_weights, traj_bandwidths, traj_neff)
            endif
            params%outvol = 'flex_pca_traj_001.mrc'
            t_blk = tic()
            call reconstruct_flex_weighted_states(params, build, pinds, traj_weights, nstates, floor_rho=.true.)
            write(logfhandle,'(A,F9.1)') '>>> FLEX_PCA STAGE states_traj seconds=', toc(t_blk)
            call write_trajectory_targets(axis_sel, nstates, ncomp, traj_targets, traj_bandwidths, traj_neff)
            deallocate(traj_weights, traj_targets, traj_bandwidths, traj_neff, traj_labels)
            if( allocated(tdist) ) deallocate(tdist, tfloor)
        endif
        ! nonuniform filtering LAST, so every delivered map is filtered the same way
        if( trim(params%nufilt) == 'yes' ) call apply_consensus_nu_filter(params, nstates)
        call write_covariance_manifest(params, nptcls, ncomp, nstates, state_axis, min_neff, sigma_loaded)

        ! in resume mode neither the mean reconstructor nor the basis reconstructors were
        ! ever built -- only the cached embedding was read
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
    end subroutine run_flex_pca

    subroutine validate_covariance_inputs( params, build, cline, pinds, nptcls )
        type(parameters), intent(inout) :: params
        type(builder),    intent(inout) :: build
        class(cmdline),   intent(in)    :: cline
        integer, allocatable, intent(out) :: pinds(:)
        integer, intent(out) :: nptcls
        integer :: q
        if( trim(params%oritype) /= 'ptcl3D' ) THROW_HARD('flex_pca requires oritype=ptcl3D')
        if( .not. cline%defined('vol1') ) THROW_HARD('flex_pca requires vol1=<consensus mean map>')
        if( trim(params%ptcl_src) /= 'raw' ) THROW_HARD('flex_pca currently requires ptcl_src=raw')
        if( params%nstates /= 1 ) THROW_HARD('flex_pca requires a one-state input project')
        call build%spproj_field%sample4rec([params%fromp,params%top],nptcls,pinds)
        if( nptcls < 100 ) THROW_HARD('flex_pca requires at least 100 active particles')
        if( build%spproj%os_ptcl2D%get_noris() > 0 .and. &
            &build%spproj%os_ptcl2D%get_noris() /= build%spproj%os_ptcl3D%get_noris() ) &
            &THROW_HARD('flex_pca requires matching ptcl2D and ptcl3D rows')
        ! assign gold-standard halfsets by index parity if the project's eo split is degenerate.
        ! MASTER-ONLY and must be persisted: a worker counts over its own range so it could reach a
        ! different verdict, and partition_eo mutates the field in memory only, so a worker re-reading
        ! the project would silently build a different -- equally valid -- split.
        if( flex_pca_is_worker() )then
            if( count([(build%spproj_field%get_eo(pinds(q))==0,q=1,nptcls)]) < 1 .or. &
                &count([(build%spproj_field%get_eo(pinds(q))==1,q=1,nptcls)]) < 1 )then
                THROW_HARD('flex_pca worker sees a degenerate eo split; the master did not persist its repair')
            endif
        else if( count([(build%spproj_field%get_eo(pinds(q))==0,q=1,nptcls)]) < 20 .or. &
            &count([(build%spproj_field%get_eo(pinds(q))==1,q=1,nptcls)]) < 20 )then
            write(logfhandle,'(A)') '>>> FLEX_PCA assigning alternating even/odd halfsets (project eo was degenerate)'
            call build%spproj_field%partition_eo
            if( flex_pca_is_master() )then
                call build%spproj%write_segment_inside(params%oritype, params%projfile)
                write(logfhandle,'(A)') '>>> FLEX_PCA persisted the repaired eo split for the workers'
            endif
        endif
        if( count([(build%spproj_field%get_eo(pinds(q))==0,q=1,nptcls)]) < 20 .or. &
            &count([(build%spproj_field%get_eo(pinds(q))==1,q=1,nptcls)]) < 20 ) &
            &THROW_HARD('flex_pca requires populated even and odd halfsets')
    end subroutine validate_covariance_inputs

    !> Report what the requested state count will cost in resident reconstructors.
    !!
    !! reconstruct_flex_weighted_states allocates all nstates reconstructors before the particle
    !! loop, plus a second nstates when the halfsets are fused, each carrying an expanded complex
    !! grid and its rho at the reconstruction box. The cost is therefore linear in nstates, cubic in
    !! box_rec, and paid up front -- so it can be estimated here and named alongside the knobs that
    !! control it.
    !!
    !! REPORT ONLY. This used to refuse to run above a fixed 64 GB budget, which is not a number this
    !! code can know: the figure below is an estimate of ONE allocation on a machine whose memory,
    !! other tenants and swap are all unknown here, and the nparts multiplier assumes every process is
    !! resident at once. Refusing on it blocked runs that fit and never caught the cases that did not.
    !! The operator sees the number and decides.
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
        ! every process pays this, not just the master: a worker owns a particle subset but all of
        ! the states, so the machine-wide peak is nparts+1 times the per-process figure
        nproc   = max(1, params%nparts) + 1
        gb      = gb_proc * real(nproc,dp)
        write(logfhandle,'(A,I0,A,I0,A,F8.2,A,I0,A,F8.2,A)') '>>> FLEX_PCA states=',nstates, &
            &' at reconstruction box ',box_rec,' -> approx ',gb_proc,' GB of reconstructors per &
            &process x ',nproc,' processes = ',gb,' GB machine-wide'
        write(logfhandle,'(A)') '>>> FLEX_PCA the knobs that move it are npreimages (linear) and &
            &box_crop (cubic)'
        ! the reconstructors are rarely what hurts: the reduced solve's accumulator is sized against
        ! COV_ATHR_BUDGET, which is likewise per process and predates distributed execution, so a
        ! distributed run multiplies it by nproc too. SIMPLE_COV_DTILDE is the knob that moves it.
        if( params%nparts > 1 ) write(logfhandle,'(A,I0,A)') '>>> FLEX_PCA NOTE: the reduced-solve &
            &accumulator is also per process; at nparts=',params%nparts,' it is paid that many times &
            &over. Cap it with SIMPLE_COV_DTILDE if the machine is tight.'
        call flush(logfhandle)
    end subroutine report_state_memory

    subroutine load_and_validate_sigma( params, build, cline, pinds, loaded )
        type(parameters), intent(inout) :: params
        type(builder),    intent(inout) :: build
        class(cmdline),   intent(inout) :: cline
        integer,          intent(in)    :: pinds(:)
        logical,          intent(out)   :: loaded
        integer :: i, k, iptcl, noris, kto
        call load_sigma2_groups(params, build%pftc, build%esig, build%spproj_field, cline, loaded)
        if( .not. loaded )then
            ! gen_fplane4rec is preceded by norm_noise_taper_edge_pad_fft, so a unit spectrum is the correct
            ! fallback for white synthetic noise after background normalization.
            noris = build%spproj_field%get_noris()
            kto   = max(1, fdim(params%box)-1)
            if( allocated(build%esig%sigma2_noise) ) deallocate(build%esig%sigma2_noise)
            allocate(build%esig%sigma2_noise(1:kto,1:noris), source=1.0)
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
        if( flex_pca_is_master() )then
            call write_sigma_state(loaded)
        else if( flex_pca_is_worker() )then
            call check_sigma_state(loaded)
        endif
    end subroutine load_and_validate_sigma

    !> Median of the chi-squared distribution with k degrees of freedom, via the Wilson-Hilferty cube-root
    !! normal approximation at the median (z=0): chi2_med(k) ~= k * (1 - 2/(9k))^3. Good to 0.02 % at
    !! k=20 and 3 % at k=1 -- far inside the tolerance of a bandwidth FLOOR, and it avoids pulling in an
    !! incomplete-gamma inverse.
    pure real(dp) function chi2_median( k )
        integer, intent(in) :: k
        real(dp) :: kk
        kk = real(max(k,1),dp)
        chi2_median = kk * (1.d0 - 2.d0/(9.d0*kk))**3
    end function chi2_median

    !> Binary cache of everything the state-weight and reconstruction stages need, so a different state
    !! count / bandwidth / placement can be tried without re-fitting the covariance basis and re-embedding
    !! every particle.
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

    !> Counterpart of write_embedding_cache.
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
        integer :: u, ver, nptcls_c, i
        if( .not. file_exists(fname) ) THROW_HARD('flex_pca embedding cache not found: '//trim(fname))
        open(newunit=u, file=fname, status='old', action='read', access='stream', form='unformatted')
        read(u) magic, ver
        if( magic /= COV_CACHE_MAGIC ) THROW_HARD('not a flex_pca embedding cache: '//trim(fname))
        if( ver /= COV_CACHE_VERSION ) THROW_HARD('flex_pca embedding cache version mismatch; re-run the fit')
        read(u) nptcls_c, ncomp
        if( nptcls_c /= nptcls ) THROW_HARD('flex_pca embedding cache particle count does not match the project')
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

    subroutine write_covariance_eigenvolumes( basis_recs, eigvals, ncomp )
        type(reconstructor), intent(inout) :: basis_recs(ncomp)
        real(dp),            intent(in)    :: eigvals(ncomp)
        integer,             intent(in)    :: ncomp
        integer :: q, u
        ! the eigenvolume MRCs are written by form_eigenbasis_from_reduced; only the table is written here
        call del_file('flex_pca_eigenvalues.txt')
        open(newunit=u,file='flex_pca_eigenvalues.txt',status='replace',action='write')
        write(u,'(A)') '# component eigenvalue'
        do q = 1, ncomp
            write(u,'(I6,1X,ES20.10)') q,eigvals(q)
        end do
        close(u)
    end subroutine write_covariance_eigenvolumes

    !> Kernel-regression reconstruction weights (supplement S.F): place nstates latent targets, then give
    !! every particle an Epanechnikov weight per state from its Mahalanobis distance to that target in the
    !! placement metric. Each state's bandwidth is floored so that at least min_neff particles fall
    !! inside its support.
    subroutine build_covariance_state_weights( z, nptcls, ncomp, nkern, nstates, axis, min_neff, &
        &eigvals, precision, weights, targets, bandwidths, neff, labels, dist_out, bfloor_out, targets_in, &
        &zmetric, comp_rho, views )
        integer,  intent(in) :: nptcls, ncomp, nkern, nstates, axis, min_neff
        real(dp), intent(in) :: z(nptcls,ncomp), eigvals(ncomp)
        real(dp), intent(in) :: precision(ncomp,ncomp,nptcls)   ! per-particle latent precision Pi_i
        real, allocatable, intent(out) :: weights(:,:), bandwidths(:), neff(:)
        real, allocatable, intent(out) :: targets(:,:)          ! (ncomp,nstates) latent target coordinates
        integer, allocatable, intent(out) :: labels(:)
        ! per-state kernel distances and bandwidth floors, so cv_select_bandwidths can rebuild
        ! weights at any bandwidth without redoing the nk^2 quadratic forms
        real(dp), allocatable, optional, intent(out) :: dist_out(:,:), bfloor_out(:)
        ! externally placed latent targets (ncomp,nstates).
        real(dp), optional,    intent(in) :: targets_in(ncomp,nstates)
        ! optional LOT pullback metric on the leading nk latent components; absent = identity
        real(dp), optional,    intent(in) :: zmetric(:,:)
        ! per-component reliability; enables the default reliability-ordered equal-occupancy placement
        real(dp), optional,    intent(in) :: comp_rho(ncomp)
        ! per-particle viewing AXIS; only needed for the GMM's orientation-coverage term
        real(dp), optional,    intent(in) :: views(3,nptcls)
        real,     allocatable :: sorted(:)
        real(dp), allocatable :: wcomp(:), tvec(:), tcen(:,:), dist(:), dvec(:), mvec(:)
        real(dp), allocatable :: pk(:,:,:), cfull(:,:), cblk(:,:), edges(:)
        real(dp), allocatable :: ppath(:), tpath(:)   ! per-particle / per-target coordinate along the path
        integer,  allocatable :: occ(:)
        real(dp) :: qlo, qhi, target, h, d2, u2, sumw, sumw2, best, zspread, bmin, chi2med, wsum_i
        real(dp) :: orient_lam
        integer  :: nrenorm, ispace
        integer  :: i, q, r, state, ilo, ihi, best_state, grow, nfed, occmax, ifloor, nunassigned, nsupp
        integer  :: nk, errflg
        logical  :: l_relpath, l_diffuse, l_gmm
        character(len=12) :: bwsrc
        ! Per-component weights for TARGET PLACEMENT (kmeans_latent_targets / path_latent_targets)
        nk = max(1, min(ncomp, nkern))
        allocate(wcomp(nk), tvec(nk), tcen(nk,nstates), dist(nptcls), dvec(nk), mvec(nk))
        wcomp = 1.d0
        ! STANDARDIZED PLACEMENT, on by default (opt out with SIMPLE_COV_STDZ=0). Weighting the placement
        ! metric by the eigenvalues concentrates every target along the highest-variance components, which
        ! are not the conformational ones; 1/var per component gives each retained direction equal say.
        if( .not. cov_env_flag_off('SIMPLE_COV_STDZ') )then
            do q = 1, nk
                zspread = sum(z(:,q)) / real(nptcls,dp)
                d2      = sum((z(:,q) - zspread)**2) / real(nptcls,dp)
                wcomp(q) = 1.d0 / max(d2, DTINY)
            end do
            wcomp = wcomp * real(nk,dp) / sum(wcomp)   ! keep the metric's overall scale
            write(logfhandle,'(A)') '>>> FLEX_PCA state placement metric: standardized (1/var per component)'
        endif
        ! degeneracy guard, in the same metric the targets are placed in
        zspread = 0.d0
        do q = 1, nk
            zspread = zspread + wcomp(q)*(maxval(z(:,q)) - minval(z(:,q)))**2
        end do
        if( sqrt(zspread) <= sqrt(DTINY) ) &
            &THROW_HARD('flex_pca latent embedding has zero spread; embedding collapsed')
        ! Restricting to the leading nk components MARGINALISES the precision (invert, slice, re-invert);
        ! slicing the precision matrix directly would condition on the dropped components instead.
        allocate(pk(nk,nk,nptcls))
        if( nk == ncomp )then
            pk = precision
        else
            !$omp parallel do default(shared) private(i,cfull,cblk,errflg) schedule(static)
            do i = 1, nptcls
                allocate(cfull(ncomp,ncomp), cblk(nk,nk))
                call matinv(precision(:,:,i), cfull, ncomp, errflg)
                cblk = cfull(1:nk,1:nk)
                call matinv(cblk, pk(:,:,i), nk, errflg)
                deallocate(cfull, cblk)
            end do
            !$omp end parallel do
            write(logfhandle,'(A,I0,A,I0,A)') '>>> FLEX_PCA state stage restricted to the leading ',nk, &
                &' of ',ncomp,' latent components (marginalised precision)'
        endif
        ! The reliability-ordered path is the default whenever a reliability vector reached us; opt out
        ! with SIMPLE_COV_KMEANS=1 to recover the historical k-means placement. Set only inside the
        ! branch that actually places the path, since it also selects the along-path weighting below.
        l_relpath = .false.
        allocate(weights(nptcls,nstates), targets(ncomp,nstates), bandwidths(nstates), neff(nstates), labels(nptcls))
        if( present(dist_out)   ) allocate(dist_out(nptcls,nstates))
        if( present(bfloor_out) ) allocate(bfloor_out(nstates))
        if( present(targets_in) )then
            ! ---- externally placed targets (FINCH geodesic) ----
            do state = 1, nstates
                tcen(:,state) = targets_in(1:nk,state)
            end do
            write(logfhandle,'(A,I0,A,I0)') '>>> FLEX_PCA state targets: supplied externally over ',nk, &
                &' components, points=',nstates
        else if( axis < 0 )then
            call path_latent_targets(z(:,1:nk), nptcls, nk, nstates, wcomp, tcen)
            write(logfhandle,'(A,I0,A,I0)') '>>> FLEX_PCA state targets: density-spread path over ',nk, &
                &' components, points=',nstates
        else if( axis == 0 .and. present(comp_rho) .and. .not. cov_env_flag_on('SIMPLE_COV_KMEANS') &
                &.and. .not. cov_env_flag_on('SIMPLE_COV_RELPATH') )then
            ! ---- DEFAULT: manifold-covering diffusion k-center ----
            ! Handles a continuous reaction coordinate and a branched set of compositional states
            ! with the same code and the same constants, because a curve is a degenerate graph.
            call diffusion_kcenter_targets(z(:,1:nk), nptcls, nk, nstates, wcomp, comp_rho(1:nk), &
                &tcen, l_diffuse)
            if( l_diffuse )then
                write(logfhandle,'(A,I0,A,I0)') '>>> FLEX_PCA state targets: diffusion k-center over ', &
                    &nk,' components, points=',nstates
            else
                write(logfhandle,'(A)') '>>> FLEX_PCA diffusion k-center unavailable; falling back &
                    &to the reliability-ordered path'
                allocate(ppath(nptcls), tpath(nstates))
                call reliability_path_targets(z(:,1:nk), nptcls, nk, nstates, wcomp, comp_rho(1:nk), tcen, &
                    &proj_out=ppath, tproj_out=tpath)
                l_relpath = .true.
            endif
        else if( axis == 0 .and. present(comp_rho) .and. .not. cov_env_flag_on('SIMPLE_COV_KMEANS') )then
            ! ---- 1-D reliability-ordered equal-occupancy path (SIMPLE_COV_RELPATH=1) ----
            ! Correct coverage on a genuine reaction coordinate, but it MERGES states on a branched
            ! manifold: on IgG-RL it collapsed 20 requested states onto 6 distinct ground-truth
            ! conformations. Kept for 1-D data and as the diffusion fallback, not as the default.
            allocate(ppath(nptcls), tpath(nstates))
            call reliability_path_targets(z(:,1:nk), nptcls, nk, nstates, wcomp, comp_rho(1:nk), tcen, &
                &proj_out=ppath, tproj_out=tpath)
            l_relpath = .true.
            write(logfhandle,'(A,I0,A,I0)') '>>> FLEX_PCA state targets: reliability-ordered &
                &equal-occupancy path over ',nk,' components, points=',nstates
        else if( axis == 0 )then
            ! ---- k-means centroids over the retained latent subspace ----
            call kmeans_latent_targets(z(:,1:nk), nptcls, nk, nstates, wcomp, tcen)
            write(logfhandle,'(A,I0,A,I0)') '>>> FLEX_PCA state targets: k-means over ',nk, &
                &' components, k=',nstates
        else
            ! ---- equal-occupancy slices along one component: every state gets the same particle count,
            ! and each target is the slice MEAN over all nk components, not a point on the axis ----
            if( axis > nk ) THROW_HARD('flex_pca state_axis exceeds the retained component count nkern')
            allocate(sorted(nptcls), source=real(z(:,axis)))
            call hpsort(sorted)
            allocate(edges(max(1,nstates-1)), source=0.d0)
            allocate(occ(nstates), source=0)
            do state = 1, nstates-1
                ifloor = max(1, min(nptcls, nint(real(state,dp)/real(nstates,dp)*real(nptcls,dp))))
                edges(state) = real(sorted(ifloor),dp)
            end do
            if( real(sorted(nptcls),dp)-real(sorted(1),dp) <= sqrt(DTINY) ) &
                &THROW_HARD('flex_pca state axis has zero range; embedding collapsed')
            tcen = 0.d0
            do i = 1, nptcls
                best_state = nstates
                do state = 1, nstates-1
                    if( z(i,axis) < edges(state) )then
                        best_state = state
                        exit
                    endif
                end do
                occ(best_state) = occ(best_state) + 1
                do q = 1, nk
                    tcen(q,best_state) = tcen(q,best_state) + z(i,q)
                end do
            end do
            do state = 1, nstates
                if( occ(state) > 0 ) tcen(:,state) = tcen(:,state) / real(occ(state),dp)
            end do
            write(logfhandle,'(A,I0,A,I0,A,I0,A,I0)') '>>> FLEX_PCA state targets: equal-occupancy slices along z', &
                &axis,' over ',nk,' components, points=',nstates,'  min slice occupancy=',minval(occ)
            deallocate(sorted, edges, occ)
        endif
        chi2med = chi2_median(nk)
        write(logfhandle,'(A,ES12.4)') &
            &'>>> FLEX_PCA Epanechnikov kernel regression, posterior-precision distance; chi2 median=', &
            &chi2med
        call flush(logfhandle)
        allocate(sorted(nptcls))
        do state = 1, nstates
            tvec = tcen(:,state)
            if( l_relpath )then
                ! ALONG-PATH distance. States placed on a path are weighted on that path, so the 15
                ! off-path directions -- mostly noise -- cannot push an on-path particle out of support.
                !$omp parallel do default(shared) private(i) schedule(static)
                do i = 1, nptcls
                    dist(i) = (ppath(i) - tpath(state))**2
                end do
                !$omp end parallel do
            else
            !$omp parallel do default(shared) private(i,q,r,d2,dvec,mvec) schedule(static)
            do i = 1, nptcls
                do q = 1, nk
                    dvec(q) = z(i,q) - tvec(q)
                end do
                ! optional LOT pullback metric; absent = identity
                if( present(zmetric) )then
                    do q = 1, nk
                        mvec(q) = 0.d0
                        do r = 1, nk
                            mvec(q) = mvec(q) + zmetric(q,r)*dvec(r)
                        end do
                    end do
                    dvec = mvec
                endif
                d2 = 0.d0
                do q = 1, nk
                    do r = 1, nk
                        d2 = d2 + dvec(q)*pk(q,r,i)*dvec(r)
                    end do
                end do
                dist(i) = max(d2, 0.d0)
            end do
            !$omp end parallel do
            endif
            sorted = real(dist)
            call hpsort(sorted)
            ifloor = max(1, min(nptcls, min_neff))
            if( l_relpath )then
                ! Tie the kernel to the TARGET SPACING rather than to the latent dimension. chi2(nk)
                ! grows with nk, so the dimension alone forces kernels wide enough to swallow several
                ! neighbouring targets -- which is why lowering nkern strands particles instead of
                ! sharpening states. One slice's worth of particles is the width at which the nstates
                ! kernels tile the path with modest overlap.
                ! Support spans TWO slices: one slice each side of the target. Adjacent supports then
                ! overlap, so every particle is inside at least one and none is stranded, while a state
                ! still draws predominantly on its own slice. The kernel support is dist < h^2 = 2*bmin,
                ! hence the half.
                ispace = max(1, min(nptcls, (2*nptcls)/max(nstates,1)))
                bmin   = 0.5d0*real(sorted(max(ifloor, ispace)),dp)
                bwsrc  = merge('path-spacing', 'min_neff-nn ', ispace >= ifloor)
            else
                bmin   = max(real(sorted(ifloor),dp), chi2med)
                bwsrc  = merge('chi2(ncomp) ', 'min_neff-nn ', chi2med >= real(sorted(ifloor),dp))
            endif
            h      = sqrt(2.d0*bmin)          ! their kernel arg is sqrt(d^2/(2b)) => h^2 = 2b
            if( present(dist_out)   ) dist_out(:,state) = dist
            if( present(bfloor_out) ) bfloor_out(state) = bmin
            ! log the distance quantiles: the chi2 floor is only meaningful if the posterior quadratic
            ! form really is on a chi2(nk) scale
            write(logfhandle,'(A,I3,A,ES11.3,A,ES11.3,A,ES11.3,A,A)') '>>>   state=',state, &
                &' dist: median=',real(sorted(max(1,nptcls/2)),dp),' p95=', &
                &real(sorted(max(1,nint(0.95*real(nptcls)))),dp),'  nn_floor=',real(sorted(ifloor),dp), &
                &'  bandwidth floor from ', bwsrc
            ! Widening loop. In nk dimensions the enclosed population grows like h^nk, so a 1.3x step is
            ! already a ~190x population step at nk=20 and every state ends up swallowing most of the
            ! dataset. The floor above should make this a no-op; COV_MAX_BW_GROW caps it if it does fire.
            nsupp = 0
            do grow = 0, COV_MAX_BW_GROW
                sumw  = 0.d0
                sumw2 = 0.d0
                nsupp = 0
                !$omp parallel do default(shared) private(i,u2) schedule(static) &
                !$omp& reduction(+:sumw,sumw2,nsupp)
                do i = 1, nptcls
                    u2 = dist(i) / (h*h)
                    weights(i,state) = real(max(0.d0, 1.d0 - u2))   ! Epanechnikov, compact support
                    sumw  = sumw  + real(weights(i,state),dp)
                    sumw2 = sumw2 + real(weights(i,state),dp)**2
                    if( weights(i,state) > 0. ) nsupp = nsupp + 1
                end do
                !$omp end parallel do
                if( nsupp >= min(min_neff, nptcls) ) exit
                if( grow >= COV_MAX_BW_GROW      ) exit
                h = 1.3d0*h                       ! safety only; should not fire
            end do
            if( nsupp < min(min_neff, nptcls) )then
                write(logfhandle,'(A,I3,A,I0,A,I0)') '>>>   WARNING state=',state, &
                    &' raw kernel support ',nsupp,' below min_neff after safety growth; requested ',min_neff
            endif
            if( maxval(weights(:,state)) > TINY ) weights(:,state)=weights(:,state)/maxval(weights(:,state))
            sumw  = sum(real(weights(:,state),dp))
            sumw2 = sum(real(weights(:,state),dp)**2)
            ! targets is reported over the FULL component set for the manifest
            targets(1:nk,state) = real(tcen(:,state))
            do q = nk+1, ncomp
                targets(q,state) = real(sum(z(:,q)) / real(nptcls,dp))
            end do
            bandwidths(state) = real(h)
            neff(state)       = real(sumw*sumw/max(sumw2,DTINY))
        end do
        ! ---- ASSIGNMENT RULE ----
        ! The tied-covariance mixture is the default; SIMPLE_COV_KERNEL=0 is not a thing, opt out with
        ! SIMPLE_COV_GMM=0 to recover the Epanechnikov kernel measured above. The kernel loop above
        ! still runs because dist_out feeds cv_select_bandwidths and the EXCLUSIVE path, and its
        ! distance quantiles remain the useful diagnostic of whether the chi2 scale is right.
        l_gmm = .not. cov_env_flag_off('SIMPLE_COV_GMM')
        if( l_gmm )then
            orient_lam = 0.d0
            call cov_env_dp('SIMPLE_COV_ORIENT_LAM', orient_lam)
            if( present(views) .and. orient_lam > 0.d0 )then
                call gmm_state_weights(z, nptcls, ncomp, nk, nstates, tcen, wcomp, weights, neff, &
                    &bandwidths, labels, views=views, orient_lam=orient_lam)
            else
                call gmm_state_weights(z, nptcls, ncomp, nk, nstates, tcen, wcomp, weights, neff, &
                    &bandwidths, labels)
            endif
            ! tcen now holds the FITTED means; refresh the reported target table so the manifest and
            ! flex_pca_state_targets.txt describe the maps that were delivered
            do state = 1, nstates
                targets(1:nk,state) = real(tcen(:,state))
            end do
        endif
        ! Nearest-state label = argmax weight, with 0 for a particle outside EVERY kernel support.
        ! Defaulting to state 1 instead would pile every unassigned particle onto the first state and make
        ! the occupancy report below read as a concentration failure when it is nothing of the kind.
        ! Under the mixture every responsibility is strictly positive, so this always assigns.
        nunassigned = 0
        do i = 1, nptcls
            best_state = 0
            best = 0.d0
            do state = 1, nstates
                if( real(weights(i,state),dp) > best )then
                    best = real(weights(i,state),dp)
                    best_state = state
                endif
            end do
            labels(i) = best_state
            if( best_state == 0 ) nunassigned = nunassigned + 1
        end do
        ! NORMALISED WEIGHTS (SIMPLE_COV_NORMW=1) -- rescale each particle's weights to a partition of
        ! unity over the states, so a particle in many kernels does not contribute to all of them at
        ! full strength.
        if( cov_env_flag_on('SIMPLE_COV_NORMW') .and. .not. l_gmm )then
            nrenorm = 0
            do i = 1, nptcls
                wsum_i = sum(real(weights(i,:), dp))
                if( wsum_i <= DTINY ) cycle
                weights(i,:) = real(real(weights(i,:), dp) / wsum_i)
                nrenorm = nrenorm + 1
            end do
            do state = 1, nstates
                sumw  = sum(real(weights(:,state), dp))
                sumw2 = sum(real(weights(:,state), dp)**2)
                neff(state) = real(sumw*sumw / max(sumw2, DTINY))
            end do
            write(logfhandle,'(A,I0,A)') '>>> FLEX_PCA NORMALISED weights: ', nrenorm, &
                &' particles renormalised to unit total weight across states'
            call flush(logfhandle)
        endif
        ! EXCLUSIVE ASSIGNMENT (SIMPLE_COV_EXCLUSIVE=1) -- hard-assign each particle to its nearest target
        ! and give under-populated states their own min_neff nearest particles, so no two state maps share
        ! any particle.
        if( cov_env_flag_on('SIMPLE_COV_EXCLUSIVE') .and. present(dist_out) .and. .not. l_gmm )then
            weights = 0.
            do i = 1, nptcls
                best_state = 1
                d2 = dist_out(i,1)
                do state = 2, nstates
                    if( dist_out(i,state) < d2 )then
                        d2 = dist_out(i,state); best_state = state
                    endif
                end do
                labels(i)             = best_state
                weights(i,best_state) = 1.
            end do
            do state = 1, nstates
                if( count(weights(:,state) > 0.) >= min(min_neff,nptcls) ) cycle
                sorted = real(dist_out(:,state))
                call hpsort(sorted)
                d2 = real(sorted(max(1,min(nptcls,min_neff))),dp)
                do i = 1, nptcls
                    if( dist_out(i,state) <= d2 ) weights(i,state) = 1.
                end do
                write(logfhandle,'(A,I3,A)') '>>>   state ',state,' under-populated; claimed its own nearest'
            end do
            do state = 1, nstates
                neff(state)       = real(count(weights(:,state) > 0.))
                bandwidths(state) = 0.
            end do
            nunassigned = count(labels < 1)
            write(logfhandle,'(A)') '>>> FLEX_PCA EXCLUSIVE assignment: states disjoint by construction'
            call flush(logfhandle)
        endif
        ! occupancy report
        allocate(occ(nstates), source=0)
        do i = 1, nptcls
            if( labels(i) >= 1 ) occ(labels(i)) = occ(labels(i)) + 1
        end do
        occmax = maxval(occ)
        nfed   = count(occ > 100)
        write(logfhandle,'(A,F6.2,A,I0,A,I0)') '>>> FLEX_PCA state occupancy (of assigned): max=', &
            &100.0*real(occmax)/real(max(nptcls-nunassigned,1)),'%  states with >100 particles=',nfed,' of ',nstates
        write(logfhandle,'(A,I0,A,F6.2,A)') '>>> FLEX_PCA outside every kernel support: ',nunassigned, &
            &' particles (',100.0*real(nunassigned)/real(max(nptcls,1)),'%) -- these contribute to no state map'
        do state = 1, nstates
            write(logfhandle,'(A,I3,A,I9,A,F10.1,A,ES11.3,A,F7.3)') '>>>   state=',state,'  particles=',occ(state), &
                &'  neff=',neff(state),'  bandwidth=',bandwidths(state),'  frac_in_support=', &
                &real(count(weights(:,state) > 0.))/real(max(nptcls,1))
        end do
        call flush(logfhandle)
        deallocate(wcomp, tvec, tcen, occ, dist, dvec, sorted, pk)
    end subroutine build_covariance_state_weights

    !>  Order a trajectory's targets ALONG the path they sit on, and permute everything that indexes
    !!  them to match.
    !!
    !!  The geodesic and single-axis placements both return a SET of targets on the manifold; nothing
    !!  guarantees the order they come out in is the order along the path, and the file names
    !!  flex_pca_traj_001..NNN are read by everyone as a sequence. Measured on Ribosembly: the sixteen
    !!  delivered frames matched ground-truth states 2, 8, 12, 12, 12, 10, 7, 9, 3, 11, 15, 14, 5, 4,
    !!  13, 15 in file order -- the volumes spanned the assembly series correctly (mass tracked the
    !!  matched state at r=0.999) but flipping through them showed no progression, which reads as a
    !!  far worse result than it is.
    !!
    !!  Ordering is by projection onto the leading principal direction of the TARGET SET, which is the
    !!  arclength direction for anything path-like. Sign is fixed by the projection of the first
    !!  target so repeat runs do not silently reverse. A closed loop has no well-defined order and is
    !!  left alone -- detected as a leading direction carrying less than half the target variance.
    subroutine order_trajectory_targets( weights, targets, bandwidths, neff, labels, nk, nstates )
        real,    allocatable, intent(inout) :: weights(:,:), targets(:,:), bandwidths(:), neff(:)
        integer, allocatable, intent(inout) :: labels(:)
        integer,              intent(in)    :: nk, nstates
        real(dp), allocatable :: T(:,:), C(:,:), ev(:), evec(:,:), proj(:)
        real,     allocatable :: wtmp(:,:), ttmp(:,:), btmp(:), ntmp(:)
        integer,  allocatable :: ord(:), inv(:)
        real(dp) :: mu_q, frac
        integer  :: q, r, s, nrot, i
        if( nstates < 3 ) return
        allocate(T(nk,nstates), C(nk,nk), ev(nk), evec(nk,nk), proj(nstates), ord(nstates))
        do s = 1, nstates
            do q = 1, nk
                T(q,s) = real(targets(q,s),dp)
            end do
        end do
        do q = 1, nk
            mu_q = sum(T(q,:)) / real(nstates,dp)
            T(q,:) = T(q,:) - mu_q
        end do
        C = matmul(T, transpose(T)) / real(nstates,dp)
        call jacobi(C, nk, nk, ev, evec, nrot)
        call eigsrt(ev, evec, nk, nk)
        frac = ev(1) / max(sum(ev), DTINY)
        if( frac < 0.5d0 )then
            write(logfhandle,'(A,F6.3,A)') '>>> FLEX_PCA trajectory targets are not path-like (leading &
                &direction carries ',frac,' of target variance); leaving the emitted order alone'
            deallocate(T, C, ev, evec, proj, ord)
            return
        endif
        do s = 1, nstates
            proj(s) = sum(evec(:,1)*T(:,s))
        end do
        do s = 1, nstates
            ord(s) = s
        end do
        ! insertion sort on nstates <= 20 entries; keeps the permutation explicit
        do s = 2, nstates
            i = ord(s)
            r = s - 1
            do while( r >= 1 )
                if( proj(ord(r)) <= proj(i) ) exit
                ord(r+1) = ord(r)
                r = r - 1
            end do
            ord(r+1) = i
        end do
        if( proj(ord(1)) > proj(ord(nstates)) ) ord = ord(nstates:1:-1)   ! deterministic sign
        allocate(wtmp(size(weights,1),nstates), ttmp(size(targets,1),nstates))
        allocate(btmp(nstates), ntmp(nstates), inv(nstates))
        do s = 1, nstates
            wtmp(:,s) = weights(:,ord(s))
            ttmp(:,s) = targets(:,ord(s))
            btmp(s)   = bandwidths(ord(s))
            ntmp(s)   = neff(ord(s))
            inv(ord(s)) = s
        end do
        weights = wtmp; targets = ttmp; bandwidths = btmp; neff = ntmp
        do i = 1, size(labels)
            if( labels(i) >= 1 .and. labels(i) <= nstates ) labels(i) = inv(labels(i))
        end do
        write(logfhandle,'(A,F6.3,A)') '>>> FLEX_PCA trajectory frames ordered along the path &
            &(leading direction carries ',frac,' of target variance)'
        call flush(logfhandle)
        deallocate(T, C, ev, evec, proj, ord, wtmp, ttmp, btmp, ntmp, inv)
    end subroutine order_trajectory_targets

    !>  Tied-covariance Gaussian-mixture responsibilities over the placed state targets.
    !!
    !!  Replaces the compact-support Epanechnikov kernel as the assignment rule. The kernel gives a
    !!  particle zero weight in every state whose support it falls outside, so it contributes to no
    !!  map at all -- 48 % of Ribosembly particles at the shipped settings. Widening the kernel to
    !!  fix that pulls every state toward consensus, which is the coverage/distinctness trade the
    !!  nkern experiment ran into. A mixture removes the trade: support is unbounded, so coverage is
    !!  100 % by construction, while the fitted tied covariance and the free means keep the states
    !!  apart.
    !!
    !!  Measured on Ribosembly against ground truth with the probe-refined embedding, holding target
    !!  placement fixed so the comparison isolates the assignment rule
    !!  (`scripts/state_assign_proto.py`):
    !!
    !!    rule                                    ARI    NMI   coverage  purity  distinct GT
    !!    Epanechnikov kernel (previous default) 0.532  0.745    52.0 %   0.770       10/16
    !!    hard nearest-target Voronoi            0.658  0.786   100.0 %   0.747       10/16
    !!    tied GMM, means pinned at the targets  0.781  0.902   100.0 %   0.874       10/16
    !!    tied GMM, means free (this routine)    0.900  0.939   100.0 %   0.945       13/16
    !!
    !!  The gain splits about evenly three ways: coverage, the fitted covariance shape, and letting
    !!  the means move. Plain softmax over the kernel bandwidth is NOT a substitute -- it reaches
    !!  100 % coverage but drops to ARI 0.310 / purity 0.688, because every particle then contributes
    !!  to every state and the maps blur back toward consensus.
    !!
    !!  FREE MEANS ARE SAFE ON A CONTINUUM, which is not obvious: mode / mean-shift refinement is
    !!  ruled out for exactly this reason (it collapsed IgG's 20 targets onto 1 conformation, because
    !!  a continuum has no modes for targets to climb to). A mixture is different in kind -- its
    !!  components repel through the shared responsibilities instead of each climbing the same
    !!  density gradient independently. Measured on IgG-RL, latent proxy, targets scored by nearest
    !!  ground-truth conformation: distinct 6/20 -> 20/20, reaction-coordinate span 27.3 % -> 85.9 %.
    !!  It spreads, it does not collapse.
    !!
    !!  Tied covariance because the within-state spread here is shared measurement error; a covariance
    !!  per component estimates K*nk*(nk+1)/2 parameters the data cannot support and fits noise.
    subroutine gmm_state_weights( z, nptcls, ncomp, nk, nstates, tcen, wcomp, weights, neff, &
        &bandwidths, labels, views, orient_lam )
        integer,  intent(in)    :: nptcls, ncomp, nk, nstates
        real(dp), intent(in)    :: z(nptcls,ncomp), wcomp(nk)
        !> in: the placed targets. out: the FITTED component means, so the reported target table
        !! describes the maps that were actually delivered rather than the k-center seed.
        real(dp), intent(inout) :: tcen(nk,nstates)
        real,     intent(inout) :: weights(nptcls,nstates), neff(nstates), bandwidths(nstates)
        integer,  intent(inout) :: labels(nptcls)
        !> per-particle viewing axis, enables the orientation-coverage term (see orient_lam)
        real(dp), optional, intent(in) :: views(3,nptcls)
        !> weight of the orientation-coverage penalty; 0 or absent = off
        real(dp), optional, intent(in) :: orient_lam
        real(dp), parameter :: GMM_REG = 1.d-6, GMM_TOL = 1.d-5
        !> responsibilities below this are zeroed so the reconstructor's live-state compaction works
        real(dp), parameter :: RESP_FLOOR = 1.d-3
        !> two means closer than this in the tied-covariance metric describe the same state
        real(dp), parameter :: GMM_MERGE_D2   = 1.0d0
        integer,  parameter :: GMM_MAX_RESPAWN = 8
        integer,  parameter :: GMM_MAXIT = 60
        real(dp), allocatable :: y(:,:), mu(:,:), S(:,:), Sinv(:,:), Syy(:,:), resp(:,:)
        real(dp), allocatable :: Smu(:,:), mSm(:), pival(:), nresp(:), evwork(:,:), ev(:), evec(:,:)
        real(dp), allocatable :: ybar(:)
        real(dp), allocatable :: Tdev(:,:,:), Tbar(:,:)
        real(dp) :: ySy, ySm, lmax, lsum, ll, prev_ll, logdet, lam, sumw, sumw2, vq, trS
        real(dp) :: bicval, ent, dmin, d2pair
        integer  :: nfree, nrespawn, kmin, kdrop, iworst
        integer  :: i, q, r, state, it, nrot, errflg
        integer(kind=8) :: nact_tot
        logical  :: l_orient
        lam      = 0.d0
        if( present(orient_lam) ) lam = orient_lam
        l_orient = present(views) .and. lam > 0.d0
        allocate(y(nptcls,nk), mu(nk,nstates), S(nk,nk), Sinv(nk,nk), Syy(nk,nk))
        allocate(resp(nptcls,nstates), Smu(nk,nstates), mSm(nstates), pival(nstates), nresp(nstates))
        allocate(evwork(nk,nk), ev(nk), evec(nk,nk))
        ! Standardised frame. wcomp is 1/var per component, so sqrt(wcomp) IS the 1/sd
        ! standardisation -- do not divide by sd as well, that standardises twice.
        !$omp parallel do default(shared) private(i,q) schedule(static)
        do i = 1, nptcls
            do q = 1, nk
                y(i,q) = z(i,q) * sqrt(wcomp(q))
            end do
        end do
        !$omp end parallel do
        do state = 1, nstates
            do q = 1, nk
                mu(q,state) = tcen(q,state) * sqrt(wcomp(q))
            end do
        end do
        ! Syy is fixed, so the tied-covariance M step below is a rank-nstates correction to it rather
        ! than a second pass over the particles: sum_k sum_i R_ik (y_i-mu_k)(y_i-mu_k)'
        ! = sum_i y_i y_i' - sum_k N_k mu_k mu_k', because sum_k R_ik = 1.
        Syy = matmul(transpose(y), y)
        S   = Syy / real(nptcls,dp)
        allocate(ybar(nk))
        do q = 1, nk
            ybar(q) = sum(y(:,q)) / real(nptcls,dp)
        end do
        do q = 1, nk
            do r = 1, nk
                S(q,r) = S(q,r) - ybar(q)*ybar(r)
            end do
        end do
        deallocate(ybar)
        do q = 1, nk
            S(q,q) = S(q,q) + GMM_REG
        end do
        pival    = 1.d0 / real(nstates,dp)
        prev_ll  = -huge(1.d0)
        nrespawn = 0
        if( l_orient )then
            allocate(Tbar(3,3), Tdev(3,3,nstates), source=0.d0)
            do i = 1, nptcls
                do q = 1, 3
                    do r = 1, 3
                        Tbar(q,r) = Tbar(q,r) + views(q,i)*views(r,i)
                    end do
                end do
            end do
            Tbar = Tbar / real(nptcls,dp)
        endif
        do it = 1, GMM_MAXIT
            call matinv(S, Sinv, nk, errflg)
            if( errflg /= 0 )then
                write(logfhandle,'(A)') '>>> FLEX_PCA GMM tied covariance singular; keeping kernel weights'
                deallocate(y, mu, S, Sinv, Syy, resp, Smu, mSm, pival, nresp, evwork, ev, evec)
                if( allocated(Tbar) ) deallocate(Tbar, Tdev)
                return
            endif
            evwork = S
            call jacobi(evwork, nk, nk, ev, evec, nrot)
            logdet = 0.d0
            do q = 1, nk
                logdet = logdet + log(max(ev(q), DTINY))
            end do
            Smu = matmul(Sinv, mu)
            do state = 1, nstates
                mSm(state) = sum(mu(:,state)*Smu(:,state))
            end do
            ll = 0.d0
            !$omp parallel do default(shared) private(i,q,r,state,ySy,ySm,lmax,lsum,vq) &
            !$omp& schedule(static) reduction(+:ll)
            do i = 1, nptcls
                ySy = 0.d0
                do q = 1, nk
                    do r = 1, nk
                        ySy = ySy + y(i,q)*Sinv(q,r)*y(i,r)
                    end do
                end do
                do state = 1, nstates
                    ySm = 0.d0
                    do q = 1, nk
                        ySm = ySm + y(i,q)*Smu(q,state)
                    end do
                    resp(i,state) = -0.5d0*(ySy - 2.d0*ySm + mSm(state)) - 0.5d0*logdet &
                        &+ log(max(pival(state), DTINY))
                    if( l_orient )then
                        ! ACK-means template with the variable swapped from class SIZE to orientation
                        ! coverage: down-weight a particle whose viewing AXIS points where this state
                        ! is already over-represented. v v', never the mean resultant -- +n and -n are
                        ! mirror projections, so bias lives in the axis and a dipole statistic is blind.
                        vq = 0.d0
                        do q = 1, 3
                            do r = 1, 3
                                vq = vq + views(q,i)*Tdev(q,r,state)*views(r,i)
                            end do
                        end do
                        resp(i,state) = resp(i,state) - lam*vq
                    endif
                end do
                lmax = maxval(resp(i,:))
                lsum = 0.d0
                do state = 1, nstates
                    resp(i,state) = exp(resp(i,state) - lmax)
                    lsum       = lsum + resp(i,state)
                end do
                resp(i,:) = resp(i,:) / max(lsum, DTINY)
                ll     = ll + (log(max(lsum, DTINY)) + lmax)
            end do
            !$omp end parallel do
            ll = ll / real(nptcls,dp)
            ! ---- M step ----
            do state = 1, nstates
                nresp(state) = sum(resp(:,state))
            end do
            pival = max(nresp, DTINY) / real(nptcls,dp)
            do state = 1, nstates
                do q = 1, nk
                    mu(q,state) = sum(resp(:,state)*y(:,q)) / max(nresp(state), DTINY)
                end do
            end do
            S = Syy
            do state = 1, nstates
                do q = 1, nk
                    do r = 1, nk
                        S(q,r) = S(q,r) - nresp(state)*mu(q,state)*mu(r,state)
                    end do
                end do
            end do
            S = S / real(nptcls,dp)
            do q = 1, nk
                S(q,q) = S(q,q) + GMM_REG
            end do
            S = 0.5d0*(S + transpose(S))
            if( l_orient )then
                do state = 1, nstates
                    do q = 1, 3
                        do r = 1, 3
                            Tdev(q,r,state) = sum(resp(:,state)*views(q,:)*views(r,:)) &
                                &/ max(nresp(state), DTINY) - Tbar(q,r)
                        end do
                    end do
                end do
            endif
            if( abs(ll - prev_ll) < GMM_TOL*abs(ll) )then
                ! ---- MERGE REDUNDANT COMPONENTS, RESPAWN WHERE COVERAGE IS WORST ----
                ! Converged EM will happily park several components on one well-populated region and
                ! leave other states with none: measured on Ribosembly, three of sixteen targets landed
                ! on ground-truth state 15 while states 0 and 11 got none -- and the embedding ceiling
                ! says all sixteen were reachable, so this is the mixture wasting targets, not a limit
                ! of the latent. Two means closer than GMM_MERGE_D2 in the tied-covariance metric
                ! describe the same state; keep one and restart the other at the worst-explained
                ! particle, which is by construction in a region no component currently owns.
                if( nrespawn < GMM_MAX_RESPAWN )then
                    kmin = 0; kdrop = 0; dmin = huge(1.d0)
                    do state = 1, nstates - 1
                        do r = state + 1, nstates
                            d2pair = 0.d0
                            do q = 1, nk
                                d2pair = d2pair + (mu(q,state) - mu(q,r))*(Smu(q,state) - Smu(q,r))
                            end do
                            if( d2pair < dmin )then
                                dmin = d2pair
                                kmin = state
                                kdrop = merge(r, state, nresp(r) < nresp(state))
                                if( kdrop == state ) kmin = r
                            endif
                        end do
                    end do
                    if( kdrop >= 1 .and. dmin < GMM_MERGE_D2 )then
                        iworst = maxloc(-maxval(resp, dim=2), dim=1)
                        write(logfhandle,'(A,I0,A,I0,A,F8.4,A,I0)') &
                            &'>>> FLEX_PCA GMM components ',kmin,' and ',kdrop, &
                            &' are redundant (separation ',real(dmin),'); respawning ',kdrop
                        call flush(logfhandle)
                        mu(:,kdrop)    = y(iworst,:)
                        pival(kdrop)   = 1.d0 / real(nstates,dp)
                        nrespawn       = nrespawn + 1
                        prev_ll        = -huge(1.d0)
                        cycle
                    endif
                endif
                exit
            endif
            prev_ll = ll
        end do
        write(logfhandle,'(A,I0,A,ES13.5,A,F7.4,A,F7.4)') &
            &'>>> FLEX_PCA GMM tied-covariance responsibilities: iters=',min(it,GMM_MAXIT), &
            &'  loglik=',ll,'  pi range ',minval(pival),' - ',maxval(pival)
        ! ---- MODEL SELECTION DIAGNOSTICS ----
        ! The kernel never had a likelihood, so the state count could only ever be supplied. A mixture
        ! does, which makes BIC/ICL available for choosing it. ICL is BIC plus the mean-entropy penalty
        ! -2*sum_ik r_ik log r_ik: BIC scores how well the density is fitted and will happily spend
        ! several components on one well-populated state, whereas ICL charges for ambiguous
        ! responsibilities and so prefers components that correspond to SEPARATED states -- which is
        ! the question actually being asked here. Reported per run so a ladder over npreimages can be
        ! read off directly; nothing is selected automatically yet.
        nfree  = nstates*nk + (nk*(nk+1))/2 + nstates - 1
        bicval = -2.d0*ll*real(nptcls,dp) + real(nfree,dp)*log(real(nptcls,dp))
        ent    = 0.d0
        !$omp parallel do default(shared) private(i,state) schedule(static) reduction(+:ent)
        do i = 1, nptcls
            do state = 1, nstates
                if( resp(i,state) > DTINY ) ent = ent - resp(i,state)*log(resp(i,state))
            end do
        end do
        !$omp end parallel do
        write(logfhandle,'(A,I0,A,I0,A,ES14.6,A,ES14.6,A,F8.4)') &
            &'>>> FLEX_PCA GMM model selection: K=',nstates,'  free params=',nfree, &
            &'  BIC=',bicval,'  ICL=',bicval + 2.d0*ent,'  mean entropy=',real(ent/real(nptcls,dp))
        ! A component whose responsibility mass is a small multiple of its own dimension cannot be
        ! estimated and is a candidate for merging; several such in one region is the signature of
        ! targets wasted on the same state.
        write(logfhandle,'(A,I0,A,I0)') '>>> FLEX_PCA GMM components below 3*nk effective mass: ', &
            &count(nresp < 3.d0*real(nk,dp)),' of ',nstates
        call flush(logfhandle)
        if( l_orient ) write(logfhandle,'(A,F7.3)') &
            &'>>> FLEX_PCA GMM orientation-coverage term active, lambda=',lam
        ! SPARSIFY before handing the responsibilities over. This is not an optimisation, it is
        ! required: insert_plane_oversamp_multi_scaled compacts the live states per particle and skips
        ! the plane entirely when none are live, and its comment states the assumption outright --
        ! "a particle typically sits inside only one or two of the nstates targets and the remaining
        ! scale pairs are EXACT zeros". Epanechnikov weights satisfy that by having compact support;
        ! mixture responsibilities are strictly positive everywhere and do not, so every particle
        ! would insert into all nstates maps. Measured consequence of leaving them dense: the state
        ! reconstruction had produced no output after 35 minutes where the kernel path takes ~1.
        ! A responsibility below RESP_FLOOR moves a 335k-particle average by <0.1 % of one particle,
        ! so this is numerically free; renormalise so each particle still carries unit total mass.
        nact_tot = 0
        !$omp parallel do default(shared) private(i,state,lsum) schedule(static) reduction(+:nact_tot)
        do i = 1, nptcls
            lsum = 0.d0
            do state = 1, nstates
                if( resp(i,state) < RESP_FLOOR )then
                    resp(i,state) = 0.d0
                else
                    lsum     = lsum + resp(i,state)
                    nact_tot = nact_tot + 1
                endif
            end do
            if( lsum > DTINY )then
                resp(i,:) = resp(i,:) / lsum
            else
                ! everything fell under the floor: keep the single best state so the particle still counts
                state           = maxloc(resp(i,:), dim=1)
                resp(i,:)       = 0.d0
                resp(i,state)   = 1.d0
                nact_tot        = nact_tot + 1
            endif
        end do
        !$omp end parallel do
        write(logfhandle,'(A,F6.2,A,ES9.2)') '>>> FLEX_PCA GMM mean active states per particle=', &
            &real(nact_tot)/real(max(nptcls,1)),'  (responsibility floor ',RESP_FLOOR
        call flush(logfhandle)
        ! responsibilities ARE the weights: the reconstructor uses them as both the data and the
        ! density scale, so each state map is a weighted average and the per-state scale divides out
        trS = 0.d0
        do q = 1, nk
            trS = trS + S(q,q)
        end do
        trS = sqrt(trS / real(nk,dp))          ! reported in the bandwidth column: the tied scale
        do state = 1, nstates
            weights(:,state) = real(resp(:,state))
            sumw             = sum(resp(:,state))
            sumw2            = sum(resp(:,state)**2)
            neff(state)      = real(sumw*sumw / max(sumw2, DTINY))
            bandwidths(state)= real(trS)
        end do
        do i = 1, nptcls
            labels(i) = maxloc(resp(i,:), dim=1)
        end do
        ! hand the fitted means back out of the standardised frame
        do state = 1, nstates
            do q = 1, nk
                tcen(q,state) = mu(q,state) / sqrt(max(wcomp(q), DTINY))
            end do
        end do
        deallocate(y, mu, S, Sinv, Syy, resp, Smu, mSm, pival, nresp, evwork, ev, evec)
        if( allocated(Tbar) ) deallocate(Tbar, Tdev)
    end subroutine gmm_state_weights

    !> Epanechnikov weights for ONE state at a GIVEN bandwidth, from precomputed distances.
    subroutine kernel_weights_at_bandwidth( dist, nptcls, h_in, min_neff, w, h_out, neff_out )
        integer,  intent(in)  :: nptcls, min_neff
        real(dp), intent(in)  :: dist(nptcls), h_in
        real,     intent(out) :: w(nptcls)
        real(dp), intent(out) :: h_out
        real,     intent(out) :: neff_out
        real(dp) :: h, u2, sumw, sumw2
        integer  :: i, grow, nsupp
        h = h_in
        ! same capped widening loop as build_covariance_state_weights
        nsupp = 0
        do grow = 0, COV_MAX_BW_GROW
            sumw = 0.d0; sumw2 = 0.d0; nsupp = 0
            !$omp parallel do default(shared) private(i,u2) schedule(static) &
            !$omp& reduction(+:sumw,sumw2,nsupp)
            do i = 1, nptcls
                u2 = dist(i) / (h*h)
                w(i) = real(max(0.d0, 1.d0 - u2))
                sumw  = sumw  + real(w(i),dp)
                sumw2 = sumw2 + real(w(i),dp)**2
                if( w(i) > 0. ) nsupp = nsupp + 1
            end do
            !$omp end parallel do
            if( nsupp >= min(min_neff, nptcls) ) exit
            if( grow >= COV_MAX_BW_GROW      ) exit
            h = 1.3d0*h
        end do
        if( maxval(w) > TINY ) w = w / maxval(w)
        sumw  = sum(real(w,dp))
        sumw2 = sum(real(w,dp)**2)
        h_out    = h
        neff_out = real(sumw*sumw/max(sumw2,DTINY))
    end subroutine kernel_weights_at_bandwidth

    !> Cross-validated bandwidth selection. For each state and each trial bandwidth, reconstruct even and
    !! odd half maps and score them against the NARROWEST bin's opposite-half map, symmetrized. Scoring on
    !! plain even/odd agreement instead would rise monotonically with bandwidth and pick the widest bin
    !! every time -- maximal smearing.
    subroutine cv_select_bandwidths( params, build, pinds, nptcls, nstates, nbins, min_neff, &
        &dist, bfloor, weights, bandwidths, neff )
        class(parameters), intent(inout) :: params
        class(builder),    intent(inout) :: build
        integer,           intent(in)    :: pinds(:), nptcls, nstates, nbins, min_neff
        real(dp),          intent(in)    :: dist(nptcls,nstates), bfloor(nstates)
        real, allocatable, intent(inout) :: weights(:,:), bandwidths(:), neff(:)
        real,     allocatable :: wbin(:,:), whalf(:,:), sorted(:)
        real(dp), allocatable :: bins(:,:), hbin(:,:), err(:,:)
        real,     allocatable :: tgt_ev(:,:), tgt_od(:,:)     ! narrow-bin targets per state
        real,     allocatable :: rmat(:,:,:)
        type(image) :: ev, od
        type(string):: fn
        real(dp) :: b_lo, b_hi, t, h_used, e1, e2
        real     :: neff_used
        integer  :: state, ib, p95, nvox, ibest
        character(len=3) :: bstr
        allocate(wbin(nptcls,nstates), whalf(nptcls,nstates), bins(nbins,nstates), &
            &hbin(nbins,nstates), err(nbins,nstates), sorted(nptcls))
        err = 0.d0
        ! --- bin grid per state: from the bandwidth floor to the 95th distance percentile, linear in sqrt
        p95 = max(1, min(nptcls, nint(0.95*real(nptcls))))
        do state = 1, nstates
            sorted = real(dist(:,state))
            call hpsort(sorted)
            b_lo = bfloor(state)
            b_hi = max(real(sorted(p95),dp), b_lo*1.0001d0)
            do ib = 1, nbins
                t = real(ib-1,dp)/real(max(1,nbins-1),dp)
                bins(ib,state) = (sqrt(b_lo) + t*(sqrt(b_hi)-sqrt(b_lo)))**2   ! linear in sqrt
            end do
        end do
        write(logfhandle,'(A,I0,A)') '>>> FLEX_PCA cross-validated bandwidth selection over ',nbins, &
            &' bins per state'
        call flush(logfhandle)
        nvox = params%box_crop**3
        allocate(tgt_ev(nvox,nstates), tgt_od(nvox,nstates), source=0.)
        do ib = 1, nbins
            do state = 1, nstates
                call kernel_weights_at_bandwidth(dist(:,state), nptcls, sqrt(2.d0*bins(ib,state)), &
                    &min_neff, wbin(:,state), h_used, neff_used)
                hbin(ib,state) = h_used
            end do
            write(bstr,'(I3.3)') ib
            whalf = wbin
            call mask_state_weights_by_half(build, pinds, 0, whalf)
            params%outvol = 'flex_pca_cv'//bstr//'_even_state_001.mrc'
            call reconstruct_flex_weighted_states(params, build, pinds, whalf, nstates, floor_rho=.true.)
            whalf = wbin
            call mask_state_weights_by_half(build, pinds, 1, whalf)
            params%outvol = 'flex_pca_cv'//bstr//'_odd_state_001.mrc'
            call reconstruct_flex_weighted_states(params, build, pinds, whalf, nstates, floor_rho=.true.)
            do state = 1, nstates
                ! trial half maps are written at box_rec; the CV score is computed at box_crop
                fn = 'flex_pca_cv'//bstr//'_even_state_'//int2str_pad(state,3)//MRC_EXT
                call ev%read_and_crop(fn, flex_rec_smpd(params), params%box_crop, params%smpd_crop)
                call del_file(fn%to_char()); call fn%kill
                fn = 'flex_pca_cv'//bstr//'_odd_state_'//int2str_pad(state,3)//MRC_EXT
                call od%read_and_crop(fn, flex_rec_smpd(params), params%box_crop, params%smpd_crop)
                call del_file(fn%to_char()); call fn%kill
                if( ib == 1 )then
                    rmat = ev%get_rmat(); tgt_ev(:,state) = reshape(rmat, [nvox])
                    rmat = od%get_rmat(); tgt_od(:,state) = reshape(rmat, [nvox])
                endif
                ! symmetrized cross-halfset error against the OPPOSITE half's narrow target
                rmat = od%get_rmat()
                e1 = sum((real(tgt_ev(:,state),dp) - real(reshape(rmat,[nvox]),dp))**2)
                rmat = ev%get_rmat()
                e2 = sum((real(reshape(rmat,[nvox]),dp) - real(tgt_od(:,state),dp))**2)
                err(ib,state) = e1 + e2
                call ev%kill; call od%kill
            end do
            write(logfhandle,'(A,I3,A,ES11.3,A,ES12.4,A,ES12.4)') '>>>   cv bin=',ib,' h(state1)=',hbin(ib,1), &
                &'  cross-halfset error: min=',minval(err(ib,:)),' max=',maxval(err(ib,:))
            call flush(logfhandle)
        end do
        ! --- pick the minimum-error bin per state and rebuild the final weights ---
        write(logfhandle,'(A)') '>>> FLEX_PCA selected bandwidths (per state, argmin cross-halfset error):'
        do state = 1, nstates
            ibest = minloc(err(:,state), dim=1)
            call kernel_weights_at_bandwidth(dist(:,state), nptcls, sqrt(2.d0*bins(ibest,state)), &
                &min_neff, weights(:,state), h_used, neff_used)
            bandwidths(state) = real(h_used)
            neff(state)       = neff_used
            write(logfhandle,'(A,I3,A,I3,A,I0,A,ES11.3,A,ES12.4,A,F10.1)') '>>>   state=',state,'  bin=',ibest,'/',nbins, &
                &'  h=',h_used,'  error=',err(ibest,state),'  neff=',neff_used
        end do
        call flush(logfhandle)
        deallocate(wbin, whalf, bins, hbin, err, sorted, tgt_ev, tgt_od)
    end subroutine cv_select_bandwidths

    !>  True only when an environment flag is set to a nonzero integer (an opt-IN switch).
    logical function cov_env_flag_on( name ) result(on)
        character(len=*), intent(in) :: name
        character(len=32) :: envval
        integer :: stat, ln, ival
        on = .false.
        call get_environment_variable(name, envval, ln, stat)
        if( stat /= 0 .or. ln < 1 ) return
        read(envval(:ln), *, iostat=stat) ival
        if( stat == 0 ) on = ival /= 0
    end function cov_env_flag_on

    !>  Read a real-valued override from the environment, leaving `val` untouched when unset or
    !!  unparseable. Companion to cov_env_int, which is integer-only and positive-only.
    subroutine cov_env_dp( name, val )
        character(len=*), intent(in)    :: name
        real(dp),         intent(inout) :: val
        character(len=32) :: envval
        integer  :: stat, ln_env
        real(dp) :: rval
        call get_environment_variable(name, envval, ln_env, stat)
        if( stat /= 0 .or. ln_env < 1 ) return
        read(envval(:ln_env), *, iostat=stat) rval
        if( stat == 0 )then
            val = rval
            write(logfhandle,'(A,A,A,ES12.4)') '>>> FLEX_PCA ',trim(name),' override: ',rval
            call flush(logfhandle)
        endif
    end subroutine cov_env_dp

    !>  True only when an environment flag is explicitly set to zero (an opt-OUT switch).
    logical function cov_env_flag_off( name ) result(off)
        character(len=*), intent(in) :: name
        character(len=32) :: envval
        integer :: stat, ln, ival
        off = .false.
        call get_environment_variable(name, envval, ln, stat)
        if( stat /= 0 .or. ln < 1 ) return
        read(envval(:ln), *, iostat=stat) ival
        if( stat == 0 ) off = ival == 0
    end function cov_env_flag_off

    !> Nonuniform local-resolution filtering of the delivered state and trajectory maps, using ONE filter
    !! derived from the CONSENSUS half maps.
    subroutine apply_consensus_nu_filter( params, nstates )
        use simple_nu_filter, only: setup_nu_dmats, optimize_nu_cutoff_finds, nu_filter_vol, &
            &cleanup_nu_filter, write_nu_local_resolution_map, get_nu_filtmap_finest_selected_lp
        class(parameters), intent(in) :: params
        integer,           intent(in) :: nstates
        type(image)  :: vol_e, vol_o, vin, vout
        type(string) :: fn, fn_out, spe, spo
        character(len=:), allocatable :: pe, po, base
        integer :: s, ldim(3), ipass
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
        do ipass = 1, 2
            do s = 1, nstates
                if( ipass == 1 )then
                    fn     = 'flex_pca_state_'//int2str_pad(s,3)//MRC_EXT
                    fn_out = 'flex_pca_state_'//int2str_pad(s,3)//'_nu'//MRC_EXT
                else
                    fn     = 'flex_pca_traj_'//int2str_pad(s,3)//MRC_EXT
                    fn_out = 'flex_pca_traj_'//int2str_pad(s,3)//'_nu'//MRC_EXT
                endif
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
        end do
        call cleanup_nu_filter()
        call vol_e%kill; call vol_o%kill
        write(logfhandle,'(A)') '>>> FLEX_PCA nonuniform-filtered maps written as *_nu.mrc &
            &(originals retained); local resolution map: flex_pca_nu_locres.mrc'
        call flush(logfhandle)
    end subroutine apply_consensus_nu_filter

    !> Rotate the covariance eigenbasis toward SPATIALLY COHERENT components: maximise the Rayleigh
    !! quotient of each component against its own smooth_lp-smoothed copy, on the molecule, as a
    !! generalized symmetric eigenproblem against the basis Gram. The latents and their second moments
    !! are transformed with it, so U' z' reproduces U z exactly and every downstream stage sees one frame.
    subroutine rotate_basis_by_smoothness( params, ncomp, smooth_lp, z, nptcls, latent_second, ok, eigdir )
        class(parameters), intent(in)    :: params
        integer,           intent(in)    :: ncomp, nptcls
        real,              intent(in)    :: smooth_lp        ! Gaussian low-pass, Angstrom
        real(dp),          intent(inout) :: z(nptcls,ncomp)
        real(dp),          intent(inout) :: latent_second(ncomp,ncomp,nptcls)
        logical,           intent(out)   :: ok
        ! a resumed run executes in a fresh job dir; the eigenvolumes sit beside the cache
        character(len=*), optional, intent(in) :: eigdir
        character(len=:), allocatable :: pcdir
        type(image),  allocatable :: eigs(:), sm(:)
        real(dp),     allocatable :: M(:,:), N(:,:), L(:,:), C(:,:), Q(:,:), Rot(:,:), Rinv(:,:)
        real(dp),     allocatable :: ev(:), tmp(:,:), zrow(:)
        real,         allocatable :: rmat(:,:,:), vsort(:)
        logical,      allocatable :: msk(:,:,:)
        type(image)  :: refvol
        type(string) :: fn
        integer :: ik, iq, ir, ii, ldim(3), nvox, ord(ncomp)
        real(dp) :: acc, thr
        ok = .false.
        if( ncomp < 2 .or. smooth_lp <= 0. ) return
        pcdir = ''
        if( present(eigdir) ) pcdir = eigdir
        if( .not. file_exists(pcdir//'flex_pca_pc001.mrc') ) return
        ! molecule mask from the consensus the run was given
        call refvol%new([params%box_crop,params%box_crop,params%box_crop], params%smpd_crop)
        call refvol%read_and_crop(params%vols(1), params%smpd, params%box_crop, params%smpd_crop)
        rmat = refvol%get_rmat(); ldim = shape(rmat); nvox = product(ldim)
        allocate(vsort(nvox)); vsort = reshape(rmat, [nvox]); call hpsort(vsort)
        thr = real(vsort(max(1,min(nvox,nint(0.92*real(nvox))))), dp)
        allocate(msk(ldim(1),ldim(2),ldim(3)))
        msk = rmat > real(thr)
        deallocate(vsort)
        call refvol%kill
        ! read the eigenvolumes and their smoothed copies
        allocate(eigs(ncomp), sm(ncomp))
        do ik = 1, ncomp
            call eigs(ik)%new([params%box_crop,params%box_crop,params%box_crop], params%smpd_crop)
            fn = pcdir//'flex_pca_pc'//int2str_pad(ik,3)//MRC_EXT
            call eigs(ik)%read(fn); call fn%kill
            call sm(ik)%copy(eigs(ik))
            call sm(ik)%fft
            call sm(ik)%bpgau3D(0., smooth_lp)
            call sm(ik)%ifft
        end do
        allocate(M(ncomp,ncomp), N(ncomp,ncomp), source=0.d0)
        do iq = 1, ncomp
            do ir = 1, ncomp
                M(iq,ir) = sum(real(pack(eigs(iq)%get_rmat(), msk),dp) * real(pack(sm(ir)%get_rmat(), msk),dp))
                N(iq,ir) = sum(real(pack(eigs(iq)%get_rmat(), msk),dp) * real(pack(eigs(ir)%get_rmat(), msk),dp))
            end do
        end do
        M = 0.5d0*(M + transpose(M)); N = 0.5d0*(N + transpose(N))
        ! generalized symmetric eigenproblem M x = lambda N x, by Cholesky reduction of N
        allocate(L(ncomp,ncomp), source=0.d0)
        do iq = 1, ncomp
            acc = N(iq,iq) - sum(L(iq,1:iq-1)**2)
            if( acc <= 0.d0 )then
                write(logfhandle,'(A)') '>>> FLEX_PCA smoothness rotation: basis Gram not positive definite, skipping'
                deallocate(M,N,L,msk); do ik=1,ncomp; call eigs(ik)%kill; call sm(ik)%kill; end do
                deallocate(eigs,sm); return
            endif
            L(iq,iq) = sqrt(acc)
            do ir = iq+1, ncomp
                L(ir,iq) = (N(ir,iq) - sum(L(ir,1:iq-1)*L(iq,1:iq-1))) / L(iq,iq)
            end do
        end do
        call invert_lower(L, ncomp)                  ! L now holds L^-1
        allocate(C(ncomp,ncomp), Q(ncomp,ncomp), ev(ncomp), tmp(ncomp,ncomp))
        tmp = matmul(L, M); C = matmul(tmp, transpose(L))
        C   = 0.5d0*(C + transpose(C))
        Q   = C
        call jacobi(Q, ncomp, ncomp, ev, tmp, ik)
        call eigsrt(ev, tmp, ncomp, ncomp)           ! descending smoothness
        allocate(Rot(ncomp,ncomp), Rinv(ncomp,ncomp))
        Rot = matmul(transpose(L), tmp)              ! R = L^-T Q
        ! z' = z * (R^-1)^T so that U' z' reproduces U z exactly; R^-1 = Q^T L
        Rinv = matmul(transpose(tmp), L)
        allocate(zrow(ncomp))
        do ii = 1, nptcls
            zrow = z(ii,:)
            do iq = 1, ncomp
                z(ii,iq) = sum(zrow * Rinv(iq,:))
            end do
        end do
        ! precision transforms as Pi' = R^T Pi R
        do ii = 1, nptcls
            latent_second(:,:,ii) = matmul(transpose(Rot), matmul(latent_second(:,:,ii), Rot))
        end do
        ! rewrite the eigenvolumes in the rotated frame
        do ik = 1, ncomp
            call sm(ik)%zero_and_unflag_ft
            do iq = 1, ncomp
                if( abs(Rot(iq,ik)) < DTINY ) cycle
                call sm(ik)%add(eigs(iq), real(Rot(iq,ik)))
            end do
        end do
        do ik = 1, ncomp
            fn = 'flex_pca_pc'//int2str_pad(ik,3)//MRC_EXT
            call sm(ik)%write(fn, del_if_exists=.true.); call fn%kill
        end do
        write(logfhandle,'(A,F6.1,A,ES11.3,A,ES11.3)') '>>> FLEX_PCA smoothness rotation applied (lp=', &
            &smooth_lp,' A); smoothness eigenvalues max=',ev(1),' min=',ev(ncomp)
        call flush(logfhandle)
        do ik = 1, ncomp
            call eigs(ik)%kill; call sm(ik)%kill
        end do
        deallocate(eigs, sm, M, N, L, C, Q, ev, tmp, Rot, Rinv, zrow, msk)
        ok = .true.
    end subroutine rotate_basis_by_smoothness

    !>  In-place inverse of a lower-triangular matrix (forward substitution per column).
    subroutine invert_lower( L, n )
        integer,  intent(in)    :: n
        real(dp), intent(inout) :: L(n,n)
        real(dp), allocatable :: X(:,:)
        integer  :: i, j, k
        real(dp) :: s
        allocate(X(n,n), source=0.d0)
        do j = 1, n
            X(j,j) = 1.d0 / L(j,j)
            do i = j+1, n
                s = 0.d0
                do k = j, i-1
                    s = s + L(i,k)*X(k,j)
                end do
                X(i,j) = -s / L(i,i)
            end do
        end do
        L = X
        deallocate(X)
    end subroutine invert_lower

    !> FINCH state placement: cluster the latents and take one representative per cluster as a state
    !! target, so the state COUNT comes from the data instead of npreimages. Clustering runs on a
    !! standardized strided subsample -- the kd-tree first-neighbour search is the expensive part, and
    !! standardizing keeps the largest-variance component from dominating the neighbour graph.
    subroutine finch_state_targets( z, nptcls, ncomp, nstates_max, centroids, nfound, ok )
        integer,  intent(in)  :: nptcls, ncomp, nstates_max
        real(dp), intent(in)  :: z(nptcls,ncomp)
        real(dp), allocatable, intent(out) :: centroids(:,:)
        integer,  intent(out) :: nfound
        logical,  intent(out) :: ok
        integer, parameter :: NFINCH_MAX = 20000
        type(finch_hierarchy) :: hier
        real,    allocatable  :: feats(:,:)
        integer, allocatable  :: sub(:), labels(:), reps(:)
        real(dp), allocatable :: sd(:)
        real(dp) :: mu
        integer  :: nsub, i, q, lev, k
        ok = .false.; nfound = 0
        if( nptcls < 100 .or. ncomp < 1 ) return
        nsub = min(nptcls, NFINCH_MAX)
        allocate(sub(nsub))
        do i = 1, nsub
            sub(i) = 1 + int(real(i-1,dp)*real(nptcls-1,dp)/real(max(1,nsub-1),dp))
        end do
        allocate(feats(ncomp,nsub), sd(ncomp))
        do q = 1, ncomp
            mu    = sum(z(:,q)) / real(nptcls,dp)
            sd(q) = sqrt(max(sum((z(:,q)-mu)**2)/real(nptcls,dp), DTINY))
            do i = 1, nsub
                feats(q,i) = real((z(sub(i),q) - mu) / sd(q))
            end do
        end do
        call fit_finch(feats, hier)
        if( hier%get_nlevels() < 1 )then
            call hier%kill; deallocate(feats, sd, sub); return
        endif
        ! report the whole hierarchy, so the chosen level can be judged against the alternatives
        write(logfhandle,'(A)',advance='no') '>>> FLEX_PCA FINCH hierarchy (clusters per level): '
        do i = 1, hier%get_nlevels()
            write(logfhandle,'(I0,1X)',advance='no') hier%get_nclusters(i)
        end do
        write(logfhandle,*)
        ! SIMPLE_COV_FINCH_LEVEL=N takes level N's NATIVE partition directly, no Ward merge and no cap, so
        ! the cluster count is genuinely FINCH's own.
        lev = 0
        call cov_env_int_pub('SIMPLE_COV_FINCH_LEVEL', lev)
        if( lev >= 1 .and. lev <= hier%get_nlevels() )then
            call hier%get_labels(lev, labels)
            k = maxval(labels)
            write(logfhandle,'(A,I0,A,I0,A)') '>>> FLEX_PCA FINCH level ',lev, &
                &' taken natively: ',k,' clusters (no Ward merge, no cap)'
        else
            call select_finch_level(feats, hier, lev)
            call refine_finch_level(feats, hier, lev, labels, k, max_clusters=nstates_max)
        endif
        if( k < 2 )then
            call hier%kill; deallocate(feats, sd, sub); if( allocated(labels) ) deallocate(labels); return
        endif
        call finch_representatives(feats, labels, reps)
        nfound = min(k, size(reps))
        allocate(centroids(ncomp,nfound))
        do i = 1, nfound
            do q = 1, ncomp
                centroids(q,i) = z(sub(reps(i)),q)
            end do
        end do
        write(logfhandle,'(A,I0,A,I0,A,I0,A,I0)') '>>> FLEX_PCA FINCH state placement: levels=', &
            &hier%get_nlevels(),' selected level=',lev,' AUTOMATIC cluster count=',nfound, &
            &' from a subsample of ',nsub
        call flush(logfhandle)
        ok = .true.
        call hier%kill
        deallocate(feats, sd, sub, labels, reps)
    end subroutine finch_state_targets

    !> Parse a comma-separated component list from an environment variable (e.g. "4,3,9"), falling back
    !! to a single default axis when unset or unparseable. Out-of-range entries are dropped.
    subroutine cov_parse_axes( name, ncomp, axis_dflt, axes, naxes )
        character(len=*),     intent(in)  :: name
        integer,              intent(in)  :: ncomp, axis_dflt
        integer, allocatable, intent(out) :: axes(:)
        integer,              intent(out) :: naxes
        character(len=128) :: envval
        integer :: stat, ln, i, j, v, tmp(32)
        naxes = 0
        call get_environment_variable(name, envval, ln, stat)
        if( stat == 0 .and. ln > 0 )then
            i = 1
            do while( i <= ln .and. naxes < 32 )
                j = index(envval(i:ln), ',')
                if( j == 0 )then
                    j = ln + 1
                else
                    j = i + j - 1
                endif
                read(envval(i:j-1), *, iostat=stat) v
                if( stat == 0 .and. v >= 1 .and. v <= ncomp )then
                    naxes = naxes + 1; tmp(naxes) = v
                endif
                i = j + 1
            end do
        endif
        if( naxes < 1 )then
            naxes = 1; tmp(1) = max(1, min(ncomp, axis_dflt))
        endif
        allocate(axes(naxes), source=tmp(1:naxes))
    end subroutine cov_parse_axes

    !> Anchor direction for the trajectory: the leading standardized-correlation eigenvector over the
    !! given axes, expressed back in raw latent units. A motion is typically spread over several
    !! components, so a multi-axis span reaches further along it than any single axis can.
    subroutine build_anchor_direction( z, nptcls, ncomp, axes, naxes, wdir )
        integer,  intent(in)  :: nptcls, ncomp, naxes, axes(naxes)
        real(dp), intent(in)  :: z(nptcls,ncomp)
        real(dp), intent(out) :: wdir(ncomp)
        real(dp), allocatable :: C(:,:), evec(:,:), eval(:), sdv(:)
        real(dp) :: mu_a, mu_b, acc
        integer  :: ia, ib, i, nrot
        wdir = 0.d0
        if( naxes < 1 ) return
        if( naxes == 1 )then
            wdir(axes(1)) = 1.d0
            return
        endif
        allocate(C(naxes,naxes), evec(naxes,naxes), eval(naxes), sdv(naxes))
        do ia = 1, naxes
            mu_a = sum(z(:,axes(ia))) / real(nptcls,dp)
            sdv(ia) = sqrt(max(sum((z(:,axes(ia))-mu_a)**2)/real(nptcls,dp), DTINY))
        end do
        do ia = 1, naxes
            mu_a = sum(z(:,axes(ia))) / real(nptcls,dp)
            do ib = ia, naxes
                mu_b = sum(z(:,axes(ib))) / real(nptcls,dp)
                acc  = 0.d0
                do i = 1, nptcls
                    acc = acc + ((z(i,axes(ia))-mu_a)/sdv(ia)) * ((z(i,axes(ib))-mu_b)/sdv(ib))
                end do
                C(ia,ib) = acc / real(nptcls,dp)
                C(ib,ia) = C(ia,ib)
            end do
        end do
        call jacobi(C, naxes, naxes, eval, evec, nrot)
        call eigsrt(eval, evec, naxes, naxes)
        ! leading direction, expressed back in the RAW latent units the kernel works in
        do ia = 1, naxes
            wdir(axes(ia)) = evec(ia,1) / sdv(ia)
        end do
        acc = sqrt(sum(wdir*wdir))
        if( acc > DTINY ) wdir = wdir / acc
        deallocate(C, evec, eval, sdv)
    end subroutine build_anchor_direction

    subroutine finch_geodesic_targets( z, nptcls, ncomp, nstates, comp_rank, ngeo, axis_sel, centroids, ok )
        integer,  intent(in)  :: nptcls, ncomp, nstates, ngeo
        real(dp), intent(in)  :: z(nptcls,ncomp)
        integer,  intent(in)  :: comp_rank(ncomp)     ! components, most reproducible first
        integer,  intent(in)  :: axis_sel             ! component defining the trajectory ENDPOINTS
        real(dp), intent(out) :: centroids(ncomp,nstates)
        logical,  intent(out) :: ok
        integer, parameter :: NFINCH_MAX = 20000      ! kd-tree 1-NN cost; a subsample fixes the manifold
        integer, parameter :: NFINCH_NODE_CAP = 1200  ! graph work is O(nrep^2); this bounds it
        type(finch_hierarchy) :: hier
        real,     allocatable :: feats(:,:)
        integer,  allocatable :: labels(:), reps(:), sub(:), path(:), prev(:), keep(:)
        real(dp), allocatable :: rz(:,:), dist(:), sd(:), arc(:)
        real,     allocatable :: zax(:)
        real(dp), allocatable :: wdir(:), anchor(:)
        integer,  allocatable :: axes(:)
        integer :: naxes
        logical,  allocatable :: done(:)
        integer :: nsub, nd, i, j, q, s, k, knn, nrep, lev, a, b, u, v, ncon, npath, best_i, nkeep
        real(dp) :: d2, dmax, w, target_arc, t, mu, qlo, qhi
        ok = .false.
        nd = max(1, min(ncomp, ngeo))
        if( nstates < 2 .or. nptcls < 50 ) return
        ! ---- subsample + standardized feature table over the retained components ----
        nsub = min(nptcls, NFINCH_MAX)
        allocate(sub(nsub))
        do i = 1, nsub
            sub(i) = 1 + int(real(i-1,dp)*real(nptcls-1,dp)/real(max(1,nsub-1),dp))
        end do
        allocate(feats(nd,nsub), sd(nd))
        do q = 1, nd
            k  = comp_rank(q)
            mu = sum(z(:,k)) / real(nptcls,dp)
            sd(q) = sqrt(max(sum((z(:,k)-mu)**2)/real(nptcls,dp), DTINY))
            do i = 1, nsub
                feats(q,i) = real((z(sub(i),k) - mu) / sd(q))
            end do
        end do
        ! ---- FINCH -> cluster medoids ----
        call fit_finch(feats, hier)
        if( hier%get_nlevels() < 1 )then
            call hier%kill; deallocate(feats, sd, sub); return
        endif
        ! FINEST level that still fits the O(nrep^2) graph work.
        lev = hier%get_nlevels()
        do i = hier%get_nlevels(), 1, -1
            if( hier%get_nclusters(i) <= NFINCH_NODE_CAP ) lev = i
        end do
        if( hier%get_nclusters(lev) < max(3*nstates, 20) )then
            ! every level below the cap is too coarse; take the finest available
            lev = 1
        endif
        call hier%get_labels(lev, labels)
        call finch_representatives(feats, labels, reps)
        nrep = size(reps)
        write(logfhandle,'(A,I0,A,I0,A,I0,A,I0)') '>>> FLEX_PCA FINCH levels=',hier%get_nlevels(), &
            &' level=',lev,' nodes=',nrep,' over components=',nd
        if( nrep < nstates )then
            call hier%kill; deallocate(feats, sd, sub, labels, reps); return
        endif
        ! node coordinates in the standardized subspace
        allocate(rz(nd,nrep))
        do i = 1, nrep
            rz(:,i) = real(feats(:,reps(i)), dp)
        end do
        ! ---- kNN graph, grown until connected ----
        allocate(dist(nrep), prev(nrep), done(nrep), arc(nrep))
        knn = max(3, min(nrep-1, 6))
        do
            call dijkstra_knn(rz, nd, nrep, knn, 1, dist, prev)
            ncon = count(dist < huge(1.d0)*0.5d0)
            if( ncon == nrep .or. knn >= nrep-1 ) exit
            knn = min(nrep-1, knn*2)
        end do
        if( ncon < nrep )then
            write(logfhandle,'(A,I0,A,I0,A)') '>>> FLEX_PCA FINCH graph stayed disconnected (',ncon, &
                &' of ',nrep,' nodes); falling back to axis placement'
            call hier%kill; deallocate(feats,sd,sub,labels,reps,rz,dist,prev,done,arc); return
        endif
        ! Endpoints come from an ANCHOR direction, not from the longest geodesic: in a multi-component
        ! subspace the longest geodesic runs along whichever direction happens to be most extended, which
        ! need not be the motion, and the path then doubles back on it.
        call cov_parse_axes('SIMPLE_COV_TRAJAXES', ncomp, axis_sel, axes, naxes)
        allocate(wdir(ncomp), anchor(nptcls))
        call build_anchor_direction(z, nptcls, ncomp, axes, naxes, wdir)
        do i = 1, nptcls
            anchor(i) = sum(z(i,:)*wdir)
        end do
        if( naxes > 1 )then
            write(logfhandle,'(A,20(I0,1X))') '>>> FLEX_PCA trajectory anchored on the leading direction of components: ', &
                &(axes(i), i=1,naxes)
        endif
        allocate(zax(nptcls))
        zax = real(anchor)
        call hpsort(zax)
        qlo = real(zax(max(1, min(nptcls, nint(0.05*real(nptcls))))), dp)
        qhi = real(zax(max(1, min(nptcls, nint(0.95*real(nptcls))))), dp)
        deallocate(zax)
        a = 1; b = 1; dmax = huge(1.d0); w = huge(1.d0)
        do i = 1, nrep
            d2 = abs(anchor(sub(reps(i))) - qlo)
            if( d2 < dmax )then
                dmax = d2; a = i
            endif
            d2 = abs(anchor(sub(reps(i))) - qhi)
            if( d2 < w )then
                w = d2; b = i
            endif
        end do
        if( a == b )then
            call hier%kill; deallocate(feats,sd,sub,labels,reps,rz,dist,prev,done,arc); return
        endif
        call dijkstra_knn(rz, nd, nrep, knn, a, dist, prev)
        if( dist(b) >= huge(1.d0)*0.5d0 )then
            call hier%kill; deallocate(feats,sd,sub,labels,reps,rz,dist,prev,done,arc); return
        endif
        ! ---- extract the a->b path ----
        npath = 0; v = b
        do
            npath = npath + 1
            if( v == a .or. npath > nrep ) exit
            v = prev(v)
        end do
        if( npath < 2 )then
            call hier%kill; deallocate(feats,sd,sub,labels,reps,rz,dist,prev,done,arc); return
        endif
        allocate(path(npath))
        v = b
        do i = npath, 1, -1
            path(i) = v
            if( v == a ) exit
            v = prev(v)
        end do
        ! ---- resample the path to nstates targets by arc length ----
        arc(1) = 0.d0
        do i = 2, npath
            d2 = 0.d0
            do q = 1, nd
                d2 = d2 + (rz(q,path(i)) - rz(q,path(i-1)))**2
            end do
            arc(i) = arc(i-1) + sqrt(d2)
        end do
        if( arc(npath) <= DTINY )then
            call hier%kill; deallocate(feats,sd,sub,labels,reps,rz,dist,prev,done,arc,path); return
        endif
        ! Keep only the path nodes that ADVANCE monotonically along the anchor, so a path that doubles
        ! back does not place several states on the same stretch, then INTERPOLATE the targets between
        ! consecutive kept nodes in FULL latent coordinates. Snapping each state to its nearest node
        ! instead would assign several states to the SAME particle, which reconstruct to identical volumes.
        allocate(keep(npath))
        nkeep = 1; keep(1) = path(1)
        w = anchor(sub(reps(path(1))))
        do i = 2, npath
            d2 = anchor(sub(reps(path(i))))
            if( (qhi > qlo .and. d2 > w) .or. (qhi <= qlo .and. d2 < w) )then
                nkeep = nkeep + 1; keep(nkeep) = path(i); w = d2
            endif
        end do
        if( nkeep < 2 )then
            call hier%kill; deallocate(feats,sd,sub,labels,reps,rz,dist,prev,done,arc,path,keep); return
        endif
        centroids = 0.d0
        do s = 1, nstates
            target_arc = qlo + (qhi - qlo) * real(s-1,dp) / real(max(1,nstates-1),dp)
            best_i = 1
            do i = 1, nkeep-1
                if( (qhi > qlo .and. anchor(sub(reps(keep(i)))) <= target_arc) .or. &
                    (qhi <= qlo .and. anchor(sub(reps(keep(i)))) >= target_arc) ) best_i = i
            end do
            j = sub(reps(keep(best_i)))
            k = sub(reps(keep(min(nkeep, best_i+1))))
            w = anchor(k) - anchor(j)
            if( abs(w) > DTINY )then
                t = (target_arc - anchor(j)) / w
                t = max(0.d0, min(1.d0, t))
            else
                t = 0.d0
            endif
            do q = 1, ncomp
                centroids(q,s) = (1.d0-t)*z(j,q) + t*z(k,q)
            end do
        end do
        deallocate(keep)
        write(logfhandle,'(A,I0,A,I0,A,I0,A,F8.3)') '>>> FLEX_PCA FINCH geodesic: knn=',knn, &
            &' endpoints=',a,'/',b,'  path arc length=',arc(npath)
        call flush(logfhandle)
        ok = .true.
        call hier%kill
        deallocate(feats, sd, sub, labels, reps, rz, dist, prev, done, arc, path)
    end subroutine finch_geodesic_targets

    !> Dijkstra from src over the symmetric kNN graph of the nrep columns of rz.
    subroutine dijkstra_knn( rz, nd, nrep, knn, src, dist, prev )
        integer,  intent(in)  :: nd, nrep, knn, src
        real(dp), intent(in)  :: rz(nd,nrep)
        real(dp), intent(out) :: dist(nrep)
        integer,  intent(out) :: prev(nrep)
        real(dp), allocatable :: d2mat(:,:), w(:)
        integer,  allocatable :: nbr(:,:), ord(:)
        logical,  allocatable :: done(:)
        real(dp) :: d2, alt, dbest
        integer  :: i, j, q, k, u, ibest
        allocate(d2mat(nrep,nrep), source=0.d0)
        do i = 1, nrep
            do j = i+1, nrep
                d2 = 0.d0
                do q = 1, nd
                    d2 = d2 + (rz(q,i)-rz(q,j))**2
                end do
                d2mat(i,j) = sqrt(d2); d2mat(j,i) = d2mat(i,j)
            end do
            d2mat(i,i) = huge(1.d0)
        end do
        ! kNN adjacency (symmetrized implicitly by scanning both directions below)
        allocate(nbr(knn,nrep), w(nrep), ord(nrep))
        do i = 1, nrep
            w = d2mat(:,i)
            do k = 1, knn
                ibest = minloc(w, dim=1)
                nbr(k,i) = ibest
                w(ibest) = huge(1.d0)
            end do
        end do
        allocate(done(nrep), source=.false.)
        dist = huge(1.d0); prev = 0; dist(src) = 0.d0
        do
            ibest = 0; dbest = huge(1.d0)
            do i = 1, nrep
                if( .not.done(i) .and. dist(i) < dbest )then
                    dbest = dist(i); ibest = i
                endif
            end do
            if( ibest == 0 ) exit
            u = ibest; done(u) = .true.
            ! edges out of u, plus edges INTO u from nodes listing u as a neighbour (symmetry)
            do k = 1, knn
                j = nbr(k,u)
                alt = dist(u) + d2mat(u,j)
                if( alt < dist(j) )then
                    dist(j) = alt; prev(j) = u
                endif
            end do
            do i = 1, nrep
                if( done(i) ) cycle
                if( any(nbr(:,i) == u) )then
                    alt = dist(u) + d2mat(u,i)
                    if( alt < dist(i) )then
                        dist(i) = alt; prev(i) = u
                    endif
                endif
            end do
        end do
        deallocate(d2mat, nbr, w, ord, done)
    end subroutine dijkstra_knn

    !> Reliability proxy computable from a cached embedding alone: observed spread over mean posterior
    !! variance, mapped through r/(1+r) onto the same "higher is better measured" scale as the
    !! split-half rho. A component whose spread merely matches its own posterior width scores ~0.5;
    !! one whose spread greatly exceeds it approaches 1. The posterior variances need a matrix inverse
    !! per particle, so they are averaged over a stride sample -- only the mean enters.
    subroutine component_reliability_proxy( z, precision, nptcls, ncomp, rho )
        integer,  intent(in)  :: nptcls, ncomp
        real(dp), intent(in)  :: z(nptcls,ncomp), precision(ncomp,ncomp,nptcls)
        real(dp), intent(out) :: rho(ncomp)
        integer,  parameter   :: NSAMPLE = 2000
        real(dp), allocatable :: cfull(:,:), pvar(:)
        real(dp) :: zbar, spread, ratio
        integer  :: i, q, errflg, stride, nused
        allocate(cfull(ncomp,ncomp), pvar(ncomp), source=0.d0)
        stride = max(1, nptcls/NSAMPLE)
        nused  = 0
        do i = 1, nptcls, stride
            call matinv(precision(:,:,i), cfull, ncomp, errflg)
            if( errflg /= 0 ) cycle
            do q = 1, ncomp
                pvar(q) = pvar(q) + max(cfull(q,q), 0.d0)
            end do
            nused = nused + 1
        end do
        if( nused > 0 ) pvar = pvar / real(nused,dp)
        do q = 1, ncomp
            zbar   = sum(z(:,q)) / real(nptcls,dp)
            spread = sum((z(:,q)-zbar)**2) / real(nptcls,dp)
            ratio  = spread / max(pvar(q), DTINY)
            rho(q) = ratio / (1.d0 + ratio)
        end do
        deallocate(cfull, pvar)
    end subroutine component_reliability_proxy

    !> Manifold-covering state targets by diffusion-map k-center.
    !!
    !! This is the placement the two validation datasets jointly select, and it is deliberately
    !! agnostic: it never asks whether the heterogeneity is continuous or compositional.
    !!
    !! k-means allocates centroids in proportion to DENSITY, so it over-samples whatever is
    !! abundant and misses sparse states. A 1-D ordered path fixes coverage on a genuine reaction
    !! coordinate but merges states on a BRANCHED manifold -- measured on IgG-RL it collapsed 20
    !! requested states onto 6 distinct ground-truth conformations. A diffusion embedding turns
    !! geodesic structure on the manifold into ordinary Euclidean distance, after which greedy
    !! farthest-point (k-center) covers whatever shape the manifold actually has: a chain gives
    !! evenly spaced path targets, a branched graph spreads targets across the branches.
    !!
    !! Measured against the alternatives at equal settings:
    !!   IgG-RL (20 states)   distinct GT conformations 20/20, dihedral span 83.0%
    !!                        (k-means 20/20 and 68.2%; 1-D path 6/20 and 53.7%; RECOVAR 68.7%)
    !!   Ribosembly (16)      6/16 distinct GT volumes matched (k-means 5, 1-D path 4)
    !!
    !! Reliability weighting is essential: on both datasets z1 carries the largest eigenvalue and
    !! the WORST split-half reliability (0.24 and 0.00), so an unweighted metric follows noise.
    subroutine diffusion_kcenter_targets( z, nptcls, ncomp, nstates, wcomp, rho, centroids, ok )
        integer,  intent(in)  :: nptcls, ncomp, nstates
        real(dp), intent(in)  :: z(nptcls,ncomp), wcomp(ncomp), rho(ncomp)
        real(dp), intent(out) :: centroids(ncomp,nstates)
        logical,  intent(out) :: ok
        integer,  parameter :: NNODE_MAX = 20000   ! graph nodes; kd-tree kNN keeps this affordable
        integer,  parameter :: KNN       = 24      ! neighbours per node in the affinity graph
        integer,  parameter :: NDIFF     = 2       ! non-trivial diffusion coordinates retained
        integer,  parameter :: NPOW      = 300     ! orthogonal-iteration sweeps
        real(dp), parameter :: RHO_FLOOR = 0.1d0
        type(kd_tree)   :: tree
        type(knn_table) :: knntab
        real,     allocatable :: feats(:,:)
        integer,  allocatable :: nodes(:), er(:), ec(:), sel(:), cell(:), ccnt(:)
        real(dp), allocatable :: ev(:), sig(:), qdeg(:), ddeg(:), V(:,:), W(:,:), lam(:)
        real(dp), allocatable :: psi(:,:), dmin(:), rw(:), zbar(:), sdv(:)
        real(dp) :: d2, s, nrm, dot, best, wk1
        integer  :: nnode, i, j, q, e, nedge, it, m, c, ibest, ni, nj
        ok = .false.
        if( nstates < 2 .or. nptcls < 100 ) return
        nnode = min(nptcls, NNODE_MAX)
        if( nnode <= KNN + 2 ) return
        m = NDIFF + 1                                  ! + the trivial eigenvector
        allocate(nodes(nnode), rw(ncomp), zbar(ncomp), sdv(ncomp))
        ! deterministic stride subsample: a uniform subsample already carries the data's density,
        ! so no separate quantiser is needed and the node set is reproducible run to run
        do i = 1, nnode
            nodes(i) = 1 + int(real(i-1,dp)*real(nptcls-1,dp)/real(max(1,nnode-1),dp))
        end do
        do q = 1, ncomp
            zbar(q) = sum(z(:,q)) / real(nptcls,dp)
            sdv(q)  = sqrt(max(sum((z(:,q)-zbar(q))**2)/real(nptcls,dp), DTINY))
        end do
        s = maxval(rho)
        if( s <= DTINY ) s = 1.d0
        do q = 1, ncomp
            ! wcomp is 1/var per component, so sqrt(wcomp) ALREADY carries the 1/sd standardisation.
            ! Dividing by sdv again standardises twice and flattens the geometry the graph reads.
            rw(q) = sqrt(max(wcomp(q),0.d0)) * max(rho(q)/s, RHO_FLOOR)
        end do
        allocate(feats(ncomp,nnode))
        do i = 1, nnode
            do q = 1, ncomp
                feats(q,i) = real((z(nodes(i),q) - zbar(q)) * rw(q))
            end do
        end do
        call tree%build(feats)
        call tree%query_all(feats, KNN, knntab)
        ! self-tuning bandwidth: sigma_i is the distance to the K-th neighbour, so the affinity
        ! adapts to local sampling density instead of imposing one global scale
        allocate(sig(nnode))
        do i = 1, nnode
            sig(i) = sqrt(max(real(knntab%distance2(KNN,i),dp), DTINY))
        end do
        ! symmetric affinity W + W^T stored as an edge list; both directions are emitted
        nedge = 2*nnode*KNN
        allocate(er(nedge), ec(nedge), ev(nedge))
        e = 0
        do i = 1, nnode
            do j = 1, KNN
                ni = knntab%neighbor(j,i)
                if( ni < 1 .or. ni > nnode .or. ni == i ) cycle
                d2 = real(knntab%distance2(j,i),dp)
                s  = exp(-d2/max(sig(i)*sig(ni), DTINY))
                e = e + 1; er(e) = i;  ec(e) = ni; ev(e) = s
                e = e + 1; er(e) = ni; ec(e) = i;  ev(e) = s
            end do
        end do
        nedge = e
        if( nedge < nnode ) return
        ! alpha = 1 (Laplace-Beltrami): divide out the sampling density so the embedding reflects
        ! manifold GEOMETRY rather than how heavily each region happens to be populated. This is
        ! what stops abundant states from dominating, the exact failure mode of k-means here.
        allocate(qdeg(nnode), source=0.d0)
        do e = 1, nedge
            qdeg(er(e)) = qdeg(er(e)) + ev(e)
        end do
        do e = 1, nedge
            ev(e) = ev(e) / max(qdeg(er(e))*qdeg(ec(e)), DTINY)
        end do
        allocate(ddeg(nnode), source=0.d0)
        do e = 1, nedge
            ddeg(er(e)) = ddeg(er(e)) + ev(e)
        end do
        do i = 1, nnode
            ddeg(i) = 1.d0 / sqrt(max(ddeg(i), DTINY))
        end do
        do e = 1, nedge
            ev(e) = ev(e) * ddeg(er(e)) * ddeg(ec(e))       ! S = D^-1/2 W D^-1/2, symmetric
        end do
        ! leading eigenvectors of S by orthogonal iteration (deterministic start)
        allocate(V(nnode,m), W(nnode,m), lam(m))
        do j = 1, m
            do i = 1, nnode
                V(i,j) = sin(real(i*j,dp)*0.7717d0) + 0.1d0*real(j,dp)
            end do
        end do
        call orth_block(V, nnode, m)
        do it = 1, NPOW
            W = 0.d0
            do e = 1, nedge
                do j = 1, m
                    W(er(e),j) = W(er(e),j) + ev(e)*V(ec(e),j)
                end do
            end do
            V = W
            call orth_block(V, nnode, m)
        end do
        ! Rayleigh quotients
        W = 0.d0
        do e = 1, nedge
            do j = 1, m
                W(er(e),j) = W(er(e),j) + ev(e)*V(ec(e),j)
            end do
        end do
        do j = 1, m
            lam(j) = sum(V(:,j)*W(:,j))
        end do
        call sort_block_desc(V, lam, nnode, m)
        ! psi = D^-1/2 V, dropping the trivial leading eigenvector; commute-time scaling
        ! 1/sqrt(1-lambda) puts the coordinates on a metric where Euclidean distance approximates
        ! diffusion distance, which is what makes plain k-center meaningful in this space
        allocate(psi(nnode,NDIFF))
        wk1 = 1.d0 / sqrt(max(1.d0 - min(lam(2), 1.d0-1.d-9), 1.d-9))
        do j = 1, NDIFF
            s = 1.d0 / sqrt(max(1.d0 - min(lam(j+1), 1.d0-1.d-9), 1.d-9))
            nrm = sqrt(max(sum((V(:,j+1)*ddeg)**2), DTINY))
            do i = 1, nnode
                psi(i,j) = V(i,j+1)*ddeg(i)/nrm * (s/wk1)
            end do
        end do
        ! greedy k-center: seed at the node farthest from the embedding's centroid, then repeatedly
        ! take the node farthest from everything already chosen -- coverage, not density
        allocate(sel(nstates), dmin(nnode))
        do j = 1, NDIFF
            s = sum(psi(:,j))/real(nnode,dp)
            psi(:,j) = psi(:,j) - s
        end do
        best = -1.d0; ibest = 1
        do i = 1, nnode
            d2 = sum(psi(i,:)**2)
            if( d2 > best )then
                best = d2; ibest = i
            endif
        end do
        sel(1) = ibest
        do i = 1, nnode
            dmin(i) = sum((psi(i,:)-psi(ibest,:))**2)
        end do
        do c = 2, nstates
            best = -1.d0; ibest = 1
            do i = 1, nnode
                if( dmin(i) > best )then
                    best = dmin(i); ibest = i
                endif
            end do
            sel(c) = ibest
            do i = 1, nnode
                dmin(i) = min(dmin(i), sum((psi(i,:)-psi(ibest,:))**2))
            end do
        end do
        ! Every node inherits the CELL of the nearest selected node in the diffusion embedding, and a
        ! target is the mean RAW latent coordinate of its whole cell. k-center selects RIM nodes on
        ! purpose -- that is exactly what buys coverage -- so the selected node, and any small
        ! neighbourhood around it, is a rim sample rather than a state. Averaging the full cell
        ! regresses the target onto the state that cell actually covers. Measured on Ribosembly:
        ! rim-node target 7/16 distinct states, cell mean 11/16.
        allocate(cell(nnode), ccnt(nstates), source=0)
        do i = 1, nnode
            best = huge(1.d0); ibest = 1
            do c = 1, nstates
                d2 = sum((psi(i,:) - psi(sel(c),:))**2)
                if( d2 < best )then
                    best = d2; ibest = c
                endif
            end do
            cell(i) = ibest
        end do
        centroids = 0.d0
        do i = 1, nnode
            c = cell(i)
            ccnt(c) = ccnt(c) + 1
            centroids(:,c) = centroids(:,c) + z(nodes(i),:)
        end do
        do c = 1, nstates
            if( ccnt(c) > 0 )then
                centroids(:,c) = centroids(:,c) / real(ccnt(c),dp)
            else
                centroids(:,c) = z(nodes(sel(c)),:)     ! degenerate cell: keep the node itself
            endif
        end do
        write(logfhandle,'(A,I0,A,I0)') '>>> FLEX_PCA diffusion cell occupancy: min=',minval(ccnt), &
            &' max=',maxval(ccnt)
        write(logfhandle,'(A,I0,A,I0,A,F7.4,A,F7.4)') '>>> FLEX_PCA diffusion k-center: nodes=',nnode, &
            &' knn=',KNN,' lambda2=',lam(2),' lambda3=',lam(3)
        ok = .true.
        call tree%kill; call knntab%kill
        deallocate(feats, nodes, er, ec, ev, sig, qdeg, ddeg, V, W, lam, psi, dmin, sel, rw, zbar, sdv, cell, ccnt)
    end subroutine diffusion_kcenter_targets

    !> Modified Gram-Schmidt orthonormalisation of an (n,m) block, m small.
    subroutine orth_block( V, n, m )
        integer,  intent(in)    :: n, m
        real(dp), intent(inout) :: V(n,m)
        real(dp) :: dot, nrm
        integer  :: j, k
        do j = 1, m
            do k = 1, j-1
                dot = sum(V(:,k)*V(:,j))
                V(:,j) = V(:,j) - dot*V(:,k)
            end do
            nrm = sqrt(max(sum(V(:,j)**2), DTINY))
            V(:,j) = V(:,j) / nrm
        end do
    end subroutine orth_block

    !> Sort an eigenvector block by descending eigenvalue.
    subroutine sort_block_desc( V, lam, n, m )
        integer,  intent(in)    :: n, m
        real(dp), intent(inout) :: V(n,m), lam(m)
        real(dp) :: tl
        real(dp), allocatable :: tv(:)
        integer  :: i, j, imax
        allocate(tv(n))
        do i = 1, m-1
            imax = i
            do j = i+1, m
                if( lam(j) > lam(imax) ) imax = j
            end do
            if( imax /= i )then
                tl = lam(i); lam(i) = lam(imax); lam(imax) = tl
                tv = V(:,i); V(:,i) = V(:,imax); V(:,imax) = tv
            endif
        end do
        deallocate(tv)
    end subroutine sort_block_desc

    !> Equal-occupancy targets along a reliability-ordered principal direction.
    !!
    !! Two departures from k-means, aimed at the same failure. (1) The ordering direction is the leading
    !! eigenvector of the latent covariance weighted by per-component RELIABILITY, so a high-variance but
    !! poorly measured nuisance component cannot set the direction. nkern cannot achieve this because it
    !! slices components by INDEX and the nuisance mode is normally component 1, the largest eigenvalue.
    !! (2) Slices carry equal PARTICLE COUNTS rather than equal width, so targets spread over the
    !! populated manifold instead of concentrating wherever the density happens to be highest -- k-means
    !! allocates centroids in proportion to density, which is why it splits abundant states and never
    !! isolates sparse ones.
    !! Each target is the slice MEAN over all retained components, so the polyline of targets follows the
    !! curvature of the manifold rather than cutting a straight chord across it.
    subroutine reliability_path_targets( z, nptcls, ncomp, nstates, wcomp, rho, centroids, proj_out, tproj_out )
        integer,  intent(in)  :: nptcls, ncomp, nstates
        real(dp), intent(in)  :: z(nptcls,ncomp), wcomp(ncomp), rho(ncomp)
        real(dp), intent(out) :: centroids(ncomp,nstates)
        ! per-particle and per-target coordinate ALONG the path. States placed on a path must also be
        ! weighted along it: a full-rank Mahalanobis kernel measures the 15 off-path directions too, so a
        ! particle sitting exactly on-path but noisy elsewhere falls outside every state's support. That
        ! is what strands the majority of the dataset once the kernels are narrow enough to be distinct.
        real(dp), optional, intent(out) :: proj_out(nptcls), tproj_out(nstates)
        integer,  parameter   :: NPOWER = 128
        real(dp), parameter   :: RHO_PATH_FLOOR = 0.1d0
        real(dp), allocatable :: zbar(:), u(:), unew(:), cov(:,:), proj(:), rw(:), zc(:,:), pedge(:)
        real,     allocatable :: sproj(:)
        integer,  allocatable :: cnt(:)
        real(dp) :: rmax, nrm, d
        integer  :: i, q, s, it, ilo, islot, slo, shi
        allocate(zbar(ncomp), u(ncomp), unew(ncomp), cov(ncomp,ncomp), proj(nptcls), &
            &sproj(nptcls), cnt(nstates), rw(ncomp), zc(nptcls,ncomp), pedge(max(1,nstates-1)))
        do q = 1, ncomp
            zbar(q) = sum(z(:,q)) / real(nptcls,dp)
        end do
        ! Reliability weights, taken RELATIVE to the best-measured component and floored so no direction
        ! is removed outright. wcomp carries the same 1/var standardisation the kernel metric uses.
        rmax = maxval(rho)
        if( rmax <= DTINY ) rmax = 1.d0
        do q = 1, ncomp
            rw(q) = sqrt(max(wcomp(q),0.d0)) * max(rho(q)/rmax, RHO_PATH_FLOOR)
        end do
        !$omp parallel do default(shared) private(i,q) schedule(static)
        do i = 1, nptcls
            do q = 1, ncomp
                zc(i,q) = rw(q)*(z(i,q)-zbar(q))
            end do
        end do
        !$omp end parallel do
        cov = matmul(transpose(zc), zc)
        ! leading eigenvector by power iteration: deterministic, and ncomp is at most a few dozen
        u = 1.d0 / sqrt(real(ncomp,dp))
        do it = 1, NPOWER
            unew = matmul(cov, u)
            nrm  = sqrt(sum(unew*unew))
            if( nrm <= DTINY ) exit
            unew = unew / nrm
            if( sum(abs(unew-u)) < 1.d-12 )then
                u = unew
                exit
            endif
            u = unew
        end do
        !$omp parallel do default(shared) private(i,q,d) schedule(static)
        do i = 1, nptcls
            d = 0.d0
            do q = 1, ncomp
                d = d + zc(i,q)*u(q)
            end do
            proj(i) = d
        end do
        !$omp end parallel do
        ! EQUAL-OCCUPANCY edges: every slice carries the same particle count by construction
        sproj = real(proj)
        call hpsort(sproj)
        do s = 1, nstates-1
            ilo = max(1, min(nptcls, nint(real(s,dp)/real(nstates,dp)*real(nptcls,dp))))
            pedge(s) = real(sproj(ilo),dp)
        end do
        centroids = 0.d0; cnt = 0
        do i = 1, nptcls
            islot = nstates
            do s = 1, nstates-1
                if( proj(i) < pedge(s) )then
                    islot = s
                    exit
                endif
            end do
            cnt(islot) = cnt(islot) + 1
            do q = 1, ncomp
                centroids(q,islot) = centroids(q,islot) + z(i,q)
            end do
        end do
        do s = 1, nstates
            if( cnt(s) > 0 )then
                centroids(:,s) = centroids(:,s) / real(cnt(s),dp)
            endif
        end do
        ! Ties in the projection can leave a slice empty. Interpolate it between its nearest occupied
        ! neighbours so the target polyline stays ordered and no state collapses onto the global mean.
        do s = 1, nstates
            if( cnt(s) > 0 ) cycle
            slo = 0; shi = 0
            do i = s-1, 1, -1
                if( cnt(i) > 0 )then
                    slo = i; exit
                endif
            end do
            do i = s+1, nstates
                if( cnt(i) > 0 )then
                    shi = i; exit
                endif
            end do
            if( slo > 0 .and. shi > 0 )then
                d = real(s-slo,dp)/real(shi-slo,dp)
                centroids(:,s) = (1.d0-d)*centroids(:,slo) + d*centroids(:,shi)
            else if( slo > 0 )then
                centroids(:,s) = centroids(:,slo)
            else if( shi > 0 )then
                centroids(:,s) = centroids(:,shi)
            else
                centroids(:,s) = zbar
            endif
        end do
        if( present(proj_out) ) proj_out = proj
        if( present(tproj_out) )then
            ! each target's own path coordinate = the mean projection of the particles in its slice
            tproj_out = 0.d0
            do i = 1, nptcls
                islot = nstates
                do s = 1, nstates-1
                    if( proj(i) < pedge(s) )then
                        islot = s
                        exit
                    endif
                end do
                tproj_out(islot) = tproj_out(islot) + proj(i)
            end do
            do s = 1, nstates
                if( cnt(s) > 0 )then
                    tproj_out(s) = tproj_out(s) / real(cnt(s),dp)
                else if( s > 1 )then
                    tproj_out(s) = tproj_out(s-1)
                endif
            end do
        endif
        write(logfhandle,'(A,I0,A,I0,A,F6.3)') '>>> FLEX_PCA path ordering direction over ',ncomp, &
            &' components; slice occupancies min=',minval(cnt),'  leading |u| on z1=',abs(u(1))
        deallocate(zbar, u, unew, cov, proj, sproj, cnt, rw, zc, pedge)
    end subroutine reliability_path_targets

    subroutine path_latent_targets( z, nptcls, ncomp, nstates, wcomp, centroids )
        integer,  intent(in)  :: nptcls, ncomp, nstates
        real(dp), intent(in)  :: z(nptcls,ncomp), wcomp(ncomp)
        real(dp), intent(out) :: centroids(ncomp,nstates)
        real(dp), allocatable :: zbar(:), pa(:), pb(:), dirv(:), proj(:)
        real,     allocatable :: sproj(:)
        integer,  allocatable :: cnt(:)
        real(dp) :: d2, best, dmax, dnorm, lo, hi, t
        integer  :: i, q, s, ia, ib, islot
        allocate(zbar(ncomp), pa(ncomp), pb(ncomp), dirv(ncomp), proj(nptcls), &
            &sproj(nptcls), cnt(nstates))
        ! endpoint 1: farthest particle from the latent mean
        do q = 1, ncomp
            zbar(q) = sum(z(:,q)) / real(nptcls,dp)
        end do
        dmax = -1.d0; ia = 1
        do i = 1, nptcls
            d2 = 0.d0
            do q = 1, ncomp
                d2 = d2 + wcomp(q)*(z(i,q)-zbar(q))**2
            end do
            if( d2 > dmax )then
                dmax = d2; ia = i
            endif
        end do
        pa = z(ia,:)
        ! endpoint 2: farthest particle from endpoint 1
        dmax = -1.d0; ib = 1
        do i = 1, nptcls
            d2 = 0.d0
            do q = 1, ncomp
                d2 = d2 + wcomp(q)*(z(i,q)-pa(q))**2
            end do
            if( d2 > dmax )then
                dmax = d2; ib = i
            endif
        end do
        pb = z(ib,:)
        dirv = pb - pa
        dnorm = sum(wcomp*dirv*dirv)
        if( dnorm <= DTINY )then
            do s = 1, nstates
                centroids(:,s) = zbar
            end do
            deallocate(zbar, pa, pb, dirv, proj, sproj, cnt)
            return
        endif
        ! project every particle onto the segment, in the same metric the kernel uses
        !$omp parallel do default(shared) private(i,q,d2) schedule(static)
        do i = 1, nptcls
            d2 = 0.d0
            do q = 1, ncomp
                d2 = d2 + wcomp(q)*(z(i,q)-pa(q))*dirv(q)
            end do
            proj(i) = d2 / dnorm
        end do
        !$omp end parallel do
        ! equal-WIDTH slices between the 0.1 % and 99.9 % projections, each target the slice mean
        sproj = real(proj)
        call hpsort(sproj)
        centroids = 0.d0; cnt = 0
        do i = 1, nptcls
            lo = real(sproj(max(1,min(nptcls,1+int(real(nptcls,dp)*0.001d0)))),dp)
            hi = real(sproj(max(1,min(nptcls,nptcls-int(real(nptcls,dp)*0.001d0)))),dp)
            if( hi-lo <= DTINY )then
                islot = 1
            else
                t     = (proj(i)-lo)/(hi-lo)
                islot = 1 + int(t*real(nstates,dp))
                islot = max(1, min(nstates, islot))
            endif
            cnt(islot) = cnt(islot) + 1
            do q = 1, ncomp
                centroids(q,islot) = centroids(q,islot) + z(i,q)
            end do
        end do
        do s = 1, nstates
            if( cnt(s) > 0 )then
                centroids(:,s) = centroids(:,s) / real(cnt(s),dp)
            else
                ! empty slice: interpolate along the segment so the path stays ordered
                t = real(s-1,dp)/real(max(1,nstates-1),dp)
                centroids(:,s) = pa + t*dirv
            endif
        end do
        write(logfhandle,'(A,I0,A,I0,A,I0)') '>>> FLEX_PCA path endpoints from particles ',ia,' and ',ib, &
            &'; slice occupancies min=',minval(cnt)
        deallocate(zbar, pa, pb, dirv, proj, sproj, cnt)
    end subroutine path_latent_targets

    !> Deterministic k-means over the retained latent coordinates, in the SAME per-component metric wcomp
    !! that the Epanechnikov kernel uses, so target placement and particle weighting agree. Seeded
    !! farthest-point from the particle nearest the latent mean, so the result does not depend on an RNG.
    subroutine kmeans_latent_targets( z, nptcls, ncomp, nstates, wcomp, centroids )
        integer,  intent(in)  :: nptcls, ncomp, nstates
        real(dp), intent(in)  :: z(nptcls,ncomp), wcomp(ncomp)
        real(dp), intent(out) :: centroids(ncomp,nstates)
        integer, parameter    :: MAXIT = 50
        real(dp), allocatable :: mind(:), csum(:,:), zbar(:)
        integer,  allocatable :: cnt(:), memb(:)
        real(dp) :: d2, best, dmax
        integer  :: i, q, s, it, ibest, iseed, nchanged
        logical  :: l_reseed
        allocate(mind(nptcls), csum(ncomp,nstates), zbar(ncomp), cnt(nstates), memb(nptcls))
        ! seed 1: the particle closest to the latent mean
        do q = 1, ncomp
            zbar(q) = sum(z(:,q)) / real(nptcls,dp)
        end do
        best = huge(1.d0); iseed = 1
        do i = 1, nptcls
            d2 = 0.d0
            do q = 1, ncomp
                d2 = d2 + wcomp(q)*(z(i,q)-zbar(q))**2
            end do
            if( d2 < best )then
                best  = d2
                iseed = i
            endif
        end do
        centroids(:,1) = z(iseed,:)
        ! seeds 2..nstates: farthest point from the already-chosen set
        mind = huge(1.d0)
        do s = 2, nstates
            dmax = -1.d0; iseed = 1
            do i = 1, nptcls
                d2 = 0.d0
                do q = 1, ncomp
                    d2 = d2 + wcomp(q)*(z(i,q)-centroids(q,s-1))**2
                end do
                mind(i) = min(mind(i), d2)
                if( mind(i) > dmax )then
                    dmax  = mind(i)
                    iseed = i
                endif
            end do
            centroids(:,s) = z(iseed,:)
        end do
        ! Lloyd
        memb = 0
        do it = 1, MAXIT
            nchanged = 0
            !$omp parallel do default(shared) private(i,q,s,d2,best,ibest) schedule(static) reduction(+:nchanged)
            do i = 1, nptcls
                best = huge(1.d0); ibest = 1
                do s = 1, nstates
                    d2 = 0.d0
                    do q = 1, ncomp
                        d2 = d2 + wcomp(q)*(z(i,q)-centroids(q,s))**2
                    end do
                    if( d2 < best )then
                        best  = d2
                        ibest = s
                    endif
                end do
                if( memb(i) /= ibest ) nchanged = nchanged + 1
                memb(i) = ibest
                mind(i)   = best
            end do
            !$omp end parallel do
            csum = 0.d0; cnt = 0
            do i = 1, nptcls
                cnt(memb(i))    = cnt(memb(i)) + 1
                csum(:,memb(i)) = csum(:,memb(i)) + z(i,:)
            end do
            l_reseed = .false.
            do s = 1, nstates
                if( cnt(s) > 0 )then
                    centroids(:,s) = csum(:,s) / real(cnt(s),dp)
                else
                    ! empty cluster: reseed on the worst-fitted particle
                    iseed          = maxloc(mind, dim=1)
                    centroids(:,s) = z(iseed,:)
                    mind(iseed)    = -1.d0
                    l_reseed       = .true.
                endif
            end do
            if( nchanged == 0 .and. .not. l_reseed ) exit
        end do
        write(logfhandle,'(A,I0,A,I0)') '>>> FLEX_PCA k-means iterations=',min(it,MAXIT), &
            &'  reassigned on last pass=',nchanged
        deallocate(mind, csum, zbar, cnt, memb)
    end subroutine kmeans_latent_targets

    !> Order the reconstructed state volumes along the dominant conformational motion.
    subroutine order_states_by_volume_trajectory( params, nstates, ncomp, eigdir, axis_sel, comp_rank )
        class(parameters), intent(in) :: params
        integer,           intent(in) :: nstates, ncomp
        ! directory holding flex_pca_pc###.mrc.
        character(len=*), optional, intent(in) :: eigdir
        ! latent component carrying the reproducible motion; 0 if none could be selected
        integer,          optional, intent(out) :: axis_sel
        ! components ordered by half-set reproducibility, most reproducible first.
        integer,          optional, intent(out) :: comp_rank(ncomp)
        character(len=:), allocatable :: pcdir
        type(image),  allocatable :: vols(:), eigs(:), evols(:), ovols(:)
        real(dp),     allocatable :: X(:,:), vbar(:), coord(:), coord_best(:), U(:), diff(:)
        real(dp),     allocatable :: Xe(:,:), Xo(:,:), deven(:), dodd(:)
        integer,      allocatable :: order(:), order_best(:), mskidx(:)
        real,         allocatable :: rmat(:,:,:), vsort(:)
        type(string) :: fn
        integer  :: s, k, nvox, ldim(3), iu, kbest, ipair, jpair, nmsk, imsk, ldim_pc(3), ifoo
        real(dp) :: sc, sc_best, sep, cons, sep_best, cons_best, sep_pair
        real(dp) :: rep, rep_best, thresh, amp, amp_best, score, score_best, mean_msk_rms
        real(dp), allocatable :: rep_all(:)
        integer,  allocatable :: rep_ord(:)
        logical  :: l_halves
        if( nstates < 2 ) return
        if( present(axis_sel) ) axis_sel = 0
        pcdir = ''
        if( present(eigdir) ) pcdir = eigdir
        ! The eigenvolumes are written at the box_crop of the run that produced them.
        fn = pcdir//'flex_pca_pc001'//MRC_EXT
        call find_ldim_nptcls(fn, ldim_pc, ifoo)
        call fn%kill
        if( ldim_pc(1) /= params%box_crop )then
            write(logfhandle,'(A,I0,A,I0,A)') '>>> FLEX_PCA skipping volume-space trajectory ordering: &
                &cached eigenvolumes are box ',ldim_pc(1),' but this run reconstructs at box ',params%box_crop, &
                &'. The state maps are unaffected; re-run without infile= to order the trajectory'
            call flush(logfhandle)
            return
        endif
        ! Read the combined state volumes AT box_crop.
        allocate(vols(nstates))
        do s = 1, nstates
            fn = 'flex_pca_state_'//int2str_pad(s,3)//MRC_EXT
            call vols(s)%read_and_crop(fn, flex_rec_smpd(params), params%box_crop, params%smpd_crop)
            call fn%kill
        end do
        rmat = vols(1)%get_rmat(); ldim = shape(rmat); nvox = ldim(1)*ldim(2)*ldim(3)
        allocate(X(nvox,nstates), vbar(nvox), source=0.d0)
        do s = 1, nstates
            X(:,s) = reshape(real(vols(s)%get_rmat(),dp), [nvox])
        end do
        do k = 1, nvox
            vbar(k) = sum(X(k,:)) / real(nstates,dp)
        end do
        do s = 1, nstates
            X(:,s) = X(:,s) - vbar
        end do
        ! ---- MOLECULE MASK (top 8 % of the state mean) ---- The axis score below is a half-set
        ! reproducibility and MUST be confined to the molecule: the solvent baseline and the global level
        ! of a state map reproduce trivially between halves and would dominate an unmasked correlation.
        allocate(vsort(nvox))
        vsort = real(vbar)
        call hpsort(vsort)
        thresh = real(vsort(max(1, min(nvox, nint(0.92*real(nvox))))), dp)
        nmsk = count(vbar > thresh)
        if( nmsk < 64 )then   ! degenerate mean: fall back to the whole box
            nmsk = nvox
            allocate(mskidx(nmsk))
            do k = 1, nvox
                mskidx(k) = k
            end do
        else
            allocate(mskidx(nmsk))
            imsk = 0
            do k = 1, nvox
                if( vbar(k) > thresh )then
                    imsk = imsk + 1; mskidx(imsk) = k
                endif
            end do
        endif
        deallocate(vsort)
        ! ---- HALF MAPS ---- Reproducibility needs the disjoint-particle reconstructions.
        l_halves = file_exists('flex_pca_even_state_'//int2str_pad(1,3)//MRC_EXT) .and. &
            &      file_exists('flex_pca_odd_state_' //int2str_pad(1,3)//MRC_EXT)
        if( l_halves )then
            allocate(evols(nstates), ovols(nstates))
            allocate(Xe(nvox,nstates), Xo(nvox,nstates), source=0.d0)
            allocate(deven(nmsk), dodd(nmsk), source=0.d0)
            do s = 1, nstates
                fn = 'flex_pca_even_state_'//int2str_pad(s,3)//MRC_EXT
                call evols(s)%read_and_crop(fn, flex_rec_smpd(params), params%box_crop, params%smpd_crop)
                call fn%kill
                Xe(:,s) = reshape(real(evols(s)%get_rmat(),dp), [nvox])
                call evols(s)%kill
                fn = 'flex_pca_odd_state_'//int2str_pad(s,3)//MRC_EXT
                call ovols(s)%read_and_crop(fn, flex_rec_smpd(params), params%box_crop, params%smpd_crop)
                call fn%kill
                Xo(:,s) = reshape(real(ovols(s)%get_rmat(),dp), [nvox])
                call ovols(s)%kill
            end do
            deallocate(evols, ovols)
        endif
        ! read the covariance eigenvolumes (candidate conformational directions)
        allocate(eigs(ncomp), U(nvox), diff(nvox), coord(nstates), coord_best(nstates), &
            &order(nstates), order_best(nstates))
        do s = 1, nstates
            order_best(s) = s
        end do
        ! AXIS SCORE = half-set reproducibility x endpoint amplitude. Reproducibility alone is not enough:
        ! an axis can order the states coherently and reproduce well while its two extreme states are
        ! nearly the same map, i.e. a collapsed trajectory.
        kbest = 0; sc_best = -1.d0; sep_best = 0.d0; cons_best = 0.d0; rep_best = -2.d0
        amp_best = 0.d0; score_best = -2.d0
        ! RMS of the state mean on the molecule, so the amplitude below is a relative change
        mean_msk_rms = sqrt(sum(U(mskidx)**2)/real(nmsk,dp))
        allocate(rep_all(ncomp), source=-2.d0)
        do k = 1, ncomp
            fn = pcdir//'flex_pca_pc'//int2str_pad(k,3)//MRC_EXT
            call eigs(k)%new([params%box_crop,params%box_crop,params%box_crop], params%smpd_crop)
            call eigs(k)%read(fn); call fn%kill
            U = reshape(real(eigs(k)%get_rmat(),dp), [nvox]); U = U - sum(U)/real(nvox,dp)
            do s = 1, nstates
                coord(s) = sum(X(:,s)*U)                     ! project each state onto U_k
            end do
            call argsort_ascending(coord, order, nstates)
            diff = X(:,order(nstates)) - X(:,order(1))       ! trajectory endpoint difference
            cons = abs(pearson_dp(diff, U, nvox))            ! self-consistency
            sep  = max(0.d0, 1.d0 - pearson_dp(X(:,order(nstates)), X(:,order(1)), nvox))
            if( l_halves )then
                ! HALF-SET REPRODUCIBILITY of the endpoint difference, on the molecule.
                deven = Xe(mskidx,order(nstates)) - Xe(mskidx,order(1))
                dodd  = Xo(mskidx,order(nstates)) - Xo(mskidx,order(1))
                rep   = pearson_dp(deven, dodd, nmsk)
            else
                rep = sep * cons                             ! legacy fallback, no half maps
            endif
            sc = sep * cons
            ! ENDPOINT AMPLITUDE on the molecule, relative to the state mean
            amp = sqrt(sum(diff(mskidx)**2)/real(nmsk,dp)) / max(mean_msk_rms, DTINY)
            score = rep*amp
            rep_all(k) = rep
            if( score > score_best )then
                score_best = score; rep_best = rep; sc_best = sc; kbest = k
                coord_best = coord; order_best = order
                sep_best = sep; cons_best = cons; amp_best = amp
            endif
            write(logfhandle,'(A,I3,A,F7.4,A,F7.4,A,F8.4,A,F7.4)') '>>>   axis pc',k, &
                &'  repro=',rep,'  amp=',amp,'  score=',score,'  self-consistency=',cons
            call eigs(k)%kill
        end do
        if( present(axis_sel) .and. l_halves ) axis_sel = kbest
        ! SIMPLE_COV_TRAJAXIS forces the trajectory axis. The automatic pick answers "which component is
        ! most reproducible", which is not the same question as "which physical motion": a dominant
        ! density mode can reproduce better than the motion of interest and be near-orthogonal to it.
        if( present(axis_sel) )then
            call cov_env_int_pub('SIMPLE_COV_TRAJAXIS', axis_sel)
            if( axis_sel > ncomp ) axis_sel = ncomp
        endif
        if( present(comp_rank) )then
            allocate(rep_ord(ncomp))
            call argsort_ascending(-rep_all, rep_ord, ncomp)   ! descending reproducibility
            comp_rank = rep_ord
            deallocate(rep_ord)
        endif
        ! most separated state pair anywhere in the set, reported alongside the chosen endpoints
        sep_pair = -1.d0; ipair = 1; jpair = 1
        do s = 1, nstates - 1
            do k = s + 1, nstates
                sc = 1.d0 - pearson_dp(X(:,s), X(:,k), nvox)
                if( sc > sep_pair )then
                    sep_pair = sc; ipair = s; jpair = k
                endif
            end do
        end do
        ! write trajectory-ordered volumes and the mapping
        do s = 1, nstates
            fn = 'flex_pca_traj_'//int2str_pad(s,3)//MRC_EXT
            call vols(order_best(s))%write(fn, del_if_exists=.true.); call fn%kill
        end do
        call del_file('flex_pca_trajectory_order.txt')
        open(newunit=iu,file='flex_pca_trajectory_order.txt',status='replace',action='write')
        write(iu,'(A,I0,A,F8.4,A,F8.4,A,F8.4,A,F8.4,A)') '# volume-space trajectory: ordering eigenvolume pc',kbest, &
            &' (half-set reproducibility=',rep_best,'; separation ',sep_best,' x self-consistency ',cons_best, &
            &' = ',sc_best,')'
        write(iu,'(A,I0,A,I0,A,F8.4)') '# most separated pair anywhere in the set: states ',ipair,' and ',jpair, &
            &'  separation=',sep_pair
        write(iu,'(A)') '# traj_position  original_state  projection_coord'
        do s = 1, nstates
            write(iu,'(I5,1X,I5,1X,ES16.8)') s, order_best(s), coord_best(order_best(s))
        end do
        close(iu)
        write(logfhandle,'(A,I0,A,I0,A,F7.4,A,F7.4,A,F7.4,A)') '>>> FLEX_PCA volume-space state trajectory: ',nstates, &
            &' states ordered by eigenvolume pc',kbest,' (half-set reproducibility=',rep_best,' separation=',sep_best, &
            &' self-consistency=',cons_best,') -> flex_pca_traj_###.mrc'
        write(logfhandle,'(A,F7.4,A,I0,A,I0,A,F7.4,A)') '>>> FLEX_PCA trajectory endpoint separation=',sep_best, &
            &'; best separated pair anywhere is states ',ipair,'/',jpair,' at ',sep_pair, &
            &' (a large shortfall means the ordering axis is not spanning the motion)'
        call flush(logfhandle)
        do s = 1, nstates
            call vols(s)%kill
        end do
        deallocate(vols, eigs, X, vbar, U, diff, coord, coord_best, order, order_best, mskidx, rep_all)
        if( allocated(Xe) ) deallocate(Xe, Xo, deven, dodd)
    end subroutine order_states_by_volume_trajectory

    !>  Ascending argsort (selection sort; n is small) returning the index permutation.
    subroutine argsort_ascending( vals, order, n )
        integer,  intent(in)  :: n
        real(dp), intent(in)  :: vals(n)
        integer,  intent(out) :: order(n)
        integer :: i, j, best, tmp
        do i = 1, n
            order(i) = i
        end do
        do i = 1, n-1
            best = i
            do j = i+1, n
                if( vals(order(j)) < vals(order(best)) ) best = j
            end do
            tmp = order(i); order(i) = order(best); order(best) = tmp
        end do
    end subroutine argsort_ascending

    !>  Pearson correlation of two length-n double vectors.
    real(dp) function pearson_dp( a, b, n )
        integer,  intent(in) :: n
        real(dp), intent(in) :: a(n), b(n)
        real(dp) :: ma, mb, sa, sb, sab
        ma = sum(a)/real(n,dp); mb = sum(b)/real(n,dp)
        sa = sum((a-ma)**2); sb = sum((b-mb)**2); sab = sum((a-ma)*(b-mb))
        if( sa > 0.d0 .and. sb > 0.d0 )then
            pearson_dp = sab / sqrt(sa*sb)
        else
            pearson_dp = 0.d0
        endif
    end function pearson_dp

    !> Held-out (cross-halfset) embedding.
    subroutine heldout_embedding( params, build, mean_rec, pinds, nptcls, ncols_req, col_sep, neigs_req, &
        &basis_recs, eigvals, ncomp, sig2_eff, z, contrast, latent_second, resid_energy, resid_mean_energy )
        type(parameters),    intent(inout) :: params
        type(builder),       intent(inout) :: build
        type(reconstructor), intent(inout) :: mean_rec
        integer,             intent(in)    :: pinds(:), nptcls, ncols_req, col_sep, neigs_req
        type(reconstructor), allocatable, intent(out) :: basis_recs(:)
        real(dp),            allocatable, intent(out) :: eigvals(:)
        integer,             intent(out)   :: ncomp
        real(dp),            intent(out)   :: sig2_eff
        real(dp),            allocatable, intent(out) :: z(:,:), contrast(:), latent_second(:,:,:)
        real(dp),            allocatable, intent(out) :: resid_energy(:), resid_mean_energy(:)
        type(reconstructor), allocatable :: recs_a(:), recs_b(:)
        type(image),         allocatable :: imgs_a(:), imgs_b(:)
        real(dp),            allocatable :: eig_a(:), eig_b(:), M(:,:), svals(:)
        real(dp),            allocatable :: z_a(:,:), z_b(:,:), c_a(:), c_b(:)
        real(dp),            allocatable :: p_a(:,:,:), p_b(:,:,:), re_a(:), re_b(:), rm_a(:), rm_b(:)
        integer,             allocatable :: pind_a(:), pind_b(:), row_a(:), row_b(:)
        integer  :: i, q, na, nb, ncomp_a, ncomp_b, ia, ib
        real(dp) :: sig2_a, sig2_b
        ! partition the sampled particles by halfset
        na = 0; nb = 0
        do i = 1, nptcls
            if( build%spproj_field%get_eo(pinds(i)) == 0 )then
                na = na + 1
            else
                nb = nb + 1
            endif
        end do
        if( na < 50 .or. nb < 50 ) THROW_HARD('flex_pca heldout=yes requires >=50 particles per halfset')
        allocate(pind_a(na), pind_b(nb), row_a(na), row_b(nb))
        ia = 0; ib = 0
        do i = 1, nptcls
            if( build%spproj_field%get_eo(pinds(i)) == 0 )then
                ia = ia + 1; pind_a(ia) = pinds(i); row_a(ia) = i
            else
                ib = ib + 1; pind_b(ib) = pinds(i); row_b(ib) = i
            endif
        end do
        write(logfhandle,'(A,I0,A,I0)') '>>> FLEX_PCA HELD-OUT EMBEDDING halfset_A=',na,' halfset_B=',nb
        call flush(logfhandle)
        ! basis A from halfset A (reference frame), basis B from halfset B
        call build_covariance_eigenbasis(params, build, mean_rec, pind_a, na, &
            &ncols_req, col_sep, neigs_req, recs_a, eig_a, ncomp_a, sig2_a, &
            &basis_imgs=imgs_a, fprefix='flex_pca_pc')
        call build_covariance_eigenbasis(params, build, mean_rec, pind_b, nb, &
            &ncols_req, col_sep, neigs_req, recs_b, eig_b, ncomp_b, sig2_b, &
            &basis_imgs=imgs_b, fprefix='flex_pca_heldout_pcB')
        ! the reference basis (A) defines the output frame and rank
        ncomp = ncomp_a
        if( ncomp < 1 .or. ncomp_b < 1 ) THROW_HARD('flex_pca heldout produced no retained components')
        write(logfhandle,'(A,I0,A,I0)') '>>> FLEX_PCA heldout components A(reference)=',ncomp_a, &
            &' B=',ncomp_b
        ! cross-halfset subspace agreement (principal-angle cosines)
        call align_basis_to_reference(imgs_a, ncomp_a, imgs_b, ncomp_b, M, svals)
        write(logfhandle,'(A)') '>>> FLEX_PCA halfset subspace principal-angle cosines (independent bases):'
        do q = 1, size(svals)
            write(logfhandle,'(A,I3,A,F8.4)') '>>>   k=',q,'  cos(theta)=',real(svals(q))
        end do
        write(logfhandle,'(A,F8.4)') '>>> FLEX_PCA mean principal-angle cosine=', &
            &real(sum(svals)/real(max(1,size(svals)),dp))
        write(logfhandle,'(A)') '>>> FLEX_PCA (these are fitted on DISJOINT particles: unlike per-component'
        write(logfhandle,'(A)') '>>> FLEX_PCA  FSC they cannot be satisfied by a shared, self-confirming basis)'
        call flush(logfhandle)
        ! embed each halfset with the OTHER halfset's basis
        allocate(z_b(nb,ncomp_a), c_b(nb), p_b(ncomp_a,ncomp_a,nb), re_b(nb), rm_b(nb))
        call embed_latents_with_contrast(params, build, mean_rec, recs_a, ncomp_a, eig_a(:ncomp_a), sig2_a, &
            &pind_b, nb, z_b, c_b, p_b, re_b, rm_b)
        allocate(z_a(na,ncomp_b), c_a(na), p_a(ncomp_b,ncomp_b,na), re_a(na), rm_a(na))
        call embed_latents_with_contrast(params, build, mean_rec, recs_b, ncomp_b, eig_b(:ncomp_b), sig2_b, &
            &pind_a, na, z_a, c_a, p_a, re_a, rm_a)
        ! merge into the caller's particle order, rotating the B-frame latents (halfset A,
        ! embedded with basis B) into the reference frame of basis A
        allocate(z(nptcls,ncomp), contrast(nptcls), latent_second(ncomp,ncomp,nptcls))
        allocate(resid_energy(nptcls), resid_mean_energy(nptcls))
        z = 0.d0; contrast = 1.d0; latent_second = 0.d0; resid_energy = 0.d0; resid_mean_energy = 0.d0
        do i = 1, nb
            z(row_b(i),:)              = z_b(i,:)
            contrast(row_b(i))         = c_b(i)
            latent_second(:,:,row_b(i))= p_b(:,:,i)
            resid_energy(row_b(i))     = re_b(i)
            resid_mean_energy(row_b(i))= rm_b(i)
        end do
        do i = 1, na
            do q = 1, ncomp
                z(row_a(i),q) = sum(M(q,1:ncomp_b)*z_a(i,1:ncomp_b))
            end do
            contrast(row_a(i))          = c_a(i)
            ! the halfset-A precision was fitted in the B frame and has rank ncomp_b
            resid_energy(row_a(i))      = re_a(i)
            resid_mean_energy(row_a(i)) = rm_a(i)
        end do
        ! the reference basis (A) is handed straight over -- never copy a reconstructor that
        ! has been expand_exp'd; move the allocation instead
        allocate(eigvals(ncomp), source=eig_a(:ncomp))
        sig2_eff = sig2_a
        call move_alloc(recs_a, basis_recs)
        do q = 1, ncomp_a
            call imgs_a(q)%kill
        end do
        do q = 1, ncomp_b
            call recs_b(q)%dealloc_rho; call recs_b(q)%kill; call imgs_b(q)%kill
        end do
        deallocate(recs_b, imgs_a, imgs_b, eig_a, eig_b, M, svals)
        deallocate(z_a, z_b, c_a, c_b, p_a, p_b, re_a, re_b, rm_a, rm_b)
        deallocate(pind_a, pind_b, row_a, row_b)
    end subroutine heldout_embedding

    !> Environment-driven wrapper around probe_external_basis.
    subroutine run_external_basis_probe( params, build, mean_rec, pinds, nptcls, ncomp, eigvals, &
        &sig2_eff, l_need_mean )
        class(parameters),   intent(inout) :: params
        type(builder),       intent(inout) :: build
        type(reconstructor), intent(inout) :: mean_rec
        integer,             intent(in)    :: pinds(:), nptcls, ncomp
        real(dp),            intent(in)    :: eigvals(:), sig2_eff
        logical,             intent(in)    :: l_need_mean
        character(len=STDLEN) :: prefix, eigdir, envval
        integer :: nprobe, ln, stat, ival
        call get_environment_variable('SIMPLE_COV_PROBEBASIS', prefix, ln, stat)
        if( stat /= 0 .or. ln < 1 ) return
        nprobe = 0
        call get_environment_variable('SIMPLE_COV_PROBEN', envval, ln, stat)
        if( stat == 0 .and. ln > 0 )then
            read(envval(:ln),*,iostat=stat) ival
            if( stat == 0 ) nprobe = ival
        endif
        if( nprobe < 0 )then
            write(logfhandle,'(A)') '>>> FLEX_PCA probe basis requested but SIMPLE_COV_PROBEN invalid'
            return
        endif
        eigdir = ''
        call get_environment_variable('SIMPLE_COV_PROBEEIG', envval, ln, stat)
        if( stat == 0 .and. ln > 0 ) eigdir = envval(:ln)
        if( l_need_mean ) call estimate_covariance_mean(params, build, mean_rec, pinds, nptcls)
        call probe_external_basis(params, build, mean_rec, pinds, nptcls, trim(eigdir), ncomp, &
            &eigvals, sig2_eff, trim(prefix), nprobe)
    end subroutine run_external_basis_probe

    !> Read externally placed latent targets from the file named by SIMPLE_COV_TARGETFILE.
    subroutine read_external_targets( ncomp, nstates, targets, ok )
        integer,  intent(in)  :: ncomp, nstates
        real(dp), allocatable, intent(out) :: targets(:,:)
        logical,  intent(out) :: ok
        character(len=STDLEN)  :: fname
        character(len=LONGSTRLEN) :: line
        integer :: u, s, q, ln, stat, nread
        ok = .false.
        call get_environment_variable('SIMPLE_COV_TARGETFILE', fname, ln, stat)
        if( stat /= 0 .or. ln < 1 ) return
        if( .not. file_exists(trim(fname)) )then
            write(logfhandle,'(A,A)') '>>> FLEX_PCA SIMPLE_COV_TARGETFILE not found: ', trim(fname)
            return
        endif
        allocate(targets(ncomp,nstates), source=0._dp)
        open(newunit=u, file=trim(fname), status='old', action='read', iostat=stat)
        if( stat /= 0 )then
            deallocate(targets); return
        endif
        nread = 0
        do
            read(u,'(A)',iostat=stat) line
            if( stat /= 0 ) exit
            if( len_trim(line) == 0 ) cycle
            if( index(adjustl(line), '#') == 1 ) cycle
            nread = nread + 1
            if( nread > nstates ) exit
            read(line,*,iostat=stat) (targets(q,nread), q=1,ncomp)
            if( stat /= 0 )then
                write(logfhandle,'(A,I0)') '>>> FLEX_PCA malformed target row ', nread
                close(u); deallocate(targets); return
            endif
        end do
        close(u)
        if( nread /= nstates )then
            write(logfhandle,'(A,I0,A,I0)') '>>> FLEX_PCA target file has ',nread, &
                &' rows, expected ',nstates
            deallocate(targets); return
        endif
        ok = .true.
    end subroutine read_external_targets

    subroutine mask_state_weights_by_half( build, pinds, wanted_eo, weights )
        type(builder), intent(inout) :: build
        integer,       intent(in)    :: pinds(:), wanted_eo
        real,          intent(inout) :: weights(:,:)
        integer :: i
        do i = 1, size(pinds)
            if( build%spproj_field%get_eo(pinds(i)) /= wanted_eo ) weights(i,:) = 0.
        end do
    end subroutine mask_state_weights_by_half

    subroutine write_covariance_tables( build, pinds, z, eigvals, prior_precision, weights, labels, &
        &targets, bandwidths, neff, resid_energy, resid_mean_energy )
        type(builder), intent(inout) :: build
        integer,       intent(in) :: pinds(:), labels(:)
        real(dp),      intent(in) :: z(:,:), eigvals(:), prior_precision(:)
        real,          intent(in) :: weights(:,:), targets(:,:), bandwidths(:), neff(:)
        real(dp),      intent(in) :: resid_energy(:), resid_mean_energy(:)
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
        call del_file('flex_pca_state_targets.txt')
        open(newunit=u,file='flex_pca_state_targets.txt',status='replace',action='write')
        ! targets are full latent-space points (k-means centroids by default), so every coordinate is
        ! written
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

    !>  Write the hard state assignment as a runnable project: a copy of the input project in which
    !!  ptcl3D/state carries the state label of every embedded particle and 0 everywhere else.
    !!  This decouples the embedding and its state assignment from the kernel-weighted reconstruction
    !!  backend, so the clusters can be judged with a plain
    !!      simple_exec prg=reconstruct3D projfile=<outfile> nstates=<nstates>
    !!  Particles the state stage left unassigned (label 0, see build_covariance_state_weights) stay at
    !!  state 0 and are excluded by reconstruct3D, exactly as they are excluded from the kernel states.
    subroutine write_discrete_state_project( spproj, pinds, labels, nstates, outfile )
        type(sp_project), intent(inout) :: spproj
        integer,          intent(in)    :: pinds(:), labels(:), nstates
        type(string),     intent(in)    :: outfile
        type(sp_project)     :: outproj
        logical, allocatable :: assigned(:)
        integer :: i, iptcl, state, nptcls, nexcluded
        if( size(pinds) < 1 .or. size(labels) /= size(pinds) .or. nstates < 2 ) &
            &THROW_HARD('invalid flex_pca discrete-state assignment')
        if( len_trim(outfile%to_char()) == 0 ) THROW_HARD('flex_pca discrete-state output project is empty')
        nptcls = spproj%os_ptcl3D%get_noris()
        allocate(assigned(nptcls), source=.false.)
        call outproj%copy(spproj)
        call outproj%update_projinfo(outfile)
        do iptcl = 1,nptcls
            call outproj%os_ptcl3D%set_state(iptcl,0)
        end do
        nexcluded = 0
        do i = 1,size(pinds)
            iptcl = pinds(i)
            state = labels(i)
            if( iptcl < 1 .or. iptcl > nptcls ) THROW_HARD('flex_pca discrete-state particle index outside project')
            if( assigned(iptcl) ) THROW_HARD('duplicate particle in flex_pca discrete-state assignment')
            if( state > nstates ) THROW_HARD('flex_pca discrete-state label outside state range')
            if( state < 1 )then
                nexcluded = nexcluded + 1
            else
                call outproj%os_ptcl3D%set_state(iptcl,state)
            endif
            assigned(iptcl) = .true.
        end do
        call outproj%write(outfile)
        write(logfhandle,'(A,A)') '>>> FLEX_PCA DISCRETE-STATE PROJECT: ',outfile%to_char()
        do state = 1,nstates
            write(logfhandle,'(A,I0,A,I0)') '>>> FLEX_PCA DISCRETE-STATE state=',state, &
                &' population=',count(labels==state)
        end do
        if( nexcluded > 0 ) write(logfhandle,'(A,I0)') &
            &'>>> FLEX_PCA DISCRETE-STATE unassigned particles left at state=0: ',nexcluded
        write(logfhandle,'(A,A,A,I0)') '>>> RECONSTRUCT WITH: simple_exec prg=reconstruct3D projfile=', &
            &outfile%to_char(),' nstates=',nstates
        call flush(logfhandle)
        call outproj%kill
        deallocate(assigned)
    end subroutine write_discrete_state_project

    !>  Latent targets of the automatically-placed conformational trajectory, so the sweep the
    !!  flex_pca_traj_###.mrc volumes represent can be read off without re-running.
    subroutine write_trajectory_targets( axis, nstates, ncomp, targets, bandwidths, neff )
        integer, intent(in) :: axis, nstates, ncomp
        real,    intent(in) :: targets(ncomp,nstates), bandwidths(nstates), neff(nstates)
        integer :: s, q, iu
        call del_file('flex_pca_trajectory_targets.txt')
        open(newunit=iu, file='flex_pca_trajectory_targets.txt', status='replace', action='write')
        write(iu,'(A,I0,A)') '# conformational trajectory placed along latent component ',axis, &
            &', selected by half-set reproducibility of the endpoint difference'
        write(iu,'(A)') '# state  bandwidth  effective_particles  axis_coord  t1..tncomp'
        do s = 1, nstates
            write(iu,'(I5,1X,ES16.8,1X,ES16.8,1X,ES16.8)',advance='no') s, bandwidths(s), neff(s), targets(axis,s)
            do q = 1, ncomp
                write(iu,'(1X,ES16.8)',advance='no') targets(q,s)
            end do
            write(iu,*)
        end do
        close(iu)
    end subroutine write_trajectory_targets

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
        write(u,'(A,I0)') 'selected_columns=',params%ncols
        write(u,'(A,I0)') 'column_separation=',params%column_separation
        write(u,'(A,I0)') 'probe_iters=',params%n_probe_iters
        write(u,'(A,L1)') 'heldout_embedding=',params%l_heldout
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
        close(u)
    end subroutine write_covariance_manifest

    ! ============ SELF-CONTAINED TESTS ============ No project, no images, no data files.

    !> Embedding-cache round trip.
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

    !> Epanechnikov kernel contract: compact support, peak normalised to 1, neff bounded by the raw
    !! support count, and widening that fires only when the support falls short and stops at the cap.
    subroutine test_flex_pca_kernel_bandwidth()
        integer,  parameter :: NP = 400
        integer  :: i, nsupp
        real(dp) :: dist(NP), h_out, h_wide
        real     :: w(NP), neff, neff_wide
        write(logfhandle,'(A)') '>>> TEST flex_pca kernel weights at bandwidth'
        do i = 1, NP
            dist(i) = real(i,dp)                    ! squared distances 1..NP
        end do
        ! h^2 = 100 -> exactly 99 particles strictly inside the support (dist < 100)
        call kernel_weights_at_bandwidth(dist, NP, 10.d0, 1, w, h_out, neff)
        nsupp = count(w > 0.)
        if( nsupp /= 99 ) THROW_HARD('kernel support count wrong for h^2=100')
        if( abs(h_out - 10.d0) > 1.d-12 ) THROW_HARD('bandwidth grew when support already sufficed')
        if( any(w(100:) /= 0.) ) THROW_HARD('kernel is not compactly supported')
        if( abs(maxval(w) - 1.) > 1.e-6 ) THROW_HARD('kernel peak is not normalised to 1')
        if( neff > real(nsupp) ) THROW_HARD('neff exceeds raw support count')
        if( neff < 1. ) THROW_HARD('neff below 1 with non-empty support')
        ! demand more support than h admits: the bandwidth must widen, bounded by COV_MAX_BW_GROW
        call kernel_weights_at_bandwidth(dist, NP, 10.d0, 200, w, h_wide, neff_wide)
        if( h_wide <= 10.d0 ) THROW_HARD('bandwidth did not widen when support fell short of min_neff')
        if( h_wide > 10.d0*1.3d0**COV_MAX_BW_GROW + 1.d-9 ) THROW_HARD('bandwidth widening exceeded its cap')
        if( count(w > 0.) <= nsupp ) THROW_HARD('widening did not increase support')
        if( neff_wide <= neff ) THROW_HARD('neff did not increase with bandwidth')
        write(logfhandle,'(A)') '>>>   PASSED (support, normalisation, neff bounds, bounded widening)'
    end subroutine test_flex_pca_kernel_bandwidth

    !> State-weight assembly on a synthetic two-cluster embedding.
    subroutine test_flex_pca_state_weights()
        integer,  parameter :: NPC = 150, NP = 2*NPC, NC = 2, NST = 2
        integer  :: i, q, r, state, nlab(NST)
        real(dp) :: z(NP,NC), eigvals(NC), prec(NC,NC,NP)
        real,     allocatable :: weights(:,:), targets(:,:), bandwidths(:), neff(:)
        integer,  allocatable :: labels(:)
        write(logfhandle,'(A)') '>>> TEST flex_pca state weights on a two-cluster embedding'
        ! two tight, well-separated blobs along component 1
        do i = 1, NPC
            z(i,      1) = -5.d0 + 0.01d0*real(mod(i,7),dp)
            z(i,      2) =  0.02d0*real(mod(i,5),dp)
            z(NPC+i,  1) =  5.d0 + 0.01d0*real(mod(i,7),dp)
            z(NPC+i,  2) =  0.02d0*real(mod(i,5),dp)
        end do
        do q = 1, NC
            eigvals(q) = 1.d0
        end do
        prec = 0.d0
        do i = 1, NP
            do q = 1, NC
                prec(q,q,i) = 1.d0                  ! identity posterior precision
            end do
        end do
        call build_covariance_state_weights(z, NP, NC, NC, NST, 0, 10, eigvals, prec, &
            &weights, targets, bandwidths, neff, labels)
        if( size(weights,1) /= NP .or. size(weights,2) /= NST ) THROW_HARD('weights shape wrong')
        if( size(labels)    /= NP ) THROW_HARD('labels shape wrong')
        if( any(labels < 0) .or. any(labels > NST) ) THROW_HARD('label outside 0..nstates')
        if( any(weights < 0.) ) THROW_HARD('negative kernel weight')
        do q = 1, NC
            do state = 1, NST
                if( real(targets(q,state),dp) < minval(z(:,q)) - 1.d-6 .or. &
                    &real(targets(q,state),dp) > maxval(z(:,q)) + 1.d-6 ) &
                    &THROW_HARD('state target lies outside the occupied latent range')
            end do
        end do
        nlab = 0
        do i = 1, NP
            if( labels(i) >= 1 ) nlab(labels(i)) = nlab(labels(i)) + 1
        end do
        do state = 1, NST
            if( nlab(state) < 1 ) THROW_HARD('a state drew no particles from a clearly bimodal embedding')
            if( neff(state) < 1. ) THROW_HARD('state has non-positive effective count')
        end do
        ! the two targets must separate along the component that carries the separation
        if( abs(real(targets(1,1),dp) - real(targets(1,2),dp)) < 5.d0 ) &
            &THROW_HARD('state targets did not separate along the bimodal component')
        deallocate(weights, targets, bandwidths, neff, labels)
        write(logfhandle,'(A)') '>>>   PASSED (shapes, labels, target hull, cluster separation)'
    end subroutine test_flex_pca_state_weights
end module simple_flex_pca_model
