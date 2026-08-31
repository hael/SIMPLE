!@descr: Standalone projection-aware low-rank covariance workflow for heterogeneous SPA data
module simple_flex_pca_model
use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
use simple_core_module_api
use simple_builder,                    only: builder
use simple_cmdline,                    only: cmdline
use simple_flex_pca_rec3D,         only: reconstruct_flex_weighted_states, flex_rec_box, flex_rec_smpd
use simple_flex_pca_em,         only: cov_env_int_pub, cov_perturb_project_poses, &
    &build_covariance_eigenbasis, embed_latents_with_contrast, &
    &estimate_covariance_mean, probe_subspace_iteration, align_basis_to_reference, probe_external_basis, &
    &save_probe_state, bag_basis_pool, basis_recs_from_images
use simple_image,                      only: image
use simple_parameters,                 only: parameters
use simple_reconstructor,              only: reconstructor
use simple_sigma2_files,               only: load_sigma2_groups
use simple_sp_project,                 only: sp_project
use simple_srch_sort_loc,              only: hpsort
use simple_flex_pca_distr,             only: flex_pca_is_worker, flex_pca_is_master, flex_pca_nparts, flex_pca_run_stage, &
    &flex_pca_set_fit, flex_pca_fit, FLEX_FIT_ALL, FLEX_FIT_A, FLEX_FIT_B
use simple_flex_pca_parts,             only: write_sigma_state, check_sigma_state, &
    &read_state_weights_round
use simple_flex_pca_distr,             only: PCA_STAGE_STATES, PCA_STAGE_EMBED
use simple_finch,                      only: finch_hierarchy, fit_finch, finch_representatives, &
    &                                          select_finch_level, refine_finch_level
use simple_kd_tree,                    only: kd_tree, knn_table
use simple_linalg,                     only: jacobi, eigsrt, matinv
use simple_flex_pca_merge,             only: flex_pca_merge_enabled, two_gate_state_merge, &
    &flex_pca_merge_force_on
implicit none
private
#include "simple_local_flags.inc"

public :: run_flex_pca
public :: test_flex_pca_embedding_cache_io, test_flex_pca_kernel_bandwidth
    public :: test_flex_pca_state_weights
public :: test_flex_pca_auto_settings
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
integer,          parameter :: COV_CACHE_VERSION = 3
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
        integer :: nptcls, ncomp, nstates, min_neff, state_axis, col_sep, neigs_req, nkern
        integer :: q, i, r, s, nfinch
        integer :: nstates_merged, kmix
        integer, allocatable :: merge_label(:)
        real,    allocatable :: merged_weights(:,:), merged_targets(:,:), merged_bw(:)
        real(dp), allocatable :: state_mass(:), merged_mass(:)
        real(dp) :: sumw_s, sumw2_s
        real(dp), allocatable :: finch_targets(:,:)
        logical :: l_finch_states
        logical :: l_rot
        real(dp), allocatable :: kdist(:,:), kfloor(:)
        real(dp), allocatable :: comp_rho(:)     ! per-component reliability, drives state-target ordering
        ! ---- B1: responsibility-delivered states from the probe-fitted mixture ----
        real(dp) :: aw_b2
        real(dp), allocatable :: pviews(:,:)     ! per-particle viewing AXIS, for the GMM orientation term
        character(len=:), allocatable :: cachedir, cachestr
        real(dp) :: sig2_eff, snr_best
        logical :: sigma_loaded, l_resume, l_split_eo
        integer(timer_int_kind) :: t_blk

        call validate_covariance_inputs(params, build, cline, pinds, nptcls)
        ! optional pose degradation, before anything reads the orientations (see the routine)
        call cov_perturb_project_poses(build, pinds, nptcls)
        call load_and_validate_sigma(params, build, cline, pinds, sigma_loaded)

        ! STATE-RECONSTRUCTION WORKER. The master already produced the weight table, so replicate its
        ! operator setup (ml_reg off, reconstruction Fourier band) and go straight to accumulation.
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
        ! No ceiling on the state count: over-provisioning is the only regime in which the merge recovers K.
        ! npreimages=0 selects the automatic ceiling; anything else is taken as the requested ceiling.
        ! npreimages is the ceiling; preimage_auto raises it and enables the collapse that makes a
        ! ceiling meaningful. An explicit npreimages alongside preimage_auto=yes is honoured as the
        ! ceiling -- auto then contributes only the merge.
        nstates = max(MIN_NSTATES, params%npreimages)
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
        ! env-only: inert by default (the GMM replaces it), live on the nbins>1 and SIMPLE_COV_GMM=0 opt-outs
        min_neff = params%min_neff
        call cov_env_int_pub('SIMPLE_COV_MIN_NEFF', min_neff)
        min_neff = max(20, min(nptcls, min_neff))
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

        call estimate_covariance_mean(params, build, mean_rec, pinds, nptcls)

        if( params%l_heldout )then
            call heldout_embedding(params, build, mean_rec, pinds, nptcls, col_sep, neigs_req, &
                &basis_recs, eigvals, ncomp, sig2_eff, z, contrast, latent_second, resid_energy, resid_mean_energy)
            if( state_axis > 0 ) state_axis = min(state_axis, min(ncomp, nkern))
            allocate(prior_precision(ncomp))
            do q = 1, ncomp
                prior_precision(q) = 1.d0 / max(eigvals(q), DTINY)
            end do
            write(logfhandle,'(A,I0)') '>>> FLEX_PCA held-out retained covariance components=',ncomp
            call flush(logfhandle)
        else if( cov_env_flag_on('SIMPLE_COV_BAG') )then
            ! two disjoint fits pooled into one rank-matched basis, then everything embedded through it
            call bagged_embedding(params, build, mean_rec, pinds, nptcls, col_sep, neigs_req, &
                &basis_recs, eigvals, ncomp, sig2_eff, z, contrast, latent_second, resid_energy, resid_mean_energy)
            if( flex_pca_is_worker() )then
                ! this worker's part file for the requested fit/stage is on disk; the master reduces it
                call qsys_job_finished(params, string('simple_flex_pca_model :: run_flex_pca'))
                call mean_rec%dealloc_rho; call mean_rec%kill
                deallocate(pinds)
                return
            endif
            if( state_axis > 0 ) state_axis = min(state_axis, min(ncomp, nkern))
            allocate(prior_precision(ncomp))
            do q = 1, ncomp
                prior_precision(q) = 1.d0 / max(eigvals(q), DTINY)
            end do
            write(logfhandle,'(A,I0)') '>>> FLEX_PCA bagged retained covariance components=',ncomp
            call flush(logfhandle)
        else

            call build_covariance_eigenbasis(params, build, mean_rec, pinds, nptcls, &
                &col_sep, neigs_req, basis_recs, eigvals, ncomp, sig2_eff)
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

            ! Optional Wiener E-step / weighted-backprojection M-step, to clean the noisy column directions.
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
        ! Resuming skips the split-half solve, so fall back to the spread-over-posterior-variance proxy.
        if( .not. allocated(comp_rho) )then
            allocate(comp_rho(ncomp))
            call component_reliability_proxy(z, latent_second, nptcls, ncomp, comp_rho)
            write(logfhandle,'(A)') '>>> FLEX_PCA resumed embedding: component reliability from the &
                &spread/posterior-variance proxy (no cached split-half rho)'
        endif
        call run_external_basis_probe(params, build, mean_rec, pinds, nptcls, ncomp, eigvals, &
            &sig2_eff, l_resume)


        ! Posterior latent shrinkage; SIMPLE_COV_ZSHRINK=0 restores the raw MAP z. Must stay
        ! after the cache write: the cache stores raw z, so a resume never shrinks twice.
        if( .not. cov_env_flag_off('SIMPLE_COV_ZSHRINK') )then
                call shrink_latents_gaussian(z, latent_second, nptcls, ncomp)
        endif

        ! Rotation runs BEFORE state placement, so every downstream stage sees one consistent frame.
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

        ! Sphere the latent cloud before state placement. A component's influence on any
        ! Euclidean or graph geometry goes as its variance, but its usefulness for telling
        ! states apart does not: measured on EMPIAR-10076, component 4 carries half the
        ! class information of component 1 and one seventy-fourth of the distance. The
        ! latent spectrum spans 1800:1 with 95% of the variance in two directions and
        ! off-diagonal correlations up to 0.74, so neither the diffusion graph's diagonal
        ! standardisation nor a k-means initialisation reads the geometry that carries the
        ! states. This is a change of frame only -- U' z' = U z holds exactly -- so the
        ! model, the reconstructions and the residuals are untouched.
        ! OPT-IN, not default. Measured at K=20 on the 87k clean particles, EM arm: whitening
        ! alone is +0.004 ARI, which is inside the noise, because the diffusion graph already
        ! standardises per-component and only the decorrelation is new. Worse, it does not
        ! compose with the mean-shaped deflation that actually fixes this -- deflation alone is
        ! 0.176 and deflation plus whitening is 0.154. That makes sense: whitening is a
        ! workaround for latents whose largest directions are nuisance, so once the nuisance is
        ! removed at source, equalising the axes just promotes the noise directions.
        if( cov_env_flag_on('SIMPLE_COV_ZWHITEN') )then
            if( l_resume .and. .not. file_exists('flex_pca_pc001.mrc') )then
                call whiten_latent_frame(params, ncomp, z, nptcls, latent_second, l_rot, eigdir=cachedir)
            else
                call whiten_latent_frame(params, ncomp, z, nptcls, latent_second, l_rot)
            endif
            ! comp_rho is indexed BY COMPONENT and the whitening remixes components, so the
            ! split-half reliabilities measured in the old frame no longer name the same
            ! directions. They cannot be carried through a rotation, so recompute the proxy
            ! in the frame the state placement will actually read.
            if( l_rot )then
                call component_reliability_proxy(z, latent_second, nptcls, ncomp, comp_rho)
                write(logfhandle,'(A)') '>>> FLEX_PCA component reliability recomputed in the whitened frame'
            endif
        endif

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
        ! Per-particle viewing AXIS, folded antipodally: +n and -n are mirror projections, so orientation
        ! bias lives in the axis and never in the mean resultant. Only the GMM's coverage term reads it.
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
        ! ---- B2: amplitude soft-weighting (SIMPLE_COV_AWEIGHT=1) ----
        ! The fitted per-particle amplitude separates junk at AUC 0.784 where the residual is
        ! blind (0.540): a featureless particle is fit perfectly by a -> 0. A soft gate on a_i
        ! downweights those particles in every state reconstruction without the hard prune's
        ! clean-particle loss. Sigmoid midpoint 0.65, width 0.10 -- between the measured junk
        ! mean (0.79) and clean mean (1.14), biting mostly below the junk mode.
        if( cov_env_flag_on('SIMPLE_COV_AWEIGHT') )then
            i = 0
            do q = 1, nptcls
                aw_b2 = 1.d0 / (1.d0 + exp(-(contrast(q) - 0.65d0)/0.10d0))
                state_weights(q,:) = state_weights(q,:) * real(aw_b2)
                if( aw_b2 < 0.05d0 ) i = i + 1
            end do
            write(logfhandle,'(A,I0,A,I0,A)') '>>> FLEX_PCA AWEIGHT: ', i, ' of ', nptcls, &
                &' particles effectively excluded (amplitude gate)'
            call flush(logfhandle)
        endif

        ! ---- AUTO-K: reproducible-state count (SIMPLE_COV_AUTO_K=1) ----
        ! Point-estimating K is measured broken on real data (BIC/ICL/FINCH/merge all fail there);
        ! instead: refit the state placement INDEPENDENTLY on the even and odd halfsets against the
        ! same basis, then ask, per full-fit state, whether its members re-cluster coherently under
        ! BOTH half fits. Two complementary criteria, either passes: MODAL coherence (a majority of
        ! the state's members share one half-fit cluster -- protects small compact states) or
        ! territory OWNERSHIP (a majority land in half-clusters this state plurality-owns --
        ! protects large basins the half fits merely subdivide). Measured 2026-08-16: the union
        ! keeps 13/13 substantial states on pristine Ribosembly, kills 16/16 phantom states fitted
        ! on structureless noise, and on cleaned EMPIAR-10076 prunes exactly the diffuse
        ! non-recurring basins. Mutual-nearest target matching with a distance threshold is
        ! REFUTED: real broad basins drift ~1 sd between half fits, the same scale as the
        ! inter-target spacing, so any tau loses genuine states (3 of 16 on Ribosembly).
        ! Design + validation targets: doc/implementation_notes/flex_pca_automation_plan.md.
        block
            integer :: vak, nE, nO, i2, s, a, j
            integer, allocatable :: idxE(:), idxO(:), labE(:), labO(:), labA(:), labB(:)
            integer, allocatable :: pop(:), cmE(:,:), cmO(:,:), ownE(:), ownO(:)
            real(dp), allocatable :: zE(:,:), zO(:,:), pE(:,:,:), pO(:,:,:), vwE(:,:), vwO(:,:)
            real,     allocatable :: wE(:,:), wO(:,:), tE(:,:), tO(:,:), bwE(:), bwO(:), nfE(:), nfO(:)
            real(dp), allocatable :: sdv(:)
            real(dp) :: d2, dbest, fmE, fmO, foE, foO
            logical,  allocatable :: reproducible(:)
            real(dp), parameter :: AUTOK_COH = 0.45d0
            vak = 0
            call cov_env_int_pub('SIMPLE_COV_AUTO_K', vak)
            if( vak > 0 )then
                nE = 0
                do i = 1, nptcls
                    if( build%spproj_field%get_eo(pinds(i)) == 0 ) nE = nE + 1
                end do
                nO = nptcls - nE
                allocate(idxE(nE), idxO(nO))
                nE = 0; nO = 0
                do i = 1, nptcls
                    if( build%spproj_field%get_eo(pinds(i)) == 0 )then
                        nE = nE + 1; idxE(nE) = i
                    else
                        nO = nO + 1; idxO(nO) = i
                    endif
                end do
                allocate(zE(nE,ncomp), pE(ncomp,ncomp,nE), vwE(3,nE))
                allocate(zO(nO,ncomp), pO(ncomp,ncomp,nO), vwO(3,nO))
                do i = 1, nE
                    zE(i,:) = z(idxE(i),:); pE(:,:,i) = latent_second(:,:,idxE(i)); vwE(:,i) = pviews(:,idxE(i))
                end do
                do i = 1, nO
                    zO(i,:) = z(idxO(i),:); pO(:,:,i) = latent_second(:,:,idxO(i)); vwO(:,i) = pviews(:,idxO(i))
                end do
                call build_covariance_state_weights(zE, nE, ncomp, nkern, nstates, state_axis, &
                    &max(20,min_neff/2), eigvals, pE, wE, tE, bwE, nfE, labE, comp_rho=comp_rho, views=vwE)
                call build_covariance_state_weights(zO, nO, ncomp, nkern, nstates, state_axis, &
                    &max(20,min_neff/2), eigvals, pO, wO, tO, bwO, nfO, labO, comp_rho=comp_rho, views=vwO)
                ! standardized metric from the FULL latent spread
                allocate(sdv(ncomp))
                do i2 = 1, ncomp
                    sdv(i2) = max(sqrt(sum((z(:,i2)-sum(z(:,i2))/nptcls)**2)/real(nptcls,dp)), 1.d-12)
                end do
                ! assign EVERY particle to its nearest half-fit target in the standardized metric
                allocate(labA(nptcls), labB(nptcls))
                !$omp parallel do default(shared) private(i,a,d2,dbest) schedule(static)
                do i = 1, nptcls
                    dbest = huge(1.d0)
                    do a = 1, nstates
                        d2 = sum(((z(i,:)-real(tE(:,a),dp))/sdv)**2)
                        if( d2 < dbest )then
                            dbest = d2; labA(i) = a
                        endif
                    end do
                    dbest = huge(1.d0)
                    do a = 1, nstates
                        d2 = sum(((z(i,:)-real(tO(:,a),dp))/sdv)**2)
                        if( d2 < dbest )then
                            dbest = d2; labB(i) = a
                        endif
                    end do
                end do
                !$omp end parallel do
                ! cross-tabulate full-fit membership against each half labeling
                allocate(pop(nstates), cmE(nstates,nstates), cmO(nstates,nstates), &
                    &ownE(nstates), ownO(nstates))
                pop = 0; cmE = 0; cmO = 0
                do i = 1, nptcls
                    s = labels(i)
                    if( s < 1 .or. s > nstates ) cycle
                    pop(s) = pop(s) + 1
                    cmE(s,labA(i)) = cmE(s,labA(i)) + 1
                    cmO(s,labB(i)) = cmO(s,labB(i)) + 1
                end do
                do j = 1, nstates
                    ownE(j) = 0; ownO(j) = 0
                    if( sum(cmE(:,j)) > 0 ) ownE(j) = maxloc(cmE(:,j), dim=1)
                    if( sum(cmO(:,j)) > 0 ) ownO(j) = maxloc(cmO(:,j), dim=1)
                end do
                allocate(reproducible(nstates), source=.false.)
                do s = 1, nstates
                    if( pop(s) < 1 ) cycle
                    fmE = real(maxval(cmE(s,:)),dp) / real(pop(s),dp)
                    fmO = real(maxval(cmO(s,:)),dp) / real(pop(s),dp)
                    foE = 0.d0; foO = 0.d0
                    do j = 1, nstates
                        if( ownE(j) == s ) foE = foE + real(cmE(s,j),dp)
                        if( ownO(j) == s ) foO = foO + real(cmO(s,j),dp)
                    end do
                    foE = foE / real(pop(s),dp); foO = foO / real(pop(s),dp)
                    reproducible(s) = ( fmE >= AUTOK_COH .and. fmO >= AUTOK_COH ) .or. &
                                     &( foE >= AUTOK_COH .and. foO >= AUTOK_COH )
                    write(logfhandle,'(A,I0,A,I0,4(A,F5.2))') '>>> FLEX_PCA AUTO-K: state ', s, &
                        &' pop=', pop(s), ' modal=', fmE, '/', fmO, ' owned=', foE, '/', foO
                end do
                write(logfhandle,'(A,I0,A,I0,A)') '>>> FLEX_PCA AUTO-K: ', count(reproducible), &
                    &' of ', nstates, ' states reproduce across independent halfset fits'
                do s = 1, nstates
                    if( .not. reproducible(s) )then
                        write(logfhandle,'(A,I0,A)') '>>> FLEX_PCA AUTO-K: state ', s, &
                            &' NOT reproducible -- pruned from delivery'
                        state_weights(:,s) = 0.
                        where( labels == s ) labels = 0
                    endif
                end do
                call flush(logfhandle)
                deallocate(idxE, idxO, zE, zO, pE, pO, vwE, vwO, wE, wO, tE, tO, bwE, bwO, nfE, nfO, &
                    &labE, labO, labA, labB, pop, cmE, cmO, ownE, ownO, sdv, reproducible)
            endif
        end block
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

        if( params%nbins > 1 )then
            call cv_select_bandwidths(params, build, pinds, nptcls, nstates, params%nbins, min_neff, &
                &kdist, kfloor, state_weights, bandwidths, neff)
        endif
        ! combined states and both halfsets in ONE pass; combined == even + odd exactly
        params%outvol = 'flex_pca_state_001.mrc'
        t_blk = tic()
        call reconstruct_flex_weighted_states(params, build, pinds, state_weights, nstates, &
            &floor_rho=.true., outvol_even=string('flex_pca_even_state_001.mrc'), &
            &outvol_odd=string('flex_pca_odd_state_001.mrc'))
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
                nstates = nstates_merged
                call reconstruct_flex_weighted_states(params, build, pinds, state_weights, nstates, &
                    &floor_rho=.true., outvol_even=string('flex_pca_even_state_001.mrc'), &
                    &outvol_odd=string('flex_pca_odd_state_001.mrc'))
            endif
            deallocate(merge_label)
            write(logfhandle,'(A,F9.1)') '>>> FLEX_PCA STAGE two_gate_merge seconds=', toc(t_blk)
        endif
        call write_covariance_tables(build, pinds, z, eigvals, prior_precision, state_weights, labels, &
            &targets, bandwidths, neff, resid_energy, resid_mean_energy, contrast)
        ! Hard labels into the project itself, so the assignment can be judged by an INDEPENDENT
        ! reconstructor. Master-only: a worker shares the master's projfile, and the worker return
        ! above covers the default path but NOT l_heldout, which falls through to here.
        if( flex_pca_is_master() )then
            call write_discrete_state_project(build%spproj, pinds, labels, nstates, params%projfile)
        endif
        allocate(half_weights(nptcls,nstates), source=state_weights)
        ! nonuniform filtering LAST, so every delivered map is filtered the same way
        if( trim(params%nufilt) == 'yes' ) call apply_consensus_nu_filter(params, nstates)
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
    end subroutine run_flex_pca

    subroutine validate_covariance_inputs( params, build, cline, pinds, nptcls )
        type(parameters), intent(inout) :: params
        type(builder),    intent(inout) :: build
        class(cmdline),   intent(in)    :: cline
        integer, allocatable, intent(out) :: pinds(:)
        integer, intent(out) :: nptcls
        integer :: q, i, cnt
        integer, allocatable :: sel(:)
        if( trim(params%oritype) /= 'ptcl3D' ) THROW_HARD('flex_pca requires oritype=ptcl3D')
        if( .not. cline%defined('vol1') ) THROW_HARD('flex_pca requires vol1=<consensus mean map>')
        if( trim(params%ptcl_src) /= 'raw' ) THROW_HARD('flex_pca currently requires ptcl_src=raw')
        if( build%spproj_field%get_n('state') /= 1 ) THROW_HARD('flex_pca requires a project with an eo split')
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
        if( nptcls < 100 ) THROW_HARD('flex_pca requires at least 100 active particles')
        if( build%spproj%os_ptcl2D%get_noris() > 0 .and. &
            &build%spproj%os_ptcl2D%get_noris() /= build%spproj%os_ptcl3D%get_noris() ) &
            &THROW_HARD('flex_pca requires matching ptcl2D and ptcl3D rows')
        ! eo split by index parity when the project's is degenerate. MASTER-ONLY and must be persisted: a
        ! worker counts over its own range, so it would silently build a different -- equally valid -- split.
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

    subroutine load_and_validate_sigma( params, build, cline, pinds, loaded )
        type(parameters), intent(inout) :: params
        type(builder),    intent(inout) :: build
        class(cmdline),   intent(inout) :: cline
        integer,          intent(in)    :: pinds(:)
        logical,          intent(out)   :: loaded
        integer :: i, k, iptcl, noris, kto
        call load_sigma2_groups(params, build%pftc, build%esig, build%spproj_field, cline, loaded)
        if( .not. loaded )then
            ! gen_fplane4rec follows norm_noise_taper_edge_pad_fft, so a unit spectrum is the correct fallback
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

    !> Median of chi-squared with k dof, Wilson-Hilferty: k*(1 - 2/(9k))^3. Good to 3 % at k=1, which is
    !! far inside the tolerance of a bandwidth FLOOR and needs no gamma inverse.
    pure real(dp) function chi2_median( k )
        integer, intent(in) :: k
        real(dp) :: kk
        kk = real(max(k,1),dp)
        chi2_median = kk * (1.d0 - 2.d0/(9.d0*kk))**3
    end function chi2_median


    !> Posterior mean under a Gaussian population prior with S = Cov(z) - mean_i(P_i^-1)
    !! (the observed scatter carries each particle's posterior covariance once), clipped PSD.
    subroutine shrink_latents_gaussian( z, precision, nptcls, ncomp )
        integer,  intent(in)    :: nptcls, ncomp
        real(dp), intent(inout) :: z(nptcls,ncomp)
        real(dp), intent(in)    :: precision(ncomp,ncomp,nptcls)
        real(dp) :: mu(ncomp), C_obs(ncomp,ncomp), Nbar(ncomp,ncomp), S(ncomp,ncomp)
        real(dp) :: Sinv(ncomp,ncomp), evec(ncomp,ncomp), ev(ncomp), rhsmu(ncomp)
        real(dp) :: Pinv(ncomp,ncomp), A(ncomp,ncomp), Ainv(ncomp,ncomp), rhs(ncomp)
        real(dp) :: floor_ev, dz2, z2
        integer  :: i, q, nrot, errflg, nbad
        mu = sum(z, dim=1) / real(nptcls,dp)
        C_obs = 0.d0
        do i = 1, nptcls
            do q = 1, ncomp
                C_obs(:,q) = C_obs(:,q) + (z(i,:) - mu) * (z(i,q) - mu(q))
            end do
        end do
        C_obs = C_obs / real(nptcls,dp)
        Nbar = 0.d0
        nbad = 0
        !$omp parallel do schedule(static) default(shared) private(i,Pinv,errflg) &
        !$omp reduction(+:Nbar,nbad)
        do i = 1, nptcls
            call matinv(precision(:,:,i), Pinv, ncomp, errflg)
            if( errflg /= 0 )then
                nbad = nbad + 1
            else
                Nbar = Nbar + Pinv
            endif
        end do
        !$omp end parallel do
        if( nbad >= nptcls ) THROW_HARD('shrink_latents_gaussian: no invertible posterior precisions')
        Nbar = Nbar / real(nptcls - nbad,dp)
        ! deconvolved population covariance, PSD-clipped (jacobi overwrites its input; matinv does not)
        S = C_obs - Nbar
        call jacobi(S, ncomp, ncomp, ev, evec, nrot)
        if( maxval(ev) <= 0.d0 )then
            write(logfhandle,'(A)') '>>> FLEX_PCA LATENT SHRINKAGE skipped: deconvolved population &
                &covariance has no positive spectrum'
            return
        endif
        floor_ev = 1.d-6 * maxval(ev)
        ev = max(ev, floor_ev)
        S    = matmul(evec, matmul(diagmat(ev, ncomp),      transpose(evec)))
        Sinv = matmul(evec, matmul(diagmat(1.d0/ev, ncomp), transpose(evec)))
        rhsmu = matmul(Sinv, mu)
        dz2 = 0.d0
        z2  = 0.d0
        !$omp parallel do schedule(static) default(shared) private(i,A,Ainv,rhs,errflg) &
        !$omp reduction(+:dz2,z2)
        do i = 1, nptcls
            A   = Sinv + precision(:,:,i)
            rhs = rhsmu + matmul(precision(:,:,i), z(i,:))
            call matinv(A, Ainv, ncomp, errflg)
            if( errflg /= 0 ) cycle          ! leave this particle at its raw MAP
            rhs = matmul(Ainv, rhs)          ! reuse rhs as the shrunk latent
            dz2 = dz2 + sum((rhs - z(i,:))**2)
            z2  = z2  + sum(z(i,:)**2)
            z(i,:) = rhs
        end do
        !$omp end parallel do
        write(logfhandle,'(A,I0,A,ES10.3,A,ES10.3,A,F6.3)') &
            &'>>> FLEX_PCA LATENT SHRINKAGE (Gaussian population prior): n=',nptcls, &
            &'  prior spectrum max=',maxval(ev),' min=',minval(ev), &
            &'  mean |dz|/|z|=',real(sqrt(dz2/max(z2,DTINY)))
        if( nbad > 0 ) write(logfhandle,'(A,I0,A)') '>>> FLEX_PCA LATENT SHRINKAGE: ',nbad, &
            &' particles kept raw (posterior precision not invertible)'
        call flush(logfhandle)

    contains

        pure function diagmat( d, n ) result( D2 )
            integer,  intent(in) :: n
            real(dp), intent(in) :: d(n)
            real(dp) :: D2(n,n)
            integer  :: j
            D2 = 0.d0
            do j = 1, n
                D2(j,j) = d(j)
            end do
        end function diagmat

    end subroutine shrink_latents_gaussian


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
    !! every particle an Epanechnikov weight per state from its Mahalanobis distance to that target. Each
    !! bandwidth is floored so at least min_neff particles fall inside its support.
    subroutine build_covariance_state_weights( z, nptcls, ncomp, nkern, nstates, axis, min_neff, &
        &eigvals, precision, weights, targets, bandwidths, neff, labels, dist_out, bfloor_out, targets_in, &
        &zmetric, comp_rho, views )
        integer,  intent(in) :: nptcls, ncomp, nkern, nstates, axis, min_neff
        real(dp), intent(in) :: z(nptcls,ncomp), eigvals(ncomp)
        real(dp), intent(in) :: precision(ncomp,ncomp,nptcls)   ! per-particle latent precision Pi_i
        real, allocatable, intent(out) :: weights(:,:), bandwidths(:), neff(:)
        real, allocatable, intent(out) :: targets(:,:)          ! (ncomp,nstates) latent target coordinates
        integer, allocatable, intent(out) :: labels(:)
        ! distances and floors, so cv_select_bandwidths can rebuild weights without redoing the nk^2 forms
        real(dp), allocatable, optional, intent(out) :: dist_out(:,:), bfloor_out(:)
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
        real(dp) :: h, d2, u2, sumw, sumw2, best, zspread, bmin, chi2med, wsum_i
        real(dp) :: orient_lam, cdeconv
        integer  :: nrenorm, ispace
        integer  :: i, q, r, state, best_state, grow, nfed, occmax, ifloor, nunassigned, nsupp
        integer  :: nk, errflg
        logical  :: l_relpath, l_diffuse, l_gmm
        character(len=12) :: bwsrc
        nk = max(1, min(ncomp, nkern))
        allocate(wcomp(nk), tvec(nk), tcen(nk,nstates), dist(nptcls), dvec(nk), mvec(nk))
        wcomp = 1.d0
        ! STANDARDIZED PLACEMENT (SIMPLE_COV_STDZ=0 opts out). Eigenvalue weighting would concentrate every
        ! target along the highest-variance components, which are not the conformational ones.
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
        ! Restricting to nk components MARGINALISES the precision (invert, slice, re-invert); slicing the
        ! precision directly would condition on the dropped components instead.
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
        ! Reliability-ordered path is the default when a reliability vector arrived; SIMPLE_COV_KMEANS=1
        ! recovers k-means. Set inside the placing branch, as it also selects the along-path weighting.
        l_relpath = .false.
        allocate(weights(nptcls,nstates), targets(ncomp,nstates), bandwidths(nstates), neff(nstates), labels(nptcls))
        if( present(dist_out)   ) allocate(dist_out(nptcls,nstates))
        if( present(bfloor_out) ) allocate(bfloor_out(nstates))
        if( present(targets_in) )then
            do state = 1, nstates
                tcen(:,state) = targets_in(1:nk,state)
            end do
            if( nstates >= 2 )then
                ! ORDERED CURVE: external targets are read as a polyline through latent space (the
                ! supplier's row order IS the path order). Each particle takes the arc-length
                ! coordinate of its projection onto the polyline, so distance, bandwidth and the
                ! frames all live ON the curve: off-curve directions (noise, and any motion
                ! components not in the curve) cannot strand particles, the bandwidth derives from
                ! target spacing, and the GMM refit is skipped exactly as for the axis path. This is
                ! what lets a high-amplitude motion be rendered along its winding trajectory through
                ! the eigencomponent ladder instead of faded along one linear axis.
                allocate(ppath(nptcls), tpath(nstates))
                call project_onto_target_polyline(z(:,1:nk), nptcls, nk, tcen, nstates, ppath, tpath)
                l_relpath = .true.
                write(logfhandle,'(A,I0,A,I0,A)') '>>> FLEX_PCA state targets: supplied externally over ', &
                    &nk,' components, points=',nstates,' -- treated as an ORDERED CURVE (arc-length kernel)'
            else
                write(logfhandle,'(A,I0,A,I0)') '>>> FLEX_PCA state targets: supplied externally over ',nk, &
                &' components, points=',nstates
            endif
        else if( axis < 0 )then
            call path_latent_targets(z(:,1:nk), nptcls, nk, nstates, wcomp, tcen)
            write(logfhandle,'(A,I0,A,I0)') '>>> FLEX_PCA state targets: density-spread path over ',nk, &
                &' components, points=',nstates
        else if( axis == 0 .and. present(comp_rho) .and. .not. cov_env_flag_on('SIMPLE_COV_KMEANS') &
                &.and. .not. cov_env_flag_on('SIMPLE_COV_RELPATH') )then
            ! DEFAULT: diffusion k-center. Handles a continuous reaction coordinate and branched compositional
            ! states with the same constants, because a curve is a degenerate graph.
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
            ! 1-D equal-occupancy path (SIMPLE_COV_RELPATH=1). Correct on a genuine reaction coordinate, but it
            ! MERGES states on a branched manifold -- kept for 1-D data and as the diffusion fallback.
            allocate(ppath(nptcls), tpath(nstates))
            call reliability_path_targets(z(:,1:nk), nptcls, nk, nstates, wcomp, comp_rho(1:nk), tcen, &
                &proj_out=ppath, tproj_out=tpath)
            l_relpath = .true.
            write(logfhandle,'(A,I0,A,I0)') '>>> FLEX_PCA state targets: reliability-ordered &
                &equal-occupancy path over ',nk,' components, points=',nstates
        else if( axis == 0 )then
            call kmeans_latent_targets(z(:,1:nk), nptcls, nk, nstates, wcomp, tcen)
            write(logfhandle,'(A,I0,A,I0)') '>>> FLEX_PCA state targets: k-means over ',nk, &
                &' components, k=',nstates
        else
            ! equal-occupancy slices along one component; each target is the slice MEAN over all nk components
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
            ! Score this placement ALONG THE AXIS, exactly as the reliability path does. The full-nk
            ! posterior quadratic form measures distance in every direction, and on a continuum the
            ! off-axis directions are noise -- they strand on-axis particles and hand the bandwidth
            ! rule a chi2(nk) noise scale instead of the target spacing. The slices were cut on this
            ! coordinate, so this is the coordinate the kernel must live on.
            allocate(ppath(nptcls), tpath(nstates))
            do i = 1, nptcls
                ppath(i) = z(i,axis)
            end do
            do state = 1, nstates
                tpath(state) = tcen(axis,state)
            end do
            l_relpath = .true.
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
                ! ALONG-PATH distance, so the off-path directions -- mostly noise -- cannot strand an on-path particle
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
                ! Tie the kernel to TARGET SPACING, not the latent dimension: chi2(nk) grows with nk, so the dimension
                ! alone forces kernels wide enough to swallow neighbouring targets -- which is why lowering nkern
                ! strands particles. Support is dist < h^2 = 2*bmin, hence the half.
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
            ! the chi2 floor is only meaningful if the posterior quadratic form really is on a chi2(nk) scale
            write(logfhandle,'(A,I3,A,ES11.3,A,ES11.3,A,ES11.3,A,A)') '>>>   state=',state, &
                &' dist: median=',real(sorted(max(1,nptcls/2)),dp),' p95=', &
                &real(sorted(max(1,nint(0.95*real(nptcls)))),dp),'  nn_floor=',real(sorted(ifloor),dp), &
                &'  bandwidth floor from ', bwsrc
            ! Enclosed population grows like h^nk (a 1.3x step is ~190x at nk=20); the floor should make this a no-op
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
        ! Tied-covariance mixture by default; SIMPLE_COV_GMM=0 recovers the kernel. The kernel loop above
        ! still runs: dist_out feeds cv_select_bandwidths and its quantiles diagnose the chi2 scale.
        l_gmm = .not. cov_env_flag_off('SIMPLE_COV_GMM')
        ! EQUAL-MASS PLACEMENT IS NOT A GMM INITIALISATION. The tied-covariance mixture is a
        ! discrete-state model: it re-fits the means, and on a continuum with one dense mode every
        ! component slides into that mode -- which silently UNDOES the equal-occupancy placement that
        ! was just constructed (measured on the RNA data: sextiles in, one state holding 88% of the
        ! particles out). Where the targets carry equal mass by construction, keep them and let the
        ! along-path kernel deliver the frames. SIMPLE_COV_GMM=1 forces the refit back on for A/B.
        if( l_relpath .and. .not. cov_env_flag_on('SIMPLE_COV_GMM') )then
            if( l_gmm ) write(logfhandle,'(A)') '>>> FLEX_PCA equal-mass targets: GMM refit SKIPPED &
                &(it would re-fit the means onto the dominant mode); along-path kernel weights kept'
            l_gmm = .false.
        endif
        if( l_gmm )then
            orient_lam = 0.d0
            call cov_env_dp('SIMPLE_COV_ORIENT_LAM', orient_lam)
            ! Measurement-error deconvolution scale. Off by default: c is a property of the DATASET's
            ! noise level, not a universal constant, so it has to be calibrated per project (see the
            ! header of gmm_state_weights for the ground-truth-free procedure) rather than shipped.
            cdeconv = 0.d0
            call cov_env_dp('SIMPLE_COV_DECONV', cdeconv)
            if( present(views) .and. orient_lam > 0.d0 )then
                call gmm_state_weights(z, nptcls, ncomp, nk, nstates, tcen, wcomp, weights, neff, &
                    &bandwidths, labels, views=views, orient_lam=orient_lam, &
                    &precision=precision, cdeconv=cdeconv)
            else
                call gmm_state_weights(z, nptcls, ncomp, nk, nstates, tcen, wcomp, weights, neff, &
                    &bandwidths, labels, precision=precision, cdeconv=cdeconv)
            endif
            ! tcen now holds the FITTED means; refresh the reported targets to describe the delivered maps
            do state = 1, nstates
                targets(1:nk,state) = real(tcen(:,state))
            end do
        endif
        ! argmax weight, 0 outside EVERY kernel support. Defaulting to state 1 would pile the unassigned
        ! onto the first state and fake a concentration failure.
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
        ! SIMPLE_COV_NORMW=1: partition of unity, so a particle in many kernels is not counted at full strength
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

    !> Arc-length coordinate of every particle on the polyline through the supplied targets.
    !! Projection is Euclidean in the raw latent (the units the targets are given in); each particle
    !! takes the closest point over all segments, clamped to the segment ends, so particles beyond
    !! either terminus map to the terminus and the end frames collect the tails.
    subroutine project_onto_target_polyline( z, nptcls, nk, tcen, nstates, ppath, tpath )
        integer,  intent(in)  :: nptcls, nk, nstates
        real(dp), intent(in)  :: z(nptcls,nk), tcen(nk,nstates)
        real(dp), intent(out) :: ppath(nptcls), tpath(nstates)
        real(dp) :: seg(nk,max(1,nstates-1)), sl2(max(1,nstates-1)), seglen(max(1,nstates-1))
        real(dp) :: dz(nk), t, d2, best_d2, best_c
        integer  :: i, s, q
        tpath(1) = 0.d0
        do s = 1, nstates-1
            seg(:,s)   = tcen(:,s+1) - tcen(:,s)
            sl2(s)     = sum(seg(:,s)**2)
            seglen(s)  = sqrt(sl2(s))
            tpath(s+1) = tpath(s) + seglen(s)
        end do
        !$omp parallel do default(shared) private(i,s,q,dz,t,d2,best_d2,best_c) &
        !$omp& schedule(static) proc_bind(close)
        do i = 1, nptcls
            best_d2 = huge(0.d0)
            best_c  = 0.d0
            do s = 1, nstates-1
                if( sl2(s) <= DTINY ) cycle
                do q = 1, nk
                    dz(q) = z(i,q) - tcen(q,s)
                end do
                t  = max(0.d0, min(1.d0, sum(dz*seg(:,s))/sl2(s)))
                d2 = 0.d0
                do q = 1, nk
                    d2 = d2 + (dz(q) - t*seg(q,s))**2
                end do
                if( d2 < best_d2 )then
                    best_d2 = d2
                    best_c  = tpath(s) + t*seglen(s)
                endif
            end do
            ppath(i) = best_c
        end do
        !$omp end parallel do
    end subroutine project_onto_target_polyline

    !>  Tied-covariance Gaussian-mixture responsibilities over the placed state targets. Replaces the
    !!  Epanechnikov kernel, whose compact support left many particles in no map at all; softmax over the
    !!  bandwidth is no substitute, it blurs every state back toward consensus. Tied covariance because
    !!  within-state spread is shared measurement error.
    !!  Measurements: doc/implementation_notes/flex_pca_state_placement_measurements.md
    subroutine gmm_state_weights( z, nptcls, ncomp, nk, nstates, tcen, wcomp, weights, neff, &
        &bandwidths, labels, views, orient_lam, precision, cdeconv )
        integer,  intent(in)    :: nptcls, ncomp, nk, nstates
        real(dp), intent(in)    :: z(nptcls,ncomp), wcomp(nk)
        !> in: placed targets. out: FITTED means, so the reported table describes the delivered maps
        real(dp), intent(inout) :: tcen(nk,nstates)
        real,     intent(inout) :: weights(nptcls,nstates), neff(nstates), bandwidths(nstates)
        integer,  intent(inout) :: labels(nptcls)
        !> per-particle viewing axis, enables the orientation-coverage term (see orient_lam)
        real(dp), optional, intent(in) :: views(3,nptcls)
        !> weight of the orientation-coverage penalty; 0 or absent = off
        real(dp), optional, intent(in) :: orient_lam
        !> per-particle latent PRECISION Pi_i; its inverse is the measurement covariance Ni
        real(dp), optional, intent(in) :: precision(ncomp,ncomp,nptcls)
        !> deconvolution scale c; <= 0 or absent = plain tied-covariance EM (shipped behaviour)
        real(dp), optional, intent(in) :: cdeconv
        real(dp), parameter :: GMM_REG = 1.d-6, GMM_TOL = 1.d-5
        !> responsibilities below this are zeroed so the reconstructor's live-state compaction works
        real(dp), parameter :: RESP_FLOOR = 1.d-3
        !> two means closer than this in the tied-covariance metric describe the same state
        real(dp), parameter :: GMM_MERGE_D2   = 1.0d0
        integer,  parameter :: GMM_MAX_RESPAWN = 8
        integer,  parameter :: GMM_MAXIT = 60
        !> ceiling on the deconvolution's per-particle covariance workspace
        real(dp), parameter :: GMM_DECONV_MAXGB = 4.d0
        !> the deconvolved fit needs far more iterations than the plain one; see the moment start
        integer,  parameter :: GMM_DECONV_MAXIT = 2000
        real(dp), allocatable :: y(:,:), mu(:,:), S(:,:), Sinv(:,:), Syy(:,:), resp(:,:)
        real(dp), allocatable :: Smu(:,:), mSm(:), pival(:), nresp(:), evwork(:,:), ev(:), evec(:,:)
        real(dp), allocatable :: ybar(:)
        real(dp), allocatable :: Tdev(:,:,:), Tbar(:,:)
        ! deconvolution workspace: Ni in the standardised frame, and the per-particle T_i^-1 / logdet
        ! rebuilt each iteration because they track S
        real(dp), allocatable :: Nz(:,:,:), Tinv_all(:,:,:), logdet_all(:)
        ! automatic, NOT allocatable: OpenMP-private allocatables arrive unallocated per thread, and
        ! allocatables are not portable array-reduction targets
        real(dp) :: cwork(ncomp,ncomp), Ti(nk,nk), Wi(nk,nk), Bi(nk,nk)
        real(dp) :: bik(nk), dev(nk)
        real(dp) :: Bsum(nk,nk), scat(nk,nk), mnum(nk,nstates), Nbar(nk,nk)
        real(dp) :: ySy, ySm, lmax, lsum, ll, prev_ll, logdet, lam, sumw, sumw2, vq, trS
        real(dp) :: bicval, ent, dmin, d2pair, cdec, gbytes, ldet, evfloor
        real(dp) :: trS_now, prev_trS
        integer  :: nfree, nrespawn, kmin, kdrop, iworst, ndegen, maxit_eff
        integer  :: i, q, r, state, it, nrot, errflg
        integer(kind=8) :: nact_tot
        logical  :: l_orient, l_deconv, l_conv
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
        ! Syy is fixed, so the tied-covariance M step is a rank-nstates correction, not a second pass:
        ! sum_k sum_i R_ik (y_i-mu_k)(y_i-mu_k)' = sum_i y_i y_i' - sum_k N_k mu_k mu_k', as sum_k R_ik = 1.
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
        ! ---- MEASUREMENT-ERROR DECONVOLUTION SETUP ----
        cdec     = 0.d0
        if( present(cdeconv) ) cdec = cdeconv
        l_deconv = cdec > 0.d0 .and. present(precision)
        if( l_deconv )then
            ! two nk x nk matrices per particle are carried across the EM loop
            gbytes = 2.d0*real(nk,dp)*real(nk,dp)*real(nptcls,dp)*8.d0 / 1.073741824d9
            if( gbytes > GMM_DECONV_MAXGB )then
                write(logfhandle,'(A,F7.2,A,F5.1,A)') &
                    &'>>> FLEX_PCA GMM deconvolution workspace ',gbytes,' GB exceeds the ', &
                    &GMM_DECONV_MAXGB,' GB budget'
                THROW_HARD('flex_pca deconvolution over budget; lower nkern or unset SIMPLE_COV_DECONV')
            endif
            allocate(Nz(nk,nk,nptcls), Tinv_all(nk,nk,nptcls), logdet_all(nptcls))
            ndegen = 0
            Nbar   = 0.d0
            !$omp parallel do default(shared) private(i,q,r,cwork,errflg) schedule(static) &
            !$omp& reduction(+:ndegen,Nbar)
            do i = 1, nptcls
                ! invert the FULL precision then slice: that is the MARGINAL posterior covariance.
                ! Inverting the nk block instead gives the conditional, which is systematically smaller.
                call matinv(precision(:,:,i), cwork, ncomp, errflg)
                if( errflg /= 0 )then
                    ! no usable noise model for this particle; contribute as if noiseless rather than
                    ! dropping it, and report the count so a wholesale failure cannot pass silently
                    Nz(:,:,i) = 0.d0
                    ndegen    = ndegen + 1
                    cycle
                endif
                do q = 1, nk
                    do r = 1, nk
                        Nz(q,r,i) = cwork(q,r) * sqrt(wcomp(q)*wcomp(r))
                        Nbar(q,r) = Nbar(q,r) + Nz(q,r,i)
                    end do
                end do
            end do
            !$omp end parallel do
            Nbar = Nbar / real(nptcls,dp)
            ! heldout=yes leaves latent_second EXACTLY zero for halfset A, so a run resumed from that
            ! cache would deconvolve half the data against a null noise model. Refuse.
            if( ndegen > nptcls/100 )then
                write(logfhandle,'(A,I0,A,I0,A)') '>>> FLEX_PCA GMM deconvolution: ',ndegen,' of ', &
                    &nptcls,' particles have no invertible latent precision'
                THROW_HARD('flex_pca deconvolution needs a normal-path precision cache (heldout=yes &
                    &leaves halfset A zero)')
            endif
            if( ndegen > 0 ) write(logfhandle,'(A,I0,A)') &
                &'>>> FLEX_PCA GMM deconvolution: ',ndegen,' particles fell back to a noiseless model'
            ! METHOD-OF-MOMENTS START. EM here contracts at a rate set by the noise-to-signal ratio
            ! and is unusably slow from the observed covariance -- at N/S = 44 it is still 2x from
            ! the fixed point after 200 iterations, and the likelihood plateaus long before S does,
            ! so it would look converged. Subtracting the mean noise first lands S on the answer for
            ! homoscedastic noise and near it otherwise.
            S = S - cdec*Nbar
            ! ... which can leave S indefinite when c overstates the noise, so floor its spectrum
            evwork = S
            call jacobi(evwork, nk, nk, ev, evec, nrot)
            call eigsrt(ev, evec, nk, nk)
            evfloor = max(GMM_REG, 1.d-3*max(ev(1), DTINY))
            do q = 1, nk
                ev(q) = max(ev(q), evfloor)
            end do
            do q = 1, nk
                do r = 1, nk
                    S(q,r) = sum(evec(q,:)*ev(:)*evec(r,:))
                end do
            end do
            S = 0.5d0*(S + transpose(S))
            write(logfhandle,'(A,F8.3,A,ES12.4)') &
                &'>>> FLEX_PCA GMM measurement-error deconvolution ACTIVE, c=',cdec, &
                &'  moment-matched initial tr(S)/nk=',sum([(S(q,q), q=1,nk)])/real(nk,dp)
            call flush(logfhandle)
        endif
        prev_trS  = -huge(1.d0)
        maxit_eff = merge(GMM_DECONV_MAXIT, GMM_MAXIT, l_deconv)
        do it = 1, maxit_eff
            call matinv(S, Sinv, nk, errflg)
            if( errflg /= 0 )then
                write(logfhandle,'(A)') '>>> FLEX_PCA GMM tied covariance singular; keeping kernel weights'
                deallocate(y, mu, S, Sinv, Syy, resp, Smu, mSm, pival, nresp, evwork, ev, evec)
                if( allocated(Tbar) ) deallocate(Tbar, Tdev)
                if( allocated(Nz)   ) deallocate(Nz, Tinv_all, logdet_all)
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
            if( l_deconv )then
                ! ---- E STEP (deconvolved) ----
                ! T_i = S + c*Ni does not depend on the component, so one SPD factorisation per
                ! particle serves all nstates. Cached in Tinv_all for the M step below rather than
                ! refactorised there.
                !$omp parallel do default(shared) schedule(static) reduction(+:ll) &
                !$omp& private(i,q,r,state,Ti,ldet,dev,ySy,lmax,lsum,vq,errflg)
                do i = 1, nptcls
                    Ti = S + cdec*Nz(:,:,i)
                    call dpotrf('U', nk, Ti, nk, errflg)
                    if( errflg /= 0 )then
                        ! T_i lost positive-definiteness; drop this particle's noise term for the
                        ! iteration rather than abandoning the fit
                        Ti = S
                        call dpotrf('U', nk, Ti, nk, errflg)
                    endif
                    ldet = 0.d0
                    do q = 1, nk
                        ldet = ldet + log(max(Ti(q,q), DTINY))
                    end do
                    ldet = 2.d0*ldet                        ! log|T| = 2 sum log diag(chol)
                    call dpotri('U', nk, Ti, nk, errflg)
                    do q = 1, nk                            ! dpotri fills one triangle only
                        do r = q+1, nk
                            Ti(r,q) = Ti(q,r)
                        end do
                    end do
                    Tinv_all(:,:,i) = Ti
                    logdet_all(i)   = ldet
                    do state = 1, nstates
                        do q = 1, nk
                            dev(q) = y(i,q) - mu(q,state)
                        end do
                        ySy = 0.d0
                        do q = 1, nk
                            do r = 1, nk
                                ySy = ySy + dev(q)*Ti(q,r)*dev(r)
                            end do
                        end do
                        resp(i,state) = -0.5d0*ySy - 0.5d0*ldet + log(max(pival(state), DTINY))
                        if( l_orient )then
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
                        lsum          = lsum + resp(i,state)
                    end do
                    resp(i,:) = resp(i,:) / max(lsum, DTINY)
                    ll = ll + (log(max(lsum, DTINY)) + lmax)
                end do
                !$omp end parallel do
            else
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
                        ! orientation coverage via v v', never the mean resultant -- a dipole statistic is blind to axis bias
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
            endif
            ll = ll / real(nptcls,dp)
            do state = 1, nstates
                nresp(state) = sum(resp(:,state))
            end do
            pival = max(nresp, DTINY) / real(nptcls,dp)
            if( l_deconv )then
                ! XD M step, Sigma tied: W_i and B_i are component-independent, and one pass suffices
                ! for the covariance by the same identity the Syy path uses,
                !   sum_ik R_ik (b_ik-mu_k)(b_ik-mu_k)' = Sbb - sum_k N_k mu_k mu_k'
                ! so the scatter never needs the updated means inside the loop.
                mnum = 0.d0; Bsum = 0.d0; scat = 0.d0
                !$omp parallel do default(shared) schedule(static) &
                !$omp& private(i,q,r,state,Wi,Bi,dev,bik) reduction(+:mnum,Bsum,scat)
                do i = 1, nptcls
                    Wi   = matmul(S, Tinv_all(:,:,i))
                    Bi   = S - matmul(Wi, S)
                    Bsum = Bsum + Bi
                    do state = 1, nstates
                        do q = 1, nk
                            dev(q) = y(i,q) - mu(q,state)
                        end do
                        do q = 1, nk
                            bik(q) = mu(q,state) + sum(Wi(q,:)*dev(:))
                        end do
                        do q = 1, nk
                            mnum(q,state) = mnum(q,state) + resp(i,state)*bik(q)
                            do r = 1, nk
                                scat(q,r) = scat(q,r) + resp(i,state)*bik(q)*bik(r)
                            end do
                        end do
                    end do
                end do
                !$omp end parallel do
                do state = 1, nstates
                    do q = 1, nk
                        mu(q,state) = mnum(q,state) / max(nresp(state), DTINY)
                    end do
                end do
                S = Bsum + scat
                do state = 1, nstates
                    do q = 1, nk
                        do r = 1, nk
                            S(q,r) = S(q,r) - nresp(state)*mu(q,state)*mu(r,state)
                        end do
                    end do
                end do
            else
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
            endif
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
                ! Converged EM parks components on one region and starves others. Two means within GMM_MERGE_D2 are
                ! the same state; keep one, restart the other at the worst-explained particle.
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
            &'>>> FLEX_PCA GMM tied-covariance responsibilities: iters=',min(it,maxit_eff), &
            &'  loglik=',ll,'  pi range ',minval(pival),' - ',maxval(pival)
        ! BIC will spend components on one populated state; ICL adds the entropy penalty and so prefers
        ! SEPARATED states. Reported per run; nothing is selected automatically yet.
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
        ! responsibility mass near a component's own dimension cannot be estimated -- a merge candidate
        write(logfhandle,'(A,I0,A,I0)') '>>> FLEX_PCA GMM components below 3*nk effective mass: ', &
            &count(nresp < 3.d0*real(nk,dp)),' of ',nstates
        call flush(logfhandle)
        if( l_orient ) write(logfhandle,'(A,F7.3)') &
            &'>>> FLEX_PCA GMM orientation-coverage term active, lambda=',lam
        ! SPARSIFY before handing over -- REQUIRED, not an optimisation. insert_plane_oversamp_multi_scaled
        ! assumes most scale pairs are EXACT zeros; left dense the state reconstruction does not finish.
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
        ! responsibilities ARE the weights: data and density scale both, so the per-state scale divides out
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
        do state = 1, nstates
            do q = 1, nk
                tcen(q,state) = mu(q,state) / sqrt(max(wcomp(q), DTINY))
            end do
        end do
        deallocate(y, mu, S, Sinv, Syy, resp, Smu, mSm, pival, nresp, evwork, ev, evec)
        if( allocated(Tbar) ) deallocate(Tbar, Tdev)
        if( allocated(Nz)   ) deallocate(Nz, Tinv_all, logdet_all)
    end subroutine gmm_state_weights

    subroutine kernel_weights_at_bandwidth( dist, nptcls, h_in, min_neff, w, h_out, neff_out )
        integer,  intent(in)  :: nptcls, min_neff
        real(dp), intent(in)  :: dist(nptcls), h_in
        real,     intent(out) :: w(nptcls)
        real(dp), intent(out) :: h_out
        real,     intent(out) :: neff_out
        real(dp) :: h, u2, sumw, sumw2
        integer  :: i, grow, nsupp
        h = h_in
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

    !> Cross-validated bandwidth selection, scored against the NARROWEST bin's opposite-half map.
    !! Plain even/odd agreement would rise monotonically with bandwidth and pick maximal smearing.
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

    !>  Real-valued environment override, leaving `val` untouched when unset. Companion to cov_env_int.
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

    !> Rotate the eigenbasis toward SPATIALLY COHERENT components: maximise each component's Rayleigh
    !! quotient against its smooth_lp-smoothed copy, as a generalized eigenproblem against the basis Gram.
    !! Latents and second moments are transformed with it, so U' z' reproduces U z exactly.
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
        integer :: ik, iq, ir, ii, ldim(3), nvox
        real(dp) :: acc, thr
        ok = .false.
        if( ncomp < 2 .or. smooth_lp <= 0. ) return
        pcdir = ''
        if( present(eigdir) ) pcdir = eigdir
        if( .not. file_exists(pcdir//'flex_pca_pc001.mrc') ) return
        call refvol%new([params%box_crop,params%box_crop,params%box_crop], params%smpd_crop)
        call refvol%read_and_crop(params%vols(1), params%smpd, params%box_crop, params%smpd_crop)
        rmat = refvol%get_rmat(); ldim = shape(rmat); nvox = product(ldim)
        allocate(vsort(nvox)); vsort = reshape(rmat, [nvox]); call hpsort(vsort)
        thr = real(vsort(max(1,min(nvox,nint(0.92*real(nvox))))), dp)
        allocate(msk(ldim(1),ldim(2),ldim(3)))
        msk = rmat > real(thr)
        deallocate(vsort)
        call refvol%kill
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

    !> Change the latent frame to one in which the population cloud is white.
    !!
    !! z' = Lam^-1/2 V^T z with V, Lam the eigenpairs of Cov(z). The basis absorbs the
    !! inverse, U' = U V Lam^1/2, so U' z' = U z exactly: this rotates and rescales the
    !! coordinates in which states are placed and leaves the model alone. Downstream that
    !! means every component contributes to distance in proportion to how much it separates
    !! particles rather than how much variance it happens to carry, which is what the
    !! diffusion graph and the mixture assignment both assume and neither currently gets.
    !!
    !! The transform is applied about the origin, not the latent mean, so the identity is
    !! exact. The PPCA prior is zero-mean and the measured offset is a few percent of a
    !! standard deviation; the mean is still used for the covariance, which is the geometry
    !! that matters. The eigenvolume amplitudes carry the variance after this, so pc001 no
    !! longer looks the same size as pc016 when it contributes forty times more.
    subroutine whiten_latent_frame( params, ncomp, z, nptcls, latent_second, ok, eigdir )
        class(parameters), intent(in)    :: params
        integer,           intent(in)    :: ncomp, nptcls
        real(dp),          intent(inout) :: z(nptcls,ncomp)
        real(dp),          intent(inout) :: latent_second(ncomp,ncomp,nptcls)
        logical,           intent(out)   :: ok
        character(len=*), optional, intent(in) :: eigdir
        character(len=:), allocatable :: pcdir
        type(image),  allocatable :: eigs(:), sm(:)
        real(dp),     allocatable :: C(:,:), V(:,:), ev(:), Rot(:,:), Rinv(:,:), zrow(:), mu(:)
        type(string) :: fn
        real(dp) :: floor_ev, cond_before, off_max
        integer  :: ii, iq, ir, ik, nrot
        ok = .false.
        if( ncomp < 2 .or. nptcls < 10 ) return
        allocate(mu(ncomp), C(ncomp,ncomp), V(ncomp,ncomp), ev(ncomp))
        do iq = 1, ncomp
            mu(iq) = sum(z(:,iq)) / real(nptcls,dp)
        end do
        C = 0.d0
        do ii = 1, nptcls
            do iq = 1, ncomp
                C(:,iq) = C(:,iq) + (z(ii,:) - mu) * (z(ii,iq) - mu(iq))
            end do
        end do
        C = C / real(nptcls,dp)
        C = 0.5d0*(C + transpose(C))
        ! report the conditioning being removed, so a run where this is a no-op says so
        off_max = 0.d0
        do iq = 1, ncomp
            do ir = 1, ncomp
                if( iq == ir ) cycle
                off_max = max(off_max, abs(C(iq,ir))/sqrt(max(C(iq,iq)*C(ir,ir),DTINY)))
            end do
        end do
        call jacobi(C, ncomp, ncomp, ev, V, nrot)
        call eigsrt(ev, V, ncomp, ncomp)             ! descending variance
        if( ev(1) <= DTINY )then
            write(logfhandle,'(A)') '>>> FLEX_PCA LATENT WHITENING skipped: latent covariance has no &
                &positive spectrum'
            return
        endif
        ! a component the fit never populated must not be inflated into the geometry
        floor_ev    = 1.d-8 * ev(1)
        cond_before = ev(1) / max(ev(ncomp), floor_ev)
        ev = max(ev, floor_ev)
        allocate(Rot(ncomp,ncomp), Rinv(ncomp,ncomp))
        do iq = 1, ncomp
            Rot(:,iq)  = V(:,iq) * sqrt(ev(iq))      ! R   = V Lam^1/2
            Rinv(iq,:) = V(:,iq) / sqrt(ev(iq))      ! R^-1 = Lam^-1/2 V^T
        end do
        allocate(zrow(ncomp))
        do ii = 1, nptcls
            zrow = z(ii,:)
            do iq = 1, ncomp
                z(ii,iq) = sum(zrow * Rinv(iq,:))
            end do
        end do
        ! precision transforms as Pi' = R^T Pi R, matching the smoothness rotation
        !$omp parallel do schedule(static) default(shared) private(ii)
        do ii = 1, nptcls
            latent_second(:,:,ii) = matmul(transpose(Rot), matmul(latent_second(:,:,ii), Rot))
        end do
        !$omp end parallel do
        ! the eigenvolumes absorb R, keeping U' z' = U z; skipped when they are not on disk,
        ! in which case the latent frame still changes and only the written pc maps go stale
        pcdir = ''
        if( present(eigdir) ) pcdir = eigdir
        if( file_exists(pcdir//'flex_pca_pc001.mrc') )then
            allocate(eigs(ncomp), sm(ncomp))
            do ik = 1, ncomp
                call eigs(ik)%new([params%box_crop,params%box_crop,params%box_crop], params%smpd_crop)
                fn = pcdir//'flex_pca_pc'//int2str_pad(ik,3)//MRC_EXT
                call eigs(ik)%read(fn); call fn%kill
                call sm(ik)%copy(eigs(ik))
            end do
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
                call eigs(ik)%kill; call sm(ik)%kill
            end do
            deallocate(eigs, sm)
        else
            write(logfhandle,'(A)') '>>> FLEX_PCA LATENT WHITENING: eigenvolumes not found, pc maps &
                &left in the old frame (state placement is unaffected)'
        endif
        write(logfhandle,'(A,ES10.3,A,F6.3,A,ES10.3,A,ES10.3)') &
            &'>>> FLEX_PCA LATENT WHITENING applied: condition number ',cond_before, &
            &' -> 1, max |corr| ',off_max,'; latent variance max=',ev(1),' min=',ev(ncomp)
        call flush(logfhandle)
        deallocate(C, V, ev, Rot, Rinv, zrow, mu)
        ok = .true.
    end subroutine whiten_latent_frame

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

    !> FINCH state placement: one representative per cluster, so the state COUNT comes from the data.
    !! Runs on a standardized strided subsample.
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
        write(logfhandle,'(A)',advance='no') '>>> FLEX_PCA FINCH hierarchy (clusters per level): '
        do i = 1, hier%get_nlevels()
            write(logfhandle,'(I0,1X)',advance='no') hier%get_nclusters(i)
        end do
        write(logfhandle,*)
        ! SIMPLE_COV_FINCH_LEVEL=N takes that level natively -- no Ward merge, no cap, FINCH's own count
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

    !> Reliability proxy from a cached embedding: observed spread over mean posterior variance, mapped
    !! through r/(1+r) onto the split-half rho scale. Posterior variances are stride-sampled.
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

    !> Manifold-covering state targets by diffusion-map k-center. k-means allocates by DENSITY and misses
    !! sparse states; a 1-D path merges states on a BRANCHED manifold. A diffusion embedding makes geodesic
    !! structure Euclidean, so greedy farthest-point covers any shape. Reliability weighting is essential:
    !! the largest-eigenvalue component is typically the worst measured.
    !! Measurements: doc/implementation_notes/flex_pca_state_placement_measurements.md
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
        real(dp) :: d2, s, nrm, best, wk1
        integer  :: nnode, i, j, q, e, nedge, it, m, c, ibest, ni
        ok = .false.
        if( nstates < 2 .or. nptcls < 100 ) return
        nnode = min(nptcls, NNODE_MAX)
        if( nnode <= KNN + 2 ) return
        m = NDIFF + 1                                  ! + the trivial eigenvector
        allocate(nodes(nnode), rw(ncomp), zbar(ncomp), sdv(ncomp))
        ! deterministic stride subsample: it already carries the data's density, and is reproducible run to run
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
        ! self-tuning bandwidth: sigma_i is the K-th neighbour distance, so affinity adapts to local density
        allocate(sig(nnode))
        do i = 1, nnode
            sig(i) = sqrt(max(real(knntab%distance2(KNN,i),dp), DTINY))
        end do
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
        ! alpha = 1 (Laplace-Beltrami): divide out the sampling density so the embedding reflects manifold
        ! GEOMETRY, not how heavily each region is populated -- the exact failure mode of k-means here.
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
        ! psi = D^-1/2 V, dropping the trivial leading eigenvector; the 1/sqrt(1-lambda) commute-time
        ! scaling puts Euclidean distance on the diffusion metric, which is what makes k-center meaningful.
        allocate(psi(nnode,NDIFF))
        wk1 = 1.d0 / sqrt(max(1.d0 - min(lam(2), 1.d0-1.d-9), 1.d-9))
        do j = 1, NDIFF
            s = 1.d0 / sqrt(max(1.d0 - min(lam(j+1), 1.d0-1.d-9), 1.d-9))
            nrm = sqrt(max(sum((V(:,j+1)*ddeg)**2), DTINY))
            do i = 1, nnode
                psi(i,j) = V(i,j+1)*ddeg(i)/nrm * (s/wk1)
            end do
        end do
        ! greedy k-center: seed farthest from the centroid, then farthest from everything chosen -- coverage
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
        ! each node takes the nearest selected node's CELL; the target is the cell's mean RAW latent, which
        ! regresses the RIM nodes k-center picks by design back onto the state they cover.
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

    !> Equal-occupancy targets along a reliability-ordered principal direction. Two departures from
    !! k-means: the direction is the RELIABILITY-weighted leading eigenvector, so a high-variance but
    !! poorly measured nuisance mode cannot set it; and slices carry equal PARTICLE COUNTS, not width.
    subroutine reliability_path_targets( z, nptcls, ncomp, nstates, wcomp, rho, centroids, proj_out, tproj_out )
        integer,  intent(in)  :: nptcls, ncomp, nstates
        real(dp), intent(in)  :: z(nptcls,ncomp), wcomp(ncomp), rho(ncomp)
        real(dp), intent(out) :: centroids(ncomp,nstates)
        ! coordinate ALONG the path. A full-rank Mahalanobis kernel also measures the off-path directions, so
        ! an on-path but noisy particle falls outside every support -- that is what strands the dataset.
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
        ! reliability RELATIVE to the best-measured component, floored so no direction is removed outright
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
        ! ties can leave a slice empty; interpolate so the polyline stays ordered and no state hits the mean
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
        real(dp) :: d2, dmax, dnorm, lo, hi, t
        integer  :: i, q, s, ia, ib, islot
        allocate(zbar(ncomp), pa(ncomp), pb(ncomp), dirv(ncomp), proj(nptcls), &
            &sproj(nptcls), cnt(nstates))
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
                t = real(s-1,dp)/real(max(1,nstates-1),dp)
                centroids(:,s) = pa + t*dirv
            endif
        end do
        write(logfhandle,'(A,I0,A,I0,A,I0)') '>>> FLEX_PCA path endpoints from particles ',ia,' and ',ib, &
            &'; slice occupancies min=',minval(cnt)
        deallocate(zbar, pa, pb, dirv, proj, sproj, cnt)
    end subroutine path_latent_targets

    !> Deterministic k-means in the SAME wcomp metric the kernel uses, so placement and weighting agree.
    !! Seeded farthest-point from the particle nearest the latent mean, so no RNG is involved.
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

    subroutine heldout_embedding( params, build, mean_rec, pinds, nptcls, col_sep, neigs_req, &
        &basis_recs, eigvals, ncomp, sig2_eff, z, contrast, latent_second, resid_energy, resid_mean_energy )
        type(parameters),    intent(inout) :: params
        type(builder),       intent(inout) :: build
        type(reconstructor), intent(inout) :: mean_rec
        integer,             intent(in)    :: pinds(:), nptcls, col_sep, neigs_req
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
            &col_sep, neigs_req, recs_a, eig_a, ncomp_a, sig2_a, &
            &basis_imgs=imgs_a, fprefix='flex_pca_pc')
        call build_covariance_eigenbasis(params, build, mean_rec, pind_b, nb, &
            &col_sep, neigs_req, recs_b, eig_b, ncomp_b, sig2_b, &
            &basis_imgs=imgs_b, fprefix='flex_pca_heldout_pcB')
        ncomp = ncomp_a
        if( ncomp < 1 .or. ncomp_b < 1 ) THROW_HARD('flex_pca heldout produced no retained components')
        write(logfhandle,'(A,I0,A,I0)') '>>> FLEX_PCA heldout components A(reference)=',ncomp_a, &
            &' B=',ncomp_b
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
        ! merge into the caller's order, rotating the B-frame latents into basis A's reference frame
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

    !> Bagged-basis embedding: fit two bases on disjoint halfsets, pool them, embed everything.
    !!
    !!  NOT a held-out estimate: the pooled basis is fitted on ALL particles, so heldout=yes remains
    !!  the route for unbiased cross-halfset diagnostics. Use heldout to MEASURE, this to PRODUCE.
    subroutine bagged_embedding( params, build, mean_rec, pinds, nptcls, col_sep, neigs_req, &
        &basis_recs, eigvals, ncomp, sig2_eff, z, contrast, latent_second, resid_energy, resid_mean_energy )
        type(parameters),    intent(inout) :: params
        type(builder),       intent(inout) :: build
        type(reconstructor), intent(inout) :: mean_rec
        integer,             intent(in)    :: pinds(:), nptcls, col_sep, neigs_req
        type(reconstructor), allocatable, intent(out) :: basis_recs(:)
        real(dp),            allocatable, intent(out) :: eigvals(:)
        integer,             intent(out)   :: ncomp
        real(dp),            intent(out)   :: sig2_eff
        real(dp),            allocatable, intent(out) :: z(:,:), contrast(:), latent_second(:,:,:)
        real(dp),            allocatable, intent(out) :: resid_energy(:), resid_mean_energy(:)
        type(reconstructor), allocatable :: recs_a(:), recs_b(:)
        type(image),         allocatable :: imgs_a(:), imgs_b(:), pooled(:)
        real(dp),            allocatable :: eig_a(:), eig_b(:), M(:,:), svals(:)
        integer,             allocatable :: pind_a(:), pind_b(:)
        type(string) :: fname
        integer  :: i, q, na, nb, ncomp_a, ncomp_b, ia, ib
        real(dp) :: sig2_a, sig2_b
        ! DISTRIBUTION. build_covariance_eigenbasis distributes internally -- as master it spawns a
        ! qsys round and reduces the part files, as worker it runs its one assigned stage. Bagging
        ! calls it twice, so the stage alone is ambiguous: it says WHAT to compute but not over which
        ! particles. flex_pca_set_fit stamps each round with the fit it serves, the worker reads it
        ! back from params%pcafit, and the select case below routes it to the matching halfset.
        ! Every process partitions identically because get_eo is a project property, so master and
        ! worker agree on the split without any of it being communicated.
        ! partition the sampled particles by halfset
        na = 0; nb = 0
        do i = 1, nptcls
            if( build%spproj_field%get_eo(pinds(i)) == 0 )then
                na = na + 1
            else
                nb = nb + 1
            endif
        end do
        if( na < 50 .or. nb < 50 ) THROW_HARD('flex_pca bagging requires >=50 particles per halfset')
        allocate(pind_a(na), pind_b(nb))
        ia = 0; ib = 0
        do i = 1, nptcls
            if( build%spproj_field%get_eo(pinds(i)) == 0 )then
                ia = ia + 1; pind_a(ia) = pinds(i)
            else
                ib = ib + 1; pind_b(ib) = pinds(i)
            endif
        end do
        write(logfhandle,'(A,I0,A,I0)') '>>> FLEX_PCA BAGGED BASIS halfset_A=',na,' halfset_B=',nb
        call flush(logfhandle)
        ! ---- WORKER: serve the one round the master scheduled us for ----
        if( flex_pca_is_worker() )then
            select case( flex_pca_fit() )
                case( FLEX_FIT_A )
                    call build_covariance_eigenbasis(params, build, mean_rec, pind_a, na, &
                        &col_sep, neigs_req, recs_a, eig_a, ncomp_a, sig2_a)
                    call release_recs(recs_a)
                case( FLEX_FIT_B )
                    call build_covariance_eigenbasis(params, build, mean_rec, pind_b, nb, &
                        &col_sep, neigs_req, recs_b, eig_b, ncomp_b, sig2_b)
                    call release_recs(recs_b)
                case default
                    ! FLEX_FIT_ALL: pooled-basis rounds. The basis is on disk under the normal
                    ! flex_pca_pc* names by then, so the full particle range applies.
                    call build_covariance_eigenbasis(params, build, mean_rec, pinds, nptcls, &
                        &col_sep, neigs_req, recs_a, eig_a, ncomp_a, sig2_a)
                    call release_recs(recs_a)
            end select
            ! the part file is on disk; run_flex_pca calls qsys_job_finished and exits
            ncomp    = 0
            sig2_eff = 0.d0
            allocate(basis_recs(0), eigvals(0))
            allocate(z(0,0), contrast(0), latent_second(0,0,0))
            allocate(resid_energy(0), resid_mean_energy(0))
            if( allocated(eig_a) ) deallocate(eig_a)
            if( allocated(eig_b) ) deallocate(eig_b)
            deallocate(pind_a, pind_b)
            return
        endif
        ! ---- MASTER: two fits, each stamped so its workers know which halfset they serve ----
        call flex_pca_set_fit(FLEX_FIT_A)
        call build_covariance_eigenbasis(params, build, mean_rec, pind_a, na, &
            &col_sep, neigs_req, recs_a, eig_a, ncomp_a, sig2_a, &
            &basis_imgs=imgs_a, fprefix='flex_pca_bagA_pc')
        call flex_pca_set_fit(FLEX_FIT_B)
        call build_covariance_eigenbasis(params, build, mean_rec, pind_b, nb, &
            &col_sep, neigs_req, recs_b, eig_b, ncomp_b, sig2_b, &
            &basis_imgs=imgs_b, fprefix='flex_pca_bagB_pc')
        ! back to the full range before anything downstream schedules a round
        call flex_pca_set_fit(FLEX_FIT_ALL)
        if( ncomp_a < 1 .or. ncomp_b < 1 ) THROW_HARD('flex_pca bagging produced no retained components')
        ncomp = ncomp_a
        write(logfhandle,'(A,I0,A,I0,A,I0)') '>>> FLEX_PCA bag components A=',ncomp_a,' B=',ncomp_b, &
            &' pooled rank=',ncomp
        ! how much agreement there was to bag, on the same statistic heldout=yes reports
        call align_basis_to_reference(imgs_a, ncomp_a, imgs_b, ncomp_b, M, svals)
        write(logfhandle,'(A)') '>>> FLEX_PCA bag input principal-angle cosines (disjoint fits):'
        do q = 1, size(svals)
            write(logfhandle,'(A,I3,A,F8.4)') '>>>   k=',q,'  cos(theta)=',real(svals(q))
        end do
        call flush(logfhandle)
        call bag_basis_pool(imgs_a, ncomp_a, eig_a(:ncomp_a), imgs_b, ncomp_b, eig_b(:ncomp_b), &
            &ncomp, pooled, eigvals)
        ! the delivered basis takes over the normal-path on-disk names, so probe_external_basis,
        ! load_probe_basis and the analysis scripts all keep working unchanged
        do q = 1, ncomp
            fname = string('flex_pca_pc')//int2str_pad(q,3)//MRC_EXT
            call pooled(q)%write(fname)
            call fname%kill
        end do
        call basis_recs_from_images(params, build, pooled, ncomp, basis_recs)
        sig2_eff = 0.5d0*(sig2_a + sig2_b)
        ! every particle through the pooled basis
        allocate(z(nptcls,ncomp), contrast(nptcls), latent_second(ncomp,ncomp,nptcls))
        allocate(resid_energy(nptcls), resid_mean_energy(nptcls))
        if( flex_pca_is_master() .and. flex_pca_nparts() > 1 )then
            ! workers rebuild the pooled basis from the volumes written above and ship statistics
            call save_probe_state(ncomp, eigvals, sig2_eff)
            call flex_pca_run_stage(PCA_STAGE_EMBED, 'bagged embedding')
            call embed_latents_with_contrast(params, build, mean_rec, basis_recs, ncomp, eigvals, &
                &sig2_eff, pinds, nptcls, z, contrast, latent_second, resid_energy, &
                &resid_mean_energy, from_parts=.true.)
        else
            call embed_latents_with_contrast(params, build, mean_rec, basis_recs, ncomp, eigvals, &
                &sig2_eff, pinds, nptcls, z, contrast, latent_second, resid_energy, resid_mean_energy)
        endif
        do q = 1, ncomp_a
            call recs_a(q)%dealloc_rho; call recs_a(q)%kill; call imgs_a(q)%kill
        end do
        do q = 1, ncomp_b
            call recs_b(q)%dealloc_rho; call recs_b(q)%kill; call imgs_b(q)%kill
        end do
        do q = 1, ncomp
            call pooled(q)%kill
        end do
        deallocate(recs_a, recs_b, imgs_a, imgs_b, pooled, eig_a, eig_b, M, svals)
        deallocate(pind_a, pind_b)
      contains
        !> a worker's build_covariance_eigenbasis returns an empty or short-lived basis; free it
        subroutine release_recs( recs )
            type(reconstructor), allocatable, intent(inout) :: recs(:)
            integer :: s
            if( .not. allocated(recs) ) return
            do s = 1, size(recs)
                call recs(s)%dealloc_rho; call recs(s)%kill
            end do
            deallocate(recs)
        end subroutine release_recs
    end subroutine bagged_embedding

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
        integer :: u, q, ln, stat, nread
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

    !> Kernel contract: compact support, peak 1, neff bounded by the support count, capped widening.
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
        call kernel_weights_at_bandwidth(dist, NP, 10.d0, 200, w, h_wide, neff_wide)
        if( h_wide <= 10.d0 ) THROW_HARD('bandwidth did not widen when support fell short of min_neff')
        if( h_wide > 10.d0*1.3d0**COV_MAX_BW_GROW + 1.d-9 ) THROW_HARD('bandwidth widening exceeded its cap')
        if( count(w > 0.) <= nsupp ) THROW_HARD('widening did not increase support')
        if( neff_wide <= neff ) THROW_HARD('neff did not increase with bandwidth')
        write(logfhandle,'(A)') '>>>   PASSED (support, normalisation, neff bounds, bounded widening)'
    end subroutine test_flex_pca_kernel_bandwidth

    subroutine test_flex_pca_state_weights()
        integer,  parameter :: NPC = 150, NP = 2*NPC, NC = 2, NST = 2
        integer  :: i, q, state, nlab(NST)
        real(dp) :: z(NP,NC), eigvals(NC), prec(NC,NC,NP)
        real,     allocatable :: weights(:,:), targets(:,:), bandwidths(:), neff(:)
        integer,  allocatable :: labels(:)
        write(logfhandle,'(A)') '>>> TEST flex_pca state weights on a two-cluster embedding'
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

    !> deterministic centred unit-variance pseudo-random draw, so the tests do not depend on an RNG
    real(dp) function punit( n ) result( u )
        integer, intent(in) :: n
        real(dp) :: t
        t = sin(real(n,dp)*12.9898d0)*43758.5453d0
        u = t - floor(t)                               ! uniform(0,1)
        u = (2.d0*u - 1.d0)*sqrt(3.d0)                 ! centred, unit variance
    end function punit
end module simple_flex_pca_model
