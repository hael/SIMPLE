!@descr: Standalone projection-aware low-rank covariance workflow for heterogeneous SPA data
module simple_flex_pca_model
use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
use simple_core_module_api
use simple_builder,                    only: builder
use simple_cmdline,                    only: cmdline
use simple_flex_diffmap_rec3D,         only: reconstruct_flex_diffmap_weighted_states, flex_rec_box, flex_rec_smpd
use simple_flex_pca_columns,    only: cov_env_int_pub, build_covariance_eigenbasis, embed_latents_with_contrast, &
    &estimate_covariance_mean, probe_subspace_iteration, align_basis_to_reference, probe_external_basis
use simple_image,                      only: image
use simple_parameters,                 only: parameters
use simple_reconstructor,              only: reconstructor
use simple_sigma2_files,               only: load_sigma2_groups
use simple_srch_sort_loc,              only: hpsort
use simple_finch,                      only: finch_hierarchy, fit_finch, finch_representatives, &
    &                                          select_finch_level, refine_finch_level
use simple_linalg,                     only: jacobi, eigsrt, matinv
implicit none
private
#include "simple_local_flags.inc"

public :: run_flex_pca
public :: test_flex_pca_embedding_cache_io, test_flex_pca_kernel_bandwidth
public :: test_flex_pca_state_weights

character(len=8), parameter :: COV_CACHE_MAGIC   = 'SIMPLFXC'
! A*Gtil^+*A, not the posterior precision A.
integer,          parameter :: COV_CACHE_VERSION = 3
! Safety cap on the bandwidth widening loop.
integer,          parameter :: COV_MAX_BW_GROW   = 4

contains

    subroutine run_flex_pca( params, build, cline )
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
        integer :: q, i, r, metric_valid_count, axis_sel, nfinch
        real(dp), allocatable :: finch_targets(:,:)
        logical :: l_finch_states
        logical :: l_geo, l_rot
        real(dp), allocatable :: kdist(:,:), kfloor(:)
        character(len=:), allocatable :: cachedir, cachestr
        real(dp) :: sig2_eff
        logical :: sigma_loaded, l_resume

        call validate_covariance_inputs(params, build, cline, pinds, nptcls)
        call load_and_validate_sigma(params, build, cline, pinds, sigma_loaded)

        neigs_req  = max(1, min(48, params%neigs))
        nstates    = max(3, min(20, params%npreimages))
        min_neff   = max(20, min(nptcls, params%min_neff))
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
        ! A low-pass finer than the working Nyquist (2*smpd_crop) cannot be honoured:
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

        ! RESUME MODE. The covariance basis and the per-particle embedding are ~77 % of the runtime (SNR 134
        ! s + columns 584 s + reduced solve 117 s + embed 35 s of 1126 s at box 128 / 100 k) and are
        ! COMPLETELY INDEPENDENT of the state count, the kernel bandwidth and the target placement.
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

        ! Mean of eq.
        call estimate_covariance_mean(params, build, mean_rec, pinds, nptcls)

        if( params%l_heldout )then
            ! Cross-halfset embedding:
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
            if( state_axis > 0 ) state_axis = min(state_axis, min(ncomp, nkern))
            write(logfhandle,'(A,I0)') '>>> FLEX_PCA retained covariance components=',ncomp
            call flush(logfhandle)

            ! Phase 3 (proposal 4.1/4.2): alternating the Wiener E-step (per-particle latents, with the
            ! anisotropic posterior second moment) and the weighted-backprojection M-step cleans the
            ! per-particle projection directions -- which is the one thing that moves the state maps off the
            ! null (the covariance columns give per-voxel-weighted directions with ~99% noise energy
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

            ! MAP embedding with the covariance prior:
            do q = 1, ncomp
                prior_precision(q) = 1.d0 / max(eigvals(q), DTINY)
            end do
            write(logfhandle,'(A,ES12.4,A,ES12.4)') '>>> FLEX_PCA covariance eigenvalues: max=', &
                &maxval(eigvals),' min=',minval(eigvals)
            call flush(logfhandle)
            ! Contrast-aware MAP embedding (S.D/S.E):
            allocate(contrast(nptcls))
            call embed_latents_with_contrast(params, build, mean_rec, basis_recs, ncomp, eigvals, sig2_eff, &
                &pinds, nptcls, z, contrast, latent_second, resid_energy, resid_mean_energy)
        endif
        write(logfhandle,'(A,F7.3,A,F7.3)') '>>> FLEX_PCA per-particle contrast: mean=', &
            &real(sum(contrast)/real(nptcls,dp)),' sd=', &
            &real(sqrt(max(sum((contrast-sum(contrast)/nptcls)**2)/real(nptcls,dp),DTINY)))
        call flush(logfhandle)
        ! EXTERNAL-BASIS PROBE (SIMPLE_COV_PROBEBASIS=<prefix>, SIMPLE_COV_PROBEN=<count>, optionally
        ! SIMPLE_COV_PROBEEIG=<dir> when resuming and the eigenvolumes live elsewhere).


        ! NO latent standardization. the reference reports z in the physical units the MAP solve (eq.

        call write_covariance_eigenvolumes(basis_recs, eigvals, ncomp)
        call write_embedding_cache('flex_pca_embedding.bin', pinds, nptcls, ncomp, z, eigvals, &
            &contrast, resid_energy, resid_mean_energy, latent_second, sig2_eff)

        endif   ! .not. l_resume
        ! EXTERNAL-BASIS PROBE (SIMPLE_COV_PROBEBASIS=<prefix>, SIMPLE_COV_PROBEN=<count>, optionally
        ! SIMPLE_COV_PROBEEIG=<dir>).
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
        ! FINCH state placement (SIMPLE_COV_FINCH_STATES=1):
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
        if( l_finch_states )then
            call build_covariance_state_weights(z, nptcls, ncomp, nkern, nstates, state_axis, min_neff, &
                &eigvals, latent_second, state_weights, targets, bandwidths, neff, labels, &
                &dist_out=kdist, bfloor_out=kfloor, targets_in=finch_targets)
            deallocate(finch_targets)
        else
            call build_covariance_state_weights(z, nptcls, ncomp, nkern, nstates, state_axis, min_neff, &
                &eigvals, latent_second, state_weights, targets, bandwidths, neff, labels, &
                &dist_out=kdist, bfloor_out=kfloor)
        endif
        ! ---- STATE-RECONSTRUCTION OPERATOR SETUP ---- This MUST precede cv_select_bandwidths, not follow
        ! it.
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

        ! Reuse the ordinary gridding reconstructor three times.
        params%outvol = 'flex_pca_state_001.mrc'
        call reconstruct_flex_diffmap_weighted_states(params, build, pinds, state_weights, nstates, floor_rho=.true.)
        allocate(half_weights(nptcls,nstates), source=state_weights)
        call mask_state_weights_by_half(build, pinds, 0, half_weights)
        params%outvol = 'flex_pca_even_state_001.mrc'
        call reconstruct_flex_diffmap_weighted_states(params, build, pinds, half_weights, nstates, floor_rho=.true.)
        half_weights = state_weights
        call mask_state_weights_by_half(build, pinds, 1, half_weights)
        params%outvol = 'flex_pca_odd_state_001.mrc'
        call reconstruct_flex_diffmap_weighted_states(params, build, pinds, half_weights, nstates, floor_rho=.true.)
        ! NOTE (transport geometry, measured on IgG-RL): the linearised optimal-transport embedding of
        ! this state set needs 3 modes where density-space PCA needs 13, and its second mode correlates
        ! 0.927 with the true Fab swing against 0.556 for the best density component. The descriptive
        ! pass that produced those numbers is retired; simple_flex_lot now exposes only
        ! lot_pullback_metric, which can supply the optional zmetric of build_covariance_state_weights.
        ! Order the reconstructed states along the dominant conformational motion in VOLUME space so the
        ! trajectory endpoints span the actual change (the covariance state axis need not order the
        ! states along the motion).
        axis_sel = 0
        allocate(comp_rank(ncomp))
        do q = 1, ncomp
            comp_rank(q) = q
        end do
        if( file_exists('flex_pca_pc001.mrc') )then
            call order_states_by_volume_trajectory(params, nstates, ncomp, axis_sel=axis_sel, comp_rank=comp_rank)
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
            ! Preferred:
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
                    &targets_in=geo_targets)
                deallocate(geo_targets)
            else
                write(logfhandle,'(A,I0,A)') '>>> FLEX_PCA placing trajectory targets along component ', &
                    &axis_sel,' (selected by half-set reproducibility) -> flex_pca_traj_###.mrc'
                call build_covariance_state_weights(z, nptcls, ncomp, nkern, nstates, axis_sel, min_neff, &
                    &eigvals, latent_second, traj_weights, traj_targets, traj_bandwidths, traj_neff, traj_labels)
            endif
            params%outvol = 'flex_pca_traj_001.mrc'
            call reconstruct_flex_diffmap_weighted_states(params, build, pinds, traj_weights, nstates, floor_rho=.true.)
            call write_trajectory_targets(axis_sel, nstates, ncomp, traj_targets, traj_bandwidths, traj_neff)
            deallocate(traj_weights, traj_targets, traj_bandwidths, traj_neff, traj_labels)
        endif
        ! Nonuniform filtering LAST:
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
        ! Assign gold-standard halfsets by index parity if the project's eo split is degenerate (e.g.
        if( count([(build%spproj_field%get_eo(pinds(q))==0,q=1,nptcls)]) < 20 .or. &
            &count([(build%spproj_field%get_eo(pinds(q))==1,q=1,nptcls)]) < 20 )then
            write(logfhandle,'(A)') '>>> FLEX_PCA assigning alternating even/odd halfsets (project eo was degenerate)'
            call build%spproj_field%partition_eo
        endif
        if( count([(build%spproj_field%get_eo(pinds(q))==0,q=1,nptcls)]) < 20 .or. &
            &count([(build%spproj_field%get_eo(pinds(q))==1,q=1,nptcls)]) < 20 ) &
            &THROW_HARD('flex_pca requires populated even and odd halfsets')
    end subroutine validate_covariance_inputs

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
    end subroutine load_and_validate_sigma

    !> Median of the chi-squared distribution with k degrees of freedom, via the Wilson-Hilferty cube-root
    !! normal approximation evaluated at the median (z=0): chi2_med(k) ~= k * (1 - 2/(9k))^3 calls scipy's
    !! chi2.ppf(0.5, df=zdim) for the same quantity (latent_density.get_log_likelihood_threshold). The
    !! approximation is exact to 0.02 % at k=20, 0.7 % at k=3 and 3 % at k=1 -- far inside the tolerance of
    !! a bandwidth FLOOR, and it avoids pulling in an incomplete-gamma inverse.
    pure real(dp) function chi2_median( k )
        integer, intent(in) :: k
        real(dp) :: kk
        kk = real(max(k,1),dp)
        chi2_median = kk * (1.d0 - 2.d0/(9.d0*kk))**3
    end function chi2_median

    !> Binary cache of everything the state-weight and reconstruction stages need, so a different state
    !! count / bandwidth / placement can be tried without re-fitting the covariance basis and re-embedding
    !! 100 k particles.
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
        ! NOTE:
        call del_file('flex_pca_eigenvalues.txt')
        open(newunit=u,file='flex_pca_eigenvalues.txt',status='replace',action='write')
        write(u,'(A)') '# component eigenvalue'
        do q = 1, ncomp
            write(u,'(I6,1X,ES20.10)') q,eigvals(q)
        end do
        close(u)
    end subroutine write_covariance_eigenvolumes

    !> reference-style kernel-regression reconstruction weights (supplement S.F). Measured on IgG-RL with 20
    !! informative components it left 92-96% of the particles on a single state and gave >100 particles to
    !! only 8 of 15 states, and choosing a better axis does not repair it: component 1 correlated 0.006 with
    !! the true motion while component 2 gave 0.927, yet forcing axis=2 moved max occupancy only 92.7% ->
    !! 96.3%.
    subroutine build_covariance_state_weights( z, nptcls, ncomp, nkern, nstates, axis, min_neff, &
        &eigvals, precision, weights, targets, bandwidths, neff, labels, dist_out, bfloor_out, targets_in, zmetric )
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
        real,     allocatable :: sorted(:)
        real(dp), allocatable :: wcomp(:), tvec(:), tcen(:,:), dist(:), dvec(:), mvec(:)
        real(dp), allocatable :: pk(:,:,:), cfull(:,:), cblk(:,:), edges(:)
        integer,  allocatable :: occ(:)
        real(dp) :: qlo, qhi, target, h, d2, u2, sumw, sumw2, best, zspread, bmin, chi2med, wsum_i
        integer  :: nrenorm
        integer  :: i, q, r, state, ilo, ihi, best_state, grow, nfed, occmax, ifloor, nunassigned, nsupp
        integer  :: nk, errflg
        ! Per-component weights for TARGET PLACEMENT (kmeans_latent_targets / path_latent_targets). comp
        ! share of the eigenvalue metric conformational information z1 57.8 % 3.56x (rank 19 of 20) z2 18.9
        ! % 2.86x (rank 20 of 20) z4 3.2 % 14.54x (rank 1 of 20) z1+z2 held 76.7 % of the placement metric
        ! while ranking LAST, and corr(log eigenvalue, conformational information) = -0.147. Targets were
        ! therefore placed along the two least informative directions, so every state drew a
        ! conformationally near-random particle mixture and all state maps collapsed onto the consensus --
        ! independently of kernel bandwidth, which is why sweeping the bandwidth across a 4x range (runs
        ! 80/81/82) never moved on-molecule state correlation below 0.989.
        nk = max(1, min(ncomp, nkern))
        allocate(wcomp(nk), tvec(nk), tcen(nk,nstates), dist(nptcls), dvec(nk), mvec(nk))
        wcomp = 1.d0
        ! STANDARDIZED PLACEMENT (SIMPLE_COV_STDZ=1). Measured on this data, SIMPLE's component 1 carries
        ! 53.4% of the total latent variance against 20.7% (per-component sd ratio 13.0 vs 5.2) -- and z1 is
        ! NOT conformational here (0.20 with the dominant density mode, 0.25 with the swing). distinct
        ! ground-truth conformations 14 -> 20 of 20, Fab-swing coverage 55.9% -> 74.5%, amplitude 0.631 ->
        ! 0.713, pairwise state correlation 0.901 -> 0.845.
        if( .not. cov_env_flag_off('SIMPLE_COV_STDZ') )then
            do q = 1, nk
                zspread = sum(z(:,q)) / real(nptcls,dp)
                d2      = sum((z(:,q) - zspread)**2) / real(nptcls,dp)
                wcomp(q) = 1.d0 / max(d2, DTINY)
            end do
            wcomp = wcomp * real(nk,dp) / sum(wcomp)   ! keep the metric's overall scale
            write(logfhandle,'(A)') '>>> FLEX_PCA state placement metric: standardized (1/var per component)'
        endif
        ! degeneracy guard, in the same metric the targets are placed in:
        zspread = 0.d0
        do q = 1, nk
            zspread = zspread + wcomp(q)*(maxval(z(:,q)) - minval(z(:,q)))**2
        end do
        if( sqrt(zspread) <= sqrt(DTINY) ) &
            &THROW_HARD('flex_pca latent embedding has zero spread; embedding collapsed')
        ! The retained precision must be MARGINALISED, not sliced.
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
        else if( axis == 0 )then
            ! ---- k-means centroids over the retained latent subspace ----
            call kmeans_latent_targets(z(:,1:nk), nptcls, nk, nstates, wcomp, tcen)
            write(logfhandle,'(A,I0,A,I0)') '>>> FLEX_PCA state targets: k-means over ',nk, &
                &' components, k=',nstates
        else
            ! ---- equal-occupancy slices along one component ---- Two changes from the former
            ! quantile-SPACED placement.
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
            !$omp parallel do default(shared) private(i,q,r,d2,dvec,mvec) schedule(static)
            do i = 1, nptcls
                do q = 1, nk
                    dvec(q) = z(i,q) - tvec(q)
                end do
                ! Optional LOT metric:
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
            sorted = real(dist)
            call hpsort(sorted)
            ifloor = max(1, min(nptcls, min_neff))
            bmin   = max(real(sorted(ifloor),dp), chi2med)
            h      = sqrt(2.d0*bmin)          ! their kernel arg is sqrt(d^2/(2b)) => h^2 = 2b
            if( present(dist_out)   ) dist_out(:,state) = dist
            if( present(bfloor_out) ) bfloor_out(state) = bmin
            ! DIAGNOSTIC -- the chi-squared floor is only meaningful if the posterior quadratic form is
            ! actually on a chi2(ncomp) scale here.
            write(logfhandle,'(A,I3,A,ES11.3,A,ES11.3,A,ES11.3,A,A)') '>>>   state=',state, &
                &' dist: median=',real(sorted(max(1,nptcls/2)),dp),' p95=', &
                &real(sorted(max(1,nint(0.95*real(nptcls)))),dp),'  nn_floor=',real(sorted(ifloor),dp), &
                &'  bandwidth floor from ', merge('chi2(ncomp)', 'min_neff-nn', chi2med >= real(sorted(ifloor),dp))
            ! With h^2 = 2b in ncomp dimensions the enclosed population grows ~h^ncomp, so each 1.3x step
            ! multiplied it by ~1.3^20 ~ 190 at the ncomp=20 these runs actually use. states 7-10 each
            ! swallowed 69-73 % of the WHOLE dataset, 6 of 105 state pairs had kernel overlap > 0.95, only
            ! 3.11 kernels were effectively independent, and the delivered state maps correlated 0.9961 on
            ! the molecule -- i.e.
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
        ! Nearest-state label = argmax weight. Defaulting it to state 1 instead (the natural argmax
        ! initialisation) silently piles every unassigned particle onto the first state and makes the
        ! occupancy report below read as a concentration failure when it is nothing of the kind -- measured
        ! 78.4 % for state 1 whose actual effective support was 18.7 k of 100 k.
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
        ! NORMALISED WEIGHTS (SIMPLE_COV_NORMW=1) -- a partition of unity over the states. the last eight
        ! trajectory volumes plateau at 75-77 deg of a 53-106 deg range, and only 10 of 20 volumes are
        ! distinct conformations).
        if( cov_env_flag_on('SIMPLE_COV_NORMW') )then
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
        ! EXCLUSIVE ASSIGNMENT (SIMPLE_COV_EXCLUSIVE=1). Measured on the final run, nn_floor grows from 35.7
        ! at the central targets to 148.7 at the extremes: 8 of the last 9 trajectory volumes fell on a
        ! single ground-truth conformation, 80.0 deg, and the swing span stalled at 42 %).
        if( cov_env_flag_on('SIMPLE_COV_EXCLUSIVE') .and. present(dist_out) )then
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
        ! Occupancy report -- this is the number Defect 2 is about.
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
        ! Same change as in build_covariance_state_weights:
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

    !> even/odd correlation 0.49-0.99 at the floor bin rising to 0.98-0.997 at the next bin, which would
    !! have selected the widest bandwidth every time -- maximal smearing.
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
        ! --- bin grid per state, exactly pick_heterogeneity_bins2 ---
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
            call reconstruct_flex_diffmap_weighted_states(params, build, pinds, whalf, nstates, floor_rho=.true.)
            whalf = wbin
            call mask_state_weights_by_half(build, pinds, 1, whalf)
            params%outvol = 'flex_pca_cv'//bstr//'_odd_state_001.mrc'
            call reconstruct_flex_diffmap_weighted_states(params, build, pinds, whalf, nstates, floor_rho=.true.)
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

    !> Density-spread path targets: 37 % of particles fell outside every kernel support (94 % of the top |z|
    !! decile), and the median pairwise state correlation was 0.98 against 0.79.
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
        type(image)  :: vol_e, vol_o, vin, vout, vmsk
        type(string) :: fn, fn_out, spe, spo
        character(len=:), allocatable :: pe, po, base
        logical, allocatable :: l_mask(:,:,:)
        integer :: s, ldim(3), ipass
        real    :: mskrad_px
        ! consensus half maps:
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
        ! spherical support, as the rec_distr path does when no automask is available
        mskrad_px = 0.5 * params%mskdiam / params%smpd_crop
        call vmsk%disc(ldim, params%smpd_crop, mskrad_px, l_mask)
        call setup_nu_dmats(vol_e, vol_o, l_mask, [real ::])
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
        call vol_e%kill; call vol_o%kill; call vmsk%kill
        if( allocated(l_mask) ) deallocate(l_mask)
        write(logfhandle,'(A)') '>>> FLEX_PCA nonuniform-filtered maps written as *_nu.mrc &
            &(originals retained); local resolution map: flex_pca_nu_locres.mrc'
        call flush(logfhandle)
    end subroutine apply_consensus_nu_filter

    !> Rotate the covariance eigenbasis toward SPATIALLY COHERENT components. z4 0.72, z3 0.46, z9 0.35), so
    !! no single-component trajectory follows the rotation. VARIMAX IS THE WRONG ROTATION HERE and was
    !! measured to make it worse (best single 0.718 -> 0.693, concentration 0.350 -> 0.247).
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
        ! generalized symmetric eigenproblem by Cholesky reduction:
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

    !> Conformational trajectory as a DISCRETE GEODESIC over FINCH cluster representatives. raw z is
    !! dominated by the largest-variance component, which is measurably NOT the conformational one (z1
    !! correlates 0.371 with the true motion against z2's 0.931), and including all 20 components buries the
    !! manifold in noise dimensions where nearest neighbours are meaningless. the kd-tree first-neighbour
    !! search is the expensive part, so this clusters a strided subsample and assigns the rest by nearest
    !! representative -- sublinear in particle count, where Lloyd is linear with a 50x iteration constant.
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
        ! report the WHOLE hierarchy:
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

    !> Parse a comma-separated component list from an environment variable, e.g.
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

    !> Anchor direction for the trajectory: z4 alone reaches R2 0.512, +z3 0.680, +z9 0.776. The ground
    !! truth is equally fragmented (its own leading density PC explains R2 0.003 of the swing, three PCs
    !! 0.744), so this is a property of the data, not of the estimator -- and a single-axis trajectory is
    !! capped near 0.51 by construction while the 3-D span reaches 0.78.
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
        ! in a 6-component subspace the longest geodesic runs along whichever direction happens to be most
        ! extended, not along the motion, and the measured trajectory duly descended the true motion for 9
        ! states and then turned back (corr -0.788 against -0.979 for plain axis placement).
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
        ! Targets are the FULL latent coordinates, not the standardized subspace ones: The geodesic is
        ! typically only a handful of hops even on a fine level, so snapping assigns several states to the
        ! SAME particle and they reconstruct to identical volumes -- measured directly as blocks of repeated
        ! values in the trajectory projection (+0.07 +0.07 +0.07 +0.07 +0.47 ...) before this was fixed. the
        ! path folded back on itself for a quarter of its length, and the reconstructed swing duly went 92
        ! 89 81 75 97 76 71 deg instead of sweeping -- even though the anchor itself tracks the swing at
        ! corr -0.881.
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
        ! equal-OCCUPANCY slices:
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

    !> Deterministic k-means over the retained latent coordinates, in the SAME eigenvalue-weighted metric
    !! the Epanechnikov kernel uses, so target placement and particle weighting agree. latents carry
    !! near-equal variance (max/min 2.54), so there is no strong ordering to restore, and the eigenvalue
    !! rank is anti-correlated with conformational information on real data (-0.147), so the weighting
    !! placed every target along the two LEAST informative components.
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
                    ! empty cluster:
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
        ! ---- MOLECULE MASK ---- The axis score below is a half-set reproducibility, and it MUST be
        ! confined to the molecule. scored over the whole box the winner was the largest-variance component
        ! (pc1, corr to the ground-truth motion -0.28), because the solvent baseline and the global level of
        ! a state map reproduce trivially between halves and they dominate an unmasked correlation. Confined
        ! to the top 8% of the state mean the winner is pc2 (corr +0.85) -- the actual motion.
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
        ! AXIS SCORE = separation x self-consistency. an axis can order the states coherently while its two
        ! extreme states are nearly the same map, which is exactly what happened at k=15 (chosen pc2 scored
        ! 0.322 self-consistency and gave endpoints correlated 0.86, while a state pair at 0.79 existed in
        ! the set).
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
            ! ENDPOINT AMPLITUDE on the molecule, relative to the state mean. at ncols=256 five to seven
            ! components are reproducible and the criterion picked pc8 (reproducibility 0.859, but
            ! self-consistency 0.316 and endpoint amplitude 0.087 -- the trajectory collapsed).
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
        ! SIMPLE_COV_TRAJAXIS forces the trajectory axis. The automatic pick is the most REPRODUCIBLE
        ! component, which is not the same question as "which physical motion" -- on IgG-RL the most
        ! reproducible component (z2) is the dominant density mode while the Fab rotation lives in z4, and
        ! the two are near-orthogonal (GT corr 0.009).
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
        ! Diagnostic:
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

    ! ============================ SELF-CONTAINED TESTS ============================ No project, no images,
    ! no data files:

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

    !> Epanechnikov kernel contract at a fixed bandwidth:
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
