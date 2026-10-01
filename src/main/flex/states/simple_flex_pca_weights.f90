!@descr: flex_pca state weights: kernel/equal-mass/on-axis placement, bandwidth selection, half masks
module simple_flex_pca_weights
use simple_core_module_api
use simple_flex_pca_records, only: flex_latent, flex_selection, flex_state_set
use simple_builder, only: builder
use simple_flex_pca_rec3D, only: reconstruct_flex_weighted_states, flex_rec_smpd
use simple_flex_pca_rounds, only: flex_pca_rounds
use simple_image, only: image
use simple_parameters, only: parameters
use simple_srch_sort_loc, only: hpsort
use simple_linalg, only: matinv
use simple_flex_pca_run_types, only: flex_run_settings
use simple_flex_pca_util, only: chi2_median, &
    &kernel_weights_at_bandwidth, project_onto_target_polyline, COV_MAX_BW_GROW
use simple_flex_pca_gmm, only: gmm_state_weights, gmm_auto_state_weights
use simple_flex_pca_targets, only: diffusion_kcenter_targets, kmeans_latent_targets, path_latent_targets, reliability_path_targets
implicit none
private
#include "simple_local_flags.inc"

public :: build_covariance_state_weights, kernel_weights_at_bandwidth, cv_select_bandwidths, project_onto_target_polyline, mask_state_weights_by_half


contains

    !> Kernel-regression reconstruction weights (supplement S.F): place nstates latent targets, then give
    !! every particle an Epanechnikov weight per state from its Mahalanobis distance to that target. Each
    !! bandwidth is floored so at least min_neff particles fall inside its support.
    subroutine build_covariance_state_weights( latent, nkern, axis, min_neff, states, macro_in, equal_occ )
        type(flex_latent),    intent(in)    :: latent  !< z, precision and (when allocated) comp_rho
        type(flex_state_set), intent(inout) :: states  !< nstates in: the requested count; weights, targets, bandwidths, neff, labels, kdist, kfloor out
        integer :: nptcls, ncomp, nstates
        integer,  intent(in) :: nkern, axis, min_neff
        ! macro-cluster label per particle from the latent deconvolution's mixture; replaces GMM AUTO's discovery fit
        integer,  optional,    intent(in) :: macro_in(:)
        ! per-particle viewing AXIS; only needed for the GMM's orientation-coverage term
        real,     allocatable :: sorted(:)
        real(dp), allocatable :: wcomp(:), tvec(:), tcen(:,:), dist(:), dvec(:), mvec(:)
        real(dp), allocatable :: pk(:,:,:), cfull(:,:), cblk(:,:), edges(:)
        real(dp), allocatable :: ppath(:), tpath(:)   ! per-particle / per-target coordinate along the path
        integer,  allocatable :: occ(:)
        real(dp) :: h, d2, u2, sumw, sumw2, best, zspread, bmin, chi2med
        integer  :: ispace
        integer  :: i, q, r, state, best_state, grow, nfed, occmax, ifloor, nunassigned, nsupp
        integer  :: nk, errflg
        !> state_placement=equal_occ: at axis=0 the reliability-ordered equal-occupancy path replaces diffusion k-center
        logical,  optional, intent(in) :: equal_occ
        logical  :: l_relpath, l_diffuse, l_gmm, l_gmm_auto, l_equal_occ
        character(len=12) :: bwsrc
        nptcls  = size(latent%z,1)
        ncomp   = size(latent%z,2)
        nstates = states%nstates
        if( allocated(states%weights) )    deallocate(states%weights)
        if( allocated(states%targets) )    deallocate(states%targets)
        if( allocated(states%bandwidths) ) deallocate(states%bandwidths)
        if( allocated(states%neff) )       deallocate(states%neff)
        if( allocated(states%labels) )     deallocate(states%labels)
        if( allocated(states%kdist) )      deallocate(states%kdist)
        if( allocated(states%kfloor) )     deallocate(states%kfloor)
        nk = max(1, min(ncomp, nkern))
        l_equal_occ = .false.
        if( present(equal_occ) ) l_equal_occ = equal_occ
        allocate(wcomp(nk), tvec(nk), tcen(nk,nstates), dist(nptcls), dvec(nk), mvec(nk))
        wcomp = 1.d0
        ! STANDARDIZED PLACEMENT. Eigenvalue weighting would concentrate every target along the
        ! highest-variance components, which are not the conformational ones.
        do q = 1, nk
            zspread = sum(latent%z(:,q)) / real(nptcls,dp)
            d2      = sum((latent%z(:,q) - zspread)**2) / real(nptcls,dp)
            wcomp(q) = 1.d0 / max(d2, DTINY)
        end do
        wcomp = wcomp * real(nk,dp) / sum(wcomp)   ! keep the metric's overall scale
        write(logfhandle,'(A)') '>>> FLEX_PCA state placement metric: standardized (1/var per component)'
        ! degeneracy guard, in the same metric the targets are placed in
        zspread = 0.d0
        do q = 1, nk
            zspread = zspread + wcomp(q)*(maxval(latent%z(:,q)) - minval(latent%z(:,q)))**2
        end do
        if( sqrt(zspread) <= sqrt(DTINY) ) &
            &THROW_HARD('flex_pca latent embedding has zero spread; embedding collapsed')
        ! Restricting to nk components MARGINALISES the precision (invert, slice, re-invert); slicing the
        ! precision directly would condition on the dropped components instead.
        allocate(pk(nk,nk,nptcls))
        if( nk == ncomp )then
            pk = latent%precision
        else
            !$omp parallel do default(shared) private(i,cfull,cblk,errflg) schedule(static)
            do i = 1, nptcls
                allocate(cfull(ncomp,ncomp), cblk(nk,nk))
                call matinv(latent%precision(:,:,i), cfull, ncomp, errflg)
                cblk = cfull(1:nk,1:nk)
                call matinv(cblk, pk(:,:,i), nk, errflg)
                deallocate(cfull, cblk)
            end do
            !$omp end parallel do
            write(logfhandle,'(A,I0,A,I0,A)') '>>> FLEX_PCA state stage restricted to the leading ',nk, &
                &' of ',ncomp,' latent components (marginalised precision)'
        endif
        ! Set inside the placing branch: along-path targets (external curve, reliability path) also select
        ! the along-path weighting.
        l_relpath = .false.
        allocate(states%weights(nptcls,nstates), states%targets(ncomp,nstates), states%bandwidths(nstates), states%neff(nstates), states%labels(nptcls))
        if( .true.   ) allocate(states%kdist(nptcls,nstates))
        allocate(states%kfloor(nstates))
        if( axis < 0 )then
            call path_latent_targets(latent%z(:,1:nk), nptcls, nk, nstates, wcomp, tcen)
            write(logfhandle,'(A,I0,A,I0)') '>>> FLEX_PCA state targets: density-spread path over ',nk, &
                &' components, points=',nstates
        else if( axis == 0 .and. allocated(latent%comp_rho) )then
            ! DEFAULT: diffusion k-center. Handles a continuous reaction coordinate and branched compositional
            ! states with the same constants, because a curve is a degenerate graph.
            ! state_placement=equal_occ skips it and takes the reliability-path fallback (equal-count slices)
            l_diffuse = .false.
            if( .not. l_equal_occ ) call diffusion_kcenter_targets(latent%z(:,1:nk), nptcls, nk, nstates, wcomp, &
                &latent%comp_rho(1:nk), tcen, l_diffuse)
            if( l_diffuse )then
                write(logfhandle,'(A,I0,A,I0)') '>>> FLEX_PCA state targets: diffusion k-center over ', &
                    &nk,' components, points=',nstates
            else
                if( l_equal_occ )then
                    write(logfhandle,'(A)') '>>> FLEX_PCA state_placement=equal_occ: reliability-ordered &
                        &equal-occupancy path'
                else
                    write(logfhandle,'(A)') '>>> FLEX_PCA diffusion k-center unavailable; falling back &
                        &to the reliability-ordered path'
                endif
                allocate(ppath(nptcls), tpath(nstates))
                call reliability_path_targets(latent%z(:,1:nk), nptcls, nk, nstates, wcomp, latent%comp_rho(1:nk), tcen, &
                    &proj_out=ppath, tproj_out=tpath)
                l_relpath = .true.
            endif
        else if( axis == 0 )then
            call kmeans_latent_targets(latent%z(:,1:nk), nptcls, nk, nstates, wcomp, tcen)
            write(logfhandle,'(A,I0,A,I0)') '>>> FLEX_PCA state targets: k-means over ',nk, &
                &' components, k=',nstates
        else
            ! equal-occupancy slices along one component; each target is the slice MEAN over all nk components
            if( axis > nk ) THROW_HARD('flex_pca state_axis exceeds the retained component count nkern')
            allocate(sorted(nptcls), source=real(latent%z(:,axis)))
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
                    if( latent%z(i,axis) < edges(state) )then
                        best_state = state
                        exit
                    endif
                end do
                occ(best_state) = occ(best_state) + 1
                do q = 1, nk
                    tcen(q,best_state) = tcen(q,best_state) + latent%z(i,q)
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
                ppath(i) = latent%z(i,axis)
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
                    dvec(q) = latent%z(i,q) - tvec(q)
                end do
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
            if( .true.   ) states%kdist(:,state) = dist
            states%kfloor(state) = bmin
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
                    states%weights(i,state) = real(max(0.d0, 1.d0 - u2))   ! Epanechnikov, compact support
                    sumw  = sumw  + real(states%weights(i,state),dp)
                    sumw2 = sumw2 + real(states%weights(i,state),dp)**2
                    if( states%weights(i,state) > 0. ) nsupp = nsupp + 1
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
            if( maxval(states%weights(:,state)) > TINY ) states%weights(:,state)=states%weights(:,state)/maxval(states%weights(:,state))
            sumw  = sum(real(states%weights(:,state),dp))
            sumw2 = sum(real(states%weights(:,state),dp)**2)
            ! targets is reported over the FULL component set for the manifest
            states%targets(1:nk,state) = real(tcen(:,state))
            do q = nk+1, ncomp
                states%targets(q,state) = real(sum(latent%z(:,q)) / real(nptcls,dp))
            end do
            states%bandwidths(state) = real(h)
            states%neff(state)       = real(sumw*sumw/max(sumw2,DTINY))
        end do
        ! Tied-covariance mixture by default; SIMPLE_COV_GMM=0 recovers the kernel. The kernel loop above
        ! still runs: dist_out feeds cv_select_bandwidths and its quantiles diagnose the chi2 scale.
        l_gmm = .true.
        ! EQUAL-MASS PLACEMENT IS NOT A GMM INITIALISATION. The tied-covariance mixture is a
        ! discrete-state model: it re-fits the means, and on a continuum with one dense mode every
        ! component slides into that mode -- which silently UNDOES the equal-occupancy placement that
        ! was just constructed (measured on the RNA data: sextiles in, one state holding 88% of the
        ! particles out). Where the targets carry equal mass by construction, keep them and let the
        ! along-path kernel deliver the frames. SIMPLE_COV_GMM=1 forces the refit back on for A/B.
        if( l_relpath )then
            write(logfhandle,'(A)') '>>> FLEX_PCA equal-mass states%targets: GMM refit SKIPPED &
                &(it would re-fit the means onto the dominant mode); along-path kernel states%weights kept'
            l_gmm = .false.
        endif
        if( l_gmm )then
            ! Hierarchical placement (default ON, SIMPLE_COV_GMM_AUTO=0 opts out): detect discrete
            ! islands vs continuum in the mixture itself and give each its own share of the budget.
            l_gmm_auto = .true.
            if( l_gmm_auto .and. nstates >= 3 )then
                call gmm_auto_state_weights(latent%z, nptcls, ncomp, nk, nstates, tcen, wcomp, min_neff, states%weights, &
                    &states%neff, states%bandwidths, states%labels, macro_in=macro_in)
            else
                call gmm_state_weights(latent%z, nptcls, ncomp, nk, nstates, tcen, wcomp, states%weights, states%neff, &
                    &states%bandwidths, states%labels)
            endif
            ! tcen now holds the FITTED means; refresh the reported targets to describe the delivered maps
            do state = 1, nstates
                states%targets(1:nk,state) = real(tcen(:,state))
            end do
        endif
        ! argmax weight, 0 outside EVERY kernel support. Defaulting to state 1 would pile the unassigned
        ! onto the first state and fake a concentration failure.
        nunassigned = 0
        do i = 1, nptcls
            best_state = 0
            best = 0.d0
            do state = 1, nstates
                if( real(states%weights(i,state),dp) > best )then
                    best = real(states%weights(i,state),dp)
                    best_state = state
                endif
            end do
            states%labels(i) = best_state
            if( best_state == 0 ) nunassigned = nunassigned + 1
        end do
        allocate(occ(nstates), source=0)
        do i = 1, nptcls
            if( states%labels(i) >= 1 ) occ(states%labels(i)) = occ(states%labels(i)) + 1
        end do
        occmax = maxval(occ)
        nfed   = count(occ > 100)
        write(logfhandle,'(A,F6.2,A,I0,A,I0)') '>>> FLEX_PCA state occupancy (of assigned): max=', &
            &100.0*real(occmax)/real(max(nptcls-nunassigned,1)),'%  states with >100 particles=',nfed,' of ',nstates
        write(logfhandle,'(A,I0,A,F6.2,A)') '>>> FLEX_PCA outside every kernel support: ',nunassigned, &
            &' particles (',100.0*real(nunassigned)/real(max(nptcls,1)),'%) -- these contribute to no state map'
        do state = 1, nstates
            write(logfhandle,'(A,I3,A,I9,A,F10.1,A,ES11.3,A,F7.3)') '>>>   state=',state,'  particles=',occ(state), &
                &'  neff=',states%neff(state),'  bandwidth=',states%bandwidths(state),'  frac_in_support=', &
                &real(count(states%weights(:,state) > 0.))/real(max(nptcls,1))
        end do
        call flush(logfhandle)
        deallocate(wcomp, tvec, tcen, occ, dist, dvec, sorted, pk)
    end subroutine build_covariance_state_weights


    !> Cross-validated bandwidth selection, scored against the NARROWEST bin's opposite-half map.
    !! Plain even/odd agreement would rise monotonically with bandwidth and pick maximal smearing.
    subroutine cv_select_bandwidths( params, cfg, build, sel, nbins, min_neff, states , rounds)
        type(flex_selection), intent(in)    :: sel
        type(flex_state_set), intent(inout) :: states  !< kdist/kfloor in; weights, bandwidths, neff re-selected
        class(flex_pca_rounds), intent(inout) :: rounds
        class(parameters), intent(inout) :: params
        type(flex_run_settings), intent(in)    :: cfg
        class(builder),    intent(inout) :: build
        integer,           intent(in) :: nbins, min_neff
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
        allocate(wbin(sel%nptcls,states%nstates), whalf(sel%nptcls,states%nstates), bins(nbins,states%nstates), &
            &hbin(nbins,states%nstates), err(nbins,states%nstates), sorted(sel%nptcls))
        err = 0.d0
        ! --- bin grid per state: from the bandwidth floor to the 95th distance percentile, linear in sqrt
        p95 = max(1, min(sel%nptcls, nint(0.95*real(sel%nptcls))))
        do state = 1, states%nstates
            sorted = real(states%kdist(:,state))
            call hpsort(sorted)
            b_lo = states%kfloor(state)
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
        allocate(tgt_ev(nvox,states%nstates), tgt_od(nvox,states%nstates), source=0.)
        do ib = 1, nbins
            do state = 1, states%nstates
                call kernel_weights_at_bandwidth(states%kdist(:,state), sel%nptcls, sqrt(2.d0*bins(ib,state)), &
                    &min_neff, wbin(:,state), h_used, neff_used)
                hbin(ib,state) = h_used
            end do
            write(bstr,'(I3.3)') ib
            whalf = wbin
            call mask_state_weights_by_half(build, sel%pinds, 0, whalf)
            params%outvol = 'flex_pca_cv'//bstr//'_even_state_001.mrc'
            call reconstruct_flex_weighted_states(params, cfg, build, sel%pinds, whalf, states%nstates, floor_rho=.true., rounds=rounds)
            whalf = wbin
            call mask_state_weights_by_half(build, sel%pinds, 1, whalf)
            params%outvol = 'flex_pca_cv'//bstr//'_odd_state_001.mrc'
            call reconstruct_flex_weighted_states(params, cfg, build, sel%pinds, whalf, states%nstates, floor_rho=.true., rounds=rounds)
            do state = 1, states%nstates
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
        do state = 1, states%nstates
            ibest = minloc(err(:,state), dim=1)
            call kernel_weights_at_bandwidth(states%kdist(:,state), sel%nptcls, sqrt(2.d0*bins(ibest,state)), &
                &min_neff, states%weights(:,state), h_used, neff_used)
            states%bandwidths(state) = real(h_used)
            states%neff(state)       = neff_used
            write(logfhandle,'(A,I3,A,I3,A,I0,A,ES11.3,A,ES12.4,A,F10.1)') '>>>   state=',state,'  bin=',ibest,'/',nbins, &
                &'  h=',h_used,'  error=',err(ibest,state),'  neff=',neff_used
        end do
        call flush(logfhandle)
        deallocate(wbin, whalf, bins, hbin, err, sorted, tgt_ev, tgt_od)
    end subroutine cv_select_bandwidths

    !> Arc-length coordinate of every particle on the polyline through the supplied targets.
    !! Projection is Euclidean in the raw latent (the units the targets are given in); each particle
    !! takes the closest point over all segments, clamped to the segment ends, so particles beyond
    !! either terminus map to the terminus and the end frames collect the tails.

    subroutine mask_state_weights_by_half( build, pinds, wanted_eo, weights )
        type(builder), intent(inout) :: build
        integer,       intent(in)    :: pinds(:), wanted_eo
        real,          intent(inout) :: weights(:,:)
        integer :: i
        do i = 1, size(pinds)
            if( build%spproj_field%get_eo(pinds(i)) /= wanted_eo ) weights(i,:) = 0.
        end do
    end subroutine mask_state_weights_by_half

end module simple_flex_pca_weights
