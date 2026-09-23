!@descr: flex_pca state weights: kernel/equal-mass/on-axis placement, bandwidth selection, half masks
module simple_flex_pca_weights
use simple_core_module_api
use simple_builder, only: builder
use simple_flex_pca_rec3D, only: reconstruct_flex_weighted_states, flex_rec_smpd
use simple_flex_pca_rounds, only: flex_pca_rounds
use simple_image, only: image
use simple_parameters, only: parameters
use simple_srch_sort_loc, only: hpsort
use simple_linalg, only: matinv
use simple_flex_pca_util, only: cov_env_flag_on, cov_env_flag_off, cov_env_dp, chi2_median, &
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
    subroutine build_covariance_state_weights( z, nptcls, ncomp, nkern, nstates, axis, min_neff, &
        &eigvals, precision, weights, targets, bandwidths, neff, labels, dist_out, bfloor_out, targets_in, &
        &zmetric, comp_rho, macro_in )
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
        logical  :: l_relpath, l_diffuse, l_gmm, l_gmm_auto
        character(len=12) :: bwsrc
        nk = max(1, min(ncomp, nkern))
        allocate(wcomp(nk), tvec(nk), tcen(nk,nstates), dist(nptcls), dvec(nk), mvec(nk))
        wcomp = 1.d0
        ! STANDARDIZED PLACEMENT. Eigenvalue weighting would concentrate every target along the
        ! highest-variance components, which are not the conformational ones.
        if( .true. )then
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
        else if( axis == 0 .and. present(comp_rho) )then
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
        else if( axis == 0 .and. present(comp_rho) )then
            ! 1-D equal-occupancy path. Correct on a genuine reaction coordinate, but it
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
        l_gmm = .true.
        ! EQUAL-MASS PLACEMENT IS NOT A GMM INITIALISATION. The tied-covariance mixture is a
        ! discrete-state model: it re-fits the means, and on a continuum with one dense mode every
        ! component slides into that mode -- which silently UNDOES the equal-occupancy placement that
        ! was just constructed (measured on the RNA data: sextiles in, one state holding 88% of the
        ! particles out). Where the targets carry equal mass by construction, keep them and let the
        ! along-path kernel deliver the frames. SIMPLE_COV_GMM=1 forces the refit back on for A/B.
        if( l_relpath )then
            write(logfhandle,'(A)') '>>> FLEX_PCA equal-mass targets: GMM refit SKIPPED &
                &(it would re-fit the means onto the dominant mode); along-path kernel weights kept'
            l_gmm = .false.
        endif
        if( l_gmm )then
            ! Hierarchical placement (default ON, SIMPLE_COV_GMM_AUTO=0 opts out): detect discrete
            ! islands vs continuum in the mixture itself and give each its own share of the budget.
            l_gmm_auto = .true.
            if( l_gmm_auto .and. nstates >= 3 )then
                call gmm_auto_state_weights(z, nptcls, ncomp, nk, nstates, tcen, wcomp, min_neff, weights, &
                    &neff, bandwidths, labels, macro_in=macro_in)
            else
                call gmm_state_weights(z, nptcls, ncomp, nk, nstates, tcen, wcomp, weights, neff, &
                    &bandwidths, labels)
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


    !> Cross-validated bandwidth selection, scored against the NARROWEST bin's opposite-half map.
    !! Plain even/odd agreement would rise monotonically with bandwidth and pick maximal smearing.
    subroutine cv_select_bandwidths( params, build, pinds, nptcls, nstates, nbins, min_neff, &
        &dist, bfloor, weights, bandwidths, neff , rounds)
        class(flex_pca_rounds), intent(inout) :: rounds
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
            call reconstruct_flex_weighted_states(params, build, pinds, whalf, nstates, floor_rho=.true., rounds=rounds)
            whalf = wbin
            call mask_state_weights_by_half(build, pinds, 1, whalf)
            params%outvol = 'flex_pca_cv'//bstr//'_odd_state_001.mrc'
            call reconstruct_flex_weighted_states(params, build, pinds, whalf, nstates, floor_rho=.true., rounds=rounds)
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
