!@descr: flex_pca state placement: tied-covariance GMM and the hierarchical GMM AUTO weights
module simple_flex_pca_gmm
use simple_core_module_api,  only: dp, dtiny, logfhandle
use simple_stat,             only: kish_ess
use simple_gmm,              only: gmm
use simple_clustering_utils, only: equal_mass_quantile_start
use simple_flex_pca_util,    only: two_gauss_unimodal
use simple_flex_pca_targets, only: kmeans_latent_targets
implicit none
private
#include "simple_local_flags.inc"

public :: gmm_state_weights, gmm_auto_state_weights


contains

    !>  Tied-covariance Gaussian-mixture responsibilities over the placed state targets (simple_gmm, in the
    !!  standardised frame); unlike the compact Epanechnikov kernel they leave no particle outside every map.
    !!  Tied covariance because within-state spread is shared measurement error. States come back ordered
    !!  by decreasing population.
    subroutine gmm_state_weights( z, nptcls, ncomp, nk, nstates, tcen, wcomp, weights, neff, &
        &bandwidths, labels, pairsep, piout, maxit, respawn, pimin )
        integer,            intent(in)    :: nptcls, ncomp, nk, nstates
        real(dp),           intent(in)    :: z(nptcls,ncomp), wcomp(nk)
        !> in: placed targets. out: FITTED means, so the reported table describes the delivered maps
        real(dp),           intent(inout) :: tcen(nk,nstates)
        real,               intent(inout) :: weights(nptcls,nstates), neff(nstates), bandwidths(nstates)
        integer,            intent(inout) :: labels(nptcls)
        !> pairwise Mahalanobis separation of the FITTED means under the tied covariance;
        !! left at the -1 sentinel when the fit bails out on a singular covariance
        real(dp), optional, intent(out)   :: pairsep(nstates,nstates)
        !> fitted mixing proportions
        real(dp), optional, intent(out)   :: piout(nstates)
        !> EM iteration cap override; the tolerance still exits early on convergence
        integer,  optional, intent(in)    :: maxit
        !> disable the redundant-pair respawn (the discovery fit must: a respawned component
        !! lands on the worst-explained particle, an outlier in the gap BETWEEN clusters, and
        !! bridges the unimodality merge so everything chains into one macro-cluster)
        logical,  optional, intent(in)    :: respawn
        !> mixing-weight floor (constrained MLE); keeps components from starving below a
        !! deliverable occupancy. Inactive when the unconstrained fit already exceeds it.
        real(dp), optional, intent(in)    :: pimin
        real(dp), parameter :: GMM_REG = 1.d-6, GMM_TOL = 1.d-5
        !> responsibilities below this are zeroed so the reconstructor's live-state compaction works
        real(dp), parameter :: RESP_FLOOR = 1.d-3
        integer,  parameter :: GMM_MAXIT  = 60
        type(gmm) :: gm
        real(dp), allocatable :: y(:,:), mu(:,:), resp(:,:), S(:,:), pival(:), nresp(:)
        real(dp) :: bicval, icl, trS
        integer  :: i, q, state, maxit_eff
        integer(kind=8) :: nact_tot
        logical  :: l_respawn, ok
        l_respawn = .true.
        if( present(respawn) ) l_respawn = respawn
        if( present(pairsep) ) pairsep = -1.d0
        if( present(piout)   ) piout   = 0.d0
        maxit_eff = GMM_MAXIT
        if( present(maxit) ) maxit_eff = maxit
        ! Standardised frame. wcomp is 1/var per component, so sqrt(wcomp) IS the 1/sd
        ! standardisation -- do not divide by sd as well, that standardises twice.
        allocate(y(nptcls,nk), mu(nk,nstates))
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
        if( present(pimin) )then
            call gm%new(y, nstates, mu, reg=GMM_REG, tol=GMM_TOL, maxits=maxit_eff, pi_floor=pimin, &
                &resp_floor=RESP_FLOOR, respawn=l_respawn)
        else
            call gm%new(y, nstates, mu, reg=GMM_REG, tol=GMM_TOL, maxits=maxit_eff, resp_floor=RESP_FLOOR, &
                &respawn=l_respawn)
        endif
        deallocate(y)
        call gm%fit(ok)
        if( .not. ok )then
            write(logfhandle,'(A)') '>>> FLEX_PCA GMM tied covariance singular; keeping kernel weights'
            call gm%kill
            deallocate(mu)
            return
        endif
        allocate(resp(nptcls,nstates), S(nk,nk), pival(nstates), nresp(nstates))
        call gm%get_means(mu)
        call gm%get_cov(S)
        call gm%get_pi(pival)
        call gm%get_mass(nresp)
        if( present(piout)   ) piout = pival
        if( present(pairsep) ) call gm%get_pairsep(pairsep, ok)
        write(logfhandle,'(A,I0,A,ES13.5,A,F7.4,A,F7.4)') &
            &'>>> FLEX_PCA GMM tied-covariance responsibilities: iters=',gm%get_niters(), &
            &'  loglik=',gm%get_loglik(),'  pi range ',minval(pival),' - ',maxval(pival)
        ! BIC will spend components on one populated state; ICL adds the entropy penalty and so prefers
        ! SEPARATED states. Reported per run; nothing is selected automatically yet.
        bicval = gm%get_bic()
        icl    = gm%get_icl()
        write(logfhandle,'(A,I0,A,I0,A,ES14.6,A,ES14.6,A,F8.4)') &
            &'>>> FLEX_PCA GMM model selection: K=',nstates,'  free params=', &
            &nstates*nk + (nk*(nk+1))/2 + nstates - 1,'  BIC=',bicval,'  ICL=',icl, &
            &'  mean entropy=',real((icl - bicval)/(2.d0*real(nptcls,dp)))
        ! responsibility mass near a component's own dimension cannot be estimated -- a merge candidate
        write(logfhandle,'(A,I0,A,I0)') '>>> FLEX_PCA GMM components below 3*nk effective mass: ', &
            &count(nresp < 3.d0*real(nk,dp)),' of ',nstates
        call flush(logfhandle)
        ! SPARSE responsibilities (floored and renormalized) -- REQUIRED, not an optimisation: the multi-state
        ! insertion assumes most scale pairs are EXACT zeros; left dense the state reconstruction does not finish.
        call gm%get_resp(resp)
        call gm%get_labels(labels)
        call gm%kill
        nact_tot = count(resp > 0.d0, kind=8)
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
            weights(:,state)  = real(resp(:,state))
            neff(state)       = real(kish_ess(resp(:,state)))
            bandwidths(state) = real(trS)
        end do
        do state = 1, nstates
            do q = 1, nk
                tcen(q,state) = mu(q,state) / sqrt(max(wcomp(q), DTINY))
            end do
        end do
        deallocate(mu, resp, S, pival, nresp)
    end subroutine gmm_state_weights

    !> Hierarchical mixture placement: over-fit a tied-covariance GMM, merge components whose
    !! pairwise density is unimodal (chains of continuum tiles connect, discrete islands stay
    !! separate), apportion the state budget over the macro-clusters by mass with a floor of one,
    !! and refit a GMM within each macro-cluster. Fixes the factor-mixing failure where a joint
    !! GMM cuts across composition x conformation: the island (e.g. a missing-domain minority)
    !! gets exactly one state and the continuum keeps the rest of the budget.
    subroutine gmm_auto_state_weights( z, nptcls, ncomp, nk, nstates, tcen, wcomp, min_neff, weights, &
        &neff, bandwidths, labels, macro_in )
        integer,           intent(in)    :: nptcls, ncomp, nk, nstates, min_neff
        !> macro-cluster per particle supplied by the latent deconvolution's mixture (full
        !! covariances, per-particle noise, held-out K): the discovery fit and the tied-covariance
        !! unimodality merge are skipped; too-small clusters still fold into their nearest, and
        !! clusters beyond the state budget fold smallest-first
        integer, optional, intent(in)    :: macro_in(:)
        real(dp),          intent(in)    :: z(nptcls,ncomp), wcomp(nk)
        real(dp),          intent(inout) :: tcen(nk,nstates)
        real,              intent(inout) :: weights(nptcls,nstates), neff(nstates), bandwidths(nstates)
        integer,           intent(inout) :: labels(nptcls)
        integer, parameter :: KFIT_MAX    = 24
        !> minimum deliverable state occupancy: below this a map is noise, so no macro-cluster
        !! or seat allocation may create one
        integer, parameter :: GMM_MIN_OCC = 5000
        real(dp), allocatable :: tcen_d(:,:), sep(:,:), pifit(:), mass(:)
        real(dp), allocatable :: zsub(:,:), tcen_m(:,:)
        real,     allocatable :: w_d(:,:), neff_d(:), bw_d(:), w_m(:,:), neff_m(:), bw_m(:)
        integer,  allocatable :: lab_d(:), lab_m(:), macro(:), idx(:), budget(:), cnt(:)
        real(dp) :: best, share
        real(dp), allocatable :: mn(:,:)
        integer  :: kfit, nmac, minisl, nseat, gstate, nm, bm, minocc
        integer  :: i, j, s, q, m, jbest, msmall, mbest
        logical  :: changed
        minocc = GMM_MIN_OCC
        kfit = min(KFIT_MAX, 2*nstates)
        kfit = min(kfit, max(2, nptcls/(3*max(nk, 1))))
        if( kfit <= nstates .or. nptcls < 10*nstates )then
            call gmm_state_weights(z, nptcls, ncomp, nk, nstates, tcen, wcomp, weights, neff, &
                &bandwidths, labels)
            return
        endif
        if( present(macro_in) )then
            ! ---- macro-clusters GIVEN: one per mixture component, means in the standardised metric
            ! for the nearest-neighbour folds below ----
            kfit = maxval(macro_in(1:nptcls))
            allocate(lab_d(nptcls), sep(kfit,kfit), pifit(kfit), mn(nk,kfit), macro(kfit))
            lab_d = macro_in(1:nptcls)
            mn = 0.d0; pifit = 0.d0
            do i = 1, nptcls
                mn(:,lab_d(i)) = mn(:,lab_d(i)) + z(i,1:nk)*sqrt(wcomp(:))
                pifit(lab_d(i)) = pifit(lab_d(i)) + 1.d0
            end do
            do s = 1, kfit
                if( pifit(s) > 0.d0 ) mn(:,s) = mn(:,s)/pifit(s)
                macro(s) = s
            end do
            pifit = pifit/real(nptcls,dp)
            do s = 1, kfit
                do j = 1, kfit
                    sep(s,j) = sqrt(sum((mn(:,s) - mn(:,j))**2))
                end do
            end do
            nmac = kfit
            write(logfhandle,'(A,I0,A)') '>>> FLEX_PCA GMM AUTO: macro-clusters from the deconvolved mixture (K=', &
                &kfit, '); discovery fit and tied-covariance merge skipped'
            call flush(logfhandle)
            deallocate(mn)
        else
            ! ---- DISCOVERY FIT at kfit components ----
            ! init quantiles along the axis with the largest OBSERVED standardised variance:
            ! component-1 init fails at some kfit when the delivered variance ordering disagrees
            ! with the eigenvalue ordering (measured on 10028: the amplitude anomaly is per-axis)
            best = -huge(0.d0)
            jbest = 1
            do q = 1, nk
                share = 0.d0
                do i = 1, nptcls
                    share = share + z(i,q)*z(i,q)
                end do
                share = share*wcomp(q)/real(nptcls,dp)
                if( share > best )then
                    best  = share
                    jbest = q
                endif
            end do
            write(logfhandle,'(A,I0)') '>>> FLEX_PCA GMM AUTO discovery init axis (max observed var): ', jbest
            allocate(tcen_d(nk,kfit))
            call equal_mass_quantile_start(z(:,1:nk), z(:,jbest), kfit, tcen_d)
            allocate(w_d(nptcls,kfit), neff_d(kfit), bw_d(kfit), lab_d(nptcls), &
                &sep(kfit,kfit), pifit(kfit))
            w_d = 0.; neff_d = 0.; bw_d = 0.; lab_d = 0
            ! measured on 10028: the island needs ~60 CLEAN iterations to separate; respawn
            ! resets eat a 60 cap, so give the discovery fit real convergence room
            call gmm_state_weights(z, nptcls, ncomp, nk, kfit, tcen_d, wcomp, w_d, neff_d, bw_d, &
                &lab_d, pairsep=sep, piout=pifit, maxit=300, respawn=.false.)
            if( any(sep < 0.d0) )then
                write(logfhandle,'(A)') '>>> FLEX_PCA GMM AUTO: discovery fit degenerate; plain GMM'
                call gmm_state_weights(z, nptcls, ncomp, nk, nstates, tcen, wcomp, weights, neff, &
                    &bandwidths, labels)
                return
            endif
            best = huge(0.d0)
            do s = 1, kfit - 1
                do j = s + 1, kfit
                    best = min(best, sep(s,j))
                end do
            end do
            write(logfhandle,'(A,F7.2,A,F7.2,A,F8.4)') '>>> FLEX_PCA GMM AUTO discovery separations: min=', &
                &real(best),'  max=',real(maxval(sep)),'  min mixing weight=',real(minval(pifit))
            ! ---- MERGE unimodal pairs into macro-clusters (transitive closure) ----
            ! Trivial-mass components carry no density evidence and can sit in the gap BETWEEN
            ! clusters, where transitive closure would chain through them; they get no vote in the
            ! graph and are attached to their nearest macro-cluster afterwards.
            allocate(macro(kfit), source=0)
            do s = 1, kfit
                if( pifit(s) < 0.25d0/real(kfit,dp) ) macro(s) = -1
            end do
            if( count(macro == 0) == 0 ) macro = 0        ! all trivial: no exclusion possible
            nmac = 0
            do s = 1, kfit
                if( macro(s) /= 0 ) cycle
                nmac     = nmac + 1
                macro(s) = nmac
                changed  = .true.
                do while( changed )
                    changed = .false.
                    do j = 1, kfit
                        if( macro(j) /= nmac ) cycle
                        do q = 1, kfit
                            if( macro(q) /= 0 ) cycle
                            if( two_gauss_unimodal(sep(j,q), pifit(j), pifit(q)) )then
                                macro(q) = nmac
                                changed  = .true.
                            endif
                        end do
                    end do
                end do
            end do
            ! attach excluded trivial components to the macro-cluster of their nearest peer
            do s = 1, kfit
                if( macro(s) /= -1 ) cycle
                best  = huge(0.d0)
                mbest = 0
                do j = 1, kfit
                    if( macro(j) <= 0 ) cycle
                    if( sep(s,j) < best )then
                        best  = sep(s,j)
                        mbest = macro(j)
                    endif
                end do
                macro(s) = max(mbest, 1)
            end do
        endif
        allocate(mass(nmac), cnt(nmac))
        mass = 0.d0
        cnt  = 0
        do s = 1, kfit
            mass(macro(s)) = mass(macro(s)) + pifit(s)
        end do
        do i = 1, nptcls
            cnt(macro(lab_d(i))) = cnt(macro(lab_d(i))) + 1
        end do
        ! macro-clusters too small to earn a map fold into their nearest neighbour. A mixture-defined
        ! cluster is already a population with its own covariance, so its floor is min_neff (a 3-4k
        ! population is a 7-8 A map), not the discovery fit's 5000.
        if( present(macro_in) )then
            minisl = max(min_neff, nptcls/200)
        else
            minisl = max(minocc, nptcls/200)
        endif
        do
            msmall = 0
            do m = 1, nmac
                if( cnt(m) < minisl .and. nmac > 1 )then
                    if( msmall == 0 )then
                        msmall = m
                    else if( cnt(m) < cnt(msmall) )then
                        msmall = m
                    endif
                endif
            end do
            if( msmall == 0 ) exit
            best  = huge(0.d0)
            mbest = 0
            do s = 1, kfit
                if( macro(s) /= msmall ) cycle
                do j = 1, kfit
                    if( macro(j) == msmall ) cycle
                    if( sep(s,j) < best )then
                        best  = sep(s,j)
                        mbest = macro(j)
                    endif
                end do
            end do
            if( mbest == 0 ) exit
            do s = 1, kfit
                if( macro(s) == msmall ) macro(s) = mbest
            end do
            do s = 1, kfit
                if( macro(s) > msmall ) macro(s) = macro(s) - 1
            end do
            nmac = nmac - 1
            mass(1:nmac) = 0.d0
            cnt(1:nmac)  = 0
            do s = 1, kfit
                mass(macro(s)) = mass(macro(s)) + pifit(s)
            end do
            do i = 1, nptcls
                cnt(macro(lab_d(i))) = cnt(macro(lab_d(i))) + 1
            end do
        end do
        ! given macro-clusters beyond the state budget: fold the smallest into its nearest until
        ! every remaining cluster can hold at least one seat
        if( present(macro_in) )then
            do while( nmac > nstates .and. nmac > 1 )
                msmall = minloc(cnt(1:nmac), dim=1)
                best  = huge(0.d0)
                mbest = 0
                do s = 1, kfit
                    if( macro(s) /= msmall ) cycle
                    do j = 1, kfit
                        if( macro(j) == msmall ) cycle
                        if( sep(s,j) < best )then
                            best  = sep(s,j)
                            mbest = macro(j)
                        endif
                    end do
                end do
                if( mbest == 0 ) exit
                do s = 1, kfit
                    if( macro(s) == msmall ) macro(s) = mbest
                end do
                do s = 1, kfit
                    if( macro(s) > msmall ) macro(s) = macro(s) - 1
                end do
                nmac = nmac - 1
                mass(1:nmac) = 0.d0
                cnt(1:nmac)  = 0
                do s = 1, kfit
                    mass(macro(s)) = mass(macro(s)) + pifit(s)
                end do
                do i = 1, nptcls
                    cnt(macro(lab_d(i))) = cnt(macro(lab_d(i))) + 1
                end do
            end do
        endif
        write(logfhandle,'(A,I0,A,I0,A)',advance='no') '>>> FLEX_PCA GMM AUTO: kfit=',kfit, &
            &'  macro-clusters=',nmac,'  (particles:'
        do m = 1, nmac
            write(logfhandle,'(A,I0)',advance='no') ' ',cnt(m)
        end do
        write(logfhandle,'(A)') ' )'
        call flush(logfhandle)
        ! mixture-defined clusters never fall back to the tied-covariance mixture: one cluster is
        ! one continuum that takes every seat as regions, and clusters were folded to the budget above
        if( .not. present(macro_in) .and. (nmac <= 1 .or. nmac >= nstates) )then
            if( nmac >= nstates ) write(logfhandle,'(A)') &
                &'>>> FLEX_PCA GMM AUTO: at least as many discrete clusters as states; plain GMM'
            call gmm_state_weights(z, nptcls, ncomp, nk, nstates, tcen, wcomp, weights, neff, &
                &bandwidths, labels)
            return
        endif
        ! ---- APPORTION the state budget: one state per macro-cluster, remaining seats to the
        ! largest unmet mass-per-seat, capacity-limited so a state always has particles behind it
        allocate(budget(nmac), source=1)
        nseat = nstates - nmac
        do while( nseat > 0 )
            jbest = 0
            best  = -1.d0
            do m = 1, nmac
                share = mass(m)/real(budget(m),dp)
                if( share > best .and. cnt(m) >= minocc*(budget(m) + 1) )then
                    best  = share
                    jbest = m
                endif
            end do
            if( jbest == 0 ) jbest = maxloc(cnt, dim=1)
            budget(jbest) = budget(jbest) + 1
            nseat = nseat - 1
        end do
        ! ---- REFIT a GMM within each macro-cluster; compose the global state arrays ----
        weights = 0.
        labels  = 0
        gstate  = 0
        do m = 1, nmac
            nm = cnt(m)
            bm = budget(m)
            allocate(idx(nm))
            j = 0
            do i = 1, nptcls
                if( macro(lab_d(i)) == m )then
                    j      = j + 1
                    idx(j) = i
                endif
            end do
            allocate(zsub(nm,ncomp), tcen_m(nk,bm), w_m(nm,bm), neff_m(bm), bw_m(bm), lab_m(nm))
            do j = 1, nm
                zsub(j,:) = z(idx(j),:)
            end do
            ! one seat: the macro mean starts the mixture; several seats: continuum_region_weights places
            ! its own k-means targets
            if( bm == 1 )then
                do q = 1, nk
                    tcen_m(q,1) = sum(zsub(:,q))/real(nm,dp)
                end do
            endif
            w_m = 0.; neff_m = 0.; bw_m = 0.; lab_m = 0
            if( bm == 1 )then
                call gmm_state_weights(zsub, nm, ncomp, nk, bm, tcen_m, wcomp, w_m, neff_m, bw_m, &
                    &lab_m, maxit=300, respawn=.false., &
                    &pimin=min(real(minocc,dp)/real(nm,dp), 0.5d0/real(bm,dp)))
            else
                ! A multi-seat macro-cluster is a continuum (one deconvolution Gaussian) where a likelihood fit
                ! collapses; its sections are k-means regions in the full standardised latent, so no axis dictates the cut.
                call continuum_region_weights(zsub(:,1:nk), nm, nk, bm, tcen_m, wcomp, min_neff, &
                    &w_m, neff_m, bw_m, lab_m)
            endif
            do s = 1, bm
                tcen(:,gstate+s)     = tcen_m(:,s)
                neff(gstate+s)       = neff_m(s)
                bandwidths(gstate+s) = bw_m(s)
                do j = 1, nm
                    weights(idx(j),gstate+s) = w_m(j,s)
                end do
            end do
            do j = 1, nm
                if( lab_m(j) >= 1 ) labels(idx(j)) = gstate + lab_m(j)
            end do
            write(logfhandle,'(A,I0,A,I0,A,I0,A,I0)') '>>> FLEX_PCA GMM AUTO macro=',m, &
                &'  particles=',nm,'  states=',bm,'  first state=',gstate + 1
            call flush(logfhandle)
            gstate = gstate + bm
            deallocate(idx, zsub, tcen_m, w_m, neff_m, bw_m, lab_m)
        end do
    end subroutine gmm_auto_state_weights

    !> State weights for one continuum macro-cluster: k-means targets on the denoised coordinates in
    !! the standardised metric (all nk dimensions, every axis at unit variance), and each particle
    !! assigned to its nearest target with unit weight. The regions partition the cluster; a compact
    !! kernel of half the target spacing was tried first and in 17 dimensions held only the core of
    !! each cell (5k of 51k particles in any support, the rest feeding no map).
    subroutine continuum_region_weights( zsub, nm, nk, bm, tcen_m, wcomp, min_neff, w_m, neff_m, bw_m, lab_m )
        integer,  intent(in)  :: nm, nk, bm, min_neff
        real(dp), intent(in)  :: zsub(nm,nk), wcomp(nk)
        real(dp), intent(out) :: tcen_m(nk,bm)
        real,     intent(out) :: w_m(nm,bm), neff_m(bm), bw_m(bm)
        integer,  intent(out) :: lab_m(nm)
        real(dp), allocatable :: dsum(:)
        real(dp) :: d2, dbest
        integer  :: j, s, q, sbest
        call kmeans_latent_targets(zsub, nm, nk, bm, wcomp, tcen_m)
        allocate(dsum(bm), source=0.d0)
        w_m = 0.
        do j = 1, nm
            sbest = 1
            dbest = huge(0.d0)
            do s = 1, bm
                d2 = 0.d0
                do q = 1, nk
                    d2 = d2 + wcomp(q)*(zsub(j,q) - tcen_m(q,s))**2
                end do
                if( d2 < dbest )then
                    dbest = d2
                    sbest = s
                endif
            end do
            lab_m(j)      = sbest
            w_m(j,sbest)  = 1.
            dsum(sbest)   = dsum(sbest) + sqrt(dbest)
        end do
        do s = 1, bm
            neff_m(s) = real(count(lab_m == s))
            bw_m(s)   = real(dsum(s)/real(max(1, count(lab_m == s)),dp))   ! mean distance to the target
        end do
        write(logfhandle,'(A,I0,A)') '>>> FLEX_PCA GMM AUTO continuum macro-cluster: ', bm, &
            &' k-means regions in the standardised latent (nearest target, unit weight)'
        do s = 1, bm
            write(logfhandle,'(A,I2,A,I0,A,ES10.3)') '>>>   region ', s, '  particles=', count(lab_m == s), &
                &'  mean distance to target=', bw_m(s)
        end do
        if( any(neff_m < real(min_neff)) ) write(logfhandle,'(A)') '>>>   WARNING: a region holds fewer than &
            &min_neff particles'
        call flush(logfhandle)
        deallocate(dsum)
    end subroutine continuum_region_weights

end module simple_flex_pca_gmm
