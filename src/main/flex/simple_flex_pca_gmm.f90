!@descr: flex_pca state placement: tied-covariance GMM and the hierarchical GMM AUTO weights
module simple_flex_pca_gmm
use simple_core_module_api
use simple_flex_pca_em, only: cov_env_int_pub
use simple_reconstructor, only: reconstructor
use simple_srch_sort_loc, only: hpsort
use simple_linalg, only: jacobi, eigsrt, matinv
use simple_flex_pca_util, only: two_gauss_unimodal
use simple_flex_pca_targets, only: kmeans_latent_targets
implicit none
private
#include "simple_local_flags.inc"

public :: gmm_state_weights, gmm_auto_state_weights


contains

    !>  Tied-covariance Gaussian-mixture responsibilities over the placed state targets. Replaces the
    !!  Epanechnikov kernel, whose compact support left many particles in no map at all; softmax over the
    !!  bandwidth is no substitute, it blurs every state back toward consensus. Tied covariance because
    !!  within-state spread is shared measurement error.
    !!  Measurements: doc/implementation_notes/flex_pca_state_placement_measurements.md
    subroutine gmm_state_weights( z, nptcls, ncomp, nk, nstates, tcen, wcomp, weights, neff, &
        &bandwidths, labels, pairsep, piout, maxit, respawn, pimin )
        integer,  intent(in)    :: nptcls, ncomp, nk, nstates
        real(dp), intent(in)    :: z(nptcls,ncomp), wcomp(nk)
        !> in: placed targets. out: FITTED means, so the reported table describes the delivered maps
        real(dp), intent(inout) :: tcen(nk,nstates)
        real,     intent(inout) :: weights(nptcls,nstates), neff(nstates), bandwidths(nstates)
        integer,  intent(inout) :: labels(nptcls)
        !> pairwise Mahalanobis separation of the FITTED means under the tied covariance;
        !! left at the -1 sentinel when the fit bails out on a singular covariance
        real(dp), optional, intent(out) :: pairsep(nstates,nstates)
        !> fitted mixing proportions
        real(dp), optional, intent(out) :: piout(nstates)
        !> EM iteration cap override; the tolerance still exits early on convergence
        integer,  optional, intent(in)  :: maxit
        !> disable the redundant-pair respawn (the discovery fit must: a respawned component
        !! lands on the worst-explained particle, an outlier in the gap BETWEEN clusters, and
        !! bridges the unimodality merge so everything chains into one macro-cluster)
        logical,  optional, intent(in)  :: respawn
        !> mixing-weight floor (constrained MLE); keeps components from starving below a
        !! deliverable occupancy. Inactive when the unconstrained fit already exceeds it.
        real(dp), optional, intent(in)  :: pimin
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
        real(dp) :: dev(nk)
        real(dp) :: ySy, ySm, lmax, lsum, ll, prev_ll, logdet, sumw, sumw2, trS
        real(dp) :: bicval, ent, dmin, d2pair
        real(dp) :: prev_trS
        integer  :: nfree, nrespawn, kmin, kdrop, iworst, maxit_eff
        integer  :: i, q, r, state, it, nrot, errflg
        integer(kind=8) :: nact_tot
        logical  :: l_respawn
        l_respawn = .true.
        if( present(respawn) ) l_respawn = respawn
        if( present(pairsep) ) pairsep = -1.d0
        if( present(piout)   ) piout   = 0.d0
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
        prev_trS  = -huge(1.d0)
        maxit_eff = GMM_MAXIT
        if( present(maxit) ) maxit_eff = maxit
        do it = 1, maxit_eff
            call matinv(S, Sinv, nk, errflg)
            if( errflg /= 0 )then
                write(logfhandle,'(A)') '>>> FLEX_PCA GMM tied covariance singular; keeping kernel weights'
                deallocate(y, mu, S, Sinv, Syy, resp, Smu, mSm, pival, nresp, evwork, ev, evec)
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
        !$omp parallel do default(shared) private(i,q,r,state,ySy,ySm,lmax,lsum) &
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
            do state = 1, nstates
                nresp(state) = sum(resp(:,state))
            end do
            pival = max(nresp, DTINY) / real(nptcls,dp)
            if( present(pimin) )then
                pival = max(pival, pimin)
                pival = pival / sum(pival)
            endif
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
            if( abs(ll - prev_ll) < GMM_TOL*abs(ll) )then
                ! Converged EM parks components on one region and starves others. Two means within GMM_MERGE_D2 are
                ! the same state; keep one, restart the other at the worst-explained particle.
                if( nrespawn < GMM_MAX_RESPAWN .and. l_respawn )then
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
        if( present(piout) ) piout = pival
        if( present(pairsep) )then
            call matinv(S, Sinv, nk, errflg)
            if( errflg == 0 )then
                pairsep = 0.d0
                do state = 1, nstates - 1
                    do r = state + 1, nstates
                        dev = mu(:,state) - mu(:,r)
                        d2pair = 0.d0
                        do q = 1, nk
                            d2pair = d2pair + dev(q)*sum(Sinv(q,:)*dev(:))
                        end do
                        pairsep(state,r) = sqrt(max(d2pair, 0.d0))
                        pairsep(r,state) = pairsep(state,r)
                    end do
                end do
            endif
        endif
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
    end subroutine gmm_state_weights

    !> Hierarchical mixture placement: over-fit a tied-covariance GMM, merge components whose
    !! pairwise density is unimodal (chains of continuum tiles connect, discrete islands stay
    !! separate), apportion the state budget over the macro-clusters by mass with a floor of one,
    !! and refit a GMM within each macro-cluster. Fixes the factor-mixing failure where a joint
    !! GMM cuts across composition x conformation: the island (e.g. a missing-domain minority)
    !! gets exactly one state and the continuum keeps the rest of the budget.
    subroutine gmm_auto_state_weights( z, nptcls, ncomp, nk, nstates, tcen, wcomp, min_neff, weights, &
        &neff, bandwidths, labels, macro_in )
        integer,  intent(in)    :: nptcls, ncomp, nk, nstates, min_neff
        !> macro-cluster per particle supplied by the latent deconvolution's mixture (full
        !! covariances, per-particle noise, held-out K): the discovery fit and the tied-covariance
        !! unimodality merge are skipped; too-small clusters still fold into their nearest, and
        !! clusters beyond the state budget fold smallest-first
        integer,  optional, intent(in) :: macro_in(:)
        real(dp), intent(in)    :: z(nptcls,ncomp), wcomp(nk)
        real(dp), intent(inout) :: tcen(nk,nstates)
        real,     intent(inout) :: weights(nptcls,nstates), neff(nstates), bandwidths(nstates)
        integer,  intent(inout) :: labels(nptcls)
        integer,  parameter   :: KFIT_MAX = 24
        !> minimum deliverable state occupancy: below this a map is noise, so no macro-cluster
        !! or seat allocation may create one. SIMPLE_COV_MIN_STATE overrides.
        integer,  parameter   :: GMM_MIN_OCC = 5000
        real(dp), allocatable :: tcen_d(:,:), sep(:,:), pifit(:), mass(:)
        real(dp), allocatable :: zsub(:,:), tcen_m(:,:), ysub(:,:), C(:,:), ev(:), evec(:,:)
        real,     allocatable :: w_d(:,:), neff_d(:), bw_d(:), w_m(:,:), neff_m(:), bw_m(:)
        real,     allocatable :: key(:), keym(:)
        integer,  allocatable :: lab_d(:), lab_m(:), macro(:), idx(:), budget(:), ord(:)
        integer,  allocatable :: ordm(:), cnt(:)
        real(dp) :: best, share, ymean(nk)
        real(dp), allocatable :: mn(:,:)
        integer  :: kfit, nmac, minisl, nseat, gstate, nm, bm, minocc
        integer  :: i, j, s, q, m, i0, i1, jbest, msmall, mbest, nrot
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
            goto 100
        endif
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
        allocate(key(nptcls), ord(nptcls))
        do i = 1, nptcls
            key(i) = real(z(i,jbest))
            ord(i) = i
        end do
        call hpsort(key, ord)
        allocate(tcen_d(nk,kfit), source=0.d0)
        do s = 1, kfit
            i0 = nint(real(s-1,dp)*real(nptcls,dp)/real(kfit,dp)) + 1
            i1 = max(nint(real(s,dp)*real(nptcls,dp)/real(kfit,dp)), i0)
            do j = i0, i1
                tcen_d(:,s) = tcen_d(:,s) + z(ord(j),1:nk)
            end do
            tcen_d(:,s) = tcen_d(:,s) / real(i1 - i0 + 1,dp)
        end do
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
100     continue
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
            if( bm == 1 )then
                do q = 1, nk
                    tcen_m(q,1) = sum(zsub(:,q))/real(nm,dp)
                end do
            else
                ! equal-mass quantile init along the macro's dominant internal axis
                allocate(ysub(nm,nk), C(nk,nk), ev(nk), evec(nk,nk), keym(nm), ordm(nm))
                do j = 1, nm
                    ysub(j,:) = zsub(j,1:nk)*sqrt(wcomp(:))
                end do
                do q = 1, nk
                    ymean(q) = sum(ysub(:,q))/real(nm,dp)
                end do
                do q = 1, nk
                    do s = 1, nk
                        C(q,s) = sum((ysub(:,q) - ymean(q))*(ysub(:,s) - ymean(s)))/real(nm,dp)
                    end do
                end do
                call jacobi(C, nk, nk, ev, evec, nrot)
                call eigsrt(ev, evec, nk, nk)
                do j = 1, nm
                    keym(j) = real(sum(ysub(j,:)*evec(:,1)))
                    ordm(j) = j
                end do
                call hpsort(keym, ordm)
                tcen_m = 0.d0
                do s = 1, bm
                    i0 = nint(real(s-1,dp)*real(nm,dp)/real(bm,dp)) + 1
                    i1 = max(nint(real(s,dp)*real(nm,dp)/real(bm,dp)), i0)
                    do j = i0, i1
                        tcen_m(:,s) = tcen_m(:,s) + zsub(ordm(j),1:nk)
                    end do
                    tcen_m(:,s) = tcen_m(:,s)/real(i1 - i0 + 1,dp)
                end do
                deallocate(ysub, C, ev, evec, keym, ordm)
            endif
            w_m = 0.; neff_m = 0.; bw_m = 0.; lab_m = 0
            if( bm == 1 )then
                call gmm_state_weights(zsub, nm, ncomp, nk, bm, tcen_m, wcomp, w_m, neff_m, bw_m, &
                    &lab_m, maxit=300, respawn=.false., &
                    &pimin=min(real(minocc,dp)/real(nm,dp), 0.5d0/real(bm,dp)))
            else
                ! A macro-cluster with several seats is a continuum (the deconvolution found it to be
                ! one Gaussian: a likelihood fit inside it has no modes to find and collapses, measured
                ! twice on 2026-09-08). Its sections are GEOMETRIC: k-means regions in the full
                ! standardised latent, so a high-variance nuisance axis (breathing) does not dictate
                ! the cut the way slices along the dominant axis did (v11: sections 2-4 identical maps).
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
