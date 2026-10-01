!@descr: flex_pca posterior inference: the per-particle latent posterior from sufficient statistics
!!
!! Consumes what an E-step formulation produces for one particle -- the projected Gram G, the
!! data and mean rows b and c, the mean/data scalars -- and returns the posterior: the MAP
!! latent (or the mixture-weighted mean), its precision/covariance, log det, likelihood
!! contribution and the fitted contrast (ECM/MCFA), plus the mixture M-step and its
!! initialisation and the estimator precision of the MAP embedding. Dense algebra on arrays;
!! no projection code, no fit state, no file.
module simple_flex_pca_posterior
use simple_core_module_api
use simple_linalg, only: jacobi, eigsrt
use simple_flex_pca_fit_types, only: flex_fit
implicit none
private

public :: spd_logdet_dp, quad_form, spd_solve_dp, spd_inv_dp
public :: probe_solve_ecm, probe_solve_mix
public :: mcfa_init, mcfa_condition, mcfa_mstep
public :: map_sampling_precision

real(dp), parameter :: COV_PINV_RCOND = 1.0d-6

contains

    !> log(det A) for symmetric positive-definite A, by Cholesky on a private copy.
    !!
    !! This is the term the probe's `resid_energy` has always been missing. resid_energy is the JOINT
    !! MAP objective ||.||^2/sig2 + z'Gamma^-1 z evaluated at zhat, which decreases by construction and
    !! therefore says nothing about convergence. The MARGINAL likelihood needs log det A_i and
    !! log det Gamma as well, and without them the EM has no objective to watch -- which is why the
    !! iteration count was a tuned constant standing in for a stopping rule.
    pure module subroutine spd_logdet_dp( A, n, logdet, ok )
        integer,  intent(in)  :: n
        real(dp), intent(in)  :: A(n,n)
        real(dp), intent(out) :: logdet
        logical,  intent(out) :: ok
        real(dp) :: L(n,n), s
        integer  :: i, j
        L      = 0.d0
        logdet = 0.d0
        ok     = .false.
        do j = 1, n
            s = A(j,j) - sum(L(j,1:j-1)**2)
            if( s <= 0.d0 ) return
            L(j,j) = sqrt(s)
            logdet = logdet + 2.d0*log(L(j,j))
            do i = j+1, n
                L(i,j) = (A(i,j) - sum(L(i,1:j-1)*L(j,1:j-1))) / L(j,j)
            end do
        end do
        ok = .true.
    end subroutine spd_logdet_dp


    !>  z' M z for symmetric M.
    pure module function quad_form( M, z, n ) result( val )
        integer,  intent(in) :: n
        real(dp), intent(in) :: M(n,n), z(n)
        real(dp) :: val
        integer  :: q, r
        val = 0.d0
        do r = 1, n
            do q = 1, n
                val = val + z(q)*M(q,r)*z(r)
            end do
        end do
    end function quad_form


    !> In-place symmetric positive-definite solve A x = b (b overwritten by x) via Cholesky. A is first
    !! scaled by its mean diagonal, so the retry ridge is RELATIVE: an absolute ridge either swamps a
    !! small-diagonal system or fails to rescue a large one, and the b=0 fallback then collapses the
    !! latents of essentially every particle.
    subroutine spd_solve_dp( A, b, n )
        integer,  intent(in)    :: n
        real(dp), intent(inout) :: A(n,n), b(n)
        real(dp) :: L(n,n), s, y(n), ridge, dscale
        integer  :: i, j, attempt
        dscale = 0.d0
        do i = 1, n
            dscale = dscale + abs(A(i,i))
        end do
        dscale = dscale / real(n,dp)
        if( dscale > 0.d0 )then
            A = A / dscale
            b = b / dscale
        endif
        do attempt = 1, 3
            L = 0.d0
            do j = 1, n
                s = A(j,j) - sum(L(j,1:j-1)**2)
                if( s <= 0.d0 ) exit
                L(j,j) = sqrt(s)
                do i = j+1, n
                    L(i,j) = (A(i,j) - sum(L(i,1:j-1)*L(j,1:j-1))) / L(j,j)
                end do
            end do
            if( j > n )then
                ! forward/back substitution
                do i = 1, n
                    y(i) = (b(i) - sum(L(i,1:i-1)*y(1:i-1))) / L(i,i)
                end do
                do i = n, 1, -1
                    b(i) = (y(i) - sum(L(i+1:n,i)*b(i+1:n))) / L(i,i)
                end do
                return
            endif
            ridge = 1.d-8 * (abs(A(1,1))+1.d0) * (10.d0**(attempt-1))
            do i = 1, n
                A(i,i) = A(i,i) + ridge
            end do
        end do
        b = 0.d0
    end subroutine spd_solve_dp


    subroutine spd_inv_dp( A, Ainv, n )
        integer,  intent(in)    :: n
        real(dp), intent(inout) :: A(n,n)
        real(dp), intent(out)   :: Ainv(n,n)
        real(dp) :: L(n,n), Linv(n,n), s, ridge, dscale
        integer  :: i, j, attempt
        Ainv   = 0.d0
        dscale = 0.d0
        do i = 1, n
            dscale = dscale + abs(A(i,i))
        end do
        dscale = dscale / real(n,dp)
        if( dscale > 0.d0 ) A = A / dscale
        do attempt = 1, 3
            L = 0.d0
            do j = 1, n
                s = A(j,j) - sum(L(j,1:j-1)**2)
                if( s <= 0.d0 ) exit
                L(j,j) = sqrt(s)
                do i = j+1, n
                    L(i,j) = (A(i,j) - sum(L(i,1:j-1)*L(j,1:j-1))) / L(j,j)
                end do
            end do
            if( j > n )then
                ! Linv = L^-1 by forward substitution on the identity, column by column
                Linv = 0.d0
                do j = 1, n
                    Linv(j,j) = 1.d0 / L(j,j)
                    do i = j+1, n
                        Linv(i,j) = -sum(L(i,j:i-1)*Linv(j:i-1,j)) / L(i,i)
                    end do
                end do
                ! A = L L' => A^-1 = (L^-1)' (L^-1), lower-triangular so the sum starts at max(i,j)
                do i = 1, n
                    do j = 1, i
                        Ainv(i,j) = sum(Linv(i:n,i)*Linv(i:n,j))
                        Ainv(j,i) = Ainv(i,j)
                    end do
                end do
                ! undo the rescaling: A_orig = dscale*A_scaled, so A_orig^-1 = A_scaled^-1 / dscale
                if( dscale > 0.d0 ) Ainv = Ainv / dscale
                return
            endif
            ridge = 1.d-8 * (abs(A(1,1))+1.d0) * (10.d0**(attempt-1))
            do i = 1, n
                A(i,i) = A(i,i) + ridge
            end do
        end do
    end subroutine spd_inv_dp


    !> One particle's MAP solve; nml>0 adds ECM contrast updates against the current basis,
    !!   a <- (m'y + b'z) / (||m||^2 + 2c'z + z'Gz + tr(G A^-1)),  clamped to [0.1, 5].
    !! tr(G A^-1) is the posterior variance; dropping it biases a high.
    !> One particle's ECM posterior solve under the plain 1/Gamma prior: the fit's E-step
    !! statistics of thread ithr in, the latent, its posterior covariance and the contrast out.
    subroutine probe_solve_ecm( fit, ithr, myv, e_mm, a, ldA, lok, quad )
        class(flex_fit), intent(inout) :: fit
        integer,         intent(in)    :: ithr
        real(dp),        intent(in)    :: myv, e_mm
        real(dp),        intent(inout) :: a
        real(dp),        intent(out)   :: ldA
        logical,         intent(out)   :: lok
        real(dp),        intent(out)   :: quad
        real(dp) :: Amat(fit%model%ncomp,fit%model%ncomp), Acp(fit%model%ncomp,fit%model%ncomp), h(fit%model%ncomp), aa, a_num, a_den
        integer  :: q, icm
        associate( n => fit%model%ncomp, nml => fit%spec%nml_plain, sig2 => fit%model%sig2, G => fit%iter%Gth(:,:,ithr), b => fit%iter%bth(:,ithr), &
            &c => fit%iter%cth(:,ithr), prior_ => fit%iter%prior, z_ => fit%iter%zth(:,ithr), Ainv_ => fit%iter%Ainvth(:,:,ithr) )
        do icm = 0, max(0, nml)
            aa   = a*a
            Amat = (aa/sig2)*G
            do q = 1, n
                Amat(q,q) = Amat(q,q) + prior_(q)
                z_(q)     = (a*b(q) - aa*c(q))/sig2
            end do
            Acp = Amat
            call spd_logdet_dp(Amat, n, ldA, lok)
            h = z_
            call spd_inv_dp(Acp, Ainv_, n)
            call spd_solve_dp(Amat, z_, n)
            quad = dot_product(h, z_)
            if( icm >= max(0, nml) ) exit
            a_num = myv + dot_product(b, z_)
            a_den = e_mm + 2.d0*dot_product(c, z_) + quad_form(G, z_, n) + sum(G*Ainv_)
            if( a_den > DTINY ) a = min(5.0d0, max(0.1d0, a_num/a_den))
        end do
        end associate
    end subroutine probe_solve_ecm


    !> Per-particle MCFA posterior from fetched sufficient statistics (G, b, c) -- the ONE
    !! source for the mixture E-step, shared by the plain CPU body and the fused device
    !! body (the device computes the same G/b/c on card; only this host solve differs from
    !! the single-Gaussian path). Updates the caller's thread-slice accumulators exactly as
    !! the historical in-line block did.
    !> One particle's posterior solve under the MCFA mixture prior: the fit's E-step statistics of
    !! thread ithr in; the latent mean, the second moment of batch row i, the thread's mixture
    !! accumulators and the likelihood term out.
    subroutine probe_solve_mix( fit, ithr, i, myv, e_mm, a, ldA, lok, nll_add )
        class(flex_fit), intent(inout) :: fit
        integer,         intent(in)    :: ithr, i
        real(dp),        intent(in)    :: myv, e_mm
        real(dp),        intent(inout) :: a
        real(dp),        intent(out)   :: ldA, nll_add
        logical,         intent(out)   :: lok
        real(dp) :: Amat(fit%model%ncomp,fit%model%ncomp), Acp(fit%model%ncomp,fit%model%ncomp), rhs0(fit%model%ncomp), mk(fit%model%ncomp,fit%spec%kmix), lw(fit%spec%kmix), rk(fit%spec%kmix)
        real(dp) :: aa, lwm, wsm, a_num, a_den
        integer  :: q, r, kk, icm
        associate( n => fit%model%ncomp, kmix => fit%spec%kmix, nml => fit%spec%nml_plain, sig2 => fit%model%sig2, &
            &G => fit%iter%Gth(:,:,ithr), b => fit%iter%bth(:,ithr), c => fit%iter%cth(:,ithr), &
            &Ominv => fit%history%mix_Ominv, Omxi => fit%history%mix_Omxi, lpi => fit%history%mix_lpi, xiOx => fit%history%mix_xiOx, &
            &zbar => fit%iter%zth(:,ithr), Ainv_ => fit%iter%Ainvth(:,:,ithr), dens_ => fit%iter%dens(:,:,i), &
            &sr_acc => fit%history%mxa_sr(:,ithr), sm_acc => fit%history%mxa_sm(:,:,ithr), smm_acc => fit%history%mxa_smm(:,:,:,ithr), sainv_acc => fit%history%mxa_sainv(:,:,ithr) )
        do icm = 0, max(0, nml)
            aa   = a*a
            Amat = (aa/sig2)*G + Ominv
            Acp  = Amat
            call spd_logdet_dp(Amat, n, ldA, lok)
            call spd_inv_dp(Acp, Ainv_, n)
            do q = 1, n
                rhs0(q) = (a*b(q) - aa*c(q))/sig2
            end do
            do kk = 1, kmix
                mk(:,kk) = matmul(Ainv_, rhs0 + Omxi(:,kk))
                lw(kk)   = lpi(kk) - 0.5d0*xiOx(kk) + 0.5d0*dot_product(rhs0 + Omxi(:,kk), mk(:,kk))
            end do
            lwm = maxval(lw)
            wsm = 0.d0
            do kk = 1, kmix
                rk(kk) = exp(max(-7.d2, lw(kk) - lwm))
                wsm    = wsm + rk(kk)
            end do
            rk = rk / max(wsm, DTINY)
            zbar = 0.d0
            do kk = 1, kmix
                zbar = zbar + rk(kk)*mk(:,kk)
            end do
            ! E[zz'|y] = A^-1 + sum_k r_k m_k m_k' (between-component spread is real variance)
            dens_ = Ainv_
            do kk = 1, kmix
                do r = 1, n
                    do q = 1, n
                        dens_(q,r) = dens_(q,r) + rk(kk)*mk(q,kk)*mk(r,kk)
                    end do
                end do
            end do
            if( icm >= max(0, nml) ) exit
            ! ECM amplitude update under the mixture: the plain path's a-update with the
            ! mixture posterior moments in place of the single-Gaussian ones. This is the
            ! lever that separates on-particle amplitude (ice) from conformation.
            a_num = myv + dot_product(b, zbar)
            a_den = e_mm + 2.d0*dot_product(c, zbar) + sum(G*dens_)
            if( a_den > DTINY ) a = min(5.0d0, max(0.1d0, a_num/a_den))
        end do
        ! mixture marginal, varying part; the N*logdet(Omega) term is added once globally
        nll_add = ldA - 2.d0*(lwm + log(max(wsm, DTINY)))
        sainv_acc = sainv_acc + Ainv_
        do kk = 1, kmix
            sr_acc(kk)   = sr_acc(kk)   + rk(kk)
            sm_acc(:,kk) = sm_acc(:,kk) + rk(kk)*mk(:,kk)
            do r = 1, n
                smm_acc(:,r,kk) = smm_acc(:,r,kk) + rk(kk)*mk(:,kk)*mk(r,kk)
            end do
        end do
        end associate
    end subroutine probe_solve_mix


    !> MCFA initialisation. K=1 pins the single component at the origin over the current Gamma
    !! diagonal -- with the diagonal constraint in mcfa_condition this makes the mixture path
    !! reproduce the plain PPCA EM exactly (the regression test). K>=2 seeds by deterministic
    !! farthest-point selection on the current latents and polishes with Lloyd iterations;
    !! nothing here is random, so a rerun reproduces bit-for-bit.
    subroutine mcfa_init( z, nptcls, ncomp, kmix, gam_sum, nval, xi, ppi, Om )
        integer,  intent(in)  :: nptcls, ncomp, kmix, nval
        real(dp), intent(in)  :: z(nptcls,ncomp), gam_sum(ncomp)
        real(dp), intent(out) :: xi(ncomp,kmix), ppi(kmix), Om(ncomp,ncomp)
        integer,  allocatable :: rows(:), lab(:), cnt(:)
        real(dp), allocatable :: d2min(:)
        real(dp) :: best, dd
        integer  :: n, i, j, k, kbest, it2
        if( kmix == 1 )then
            xi  = 0.d0
            ppi = 1.d0
            Om  = 0.d0
            do i = 1, ncomp
                Om(i,i) = max(gam_sum(i)/real(max(1,nval),dp), DTINY)
            end do
            return
        endif
        ! particles the E-step skipped (state zero) keep z = 0 and must not seed a component
        allocate(rows(nptcls))
        n = 0
        do i = 1, nptcls
            if( any(z(i,1:ncomp) /= 0.d0) )then
                n = n + 1
                rows(n) = i
            endif
        end do
        allocate(lab(n), cnt(kmix), d2min(n))
        best = -1.d0
        j = 1
        do i = 1, n
            dd = sum(z(rows(i),1:ncomp)**2)
            if( dd > best )then
                best = dd
                j = i
            endif
        end do
        xi(:,1) = z(rows(j),1:ncomp)
        d2min = huge(0.d0)
        do k = 2, kmix
            best  = -1.d0
            kbest = 1
            do i = 1, n
                dd = sum((z(rows(i),1:ncomp) - xi(:,k-1))**2)
                if( dd < d2min(i) ) d2min(i) = dd
                if( d2min(i) > best )then
                    best  = d2min(i)
                    kbest = i
                endif
            end do
            xi(:,k) = z(rows(kbest),1:ncomp)
        end do
        do it2 = 1, 12
            cnt = 0
            do i = 1, n
                best  = huge(0.d0)
                kbest = 1
                do k = 1, kmix
                    dd = sum((z(rows(i),1:ncomp) - xi(:,k))**2)
                    if( dd < best )then
                        best  = dd
                        kbest = k
                    endif
                end do
                lab(i) = kbest
            end do
            xi = 0.d0
            do i = 1, n
                xi(:,lab(i)) = xi(:,lab(i)) + z(rows(i),1:ncomp)
                cnt(lab(i))  = cnt(lab(i)) + 1
            end do
            do k = 1, kmix
                if( cnt(k) > 0 ) xi(:,k) = xi(:,k)/real(cnt(k),dp)
            end do
        end do
        ! mixing proportions floored so no component is born dead
        do k = 1, kmix
            ppi(k) = max(real(cnt(k),dp)/real(max(1,n),dp), 0.25d0/real(kmix,dp))
        end do
        ppi = ppi / sum(ppi)
        ! tied Omega from the pooled within-cluster scatter of the POINT latents; the M-step
        ! adds the posterior covariance back from its own accumulators, so this only has to
        ! carry a sane scale, not be unbiased
        Om = 0.d0
        do i = 1, n
            do j = 1, ncomp
                Om(:,j) = Om(:,j) + (z(rows(i),1:ncomp) - xi(:,lab(i)))*(z(rows(i),j) - xi(j,lab(i)))
            end do
        end do
        Om = Om / real(max(1,n),dp)
        deallocate(rows, lab, cnt, d2min)
    end subroutine mcfa_init


    !> Symmetrise, eigen-floor and invert the tied prior covariance. diag_only enforces the
    !! K=1 reduction-test constraint (matching the plain path's diagonal Gamma exactly).
    subroutine mcfa_condition( ncomp, diag_only, Om, Ominv, ldOm )
        integer,  intent(in)    :: ncomp
        logical,  intent(in)    :: diag_only
        real(dp), intent(inout) :: Om(ncomp,ncomp)
        real(dp), intent(out)   :: Ominv(ncomp,ncomp), ldOm
        real(dp) :: ev(ncomp), evec(ncomp,ncomp), work(ncomp,ncomp), floor_ev
        integer  :: q, r2, nrot
        if( diag_only )then
            do q = 1, ncomp
                do r2 = 1, ncomp
                    if( q /= r2 ) Om(q,r2) = 0.d0
                end do
            end do
            Ominv = 0.d0
            ldOm  = 0.d0
            do q = 1, ncomp
                Om(q,q)    = max(Om(q,q), DTINY)
                Ominv(q,q) = 1.d0/Om(q,q)
                ldOm       = ldOm + log(Om(q,q))
            end do
            return
        endif
        Om   = 0.5d0*(Om + transpose(Om))
        work = Om
        call jacobi(work, ncomp, ncomp, ev, evec, nrot)
        floor_ev = max(1.d-6*maxval(ev), DTINY)
        ldOm = 0.d0
        do q = 1, ncomp
            ev(q) = max(ev(q), floor_ev)
            ldOm  = ldOm + log(ev(q))
        end do
        do q = 1, ncomp
            do r2 = 1, ncomp
                Om(q,r2)    = sum(evec(q,:)*ev(:)*evec(r2,:))
                Ominv(q,r2) = sum(evec(q,:)/ev(:)*evec(r2,:))
            end do
        end do
    end subroutine mcfa_condition


    !> MCFA M-step from REDUCED sufficient statistics (the thread reduction happens at
    !! the call site, so a rotated running-averaged history can be substituted for the
    !! this-iteration statistics transparently):
    !!   pi_k = Sr_k/N,  xi_k = Sm_k/Sr_k,
    !!   Omega = (1/N)[ sum_i A_i^-1 + sum_k (Smm_k - xi Sm' - Sm xi' + Sr xi xi') ]
    !! pin_origin (K=1) keeps xi at 0 and Omega diagonal -- the plain-PPCA reduction.
    !> The MCFA M-step from the reduced sufficient statistics (sr, sm, smm, sai) into the fit's
    !! mixture (weights, means, the tied covariance and its inverse and log-determinant).
    subroutine mcfa_mstep( fit, sr, sm, smm, sai )
        class(flex_fit), intent(inout) :: fit
        real(dp),        intent(in)    :: sr(fit%spec%kmix), sm(fit%model%ncomp,fit%spec%kmix)
        real(dp),        intent(in)    :: smm(fit%model%ncomp,fit%model%ncomp,fit%spec%kmix), sai(fit%model%ncomp,fit%model%ncomp)
        integer  :: k, r2
        associate( ncomp => fit%model%ncomp, kmix => fit%spec%kmix, nval => fit%iter%nval, pin_origin => (fit%spec%kmix == 1), &
            &ppi => fit%history%mix_pi, xi => fit%history%mix_xi, Om => fit%history%mix_Om, Ominv => fit%history%mix_Ominv, ldOm => fit%history%ldOm_mix )
        do k = 1, kmix
            ppi(k) = max(sr(k)/real(max(1,nval),dp), 1.d-4)
        end do
        ppi = ppi / sum(ppi)
        if( .not. pin_origin )then
            do k = 1, kmix
                ! a starved component keeps its old mean rather than dividing by ~0
                if( sr(k) > 1.d-8*real(max(1,nval),dp) ) xi(:,k) = sm(:,k)/sr(k)
            end do
        endif
        Om = sai
        do k = 1, kmix
            do r2 = 1, ncomp
                Om(:,r2) = Om(:,r2) + smm(:,r2,k) - xi(:,k)*sm(r2,k) - sm(:,k)*xi(r2,k) &
                    &+ sr(k)*xi(:,k)*xi(r2,k)
            end do
        end do
        Om = Om / real(max(1,nval),dp)
        call mcfa_condition(ncomp, pin_origin, Om, Ominv, ldOm)
        end associate
    end subroutine mcfa_mstep



    !> Sampling precision of the MAP latent estimate, Q = A*Gtil^+*A with A = Gtil + diag(prior). This is
    !! the precision of the ESTIMATOR z_hat, not the posterior precision A, so distances measured with it
    !! reflect how well each component was actually determined for the particle.
    subroutine map_sampling_precision( Gtil, prior, n, Qout )
        integer,  intent(in)  :: n
        real(dp), intent(in)  :: Gtil(n,n), prior(n)
        real(dp), intent(out) :: Qout(n,n)
        real(dp) :: Amat(n,n), Gpinv(n,n), Vmat(n,n), Awork(n,n), ev(n), thresh
        integer  :: ii, jj, kk, nrot
        Amat = Gtil
        do ii = 1, n
            Amat(ii,ii) = Amat(ii,ii) + prior(ii)
        end do
        Awork = Gtil
        call jacobi(Awork, n, n, ev, Vmat, nrot)   ! symmetric eigendecomposition (LAPACK dsyev)
        thresh = COV_PINV_RCOND * maxval(abs(ev))
        Gpinv  = 0.d0
        do kk = 1, n
            if( abs(ev(kk)) <= thresh ) cycle      ! drop the null space, as pinv does
            do jj = 1, n
                do ii = 1, n
                    Gpinv(ii,jj) = Gpinv(ii,jj) + Vmat(ii,kk)*Vmat(jj,kk)/ev(kk)
                end do
            end do
        end do
        Qout = matmul(Amat, matmul(Gpinv, Amat))
        Qout = 0.5d0*(Qout + transpose(Qout))      ! symmetrise away round-off
    end subroutine map_sampling_precision


end module simple_flex_pca_posterior
