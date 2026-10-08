!@descr: flex_pca posterior inference: the per-particle latent posterior from sufficient statistics
!!
!! Consumes what an E-step formulation produces for one particle -- the projected Gram G, the
!! data and mean rows b and c, the mean/data scalars -- and returns the posterior: the MAP
!! latent (or the mixture-weighted mean), its precision/covariance, log det, likelihood
!! contribution at the fixed per-particle contrast, plus the mixture M-step and its
!! initialisation and the estimator precision of the MAP embedding. Dense algebra on arrays;
!! no projection code, no fit state, no file.
module simple_flex_pca_posterior
use simple_core_module_api,    only: dp, dtiny, simple_exception
use simple_linalg,             only: jacobi, cholesky, chol_forward, chol_backward, spd_inverse, spd_logdet
use simple_flex_pca_fit_types, only: flex_fit
use simple_kmeans,             only: kmeans
implicit none
private
#include "simple_local_flags.inc"

public :: quad_form, spd_solve_dp, spd_inv_dp
public :: probe_solve_plain, probe_solve_mix
public :: mcfa_init, mcfa_condition, mcfa_mstep
public :: map_sampling_precision

real(dp), parameter :: COV_PINV_RCOND = 1.0d-6

contains

    !>  z' M z for symmetric M.
    pure function quad_form( M, z, n ) result( val )
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

    !> In-place symmetric positive-definite solve A x = b (b overwritten by x) via Cholesky (simple_linalg).
    !! A is first scaled by its mean diagonal, so the retry ridge is RELATIVE: an absolute ridge either
    !! swamps a small-diagonal system or fails to rescue a large one, and the b=0 fallback then collapses
    !! the latents of essentially every particle. A is overwritten.
    subroutine spd_solve_dp( A, b, n )
        integer,  intent(in)    :: n
        real(dp), intent(inout) :: A(n,n), b(n)
        real(dp) :: L(n,n), y(n), ridge, dscale
        integer  :: i, attempt
        logical  :: ok
        dscale = mean_abs_diag(A, n)
        if( dscale > 0.d0 )then
            A = A / dscale
            b = b / dscale
        endif
        do attempt = 1, 3
            call cholesky(A, L, n, ok)
            if( ok )then
                call chol_forward(L, b, y, n)
                call chol_backward(L, y, b, n)
                return
            endif
            ridge = 1.d-8 * (abs(A(1,1))+1.d0) * (10.d0**(attempt-1))
            do i = 1, n
                A(i,i) = A(i,i) + ridge
            end do
        end do
        b = 0.d0
    end subroutine spd_solve_dp

    !> SPD inverse by Cholesky (simple_linalg), same rescaling and ridge escalation as spd_solve_dp;
    !! zeros if all attempts fail. A is overwritten.
    subroutine spd_inv_dp( A, Ainv, n )
        integer,  intent(in)    :: n
        real(dp), intent(inout) :: A(n,n)
        real(dp), intent(out)   :: Ainv(n,n)
        real(dp) :: ridge, dscale
        integer  :: i, attempt
        logical  :: ok
        Ainv   = 0.d0
        dscale = mean_abs_diag(A, n)
        if( dscale > 0.d0 ) A = A / dscale
        do attempt = 1, 3
            call spd_inverse(A, Ainv, n, ok)
            if( ok )then
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

    pure real(dp) function mean_abs_diag( A, n ) result( dscale )
        integer,  intent(in) :: n
        real(dp), intent(in) :: A(n,n)
        integer :: i
        dscale = 0.d0
        do i = 1, n
            dscale = dscale + abs(A(i,i))
        end do
        dscale = dscale / real(n,dp)
    end function mean_abs_diag

    !> One particle's posterior solve under the plain 1/Gamma prior at fixed contrast: the
    !! fit's E-step statistics of thread ithr in, the latent and its posterior covariance out.
    subroutine probe_solve_plain( fit, ithr, a, ldA, lok, quad )
        class(flex_fit), intent(inout) :: fit
        integer,         intent(in)    :: ithr
        real(dp),        intent(in)    :: a
        real(dp),        intent(out)   :: ldA
        logical,         intent(out)   :: lok
        real(dp),        intent(out)   :: quad
        real(dp) :: Amat(fit%model%ncomp,fit%model%ncomp)
        real(dp) :: Acp(fit%model%ncomp,fit%model%ncomp), h(fit%model%ncomp), aa
        integer  :: q
        associate( n => fit%model%ncomp, sig2 => fit%model%sig2, G => fit%iter%Gth(:,:,ithr), b => fit%iter%bth(:,ithr), &
            &c => fit%iter%cth(:,ithr), prior_ => fit%iter%prior, z_ => fit%iter%zth(:,ithr), Ainv_ => fit%iter%Ainvth(:,:,ithr) )
            aa   = a*a
            Amat = (aa/sig2)*G
            do q = 1, n
                Amat(q,q) = Amat(q,q) + prior_(q)
                z_(q)     = (a*b(q) - aa*c(q))/sig2
            end do
            Acp = Amat
            call spd_logdet(Amat, n, ldA, lok)
            h = z_
            call spd_inv_dp(Acp, Ainv_, n)
            call spd_solve_dp(Amat, z_, n)
            quad = dot_product(h, z_)
        end associate
    end subroutine probe_solve_plain

    !> Per-particle MCFA posterior from fetched sufficient statistics (G, b, c) -- the ONE
    !! source for the mixture E-step, shared by the plain CPU body and the fused device
    !! body (the device computes the same G/b/c on card; only this host solve differs from
    !! the single-Gaussian path). Updates the caller's thread-slice accumulators exactly as
    !! the historical in-line block did.
    !> One particle's posterior solve under the MCFA mixture prior: the fit's E-step statistics of
    !! thread ithr in; the latent mean, the second moment of batch row i, the thread's mixture
    !! accumulators and the likelihood term out.
    subroutine probe_solve_mix( fit, ithr, i, a, ldA, lok, nll_add )
        class(flex_fit), intent(inout) :: fit
        integer,         intent(in)    :: ithr, i
        real(dp),        intent(in)    :: a
        real(dp),        intent(out)   :: ldA, nll_add
        logical,         intent(out)   :: lok
        real(dp) :: Amat(fit%model%ncomp,fit%model%ncomp), Acp(fit%model%ncomp,fit%model%ncomp)
        real(dp) :: rhs0(fit%model%ncomp), mk(fit%model%ncomp,fit%spec%kmix)
        real(dp) :: lw(fit%spec%kmix), rk(fit%spec%kmix)
        real(dp) :: aa, lwm, wsm
        integer  :: q, r, kk
        associate( n => fit%model%ncomp, kmix => fit%spec%kmix, sig2 => fit%model%sig2, &
            &G => fit%iter%Gth(:,:,ithr), b => fit%iter%bth(:,ithr), c => fit%iter%cth(:,ithr), &
            &Ominv => fit%history%mix_Ominv, Omxi => fit%history%mix_Omxi, &
            &lpi => fit%history%mix_lpi, xiOx => fit%history%mix_xiOx, &
            &zbar => fit%iter%zth(:,ithr), Ainv_ => fit%iter%Ainvth(:,:,ithr), dens_ => fit%iter%dens(:,:,i), &
            &sr_acc => fit%history%mxa_sr(:,ithr), sm_acc => fit%history%mxa_sm(:,:,ithr), &
            &smm_acc => fit%history%mxa_smm(:,:,:,ithr), &
            &sainv_acc => fit%history%mxa_sainv(:,:,ithr) )
            aa   = a*a
            Amat = (aa/sig2)*G + Ominv
            Acp  = Amat
            call spd_logdet(Amat, n, ldA, lok)
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
    !! reproduce the plain PPCA EM exactly (the regression test). K>=2 clusters the current latents
    !! by deterministic k-means (simple_kmeans: seeded at the point nearest the mean, then
    !! farthest-point, Lloyd to convergence, empty clusters recovered); nothing here is random, so a
    !! rerun reproduces bit-for-bit.
    subroutine mcfa_init( z, nptcls, ncomp, kmix, gam_sum, nval, xi, ppi, Om )
        integer,  intent(in)  :: nptcls, ncomp, kmix, nval
        real(dp), intent(in)  :: z(nptcls,ncomp), gam_sum(ncomp)
        real(dp), intent(out) :: xi(ncomp,kmix), ppi(kmix), Om(ncomp,ncomp)
        type(kmeans)          :: km
        real(dp), allocatable :: zrows(:,:)
        integer,  allocatable :: lab(:), cnt(:)
        integer  :: n, i, j, k
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
        n = 0
        do i = 1, nptcls
            if( any(z(i,1:ncomp) /= 0.d0) ) n = n + 1
        end do
        if( n < kmix ) THROW_HARD('fewer latents with signal than mixture components; mcfa_init')
        allocate(zrows(n,ncomp), lab(n), cnt(kmix))
        n = 0
        do i = 1, nptcls
            if( any(z(i,1:ncomp) /= 0.d0) )then
                n = n + 1
                zrows(n,:) = z(i,1:ncomp)
            endif
        end do
        call km%new(zrows, kmix)
        call km%cluster(lab, xi)
        call km%kill
        cnt = 0
        do i = 1, n
            cnt(lab(i)) = cnt(lab(i)) + 1
        end do
        ! mixing proportions floored so no component is born dead
        do k = 1, kmix
            ppi(k) = max(real(cnt(k),dp)/real(n,dp), 0.25d0/real(kmix,dp))
        end do
        ppi = ppi / sum(ppi)
        ! tied Omega from the pooled within-cluster scatter of the POINT latents; the M-step
        ! adds the posterior covariance back from its own accumulators, so this only has to
        ! carry a sane scale, not be unbiased
        Om = 0.d0
        do i = 1, n
            do j = 1, ncomp
                Om(:,j) = Om(:,j) + (zrows(i,:) - xi(:,lab(i)))*(zrows(i,j) - xi(j,lab(i)))
            end do
        end do
        Om = Om / real(n,dp)
        deallocate(zrows, lab, cnt)
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
