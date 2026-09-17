!@descr: flex_pca latent deconvolution: calibrated per-particle noise + an empirical-Bayes mixture prior fitted through it
!!
!! The delivered latent z_i is a MAP solve: z_i = A_i^-1 (D_i x_i + e_i) with D_i the data precision
!! (the whitened projected Gram over sig2), A_i = D_i + P the posterior precision (P = diag prior),
!! and e_i whitened image noise of covariance D_i. So E[z_i|x_i] = R_i x_i with R_i = A_i^-1 D_i and
!! Cov[z_i|x_i] = A_i^-1 D_i A_i^-1 =: N_i. The stored precision A_i understates the real noise
!! (basis error, model misfit), so N_i is scaled by ONE factor a, calibrated from the even/odd
!! Fourier-half solutions of every particle: d_i = z_i^even - z_i^odd has predicted covariance
!! A_h^-1 D_i A_h^-1 with A_h = P + D_i/2, and a = sum|d_i|^2 / sum tr(predicted). No likelihood
!! chooses a.
!!
!! The population prior G is a K-component Gaussian mixture fitted by extreme deconvolution
!! (Bovy, Hogg & Roweis 2011): z_i ~ sum_k pi_k N(R_i mu_k, R_i Sigma_k R_i^T + a N_i). K is chosen
!! by particle-half cross-validation (fit on one half of the particles, score the held-out
!! log-likelihood of the other, both ways). The delivered coordinates become the posterior means
!! E[x_i|z_i] under G, with their posterior covariances, so every downstream consumer (UMAP, state
!! placement, kernels) works on the population instead of the noise ellipsoid.
module simple_flex_pca_deconv
use simple_core_module_api
use simple_srch_sort_loc, only: hpsort
implicit none
private
#include "simple_local_flags.inc"

public :: calibrate_noise_scale, deconvolve_latent, test_flex_pca_deconv
public :: noise_and_projection, signal_subspace

integer,  parameter :: XD_MAXIT   = 150
real(dp), parameter :: XD_TOL     = 1.d-6      !< relative log-likelihood change
real(dp), parameter :: XD_CV_MAX  = 20000       !< particles in the K-selection ladder (both halves together)
real(dp), parameter :: XD_RIDGE   = 1.d-6      !< Sigma_k ridge, relative to tr(C_obs)/d
real(dp), parameter :: XD_PI_MIN  = 1.d-4
integer,  parameter :: XD_KMAX    = 16

contains

    ! ======================================================================= calibration

    !> a = sum_i |z_i^even - z_i^odd|^2 / sum_i tr( A_h^-1 D_i A_h^-1 ),  A_h = P + D_i/2, D_i = A_i - P
    subroutine calibrate_noise_scale( z, zhalf, precision, prior, nptcls, ncomp, a, a_comp )
        integer,  intent(in)  :: nptcls, ncomp
        real(dp), intent(in)  :: z(nptcls,ncomp), zhalf(nptcls,ncomp,2)
        real(dp), intent(in)  :: precision(ncomp,ncomp,nptcls), prior(ncomp)
        real(dp), intent(out) :: a, a_comp(ncomp)
        real(dp) :: Ah(ncomp,ncomp), Ahinv(ncomp,ncomp), Dm(ncomp,ncomp), C(ncomp,ncomp), dvec(ncomp)
        real(dp) :: num, den, num_q(ncomp), den_q(ncomp)
        integer  :: i, q
        logical  :: ok
        num = 0.d0; den = 0.d0; num_q = 0.d0; den_q = 0.d0
        !$omp parallel do default(shared) private(i,q,Ah,Ahinv,Dm,C,dvec,ok) schedule(static) &
        !$omp& reduction(+:num,den,num_q,den_q)
        do i = 1, nptcls
            dvec = zhalf(i,:,1) - zhalf(i,:,2)
            if( all(abs(dvec) < DTINY) ) cycle          ! no half solve for this particle
            Dm = precision(:,:,i)
            do q = 1, ncomp
                Dm(q,q) = Dm(q,q) - prior(q)
            end do
            Ah = 0.5d0*Dm
            do q = 1, ncomp
                Ah(q,q) = Ah(q,q) + prior(q)
            end do
            call spd_inverse(Ah, Ahinv, ncomp, ok)
            if( .not. ok ) cycle
            C = matmul(Ahinv, matmul(Dm, Ahinv))
            num = num + sum(dvec*dvec)
            do q = 1, ncomp
                den      = den + C(q,q)
                num_q(q) = num_q(q) + dvec(q)*dvec(q)
                den_q(q) = den_q(q) + C(q,q)
            end do
        end do
        !$omp end parallel do
        a = 1.d0
        if( den > DTINY ) a = num/den
        do q = 1, ncomp
            a_comp(q) = 1.d0
            if( den_q(q) > DTINY ) a_comp(q) = num_q(q)/den_q(q)
        end do
        write(logfhandle,'(A,F8.3)') '>>> FLEX_PCA NOISE CALIBRATION (even/odd half solves vs the model): a=', a
        write(logfhandle,'(A)',advance='no') '>>>   per component:'
        do q = 1, ncomp
            write(logfhandle,'(1X,F6.2)',advance='no') a_comp(q)
        end do
        write(logfhandle,'(A)') ''
        call flush(logfhandle)
    end subroutine calibrate_noise_scale

    ! ======================================================================= deconvolution

    !> Fit the mixture prior through the calibrated noise, choose K by particle-half cross-validation,
    !! and replace z by the posterior means and precision by the posterior precisions.
    subroutine deconvolve_latent( z, precision, prior, nptcls, ncomp, a, kmax_in, k_out, prior_fname, labels_fname, pinds, labels_out )
        integer,  intent(in)    :: nptcls, ncomp
        real(dp), intent(inout) :: z(nptcls,ncomp)
        real(dp), intent(inout) :: precision(ncomp,ncomp,nptcls)
        real(dp), intent(in)    :: prior(ncomp), a
        integer,  intent(in)    :: kmax_in
        integer,  intent(out)   :: k_out
        character(len=*), optional, intent(in) :: prior_fname
        character(len=*), optional, intent(in) :: labels_fname   !< per-particle argmax component + its responsibility
        integer,          optional, intent(in) :: pinds(:)        !< project rows for that file
        integer, allocatable, optional, intent(out) :: labels_out(:)  !< argmax component per particle
        real(dp), allocatable :: R(:,:,:), Nz(:,:,:), mu(:,:), Sig(:,:,:), pik(:), xhat(:,:), xcov(:,:,:)
        real(dp), allocatable :: mu_h(:,:), Sig_h(:,:,:), pik_h(:), score(:), resp(:,:)
        real(dp) :: ll, ll_h, cobs(ncomp,ncomp), zmean(ncomp), vpop(ncomp), vobs(ncomp), mbar(ncomp)
        integer  :: kmax, kcv, i, q, k, nh, u, ndrop
        logical, allocatable :: inA(:), inB(:)
        integer :: stride
        integer(timer_int_kind) :: t0
        t0 = tic()
        allocate(R(ncomp,ncomp,nptcls), Nz(ncomp,ncomp,nptcls))
        call noise_and_projection(precision, prior, nptcls, ncomp, a, R, Nz)
        ! observed moments (for the init, the ridge and the report)
        do q = 1, ncomp
            zmean(q) = sum(z(:,q))/real(nptcls,dp)
        end do
        do q = 1, ncomp
            do k = 1, ncomp
                cobs(q,k) = sum((z(:,q)-zmean(q))*(z(:,k)-zmean(k)))/real(nptcls,dp)
            end do
        end do
        ! ---- K by particle-half cross-validation ----
        kmax = max(1, min(kmax_in, XD_KMAX, nptcls/2000))
        ! K is a coarse choice: select it on a strided subsample (measured 1178 s for the ladder on
        ! 105k particles in 17 dimensions; the final fit below still sees every particle)
        allocate(inA(nptcls), inB(nptcls), score(kmax))
        stride = max(1, nint(real(nptcls,dp)/XD_CV_MAX))
        do i = 1, nptcls
            inA(i) = mod(i, 2*stride) == 1
            inB(i) = mod(i, 2*stride) == mod(stride + 1, 2*stride)
        end do
        nh = count(inA)
        write(logfhandle,'(A,I0,A,I0,A)') '>>> FLEX_PCA DECONV K selection on ', count(inA) + count(inB), &
            &' of ', nptcls, ' particles (two disjoint halves); final fit on all'
        score = -huge(0.d0)
        ndrop = 0
        do kcv = 1, kmax
            score(kcv) = 0.d0
            call fit_half(kcv, inA, inB, ll_h);  score(kcv) = score(kcv) + ll_h
            call fit_half(kcv, inB, inA, ll_h);  score(kcv) = score(kcv) + ll_h
            score(kcv) = score(kcv)/real(count(inA) + count(inB),dp)
            write(logfhandle,'(A,I2,A,F14.6)') '>>> FLEX_PCA DECONV CV K=', kcv, '  held-out loglik/particle=', score(kcv)
            call flush(logfhandle)
            if( kcv >= 2 )then
                if( score(kcv) < score(kcv-1) )then
                    ndrop = ndrop + 1
                else
                    ndrop = 0
                endif
                if( ndrop >= 2 ) exit
            endif
        end do
        k_out = maxloc(score(1:min(kcv,kmax)), dim=1)
        ! ---- final fit on all particles ----
        call xd_ladder(z, R, Nz, nptcls, ncomp, k_out, cobs, zmean, mu, Sig, pik, ll)
        allocate(xhat(nptcls,ncomp), xcov(ncomp,ncomp,nptcls), resp(nptcls,k_out))
        call xd_posterior(z, R, Nz, nptcls, ncomp, k_out, mu, Sig, pik, xhat, xcov, resp)
        ! ---- report: population variance per axis vs observed ----
        mbar = 0.d0
        do k = 1, k_out
            mbar = mbar + pik(k)*mu(:,k)
        end do
        do q = 1, ncomp
            vpop(q) = 0.d0
            do k = 1, k_out
                vpop(q) = vpop(q) + pik(k)*(Sig(q,q,k) + (mu(q,k)-mbar(q))**2)
            end do
            vobs(q) = cobs(q,q)
        end do
        write(logfhandle,'(A,I0,A,F14.6,A,F8.1)') '>>> FLEX_PCA DECONV: K=', k_out, '  loglik/particle=', &
            &ll/real(nptcls,dp), '  seconds=', toc(t0)
        write(logfhandle,'(A)') '>>> FLEX_PCA DECONV population variance per axis (deconvolved / observed = signal fraction):'
        do q = 1, ncomp
            write(logfhandle,'(A,I3,A,ES11.3,A,ES11.3,A,F7.3)') '>>>   z', q, '  pop=', vpop(q), &
                &'  obs=', vobs(q), '  fraction=', vpop(q)/max(vobs(q), DTINY)
        end do
        do k = 1, k_out
            write(logfhandle,'(A,I2,A,F7.4,A,ES11.3)') '>>>   component ', k, '  weight=', pik(k), &
                &'  tr(Sigma)=', sum([(Sig(q,q,k), q=1,ncomp)])
        end do
        call flush(logfhandle)
        if( present(prior_fname) )then
            open(newunit=u, file=prior_fname, status='replace', action='write')
            write(u,'(A,I0,A,I0,A,F10.4)') '# flex_pca deconvolved prior: K=', k_out, ' ncomp=', ncomp, ' noise_scale=', a
            write(u,'(A)') '# rows: k weight mu(1:ncomp) diag(Sigma)(1:ncomp)'
            do k = 1, k_out
                write(u,'(I3,1X,F8.5,*(1X,ES14.6))') k, pik(k), mu(:,k), [(Sig(q,q,k), q=1,ncomp)]
            end do
            close(u)
        endif
        if( present(labels_fname) )then
            open(newunit=u, file=labels_fname, status='replace', action='write')
            write(u,'(A)') '# particle  component(argmax responsibility)  responsibility'
            do i = 1, nptcls
                k = maxloc(resp(i,:), dim=1)
                if( present(pinds) )then
                    write(u,'(I10,1X,I3,1X,F8.5)') pinds(i), k, resp(i,k)
                else
                    write(u,'(I10,1X,I3,1X,F8.5)') i, k, resp(i,k)
                endif
            end do
            close(u)
        endif
        ! ---- hand back ----
        if( present(labels_out) )then
            if( allocated(labels_out) ) deallocate(labels_out)
            allocate(labels_out(nptcls))
            do i = 1, nptcls
                labels_out(i) = maxloc(resp(i,:), dim=1)
            end do
        endif
        z = xhat
        do i = 1, nptcls
            call spd_inverse(xcov(:,:,i), precision(:,:,i), ncomp)
        end do
        deallocate(R, Nz, mu, Sig, pik, xhat, xcov, inA, inB, score, resp)

    contains

        !> fit nk components on the particles flagged in sel, return the held-out log-likelihood of those in tst
        subroutine fit_half( nk, sel, tst, ll_held )
            integer, intent(in)  :: nk
            logical, intent(in)  :: sel(nptcls), tst(nptcls)
            real(dp), intent(out) :: ll_held
            real(dp), allocatable :: zs(:,:), Rs(:,:,:), Nsel(:,:,:), zt(:,:), Rt(:,:,:), Ntst(:,:,:)
            real(dp) :: llfit, cs(ncomp,ncomp), ms(ncomp)
            integer  :: ns, nt, j, jj, qq, kk
            ns = count(sel); nt = count(tst)
            allocate(zs(ns,ncomp), Rs(ncomp,ncomp,ns), Nsel(ncomp,ncomp,ns), zt(nt,ncomp), Rt(ncomp,ncomp,nt), Ntst(ncomp,ncomp,nt))
            j = 0; jj = 0
            do i = 1, nptcls
                if( sel(i) )then
                    j = j + 1; zs(j,:) = z(i,:); Rs(:,:,j) = R(:,:,i); Nsel(:,:,j) = Nz(:,:,i)
                else if( tst(i) )then
                    jj = jj + 1; zt(jj,:) = z(i,:); Rt(:,:,jj) = R(:,:,i); Ntst(:,:,jj) = Nz(:,:,i)
                endif
            end do
            do qq = 1, ncomp
                ms(qq) = sum(zs(:,qq))/real(ns,dp)
            end do
            do qq = 1, ncomp
                do kk = 1, ncomp
                    cs(qq,kk) = sum((zs(:,qq)-ms(qq))*(zs(:,kk)-ms(kk)))/real(ns,dp)
                end do
            end do
            allocate(mu_h(ncomp,nk), Sig_h(ncomp,ncomp,nk), pik_h(nk))
            call xd_init(zs, ns, ncomp, nk, cs, ms, Nsel, mu_h, Sig_h, pik_h)
            call xd_fit(zs, Rs, Nsel, ns, ncomp, nk, mu_h, Sig_h, pik_h, llfit)
            ll_held = xd_loglik(zt, Rt, Ntst, nt, ncomp, nk, mu_h, Sig_h, pik_h)
            deallocate(zs, Rs, Nsel, zt, Rt, Ntst, mu_h, Sig_h, pik_h)
        end subroutine fit_half

    end subroutine deconvolve_latent

    !> R_i = A_i^-1 D_i and N_i = a A_i^-1 D_i A_i^-1 from the posterior precision A_i and the prior P
    subroutine noise_and_projection( precision, prior, nptcls, ncomp, a, R, N )
        integer,  intent(in)  :: nptcls, ncomp
        real(dp), intent(in)  :: precision(ncomp,ncomp,nptcls), prior(ncomp), a
        real(dp), intent(out) :: R(ncomp,ncomp,nptcls), N(ncomp,ncomp,nptcls)
        real(dp) :: Ainv(ncomp,ncomp), D(ncomp,ncomp)
        integer  :: i, q
        logical  :: ok
        !$omp parallel do default(shared) private(i,q,Ainv,D,ok) schedule(static)
        do i = 1, nptcls
            call spd_inverse(precision(:,:,i), Ainv, ncomp, ok)
            D = precision(:,:,i)
            do q = 1, ncomp
                D(q,q) = D(q,q) - prior(q)
            end do
            if( ok )then
                R(:,:,i) = matmul(Ainv, D)
                N(:,:,i) = a*matmul(R(:,:,i), Ainv)
                N(:,:,i) = 0.5d0*(N(:,:,i) + transpose(N(:,:,i)))
            else
                R(:,:,i) = 0.d0
                N(:,:,i) = 0.d0
                do q = 1, ncomp
                    R(q,q,i) = 1.d0
                    N(q,q,i) = a/max(precision(q,q,i), DTINY)
                end do
            endif
        end do
        !$omp end parallel do
    end subroutine noise_and_projection

    !> fit nk components from the equal-mass quantile init (seeding new components on the
    !! worst-explained particles was tried on 2026-09-08 and catches noise outliers, not compact
    !! minority populations: five of eight components tiny and 20x wider than the body)
    subroutine xd_ladder( z, R, Nz, n, d, nk, cobs, zmean, mu, Sig, pik, ll )
        integer,  intent(in)  :: n, d, nk
        real(dp), intent(in)  :: z(n,d), R(d,d,n), Nz(d,d,n), cobs(d,d), zmean(d)
        real(dp), allocatable, intent(out) :: mu(:,:), Sig(:,:,:), pik(:)
        real(dp), intent(out) :: ll
        allocate(mu(d,nk), Sig(d,d,nk), pik(nk))
        call xd_init(z, n, d, nk, cobs, zmean, Nz, mu, Sig, pik)
        call xd_fit(z, R, Nz, n, d, nk, mu, Sig, pik, ll)
    end subroutine xd_ladder

    !> equal-mass quantile means along the axis of largest observed variance; Sigma_k from the
    !! deconvolved observed covariance (PSD-clipped), pi uniform
    subroutine xd_init( z, n, d, nk, cobs, zmean, Nz, mu, Sig, pik )
        integer,  intent(in)  :: n, d, nk
        real(dp), intent(in)  :: z(n,d), cobs(d,d), zmean(d), Nz(d,d,n)
        real(dp), intent(out) :: mu(d,nk), Sig(d,d,nk), pik(nk)
        real,     allocatable :: key(:)
        integer,  allocatable :: ord(:)
        real(dp) :: S(d,d), Nbar(d,d), ridge
        integer  :: i, k, i0, i1, q, jbest
        jbest = maxloc([(cobs(q,q), q=1,d)], dim=1)
        allocate(key(n), ord(n))
        do i = 1, n
            key(i) = real(z(i,jbest)); ord(i) = i
        end do
        call hpsort(key, ord)
        mu = 0.d0
        do k = 1, nk
            i0 = nint(real(k-1,dp)*real(n,dp)/real(nk,dp)) + 1
            i1 = max(nint(real(k,dp)*real(n,dp)/real(nk,dp)), i0)
            do i = i0, i1
                mu(:,k) = mu(:,k) + z(ord(i),:)
            end do
            mu(:,k) = mu(:,k)/real(i1-i0+1,dp)
        end do
        Nbar = 0.d0
        do i = 1, n
            Nbar = Nbar + Nz(:,:,i)
        end do
        Nbar = Nbar/real(n,dp)
        S = cobs - Nbar
        ridge = XD_RIDGE*max(sum([(cobs(q,q), q=1,d)])/real(d,dp), DTINY)
        call psd_clip(S, d, ridge)
        do k = 1, nk
            Sig(:,:,k) = S
            pik(k)     = 1.d0/real(nk,dp)
        end do
        deallocate(key, ord)
    end subroutine xd_init

    !> extreme-deconvolution EM with per-particle projection R_i and noise N_i
    subroutine xd_fit( z, R, Nz, n, d, nk, mu, Sig, pik, ll )
        integer,  intent(in)    :: n, d, nk
        real(dp), intent(in)    :: z(n,d), R(d,d,n), Nz(d,d,n)
        real(dp), intent(inout) :: mu(d,nk), Sig(d,d,nk), pik(nk)
        real(dp), intent(out)   :: ll
        real(dp), allocatable :: sw(:,:), sb(:,:,:), sBB(:,:,:,:)
        real(dp) :: T(d,d), Tinv(d,d), SR(d,d), b(d), Bk(d,d), rz(d), logp(nk), lmax, lse, rk(nk)
        real(dp) :: ll_prev, ridge, trc
        integer  :: it, i, k, ithr, nthr, q
        logical  :: ok
        nthr = omp_get_max_threads()
        allocate(sw(nk,nthr), sb(d,nk,nthr), sBB(d,d,nk,nthr))
        trc = 0.d0
        do k = 1, nk
            trc = trc + sum([(Sig(q,q,k), q=1,d)])
        end do
        ridge = XD_RIDGE*max(trc/real(d*nk,dp), DTINY)
        ll_prev = -huge(0.d0)
        do it = 1, XD_MAXIT
            sw = 0.d0; sb = 0.d0; sBB = 0.d0
            ll = 0.d0
            !$omp parallel do default(shared) private(i,k,ithr,T,Tinv,SR,b,Bk,rz,logp,lmax,lse,rk,ok,q) &
            !$omp& schedule(static) reduction(+:ll)
            do i = 1, n
                ithr = omp_get_thread_num() + 1
                do k = 1, nk
                    call component_terms(z(i,:), R(:,:,i), Nz(:,:,i), mu(:,k), Sig(:,:,k), d, logp(k), ok)
                    logp(k) = logp(k) + log(max(pik(k), XD_PI_MIN))
                    if( .not. ok ) logp(k) = -huge(0.d0)/2
                end do
                lmax = maxval(logp)
                lse  = lmax + log(sum(exp(logp - lmax)))
                ll   = ll + lse
                rk   = exp(logp - lse)
                do k = 1, nk
                    if( rk(k) < 1.d-12 ) cycle
                    ! b_ik, B_ik
                    SR = matmul(Sig(:,:,k), transpose(R(:,:,i)))                  ! Sigma R^T
                    T  = matmul(R(:,:,i), SR) + Nz(:,:,i)                          ! R Sigma R^T + Nz
                    call spd_inverse(T, Tinv, d, ok)
                    if( .not. ok ) cycle
                    rz = z(i,:) - matmul(R(:,:,i), mu(:,k))
                    b  = mu(:,k) + matmul(SR, matmul(Tinv, rz))
                    Bk = Sig(:,:,k) - matmul(SR, matmul(Tinv, transpose(SR)))
                    sw(k,ithr)     = sw(k,ithr) + rk(k)
                    sb(:,k,ithr)   = sb(:,k,ithr) + rk(k)*b
                    sBB(:,:,k,ithr) = sBB(:,:,k,ithr) + rk(k)*(Bk + outer(b, b, d))
                end do
            end do
            !$omp end parallel do
            ! M-step
            do k = 1, nk
                if( sum(sw(k,:)) <= XD_PI_MIN*real(n,dp) )then
                    pik(k) = XD_PI_MIN
                    cycle
                endif
                pik(k)  = sum(sw(k,:))/real(n,dp)
                mu(:,k) = sum(sb(:,k,:), dim=2)/sum(sw(k,:))
                Sig(:,:,k) = sum(sBB(:,:,k,:), dim=3)/sum(sw(k,:)) - outer(mu(:,k), mu(:,k), d)
                Sig(:,:,k) = 0.5d0*(Sig(:,:,k) + transpose(Sig(:,:,k)))
                call psd_clip(Sig(:,:,k), d, ridge)
            end do
            pik = pik/sum(pik)
            if( abs(ll - ll_prev) <= XD_TOL*abs(ll) ) exit
            ll_prev = ll
        end do
        deallocate(sw, sb, sBB)
    end subroutine xd_fit

    !> log-likelihood of a set under a fitted prior
    function xd_loglik( z, R, Nz, n, d, nk, mu, Sig, pik ) result( ll )
        integer,  intent(in) :: n, d, nk
        real(dp), intent(in) :: z(n,d), R(d,d,n), Nz(d,d,n), mu(d,nk), Sig(d,d,nk), pik(nk)
        real(dp) :: ll, logp(nk), lmax
        integer  :: i, k
        logical  :: ok
        ll = 0.d0
        !$omp parallel do default(shared) private(i,k,logp,lmax,ok) schedule(static) reduction(+:ll)
        do i = 1, n
            do k = 1, nk
                call component_terms(z(i,:), R(:,:,i), Nz(:,:,i), mu(:,k), Sig(:,:,k), d, logp(k), ok)
                logp(k) = logp(k) + log(max(pik(k), XD_PI_MIN))
                if( .not. ok ) logp(k) = -huge(0.d0)/2
            end do
            lmax = maxval(logp)
            ll   = ll + lmax + log(sum(exp(logp - lmax)))
        end do
        !$omp end parallel do
    end function xd_loglik

    !> posterior mean and covariance of every particle under the fitted prior
    subroutine xd_posterior( z, R, Nz, n, d, nk, mu, Sig, pik, xhat, xcov, resp )
        integer,  intent(in)  :: n, d, nk
        real(dp), intent(in)  :: z(n,d), R(d,d,n), Nz(d,d,n), mu(d,nk), Sig(d,d,nk), pik(nk)
        real(dp), intent(out) :: xhat(n,d), xcov(d,d,n)
        real(dp), optional, intent(out) :: resp(n,nk)
        real(dp) :: T(d,d), Tinv(d,d), SR(d,d), b(d,nk), Bk(d,d,nk), rz(d), logp(nk), lmax, lse, rk(nk)
        integer  :: i, k
        logical  :: ok
        !$omp parallel do default(shared) private(i,k,T,Tinv,SR,b,Bk,rz,logp,lmax,lse,rk,ok) schedule(static)
        do i = 1, n
            do k = 1, nk
                call component_terms(z(i,:), R(:,:,i), Nz(:,:,i), mu(:,k), Sig(:,:,k), d, logp(k), ok)
                logp(k) = logp(k) + log(max(pik(k), XD_PI_MIN))
                if( .not. ok ) logp(k) = -huge(0.d0)/2
                SR = matmul(Sig(:,:,k), transpose(R(:,:,i)))
                T  = matmul(R(:,:,i), SR) + Nz(:,:,i)
                call spd_inverse(T, Tinv, d, ok)
                rz = z(i,:) - matmul(R(:,:,i), mu(:,k))
                b(:,k)    = mu(:,k) + matmul(SR, matmul(Tinv, rz))
                Bk(:,:,k) = Sig(:,:,k) - matmul(SR, matmul(Tinv, transpose(SR)))
            end do
            lmax = maxval(logp)
            lse  = lmax + log(sum(exp(logp - lmax)))
            rk   = exp(logp - lse)
            if( present(resp) ) resp(i,:) = rk
            xhat(i,:) = 0.d0
            do k = 1, nk
                xhat(i,:) = xhat(i,:) + rk(k)*b(:,k)
            end do
            xcov(:,:,i) = 0.d0
            do k = 1, nk
                xcov(:,:,i) = xcov(:,:,i) + rk(k)*(Bk(:,:,k) + outer(b(:,k)-xhat(i,:), b(:,k)-xhat(i,:), d))
            end do
        end do
        !$omp end parallel do
    end subroutine xd_posterior

    !> log N(z; R mu, R Sigma R^T + N)
    subroutine component_terms( zi, Ri, Ni, muk, Sigk, d, logp, ok )
        integer,  intent(in)  :: d
        real(dp), intent(in)  :: zi(d), Ri(d,d), Ni(d,d), muk(d), Sigk(d,d)
        real(dp), intent(out) :: logp
        logical,  intent(out) :: ok
        real(dp) :: T(d,d), L(d,d), rz(d), y(d), logdet
        integer  :: q
        T  = matmul(Ri, matmul(Sigk, transpose(Ri))) + Ni
        call cholesky(T, L, d, ok)
        if( .not. ok )then
            logp = -huge(0.d0)/2
            return
        endif
        rz = zi - matmul(Ri, muk)
        call chol_forward(L, rz, y, d)
        logdet = 0.d0
        do q = 1, d
            logdet = logdet + 2.d0*log(L(q,q))
        end do
        logp = -0.5d0*(sum(y*y) + logdet + real(d,dp)*log(2.d0*DPI))
    end subroutine component_terms

    ! ======================================================================= small dense algebra

    pure function outer( u, v, d ) result( M )
        integer,  intent(in) :: d
        real(dp), intent(in) :: u(d), v(d)
        real(dp) :: M(d,d)
        integer  :: q
        do q = 1, d
            M(:,q) = u*v(q)
        end do
    end function outer

    subroutine cholesky( A, L, d, ok )
        integer,  intent(in)  :: d
        real(dp), intent(in)  :: A(d,d)
        real(dp), intent(out) :: L(d,d)
        logical,  intent(out) :: ok
        real(dp) :: s
        integer  :: i, j, k
        L  = 0.d0
        ok = .true.
        do j = 1, d
            s = A(j,j)
            do k = 1, j-1
                s = s - L(j,k)*L(j,k)
            end do
            if( s <= 0.d0 )then
                ok = .false.
                return
            endif
            L(j,j) = sqrt(s)
            do i = j+1, d
                s = A(i,j)
                do k = 1, j-1
                    s = s - L(i,k)*L(j,k)
                end do
                L(i,j) = s/L(j,j)
            end do
        end do
    end subroutine cholesky

    !> y = L^-1 r
    subroutine chol_forward( L, r, y, d )
        integer,  intent(in)  :: d
        real(dp), intent(in)  :: L(d,d), r(d)
        real(dp), intent(out) :: y(d)
        integer :: i, k
        do i = 1, d
            y(i) = r(i)
            do k = 1, i-1
                y(i) = y(i) - L(i,k)*y(k)
            end do
            y(i) = y(i)/L(i,i)
        end do
    end subroutine chol_forward

    !> inverse of an SPD matrix by Cholesky; ok=.false. (and identity-scaled fallback) when not SPD
    subroutine spd_inverse( A, Ainv, d, ok )
        integer,  intent(in)  :: d
        real(dp), intent(in)  :: A(d,d)
        real(dp), intent(out) :: Ainv(d,d)
        logical, optional, intent(out) :: ok
        real(dp) :: L(d,d), Linv(d,d), s
        logical  :: lok
        integer  :: i, j, k
        call cholesky(A, L, d, lok)
        if( .not. lok )then
            Ainv = 0.d0
            do i = 1, d
                Ainv(i,i) = 1.d0/max(A(i,i), DTINY)
            end do
            if( present(ok) ) ok = .false.
            return
        endif
        ! Linv = L^-1 (lower triangular)
        Linv = 0.d0
        do j = 1, d
            Linv(j,j) = 1.d0/L(j,j)
            do i = j+1, d
                s = 0.d0
                do k = j, i-1
                    s = s - L(i,k)*Linv(k,j)
                end do
                Linv(i,j) = s/L(i,i)
            end do
        end do
        Ainv = matmul(transpose(Linv), Linv)
        if( present(ok) ) ok = .true.
    end subroutine spd_inverse

    !> Directions of the deconvolved population that carry signal: S_pop = sum_k pi_k (Sig_k + mu_k mu_k') - m m'
    !! against the mean measurement noise Nbar; whiten by Nbar^-1/2, diagonalise, keep eigenvalues > thresh
    !! (default 1: population variance above noise). W = U_kept' Nbar^-1/2, so W z has unit noise per axis.
    subroutine signal_subspace( mu, Sig, pik, nk, Nz, nptcls, d, thresh, W, k )
        use simple_linalg, only: jacobi, eigsrt
        integer,  intent(in)  :: nk, nptcls, d
        real(dp), intent(in)  :: mu(d,nk), Sig(d,d,nk), pik(nk), Nz(d,d,nptcls)
        real(dp), optional, intent(in) :: thresh
        real(dp), allocatable, intent(out) :: W(:,:)
        integer,  intent(out) :: k
        real(dp) :: S(d,d), Nbar(d,d), Nh(d,d), Mw(d,d), mvec(d), ev(d), evec(d,d), evn(d), evecn(d,d), thr, floor_ev
        integer  :: i, q, kk, nrot
        thr = 1.d0
        if( present(thresh) ) thr = thresh
        mvec = 0.d0
        do kk = 1, nk
            mvec = mvec + pik(kk)*mu(:,kk)
        end do
        S = 0.d0
        do kk = 1, nk
            S = S + pik(kk)*(Sig(:,:,kk) + outer(mu(:,kk), mu(:,kk), d))
        end do
        S = S - outer(mvec, mvec, d)
        S = 0.5d0*(S + transpose(S))
        Nbar = 0.d0
        do i = 1, nptcls
            Nbar = Nbar + Nz(:,:,i)
        end do
        Nbar = Nbar/real(nptcls,dp)
        Nbar = 0.5d0*(Nbar + transpose(Nbar))
        call jacobi(Nbar, d, d, evn, evecn, nrot)
        floor_ev = 1.d-8*max(maxval(evn), DTINY)
        Nh = 0.d0
        do q = 1, d
            Nh = Nh + outer(evecn(:,q), evecn(:,q), d)/sqrt(max(evn(q), floor_ev))
        end do
        Mw = matmul(Nh, matmul(S, Nh))
        Mw = 0.5d0*(Mw + transpose(Mw))
        call jacobi(Mw, d, d, ev, evec, nrot)
        call eigsrt(ev, evec, d, d)
        k = count(ev > thr)
        k = max(1, k)
        allocate(W(k,d))
        do q = 1, k
            W(q,:) = matmul(evec(:,q), Nh)
        end do
        write(logfhandle,'(A,I0,A,I0,A,F5.2,A)') '>>> FLEX_PCA SIGNAL SUBSPACE: ', k, ' of ', d, &
            &' noise-whitened directions carry population variance above ', thr, ' x noise'
        write(logfhandle,'(A,20(1X,F8.2))') '>>>   signal/noise per direction:', (ev(q), q=1,min(d,20))
        call flush(logfhandle)
    end subroutine signal_subspace

    !> clip a symmetric matrix to eigenvalues >= ridge (Jacobi)
    subroutine psd_clip( S, d, ridge )
        use simple_linalg, only: jacobi
        integer,  intent(in)    :: d
        real(dp), intent(inout) :: S(d,d)
        real(dp), intent(in)    :: ridge
        real(dp) :: A(d,d), ev(d), evec(d,d)
        integer  :: nrot, q
        A = S
        call jacobi(A, d, d, ev, evec, nrot)
        do q = 1, d
            ev(q) = max(ev(q), ridge)
        end do
        S = 0.d0
        do q = 1, d
            S = S + ev(q)*outer(evec(:,q), evec(:,q), d)
        end do
        S = 0.5d0*(S + transpose(S))
    end subroutine psd_clip

    ! ======================================================================= self-test

    !> Synthetic check: a 2-component population in d=4 under heavy heteroscedastic noise. The
    !! held-out rule must pick K=2, the means must be recovered, and the posterior means must be
    !! closer to the truth than the observations.
    subroutine test_flex_pca_deconv()
        integer,  parameter :: n = 20000, d = 4
        real(dp), allocatable :: x(:,:), z(:,:), prec(:,:,:), zhalf(:,:,:), z0(:,:)
        real(dp) :: prior(d), a, a_comp(d), s, g(d), mse_z, mse_x, mu_true(d,2)
        integer  :: i, q, k_out, kt
        real(dp) :: u
        allocate(x(n,d), z(n,d), prec(d,d,n), zhalf(n,d,2), z0(n,d))
        mu_true = 0.d0
        mu_true(1,1) = -1.5d0; mu_true(1,2) = 1.5d0
        mu_true(2,1) =  0.5d0; mu_true(2,2) = -0.5d0
        prior = 1.d0/4.d0                              ! prior variance 4 per axis (weak)
        call random_seed()
        do i = 1, n
            kt = merge(1, 2, mod(i,3) == 0)            ! weights 1/3, 2/3
            do q = 1, d
                x(i,q) = mu_true(q,kt) + 0.3d0*gauss()
            end do
            call random_number(u)
            s = 2.d0 + 18.d0*u                         ! noise variance 2..20 per axis (signal ~0.09 + means)
            prec(:,:,i) = 0.d0
            do q = 1, d
                prec(q,q,i)  = 1.d0/s + prior(q)     ! posterior precision = data + prior
                g(q)         = gauss()*sqrt(s)
            end do
            ! z = A^-1 (D x + e), with D = 1/s I:  R = D/A, N = D/A^2
            do q = 1, d
                z(i,q) = ((x(i,q) + g(q))/s)/prec(q,q,i)
                ! halves: each with half the data precision and independent noise
                zhalf(i,q,1) = ((x(i,q) + gauss()*sqrt(2.d0*s))/(2.d0*s))/(1.d0/(2.d0*s) + prior(q))
                zhalf(i,q,2) = ((x(i,q) + gauss()*sqrt(2.d0*s))/(2.d0*s))/(1.d0/(2.d0*s) + prior(q))
            end do
        end do
        z0 = z
        call calibrate_noise_scale(z, zhalf, prec, prior, n, d, a, a_comp)
        if( abs(a - 1.d0) > 0.15d0 ) THROW_HARD('test_flex_pca_deconv: calibration off (expected 1)')
        call deconvolve_latent(z, prec, prior, n, d, a, 4, k_out)
        if( k_out /= 2 ) THROW_HARD('test_flex_pca_deconv: held-out rule did not pick K=2')
        mse_z = sum((z0 - x)**2)/real(n*d,dp)
        mse_x = sum((z  - x)**2)/real(n*d,dp)
        write(logfhandle,'(A,F8.4,A,F8.4)') '>>> test_flex_pca_deconv: mse(z)=', mse_z, '  mse(xhat)=', mse_x
        if( mse_x > 0.6d0*mse_z ) THROW_HARD('test_flex_pca_deconv: posterior means not closer to the truth')
        write(logfhandle,'(A)') '>>> test_flex_pca_deconv PASSED'
        deallocate(x, z, prec, zhalf, z0)

    contains

        real(dp) function gauss()
            real(dp) :: u1, u2
            call random_number(u1); call random_number(u2)
            gauss = sqrt(-2.d0*log(max(u1, 1.d-12)))*cos(2.d0*DPI*u2)
        end function gauss

    end subroutine test_flex_pca_deconv

end module simple_flex_pca_deconv
