!@descr: flex_pca latent deconvolution: calibrated per-particle noise + an empirical-Bayes mixture prior fitted through it
!! MAP latents z_i = A_i^-1(D_i x_i + e_i), D_i = A_i - P: E[z|x] = R_i x, R_i = A_i^-1 D_i, noise a*A_i^-1 D_i A_i^-1,
!! with the scalar a calibrated from even/odd half solutions. The prior is a K-Gaussian mixture fitted by extreme
!! deconvolution (simple_xd_gmm). K is chosen by particle-half held-out log-likelihood on a strided
!! subsample of ~XD_CV_MAX particles; the final fit uses all. z/precision become posterior means/precisions.
module simple_flex_pca_deconv
use simple_core_module_api, only: dp, dtiny, logfhandle, tic, timer_int_kind, toc
use simple_linalg,          only: spd_inverse
use simple_xd_gmm,          only: xd_gmm, xd_select_k
implicit none
private
#include "simple_local_flags.inc"

public :: calibrate_noise_scale, deconvolve_latent

integer,  parameter :: XD_MAXIT  = 150
real(dp), parameter :: XD_TOL    = 1.d-6      !< relative log-likelihood change
integer,  parameter :: XD_CV_MAX = 20000      !< particles in the K-selection ladder (both halves together)
real(dp), parameter :: XD_PI_MIN = 1.d-4
integer,  parameter :: XD_KMAX   = 16

contains

    ! ======================================================================= calibration

    !> a = sum_i |z_i^even - z_i^odd|^2 / sum_i tr( A_h^-1 D_i A_h^-1 ),  A_h = P + D_i/2, D_i = A_i - P
    subroutine calibrate_noise_scale( zhalf, precision, prior, nptcls, ncomp, a, a_comp )
        integer,  intent(in)  :: nptcls, ncomp
        real(dp), intent(in)  :: zhalf(nptcls,ncomp,2)
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
            call spd_inverse_or_diag(Ah, Ahinv, ncomp, ok)
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
        integer,                    intent(in)    :: nptcls, ncomp
        real(dp),                   intent(inout) :: z(nptcls,ncomp)
        real(dp),                   intent(inout) :: precision(ncomp,ncomp,nptcls)
        real(dp),                   intent(in)    :: prior(ncomp), a
        integer,                    intent(in)    :: kmax_in
        integer,                    intent(out)   :: k_out
        character(len=*), optional, intent(in)    :: prior_fname
        character(len=*), optional, intent(in)    :: labels_fname   !< per-particle argmax component + its responsibility
        integer,          optional, intent(in)    :: pinds(:)        !< project rows for that file
        integer, allocatable, optional, intent(out)   :: labels_out(:)  !< argmax component per particle
        real(dp), allocatable :: R(:,:,:), Nz(:,:,:), mu(:,:), Sig(:,:,:), pik(:), xhat(:,:), xcov(:,:,:)
        real(dp), allocatable :: score(:), resp(:,:)
        type(xd_gmm) :: xd
        real(dp) :: ll, cobs(ncomp,ncomp), zmean(ncomp), vpop(ncomp), vobs(ncomp), mbar(ncomp)
        integer  :: kmax, i, q, k, u, stride, nsub
        integer(timer_int_kind) :: t0
        t0 = tic()
        allocate(R(ncomp,ncomp,nptcls), Nz(ncomp,ncomp,nptcls))
        call noise_and_projection(precision, prior, nptcls, ncomp, a, R, Nz)
        ! observed moments (for the report)
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
        ! K is a coarse choice: select it on a strided subsample of ~XD_CV_MAX particles; the final
        ! fit below still sees every particle
        stride = max(1, nint(real(nptcls,dp)/real(XD_CV_MAX,dp)))
        nsub   = count([(mod(i, 2*stride) == 1 .or. mod(i, 2*stride) == mod(stride + 1, 2*stride), i=1,nptcls)])
        write(logfhandle,'(A,I0,A,I0,A)') '>>> FLEX_PCA DECONV K selection on ', nsub, &
            &' of ', nptcls, ' particles (two disjoint halves); final fit on all'
        allocate(score(kmax))
        call xd_select_k(z, R, Nz, kmax, XD_CV_MAX, XD_MAXIT, XD_TOL, XD_PI_MIN, k_out, score)
        ! ---- final fit on all particles ----
        call xd%new(ncomp, k_out, XD_MAXIT, XD_TOL, XD_PI_MIN)
        call xd%init(z, Nz)
        call xd%fit(z, R, Nz)
        allocate(mu(ncomp,k_out), Sig(ncomp,ncomp,k_out), pik(k_out))
        call xd%get_means(mu)
        call xd%get_covs(Sig)
        call xd%get_pi(pik)
        ll = xd%get_loglik()
        allocate(xhat(nptcls,ncomp), xcov(ncomp,ncomp,nptcls), resp(nptcls,k_out))
        call xd%posterior(z, R, Nz, xhat, xcov, resp)
        call xd%kill
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
            call spd_inverse_or_diag(xcov(:,:,i), precision(:,:,i), ncomp)
        end do
        deallocate(R, Nz, mu, Sig, pik, xhat, xcov, score, resp)
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
            call spd_inverse_or_diag(precision(:,:,i), Ainv, ncomp, ok)
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

    ! ======================================================================= small dense algebra

    !> spd_inverse (simple_linalg) with the deconvolution's fallback: the inverse diagonal when A is not SPD
    subroutine spd_inverse_or_diag( A, Ainv, d, ok )
        integer,           intent(in)  :: d
        real(dp),          intent(in)  :: A(d,d)
        real(dp),          intent(out) :: Ainv(d,d)
        logical, optional, intent(out) :: ok
        logical :: lok
        integer :: i
        call spd_inverse(A, Ainv, d, lok)
        if( .not. lok )then
            Ainv = 0.d0
            do i = 1, d
                Ainv(i,i) = 1.d0/max(A(i,i), DTINY)
            end do
        endif
        if( present(ok) ) ok = lok
    end subroutine spd_inverse_or_diag

end module simple_flex_pca_deconv
