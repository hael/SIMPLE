!@descr: extreme deconvolution: a full-covariance Gaussian mixture fitted to noisy, projected observations
! Bovy, Hogg & Roweis (2011). Point i is observed as x_i = R_i v_i + e_i with known projection R_i and
! noise covariance N_i, and v_i drawn from sum_k pi_k N(mu_k, Sigma_k). The model is fitted to the
! underlying v by EM: per point and component one Cholesky factor of T_ik = R_i Sigma_k R_i' + N_i serves
! the responsibility (log-sum-exp, simple_stat) and the conditional moments of v. Start: equal-mass
! quantile means along the axis of largest observed variance (simple_clustering_utils), every Sigma_k
! the observed covariance minus the mean noise, clipped to a ridge. Components below the mixing floor
! freeze. After a fit the components are ordered by decreasing mixing proportion (ties by lower index).
! The data are arguments of every call, not copied: R and N are d x d per point. xd_select_k chooses the
! number of components by held-out log-likelihood on two disjoint strided halves.
module simple_xd_gmm
!$ use omp_lib, only: omp_get_max_threads, omp_get_thread_num
use simple_core_module_api,  only: dp, dpi, dtiny, logfhandle, simple_exception
use simple_linalg,           only: cholesky, chol_forward, jacobi
use simple_stat,             only: logsumexp
use simple_clustering_utils, only: equal_mass_quantile_start
implicit none

public :: xd_gmm, xd_select_k
private
#include "simple_local_flags.inc"

real(dp), parameter :: XD_RIDGE_REL = 1.d-6 !< Sigma_k ridge, relative to the mean variance

type xd_gmm
    private
    integer               :: d      = 0
    integer               :: k      = 0
    integer               :: maxits = 0
    real(dp)              :: tol    = 0._dp   !< relative log-likelihood change of convergence
    real(dp)              :: pi_min = 0._dp   !< mixing floor: a component below it freezes
    real(dp)              :: ll     = 0._dp   !< log-likelihood of the last fit (sum over points)
    real(dp), allocatable :: mu(:,:)          !< (d,k)
    real(dp), allocatable :: Sig(:,:,:)       !< (d,d,k)
    real(dp), allocatable :: pik(:)           !< (k)
    logical               :: exists = .false.
  contains
    procedure :: new
    procedure :: init
    procedure :: fit
    procedure :: loglik
    procedure :: posterior
    procedure :: get_means
    procedure :: get_covs
    procedure :: get_pi
    procedure :: get_loglik
    procedure :: kill
    procedure, private :: order_by_weight
end type xd_gmm

contains

    !>  \brief  is a constructor
    subroutine new( self, d, k, maxits, tol, pi_min )
        class(xd_gmm), intent(inout) :: self
        integer,       intent(in)    :: d, k, maxits
        real(dp),      intent(in)    :: tol, pi_min
        call self%kill
        if( d < 1 .or. k < 1 ) THROW_HARD('d and k must be >= 1; new')
        self%d      = d
        self%k      = k
        self%maxits = max(1, maxits)
        self%tol    = tol
        self%pi_min = pi_min
        allocate(self%mu(d,k), self%Sig(d,d,k), self%pik(k), source=0._dp)
        self%exists = .true.
    end subroutine new

    !>  \brief  start: equal-mass quantile means along the axis of largest observed variance (or the means the
    !!          caller supplies), every Sigma_k the observed covariance minus the mean noise (PSD-clipped),
    !!          uniform mixing proportions
    subroutine init( self, X, Nz, means0 )
        class(xd_gmm),      intent(inout) :: self
        real(dp),           intent(in)    :: X(:,:)      !< (n,d) observations
        real(dp),           intent(in)    :: Nz(:,:,:)   !< (d,d,n) noise covariances
        real(dp), optional, intent(in)    :: means0(:,:) !< (d,k) initial means
        real(dp) :: cobs(self%d,self%d), xmean(self%d), S(self%d,self%d), Nbar(self%d,self%d), ridge
        integer  :: n, i, q, r, jbest
        n = size(X,1)
        if( size(X,2) /= self%d .or. size(Nz,3) /= n ) THROW_HARD('inconsistent data dimensions; init')
        if( n < self%k ) THROW_HARD('fewer points than components; init')
        do q = 1, self%d
            xmean(q) = sum(X(:,q))/real(n,dp)
        end do
        do q = 1, self%d
            do r = 1, self%d
                cobs(q,r) = sum((X(:,q)-xmean(q))*(X(:,r)-xmean(r)))/real(n,dp)
            end do
        end do
        if( present(means0) )then
            if( size(means0,1) /= self%d .or. size(means0,2) /= self%k ) THROW_HARD('means0 must be (d,k); init')
            self%mu = means0
        else
            jbest = maxloc([(cobs(q,q), q=1,self%d)], dim=1)
            call equal_mass_quantile_start(X, X(:,jbest), self%k, self%mu)
        endif
        Nbar = 0._dp
        do i = 1, n
            Nbar = Nbar + Nz(:,:,i)
        end do
        Nbar  = Nbar/real(n,dp)
        S     = cobs - Nbar
        ridge = XD_RIDGE_REL*max(sum([(cobs(q,q), q=1,self%d)])/real(self%d,dp), DTINY)
        call psd_clip(S, self%d, ridge)
        do q = 1, self%k
            self%Sig(:,:,q) = S
            self%pik(q)     = 1._dp/real(self%k,dp)
        end do
    end subroutine init

    !>  \brief  extreme-deconvolution EM from the current model
    subroutine fit( self, X, R, Nz )
        class(xd_gmm), intent(inout) :: self
        real(dp),      intent(in)    :: X(:,:)     !< (n,d) observations
        real(dp),      intent(in)    :: R(:,:,:)   !< (d,d,n) projections
        real(dp),      intent(in)    :: Nz(:,:,:)  !< (d,d,n) noise covariances
        real(dp), allocatable :: sw(:,:), sb(:,:,:), sBB(:,:,:,:)
        real(dp) :: Lk(self%d,self%d,self%k), RSk(self%d,self%d,self%k), yk(self%d,self%k), b(self%d)
        real(dp) :: Bk(self%d,self%d), logp(self%k), lse, rk(self%k), xi(self%d)
        real(dp) :: ll, ll_prev, ridge, trc, wk
        integer  :: n, d, nk, it, i, c, ithr, nthr, q
        logical  :: okk(self%k)
        n  = size(X,1)
        d  = self%d
        nk = self%k
        if( size(X,2) /= d .or. size(R,3) /= n .or. size(Nz,3) /= n ) THROW_HARD('inconsistent data dimensions; fit')
        nthr = 1
        !$ nthr = omp_get_max_threads()
        allocate(sw(nk,nthr), sb(d,nk,nthr), sBB(d,d,nk,nthr))
        trc = 0._dp
        do c = 1, nk
            trc = trc + sum([(self%Sig(q,q,c), q=1,d)])
        end do
        ridge   = XD_RIDGE_REL*max(trc/real(d*nk,dp), DTINY)
        ll      = 0._dp
        ll_prev = -huge(0._dp)
        do it = 1, self%maxits
            sw = 0._dp; sb = 0._dp; sBB = 0._dp
            ll = 0._dp
            !$omp parallel do default(shared) private(i,c,ithr,Lk,RSk,yk,b,Bk,logp,lse,rk,okk,xi) &
            !$omp& schedule(static) reduction(+:ll)
            do i = 1, n
                ithr = 1
                !$ ithr = omp_get_thread_num() + 1
                xi = X(i,:)                ! the row is strided in X(n,d)
                do c = 1, nk
                    call component_factor(xi, R(:,:,i), Nz(:,:,i), self%mu(:,c), self%Sig(:,:,c), d, &
                        &logp(c), Lk(:,:,c), RSk(:,:,c), yk(:,c), okk(c))
                    logp(c) = logp(c) + log(max(self%pik(c), self%pi_min))
                    if( .not. okk(c) ) logp(c) = -huge(0._dp)/2
                end do
                lse = logsumexp(logp, rk)
                ll  = ll + lse
                do c = 1, nk
                    if( rk(c) < 1.d-12 .or. .not. okk(c) ) cycle
                    call component_moments(Lk(:,:,c), RSk(:,:,c), yk(:,c), self%mu(:,c), self%Sig(:,:,c), d, b, Bk)
                    sw(c,ithr)      = sw(c,ithr) + rk(c)
                    sb(:,c,ithr)    = sb(:,c,ithr) + rk(c)*b
                    sBB(:,:,c,ithr) = sBB(:,:,c,ithr) + rk(c)*(Bk + outer(b, b, d))
                end do
            end do
            !$omp end parallel do
            ! M step
            do c = 1, nk
                wk = sum(sw(c,:))
                if( wk <= self%pi_min*real(n,dp) )then
                    self%pik(c) = self%pi_min
                    cycle
                endif
                self%pik(c)     = wk/real(n,dp)
                self%mu(:,c)    = sum(sb(:,c,:), dim=2)/wk
                self%Sig(:,:,c) = sum(sBB(:,:,c,:), dim=3)/wk - outer(self%mu(:,c), self%mu(:,c), d)
                self%Sig(:,:,c) = 0.5_dp*(self%Sig(:,:,c) + transpose(self%Sig(:,:,c)))
                call psd_clip(self%Sig(:,:,c), d, ridge)
            end do
            self%pik = self%pik/sum(self%pik)
            if( abs(ll - ll_prev) <= self%tol*abs(ll) ) exit
            ll_prev = ll
        end do
        self%ll = ll
        deallocate(sw, sb, sBB)
        call self%order_by_weight
    end subroutine fit

    !>  \brief  orders the components by decreasing mixing proportion, ties by lower index (stable)
    subroutine order_by_weight( self )
        class(xd_gmm), intent(inout) :: self
        integer,  allocatable :: order(:)
        real(dp), allocatable :: tmu(:,:), tSig(:,:,:), tpik(:)
        integer :: i, j, tmp
        allocate(order(self%k))
        order = [(i, i=1,self%k)]
        do i = 2, self%k
            j = i
            do while( j > 1 )
                if( self%pik(order(j-1)) >= self%pik(order(j)) ) exit
                tmp        = order(j-1)
                order(j-1) = order(j)
                order(j)   = tmp
                j = j - 1
            end do
        end do
        tmu  = self%mu(:,order)
        tSig = self%Sig(:,:,order)
        tpik = self%pik(order)
        self%mu  = tmu
        self%Sig = tSig
        self%pik = tpik
        deallocate(order, tmu, tSig, tpik)
    end subroutine order_by_weight

    !>  \brief  log-likelihood (sum over points) of an observation set under the model
    function loglik( self, X, R, Nz ) result( ll )
        class(xd_gmm), intent(in) :: self
        real(dp),      intent(in) :: X(:,:), R(:,:,:), Nz(:,:,:)
        real(dp) :: ll, logp(self%k), L(self%d,self%d), RS(self%d,self%d), y(self%d), xi(self%d)
        integer  :: i, c
        logical  :: ok
        ll = 0._dp
        !$omp parallel do default(shared) private(i,c,logp,L,RS,y,ok,xi) schedule(static) reduction(+:ll)
        do i = 1, size(X,1)
            xi = X(i,:)
            do c = 1, self%k
                call component_factor(xi, R(:,:,i), Nz(:,:,i), self%mu(:,c), self%Sig(:,:,c), self%d, logp(c), L, RS, y, ok)
                logp(c) = logp(c) + log(max(self%pik(c), self%pi_min))
                if( .not. ok ) logp(c) = -huge(0._dp)/2
            end do
            ll = ll + logsumexp(logp)
        end do
        !$omp end parallel do
    end function loglik

    !>  \brief  posterior mean and covariance of every underlying point, and its component responsibilities
    subroutine posterior( self, X, R, Nz, xhat, xcov, resp )
        class(xd_gmm),      intent(in)  :: self
        real(dp),           intent(in)  :: X(:,:), R(:,:,:), Nz(:,:,:)
        real(dp),           intent(out) :: xhat(:,:)    !< (n,d)
        real(dp),           intent(out) :: xcov(:,:,:)  !< (d,d,n)
        real(dp), optional, intent(out) :: resp(:,:)    !< (n,k)
        real(dp) :: L(self%d,self%d), RS(self%d,self%d), y(self%d), b(self%d,self%k), Bk(self%d,self%d,self%k)
        real(dp) :: logp(self%k), lse, rk(self%k), xi(self%d)
        integer  :: i, c
        logical  :: ok
        !$omp parallel do default(shared) private(i,c,L,RS,y,b,Bk,logp,lse,rk,ok,xi) schedule(static)
        do i = 1, size(X,1)
            xi = X(i,:)
            do c = 1, self%k
                call component_factor(xi, R(:,:,i), Nz(:,:,i), self%mu(:,c), self%Sig(:,:,c), self%d, logp(c), L, RS, y, ok)
                logp(c) = logp(c) + log(max(self%pik(c), self%pi_min))
                if( .not. ok ) logp(c) = -huge(0._dp)/2
                call component_moments(L, RS, y, self%mu(:,c), self%Sig(:,:,c), self%d, b(:,c), Bk(:,:,c))
            end do
            lse = logsumexp(logp, rk)
            if( present(resp) ) resp(i,:) = rk
            xhat(i,:) = 0._dp
            do c = 1, self%k
                xhat(i,:) = xhat(i,:) + rk(c)*b(:,c)
            end do
            xcov(:,:,i) = 0._dp
            do c = 1, self%k
                xcov(:,:,i) = xcov(:,:,i) + rk(c)*(Bk(:,:,c) + outer(b(:,c)-xhat(i,:), b(:,c)-xhat(i,:), self%d))
            end do
        end do
        !$omp end parallel do
    end subroutine posterior

    subroutine get_means( self, mu )
        class(xd_gmm), intent(in)  :: self
        real(dp),      intent(out) :: mu(:,:) !< (d,k)
        mu = self%mu
    end subroutine get_means

    subroutine get_covs( self, Sig )
        class(xd_gmm), intent(in)  :: self
        real(dp),      intent(out) :: Sig(:,:,:) !< (d,d,k)
        Sig = self%Sig
    end subroutine get_covs

    subroutine get_pi( self, pik )
        class(xd_gmm), intent(in)  :: self
        real(dp),      intent(out) :: pik(:) !< (k)
        pik = self%pik
    end subroutine get_pi

    !>  \brief  log-likelihood (sum over points) of the last fit
    pure real(dp) function get_loglik( self )
        class(xd_gmm), intent(in) :: self
        get_loglik = self%ll
    end function get_loglik

    subroutine kill( self )
        class(xd_gmm), intent(inout) :: self
        if( self%exists )then
            deallocate(self%mu, self%Sig, self%pik)
            self%d = 0; self%k = 0; self%maxits = 0
            self%ll = 0._dp
            self%exists = .false.
        endif
    end subroutine kill

    ! K SELECTION

    !>  \brief  the number of components by held-out log-likelihood: two disjoint strided halves of about
    !!          nsub points together (every point when nsub >= n); for K = 1, 2, ... each half is fitted and
    !!          scored on the other, and the ladder stops after two consecutive decreases or at kmax. kbest
    !!          maximizes the score (mean held-out log-likelihood per point); scores beyond the stop are -huge
    subroutine xd_select_k( X, R, Nz, kmax, nsub, maxits, tol, pi_min, kbest, scores )
        real(dp), intent(in)  :: X(:,:), R(:,:,:), Nz(:,:,:)
        integer,  intent(in)  :: kmax, nsub, maxits
        real(dp), intent(in)  :: tol, pi_min
        integer,  intent(out) :: kbest
        real(dp), intent(out) :: scores(:)   !< (kmax)
        logical, allocatable :: inA(:), inB(:)
        integer :: n, i, kcv, ndrop, stride
        n = size(X,1)
        if( kmax < 1 .or. size(scores) < kmax ) THROW_HARD('kmax must be >= 1 and scores hold kmax entries; xd_select_k')
        allocate(inA(n), inB(n))
        stride = max(1, nint(real(n,dp)/real(max(1,nsub),dp)))
        do i = 1, n
            inA(i) = mod(i, 2*stride) == 1
            inB(i) = mod(i, 2*stride) == mod(stride + 1, 2*stride)
        end do
        scores = -huge(0._dp)
        ndrop  = 0
        do kcv = 1, kmax
            scores(kcv) = (fit_half(kcv, inA, inB) + fit_half(kcv, inB, inA)) / real(count(inA) + count(inB),dp)
            write(logfhandle,'(A,I2,A,F14.6)') '>>> XD K=', kcv, '  held-out loglik/point=', scores(kcv)
            call flush(logfhandle)
            if( kcv >= 2 )then
                if( scores(kcv) < scores(kcv-1) )then
                    ndrop = ndrop + 1
                else
                    ndrop = 0
                endif
                if( ndrop >= 2 ) exit
            endif
        end do
        kbest = maxloc(scores(1:min(kcv,kmax)), dim=1)
        deallocate(inA, inB)

    contains

        !> fit nk components on the points flagged in sel, return the log-likelihood of those in tst
        real(dp) function fit_half( nk, sel, tst ) result( ll_held )
            integer, intent(in) :: nk
            logical, intent(in) :: sel(n), tst(n)
            real(dp), allocatable :: Xsel(:,:), Rsel(:,:,:), Nsel(:,:,:), Xtst(:,:), Rtst(:,:,:), Ntst(:,:,:)
            type(xd_gmm) :: xd
            integer :: n_sel, n_tst, j, jj, ii, d
            d     = size(X,2)
            n_sel = count(sel)
            n_tst = count(tst)
            allocate(Xsel(n_sel,d), Rsel(d,d,n_sel), Nsel(d,d,n_sel), Xtst(n_tst,d), Rtst(d,d,n_tst), Ntst(d,d,n_tst))
            j = 0; jj = 0
            do ii = 1, n
                if( sel(ii) )then
                    j = j + 1; Xsel(j,:) = X(ii,:); Rsel(:,:,j) = R(:,:,ii); Nsel(:,:,j) = Nz(:,:,ii)
                else if( tst(ii) )then
                    jj = jj + 1; Xtst(jj,:) = X(ii,:); Rtst(:,:,jj) = R(:,:,ii); Ntst(:,:,jj) = Nz(:,:,ii)
                endif
            end do
            call xd%new(d, nk, maxits, tol, pi_min)
            call xd%init(Xsel, Nsel)
            call xd%fit(Xsel, Rsel, Nsel)
            ll_held = xd%loglik(Xtst, Rtst, Ntst)
            call xd%kill
            deallocate(Xsel, Rsel, Nsel, Xtst, Rtst, Ntst)
        end function fit_half

    end subroutine xd_select_k

    ! NUMERICS

    !>  \brief  log N(x; R mu, T), T = R Sigma R' + N, with what component_moments needs: the lower Cholesky
    !!          factor L of T, RS = R Sigma and y = L^-1 (x - R mu). When T is not SPD: ok=.false.,
    !!          logp=-huge/2 and L = diag(T)^(1/2), so the moments fall back to diag(T)^-1
    subroutine component_factor( xi, Ri, Ni, muk, Sigk, d, logp, L, RS, y, ok )
        integer,  intent(in)  :: d
        real(dp), intent(in)  :: xi(d), Ri(d,d), Ni(d,d), muk(d), Sigk(d,d)
        real(dp), intent(out) :: logp, L(d,d), RS(d,d), y(d)
        logical,  intent(out) :: ok
        real(dp) :: T(d,d), rx(d), logdet
        integer  :: q
        RS = matmul(Ri, Sigk)
        T  = matmul(RS, transpose(Ri)) + Ni
        call cholesky(T, L, d, ok)
        if( .not. ok )then
            L = 0._dp
            do q = 1, d
                L(q,q) = sqrt(max(T(q,q), DTINY))
            end do
        endif
        rx = xi - matmul(Ri, muk)
        call chol_forward(L, rx, y, d)
        if( .not. ok )then
            logp = -huge(0._dp)/2
            return
        endif
        logdet = 0._dp
        do q = 1, d
            logdet = logdet + 2._dp*log(L(q,q))
        end do
        logp = -0.5_dp*(sum(y*y) + logdet + real(d,dp)*log(2._dp*DPI))
    end subroutine component_factor

    !>  \brief  conditional moments of v given x under one component from component_factor: with
    !!          W = L^-1 R Sigma, bvec = mu + W' y = mu + Sigma R' T^-1 (x - R mu) and
    !!          Bcov = Sigma - W' W = Sigma - Sigma R' T^-1 R Sigma
    subroutine component_moments( L, RS, y, muk, Sigk, d, bvec, Bcov )
        integer,  intent(in)  :: d
        real(dp), intent(in)  :: L(d,d), RS(d,d), y(d), muk(d), Sigk(d,d)
        real(dp), intent(out) :: bvec(d), Bcov(d,d)
        real(dp) :: W(d,d)
        integer  :: q
        do q = 1, d
            call chol_forward(L, RS(:,q), W(:,q), d)
        end do
        bvec = muk + matmul(y, W)
        Bcov = Sigk - matmul(transpose(W), W)
    end subroutine component_moments

    pure function outer( u, v, d ) result( M )
        integer,  intent(in) :: d
        real(dp), intent(in) :: u(d), v(d)
        real(dp) :: M(d,d)
        integer  :: q
        do q = 1, d
            M(:,q) = u*v(q)
        end do
    end function outer

    !>  \brief  clips a symmetric matrix to eigenvalues >= ridge (Jacobi)
    subroutine psd_clip( S, d, ridge )
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
        S = 0._dp
        do q = 1, d
            S = S + ev(q)*outer(evec(:,q), evec(:,q), d)
        end do
        S = 0.5_dp*(S + transpose(S))
    end subroutine psd_clip

end module simple_xd_gmm
