!@descr: Gaussian mixture with a tied (shared) covariance fitted by expectation maximization, deterministic
! Points are the rows of X(n,d); the caller supplies the initial means (e.g. from simple_kmeans or
! equal_mass_quantile_start in simple_clustering_utils). Each iteration: E step by log-sum-exp
! (simple_stat), mixing proportions (optionally floored, a constrained maximum likelihood that keeps a
! component from starving), means, and the tied covariance as the rank-k correction of the fixed second
! moment, sum_k sum_i r_ik (x_i - mu_k)(x_i - mu_k)' = sum_i x_i x_i' - sum_k N_k mu_k mu_k', plus a ridge.
! Converged EM can park two components on one region: with respawn on, a pair closer than one standard
! deviation (squared Mahalanobis distance < 1) keeps its larger member and restarts the other at the
! worst-explained point, at most GMM_MAX_RESPAWN times. Components are labelled by decreasing population
! of the hard (argmax) assignment, ties by smallest member index. Responsibilities are returned floored
! at resp_floor and renormalized (exact zeros below the floor, the best component kept when all fall
! below), the hard labels as their argmax.
module simple_gmm
use simple_core_module_api, only: dp, dtiny, logfhandle, simple_exception
use simple_linalg,          only: jacobi, matinv
use simple_stat,            only: logsumexp
implicit none

public :: gmm
private
#include "simple_local_flags.inc"

integer,  parameter :: GMM_MAX_RESPAWN = 8      !< respawns per fit
real(dp), parameter :: GMM_MERGE_D2    = 1.0_dp !< squared Mahalanobis separation below which two means are one

type gmm
    private
    integer               :: n = 0, d = 0, k = 0
    integer               :: maxits     = 0
    integer               :: niters     = 0       !< EM iterations of the last fit
    integer               :: nrespawn   = 0       !< respawns of the last fit
    real(dp)              :: reg        = 0._dp   !< ridge added to the tied covariance
    real(dp)              :: tol        = 0._dp   !< relative log-likelihood change of convergence
    real(dp)              :: pi_floor   = 0._dp   !< mixing-proportion floor (0: none)
    real(dp)              :: resp_floor = 0._dp   !< responsibility floor of get_resp (0: none)
    real(dp)              :: ll         = 0._dp   !< mean log-likelihood per point of the last fit
    logical               :: l_respawn  = .true.
    logical               :: l_fitted   = .false.
    real(dp), allocatable :: X(:,:)               !< (n,d) points
    real(dp), allocatable :: mu(:,:)              !< (d,k) component means
    real(dp), allocatable :: S(:,:)               !< (d,d) tied covariance
    real(dp), allocatable :: pival(:)             !< (k) mixing proportions
    real(dp), allocatable :: resp(:,:)            !< (n,k) responsibilities
    real(dp), allocatable :: nresp(:)             !< (k) responsibility mass
    logical               :: exists = .false.
  contains
    procedure :: new
    procedure :: fit
    procedure :: get_resp
    procedure :: get_labels
    procedure :: get_means
    procedure :: get_cov
    procedure :: get_pi
    procedure :: get_mass
    procedure :: get_loglik
    procedure :: get_niters
    procedure :: get_bic
    procedure :: get_icl
    procedure :: get_pairsep
    procedure :: kill
    procedure, private :: order_by_population
end type gmm

contains

    !>  \brief  is a constructor
    subroutine new( self, X, k, means0, reg, tol, maxits, pi_floor, resp_floor, respawn )
        class(gmm),         intent(inout) :: self
        real(dp),           intent(in)    :: X(:,:)       !< (n,d) points
        integer,            intent(in)    :: k            !< nr of components
        real(dp),           intent(in)    :: means0(:,:)  !< (d,k) initial means
        real(dp),           intent(in)    :: reg          !< ridge of the tied covariance
        real(dp),           intent(in)    :: tol          !< relative log-likelihood change of convergence
        integer,            intent(in)    :: maxits       !< EM iteration cap
        real(dp), optional, intent(in)    :: pi_floor     !< mixing-proportion floor (default none)
        real(dp), optional, intent(in)    :: resp_floor   !< responsibility floor of get_resp (default none)
        logical,  optional, intent(in)    :: respawn      !< respawn redundant components (default yes)
        call self%kill
        self%n = size(X,1)
        self%d = size(X,2)
        if( self%n < 1 .or. self%d < 1 ) THROW_HARD('empty point set; new')
        if( k < 1 ) THROW_HARD('k must be >= 1; new')
        if( size(means0,1) /= self%d .or. size(means0,2) /= k ) THROW_HARD('means0 must be (d,k); new')
        self%k      = k
        self%reg    = reg
        self%tol    = tol
        self%maxits = max(1, maxits)
        self%pi_floor   = 0._dp
        self%resp_floor = 0._dp
        self%l_respawn  = .true.
        if( present(pi_floor)   ) self%pi_floor   = pi_floor
        if( present(resp_floor) ) self%resp_floor = resp_floor
        if( present(respawn)    ) self%l_respawn  = respawn
        allocate(self%X(self%n,self%d), source=X)
        allocate(self%mu(self%d,k), source=means0)
        allocate(self%S(self%d,self%d), self%pival(k), self%resp(self%n,k), self%nresp(k), source=0._dp)
        self%exists = .true.
    end subroutine new

    !>  \brief  expectation maximization from the initial means; ok is .false. when the tied covariance
    !!          turns singular, and the object then holds no fit
    subroutine fit( self, ok )
        class(gmm), intent(inout) :: self
        logical,    intent(out)   :: ok
        real(dp), allocatable :: Sinv(:,:), Sxx(:,:), Smu(:,:), mSm(:), xbar(:), evwork(:,:), ev(:), evec(:,:)
        real(dp) :: lp(self%k), xSx, xSm, lse, ll, prev_ll, logdet, dmin, d2pair
        integer  :: i, q, r, c, it, nrot, errflg, kmin, kdrop, iworst
        if( .not. self%exists ) THROW_HARD('object does not exist; fit')
        ok = .false.
        self%l_fitted = .false.
        allocate(Sinv(self%d,self%d), Sxx(self%d,self%d), Smu(self%d,self%k), mSm(self%k), xbar(self%d), &
            &evwork(self%d,self%d), ev(self%d), evec(self%d,self%d))
        ! the second moment is fixed, so the tied-covariance M step is a rank-k correction of it
        Sxx = matmul(transpose(self%X), self%X)
        self%S = Sxx / real(self%n,dp)
        do q = 1, self%d
            xbar(q) = sum(self%X(:,q)) / real(self%n,dp)
        end do
        do q = 1, self%d
            do r = 1, self%d
                self%S(q,r) = self%S(q,r) - xbar(q)*xbar(r)
            end do
        end do
        do q = 1, self%d
            self%S(q,q) = self%S(q,q) + self%reg
        end do
        self%pival    = 1._dp / real(self%k,dp)
        prev_ll       = -huge(1._dp)
        ll            = 0._dp
        self%nrespawn = 0
        do it = 1, self%maxits
            call matinv(self%S, Sinv, self%d, errflg)
            if( errflg /= 0 )then
                deallocate(Sinv, Sxx, Smu, mSm, xbar, evwork, ev, evec)
                return
            endif
            evwork = self%S
            call jacobi(evwork, self%d, self%d, ev, evec, nrot)
            logdet = 0._dp
            do q = 1, self%d
                logdet = logdet + log(max(ev(q), DTINY))
            end do
            Smu = matmul(Sinv, self%mu)
            do c = 1, self%k
                mSm(c) = sum(self%mu(:,c)*Smu(:,c))
            end do
            ! E step
            ll = 0._dp
            !$omp parallel do default(shared) private(i,q,r,c,xSx,xSm,lp,lse) schedule(static) reduction(+:ll)
            do i = 1, self%n
                xSx = 0._dp
                do q = 1, self%d
                    do r = 1, self%d
                        xSx = xSx + self%X(i,q)*Sinv(q,r)*self%X(i,r)
                    end do
                end do
                do c = 1, self%k
                    xSm = 0._dp
                    do q = 1, self%d
                        xSm = xSm + self%X(i,q)*Smu(q,c)
                    end do
                    lp(c) = -0.5_dp*(xSx - 2._dp*xSm + mSm(c)) - 0.5_dp*logdet + log(max(self%pival(c), DTINY))
                end do
                lse = logsumexp(lp, self%resp(i,:))
                ll  = ll + lse
            end do
            !$omp end parallel do
            ll = ll / real(self%n,dp)
            ! M step
            do c = 1, self%k
                self%nresp(c) = sum(self%resp(:,c))
            end do
            self%pival = max(self%nresp, DTINY) / real(self%n,dp)
            if( self%pi_floor > 0._dp )then
                self%pival = max(self%pival, self%pi_floor)
                self%pival = self%pival / sum(self%pival)
            endif
            do c = 1, self%k
                do q = 1, self%d
                    self%mu(q,c) = sum(self%resp(:,c)*self%X(:,q)) / max(self%nresp(c), DTINY)
                end do
            end do
            self%S = Sxx
            do c = 1, self%k
                do q = 1, self%d
                    do r = 1, self%d
                        self%S(q,r) = self%S(q,r) - self%nresp(c)*self%mu(q,c)*self%mu(r,c)
                    end do
                end do
            end do
            self%S = self%S / real(self%n,dp)
            do q = 1, self%d
                self%S(q,q) = self%S(q,q) + self%reg
            end do
            self%S = 0.5_dp*(self%S + transpose(self%S))
            if( abs(ll - prev_ll) < self%tol*abs(ll) )then
                if( self%l_respawn .and. self%nrespawn < GMM_MAX_RESPAWN )then
                    ! the closest pair in the metric of this iteration's precision; the less populated respawns
                    kmin = 0; kdrop = 0; dmin = huge(1._dp)
                    do c = 1, self%k - 1
                        do r = c + 1, self%k
                            d2pair = 0._dp
                            do q = 1, self%d
                                d2pair = d2pair + (self%mu(q,c) - self%mu(q,r))*(Smu(q,c) - Smu(q,r))
                            end do
                            if( d2pair < dmin )then
                                dmin  = d2pair
                                kmin  = c
                                kdrop = merge(r, c, self%nresp(r) < self%nresp(c))
                                if( kdrop == c ) kmin = r
                            endif
                        end do
                    end do
                    if( kdrop >= 1 .and. dmin < GMM_MERGE_D2 )then
                        iworst = maxloc(-maxval(self%resp, dim=2), dim=1)
                        write(logfhandle,'(A,I0,A,I0,A,F8.4,A,I0)') '>>> GMM components ', kmin, ' and ', kdrop, &
                            &' are redundant (separation ', real(dmin), '); respawning ', kdrop
                        call flush(logfhandle)
                        self%mu(:,kdrop)  = self%X(iworst,:)
                        self%pival(kdrop) = 1._dp / real(self%k,dp)
                        self%nrespawn     = self%nrespawn + 1
                        prev_ll           = -huge(1._dp)
                        cycle
                    endif
                endif
                exit
            endif
            prev_ll = ll
        end do
        self%niters = min(it, self%maxits)
        self%ll     = ll
        deallocate(Sinv, Sxx, Smu, mSm, xbar, evwork, ev, evec)
        call self%order_by_population
        self%l_fitted = .true.
        ok = .true.
    end subroutine fit

    !>  \brief  relabels the components by decreasing hard-assignment population, ties by smallest member
    subroutine order_by_population( self )
        class(gmm), intent(inout) :: self
        integer,  allocatable :: cnt(:), first_member(:), order(:)
        real(dp), allocatable :: tmp2(:,:), tmp1(:)
        integer :: i, j, c, tmp
        allocate(cnt(self%k), source=0)
        allocate(first_member(self%k), source=huge(1))
        do i = self%n, 1, -1
            c = maxloc(self%resp(i,:), dim=1)
            cnt(c) = cnt(c) + 1
            first_member(c) = i
        end do
        allocate(order(self%k))
        order = [(c, c=1,self%k)]
        do i = 2, self%k
            j = i
            do while( j > 1 )
                if( cnt(order(j-1)) > cnt(order(j)) ) exit
                if( cnt(order(j-1)) == cnt(order(j)) .and. first_member(order(j-1)) < first_member(order(j)) ) exit
                tmp        = order(j-1)
                order(j-1) = order(j)
                order(j)   = tmp
                j = j - 1
            end do
        end do
        tmp2 = self%mu(:,order)
        self%mu = tmp2
        tmp1 = self%pival(order)
        self%pival = tmp1
        tmp1 = self%nresp(order)
        self%nresp = tmp1
        tmp2 = self%resp(:,order)
        self%resp = tmp2
        deallocate(cnt, first_member, order, tmp2, tmp1)
    end subroutine order_by_population

    !>  \brief  responsibilities floored at resp_floor and renormalized: entries below the floor become exact
    !!          zeros, and a point whose every entry falls below keeps its best component with weight one
    subroutine get_resp( self, resp )
        class(gmm), intent(in)  :: self
        real(dp),   intent(out) :: resp(:,:) !< (n,k)
        real(dp) :: lsum
        integer  :: i, c
        if( .not. self%l_fitted ) THROW_HARD('no fit; get_resp')
        resp = self%resp
        if( self%resp_floor <= 0._dp ) return
        !$omp parallel do default(shared) private(i,c,lsum) schedule(static)
        do i = 1, self%n
            lsum = 0._dp
            do c = 1, self%k
                if( resp(i,c) < self%resp_floor )then
                    resp(i,c) = 0._dp
                else
                    lsum = lsum + resp(i,c)
                endif
            end do
            if( lsum > DTINY )then
                resp(i,:) = resp(i,:) / lsum
            else
                c         = maxloc(self%resp(i,:), dim=1)
                resp(i,:) = 0._dp
                resp(i,c) = 1._dp
            endif
        end do
        !$omp end parallel do
    end subroutine get_resp

    !>  \brief  hard labels, the argmax of the responsibilities
    subroutine get_labels( self, labels )
        class(gmm), intent(in)  :: self
        integer,    intent(out) :: labels(:) !< (n)
        integer :: i
        if( .not. self%l_fitted ) THROW_HARD('no fit; get_labels')
        do i = 1, self%n
            labels(i) = maxloc(self%resp(i,:), dim=1)
        end do
    end subroutine get_labels

    subroutine get_means( self, means )
        class(gmm), intent(in)  :: self
        real(dp),   intent(out) :: means(:,:) !< (d,k)
        means = self%mu
    end subroutine get_means

    subroutine get_cov( self, cov )
        class(gmm), intent(in)  :: self
        real(dp),   intent(out) :: cov(:,:) !< (d,d)
        cov = self%S
    end subroutine get_cov

    subroutine get_pi( self, pival )
        class(gmm), intent(in)  :: self
        real(dp),   intent(out) :: pival(:) !< (k)
        pival = self%pival
    end subroutine get_pi

    !>  \brief  responsibility mass of every component (its effective population)
    subroutine get_mass( self, mass )
        class(gmm), intent(in)  :: self
        real(dp),   intent(out) :: mass(:) !< (k)
        mass = self%nresp
    end subroutine get_mass

    !>  \brief  mean log-likelihood per point of the last fit
    pure real(dp) function get_loglik( self )
        class(gmm), intent(in) :: self
        get_loglik = self%ll
    end function get_loglik

    !>  \brief  EM iterations of the last fit
    pure integer function get_niters( self )
        class(gmm), intent(in) :: self
        get_niters = self%niters
    end function get_niters

    !>  \brief  Bayesian information criterion of the tied model, k*d + d(d+1)/2 + k - 1 free parameters
    pure real(dp) function get_bic( self )
        class(gmm), intent(in) :: self
        integer :: nfree
        nfree    = self%k*self%d + (self%d*(self%d+1))/2 + self%k - 1
        get_bic  = -2._dp*self%ll*real(self%n,dp) + real(nfree,dp)*log(real(self%n,dp))
    end function get_bic

    !>  \brief  integrated completed likelihood: the BIC plus twice the assignment entropy, which prefers
    !!          separated components where the BIC spends components on one populated region
    real(dp) function get_icl( self )
        class(gmm), intent(in) :: self
        real(dp) :: ent
        integer  :: i, c
        ent = 0._dp
        !$omp parallel do default(shared) private(i,c) schedule(static) reduction(+:ent)
        do i = 1, self%n
            do c = 1, self%k
                if( self%resp(i,c) > DTINY ) ent = ent - self%resp(i,c)*log(self%resp(i,c))
            end do
        end do
        !$omp end parallel do
        get_icl = self%get_bic() + 2._dp*ent
    end function get_icl

    !>  \brief  pairwise Mahalanobis separation of the fitted means under the tied covariance; ok is
    !!          .false. (and sep left at -1) when the covariance is singular
    subroutine get_pairsep( self, sep, ok )
        class(gmm), intent(in)  :: self
        real(dp),   intent(out) :: sep(:,:) !< (k,k)
        logical,    intent(out) :: ok
        real(dp), allocatable :: Sinv(:,:), dev(:)
        real(dp) :: d2pair
        integer  :: c, r, q, errflg
        sep = -1._dp
        ok  = .false.
        if( .not. self%l_fitted ) return
        allocate(Sinv(self%d,self%d), dev(self%d))
        call matinv(self%S, Sinv, self%d, errflg)
        if( errflg == 0 )then
            sep = 0._dp
            do c = 1, self%k - 1
                do r = c + 1, self%k
                    dev    = self%mu(:,c) - self%mu(:,r)
                    d2pair = 0._dp
                    do q = 1, self%d
                        d2pair = d2pair + dev(q)*sum(Sinv(q,:)*dev(:))
                    end do
                    sep(c,r) = sqrt(max(d2pair, 0._dp))
                    sep(r,c) = sep(c,r)
                end do
            end do
            ok = .true.
        endif
        deallocate(Sinv, dev)
    end subroutine get_pairsep

    subroutine kill( self )
        class(gmm), intent(inout) :: self
        if( self%exists )then
            deallocate(self%X, self%mu, self%S, self%pival, self%resp, self%nresp)
            self%n = 0; self%d = 0; self%k = 0; self%maxits = 0; self%niters = 0; self%nrespawn = 0
            self%ll = 0._dp
            self%l_fitted = .false.
            self%exists   = .false.
        endif
    end subroutine kill

end module simple_gmm
