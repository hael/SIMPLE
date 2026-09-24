!@descr: clustering of a similartity matrix using affinity propagation
module simple_aff_prop
use simple_core_module_api
implicit none

public :: aff_prop
private
#include "simple_local_flags.inc"

integer, parameter :: EXEMPLAR_STABLE_CHECKS = 10

type aff_prop
    private
    integer              :: N                    !< nr of data entries
    real,    allocatable :: A(:,:), R(:,:)       !< A is affinities & R responsibilities
    real,    allocatable :: Aold(:,:), Rold(:,:)
    real,    allocatable :: Rp(:,:), AS(:,:)
    real,    allocatable :: S(:,:)               !< private similarity matrix
    integer              :: maxits=500           !< maximum number of iterations
    real                 :: ftol=1e-9            !< fractional convergence tolerance score
    real                 :: lam=0.5              !< dampening factor
    real                 :: Smin, Smax           !< similarity min/max
    logical              :: exists=.false.       !< to indicate existence
  contains
    procedure :: new
    procedure :: propagate
    procedure :: kill
end type aff_prop

contains

    ! CONSTRUCTOR

    !>  \brief  is a constructor
    subroutine new( self, N, S, ftol, lam, pref, maxits )
        class(aff_prop),   intent(inout) :: self
        integer,           intent(in)    :: N               ! # data entries
        real,              intent(in)    :: S(N,N)          ! similarity matrix
        real,    optional, intent(in)    :: ftol, lam, pref ! tolerance, dampening factor, preference thres
        integer, optional, intent(in)    :: maxits
        integer :: i, j
        real    :: ppref, diff
        call self%kill
        ! set constants
        self%N = N
        if( present(ftol)    ) self%ftol    = ftol
        if( present(lam)     ) self%lam     = lam
        if( present(maxits)  ) self%maxits  = maxits
        ! calculate similarity characteristics
        allocate(self%S(N,N), source=S)
        call analyze_smat(self%S, .false., self%Smin, self%Smax)
        ! low pref [0.1,1] leads to small numbers of clusters and
        ! high pref leads to large numbers of clusters. Frey suggests Smin - (Smax - Smin) as a threshold for fewer cluster
        ppref = self%Smin
        if( present(pref) ) ppref = pref
        ! remove degeneracies reproducibly; volcluster/AP restarts must be deterministic
        do i=1,N
            self%S(i,i) = ppref
        end do
        diff = (self%Smax - self%Smin) * self%ftol
        do i=1,self%N
            do j=1,self%N
                self%S(i,j) = self%S(i,j) + deterministic_tiebreak(i, j) * diff
            end do
        end do
        ! allocate
        allocate( self%A(N,N), self%R(N,N), self%Aold(N,N), self%Rp(N,N),&
        self%Rold(N,N), self%AS(N,N) )
        self%exists = .true.
    end subroutine new

    ! PROPAGATOR

    !>  \brief  is the message passing algorithm
    subroutine propagate( self, centers, labels, simsum )
        class(aff_prop),      intent(inout) :: self
        integer, allocatable, intent(inout) :: centers(:) !< cluster centres
        integer, allocatable, intent(inout) :: labels(:)  !< cluster labels
        real,                 intent(inout) :: simsum     !< similarity sum
        real, allocatable :: similarities(:)
        logical, allocatable :: exemplars(:), exemplars_prev(:)
        real              :: best, second, avail, colsum, val, pseudo, exemplar_tol, realmax
        integer           :: i, j, k, ncls, nborder, fallback_center, nstable, best_idx
        logical           :: converged, use_fallback_center
        ! initialize
        realmax   = huge(realmax)
        converged = .false.
        nstable   = 0
        allocate(exemplars(self%N), exemplars_prev(self%N), source=.false.)
        self%A    = 0.
        self%R    = 0.
        self%Aold = 0.
        self%Rp   = 0.
        self%Rold = 0.
        self%AS   = 0.
        ! iterate
        do i=1,self%maxits
            !------------------------
            ! RESPONSIBILITIES
            !------------------------
            !$omp parallel do default(shared) private(j,k) schedule(static)
            do k = 1, self%N
                do j = 1, self%N
                    self%Rold(j,k) = self%R(j,k)
                    self%AS(j,k)   = self%A(j,k) + self%S(j,k)
                end do
            end do
            !$omp end parallel do
            !$omp parallel do default(shared) private(j,k,best,second,best_idx,val) schedule(static)
            do j = 1, self%N
                best_idx   = 1
                best       = self%AS(j,1)
                second     = -huge(second)
                do k = 2, self%N
                    val = self%AS(j,k)
                    if( val > best )then
                        second     = best
                        best       = val
                        best_idx   = k
                    else if( val > second )then
                        second     = val
                    endif
                end do
                do k = 1, self%N
                    if( k == best_idx )then
                        self%R(j,k) = self%S(j,k) - second
                    else
                        self%R(j,k) = self%S(j,k) - best
                    endif
                end do
            end do
            !$omp end parallel do
            ! update responsibilities (in a dampened fashion)
            !$omp parallel do default(shared) private(j,k) schedule(static)
            do k = 1, self%N
                do j = 1, self%N
                    self%R(j,k) = (1. - self%lam) * self%R(j,k) + self%lam * self%Rold(j,k)
                end do
            end do
            !$omp end parallel do
            !------------------------
            ! AVAILABILITIES
            !------------------------
            !$omp parallel do default(shared) private(j,k) schedule(static)
            do k = 1, self%N
                do j = 1, self%N
                    self%Aold(j,k) = self%A(j,k)
                    self%Rp(j,k)   = max(0., self%R(j,k))
                end do
            end do
            !$omp end parallel do
            !$omp parallel do default(shared) private(j,k,colsum,avail) schedule(static)
            do k = 1, self%N
                colsum = 0.
                do j = 1, self%N
                    if( j /= k ) colsum = colsum + self%Rp(j,k)
                end do
                do j = 1, self%N
                    if( j == k )then
                        self%A(j,k) = colsum
                    else
                        avail = self%R(k,k) + colsum
                        if( self%R(j,k) > 0. ) avail = avail - self%R(j,k)
                        self%A(j,k) = min(0., avail)
                    endif
                end do
            end do
            !$omp end parallel do
            ! update availabilities (in a dampened fashion)
            !$omp parallel do default(shared) private(j,k) schedule(static)
            do k = 1, self%N
                do j = 1, self%N
                    self%A(j,k) = (1. - self%lam) * self%A(j,k) + self%lam * self%Aold(j,k)
                end do
            end do
            !$omp end parallel do
            !========================
            ! Convergence check each 5th iter
            !========================
            if (mod(i,5) == 0) then
                call calc_exemplar_mask(exemplars, exemplar_tol)
                if( any(exemplars) .and. all(exemplars .eqv. exemplars_prev) )then
                    nstable = nstable + 1
                else
                    nstable = 0
                    exemplars_prev = exemplars
                endif
                if( nstable >= EXEMPLAR_STABLE_CHECKS )then
                    write(logfhandle,'(a,i6,1x,a,i0)') 'aff_prop exemplar set stable at iter=', i, ' checks=', nstable
                    converged = .true.
                    exit
                end if
            endif
        end do
        if( .not. converged ) write(logfhandle,'(a,i6)') 'aff_prop WARNING: reached maxits without convergence, iter=', self%maxits
        !$omp parallel do default(shared) private(j,k) schedule(static)
        do k = 1, self%N
            do j = 1, self%N
                self%AS(j,k) = self%A(j,k) + self%R(j,k) ! pseudomarginals
            end do
        end do
        !$omp end parallel do
        call calc_exemplar_mask(exemplars, exemplar_tol)
        nborder = 0
        fallback_center = 1
        ! count the number of clusters
        ncls = 0
        do j=1,self%N
            if( self%AS(j,j) > self%AS(fallback_center,fallback_center) ) fallback_center = j
            if( exemplars(j) )then
                ncls = ncls + 1
            else if( abs(self%AS(j,j)) <= exemplar_tol )then
                nborder = nborder + 1
            endif
        end do
        use_fallback_center = ncls == 0
        if( use_fallback_center ) ncls = 1
        if( nborder > 0 ) write(logfhandle,'(a,i0,1x,es12.4)') 'aff_prop near-zero exemplar candidates ignored: ', nborder, exemplar_tol
        if( allocated(centers) ) deallocate(centers)
        if( allocated(labels) )  deallocate(labels)
        allocate( centers(ncls), similarities(ncls), labels(self%N) )
        ! set the cluster centers
        if( use_fallback_center )then
            centers(1) = fallback_center
        else
            ncls = 0
            do j=1,self%N
                if( exemplars(j) )then
                    ncls = ncls + 1
                    centers(ncls) = j
                endif
            end do
        endif
        ! report back the labeling
        simsum = 0.
        do j=1,self%N
            do k=1,ncls
                if( j .ne. centers(k) )then
                    similarities(k) = self%S(centers(k),j)
                else
                    similarities(k) = realmax
                endif
            end do
            labels(j) = 1
            do k=2,ncls
                if( similarities(k) > similarities(labels(j)) ) labels(j) = k
            end do
            if( j .ne. centers(labels(j)) ) simsum = simsum + similarities(labels(j))
        end do
        simsum = simsum / real(self%N)
        deallocate(similarities, exemplars, exemplars_prev)

        contains

            subroutine calc_exemplar_mask( mask, tol )
                logical, intent(inout) :: mask(:)
                real,    intent(out)   :: tol
                real :: pseudo_absmax_here
                integer :: ii
                pseudo_absmax_here = 0.
                do ii=1,self%N
                    pseudo = self%A(ii,ii) + self%R(ii,ii)
                    pseudo_absmax_here = max(pseudo_absmax_here, abs(pseudo))
                end do
                tol = max(10. * epsilon(tol), self%ftol) * max(1., pseudo_absmax_here)
                do ii=1,self%N
                    mask(ii) = (self%A(ii,ii) + self%R(ii,ii)) > tol
                end do
            end subroutine calc_exemplar_mask
    end subroutine propagate

    pure real function deterministic_tiebreak( i, j ) result( val )
        integer, intent(in) :: i, j
        integer, parameter :: HASH_MOD = 8191, HASH_A = 37, HASH_B = 1009
        integer :: ii, jj, h
        ii  = mod(i - 1, HASH_MOD)
        jj  = mod(j - 1, HASH_MOD)
        h   = mod(HASH_A * ii + HASH_B * jj + mod(ii * jj, HASH_MOD), HASH_MOD)
        val = real(h) / real(HASH_MOD)
    end function deterministic_tiebreak

    ! DESTRUCTOR

    subroutine kill( self )
        class(aff_prop), intent(inout) :: self
        if( self%exists )then
            deallocate( self%S, self%A, self%R, self%Aold, self%Rp,&
            self%Rold, self%AS )
            self%exists = .false.
        endif
    end subroutine kill

end module simple_aff_prop
