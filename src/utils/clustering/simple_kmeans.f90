!@descr: k-means clustering of feature vectors with deterministic farthest-point seeding and empty-cluster recovery
! Points are the rows of X(n,d); the squared distance takes optional per-dimension weights. Seeding: the
! point nearest the mean, then farthest-point (simple_clustering_utils), so no random draw is involved.
! Lloyd iterations run to convergence (no reassignment) or maxits; an empty cluster is reseeded at the
! worst-fitted point. Labels by decreasing population (ties by smallest member index).
module simple_kmeans
use simple_core_module_api,  only: dp, simple_exception
use simple_clustering_utils, only: farthest_point_seeds
implicit none

public :: kmeans
private
#include "simple_local_flags.inc"

integer, parameter :: KMEANS_MAXITS = 100

type kmeans
    private
    integer               :: n        = 0       !< nr of points
    integer               :: d        = 0       !< nr of dimensions
    integer               :: k        = 0       !< nr of clusters
    integer               :: maxits   = 0       !< Lloyd iteration cap
    integer               :: niters   = 0       !< Lloyd iterations of the last cluster call
    integer               :: nchanged = 0       !< reassignments on the last iteration of the last cluster call
    real(dp), allocatable :: X(:,:)             !< (n,d) points
    real(dp), allocatable :: wdim(:)            !< (d) per-dimension weights of the squared distance
    logical               :: exists = .false.
  contains
    procedure :: new
    procedure :: cluster
    procedure :: get_niters
    procedure :: get_nchanged
    procedure :: kill
end type kmeans

contains

    !>  \brief  is a constructor
    subroutine new( self, X, k, wdim, maxits )
        class(kmeans),      intent(inout) :: self
        real(dp),           intent(in)    :: X(:,:)  !< (n,d) points
        integer,            intent(in)    :: k       !< nr of clusters
        real(dp), optional, intent(in)    :: wdim(:) !< (d) per-dimension weights (default 1)
        integer,  optional, intent(in)    :: maxits  !< Lloyd iteration cap (default KMEANS_MAXITS)
        call self%kill
        self%n = size(X,1)
        self%d = size(X,2)
        if( self%n < 1 .or. self%d < 1 ) THROW_HARD('empty point set; new')
        if( k < 1 .or. k > self%n )      THROW_HARD('k must be in [1,n]; new')
        self%k      = k
        self%maxits = KMEANS_MAXITS
        if( present(maxits) ) self%maxits = max(1, maxits)
        allocate(self%X(self%n,self%d), source=X)
        allocate(self%wdim(self%d), source=1._dp)
        if( present(wdim) )then
            if( size(wdim) /= self%d ) THROW_HARD('wdim must have one weight per dimension; new')
            self%wdim = wdim
        endif
        self%exists = .true.
    end subroutine new

    !>  \brief  seeds, iterates and labels; centroids(:,c) is the mean of the points labelled c
    subroutine cluster( self, labels, centroids )
        class(kmeans), intent(inout) :: self
        integer,       intent(out)   :: labels(:)      !< (n) cluster labels, 1 the most populated
        real(dp),      intent(out)   :: centroids(:,:) !< (d,k)
        real(dp), allocatable :: mind(:), csum(:,:), xbar(:), cen(:,:)
        integer,  allocatable :: cnt(:), memb(:), seeds(:), first_member(:), order(:), newlab(:)
        real(dp) :: d2, best
        integer  :: i, q, s, it, ibest, iseed, nchanged, j, tmp
        logical  :: l_reseed
        if( .not. self%exists ) THROW_HARD('object does not exist; cluster')
        if( size(labels) /= self%n ) THROW_HARD('labels must have one entry per point; cluster')
        if( size(centroids,1) /= self%d .or. size(centroids,2) /= self%k ) THROW_HARD('centroids must be (d,k); cluster')
        allocate(mind(self%n), csum(self%d,self%k), xbar(self%d), cen(self%d,self%k), cnt(self%k), memb(self%n), seeds(self%k))
        ! first seed: the point nearest the mean in the weighted metric
        do q = 1, self%d
            xbar(q) = sum(self%X(:,q)) / real(self%n,dp)
        end do
        best = huge(1._dp); iseed = 1
        do i = 1, self%n
            d2 = 0._dp
            do q = 1, self%d
                d2 = d2 + self%wdim(q)*(self%X(i,q) - xbar(q))**2
            end do
            if( d2 < best )then
                best  = d2
                iseed = i
            endif
        end do
        call farthest_point_seeds(self%X, self%k, iseed, seeds, wdim=self%wdim)
        do s = 1, self%k
            cen(:,s) = self%X(seeds(s),:)
        end do
        ! Lloyd iterations
        memb     = 0
        nchanged = 0
        do it = 1, self%maxits
            nchanged = 0
            !$omp parallel do default(shared) private(i,q,s,d2,best,ibest) schedule(static) reduction(+:nchanged)
            do i = 1, self%n
                best = huge(1._dp); ibest = 1
                do s = 1, self%k
                    d2 = 0._dp
                    do q = 1, self%d
                        d2 = d2 + self%wdim(q)*(self%X(i,q) - cen(q,s))**2
                    end do
                    if( d2 < best )then
                        best  = d2
                        ibest = s
                    endif
                end do
                if( memb(i) /= ibest ) nchanged = nchanged + 1
                memb(i) = ibest
                mind(i) = best
            end do
            !$omp end parallel do
            csum = 0._dp; cnt = 0
            do i = 1, self%n
                cnt(memb(i))    = cnt(memb(i)) + 1
                csum(:,memb(i)) = csum(:,memb(i)) + self%X(i,:)
            end do
            l_reseed = .false.
            do s = 1, self%k
                if( cnt(s) > 0 )then
                    cen(:,s) = csum(:,s) / real(cnt(s),dp)
                else
                    ! empty cluster: reseed on the worst-fitted point
                    iseed       = maxloc(mind, dim=1)
                    cen(:,s)    = self%X(iseed,:)
                    mind(iseed) = -1._dp
                    l_reseed    = .true.
                endif
            end do
            if( nchanged == 0 .and. .not. l_reseed ) exit
        end do
        self%niters   = min(it, self%maxits)
        self%nchanged = nchanged
        ! order the clusters by decreasing population, ties by smallest member index (insertion sort)
        allocate(first_member(self%k), source=huge(1))
        do i = self%n, 1, -1
            first_member(memb(i)) = i
        end do
        allocate(order(self%k))
        order = [(s, s=1,self%k)]
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
        allocate(newlab(self%k))
        do s = 1, self%k
            newlab(order(s)) = s
            centroids(:,s)   = cen(:,order(s))
        end do
        do i = 1, self%n
            labels(i) = newlab(memb(i))
        end do
        deallocate(mind, csum, xbar, cen, cnt, memb, seeds, first_member, order, newlab)
    end subroutine cluster

    !>  \brief  Lloyd iterations of the last cluster call
    pure integer function get_niters( self )
        class(kmeans), intent(in) :: self
        get_niters = self%niters
    end function get_niters

    !>  \brief  reassignments on the last Lloyd iteration of the last cluster call (0: converged)
    pure integer function get_nchanged( self )
        class(kmeans), intent(in) :: self
        get_nchanged = self%nchanged
    end function get_nchanged

    subroutine kill( self )
        class(kmeans), intent(inout) :: self
        if( self%exists )then
            deallocate(self%X, self%wdim)
            self%n = 0; self%d = 0; self%k = 0; self%maxits = 0; self%niters = 0; self%nchanged = 0
            self%exists = .false.
        endif
    end subroutine kill

end module simple_kmeans
