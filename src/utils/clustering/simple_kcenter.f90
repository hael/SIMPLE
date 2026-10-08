!@descr: greedy farthest-point k-center clustering of feature vectors (Gonzalez), deterministic
! Points are the rows of X(n,d). The first center is the point farthest from the mean, every next one the
! point farthest from all centers chosen so far (simple_clustering_utils), the first index on a tie; each
! point joins its nearest center. The centers cover the point set within twice the optimal radius, so a
! sparse region gets a center that k-means would spend on a dense one. Labels by decreasing population
! (ties by smallest member index).
module simple_kcenter
use simple_core_module_api,  only: dp, simple_exception
use simple_clustering_utils, only: farthest_point_seeds
implicit none

public :: kcenter
private
#include "simple_local_flags.inc"

type kcenter
    private
    integer               :: n = 0, d = 0, k = 0
    real(dp), allocatable :: X(:,:)          !< (n,d) points
    logical               :: exists = .false.
  contains
    procedure :: new
    procedure :: cluster
    procedure :: kill
end type kcenter

contains

    !>  \brief  is a constructor
    subroutine new( self, X, k )
        class(kcenter), intent(inout) :: self
        real(dp),       intent(in)    :: X(:,:) !< (n,d) points
        integer,        intent(in)    :: k      !< nr of centers
        call self%kill
        self%n = size(X,1)
        self%d = size(X,2)
        if( self%n < 1 .or. self%d < 1 ) THROW_HARD('empty point set; new')
        if( k < 1 .or. k > self%n )      THROW_HARD('k must be in [1,n]; new')
        self%k = k
        allocate(self%X(self%n,self%d), source=X)
        self%exists = .true.
    end subroutine new

    !>  \brief  chooses the centers and assigns every point to its nearest one
    subroutine cluster( self, centers, labels )
        class(kcenter), intent(inout) :: self
        integer,        intent(out)   :: centers(:) !< (k) point index of every center, by label
        integer,        intent(out)   :: labels(:)  !< (n) cluster labels, 1 the most populated
        real(dp), allocatable :: xbar(:)
        integer,  allocatable :: seeds(:), memb(:), cnt(:), first_member(:), order(:), newlab(:)
        real(dp) :: d2, best
        integer  :: i, q, c, ibest, j, tmp
        if( .not. self%exists ) THROW_HARD('object does not exist; cluster')
        if( size(centers) /= self%k ) THROW_HARD('centers must have k entries; cluster')
        if( size(labels)  /= self%n ) THROW_HARD('labels must have one entry per point; cluster')
        allocate(xbar(self%d), seeds(self%k), memb(self%n), cnt(self%k))
        ! first center: the point farthest from the mean
        do q = 1, self%d
            xbar(q) = sum(self%X(:,q)) / real(self%n,dp)
        end do
        best = -1._dp; ibest = 1
        do i = 1, self%n
            d2 = 0._dp
            do q = 1, self%d
                d2 = d2 + (self%X(i,q) - xbar(q))**2
            end do
            if( d2 > best )then
                best  = d2
                ibest = i
            endif
        end do
        call farthest_point_seeds(self%X, self%k, ibest, seeds)
        ! nearest center, the first on a tie
        cnt = 0
        do i = 1, self%n
            best = huge(1._dp); ibest = 1
            do c = 1, self%k
                d2 = 0._dp
                do q = 1, self%d
                    d2 = d2 + (self%X(i,q) - self%X(seeds(c),q))**2
                end do
                if( d2 < best )then
                    best  = d2
                    ibest = c
                endif
            end do
            memb(i)    = ibest
            cnt(ibest) = cnt(ibest) + 1
        end do
        ! order the clusters by decreasing population, ties by smallest member index (insertion sort)
        allocate(first_member(self%k), source=huge(1))
        do i = self%n, 1, -1
            first_member(memb(i)) = i
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
        allocate(newlab(self%k))
        do c = 1, self%k
            newlab(order(c)) = c
            centers(c)       = seeds(order(c))
        end do
        do i = 1, self%n
            labels(i) = newlab(memb(i))
        end do
        deallocate(xbar, seeds, memb, cnt, first_member, order, newlab)
    end subroutine cluster

    subroutine kill( self )
        class(kcenter), intent(inout) :: self
        if( self%exists )then
            deallocate(self%X)
            self%n = 0; self%d = 0; self%k = 0
            self%exists = .false.
        endif
    end subroutine kill

end module simple_kcenter
