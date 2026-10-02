!@descr: clustering of a distance matrix by agglomerative hierarchical clustering with average linkage (UPGMA)
! Labels by decreasing population (ties by smallest member index); medoid = minimal summed distance.
module simple_avg_linkage
use simple_core_module_api
implicit none

public :: avg_linkage
private
#include "simple_local_flags.inc"

type avg_linkage
    private
    integer               :: N = 0             !< nr of data entries
    integer               :: nclust = 0        !< nr of clusters at which merging stops
    real,     allocatable :: D(:,:)            !< private distance matrix
    real(dp), allocatable :: dlink(:,:)        !< mean pairwise distance between clusters (active rows/columns)
    integer,  allocatable :: rep(:)            !< cluster of every entry, named by its smallest member index
    integer,  allocatable :: csize(:)          !< population of the cluster named by an entry
    logical,  allocatable :: active(:)         !< entry names a cluster that still exists
    logical               :: exists = .false.  !< to indicate existence
  contains
    procedure :: new
    procedure :: cluster
    procedure :: kill
    procedure, private :: merge2nclust
    procedure, private :: label_clusters
    procedure, private :: find_medoids
end type avg_linkage

contains

    ! CONSTRUCTOR

    !>  \brief  is a constructor
    subroutine new( self, N, D, nclust )
        class(avg_linkage), intent(inout) :: self
        integer,            intent(in)    :: N       ! # data entries
        real,               intent(in)    :: D(N,N)  ! distance matrix
        integer,            intent(in)    :: nclust  ! # clusters to form
        call self%kill
        if( N < 1 )                       THROW_HARD('N must be >= 1; new')
        if( nclust < 1 .or. nclust > N ) THROW_HARD('nclust must be in [1,N]; new')
        self%N      = N
        self%nclust = nclust
        allocate(self%D(N,N), source=D)
        allocate(self%dlink(N,N), self%rep(N), self%csize(N), self%active(N))
        self%exists = .true.
    end subroutine new

    ! CLUSTERING

    !>  \brief  merges the closest clusters until nclust remain, then labels them and finds their medoids
    subroutine cluster( self, medoids, labels )
        class(avg_linkage),   intent(inout) :: self
        integer, allocatable, intent(inout) :: medoids(:) !< medoid of every cluster
        integer, allocatable, intent(inout) :: labels(:)  !< cluster labels, 1 the most populated
        if( .not. self%exists ) THROW_HARD('object does not exist; cluster')
        if( allocated(medoids) ) deallocate(medoids)
        if( allocated(labels)  ) deallocate(labels)
        call self%merge2nclust
        call self%label_clusters(labels)
        call self%find_medoids(labels, medoids)
    end subroutine cluster

    !>  \brief  agglomerates from singletons until nclust clusters are active
    subroutine merge2nclust( self )
        class(avg_linkage), intent(inout) :: self
        real(dp) :: dmin
        integer  :: nactive, i, j, k, imin, jmin
        self%dlink  = real(self%D, dp)
        self%rep    = [(i, i=1,self%N)]
        self%csize  = 1
        self%active = .true.
        nactive     = self%N
        do while( nactive > self%nclust )
            ! closest pair of active clusters, the first one on a tie (i < j)
            dmin = huge(dmin)
            imin = 0
            jmin = 0
            do j = 2, self%N
                if( .not. self%active(j) ) cycle
                do i = 1, j - 1
                    if( .not. self%active(i) ) cycle
                    if( self%dlink(i,j) < dmin )then
                        dmin = self%dlink(i,j)
                        imin = i
                        jmin = j
                    endif
                enddo
            enddo
            if( imin == 0 ) THROW_HARD('non-finite distances; merge2nclust')
            ! merge jmin into imin; the mean distance to every other cluster
            ! is the population-weighted mean of the two
            do k = 1, self%N
                if( .not. self%active(k) .or. k == imin .or. k == jmin ) cycle
                self%dlink(imin,k) = (real(self%csize(imin),dp) * self%dlink(imin,k) + &
                    &real(self%csize(jmin),dp) * self%dlink(jmin,k)) / real(self%csize(imin) + self%csize(jmin), dp)
                self%dlink(k,imin) = self%dlink(imin,k)
            enddo
            self%csize(imin)  = self%csize(imin) + self%csize(jmin)
            self%active(jmin) = .false.
            where( self%rep == jmin ) self%rep = imin
            nactive = nactive - 1
        enddo
    end subroutine merge2nclust

    !>  \brief  labels the active clusters by decreasing population, ties by representative
    subroutine label_clusters( self, labels )
        class(avg_linkage),   intent(in)    :: self
        integer, allocatable, intent(inout) :: labels(:)
        integer, allocatable :: reps(:)
        integer :: i, j, ic, tmp
        reps = pack([(i, i=1,self%N)], mask=self%active)
        ! insertion sort, stable on the representative order of pack
        do i = 2, self%nclust
            j = i
            do while( j > 1 )
                if( self%csize(reps(j-1)) >= self%csize(reps(j)) ) exit
                tmp       = reps(j-1)
                reps(j-1) = reps(j)
                reps(j)   = tmp
                j = j - 1
            enddo
        enddo
        allocate(labels(self%N), source=0)
        do ic = 1, self%nclust
            where( self%rep == reps(ic) ) labels = ic
        enddo
        deallocate(reps)
    end subroutine label_clusters

    !>  \brief  the member of every cluster with the smallest summed distance to the others, first on a tie
    subroutine find_medoids( self, labels, medoids )
        class(avg_linkage),   intent(in)    :: self
        integer,              intent(in)    :: labels(:)
        integer, allocatable, intent(inout) :: medoids(:)
        integer, allocatable :: members(:)
        real,    allocatable :: dsums(:)
        integer :: i, ic, loc(1)
        allocate(medoids(self%nclust), source=0)
        do ic = 1, self%nclust
            members = pack([(i, i=1,self%N)], mask=labels == ic)
            allocate(dsums(size(members)), source=0.)
            do i = 1, size(members)
                dsums(i) = sum(self%D(members(i),members))
            enddo
            loc = minloc(dsums)
            medoids(ic) = members(loc(1))
            deallocate(dsums, members)
        enddo
    end subroutine find_medoids

    ! DESTRUCTOR

    subroutine kill( self )
        class(avg_linkage), intent(inout) :: self
        if( self%exists )then
            deallocate(self%D, self%dlink, self%rep, self%csize, self%active)
            self%N      = 0
            self%nclust = 0
            self%exists = .false.
        endif
    end subroutine kill

end module simple_avg_linkage
