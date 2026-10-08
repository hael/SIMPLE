!@descr: agglomerative hierarchical clustering of a distance matrix, average (UPGMA) or complete linkage
! Merging stops at a cluster count or at a distance threshold (exactly one is given). Labels by decreasing
! population (ties by smallest member index); medoid = minimal summed distance. Entries outside the
! optional mask never merge and stay singletons. The merge history (the two clusters joined, named by
! their smallest member index, and the linkage distance) is kept for callers that report it.
module simple_hac
use simple_core_module_api, only: dp, simple_exception
implicit none

public :: hac
private
#include "simple_local_flags.inc"

type hac
    private
    integer               :: N          = 0       !< nr of data entries
    integer               :: nclust     = 0       !< nr of clusters at which merging stops (0: threshold stop)
    integer               :: nclust_out = 0       !< nr of clusters the last cluster call formed
    integer               :: nmerge     = 0       !< nr of merges of the last cluster call
    real(dp)              :: thres      = 0._dp   !< largest linkage distance that merges (threshold stop)
    logical               :: l_complete = .false. !< complete linkage, else average
    real(dp), allocatable :: D(:,:)               !< private distance matrix
    real(dp), allocatable :: dlink(:,:)           !< linkage distance between clusters (active rows/columns)
    integer,  allocatable :: rep(:)               !< cluster of every entry, named by its smallest member index
    integer,  allocatable :: csize(:)             !< population of the cluster named by an entry
    logical,  allocatable :: active(:)            !< entry names a cluster that still exists
    logical,  allocatable :: incl(:)              !< entry may merge (the mask)
    integer,  allocatable :: hpairs(:,:)          !< (2,N-1) clusters joined at each merge, smaller name first
    real(dp), allocatable :: heights(:)           !< (N-1) linkage distance of each merge
    logical               :: exists = .false.     !< to indicate existence
  contains
    procedure :: new
    procedure :: cluster
    procedure :: get_nclust
    procedure :: get_history
    procedure :: kill
    procedure, private :: merge_clusters
    procedure, private :: label_clusters
    procedure, private :: find_medoids
end type hac

contains

    ! CONSTRUCTOR

    !>  \brief  is a constructor; linkage is 'average' or 'complete', and exactly one of nclust and thres
    !!          ends the merging
    subroutine new( self, N, D, linkage, nclust, thres, mask )
        class(hac),         intent(inout) :: self
        integer,            intent(in)    :: N       ! # data entries
        real(dp),           intent(in)    :: D(N,N)  ! distance matrix
        character(len=*),   intent(in)    :: linkage ! 'average' | 'complete'
        integer,  optional, intent(in)    :: nclust  ! # clusters to form
        real(dp), optional, intent(in)    :: thres   ! largest linkage distance that merges
        logical,  optional, intent(in)    :: mask(N) ! .false. entries stay singletons
        call self%kill
        if( N < 1 ) THROW_HARD('N must be >= 1; new')
        if( present(nclust) .eqv. present(thres) ) THROW_HARD('give exactly one of nclust and thres; new')
        select case(trim(linkage))
            case('average')
                self%l_complete = .false.
            case('complete')
                self%l_complete = .true.
            case DEFAULT
                THROW_HARD('unsupported linkage: '//trim(linkage)//'; new')
        end select
        self%N = N
        if( present(nclust) )then
            if( nclust < 1 .or. nclust > N ) THROW_HARD('nclust must be in [1,N]; new')
            self%nclust = nclust
        else
            self%nclust = 0
            self%thres  = thres
        endif
        allocate(self%D(N,N), source=D)
        allocate(self%dlink(N,N), self%rep(N), self%csize(N), self%active(N), self%incl(N))
        allocate(self%hpairs(2,max(1,N-1)), self%heights(max(1,N-1)))
        self%incl = .true.
        if( present(mask) ) self%incl = mask
        self%exists = .true.
    end subroutine new

    ! CLUSTERING

    !>  \brief  merges the closest clusters until the stop rule holds, then labels them and finds their medoids
    subroutine cluster( self, medoids, labels )
        class(hac),           intent(inout) :: self
        integer, allocatable, intent(inout) :: medoids(:) !< medoid of every cluster
        integer, allocatable, intent(inout) :: labels(:)  !< cluster labels, 1 the most populated
        if( .not. self%exists ) THROW_HARD('object does not exist; cluster')
        if( allocated(medoids) ) deallocate(medoids)
        if( allocated(labels)  ) deallocate(labels)
        call self%merge_clusters
        call self%label_clusters(labels)
        call self%find_medoids(labels, medoids)
    end subroutine cluster

    !>  \brief  agglomerates from singletons: the closest pair of active, unmasked clusters merges while
    !!          more than nclust clusters remain (count stop) or while its distance is <= thres
    !!          (threshold stop). Clusters are a union-find flattened into rep: every entry names its
    !!          cluster's smallest member index, and a merge renames the larger name to the smaller.
    subroutine merge_clusters( self )
        class(hac), intent(inout) :: self
        real(dp) :: dmin
        integer  :: nactive, i, j, k, imin, jmin
        self%dlink  = self%D
        self%rep    = [(i, i=1,self%N)]
        self%csize  = 1
        self%active = .true.
        self%nmerge = 0
        nactive     = self%N
        do
            if( self%nclust > 0 )then
                if( nactive <= self%nclust ) exit
            endif
            ! closest pair of active, unmasked clusters, the first one on a tie (j, then i < j)
            dmin = huge(dmin)
            imin = 0
            jmin = 0
            do j = 2, self%N
                if( .not. self%active(j) ) cycle
                if( .not. self%incl(j)   ) cycle
                do i = 1, j - 1
                    if( .not. self%active(i) ) cycle
                    if( .not. self%incl(i)   ) cycle
                    if( self%dlink(i,j) < dmin )then
                        dmin = self%dlink(i,j)
                        imin = i
                        jmin = j
                    endif
                enddo
            enddo
            if( imin == 0 ) exit                                  ! no mergeable pair left
            if( self%nclust == 0 )then
                if( dmin > self%thres ) exit
            endif
            ! merge jmin into imin; the linkage distance to every other cluster is the
            ! population-weighted mean of the two (average) or the larger of the two (complete)
            do k = 1, self%N
                if( .not. self%active(k) .or. k == imin .or. k == jmin ) cycle
                if( self%l_complete )then
                    self%dlink(imin,k) = max(self%dlink(imin,k), self%dlink(jmin,k))
                else
                    self%dlink(imin,k) = (real(self%csize(imin),dp) * self%dlink(imin,k) + &
                        &real(self%csize(jmin),dp) * self%dlink(jmin,k)) / real(self%csize(imin) + self%csize(jmin), dp)
                endif
                self%dlink(k,imin) = self%dlink(imin,k)
            enddo
            self%csize(imin)  = self%csize(imin) + self%csize(jmin)
            self%active(jmin) = .false.
            where( self%rep == jmin ) self%rep = imin
            nactive     = nactive - 1
            self%nmerge = self%nmerge + 1
            self%hpairs(:,self%nmerge) = [imin, jmin]
            self%heights(self%nmerge)  = dmin
        enddo
        self%nclust_out = nactive
    end subroutine merge_clusters

    !>  \brief  labels the active clusters by decreasing population, ties by representative
    subroutine label_clusters( self, labels )
        class(hac),           intent(in)    :: self
        integer, allocatable, intent(inout) :: labels(:)
        integer, allocatable :: reps(:)
        integer :: i, j, ic, tmp
        reps = pack([(i, i=1,self%N)], mask=self%active)
        ! insertion sort, stable on the representative order of pack
        do i = 2, self%nclust_out
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
        do ic = 1, self%nclust_out
            where( self%rep == reps(ic) ) labels = ic
        enddo
        deallocate(reps)
    end subroutine label_clusters

    !>  \brief  the member of every cluster with the smallest summed distance to the others, first on a tie
    subroutine find_medoids( self, labels, medoids )
        class(hac),           intent(in)    :: self
        integer,              intent(in)    :: labels(:)
        integer, allocatable, intent(inout) :: medoids(:)
        integer,  allocatable :: members(:)
        real(dp), allocatable :: dsums(:)
        integer :: i, ic, loc(1)
        allocate(medoids(self%nclust_out), source=0)
        do ic = 1, self%nclust_out
            members = pack([(i, i=1,self%N)], mask=labels == ic)
            allocate(dsums(size(members)), source=0._dp)
            do i = 1, size(members)
                dsums(i) = sum(self%D(members(i),members))
            enddo
            loc = minloc(dsums)
            medoids(ic) = members(loc(1))
            deallocate(dsums, members)
        enddo
    end subroutine find_medoids

    ! GETTERS

    !>  \brief  nr of clusters the last cluster call formed
    pure integer function get_nclust( self )
        class(hac), intent(in) :: self
        get_nclust = self%nclust_out
    end function get_nclust

    !>  \brief  the merges of the last cluster call in order: pairs(:,m) are the two clusters joined (named
    !!          by their smallest member index, smaller first) and heights(m) their linkage distance
    subroutine get_history( self, pairs, heights )
        class(hac),            intent(in)    :: self
        integer,  allocatable, intent(inout) :: pairs(:,:)
        real(dp), allocatable, intent(inout) :: heights(:)
        if( allocated(pairs)   ) deallocate(pairs)
        if( allocated(heights) ) deallocate(heights)
        allocate(pairs(2,self%nmerge), heights(self%nmerge))
        if( self%nmerge > 0 )then
            pairs   = self%hpairs(:,1:self%nmerge)
            heights = self%heights(1:self%nmerge)
        endif
    end subroutine get_history

    ! DESTRUCTOR

    subroutine kill( self )
        class(hac), intent(inout) :: self
        if( self%exists )then
            deallocate(self%D, self%dlink, self%rep, self%csize, self%active, self%incl, self%hpairs, self%heights)
            self%N          = 0
            self%nclust     = 0
            self%nclust_out = 0
            self%nmerge     = 0
            self%thres      = 0._dp
            self%l_complete = .false.
            self%exists     = .false.
        endif
    end subroutine kill

end module simple_hac
