!@descr: agglomerative hierarchical clustering of a distance matrix using nearest-neighbor chain
module simple_hclust
use simple_core_module_api
implicit none

public :: hclust, test_hclust
private
#include "simple_local_flags.inc"

integer, parameter :: LINK_SINGLE   = 1
integer, parameter :: LINK_COMPLETE = 2
integer, parameter :: LINK_AVERAGE  = 3

type hclust
    private
    integer              :: N = 0
    integer              :: linkage = LINK_AVERAGE
    real,    allocatable :: dmat(:,:) !< working distances
    integer, allocatable :: size(:)
    integer, allocatable :: active(:)
    integer, allocatable :: map(:)
    integer, allocatable :: chain(:)
    integer              :: step = 0
    logical              :: exists = .false.
  contains
    procedure :: new
    procedure :: cluster
    procedure :: kill
end type hclust

contains

    !===============================================================================
    ! CONSTRUCTOR
    !===============================================================================

    subroutine new(self, N, dmat, linkage)
        class(hclust),     intent(inout) :: self
        integer,           intent(in)    :: N
        real,    target,   intent(in)    :: dmat(N,N)
        integer, optional, intent(in)    :: linkage
        integer :: i
        call self%kill
        self%N = N
        if(present(linkage)) self%linkage = linkage
        allocate(self%dmat(N,N), source=dmat)
        allocate(self%size(N), self%active(N), self%map(N), self%chain(N))
        self%size   = 1
        self%active = 1
        self%map    = [(i,i=1,N)]
        self%step   = 0
        self%exists = .true.
    end subroutine new

    !===============================================================================
    ! CLUSTERING (Nearest-Neighbor Chain)
    !===============================================================================

    subroutine cluster(self, merge_mat, height, labels, n_clusters)
        class(hclust), intent(inout)         :: self
        integer,       intent(out)           :: merge_mat(:,:)
        real,          intent(out)           :: height(:)
        integer,       intent(out), optional :: labels(:)
        integer,       intent(in),  optional :: n_clusters
        integer :: i, j, k, x, y, nn, na, nb, nc, clen
        real    :: dmin, dnew, val
        logical :: l_msk(self%N)
        if (.not. self%exists) THROW_HARD("hclust not initialized")
        do while (count(self%active == 1) > 1)
            !-----------------------------------------------------------
            ! start chain
            !-----------------------------------------------------------
            do i = 1, self%N
                if (self%active(i) == 1) then
                    self%chain(1) = i
                    clen = 1
                    exit
                endif
            end do
            !-----------------------------------------------------------
            ! grow nearest-neighbor chain
            !-----------------------------------------------------------
            do
                x    = self%chain(clen)
                nn   = -1
                dmin = huge(dmin)
                l_msk = (self%active == 1)
                l_msk(x) = .false. 
                l_msk = l_msk(self%N:1:-1)
                nn = self%N - minloc(self%dmat(x,self%N:1:-1), dim = 1, mask = l_msk) + 1
                if (clen > 1) then
                    y = self%chain(clen-1)
                    if (nn == y) exit
                endif
                clen = clen + 1
                self%chain(clen) = nn
            end do
            !-----------------------------------------------------------
            ! merge
            !-----------------------------------------------------------
            y                      = self%chain(clen-1)
            x                      = self%chain(clen)
            self%step              = self%step + 1
            merge_mat(1,self%step) = self%map(y)
            merge_mat(2,self%step) = self%map(x)
            height(self%step)      = self%dmat(y,x)
            na                     = self%size(y)
            nb                     = self%size(x)
            !-----------------------------------------------------------
            ! distance update (dominant cost)
            !-----------------------------------------------------------
            !$omp parallel do default(shared) private(k,nc,dnew)
            do k = 1, self%N
                if (self%active(k) == 0 .or. k == x .or. k == y) cycle
                nc = self%size(k)
                select case (self%linkage)
                    case (LINK_SINGLE)
                        dnew = min(self%dmat(y,k), self%dmat(x,k))
                    case (LINK_COMPLETE)
                        dnew = max(self%dmat(y,k), self%dmat(x,k))
                    case (LINK_AVERAGE)
                        dnew = (real(na)*self%dmat(y,k) + real(nb)*self%dmat(x,k)) / real(na)+real(nb)
                end select
                self%dmat(y,k) = dnew
                self%dmat(k,y) = dnew
            end do
            !$omp end parallel do
            self%active(x) = 0
            self%size(y)   = na + nb
            self%map(y)    = self%N + self%step
        end do
        if (present(n_clusters)) then
            if( .not. present(labels) ) THROW_HARD('Need labels optional dummy input when n_clusters are given')
            call cut_tree(self%N, merge_mat, n_clusters, labels)
        endif
    end subroutine cluster

    !===============================================================================
    ! TREE CUT
    !===============================================================================

    subroutine cut_tree(N, merge_mat, k, labels)
        integer, intent(in)  :: N, k
        integer, intent(in)  :: merge_mat(2,N-1)
        integer, intent(out) :: labels(N)
        integer :: parent(2*N), node_label(2*N - 1)
        integer :: i, s, c, lab
        parent = [(i,i=1,2*N)]
        node_label = 0
        do s = 1, N-k
            parent(merge_mat(1,s)) = N + s
            parent(merge_mat(2,s)) = N + s
        end do 
        labels = 0
        lab = 0
        do i = 1, N
            c = i
            do while (parent(c) /= c)
                c = parent(c)
            end do
            if (node_label(c) == 0) then
                lab = lab + 1
                node_label(c) = lab
            end if
            labels(i) = node_label(c)
        end do
    end subroutine cut_tree

    !===============================================================================
    ! DESTRUCTOR
    !===============================================================================

    subroutine kill(self)
        class(hclust), intent(inout) :: self
        if (self%exists) then
            deallocate(self%dmat, self%size, self%active, self%map, self%chain)
            self%exists = .false.
        endif
    end subroutine kill

    !===============================================================================
    ! UNIT TEST
    !===============================================================================

    subroutine test_hclust
        real         :: data(300,2), D(300,300)
        type(hclust) :: hc
        integer :: merge_mat(2,299)
        real    :: height(299)
        integer :: labels(300)
        integer :: i, j, nerr
        do i=1,100; data(i,:) = 0.; enddo
        do i=101,200; data(i,:) = 5.; enddo
        do i=201,300; data(i,:) = 10.; enddo
        do i=1,300
            do j=i,300
                D(i,j) = euclid(data(i,:),data(j,:))
                D(j,i) = D(i,j)
            end do
        end do
        call hc%new(300, D, LINK_AVERAGE)
        call hc%cluster(merge_mat, height, labels, 3)
        nerr = 0
        do i=1,99
            do j=i+1,100
                if( labels(i) /= labels(j) ) nerr = nerr+1
            end do
        end do
        do i=101,199
            do j=i+1,200
                if( labels(i) /= labels(j) ) nerr = nerr+1
            end do
        end do
        do i=201,299
            do j=i+1,300
                if( labels(i) /= labels(j) ) nerr = nerr+1
            end do
        end do
        write(logfhandle,*) 'NR OF ASSIGNMENT ERRORS:', nerr
        if( nerr == 0 )then
            write(logfhandle,'(a)') 'SIMPLE_HCLUST_UNIT_TEST COMPLETED ;-)'
        else
            write(logfhandle,'(a)') 'SIMPLE_HCLUST_UNIT_TEST FAILED!'
        endif
    end subroutine test_hclust

end module simple_hclust
