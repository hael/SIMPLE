! k-medoids clustering of a similarity matrix
module simple_kmedoids_online
include 'simple_lib.f08'
implicit none

public  :: kmedoids_online, distfun
private
#include "simple_local_flags.inc"

type kmedoids_online
    private
    procedure(distfun), pointer, nopass :: distfun => null() !< defines the distance function
    real,               allocatable     :: dists2meds(:,:)
    integer,            allocatable     :: i_medoids(:), cls_labels(:), cls_pops(:)
    integer :: n=0, ncls=0
    logical :: exists = .false.
  contains
    ! init
    procedure           :: new
    procedure           :: assign_labels
    ! refinement
    procedure           :: refine_clustering
    procedure, private  :: find_medoids
    procedure, private  :: find_medoid
    ! getters
    procedure           :: get_labels
    procedure           :: get_medoids
    ! destructor
    procedure           :: kill
end type kmedoids_online

!>  \brief  defines the distance function interface
abstract interface
    function distfun( i, j ) result( dist )
        integer,  intent(in) :: i, j
        real                 :: dist
    end function
end interface

contains

    subroutine new( self, n, i_medoids, distfun )
        class(kmedoids_online), intent(inout) :: self
        integer,                intent(in)    :: n
        integer,                intent(in)    :: i_medoids(:) ! need an initial estimate of medoids for online refinement
        interface
            function distfun( i, j ) result( dist )
                integer,  intent(in) :: i, j
                real                 :: dist
            end function
        end interface
        call self%kill
        self%distfun => distfun
        self%n       =  n
        self%ncls    =  size(i_medoids)
        allocate(self%i_medoids(self%ncls), source=i_medoids)
        allocate(self%cls_labels(self%n), self%cls_pops(self%ncls), source=0)
        allocate(self%dists2meds(self%n,self%ncls), source=0.)
        self%exists = .true.
    end subroutine new

    subroutine assign_labels( self, nchanges )
        class(kmedoids_online), intent(inout) :: self
        integer,                intent(out)   :: nchanges
        integer :: i, icls, loc(self%n)
        real :: x, rhuge
        ! extract distances to medoids
        rhuge = huge(x)
        do i = 1, self%n
            do icls = 1, self%ncls
                if( self%i_medoids(icls) == 0 )then
                    self%dists2meds(i,icls) = rhuge
                else
                    self%dists2meds(i,icls) = self%distfun(self%i_medoids(icls), i)
                endif
            end do
        end do
        ! assign clusters
        loc = minloc(self%dists2meds, dim=2)
        nchanges = count(self%cls_labels /= loc)
        self%cls_labels = loc
    end subroutine assign_labels

    subroutine refine_clustering( self )
        class(kmedoids_online), intent(inout) :: self
        integer, parameter :: MAXITS = 10
        integer :: iter, nchanges, i
        real    :: dist
        if( all(self%cls_labels == 0) ) THROW_HARD('Cluster labels must be initialized somehow')
        write(logfhandle,'(A)') 'ONLINE K-MEDOIDS CLUSTERING REFINEMENT'
        iter = 0
        do
            call self%find_medoids
            call self%assign_labels(nchanges)
            ! update iteration counter
            iter = iter + 1
            ! report joint distance
            dist = 0.
            do i = 1, self%n
                dist = dist + self%dists2meds(i,self%cls_labels(i))
            end do
            write(logfhandle,'(a,1x,f8.2)') 'ITER: '//int2str_pad(iter,2)//' DIST: ', dist
            ! set l_converged flag
            if( nchanges == 0 .or. iter == MAXITS) exit
        end do
    end subroutine refine_clustering

    subroutine find_medoids( self )
        class(kmedoids_online), intent(inout) :: self
        integer :: icls
        !$omp parallel do default(shared) private(icls) proc_bind(close)
        do icls = 1, self%ncls
            call self%find_medoid(icls)
        enddo
        !$omp end parallel do
    end subroutine find_medoids

    subroutine find_medoid( self, icls )
        class(kmedoids_online), intent(inout) :: self
        integer,                intent(in)    :: icls
        real    :: dists(self%n)
        integer :: i, j, loc(1)
        self%cls_pops(icls) = count(self%cls_labels == icls)
        if( self%cls_pops(icls) == 0 )then
            self%i_medoids(icls) = 0
        else if( self%cls_pops(icls) == 1 )then
            do i = 1, self%n
                if( self%cls_labels(i) == icls )then
                    self%i_medoids(icls) = i
                    exit
                endif
            enddo
        else
            do i = 1, self%n
                dists(i) = 0.                        
                do j = 1, self%n
                    if( self%cls_labels(i) == icls .and. self%cls_labels(j) == icls )then
                        dists(i) = dists(i) + self%distfun(i, j) 
                    endif
                enddo
            end do
            loc = minloc(dists, mask=self%cls_labels == icls)
            self%i_medoids(icls) = loc(1)
        endif
    end subroutine find_medoid

    subroutine get_labels( self, cls_labels )
        class(kmedoids_online), intent(inout) :: self
        integer,                intent(out)   :: cls_labels(self%n)
        cls_labels = self%cls_labels
    end subroutine get_labels

    subroutine get_medoids( self, i_medoids )
        class(kmedoids_online), intent(inout) :: self
        integer,                intent(out)   :: i_medoids(self%ncls)
        i_medoids = self%i_medoids
    end subroutine get_medoids

    subroutine kill( self )
        class(kmedoids_online), intent(inout) :: self
        if( self%exists )then
            self%distfun => null()
            self%n       =  0
            self%ncls    =  0
            deallocate(self%i_medoids, self%cls_labels, self%cls_pops, self%dists2meds)
            self%exists   =  .false.
        endif
    end subroutine kill

end module simple_kmedoids_online
