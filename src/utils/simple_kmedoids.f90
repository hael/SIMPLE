! k-medoids clustering of a similarity matrix
module simple_kmedoids
include 'simple_lib.f08'
implicit none

public  :: kmedoids
private
#include "simple_local_flags.inc"

type kmedoids
    private
    real, pointer        :: ptr_dmat(:,:) => null()
    integer              :: n=0, ncls=0
    integer              :: i_medoid_glob = 0
    integer, allocatable :: i_medoids(:), cls_labels(:), cls_pops(:)
    logical              :: exists = .false.
  contains
    procedure            :: new
    procedure            :: assign_cls_labels
    procedure            :: init
    procedure            :: cluster
    procedure            :: find_medoids
    procedure            :: kill
end type kmedoids

contains

    subroutine new( self, n, dmat, ncls )
        class(kmedoids), intent(inout) :: self
        integer,         intent(in)    :: n, ncls
        real, target,    intent(in)    :: dmat(n,n)
        call self%kill
        self%n        =  n
        self%ptr_dmat => dmat
        self%ncls     =  ncls
        allocate(self%i_medoids(self%ncls),self%cls_labels(self%n), self%cls_pops(self%ncls), source=0)
        self%exists   = .true.
    end subroutine new

    subroutine set_cls_labels( self, cls_labels )
        class(kmedoids), intent(inout) :: self
        integer,         intent(in)    :: cls_labels(self%ncls)
        self%cls_labels = cls_labels
    end subroutine set_cls_labels

    ! initialize throuhg distance to medoid analysis
    subroutine init( self )
        class(kmedoids), intent(inout) :: self
        integer, allocatable :: parts(:,:)
        real    :: medoid_dists(self%n)
        integer :: order(self%n), i, icls
        call medoid_from_dmat(self%ptr_dmat, self%i_medoid_glob)
        ! order according to similarity to medoid
        !$omp parallel do default(shared) private(i) proc_bind(close)
        do i = 1, self%n
            medoid_dists(i) = self%ptr_dmat(i,self%i_medoid_glob)
        end do
        !$omp end parallel do
        order = (/(i,i=1,self%n)/)
        call hpsort(medoid_dists, order)
        parts = split_nobjs_even(self%n, self%ncls)
        ! assign classes based on similarity to medoid
        do icls = 1, self%ncls
            do i = parts(icls,1), parts(icls,2)
                self%cls_labels(order(i)) = icls
            end do
        end do
        ! set class populations
        do icls = 1, self%ncls
            self%cls_pops(icls) = count(self%cls_labels == icls) 
        end do
        deallocate(parts)
    end subroutine init

    subroutine cluster( self )
        class(kmedoids), intent(inout) :: self
        integer, parameter :: MAXITS = 10
        integer :: iter, nchanges
        if( all(self%cls_labels == 0) ) THROW_HARD('Cluster labels must be initialized somehow')
        write(logfhandle,'(A)') 'K-MEDOIDS CLUSTERING OF DISTANCE MATRIX'
        iter = 0
        do
            call self%find_medoids
            call self%assign_cls_labels(nchanges)
            if( nchanges == 0 .or. iter == MAXITS) exit
        end do
    end subroutine cluster

    subroutine find_medoids( self )
        class(kmedoids), intent(inout) :: self
        real    :: dists(self%n)
        integer :: icls, i, j, loc(1)
        !$omp parallel do default(shared) private(icls,i,j,dists,loc) proc_bind(close)
        do icls = 1, self%ncls
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
                            dists(i) = dists(i) + self%ptr_dmat(i,j) 
                        endif
                    enddo
                end do
                loc = minloc(dists, mask=self%cls_labels == icls)
                self%i_medoids(icls) = loc(1)
            endif
        enddo
        !$omp end parallel do
    end subroutine find_medoids

    subroutine assign_cls_labels( self, nchanges )
        class(kmedoids), intent(inout) :: self
        integer,         intent(out)   :: nchanges
        integer :: i, icls, loc(self%n)
        real :: x, dists2meds(self%n,self%ncls), rhuge
        ! extract distances to medoids
        rhuge = huge(x)
        do i = 1, self%n
            do icls = 1, self%ncls
                if( self%i_medoids(icls) == 0 )then
                    dists2meds(i,icls) = rhuge
                else
                    dists2meds(i,icls) = self%ptr_dmat(i,self%i_medoids(icls))
                endif
            end do
        end do
        ! assign clusters
        loc = minloc(dists2meds, dim=2)
        nchanges = count(self%cls_labels /= loc)
        self%cls_labels = loc
    end subroutine assign_cls_labels

    subroutine kill( self )
        class(kmedoids), intent(inout) :: self
        if( self%exists )then
            self%ptr_dmat => null()
            self%n    = 0
            self%ncls = 0
            deallocate(self%i_medoids, self%cls_labels, self%cls_pops)
        endif
    end subroutine kill

end module simple_kmedoids