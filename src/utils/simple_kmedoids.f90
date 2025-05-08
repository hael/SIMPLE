! k-medoids clustering of a similarity matrix
module simple_kmedoids
include 'simple_lib.f08'
implicit none

public  :: kmedoids
private
#include "simple_local_flags.inc"

type kmedoids
    private
    real,    pointer     :: ptr_dmat(:,:) => null()
    real,    allocatable :: dists2meds(:,:)
    integer              :: n=0, ncls=0
    integer              :: i_medoid_glob=0
    integer, allocatable :: i_medoids(:), cls_labels(:), cls_pops(:)
    logical              :: exists = .false.
  contains
    procedure, private   :: new_1, new_2
    generic              :: new => new_1, new_2
    procedure            :: get_labels
    procedure            :: init
    procedure            :: cluster
    procedure            :: find_medoids
    procedure            :: assign_labels
    procedure            :: merge
    procedure            :: kill
end type kmedoids

contains

    subroutine new_1( self, n, dmat, ncls )
        class(kmedoids), intent(inout) :: self
        integer,         intent(in)    :: n, ncls
        real, target,    intent(in)    :: dmat(n,n)
        call self%kill
        self%n        =  n
        self%ptr_dmat => dmat
        self%ncls     =  ncls
        allocate(self%i_medoids(self%ncls), self%cls_labels(self%n), self%cls_pops(self%ncls), source=0)
        allocate(self%dists2meds(self%n,self%ncls), source=0.)
        self%exists   = .true.
    end subroutine new_1

    subroutine new_2( self, cls_labels, dmat )
        class(kmedoids), intent(inout) :: self
        integer,         intent(in)    :: cls_labels(:)
        real, target,    intent(in)    :: dmat(size(cls_labels),size(cls_labels))
        call self%new_1(size(cls_labels), dmat, maxval(cls_labels))
        self%cls_labels = cls_labels
        call self%find_medoids
    end subroutine new_2

    subroutine get_labels( self, cls_labels )
        class(kmedoids), intent(inout) :: self
        integer,         intent(out)   :: cls_labels(self%n)
        cls_labels = self%cls_labels
    end subroutine get_labels

    ! initialize through distance to medoid analysis
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
        ! assign clusters based on similarity to medoid
        do icls = 1, self%ncls
            do i = parts(icls,1), parts(icls,2)
                self%cls_labels(order(i)) = icls
            end do
        end do
        ! set cluster populations
        do icls = 1, self%ncls
            self%cls_pops(icls) = count(self%cls_labels == icls) 
        end do
        deallocate(parts)
    end subroutine init

    subroutine cluster( self )
        class(kmedoids), intent(inout) :: self
        integer, parameter :: MAXITS = 10
        integer :: iter, nchanges, i
        real    :: dist
        if( all(self%cls_labels == 0) ) THROW_HARD('Cluster labels must be initialized somehow')
        write(logfhandle,'(A)') 'K-MEDOIDS CLUSTERING OF DISTANCE MATRIX'
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
            write(logfhandle,'(a,1x,f8.2)') 'ITER: '//int2str(iter)//' DIST: ', dist
            ! set l_converged flag
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

    subroutine assign_labels( self, nchanges )
        class(kmedoids), intent(inout) :: self
        integer,         intent(out)   :: nchanges
        integer :: i, icls, loc(self%n)
        real :: x, rhuge
        ! extract distances to medoids
        rhuge = huge(x)
        do i = 1, self%n
            do icls = 1, self%ncls
                if( self%i_medoids(icls) == 0 )then
                    self%dists2meds(i,icls) = rhuge
                else
                    self%dists2meds(i,icls) = self%ptr_dmat(i,self%i_medoids(icls))
                endif
            end do
        end do
        ! assign clusters
        loc = minloc(self%dists2meds, dim=2)
        nchanges = count(self%cls_labels /= loc)
        self%cls_labels = loc
    end subroutine assign_labels

    subroutine merge( self, ncls_new )
        class(kmedoids), intent(inout) :: self
        integer,         intent(in)    :: ncls_new
        real,    allocatable :: dists_btw_meds(:,:)
        logical, allocatable :: med_pair_mask(:,:)
        integer              :: loc(2)
        if( ncls_new >= self%ncls ) THROW_HARD('New # clusters (ncls_new) must be less than old # clusters')
        do
            call calc_med_dists
            loc = minloc(dists_btw_meds, mask=med_pair_mask)
            call merge_clusters(self%i_medoids(loc(1)), self%i_medoids(loc(2)))
            if( self%ncls == ncls_new ) exit
        end do
        if( allocated(dists_btw_meds) ) deallocate(dists_btw_meds)
        if( allocated(med_pair_mask)  ) deallocate(med_pair_mask)

    contains

        subroutine calc_med_dists
            integer :: imed, jmed
            ! extract distances between medoids
            if( allocated(dists_btw_meds) ) deallocate(dists_btw_meds)
            allocate(dists_btw_meds(self%ncls,self%ncls), source=0.)
            if( allocated(med_pair_mask)  ) deallocate(med_pair_mask)
            allocate(med_pair_mask(self%ncls,self%ncls),  source=.false.)
            do imed = 1, self%ncls - 1
                do jmed = imed + 1, self%ncls
                    dists_btw_meds(imed,jmed) = self%ptr_dmat(imed,jmed)
                    dists_btw_meds(jmed,imed) = dists_btw_meds(imed,jmed)
                    med_pair_mask(imed,jmed)  = .true.
                end do
            end do
        end subroutine calc_med_dists

        subroutine merge_clusters( i_cls, j_cls )
            integer, intent(in) :: i_cls, j_cls
            integer :: new_cls_labels(self%n), cnt, icls
            where(self%cls_labels == j_cls) self%cls_labels = i_cls
            ! reorder labels
            cnt            = 0
            new_cls_labels = 0
            do icls = 1, self%ncls
                if( any(self%cls_labels == icls) )then 
                    cnt = cnt + 1
                    where(self%cls_labels == icls) new_cls_labels = cnt
                endif
            enddo
            self%cls_labels = new_cls_labels
            ! update # clusters
            self%ncls = maxval(self%cls_labels)
            ! update medoids
            deallocate(self%i_medoids, self%cls_pops, self%dists2meds)
            allocate(self%i_medoids(self%ncls), self%cls_pops(self%ncls), source=0)
            allocate(self%dists2meds(self%n,self%ncls), source=0.)
            call self%find_medoids
        end subroutine merge_clusters

    end subroutine merge

    subroutine kill( self )
        class(kmedoids), intent(inout) :: self
        if( self%exists )then
            self%ptr_dmat => null()
            self%n    = 0
            self%ncls = 0
            deallocate(self%i_medoids, self%cls_labels, self%cls_pops, self%dists2meds)
        endif
    end subroutine kill

end module simple_kmedoids