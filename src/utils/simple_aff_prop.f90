! clustering based on a similartity matrix using affinity propagation
module simple_aff_prop
include 'simple_lib.f08'
implicit none

public :: aff_prop, test_aff_prop
public :: aggregate, silhouette_score
private
#include "simple_local_flags.inc"

type aff_prop
    private
    integer              :: N                    !< nr of data entries
    real,    allocatable :: A(:,:), R(:,:)       !< A is affinities & R responsibilities
    real,    allocatable :: Aold(:,:), Rold(:,:)
    real,    allocatable :: Rp(:,:), tmp(:)
    real,    allocatable :: AS(:,:), dA(:)
    real,    pointer     :: S(:,:)               !< pointer to similarity matrix
    real,    allocatable :: Y(:), Y2(:)          !< maxvals
    integer, allocatable :: I(:), I2(:)          !< index arrays
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
        real,    target,   intent(inout) :: S(N,N)          ! similarity matrix
        real,    optional, intent(in)    :: ftol, lam, pref ! tolerance, dampening factor, preference thres
        integer, optional, intent(in)    :: maxits
        integer :: i, j
        real    :: ppref, diff
        call self%kill
        ! set constants
        self%S => S
        self%N = N
        if( present(ftol)    ) self%ftol    = ftol
        if( present(lam)     ) self%lam     = lam
        if( present(maxits)  ) self%maxits  = maxits
        ! calculate similarity characteristics
        call analyze_smat(S, .false., self%Smin, self%Smax)
        ! low pref [0.1,1] leads to small numbers of clusters and
        ! high pref leads to large numbers of clusters. Frey suggests Smin - (Smax - Smin) as a threshold for fewer cluster
        ppref = self%Smin
        if( present(pref) ) ppref = pref
        ! remove degeneracies
        forall(i=1:N) self%S(i,i) = ppref
        diff = (self%Smax - self%Smin) * self%ftol
        do i=1,self%N
            do j=1,self%N
                self%S(i,j) = self%S(i,j) + ran3() * diff
            end do
        end do
        ! allocate
        allocate( self%A(N,N), self%R(N,N), self%Aold(N,N), self%Rp(N,N),&
        self%Rold(N,N), self%AS(N,N), self%Y(N), self%Y2(N), self%tmp(N),&
        self%I(N), self%I2(N), self%dA(N) )
        self%A      = 0.
        self%R      = 0.
        self%Aold   = 0.
        self%Rp     = 0.
        self%Rold   = 0.
        self%AS     = 0.
        self%Y      = 0.
        self%Y2     = 0.
        self%tmp    = 0.
        self%I      = 0
        self%I2     = 0
        self%dA     = 0.
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
        real              :: x, realmax
        integer           :: i, j, k, ncls
        ! initialize
        realmax   = huge(x)
        self%A    = 0.
        self%R    = 0.
        self%Aold = 0.
        self%Rp   = 0.
        self%Rold = 0.
        self%AS   = 0.
        self%Y    = 0.
        self%Y2   = 0.
        self%tmp  = 0.
        self%I    = 0
        self%I2   = 0
        self%dA   = 0.
        ! iterate
        do i=1,self%maxits
            ! FIRST, COMPUTE THE RESPONSIBILITIES
            self%Rold = self%R
            self%AS   = self%A + self%S
            self%I    = maxloc(self%AS, dim=2)
            forall(j=1:self%N) self%Y(j) = self%AS(j,self%I(j))
            do j=1,self%N
                self%AS(j,self%I(j)) = -realmax
            end do
            self%I2 = maxloc(self%AS, dim=2)
            forall(j=1:self%N) self%Y2(j) = self%AS(j,self%I2(j))
            self%R = self%S
            forall(j=1:self%N) self%R(j,:) = self%R(j,:) - self%Y
            do j=1,self%N
               self%R(j,self%I(j)) = self%S(j,self%I(j)) - self%Y2(j)
            end do
            ! update responsibilities (in a dampened fashion)
            self%R = (1. - self%lam) * self%R + self%lam * self%Rold
            ! THEN, COMPUTE THE AVAILABILITIES
            self%Aold = self%A
            where(self%R > 0.)
                self%Rp = self%R
            elsewhere
                self%Rp = 0.
            end where
            forall(k=1:self%N) self%Rp(k,k) = self%R(k,k)
            forall(k=1:self%N) self%tmp(k)  = sum(self%Rp(:,k))
            self%A = -self%Rp
            forall(j=1:self%N) self%A(j,:)  = self%A(j,:) + self%tmp
            forall(k=1:self%N) self%dA(k)   = self%A(k,k)
            where(self%A > 0.) self%A       = 0.
            forall(k=1:self%N) self%A(k,k)  = self%dA(k)
            ! update availabilities (in a dampened fashion)
            self%A = (1. - self%lam) * self%A + self%lam * self%Aold
        end do
        self%R = self%R + self%A ! pseudomarginals
        ! count the number of clusters
        ncls = 0
        do j=1,self%N
            if( self%R(j,j) > 0. ) ncls = ncls + 1
        end do
        if( allocated(centers) ) deallocate(centers)
        if( allocated(labels) )  deallocate(labels)
        allocate( centers(ncls), similarities(ncls), labels(self%N) )
        ! set the cluster centers
        ncls = 0
        do j=1,self%N
            if( self%R(j,j) > 0. )then
                ncls = ncls + 1
                centers(ncls) = j
            endif
        end do
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
            labels(j) = maxloc(similarities,dim=1)
            if( j .ne. centers(labels(j)) ) simsum = simsum + similarities(labels(j))
        end do
        simsum = simsum / real(self%N)
        deallocate(similarities)
    end subroutine propagate

    ! UNIT TEST

    !>  \brief  is the aff_prop unit test
    subroutine test_aff_prop
        real                 :: datavecs(900,5)
        type(aff_prop)       :: apcls
        real                 :: simmat(900,900), simsum
        integer, allocatable :: centers(:), labels(:)
        integer              :: i, j, ncls, nerr
        write(logfhandle,'(a)') '**info(simple_aff_prop_unit_test): testing all functionality'
        ! make data
        do i=1,300
            datavecs(i,:) = 1.
        end do
        do i=301,600
            datavecs(i,:) = 5.
        end do
        do i=601,900
            datavecs(i,:) = 10.
        end do
        do i=1,900-1
            simmat(i,i) = 0.
            do j=i+1,900
                simmat(i,j) = -euclid(datavecs(i,:),datavecs(j,:))
                simmat(j,i) = simmat(i,j)
            end do
        end do
        simmat(900,900) = 0.
        call apcls%new(900, simmat)
        call apcls%propagate(centers, labels, simsum)
        ncls = size(centers)
        nerr = 0
        do i=1,299
            do j=i+1,300
                if( labels(i) /= labels(j) ) nerr = nerr+1
            end do
        end do
        do i=301,599
            do j=i+1,600
                if( labels(i) /= labels(j) ) nerr = nerr+1
            end do
        end do
        do i=601,899
            do j=i+1,900
                if( labels(i) /= labels(j) ) nerr = nerr+1
            end do
        end do
        write(logfhandle,*) 'NR OF CLUSTERS FOUND:', ncls
        write(logfhandle,*) 'NR OF ASSIGNMENT ERRORS:', nerr
        write(logfhandle,*) 'CENTERS'
        do i=1,size(centers)
            write(logfhandle,*) datavecs(centers(i),:)
        end do
        if( ncls == 3 .and. nerr == 0 )then
            write(logfhandle,'(a)') 'SIMPLE_AFF_PROP_UNIT_TEST COMPLETED ;-)'
        else
            write(logfhandle,'(a)') 'SIMPLE_AFF_PROP_UNIT_TEST FAILED!'
        endif
    end subroutine test_aff_prop

    ! DESTRUCTOR

    subroutine kill( self )
        class(aff_prop), intent(inout) :: self
        if( self%exists )then
            self%S => null()
            deallocate( self%A, self%R, self%Aold, self%Rp,&
            self%Rold, self%AS, self%Y, self%Y2, self%tmp, self%I,&
            self%I2, self%dA )
            self%exists = .false.
        endif
    end subroutine kill

    ! PUBLIC METHODS

    !>  \brief Given labelling into k partitions and corresponding distance matrix iteratively
    !          combines the pair of clusters with lowest average inter-cluster distance,
    !          and produces labelling for partitions in [1;k] (1/k is for convenience)
    subroutine aggregate( labels, distmat, labelsmat )
        integer,              intent(in)    :: labels(:)
        real,                 intent(in)    :: distmat(:,:)
        integer, allocatable, intent(inout) :: labelsmat(:,:)
        real,    allocatable :: clmat(:,:)
        integer, allocatable :: iinds(:), jinds(:)
        real    :: d
        integer :: pair(2), nc, nl, cl, i,j,k, ni, aggc, oldc
        nl = size(labels)
        nc = maxval(labels) ! 1-base numbering assumed
        allocate(labelsmat(nl,nc), source=0)
        labelsmat(:,1)  = 1
        labelsmat(:,nc) = labels
        do cl = nc-1,2,-1
            ! average inter-cluster distance
            allocate(clmat(cl+1,cl+1),source=huge(d))
            do i = 1,cl
                iinds = pack((/(k,k=1,nl)/),mask=labelsmat(:,cl+1)==i)
                ni = size(iinds)
                do j = i+1, cl+1
                    jinds = pack((/(k,k=1,nl)/),mask=labelsmat(:,cl+1)==j)
                    d = 0.
                    do k = 1,ni
                        d = d + sum(distmat(iinds(k),jinds(:)))
                    enddo
                    clmat(i,j) = d / real(ni*size(jinds))
                enddo
            enddo
            ! select pait of clusters to aggregate
            pair = minloc(clmat)
            aggc = minval(pair)
            oldc = maxval(pair)
            ! update labelling
            labelsmat(:,cl) = labelsmat(:,cl+1)
            where(labelsmat(:,cl) == oldc) labelsmat(:,cl) = aggc
            where(labelsmat(:,cl) >  oldc) labelsmat(:,cl)  = labelsmat(:,cl)-1
            deallocate(clmat,iinds,jinds)
        enddo
    end subroutine aggregate

    ! Average silhouette score, for instance https://en.wikipedia.org/wiki/Silhouette_(clustering)
    ! -1 <= SC <= 1, higher SC => better clustering
    real function silhouette_score( labels, distmat )
        integer, intent(in)  :: labels(:)
        real,    intent(in)  :: distmat(:,:)
        integer, allocatable :: inds(:)
        real    :: a,b,s
        integer :: i,j,k, nl, nc, ni
        nl = size(labels)
        nc = maxval(labels)
        s  = 0.
        do i = 1,nl
            inds = pack((/(j,j=1,nl)/),mask=labels==labels(i))
            ni   = size(inds)
            if( ni == 1 ) cycle
            a    = sum(distmat(i,inds)) / real(ni-1) ! -1 because i is part of the set, and d(i,i)=0
            b    = huge(b)
            do k = 1,nc
                if( labels(i)==k ) cycle
                inds = pack((/(j,j=1,nl)/),mask=labels==k)
                if( .not.allocated(inds) ) cycle
                if( size(inds)==0 ) cycle
                b = min(b, sum(distmat(i,inds)) / real(size(inds)))
            enddo
            s = s + (b-a) / max(a,b)
        enddo
        silhouette_score = s / real(nl)
    end function silhouette_score

end module simple_aff_prop
