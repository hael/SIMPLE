! clustering based on a similartity matrix using affinity propagation
module simple_aff_prop
!$ use omp_lib
use simple_core_module_api
implicit none

public :: aff_prop, test_aff_prop
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
        real              :: x, realmax, maxdiff, val
        integer           :: i, j, k, ncls
        logical           :: converged
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
            !------------------------
            ! RESPONSIBILITIES
            !------------------------
            self%Rold = self%R
            self%AS   = self%A + self%S
            self%I    = maxloc(self%AS, dim=2)
            !$omp parallel do default(shared) private(j,k)
            do j = 1, self%N
                self%Y(j) = self%AS(j,self%I(j))
                self%AS(j,self%I(j)) = -realmax
            enddo
            !$omp end parallel do
            self%I2 = maxloc(self%AS, dim=2)
            !$omp parallel default(shared) private(j)
            !$omp do
            do j = 1, self%N
                self%Y2(j) = self%AS(j,self%I2(j))
            enddo
            !$omp end do
            !$omp workshare
            self%R = self%S
            !$omp end workshare
            !$omp do
            do j = 1, self%N
                self%R(j,:) = self%R(j,:) - self%Y
            enddo
            !$omp end do
            !$omp do 
            do j=1,self%N
                self%R(j,self%I(j)) = self%S(j,self%I(j)) - self%Y2(j)
            end do
            !$omp end do
            !$omp workshare
            ! update responsibilities (in a dampened fashion)
            self%R = (1. - self%lam) * self%R + self%lam * self%Rold
            !------------------------
            ! AVAILABILITIES
            !------------------------
            self%Aold = self%A
            where(self%R > 0.)
                self%Rp = self%R
            elsewhere
                self%Rp = 0.
            end where
            !$omp end workshare
            !$omp do
            do k = 1, self%N
                self%Rp(k,k) = self%R(k,k)
            end do
            !$omp end do
            !$omp do
            do k = 1, self%N
                self%tmp(k)  = sum(self%Rp(:,k))
            end do
            !$omp end do
            !$omp workshare
            self%A = -self%Rp
            !$omp end workshare
            !$omp do
            do j = 1, self%N
                self%A(j,:) = self%A(j,:) + self%tmp
            end do
            !$omp end do
            !$omp do
            do k = 1, self%N
                self%dA(k) = self%A(k,k)
            end do
            !$omp end do
            !$omp workshare
            where(self%A > 0.) self%A = 0.
            !$omp end workshare
            !$omp do
            do k = 1, self%N
                self%A(k,k) = self%dA(k)
            end do
            !$omp end do
            !$omp workshare
            ! update availabilities (in a dampened fashion)
            self%A = (1. - self%lam) * self%A + self%lam * self%Aold
            !$omp end workshare
            !$omp end parallel
            !========================
            ! Convergence check each 5th iter
            !========================
            if (mod(i,5) == 0) then
                maxdiff = 0.0
                !$omp parallel do default(shared) private(j,k,val) reduction(max:maxdiff)
                do j = 1, self%N
                    do k = 1, self%N
                        val = max( abs(self%A(j,k)-self%Aold(j,k)), abs(self%R(j,k)-self%Rold(j,k)) )
                        if (val > maxdiff) maxdiff = val
                    end do
                end do
                !$omp end parallel do
                if (maxdiff < self%ftol) then
                    write(logfhandle,'(a,i6,1x,es12.4)') 'aff_prop converged at iter=', i, maxdiff
                    converged = .true.
                    exit
                end if
            endif
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

end module simple_aff_prop
