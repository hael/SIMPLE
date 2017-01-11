!==Class simple_aff_prop
!
! simple_aff_prop is the SIMPLE class for clustering by affinity propagation.
! Affinity propagation is based on loopy belief propagation.
! The code is distributed with the hope that it will be useful, but _WITHOUT_ 
! _ANY_ _WARRANTY_. Redistribution or modification is regulated by the 
! GNU General Public License. *Author:* Hans Elmlund, 2011-09-03.
!
!==Changes are documented below
!
module simple_aff_prop
use simple_math
use simple_rnd
use simple_stat
implicit none

public :: aff_prop, test_aff_prop
private

type aff_prop
    private
    integer              :: N                    !< nr of data entries
    real, allocatable    :: A(:,:), R(:,:)       !< affinities & responsibilities
    real, allocatable    :: Aold(:,:), Rold(:,:) !< old affinities & responsibilities
    real, allocatable    :: Rp(:,:), tmp(:)      !< other stufff needed
    real, allocatable    :: AS(:,:), dA(:)       !< A+S & diag(A)
    real, pointer        :: S(:,:)               !< pointer to similarity matrix
    real, allocatable    :: Y(:), Y2(:)          !< maxvals
    integer, allocatable :: I(:), I2(:)          !< index arrays
    integer              :: maxits=200           !< maximum number of iterations
    integer              :: convits=15           !< convergence factor
    real                 :: ftol=1e-9            !< fractional convergence tolerance score
    real                 :: lam=0.5              !< dampening factor
    real                 :: Smin, Smax           !< similarity characteristics
    logical              :: exists=.false.       !< to indicate existence
  contains
    ! CONSTRUCTOR
    procedure :: new
    ! PROPAGATOR
    procedure :: propagate
    ! GETTERS
    procedure :: get_resp
    ! DESTRUCTOR
    procedure :: kill
end type

contains

    ! CONSTRUCTOR

    !>  \brief  is a constructor
    subroutine new( self, N, S, ftol, lam, pref, maxits, convits )
        class(aff_prop),   intent(inout) :: self
        integer,           intent(in)    :: N
        real,    target,   intent(inout) :: S(N,N)
        real,    optional, intent(in)    :: ftol, lam, pref
        integer, optional, intent(in)    :: maxits, convits
        integer :: alloc_stat, i, j
        real    :: ppref, diff
        call self%kill
        ! set constants
        self%S => S
        self%N = N
        if( present(ftol) )    self%ftol    = ftol
        if( present(lam) )     self%lam     = lam
        if( present(maxits) )  self%maxits  = maxits
        if( present(convits) ) self%convits = convits
        ! calculate similarity characteristics
        call analyze_smat(s, .false., self%Smin, self%Smax)
        ! low pref [0.1,1] leads to small numbers of clusters and 
        ! high pref leads to large numbers of clusters
        ppref = self%Smin
        if( present(pref) ) ppref = pref
        ! remove degeneracies
        forall(i=1:N) self%S(i,i) = ppref
        diff = (self%Smax-self%Smin)*self%ftol
        do i=1,self%N
            do j=1,self%N
                self%S(i,j) = self%S(i,j)+ran3()*diff
            end do
        end do
        ! allocate
        allocate( self%A(N,N), self%R(N,N), self%Aold(N,N), self%Rp(N,N),&
        self%Rold(N,N), self%AS(N,N), self%Y(N), self%Y2(N), self%tmp(N),&
        self%I(N), self%I2(N), self%dA(N), stat=alloc_stat )
        call alloc_err('new; simple_aff_prop, 2', alloc_stat)
        ! initialize
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
        integer, allocatable, intent(out)   :: centers(:), labels(:)
        real,                 intent(out)   :: simsum
        real, allocatable :: similarities(:)
        real              :: x, realmax
        integer           :: i, j, k, alloc_stat, ncls, loc(1), convits, ncls_prev
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
        convits   = 0
        ncls      = 0
        ! iterate        
        do i=1,self%maxits
            ! FIRST, COMPUTE THE RESPONSIBILITIES
            self%Rold = self%R
            self%AS   = self%A+self%S
            self%I    = maxloc(self%AS, dim=2)
            forall(j=1:self%N) self%Y(j) = self%AS(j,self%I(j))
            do j=1,self%N
                self%AS(j,self%I(j)) = -realmax
            end do
            self%I2 = maxloc(self%AS, dim=2)
            forall(j=1:self%N) self%Y2(j) = self%AS(j,self%I2(j))
            self%R = self%S
            forall(j=1:self%N) self%R(j,:) = self%R(j,:)-self%Y
            do j=1,self%N
               self%R(j,self%I(j)) = self%S(j,self%I(j))-self%Y2(j)
            end do
            ! update responsibilities (in a dampened fashion)
            self%R = (1.-self%lam)*self%R+self%lam*self%Rold
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
            forall(j=1:self%N) self%A(j,:)  = self%A(j,:)+self%tmp
            forall(k=1:self%N) self%dA(k)   = self%A(k,k)
            where(self%A > 0.) self%A       = 0.
            forall(k=1:self%N) self%A(k,k)  = self%dA(k)
            ! update availabilities (in a dampened fashion)
            self%A = (1.-self%lam)*self%A+self%lam*self%Aold
            ! count the number of clusters
            ncls_prev = ncls
            ncls = 0
            do j=1,self%N
                if( self%R(j,j)+self%A(j,j) > 0. )then
                    ncls = ncls+1
                endif
            end do
            if( ncls_prev == ncls )then
                convits = convits+1
            else
                convits = 0
            endif
            if( convits == self%convits ) exit
        end do
        self%R = self%R+self%A ! pseudomarginals
        if( ncls==0 )ncls=1
        if( allocated(centers) ) deallocate(centers)
        if( allocated(labels) )  deallocate(labels)
        allocate( centers(ncls), similarities(ncls), labels(self%N), stat=alloc_stat )
        call alloc_err("In: propagate; simple_aff_prop", alloc_stat)
        ncls = 0
        do j=1,self%N
            if( self%R(j,j) > 0. )then
                ncls = ncls+1
                centers(ncls) = j
            endif
        end do
        ! report back the labeling
        simsum = 0.
        do j=1,self%N
            do k=1,ncls
                if( j.ne.centers(k) )then
                    similarities(k) = self%S(centers(k),j)
                else
                    similarities(k) = realmax
                endif
            end do
            loc = maxloc(similarities)
            labels(j) = loc(1)
            if( j.ne.centers(labels(j)) ) simsum = simsum+similarities(labels(j))
        end do
        simsum = simsum/real(self%N)
        deallocate(similarities)
    end subroutine propagate

    ! GETTERS

    function get_resp( self )result( resp )
        class(aff_prop), intent(inout) :: self
        real :: resp( self%N,self%N )
        resp = self%R
    end function get_resp

    ! UNIT TEST
    
    !>  \brief  is the aff_prop unit test
    subroutine test_aff_prop
        real                 :: datavecs(900,5)
        type(aff_prop)       :: apcls
        real                 :: simmat(900,900), simsum
        integer, allocatable :: centers(:), labels(:)
        integer              :: i, j, ncls, nerr
        write(*,'(a)') '**info(simple_aff_prop_unit_test): testing all functionality'
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
            do j=i+1,900
                simmat(i,j) = -euclid(datavecs(i,:),datavecs(j,:))
                simmat(j,i) = simmat(i,j)
            end do
        end do
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
        write(*,*) 'NR OF CLUSTERS FOUND:', ncls
        write(*,*) 'NR OF ASSIGNMENT ERRORS:', nerr
        write(*,*) 'CENTERS'
        do i=1,size(centers)
            write(*,*) datavecs(centers(i),:)
        end do
        if( ncls == 3 .and. nerr == 0 )then
            write(*,'(a)') 'SIMPLE_AFF_PROP_UNIT_TEST COMPLETED ;-)'
        else
            write(*,'(a)') 'SIMPLE_AFF_PROP_UNIT_TEST FAILED!'
        endif
    end subroutine test_aff_prop
    
    ! DESTRUCTOR
    
    subroutine kill( self )
        class(aff_prop), intent(inout) :: self
        if( self%exists )then
            deallocate( self%A, self%R, self%Aold, self%Rp,&
            self%Rold, self%AS, self%Y, self%Y2, self%tmp, self%I,&
            self%I2, self%dA )
            self%exists = .false.
        endif
    end subroutine kill

end module simple_aff_prop
