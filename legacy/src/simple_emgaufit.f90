!==Class simple_emgaufit
!
! simple_emgaufit is a SIMPLE class for fitting K probability distributions to one-
! dimensional data. This is just to test my knowledge about the EM algorithm.,
! The code is distributed with the hope that it will be useful, but _WITHOUT_ _ANY_ _WARRANTY_. 
! Redistribution or modification is regulated by the GNU General Public License. 
! *Author:* Hans Elmlund, 2012-02-02.
!
!==Changes are documented below
!
module simple_emgaufit
use simple_defs  ! singleton
implicit none

public :: emgaufit, test_emgaufit
private

type emgaufit
    private
    integer           :: N    ! number of data
    integer           :: K    ! number of components
    real, allocatable :: avgs(:), vars(:), mix(:), gammas(:,:)
    real, pointer     :: dat(:)=>null()
    logical           :: exists=.false.
  contains
    ! CONSTRUCTOR
    procedure :: new
    ! GETTERS/SETTERS
    procedure :: set_data
    procedure :: get_avgs
    procedure :: get_vars
    ! MACHINERY
    procedure :: fit
    procedure, private :: init
    procedure, private :: logL
    procedure, private :: estep
    procedure, private :: mstep
    ! DESTRUCTOR
    procedure :: kill
end type

contains

    ! CONSTRUCTOR

    !>  \brief  is a constructor
    subroutine new( self, N, K )
        use simple_jiffys, only: alloc_err
        class(emgaufit), intent(inout) :: self
        integer, intent(in)            :: N, K
        integer :: alloc_stat
        if( self%exists ) call self%kill
        self%N = N
        self%K = K
        allocate( self%avgs(K), self%vars(K), self%mix(K), self%gammas(N,K), stat=alloc_stat )
        call alloc_err('In: new; simple_emgaufit', alloc_stat )
    end subroutine
    
    ! GETTERS/SETTERS
    
    !>  \brief  is for setting pointer to data
    subroutine set_data( self, dat )
        class(emgaufit), intent(inout) :: self
        real, intent(in), target       :: dat(:)
        self%dat => dat
    end subroutine
    
    !>  \brief  is for getting the averages
    function get_avgs( self ) result( avgs )
        class(emgaufit), intent(in) :: self
        real :: avgs(self%K)
        avgs = self%avgs
    end function
    
    !>  \brief  is for getting the variances
    function get_vars( self ) result( vars )
        class(emgaufit), intent(in) :: self
        real :: vars(self%K)
        vars = self%vars
    end function
    
    ! MACHINERY
    
    !>  \brief  is the master fitting subroutine
    subroutine fit( self, mits )
        class(emgaufit), intent(inout) :: self
        integer, intent(in) :: mits
        integer             :: k
        real                :: L
        ! initialize
        call self%init
        ! calculate initial log-likelihood
        L = self%logL()
        print *, 'initial log(L):', L
        ! iterate
        do k=1,mits
            call self%estep
            call self%mstep
            L = self%logL()
            if( k == 1 .or. mod(k,5) == 0 )then
                write(*,"(1X,A,1X,I3,1X,A,1X,F10.0)") 'Iteration:', k, 'Log-likelihood:', L
            endif
        end do
    end subroutine
    
    !>  \brief  is for random initialization of the gammas
    subroutine init( self )
        use simple_rnd, only: ran3
        class(emgaufit), intent(inout) :: self
        integer :: i, j
        real    :: P(self%K), S
        do i=1,self%N
            S = 0.
            do j=1,self%K
                P(j) = ran3()
                S = S+P(j)
            end do
            do j=1,self%K
                self%gammas(i,j) = P(j)/S
            end do
        end do
        call self%mstep
    end subroutine
    
    !>  \brief  calculates the log likelihood
    function logL( self ) result ( L )
        use simple_math, only: gaussian
        class(emgaufit), intent(inout) :: self
        integer :: i, j
        real    :: L, S
        L = 0.
        do i=1,self%N
            S = 0.
            do j=1,self%K
                S = S+self%mix(j)*gaussian(self%dat(i),self%avgs(j),sqrt(self%vars(j)))
            end do
            L = L+log(S)
        end do    
    end function
    
    !>  \brief  evaluates the responsibilities (gammas) using the current parameter values
    subroutine estep( self )
        use simple_math, only: gaussian
        class(emgaufit), intent(inout) :: self
        real    :: P(self%K), S
        integer :: i, j
        do i=1,self%N
            S = 0.
            do j=1,self%K
                P(j) = self%mix(j)*gaussian(self%dat(i),self%avgs(j),sqrt(self%vars(j)))
                S = S+P(j)
            end do
            do j=1,self%K
                self%gammas(i,j) = P(j)/S
            end do
        end do
    end subroutine
    
    !>  \brief  updates the statistics based on the new responsibilities
    subroutine mstep( self )
        class(emgaufit), intent(inout) :: self
        real    :: NK, dev
        integer :: i, j
        do j=1,self%K
            ! calculate average
            self%avgs(j) = 0.
            NK = 0.
            do i=1,self%N
                self%avgs(j) = self%avgs(j)+self%dat(i)*self%gammas(i,j)
                NK = NK+self%gammas(i,j)
            end do
            self%avgs(j) = self%avgs(j)/NK
            ! calculate variance
            self%vars(j) = 0.
            do i=1,self%N
                dev = self%dat(i)-self%avgs(j)
                self%vars(j) = self%vars(j)+self%gammas(i,j)*dev*dev
            end do
            self%vars(j) = self%vars(j)/NK
            ! calculate mixing coefficients
            self%mix(j) = NK/real(self%N)
        end do
    end subroutine
    
    ! DESTRUCTOR
    
    !>  \brief  is a destructor
    subroutine kill( self )
        class(emgaufit), intent(inout) :: self
        if( self%exists )then
            deallocate(self%avgs, self%vars, self%mix, self%gammas)
            self%exists = .false.
        endif
    end subroutine
    
    ! UNIT TEST
    
    !>  \brief  is the unit test for the class
    subroutine test_emgaufit
        use simple_rnd,  only: gasdev
        use simple_stat, only: moment
        real           :: mydata(1000), avgs(2), vars(2)
        real           :: ave, sdev, var
        integer        :: i
        logical        :: err
        type(emgaufit) :: emfit
        ! generate the data
        do i=1,1000
            if( i <= 500 )then
                mydata(i) = gasdev(1.,1.)
            else
                mydata(i) = gasdev(5.,2.)
            endif
        end do
        call moment(mydata(:500), ave, sdev, var, err)
        write(*,*) 'FACIT AVG/VAR 1:', ave, var
        call moment(mydata(501:), ave, sdev, var, err)
        write(*,*) 'FACIT AVG/VAR 2:', ave, var
        ! fit
        call emfit%new(1000,2)
        call emfit%set_data(mydata)
        call emfit%fit(50)
        avgs = emfit%get_avgs()
        vars = emfit%get_vars()
        write(*,*) 'AVG/VAR 1:', avgs(1), vars(1)
        write(*,*) 'AVG/VAR 2:', avgs(2), vars(2)
    end subroutine
    
end module simple_emgaufit
