module simple_nanoML
use simple_defs
use simple_math
implicit none

public :: nanoML, test_nanoML
private

type nanoML
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
    procedure :: set_avgs
    procedure :: get_avgs
    procedure :: get_vars
    procedure :: get_gammas
    ! MACHINERY
    procedure :: kmeans_biased2classes
    procedure :: fit
    procedure, private :: init, random_init
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
        class(nanoML), intent(inout) :: self
        integer, intent(in)          :: N, K
        if( self%exists ) call self%kill
        self%N = N
        self%K = K
        allocate( self%avgs(K), self%vars(K), self%mix(K), self%gammas(N,K) )
    end subroutine

    ! GETTERS/SETTERS

    !>  \brief  is for setting pointer to data
    subroutine set_data( self, dat )
        class(nanoML), intent(inout) :: self
        real, intent(in), target     :: dat(:)
        self%dat => dat
    end subroutine

    !>  \brief  is for setting the averages
    subroutine set_avgs( self, avgs)
        class(nanoML), intent(inout) :: self
        real,          intent(in)    :: avgs(self%K)
        self%avgs = avgs
    end subroutine

    !>  \brief  is for getting the averages
    function get_avgs( self ) result( avgs )
        class(nanoML), intent(in) :: self
        real :: avgs(self%K)
        avgs = self%avgs
    end function

    !>  \brief  is for getting the variances
    function get_vars( self ) result( vars )
        class(nanoML), intent(in) :: self
        real :: vars(self%K)
        vars = self%vars
    end function

    !>  \brief  is for getting the gammas TO REMOVEE
    function get_gammas( self ) result( gammas )
        class(nanoML), intent(in) :: self
        real :: gammas(self%N,self%K)
        gammas = self%gammas
    end function

    ! MACHINERY

    !>  \brief  is the master fitting subroutine
    subroutine fit( self, mits, avgs)
        class(nanoML), intent(inout) :: self
        integer,         intent(in)  :: mits
        real   ,         intent(in)  :: avgs(self%K)
        integer :: k
        real    :: L
        write(logfhandle,*) '****ML with k-means fit, init'
        ! initialize
        call self%init(avgs)
        ! call self%random_init()
        ! calculate initial log-likelihood
        L = self%logL()
        ! iterate
        do k=1,mits
            call self%estep
            call self%mstep
            L = self%logL()
            if( k == 1 .or. mod(k,5) == 0 )then
                write(*,"(1X,A,1X,I3,1X,A,1X,F10.0)") 'Iteration:', k, 'Log-likelihood:', L
            endif
        end do
        write(logfhandle,*) '****ML with k-means fit, completed'
    end subroutine

    !>  \brief  is for random initialization of the gammas
    subroutine random_init( self )
        use simple_rnd, only: ran3
        class(nanoML), intent(inout) :: self
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

    !>  \brief  is for initialization with centers
    !   identified by kmeans
    !   mixing coeffs are initialised uniformly
    !   variances are initialised randomly
    !   gammas are consequentely calculated
    subroutine init( self, avgs )
        use simple_rnd, only: ran3
        class(nanoML), intent(inout) :: self
        real,            intent(in)  :: avgs(self%K) ! input avgs are kmeans output
        integer :: i, j
        real    :: P(self%K), S
        ! set averages
        call self%set_avgs(avgs)
        ! set mixing coeffs
        self%mix = 1./real(self%K)
        ! set variances
        do i = 1, self%K
            self%vars(i) = 10.*ran3()
        enddo
        ! set gammas
        call self%estep
    end subroutine

    !>  \brief  calculates the log likelihood
    function logL( self ) result ( L )
        class(nanoML), intent(inout) :: self
        integer :: i, j
        real    :: L, S
        L = 0.
        do i=1,self%N
            S = 0.
            do j=1,self%K
                S = S+self%mix(j)*gaussian1D(self%dat(i),self%avgs(j),sqrt(self%vars(j)))
            end do
            L = L+log(S)
        end do
    end function

    !>  \brief  evaluates the responsibilities (gammas) using the current parameter values
    subroutine estep( self )
        class(nanoML), intent(inout) :: self
        real    :: P(self%K), S
        integer :: i, j
        do i=1,self%N
            S = 0.
            do j=1,self%K
                P(j) = self%mix(j)*gaussian1D(self%dat(i),self%avgs(j),sqrt(self%vars(j)))
                S = S+P(j)
            end do
            do j=1,self%K
                self%gammas(i,j) = P(j)/S
            end do
        end do
    end subroutine

    !>  \brief  updates the statistics based on the new responsibilities
    subroutine mstep( self )
        class(nanoML), intent(inout) :: self
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

    ! This subroutine clusters the atoms with respect to the maximum intensity
    ! using kmean algorithm for 2 classes. The initial guess fo the centers
    ! is intentionally biased. It supposes there are two distinguished classes
    ! with different avgs (proved with simulated data).
    ! It returns the centers of the classes, which is going to be the
    ! input for ML algorithm
    subroutine kmeans_biased2classes(self,data, centers)
        class(nanoML), intent(inout) :: self
        real,          intent(inout) :: data(:)
        real,          intent(inout) :: centers(2) ! centers of the classes, output
        real    :: data_copy(size(data))
        integer :: i, cnt1, cnt2
        integer :: val1, val2
        logical :: converged
        integer :: N_dat
        N_dat = size(data)
        write(logfhandle,*) '****kmeans_biased2classes, init'
        data_copy = data
        ! Initialise
        call initialise_centers(data_copy,centers(1),centers(2))
        converged = .false.
        do i = 1, N_dat  ! maximum number of iterations = N_dat
            if(.not. converged) then
                call update_centers(centers(1),centers(2),converged,val1,val2)
            else
                exit
            endif
        enddo
        write(logfhandle,*) '****kmeans_biased2classes, completed'

    contains

        subroutine initialise_centers(data,cen1,cen2)
            real, intent(inout) :: cen1,cen2
            real, intent(inout) :: data(:)
            !>   rheapsort from numerical recepies (largest last)
            call hpsort(data)
            cen1 = sum(data(1:N_dat/2))/real(N_dat/2)
            cen2 = sum(data(N_dat/2+1:size(data)))/real(N_dat/2)
        end subroutine initialise_centers

        subroutine update_centers(cen1,cen2,converged,val1,val2)
            real,    intent(inout) :: cen1,cen2
            logical, intent(inout) :: converged
            integer, intent(inout) :: val1, val2
            integer :: i
            integer :: cnt1, cnt2
            real    :: sum1, sum2
            real :: cen1_new, cen2_new
            sum1 = 0.
            cnt1 = 0
            do i=1,N_dat
                if( (cen1-data(i))**2. < (cen2-data(i))**2. )then
                    cnt1 = cnt1 + 1 ! number of elements in cluster 1
                    sum1 = sum1 + data(i)
                endif
            end do
            cnt2 = N_dat - cnt1       ! number of elements in cluster 2
            sum2 = sum(data)- sum1
            cen1_new = sum1 / real(cnt1)
            cen2_new = sum2 / real(cnt2)
            if(abs(cen1_new - cen1) < TINY .and. abs(cen2_new - cen2) < TINY) then
                converged = .true.
            else
                converged = .false.
            endif
            ! update
            cen1 = cen1_new
            cen2 = cen2_new
            ! assign values to the centers
            if( cen1 > cen2 )then
                val1           = 1
                val2           = 0
            else
                val1           = 0
                val2           = 1
            endif
        end subroutine update_centers
    end subroutine kmeans_biased2classes

    ! DESTRUCTOR

    !>  \brief  is a destructor
    subroutine kill( self )
        class(nanoML), intent(inout) :: self
        if( self%exists )then
            deallocate(self%avgs, self%vars, self%mix, self%gammas)
            self%exists = .false.
        endif
    end subroutine

    ! UNIT TEST

    !>  \brief  is the unit test for the class
    subroutine test_nanoML
        use simple_rnd,  only: gasdev
        use simple_stat, only: moment
        real           :: mydata(1000), avgs(2), vars(2)
        real           :: ave, sdev, var
        integer        :: i
        logical        :: err
        type(nanoML)   :: emfit
        real           :: centers_kmeans(2)
        ! generate the data
        do i=1,1000
            if( i <= 500 )then
                mydata(i) = gasdev(1.,1.)
            else
                mydata(i) = gasdev(5.,2.)
            endif
        end do
        ! kmeans
        call emfit%kmeans_biased2classes(mydata, centers_kmeans)
        ! expected output
        call moment(mydata(:500), ave, sdev, var, err)
        write(*,*) 'FACIT AVG/VAR 1:', ave, var
        call moment(mydata(501:), ave, sdev, var, err)
        write(*,*) 'FACIT AVG/VAR 2:', ave, var
        ! ML, dfit
        call emfit%new(1000,2)
        call emfit%set_data(mydata)
        call emfit%fit(50,centers_kmeans) !50: max iterations
        avgs = emfit%get_avgs()
        vars = emfit%get_vars()
        write(*,*) 'AVG/VAR 1:', avgs(1), vars(1)
        write(*,*) 'AVG/VAR 2:', avgs(2), vars(2)
    end subroutine
end module simple_nanoML
