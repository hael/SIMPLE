!> simple optimisation module: simulated annealing
!!
!! simple_sa_opt performs cobinatorial minimization of an externally defined
!! function by simulated annealing
! The code is distributed with the hope that it will be useful, but _WITHOUT_ _ANY_ _WARRANTY_.
! Redistribution or modification is regulated by the GNU General Public License.
! *Author:* Hans Elmlund, 2009-05-25.
!
!==Changes are documented below
!
!* incorporated in the _SIMPLE_ library, HE 2009-06-25
!
module simple_sa_opt
use simple_optimizer, only: optimizer
use simple_ran_tabu,  only: ran_tabu
use simple_jiffys,    only: alloc_err
implicit none

public :: sa_opt
private

type sa_opt
    private
    procedure(costfun),   pointer, nopass :: costfun =>null()   !< defines cost function
    procedure(acceptfun), pointer, nopass :: acceptfun =>null() !< defines acceptance function
    procedure(statefun),  pointer, nopass :: statefun=>null()   !< defines state function
    type(ran_tabu), allocatable :: rt(:)                        !< for constraint satisfaction
    integer, allocatable        :: vec(:,:)                     !< solution vector
    real, allocatable           :: cost(:)                      !< costs (one per label)
    integer, allocatable        :: params(:)                    !< nr of params per discrete variable
    logical, allocatable        :: nevars(:)                    !< indicates no two equal constraint
    integer                     :: N=0                          !< nr of particles
    integer                     :: L=0                          !< number of label series (dimensions)
    integer                     :: nprobs                       !< number of probabilities calculated per T level
    real                        :: avgprob                      !< average probability at one temperature level
    logical                     :: existence=.false.            !< indicates existence
  contains
    procedure          :: new
    procedure          :: set_statefun
    procedure          :: minimize
    procedure          :: minimize_climb
    procedure, private :: minimize_one_chain
    procedure, private :: minimize_one_chain_climb
    procedure          :: kill
end type sa_opt

!>  \brief  defines the function interfaces
abstract interface
    !> Cost minimisation/maximisation function of optimisation
    function costfun( vec, i, N, L ) result( cost )
        integer, intent(in) :: N, L, vec(N,L), i
        real                :: cost
    end function costfun
    !> Definition of completion function of optimisation
    subroutine acceptfun( vec, i, L )
        integer, intent(in) :: L, vec(L), i
    end subroutine acceptfun
    !< Current state of optimisation
    subroutine statefun( T, Tinit )
        real, intent(in) :: T, Tinit
    end subroutine
end interface

interface sa_opt
    module procedure constructor
end interface

contains

    !> \brief is a constructor
    !! \param costfun Cost minimisation/maximisation function
    !! \param acceptfun Definition of completion function
    function constructor( params, nevars, N, L, costfun, acceptfun  ) result( self )
        integer, intent(in) :: N         !< nr of individuals to optimize (nr of ptcls)
        integer, intent(in) :: L         !< vector length (problem dimensionality)
        integer, intent(in) :: params(L) !< array with nr of params per discrete variable
        logical, intent(in) :: nevars(L) !< indicates which variables that are not allowed to be equal to each other
        interface
            function costfun( vec, i, N, L ) result( cost )
                integer, intent(in) :: N, L, vec(N,L), i
                real                :: cost
            end function costfun
            subroutine acceptfun( vec, i, L )
                integer, intent(in) :: L, vec(L), i
            end subroutine
        end interface
        type(sa_opt)        :: self
        call self%new( params, nevars, N, L, costfun, acceptfun )
    end function

    !> \brief is a constructor
    !> \param costfun Cost minimisation/maximisation function
    !> \param acceptfun Definition of completion function
    subroutine new( self, params, nevars, N, L, costfun, acceptfun )
        class(sa_opt), intent(inout) :: self      !< object
        integer, intent(in)          :: N         !< nr of individuals to optimize (nr of ptcls)
        integer, intent(in)          :: L         !< vector length (problem dimensionality)
        integer, intent(in)          :: params(L) !< array with nr of params per discrete variable
        logical, intent(in)          :: nevars(L) !< variables that are not allowed to be equal
        interface
             !> Cost minimisation/maximisation function
            function costfun( vec, i, N, L ) result( cost )
                integer, intent(in) :: N, L, vec(N,L), i
                real                :: cost
            end function costfun
            !> Definition of completion function
            subroutine acceptfun( vec, i, L )
                integer, intent(in) :: L, vec(L), i
            end subroutine
        end interface
        integer                      :: alloc_stat, j
        call self%kill
        allocate(self%vec(N,L), self%cost(N), self%params(L),&
        self%nevars(L), self%rt(L), stat=alloc_stat)
        call alloc_err('new; simple_sa_opt', alloc_stat)
        self%params = params
        self%nevars = nevars
        self%N      = N
        self%L      = L
        do j=1,self%L
           if( self%nevars(j) ) self%rt(j) = ran_tabu(self%params(j))
        end do
        ! set function pointers
        self%costfun => costfun
        self%acceptfun => acceptfun
        self%existence= .true.
    end subroutine

    !> \brief  communicates statefun to the object
    subroutine set_statefun( self, statefun )
        class(sa_opt), intent(inout) :: self
        interface
            subroutine statefun( T, Tinit )
                real, intent(in) :: T, Tinit
            end subroutine
        end interface
        self%statefun => statefun
    end subroutine

    !> \brief is the combinatorial minimization. I have found the params defined by Kirkpatrick to work well:
    !!        T_init=1000000, TC=0.9, max_rearr=10**L*N, and T_lowlim=0.01. The parameter max_rearr (maximum
    !!        number of rearrangements per T_-level) affects performance and CPU time most, requires testing
    !!        for specific problems.
    subroutine minimize( self, cost_err, T_init, TC, max_rearr, T_lowlim, out_solution, cost )
        use simple_rnd,    only: irnd_uni
        use simple_jiffys, only: progress
        class(sa_opt), intent(inout)  :: self
        real, intent(in)              :: cost_err    !< Convergence factor, range of error
        real, intent(in)              :: T_init, TC, T_lowlim
        integer, intent(in)           :: max_rearr   !< maximum number of rearrangements per T_-level
        integer, intent(out)          :: out_solution(self%N,self%L) !< final output of params
        real, intent(out)             :: cost                    !< final cost value of solution
        real                          :: joint_cost, T, costs(2)
        integer                       :: it, convergence, mits, i, j
        logical                       :: frozen
        ! Initialization of solution
        do j=1,self%L
            if( self%nevars(j) ) call self%rt(j)%reset
            do i=1,self%N
                if( self%nevars(j) ) then
                    self%vec(i,j) = self%rt(j)%irnd()
                    call self%rt(j)%insert(self%vec(i,j))
                else
                    self%vec(i,j) = irnd_uni(self%params(j))
                endif
            end do
        end do
        ! communicate the init solution via the acceptor function
        do i=1,self%N
            call self%acceptfun( self%vec(i,:), i, self%L )
        end do
        ! Initialization of variables
        T           = T_init
        frozen      = .false. ! used to minitor convergence
        it          = 0       ! iteration counter
        convergence = 0       ! used to minitor convergence
        costs       = 99999.  ! used to minitor convergence
        mits        = nits()  ! maximum nr of iterations
        write(*,'(A)') '>>> ANNEALING_STARTS'
        do while (frozen .neqv. .true.)
            it = it+1
            call progress(it, mits)
            if( associated(self%statefun) ) then
                ! use the state function to perform temerature dependent external changes
                call self%statefun(T,T_init)
            endif
            self%nprobs = 0
            self%avgprob = 0.
            call self%minimize_one_chain(T, max_rearr, joint_cost)
            self%avgprob = self%avgprob/real(self%nprobs)
!            write(*,'(A)') '***************************'
!            write(*,'(A,14XF10.2)') 'T: ', T
!            write(*,'(A,15XF6.2)') 'P(T): ', self%avgprob
!            write(*,'(A,13XF8.4)') 'COST: ', joint_cost
            T           = TC*T     ! Temperature update
            costs(1)    = costs(2) ! Cost update
            costs(2)    = joint_cost
            convergence = 0
            if( abs(costs(1)-costs(2)) <= cost_err*0.5 ) convergence = convergence+1
            ! check if system frozen
            frozen = .false.
            if( T <= T_lowlim .or. convergence >= 10 ) frozen = .true.
            if( frozen ) then
!                write(*,'(A)') 'System frozen'
!                write(*,'(A,4XF10.2)') 'T: ', T
!                write(*,'(A,5XF6.2)') 'P(T): ', self%avgprob
!                write(*,'(A,3XF8.4)') 'COST: ', joint_cost
            end if
        end do
        out_solution = self%vec
        cost = joint_cost

        contains

            function nits( ) result( n )
                integer :: n
                real :: tmp
                tmp = T_init
                n = 0
                do while( tmp > T_lowlim )
                    n = n+1
                    tmp = tmp*TC
                end do
            end function

    end subroutine

    !> \brief  provides the core functionality of the simulated annealing
    !!         and parallel tempering algorithms
    subroutine minimize_one_chain( self, T, max_rearr, cost )
        use simple_rnd, only: irnd_uni
        class(sa_opt), intent(inout) :: self
        real, intent(in)             :: T
        integer, intent(in)          :: max_rearr  !< maximum number of rearrangements per T_-level
        real, intent(out)            :: cost
        integer                      :: i, j, nr_rearr, olds(self%L)
        real                         :: cost_perm
        ! initialize costs
        do i=1,self%N
            self%cost(i) = self%costfun( self%vec, i, self%N, self%L )
        end do
        nr_rearr  = 0 ! nr of rearrangements at given T
        do while(nr_rearr < max_rearr)
            do i=1,self%N
                ! randomly permute solution
                do j=1,self%L
                    olds(j) = self%vec(i,j)
                    if( self%nevars(j) ) then
                        self%vec(i,j) = self%rt(j)%irnd()
                    else
                        do while( self%vec(i,j) == olds(j) )
                            self%vec(i,j) = irnd_uni(self%params(j))
                        end do
                    endif
                end do
                ! calculate cost of new configuration
                cost_perm = self%costfun( self%vec, i, self%N, self%L )
                ! Metropolis-Hastings acceptance criterion
                call metropolis_accept
            end do
            ! count the number of rearrangements
            nr_rearr = nr_rearr+1
        end do
        cost = sum(self%cost)/real(self%N)

        contains

            !> \brief  is for temperature-weighted probabilistic acceptance of
            !!         trial solutions according to the Metropolis-Hastings criterion
            subroutine metropolis_accept
                use simple_rnd, only: ran3
                real            :: prob
                logical         :: perm_accepted
                integer         :: j
                real, parameter :: bolzmann=10000.
                perm_accepted = .false.
                if( cost_perm < self%cost(i) ) then
                    ! always accept the downhill climb
                    perm_accepted = .true.
                else
                    prob         = min(1.,exp(((-bolzmann)*(cost_perm-self%cost(i)))/T))
                    self%nprobs  = self%nprobs+1
                    self%avgprob = self%avgprob+prob
                    if( ran3() <= prob ) then
                        ! sometimes accept an uphill climb
                        perm_accepted = .true.
                    endif
                endif
                if( perm_accepted )then
                    self%cost(i) = cost_perm
                    ! take care of the constraint satisfaction
                    do j=1,self%L
                        if( self%nevars(j) ) then
                            call self%rt(j)%remove(olds(j))
                            call self%rt(j)%insert(self%vec(i,j))
                        endif
                    end do
                    ! use the observer function to notify the change
                    call self%acceptfun(self%vec(i,:), i, self%L)
                else
                    do j=1,self%L
                        self%vec(i,j) = olds(j)
                    end do
                endif
            end subroutine

    end subroutine

    !> \brief combinatorial minimization by simple hill climbing
    subroutine minimize_climb( self, max_rearr, in_solution, out_solution, cost )
        class(sa_opt), intent(inout) :: self
        integer, intent(in)          :: max_rearr
        integer, intent(in)          :: in_solution(self%N,self%L)
        integer, intent(out)         :: out_solution(self%N,self%L)
        real, intent(out)            :: cost
        real                         :: joint_cost
        integer                      :: i, j
        ! Initialization based on input
        do j=1,self%L
            if( self%nevars(j) ) call self%rt(j)%reset
            do i=1,self%N
                if( self%nevars(j) ) then
                    self%vec(i,j) = in_solution(i,j)
                    call self%rt(j)%insert(self%vec(i,j))
                else
                    self%vec(i,j) = in_solution(i,j)
                endif
            end do
        end do
        ! communicate the init solution via the acceptor function
        do i=1,self%N
            call self%acceptfun( self%vec(i,:), i, self%L )
        end do
        call self%minimize_one_chain_climb(max_rearr, joint_cost)
        out_solution = self%vec
        cost = joint_cost
    end subroutine

    !> \brief  provides local search functionality
    subroutine minimize_one_chain_climb( self, max_rearr, cost )
        class(sa_opt), intent(inout) :: self
        integer, intent(in)          :: max_rearr  !< maximum rearrangements
        real, intent(out)            :: cost       !< final cost value
        integer                      :: i, j, k, nr_rearr, accepted, no_change, olds(self%L), alloc_stat
        logical                      :: steady
        real, allocatable            :: costs(:)
        real                         :: x
        write(*,'(A)') '>>> HILL_CLIMBING_STARTS'
        ! initialize costs
        do i=1,self%N
            self%cost(i) = self%costfun( self%vec, i, self%N, self%L )
        end do
        nr_rearr  = 0       ! nr of maximum rearrangements
        steady    = .false. ! to monitor steady state
        no_change = 0       ! nr of no changes
        do while(steady .neqv. .true.)
            accepted  = 0  ! nr of accepted permutations
            do i=1,self%N
                ! evaluate the entire available neighborhood
                do j=1,self%L
                    olds(j) = self%vec(i,j)
                    allocate( costs(self%params(j)), stat=alloc_stat )
                    call alloc_err( 'minimize_one_chain_climb; simple_sa_opt', alloc_stat )
                    do k=1,self%params(j)
                        self%vec(i,j) = k
                        if( self%nevars(j) ) then
                            if( self%rt(j)%is(k) )then ! element is tabu
                                 costs(k) = huge(x)
                            else
                                ! calculate cost of new configuration
                                costs(k) = self%costfun( self%vec, i, self%N, self%L )
                            endif
                        else
                            ! calculate cost of new configuration
                            costs(k) = self%costfun( self%vec, i, self%N, self%L )
                        endif
                    end do
                    call climb_accept
                    deallocate( costs )
                end do
            end do
            ! count the changes
            if( accepted == 0 ) then
                no_change = no_change+1
            else
                no_change = 0
            end if
            ! count the number of rearrangements
            nr_rearr = nr_rearr+1
            ! check if steady state
            steady = .false.
            if( no_change >= 2 .or. nr_rearr >= max_rearr ) steady = .true.
            write(*,'(A)') '***************************'
            write(*,'(A,13X,F8.4)') 'COST: ', sum(self%cost)/real(self%N)
        end do
        cost = sum(self%cost)/real(self%N)

        contains

            subroutine climb_accept
                integer :: loc(1)
                loc = minloc(costs)
                ! always accept a new downhill climb
                if( costs(loc(1)) < self%cost(i) .and. olds(j) /= loc(1) )then
                    accepted = accepted+1
                    self%vec(i,j) = loc(1)
                    self%cost(i)  = costs(loc(1))
                    ! take care of the constraint satisfaction
                    if( self%nevars(j) ) then
                        call self%rt(j)%remove(olds(j))
                        call self%rt(j)%insert(self%vec(i,j))
                    endif
                    ! use the observer function to notify the change
                    call self%acceptfun(self%vec(i,:), i, self%L)
                else
                    self%vec(i,j) = olds(j)
                endif
            end subroutine

    end subroutine

    !> \brief is a destructor
    subroutine kill( self )
        class(sa_opt), intent(inout) :: self
        integer                      :: j
        if( self%existence )then
            do j=1,self%L
                if( self%nevars(j) ) then
                    call self%rt(j)%kill
                endif
            end do
            deallocate( self%vec, self%cost, self%params, self%nevars, self%rt )
            self%existence= .false.
        endif
    end subroutine

end module simple_sa_opt
