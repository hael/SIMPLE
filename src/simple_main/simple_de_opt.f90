module simple_de_opt
use simple_optimizer, only: optimizer
use simple_opt_spec,  only: opt_spec
use simple_rnd,       only: ran3, irnd_uni
use simple_defs       ! singleton
implicit none

public :: de_opt
private

integer, parameter :: N_general = 136,    N_valley = 126,    N_multimodal = 103,    N_flat = 106
real, parameter    :: F_general = 0.2790, F_valley = 0.4027, F_multimodal = 0.3976, F_flat = 0.5860
real, parameter    :: X_general = 0.9813, X_valley = 0.9211, X_multimodal = 0.9794, X_flat = 0.3345
logical, parameter :: debug = .false.

type, extends(optimizer) :: de_opt
    private
    real, allocatable :: pop(:,:)          !< solution vector population
    real, allocatable :: costs(:)          !< costs
    integer           :: best              !< index of best population member
    real              :: F  = F_multimodal !< amplification factor
    real              :: CR = X_multimodal !< cross-over rate
    logical           :: exists = .false.  !< to indicate existence
  contains
    procedure :: new      => new_de
    procedure :: minimize => de_minimize
    procedure :: kill     => kill_de
end type de_opt

contains
    
    !> \brief  is a constructor
    subroutine new_de( self, spec )
        use simple_jiffys, only: alloc_err
        class(de_opt),   intent(inout) :: self !< instance
        class(opt_spec), intent(inout) :: spec !< specification
        integer :: alloc_stat
        ! destruct if exists
        call self%kill
        ! adjust control parameters according to mode
        select case(trim(spec%str_mode))
            case('general')
                spec%npop = N_general
                self%F    = F_general
                self%CR   = X_general
            case('valley')
                spec%npop = N_valley
                self%F    = F_valley
                self%CR   = X_valley
            case('multimodal')
                spec%npop = N_multimodal
                self%F    = F_multimodal
                self%CR   = X_multimodal
            case('flat')
                spec%npop = N_flat
                self%F    = F_flat
                self%CR   = X_flat
            case DEFAULT
                spec%npop = N_multimodal
                self%F    = F_multimodal
                self%CR   = X_multimodal
        end select
        ! allocate
        allocate(self%pop(spec%npop,spec%ndim), self%costs(spec%npop), stat=alloc_stat)
        call alloc_err("In: new_de; simple_de_opt", alloc_stat)
        self%exists = .true. ! indicates existence
        if( spec%debug ) write(*,*) 'created new differential evolution population'
    end subroutine new_de
    
    !> \brief  is the particle swarm minimize minimization routine
    subroutine de_minimize( self, spec, lowest_cost )
        class(de_opt), intent(inout)   :: self        !< instance     
        class(opt_spec), intent(inout) :: spec        !< specification
        real, intent(out)              :: lowest_cost !< lowest cost
        integer :: t                ! iteration counter
        integer :: npeaks           ! number of local minima
        real    :: rtol             ! relative tolerance
        integer :: loc(1), nworse, X
        if( .not. associated(spec%costfun) )then
            stop 'cost function not associated in opt_spec; de_minimize; simple_de_opt'
        endif
        ! initialize
        call init
        spec%nevals = 0
        npeaks      = 0
        nworse      = 0
        do t=1,spec%maxits ! generations loop
            ! select solution to modify
            X  = irnd_uni(spec%npop)
            call update_agent( X )
            if( spec%npeaks > 0 )then
                ! exit when we have identified spec%npeaks local optima
                if( npeaks == spec%npeaks ) exit
            endif
            if( nworse == spec%npop ) exit ! if no better solutions have been found in the last 4*ndim iterations
        end do
        lowest_cost = self%costs(self%best)
        spec%x = self%pop(self%best,:)

      contains

        !> \brief  initialize the population & set best
        subroutine init
            integer :: i, j
            real :: L
            ! obtain initial solutions by randomized bounds
            do i=1,spec%npop
                do j=1,spec%ndim
                    L = spec%limits(j,2)-spec%limits(j,1)
                    self%pop(i,j) = spec%limits(j,1)+ran3()*L
                end do
                if( i == 1 )then
                    if( .not. all(spec%x == 0.) )then ! set one point in the swarm to best point in spec
                        self%pop(1,:)= spec%x
                    endif
                endif
            end do
            ! calculate initial costs
            do i=1,spec%npop
                self%costs(i) = spec%costfun(self%pop(i,:), spec%ndim)
                spec%nevals = spec%nevals+1
            end do
            loc       = minloc(self%costs)
            self%best = loc(1)
        end subroutine init
        
        subroutine update_agent( X )
            integer, intent(in) :: X
            integer :: a, rb, b, i
            real :: trial(spec%ndim), cost_trial, L
            ! select random disjoint pair
            a  = irnd_uni(spec%npop)
            rb = irnd_uni(spec%npop-1)
            b  = a+rb 
            if( b <= spec%npop )then
            else
                b = a+rb-spec%npop
            endif
            ! create a trial solution
            do i=1,spec%ndim
                if( i == X .or. ran3() < self%CR )then
                    trial(i) = self%pop(self%best,i)+self%F*(self%pop(a,i)-self%pop(b,i))
                    ! enforce limits 
                    trial(i) = min(spec%limits(i,2),trial(i))
                    trial(i) = max(spec%limits(i,1),trial(i))
                else
                    trial(i) = self%pop(X,i) 
                endif
            end do
            ! calculate cost 
            cost_trial = spec%costfun(trial, spec%ndim)
            spec%nevals = spec%nevals+1
            ! update pop if better solution is found
            if( cost_trial < self%costs(X) )then
                nworse = 0
                self%pop(X,:) = trial
                self%costs(X) = cost_trial
                ! update global best if needed
                if( cost_trial < self%costs(self%best) ) self%best = X
                if( spec%npeaks > 0 )then
                    ! we will output local optima
                    if( cost_trial < spec%yb )then
                        npeaks = npeaks+1
                        spec%peaks(npeaks,:spec%ndim)  = trial
                        spec%peaks(npeaks,spec%ndim+1) = cost_trial
                    endif
                endif
            else
                nworse = nworse+1
            endif
        end subroutine update_agent
        
    end subroutine de_minimize

    !> \brief  is a destructor
    subroutine kill_de( self )
        class(de_opt), intent(inout) :: self !< instance
        integer :: i
        if( self%exists )then
            deallocate(self%pop,self%costs)
            self%exists = .false.
        endif
    end subroutine kill_de

end module simple_de_opt
