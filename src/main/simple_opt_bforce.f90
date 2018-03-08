! brute force function minimisation
module simple_opt_bforce
include 'simple_lib.f08'
use simple_optimizer, only: optimizer
implicit none

public :: opt_bforce
private
#include "simple_local_flags.inc"

type, extends(optimizer) :: opt_bforce
    private
    real, allocatable :: pb(:)          !< best point
    real, allocatable :: pc(:)          !< current point
    real              :: yb=0.          !< best cost function value
    logical           :: exists=.false. !< to indicate existence
  contains
    procedure :: new      => new_opt_bforce
    procedure :: minimize => bforce_minimize
    procedure :: kill     => kill_opt_bforce
end type

contains

    !> \brief  is a constructor
    subroutine new_opt_bforce( self, spec )
        use simple_opt_spec, only: opt_spec
        class(opt_bforce), intent(inout) :: self !< instance
        class(opt_spec), intent(inout)   :: spec !< specification
        integer                          :: i
        real                             :: x
        call self%kill
        allocate(self%pb(spec%ndim), self%pc(spec%ndim), stat=alloc_stat)
        if(alloc_stat.ne.0)call allocchk("In: new_opt_bforce", alloc_stat)
        self%pb = spec%limits(:,1)
        self%pc = spec%limits(:,1)
        if( all(spec%stepsz == 0.) ) stop 'step size (stepsz) not set in&
        &specification (opt_spec); new_opt_bforce; simple_opt_bforce'
        ! initialize best cost to huge number
        self%yb = huge(x)
        self%exists = .true. ! indicates existence
        DebugPrint  'created new opt_bforce object'
    end subroutine new_opt_bforce
    
    !> \brief  brute force minimization
    subroutine bforce_minimize( self, spec, fun_self, lowest_cost )
        use simple_opt_spec, only: opt_spec
        class(opt_bforce), intent(inout) :: self        !< instance
        class(opt_spec),   intent(inout) :: spec        !< specification
        class(*),          intent(inout) :: fun_self    !< self-pointer for cost function
        real, intent(out)                :: lowest_cost !< lowest cost
        real :: y
        if( .not. associated(spec%costfun) )then
            stop 'cost function not associated in opt_spec; bforce_minimize; simple_opt_bforce'
        endif
        ! init nevals counter
        spec%nevals = 0
        ! generate initial vector (lower bounds)
        spec%x = spec%limits(:,1)
        DebugPrint  'generated initial vector'
        ! set best and current point to best point in spec
        self%pb = spec%x
        self%pc = spec%x
        DebugPrint  'did set best and current point'
        ! set best cost
        self%yb     = spec%costfun(fun_self, self%pb, spec%ndim)
        if( debug ) write(*,'(a,1x,f7.3)') 'Initial cost:', self%yb 
        spec%nevals = spec%nevals+1
        ! search: we will start at the lowest value for each dimension, then 
        ! go in steps of stepsz until we get to the upper bounds
        spec%niter = 0
        DebugPrint  'starting brute force search'
        do while( srch_not_done() )
            y = spec%costfun(fun_self, self%pc, spec%ndim)
            spec%nevals = spec%nevals+1
            spec%niter  = spec%niter+1
            if( y <= self%yb )then
                self%yb = y       ! updating the best cost 
                self%pb = self%pc ! updating the best solution
                DebugPrint  'Found better best, cost:', self%yb 
            endif
        end do
        spec%x = self%pb
        lowest_cost = self%yb
        
        contains
            
            function srch_not_done() result( snd )
                integer :: i
                logical :: snd
                snd = .false.
                do i=1,spec%ndim
                    ! if we are still below the upper bound, increment this dim and 
                    ! set all other dims to the starting point (lower bound)
                    if( self%pc(i) < spec%limits(i,2) )then
                        ! if we got here, the search is not over
                        snd = .true.
                        ! increment the ith dimension
                        self%pc(i) = self%pc(i)+spec%stepsz(i)
                        ! reset all previous dimensions to the lower bound
                        if( i > 1 ) self%pc(1:i-1) = spec%limits(1:i-1,1)
                        ! if the ith dimension has reached or gone over its 
                        ! upper bound, set it to the upper bound 
                        self%pc(i) = min(self%pc(i),spec%limits(i,2))
                        exit
                    endif
                end do
                DebugPrint  'New configuration:', self%pc(:)
            end function
            
    end subroutine bforce_minimize

    !> \brief  is a destructor
    subroutine kill_opt_bforce( self )
        class(opt_bforce), intent(inout) :: self
        if( self%exists )then
            deallocate(self%pb, self%pc)
            self%exists = .false.
        endif
    end subroutine kill_opt_bforce
    
end module simple_opt_bforce
