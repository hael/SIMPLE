!==Class simple_bforce_opt
!!
!! This class is based on a class used in CTFFIND4, developed by Alexis Rohou
!! and Nikolaus Grigorieff at Janelia Farm. The below copyright statement therefore 
!! needs to be included here:
!! Copyright 2014 Howard Hughes Medical Institute
!! All rights reserved
!! Use is subject to Janelia Farm Research Campus Software Copyright 1.1
!! license terms ( http://license.janelia.org/license/jfrc_copyright_1_1.html )
!!
module simple_bforce_opt
use simple_optimizer, only: optimizer
implicit none

public :: bforce_opt
private

logical :: debug=.false.

type, extends(optimizer) :: bforce_opt
    private
    real, allocatable    :: pb(:)          !< best point
    real, allocatable    :: pc(:)          !< current point
    real                 :: yb=0.          !< best cost function value
    logical              :: exists=.false. !< to indicate existence
  contains
    procedure :: new          => new_bforce_opt
    procedure :: minimize     => bforce_minimize
    procedure :: get_vertices => bforce_get_vertices
    procedure :: kill         => kill_bforce_opt
end type

contains

    !> \brief  is a constructor
    subroutine new_bforce_opt( self, spec )
        use simple_opt_spec, only: opt_spec
        use simple_jiffys,   only: alloc_err
        class(bforce_opt), intent(inout) :: self !< instance
        class(opt_spec), intent(inout)   :: spec !< specification
        integer                          :: alloc_stat, i
        real                             :: x
        call self%kill
        allocate(self%pb(spec%ndim), self%pc(spec%ndim), stat=alloc_stat)
        call alloc_err("In: new_bforce_opt", alloc_stat)
        self%pb = spec%limits(:,1)
        self%pc = spec%limits(:,1)
        if( all(spec%stepsz == 0.) ) stop 'step size (stepsz) not set in&
        specification (opt_spec); new_bforce_opt; simple_bforce_opt'
        ! initialize best cost to huge number
        self%yb = huge(x)
        self%exists = .true. ! indicates existence
        if( debug ) write(*,'(a)') 'created new bforce_opt object'
    end subroutine
    
    !> \brief  brute force minimization
    subroutine bforce_minimize( self, spec, lowest_cost )
        use simple_opt_spec, only: opt_spec
        class(bforce_opt), intent(inout) :: self        !< instance
        class(opt_spec), intent(inout)   :: spec        !< specification
        real, intent(out)                :: lowest_cost !< lowest cost
        real :: y
        if( .not. associated(spec%costfun) )then
            stop 'cost function not associated in opt_spec; bforce_minimize; simple_bforce_opt'
        endif
        ! generate initial vector (lower bounds)
        spec%x = spec%limits(:,1)
        if( debug ) write(*,'(a)') 'generated initial vector'
        ! set best & current point to best point in spec
        self%pb = spec%x
        self%pc = spec%x
        if( debug ) write(*,'(a)') 'did set best & current point'
        ! set best cost
        spec%nevals = 0
        self%yb     = spec%costfun(self%pb, spec%ndim)
        if( debug ) write(*,'(a,1x,f7.3)') 'Initial cost:', self%yb 
        spec%nevals = spec%nevals+1
        ! search: we will start at the lowest value for each dimension, then 
        ! go in steps of stepsz until we get to the upper bounds
        spec%niter = 0
        if( debug ) write(*,'(a)') 'starting brute force search'
        do while( srch_not_done() )
            y = spec%costfun(self%pc, spec%ndim)
            spec%nevals = spec%nevals+1
            spec%niter  = spec%niter+1
            if( y <= self%yb )then
                self%yb = y       ! updating the best cost 
                self%pb = self%pc ! updating the best solution
                if( debug ) write(*,*) 'Found better best, cost:', self%yb 
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
                if( debug ) write(*,*) 'New configuration:', self%pc(:)
            end function
            
    end subroutine

    !> \brief  dummy procedure, only defined in simplex
    subroutine bforce_get_vertices( self, spec, vertices, costs )
        use simple_opt_spec,           only: opt_spec
        class(bforce_opt), intent(inout) :: self
        class(opt_spec),    intent(inout) :: spec
        real, allocatable,  intent(inout) :: vertices(:,:), costs(:)
        stop 'procedure only defined for simplex optimisation'
    end subroutine bforce_get_vertices

    !> \brief  is a destructor
    subroutine kill_bforce_opt( self )
        class(bforce_opt), intent(inout) :: self
        if( self%exists )then
            deallocate(self%pb, self%pc)
            self%exists = .false.
        endif
    end subroutine
    
end module simple_bforce_opt
