! abstract optimiser defining an interface for the extending optimisation classes
module simple_optimizer
implicit none

public :: optimizer
private

type, abstract :: optimizer
  contains
    procedure(generic_new),          deferred :: new
    procedure(generic_minimize),     deferred :: minimize
    procedure(generic_kill),         deferred :: kill
end type optimizer

abstract interface

    !>  \brief  is a constructor
    !! must define specifications for a new optimizer
    !! \param spec specifications object
    subroutine generic_new( self, spec )
        use simple_opt_spec, only: opt_spec
        import :: optimizer
        class(optimizer), intent(inout) :: self
        class(opt_spec), intent(inout)  :: spec
    end subroutine generic_new

    !> \brief  minimization of the cost-function
    !! minimization requires specifications and a return param
    !! \param spec specifications of opt
    !! \param lowest_cost return param for best result at end of optimisation
    subroutine generic_minimize( self, spec, fun_self, lowest_cost )
        use simple_opt_spec, only: opt_spec
        import :: optimizer
        class(optimizer), intent(inout) :: self
        class(opt_spec),  intent(inout) :: spec
        class(*),         intent(inout) :: fun_self        
        real, intent(out)               :: lowest_cost
    end subroutine generic_minimize

    !>  \brief  is a destructor
    subroutine generic_kill( self )
        import :: optimizer
        class(optimizer), intent(inout) :: self
    end subroutine generic_kill

end interface

end module simple_optimizer
