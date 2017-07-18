!------------------------------------------------------------------------------!
! SIMPLE v2.5         Elmlund & Elmlund Lab          simplecryoem.com          !
!------------------------------------------------------------------------------!
!> Simple optimizer module: the interface to the simple optimizer class.
!
! The code is distributed with the hope that it will be useful, but _WITHOUT_
! _ANY_ _WARRANTY_. Redistribution or modification is regulated by the GNU
! General Public License. *Author:* Hans Elmlund, 2014-01-07

module simple_optimizer
implicit none

public :: optimizer
private

type, abstract :: optimizer
  contains
    procedure(generic_new),          deferred :: new
    procedure(generic_minimize),     deferred :: minimize
    procedure(generic_get_vertices), deferred :: get_vertices
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
    subroutine generic_minimize( self, spec, lowest_cost )
        use simple_opt_spec, only: opt_spec
        import :: optimizer
        class(optimizer), intent(inout) :: self
        class(opt_spec), intent(inout)  :: spec
        real, intent(out)               :: lowest_cost
    end subroutine generic_minimize

    !> \brief  getter of vertices and their cost for simplex only
    subroutine generic_get_vertices( self, spec, vertices, costs )
        use simple_opt_spec, only: opt_spec
        import :: optimizer
        class(optimizer),  intent(inout) :: self
        class(opt_spec),   intent(inout) :: spec
        real, allocatable ,intent(inout) :: vertices(:,:), costs(:)
    end subroutine generic_get_vertices

    !>  \brief  is a destructor
    subroutine generic_kill( self )
        import :: optimizer
        class(optimizer), intent(inout) :: self
    end subroutine generic_kill

end interface

end module simple_optimizer
