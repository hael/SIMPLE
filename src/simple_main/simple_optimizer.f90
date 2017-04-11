!==Class simple_optimizer
!
! Is the abstract simple optimizer class.
! The code is distributed with the hope that it will be useful, but _WITHOUT_ _ANY_ _WARRANTY_. 
! Redistribution or modification is regulated by the GNU General Public License. 
! *Author:* Hans Elmlund, 2014-01-07

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
end type

abstract interface

    !>  \brief  is a constructor
    subroutine generic_new( self, spec ) 
        use simple_opt_spec, only: opt_spec
        import :: optimizer
        class(optimizer), intent(inout) :: self
        class(opt_spec), intent(inout)  :: spec
    end subroutine generic_new
    
    !> \brief  minimization of the costfunction
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
