!@descr: abstract base class defining the common strategy3D interface
module simple_strategy3D
use simple_strategy3D_srch, only: strategy3D_srch, strategy3D_spec
use simple_oris,            only: oris
implicit none

public :: strategy3D
private

type, abstract :: strategy3D
    type(strategy3D_srch) :: s
    type(strategy3D_spec) :: spec
contains
    procedure(generic_new),         deferred :: new
    procedure(generic_srch),        deferred :: srch
    procedure(generic_oris_assign), deferred :: oris_assign
    procedure(generic_kill),        deferred :: kill
end type strategy3D

abstract interface

    subroutine generic_new( self, spec, build )
        use simple_strategy3D_srch, only: strategy3D_spec
        use simple_builder,         only: builder
        import :: strategy3D
        class(strategy3D),      intent(inout) :: self
        class(strategy3D_spec), intent(inout) :: spec
        class(builder), target, intent(inout) :: build
    end subroutine generic_new

    subroutine generic_srch( self, os, ithr )
        import :: strategy3D
        import :: oris
        class(strategy3D), intent(inout) :: self
        class(oris),       intent(inout) :: os
        integer,           intent(in)    :: ithr
    end subroutine generic_srch

    subroutine generic_oris_assign( self )
        import :: strategy3D
        class(strategy3D), intent(inout) :: self
    end subroutine generic_oris_assign

    subroutine generic_kill( self )
        import :: strategy3D
        class(strategy3D), intent(inout) :: self
    end subroutine generic_kill

end interface

end module simple_strategy3D
