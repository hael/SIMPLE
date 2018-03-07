! abstract strategy3D base class
module simple_strategy3D
implicit none

public :: strategy3D
private

type, abstract :: strategy3D
  contains
    procedure(generic_new),         deferred :: new
    procedure(generic_srch),        deferred :: srch
    procedure(generic_oris_assign), deferred :: oris_assign
    procedure(generic_kill),        deferred :: kill
end type strategy3D

abstract interface

    subroutine generic_new( self, spec, npeaks )
        use simple_strategy3D_srch, only: strategy3D_spec
        import :: strategy3D
        class(strategy3D),      intent(inout) :: self
        class(strategy3D_spec), intent(inout) :: spec
        integer,                intent(in)    :: npeaks
    end subroutine generic_new

    subroutine generic_srch( self )
        import :: strategy3D
        class(strategy3D), intent(inout) :: self
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
