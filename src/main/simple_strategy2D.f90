module simple_strategy2D
implicit none

public :: strategy2D, strategy2D_per_ptcl
private

! abstract strategy2D base class
type, abstract :: strategy2D
  contains
    procedure(generic_new),  deferred :: new
    procedure(generic_srch), deferred :: srch
    procedure(generic_kill), deferred :: kill
end type strategy2D

abstract interface

    subroutine generic_new( self, spec )
        use simple_strategy2D_srch, only: strategy2D_spec
        import :: strategy2D
        class(strategy2D),      intent(inout) :: self
        class(strategy2D_spec), intent(inout) :: spec
    end subroutine generic_new

    subroutine generic_srch( self )
        import :: strategy2D
        class(strategy2D), intent(inout) :: self
    end subroutine generic_srch

    subroutine generic_kill( self )
        import :: strategy2D
        class(strategy2D), intent(inout) :: self
    end subroutine generic_kill

end interface

! This type is to allow particle-dependent decision about which 2D strategy to
! use hybrid or combined search strategies can then be implemented as extensions
! of the relevant strategy2D base class
type :: strategy2D_per_ptcl
    class(strategy2D), pointer :: ptr => null()
end type strategy2D_per_ptcl

end module simple_strategy2D
