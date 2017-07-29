!------------------------------------------------------------------------------!
! SIMPLE v2.5         Elmlund & Elmlund Lab          simplecryoem.com          !
!------------------------------------------------------------------------------!
!> \brief Abstract module for commander classes
module simple_commander_base
implicit none

public :: commander_base
private

type, abstract :: commander_base
  contains
    procedure(generic_execute), deferred :: execute
end type commander_base

abstract interface

    !>  \brief  executes the commander
    subroutine generic_execute( self, cline )
        use simple_cmdline, only: cmdline
        import :: commander_base
        class(commander_base), intent(inout) :: self
        class(cmdline),        intent(inout) :: cline
    end subroutine generic_execute

end interface

end module simple_commander_base
