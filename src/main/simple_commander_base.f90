! abstract commander
module simple_commander_base
include 'simple_lib.f08'
implicit none

public :: commander_base, execute_commander
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

contains

    subroutine execute_commander( self, cline )
        use simple_cmdline, only: cmdline
        class(commander_base), intent(inout) :: self
        class(cmdline),        intent(inout) :: cline
        type(cmdline)         :: cline_original
        character(len=STDLEN) :: cwd
        integer :: nrestarts, i
        nrestarts = 1
        if( cline%defined('nrestarts') ) nrestarts = nint(cline%get_rarg('nrestarts'))
        if( nrestarts > 1 ) call cline%set('mkdir', 'yes') ! need to make unique execution dir for multiple restarts
        ! update original CWD global in defs
        call simple_getcwd(cwd)
        if( allocated(CWD_GLOB_ORIGINAL) ) deallocate(CWD_GLOB_ORIGINAL)
        allocate(CWD_GLOB_ORIGINAL, source=trim(cwd))
        ! save original command line
        cline_original = cline
        do i=1,nrestarts
            call self%execute(cline)
            ! go back to original working directory
            print *, 'IN COMMANDER_BASE, CHANGING BACK TO DIR: ', CWD_GLOB_ORIGINAL
            call simple_chdir(CWD_GLOB_ORIGINAL)
            ! put back the original command line
            cline = cline_original
        end do
    end subroutine execute_commander

end module simple_commander_base
