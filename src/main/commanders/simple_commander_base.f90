! abstract commander
module simple_commander_base
include 'simple_lib.f08'
implicit none

public :: commander_base
private

type, abstract :: commander_base
  contains
  procedure(generic_execute), deferred :: execute
  procedure                            :: execute_safe
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

    subroutine execute_safe( self, cline )
        use simple_cmdline,    only: cmdline
        use simple_parameters, only: parameters, params_glob
        use simple_builder,    only: builder, build_glob
        class(commander_base), intent(inout) :: self
        class(cmdline),        intent(inout) :: cline
        class(parameters), pointer :: params_ptr
        class(builder),    pointer :: build_ptr
        params_ptr => params_glob
        build_ptr  => build_glob
        nullify(build_glob,params_glob)
        call self%execute(cline)
        params_glob => params_ptr
        build_glob  => build_ptr
        nullify(build_ptr,params_ptr)
    end subroutine execute_safe

end module simple_commander_base
