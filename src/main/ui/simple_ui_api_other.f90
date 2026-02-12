module simple_ui_api_other
use simple_ui_program, only: ui_program
use simple_ui_hash,    only: ui_hash
implicit none
type(ui_program), target :: match_stacks
type(ui_program), target :: mkdir_
type(ui_program), target :: normalize_
type(ui_program), target :: print_ui_stream
type(ui_program), target :: split_
type(ui_program), target :: split_stack

contains

    subroutine register_simple_ui_other( prgtab )
        class(ui_hash), intent(inout) :: prgtab
        call prgtab%set_ui_program('match_stacks',    match_stacks)
        call prgtab%set_ui_program('mkdir_',          mkdir_)
        call prgtab%set_ui_program('normalize_',      normalize_)
        call prgtab%set_ui_program('print_ui_stream', print_ui_stream)
        call prgtab%set_ui_program('split_',          split_)
        call prgtab%set_ui_program('split_stack',     split_stack)
    end subroutine register_simple_ui_other

end module simple_ui_api_other
