module other_ui_api
use simple_ui_program, only: ui_program
use simple_ui_hash,    only: ui_hash
implicit none
type(ui_program), target :: mkdir_
type(ui_program), target :: normalize_
type(ui_program), target :: print_ui_stream
type(ui_program), target :: match_stacks
type(ui_program), target :: split_

contains

    subroutine add_other_ui_api2prgtab( prgtab )
        class(ui_hash), intent(inout) :: prgtab
        call prgtab%set_ui_program('mkdir_',          mkdir_)
        call prgtab%set_ui_program('normalize_',      normalize_)
        call prgtab%set_ui_program('print_ui_stream', print_ui_stream)
        call prgtab%set_ui_program('match_stacks',    match_stacks)
        call prgtab%set_ui_program('split_',          split_)
    end subroutine add_other_ui_api2prgtab

end module other_ui_api
