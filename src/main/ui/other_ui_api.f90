module other_ui_api
use simple_ui_program, only: ui_program
implicit none
type(ui_program), target :: mkdir_
type(ui_program), target :: normalize_
type(ui_program), target :: print_ui_stream
type(ui_program), target :: match_stacks
type(ui_program), target :: split_
end module other_ui_api
