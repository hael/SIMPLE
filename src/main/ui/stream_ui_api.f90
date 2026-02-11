module stream_ui_api
use simple_ui_program, only: ui_program
implicit none
type(ui_program), target :: abinitio2D_stream
type(ui_program), target :: assign_optics
type(ui_program), target :: cluster2D_stream
type(ui_program), target :: gen_pickrefs
type(ui_program), target :: pick_extract
type(ui_program), target :: preproc
type(ui_program), target :: sieve_cavgs
end module stream_ui_api
