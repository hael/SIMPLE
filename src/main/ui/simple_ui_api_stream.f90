!@descr: "ui_api_stream" UI api (concrete implementation)
module simple_ui_api_stream
use simple_ui_program, only: ui_program
use simple_ui_hash,    only: ui_hash
use simple_ui_utils,   only: add_ui_program
implicit none

type(ui_program), target :: abinitio2D_stream
type(ui_program), target :: assign_optics
type(ui_program), target :: cluster2D_stream
type(ui_program), target :: gen_pickrefs
type(ui_program), target :: pick_extract
type(ui_program), target :: preproc
type(ui_program), target :: sieve_cavgs

contains

    subroutine register_simple_ui_stream(prgtab)
        class(ui_hash), intent(inout) :: prgtab
        call add_ui_program('abinitio2D_stream', abinitio2D_stream, prgtab)
        call add_ui_program('assign_optics',     assign_optics,     prgtab)
        call add_ui_program('cluster2D_stream',  cluster2D_stream,  prgtab)
        call add_ui_program('gen_pickrefs',      gen_pickrefs,      prgtab)
        call add_ui_program('pick_extract',      pick_extract,      prgtab)
        call add_ui_program('preproc',           preproc,           prgtab)
        call add_ui_program('sieve_cavgs',       sieve_cavgs,       prgtab)
    end subroutine register_simple_ui_stream

end module simple_ui_api_stream
