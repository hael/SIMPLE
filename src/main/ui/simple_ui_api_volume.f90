!@descr: "volume" UI api (concrete implementation)
module simple_ui_api_volume
use simple_ui_program, only: ui_program
use simple_ui_hash,    only: ui_hash
use simple_ui_utils,   only: add_ui_program
implicit none
public :: register_ui_volume

type(ui_program), target :: center
type(ui_program), target :: reproject
type(ui_program), target :: volops

contains

    subroutine register_ui_volume(prgtab)
        class(ui_hash), intent(inout) :: prgtab
        call add_ui_program('center',    center,    prgtab)
        call add_ui_program('reproject', reproject, prgtab)
        call add_ui_program('volops',    volops,    prgtab)
    end subroutine register_ui_volume

end module simple_ui_api_volume
