!@descr: "single_ui_api_tseries" UI api (concrete implementation)
module single_ui_api_tseries
use simple_ui_program, only: ui_program
use simple_ui_hash,    only: ui_hash
use simple_ui_utils,   only: add_ui_program
implicit none

type(ui_program), target :: track_particles
type(ui_program), target :: tseries_import
type(ui_program), target :: tseries_make_pickavg
type(ui_program), target :: tseries_motion_correct

contains

    subroutine register_single_ui_tseries(prgtab)
        class(ui_hash), intent(inout) :: prgtab
        call add_ui_program('track_particles', track_particles, prgtab)
        call add_ui_program('tseries_import', tseries_import, prgtab)
        call add_ui_program('tseries_make_pickavg', tseries_make_pickavg, prgtab)
        call add_ui_program('tseries_motion_correct', tseries_motion_correct, prgtab)
    end subroutine register_single_ui_tseries

end module single_ui_api_tseries
