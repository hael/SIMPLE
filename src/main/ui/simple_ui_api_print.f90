!@descr: "print" UI api (concrete implementation)
module simple_ui_api_print
use simple_ui_program, only: ui_program
use simple_ui_hash,    only: ui_hash
use simple_ui_utils,   only: add_ui_program
implicit none

type(ui_program), target :: info_image
type(ui_program), target :: info_stktab
type(ui_program), target :: print_dose_weights
type(ui_program), target :: print_fsc
type(ui_program), target :: print_magic_boxes

contains

    subroutine register_ui_print(prgtab)
        class(ui_hash), intent(inout) :: prgtab
        call add_ui_program('info_image',         info_image,         prgtab)
        call add_ui_program('info_stktab',        info_stktab,        prgtab)
        call add_ui_program('print_dose_weights', print_dose_weights, prgtab)
        call add_ui_program('print_fsc',          print_fsc,          prgtab)
        call add_ui_program('print_magic_boxes',  print_magic_boxes,  prgtab)
    end subroutine register_ui_print

end module simple_ui_api_print
