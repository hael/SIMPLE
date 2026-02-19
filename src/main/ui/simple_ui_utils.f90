!@descr: module defining utility procedures for the simple_ui modules
module simple_ui_utils
use simple_ui_program, only: ui_program
use simple_ui_hash,    only: ui_hash
use simple_error,      only: simple_exception
implicit none
#include "simple_local_flags.inc"

contains

    subroutine add_ui_program( prg_name, ui_prg_instance, prgtab )
        character(len=*),          intent(in)    :: prg_name
        class(ui_program), target, intent(inout) :: ui_prg_instance
        class(ui_hash),            intent(inout) :: prgtab
        character(len=:), allocatable :: k
        k = trim(adjustl(prg_name))
        if( prgtab%has_key(k) )then
            THROW_HARD('Key: '//k//' already in ui_hash')
        endif
        call prgtab%set_ui_program(k, ui_prg_instance)
    end subroutine add_ui_program

end module simple_ui_utils
