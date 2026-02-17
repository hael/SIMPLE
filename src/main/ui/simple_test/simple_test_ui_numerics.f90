!@ descr: module defining the user interfaces for numerics programs in the simple_test_exec suite
module simple_test_ui_numerics
use simple_ui_modules
implicit none

type(ui_program), target :: simple_test_eigh
type(ui_program), target :: simple_test_kbinterpol_fast
type(ui_program), target :: simple_test_neigh
type(ui_program), target :: simple_test_maxnloc

contains

    subroutine construct_numerics_programs(prgtab)
        class(ui_hash), intent(inout) :: prgtab
        call new_simple_test_eigh(prgtab)
        call new_simple_test_kbinterpol_fast(prgtab)
        call new_simple_test_neigh(prgtab)
        call new_simple_test_maxnloc(prgtab)
    end subroutine construct_numerics_programs

    subroutine print_numerics_programs(logfhandle)
        integer, intent(in) :: logfhandle
        write(logfhandle,'(A)') format_str('NUMERICS:', C_UNDERLINED)
        write(logfhandle,'(A)') simple_test_eigh%name%to_char()
        write(logfhandle,'(A)') simple_test_kbinterpol_fast%name%to_char()
        write(logfhandle,'(A)') simple_test_neigh%name%to_char()
        write(logfhandle,'(A)') simple_test_maxnloc%name%to_char()
        write(logfhandle,'(A)') ''
    end subroutine print_numerics_programs

    subroutine new_simple_test_eigh( prgtab )
        class(ui_hash), intent(inout) :: prgtab
        ! add to ui_hash
        call add_ui_program('simple_test_eigh', simple_test_eigh, prgtab)
    end subroutine new_simple_test_eigh

    subroutine new_simple_test_kbinterpol_fast( prgtab )
        class(ui_hash), intent(inout) :: prgtab
        ! add to ui_hash
        call add_ui_program('simple_test_kbinterpol_fast', simple_test_kbinterpol_fast, prgtab)
    end subroutine new_simple_test_kbinterpol_fast

    subroutine new_simple_test_neigh( prgtab )
        class(ui_hash), intent(inout) :: prgtab
        ! add to ui_hash
        call add_ui_program('simple_test_neigh', simple_test_neigh, prgtab)
    end subroutine new_simple_test_neigh

    subroutine new_simple_test_maxnloc( prgtab )
        class(ui_hash), intent(inout) :: prgtab
        ! add to ui_hash
        call add_ui_program('simple_test_maxnloc', simple_test_maxnloc, prgtab)
    end subroutine new_simple_test_maxnloc

end module simple_test_ui_numerics
