!@ descr: module defining the user interfaces for optimize programs in the simple_test_exec suite
module simple_test_ui_optimize
use simple_ui_modules
implicit none

type(ui_program), target :: simple_test_lbfgsb
type(ui_program), target :: simple_test_lbfgsb_cosine
type(ui_program), target :: simple_test_opt_lp
type(ui_program), target :: simple_test_lplims
type(ui_program), target :: simple_test_lpstages
type(ui_program), target :: simple_test_tree_srch

contains

    subroutine construct_optimize_programs(prgtab)
        class(ui_hash), intent(inout) :: prgtab
        call new_simple_test_lbfgsb(prgtab)
        call new_simple_test_lbfgsb_cosine(prgtab)
        call new_simple_test_opt_lp(prgtab)
        call new_simple_test_lplims(prgtab)
        call new_simple_test_lpstages(prgtab)
        call new_simple_test_tree_srch(prgtab)
    end subroutine construct_optimize_programs

    subroutine print_optimize_programs(logfhandle)
        integer, intent(in) :: logfhandle
        write(logfhandle,'(A)') format_str('OPTMIZE:', C_UNDERLINED)
        write(logfhandle,'(A)') simple_test_lbfgsb%name%to_char()
        write(logfhandle,'(A)') simple_test_lbfgsb_cosine%name%to_char()
        write(logfhandle,'(A)') simple_test_opt_lp%name%to_char()
        write(logfhandle,'(A)') simple_test_lplims%name%to_char()
        write(logfhandle,'(A)') simple_test_lpstages%name%to_char()
        write(logfhandle,'(A)') simple_test_tree_srch%name%to_char()
        write(logfhandle,'(A)') ''
    end subroutine print_optimize_programs

    subroutine new_simple_test_lbfgsb( prgtab )
        class(ui_hash), intent(inout) :: prgtab
        ! add to ui_hash
        call add_ui_program('simple_test_lbfgsb', simple_test_lbfgsb, prgtab)
    end subroutine new_simple_test_lbfgsb

    subroutine new_simple_test_lbfgsb_cosine( prgtab )
        class(ui_hash), intent(inout) :: prgtab
        ! add to ui_hash
        call add_ui_program('simple_test_lbfgsb_cosine', simple_test_lbfgsb_cosine, prgtab)
    end subroutine new_simple_test_lbfgsb_cosine

    subroutine new_simple_test_opt_lp( prgtab )
        class(ui_hash), intent(inout) :: prgtab
        ! add to ui_hash
        call add_ui_program('simple_test_opt_lp', simple_test_opt_lp, prgtab)
    end subroutine new_simple_test_opt_lp

    subroutine new_simple_test_lplims( prgtab )
        class(ui_hash), intent(inout) :: prgtab
        ! add to ui_hash
        call add_ui_program('simple_test_lplims', simple_test_lplims, prgtab)
    end subroutine new_simple_test_lplims

    subroutine new_simple_test_lpstages( prgtab )
        class(ui_hash), intent(inout) :: prgtab
        ! add to ui_hash
        call add_ui_program('simple_test_lpstages', simple_test_lpstages, prgtab)
    end subroutine new_simple_test_lpstages

    subroutine new_simple_test_tree_srch( prgtab )
        class(ui_hash), intent(inout) :: prgtab
        ! add to ui_hash
        call add_ui_program('simple_test_tree_srch', simple_test_tree_srch, prgtab)
    end subroutine new_simple_test_tree_srch

end module simple_test_ui_optimize
