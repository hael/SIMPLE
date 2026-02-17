!@ descr: module defining the user interfaces for geometry programs in the simple_test_exec suite
module simple_test_ui_geometry
use simple_ui_modules
implicit none

type(ui_program), target :: simple_test_angres
type(ui_program), target :: simple_test_ori
type(ui_program), target :: simple_test_oris
type(ui_program), target :: simple_test_uniform_euler
type(ui_program), target :: simple_test_uniform_rot
type(ui_program), target :: simple_test_sym

contains

    subroutine construct_geometry_programs(prgtab)
        class(ui_hash), intent(inout) :: prgtab
        call new_simple_test_angres(prgtab)
        call new_simple_test_ori(prgtab)
        call new_simple_test_oris(prgtab)
        call new_simple_test_uniform_euler(prgtab)
        call new_simple_test_uniform_rot(prgtab)
        call new_simple_test_sym(prgtab)
    end subroutine construct_geometry_programs

    subroutine print_geometry_programs(logfhandle)
        integer, intent(in) :: logfhandle
        write(logfhandle,'(A)') format_str('GEOMETRY:', C_UNDERLINED)
        write(logfhandle,'(A)') simple_test_angres%name%to_char()
        write(logfhandle,'(A)') simple_test_ori%name%to_char()
        write(logfhandle,'(A)') simple_test_oris%name%to_char()
        write(logfhandle,'(A)') simple_test_uniform_euler%name%to_char()
        write(logfhandle,'(A)') simple_test_uniform_rot%name%to_char()
        write(logfhandle,'(A)') simple_test_sym%name%to_char()
        write(logfhandle,'(A)') ''
    end subroutine print_geometry_programs

    subroutine new_simple_test_angres( prgtab )
        class(ui_hash), intent(inout) :: prgtab
        ! add to ui_hash
        call add_ui_program('simple_test_angres', simple_test_angres, prgtab)
    end subroutine new_simple_test_angres

    subroutine new_simple_test_ori( prgtab )
        class(ui_hash), intent(inout) :: prgtab
        ! add to ui_hash
        call add_ui_program('simple_test_ori', simple_test_ori, prgtab)
    end subroutine new_simple_test_ori

    subroutine new_simple_test_oris( prgtab )
        class(ui_hash), intent(inout) :: prgtab
        ! add to ui_hash
        call add_ui_program('simple_test_oris', simple_test_oris, prgtab)
    end subroutine new_simple_test_oris

    subroutine new_simple_test_uniform_euler( prgtab )
        class(ui_hash), intent(inout) :: prgtab
        ! add to ui_hash
        call add_ui_program('simple_test_uniform_euler', simple_test_uniform_euler, prgtab)
    end subroutine new_simple_test_uniform_euler

    subroutine new_simple_test_uniform_rot( prgtab )
        class(ui_hash), intent(inout) :: prgtab
        ! add to ui_hash
        call add_ui_program('simple_test_uniform_rot', simple_test_uniform_rot, prgtab)
    end subroutine new_simple_test_uniform_rot

    subroutine new_simple_test_sym( prgtab )
        class(ui_hash), intent(inout) :: prgtab
        ! add to ui_hash
        call add_ui_program('simple_test_sym', simple_test_sym, prgtab)
    end subroutine new_simple_test_sym

end module simple_test_ui_geometry
