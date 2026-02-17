!@ descr: module defining the user interfaces for utils programs in the simple_test_exec suite
module simple_test_ui_utils
use simple_ui_modules
implicit none

type(ui_program), target :: simple_test_cmdline
type(ui_program), target :: simple_test_stringmatch
type(ui_program), target :: simple_test_ansi_colors
type(ui_program), target :: simple_test_units
type(ui_program), target :: simple_test_serialize
type(ui_program), target :: simple_test_install
type(ui_program), target :: simple_test_nice
type(ui_program), target :: simple_test_binoris_io
type(ui_program), target :: simple_test_binoris

contains

    subroutine construct_utils_programs(prgtab)
        class(ui_hash), intent(inout) :: prgtab
        call new_simple_test_cmdline(prgtab)
        call new_simple_test_stringmatch(prgtab)
        call new_simple_test_ansi_colors(prgtab)
        call new_simple_test_units(prgtab)
        call new_simple_test_serialize(prgtab)
        call new_simple_test_install(prgtab)
        call new_simple_test_nice(prgtab)
        call new_simple_test_binoris_io(prgtab)
        call new_simple_test_binoris(prgtab)
    end subroutine construct_utils_programs

    subroutine print_utils_programs(logfhandle)
        integer, intent(in) :: logfhandle
        write(logfhandle,'(A)') format_str('UTILS:', C_UNDERLINED)
        write(logfhandle,'(A)') simple_test_cmdline%name%to_char()
        write(logfhandle,'(A)') simple_test_stringmatch%name%to_char()
        write(logfhandle,'(A)') simple_test_ansi_colors%name%to_char()
        write(logfhandle,'(A)') simple_test_units%name%to_char()
        write(logfhandle,'(A)') simple_test_serialize%name%to_char()
        write(logfhandle,'(A)') simple_test_install%name%to_char()
        write(logfhandle,'(A)') simple_test_nice%name%to_char()
        write(logfhandle,'(A)') simple_test_binoris_io%name%to_char()
        write(logfhandle,'(A)') simple_test_binoris%name%to_char()
        write(logfhandle,'(A)') ''
    end subroutine print_utils_programs

    subroutine new_simple_test_cmdline( prgtab )
        class(ui_hash), intent(inout) :: prgtab
        ! add to ui_hash
        call add_ui_program('simple_test_cmdline', simple_test_cmdline, prgtab)
    end subroutine new_simple_test_cmdline

    subroutine new_simple_test_stringmatch( prgtab )
        class(ui_hash), intent(inout) :: prgtab
        ! add to ui_hash
        call add_ui_program('simple_test_stringmatch', simple_test_stringmatch, prgtab)
    end subroutine new_simple_test_stringmatch

    subroutine new_simple_test_ansi_colors( prgtab )
        class(ui_hash), intent(inout) :: prgtab
        ! add to ui_hash
        call add_ui_program('simple_test_ansi_colors', simple_test_ansi_colors, prgtab)
    end subroutine new_simple_test_ansi_colors

    subroutine new_simple_test_units( prgtab )
        class(ui_hash), intent(inout) :: prgtab
        ! add to ui_hash
        call add_ui_program('simple_test_units', simple_test_units, prgtab)
    end subroutine new_simple_test_units

    subroutine new_simple_test_serialize( prgtab )
        class(ui_hash), intent(inout) :: prgtab
        ! add to ui_hash
        call add_ui_program('simple_test_serialize', simple_test_serialize, prgtab)
    end subroutine new_simple_test_serialize

    subroutine new_simple_test_install( prgtab )
        class(ui_hash), intent(inout) :: prgtab
        ! add to ui_hash
        call add_ui_program('simple_test_install', simple_test_install, prgtab)
    end subroutine new_simple_test_install

    subroutine new_simple_test_nice( prgtab )
        class(ui_hash), intent(inout) :: prgtab
        ! add to ui_hash
        call add_ui_program('simple_test_nice', simple_test_nice, prgtab)
    end subroutine new_simple_test_nice

    subroutine new_simple_test_binoris_io( prgtab )
        class(ui_hash), intent(inout) :: prgtab
        ! add to ui_hash
        call add_ui_program('simple_test_binoris_io', simple_test_binoris_io, prgtab)
    end subroutine new_simple_test_binoris_io

    subroutine new_simple_test_binoris( prgtab )
        class(ui_hash), intent(inout) :: prgtab
        ! add to ui_hash
        call add_ui_program('simple_test_binoris', simple_test_binoris, prgtab)
    end subroutine new_simple_test_binoris

end module simple_test_ui_utils
