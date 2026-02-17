!@ descr: module defining the user interfaces for io programs in the simple_test_exec suite
module simple_test_ui_io
use simple_ui_modules
implicit none

type(ui_program), target :: simple_test_imgfile
type(ui_program), target :: simple_test_io
type(ui_program), target :: simple_test_io_parallel
type(ui_program), target :: simple_test_stack_io
type(ui_program), target :: simple_test_mrc_validation
type(ui_program), target :: simple_test_mrc2jpeg
type(ui_program), target :: simple_test_starfile
type(ui_program), target :: simple_test_star_export
type(ui_program), target :: simple_test_inside_write

contains

    subroutine construct_io_programs(prgtab)
        class(ui_hash), intent(inout) :: prgtab
        call new_simple_test_imgfile(prgtab)
        call new_simple_test_io(prgtab)
        call new_simple_test_io_parallel(prgtab)
        call new_simple_test_stack_io(prgtab)
        call new_simple_test_mrc_validation(prgtab)
        call new_simple_test_mrc2jpeg(prgtab)
        call new_simple_test_starfile(prgtab)
        call new_simple_test_star_export(prgtab)
        call new_simple_test_inside_write(prgtab)
    end subroutine construct_io_programs

    subroutine print_io_programs(logfhandle)
        integer, intent(in) :: logfhandle
        write(logfhandle,'(A)') format_str('INPUT/OUTPUT:', C_UNDERLINED)
        write(logfhandle,'(A)') simple_test_imgfile%name%to_char()
        write(logfhandle,'(A)') simple_test_io%name%to_char()
        write(logfhandle,'(A)') simple_test_io_parallel%name%to_char()
        write(logfhandle,'(A)') simple_test_stack_io%name%to_char()
        write(logfhandle,'(A)') simple_test_mrc_validation%name%to_char()
        write(logfhandle,'(A)') simple_test_mrc2jpeg%name%to_char()
        write(logfhandle,'(A)') simple_test_starfile%name%to_char()
        write(logfhandle,'(A)') simple_test_star_export%name%to_char()
        write(logfhandle,'(A)') simple_test_inside_write%name%to_char()
        write(logfhandle,'(A)') ''
    end subroutine print_io_programs

    subroutine new_simple_test_imgfile( prgtab )
        class(ui_hash), intent(inout) :: prgtab
        ! add to ui_hash
        call add_ui_program('simple_test_imgfile', simple_test_imgfile, prgtab)
    end subroutine new_simple_test_imgfile

    subroutine new_simple_test_io( prgtab )
        class(ui_hash), intent(inout) :: prgtab
        ! add to ui_hash
        call add_ui_program('simple_test_io', simple_test_io, prgtab)
    end subroutine new_simple_test_io

    subroutine new_simple_test_io_parallel( prgtab )
        class(ui_hash), intent(inout) :: prgtab
        ! add to ui_hash
        call add_ui_program('simple_test_io_parallel', simple_test_io_parallel, prgtab)
    end subroutine new_simple_test_io_parallel

    subroutine new_simple_test_stack_io( prgtab )
        class(ui_hash), intent(inout) :: prgtab
        ! add to ui_hash
        call add_ui_program('simple_test_stack_io', simple_test_stack_io, prgtab)
    end subroutine new_simple_test_stack_io

    subroutine new_simple_test_mrc_validation( prgtab )
        class(ui_hash), intent(inout) :: prgtab
        ! add to ui_hash
        call add_ui_program('simple_test_mrc_validation', simple_test_mrc_validation, prgtab)
    end subroutine new_simple_test_mrc_validation

    subroutine new_simple_test_mrc2jpeg( prgtab )
        class(ui_hash), intent(inout) :: prgtab
        ! add to ui_hash
        call add_ui_program('simple_test_mrc2jpeg', simple_test_mrc2jpeg, prgtab)
    end subroutine new_simple_test_mrc2jpeg

    subroutine new_simple_test_starfile( prgtab )
        class(ui_hash), intent(inout) :: prgtab
        ! add to ui_hash
        call add_ui_program('simple_test_starfile', simple_test_starfile, prgtab)
    end subroutine new_simple_test_starfile

    subroutine new_simple_test_star_export( prgtab )
        class(ui_hash), intent(inout) :: prgtab
        ! add to ui_hash
        call add_ui_program('simple_test_star_export', simple_test_star_export, prgtab)
    end subroutine new_simple_test_star_export

    subroutine new_simple_test_inside_write( prgtab )
        class(ui_hash), intent(inout) :: prgtab
        ! add to ui_hash
        call add_ui_program('simple_test_inside_write', simple_test_inside_write, prgtab)
    end subroutine new_simple_test_inside_write

end module simple_test_ui_io
