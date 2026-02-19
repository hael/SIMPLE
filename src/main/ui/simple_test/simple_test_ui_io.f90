!@descr: module defining the user interfaces for io  programs in the simple_test_exec suite
module simple_test_ui_io
use simple_ui_modules
implicit none

type(ui_program), target :: imgfile
type(ui_program), target :: inside_write
type(ui_program), target :: io
type(ui_program), target :: io_parallel
type(ui_program), target :: mrc2jpeg
type(ui_program), target :: mrc_validation
type(ui_program), target :: stack_io
type(ui_program), target :: star_export
type(ui_program), target :: starfile_test

contains

    subroutine construct_io_programs(prgtab)
        class(ui_hash), intent(inout) :: prgtab
        call new_imgfile(prgtab)
        call new_inside_write(prgtab)
        call new_io(prgtab)
        call new_io_parallel(prgtab)
        call new_mrc2jpeg(prgtab)
        call new_mrc_validation(prgtab)
        call new_stack_io(prgtab)
        call new_star_export(prgtab)
        call new_starfile_test(prgtab)
    end subroutine construct_io_programs

    subroutine print_io_programs(logfhandle)
        integer, intent(in) :: logfhandle
        write(logfhandle,'(A)') format_str('INPUT/OUTPUT:', C_UNDERLINED)
        write(logfhandle,'(A)') imgfile%name%to_char()
        write(logfhandle,'(A)') inside_write%name%to_char()
        write(logfhandle,'(A)') io%name%to_char()
        write(logfhandle,'(A)') io_parallel%name%to_char()
        write(logfhandle,'(A)') mrc2jpeg%name%to_char()
        write(logfhandle,'(A)') mrc_validation%name%to_char()
        write(logfhandle,'(A)') stack_io%name%to_char()
        write(logfhandle,'(A)') star_export%name%to_char()
        write(logfhandle,'(A)') starfile_test%name%to_char()
        write(logfhandle,'(A)') ''
    end subroutine print_io_programs

    subroutine new_imgfile( prgtab )
        class(ui_hash), intent(inout) :: prgtab
        ! PROGRAM SPECIFICATION
        call imgfile%new(&
        &'imgfile',&                           ! name
        &'imgfile ',&                          ! descr_short
        &'is a test program for ',&
        &'simple_test_exec',&                  ! executable
        &.false.)                              ! requires sp_project
        ! INPUT PARAMETER SPECIFICATIONS
        ! image input/output
        !call imgfile%add_input(UI_IO, )
        ! parameter input/output
        !call imgfile%add_input(UI_IMG, )
        ! alternative inputs
        !call imgfile%add_input(UI_PARM, )
        ! search controls
        !call imgfile%add_input(UI_SRCH, )
        ! filter controls
        !call imgfile%add_input(UI_FILT, )
        ! mask controls
        !call imgfile%add_input(UI_MASK, )
        ! computer controls
        !call imgfile%add_input(UI_COMP, )
        ! add to ui_hash
        call add_ui_program('imgfile', imgfile, prgtab)
    end subroutine new_imgfile

    subroutine new_inside_write( prgtab )
        class(ui_hash), intent(inout) :: prgtab
        ! PROGRAM SPECIFICATION
        call inside_write%new(&
        &'inside_write',&                      ! name
        &'inside_write ',&                     ! descr_short
        &'is a test program for ',&
        &'simple_test_exec',&                  ! executable
        &.false.)                              ! requires sp_project
        ! INPUT PARAMETER SPECIFICATIONS
        ! image input/output
        !call inside_write%add_input(UI_IO, )
        ! parameter input/output
        !call inside_write%add_input(UI_IMG, )
        ! alternative inputs
        !call inside_write%add_input(UI_PARM, )
        ! search controls
        !call inside_write%add_input(UI_SRCH, )
        ! filter controls
        !call inside_write%add_input(UI_FILT, )
        ! mask controls
        !call inside_write%add_input(UI_MASK, )
        ! computer controls
        !call inside_write%add_input(UI_COMP, )
        ! add to ui_hash
        call add_ui_program('inside_write', inside_write, prgtab)
    end subroutine new_inside_write

    subroutine new_io( prgtab )
        class(ui_hash), intent(inout) :: prgtab
        ! PROGRAM SPECIFICATION
        call io%new(&
        &'io',&                                ! name
        &'io ',&                               ! descr_short
        &'is a test program for ',&
        &'simple_test_exec',&                  ! executable
        &.false.)                              ! requires sp_project
        ! INPUT PARAMETER SPECIFICATIONS
        ! image input/output
        !call io%add_input(UI_IO, )
        ! parameter input/output
        !call io%add_input(UI_IMG, )
        ! alternative inputs
        !call io%add_input(UI_PARM, )
        ! search controls
        !call io%add_input(UI_SRCH, )
        ! filter controls
        !call io%add_input(UI_FILT, )
        ! mask controls
        !call io%add_input(UI_MASK, )
        ! computer controls
        !call io%add_input(UI_COMP, )
        ! add to ui_hash
        call add_ui_program('io', io, prgtab)
    end subroutine new_io

    subroutine new_io_parallel( prgtab )
        class(ui_hash), intent(inout) :: prgtab
        ! PROGRAM SPECIFICATION
        call io_parallel%new(&
        &'io_parallel',&                       ! name
        &'io_parallel ',&                      ! descr_short
        &'is a test program for ',&
        &'simple_test_exec',&                  ! executable
        &.false.)                              ! requires sp_project
        ! INPUT PARAMETER SPECIFICATIONS
        ! image input/output
        !call io_parallel%add_input(UI_IO, )
        ! parameter input/output
        !call io_parallel%add_input(UI_IMG, )
        ! alternative inputs
        !call io_parallel%add_input(UI_PARM, )
        ! search controls
        !call io_parallel%add_input(UI_SRCH, )
        ! filter controls
        !call io_parallel%add_input(UI_FILT, )
        ! mask controls
        !call io_parallel%add_input(UI_MASK, )
        ! computer controls
        !call io_parallel%add_input(UI_COMP, )
        ! add to ui_hash
        call add_ui_program('io_parallel', io_parallel, prgtab)
    end subroutine new_io_parallel

    subroutine new_mrc2jpeg( prgtab )
        class(ui_hash), intent(inout) :: prgtab
        ! PROGRAM SPECIFICATION
        call mrc2jpeg%new(&
        &'mrc2jpeg',&                          ! name
        &'mrc2jpeg ',&                         ! descr_short
        &'is a test program for ',&
        &'simple_test_exec',&                  ! executable
        &.false.)                              ! requires sp_project
        ! INPUT PARAMETER SPECIFICATIONS
        ! image input/output
        !call mrc2jpeg%add_input(UI_IO, )
        ! parameter input/output
        !call mrc2jpeg%add_input(UI_IMG, )
        ! alternative inputs
        !call mrc2jpeg%add_input(UI_PARM, )
        ! search controls
        !call mrc2jpeg%add_input(UI_SRCH, )
        ! filter controls
        !call mrc2jpeg%add_input(UI_FILT, )
        ! mask controls
        !call mrc2jpeg%add_input(UI_MASK, )
        ! computer controls
        !call mrc2jpeg%add_input(UI_COMP, )
        ! add to ui_hash
        call add_ui_program('mrc2jpeg', mrc2jpeg, prgtab)
    end subroutine new_mrc2jpeg

    subroutine new_mrc_validation( prgtab )
        class(ui_hash), intent(inout) :: prgtab
        ! PROGRAM SPECIFICATION
        call mrc_validation%new(&
        &'mrc_validation',&                    ! name
        &'mrc_validation ',&                   ! descr_short
        &'is a test program for ',&
        &'simple_test_exec',&                  ! executable
        &.false.)                              ! requires sp_project
        ! INPUT PARAMETER SPECIFICATIONS
        ! image input/output
        !call mrc_validation%add_input(UI_IO, )
        ! parameter input/output
        !call mrc_validation%add_input(UI_IMG, )
        ! alternative inputs
        !call mrc_validation%add_input(UI_PARM, )
        ! search controls
        !call mrc_validation%add_input(UI_SRCH, )
        ! filter controls
        !call mrc_validation%add_input(UI_FILT, )
        ! mask controls
        !call mrc_validation%add_input(UI_MASK, )
        ! computer controls
        !call mrc_validation%add_input(UI_COMP, )
        ! add to ui_hash
        call add_ui_program('mrc_validation', mrc_validation, prgtab)
    end subroutine new_mrc_validation

    subroutine new_stack_io( prgtab )
        class(ui_hash), intent(inout) :: prgtab
        ! PROGRAM SPECIFICATION
        call stack_io%new(&
        &'stack_io',&                          ! name
        &'stack_io ',&                         ! descr_short
        &'is a test program for ',&
        &'simple_test_exec',&                  ! executable
        &.false.)                              ! requires sp_project
        ! INPUT PARAMETER SPECIFICATIONS
        ! image input/output
        !call stack_io%add_input(UI_IO, )
        ! parameter input/output
        !call stack_io%add_input(UI_IMG, )
        ! alternative inputs
        !call stack_io%add_input(UI_PARM, )
        ! search controls
        !call stack_io%add_input(UI_SRCH, )
        ! filter controls
        !call stack_io%add_input(UI_FILT, )
        ! mask controls
        !call stack_io%add_input(UI_MASK, )
        ! computer controls
        !call stack_io%add_input(UI_COMP, )
        ! add to ui_hash
        call add_ui_program('stack_io', stack_io, prgtab)
    end subroutine new_stack_io

    subroutine new_star_export( prgtab )
        class(ui_hash), intent(inout) :: prgtab
        ! PROGRAM SPECIFICATION
        call star_export%new(&
        &'star_export',&                       ! name
        &'star_export ',&                      ! descr_short
        &'is a test program for ',&
        &'simple_test_exec',&                  ! executable
        &.false.)                              ! requires sp_project
        ! INPUT PARAMETER SPECIFICATIONS
        ! image input/output
        !call star_export%add_input(UI_IO, )
        ! parameter input/output
        !call star_export%add_input(UI_IMG, )
        ! alternative inputs
        !call star_export%add_input(UI_PARM, )
        ! search controls
        !call star_export%add_input(UI_SRCH, )
        ! filter controls
        !call star_export%add_input(UI_FILT, )
        ! mask controls
        !call star_export%add_input(UI_MASK, )
        ! computer controls
        !call star_export%add_input(UI_COMP, )
        ! add to ui_hash
        call add_ui_program('star_export', star_export, prgtab)
    end subroutine new_star_export

    subroutine new_starfile_test( prgtab )
        class(ui_hash), intent(inout) :: prgtab
        ! PROGRAM SPECIFICATION
        call starfile_test%new(&
        &'starfile_test',&                     ! name
        &'starfile_test ',&                    ! descr_short
        &'is a test program for ',&
        &'simple_test_exec',&                  ! executable
        &.false.)                              ! requires sp_project
        ! INPUT PARAMETER SPECIFICATIONS
        ! image input/output
        !call starfile_test%add_input(UI_IO, )
        ! parameter input/output
        !call starfile_test%add_input(UI_IMG, )
        ! alternative inputs
        !call starfile_test%add_input(UI_PARM, )
        ! search controls
        !call starfile_test%add_input(UI_SRCH, )
        ! filter controls
        !call starfile_test%add_input(UI_FILT, )
        ! mask controls
        !call starfile_test%add_input(UI_MASK, )
        ! computer controls
        !call starfile_test%add_input(UI_COMP, )
        ! add to ui_hash
        call add_ui_program('starfile_test', starfile_test, prgtab)
    end subroutine new_starfile_test

end module simple_test_ui_io
