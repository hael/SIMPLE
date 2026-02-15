module simple_ui_api_other
use simple_ui_api_modules
implicit none

type(ui_program), target :: mkdir_

type(ui_program), target :: split_
type(ui_program), target :: split_stack

contains

    subroutine construct_other_programs(prgtab)
        class(ui_hash), intent(inout) :: prgtab
        call new_mkdir_(prgtab)
        call new_split_(prgtab)
        call new_split_stack(prgtab)
    end subroutine construct_other_programs

    subroutine print_other_programs(logfhandle)
        integer, intent(in) :: logfhandle
        write(logfhandle,'(A)') format_str('OTHER UTILITIES:', C_UNDERLINED)
        write(logfhandle,'(A)') mkdir_%name%to_char()
        write(logfhandle,'(A)') split_%name%to_char()
        write(logfhandle,'(A)') split_stack%name%to_char()
        write(logfhandle,'(A)') ''
    end subroutine print_other_programs

    subroutine new_mkdir_( prgtab )
        class(ui_hash), intent(inout) :: prgtab
        ! PROGRAM SPECIFICATION
        call mkdir_%new(&
        &'mkdir',&                                                       ! name
        &'Make directory',&                                              ! descr_short
        &'is a program for making an automatically numbered directory',& ! descr_long
        &'simple_exec',&                                                 ! executable
        &.false.)                                                        ! requires sp_project
        ! INPUT PARAMETER SPECIFICATIONS
        ! image input/output
        ! <empty>
        ! parameter input/output
        call mkdir_%add_input(UI_PARM, 'dir', 'dir', 'Name of directory to create', 'Name of directory name to create(e.g. X_dir/)', 'e.g. X_dir/', .true., ' ')
        ! alternative inputs
        ! <empty>
        ! search controls
        ! <empty>
        ! filter controls
        ! <empty>
        ! mask controls
        ! <empty>
        ! computer controls
        ! <empty>
        call add_ui_program('mkdir_', mkdir_, prgtab)
    end subroutine new_mkdir_

    subroutine new_split_( prgtab )
        class(ui_hash), intent(inout) :: prgtab
        call split_%new(&
        &'split',&                                   ! name
        &'Split stack into substacks',&              ! descr_short
        &'is a program for splitting a stack into evenly partitioned substacks',& ! descr_long
        &'simple_exec',&                             ! executable
        &.false.)                                    ! requires sp_project
        ! INPUT PARAMETER SPECIFICATIONS
        ! image input/output
        call split_%add_input(UI_IMG, stk, required_override=.true.)
        ! parameter input/output
        call split_%add_input(UI_PARM, smpd)
        ! computer controls
        call split_%add_input(UI_COMP, nparts)
        ! alternative inputs
        ! <empty>
        ! search controls
        ! <empty>               
        ! filter controls
        ! <empty>
        ! mask controls
        ! <empty>
        ! add to ui_hash
        call add_ui_program('split', split_, prgtab)
    end subroutine new_split_

    subroutine new_split_stack( prgtab )
        class(ui_hash), intent(inout) :: prgtab
        ! PROGRAM SPECIFICATION
        call split_stack%new(&
        &'split_stack',&                                              ! name
        &'split stack in project',&                                   ! descr_short
        &'is a program for splitting a stack into nparts substacks',& ! descr_long
        &'simple_exec',&                                              ! executable
        &.true.)                                                      ! requires sp_project
        ! INPUT PARAMETER SPECIFICATIONS
        ! image input/output
        ! <empty>
        ! parameter input/output
        call split_stack%add_input(UI_PARM, 'nparts', 'num', 'Number of parts balanced splitting of the stack', '# parts', '# parts', .true., 1.0)
        ! alternative inputs
        ! <empty>
        ! search controls
        ! <empty>
        ! filter controls
        ! <empty>
        ! mask controls
        ! <empty>
        ! computer controls
        ! <empty>
        ! add to ui_hash
        call add_ui_program('split_stack', split_stack, prgtab)
    end subroutine new_split_stack

end module simple_ui_api_other
