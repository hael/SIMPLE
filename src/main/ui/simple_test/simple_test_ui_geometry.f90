!@ descr: module defining the user interfaces for geometry  programs in the simple_test_exec suite
module simple_test_ui_geometry
use simple_ui_modules
implicit none

type(ui_program), target :: angres
type(ui_program), target :: ori_test
type(ui_program), target :: oris_test
type(ui_program), target :: sym_test
type(ui_program), target :: uniform_euler
type(ui_program), target :: uniform_rot

contains

    subroutine construct_geometry_programs(prgtab)
        class(ui_hash), intent(inout) :: prgtab
        call new_angres(prgtab)
        call new_ori_test(prgtab)
        call new_oris_test(prgtab)
        call new_sym_test(prgtab)
        call new_uniform_euler(prgtab)
        call new_uniform_rot(prgtab)
    end subroutine construct_geometry_programs

    subroutine print_geometry_programs(logfhandle)
        integer, intent(in) :: logfhandle
        write(logfhandle,'(A)') format_str('GEOMETRY:', C_UNDERLINED)
        write(logfhandle,'(A)') angres%name%to_char()
        write(logfhandle,'(A)') ori_test%name%to_char()
        write(logfhandle,'(A)') oris_test%name%to_char()
        write(logfhandle,'(A)') sym_test%name%to_char()
        write(logfhandle,'(A)') uniform_euler%name%to_char()
        write(logfhandle,'(A)') uniform_rot%name%to_char()
        write(logfhandle,'(A)') ''
    end subroutine print_geometry_programs

    subroutine new_angres( prgtab )
        class(ui_hash), intent(inout) :: prgtab
        ! PROGRAM SPECIFICATION
        call angres%new(&
        &'angres',&                         ! name
        &'angres ',&                        ! descr_short
        &'is a test program for ',&
        &'simple_test_exec',&                  ! executable
        &.false.)                              ! requires sp_project
        ! INPUT PARAMETER SPECIFICATIONS
        ! image input/output
        !call angres%add_input(UI_IO, )
        ! parameter input/output
        !call angres%add_input(UI_IMG, )
        ! alternative inputs
        !call angres%add_input(UI_PARM, )
        ! search controls
        !call angres%add_input(UI_SRCH, )
        ! filter controls
        !call angres%add_input(UI_FILT, )
        ! mask controls
        !call angres%add_input(UI_MASK, )
        ! computer controls
        !call angres%add_input(UI_COMP, )
        ! add to ui_hash
        call add_ui_program('angres', angres, prgtab)
    end subroutine new_angres

    subroutine new_ori_test( prgtab )
        class(ui_hash), intent(inout) :: prgtab
        ! PROGRAM SPECIFICATION
        call ori_test%new(&
        &'ori_test',&                         ! name
        &'ori_test ',&                        ! descr_short
        &'is a test program for ',&
        &'simple_test_exec',&                  ! executable
        &.false.)                              ! requires sp_project
        ! INPUT PARAMETER SPECIFICATIONS
        ! image input/output
        !call ori_test%add_input(UI_IO, )
        ! parameter input/output
        !call ori_test%add_input(UI_IMG, )
        ! alternative inputs
        !call ori_test%add_input(UI_PARM, )
        ! search controls
        !call ori_test%add_input(UI_SRCH, )
        ! filter controls
        !call ori_test%add_input(UI_FILT, )
        ! mask controls
        !call ori_test%add_input(UI_MASK, )
        ! computer controls
        !call ori_test%add_input(UI_COMP, )
        ! add to ui_hash
        call add_ui_program('ori_test', ori_test, prgtab)
    end subroutine new_ori_test

    subroutine new_oris_test( prgtab )
        class(ui_hash), intent(inout) :: prgtab
        ! PROGRAM SPECIFICATION
        call oris_test%new(&
        &'oris_test',&                         ! name
        &'oris_test ',&                        ! descr_short
        &'is a test program for ',&
        &'simple_test_exec',&                  ! executable
        &.false.)                              ! requires sp_project
        ! INPUT PARAMETER SPECIFICATIONS
        ! image input/output
        !call oris_test%add_input(UI_IO, )
        ! parameter input/output
        !call oris_test%add_input(UI_IMG, )
        ! alternative inputs
        !call oris_test%add_input(UI_PARM, )
        ! search controls
        !call oris_test%add_input(UI_SRCH, )
        ! filter controls
        !call oris_test%add_input(UI_FILT, )
        ! mask controls
        !call oris_test%add_input(UI_MASK, )
        ! computer controls
        !call oris_test%add_input(UI_COMP, )
        ! add to ui_hash
        call add_ui_program('oris_test', oris_test, prgtab)
    end subroutine new_oris_test

    subroutine new_sym_test( prgtab )
        class(ui_hash), intent(inout) :: prgtab
        ! PROGRAM SPECIFICATION
        call sym_test%new(&
        &'sym_test',&                         ! name
        &'sym_test ',&                        ! descr_short
        &'is a test program for ',&
        &'simple_test_exec',&                  ! executable
        &.false.)                              ! requires sp_project
        ! INPUT PARAMETER SPECIFICATIONS
        ! image input/output
        !call sym_test%add_input(UI_IO, )
        ! parameter input/output
        !call sym_test%add_input(UI_IMG, )
        ! alternative inputs
        !call sym_test%add_input(UI_PARM, )
        ! search controls
        !call sym_test%add_input(UI_SRCH, )
        ! filter controls
        !call sym_test%add_input(UI_FILT, )
        ! mask controls
        !call sym_test%add_input(UI_MASK, )
        ! computer controls
        !call sym_test%add_input(UI_COMP, )
        ! add to ui_hash
        call add_ui_program('sym_test', sym_test, prgtab)
    end subroutine new_sym_test

    subroutine new_uniform_euler( prgtab )
        class(ui_hash), intent(inout) :: prgtab
        ! PROGRAM SPECIFICATION
        call uniform_euler%new(&
        &'uniform_euler',&                         ! name
        &'uniform_euler ',&                        ! descr_short
        &'is a test program for ',&
        &'simple_test_exec',&                  ! executable
        &.false.)                              ! requires sp_project
        ! INPUT PARAMETER SPECIFICATIONS
        ! image input/output
        !call uniform_euler%add_input(UI_IO, )
        ! parameter input/output
        !call uniform_euler%add_input(UI_IMG, )
        ! alternative inputs
        !call uniform_euler%add_input(UI_PARM, )
        ! search controls
        !call uniform_euler%add_input(UI_SRCH, )
        ! filter controls
        !call uniform_euler%add_input(UI_FILT, )
        ! mask controls
        !call uniform_euler%add_input(UI_MASK, )
        ! computer controls
        !call uniform_euler%add_input(UI_COMP, )
        ! add to ui_hash
        call add_ui_program('uniform_euler', uniform_euler, prgtab)
    end subroutine new_uniform_euler

    subroutine new_uniform_rot( prgtab )
        class(ui_hash), intent(inout) :: prgtab
        ! PROGRAM SPECIFICATION
        call uniform_rot%new(&
        &'uniform_rot',&                         ! name
        &'uniform_rot ',&                        ! descr_short
        &'is a test program for ',&
        &'simple_test_exec',&                  ! executable
        &.false.)                              ! requires sp_project
        ! INPUT PARAMETER SPECIFICATIONS
        ! image input/output
        !call uniform_rot%add_input(UI_IO, )
        ! parameter input/output
        !call uniform_rot%add_input(UI_IMG, )
        ! alternative inputs
        !call uniform_rot%add_input(UI_PARM, )
        ! search controls
        !call uniform_rot%add_input(UI_SRCH, )
        ! filter controls
        !call uniform_rot%add_input(UI_FILT, )
        ! mask controls
        !call uniform_rot%add_input(UI_MASK, )
        ! computer controls
        !call uniform_rot%add_input(UI_COMP, )
        ! add to ui_hash
        call add_ui_program('uniform_rot', uniform_rot, prgtab)
    end subroutine new_uniform_rot

end module simple_test_ui_geometry
