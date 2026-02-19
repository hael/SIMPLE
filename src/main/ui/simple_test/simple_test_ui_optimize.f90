!@descr: module defining the user interfaces for optimize test programs in the simple_test_exec suite
module simple_test_ui_optimize
use simple_ui_modules
implicit none

type(ui_program), target :: lbfgsb
type(ui_program), target :: lbfgsb_cosine
type(ui_program), target :: lplims
type(ui_program), target :: lpstages_test
type(ui_program), target :: opt_lp
type(ui_program), target :: tree_srch

contains

    subroutine construct_test_optimize_programs( tsttab )
        class(ui_hash), intent(inout) :: tsttab
        call new_lbfgsb(tsttab)
        call new_lbfgsb_cosine(tsttab)
        call new_lplims(tsttab)
        call new_lpstages_test(tsttab)
        call new_opt_lp(tsttab)
        call new_tree_srch(tsttab)
    end subroutine construct_test_optimize_programs

    subroutine print_test_optimize_programs(logfhandle)
        integer, intent(in) :: logfhandle
        write(logfhandle,'(A)') format_str('OPTIMIZE:', C_UNDERLINED)
        write(logfhandle,'(A)') lbfgsb%name%to_char()
        write(logfhandle,'(A)') lbfgsb_cosine%name%to_char()
        write(logfhandle,'(A)') lplims%name%to_char()
        write(logfhandle,'(A)') lpstages_test%name%to_char()
        write(logfhandle,'(A)') opt_lp%name%to_char()
        write(logfhandle,'(A)') tree_srch%name%to_char()
        write(logfhandle,'(A)') ''
    end subroutine print_test_optimize_programs

    subroutine new_lbfgsb( tsttab )
        class(ui_hash), intent(inout) :: tsttab
        ! PROGRAM SPECIFICATION
        call lbfgsb%new(&
        &'lbfgsb',&                         ! name
        &'lbfgsb ',&                        ! descr_short
        &'is a test program for ',&
        &'simple_test_exec',&               ! executable
        &.false.)                           ! requires sp_project
        ! INPUT PARAMETER SPECIFICATIONS
        ! image input/output
        !call lbfgsb%add_input(UI_IO, )
        ! parameter input/output
        !call lbfgsb%add_input(UI_IMG, )
        ! alternative inputs
        !call lbfgsb%add_input(UI_PARM, )
        ! search controls
        !call lbfgsb%add_input(UI_SRCH, )
        ! filter controls
        !call lbfgsb%add_input(UI_FILT, )
        ! mask controls
        !call lbfgsb%add_input(UI_MASK, )
        ! computer controls
        !call lbfgsb%add_input(UI_COMP, )
        ! add to ui_hash
        call add_ui_program('lbfgsb', lbfgsb, tsttab)
    end subroutine new_lbfgsb

    subroutine new_lbfgsb_cosine( tsttab )
        class(ui_hash), intent(inout) :: tsttab
        ! PROGRAM SPECIFICATION
        call lbfgsb_cosine%new(&
        &'lbfgsb_cosine',&                     ! name
        &'lbfgsb_cosine ',&                    ! descr_short
        &'is a test program for ',&
        &'simple_test_exec',&                  ! executable
        &.false.)                              ! requires sp_project
        ! INPUT PARAMETER SPECIFICATIONS
        ! image input/output
        !call lbfgsb_cosine%add_input(UI_IO, )
        ! parameter input/output
        !call lbfgsb_cosine%add_input(UI_IMG, )
        ! alternative inputs
        !call lbfgsb_cosine%add_input(UI_PARM, )
        ! search controls
        !call lbfgsb_cosine%add_input(UI_SRCH, )
        ! filter controls
        !call lbfgsb_cosine%add_input(UI_FILT, )
        ! mask controls
        !call lbfgsb_cosine%add_input(UI_MASK, )
        ! computer controls
        !call lbfgsb_cosine%add_input(UI_COMP, )
        ! add to ui_hash
        call add_ui_program('lbfgsb_cosine', lbfgsb_cosine, tsttab)
    end subroutine new_lbfgsb_cosine

    subroutine new_lplims( tsttab )
        class(ui_hash), intent(inout) :: tsttab
        ! PROGRAM SPECIFICATION
        call lplims%new(&
        &'lplims',&                            ! name
        &'lplims ',&                           ! descr_short
        &'is a test program for ',&
        &'simple_test_exec',&                  ! executable
        &.false.)                              ! requires sp_project
        ! INPUT PARAMETER SPECIFICATIONS
        ! image input/output
        !call lplims%add_input(UI_IO, )
        ! parameter input/output
        !call lplims%add_input(UI_IMG, )
        ! alternative inputs
        !call lplims%add_input(UI_PARM, )
        ! search controls
        !call lplims%add_input(UI_SRCH, )
        ! filter controls
        !call lplims%add_input(UI_FILT, )
        ! mask controls
        !call lplims%add_input(UI_MASK, )
        ! computer controls
        !call lplims%add_input(UI_COMP, )
        ! add to ui_hash
        call add_ui_program('lplims', lplims, tsttab)
    end subroutine new_lplims

    subroutine new_lpstages_test( tsttab )
        class(ui_hash), intent(inout) :: tsttab
        ! PROGRAM SPECIFICATION
        call lpstages_test%new(&
        &'lpstages_test',&                     ! name
        &'lpstages_test ',&                    ! descr_short
        &'is a test program for ',&
        &'simple_test_exec',&                  ! executable
        &.false.)                              ! requires sp_project
        ! INPUT PARAMETER SPECIFICATIONS
        ! image input/output
        !call lpstages_test%add_input(UI_IO, )
        ! parameter input/output
        !call lpstages_test%add_input(UI_IMG, )
        ! alternative inputs
        !call lpstages_test%add_input(UI_PARM, )
        ! search controls
        !call lpstages_test%add_input(UI_SRCH, )
        ! filter controls
        !call lpstages_test%add_input(UI_FILT, )
        ! mask controls
        !call lpstages_test%add_input(UI_MASK, )
        ! computer controls
        !call lpstages_test%add_input(UI_COMP, )
        ! add to ui_hash
        call add_ui_program('lpstages_test', lpstages_test, tsttab)
    end subroutine new_lpstages_test

    subroutine new_opt_lp( tsttab )
        class(ui_hash), intent(inout) :: tsttab
        ! PROGRAM SPECIFICATION
        call opt_lp%new(&
        &'opt_lp',&                            ! name
        &'opt_lp ',&                           ! descr_short
        &'is a test program for ',&
        &'simple_test_exec',&                  ! executable
        &.false.)                              ! requires sp_project
        ! INPUT PARAMETER SPECIFICATIONS
        ! image input/output
        !call opt_lp%add_input(UI_IO, )
        ! parameter input/output
        !call opt_lp%add_input(UI_IMG, )
        ! alternative inputs
        !call opt_lp%add_input(UI_PARM, )
        ! search controls
        !call opt_lp%add_input(UI_SRCH, )
        ! filter controls
        !call opt_lp%add_input(UI_FILT, )
        ! mask controls
        !call opt_lp%add_input(UI_MASK, )
        ! computer controls
        !call opt_lp%add_input(UI_COMP, )
        ! add to ui_hash
        call add_ui_program('opt_lp', opt_lp, tsttab)
    end subroutine new_opt_lp

    subroutine new_tree_srch( tsttab )
        class(ui_hash), intent(inout) :: tsttab
        ! PROGRAM SPECIFICATION
        call tree_srch%new(&
        &'tree_srch',&                         ! name
        &'tree_srch ',&                        ! descr_short
        &'is a test program for ',&
        &'simple_test_exec',&                  ! executable
        &.false.)                              ! requires sp_project
        ! INPUT PARAMETER SPECIFICATIONS
        ! image input/output
        !call tree_srch%add_input(UI_IO, )
        ! parameter input/output
        !call tree_srch%add_input(UI_IMG, )
        ! alternative inputs
        !call tree_srch%add_input(UI_PARM, )
        ! search controls
        !call tree_srch%add_input(UI_SRCH, )
        ! filter controls
        !call tree_srch%add_input(UI_FILT, )
        ! mask controls
        !call tree_srch%add_input(UI_MASK, )
        ! computer controls
        !call tree_srch%add_input(UI_COMP, )
        ! add to ui_hash
        call add_ui_program('tree_srch', tree_srch, tsttab)
    end subroutine new_tree_srch

end module simple_test_ui_optimize
