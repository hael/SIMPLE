!@ descr: module defining the user interfaces for numerics  programs in the simple_test_exec suite
module simple_test_ui_numerics
use simple_ui_modules
implicit none

type(ui_program), target :: eigh_test
type(ui_program), target :: kbinterpol_fast
type(ui_program), target :: maxnloc_test
type(ui_program), target :: neigh

contains

    subroutine construct_numerics_programs(prgtab)
        class(ui_hash), intent(inout) :: prgtab
        call new_eigh_test(prgtab)
        call new_kbinterpol_fast(prgtab)
        call new_maxnloc_test(prgtab)
        call new_neigh(prgtab)
    end subroutine construct_numerics_programs

    subroutine print_numerics_programs(logfhandle)
        integer, intent(in) :: logfhandle
        write(logfhandle,'(A)') format_str('NUMERICS:', C_UNDERLINED)
        write(logfhandle,'(A)') eigh_test%name%to_char()
        write(logfhandle,'(A)') kbinterpol_fast%name%to_char()
        write(logfhandle,'(A)') maxnloc_test%name%to_char()
        write(logfhandle,'(A)') neigh%name%to_char()
        write(logfhandle,'(A)') ''
    end subroutine print_numerics_programs

    subroutine new_eigh_test( prgtab )
        class(ui_hash), intent(inout) :: prgtab
        ! PROGRAM SPECIFICATION
        call eigh_test%new(&
        &'eigh_test',&                         ! name
        &'eigh_test ',&                        ! descr_short
        &'is a test program for ',&
        &'simple_test_exec',&                  ! executable
        &.false.)                              ! requires sp_project
        ! INPUT PARAMETER SPECIFICATIONS
        ! image input/output
        !call eigh_test%add_input(UI_IO, )
        ! parameter input/output
        !call eigh_test%add_input(UI_IMG, )
        ! alternative inputs
        !call eigh_test%add_input(UI_PARM, )
        ! search controls
        !call eigh_test%add_input(UI_SRCH, )
        ! filter controls
        !call eigh_test%add_input(UI_FILT, )
        ! mask controls
        !call eigh_test%add_input(UI_MASK, )
        ! computer controls
        !call eigh_test%add_input(UI_COMP, )
        ! add to ui_hash
        call add_ui_program('eigh_test', eigh_test, prgtab)
    end subroutine new_eigh_test

    subroutine new_kbinterpol_fast( prgtab )
        class(ui_hash), intent(inout) :: prgtab
        ! PROGRAM SPECIFICATION
        call kbinterpol_fast%new(&
        &'kbinterpol_fast',&                         ! name
        &'kbinterpol_fast ',&                        ! descr_short
        &'is a test program for ',&
        &'simple_test_exec',&                  ! executable
        &.false.)                              ! requires sp_project
        ! INPUT PARAMETER SPECIFICATIONS
        ! image input/output
        !call kbinterpol_fast%add_input(UI_IO, )
        ! parameter input/output
        !call kbinterpol_fast%add_input(UI_IMG, )
        ! alternative inputs
        !call kbinterpol_fast%add_input(UI_PARM, )
        ! search controls
        !call kbinterpol_fast%add_input(UI_SRCH, )
        ! filter controls
        !call kbinterpol_fast%add_input(UI_FILT, )
        ! mask controls
        !call kbinterpol_fast%add_input(UI_MASK, )
        ! computer controls
        !call kbinterpol_fast%add_input(UI_COMP, )
        ! add to ui_hash
        call add_ui_program('kbinterpol_fast', kbinterpol_fast, prgtab)
    end subroutine new_kbinterpol_fast

    subroutine new_maxnloc_test( prgtab )
        class(ui_hash), intent(inout) :: prgtab
        ! PROGRAM SPECIFICATION
        call maxnloc_test%new(&
        &'maxnloc_test',&                         ! name
        &'maxnloc_test ',&                        ! descr_short
        &'is a test program for ',&
        &'simple_test_exec',&                  ! executable
        &.false.)                              ! requires sp_project
        ! INPUT PARAMETER SPECIFICATIONS
        ! image input/output
        !call maxnloc_test%add_input(UI_IO, )
        ! parameter input/output
        !call maxnloc_test%add_input(UI_IMG, )
        ! alternative inputs
        !call maxnloc_test%add_input(UI_PARM, )
        ! search controls
        !call maxnloc_test%add_input(UI_SRCH, )
        ! filter controls
        !call maxnloc_test%add_input(UI_FILT, )
        ! mask controls
        !call maxnloc_test%add_input(UI_MASK, )
        ! computer controls
        !call maxnloc_test%add_input(UI_COMP, )
        ! add to ui_hash
        call add_ui_program('maxnloc_test', maxnloc_test, prgtab)
    end subroutine new_maxnloc_test

    subroutine new_neigh( prgtab )
        class(ui_hash), intent(inout) :: prgtab
        ! PROGRAM SPECIFICATION
        call neigh%new(&
        &'neigh',&                         ! name
        &'neigh ',&                        ! descr_short
        &'is a test program for ',&
        &'simple_test_exec',&                  ! executable
        &.false.)                              ! requires sp_project
        ! INPUT PARAMETER SPECIFICATIONS
        ! image input/output
        !call neigh%add_input(UI_IO, )
        ! parameter input/output
        !call neigh%add_input(UI_IMG, )
        ! alternative inputs
        !call neigh%add_input(UI_PARM, )
        ! search controls
        !call neigh%add_input(UI_SRCH, )
        ! filter controls
        !call neigh%add_input(UI_FILT, )
        ! mask controls
        !call neigh%add_input(UI_MASK, )
        ! computer controls
        !call neigh%add_input(UI_COMP, )
        ! add to ui_hash
        call add_ui_program('neigh', neigh, prgtab)
    end subroutine new_neigh

end module simple_test_ui_numerics
