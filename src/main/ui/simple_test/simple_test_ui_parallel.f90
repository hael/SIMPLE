!@descr: module defining the user interfaces for parallel  programs in the simple_test_exec suite
module simple_test_ui_parallel
use simple_ui_modules
implicit none

type(ui_program), target :: coarrays
type(ui_program), target :: openacc
type(ui_program), target :: openmp
type(ui_program), target :: simd

contains

    subroutine construct_parallel_programs(prgtab)
        class(ui_hash), intent(inout) :: prgtab
        call new_coarrays(prgtab)
        call new_openacc(prgtab)
        call new_openmp(prgtab)
        call new_simd(prgtab)
    end subroutine construct_parallel_programs

    subroutine print_parallel_programs(logfhandle)
        integer, intent(in) :: logfhandle
        write(logfhandle,'(A)') format_str('PARALLEL:', C_UNDERLINED)
        write(logfhandle,'(A)') coarrays%name%to_char()
        write(logfhandle,'(A)') openacc%name%to_char()
        write(logfhandle,'(A)') openmp%name%to_char()
        write(logfhandle,'(A)') simd%name%to_char()
        write(logfhandle,'(A)') ''
    end subroutine print_parallel_programs

    subroutine new_coarrays( prgtab )
        class(ui_hash), intent(inout) :: prgtab
        ! PROGRAM SPECIFICATION
        call coarrays%new(&
        &'coarrays',&                          ! name
        &'coarrays ',&                         ! descr_short
        &'is a test program for ',&
        &'simple_test_exec',&                  ! executable
        &.false.)                              ! requires sp_project
        ! INPUT PARAMETER SPECIFICATIONS
        ! image input/output
        !call coarrays%add_input(UI_IO, )
        ! parameter input/output
        !call coarrays%add_input(UI_IMG, )
        ! alternative inputs
        !call coarrays%add_input(UI_PARM, )
        ! search controls
        !call coarrays%add_input(UI_SRCH, )
        ! filter controls
        !call coarrays%add_input(UI_FILT, )
        ! mask controls
        !call coarrays%add_input(UI_MASK, )
        ! computer controls
        !call coarrays%add_input(UI_COMP, )
        ! add to ui_hash
        call add_ui_program('coarrays', coarrays, prgtab)
    end subroutine new_coarrays

    subroutine new_openacc( prgtab )
        class(ui_hash), intent(inout) :: prgtab
        ! PROGRAM SPECIFICATION
        call openacc%new(&
        &'openacc',&                           ! name
        &'openacc ',&                          ! descr_short
        &'is a test program for ',&
        &'simple_test_exec',&                  ! executable
        &.false.)                              ! requires sp_project
        ! INPUT PARAMETER SPECIFICATIONS
        ! image input/output
        !call openacc%add_input(UI_IO, )
        ! parameter input/output
        !call openacc%add_input(UI_IMG, )
        ! alternative inputs
        !call openacc%add_input(UI_PARM, )
        ! search controls
        !call openacc%add_input(UI_SRCH, )
        ! filter controls
        !call openacc%add_input(UI_FILT, )
        ! mask controls
        !call openacc%add_input(UI_MASK, )
        ! computer controls
        !call openacc%add_input(UI_COMP, )
        ! add to ui_hash
        call add_ui_program('openacc', openacc, prgtab)
    end subroutine new_openacc

    subroutine new_openmp( prgtab )
        class(ui_hash), intent(inout) :: prgtab
        ! PROGRAM SPECIFICATION
        call openmp%new(&
        &'openmp',&                            ! name
        &'openmp ',&                           ! descr_short
        &'is a test program for ',&
        &'simple_test_exec',&                  ! executable
        &.false.)                              ! requires sp_project
        ! INPUT PARAMETER SPECIFICATIONS
        ! image input/output
        !call openmp%add_input(UI_IO, )
        ! parameter input/output
        !call openmp%add_input(UI_IMG, )
        ! alternative inputs
        !call openmp%add_input(UI_PARM, )
        ! search controls
        !call openmp%add_input(UI_SRCH, )
        ! filter controls
        !call openmp%add_input(UI_FILT, )
        ! mask controls
        !call openmp%add_input(UI_MASK, )
        ! computer controls
        !call openmp%add_input(UI_COMP, )
        ! add to ui_hash
        call add_ui_program('openmp', openmp, prgtab)
    end subroutine new_openmp

    subroutine new_simd( prgtab )
        class(ui_hash), intent(inout) :: prgtab
        ! PROGRAM SPECIFICATION
        call simd%new(&
        &'simd',&                              ! name
        &'simd ',&                             ! descr_short
        &'is a test program for ',&
        &'simple_test_exec',&                  ! executable
        &.false.)                              ! requires sp_project
        ! INPUT PARAMETER SPECIFICATIONS
        ! image input/output
        !call simd%add_input(UI_IO, )
        ! parameter input/output
        !call simd%add_input(UI_IMG, )
        ! alternative inputs
        !call simd%add_input(UI_PARM, )
        ! search controls
        !call simd%add_input(UI_SRCH, )
        ! filter controls
        !call simd%add_input(UI_FILT, )
        ! mask controls
        !call simd%add_input(UI_MASK, )
        ! computer controls
        !call simd%add_input(UI_COMP, )
        ! add to ui_hash
        call add_ui_program('simd', simd, prgtab)
    end subroutine new_simd

end module simple_test_ui_parallel
