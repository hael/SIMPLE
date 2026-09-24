!@descr: module defining the user interfaces for parallel test programs in the simple_test_exec suite
module simple_test_ui_parallel
use simple_ui_modules
implicit none

type(category_descriptor), parameter :: UI_CATEGORY = category_descriptor('parallel', 'Parallel', 100)
type(ui_program), target :: coarrays

contains

    subroutine construct_test_parallel_programs( tsttab )
        class(ui_hash), intent(inout) :: tsttab
        call new_coarrays(tsttab)
    end subroutine construct_test_parallel_programs

subroutine new_coarrays( tsttab )
        class(ui_hash), intent(inout) :: tsttab
        ! PROGRAM SPECIFICATION
        call coarrays%new(&
        &'coarrays',&                          ! name
        &'coarrays ',&                         ! summary
        &'is a test program for ',&
        &'simple_test_exec',&                  ! executable
        &.false.)                              ! requires sp_project
        ! INPUT PARAMETER SPECIFICATIONS
        ! image input/output
        !call coarrays%add_input(UI_IO, )
        ! parameter input/output
        !call coarrays%add_input(UI_IMG, )
        ! <no additional inputs>
        !call coarrays%add_input(UI_PARM, )
        ! search controls
        !call coarrays%add_input(UI_SRCH, )
        ! filter controls
        !call coarrays%add_input(UI_FILT, )
        ! mask controls
        !call coarrays%add_input(UI_MASK, )
        ! computer controls
        call coarrays%add_input(UI_COMP, nparts, required_override=.false.)
        call coarrays%add_input(UI_COMP, 'ncunits', 'num', 'Number of coarray images',&
        &'Number of coarray images to launch concurrently; defaults to nparts',&
        &'# coarray images', .false., 0.)
        ! add to ui_hash
        call add_ui_program('coarrays', coarrays, tsttab, UI_CATEGORY)
    end subroutine new_coarrays

end module simple_test_ui_parallel
