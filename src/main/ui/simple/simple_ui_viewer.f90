!@descr: module defining the user interfaces of the simple_viewer programs (select2D)
module simple_ui_viewer
use simple_ui_modules
implicit none

type(category_descriptor), parameter :: UI_CATEGORY = category_descriptor('viewer', 'Viewer', 10)
type(ui_program),          target    :: select2D

contains

    subroutine construct_viewer_programs(prgtab)
        class(ui_hash), intent(inout) :: prgtab
        call new_select2D(prgtab)
    end subroutine construct_viewer_programs

    subroutine new_select2D( prgtab )
        class(ui_hash), intent(inout) :: prgtab
        ! PROGRAM SPECIFICATION
        call select2D%new(&
        &'select2D',&                                                                     ! name
        &'Select class averages interactively and record the selection in the project',&  ! summary
        &'is a program that shows the class averages of a project in a window, where classes are selected '//&
        &'by clicking and sorted by resolution or population, with adjustable display. Nothing is written '//&
        &'until Save selection, which makes <n>_selection/ beside the numbered run folder that holds the input '//&
        &'project (or inside the folder that holds it, when that is not a run folder) with the project and the '//&
        &'selection applied, as program selection does (each particle takes the state of its class); the input '//&
        &'project is not changed. '//&
        &'Run as simple_viewer prg=select2D projfile=<project>.simple; needs a desktop session',&  ! help
        &'simple_viewer',&                                                                ! executable
        &.true., &
        &visibility=UI_VIS_STANDARD, &
        &display_name='Select Class Averages')                                            ! requires sp_project
        ! INPUT PARAMETER SPECIFICATIONS
        ! image input/output
        ! <empty>
        ! parameter input/output
        call select2D%add_input(UI_FILE, projfile, &
        &visibility=UI_VIS_STANDARD)
        ! search controls
        ! <empty>
        ! filter controls
        ! <empty>
        ! mask controls
        ! <empty>
        ! computer controls
        ! <empty>
        ! add to ui_hash
        call add_ui_program('select2D', select2D, prgtab, UI_CATEGORY)
    end subroutine new_select2D

end module simple_ui_viewer
