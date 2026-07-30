!@descr: module defining the user interfaces for printing programs in the simple_exec suite
module simple_ui_print
use simple_ui_modules
implicit none

type(category_descriptor), parameter :: UI_CATEGORY = category_descriptor('print', 'Print Information', 120)
type(ui_program), target :: info_image
type(ui_program), target :: info_stktab
type(ui_program), target :: print_dose_weights
type(ui_program), target :: print_fsc
type(ui_program), target :: print_magic_boxes

contains

    subroutine construct_print_programs(prgtab)
        class(ui_hash), intent(inout) :: prgtab
        call new_info_image(prgtab)
        call new_info_stktab(prgtab)
        call new_print_dose_weights(prgtab)
        call new_print_fsc(prgtab)
        call new_print_magic_boxes(prgtab)
    end subroutine construct_print_programs

subroutine new_info_image( prgtab )
        class(ui_hash), intent(inout) :: prgtab
        ! PROGRAM SPECIFICATION
        call info_image%new(&
        &'info_image', &                                                                       ! name
        &'Display header metadata for an MRC or SPIDER image file',& ! summary
        &'is a program for printing header information in MRC and SPIDER stacks and volumes',& ! help
        &'simple_exec',&                                                                       ! executable
        &.false., &
        &visibility=UI_VIS_STANDARD, display_name='Inspect Image Metadata') ! requires sp_project
        ! TEMPLATE
        ! INPUT PARAMETER SPECIFICATIONS
        ! image input/output
        call info_image%add_input(UI_FILE, 'fname', 'file', 'Name of image file', 'Name of image file', 'e.g. image.mrc', .true., '', &
        &visibility=UI_VIS_STANDARD)
        ! parameter input/output
        call info_image%add_input(UI_PARM, 'stats', 'binary', 'Output statistics', 'Output statistics(yes|no){no}','', .false., 'no', &
        &choices=ui_choices([character(len=3) :: 'yes', 'no']), &
        &visibility=UI_VIS_ADVANCED)
        call info_image%add_input(UI_PARM, 'vis',   'binary', 'Visualize image',   'Visualize image with gnuplot(yes|no){yes}','', .false., 'no', &
        &choices=ui_choices([character(len=3) :: 'yes', 'no']), &
        &visibility=UI_VIS_ADVANCED)
        ! <no additional inputs>
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
        call add_ui_program('info_image', info_image, prgtab, UI_CATEGORY)
    end subroutine new_info_image

    subroutine new_info_stktab( prgtab )
        class(ui_hash), intent(inout) :: prgtab
        ! PROGRAM SPECIFICATION
        call info_stktab%new(&
        &'info_stktab', &                                                        ! name
        &'Display the particle-stack entries listed in a stack-table file',& ! summary
        &'is a program for printing information about stktab (list of stacks)',& ! help
        &'simple_exec',&                                                         ! executable
        &.false., &
        &visibility=UI_VIS_STANDARD, display_name='Inspect Stack Table') ! requires sp_project
        ! TEMPLATE
        ! INPUT PARAMETER SPECIFICATIONS
        ! image input/output
        stktab%required = .true.
        call info_stktab%add_input(UI_FILE, stktab, &
        &visibility=UI_VIS_STANDARD)
        ! parameter input/output
        ! <empty>
        ! <no additional inputs>
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
        call add_ui_program('info_stktab', info_stktab, prgtab, UI_CATEGORY)
    end subroutine new_info_stktab

    subroutine new_print_dose_weights( prgtab )
        class(ui_hash), intent(inout) :: prgtab
        ! PROGRAM SPECIFICATION
        call print_dose_weights%new(&
        &'print_dose_weights', &                                                  ! name
        &'Print dose weights used in motion correction',&                         ! summary
        &'is a program for printing the dose weights used in motion correction',& ! help
        &'simple_exec',&                                                          ! executable
        &.false., &
        &visibility=UI_VIS_ADVANCED)                                                                 ! requires sp_project
        ! INPUT PARAMETER SPECIFICATIONS
        ! image input/output
        ! <empty>
        ! parameter input/output
        call print_dose_weights%add_input(UI_PARM, smpd, &
        &visibility=UI_VIS_STANDARD)
        call print_dose_weights%add_input(UI_PARM, box, &
        &visibility=UI_VIS_STANDARD)
        call print_dose_weights%add_input(UI_PARM, 'nframes',   'num', 'Number of frames', 'Number of movie frames', '# frames', .true., 0., &
        &visibility=UI_VIS_STANDARD)
        call print_dose_weights%add_input(UI_PARM, kv, &
        &visibility=UI_VIS_ADVANCED)
        call print_dose_weights%add_input(UI_PARM, total_dose, required_override=.true., &
        &visibility=UI_VIS_STANDARD)
        ! <no additional inputs>
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
        call add_ui_program('print_dose_weights', print_dose_weights, prgtab, UI_CATEGORY)
    end subroutine new_print_dose_weights

    subroutine new_print_fsc( prgtab )
        class(ui_hash), intent(inout) :: prgtab
        ! PROGRAM SPECIFICATION
        call print_fsc%new(&
        &'print_fsc', &                                                          ! name
        &'Print FSC file produced by REFINE3D',&                                 ! summary
        &'is a program for printing the binary FSC files produced by REFINE3D',& ! help
        &'simple_exec',&                                                         ! executable
        &.false., &
        &visibility=UI_VIS_ADVANCED)                                                                ! requires sp_project
        ! INPUT PARAMETER SPECIFICATIONS
        ! image input/output
        ! <empty>
        ! parameter input/output
        call print_fsc%add_input(UI_PARM, smpd, required_override=.false., &
        &visibility=UI_VIS_ADVANCED)
        call print_fsc%add_input(UI_PARM, box,  required_override=.false., &
        &visibility=UI_VIS_ADVANCED)
        call print_fsc%add_input(UI_FILE, 'fsc', 'file', 'FSC file', 'Binary file with FSC info',&
        'input binary file e.g. fsc_state01.bin', .false., 'fsc_state01.bin', &
        &visibility=UI_VIS_ADVANCED)
        call print_fsc%add_input(UI_PARM, frcs, &
        &visibility=UI_VIS_ADVANCED)
        ! <no additional inputs>
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
        call add_ui_program('print_fsc', print_fsc, prgtab, UI_CATEGORY)
    end subroutine new_print_fsc

    subroutine new_print_magic_boxes( prgtab )
        class(ui_hash), intent(inout) :: prgtab
        ! PROGRAM SPECIFICATION
        call print_magic_boxes%new(&
        &'print_magic_boxes', &                                   ! name
        &'List FFT-efficient image sizes for a requested range',& ! summary
        &'is a program for printing magic box sizes (fast FFT)',& ! help
        &'simple_exec',&                                          ! executable
        &.false.)                                                 ! requires sp_project
        ! INPUT PARAMETER SPECIFICATIONS
        ! image input/output
        ! <empty>
        ! parameter input/output
        call print_magic_boxes%add_input(UI_PARM, smpd, &
        &visibility=UI_VIS_STANDARD)
        call print_magic_boxes%add_input(UI_PARM, box, &
        &visibility=UI_VIS_STANDARD)
        call print_magic_boxes%add_input(UI_PARM, moldiam, &
        &visibility=UI_VIS_DEVELOPER)
        ! <no additional inputs>
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
        call add_ui_program('print_magic_boxes', print_magic_boxes, prgtab, UI_CATEGORY)
    end subroutine new_print_magic_boxes

end module simple_ui_print
