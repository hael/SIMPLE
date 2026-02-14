!@descr: "print" UI api (concrete implementation)
module simple_ui_api_print
use simple_ui_api_modules
implicit none

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

    subroutine print_print_programs(logfhandle)
        integer, intent(in) :: logfhandle
        write(logfhandle,'(A)') format_str('PRINT INFO:', C_UNDERLINED)
        write(logfhandle,'(A)') info_image%name%to_char()
        write(logfhandle,'(A)') info_stktab%name%to_char()
        write(logfhandle,'(A)') print_dose_weights%name%to_char()
        write(logfhandle,'(A)') print_fsc%name%to_char()
        write(logfhandle,'(A)') print_magic_boxes%name%to_char()
        write(logfhandle,'(A)') ''
    end subroutine print_print_programs

    subroutine new_info_image( prgtab ) 
        class(ui_hash), intent(inout) :: prgtab
        ! PROGRAM SPECIFICATION
        call info_image%new(&
        &'info_image', &                                                                       ! name
        &'Print header information',&                                                          ! descr_short
        &'is a program for printing header information in MRC and SPIDER stacks and volumes',& ! descr_long
        &'simple_exec',&                                                                       ! executable
        &.false.)                                                                              ! requires sp_project
        ! TEMPLATE
        ! INPUT PARAMETER SPECIFICATIONS
        ! image input/output
        call info_image%add_input(UI_IMG, 'fname', 'file', 'Name of image file', 'Name of image file', 'xxx.mrc file', .true., '')
        ! parameter input/output
        call info_image%add_input(UI_PARM, 'stats', 'binary', 'Output statistics', 'Output statistics(yes|no){no}',             '(yes|no){no}', .false., 'no')
        call info_image%add_input(UI_PARM, 'vis',   'binary', 'Visualize image',   'Visualize image with gnuplot(yes|no){yes}', '(yes|no){no}', .false., 'no')
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
        call add_ui_program('info_image', info_image, prgtab)
    end subroutine new_info_image

    subroutine new_info_stktab( prgtab )
        class(ui_hash), intent(inout) :: prgtab
        ! PROGRAM SPECIFICATION
        call info_stktab%new(&
        &'info_stktab', &                                                        ! name
        &'Print stktab information',&                                            ! descr_short
        &'is a program for printing information about stktab (list of stacks)',& ! descr_long
        &'simple_exec',&                                                         ! executable
        &.false.)                                                                ! requires sp_project
        ! TEMPLATE
        ! INPUT PARAMETER SPECIFICATIONS
        ! image input/output
        stktab%required = .true.
        call info_stktab%add_input(UI_IMG, stktab)
        ! parameter input/output
        ! <empty>
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
        call add_ui_program('info_stktab', info_stktab, prgtab)
    end subroutine new_info_stktab


    subroutine new_print_dose_weights( prgtab )
        class(ui_hash), intent(inout) :: prgtab
        ! PROGRAM SPECIFICATION
        call print_dose_weights%new(&
        &'print_dose_weights', &                                                  ! name
        &'Print dose weights used in motion correction',&                         ! descr_short
        &'is a program for printing the dose weights used in motion correction',& ! descr_long
        &'simple_exec',&                                                          ! executable
        &.false.)                                                                 ! requires sp_project
        ! INPUT PARAMETER SPECIFICATIONS
        ! image input/output
        ! <empty>
        ! parameter input/output
        call print_dose_weights%add_input(UI_PARM, smpd)
        call print_dose_weights%add_input(UI_PARM, box)
        call print_dose_weights%add_input(UI_PARM, 'nframes',   'num', 'Number of frames', 'Number of movie frames', '# frames', .true., 0.)
        call print_dose_weights%add_input(UI_PARM, kv)
        call print_dose_weights%add_input(UI_PARM, total_dose, required_override=.true.)
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
        call add_ui_program('print_dose_weights', print_dose_weights, prgtab)
    end subroutine new_print_dose_weights

    subroutine new_print_fsc( prgtab )
        class(ui_hash), intent(inout) :: prgtab
        ! PROGRAM SPECIFICATION
        call print_fsc%new(&
        &'print_fsc', &                                                          ! name
        &'Print FSC file produced by REFINE3D',&                                 ! descr_short
        &'is a program for printing the binary FSC files produced by REFINE3D',& ! descr_long
        &'simple_exec',&                                                         ! executable
        &.false.)                                                                ! requires sp_project
        ! INPUT PARAMETER SPECIFICATIONS
        ! image input/output
        ! <empty>
        ! parameter input/output
        call print_fsc%add_input(UI_PARM, smpd, required_override=.false.)
        call print_fsc%add_input(UI_PARM, box,  required_override=.false.)
        call print_fsc%add_input(UI_PARM, 'fsc', 'file', 'FSC file', 'Binary file with FSC info',&
        'input binary file e.g. fsc_state01.bin', .false., 'fsc_state01.bin')
        call print_fsc%add_input(UI_PARM, frcs)
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
        call add_ui_program('print_fsc', print_fsc, prgtab)
    end subroutine new_print_fsc

    subroutine new_print_magic_boxes( prgtab )
        class(ui_hash), intent(inout) :: prgtab
        ! PROGRAM SPECIFICATION
        call print_magic_boxes%new(&
        &'print_magic_boxes', &                                   ! name
        &'Print magic boxes (fast FFT)',&                         ! descr_short
        &'is a program for printing magic box sizes (fast FFT)',& ! descr_long
        &'simple_exec',&                                          ! executable
        &.false.)                                                 ! requires sp_project
        ! INPUT PARAMETER SPECIFICATIONS
        ! image input/output
        ! <empty>
        ! parameter input/output
        call print_magic_boxes%add_input(UI_PARM, smpd)
        call print_magic_boxes%add_input(UI_PARM, box)
        call print_magic_boxes%add_input(UI_PARM, moldiam)
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
        call add_ui_program('print_magic_boxes', print_magic_boxes, prgtab)
    end subroutine new_print_magic_boxes

end module simple_ui_api_print
