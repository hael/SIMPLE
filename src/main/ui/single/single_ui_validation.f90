!@descr: module defining the user interfaces for validation programs in the single_exec suite
module single_ui_validate
use simple_ui_modules
implicit none

type(ui_program), target :: cavgseoproc_nano
type(ui_program), target :: cavgsproc_nano
type(ui_program), target :: ptclsproc_nano

contains

    subroutine construct_single_validate_programs(prgtab)
        class(ui_hash), intent(inout) :: prgtab
        call new_cavgseoproc_nano(prgtab)
        call new_cavgsproc_nano(prgtab)
        call new_ptclsproc_nano(prgtab)
    end subroutine construct_single_validate_programs

    subroutine print_single_validate_programs(logfhandle)
        integer, intent(in) :: logfhandle
        write(logfhandle,'(A)') format_str('VALIDATION:', C_UNDERLINED)
        write(logfhandle,'(A)') cavgseoproc_nano%name%to_char()
        write(logfhandle,'(A)') cavgsproc_nano%name%to_char()
        write(logfhandle,'(A)') ptclsproc_nano%name%to_char()
        write(logfhandle,'(A)') ''
    end subroutine print_single_validate_programs

    subroutine new_cavgseoproc_nano( prgtab )
        class(ui_hash), intent(inout) :: prgtab
        ! PROGRAM SPECIFICATION
        call cavgseoproc_nano%new(&
        &'cavgseoproc_nano',&                                           ! name
        &'Analysis of even and odd class averages along nanocrystal time-series',& ! descr_short
        &'is a program to analyze the core/surface dynamics of nanocrystals using even and odd class averages',& ! descr_long
        &'single_exec',&                                                ! executable
        &.true., gui_advanced=.false.)                                  ! requires sp_project
        ! INPUT PARAMETER SPECIFICATIONS
        ! image input/output
        call cavgseoproc_nano%add_input(UI_IMG, 'vol1', 'file', 'Volume', 'Input volume', 'input volume e.g. vol.mrc', .true., '')
        ! parameter input/output
        call cavgseoproc_nano%add_input(UI_PARM, smpd)
        ! alternative inputs
        ! <empty>
        ! search controls
        call cavgseoproc_nano%add_input(UI_SRCH, pgrp)
        ! filter controls
        ! <empty>
        ! mask controls
        call cavgseoproc_nano%add_input(UI_MASK, mskdiam)
        ! computer controls
        call cavgseoproc_nano%add_input(UI_COMP, nthr)
        call cavgseoproc_nano%add_input(UI_COMP, script)
        call add_ui_program('cavgseoproc_nano', cavgseoproc_nano, prgtab)
    end subroutine new_cavgseoproc_nano

    subroutine new_cavgsproc_nano( prgtab )
        class(ui_hash), intent(inout) :: prgtab
        ! PROGRAM SPECIFICATION
        call cavgsproc_nano%new(&
        &'cavgsproc_nano',&                                           ! name
        &'Analysis of class averages along nanocrystal time-series',& ! descr_short
        &'is a program to analyze the core/surface dynamics of nanocrystals using class averages and re-projections',& ! descr_long
        &'single_exec',&                                              ! executable
        &.true., gui_advanced=.false.)                                ! requires sp_project
        ! INPUT PARAMETER SPECIFICATIONS
        ! image input/output
        call cavgsproc_nano%add_input(UI_IMG, 'vol1', 'file', 'Volume', 'Input volume', 'input volume e.g. vol.mrc', .true., '')
        ! parameter input/output
        call cavgsproc_nano%add_input(UI_PARM, smpd)
        ! alternative inputs
        ! <empty>
        ! search controls
        call cavgsproc_nano%add_input(UI_SRCH, pgrp)
        ! filter controls
        ! <empty>
        ! mask controls
        call cavgsproc_nano%add_input(UI_MASK, mskdiam)
        ! computer controls
        call cavgsproc_nano%add_input(UI_COMP, nthr)
        call cavgsproc_nano%add_input(UI_COMP, script)
        call add_ui_program('cavgsproc_nano', cavgsproc_nano, prgtab)
    end subroutine new_cavgsproc_nano

    subroutine new_ptclsproc_nano( prgtab )
        class(ui_hash), intent(inout) :: prgtab
        ! PROGRAM SPECIFICATION
        call ptclsproc_nano%new(&
        &'ptclsproc_nano',&                                           ! name
        &'Analysis of particle images inside a class along nanocrystal time-series using radial cross-correlation',& ! descr_short
        &'is a program to analyze the core/surface dynamics of nanocrystals using particle images inside a class',& ! descr_long
        &'single_exec',&                                              ! executable
        &.true., gui_advanced=.false.)                                ! requires sp_project
        ! INPUT PARAMETER SPECIFICATIONS
        ! image input/output
        ! parameter input/output
        call ptclsproc_nano%add_input(UI_PARM, smpd)
        ! alternative inputs
        ! <empty>
        ! search controls
        call ptclsproc_nano%add_input(UI_SRCH, pgrp)
        ! filter controls
        ! <empty>
        ! mask controls
        call ptclsproc_nano%add_input(UI_MASK, mskdiam)
        ! computer controls
        call ptclsproc_nano%add_input(UI_COMP, nthr)
        call ptclsproc_nano%add_input(UI_COMP, script)
        call add_ui_program('ptclsproc_nano', ptclsproc_nano, prgtab)
    end subroutine new_ptclsproc_nano

end module single_ui_validate
