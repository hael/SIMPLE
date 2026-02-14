!@descr: "dock" UI api (concrete implementation)
module simple_ui_api_dock
use simple_ui_api_modules
implicit none

type(ui_program), target :: dock_volpair
type(ui_program), target :: volanalyze

contains

    subroutine construct_dock_programs(prgtab)
        class(ui_hash), intent(inout) :: prgtab
        call new_dock_volpair(prgtab)
        call new_volanalyze(prgtab)
    end subroutine construct_dock_programs

    subroutine print_dock_programs(logfhandle)
        integer, intent(in) :: logfhandle
        write(logfhandle,'(A)') format_str('VOLUME DOCKING:', C_UNDERLINED)
        write(logfhandle,'(A)') dock_volpair%name%to_char()
        write(logfhandle,'(A)') volanalyze%name%to_char()
        write(logfhandle,'(A)') ''
    end subroutine print_dock_programs

    subroutine new_dock_volpair( prgtab )
        class(ui_hash), intent(inout) :: prgtab
        ! PROGRAM SPECIFICATION
        call dock_volpair%new(&
        &'dock_volpair', &                              ! name
        &'Dock a pair of volumes',&                     ! descr_short
        &'is a program for docking a pair of volumes',& ! descr long
        &'simple_exec',&                                ! executable
        &.false.)                                       ! requires sp_project
        ! INPUT PARAMETER SPECIFICATIONS
        ! image input/output
        call dock_volpair%add_input(UI_IMG, 'vol1', 'file', 'Volume', 'Reference volume', &
        & 'input reference volume e.g. vol1.mrc', .true., '')
        call dock_volpair%add_input(UI_IMG, 'vol2', 'file', 'Volume', 'Target volume', &
        & 'input target volume e.g. vol2.mrc', .true., '')
        call dock_volpair%add_input(UI_IMG, outvol)
        ! parameter input/output
        call dock_volpair%add_input(UI_PARM, smpd)
        call dock_volpair%add_input(UI_PARM, outfile)
        ! alternative inputs
        ! <empty>
        ! search controls
        call dock_volpair%add_input(UI_SRCH, trs)
        ! filter controls
        call dock_volpair%add_input(UI_FILT, hp)
        call dock_volpair%add_input(UI_FILT, lp)
        ! mask controls
        call dock_volpair%add_input(UI_MASK, mskdiam)
        ! computer controls
        call dock_volpair%add_input(UI_COMP, nthr)
        ! add to ui_hash
        call add_ui_program('dock_volpair', dock_volpair, prgtab)
    end subroutine new_dock_volpair

    subroutine new_volanalyze( prgtab )
        class(ui_hash), intent(inout) :: prgtab
        ! PROGRAM SPECIFICATION
        call volanalyze%new(&
        &'volanalyze',&                                                             ! name
        &'Analyze an emsemble of ab initio volumes',&                               ! descr_short
        &'is a program for statistical analysis an ensemble of ab initio volumes',& ! descr_long
        &'simple_exec',&                                                            ! executable
        &.false.)                                                                   ! requires sp_project
        ! INPUT PARAMETER SPECIFICATIONS
        ! image input/output
        call volanalyze%add_input(UI_IMG, 'filetab', 'file', 'Volumes list',&
        &'List of volumes to analyze', 'list input e.g. voltab.txt', .true., '')
        ! parameter input/output
        call volanalyze%add_input(UI_PARM, smpd)
        call volanalyze%add_input(UI_PARM, 'ref_ind', 'num', 'Reference volume index', 'Index of volume in voltab to use as reference', 'ref idx', .false., 0.)
        ! alternative inputs
        ! <empty>
        ! search controls
        ! <empty>
        ! filter controls
        call volanalyze%add_input(UI_FILT, hp, required_override=.true.)
        call volanalyze%add_input(UI_FILT, lp, required_override=.true.)
        ! mask controls
        ! mask controls
        call volanalyze%add_input(UI_MASK, mskdiam)
        ! computer controls
        call volanalyze%add_input(UI_COMP, nthr)
        ! add to ui_hash
        call add_ui_program('volanalyze', volanalyze, prgtab)
    end subroutine new_volanalyze

end module simple_ui_api_dock
