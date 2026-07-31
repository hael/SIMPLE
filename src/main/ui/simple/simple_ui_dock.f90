!@descr: module defining the user interfaces for docking programs in the simple_exec suite
module simple_ui_dock
use simple_ui_modules
implicit none

type(category_descriptor), parameter :: UI_CATEGORY = category_descriptor('dock', 'Volume Docking', 170)
type(ui_program), target :: dock_volpair
type(ui_program), target :: volanalyze
type(ui_program), target :: volcluster

contains

    subroutine construct_dock_programs(prgtab)
        class(ui_hash), intent(inout) :: prgtab
        call new_dock_volpair(prgtab)
        call new_volanalyze(prgtab)
        call new_volcluster(prgtab)
    end subroutine construct_dock_programs

subroutine new_dock_volpair( prgtab )
        class(ui_hash), intent(inout) :: prgtab
        ! PROGRAM SPECIFICATION
        call dock_volpair%new(&
        &'dock_volpair', &                              ! name
        &'Align and dock two input volumes',& ! summary
        &'is a program for docking a pair of volumes',& ! descr long
        &'simple_exec',&                                ! executable
        &.false., &
        &visibility=UI_VIS_ADVANCED)                                       ! requires sp_project
        ! INPUT PARAMETER SPECIFICATIONS
        ! image input/output
        call dock_volpair%add_input(UI_IMG, 'vol1', 'file', 'Volume', 'Reference volume', &
        & 'input reference volume e.g. vol1.mrc', .true., '', &
        &visibility=UI_VIS_STANDARD)
        call dock_volpair%add_input(UI_IMG, 'vol2', 'file', 'Volume', 'Target volume', &
        & 'input target volume e.g. vol2.mrc', .true., '', &
        &visibility=UI_VIS_STANDARD)
        call dock_volpair%add_input(UI_IMG, outvol, &
        &visibility=UI_VIS_ADVANCED)
        ! parameter input/output
        call dock_volpair%add_input(UI_PARM, smpd, &
        &visibility=UI_VIS_STANDARD)
        call dock_volpair%add_input(UI_FILE, outfile, &
        &visibility=UI_VIS_ADVANCED)
        ! <no additional inputs>
        ! <empty>
        ! search controls
        call dock_volpair%add_input(UI_SRCH, trs, &
        &visibility=UI_VIS_ADVANCED)
        ! filter controls
        call dock_volpair%add_input(UI_FILT, hp, &
        &visibility=UI_VIS_ADVANCED)
        call dock_volpair%add_input(UI_FILT, lp, &
        &visibility=UI_VIS_ADVANCED)
        ! mask controls
        call dock_volpair%add_input(UI_MASK, mskdiam, &
        &visibility=UI_VIS_STANDARD)
        ! computer controls
        call dock_volpair%add_input(UI_COMP, nthr, &
        &visibility=UI_VIS_STANDARD)
        ! add to ui_hash
        call add_ui_program('dock_volpair', dock_volpair, prgtab, UI_CATEGORY)
    end subroutine new_dock_volpair

    subroutine new_volanalyze( prgtab )
        class(ui_hash), intent(inout) :: prgtab
        ! PROGRAM SPECIFICATION
        call volanalyze%new(&
        &'volanalyze',&                                                             ! name
        &'Analyze an emsemble of ab initio volumes',&                               ! summary
        &'is a program for statistical analysis an ensemble of ab initio volumes',& ! help
        &'simple_exec',&                                                            ! executable
        &.false., &
        &visibility=UI_VIS_ADVANCED)                                                                   ! requires sp_project
        ! INPUT PARAMETER SPECIFICATIONS
        ! image input/output
        call volanalyze%add_input(UI_FILE, 'filetab', 'file', 'Volumes list',&
        &'List of volumes to analyze', 'list input e.g. voltab.txt', .true., '', &
        &visibility=UI_VIS_STANDARD)
        ! parameter input/output
        call volanalyze%add_input(UI_PARM, smpd, &
        &visibility=UI_VIS_STANDARD)
        call volanalyze%add_input(UI_PARM, 'ref_ind', 'num', 'Reference volume index', 'Index of volume in voltab to use as reference', 'ref idx', .false., 0., &
        &visibility=UI_VIS_ADVANCED)
        ! <no additional inputs>
        ! <empty>
        ! search controls
        ! <empty>
        ! filter controls
        call volanalyze%add_input(UI_FILT, hp, required_override=.true., &
        &visibility=UI_VIS_STANDARD)
        call volanalyze%add_input(UI_FILT, lp, required_override=.true., &
        &visibility=UI_VIS_STANDARD)
        ! mask controls
        ! mask controls
        call volanalyze%add_input(UI_MASK, mskdiam, &
        &visibility=UI_VIS_STANDARD)
        ! computer controls
        call volanalyze%add_input(UI_COMP, nthr, &
        &visibility=UI_VIS_STANDARD)
        ! add to ui_hash
        call add_ui_program('volanalyze', volanalyze, prgtab, UI_CATEGORY)
    end subroutine new_volanalyze

    subroutine new_volcluster( prgtab )
        class(ui_hash), intent(inout) :: prgtab
        ! PROGRAM SPECIFICATION
        call volcluster%new(&
        &'volcluster',&                                                             ! name
        &'Cluster pre-docked volumes by affinity propagation',& ! summary
        &'is a program for affinity-propagation clustering of pre-docked volumes',& ! help
        &'simple_exec',&                                                            ! executable
        &.false., &
        &visibility=UI_VIS_ADVANCED)                                                                   ! requires sp_project
        ! INPUT PARAMETER SPECIFICATIONS
        ! image input/output
        call volcluster%add_input(UI_FILE, 'filetab', 'file', 'Volumes list',&
        &'List of docked volumes to cluster', 'list input e.g. voltab.txt', .true., '', &
        &visibility=UI_VIS_STANDARD)
        ! parameter input/output
        call volcluster%add_input(UI_PARM, smpd, &
        &visibility=UI_VIS_STANDARD)
        call volcluster%add_input(UI_FILE, outfile, &
        &visibility=UI_VIS_ADVANCED)
        call volcluster%add_input(UI_PARM, ncls, required_override=.false., &
        &visibility=UI_VIS_ADVANCED)
        ! <no additional inputs>
        ! <empty>
        ! search controls
        ! <empty>
        ! filter controls
        call volcluster%add_input(UI_FILT, hp, required_override=.false., &
        &visibility=UI_VIS_ADVANCED)
        call volcluster%add_input(UI_FILT, lp, required_override=.true., &
        &visibility=UI_VIS_STANDARD)
        ! mask controls
        call volcluster%add_input(UI_MASK, mskdiam, &
        &visibility=UI_VIS_STANDARD)
        ! computer controls
        call volcluster%add_input(UI_COMP, nthr, required_override=.false., &
        &visibility=UI_VIS_ADVANCED)
        ! add to ui_hash
        call add_ui_program('volcluster', volcluster, prgtab, UI_CATEGORY)
    end subroutine new_volcluster

end module simple_ui_dock
