!@descr: "symmetry" UI api (concrete implementation)
module simple_ui_api_symmetry
use simple_ui_api_modules
implicit none

type(ui_program), target :: symaxis_search
type(ui_program), target :: symmetrize_map
type(ui_program), target :: symmetry_test

contains

    subroutine new_symaxis_search( prgtab )
        class(ui_hash), intent(inout) :: prgtab
        ! PROGRAM SPECIFICATION
        call symaxis_search%new(&
        &'symaxis_search',&                                                                                 ! name
        &'Search for symmetry axis',&                                                                       ! descr_short
        &'is a program for searching for the principal symmetry axis of a volume reconstructed in C1. &
        &The rotational transformation is applied to the oritype field in the project and the project &
        &file is updated. If you are unsure about the point-group, use the symmetry_test program instead',& ! descr_long
        &'simple_exec',&                                                                                    ! executable
        &.false.)                                                                                           ! requires sp_project
        ! INPUT PARAMETER SPECIFICATIONS
        ! image input/output
        call symaxis_search%add_input(UI_IMG, 'vol1', 'file', 'C1 Volume to identify symmetry axis of', 'C1 Volume to identify symmetry axis of', &
        & 'input volume e.g. vol_C1.mrc', .true., '')
        ! parameter input/output
        call symaxis_search%add_input(UI_PARM, smpd)
        ! alternative inputs
        ! <empty>
        ! search controls
        call symaxis_search%add_input(UI_SRCH, pgrp)
        call symaxis_search%add_input(UI_SRCH, 'center', 'binary', 'Center input volume', 'Center input volume by its &
        &center of gravity before symmetry axis search(yes|no){yes}', '(yes|no){yes}', .false., 'yes')
        ! filter controls
        call symaxis_search%add_input(UI_FILT, lp)
        call symaxis_search%add_input(UI_FILT, hp)
        call symaxis_search%add_input(UI_FILT, 'cenlp', 'num', 'Centering low-pass limit', 'Limit for low-pass filter used in binarisation &
        &prior to determination of the center of gravity of the input volume and centering', 'centering low-pass limit in &
        &Angstroms{30}', .false., 30.)
        ! mask controls
        call symaxis_search%add_input(UI_MASK, mskdiam)
        ! computer controls
        call symaxis_search%add_input(UI_COMP, nthr)
        ! add to ui_hash
        call add_ui_program('symaxis_search', symaxis_search, prgtab)
    end subroutine new_symaxis_search

    subroutine new_symmetrize_map( prgtab )
        class(ui_hash), intent(inout) :: prgtab
        ! PROGRAM SPECIFICATION
        call symmetrize_map%new(&
        &'symmetrize_map',&                                                                                          ! name
        &'Symmetrization of density map',&                                                                           ! descr_short
        &'is a program that implements symmetrization of the input density map. &
        &Input is a volume and point-group symmetry, output is the volume aligned to the principal symmetry axis and averaged over the symmetry operations',& ! descr long
        &'simple_exec',&                                                                                             ! executable
        &.false.)                                                                                                    ! requires sp_project
        ! INPUT PARAMETER SPECIFICATIONS
        ! image input/output
        call symmetrize_map%add_input(UI_IMG, 'vol1', 'file', 'Volume to symmetrize', 'Volume to symmetrize', &
        & 'input volume e.g. vol.mrc', .true., '')
        call symmetrize_map%add_input(UI_IMG, outvol)
        ! parameter input/output
        call symmetrize_map%add_input(UI_PARM, smpd)
        ! alternative inputs
        ! <empty>
        ! search controls
        call symmetrize_map%add_input(UI_SRCH, pgrp)
        call symmetrize_map%add_input(UI_SRCH, 'center', 'binary', 'Center input volume', 'Center input volume by its &
        &center of gravity before symmetry axis search(yes|no){yes}', '(yes|no){yes}', .false., 'yes')
        ! filter controls
        call symmetrize_map%add_input(UI_FILT, lp)
        call symmetrize_map%add_input(UI_FILT, hp)
        call symmetrize_map%add_input(UI_FILT, 'cenlp', 'num', 'Centering low-pass limit', 'Limit for low-pass filter used in binarisation &
        &prior to determination of the center of gravity of the input volume and centering', 'centering low-pass limit in &
        &Angstroms{30}', .false., 20.)
        ! mask controls
        call symmetrize_map%add_input(UI_MASK, mskdiam)
        ! computer controls
        call symmetrize_map%add_input(UI_COMP, nthr)
        ! add to ui_hash
        call add_ui_program('symmetrize_map', symmetrize_map, prgtab)
    end subroutine new_symmetrize_map

    subroutine new_symmetry_test( prgtab )
        class(ui_hash), intent(inout) :: prgtab
        ! PROGRAM SPECIFICATION
        call symmetry_test%new(&
        &'symmetry_test',&                                                                                           ! name
        &'Statistical test for symmetry',&                                                                           ! descr_short
        &'is a program that implements a statistical test for point-group symmetry. &
        &Input is a volume reconstructed without symmetry (c1) and output is the most likely point-group symmetry',& ! descr long
        &'simple_exec',&                                                                                             ! executable
        &.false.)                                                                                                    ! requires sp_project
        ! INPUT PARAMETER SPECIFICATIONS
        ! image input/output
        call symmetry_test%add_input(UI_IMG, 'vol1', 'file', 'C1 Volume to identify symmetry of', 'C1 Volume to identify symmetry of', &
        & 'input volume e.g. vol_C1.mrc', .true., '')
        ! parameter input/output
        call symmetry_test%add_input(UI_PARM, smpd)
        ! alternative inputs
        ! <empty>
        ! search controls
        call symmetry_test%add_input(UI_SRCH, 'cn_stop',  'num', 'Rotational symmetry order stop index',  'Rotational symmetry order stop index',  'give stop index',  .false., 10.)
        call symmetry_test%add_input(UI_SRCH, 'center', 'binary', 'Center input volume', 'Center input volume by its &
        &center of gravity before symmetry axis search(yes|no){yes}', '(yes|no){yes}', .false., 'yes')
        call symmetry_test%add_input(UI_SRCH, 'platonic', 'binary', 'Search for Platonic symmetries', 'Search for Platonic symmetries(yes|no){yes}', '(yes|no){yes}', .false., 'yes')
        ! filter controls
        call symmetry_test%add_input(UI_FILT, lp)
        call symmetry_test%add_input(UI_FILT, hp)
        call symmetry_test%add_input(UI_FILT, 'cenlp', 'num', 'Centering low-pass limit', 'Limit for low-pass filter used in binarisation &
        &prior to determination of the center of gravity of the input volume and centering', 'centering low-pass limit in &
        &Angstroms{30}', .false., 30.)
        ! mask controls
        call symmetry_test%add_input(UI_MASK, mskdiam)
        ! computer controls
        call symmetry_test%add_input(UI_COMP, nthr)
        ! add to ui_hash
        call add_ui_program('symmetry_test', symmetry_test, prgtab)
    end subroutine new_symmetry_test

end module simple_ui_api_symmetry
