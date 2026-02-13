!@descr: "ori" UI api (concrete implementation)
module simple_ui_api_ori
use simple_ui_api_modules
implicit none

type(ui_program), target :: make_oris
type(ui_program), target :: orisops
type(ui_program), target :: oristats
type(ui_program), target :: vizoris

contains

    subroutine new_make_oris( prgtab )
        class(ui_hash), intent(inout) :: prgtab
        ! PROGRAM SPECIFICATION
        call make_oris%new(&
        &'make_oris',&                       ! name
        &'Make orientations',&               ! descr_short
        &'is a program for making SIMPLE orientation files. Random Euler angles and random origin shifts are generated.&
        & If ndiscrete is set to an integer number > 0, the orientations produced are randomly sampled from the set of&
        & ndiscrete quasi-even projection directions, and the in-plane parameters are assigned randomly. If even is set&
        & to yes, then all nptcls orientations are assigned quasi-even projection directions and  random in-plane parameters.&
        & If nstates is set to some integer number > 0, then states are assigned randomly',&
        &'simple_exec',&                     ! executable
        &.false.)                            ! requires sp_project
        ! INPUT PARAMETER SPECIFICATIONS
        ! image input/output
        ! <empty>
        ! parameter input/output
        call make_oris%add_input(UI_PARM, 'nptcls', 'num', 'Number of per-particle orientations', 'Number of per-particle orientations to produce', '# per-ptcl oris', .true., 1.0)
        call make_oris%add_input(UI_PARM, 'ncls', 'num', 'Number of random class labels', 'Number of random class labels to produce', '# classes', .false., 0.)
        call make_oris%add_input(UI_PARM, outfile)
        call make_oris%add_input(UI_PARM, 'nstates', 'num', 'Number of random state labels', 'Number of random state labels to produce', '# states', .false., 0.0)
        call make_oris%add_input(UI_PARM, pgrp, required_override=.false.)
        call make_oris%add_input(UI_PARM, sherr)
        call make_oris%add_input(UI_PARM, angerr)
        call make_oris%add_input(UI_PARM, 'even', 'binary', 'Generate even projections', 'Generate quasi-even projection directions(yes|no){no}', '(yes|no){no}', .false., 'no')
        call make_oris%add_input(UI_PARM, 'ndiscrete', 'num', 'Number of discrete projection directions', 'Number of discrete projection directions to sample from', '# discrete projs', .false., 0.)
        call make_oris%add_input(UI_PARM, oritype)
        ! alternative inputs
        ! <empty>
        ! search controls
        ! <empty>
        ! filter controls
        ! <empty>
        ! mask controls
        ! <empty>
        ! computer controls
        call make_oris%add_input(UI_COMP, nthr)
        ! add to ui_hash
        call add_ui_program('make_oris', make_oris, prgtab)
    end subroutine new_make_oris

    subroutine new_orisops( prgtab )
        class(ui_hash), intent(inout) :: prgtab
        ! PROGRAM SPECIFICATION
        call orisops%new(&
        &'orisops',&                      ! name
        &'Standard orientation editing',& ! descr_short
        &'is a program for modifying SIMPLE orientation/parameter files. If errify=yes,&
        & uniform random angular errors, and uniform origin shift errors, &
        & and uniform random defocus errors are introduced. If nstates > 1 then random states are assigned.&
        & If mirr=2d, then the Euler angles in oritab are mirrored according to the relation&
        & e1=e1, e2=180+e2, e3=-e3. If mirr=3d, then the Euler angles in oritab are mirrored according to the&
        & relation R=M(M*R), where R is the rotation matrix calculated from the Euler angle triplet and M is a&
        & 3D reflection matrix. If e1, e2, or e3 is&
        & inputted, the orientations are rotated correspondingly. If you input state,&
        & only the orientations assigned to state state are rotated. If mul is defined, the origin shifts are multiplied with mul.&
        & If zero=yes, then the shifts are zeroed',&
        &'simple_exec',&                  ! executable
        &.false.)                         ! requires sp_project
        ! INPUT PARAMETER SPECIFICATIONS
        ! image input/output
        ! <empty>
        ! parameter input/output
        call orisops%add_input(UI_PARM, oritab, required_override=.true.)
        call orisops%add_input(UI_PARM, outfile)
        call orisops%add_input(UI_PARM, e1)
        call orisops%add_input(UI_PARM, e2)
        call orisops%add_input(UI_PARM, e3)
        call orisops%add_input(UI_PARM, 'nstates', 'num', 'Number of random state labels', 'Number of random state labels to insert', '# states', .false., 0.0)
        call orisops%add_input(UI_PARM, pgrp, required_override=.false.)
        call orisops%add_input(UI_PARM, ctf,  required_override=.false.)
        call orisops%add_input(UI_PARM, angerr)
        call orisops%add_input(UI_PARM, sherr)
        call orisops%add_input(UI_PARM, dferr)
        call orisops%add_input(UI_PARM, 'zero', 'binary', 'Zero shifts', 'Zero shifts(yes|no){no}', '(yes|no){no}', .false., 'no')
        call orisops%add_input(UI_PARM, 'ndiscrete', 'num', 'Number of discrete projection directions',&
        &'Number of projection directions to use for discretization of input orientations', '# discrete projs', .false., 0.)
        call orisops%add_input(UI_PARM, 'state', 'num', 'State to modify', 'Index of state to modify', 'give state index', .false., 1.)
        call orisops%add_input(UI_PARM, 'mul', 'num', 'Shift multiplication factor',&
        &'Origin shift multiplication factor{1}','1/scale in pixels{1}', .false., 1.)
        call orisops%add_input(UI_PARM, 'mirr', 'multi', 'Mirror orientations', 'Mirror orientations(2d|3d|no){no}', '(2d|3d|no){no}', .false., 'no')
        call orisops%add_input(UI_PARM, 'symrnd', 'binary', 'Randomize over subgroubs of point-group', 'Expand orientations over entire unit sphere by &
        &permutation according to randomly selected subgroup symmetry(yes|no){no}', '(yes|no){no}', .false., 'no')
        call orisops%add_input(UI_PARM, oritype)
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
        call add_ui_program('orisops', orisops, prgtab)
    end subroutine new_orisops
    
    subroutine new_oristats( prgtab )
        class(ui_hash), intent(inout) :: prgtab
        ! PROGRAM SPECIFICATION
        call oristats%new(&
        &'oristats',&                             ! name
        &'Statistical analyses of orientations',& ! descr_short
        &'is a program for analyzing SIMPLE orientation/parameter files. If two orientation&
        & tables (oritab and oritab2) are inputted, statistics of the distances between the orientations&
        & in the two documents are provided',&
        &'simple_exec',&                          ! executable
        &.false.)                                 ! requires sp_project
        ! INPUT PARAMETER SPECIFICATIONS
        ! image input/output
        ! <empty>
        ! parameter input/output
        call oristats%add_input(UI_PARM, oritab, required_override=.true.)
        call oristats%add_input(UI_PARM, oritab2)
        call oristats%add_input(UI_PARM, pgrp, required_override=.false.)
        call oristats%add_input(UI_PARM, nspace)
        call oristats%add_input(UI_PARM, oritype)
        call oristats%add_input(UI_PARM, 'ctfstats',  'binary', 'CTF statistics',        'Provide statistics about CTF(yes|no){no}',                      '(yes|no){no}', .false., 'no')
        call oristats%add_input(UI_PARM, 'classtats', 'binary', 'Class statistics',      'Provide statistics about 2D clusters(yes|no){no}',              '(yes|no){no}', .false., 'no')
        call oristats%add_input(UI_PARM, 'projstats', 'binary', 'Projection statistics', 'Provide statistics about projection directions(yes|no){no}',    '(yes|no){no}', .false., 'no')
        call oristats%add_input(UI_PARM, 'trsstats',  'binary', 'Shift statistics',      'Provide statistics about rotational origin shifts(yes|no){no}', '(yes|no){no}', .false., 'no')
        ! alternative inputs
        ! <empty>
        ! search controls
        ! <empty>
        ! filter controls
        ! <empty>
        ! mask controls
        ! <empty>
        ! computer controls
        call oristats%add_input(UI_COMP, nthr)
        ! add to ui_hash
        call add_ui_program('oristats', oristats, prgtab)
    end subroutine new_oristats

    subroutine new_vizoris( prgtab )
        class(ui_hash), intent(inout) :: prgtab
        ! PROGRAM SPECIFICATION
        call vizoris%new(&
        &'vizoris',&                                                                                               ! name
        &'Visualization of orientation distribution',&                                                             ! descr_short
        &'is a program for extracting projection directions from orientations for visualization in UCSF Chimera',& ! descr_long
        &'all',&                                                                                                   ! executable
        &.false.)                                                                                                  ! requires sp_project
        ! INPUT PARAMETER SPECIFICATIONS
        ! image input/output
        ! <empty>
        ! parameter input/output
        call vizoris%add_input(UI_PARM, oritab, required_override=.true.)
        call vizoris%add_input(UI_PARM, nspace)
        call vizoris%add_input(UI_PARM, pgrp)
        call vizoris%add_input(UI_PARM, oritype)
        call vizoris%add_input(UI_PARM, 'tseries', 'binary', 'Time series analysis', 'Orientations originate from analysis of a time-series(yes|no){no}', '(yes|no){no}', .false., 'no')
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
        call add_ui_program('vizoris', vizoris, prgtab)
    end subroutine new_vizoris

end module simple_ui_api_ori
