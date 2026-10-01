!@descr: module defining the user interface for the Cartesian pose refinement workflow
module simple_ui_refine3D_pose_cont
use simple_ui_modules
implicit none
private

! The parent refine3D UI module supplies the category so category identity and
! ordering remain single-source while this independently evolving workflow owns
! its complete public parameter surface.
public :: construct_refine3D_pose_cont_program

type(ui_program), target :: refine3D_pose_cont

contains

    subroutine construct_refine3D_pose_cont_program(prgtab, category)
        class(ui_hash), intent(inout) :: prgtab
        type(category_descriptor), intent(in) :: category

        ! This is a top-level workflow contract, not a direct matcher contract.
        ! The commander translates pose_cont_mode into the internal refine,
        ! pose_cont, inpl_cont, and pose_cont_route settings for each child
        ! stage; exposing those controls here would permit contradictory stage
        ! schedules.
        call refine3D_pose_cont%new(&
        &'refine3D_pose_cont',&
        &'Automatically refine and continuously polish one 3D structure',&
        &'is an automated single-state 3D workflow with mandatory Cartesian pose refinement',&
        &'simple_exec',&
        &.true.,&
        &visibility=UI_VIS_DEVELOPER, display_name='Automated 3D Refinement with Cartesian Pose Polishing')

        ! Workflow input and reconstruction ownership. An optional vol1 seeds
        ! the inherited refine3D_auto startup; rec_backend applies consistently
        ! to iterative half-map assembly and the native-sampling final solve.
        call refine3D_pose_cont%add_input(UI_IMG, 'vol1', 'file', 'Starting template volume', &
        &'Starting reference volume for particle matching', 'input starting volume e.g. vol.mrc', .false., '', &
        &visibility=UI_VIS_ADVANCED)
        call refine3D_pose_cont%add_input(UI_PARM, 'rec_backend', 'multi', 'Reconstruction backend', &
        &'Reconstruction backend for per-iteration half-map assembly and final reconstruction'//&
        &'(gridding|pcg){gridding}', '', .false., 'gridding', group='search', &
        &choices=ui_choices([character(len=8) :: 'gridding', 'pcg']), visibility=UI_VIS_ADVANCED)

        call refine3D_pose_cont%add_input(UI_SRCH, 'objfun', 'multi', 'Objective function', &
        &'Outer probabilistic-matcher objective(euclid|cc){euclid}', '', .false., 'euclid', group='search', &
        &choices=ui_choices([character(len=6) :: 'euclid', 'cc']), visibility=UI_VIS_STANDARD)

        ! Cartesian placement policy:
        !   post_matcher     polishes each probabilistic winner in-place;
        !   standalone_final checkpoints the completed probabilistic state and
        !                    adds one full Cartesian pass after the main loop.
        ! The global parameter default remains off for other programs. This
        ! mandatory-Cartesian workflow defaults to the validated local-polish
        ! placement after each probabilistic matcher winner.
        call refine3D_pose_cont%add_input(UI_SRCH, 'pose_cont_mode', 'multi', 'Cartesian workflow mode', &
        &'Placement policy for mandatory Cartesian polishing'//&
        &'(post_matcher|standalone_final){post_matcher}', '', .false., 'post_matcher', group='search', &
        &choices=ui_choices([character(len=16) :: 'post_matcher', 'standalone_final']), &
        &visibility=UI_VIS_STANDARD)
        call refine3D_pose_cont%add_input(UI_SRCH, maxits, required_override=.false., group='search', &
        &visibility=UI_VIS_ADVANCED)
        call refine3D_pose_cont%add_input(UI_SRCH, 'minits', 'num', 'Minimum automatic iterations', &
        &'Minimum main iterations before convergence may stop the workflow', &
        &'minimum iterations', .false., 3., group='search', visibility=UI_VIS_DEVELOPER)
        call refine3D_pose_cont%add_input(UI_SRCH, 'nsample', 'num', 'Projection samples', &
        &'Number of projection samples used by automatic refinement', &
        &'number of samples', .false., 25000., group='search', visibility=UI_VIS_DEVELOPER)

        ! Startup reference provenance and pose scaffold. ref_pose_init=cc is
        ! reserved for an independently supplied vol1; the optional global
        ! registration pass establishes discrete basins and then applies the
        ! first mandatory joint Cartesian polish. Its FSC threshold limits
        ! that pass rather than the later local refinement objective.
        call refine3D_pose_cont%add_input(UI_SRCH, 'ref_pose_init', 'multi', 'External-reference pose initialization', &
        &'When vol1 has independent provenance, run one fixed-reference CC pose-initialization pass at 15 Angstroms '//&
        &'before Euclidean refinement(cc|none){none}', '', .false., 'none', group='search', &
        &choices=ui_choices([character(len=4) :: 'cc', 'none']), visibility=UI_VIS_ADVANCED)
        call refine3D_pose_cont%add_input(UI_SRCH, pgrp, group='search', visibility=UI_VIS_STANDARD)
        call refine3D_pose_cont%add_input(UI_SRCH, sigma_est, group='search', visibility=UI_VIS_ADVANCED)
        call refine3D_pose_cont%add_input(UI_SRCH, 'center', 'binary', 'Center reference volume(s)', &
        &'Center reference volume(s) by their center of gravity and map shifts back to the particles(yes|no){no}', &
        &'', .false., 'no', group='search', choices=ui_choices([character(len=3) :: 'yes', 'no']), &
        &visibility=UI_VIS_ADVANCED)
        call refine3D_pose_cont%add_input(UI_SRCH, 'regpass', 'binary', 'Global registration pass', &
        &'One global refine=greedy pass of all particles against the masked startup references before refinement'//&
        &'(yes|no){yes}', '', .false., 'yes', group='search', &
        &choices=ui_choices([character(len=3) :: 'yes', 'no']), visibility=UI_VIS_ADVANCED)
        call refine3D_pose_cont%add_input(UI_SRCH, 'regpass_fsc', 'num', 'Registration-pass FSC criterion', &
        &'FSC value of the startup pair whose resolution band-limits the registration pass', &
        &'FSC value in (0,1){0.143}', .false., 0.143, group='search', visibility=UI_VIS_ADVANCED)

        ! Working-grid and lifecycle controls inherited from refine3D_auto.
        ! Autoscaling changes the computational grid while project metadata
        ! remains authoritative for native sampling and final reconstruction.
        call refine3D_pose_cont%add_input(UI_SRCH, 'autoscale', 'binary', 'Automatic down-scaling', &
        &'Automatic down-scaling of images for accelerated computation(yes|no){yes}', '', .false., 'yes', &
        &group='search', choices=ui_choices([character(len=3) :: 'yes', 'no']), visibility=UI_VIS_ADVANCED)
        call refine3D_pose_cont%add_input(UI_SRCH, 'continue', 'binary', 'Continue previous refinement', &
        &'Continue previous refinement(yes|no){no}', '', .false., 'no', group='search', &
        &choices=ui_choices([character(len=3) :: 'yes', 'no']), visibility=UI_VIS_ADVANCED)

        ! Volume conditioning between pose stages. These controls shape the
        ! references consumed by subsequent iterations and the final products;
        ! they do not configure the Cartesian optimizer. PCG-only controls are
        ! activated from rec_backend so the gridding contract stays uncluttered.
        call refine3D_pose_cont%add_input(UI_FILT, 'amsklp', 'num', 'NU evidence envelope smoothing limit', &
        &'Low-pass limit for NU evidence envelope generation in Angstroms', 'low-pass limit in Angstroms', &
        &.false., 8., group='filter', visibility=UI_VIS_ADVANCED)
        call refine3D_pose_cont%add_input(UI_FILT, 'filt_mode', 'multi', 'Filtering mode', &
        &'Filtering mode(none|nonuniform|nonuniform_lpset){nonuniform}', '', .false., 'nonuniform', &
        &group='filter', choices=ui_choices([character(len=16) :: 'none', 'nonuniform', 'nonuniform_lpset']), &
        &visibility=UI_VIS_ADVANCED)
        call refine3D_pose_cont%add_input(UI_FILT, 'maxits_pcg', 'num', 'PCG maximum iterations', &
        &'Maximum kernel PCG iterations during refinement; the final cold solve uses at least 5', &
        &'iterations{2}', .false., 2., group='filter', visibility=UI_VIS_ADVANCED, &
        &activation=ui_activation_equals_any('rec_backend', [character(len=3) :: 'pcg']))
        call refine3D_pose_cont%add_input(UI_FILT, 'maxits_ml', 'num', 'Regularized-solve PCG iterations', &
        &'Coupled PCG iterations of the ML-regularized system; 0 = closed form only', 'iterations{0}', &
        &.false., 0., group='filter', visibility=UI_VIS_ADVANCED, &
        &activation=ui_activation_equals_any('rec_backend', [character(len=3) :: 'pcg']))
        call refine3D_pose_cont%add_input(UI_FILT, 'pcg_solvent', 'binary', 'PCG soft solvent prior', &
        &'Enable the evidence-derived soft solvent prior in PCG reconstruction(yes|no){no}', '', .false., 'no', &
        &group='filter', visibility=UI_VIS_ADVANCED, choices=ui_choices([character(len=3) :: 'yes', 'no']), &
        &activation=ui_activation_equals_any('rec_backend', [character(len=3) :: 'pcg']))
        call refine3D_pose_cont%add_input(UI_FILT, 'pcg_solvent_lambda', 'num', 'PCG solvent prior strength', &
        &'Ridge coefficient of the solvent prior relative to the data scale', 'coefficient{1.0}', &
        &.false., 1.0, group='filter', visibility=UI_VIS_ADVANCED, &
        &activation=ui_activation_equals_any('pcg_solvent', [character(len=3) :: 'yes']))
        call refine3D_pose_cont%add_input(UI_FILT, envfsc, group='filter', visibility=UI_VIS_ADVANCED)
        call refine3D_pose_cont%add_input(UI_FILT, envmsklp, group='filter', visibility=UI_VIS_ADVANCED)
        call refine3D_pose_cont%add_input(UI_FILT, combine_eo, group='filter', visibility=UI_VIS_ADVANCED)
        call refine3D_pose_cont%add_input(UI_FILT, 'res_target', 'num', 'Resolution target (in A)', &
        &'Resolution target in Angstroms', 'Resolution target in Angstroms', .false., 3., &
        &group='filter', visibility=UI_VIS_ADVANCED)

        ! Mask controls define the common reference support seen by startup,
        ! probabilistic matching, Cartesian polishing, and reconstruction.
        call refine3D_pose_cont%add_input(UI_MASK, mskdiam, group='mask', visibility=UI_VIS_STANDARD)
        call refine3D_pose_cont%add_input(UI_MASK, automsk_refine3D, group='mask', visibility=UI_VIS_ADVANCED)
        call refine3D_pose_cont%add_input(UI_MASK, nu_msk_sig, group='mask', visibility=UI_VIS_ADVANCED)

        ! Execution topology changes partitioning and threading only; scientific
        ! stage order remains the commander's responsibility.
        call refine3D_pose_cont%add_input(UI_COMP, nparts, group='compute', visibility=UI_VIS_STANDARD)
        call refine3D_pose_cont%add_input(UI_COMP, nthr, group='compute', visibility=UI_VIS_STANDARD)

        call add_ui_program('refine3D_pose_cont', refine3D_pose_cont, prgtab, category)
    end subroutine construct_refine3D_pose_cont_program

end module simple_ui_refine3D_pose_cont
