!@descr: module defining the user interfaces for 3D reconstruction programs in the simple_exec suite
module simple_ui_reconstruct3D
use simple_ui_modules
implicit none

type(category_descriptor), parameter :: UI_CATEGORY = &
    category_descriptor('reconstruct3d', 'Reconstruct 3D Workflows', 68)
type(ui_program), target :: reconstruct3D
type(ui_program), target :: bootstrap_rec3D

contains

    subroutine construct_reconstruct3D_programs(prgtab)
        class(ui_hash), intent(inout) :: prgtab
        call new_reconstruct3D(prgtab)
        call new_bootstrap_rec3D(prgtab)
    end subroutine construct_reconstruct3D_programs

    subroutine new_reconstruct3D( prgtab )
        class(ui_hash), intent(inout) :: prgtab
        ! PROGRAM SPECIFICATION
        call reconstruct3D%new(&
        &'reconstruct3D',&                                               ! name
        &'3D reconstruction from oriented particles',&                   ! summary
        &'is a distributed workflow for reconstructing volumes from MRC and SPIDER stacks,&
        & given input orientations and state assignments. The algorithm is based on direct Fourier inversion&
        & with a Kaiser-Bessel (KB) interpolation kernel',&
        &'simple_exec',&                                                 ! executable
        &.true., &
        &visibility=UI_VIS_ADVANCED)                                                         ! requires sp_project
        ! INPUT PARAMETER SPECIFICATIONS
        ! image input/output
        ! <empty>
        ! parameter input/output
        call reconstruct3D%add_input(UI_PARM, 'rec_backend', 'multi', 'Reconstruction backend', &
        &'Reconstruction backend; PCG runs independent kernel solves for the two halfsets(gridding|pcg){gridding}', &
        &'', .false., 'gridding', &
        &choices=ui_choices([character(len=8) :: 'gridding', 'pcg']), &
        &visibility=UI_VIS_STANDARD)
        call reconstruct3D%add_input(UI_PARM, 'box_crop', 'num', 'Reconstruction box', &
        &'Even Fourier-cropped reconstruction box; native project geometry remains authoritative', &
        &'pixels{native box}', .false., 0.0, visibility=UI_VIS_ADVANCED)
        call reconstruct3D%add_input(UI_PARM, 'euclid_diag', 'binary', 'Euclid scale diagnostics', &
        &'Per-iteration report of the reference/particle amplitude ratio per band and the euclid objective quantiles(yes|no){no}','', .false., 'no', visibility=UI_VIS_ADVANCED, &
        &choices=ui_choices([character(len=3) :: 'yes', 'no']))
        call reconstruct3D%add_input(UI_PARM, 'projrec', 'binary', 'Projection-direction reconstruction',&
        &'Assemble raw 2D Fourier numerator/CTF-squared sums by projection direction before compact 3D reconstruction(yes|no){no}','', .false., 'no', &
        &choices=ui_choices([character(len=3) :: 'yes', 'no']), &
        &visibility=UI_VIS_ADVANCED)
        ! <no additional inputs>
        ! <empty>
        ! search controls
        call reconstruct3D%add_input(UI_SRCH, trs, &
        &visibility=UI_VIS_ADVANCED)
        call reconstruct3D%add_input(UI_SRCH, pgrp, &
        &visibility=UI_VIS_STANDARD)
        call reconstruct3D%add_input(UI_SRCH, ptcl_src, &
        &visibility=UI_VIS_ADVANCED)
        ! filter controls
        call reconstruct3D%add_input(UI_FILT, envfsc, &
        &visibility=UI_VIS_ADVANCED)
        call reconstruct3D%add_input(UI_FILT, envmsklp, &
        &visibility=UI_VIS_ADVANCED)
        call reconstruct3D%add_input(UI_FILT, 'postprocess', 'binary', 'Postprocess final map',&
        &'Postprocess reconstructed volumes using the generated FSC curves','', .false., 'yes', &
        &choices=ui_choices([character(len=3) :: 'yes', 'no']), &
        &visibility=UI_VIS_ADVANCED)
        call reconstruct3D%add_input(UI_FILT, 'maxits_pcg', 'num', 'PCG maximum iterations', &
        &'Maximum kernel PCG iterations; used only when rec_backend=pcg', 'iterations{2}', .false., 2., &
        &visibility=UI_VIS_ADVANCED, &
        &activation=ui_activation_equals_any('rec_backend', [character(len=3) :: 'pcg']))
        call reconstruct3D%add_input(UI_FILT, 'maxits_ml', 'num', 'Regularized-solve PCG iterations', &
        &'Coupled PCG iterations of the ML-regularized system from the closed-form Wiener start; 0 = closed form only', 'iterations{0}', &
        &.false., 0., visibility=UI_VIS_ADVANCED, &
        &activation=ui_activation_equals_any('rec_backend', [character(len=3) :: 'pcg']))
        call reconstruct3D%add_input(UI_FILT, 'rtol', 'num', 'PCG relative residual tolerance', &
        &'Stop at this true L2 relative residual; use <=0 for exactly maxits_pcg iterations', 'tolerance{0}', &
        &.false., 0.0, visibility=UI_VIS_ADVANCED, &
        &activation=ui_activation_equals_any('rec_backend', [character(len=3) :: 'pcg']))
        call reconstruct3D%add_input(UI_FILT, 'pcg_mskfile', 'file', 'PCG support-constraint mask volume', &
        &'Real-space [0,1] mask volume installed as the hard support constraint of the PCG solve (the projected '//&
        &'system P H P; experimental focused/support mode); spherical mskdiam support when absent', &
        &'e.g. focusmask.mrc', .false., '', group="filter", visibility=UI_VIS_ADVANCED, &
        &activation=ui_activation_equals_any('rec_backend', [character(len=3) :: 'pcg']))
        call reconstruct3D%add_input(UI_MASK, mskdiam, &
        &visibility=UI_VIS_STANDARD)
        ! computer controls
        call reconstruct3D%add_input(UI_COMP, nparts, required_override=.false., &
        &visibility=UI_VIS_ADVANCED)
        call reconstruct3D%add_input(UI_COMP, nthr, &
        &visibility=UI_VIS_STANDARD)
        ! add to ui_hash
        call add_ui_program('reconstruct3D', reconstruct3D, prgtab, UI_CATEGORY)
    end subroutine new_reconstruct3D

    subroutine new_bootstrap_rec3D( prgtab )
        class(ui_hash), intent(inout) :: prgtab
        ! PROGRAM SPECIFICATION
        call bootstrap_rec3D%new(&
        &'bootstrap_rec3D',&                                             ! name
        &'bootstrap sigma2 and ML-regularized 3D reconstruction',&       ! summary
        &'complete final-reconstruction sequence for a project with 3D orientations: seeds the sigma2 estimate from particle&
        & power spectra (calc_pspec) as the given iteration, assembles a gridding ML-regularized bootstrap map on it (with the&
        & given filt_mode/automsk, since the residual sigmas depend on the reference regularization), runs one&
        & residual sigma2 pass (refine=sigma, no search) against that map, consolidates the residual groups as the next&
        & iteration and reconstructs the shipped ML-regularized map on them with the requested backend (PCG gets the cold-solve&
        & iteration budget); standalone test entry point for the final reconstruction stage of abinitio3D and refine3D_auto',&
        &'simple_exec',&                                                 ! executable
        &.true.)                                                         ! requires sp_project
        ! INPUT PARAMETER SPECIFICATIONS
        ! image input/output
        ! <empty>
        ! parameter input/output
        call bootstrap_rec3D%add_input(UI_PARM, 'rec_backend', 'multi', 'Reconstruction backend', &
        &'Reconstruction backend; PCG runs independent kernel solves for the two halfsets(gridding|pcg){gridding}', &
        &'', .false., 'gridding', &
        &choices=ui_choices([character(len=8) :: 'gridding', 'pcg']), &
        &visibility=UI_VIS_STANDARD)
        call bootstrap_rec3D%add_input(UI_PARM, 'which_iter', 'num', 'Sigma iteration index',&
        &'Iteration number given to the residual sigma pass and its iteration files{1}', 'iteration{1}', .false., 1.0, &
        &visibility=UI_VIS_DEVELOPER)
        call bootstrap_rec3D%add_input(UI_FILE, 'outfile', 'file', 'Resolution output prefix',&
        &'Optional FSC/resolution text output prefix; state tags are appended', 'e.g. resolution',&
        &.false., 'RESOLUTION_FINAL.txt', &
        &visibility=UI_VIS_DEVELOPER)
        ! <no additional inputs>
        ! <empty>
        ! search controls
        call bootstrap_rec3D%add_input(UI_SRCH, pgrp, &
        &visibility=UI_VIS_STANDARD)
        call bootstrap_rec3D%add_input(UI_SRCH, nstates, &
        &visibility=UI_VIS_DEVELOPER)
        ! filter controls
        call bootstrap_rec3D%add_input(UI_FILT, 'postprocess', 'binary', 'Postprocess final map',&
        &'Postprocess ML-regularized reconstructed volumes using the generated FSC curves','', .false., 'yes', &
        &choices=ui_choices([character(len=3) :: 'yes', 'no']), &
        &visibility=UI_VIS_DEVELOPER)
        call bootstrap_rec3D%add_input(UI_FILT, combine_eo, &
        &visibility=UI_VIS_DEVELOPER)
        call bootstrap_rec3D%add_input(UI_FILT, 'maxits_pcg', 'num', 'PCG maximum iterations', &
        &'Maximum kernel PCG iterations; used only when rec_backend=pcg', 'iterations{2}', .false., 2., &
        &visibility=UI_VIS_ADVANCED, &
        &activation=ui_activation_equals_any('rec_backend', [character(len=3) :: 'pcg']))
        call bootstrap_rec3D%add_input(UI_FILT, 'maxits_ml', 'num', 'Regularized-solve PCG iterations', &
        &'Coupled PCG iterations of the ML-regularized system from the closed-form Wiener start; 0 = closed form only', 'iterations{0}', &
        &.false., 0., visibility=UI_VIS_ADVANCED, &
        &activation=ui_activation_equals_any('rec_backend', [character(len=3) :: 'pcg']))
        call bootstrap_rec3D%add_input(UI_FILT, 'rtol', 'num', 'PCG relative residual tolerance', &
        &'Stop at this true L2 relative residual; use <=0 for exactly maxits_pcg iterations', 'tolerance{0}', &
        &.false., 0.0, visibility=UI_VIS_ADVANCED, &
        &activation=ui_activation_equals_any('rec_backend', [character(len=3) :: 'pcg']))
        ! mask controls
        call bootstrap_rec3D%add_input(UI_MASK, mskdiam, &
        &visibility=UI_VIS_STANDARD)
        ! computer controls
        call bootstrap_rec3D%add_input(UI_COMP, nparts, required_override=.false., &
        &visibility=UI_VIS_DEVELOPER)
        call bootstrap_rec3D%add_input(UI_COMP, nthr, &
        &visibility=UI_VIS_STANDARD)
        ! add to ui_hash
        call add_ui_program('bootstrap_rec3D', bootstrap_rec3D, prgtab, UI_CATEGORY)
    end subroutine new_bootstrap_rec3D

end module simple_ui_reconstruct3D
