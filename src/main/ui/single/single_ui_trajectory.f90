!@descr: module defining the user interfaces for trajectory analysis of nanoparticles in the single_exec suite
module single_ui_trajectory
use simple_ui_modules
implicit none

type(category_descriptor), parameter :: UI_CATEGORY = category_descriptor('trajectory', 'Trajectory Analysis', 20)
type(ui_program), target :: extract_substk
type(ui_program), target :: graphene_subtr
type(ui_program), target :: import_trajectory
type(ui_program), target :: trajectory_denoise
type(ui_program), target :: trajectory_make_projavgs
type(ui_program), target :: trajectory_reconstruct3D
type(ui_program), target :: trajectory_swap_stack

contains

    subroutine construct_single_trajectory_programs(prgtab)
        class(ui_hash), intent(inout) :: prgtab
        call new_extract_substk(prgtab)
        call new_graphene_subtr(prgtab)
        call new_import_trajectory(prgtab)
        call new_trajectory_denoise(prgtab)
        call new_trajectory_make_projavgs(prgtab)
        call new_trajectory_reconstruct3D(prgtab)
        call new_trajectory_swap_stack(prgtab)
    end subroutine construct_single_trajectory_programs

subroutine new_extract_substk( prgtab )
        class(ui_hash), intent(inout) :: prgtab
        ! PROGRAM SPECIFICATION
        call extract_substk%new(&
        &'extract_substk',&                                                                                            ! name
        &'extraction of a substack segment of time-series of metallic nanoparticles',&                                 ! summary
        &'is a shared-memory workflow for extraction of a substack segment of time-series of metallic nanoparticles',& ! help
        &'single_exec',&                                                                                               ! executable
        &.true., visibility=UI_VIS_ADVANCED)                                                                                 ! requires sp_project
        ! INPUT PARAMETER SPECIFICATIONS
        ! image input/output
        ! <empty>
        ! parameter input/output
        call extract_substk%add_input(UI_FILE, projfile, &
        &visibility=UI_VIS_STANDARD)
        call extract_substk%add_input(UI_PARM, 'fromp', 'num', 'From index', 'Start index for stack copy', 'start index', .false., 1.0, &
        &visibility=UI_VIS_ADVANCED)
        call extract_substk%add_input(UI_PARM, 'top',   'num', 'To index', 'Stop index for stack copy', 'stop index', .false., 1.0, &
        &visibility=UI_VIS_ADVANCED)
        call extract_substk%add_input(UI_PARM, 'state', 'num', 'State index', 'Only particles with this state are extracted{1}; use state<0 for legacy include-all behavior', 'state index', .false., 1.0, &
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
        call add_ui_program('extract_substk', extract_substk, prgtab, UI_CATEGORY)
    end subroutine new_extract_substk

    subroutine new_graphene_subtr( prgtab )
        class(ui_hash), intent(inout) :: prgtab
        ! PROGRAM SPECIFICATION
        call graphene_subtr%new(&
        &'graphene_subtr',&                                  ! name
        &'Suppress graphene Fourier peaks in nanoparticle time-series images',& ! summary
        &'Removes graphene Fourier peaks in time-series',&   ! help
        &'single_exec',&                                     ! executable
        &.false., visibility=UI_VIS_STANDARD, display_name='Remove Graphene Signal') ! requires sp_project
        ! INPUT PARAMETER SPECIFICATIONS
        ! image input/output
        call graphene_subtr%add_input(UI_IMG, stk_traj, &
        &visibility=UI_VIS_STANDARD)
        call graphene_subtr%add_input(UI_IMG, stk_backgr, &
        &visibility=UI_VIS_STANDARD)
        call graphene_subtr%add_input(UI_IMG, outstk, &
        &visibility=UI_VIS_ADVANCED)
        ! parameter input/output
        call graphene_subtr%add_input(UI_PARM, smpd, &
        &visibility=UI_VIS_STANDARD)
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
        call graphene_subtr%add_input(UI_COMP, nthr, &
        &visibility=UI_VIS_STANDARD)
        ! add to ui_hash
        call add_ui_program('graphene_subtr', graphene_subtr, prgtab, UI_CATEGORY)
    end subroutine new_graphene_subtr

    subroutine new_import_trajectory( prgtab )
        class(ui_hash), intent(inout) :: prgtab
        ! PROGRAM SPECIFICATION
        call import_trajectory%new(&
        &'import_trajectory',&                    ! name
        &'Import particle-image time-series data for trajectory analysis',& ! summary
        &'is a workflow for importing time-series data',&   ! help
        &'single_exec',&                                    ! executable
        &.true., visibility=UI_VIS_STANDARD, display_name='Import Particle Time Series') ! requires sp_project
        ! INPUT PARAMETER SPECIFICATIONS
        ! image input/output
        call import_trajectory%add_input(UI_IMG, stk, required_override=.true., &
        &visibility=UI_VIS_STANDARD)
        ! parameter input/output
        call import_trajectory%add_input(UI_PARM, smpd, &
        &visibility=UI_VIS_STANDARD)
        call import_trajectory%add_input(UI_FILE, deftab, &
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
        call add_ui_program('import_trajectory', import_trajectory, prgtab, UI_CATEGORY)
    end subroutine new_import_trajectory

    subroutine new_trajectory_denoise( prgtab )
        class(ui_hash), intent(inout) :: prgtab
        ! PROGRAM SPECIFICATION
        call trajectory_denoise%new(&
        &'trajectory_denoise',&                                       ! name
        &'Denoise time-series particle images with diffusion-map analysis',& ! summary
        &'is a program for diffusion-map denoising of an image stack',& ! help
        &'single_exec',&                                              ! executable
        &.false., visibility=UI_VIS_STANDARD, display_name='Denoise Particle Trajectories') ! requires sp_project
        ! INPUT PARAMETER SPECIFICATIONS
        ! image input/output
        call trajectory_denoise%add_input(UI_IMG, 'stk',  'file', 'Stack to denoise',  'Stack of images to denoise', 'e.g. stk.mrcs', .true., '', &
        &visibility=UI_VIS_STANDARD)
        call trajectory_denoise%add_input(UI_IMG, outstk, &
        &visibility=UI_VIS_ADVANCED)
        ! parameter input/output
        call trajectory_denoise%add_input(UI_PARM, smpd, &
        &visibility=UI_VIS_STANDARD)
        ! <no additional inputs>
        ! <empty>
        ! search controls
        ! <empty>
        ! filter controls
        call trajectory_denoise%add_input(UI_FILT, 'neigs', 'num', 'Number of diffusion-map components (0 => auto; default 0)', 'Number of diffusion-map components (0 => auto; default 0)', '# eigenvecs', .false., 0.0, &
        &visibility=UI_VIS_ADVANCED)
        call trajectory_denoise%add_input(UI_FILT, 'pca_mode', 'multi', 'PCA methods: diffusion maps, PPCA, PPCA plus residual kPCA, standard SVD PCA, or kernel PCA', 'PCA methods','', .false., 'diffusion_maps', &
        &choices=ui_choices([character(len=15) :: 'diffusion_maps', 'ppca', 'ppca_kpca_resid', 'pca_svd', 'kpca']), &
        &visibility=UI_VIS_ADVANCED)
        call trajectory_denoise%add_input(UI_FILT, 'k_nn', 'num', 'Diffusion graph neighbors (default 5; try 5-30)', 'Local nearest neighbors used for pca_mode=diffusion_maps', '# neighbors', .false., 5.0, &
        &visibility=UI_VIS_ADVANCED)
        call trajectory_denoise%add_input(UI_FILT, 'kpca_ker', 'multi', 'Kernel PCA kernel', 'Kernel PCA kernel(rbf|cosine){rbf}','', .false., 'rbf', &
        &choices=ui_choices([character(len=6) :: 'rbf', 'cosine']), &
        &visibility=UI_VIS_ADVANCED)
        call trajectory_denoise%add_input(UI_FILT, 'kpca_backend', 'multi', 'Kernel PCA backend', 'Kernel PCA backend(exact|nystrom){nystrom}','', .false., 'nystrom', &
        &choices=ui_choices([character(len=7) :: 'exact', 'nystrom']), &
        &visibility=UI_VIS_ADVANCED)
        call trajectory_denoise%add_input(UI_FILT, 'kpca_rbf_gamma', 'num', 'RBF gamma (0 => auto)', 'RBF gamma (0 => auto)', 'gamma', .false., 0.0, &
        &visibility=UI_VIS_ADVANCED)
        call trajectory_denoise%add_input(UI_FILT, 'ppca_kpca_resid_alpha', 'num', 'Residual hybrid damping (0 => PPCA only; default 0.5)', 'Residual hybrid damping (0 => PPCA only; default 0.5)', 'hybrid alpha', .false., 0.5, &
        &visibility=UI_VIS_ADVANCED)
        call trajectory_denoise%add_input(UI_FILT, 'kpca_nystrom_npts', 'num', 'Nyström landmark count (0 => auto=max(128,2*neigs), capped at 512; try 256, 512)', 'Nyström landmark count (0 => auto=max(128,2*neigs), capped at 512; try 256, 512)', '# landmarks', .false., 512.0, &
        &visibility=UI_VIS_ADVANCED)
        call trajectory_denoise%add_input(UI_FILT, 'kpca_nystrom_local_nbrs', 'num', 'Nyström max local support neighbors (default 96; try 96, 128)', 'Nyström max local support neighbors (default 96; try 96, 128)', '# max local nbrs', .false., 96.0, &
        &visibility=UI_VIS_ADVANCED)
        call trajectory_denoise%add_input(UI_FILT, hp, &
        &visibility=UI_VIS_ADVANCED)
        call trajectory_denoise%add_input(UI_FILT, lp, &
        &visibility=UI_VIS_ADVANCED)
        ! mask controls
        call trajectory_denoise%add_input(UI_MASK, mskdiam, required_override=.false., &
        &visibility=UI_VIS_ADVANCED)
        ! computer controls
        call trajectory_denoise%add_input(UI_COMP, nthr, &
        &visibility=UI_VIS_STANDARD)
        ! add to ui_hash
        call add_ui_program('trajectory_denoise', trajectory_denoise, prgtab, UI_CATEGORY)
    end subroutine new_trajectory_denoise

    subroutine new_trajectory_make_projavgs( prgtab )
        class(ui_hash), intent(inout) :: prgtab
        ! PROGRAM SPECIFICATION
        call trajectory_make_projavgs%new(&
        &'trajectory_make_projavgs',&                                                    ! name
        &'Align & average the first few frames of the time-series',&                     ! summary
        &'is a program for aligning & averaging the first few frames of the time-series&
        & to accomplish SNR enhancement for particle identification',&                   ! help
        &'single_exec',&                                                                 ! executable
        &.true., visibility=UI_VIS_DEVELOPER)                                                   ! requires sp_project
        ! INPUT PARAMETER SPECIFICATIONS
        ! image input/output
        ! <empty>
        ! parameter input/output
        call trajectory_make_projavgs%add_input(UI_PARM, nspace, &
        &visibility=UI_VIS_DEVELOPER)
        call trajectory_make_projavgs%add_input(UI_PARM, 'athres', 'num', 'Angular threshold (degrees)', 'Angular threshold (degrees)', 'in degrees{10}', .false., 10., &
        &visibility=UI_VIS_DEVELOPER)
        ! <no additional inputs>
        ! <empty>
        ! search controls
        ! <empty>
        ! filter controls
        ! <empty>
        ! mask controls
        call trajectory_make_projavgs%add_input(UI_MASK, mskdiam, &
        &visibility=UI_VIS_STANDARD)
        ! computer controls
        call trajectory_make_projavgs%add_input(UI_COMP, nthr, &
        &visibility=UI_VIS_STANDARD)
        ! add to ui_hash
        call add_ui_program('trajectory_make_projavgs', trajectory_make_projavgs, prgtab, UI_CATEGORY)
    end subroutine new_trajectory_make_projavgs

    subroutine new_trajectory_reconstruct3D( prgtab )
        class(ui_hash), intent(inout) :: prgtab
        ! PROGRAM SPECIFICATION
        call trajectory_reconstruct3D%new(&
        &'trajectory_reconstruct3D',&                                    ! name
        &'Time windowed 3D reconstruction from oriented particles',&     ! help
        &'Time windowed 3D reconstruction from oriented particles',&
        &'single_exec',&                                                 ! executable
        &.true., visibility=UI_VIS_DEVELOPER)                                   ! requires sp_project
        ! INPUT PARAMETER SPECIFICATIONS
        ! image input/output
        ! parameter input/output
        call trajectory_reconstruct3D%add_input(UI_PARM, 'stepsz',  'num', 'Time window size (# frames){500}', 'Time window size (# frames) for windowed 3D rec{500}', 'give # frames',  .false., 500., &
        &visibility=UI_VIS_DEVELOPER)
        call trajectory_reconstruct3D%add_input(UI_PARM, 'fromp', 'num', 'From particle index', 'Start index for 3D reconstruction', 'start index', .false., 1.0, &
        &visibility=UI_VIS_DEVELOPER)
        call trajectory_reconstruct3D%add_input(UI_PARM, 'top',   'num', 'To particle index', 'Stop index for 3D reconstruction', 'stop index', .false., 1.0, &
        &visibility=UI_VIS_DEVELOPER)
        call trajectory_reconstruct3D%add_input(UI_PARM, 'nchunks', 'num', 'Number of temporal chunks', &
        &'Fixed number of contiguous chunks; overrides the stepsz-derived count when positive', '# chunks', .false., 0., &
        &visibility=UI_VIS_DEVELOPER)
        ! <no additional inputs>
        ! <empty>
        ! search controls
        call trajectory_reconstruct3D%add_input(UI_SRCH, pgrp, &
        &visibility=UI_VIS_STANDARD)
        call trajectory_reconstruct3D%add_input(UI_SRCH, maxits, required_override=.false., &
        &visibility=UI_VIS_DEVELOPER)
        call trajectory_reconstruct3D%add_input(UI_SRCH, trs, &
        &visibility=UI_VIS_DEVELOPER)
        ! filter controls
        call trajectory_reconstruct3D%add_input(UI_FILT, lp, required_override=.false., &
        &visibility=UI_VIS_DEVELOPER)
        ! mask controls
        call trajectory_reconstruct3D%add_input(UI_MASK, mskdiam, &
        &visibility=UI_VIS_STANDARD)
        ! computer controls
        call trajectory_reconstruct3D%add_input(UI_COMP, nparts, required_override=.false., &
        &visibility=UI_VIS_DEVELOPER)
        call trajectory_reconstruct3D%add_input(UI_COMP, nthr, &
        &visibility=UI_VIS_STANDARD)
        ! add to ui_hash
        call add_ui_program('trajectory_reconstruct3D', trajectory_reconstruct3D, prgtab, UI_CATEGORY)
    end subroutine new_trajectory_reconstruct3D

    subroutine new_trajectory_swap_stack( prgtab )
        class(ui_hash), intent(inout) :: prgtab
        ! PROGRAM SPECIFICATION
        call trajectory_swap_stack%new(&
        &'trajectory_swap_stack',&                                        ! name
        &'Substitutes stack into an existing project',&                   ! summary
        &'is a program for substituting stack into an existing project',& ! help
        &'single_exec',&                                                  ! executable
        &.true., visibility=UI_VIS_DEVELOPER)                                    ! requires sp_project
        ! INPUT PARAMETER SPECIFICATIONS
        ! image input/output
        call trajectory_swap_stack%add_input(UI_IMG, stk, required_override=.true., &
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
        call add_ui_program('trajectory_swap_stack', trajectory_swap_stack, prgtab, UI_CATEGORY)
    end subroutine new_trajectory_swap_stack

end module single_ui_trajectory
