!@descr: module defining the user interfaces for trajectory analysis of nanoparticles in the single_exec suite
module single_ui_trajectory
use simple_ui_modules
implicit none

type(category_descriptor), parameter :: UI_CATEGORY = category_descriptor('trajectory', 'Trajectory Analysis', 20)
type(ui_program), target :: trajectory_extract_substk
type(ui_program), target :: trajectory_graphene_subtr
type(ui_program), target :: trajectory_import
type(ui_program), target :: trajectory_denoise
type(ui_program), target :: trajectory_make_projavgs
type(ui_program), target :: trajectory_reconstruct3D
type(ui_program), target :: trajectory_swap_stack
type(ui_program), target :: trajectory_solve3D

contains

    subroutine construct_single_trajectory_programs(prgtab)
        class(ui_hash), intent(inout) :: prgtab
        call new_trajectory_extract_substk(prgtab)
        call new_trajectory_graphene_subtr(prgtab)
        call new_trajectory_import(prgtab)
        call new_trajectory_denoise(prgtab)
        call new_trajectory_make_projavgs(prgtab)
        call new_trajectory_reconstruct3D(prgtab)
        call new_trajectory_swap_stack(prgtab)
        call new_trajectory_solve3D(prgtab)
    end subroutine construct_single_trajectory_programs

subroutine new_trajectory_extract_substk( prgtab )
        class(ui_hash), intent(inout) :: prgtab
        ! PROGRAM SPECIFICATION
        call trajectory_extract_substk%new(&
        &'trajectory_extract_substk',&                                                                                            ! name
        &'extraction of a substack segment of time-series of metallic nanoparticles',&                                 ! summary
        &'is a shared-memory workflow for extraction of a substack segment of time-series of metallic nanoparticles',& ! help
        &'single_exec',&                                                                                               ! executable
        &.true., visibility=UI_VIS_ADVANCED, &
        &display_name='Extract a Time-Series Substack') ! requires sp_project
        ! INPUT PARAMETER SPECIFICATIONS
        ! image input/output
        ! <empty>
        ! parameter input/output
        call trajectory_extract_substk%add_input(UI_FILE, projfile, &
        &visibility=UI_VIS_STANDARD)
        call trajectory_extract_substk%add_input(UI_PARM, 'fromp', 'num', 'From index', 'Start index for stack copy', 'start index', .false., 1.0, &
        &visibility=UI_VIS_ADVANCED)
        call trajectory_extract_substk%add_input(UI_PARM, 'top',   'num', 'To index', 'Stop index for stack copy', 'stop index', .false., 1.0, &
        &visibility=UI_VIS_ADVANCED)
        call trajectory_extract_substk%add_input(UI_PARM, 'state', 'num', 'State index', 'Only particles with this state are extracted{1}', 'state index', .false., 1.0, &
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
        call add_ui_program('trajectory_extract_substk', trajectory_extract_substk, prgtab, UI_CATEGORY)
    end subroutine new_trajectory_extract_substk

    subroutine new_trajectory_graphene_subtr( prgtab )
        class(ui_hash), intent(inout) :: prgtab
        ! PROGRAM SPECIFICATION
        call trajectory_graphene_subtr%new(&
        &'trajectory_graphene_subtr',&                                  ! name
        &'Suppress graphene Fourier peaks in nanoparticle time-series images',& ! summary
        &'Removes graphene Fourier peaks in time-series',&   ! help
        &'single_exec',&                                     ! executable
        &.false., visibility=UI_VIS_STANDARD, display_name='Remove Graphene Signal') ! requires sp_project
        ! INPUT PARAMETER SPECIFICATIONS
        ! image input/output
        call trajectory_graphene_subtr%add_input(UI_IMG, stk, label_override='Tracked NP image stack', &
        &help_override='xxx.mrc file with NP time-trajectory', placeholder_override='e.g. trajectory.mrcs', &
        &required_override=.true., visibility=UI_VIS_STANDARD)
        call trajectory_graphene_subtr%add_input(UI_IMG, stk2, label_override='Background power spectra stack', &
        &help_override='background power spectra stack, one per image of the trajectory', &
        &placeholder_override='e.g. background_pspec.mrcs', required_override=.true., visibility=UI_VIS_STANDARD)
        call trajectory_graphene_subtr%add_input(UI_IMG, outstk, &
        &visibility=UI_VIS_ADVANCED)
        ! parameter input/output
        call trajectory_graphene_subtr%add_input(UI_PARM, smpd, &
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
        call trajectory_graphene_subtr%add_input(UI_COMP, nthr, &
        &visibility=UI_VIS_STANDARD)
        ! add to ui_hash
        call add_ui_program('trajectory_graphene_subtr', trajectory_graphene_subtr, prgtab, UI_CATEGORY)
    end subroutine new_trajectory_graphene_subtr

    subroutine new_trajectory_import( prgtab )
        class(ui_hash), intent(inout) :: prgtab
        ! PROGRAM SPECIFICATION
        call trajectory_import%new(&
        &'trajectory_import',&                    ! name
        &'Import particle-image time-series data for trajectory analysis',& ! summary
        &'is a workflow for importing time-series data',&   ! help
        &'single_exec',&                                    ! executable
        &.true., visibility=UI_VIS_STANDARD, display_name='Import Particle Time Series') ! requires sp_project
        ! INPUT PARAMETER SPECIFICATIONS
        ! image input/output
        call trajectory_import%add_input(UI_IMG, stk, required_override=.true., &
        &visibility=UI_VIS_STANDARD)
        ! parameter input/output
        call trajectory_import%add_input(UI_PARM, smpd, &
        &visibility=UI_VIS_STANDARD)
        call trajectory_import%add_input(UI_FILE, deftab, &
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
        call add_ui_program('trajectory_import', trajectory_import, prgtab, UI_CATEGORY)
    end subroutine new_trajectory_import

    subroutine new_trajectory_denoise( prgtab )
        class(ui_hash), intent(inout) :: prgtab
        ! PROGRAM SPECIFICATION
        call trajectory_denoise%new(&
        &'trajectory_denoise',&                                       ! name
        &'Denoise time-series particle images',&                      ! summary
        &'is a program for denoising a particle trajectory: every frame becomes the kernel pre-image of the frames '//&
        &'that look like it, without alignment',&                     ! help
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
        call trajectory_denoise%add_input(UI_FILT, 'neigs', 'num', 'Largest latent rank (0 => 32; default 32)', &
        &'The largest latent rank the cross-validation may choose (0 => 32; default 32)', '# eigenvecs', .false., 0.0, &
        &visibility=UI_VIS_ADVANCED)
        call trajectory_denoise%add_input(UI_FILT, 'halfwin', 'num', 'Window half-width in frames (0 => the whole trajectory)', &
        &'Frames pool only within this many frames of their block; 0 => every frame may pool with every frame (default 0)', &
        &'# frames', .false., 0.0, visibility=UI_VIS_ADVANCED)
        call trajectory_denoise%add_input(UI_FILT, 'eo_pool', 'binary', 'Pool by row parity', &
        &'A frame pools only frames of its own row parity, so the even/odd halves stay independent(yes|no){yes}', '', .false., 'yes', &
        &choices=ui_choices([character(len=3) :: 'yes', 'no']), visibility=UI_VIS_ADVANCED)
        call trajectory_denoise%add_input(UI_FILT, lp, &
        &visibility=UI_VIS_ADVANCED)
        ! mask controls
        call trajectory_denoise%add_input(UI_MASK, mskdiam, required_override=.false., &
        &help_override='Mask diameter in A (shapes the frame similarity only, the output is unmasked; '//&
        &'default the inscribed circle)', visibility=UI_VIS_ADVANCED)
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
        &.true., visibility=UI_VIS_DEVELOPER, &
        &display_name='Average Initial Trajectory Frames')                                                   ! requires sp_project
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
        &.true., visibility=UI_VIS_DEVELOPER, &
        &display_name='Time-Windowed 3D Reconstruction')                                   ! requires sp_project
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
        &.true., visibility=UI_VIS_DEVELOPER, &
        &display_name='Swap a Trajectory Stack')                                    ! requires sp_project
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

    subroutine new_trajectory_solve3D( prgtab )
        class(ui_hash), intent(inout) :: prgtab
        ! PROGRAM SPECIFICATION
        call trajectory_solve3D%new(&
        &'trajectory_solve3D',&                                                                               ! name
        &'De novo time-integrated map of a tracked nanoparticle from its graphene-subtracted frames',&      ! summary
        &'is a shared-memory workflow from the graphene-subtracted frames of one tracked particle to its time-integrated '//&
        &'map. On the raw frames: import, center2D_nano on ten time-chunk averages, the diameter on the centred averages, '//&
        &'the start model simulated from it, refine2D_nano from the centering''s shifts, then two coupled gates on its '//&
        &'class averages (the lattice gate, and the model gate: the class averages aligned to the start model, the ones '//&
        &'far below the median correlation rejected) and the prune. Only the kept frames are then denoised (flex) and '//&
        &'refined in 3D by autorefine3D_nano. One numbered subdirectory per stage; '//&
        &'trajectory_solve3D.txt records each stage. Rerun with mkdir=no in the execution directory to resume after '//&
        &'the last finished stage',&                                                                          ! help
        &'single_exec',&                                                                                      ! executable
        &.true., visibility=UI_VIS_STANDARD, display_name='Solve a Nanoparticle Trajectory') ! requires sp_project
        ! INPUT PARAMETER SPECIFICATIONS
        ! image input/output
        call trajectory_solve3D%add_input(UI_IMG, 'stk', 'file', 'Graphene-subtracted frames', &
        &'Graphene-subtracted frames of the tracked particle', 'e.g. subtr_frames.mrc', .true., '', &
        &visibility=UI_VIS_STANDARD)
        ! parameter input/output
        call trajectory_solve3D%add_input(UI_PARM, smpd, &
        &visibility=UI_VIS_STANDARD)
        call trajectory_solve3D%add_input(UI_PARM, element, required_override=.false., &
        &help_override='Atom kind of the simulated start model{Pt}', visibility=UI_VIS_ADVANCED)
        call trajectory_solve3D%add_input(UI_PARM, 'moldiam', 'num', 'Particle diameter of the start model', &
        &'Particle diameter (A) of the simulated start model; default the smallest diameter estimated on the centred '//&
        &'time-chunk averages', 'diameter in A', .false., 0., visibility=UI_VIS_ADVANCED)
        call trajectory_solve3D%add_input(UI_PARM, 'thres', 'num', 'Lattice-peak threshold of the lattice gate', &
        &'Class averages whose lattice_peak (log10 of the strongest Fourier peak over its ring mean) is below this are '//&
        &'rejected by the lattice gate; the default is the grade set for averages of raw frames{0.52}', &
        &'log10 peak over ring{0.52}', .false., 0.52, visibility=UI_VIS_ADVANCED)
        call trajectory_solve3D%add_input(UI_PARM, 'model_gate_nsd', 'num', 'Model-gate cut in robust spreads', &
        &'The model gate rejects a class average whose correlation to its best reprojection of the start model lies '//&
        &'more than this many robust spreads (1.4826 x MAD) below the median of the lattice gate''s survivors{3}', &
        &'robust spreads{3}', .false., 3., visibility=UI_VIS_ADVANCED)
        ! search controls
        call trajectory_solve3D%add_input(UI_SRCH, nspace, &
        &visibility=UI_VIS_ADVANCED)
        call trajectory_solve3D%add_input(UI_SRCH, 'maxits', 'num', 'Max 3D refinement rounds', &
        &'Rounds of refinement and atom detection of autorefine3D_nano', 'Max # rounds{5}', .false., 5., &
        &visibility=UI_VIS_ADVANCED)
        call trajectory_solve3D%add_input(UI_SRCH, 'maxits_between', 'num', 'Max iterations per 3D round', &
        &'Refinement iterations between two atom detections of autorefine3D_nano', 'Max # iterations{10}', .false., 10., &
        &visibility=UI_VIS_ADVANCED)
        ! filter controls
        call trajectory_solve3D%add_input(UI_FILT, 'lp', 'num', 'Initial 3D low-pass limit', &
        &'Initial low-pass limit of autorefine3D_nano', 'low-pass limit in Angstroms{1.5}', .false., 1.5, &
        &visibility=UI_VIS_STANDARD)
        call trajectory_solve3D%add_input(UI_FILT, 'halfwin', 'num', 'Denoising window half-width in frames', &
        &'Frames pool only within this many frames of their block; 0 => the whole trajectory (default 0)', &
        &'# frames', .false., 0.0, visibility=UI_VIS_ADVANCED)
        call trajectory_solve3D%add_input(UI_FILT, 'neigs', 'num', 'Denoising components', &
        &'The largest latent rank the denoiser''s cross-validation may choose (default 32)', '# eigenvecs', .false., 0.0, &
        &visibility=UI_VIS_ADVANCED)
        call trajectory_solve3D%add_input(UI_FILT, 'eo_pool', 'binary', 'Denoise by row parity', &
        &'A frame pools only frames of its own row parity(yes|no){yes}', '', .false., 'yes', &
        &choices=ui_choices([character(len=3) :: 'yes', 'no']), visibility=UI_VIS_ADVANCED)
        ! mask controls
        call trajectory_solve3D%add_input(UI_MASK, mskdiam, required_override=.false., &
        &help_override='Mask diameter in A of every stage after the centering; default 1.5 x the largest diameter '//&
        &'estimated on the centred time-chunk averages', visibility=UI_VIS_ADVANCED)
        ! computer controls
        call trajectory_solve3D%add_input(UI_COMP, nthr, &
        &visibility=UI_VIS_STANDARD)
        ! add to ui_hash
        call add_ui_program('trajectory_solve3D', trajectory_solve3D, prgtab, UI_CATEGORY)
    end subroutine new_trajectory_solve3D

end module single_ui_trajectory
