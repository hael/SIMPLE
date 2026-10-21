!@descr: module defining the user interfaces for 2D clustering and related programs in the simple_exec suite
module simple_ui_refine2D
use simple_ui_modules
implicit none

type(category_descriptor), parameter :: UI_CATEGORY = category_descriptor('refine2d', 'Refine 2D Workflows', 30)
type(ui_program), target :: solve2D
type(ui_program), target :: refine2D
type(ui_program), target :: solve2D_chunks
type(ui_program), target :: make_cavgs
type(ui_program), target :: bootstrap_cavgs
type(ui_program), target :: unbootstrap_cavgs
type(ui_program), target :: map_cavgs_selection
type(ui_program), target :: sample_classes
type(ui_program), target :: write_classes
type(ui_program), target :: cls_expansion

contains

    subroutine construct_refine2D_programs( prgtab )
        class(ui_hash), intent(inout) :: prgtab
        call new_solve2D(prgtab)
        call new_refine2D(prgtab)
        call new_solve2D_chunks(prgtab)
        call new_make_cavgs(prgtab)
        call new_bootstrap_cavgs(prgtab)
        call new_unbootstrap_cavgs(prgtab)
        call new_map_cavgs_selection(prgtab)
        call new_sample_classes(prgtab)
        call new_write_classes(prgtab)
        call new_cls_expansion(prgtab)
    end subroutine construct_refine2D_programs

    subroutine new_solve2D( prgtab )
        class(ui_hash), intent(inout) :: prgtab
        call new_solve2D_descriptor(solve2D, 'solve2D', &
            &'De novo 2D class-average determination from particle images', &
            &'is a distributed workflow, the 2D equivalent of solve3D, that couples ab initio 2D '//&
            &'classification with initial 2D refinement of the class averages', &
            &UI_VIS_STANDARD, 'De Novo 2D Class Averages')
        call add_ui_program('solve2D', solve2D, prgtab, UI_CATEGORY)
    end subroutine new_solve2D

    !> refine2D, the 2D counterpart of refine3D: shared memory, or the distributed master with nparts (its private
    !! worker of the same name, simple_private_exec prg=refine2D, is checked against these required keys)
    subroutine new_refine2D( prgtab )
        class(ui_hash), intent(inout) :: prgtab
        ! PROGRAM SPECIFICATION
        call refine2D%new(&
        &'refine2D',&                                                                       ! name
        &'Refine 2D class averages from particle images',&                                 ! summary
        &'is a distributed workflow for 2D alignment and clustering of particle images into class averages, '//&
        &'the 2D counterpart of refine3D and the program solve2D runs its stages with',&    ! help
        &'simple_exec',&                                                                    ! executable
        &.true.,&                                                                           ! requires sp_project
        &visibility=UI_VIS_ADVANCED, display_name='Refine 2D Class Averages')
        ! INPUT PARAMETER SPECIFICATIONS
        ! image input/output
        call refine2D%add_input(UI_FILE, projfile, required_override=.true., visibility=UI_VIS_STANDARD)
        ! parameter input/output
        call refine2D%add_input(UI_PARM, 'cache', 'binary', 'Cache downscaled particles', &
         &'Write Fourier-cropped particles once and read those instead of the originals on every iteration, &
         &trading disk space for I/O(yes|no){no}', '', &
         &.false., 'no', choices=ui_choices([character(len=3) :: 'yes', 'no']), &
         &visibility=UI_VIS_DEVELOPER)
        ! search controls
        call refine2D%add_input(UI_SRCH, ncls, required_override=.false., group="search", visibility=UI_VIS_ADVANCED)
        call refine2D%add_input(UI_SRCH, 'refine', 'multi', 'Refinement mode', &
        &'Refinement mode: snhc_smpl and snhc search stochastically, prob assigns from the dense probabilistic table, &
        &prob_snhc from a sparse one(snhc_smpl|snhc|prob|prob_snhc){snhc_smpl}', '', &
        &.false., 'snhc_smpl', group="search", &
        &choices=ui_choices([character(len=9) :: 'snhc_smpl', 'snhc', 'prob', 'prob_snhc']), &
        &visibility=UI_VIS_ADVANCED)
        call refine2D%add_input(UI_SRCH, maxits, group="search", visibility=UI_VIS_ADVANCED)
        call refine2D%add_input(UI_SRCH, startit, group="search", visibility=UI_VIS_ADVANCED)
        call refine2D%add_input(UI_SRCH, update_frac, group="search", visibility=UI_VIS_ADVANCED)
        call refine2D%add_input(UI_SRCH, trs, group="search", visibility=UI_VIS_ADVANCED)
        call refine2D%add_input(UI_SRCH, objfun, group="search", visibility=UI_VIS_ADVANCED)
        call refine2D%add_input(UI_SRCH, cls_init, group="search", visibility=UI_VIS_ADVANCED)
        call refine2D%add_input(UI_SRCH, 'center', 'binary', 'Center class averages', 'Center class averages by their &
        &center of gravity and map shifts back to the particles(yes|no){yes}','', .false., 'yes', group="search", &
        &choices=ui_choices([character(len=3) :: 'yes', 'no']), &
        &visibility=UI_VIS_ADVANCED)
        call refine2D%add_input(UI_SRCH, 'autoscale', 'binary', 'Automatic down-scaling', 'Automatic down-scaling of images &
        &for accelerated computation(yes|no){yes}','', .false., 'yes', group="search", &
        &choices=ui_choices([character(len=3) :: 'yes', 'no']), &
        &visibility=UI_VIS_ADVANCED)
        ! filter controls
        call refine2D%add_input(UI_FILT, hp, group="filter", visibility=UI_VIS_ADVANCED)
        call refine2D%add_input(UI_FILT, lp, group="filter", visibility=UI_VIS_ADVANCED)
        call refine2D%add_input(UI_FILT, ml_reg, group="filter", visibility=UI_VIS_ADVANCED)
        ! mask controls
        call refine2D%add_input(UI_MASK, mskdiam, group="mask", visibility=UI_VIS_STANDARD)
        ! computer controls
        call refine2D%add_input(UI_COMP, nparts, required_override=.false., group="compute", visibility=UI_VIS_STANDARD)
        call refine2D%add_input(UI_COMP, nthr, required_override=.false., group="compute", visibility=UI_VIS_STANDARD)
        call add_ui_program('refine2D', refine2D, prgtab, UI_CATEGORY)
    end subroutine new_refine2D

    subroutine new_solve2D_descriptor( program, name, summary, help, visibility, display_name )
        type(ui_program), intent(inout) :: program
        character(len=*), intent(in) :: name, summary, help, display_name
        integer, intent(in) :: visibility
        ! PROGRAM SPECIFICATION
        call program%new(&
        &name,&                                                                         ! name
        &summary,&                                                                      ! summary
        &help,&                                                                         ! help
        &'simple_exec',&                                                               ! executable
        &.true.,&                                                                      ! requires sp_project
        &visibility=visibility, display_name=display_name)
        ! INPUT PARAMETER SPECIFICATIONS
        ! image input/output
        ! <empty>
        ! parameter input/output
        call program%add_input(UI_PARM, 'cache', 'binary', 'Cache downscaled particles', &
         &'Write Fourier-cropped particles once and read those instead of the originals on every iteration, &
         &trading disk space for I/O(yes|no){no}', '', &
         &.false., 'no', choices=ui_choices([character(len=3) :: 'yes', 'no']), &
         &visibility=UI_VIS_DEVELOPER)
        call program%add_input(UI_PARM, 'cache_dir', 'dir', 'Particle cache directory', &
         &'Where to keep the downscaled particle cache; point it at a fast local disk when the project lives &
         &on a slow one. Defaults to the execution directory', 'e.g. /scratch/ptcl_cache/', &
         &.false., '', visibility=UI_VIS_DEVELOPER)
        ! <no additional inputs>
        ! <empty>
        ! search controls
        call program%add_input(UI_SRCH, ncls, group="search", visibility=UI_VIS_STANDARD)
        call program%add_input(UI_SRCH, 'center', 'binary', 'Center class averages', 'Center class averages by their &
        &center of gravity and map shifts back to the particles(yes|no){yes}','', .false., 'yes', group="model", &
        &choices=ui_choices([character(len=3) :: 'yes', 'no']), &
        &visibility=UI_VIS_ADVANCED)
        call program%add_input(UI_SRCH, 'autoscale', 'binary', 'Automatic down-scaling', 'Automatic down-scaling of images &
        &for accelerated computation(yes|no){yes}','', .false., 'yes', group="model", &
        &choices=ui_choices([character(len=3) :: 'yes', 'no']), &
        &visibility=UI_VIS_ADVANCED)
        call program%add_input(UI_SRCH, 'refine', 'multi', 'Refinement mode',&
        &'Refinement mode of the stages from the third on (the first two always use snhc_smpl): prob_snhc scores &
        &a sparse probabilistic table and assigns the last stage from a dense one; prob uses the dense table in every &
        &stage; snhc_smpl is prob_snhc with the last stage sparse too(prob_snhc|prob|snhc_smpl){prob_snhc}','', &
        &.false., 'prob_snhc', group="search", &
        &choices=ui_choices([character(len=9) :: 'prob_snhc', 'prob', 'snhc_smpl']), &
        &visibility=UI_VIS_ADVANCED)
        call program%add_input(UI_SRCH, cls_init, group="search", &
        &visibility=UI_VIS_ADVANCED)
        call program%add_input(UI_SRCH, 'nsample', 'num', 'Number of particles to sample', 'Number of particles &
        &to sample each iteration; every active particle when there are fewer{200000}', '# particles to sample{200000}', &
        &.false., 200000., group="search", visibility=UI_VIS_STANDARD, preserve_default=.true.)
        call program%add_input(UI_SRCH, 'nstages', 'num', 'Number of stages', 'Number of solve2D stages, at most 6 &
        &(5 with eo_stage=no): the low-pass limits from lpstart to lpstop are spread over the stages before the last, &
        &and the last is the gold-standard even/odd stage{6}', '# stages{6}', .false., 6., group="search", &
        &visibility=UI_VIS_ADVANCED, preserve_default=.true.)
        call program%add_input(UI_SRCH, 'stats', 'binary', 'Write class-average statistics', &
        &'Write class-average statistics during solve2D(yes|no){no}', '', .false., 'no', group="search", &
        &choices=ui_choices([character(len=3) :: 'yes', 'no']), &
        &visibility=UI_VIS_DEVELOPER)
        call program%add_input(UI_SRCH, 'extr_lim', 'num', 'Extreme-particle iteration limit', &
        &'Iteration limit used to identify and exclude extreme particles during solve2D', &
        &'iteration limit', .false., 20., group="search", visibility=UI_VIS_DEVELOPER)
        call program%add_input(UI_SRCH, 'nits_per_stage', 'num', 'Iterations per stage', &
        &'Number of solve2D iterations added at each frequency-marching stage', &
        &'iterations per stage', .false., 5., group="search", visibility=UI_VIS_DEVELOPER)
        call program%add_input(UI_SRCH, 'eo_stage', 'binary', 'Use even/odd stage', &
         &'Run the final solve2D frequency-marching stage with even/odd class averages(yes|no){yes}', '', &
         &.false., 'yes', group="search", choices=ui_choices([character(len=3) :: 'yes', 'no']), &
         &visibility=UI_VIS_DEVELOPER)
        ! filter controls
        call program%add_input(UI_FILT, hp, group="filter", &
        &visibility=UI_VIS_ADVANCED)
        call program%add_input(UI_FILT, ml_reg, group="filter", &
        &visibility=UI_VIS_ADVANCED)
        call program%add_input(UI_FILT, 'cenlp', 'num', 'Centering low-pass limit', 'Limit for low-pass filter used in binarisation &
        &prior to determination of the center of gravity of the class averages and centering; by default the mask &
        &diameter / 6, kept within 20-30 A and at least twice the cropped pixel size', 'centering low-pass limit in &
        &Angstroms', .false., 30., group="filter", &
        &visibility=UI_VIS_ADVANCED)
        call program%add_input(UI_FILT, 'lpstart', 'num', 'Initial low-pass limit', 'Low-pass limit of the first solve2D &
            &stage; by default the mask diameter / 12, kept within 15-20 A and at least twice the cropped pixel size',&
            &'low-pass limit in Angstroms', .false., 30., group="filter", &
        &visibility=UI_VIS_ADVANCED)
        call program%add_input(UI_FILT, 'lpstop',  'num', 'Final low-pass limit', 'Low-pass limit the stages march to, &
            &reached before the gold-standard stage; by default the mask diameter / 22, kept within 6-8 A and at least &
            &twice the cropped pixel size', 'low-pass limit in Angstroms', .false., 6., group="filter", &
        &visibility=UI_VIS_ADVANCED)
        call program%add_input(UI_FILT, lp, group="filter", &
        &visibility=UI_VIS_ADVANCED)
        ! mask controls
        call program%add_input(UI_MASK, mskdiam, group="mask", visibility=UI_VIS_STANDARD)
        ! computer controls
        call program%add_input(UI_COMP, nparts, required_override=.false., group="compute", visibility=UI_VIS_STANDARD)
        call program%add_input(UI_COMP, nthr, group="compute", visibility=UI_VIS_STANDARD)
    end subroutine new_solve2D_descriptor

    subroutine new_solve2D_chunks( prgtab )
        class(ui_hash), intent(inout) :: prgtab
        ! PROGRAM SPECIFICATION
        call solve2D_chunks%new(&
        &'solve2D_chunks',&                                                                       ! name
        &'Independent solve2D analysis of particle subsets',&                                     ! summary
        &'splits a project into particle-balanced subsets and runs independent solve2D jobs',& ! help
        &'simple_exec',&                                                                          ! executable
        &.true.,&                                                                                 ! requires sp_project
        &visibility=UI_VIS_ADVANCED, display_name='Solve 2D in Particle Chunks')
        ! INPUT PARAMETER SPECIFICATIONS
        ! image input/output
        ! <empty>
        ! parameter input/output
        call solve2D_chunks%add_input(UI_PARM, 'nchunks', 'num', 'Number of chunks', &
            &'Number of particle-balanced subset projects to run with independent solve2D jobs. &
            &Omit or set to 0 to target about 100 classes per chunk.', &
            &'# of chunks (0=auto)', .false., 0., &
        &visibility=UI_VIS_ADVANCED)
        call solve2D_chunks%add_input(UI_PARM, 'nchunks_in_parallel', 'num', 'Number of chunks run in parallel', &
            &'Maximum number of chunks executed at the same time. Each chunk uses nparts parts, so the &
            &peak load is about nchunks_in_parallel*nparts processes. Defaults to 1 (one chunk at a time). &
            &Set to 0 to run all chunks in parallel.', &
            &'# of parallel chunks (0=all){1}', .false., 1., &
        &visibility=UI_VIS_ADVANCED)
        ! <no additional inputs>
        ! <empty>
        ! search controls
        call solve2D_chunks%add_input(UI_SRCH, nptcls_per_cls, placeholder_override='# of particles per cluster{500}', group="cluster 2D", visibility=UI_VIS_STANDARD)
        call solve2D_chunks%add_input(UI_SRCH, 'center', 'binary', 'Center class averages', 'Center class averages by their center of &
            &gravity and map shifts back to the particles(yes|no){yes}','', .false., 'yes', group="cluster 2D", &
        &choices=ui_choices([character(len=3) :: 'yes', 'no']), &
        &visibility=UI_VIS_ADVANCED)
        call solve2D_chunks%add_input(UI_SRCH, 'refine', 'multi', 'Refinement mode',&
        &'Refinement mode(prob_snhc|prob|snhc_smpl){prob_snhc}','', .false., 'prob_snhc', group="cluster 2D", &
        &choices=ui_choices([character(len=9) :: 'prob_snhc', 'prob', 'snhc_smpl']), &
        &visibility=UI_VIS_ADVANCED)
        ! filter controls
        call solve2D_chunks%add_input(UI_FILT, hp, group="cluster 2D", &
        &visibility=UI_VIS_ADVANCED)
        call solve2D_chunks%add_input(UI_FILT, lp, group="cluster 2D", &
        &visibility=UI_VIS_ADVANCED)
        call solve2D_chunks%add_input(UI_FILT, 'cenlp', 'num', 'Centering low-pass limit', 'Limit for low-pass filter used in binarisation &
        &prior to determination of the center of gravity of the class averages and centering', 'centering low-pass limit in &
        &Angstroms{30}', .false., 30., group="cluster 2D", &
        &visibility=UI_VIS_ADVANCED)
        call solve2D_chunks%add_input(UI_FILT, 'lpstop', 'num', 'Final low-pass limit', 'Low-pass limit that controls the degree of &
        &downsampling in the second phase. Give estimated best final resolution', 'final low-pass limit in Angstroms', .false., 8.,&
        &group="filter", visibility=UI_VIS_ADVANCED)
        ! mask controls
        call solve2D_chunks%add_input(UI_MASK, mskdiam, group="cluster 2D", visibility=UI_VIS_STANDARD)
        ! computer controls
        call solve2D_chunks%add_input(UI_COMP, nparts, required_override=.false., group="compute", visibility=UI_VIS_STANDARD)
        call solve2D_chunks%add_input(UI_COMP, nthr, group="compute", visibility=UI_VIS_STANDARD)
        call solve2D_chunks%add_input(UI_COMP, 'walltime', 'num', 'Walltime', 'Maximum execution time for job scheduling and &
        &management(29mins){1740}', 'in seconds(29mins){1740}', .false., 1740., group="compute", &
        &visibility=UI_VIS_ADVANCED)
        ! add to ui_hash
        call add_ui_program('solve2D_chunks', solve2D_chunks, prgtab, UI_CATEGORY)
    end subroutine new_solve2D_chunks

    subroutine new_make_cavgs( prgtab )
        class(ui_hash), intent(inout) :: prgtab
        ! PROGRAM SPECIFICATION
        call make_cavgs%new(&
        &'make_cavgs', &                           ! name
        &'Create class averages from aligned particle images',& ! summary
        &'is a distributed workflow for generating class averages or initial random references&
        & for refine2D execution',&               ! help
        &'simple_exec',&                           ! executable
        &.true., &
        &visibility=UI_VIS_ADVANCED, display_name='Create Class Averages')                                   ! requires sp_project
        ! INPUT PARAMETER SPECIFICATIONS
        ! image input/output
        call make_cavgs%add_input(UI_IMG, 'refs', 'file', 'Output 2D references',&
        &'Output 2D references', 'xxx.mrc file with references', .false., '', &
        &visibility=UI_VIS_ADVANCED)
        ! parameter input/output
        call make_cavgs%add_input(UI_PARM, ncls, required_override=.false., &
        &visibility=UI_VIS_ADVANCED)
        call make_cavgs%add_input(UI_PARM, 'mul', 'num', 'Shift multiplication factor',&
        &'Origin shift multiplication factor{1}','1/scale in pixels{1}', .false., 1., &
        &visibility=UI_VIS_ADVANCED)
        call make_cavgs%add_input(UI_PARM, remap_cls, &
        &visibility=UI_VIS_ADVANCED)
        call make_cavgs%add_input(UI_PARM, nspace, required_override=.false., &
        &visibility=UI_VIS_ADVANCED)
        ! <no additional inputs>
        ! <empty>
        ! search controls
        ! <empty>
        ! filter controls
        ! mask controls
        ! <empty>
        ! computer controls
        call make_cavgs%add_input(UI_COMP, nparts, &
        &visibility=UI_VIS_STANDARD)
        call make_cavgs%add_input(UI_COMP, nthr, &
        &visibility=UI_VIS_STANDARD)
        ! add to ui_hash
        call add_ui_program('make_cavgs', make_cavgs, prgtab, UI_CATEGORY)
    end subroutine new_make_cavgs

    subroutine new_bootstrap_cavgs( prgtab )
        class(ui_hash), intent(inout) :: prgtab
        ! PROGRAM SPECIFICATION
        call bootstrap_cavgs%new(&
        &'bootstrap_cavgs', &                                                   ! name
        &'Bootstrap class averages from existing 2D class memberships',& ! summary
        &'creates an oversampled class-average stack by stochastic expansion &
        &of existing 2D class memberships and then runs make_cavgs',&           ! help
        &'simple_exec',&                                                        ! executable
        &.true., &
        &display_name='Bootstrap Class Averages') ! requires sp_project
        ! INPUT PARAMETER SPECIFICATIONS
        ! image input/output
        ! <empty>
        ! parameter input/output
        call bootstrap_cavgs%add_input(UI_PARM, 'osmpl_fac', 'num', 'Oversampling factor',&
        &'Integer factor for class-average oversampling; final class count is usable parent classes times osmpl_fac{2}',&
        &'oversampling factor{2}', .false., 2., &
        &visibility=UI_VIS_DEVELOPER)
        call bootstrap_cavgs%add_input(UI_PARM, 'frac_best', 'num', 'Anchor fraction',&
        &'Fraction used with the median class population to define the objective-ranked anchor set(0-1){0.5}',&
        &'fraction{0.5}', .false., 0.5, &
        &visibility=UI_VIS_DEVELOPER)
        call bootstrap_cavgs%add_input(UI_IMG, 'refs', 'file', 'Bootstrap reference class averages', &
        &'Reference class-average stack written for bootstrap processing', 'e.g. cavgs_bootstrap_001.mrc', .false., 'cavgs_bootstrap_001.mrc', &
        &visibility=UI_VIS_DEVELOPER)
        ! <no additional inputs>
        ! <empty>
        ! search controls
        ! <empty>
        ! filter controls
        ! <empty>
        ! mask controls
        ! <empty>
        ! computer controls
        call bootstrap_cavgs%add_input(UI_COMP, nparts, required_override=.false., &
        &visibility=UI_VIS_DEVELOPER)
        call bootstrap_cavgs%add_input(UI_COMP, nthr, &
        &visibility=UI_VIS_STANDARD)
        ! add to ui_hash
        call add_ui_program('bootstrap_cavgs', bootstrap_cavgs, prgtab, UI_CATEGORY)
    end subroutine new_bootstrap_cavgs

    subroutine new_unbootstrap_cavgs( prgtab )
        class(ui_hash), intent(inout) :: prgtab
        ! PROGRAM SPECIFICATION
        call unbootstrap_cavgs%new(&
        &'unbootstrap_cavgs', &                                                ! name
        &'Map bootstrap cls3D back to original project',&                      ! summary
        &'transfers cls3D alignment from bootstrap originals (child=0) to their bootstrap_parent class indices in the original project and maps to particles',& ! help
        &'simple_exec',&                                                        ! executable
        &.true., &
        &display_name='Map Bootstrap Classes Back') ! requires sp_project
        ! INPUT PARAMETER SPECIFICATIONS
        ! image input/output
        ! <empty>
        ! parameter input/output
        ! <empty>
        ! image input/output
        call unbootstrap_cavgs%add_input(UI_FILE, 'projfile_orig', 'file', 'Original project file', &
        &'Project file that was used as input to bootstrap_cavgs and should receive the mapped cls3D/ptcl3D parameters', &
        &'e.g. original.simple', .true., '', &
        &visibility=UI_VIS_STANDARD)
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
        call add_ui_program('unbootstrap_cavgs', unbootstrap_cavgs, prgtab, UI_CATEGORY)
    end subroutine new_unbootstrap_cavgs

    subroutine new_map_cavgs_selection( prgtab )
        class(ui_hash), intent(inout) :: prgtab
        ! PROGRAM SPECIFICATION
        call map_cavgs_selection%new(&
        &'map_cavgs_selection',&                                         ! name
        &'Map class average selection to particles in project file',&    ! summary
        &'is a program for mapping selection based on class averages to the individual particles using correlation matching',& ! help
        &'simple_exec',&                                                 ! executable
        &.true., &
        &visibility=UI_VIS_ADVANCED, &
        &display_name='Map Class Selection to Particles') ! requires sp_project
        ! INPUT PARAMETER SPECIFICATIONS
        ! image input/output
        call  map_cavgs_selection%add_input(UI_IMG, 'stk', 'file', 'Stack of cavgs to select from', 'Stack of cavgs to select from', 'e.g. cavgs_iter0XX.mrc', .false., '', &
        &visibility=UI_VIS_ADVANCED)
        call  map_cavgs_selection%add_input(UI_IMG, 'stk2', 'file', 'Stack of selected cavgs', 'Stack of selected cavgs', 'e.g. selected.spi', .true., '', &
        &visibility=UI_VIS_STANDARD)
        ! parameter input/output
        call  map_cavgs_selection%add_input(UI_PARM, prune, &
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
        call add_ui_program('map_cavgs_selection', map_cavgs_selection, prgtab, UI_CATEGORY)
    end subroutine new_map_cavgs_selection

    subroutine new_sample_classes( prgtab )
        class(ui_hash), intent(inout) :: prgtab
        ! PROGRAM SPECIFICATION
        call sample_classes%new(&
        &'sample_classes',&                                                                 ! name
        &'Probabilistic sampling of particles based on class statistics',&                  ! summary
        &'is a program for probabilistic sampling of particles based on class statistics',& ! help
        &'simple_exec',&                                                                    ! executable
        &.true., &
        &visibility=UI_VIS_ADVANCED, &
        &display_name='Sample Particles by Class') ! requires sp_project
        ! INPUT PARAMETER SPECIFICATIONS
        ! image input/output
        ! <empty>
        ! parameter input/output
        call sample_classes%add_input(UI_PARM, 'nptcls_per_part', 'num',    'Number of ptcls per part to select when balancing', '# ptcls per part after balancing', '{100000}', .false., 0.0, &
        &visibility=UI_VIS_ADVANCED)
        call sample_classes%add_input(UI_PARM, 'greedy_sampling', 'binary', 'Greedy balanced selection', 'Greedy balanced selection(yes|no){yes}','', .false., 'yes', &
        &choices=ui_choices([character(len=3) :: 'yes', 'no']), &
        &visibility=UI_VIS_ADVANCED)
        call sample_classes%add_input(UI_PARM, 'nparts',          'num',    'Number of partitions in balancing', '# balanced parts', '# balanced parts', .false., 1., &
        &visibility=UI_VIS_ADVANCED)
        call sample_classes%add_input(UI_PARM, nsample, &
        &visibility=UI_VIS_ADVANCED)
        call sample_classes%add_input(UI_PARM, 'frac_best',       'num',    'Fraction of best particles to sample from', 'Fraction of best particles to sample from(0-1)', '{0.5}', .false., 0.5, &
        &visibility=UI_VIS_ADVANCED)
        call sample_classes%add_input(UI_PARM, 'frac_worst',      'num',    'Fraction of worst particles to sample from', 'Fraction of worst particles to sample from(0-1)', '{0.5}', .false., 0.5, &
        &visibility=UI_VIS_ADVANCED)
        call sample_classes%add_input(UI_PARM, prune, &
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
        call sample_classes%add_input(UI_COMP, nthr, &
        &visibility=UI_VIS_STANDARD)
        ! add to ui_hash
        call add_ui_program('sample_classes', sample_classes, prgtab, UI_CATEGORY)
    end subroutine new_sample_classes

    subroutine new_write_classes( prgtab )
        class(ui_hash), intent(inout) :: prgtab
        ! PROGRAM SPECIFICATION
        call write_classes%new(&
        &'write_classes',&                                                                                  ! name
        &'Writes the class averages and the individual (rotated and shifted) particles part of the class',& ! summary
        &'is a program for the class averages and the individual (rotated and shifted) particles part of the classto to individual stacks',& ! help
        &'simple_exec',&                                                                                    ! executable
        &.true., &
        &display_name='Write Class Members') ! requires sp_project
        ! INPUT PARAMETER SPECIFICATIONS
        ! image input/output
        ! <empty>
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
        call add_ui_program('write_classes', write_classes, prgtab, UI_CATEGORY)
    end subroutine new_write_classes

    subroutine new_cls_expansion( prgtab )
        class(ui_hash), intent(inout) :: prgtab
        call cls_expansion%new(&
        &'cls_expansion',&
        &'Split 2D classes with a per-class low-rank covariance model',&
        &'splits the 2D classes into ncls subclasses in total. The classes are grouped by class-average correlation &
        &(nclust groups); a group gets subclasses in proportion to the square root of its population, counting the &
        &classes it already has, and its classes in proportion to their populations, at most pop/minpop each; a class &
        &with fewer than 2*minpop particles (or than 3*(neigs+1), the model floor) is kept whole. Per split class, &
        &its particles masked by the box''s inscribed circle: a CTF- and noise-weighted low-rank &
        &covariance model (rank neigs, pose fixed, pose tangents and contrast as nuisance) on a data-derived band. &
        &The latent of every particle and its uncertainty weigh the other particles (every particle''s pre-image), &
        &the latent rank and the bandwidth are chosen by predicting held-out Fourier coefficients, the subclasses are &
        &the pre-images of as many particles as the budget, chosen to cover the latent distribution (preferring &
        &pre-images of at least minpop effective particles), and their averages are kernel-weighted, CTF-corrected &
        &means of the particles (cls_expansion_cavgs.mrc, each filtered by its own even/odd FRC, plus the even/odd &
        &halves; the FRCs as class FRCs); per subclass the population, the effective population (also per half), the &
        &cross-validation z-score of the pre-images over the class mean and the cross-half reproducibility &
        &(os_cls2D pop, neff, sep, repro); a subclass half without members is written inactive (os_cls2D &
        &state_even/state_odd 0)',&
        &'simple_exec',&
        &.true., &
        &visibility=UI_VIS_ADVANCED, display_name='Class Expansion')
        call cls_expansion%add_input(UI_PARM, 'class', 'num', &
            'Optional class index to split', &
            'Optional 2D class index to split; omit to process all classes', &
            'e.g. 5', .false., 0.0, &
        &visibility=UI_VIS_ADVANCED)
        call cls_expansion%add_input(UI_PARM, 'ncls', 'num', &
            'Total number of subclasses', &
            'Total number of subclasses over all classes, at least the number of classes', &
            '# subclasses', .true., 0.0, &
        &visibility=UI_VIS_STANDARD)
        call cls_expansion%add_input(UI_PARM, 'minpop', 'num', &
            'Minimum particles per subclass', &
            'Minimum particles per subclass: a class holds at most pop/minpop subclasses and is kept whole below &
            &2*minpop particles (and below 3*(neigs+1)){25}', &
            '# particles{25}', .false., 25.0, &
        &visibility=UI_VIS_ADVANCED)
        call cls_expansion%add_input(UI_PARM, 'nclust', 'num', &
            'Number of class-average groups', &
            'Number of groups of similar class averages the subclass budget is spread over, by average linkage on &
            &their aligned correlation; without class averages every class is its own group{20}', &
            '# groups{20}', .false., 20.0, &
        &visibility=UI_VIS_ADVANCED)
        call cls_expansion%add_input(UI_FILT, 'neigs', 'num', &
            'Rank of the per-class covariance model (default 16)', &
            'Rank of the per-class covariance model (default 16)', &
            '# components', .false., 16.0, &
        &visibility=UI_VIS_ADVANCED)
        call cls_expansion%add_input(UI_FILT, 'lp', 'num', &
            'Fit band limit in A (default: data-derived per class)', &
            'The covariance model is fitted up to this resolution; without it the band is derived per class from the &
            &even/odd FRCs, at most 20 A; the sub-class averages use the full band', &
            'low-pass limit(in A)', .false., 20.0, &
        &visibility=UI_VIS_ADVANCED)
        call cls_expansion%add_input(UI_MASK, mskdiam, required_override=.true., group="mask", visibility=UI_VIS_STANDARD)
        call cls_expansion%add_input(UI_COMP, nparts, required_override=.false., group="compute", visibility=UI_VIS_STANDARD)
        call cls_expansion%add_input(UI_COMP, nthr,   group="compute", visibility=UI_VIS_STANDARD)
        call add_ui_program('cls_expansion', cls_expansion, prgtab, UI_CATEGORY)
    end subroutine new_cls_expansion

end module simple_ui_refine2D
