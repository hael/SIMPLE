!@descr: module defining the user interfaces for 2D clustering and related programs in the simple_exec suite
module simple_ui_cluster2D
use simple_ui_modules
implicit none

type(category_descriptor), parameter :: UI_CATEGORY = category_descriptor('cluster2d', 'Cluster2D Workflows', 30)
type(ui_program), target :: abinitio2D
type(ui_program), target :: cluster2D
type(ui_program), target :: abinitio2D_chunks
type(ui_program), target :: make_cavgs
type(ui_program), target :: bootstrap_cavgs
type(ui_program), target :: unbootstrap_cavgs
type(ui_program), target :: map_cavgs_selection
type(ui_program), target :: sample_classes
type(ui_program), target :: write_classes
type(ui_program), target :: cls_split

contains

    subroutine construct_cluster2D_programs( prgtab )
        class(ui_hash), intent(inout) :: prgtab
        call new_abinitio2D(prgtab)
        call new_abinitio2D_chunks(prgtab)
        call new_make_cavgs(prgtab)
        call new_bootstrap_cavgs(prgtab)
        call new_unbootstrap_cavgs(prgtab)
        call new_map_cavgs_selection(prgtab)
        call new_sample_classes(prgtab)
        call new_write_classes(prgtab)
        call new_cls_split(prgtab)
    end subroutine construct_cluster2D_programs

    subroutine new_abinitio2D( prgtab )
        class(ui_hash), intent(inout) :: prgtab
        call new_abinitio2D_descriptor(abinitio2D, 'abinitio2D', &
            &'Generate reference-free 2D class averages from particle images', &
            &'is a distributed workflow for generating 2D class averages from particles', &
            &UI_VIS_STANDARD, 'Create 2D Class Averages', 'yes')
        call add_ui_program('abinitio2D', abinitio2D, prgtab, UI_CATEGORY)
    end subroutine new_abinitio2D

    subroutine new_abinitio2D_descriptor( program, name, summary, help, visibility, display_name, &
        &inpl_cont_default )
        type(ui_program), intent(inout) :: program
        character(len=*), intent(in) :: name, summary, help, display_name
        integer, intent(in) :: visibility
        character(len=*), intent(in) :: inpl_cont_default
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
        &center of gravity and map shifts back to the particles(yes|no){no}','', .false., 'no', group="model", &
        &choices=ui_choices([character(len=3) :: 'yes', 'no']), &
        &visibility=UI_VIS_ADVANCED)
        call program%add_input(UI_SRCH, 'autoscale', 'binary', 'Automatic down-scaling', 'Automatic down-scaling of images &
        &for accelerated computation(yes|no){yes}','', .false., 'yes', group="model", &
        &choices=ui_choices([character(len=3) :: 'yes', 'no']), &
        &visibility=UI_VIS_ADVANCED)
        call program%add_input(UI_SRCH, 'refine', 'multi', 'Refinement mode',&
        &'Refinement mode(prob_snhc|prob|snhc_smpl){prob_snhc}','', .false., 'prob_snhc', group="search", &
        &choices=ui_choices([character(len=9) :: 'prob_snhc', 'prob', 'snhc_smpl']), &
        &visibility=UI_VIS_ADVANCED)
        call program%add_input(UI_SRCH, 'inpl_cont', 'binary', &
        &'Continuous in-plane refinement', &
        &'Joint continuous Euclidean in-plane and shift refinement(yes|no){'//trim(inpl_cont_default)//'}', '', &
        &.false., trim(inpl_cont_default), group="search", &
        &choices=ui_choices([character(len=3) :: 'yes', 'no']), &
        &visibility=UI_VIS_ADVANCED)
        call program%add_input(UI_SRCH, 'sigma_est', 'multi', 'Sigma estimation method',&
        &'Sigma estimation method(group|global){global}','', .false., 'global', group="search", &
        &choices=ui_choices([character(len=6) :: 'group', 'global']), &
        &visibility=UI_VIS_ADVANCED)
        call program%add_input(UI_SRCH, cls_init, group="search", &
        &visibility=UI_VIS_ADVANCED)
        call program%add_input(UI_SRCH, nsample, group="search", visibility=UI_VIS_ADVANCED)
        call program%add_input(UI_SRCH, 'stats', 'binary', 'Write class-average statistics', &
        &'Write class-average statistics during ab-initio 2D analysis(yes|no){no}', '', .false., 'no', group="search", &
        &choices=ui_choices([character(len=3) :: 'yes', 'no']), &
        &visibility=UI_VIS_DEVELOPER)
        call program%add_input(UI_SRCH, 'extr_lim', 'num', 'Extreme-particle iteration limit', &
        &'Iteration limit used to identify and exclude extreme particles during ab-initio 2D analysis', &
        &'iteration limit', .false., 20., group="search", visibility=UI_VIS_DEVELOPER)
        call program%add_input(UI_SRCH, 'nits_per_stage', 'num', 'Iterations per stage', &
        &'Number of ab-initio 2D iterations added at each frequency-marching stage', &
        &'iterations per stage', .false., 5., group="search", visibility=UI_VIS_DEVELOPER)
        call program%add_input(UI_SRCH, 'eo_stage', 'binary', 'Use even/odd stage', &
         &'Run the final ab-initio 2D frequency-marching stage with even/odd class averages(yes|no){yes}', '', &
         &.false., 'yes', group="search", choices=ui_choices([character(len=3) :: 'yes', 'no']), &
         &visibility=UI_VIS_DEVELOPER)
        ! filter controls
        call program%add_input(UI_FILT, hp, group="filter", &
        &visibility=UI_VIS_ADVANCED)
        call program%add_input(UI_FILT, ml_reg, group="filter", &
        &visibility=UI_VIS_ADVANCED)
        call program%add_input(UI_FILT, 'cenlp', 'num', 'Centering low-pass limit', 'Limit for low-pass filter used in binarisation &
        &prior to determination of the center of gravity of the reference volume(s) and centering', 'centering low-pass limit in &
        &Angstroms{30}', .false., 30., group="filter", &
        &visibility=UI_VIS_ADVANCED)
        call program%add_input(UI_FILT, 'lpstart', 'num', 'Initial low-pass limit', 'Initial low-pass resolution limit for the first stage of ab-initio model generation',&
            &'low-pass limit in Angstroms', .false., 30., group="filter", &
        &visibility=UI_VIS_ADVANCED)
        call program%add_input(UI_FILT, 'lpstop',  'num', 'Final low-pass limit', 'Final low-pass limit',&
            &'low-pass limit for the second stage (no e/o cavgs refinement) in Angstroms', .false., 6., group="filter", &
        &visibility=UI_VIS_ADVANCED)
        call program%add_input(UI_FILT, lp, group="filter", &
        &visibility=UI_VIS_ADVANCED)
        ! mask controls
        call program%add_input(UI_MASK, mskdiam, group="mask", visibility=UI_VIS_STANDARD)
        ! computer controls
        call program%add_input(UI_COMP, nparts, required_override=.false., group="compute", visibility=UI_VIS_STANDARD)
        call program%add_input(UI_COMP, nthr, group="compute", visibility=UI_VIS_STANDARD)
    end subroutine new_abinitio2D_descriptor

    subroutine new_abinitio2D_chunks( prgtab )
        class(ui_hash), intent(inout) :: prgtab
        ! PROGRAM SPECIFICATION
        call abinitio2D_chunks%new(&
        &'abinitio2D_chunks',&                                                                    ! name
        &'Independent ab initio 2D analysis of particle subsets',&                                ! summary
        &'splits a project into particle-balanced subsets and runs independent abinitio2D jobs',& ! help
        &'simple_exec',&                                                                          ! executable
        &.true.,&                                                                                 ! requires sp_project
        &visibility=UI_VIS_ADVANCED)
        ! INPUT PARAMETER SPECIFICATIONS
        ! image input/output
        ! <empty>
        ! parameter input/output
        call abinitio2D_chunks%add_input(UI_PARM, 'nchunks', 'num', 'Number of chunks', &
            &'Number of particle-balanced subset projects to run with independent abinitio2D jobs. &
            &Omit or set to 0 to target about 100 classes per chunk.', &
            &'# of chunks (0=auto)', .false., 0., &
        &visibility=UI_VIS_ADVANCED)
        call abinitio2D_chunks%add_input(UI_PARM, 'nchunks_in_parallel', 'num', 'Number of chunks run in parallel', &
            &'Maximum number of chunks executed at the same time. Each chunk uses nparts parts, so the &
            &peak load is about nchunks_in_parallel*nparts processes. Defaults to 1 (one chunk at a time). &
            &Set to 0 to run all chunks in parallel.', &
            &'# of parallel chunks (0=all){1}', .false., 1., &
        &visibility=UI_VIS_ADVANCED)
        ! <no additional inputs>
        ! <empty>
        ! search controls
        call abinitio2D_chunks%add_input(UI_SRCH, nptcls_per_cls, placeholder_override='# of particles per cluster{500}', group="cluster 2D", visibility=UI_VIS_STANDARD)
        call abinitio2D_chunks%add_input(UI_SRCH, 'center', 'binary', 'Center class averages', 'Center class averages by their center of &
            &gravity and map shifts back to the particles(yes|no){yes}','', .false., 'yes', group="cluster 2D", &
        &choices=ui_choices([character(len=3) :: 'yes', 'no']), &
        &visibility=UI_VIS_ADVANCED)
        call abinitio2D_chunks%add_input(UI_SRCH, 'refine', 'multi', 'Refinement mode',&
        &'Refinement mode(prob_snhc|prob|snhc_smpl){prob_snhc}','', .false., 'prob_snhc', group="cluster 2D", &
        &choices=ui_choices([character(len=9) :: 'prob_snhc', 'prob', 'snhc_smpl']), &
        &visibility=UI_VIS_ADVANCED)
        ! filter controls
        call abinitio2D_chunks%add_input(UI_FILT, hp, group="cluster 2D", &
        &visibility=UI_VIS_ADVANCED)
        call abinitio2D_chunks%add_input(UI_FILT, lp, group="cluster 2D", &
        &visibility=UI_VIS_ADVANCED)
        call abinitio2D_chunks%add_input(UI_FILT, 'cenlp', 'num', 'Centering low-pass limit', 'Limit for low-pass filter used in binarisation &
        &prior to determination of the center of gravity of the class averages and centering', 'centering low-pass limit in &
        &Angstroms{30}', .false., 30., group="cluster 2D", &
        &visibility=UI_VIS_ADVANCED)
        call abinitio2D_chunks%add_input(UI_FILT, 'lpstop', 'num', 'Final low-pass limit', 'Low-pass limit that controls the degree of &
        &downsampling in the second phase. Give estimated best final resolution', 'final low-pass limit in Angstroms', .false., 8.,&
        &group="filter", visibility=UI_VIS_ADVANCED)
        ! mask controls
        call abinitio2D_chunks%add_input(UI_MASK, mskdiam, group="cluster 2D", visibility=UI_VIS_STANDARD)
        ! computer controls
        call abinitio2D_chunks%add_input(UI_COMP, nparts, required_override=.false., group="compute", visibility=UI_VIS_STANDARD)
        call abinitio2D_chunks%add_input(UI_COMP, nthr, group="compute", visibility=UI_VIS_STANDARD)
        call abinitio2D_chunks%add_input(UI_COMP, 'walltime', 'num', 'Walltime', 'Maximum execution time for job scheduling and &
        &management(29mins){1740}', 'in seconds(29mins){1740}', .false., 1740., group="compute", &
        &visibility=UI_VIS_ADVANCED)
        ! add to ui_hash
        call add_ui_program('abinitio2D_chunks', abinitio2D_chunks, prgtab, UI_CATEGORY)
    end subroutine new_abinitio2D_chunks

    subroutine new_make_cavgs( prgtab )
        class(ui_hash), intent(inout) :: prgtab
        ! PROGRAM SPECIFICATION
        call make_cavgs%new(&
        &'make_cavgs', &                           ! name
        &'Create class averages from aligned particle images',& ! summary
        &'is a distributed workflow for generating class averages or initial random references&
        & for cluster2D execution',&               ! help
        &'simple_exec',&                           ! executable
        &.true., &
        &visibility=UI_VIS_ADVANCED)                                   ! requires sp_project
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
        &.true.)                                                               ! requires sp_project
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
        &.true.)                                                                ! requires sp_project
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
        &visibility=UI_VIS_ADVANCED)                                                         ! requires sp_project
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
        &visibility=UI_VIS_ADVANCED)                                                                            ! requires sp_project
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
        &.true.)                                                                                            ! requires sp_project
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

    subroutine new_cls_split( prgtab )
        class(ui_hash), intent(inout) :: prgtab
        call cls_split%new(&
        &'cls_split',&
        &'Split classes with latent clustering',&
        &'splits 2D/3D particle classes into subclasses using diffusion-map or kPCA embeddings and k-medoids clustering',&
        &'simple_exec',&
        &.true., &
        &visibility=UI_VIS_ADVANCED)
        call cls_split%add_input(UI_PARM, 'class', 'num', &
            'Optional class index to split', &
            'Optional 2D class index or 3D projection/class index to split; omit to process all classes', &
            'e.g. 5', .false., 0.0, &
        &visibility=UI_VIS_ADVANCED)
        call cls_split%add_input(UI_PARM, 'ncls', 'num', &
            'Fixed number of subclasses (0 => auto)', &
            'Fixed number of subclasses (0 => auto)', &
            '# subclasses', .false., 0.0, &
        &visibility=UI_VIS_ADVANCED)
        call cls_split%add_input(UI_PARM, 'nsubcls_min', 'num', &
            'Minimum subclass trial count in auto mode (default 3)', &
            'Used only when ncls=0: optimization tries every subclass count from nsubcls_min through nsubcls_max', &
            '# min trial subclasses', .false., 3.0, &
        &visibility=UI_VIS_ADVANCED)
        call cls_split%add_input(UI_PARM, 'nsubcls_max', 'num', &
            'Maximum subclass trial count in auto mode (default 10)', &
            'Used only when ncls=0: optimization tries every subclass count from nsubcls_min through nsubcls_max', &
            '# max trial subclasses', .false., 10.0, &
        &visibility=UI_VIS_ADVANCED)
        call cls_split%add_input(UI_PARM, 'k_nn', 'num', &
            'Diffusion graph neighbors (default 10; try 5-30)', &
            'Local nearest neighbors used only for diffusion-map modes; larger values smooth the graph', &
            '# neighbors', .false., real(DIFFMAP_GRAPH_KNN_DEFAULT), &
        &visibility=UI_VIS_ADVANCED)
        call cls_split%add_input(UI_PARM, 'oritype', 'multi', &
            'Particle type to split', 'Particle type to split(ptcl2D|ptcl3D){ptcl2D}','', .false., 'ptcl2D', &
        &choices=ui_choices([character(len=6) :: 'ptcl2D', 'ptcl3D']), &
        &visibility=UI_VIS_ADVANCED)
        call cls_split%add_input(UI_FILT, 'neigs', 'num', &
            'Number of eigencomponents (0 => auto scan; default 200)', &
            'Number of eigencomponents used as the scan upper bound before ICM dimension selection', &
            '# eigenvecs', .false., real(DIFFMAP_NEIGS_SCAN_DEFAULT), &
        &visibility=UI_VIS_ADVANCED)
        call cls_split%add_input(UI_FILT, 'pca_mode', 'multi', &
            'Class split embedding method', &
            'Class split embedding method(diffusion_maps|kpca){diffusion_maps}','', .false., 'diffusion_maps', &
        &choices=ui_choices([character(len=14) :: 'diffusion_maps', 'kpca']), &
        &visibility=UI_VIS_ADVANCED)
        call cls_split%add_input(UI_FILT, 'graph', 'multi', &
            'Class split graph', 'Class split graph(euc|ori){euc}','', .false., 'euc', &
        &choices=ui_choices([character(len=3) :: 'euc', 'ori']), &
        &visibility=UI_VIS_ADVANCED)
        call cls_split%add_input(UI_FILT, 'bandwidth_mode', 'multi', &
            'Diffusion-map bandwidth mode', &
            'Kernel bandwidth policy for diffusion maps(median|ferguson){median}','', .false., 'median', &
        &choices=ui_choices([character(len=8) :: 'median', 'ferguson']), &
        &visibility=UI_VIS_ADVANCED)
        call cls_split%add_input(UI_FILT, 'bandwidth_tune', 'num', &
            'Ferguson bandwidth multiplier (default 1)', &
            'Linear multiplier of the Ferguson-optimal kernel bandwidth (1=optimum); only used when bandwidth_mode=ferguson', &
            'tune >= 0', .false., 1.0, &
        &visibility=UI_VIS_ADVANCED)
        call cls_split%add_input(UI_FILT, 'dm_alpha', 'num', &
            'Diffusion-map density normalization (default 0)', &
            'Coifman-Lafon alpha: 0=graph Laplacian, 0.5=Fokker-Planck, 1=Laplace-Beltrami (divides out sampling density)', &
            '0 <= alpha <= 1', .false., 0.0, &
        &visibility=UI_VIS_ADVANCED)
        call cls_split%add_input(UI_MASK, mskdiam, required_override=.false., group="mask", visibility=UI_VIS_STANDARD)
        call cls_split%add_input(UI_COMP, nparts, required_override=.false., group="compute", visibility=UI_VIS_STANDARD)
        call cls_split%add_input(UI_COMP, nthr,   group="compute", visibility=UI_VIS_STANDARD)
        call add_ui_program('cls_split', cls_split, prgtab, UI_CATEGORY)
    end subroutine new_cls_split

end module simple_ui_cluster2D
