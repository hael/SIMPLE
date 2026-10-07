!@descr: module defining the user interfaces for streaming workflows in the simple_exec suite
module simple_ui_stream
use simple_ui_modules
implicit none

type(category_descriptor), parameter :: UI_CATEGORY = category_descriptor('stream', 'Stream Workflows', 10)
type(ui_program), target :: pool2D
type(ui_program), target :: solve3D_stream
type(ui_program), target :: assign_optics
type(ui_program), target :: gen_pickrefs
type(ui_program), target :: master
type(ui_program), target :: pick_extract
type(ui_program), target :: preproc
type(ui_program), target :: sieve_cavgs

contains

    subroutine construct_stream_programs(prgtab)
        class(ui_hash), intent(inout) :: prgtab
        call new_pool2D(prgtab)
        call new_solve3D_stream(prgtab)
        call new_assign_optics(prgtab)
        call new_gen_pickrefs(prgtab)
        call new_master(prgtab)
        call new_pick_extract(prgtab)
        call new_preproc(prgtab)
        call new_sieve_cavgs(prgtab)
    end subroutine construct_stream_programs

subroutine new_pool2D( prgtab )
        class(ui_hash), intent(inout) :: prgtab
        ! PROGRAM SPECIFICATION
        call pool2D%new(&
        &'pool2D', &                                                             ! name
        &'Run streaming 2D analysis as new data arrive',& ! summary
        &'is a distributed workflow that executes 2D analysis'//&                ! help
        &' in streaming mode as the microscope collects the data',&
        &'simple_stream',&                                                       ! executable
        &.true.,&                                                                ! requires sp_project
        &visibility=UI_VIS_DEVELOPER, display_name='Run streaming 2D analysis as new data arrive')
        ! image input/output
        ! <empty>
        ! parameter input/output
        call pool2D%add_input(UI_FILE, 'dir_target', 'file', 'Target directory',&
        &'Directory where the sieve_cavgs application is running', 'e.g. 3_sieve_cavgs', .true., '', group="data", visibility=UI_VIS_STANDARD)
        call pool2D%add_input(UI_FILE, 'dir_exec', 'file', 'Previous run directory',&
        &'Directory where previous 2D analysis took place', 'e.g. 3_pool2D', .false., '', group="data", &
        &visibility=UI_VIS_DEVELOPER)
        ! <no additional inputs>
        ! <empty>
        ! search controls
        call pool2D%add_input(UI_SRCH, 'ncls', 'num', 'Maximum number of 2D clusters',&
        &'Maximum number of 2D class averages for the pooled particles subsets', 'Maximum # 2D clusters', .true., 200., group="cluster 2D",&
        &visibility=UI_VIS_STANDARD)
        call pool2D%add_input(UI_SRCH, 'center', 'binary', 'Center class averages', &
        &'Center class averages by their center of gravity and map shifts back to the particles(yes|no){yes}', '', .false., 'yes', group="cluster 2D", &
        &choices=ui_choices([character(len=3) :: 'yes', 'no']), &
        &visibility=UI_VIS_ADVANCED)
        ! preserve_default: the stage default lives in the p06 commander's set_pool2D_cline, which default_audit.py does not trace
        call pool2D%add_input(UI_SRCH, 'stepwise', 'binary', 'Stepwise set import', &
        &'Each import takes only as many sieved particle sets as its particles need to reach the particle threshold; &
        &the other sets wait for a later import(yes|no){yes}', '', .false., 'yes', group="cluster 2D", &
        &choices=ui_choices([character(len=3) :: 'yes', 'no']), &
        &visibility=UI_VIS_DEVELOPER, preserve_default=.true.)
        ! preserve_default: the pool's fallback (seg) lives in its utilities, which default_audit.py does not trace
        call pool2D%add_input(UI_SRCH, 'center_type', 'multi', 'Centering scheme', &
        &'How class averages are centered: by their mass, by segmentation, or from the parameters(mass|seg|params){seg}', &
        &'', .false., 'seg', group="cluster 2D", &
        &choices=ui_choices([character(len=6) :: 'mass', 'seg', 'params']), &
        &visibility=UI_VIS_DEVELOPER, preserve_default=.true.)
        ! preserve_default: the fallback is STREAM_NPTCLS_MAX in the pool utilities, which default_audit.py does not trace
        call pool2D%add_input(UI_SRCH, 'nsample_max', 'num', 'Particles before fractional updates', &
        &'Number of selected particles in the pool beyond which each 2D iteration updates only a fraction of them', &
        &'# particles', .false., real(STREAM_NPTCLS_MAX), group="cluster 2D", &
        &visibility=UI_VIS_DEVELOPER, preserve_default=.true.)
        call pool2D%add_input(UI_SRCH, update_frac, help_override='Fraction of the particles updated per iteration once &
        &the pool uses fractional updates (beyond nsample_max particles); when not given, it is derived from the particles &
        &added since the previous iteration', group="cluster 2D", visibility=UI_VIS_DEVELOPER)
        ! filter controls
        ! preserve_default: the stage default lives in the p06 commander's set_pool2D_cline, which default_audit.py does not trace
        call pool2D%add_input(UI_FILT, 'dynreslim', 'binary', 'Dynamic resolution limit', &
        &'Enlarge the working images of the pool once its resolution has stayed at their Nyquist limit, &
        &which allows a finer resolution(yes|no){yes}', '', .false., 'yes', group="cluster 2D", &
        &choices=ui_choices([character(len=3) :: 'yes', 'no']), &
        &visibility=UI_VIS_DEVELOPER, preserve_default=.true.)
        call pool2D%add_input(UI_FILT, 'lpstop', 'num', 'Hard resolution limit', &
        &'Hard resolution limit of the pool 2D analysis (in Angstroms), never finer than the Nyquist limit of the &
        &working images; when not given, that Nyquist limit', 'low-pass limit in Angstroms', .false., 8., group="cluster 2D", &
        &visibility=UI_VIS_DEVELOPER)
        call pool2D%add_input(UI_FILT, 'cenlp', 'num', 'Centering low-pass limit', &
        &'Low-pass limit (in Angstroms) applied before the class averages are centered; when not given, it is derived &
        &from the mask diameter', 'centering low-pass limit in Angstroms', .false., 20., group="cluster 2D", &
        &visibility=UI_VIS_DEVELOPER)
        ! mask controls
        call pool2D%add_input(UI_MASK, 'mskdiam', 'num', 'Mask diameter', 'Mask diameter (in A) for application of a soft-edged circular mask to &
        &remove background noise', 'mask diameter in A', .false., 0., group="cluster 2D", visibility=UI_VIS_STANDARD)
        ! computer controls
        call pool2D%add_input(UI_COMP, nparts, group="compute", visibility=UI_VIS_STANDARD)
        call pool2D%add_input(UI_COMP, nthr, group="compute", visibility=UI_VIS_STANDARD)
        call pool2D%add_input(UI_COMP, 'walltime', 'num', 'Walltime', 'Maximum execution time for job scheduling and management in seconds{1740}(29mins)',&
        &'in seconds(29mins){1740}', .false., 1740., group="compute", &
        &visibility=UI_VIS_DEVELOPER)
        ! add to ui_hash
        call add_ui_program('pool2D', pool2D, prgtab, UI_CATEGORY)
    end subroutine new_pool2D
    
    subroutine new_solve3D_stream( prgtab )
        class(ui_hash), intent(inout) :: prgtab
        ! PROGRAM SPECIFICATION
        call solve3D_stream%new(&
        &'solve3D_stream', &                                                     ! name
        &'Run streaming 3D analysis as new data arrive',& ! summary
        &'is a distributed workflow that executes 3D analysis'//&                ! help
        &' in streaming mode as the microscope collects the data',&
        &'simple_stream',&                                                       ! executable
        &.true.,&                                                                ! requires sp_project
        &visibility=UI_VIS_DEVELOPER, display_name='Run streaming 3D analysis as new data arrive')
        ! image input/output
        ! <empty>
        ! parameter input/output
        call solve3D_stream%add_input(UI_FILE, 'dir_target', 'file', 'Target directory',&
        &'Directory where the pool2D application is running', 'e.g. 4_pool2D', .true., '', group="data", visibility=UI_VIS_STANDARD)
        ! <no additional inputs>
        ! <empty>
        ! preserve_default: the stage defaults below live in the p07 commander's set_solve3D_cline,
        ! which default_audit.py does not trace
        ! search controls
        call solve3D_stream%add_input(UI_SRCH, 'nstates', 'num', 'Number of states', &
        &'Number of states reconstructed by each streaming solve3D run{3}', '# states', .false., 3., group="search", &
        &visibility=UI_VIS_DEVELOPER, preserve_default=.true.)
        call solve3D_stream%add_input(UI_SRCH, 'nstages', 'num', 'Number of solve3D stages', &
        &'Number of low-pass limit stages of each streaming solve3D run{5}', '# stages', .false., 5., group="search", &
        &visibility=UI_VIS_DEVELOPER, preserve_default=.true.)
        call solve3D_stream%add_input(UI_SRCH, 'nptcls3D_max', 'num', 'Particles of the first solve3D', &
        &'Maximum number of particles of the first streaming solve3D run, drawn class-balanced; the others go to the first solve3D_addon run{100000}', &
        &'max # particles', .false., 100000., group="search", visibility=UI_VIS_DEVELOPER, preserve_default=.true.)
        ! filter controls
        call solve3D_stream%add_input(UI_FILT, 'lpstart', 'num', 'Starting low-pass limit', &
        &'Low-pass limit of the first stage of each streaming solve3D run (in Angstroms){50}', 'low-pass limit in Angstroms', &
        &.false., 50., group="filter", visibility=UI_VIS_DEVELOPER, preserve_default=.true.)
        call solve3D_stream%add_input(UI_FILT, 'lpstop', 'num', 'Final low-pass limit', &
        &'Low-pass limit of the last stage of each streaming solve3D run (in Angstroms){10}', 'low-pass limit in Angstroms', &
        &.false., 10., group="filter", visibility=UI_VIS_DEVELOPER, preserve_default=.true.)
        ! mask controls
        ! computer controls
        call solve3D_stream%add_input(UI_COMP, nthr, group="compute", visibility=UI_VIS_STANDARD)
        call solve3D_stream%add_input(UI_COMP, 'nparts3D', 'num', 'Partitions of each 3D job', &
        &'Number of partitions of each solve3D and solve3D_addon job{8}', '# partitions', .false., 8., group="compute", &
        &visibility=UI_VIS_DEVELOPER, preserve_default=.true.)
        call solve3D_stream%add_input(UI_COMP, 'nthr3D', 'num', 'Threads of each 3D job', &
        &'Number of OpenMP threads of each solve3D and solve3D_addon job{8}', '# threads', .false., 8., group="compute", &
        &visibility=UI_VIS_DEVELOPER, preserve_default=.true.)
        call solve3D_stream%add_input(UI_COMP, 'walltime', 'num', 'Walltime', 'Maximum execution time for job scheduling and management in seconds{1740}(29mins)',&
        &'in seconds(29mins){1740}', .false., 1740., group="compute", &
        &visibility=UI_VIS_DEVELOPER)
        ! add to ui_hash
        call add_ui_program('solve3D_stream', solve3D_stream, prgtab, UI_CATEGORY)
    end subroutine new_solve3D_stream

    subroutine new_assign_optics( prgtab )
        class(ui_hash), intent(inout) :: prgtab
        ! PROGRAM SPECIFICATION
        call assign_optics%new(&
        &'assign_optics', &                                              ! name
        &'Assign optics groups from microscope metadata',& ! summary
        &'is a program to assign optics groups during streaming',&       ! descr long
        &'simple_stream',&                                               ! executable
        &.true., &
        &visibility=UI_VIS_DEVELOPER, display_name='Assign optics groups from microscope metadata')                                                         ! requires sp_project
        ! INPUT PARAMETER SPECIFICATIONS
        ! image input/output
        ! parameter input/output
        call assign_optics%add_input(UI_FILE, 'dir_target', 'file', 'Target directory',&
        &'Directory where the preprocess_stream application is running', 'e.g. 1_preproc', .true., '', &
        &visibility=UI_VIS_STANDARD)
        call assign_optics%add_input(UI_PARM, 'nmics', 'num', 'Micrographs before termination', &
        &'Number of micrographs after which the optics assignment terminates; 0 = no limit{0}', '# micrographs{0}', .false., 0., &
        &visibility=UI_VIS_DEVELOPER)
        ! <no additional inputs>
        ! <empty>
        ! search controls
        call assign_optics%add_input(UI_SRCH, 'tilt_thres', 'num', 'Beam-shift clustering threshold', &
        &'Distance threshold of the hierarchical clustering of the beam-image shifts into optics groups{0.05}', 'e.g 0.05', &
        &.false., 0.05, group="optics groups", visibility=UI_VIS_DEVELOPER)
        call assign_optics%add_input(UI_SRCH, 'beamtilt', 'binary', 'Use beam-tilt groups', &
        &'Split the micrographs by beam-tilt group before their beam-image shifts are clustered into optics groups(yes|no){no}', &
        &'', .false., 'no', group="optics groups", &
        &choices=ui_choices([character(len=3) :: 'yes', 'no']), &
        &visibility=UI_VIS_DEVELOPER)
        ! filter controls
        ! <empty>
        ! mask controls
        ! <empty>
        ! computer controls
        call assign_optics%add_input(UI_COMP, nthr, group="compute", visibility=UI_VIS_STANDARD)
        ! add to ui_hash
        call add_ui_program('assign_optics', assign_optics, prgtab, UI_CATEGORY)
    end subroutine new_assign_optics

    subroutine new_gen_pickrefs( prgtab )
        class(ui_hash), intent(inout) :: prgtab
        ! PROGRAM SPECIFICATION
        call gen_pickrefs%new(&
        &'gen_pickrefs', &                                               ! name
        &'Do a mini stream to create the opening 2D for generation of picking references',&  ! summary
        &'is a program to do a mini stream to create the opening 2D',&   ! descr long
        &'simple_stream',&                                               ! executable
        &.true., &
        &visibility=UI_VIS_ADVANCED, &
        &display_name='Do a mini stream to create the opening 2D for generation of picking references')                                                         ! requires sp_project
        ! INPUT PARAMETER SPECIFICATIONS
        ! image input/output
        ! parameter input/output
        call gen_pickrefs%add_input(UI_FILE, 'dir_target', 'file', 'Target directory',&
        &'Directory where the preprocess_stream application is running', 'e.g. 1_preproc', .true., '', &
        &visibility=UI_VIS_STANDARD)
        call gen_pickrefs%add_input(UI_PARM, pcontrast, group="picking", &
        &visibility=UI_VIS_DEVELOPER)
        call gen_pickrefs%add_input(UI_SRCH, nptcls_per_cls, group="cluster 2D", &
        &visibility=UI_VIS_DEVELOPER)
        ! <no additional inputs>
        ! <empty>
        ! search controls
        ! filter controls
        call gen_pickrefs%add_input(UI_FILT, 'amsklp', 'num', 'Automask low-pass limit', &
        &'Low-pass limit used before opening-2D automask generation', 'low-pass limit in Angstroms{20}', .false., 20., group="picking", &
        &visibility=UI_VIS_DEVELOPER)
        ! mask controls
        call gen_pickrefs%add_input(UI_MASK, 'ngrow', 'num', 'Automask growth layers', &
        &'Number of binary-image layers grown during opening-2D automasking', '# layers{3}', .false., 3., group="picking", &
        &visibility=UI_VIS_DEVELOPER)
        call gen_pickrefs%add_input(UI_MASK, 'winsz', 'num', 'Automask window size', &
        &'Window size used during opening-2D automask estimation', 'window size{5}', .false., 5., group="picking", &
        &visibility=UI_VIS_DEVELOPER)
        call gen_pickrefs%add_input(UI_MASK, 'edge', 'num', 'Automask soft edge', &
        &'Cosine edge width used to soften the opening-2D automask', '# pixels{6}', .false., 6., group="picking", &
        &visibility=UI_VIS_DEVELOPER)
        ! the 3D route to the picking references (the master forwards its own); 0 keeps the default. preserve_default: the defaults of nrestarts_collapse, lpstart_ini3D and
        ! lpstop_ini3D live in the initial analysis stage, which default_audit.py does not trace
        ! search controls
        call gen_pickrefs%add_input(UI_SRCH, 'nstates_pickrefs', 'int',         'Picking-reference 3D states', &
        &'Number of states of the solve3D_cavgs run in the 3D route of the initial analysis (picking references); &
        &0 uses its default of 3', '0', .false., '', group="3D route", visibility=UI_VIS_DEVELOPER)
        call gen_pickrefs%add_input(UI_SRCH, 'nstages_pickrefs', 'int',         'Picking-reference 3D stages', &
        &'Number of solve3D_cavgs stages in the 3D route of the initial analysis (picking references); &
        &0 uses its default of 3', '0', .false., '', group="3D route", visibility=UI_VIS_DEVELOPER)
        call gen_pickrefs%add_input(UI_SRCH, 'nspace_pickrefs', 'int',          'Picking-reference reprojections', &
        &'Number of reprojections of the 3D route''s volume used as picking references by the initial analysis; &
        &0 uses its default of 50', '0', .false., '', group="3D route", visibility=UI_VIS_DEVELOPER)
        call gen_pickrefs%add_input(UI_SRCH, 'nrestarts_collapse', 'int',       'Picking-reference 3D restarts', &
        &'Number of solve3D_cavgs restarts when states collapse, in the 3D route of the initial analysis &
        &(picking references){3}', '3', .false., 3., group="3D route", visibility=UI_VIS_DEVELOPER, preserve_default=.true.)
        ! filter controls
        call gen_pickrefs%add_input(UI_FILT, 'lpstop_pickrefs', 'float',        'Picking-reference final low-pass limit', &
        &'Final low-pass limit (in Angstroms) of the solve3D_cavgs run in the 3D route of the initial analysis &
        &(picking references); 0 uses its default of 8', '0', .false., '', group="3D route", visibility=UI_VIS_DEVELOPER)
        call gen_pickrefs%add_input(UI_FILT, 'lpstart_ini3D', 'float',          'Picking-reference initial 3D starting low-pass limit', &
        &'Starting low-pass limit (in Angstroms) of the solve3D_cavgs initial 3D model in the 3D route of the &
        &initial analysis (picking references){100}', '100', .false., 100., group="3D route", visibility=UI_VIS_DEVELOPER, preserve_default=.true.)
        call gen_pickrefs%add_input(UI_FILT, 'lpstop_ini3D', 'float',           'Picking-reference initial 3D final low-pass limit', &
        &'Final low-pass limit (in Angstroms) of the solve3D_cavgs initial 3D model in the 3D route of the &
        &initial analysis (picking references){20}', '20', .false., 20., group="3D route", visibility=UI_VIS_DEVELOPER, preserve_default=.true.)
        ! mask controls
        ! computer controls
        call gen_pickrefs%add_input(UI_COMP, 'nthr3D_pickrefs', 'int',          'Picking-reference 3D threads', &
        &'Number of OpenMP threads of the solve3D_cavgs run and the reprojection in the 3D route of the initial analysis &
        &(picking references); 0 uses its default of 16', '0', .false., '', group="3D route", visibility=UI_VIS_DEVELOPER)
        ! computer controls
        call gen_pickrefs%add_input(UI_COMP, nthr, group="compute", visibility=UI_VIS_STANDARD)
        ! add to ui_hash
        call add_ui_program('gen_pickrefs', gen_pickrefs, prgtab, UI_CATEGORY)
    end subroutine new_gen_pickrefs

    subroutine new_master( prgtab )
        class(ui_hash), intent(inout) :: prgtab
        ! PROGRAM SPECIFICATION
        call master%new(&
        &'master', &                                                                 ! name
        &'Coordinate streaming jobs, metadata, and NICE communication',& ! summary
        &'master process that forks streaming programs, collates metadata,'//&       ! help
        &'communicates with Nice and provides job control',&
        &'simple_stream',&                                                           ! executable
        &.false., display_name='Coordinate streaming jobs, metadata, and NICE communication')                                                                     ! requires sp_project
        ! please note: globally declared inputs not used as allows custom descriptions for GUI
        ! image input/output
        call master%add_input(UI_FILE, 'dir_movies', 'dir',  'Input movies directory',   'Input movies directory',   '', .true.,  '', &
        &visibility=UI_VIS_STANDARD)
        call master%add_input(UI_FILE, 'dir_meta',   'dir',  'Input metadata directory', 'Input metadata directory', '', .false., '', &
        &visibility=UI_VIS_STANDARD)
        call master%add_input(UI_FILE, 'gainref',    'file', 'Gain reference',           'Gain reference',           '', .false., '', &
        &visibility=UI_VIS_STANDARD)
        ! parameter input/output
        call master%add_input(UI_PARM, 'flipgain',       'multi',         'Gain processing', 'Gain processing(none|flip_auto|flip_x|flip_y|flip_xy|generate){none}', '', .false., 'none', &
        &choices=ui_choices([character(len=9) :: 'none', 'flip_auto', 'flip_x', 'flip_y', 'flip_xy', 'generate']), &
        &visibility=UI_VIS_STANDARD)
        call master%add_input(UI_PARM, 'cs',             'float',  'Spherical aberration (mm)',   'Spherical aberration (mm)',   '2.7',                    .true.,  '', &
        &visibility=UI_VIS_STANDARD)
        call master%add_input(UI_PARM, 'fraca',          'float',  'Amplitude contrast fraction', 'Amplitude contrast fraction', '0.1',                    .true.,  '', &
        &visibility=UI_VIS_STANDARD)
        call master%add_input(UI_PARM, 'kv',             'int',    'Acceleration voltage (kV)',   'Acceleration voltage (kV)',   '300',                    .true.,  '', &
        &visibility=UI_VIS_STANDARD)
        call master%add_input(UI_PARM, 'smpd',           'float',  'Pixel size (A)',              'Pixel size (A)',              '',                       .true.,  '', &
        &visibility=UI_VIS_STANDARD)
        call master%add_input(UI_PARM, 'fit_phshift',    'binary', 'Fit CTF phase shift', &
        &'Fit the additive phase shift during CTF estimation (yes|no){no}', '', .false., 'no', &
        &choices=ui_choices([character(len=3) :: 'yes', 'no']), visibility=UI_VIS_STANDARD)
        call master%add_input(UI_PARM, 'phshift_min',    'float',         'Minimum CTF phase shift', 'Minimum fitted additive phase shift in degrees, 0-360; a window narrower than 180 degrees fixes the sign of the fitted CTF', '0', .false., '', visibility=UI_VIS_DEVELOPER)
        call master%add_input(UI_PARM, 'phshift_max',    'float',         'Maximum CTF phase shift', 'Maximum fitted additive phase shift in degrees, 0-360; fitting is blind to a 180-degree offset, so narrow the window around the expected phase', '180', .false., '', visibility=UI_VIS_DEVELOPER)
        call master%add_input(UI_PARM, 'phshift_step',   'float',         'CTF phase-shift step', 'Initial phase-shift grid step in degrees', '10', .false., '', visibility=UI_VIS_DEVELOPER)
        call master%add_input(UI_PARM, dfmin, visibility=UI_VIS_DEVELOPER)
        call master%add_input(UI_PARM, dfmax, visibility=UI_VIS_DEVELOPER)
        call master%add_input(UI_PARM, 'smpd_downscale', 'float', 'Downscaled pixel size (A)', &
        &'Downscaled pixel size (A)', '', .true., STREAM_DEFAULT_SMPD_DOWNSCALE, &
        &visibility=UI_VIS_STANDARD, preserve_default=.true.)
        call master%add_input(UI_PARM, 'total_dose',     'float',         'Total exposure dose (e/A2)', 'Total exposure dose (e/A2)', '', .true., '', visibility=UI_VIS_STANDARD)
        call master%add_input(UI_PARM, 'pickrefs',       'file',          '2D averages for use as picking references (optional)', '2D averages for use as picking references (optional)', '', .false., '', visibility=UI_VIS_STANDARD)
        call master%add_input(UI_PARM, 'box_extract',    'int',           'Force box size (px, optional)', 'force a box size (px) eg. to match an existing dataset"', '', .false., '', visibility=UI_VIS_STANDARD)
        call master%add_input(UI_FILE, 'dir_preprocess', 'dir',           'Pre-existing preprocessing directory', 'Pre-existing preprocessing directory', '', .false., '', visibility=UI_VIS_DEVELOPER)
        call master%add_input(UI_PARM, 'nicedispid',     'int',           'Optics group offset delta multiplier', 'Optics group offset delta multiplier', '0', .false., '', visibility=UI_VIS_DEVELOPER)
        call master%add_input(UI_PARM, 'thres',          'float',         'Distance threshold for peak picking(A)', 'Distance threshold for peak picking(A)', '0', .false., '', visibility=UI_VIS_DEVELOPER)
        call master%add_input(UI_PARM, 'nmics',          'int',           'Number of micrographs', 'Number of micrographs to collect before termination', '0', .false., '', visibility=UI_VIS_DEVELOPER)
        ! optics assignment, forwarded by the master
        call master%add_input(UI_PARM, 'beamtilt',       'binary',        'Use beam-tilt groups', &
        &'Split the micrographs by beam-tilt group before their beam-image shifts are clustered into optics groups(yes|no){no}', &
        &'', .false., 'no', choices=ui_choices([character(len=3) :: 'yes', 'no']), visibility=UI_VIS_DEVELOPER)
        call master%add_input(UI_PARM, 'tilt_thres',     'float',         'Beam-shift clustering threshold', &
        &'Distance threshold of the hierarchical clustering of the beam-image shifts into optics groups', '0.05', .false., '', &
        &visibility=UI_VIS_DEVELOPER)
        ! <no additional inputs>
        ! the 3D route of the initial analysis (picking references), forwarded by the master; 0 keeps the
        ! initial analysis' default. preserve_default: the defaults of nrestarts_collapse, lpstart_ini3D and
        ! lpstop_ini3D live in the initial analysis stage, which default_audit.py does not trace
        ! search controls
        call master%add_input(UI_SRCH, 'nstates_pickrefs', 'int',         'Picking-reference 3D states', &
        &'Number of states of the solve3D_cavgs run in the 3D route of the initial analysis (picking references); &
        &0 uses its default of 3', '0', .false., '', visibility=UI_VIS_DEVELOPER)
        call master%add_input(UI_SRCH, 'nstages_pickrefs', 'int',         'Picking-reference 3D stages', &
        &'Number of solve3D_cavgs stages in the 3D route of the initial analysis (picking references); &
        &0 uses its default of 3', '0', .false., '', visibility=UI_VIS_DEVELOPER)
        call master%add_input(UI_SRCH, 'nspace_pickrefs', 'int',          'Picking-reference reprojections', &
        &'Number of reprojections of the 3D route''s volume used as picking references by the initial analysis; &
        &0 uses its default of 50', '0', .false., '', visibility=UI_VIS_DEVELOPER)
        call master%add_input(UI_SRCH, 'nrestarts_collapse', 'int',       'Picking-reference 3D restarts', &
        &'Number of solve3D_cavgs restarts when states collapse, in the 3D route of the initial analysis &
        &(picking references){3}', '3', .false., 3., visibility=UI_VIS_DEVELOPER, preserve_default=.true.)
        ! filter controls
        call master%add_input(UI_FILT, 'lpstop_pickrefs', 'float',        'Picking-reference final low-pass limit', &
        &'Final low-pass limit (in Angstroms) of the solve3D_cavgs run in the 3D route of the initial analysis &
        &(picking references); 0 uses its default of 8', '0', .false., '', visibility=UI_VIS_DEVELOPER)
        call master%add_input(UI_FILT, 'lpstart_ini3D', 'float',          'Picking-reference initial 3D starting low-pass limit', &
        &'Starting low-pass limit (in Angstroms) of the solve3D_cavgs initial 3D model in the 3D route of the &
        &initial analysis (picking references){100}', '100', .false., 100., visibility=UI_VIS_DEVELOPER, preserve_default=.true.)
        call master%add_input(UI_FILT, 'lpstop_ini3D', 'float',           'Picking-reference initial 3D final low-pass limit', &
        &'Final low-pass limit (in Angstroms) of the solve3D_cavgs initial 3D model in the 3D route of the &
        &initial analysis (picking references){20}', '20', .false., 20., visibility=UI_VIS_DEVELOPER, preserve_default=.true.)
        ! mask controls
        ! computer controls
        call master%add_input(UI_COMP, 'nthr3D_pickrefs', 'int',          'Picking-reference 3D threads', &
        &'Number of OpenMP threads of the solve3D_cavgs run and the reprojection in the 3D route of the initial analysis &
        &(picking references); 0 uses the master''s resources table (16, or SIMPLE_STREAM_REFGEN_NTHR)', '0', .false., '', &
        &visibility=UI_VIS_DEVELOPER)
        ! add to ui_hash
        call add_ui_program('master', master, prgtab, UI_CATEGORY)
    end subroutine new_master

    subroutine new_pick_extract( prgtab )
        class(ui_hash), intent(inout) :: prgtab
        ! PROGRAM SPECIFICATION
        call pick_extract%new(&
        &'pick_extract', &                                                               ! name
        &'Pick particles and extract images as new microscope data arrive',&             ! summary
        &'is a distributed workflow that executes picking and extraction'//&             ! help
        &' in streaming mode as the microscope collects the data',&
        &'simple_stream',&                                                               ! executable
        &.true.,&                                                                        ! requires sp_project
        &visibility=UI_VIS_STANDARD, display_name='Pick and Extract During Acquisition')
        ! image input/output
        call pick_extract%add_input(UI_IMG, pickrefs, group="picking", visibility=UI_VIS_STANDARD)
        call pick_extract%add_input(UI_FILE, 'dir_exec', 'file', 'Previous run directory',&
        &'Directory where a previous pick_extract application was run', 'e.g. 2_pick_extract', .false., '', group="data", &
        &visibility=UI_VIS_ADVANCED)
        ! parameter input/output
        call pick_extract%add_input(UI_PARM, pcontrast,   group="picking", &
        &visibility=UI_VIS_ADVANCED)
        call pick_extract%add_input(UI_PARM, box_extract, group="extract", &
        &visibility=UI_VIS_ADVANCED)
        call pick_extract%add_input(UI_FILE, 'dir_target', 'file', 'Target directory',&
        &'Directory where the preprocess_stream application is running', 'e.g. 1_preproc', .true., '', group="data", &
        &visibility=UI_VIS_STANDARD)
        call pick_extract%add_input(UI_FILE, 'optics_dir', 'dir', 'Optics assignment directory',&
        &'Directory where the assign_optics application publishes its optics maps; the optics groups of the newest &
        &map are applied to the picked micrographs', 'e.g. optics_assignment', .false., '', group="data", &
        &visibility=UI_VIS_DEVELOPER)
        call pick_extract%add_input(UI_PARM, backgr_subtr, group="picking", &
        &visibility=UI_VIS_DEVELOPER)
        ! <no additional inputs>
        ! <empty>
        ! search controls
        call pick_extract%add_input(UI_SRCH, pick_roi, group="picking", &
        &visibility=UI_VIS_DEVELOPER)
        call pick_extract%add_input(UI_SRCH, 'thres', 'num', 'Peak-picking distance threshold', &
        &'Distance threshold in Angstroms for peak picking; 0 uses the picker default', 'distance threshold{0}', .false., 0., group="picking", &
        &visibility=UI_VIS_DEVELOPER)
        ! filter controls
        call pick_extract%add_input(UI_FILT, lp_pick,          group="picking", &
        &visibility=UI_VIS_ADVANCED)
        ! mask controls
        ! <empty>
        ! computer controls
        call pick_extract%add_input(UI_COMP, nthr,   group="compute", visibility=UI_VIS_STANDARD)
        call pick_extract%add_input(UI_COMP, nparts, group="compute", visibility=UI_VIS_STANDARD)
        call pick_extract%add_input(UI_COMP, 'walltime', 'num', 'Walltime', 'Maximum execution time for job scheduling and management in seconds{1740}(29mins)',&
        &'in seconds(29mins){1740}', .false., 1740., group="compute", &
        &visibility=UI_VIS_ADVANCED)
        ! add to ui_hash
        call add_ui_program('pick_extract', pick_extract, prgtab, UI_CATEGORY)
    end subroutine new_pick_extract

    subroutine new_preproc( prgtab )
        class(ui_hash), intent(inout) :: prgtab
        ! PROGRAM SPECIFICATION
        call preproc%new(&
        &'preproc', &                                                                       ! name
        &'Run streaming preprocessing as new data arrive',& ! summary
        &'is a distributed workflow that executes motion_correct, ctf_estimate and pick'//& ! help
        &' in sequence',&
        &'simple_stream',&                                                                    ! executable
        &.true., display_name='Run streaming preprocessing as new data arrive') ! requires sp_project
        ! INPUT PARAMETER SPECIFICATIONS
        ! image input/output
        call preproc%add_input(UI_FILE, dir_movies, group="data", visibility=UI_VIS_STANDARD)
        call preproc%add_input(UI_FILE, gainref,    group="data", visibility=UI_VIS_STANDARD)
        call preproc%add_input(UI_FILE, 'dir_meta', 'dir', 'Directory containing per-movie metadata in XML format',&
            &'Directory containing per-movie metadata XML files from EPU', 'e.g. /dataset/metadata', .false., '', group="data", visibility=UI_VIS_STANDARD)
        ! parameter input/output
        call preproc%add_input(UI_PARM, total_dose,                      group="data",              visibility=UI_VIS_STANDARD)
        call preproc%add_input(UI_PARM, fraction_dose_target,            group="data",              visibility=UI_VIS_STANDARD)
        call preproc%add_input(UI_PARM, 'nmics', 'num', 'Micrographs before termination', &
        &'Number of micrographs after which the preprocessing terminates; 0 = no limit{0}', '# micrographs{0}', .false., 0., &
        &group="data", visibility=UI_VIS_DEVELOPER)
        call preproc%add_input(UI_PARM, 'smpd_downscale', 'num', 'Sampling distance after downscale', &
        &'Distance between neighbouring pixels in Angstroms after downscale', 'pixel size in Angstroms', &
        &.false., STREAM_DEFAULT_SMPD_DOWNSCALE, group="motion correction", visibility=UI_VIS_STANDARD, &
        &preserve_default=.true.)
        call preproc%add_input(UI_PARM, eer_fraction,                    group="motion correction", &
        &visibility=UI_VIS_DEVELOPER)
        call preproc%add_input(UI_PARM, 'eer_upsampling', 'num', 'EER up-sampling factor', &
        &'Up-sampling factor of EER movies (1 or 2): 1 renders 4K and 2 renders 8K frames{1}', 'up-sampling factor{1}', &
        &.false., 1., group="motion correction", visibility=UI_VIS_DEVELOPER)
        call preproc%add_input(UI_PARM, max_dose,                        group="motion correction", &
        &visibility=UI_VIS_DEVELOPER)
        call preproc%add_input(UI_PARM, kv,    required_override=.true., group="data",              visibility=UI_VIS_STANDARD)
        call preproc%add_input(UI_PARM, cs,    required_override=.true., group="data",              visibility=UI_VIS_STANDARD)
        call preproc%add_input(UI_PARM, fraca, required_override=.true., group="data",              visibility=UI_VIS_STANDARD)
        call preproc%add_input(UI_PARM, smpd,  required_override=.true., group="data",              visibility=UI_VIS_STANDARD)
        call preproc%add_input(UI_PARM, fit_phshift, group="CTF estimation", visibility=UI_VIS_STANDARD)
        call preproc%add_input(UI_PARM, pspecsz, group="CTF estimation", &
        &visibility=UI_VIS_DEVELOPER)
        call preproc%add_input(UI_PARM, ctfpatch, group="CTF estimation", &
        &visibility=UI_VIS_DEVELOPER)
        call preproc%add_input(UI_PARM, flipgain, label_override='Gain reference processing', &
        &help_override='Flip the gain reference along the given axis, detect the flip from the movies (flip_auto), or &
        &generate a gain reference from the movies (generate); none and flip_x|flip_y|flip_xy are accepted as aliases of &
        &no and x|y|xy(no|x|y|xy|yx|flip_auto|generate){no}', &
        &choices_override=ui_choices([character(len=9) :: 'no', 'x', 'y', 'xy', 'yx', 'flip_auto', 'generate']), &
        &group="motion correction", visibility=UI_VIS_DEVELOPER)
        call preproc%add_input(UI_PARM, algorithm, group="motion correction", &
        &visibility=UI_VIS_DEVELOPER)
        call preproc%add_input(UI_PARM, mcconvention, group="motion correction", &
        &visibility=UI_VIS_DEVELOPER)
        call preproc%add_input(UI_PARM, mcpatch, group="motion correction", &
        &visibility=UI_VIS_DEVELOPER)
        call preproc%add_input(UI_PARM, mcpatch_thres, group="motion correction", &
        &visibility=UI_VIS_DEVELOPER)
        call preproc%add_input(UI_PARM, 'ninipick', 'num', 'Number of micrographs to perform initial picking preprocessing on',&
        & 'Number of micrographs to perform initial picking preprocessing on', 'e.g 500', .false., 0.0, &
        &visibility=UI_VIS_DEVELOPER)
        ! <no additional inputs>
        ! <empty>
        ! search controls
        call preproc%add_input(UI_SRCH, trs_mc, group="motion correction", &
        &visibility=UI_VIS_ADVANCED)
        call preproc%add_input(UI_SRCH, dfmin, group="CTF estimation", &
        &visibility=UI_VIS_DEVELOPER)
        call preproc%add_input(UI_SRCH, dfmax, group="CTF estimation", &
        &visibility=UI_VIS_DEVELOPER)
        call preproc%add_input(UI_SRCH, phshift_min, group="CTF estimation", &
        &visibility=UI_VIS_DEVELOPER)
        call preproc%add_input(UI_SRCH, phshift_max, group="CTF estimation", &
        &visibility=UI_VIS_DEVELOPER)
        call preproc%add_input(UI_SRCH, phshift_step, group="CTF estimation", &
        &visibility=UI_VIS_DEVELOPER)
        ! filter controls
        call preproc%add_input(UI_FILT, 'lpstart', 'num', 'Motion-correction low-pass start', &
        &'Starting low-pass limit for motion correction', 'low-pass limit in Angstroms{8}', .false., 8., group="motion correction", &
        &visibility=UI_VIS_DEVELOPER)
        call preproc%add_input(UI_FILT, 'lpstop', 'num', 'Motion-correction low-pass stop', &
        &'Final low-pass limit for motion correction', 'low-pass limit in Angstroms{5}', .false., 5., group="motion correction", &
        &visibility=UI_VIS_DEVELOPER)
        call preproc%add_input(UI_FILT, bfac, group="motion correction", &
        &visibility=UI_VIS_DEVELOPER)
        call preproc%add_input(UI_FILT, 'hp_ctf_estimate', 'num', 'CTF estimation high-pass limit', &
        &'High-pass limit for CTF parameter estimation', 'high-pass limit in Angstroms{30}', .false., 30., group="CTF estimation", &
        &visibility=UI_VIS_DEVELOPER)
        call preproc%add_input(UI_FILT, 'lp_ctf_estimate', 'num', 'CTF estimation low-pass limit', &
        &'Low-pass limit for CTF parameter estimation', 'low-pass limit in Angstroms{5}', .false., 5., group="CTF estimation", &
        &visibility=UI_VIS_DEVELOPER)
        call preproc%add_input(UI_FILT, ctfresthreshold, group="CTF estimation", &
        &visibility=UI_VIS_DEVELOPER)
        call preproc%add_input(UI_FILT, icefracthreshold, group="CTF estimation", &
        &visibility=UI_VIS_DEVELOPER)
        call preproc%add_input(UI_FILT, astigthreshold, group="CTF estimation", &
        &visibility=UI_VIS_DEVELOPER)
        ! mask controls
        ! <empty>
        ! computer controls
        call preproc%add_input(UI_COMP, nparts, group="compute", visibility=UI_VIS_STANDARD)
        call preproc%add_input(UI_COMP, nthr,   group="compute", visibility=UI_VIS_STANDARD)
        call preproc%add_input(UI_COMP, 'walltime', 'num', 'Walltime', 'Maximum execution time for job scheduling and management in seconds{1740}(29mins)',&
        &'in seconds(29mins){1740}', .false., 1740., group="compute", &
        &visibility=UI_VIS_DEVELOPER)
        ! add to ui_hash
        call add_ui_program('preproc', preproc, prgtab, UI_CATEGORY)
    end subroutine new_preproc

    subroutine new_sieve_cavgs( prgtab )
        class(ui_hash), intent(inout) :: prgtab
        ! PROGRAM SPECIFICATION
        call sieve_cavgs%new(&
        &'sieve_cavgs', &                                                       ! name
        &'Run 2D particle analysis automatically as new data arrive',& ! summary
        &'is a distributed workflow that executes 2D analysis'//&               ! help
        &' in streaming mode as the microscope collects the data',&
        &'simple_stream',&                                                      ! executable
        &.true.,&                                                               ! requires sp_project
        &visibility=UI_VIS_STANDARD, display_name='Analyze Streaming 2D Data')
        ! image input/output
        call sieve_cavgs%add_input(UI_IMG, 'refs', 'file', 'Compatibility-model references', &
        &'Class averages that pre-train the coarse and fine size-compatibility models before sieving starts; &
        &skipped when the file does not exist', 'e.g. references.mrc', .false., '', group="data", &
        &visibility=UI_VIS_DEVELOPER)
        ! parameter input/output
        call sieve_cavgs%add_input(UI_FILE, 'dir_target', 'file', 'Target directory',&
        &'Directory where the pick_extract application is running', 'e.g. 2_pick_extract', .true., '', group="data", visibility=UI_VIS_STANDARD)
        ! <no additional inputs>
        ! <empty>
        ! search controls
        call sieve_cavgs%add_input(UI_SRCH, 'nptcls_coarse', 'num', 'Target coarse-pass particle count', &
        &'Target number of particles in each coarse sieving chunk', '# particles{5000}', .false., 5000., group="cluster 2D", &
        &visibility=UI_VIS_DEVELOPER)
        call sieve_cavgs%add_input(UI_SRCH, 'nptcls_fine', 'num', 'Target fine-pass particle count', &
        &'Target number of particles in each fine sieving chunk', '# particles{10000}', .false., 10000., group="cluster 2D", &
        &visibility=UI_VIS_DEVELOPER)
        call sieve_cavgs%add_input(UI_SRCH, 'box_coarse', 'num', 'Coarse-pass box size', &
        &'Box size used during coarse streaming sieving', '# pixels{128}', .false., 128., group="cluster 2D", &
        &visibility=UI_VIS_DEVELOPER)
        call sieve_cavgs%add_input(UI_SRCH, 'box_fine', 'num', 'Fine-pass box size', &
        &'Box size used during fine streaming sieving', '# pixels{128}', .false., 128., group="cluster 2D", &
        &visibility=UI_VIS_DEVELOPER)
        call sieve_cavgs%add_input(UI_SRCH, 'nsample_coarse', 'num', 'Coarse-pass sample count', &
        &'Number of particles sampled during coarse streaming sieving', '# particles{2000}', .false., 2000., group="cluster 2D", &
        &visibility=UI_VIS_DEVELOPER)
        call sieve_cavgs%add_input(UI_SRCH, 'nsample_fine', 'num', 'Fine-pass sample count', &
        &'Number of particles sampled during fine streaming sieving', '# particles{2000}', .false., 2000., group="cluster 2D", &
        &visibility=UI_VIS_DEVELOPER)
        call sieve_cavgs%add_input(UI_SRCH, 'ncls_coarse', 'num', 'Coarse-pass class count', &
        &'Number of 2D classes used during coarse streaming sieving', '# classes{100}', .false., 100., group="cluster 2D", &
        &visibility=UI_VIS_DEVELOPER)
        call sieve_cavgs%add_input(UI_SRCH, 'ncls_fine', 'num', 'Fine-pass class count', &
        &'Number of 2D classes used during fine streaming sieving', '# classes{100}', .false., 100., group="cluster 2D", &
        &visibility=UI_VIS_DEVELOPER)
        call sieve_cavgs%add_input(UI_SRCH, 'use_model', 'binary', 'Use class-average rejection model', &
        &'Use the class-average rejection model during streaming sieving(yes|no){yes}', '', .false., 'yes', group="cluster 2D", &
        &choices=ui_choices([character(len=3) :: 'yes', 'no']), &
        &visibility=UI_VIS_DEVELOPER)
        call sieve_cavgs%add_input(UI_SRCH, 'single_pass', 'binary', 'Run coarse sieving pass only', &
        &'Run only the coarse pass of streaming sieving(yes|no){no}', '', .false., 'no', group="cluster 2D", &
        &choices=ui_choices([character(len=3) :: 'yes', 'no']), &
        &visibility=UI_VIS_DEVELOPER)
        ! filter controls
        call sieve_cavgs%add_input(UI_FILT, 'lpstart', 'num', 'Initial sieving low-pass limit', &
        &'Low-pass limit used to initialize streaming sieving', 'low-pass limit in Angstroms{15}', .false., 15., group="cluster 2D", &
        &visibility=UI_VIS_DEVELOPER)
        call sieve_cavgs%add_input(UI_FILT, 'lpstop_coarse', 'num', 'Coarse-pass low-pass limit', &
        &'Final low-pass limit for coarse streaming sieving', 'low-pass limit in Angstroms{15}', .false., 15., group="cluster 2D", &
        &visibility=UI_VIS_DEVELOPER)
        call sieve_cavgs%add_input(UI_FILT, 'lpstop_fine', 'num', 'Fine-pass low-pass limit', &
        &'Final low-pass limit for fine streaming sieving', 'low-pass limit in Angstroms{10}', .false., 10., group="cluster 2D", &
        &visibility=UI_VIS_DEVELOPER)
        ! mask controls
        ! <empty>: the mask diameter is the one of the picking references (moldiam.txt in dir_target)
        ! computer controls
        call sieve_cavgs%add_input(UI_COMP, nchunks,                          group="compute", visibility=UI_VIS_STANDARD)
        call sieve_cavgs%add_input(UI_COMP, nparts, required_override=.true., group="compute", visibility=UI_VIS_STANDARD)
        call sieve_cavgs%add_input(UI_COMP, nthr, group="compute", visibility=UI_VIS_STANDARD)
        call sieve_cavgs%add_input(UI_COMP, 'walltime', 'num', 'Walltime', 'Maximum execution time for job scheduling and management in seconds{1740}(29mins)',&
        &'in seconds(29mins){1740}', .false., 1740., group="compute", &
        &visibility=UI_VIS_ADVANCED)
        ! add to ui_hash
        call add_ui_program('sieve_cavgs', sieve_cavgs, prgtab, UI_CATEGORY)
    end subroutine new_sieve_cavgs

end module simple_ui_stream
