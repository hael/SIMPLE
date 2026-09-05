!@descr: module defining the user interfaces for streaming workflows in the simple_exec suite
module simple_ui_stream
use simple_ui_modules
implicit none

type(category_descriptor), parameter :: UI_CATEGORY = category_descriptor('stream', 'Stream Workflows', 10)
type(ui_program), target :: abinitio2D_stream
type(ui_program), target :: assign_optics
type(ui_program), target :: gen_pickrefs
type(ui_program), target :: master
type(ui_program), target :: pick_extract
type(ui_program), target :: preproc
type(ui_program), target :: sieve_cavgs

contains

    subroutine construct_stream_programs(prgtab)
        class(ui_hash), intent(inout) :: prgtab
        call new_abinitio2D_stream(prgtab)
        call new_assign_optics(prgtab)
        call new_gen_pickrefs(prgtab)
        call new_master(prgtab)
        call new_pick_extract(prgtab)
        call new_preproc(prgtab)
        call new_sieve_cavgs(prgtab)
    end subroutine construct_stream_programs

subroutine new_abinitio2D_stream( prgtab )
        class(ui_hash), intent(inout) :: prgtab
        ! PROGRAM SPECIFICATION
        call abinitio2D_stream%new(&
        &'abinitio2D_stream', &                                                  ! name
        &'Run streaming 2D analysis as new data arrive',& ! summary
        &'is a distributed workflow that executes 2D analysis'//&                ! help
        &' in streaming mode as the microscope collects the data',&
        &'simple_stream',&                                                       ! executable
        &.true.,&                                                                ! requires sp_project
        &visibility=UI_VIS_DEVELOPER)
        ! image input/output
        ! <empty>
        ! parameter input/output
        call abinitio2D_stream%add_input(UI_PARM, sigma_store, group="cluster 2D", visibility=UI_VIS_ADVANCED)
        call abinitio2D_stream%add_input(UI_FILE, 'dir_target', 'file', 'Target directory',&
        &'Directory where the pick_extract application is running', 'e.g. 2_pick_extract', .true., '', group="data", visibility=UI_VIS_STANDARD)
        call abinitio2D_stream%add_input(UI_FILE, 'dir_exec', 'file', 'Previous run directory',&
        &'Directory where previous 2D analysis took place', 'e.g. 3_abinitio2D_stream', .false., '', group="data", &
        &visibility=UI_VIS_DEVELOPER)
        ! <no additional inputs>
        ! <empty>
        ! search controls
        call abinitio2D_stream%add_input(UI_SRCH, 'ncls', 'num', 'Maximum number of 2D clusters',&
        &'Maximum number of 2D class averages for the pooled particles subsets', 'Maximum # 2D clusters', .true., 200., group="cluster 2D",&
        &visibility=UI_VIS_STANDARD)
        call abinitio2D_stream%add_input(UI_SRCH, 'center', 'binary', 'Center class averages', &
        &'Center class averages by their center of gravity and map shifts back to the particles(yes|no){yes}', '', .false., 'yes', group="cluster 2D", &
        &choices=ui_choices([character(len=3) :: 'yes', 'no']), &
        &visibility=UI_VIS_ADVANCED)
        ! filter controls
        ! <empty>
        ! mask controls
        call abinitio2D_stream%add_input(UI_MASK, 'mskdiam', 'num', 'Mask diameter', 'Mask diameter (in A) for application of a soft-edged circular mask to &
        &remove background noise', 'mask diameter in A', .false., 0., group="cluster 2D", visibility=UI_VIS_STANDARD)
        ! computer controls
        call abinitio2D_stream%add_input(UI_COMP, nparts, group="compute", visibility=UI_VIS_STANDARD)
        call abinitio2D_stream%add_input(UI_COMP, nthr, group="compute", visibility=UI_VIS_STANDARD)
        call abinitio2D_stream%add_input(UI_COMP, 'walltime', 'num', 'Walltime', 'Maximum execution time for job scheduling and management in seconds{1740}(29mins)',&
        &'in seconds(29mins){1740}', .false., 1740., group="compute", &
        &visibility=UI_VIS_DEVELOPER)
        ! add to ui_hash
        call add_ui_program('abinitio2D_stream', abinitio2D_stream, prgtab, UI_CATEGORY)
    end subroutine new_abinitio2D_stream

    subroutine new_assign_optics( prgtab )
        class(ui_hash), intent(inout) :: prgtab
        ! PROGRAM SPECIFICATION
        call assign_optics%new(&
        &'assign_optics', &                                              ! name
        &'Assign optics groups from microscope metadata',& ! summary
        &'is a program to assign optics groups during streaming',&       ! descr long
        &'simple_stream',&                                               ! executable
        &.true., &
        &visibility=UI_VIS_DEVELOPER)                                                         ! requires sp_project
        ! INPUT PARAMETER SPECIFICATIONS
        ! image input/output
        ! parameter input/output
        call assign_optics%add_input(UI_FILE, 'dir_target', 'file', 'Target directory',&
        &'Directory where the preprocess_stream application is running', 'e.g. 1_preproc', .true., '', &
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
        &visibility=UI_VIS_ADVANCED)                                                         ! requires sp_project
        ! INPUT PARAMETER SPECIFICATIONS
        ! image input/output
        ! parameter input/output
        call gen_pickrefs%add_input(UI_FILE, 'dir_target', 'file', 'Target directory',&
        &'Directory where the preprocess_stream application is running', 'e.g. 1_preproc', .true., '', &
        &visibility=UI_VIS_STANDARD)
        call gen_pickrefs%add_input(UI_PARM, 'nmics', 'num', 'Number of micrographs to import',&
        &'Number of micrographs to import for opening 2D', 'Number micrographs', .false., 100., &
        &visibility=UI_VIS_ADVANCED)
        call gen_pickrefs%add_input(UI_SRCH, nptcls_per_cls, group="cluster 2D", &
        &visibility=UI_VIS_DEVELOPER)
        call gen_pickrefs%add_input(UI_FILE, 'optics_dir', 'dir', 'Target directory for optics import',&
        &'Directory where assign_optics application is running', 'e.g. optics_assignment', .false., '', &
        &visibility=UI_VIS_ADVANCED)
        ! <no additional inputs>
        ! <empty>
        ! search controls
        call gen_pickrefs%add_input(UI_SRCH, pick_roi, group="picking", &
        &visibility=UI_VIS_DEVELOPER)
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
        &.false.)                                                                     ! requires sp_project
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
        call master%add_input(UI_PARM, sigma_store, group="cluster 2D", visibility=UI_VIS_ADVANCED)
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
        &'Downscaled pixel size (A)', '', .false., STREAM_DEFAULT_SMPD_DOWNSCALE, &
        &visibility=UI_VIS_DEVELOPER, preserve_default=.true.)
        call master%add_input(UI_PARM, 'total_dose',     'float',         'Total exposure dose (e/A2)', 'Total exposure dose (e/A2)', '', .true., '', visibility=UI_VIS_STANDARD)
        call master%add_input(UI_FILE, 'pickrefs',       'file',          '2D averages for use as picking references (optional)', '2D averages for use as picking references (optional)', '', .false., '', visibility=UI_VIS_STANDARD)
        call master%add_input(UI_PARM, 'box_extract',    'int',           'Force box size (px, optional)', 'force a box size (px) eg. to match an existing dataset"', '', .false., '', visibility=UI_VIS_STANDARD)
        call master%add_input(UI_FILE, 'dir_preprocess', 'dir',           'Pre-existing preprocessing directory', 'Pre-existing preprocessing directory', '', .false., '', visibility=UI_VIS_DEVELOPER)
        call master%add_input(UI_PARM, 'nicedispid',     'int',           'Optics group offset delta multiplier', 'Optics group offset delta multiplier', '0', .false., '', visibility=UI_VIS_DEVELOPER)
        call master%add_input(UI_PARM, 'thres',          'float',         'Distance threshold for peak picking(A)', 'Distance threshold for peak picking(A)', '0', .false., '', visibility=UI_VIS_DEVELOPER)
        call master%add_input(UI_PARM, 'nmics',          'int',           'Number of micrographs', 'Number of micrographs to collect before termination', '0', .false., '', visibility=UI_VIS_DEVELOPER)
        ! <no additional inputs>
        ! search controls
        ! filter controls
        ! mask controls
        ! computer controls
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
        call pick_extract%add_input(UI_PARM, moldiam,     group="picking", &
        &visibility=UI_VIS_ADVANCED)
        call pick_extract%add_input(UI_FILE, 'dir_target', 'file', 'Target directory',&
        &'Directory where the preprocess_stream application is running', 'e.g. 1_preproc', .true., '', group="data", &
        &visibility=UI_VIS_STANDARD)
        call pick_extract%add_input(UI_PARM, 'nmoldiams', 'num', 'Number of molecular diameters to investigate', 'Number of molecular diameters tested',&
        &'e.g. 5', .false., 5., group="picking", &
        &visibility=UI_VIS_ADVANCED)
        call pick_extract%add_input(UI_PARM, moldiam_max, group="picking", &
        &visibility=UI_VIS_ADVANCED)
        call pick_extract%add_input(UI_PARM, backgr_subtr, group="picking", &
        &visibility=UI_VIS_DEVELOPER)
        ! <no additional inputs>
        ! <empty>
        ! search controls
        call pick_extract%add_input(UI_SRCH, pgrp, required_override=.false., group="picking", visibility=UI_VIS_STANDARD)
        call pick_extract%add_input(UI_SRCH, pick_roi, group="picking", &
        &visibility=UI_VIS_DEVELOPER)
        call pick_extract%add_input(UI_SRCH, 'thres', 'num', 'Peak-picking distance threshold', &
        &'Distance threshold in Angstroms for peak picking; 0 uses the picker default', 'distance threshold{0}', .false., 0., group="picking", &
        &visibility=UI_VIS_DEVELOPER)
        ! filter controls
        call pick_extract%add_input(UI_FILT, lp_pick,          group="picking", &
        &visibility=UI_VIS_ADVANCED)
        call pick_extract%add_input(UI_FILT, ctfresthreshold,  group="data", &
        &visibility=UI_VIS_ADVANCED)
        call pick_extract%add_input(UI_FILT, icefracthreshold, group="data", &
        &visibility=UI_VIS_ADVANCED)
        call pick_extract%add_input(UI_FILT, astigthreshold,   group="data", &
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
        &.true.)                                                                            ! requires sp_project
        ! INPUT PARAMETER SPECIFICATIONS
        ! image input/output
        call preproc%add_input(UI_FILE, dir_movies, group="data", visibility=UI_VIS_STANDARD)
        call preproc%add_input(UI_FILE, gainref,    group="data", visibility=UI_VIS_STANDARD)
        call preproc%add_input(UI_FILE, 'dir_prev', 'file', 'Previous run directory',&
            &'Directory where a previous stream application was run', 'e.g. 2_preproc', .false., '', group="data", &
        &visibility=UI_VIS_DEVELOPER)
        call preproc%add_input(UI_FILE, 'dir_meta', 'dir', 'Directory containing per-movie metadata in XML format',&
            &'Directory containing per-movie metadata XML files from EPU', 'e.g. /dataset/metadata', .false., '', group="data", visibility=UI_VIS_STANDARD)
        ! parameter input/output
        call preproc%add_input(UI_PARM, total_dose,                      group="data",              visibility=UI_VIS_STANDARD)
        call preproc%add_input(UI_PARM, fraction_dose_target,            group="data",              visibility=UI_VIS_STANDARD)
        call preproc%add_input(UI_PARM, 'smpd_downscale', 'num', 'Sampling distance after downscale', &
        &'Distance between neighbouring pixels in Angstroms after downscale', 'pixel size in Angstroms', &
        &.false., STREAM_DEFAULT_SMPD_DOWNSCALE, group="motion correction", visibility=UI_VIS_STANDARD, &
        &preserve_default=.true.)
        call preproc%add_input(UI_PARM, eer_fraction,                    group="motion correction", &
        &visibility=UI_VIS_DEVELOPER)
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
        call preproc%add_input(UI_PARM, flipgain, group="motion correction", &
        &visibility=UI_VIS_DEVELOPER)
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
        call preproc%add_input(UI_SRCH, 'tilt_thres', 'num', 'Threshold for hierarchical clustering of beamtilts',&
        & 'Threshold for hierarchical clustering of beamtilts', 'e.g 0.05', .false., 0.05, group="optics groups", &
        &visibility=UI_VIS_DEVELOPER)
        call preproc%add_input(UI_SRCH, 'beamtilt', 'binary', 'Use beamtilts in optics group assignment',&
        & 'Use beamtilt values (if found in EPU filenames) during optics group assignment(yes|no){yes}','', .false., 'no', group="optics groups", &
        &choices=ui_choices([character(len=3) :: 'yes', 'no']), &
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
        ! <empty>
        ! parameter input/output
        call sieve_cavgs%add_input(UI_FILE, 'dir_target', 'file', 'Target directory',&
        &'Directory where the pick_extract application is running', 'e.g. 2_pick_extract', .true., '', group="data", visibility=UI_VIS_STANDARD)
        call sieve_cavgs%add_input(UI_FILE, 'dir_exec', 'file', 'Previous run directory',&
        &'Directory where previous 2D analysis took place', 'e.g. 3_sieve_cavgs', .false., '', group="data", &
        &visibility=UI_VIS_ADVANCED)
        call sieve_cavgs%add_input(UI_PARM, 'nmics', 'num', 'Micrographs per sieve cycle', &
        &'Number of micrographs imported before each streaming sieving cycle', '# micrographs{100}', .false., 100., group="data", &
        &visibility=UI_VIS_DEVELOPER)
        ! <no additional inputs>
        ! <empty>
        ! search controls
        call sieve_cavgs%add_input(UI_SRCH, ncls,                                     group="cluster 2D", visibility=UI_VIS_STANDARD)
        call sieve_cavgs%add_input(UI_PARM, sigma_store, group="cluster 2D", visibility=UI_VIS_ADVANCED)
        call sieve_cavgs%add_input(UI_SRCH, nptcls_per_cls, required_override=.true., group="cluster 2D", visibility=UI_VIS_STANDARD)
        call sieve_cavgs%add_input(UI_SRCH, nchunksperset,                                                      visibility=UI_VIS_STANDARD)
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
        call sieve_cavgs%add_input(UI_SRCH, maxnchunks, group="cluster 2D", &
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
        call sieve_cavgs%add_input(UI_MASK, 'mskdiam', 'num', 'Mask diameter', 'Mask diameter (in A) for application of a soft-edged circular mask to &
        &remove background noise', 'mask diameter in A', .false., 0., group="cluster 2D", visibility=UI_VIS_STANDARD)
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
