!@descr: module defining the user interfaces for pre-processing programs in the simple_exec suite
module simple_ui_preproc
use simple_ui_modules
implicit none

type(category_descriptor), parameter :: UI_CATEGORY = category_descriptor('preproc', 'Pre-processing', 20)
type(ui_program), target :: assign_optics_groups
type(ui_program), target :: ctf_estimate
type(ui_program), target :: extract
type(ui_program), target :: gen_pspecs_and_thumbs
type(ui_program), target :: motion_correct
type(ui_program), target :: particle_sieving
type(ui_program), target :: pick
type(ui_program), target :: preprocess
type(ui_program), target :: reextract

contains

    subroutine construct_preproc_programs(prgtab)
        class(ui_hash), intent(inout) :: prgtab
        call new_assign_optics_groups(prgtab)
        call new_ctf_estimate(prgtab)
        call new_extract(prgtab)
        call new_gen_pspecs_and_thumbs(prgtab)
        call new_motion_correct(prgtab)
        call new_particle_sieving(prgtab)
        call new_pick(prgtab)
        call new_preprocess(prgtab)
        call new_reextract(prgtab)
    end subroutine construct_preproc_programs

    subroutine print_preproc_programs(logfhandle)
        integer, intent(in) :: logfhandle
        write(logfhandle,'(A)') format_str('PRE-PROCESSING:', C_UNDERLINED)
        write(logfhandle,'(A)') assign_optics_groups%name%to_char()
        write(logfhandle,'(A)') ctf_estimate%name%to_char()
        write(logfhandle,'(A)') extract%name%to_char()
        write(logfhandle,'(A)') gen_pspecs_and_thumbs%name%to_char()
        write(logfhandle,'(A)') motion_correct%name%to_char()
        write(logfhandle,'(A)') particle_sieving%name%to_char()
        write(logfhandle,'(A)') pick%name%to_char()
        write(logfhandle,'(A)') preprocess%name%to_char()
        write(logfhandle,'(A)') reextract%name%to_char()
        write(logfhandle,'(A)') ''
    end subroutine print_preproc_programs

    subroutine new_assign_optics_groups( prgtab )
        class(ui_hash), intent(inout) :: prgtab
        ! PROGRAM SPECIFICATION
        call assign_optics_groups%new(&
        &'assign_optics_groups', &                                              ! name
        &'Assign optics groups from microscope metadata',& ! summary
        &'is a program to assign optics groups',&                               ! descr long
        &'simple_exec',&                                                        ! executable
        &.true.)                                                                ! requires sp_project
        ! INPUT PARAMETER SPECIFICATIONS
        ! image input/output
        ! parameter input/output
        call assign_optics_groups%add_input(UI_FILE, 'xmldir', 'dir', 'Directory containing per movie EPU XML files',&
        & 'Directory containing per movie EPY XML files', 'e.g. /data/datasetid/xml', .false., '', &
        &visibility=UI_VIS_DEVELOPER)
        call assign_optics_groups%add_input(UI_PARM, 'maxpop', 'num', 'Maximum number of movies/micrographs/stacks in each optics group',&
        & 'Maximum number of movies/micrographs/stacks in each optics group', 'e.g. 100', .false., '', &
        &visibility=UI_VIS_DEVELOPER)
        call assign_optics_groups%add_input(UI_PARM, 'optics_offset', 'num', 'Numbering offset to apply to optics groups',&
        & 'Numbering offset to apply to optics groups. Aids with combining datasets', 'e.g. 10', .false., '', &
        &visibility=UI_VIS_DEVELOPER)
        call assign_optics_groups%add_input(UI_PARM, 'tilt_thres', 'num', 'Threshold for hierarchical clustering of beamtilts',&
        & 'Threshold for hierarchical clustering of beamtilts', 'e.g 0.05', .false., 0.05, &
        &visibility=UI_VIS_DEVELOPER)
        call assign_optics_groups%add_input(UI_PARM, 'beamtilt', 'binary', 'Use beamtilts in optics group assignment',&
        &'Use beamtilt values (if found in EPU filenames) during optics group assignment(yes|no){yes}','', .false., 'yes', &
        &choices=ui_choices([character(len=3) :: 'yes', 'no']), &
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
        ! <empty>
        ! add to ui_hash
        call add_ui_program('assign_optics_groups', assign_optics_groups, prgtab, UI_CATEGORY)
    end subroutine new_assign_optics_groups

    subroutine new_ctf_estimate( prgtab )
        class(ui_hash), intent(inout) :: prgtab
        ! PROGRAM SPECIFICATION
        call ctf_estimate%new(&
        &'ctf_estimate', &                                                  ! name
        &'Estimate CTF parameters from micrographs',& ! summary
        &'is a distributed SIMPLE workflow for CTF parameter fitting',&     ! help
        &'simple_exec',&                                                    ! executable
        &.true.,&                                                           ! requires sp_project
        &visibility=UI_VIS_STANDARD)
        ! INPUT PARAMETER SPECIFICATIONS
        ! image input/output
        ! <empty>
        ! parameter input/output
        call ctf_estimate%add_input(UI_PARM, pspecsz, group="CTF estimation", &
        &visibility=UI_VIS_ADVANCED)
        call ctf_estimate%add_input(UI_PARM, ctfpatch, group="CTF estimation", &
        &visibility=UI_VIS_ADVANCED)
        call ctf_estimate%add_input(UI_PARM, fit_phshift, group="CTF estimation", &
        &visibility=UI_VIS_ADVANCED)
        ! <no additional inputs>
        ! <empty>
        ! search controls
        call ctf_estimate%add_input(UI_SRCH, dfmin, group="CTF estimation", &
        &visibility=UI_VIS_ADVANCED)
        call ctf_estimate%add_input(UI_SRCH, dfmax, group="CTF estimation", &
        &visibility=UI_VIS_ADVANCED)
        call ctf_estimate%add_input(UI_SRCH, astigtol, group="CTF estimation", &
        &visibility=UI_VIS_ADVANCED)
        call ctf_estimate%add_input(UI_SRCH, phshift_min, group="CTF estimation", &
        &visibility=UI_VIS_ADVANCED)
        call ctf_estimate%add_input(UI_SRCH, phshift_max, group="CTF estimation", &
        &visibility=UI_VIS_ADVANCED)
        call ctf_estimate%add_input(UI_SRCH, phshift_step, group="CTF estimation", &
        &visibility=UI_VIS_ADVANCED)
        ! filter controls
        call ctf_estimate%add_input(UI_FILT, lp, required_override=.false., group="CTF estimation", &
        &visibility=UI_VIS_ADVANCED)
        call ctf_estimate%add_input(UI_FILT, hp, required_override=.false., group="CTF estimation", &
        &visibility=UI_VIS_ADVANCED)
        ! mask controls
        ! <empty>
        ! computer controls
        call ctf_estimate%add_input(UI_COMP, nparts, group="compute", &
        &visibility=UI_VIS_STANDARD)
        call ctf_estimate%add_input(UI_COMP, nthr, group="compute", &
        &visibility=UI_VIS_STANDARD)
        ! add to ui_hash
        call add_ui_program('ctf_estimate', ctf_estimate, prgtab, UI_CATEGORY)
    end subroutine new_ctf_estimate

    subroutine new_extract( prgtab )
        class(ui_hash), intent(inout) :: prgtab
        ! PROGRAM SPECIFICATION
        call extract%new(&
        &'extract', &                                                           ! name
        &'Extract particle images from integrated movies',&                     ! summary
        &'is a program for extracting particle images from integrated movies',& ! descr long
        &'simple_exec',&                                                        ! executable
        &.true.,&                                                               ! requires sp_project
        &visibility=UI_VIS_STANDARD)
        ! INPUT PARAMETER SPECIFICATIONS
        ! image input/output
        call extract%add_input(UI_FILE, 'dir_box', 'dir', 'Box files directory', 'Directory to read the box files from', 'e.g. boxes/', .false., '',&
        &group="extract", &
        &visibility=UI_VIS_ADVANCED)
        ! parameter input/output
        call extract%add_input(UI_PARM, box, required_override=.false., group="extract", visibility=UI_VIS_STANDARD)
        call extract%add_input(UI_PARM, pcontrast, group="extract", visibility=UI_VIS_STANDARD)
        call extract%add_input(UI_PARM, outside, group="extract", &
        &visibility=UI_VIS_ADVANCED)
        call extract%add_input(UI_PARM, backgr_subtr, group="extract", &
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
        call extract%add_input(UI_COMP, nparts, group="compute", visibility=UI_VIS_STANDARD)
        call extract%add_input(UI_COMP, nthr, group="compute", visibility=UI_VIS_STANDARD)
        ! add to ui_hash
        call add_ui_program('extract', extract, prgtab, UI_CATEGORY)
    end subroutine new_extract

    subroutine new_gen_pspecs_and_thumbs( prgtab )
        class(ui_hash), intent(inout) :: prgtab
        ! PROGRAM SPECIFICATION
        call gen_pspecs_and_thumbs%new(&
        &'gen_pspecs_and_thumbs', &                                              ! name
        &'Correct anisotropic motion in movie frames',& ! summary
        &'is a distributed workflow for generating power spectra and thumbnails&
        & for imported integrated movies',&                                      ! help
        &'simple_exec',&                                                         ! executable
        &.true.)                                                                 ! requires sp_project
        ! INPUT PARAMETER SPECIFICATIONS
        ! image input/output
        ! <empty>
        ! parameter input/output
        call gen_pspecs_and_thumbs%add_input(UI_PARM, pspecsz, &
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
        call gen_pspecs_and_thumbs%add_input(UI_COMP, nparts, &
        &visibility=UI_VIS_STANDARD)
        call gen_pspecs_and_thumbs%add_input(UI_COMP, nthr, &
        &visibility=UI_VIS_STANDARD)
        ! add to ui_hash
        call add_ui_program('gen_pspecs_and_thumbs', gen_pspecs_and_thumbs, prgtab, UI_CATEGORY)
    end subroutine new_gen_pspecs_and_thumbs

    subroutine new_motion_correct( prgtab )
        class(ui_hash), intent(inout) :: prgtab
        ! PROGRAM SPECIFICATION
        call motion_correct%new(&
        &'motion_correct', &                                                            ! name
        &'Anisotropic motion correction of movies',&                                    ! summary
        &'is a distributed workflow for anisotropic motion correction of movies.&
        & If then total dose is given the individual frames will be filtered accordingly&
        & (dose-weighting strategy). If scale is given, the movie will be Fourier cropped according to&
        & the down-scaling factor (for super-resolution movies).',&                     ! help
        &'simple_exec',&                                                                ! executable
        &.true.,&                                                                       ! requires sp_project
        &visibility=UI_VIS_STANDARD)
        ! INPUT PARAMETER SPECIFICATIONS
        ! image input/output
        call motion_correct%add_input(UI_FILE, gainref, group="data", visibility=UI_VIS_STANDARD)
        ! parameter input/output
        call motion_correct%add_input(UI_PARM, total_dose, group="data", visibility=UI_VIS_STANDARD)
        call motion_correct%add_input(UI_PARM, fraction_dose_target, group="data", &
        &visibility=UI_VIS_ADVANCED)
        call motion_correct%add_input(UI_PARM, max_dose, group="data", &
        &visibility=UI_VIS_ADVANCED)
        call motion_correct%add_input(UI_PARM, smpd_downscale, group="data", &
        &visibility=UI_VIS_ADVANCED)
        call motion_correct%add_input(UI_PARM, 'fbody', 'string', 'Template output micrograph name',&
        &'Template output integrated movie name', 'e.g. mic_', .false., '', group="data", &
        &visibility=UI_VIS_ADVANCED)
        call motion_correct%add_input(UI_PARM, pspecsz, group="motion correction", &
        &visibility=UI_VIS_ADVANCED)
        call motion_correct%add_input(UI_PARM, eer_fraction, group="data", &
        &visibility=UI_VIS_ADVANCED)
        call motion_correct%add_input(UI_PARM, flipgain, group="motion correction", &
        &visibility=UI_VIS_ADVANCED)
        ! <no additional inputs>
        ! <empty>
        ! search controls
        call motion_correct%add_input(UI_SRCH, trs_mc, group ="motion correction", &
        &visibility=UI_VIS_ADVANCED)
        call motion_correct%add_input(UI_SRCH, 'bfac', 'num', 'B-factor applied to frames', 'B-factor applied to frames (in Angstroms^2)', &
        &'in Angstroms^2{50}', .false., 50., group="motion correction", &
        &visibility=UI_VIS_ADVANCED)
        call motion_correct%add_input(UI_SRCH, mcpatch, group="motion correction", &
        &visibility=UI_VIS_ADVANCED)
        call motion_correct%add_input(UI_SRCH, nxpatch, group="motion correction", &
        &visibility=UI_VIS_ADVANCED)
        call motion_correct%add_input(UI_SRCH, nypatch, group="motion correction", &
        &visibility=UI_VIS_ADVANCED)
        call motion_correct%add_input(UI_SRCH, mcconvention, group="motion correction", &
        &visibility=UI_VIS_ADVANCED)
        call motion_correct%add_input(UI_SRCH, algorithm, group="motion correction", &
        &visibility=UI_VIS_ADVANCED)
        call motion_correct%add_input(UI_SRCH, mcpatch_thres, group="motion correction", &
        &visibility=UI_VIS_ADVANCED)
        ! filter controls
        call motion_correct%add_input(UI_FILT, 'lpstart', 'num', 'Initial low-pass limit', 'Low-pass limit to be applied in the first &
        &iterations of movie alignment (in Angstroms){8}', 'in Angstroms{8}', .false., 8., group="motion correction", &
        &visibility=UI_VIS_ADVANCED)
        call motion_correct%add_input(UI_FILT, 'lpstop', 'num', 'Final low-pass limit', 'Low-pass limit to be applied in the last &
        &iterations of movie alignment (in Angstroms){5}', 'in Angstroms{5}', .false., 5., group="motion correction", &
        &visibility=UI_VIS_ADVANCED)
        call motion_correct%add_input(UI_FILT, wcrit, group="motion correction", &
        &visibility=UI_VIS_ADVANCED)
        ! mask controls
        ! <empty>
        ! computer controls
        call motion_correct%add_input(UI_COMP, nparts, group="compute", &
        &visibility=UI_VIS_STANDARD)
        call motion_correct%add_input(UI_COMP, nthr, group="compute", &
        &visibility=UI_VIS_STANDARD)
        ! add to ui_hash
        call add_ui_program('motion_correct', motion_correct, prgtab, UI_CATEGORY)
    end subroutine new_motion_correct

    subroutine new_particle_sieving( prgtab )
        class(ui_hash), intent(inout) :: prgtab
        ! PROGRAM SPECIFICATION
        call particle_sieving%new(&
        &'particle_sieving', &                                               ! name
        &'Select particles automatically using a sieving workflow',& ! summary
        &'workflow for automated particle sieving',& ! help
        &'simple_exec',&                                                   ! executable
        &.true.,&                                                          ! requires sp_project
        &visibility=UI_VIS_STANDARD)
        ! INPUT PARAMETER SPECIFICATIONS
        ! image input/output
        ! <empty>
        ! parameter input/output
        ! <empty>
        ! <no additional inputs>
        ! <empty>
        ! search controls
        call particle_sieving%add_input(UI_SRCH, 'nmics', 'num', 'Max # of micrographs per chunk', &
        &'Maximum number of micrographs accumulated into one chunk{50}', &
        &'max # of micrographs per chunk{50}', .false., 50., group="search", visibility=UI_VIS_ADVANCED)
        call particle_sieving%add_input(UI_SRCH, 'nptcls_coarse', 'num', 'Target # of particles per coarse chunk', 'Target # of particles per coarse chunk{5000}',&
        &'Target # of particles per coarse chunk{5000}', .false., 5000., group="search", visibility=UI_VIS_ADVANCED)
        call particle_sieving%add_input(UI_SRCH, 'nptcls_fine', 'num', 'Target # of particles per fine chunk', 'Target # of particles per fine chunk{10000}',&
        &'Target # of particles per fine chunk{10000}', .false., 10000., group="search", visibility=UI_VIS_ADVANCED)
        call particle_sieving%add_input(UI_SRCH, 'maxnchunks', 'num', 'Max # of chunks to process', &
        &'Cap on the total number of chunks processed, 0 = no limit{0}', &
        &'max # of chunks (0=no limit){0}', .false., 0., group="search", visibility=UI_VIS_ADVANCED)
        call particle_sieving%add_input(UI_SRCH, 'box_fine', 'num', 'Box size for fine sieving', 'Box size for fine sieving{128}', 'in pixels{128}', .false., 128., group="search", visibility=UI_VIS_ADVANCED)
        call particle_sieving%add_input(UI_SRCH, 'box_coarse', 'num', 'Box size for coarse sieving', 'Box size for coarse sieving{128}', 'in pixels{128}', .false., 128., group="search", visibility=UI_VIS_ADVANCED)
        call particle_sieving%add_input(UI_SRCH, 'ncls_coarse', 'num', 'Number of clusters for coarse sieving', 'Number of clusters for coarse sieving{100}', 'in clusters{100}', .false., 100., group="search", visibility=UI_VIS_ADVANCED)
        call particle_sieving%add_input(UI_SRCH, 'ncls_fine', 'num', 'Number of clusters for fine sieving', 'Number of clusters for fine sieving{100}', 'in clusters{100}', .false., 100., group="search", visibility=UI_VIS_ADVANCED)
        ! <no additional inputs>
        call particle_sieving%add_input(UI_SRCH, 'nsample_coarse', 'num', 'Number of particles to sample in coarse sieving', 'Number of particles to sample in coarse sieving{2000}', 'in particles{2000}', .false., 2000., group="search", visibility=UI_VIS_ADVANCED)
        call particle_sieving%add_input(UI_SRCH, 'nsample_fine', 'num', 'Number of particles to sample in fine sieving', 'Number of particles to sample in fine sieving{2000}', 'in particles{2000}', .false., 2000., group="search", visibility=UI_VIS_ADVANCED)
        call particle_sieving%add_input(UI_SRCH, 'use_model', 'str', 'Use model for class-average rejection in sieving', 'Use model for class-average rejection in sieving(yes|no){yes}', 'yes|no{yes}', .false., 'yes', group="search", visibility=UI_VIS_ADVANCED)
        call particle_sieving%add_input(UI_SRCH, 'refs', 'file', 'Reference class averages to initialise size compatibility model', 'Reference class averages to initialise size compatibility model', 'e.g. refs.mrcs', .false., '', group="search", visibility=UI_VIS_ADVANCED)
        call particle_sieving%add_input(UI_SRCH, 'single_pass', 'str', 'Coarse pass only', 'Coarse pass only(yes|no){no}', 'yes|no{no}', .false., 'no', group="search", visibility=UI_VIS_ADVANCED)
        ! filter controls
        call particle_sieving%add_input(UI_FILT, 'lpstart', 'num', 'Initial low-pass limit', 'Low-pass limit to be applied in the first &
        &iterations of particle sieving (in Angstroms){15.0}', 'in Angstroms{15.0}', .false., 15., group="filter", visibility=UI_VIS_ADVANCED)
        call particle_sieving%add_input(UI_FILT, 'lpstop_coarse', 'num', 'Final low-pass limit for coarse sieving', 'Low-pass limit to be applied in the last iterations of coarse particle sieving (in Angstroms){15.0}',&
        &'in Angstroms{15.0}', .false., 15., group="filter", visibility=UI_VIS_ADVANCED)
        call particle_sieving%add_input(UI_FILT, 'lpstop_fine', 'num', 'Final low-pass limit for fine sieving', 'Low-pass limit to be applied in the last iterations of fine particle sieving (in Angstroms){10.0}',&
        &'in Angstroms{10.0}', .false., 10., group="filter", visibility=UI_VIS_ADVANCED)
        ! mask controls
        ! <empty>
        ! computer controls
        call particle_sieving%add_input(UI_COMP, nthr, group="compute", visibility=UI_VIS_STANDARD)
        call particle_sieving%add_input(UI_COMP, 'nparts', 'num', 'Number of chunks classified simultaneously', &
        &'Number of particle-subset (chunk) abinitio2D jobs run concurrently on the local machine. Each chunk job &
        &itself runs shared-memory with nthr threads (per-chunk MPI partitioning is not used in offline sieving){1}', &
        &'# of concurrent chunks{1}', .false., 1., group="compute", visibility=UI_VIS_STANDARD)
        call particle_sieving%add_input(UI_COMP, 'walltime', 'num', 'Walltime', 'Maximum execution time for job scheduling and &
        &management(29mins){1740}', 'in seconds(29mins){1740}', .false., 1740., group="compute", &
        &visibility=UI_VIS_ADVANCED)
        ! add to ui_hash
        call add_ui_program('particle_sieving', particle_sieving, prgtab, UI_CATEGORY)
    end subroutine new_particle_sieving

    subroutine new_pick( prgtab )
        class(ui_hash), intent(inout) :: prgtab
        ! PROGRAM SPECIFICATION
        call pick%new(&
        &'pick', &                                                         ! name
        &'Template-based particle picking',&                               ! summary
        &'is a distributed workflow for template-based particle picking',& ! help
        &'simple_exec',&                                                   ! executable
        &.true.,&                                                          ! requires sp_project
        &visibility=UI_VIS_STANDARD)
        ! INPUT PARAMETER SPECIFICATIONS
        ! image input/output
        call pick%add_input(UI_IMG, pickrefs, group="picking", visibility=UI_VIS_STANDARD)
        ! parameter input/output
        call pick%add_input(UI_FILE, 'dir', 'dir', 'Output directory', 'Output directory', 'e.g. pick/', .false., 'pick', group="picking", &
        &visibility=UI_VIS_ADVANCED)
        call pick%add_input(UI_PARM, pcontrast, group="picking", &
        &visibility=UI_VIS_ADVANCED)
        call pick%add_input(UI_PARM, moldiam, group="picking", &
        &visibility=UI_VIS_ADVANCED)
        call pick%add_input(UI_PARM, picker, group="picking", &
        &visibility=UI_VIS_ADVANCED)
        call pick%add_input(UI_PARM, 'nmoldiams', 'num', 'Number of molecular diameters to investigate', 'Number of molecular diameters tested', 'e.g. 5',&
        &.false., 5., group="picking", &
        &visibility=UI_VIS_ADVANCED)
        call pick%add_input(UI_PARM, moldiam_max, group="picking", &
        &visibility=UI_VIS_ADVANCED)
        call pick%add_input(UI_PARM, 'multi_moldiams', 'str', 'Comma-separated molecular diameters with which to execute multiple gaussian pick ',&
        &'Molecular diameters with which to execulte multiple gaussian pick', 'e.g. 100,150', .false., '', group="picking", &
        &visibility=UI_VIS_ADVANCED)
        ! <no additional inputs>
        ! <empty>
        ! search controls
        call pick%add_input(UI_SRCH, 'ndev', 'num', '# of sigmas for outlier detection', '# of standard deviations threshold for outlier detection{2.5}',&
        &'{2.5}', .false., 2.5, group="picking", visibility=UI_VIS_STANDARD)
        call pick%add_input(UI_SRCH, 'ncls', 'num', 'Cluster input pickrefs before template generation', &
        &'If >0, run cluster_cavgs on pickrefs and use medoids only for make_pickrefs{0}', &
        &'# clusters{0}', .false., 0., group="picking", visibility=UI_VIS_ADVANCED)
        call pick%add_input(UI_SRCH, pick_roi, group="picking", &
        &visibility=UI_VIS_ADVANCED)
        call pick%add_input(UI_SRCH, backgr_subtr, group="picking", &
        &visibility=UI_VIS_ADVANCED)
        call pick%add_input(UI_SRCH, particle_density, group="picking", &
        &visibility=UI_VIS_ADVANCED)
        call pick%add_input(UI_SRCH, 'winsz', 'num', 'Window size for sauvola', 'Window size for local sauvola binarisation', 'winsz in pixels ', .false., 32., &
        &visibility=UI_VIS_ADVANCED)
        call pick%add_input(UI_SRCH, nboxes_max, group="picking", &
        &visibility=UI_VIS_ADVANCED)
        ! filter controls
        call pick%add_input(UI_FILT, lp, group="picking", &
        &visibility=UI_VIS_ADVANCED)
        ! mask controls
        ! <empty>
        ! computer controls
        call pick%add_input(UI_COMP, nparts, group="compute", visibility=UI_VIS_STANDARD)
        call pick%add_input(UI_COMP, nthr,   group="compute", visibility=UI_VIS_STANDARD)
        ! add to ui_hash
        call add_ui_program('pick', pick, prgtab, UI_CATEGORY)
    end subroutine new_pick

    subroutine new_preprocess( prgtab )
        class(ui_hash), intent(inout) :: prgtab
        ! PROGRAM SPECIFICATION
        call preprocess%new(&
        &'preprocess', &                                                                    ! name
        &'Run motion correction, CTF estimation, and particle picking',& ! summary
        &'is a distributed workflow that executes motion_correct, ctf_estimate and pick'//& ! help
        &' in sequence',&
        &'simple_exec',&                                                                    ! executable
        &.true., &
        &visibility=UI_VIS_STANDARD)                                                                            ! requires sp_project
        ! INPUT PARAMETER SPECIFICATIONS
        ! image input/output
        call preprocess%add_input(UI_FILE, gainref, &
        &visibility=UI_VIS_ADVANCED)
        ! parameter input/output
        call preprocess%add_input(UI_PARM, total_dose, &
        &visibility=UI_VIS_ADVANCED)
        call preprocess%add_input(UI_PARM, fraction_dose_target, &
        &visibility=UI_VIS_ADVANCED)
        call preprocess%add_input(UI_PARM, max_dose, &
        &visibility=UI_VIS_ADVANCED)
        call preprocess%add_input(UI_PARM, smpd_downscale, &
        &visibility=UI_VIS_ADVANCED)
        call preprocess%add_input(UI_PARM, eer_fraction, &
        &visibility=UI_VIS_ADVANCED)
        call preprocess%add_input(UI_PARM, 'fbody', 'string', 'Template output micrograph name',&
        &'Template output integrated movie name', 'e.g. mic_', .false., 'mic_', &
        &visibility=UI_VIS_ADVANCED)
        call preprocess%add_input(UI_PARM, pspecsz, &
        &visibility=UI_VIS_ADVANCED)
        call preprocess%add_input(UI_PARM, numlen, &
        &visibility=UI_VIS_ADVANCED)
        call preprocess%add_input(UI_PARM, ctfpatch, required_override=.false., &
        &visibility=UI_VIS_ADVANCED)
        call preprocess%add_input(UI_PARM, fit_phshift, required_override=.false., &
        &visibility=UI_VIS_ADVANCED)
        call preprocess%add_input(UI_PARM, flipgain, &
        &visibility=UI_VIS_ADVANCED)
        ! <no additional inputs>
        ! <empty>
        ! search controls
        call preprocess%add_input(UI_SRCH, trs_mc, &
        &visibility=UI_VIS_ADVANCED)
        call preprocess%add_input(UI_SRCH, dfmin, &
        &visibility=UI_VIS_ADVANCED)
        call preprocess%add_input(UI_SRCH, dfmax, &
        &visibility=UI_VIS_ADVANCED)
        call preprocess%add_input(UI_SRCH, astigtol, &
        &visibility=UI_VIS_ADVANCED)
        call preprocess%add_input(UI_SRCH, phshift_min, &
        &visibility=UI_VIS_ADVANCED)
        call preprocess%add_input(UI_SRCH, phshift_max, &
        &visibility=UI_VIS_ADVANCED)
        call preprocess%add_input(UI_SRCH, phshift_step, &
        &visibility=UI_VIS_ADVANCED)
        call preprocess%add_input(UI_SRCH, 'bfac', 'num', 'B-factor applied to frames', 'B-factor applied to frames (in Angstroms^2)', 'in Angstroms^2{50}', .false., 50., &
        &visibility=UI_VIS_ADVANCED)
        call preprocess%add_input(UI_SRCH, mcpatch, &
        &visibility=UI_VIS_ADVANCED)
        call preprocess%add_input(UI_SRCH, nxpatch, &
        &visibility=UI_VIS_ADVANCED)
        call preprocess%add_input(UI_SRCH, nypatch, &
        &visibility=UI_VIS_ADVANCED)
        call preprocess%add_input(UI_SRCH, mcconvention, &
        &visibility=UI_VIS_ADVANCED)
        call preprocess%add_input(UI_SRCH, algorithm, &
        &visibility=UI_VIS_ADVANCED)
        call preprocess%add_input(UI_SRCH, mcpatch_thres, &
        &visibility=UI_VIS_ADVANCED)
        ! filter controls
        call preprocess%add_input(UI_FILT, 'lpstart', 'num', 'Initial low-pass limit for movie alignment', 'Low-pass limit to be applied in the first &
        &iterations of movie alignment(in Angstroms){8}', 'in Angstroms{8}', .false., 8., &
        &visibility=UI_VIS_ADVANCED)
        call preprocess%add_input(UI_FILT, 'lpstop', 'num', 'Final low-pass limit for movie alignment', 'Low-pass limit to be applied in the last &
        &iterations of movie alignment(in Angstroms){5}', 'in Angstroms{5}', .false., 5., &
        &visibility=UI_VIS_ADVANCED)
        call preprocess%add_input(UI_FILT, 'lp_ctf_estimate', 'num', 'Low-pass limit for CTF parameter estimation',&
        & 'Low-pass limit for CTF parameter estimation in Angstroms{5}', 'in Angstroms{5}', .false., 5., &
        &visibility=UI_VIS_ADVANCED)
        call preprocess%add_input(UI_FILT, 'hp_ctf_estimate', 'num', 'High-pass limit for CTF parameter estimation',&
        & 'High-pass limit for CTF parameter estimation  in Angstroms{30}', 'in Angstroms{30}', .false., 30., &
        &visibility=UI_VIS_ADVANCED)
        ! mask controls
        ! <empty>
        ! computer controls
        call preprocess%add_input(UI_COMP, nparts, &
        &visibility=UI_VIS_STANDARD)
        call preprocess%add_input(UI_COMP, nthr, &
        &visibility=UI_VIS_STANDARD)
        ! add to ui_hash
        call add_ui_program('preprocess', preprocess, prgtab, UI_CATEGORY)
    end subroutine new_preprocess

    subroutine new_reextract( prgtab )
        class(ui_hash), intent(inout) :: prgtab
        ! PROGRAM SPECIFICATION
        call reextract%new(&
        &'reextract', &                                                         ! name
        &'Re-extract particle images from integrated movies',&                  ! summary
        &'is a program for re-extracting particle images from integrated movies based on determined 2D/3D shifts',& ! descr long
        &'simple_exec',&                                                        ! executable
        &.true., &
        &visibility=UI_VIS_STANDARD)                                                                ! requires sp_project
        ! INPUT PARAMETER SPECIFICATIONS
        ! image input/output
        ! <empty>
        ! parameter input/output
        call reextract%add_input(UI_PARM, box,     required_override=.false., &
        &visibility=UI_VIS_ADVANCED)
        call reextract%add_input(UI_PARM, oritype, &
        &choices_override=ui_choices([character(len=6) :: 'ptcl2D', 'ptcl3D']), &
        &visibility=UI_VIS_ADVANCED)
        call reextract%add_input(UI_PARM, pcontrast, &
        &visibility=UI_VIS_ADVANCED)
        call reextract%add_input(UI_PARM, backgr_subtr, &
        &visibility=UI_VIS_ADVANCED)
        call reextract%add_input(UI_PARM, outside, &
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
        call reextract%add_input(UI_COMP, nparts, &
        &visibility=UI_VIS_STANDARD)
        call reextract%add_input(UI_COMP, nthr, &
        &visibility=UI_VIS_STANDARD)
        ! add to ui_hash
        call add_ui_program('reextract', reextract, prgtab, UI_CATEGORY)
    end subroutine new_reextract

end module simple_ui_preproc
