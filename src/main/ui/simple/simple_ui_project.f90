!@descr: module defining the user interfaces for project management programs in the simple_exec suite
module simple_ui_project
use simple_ui_modules
implicit none

type(ui_program), target :: export_relion
type(ui_program), target :: export_starproject
type(ui_program), target :: extract_subproj
type(ui_program), target :: import_boxes
type(ui_program), target :: import_cavgs
type(ui_program), target :: import_movies
type(ui_program), target :: import_particles
type(ui_program), target :: import_starproject
type(ui_program), target :: merge_projects
type(ui_program), target :: new_project
type(ui_program), target :: print_project_field
type(ui_program), target :: print_project_info
type(ui_program), target :: prune_project
type(ui_program), target :: replace_project_field
type(ui_program), target :: selection
type(ui_program), target :: update_project
type(ui_program), target :: zero_project_shifts
type(ui_program), target :: write_mic_filetab

contains

    subroutine construct_project_programs(prgtab)
        class(ui_hash), intent(inout) :: prgtab
        call new_export_relion(prgtab)
        call new_export_starproject(prgtab)
        call new_extract_subproj(prgtab)
        call new_import_boxes(prgtab)
        call new_import_cavgs(prgtab)
        call new_import_movies(prgtab)
        call new_import_particles(prgtab)
        call new_import_starproject(prgtab)
        call new_merge_projects(prgtab)
        call new_new_project(prgtab)
        call new_print_project_field(prgtab)
        call new_print_project_info(prgtab)
        call new_prune_project(prgtab)
        call new_replace_project_field(prgtab)
        call new_selection(prgtab)
        call new_update_project(prgtab)
        call new_zero_project_shifts(prgtab)
        call new_write_mic_filetab(prgtab)
    end subroutine construct_project_programs

    subroutine print_project_programs(logfhandle)
        integer, intent(in) :: logfhandle
        write(logfhandle,'(A)') format_str('PROJECT MANAGEMENT:', C_UNDERLINED)
        write(logfhandle,'(A)') export_relion%name%to_char()
        write(logfhandle,'(A)') export_starproject%name%to_char()
        write(logfhandle,'(A)') extract_subproj%name%to_char()
        write(logfhandle,'(A)') import_boxes%name%to_char()
        write(logfhandle,'(A)') import_cavgs%name%to_char()
        write(logfhandle,'(A)') import_movies%name%to_char()
        write(logfhandle,'(A)') import_particles%name%to_char()
        write(logfhandle,'(A)') import_starproject%name%to_char()
        write(logfhandle,'(A)') merge_projects%name%to_char()
        write(logfhandle,'(A)') new_project%name%to_char()
        write(logfhandle,'(A)') print_project_field%name%to_char()
        write(logfhandle,'(A)') print_project_info%name%to_char()
        write(logfhandle,'(A)') prune_project%name%to_char()
        write(logfhandle,'(A)') replace_project_field%name%to_char()
        write(logfhandle,'(A)') selection%name%to_char()
        write(logfhandle,'(A)') update_project%name%to_char()
        write(logfhandle,'(A)') zero_project_shifts%name%to_char()
        write(logfhandle,'(A)') write_mic_filetab%name%to_char()
        write(logfhandle,'(A)') ''
    end subroutine print_project_programs

    subroutine new_export_relion( prgtab )
        class(ui_hash), intent(inout) :: prgtab
        ! PROGRAM SPECIFICATION
        call export_relion%new(&
        &'export_relion',&                                              ! name
        &'Export project to relion ',&                                  ! descr_short
        &'is a program to export simple project to relion',&
        &'simple_exec',&                                                ! executable
        &.true.)                                                        ! requires sp_project
        ! INPUT PARAMETER SPECIFICATIONS
        ! image input/output
        ! <empty>
        ! parameter input/output
        call export_relion%add_input(UI_PARM, 'tiltgroupmax', 'num', 'Max movies in a tilt/optics group', &
            'Sub-divide beamtilt/optics groups', '0', .false., 0.0)
        call export_relion%add_input(UI_PARM, 'reliongroups', 'num', 'Number of Relion groups based on defocus', &
            'Divide particles into X groups based on defocus for relion', '# micrographs', .false., 0.0)
        call export_relion%add_input(UI_PARM, 'xmlloc', 'file', 'Pathname of EPU XML files',&
            'Pathname of EPU XML files ', 'e.g. /data/xml', .false., 'NONE')
        call export_relion%add_input(UI_PARM, 'tilt_thres', 'num', 'Distance threshold',&
            'Distance threshold for hierarchical clustering of beamtilt/shift groups ', '{0.05}', .false., 0.05)
        call export_relion%add_input(UI_PARM, 'optics_offset', 'num', 'Offset to apply to optics group numbering', &
            'Offset to apply to optics group numbering', '{0}', .false., 0.0)
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
        call add_ui_program('export_relion', export_relion, prgtab)
    end subroutine new_export_relion

    subroutine new_export_starproject( prgtab )
        class(ui_hash), intent(inout) :: prgtab
        ! PROGRAM SPECIFICATION
        call export_starproject%new(&
        &'export_starproject', &                                                ! name
        &'Export projectfile in star format',&                                  ! descr_short
        &'is a program to export a SIMPLE projectfile in star format',&         ! descr long
        &'simple_exec',&                                                        ! executable
        &.true.)                                                                ! requires sp_project
        ! INPUT PARAMETER SPECIFICATIONS
        ! image input/output
        ! call export_starproject%set_input('parm_ios', 1, projfile)
        ! parameter input/output
        ! alternative inputs
        ! <empty>
        ! search controls
        ! <empty>
        ! filter controls
        ! <empty>
        ! mask controls
        ! <empty>
        ! computer controls
        ! add to ui_hash
        call add_ui_program('export_starproject', export_starproject, prgtab)
    end subroutine new_export_starproject

    subroutine new_extract_subproj( prgtab )
        class(ui_hash), intent(inout) :: prgtab
        ! PROGRAM SPECIFICATION
        call extract_subproj%new(&
        &'extract_subproj',&                                                                                     ! name
        &'extraction of a subproject of time-series of metallic nanoparticles',&                                 ! descr_short
        &'is a shared-memory workflow for extraction of a subproject of time-series of metallic nanoparticles',& ! descr_long
        &'all',&                                                                                                 ! executable
        &.true., gui_advanced=.false.)                                                                           ! requires sp_project
        ! INPUT PARAMETER SPECIFICATIONS
        ! image input/output
        ! <empty>
        ! parameter input/output
        call extract_subproj%add_input(UI_PARM, projfile)
        call extract_subproj%add_input(UI_PARM, 'subprojname', 'str', 'Subproject name', 'Name of subproject to create ./myproject/myproject.simple',&
        &'e.g. to create ./myproject/myproject.simple', .true., '')
        ! alternative inputs
        call extract_subproj%add_input(UI_ALT, 'fromp',    'num', 'From index', 'Start index for extraction', 'start index', .false., 1.0)
        call extract_subproj%add_input(UI_ALT, 'top',      'num', 'To index', 'Stop index for extraction', 'stop index', .false., 1.0)
        call extract_subproj%add_input(UI_ALT, 'clustind', 'num', 'Cluster index', 'Cluster index', 'e.g. 5', .false., 0.)
        ! search controls
        ! <empty>
        ! filter controls
        ! <empty>
        ! mask controls
        ! <empty>
        ! computer controls
        ! <empty>
        ! add to ui_hash
        call add_ui_program('extract_subproj', extract_subproj, prgtab)
    end subroutine new_extract_subproj

    subroutine new_import_boxes( prgtab )
        class(ui_hash), intent(inout) :: prgtab
        ! PROGRAM SPECIFICATION
        call import_boxes%new(&
        &'import_boxes',&                                  ! name
        &'Import EMAN box coordinates to SIMPLE project',& ! descr_short
        &'is a program for importing EMAN1.9 box coordinates to the project. The *box (text) files should be listed in boxtab',&
        &'simple_exec',&                                   ! executable
        &.true.)                                           ! requires sp_project
        ! INPUT PARAMETER SPECIFICATIONS
        ! image input/output
        ! <empty>
        ! parameter input/output
        call import_boxes%add_input(UI_PARM, 'boxtab', 'file', 'List of box files', &
            'List of per-micrograph box files (*.box) to import', 'e.g. boxes.txt', .true., '')
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
        call add_ui_program('import_boxes', import_boxes, prgtab)
    end subroutine new_import_boxes

    subroutine new_import_cavgs( prgtab )
        class(ui_hash), intent(inout) :: prgtab
        ! PROGRAM SPECIFICATION
        call import_cavgs%new(&
        &'import_cavgs',&                                        ! name
        &'Import class averages to SIMPLE project',&             ! descr_short
        &'is a program for importing class averages movies to the project',&
        &'simple_exec',&                                         ! executable
        &.true.)                                                 ! requires sp_project
        call import_cavgs%add_input(UI_IMG, 'stk', 'file', 'Stack of class averages',&
        &'Stack of class average images to import', 'e.g. cavgs.mrcs', .true., '')
        ! parameter input/output
        call import_cavgs%add_input(UI_PARM, smpd)
        ! add to ui_hash
        call add_ui_program('import_cavgs', import_cavgs, prgtab)
    end subroutine new_import_cavgs

    subroutine new_import_movies( prgtab )
        class(ui_hash), intent(inout) :: prgtab
        ! PROGRAM SPECIFICATION
        call import_movies%new(&
        &'import_movies',&                                       ! name
        &'Import movies to SIMPLE project',&                     ! descr_short
        &'is a program for importing DDD movies to the project. The movies can be located in any read-only location&
        & accessible to the project. If the movies contain only a single frame, they will be interpreted as motion-corrected&
        & and integrated. Box files (in EMAN format) can be imported along with the movies',&
        &'simple_exec',&                                         ! executable
        &.true.)                                                 ! requires sp_project
        ! INPUT PARAMETER SPECIFICATIONS
        ! image input/output
        call import_movies%add_input(UI_IMG, 'filetab',    'file', 'List of movie files',    'List of movie files (*.mrcs) to import', 'e.g. movies.txt', .false., '')
        call import_movies%add_input(UI_IMG, 'dir_movies', 'dir',  'Input movies directory', 'Where the movies to process are located or will squentially appear', 'e.g. /cryodata/', .false., 'preprocess/')
        ! parameter input/output
        call import_movies%add_input(UI_PARM, smpd)
        call import_movies%add_input(UI_PARM, kv,    required_override=.true.)
        call import_movies%add_input(UI_PARM, cs,    required_override=.true.)
        call import_movies%add_input(UI_PARM, fraca, required_override=.true.)
        call import_movies%add_input(UI_PARM, ctf_yes)
        call import_movies%add_input(UI_PARM, phaseplate)
        call import_movies%add_input(UI_PARM, 'boxtab', 'file', 'List of box files', 'List of per-micrograph box files (*.box) to import', 'e.g. boxes.txt', .false., '')
        call import_movies%add_input(UI_PARM, 'deftab', 'file','Pre-determined per-micrograph CTF parameters',&
        &'List of CTF parmeters for micrographs import only', 'e.g. deftab.txt', .false., '')
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
        call add_ui_program('import_movies', import_movies, prgtab)
    end subroutine new_import_movies

    subroutine new_import_particles( prgtab )
        class(ui_hash), intent(inout) :: prgtab
        ! PROGRAM SPECIFICATION
        call import_particles%new(&
        &'import_particles',&                                       ! name
        &'Import particles to SIMPLE project',&                     ! descr_short
        &'is a program for importing extracted particle images to the project',&
        &'all',&                                                    ! executable
        &.true.)                                                    ! requires sp_project
        ! INPUT PARAMETER SPECIFICATIONS
        ! image input/output
        ! <empty>
        ! parameter input/output
        call import_particles%add_input(UI_PARM, smpd)
        call import_particles%add_input(UI_PARM, kv,    required_override=.true.)
        call import_particles%add_input(UI_PARM, cs,    required_override=.true.)
        call import_particles%add_input(UI_PARM, fraca, required_override=.true.)
        call import_particles%add_input(UI_PARM, ctf_yes)
        call import_particles%add_input(UI_PARM, phaseplate)
        call import_particles%add_input(UI_PARM, oritab)
        call import_particles%add_input(UI_PARM, deftab)
        call import_particles%add_input(UI_PARM, 'plaintexttab', 'file', 'Plain text file of input parameters',&
        'Plain text file of tabulated per-particle input parameters: dfx, dfy, angast, phshift', 'e.g. params.txt', .false., '')
        call import_particles%add_input(UI_PARM, 'dfunit', 'multi', 'Underfocus unit', 'Underfocus unit(A|microns){microns}', '(A|microns){microns}', .false., 'microns')
        call import_particles%add_input(UI_PARM, 'angastunit', 'multi', 'Angle of astigmatism unit', 'Angle of astigmatism unit(radians|degrees){degrees}', '(radians|degrees){degrees}', .false., 'degrees')
        call import_particles%add_input(UI_PARM, 'phshiftunit', 'multi', 'Phase-shift unit', 'Phase-shift unit(radians|degrees){radians}', '(radians|degrees){radians}', .false., 'degrees')
        ! alternative inputs
        call import_particles%add_input(UI_ALT, 'stktab', 'file', 'List of per-micrograph particle stacks',&
        &'List of per-micrograph particle image stacks to import', 'per-micrograph stack list; e.g. stktab.txt', .false., '')
        call import_particles%add_input(UI_ALT, 'stk', 'file', 'Stack of particles',&
        &'Stack of particle images to import', 'e.g. stk.mrcs', .false., '')
        ! search controls
        ! <empty>
        ! filter controls
        ! <empty>
        ! mask controls
        ! <empty>
        ! computer controls
        ! <empty>
        ! add to ui_hash
        call add_ui_program('import_particles', import_particles, prgtab)
    end subroutine new_import_particles

    subroutine new_import_starproject( prgtab )
        class(ui_hash), intent(inout) :: prgtab
        ! PROGRAM SPECIFICATION
        call import_starproject%new(&
        &'import_starproject', &                                                ! name
        &'Import project in in star format',&                                   ! descr_short
        &'is a program to import a SIMPLE projectfile from star format',&       ! descr long
        &'simple_exec',&                                                        ! executable
        &.false.,&                                                              ! requires sp_project
        &gui_advanced=.false., gui_submenu_list = "data")                       ! GUI
        ! INPUT PARAMETER SPECIFICATIONS
        ! image input/output
        ! parameter input/output
        call import_starproject%add_input(UI_PARM, 'import_dir',  'dir',   'Import directory', 'Directory to import data from. In auto mode should be output &
        &from an external job e.g. relion', 'e.g. MotionCorr/job001', .true., '', gui_submenu="data", gui_advanced=.false.)
        call import_starproject%add_input(UI_PARM, 'starfile',    'file',  'Metadata starfile', 'Path to starfile containing micrograph/particle metadata. Only &
        &required when not using auto mode', 'e.g. micrographs.star', .false., '', gui_submenu="data", gui_advanced=.false.)
        call import_starproject%add_input(UI_PARM, 'import_type', 'multi', 'Import type', 'Type of data contained in starfile (auto|mic|ptcl2D|ptcl3D){auto}. &
        &Auto mode (default) will attempt to determine this automatically', '(auto|mic|ptcl2D|ptcl3D){auto}', .false., 'auto', gui_submenu="data", gui_advanced=.false.)
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
        call add_ui_program('import_starproject', import_starproject, prgtab)
    end subroutine new_import_starproject

    subroutine new_merge_projects( prgtab )
        class(ui_hash), intent(inout) :: prgtab
        ! PROGRAM SPECIFICATION
        call merge_projects%new(&
        &'merge_projects', &                                            ! name
        &'Merge two projects',&                                         ! descr_short
        &'is a program to merge two projects into one', &               ! descr_long
        &'simple_exec',&                                                ! executable
        &.true.)                                                        ! requires sp_project
        ! INPUT PARAMETER SPECIFICATIONS
        ! image input/output
        ! <empty>
        ! parameter input/output
        call merge_projects%add_input(UI_PARM, projfile_target,&
        &descr_short_override       = 'Project to merge',&
        &descr_long_override        = 'Location of project file to append',&
        &gui_submenu="data", gui_advanced=.false.)
        call merge_projects%add_input(UI_PARM, oritype,&
        &descr_long_override        = 'Oritype segment in project(ptcl2D|ptcl3D){ptcl2D}',&
        &descr_placeholder_override = '(ptcl2D|ptcl3D){ptcl2D}',&
        gui_submenu="extract", gui_advanced=.false.)
        call merge_projects%add_input(UI_PARM, box, required_override=.false., gui_submenu="extract", gui_advanced=.false.)
        call merge_projects%add_input(UI_PARM, pcontrast, gui_submenu="extract")
        call merge_projects%add_input(UI_PARM, backgr_subtr, gui_submenu="extract")
        call merge_projects%add_input(UI_PARM, outside, gui_submenu="extract")
        ! alternative inputs
        ! <empty>
        ! search controls
        ! <empty>
        ! filter controls
        ! <empty>
        ! mask controls
        ! <empty>
        ! computer controls
        call merge_projects%add_input(UI_COMP, nparts, gui_submenu="compute")
        call merge_projects%add_input(UI_COMP, nthr,   gui_submenu="compute")
        ! add to ui_hash
        call add_ui_program('merge_projects', merge_projects, prgtab)
    end subroutine new_merge_projects

    subroutine new_new_project( prgtab )
        class(ui_hash), intent(inout) :: prgtab
        ! PROGRAM SPECIFICATION
        call new_project%new(&
        &'new_project',&                     ! name
        &'Create a new project',&            ! descr_short
        &'is a program for creating a new project. SIMPLE3.0 relies on a monolithic project file for controlling &
        &execution on distributed and shared-memory systems and for unified meta-data management. This program &
        &creates a directory named projname and a file projname.simple inside that directory that contains all &
        &information about the project as well as all meta data generated by the different SIMPLE programs. This &
        &file is mirrored by an abstract data type in the back-end, which manages the parameters and &
        &meta-data I/O required for execution of SIMPLE',& ! descr_longg
        &'all',&                             ! executable
        &.false.)                            ! requires sp_project
        ! INPUT PARAMETER SPECIFICATIONS
        ! image input/output
        ! <empty>
        ! parameter input/output
        call new_project%add_input(UI_PARM, user_email)
        call new_project%add_input(UI_PARM, projfile, required_override=.false.)
        ! alternative inputs
        call new_project%add_input(UI_ALT, projname, required_override=.false.)
        call new_project%add_input(UI_ALT, 'dir', 'dir', 'Project directory', 'Project directory', 'give dir', .false., '')
        ! <empty>
        ! search controls
        ! <empty>
        ! filter controls
        ! <empty>
        ! mask controls
        ! <empty>
        ! computer controls
        call new_project%add_input(UI_COMP, time_per_image)
        call new_project%add_input(UI_COMP, user_account)
        call new_project%add_input(UI_COMP, user_project)
        call new_project%add_input(UI_COMP, qsys_partition)
        call new_project%add_input(UI_COMP, qsys_qos)
        call new_project%add_input(UI_COMP, qsys_reservation)
        call new_project%add_input(UI_COMP, job_memory_per_task)
        call new_project%add_input(UI_COMP, qsys_name)
        call new_project%add_input(UI_COMP, walltime)
        ! add to ui_hash
        call add_ui_program('new_project', new_project, prgtab)
    end subroutine new_new_project

    subroutine new_print_project_field( prgtab )
        class(ui_hash), intent(inout) :: prgtab
        ! PROGRAM SPECIFICATION
        call print_project_field%new(&
        &'print_project_field', &                                             ! name
        &'Print project field',&                                              ! descr_short
        &'is a program for printing an orientation field in the project data structure (segment in *.simple project file)',&  ! descr_long
        &'all',&                                                              ! executable
        &.true.)                                                              ! requires sp_project
        ! INPUT PARAMETER SPECIFICATIONS
        ! image input/output
        ! <empty>
        ! parameter input/output
        call print_project_field%add_input(UI_PARM, oritype, required_override=.true.)
        call print_project_field%add_input(UI_PARM, 'json',     'binary', 'output in JSON format', 'output in JSON format', '(yes|no){no}', .false., 'no')
        call print_project_field%add_input(UI_PARM, 'boxes',    'binary', 'output coordinates in JSON format', 'output coordinates in JSON format', '(yes|no){no}', .false., 'no')
        call print_project_field%add_input(UI_PARM, 'sort',     'string', 'sort oris on key', 'sort oris on key', 'e.g. ctfres', .false., '')
        call print_project_field%add_input(UI_PARM, 'plot_key', 'string', 'plot plot_key on , sort on x', 'plot plot_key on , sort on x', 'e.g. dfx', .false., '')
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
        call add_ui_program('print_project_field', print_project_field, prgtab)
    end subroutine new_print_project_field

    subroutine new_print_project_info( prgtab )
        class(ui_hash), intent(inout) :: prgtab
        ! PROGRAM SPECIFICATION
        call print_project_info%new(&
        &'print_project_info', &                                           ! name
        &'Print project info',&                                            ! descr_short
        &'is a program prints information about a *.simple project file',& ! descr_long
        &'all',&                                                           ! executable
        &.true.)                                                           ! requires sp_project
        ! INPUT PARAMETER SPECIFICATIONS
        ! image input/output
        ! <empty>
        ! parameter input/output
        ! <empty>
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
        call add_ui_program('print_project_info', print_project_info, prgtab)
    end subroutine new_print_project_info

    subroutine new_prune_project( prgtab )
        class(ui_hash), intent(inout) :: prgtab
        ! PROGRAM SPECIFICATION
        call prune_project%new(&
        &'prune_project',&                            ! name
        &'discards deselected data from a project',&  ! descr_short
        &'is a program for discarding deselected data (particles,stacks) from a project',& ! descr_long
        &'all',&                                      ! executable
        &.true.)                                      ! requires sp_project
        ! INPUT PARAMETER SPECIFICATIONS
        ! image input/output
        ! <empty>
        ! parameter input/output
        call prune_project%add_input(UI_PARM, 'state', 'num', 'State index', 'Index of state to extract', 'give state index', .false., 1.)
        ! alternative inputs
        ! <empty>
        ! search controls
        ! <empty>
        ! filter controls
        ! <empty>
        ! mask controls
        ! <empty>
        ! computer controls
        call prune_project%add_input(UI_COMP, nparts)
        ! add to ui_hash
        call add_ui_program('prune_project', prune_project, prgtab)
    end subroutine new_prune_project

    subroutine new_replace_project_field( prgtab )
        class(ui_hash), intent(inout) :: prgtab
        ! PROGRAM SPECIFICATION
        call replace_project_field%new(&
        &'replace_project_field',&                    ! name
        &'hard substitution of project field',&       ! descr_short
        &'is a program for hard substitution of project field, for development purposes',& ! descr_long
        &'simple_exec',&                              ! executable
        &.false.)                                     ! requires sp_project
        ! INPUT PARAMETER SPECIFICATIONS
        ! image input/output
        ! <empty>
        ! parameter input/output
        call replace_project_field%add_input(UI_PARM, projfile)
        call replace_project_field%add_input(UI_PARM, projfile_target)
        call replace_project_field%add_input(UI_PARM, oritype, required_override=.true.)
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
        call add_ui_program('replace_project_field', replace_project_field, prgtab)
    end subroutine new_replace_project_field

    subroutine new_selection( prgtab )
        class(ui_hash), intent(inout) :: prgtab 
        ! PROGRAM SPECIFICATION
        call selection%new(&
        &'selection',&                                                                  ! name
        &'Reports external selection through state 0/1 tags to project',&               ! descr_short
        &'is a program for reporting external (GUI) selections to the SIMPLE project',& ! descr_long
        &'simple_exec',&                                                                ! executable
        &.true.)                                                                        ! requires sp_project
        ! INPUT PARAMETER SPECIFICATIONS
        ! image input/output
        ! <empty>
        ! parameter input/output
        call selection%add_input(UI_PARM, oritype)
        call selection%add_input(UI_PARM, 'state',           'num',    'State number', 'Map selection to oris with this state only', '{1}', .false., 1.0)
        call selection%add_input(UI_PARM, prune)
        call selection%add_input(UI_PARM, 'append',          'binary', 'Append selection to existing', 'Previously deselected particles will stay deselected(yes|no){no}', '(yes|no){no}', .false., 'no')
        call selection%add_input(UI_PARM, 'balance',         'binary', 'Balanced selection of particles across classes', 'Balanced selection(yes|no){no}', '(yes|no){no}', .false., 'no')
        call selection%add_input(UI_PARM, 'nptcls_per_part', 'num',    'Number of ptcls per part to select when balancing', '# ptcls per part after balancing', '{100000}', .false., 0.0)
        call selection%add_input(UI_PARM, 'greedy_sampling', 'binary', 'Greedy balanced selection', 'Greedy balanced selection(yes|no){yes}', '(yes|no){yes}', .false., 'yes')
        call selection%add_input(UI_PARM, 'nparts',          'num',    'Number of partitions in balancing', '# balanced parts', '# balanced parts', .false., 1.)
        call selection%add_input(UI_PARM, dfmin)
        ! alternative inputs
        call selection%add_input(UI_ALT, 'infile', 'file', 'File with selection state (0/1) flags', 'Plain text file (.txt) with selection state (0/1) flags',&
        &'give .txt selection file', .false., '')
        call selection%add_input(UI_ALT, 'deselfile', 'file', 'File with deselection indices', 'Plain text file (.txt) with deselection indices',&
        &'give .txt deselection file', .false., '')
        call selection%add_input(UI_ALT, nran)
        call selection%add_input(UI_ALT, ctfresthreshold)
        call selection%add_input(UI_ALT, icefracthreshold)
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
        call add_ui_program('selection', selection, prgtab)
    end subroutine new_selection

    subroutine new_update_project( prgtab )
        class(ui_hash), intent(inout) :: prgtab
        ! PROGRAM SPECIFICATION
        call update_project%new(&
        &'update_project',&                  ! name
        &'Update an existing project',&      ! descr_short
        &'is a program for updating an existing project: changing the name/user_email/computer controls',& ! descr_long
        &'all',&                             ! executable
        &.true.)                             ! requires sp_project
        ! INPUT PARAMETER SPECIFICATIONS
        ! image input/output
        ! <empty>
        ! parameter input/output
        call update_project%add_input(UI_PARM, projfile)
        call update_project%add_input(UI_PARM, user_email)
        ! alternative inputs
        ! <empty>
        ! search controls
        ! <empty>
        ! filter controls
        ! <empty>
        ! mask controls
        ! <empty>
        ! computer controls
        call update_project%add_input(UI_COMP, time_per_image)
        call update_project%add_input(UI_COMP, user_account)
        call update_project%add_input(UI_COMP, user_project)
        call update_project%add_input(UI_COMP, qsys_partition)
        call update_project%add_input(UI_COMP, qsys_qos)
        call update_project%add_input(UI_COMP, qsys_reservation)
        call update_project%add_input(UI_COMP, job_memory_per_task)
        call update_project%add_input(UI_COMP, qsys_name)
        call update_project%add_input(UI_COMP, walltime)
        ! add to ui_hash
        call add_ui_program('update_project', update_project, prgtab)
    end subroutine new_update_project

    subroutine new_zero_project_shifts( prgtab )
        class(ui_hash), intent(inout) :: prgtab
        ! PROGRAM SPECIFICATION
        call zero_project_shifts%new(&
        &'zero_project_shifts', &                                       ! name
        &'zero_project_shifts',&                                        ! descr_short
        &'is a program that zeroes the shifts in the ptcl2D/ptcl3D fields in the project',& ! descr_long
        &'simple_exec',&                                                ! executable
        &.true.)                                                        ! requires sp_project
        ! INPUT PARAMETER SPECIFICATIONS
        ! image input/output
        ! <empty>
        ! parameter input/output
        ! <empty>
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
        call add_ui_program('zero_project_shifts', zero_project_shifts, prgtab)
    end subroutine new_zero_project_shifts

    subroutine new_write_mic_filetab( prgtab )
        class(ui_hash), intent(inout) :: prgtab
        ! PROGRAM SPECIFICATION
        call write_mic_filetab%new(&
        &'write_mic_filetab',&                                            ! name
        &'Writes a filetable of state > 0 micrographs',&                  ! descr_short
        &'is a program for writing a filetable of selected micrographs',& ! descr_long
        &'simple_exec',&                                                  ! executable
        &.true.)                                                          ! requires sp_project
        ! INPUT PARAMETER SPECIFICATIONS
        ! image input/output
        call write_mic_filetab%add_input(UI_IMG, 'fname', 'file', 'Filename micrograph list', 'Filename for list of micrograph files (*.mrc)', 'e.g. mics.txt', .true., '')
        ! parameter input/output
        call write_mic_filetab%add_input(UI_PARM, projfile)
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
        call add_ui_program('write_mic_filetab', write_mic_filetab, prgtab)
    end subroutine new_write_mic_filetab

end module simple_ui_project
