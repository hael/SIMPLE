!@descr: state and steps of stream task 4 (reference-based picking): pick and extract the accepted micrographs of each preprocessed batch
!==============================================================================
! MODULE: simple_stream_stage_refpick
!
! PURPOSE:
!   The body of stream p04 as a type. Each completed preprocessing project
!   becomes one job set: its accepted micrographs, picked with the picking
!   references and extracted by one pick_extract job. Finished sets are moved
!   to the completed folder; their micrographs join the stage's project, and
!   the stacks and particles of every set are assembled into it when the
!   project is written. The commander (simple_commanders_stream_p04_refpick_extract) only
!   normalises the command line and loops over iterate() until finished().
!
!   Micrographs are not thresholded here: preprocessing has already rejected
!   them in the projects it hands on (and follows the GUI's threshold updates,
!   which this stage does not receive).
!
!   make_pickrefs writes moldiam.txt in the stage directory before any set is
!   submitted; the particle-sieving stage reads its mask diameter there.
!
!   What happens is delegated:
!     - job sets (naming, submission, collection, restart) -> simple_stream_job_sets
!     - picking references                  -> make_pickrefs (once, in this process)
!     - micrograph import                   -> simple_mic_import
!     - STAR files                          -> starproject_stream
!     - optics groups                       -> the newest map of optics assignment (simple_optics_maps)
!     - GUI                                 -> simple_stream_pipe, simple_stream_gui_senders
!
! LIFECYCLE:
!   new(cline) -> { iterate() } until finished() -> finalize() -> kill()
!
! RESTART:
!   The completed sets are imported again, and the upstream projects they were
!   made from are put in the watcher history, so none is picked twice.
!==============================================================================
module simple_stream_stage_refpick
use simple_defs,                        only: logfhandle, STDLEN, PATH_HERE, PATH_PARENT
use simple_defs_fname,                  only: TERM_STREAM, PICKREFS_FBODY, STREAM_MOLDIAM, DIR_PICKER,&
                                             &DIR_EXTRACT
use simple_defs_stream,                 only: DIR_STREAM, DIR_STREAM_COMPLETED, SHORTWAIT, WAITTIME, INACTIVE_TIME
use simple_defs_environment,            only: SIMPLE_STREAM_PICK_PARTITION
use simple_error,                       only: simple_exception
use simple_string,                      only: string
use simple_fileio,                      only: del_file, file_exists, simple_abspath, simple_getcwd, simple_rmdir
use simple_syslib,                      only: dir_exists, simple_mkdir
use simple_timer,                       only: simple_gettime, cast_time_char
use simple_cmdline,                     only: cmdline
use simple_oris,                        only: oris
use simple_parameters,                  only: parameters
use simple_sp_project,                  only: sp_project
use simple_qsys_env,                    only: qsys_env
use simple_qsys_funs,                   only: qsys_cleanup
use simple_starproject_stream,          only: starproject_stream
use simple_stream_watcher,              only: stream_watcher
use simple_stream_state,                only: ipc_pipe_refpick_in
use simple_stream_utils,                only: create_stream_project
use simple_commanders_pick,             only: commander_make_pickrefs
use simple_ptcl_sieve,                  only: DEFAULT_COARSE_BOX
use simple_gui_utils,                   only: mrc2jpeg_tiled
use simple_gui_metadata_utils,          only: max_metadata_size
use simple_gui_metadata_types,          only: GUI_METADATA_STREAM_REFERENCE_PICKING_TYPE,&
                                             &GUI_METADATA_STREAM_REFERENCE_PICKING_MICROGRAPH_TYPE,&
                                             &GUI_METADATA_STREAM_REFERENCE_PICKING_CLS2D_TYPE
use simple_gui_metadata_micrograph,     only: gui_metadata_micrograph
use simple_gui_metadata_cavg2D,         only: gui_metadata_cavg2D
use simple_gui_metadata_stream_picking, only: gui_metadata_stream_picking
use simple_stream_pipe,                 only: stream_pipe
use simple_stream_gui_senders,          only: send_cavgs, send_recent_micrographs
use simple_mic_import,                  only: append_mics_from_projects
use simple_optics_maps,                 only: import_latest_optics_map
use simple_stream_job_sets,             only: stream_job_sets
implicit none

public :: stream_stage_refpick
private
#include "simple_local_flags.inc"

integer, parameter :: NTHUMB_MAX          = 20                 ! most recent micrograph thumbnails sent to the GUI
integer, parameter :: DEFAULT_EXTRACT_BOX = DEFAULT_COARSE_BOX ! smallest extraction box
integer, parameter :: STAR_EVERY_NMICS    = 1000 ! below this many micrographs the STAR file is rewritten on every import...
integer, parameter :: STAR_STEP_NMICS     = 100  ! ...above it, every this many new micrographs

! Components and steps are public so simple_stream_stage_refpick_tester can assemble a stage
! and run one step at a time; production code uses new/iterate/finished/finalize/kill.
type :: stream_stage_refpick
    type(parameters)                  :: params
    type(cmdline)                     :: cline_exec       ! command line of the pick_extract worker jobs
    type(cmdline)                     :: cline_pickrefs   ! make_pickrefs, run once the pixel size is known
    type(qsys_env)                    :: qenv
    type(sp_project)                  :: spproj           ! every imported micrograph; stacks and particles when written
    type(stream_job_sets)             :: sets             ! one per upstream project, picked and extracted by one job
    type(stream_watcher)              :: project_buff     ! completed preprocessing projects
    type(starproject_stream)          :: starproj_stream
    type(stream_pipe)                 :: pipe             ! to the master
    type(gui_metadata_stream_picking) :: meta_status
    type(gui_metadata_micrograph)     :: meta_micrograph
    type(gui_metadata_cavg2D)         :: meta_pickrefs
    type(string), allocatable         :: set_projects(:)  ! imported sets (completed folder), in import order
    type(string), allocatable         :: restored_sources(:) ! restart: upstream projects for the watcher history
    type(string)                      :: cwd              ! absolute stage directory
    integer :: n_mics_submitted = 0       ! micrographs sent to picking
    integer :: nptcls_glob      = 0       ! particles of the imported micrographs
    integer :: n_failed_jobs    = 0
    integer :: nmic_star        = 0       ! micrographs in the last STAR snapshot (above STAR_EVERY_NMICS)
    integer :: prev_stacksz     = 0
    integer :: last_injection   = 0       ! time of the last import
    logical :: l_attached       = .false. ! the upstream completed-projects folder exists and is watched
    logical :: l_waiting_logged = .false.
    logical :: l_pickrefs_found = .false. ! the input picking references exist
    logical :: l_pickrefs_ready = .false. ! make_pickrefs has run; the workers have their references
    logical :: l_projects_left  = .false. ! the last watch was capped, more upstream projects may be waiting
    logical :: l_haschanged     = .false. ! imports since the last idle project write
    logical :: l_restart        = .false.
    logical :: l_exists         = .false.
    ! waits (s); tests set them to 0, and settle_s to -1 to take files written in the same second
    integer :: settle_s         = SHORTWAIT ! an upstream project is taken once untouched longer than this
    integer :: wait_s           = WAITTIME  ! pause at the end of a pass
contains
    procedure :: new
    procedure :: iterate
    procedure :: finished
    procedure :: finalize
    procedure :: kill
    ! the steps of new(), in order
    procedure :: init_params
    procedure :: init_job_dirs
    procedure :: resume_previous_run
    procedure :: init_queue
    procedure :: build_worker_cline
    procedure :: init_gui
    ! the steps of iterate() and their helpers
    procedure :: attach_upstream
    procedure :: pickrefs_available
    procedure :: submit_new_projects
    procedure :: create_set_project
    procedure :: prepare_pickrefs
    procedure :: schedule_jobs
    procedure :: collect_jobs
    procedure :: import_finished_sets
    procedure :: import_sets
    procedure :: import_previous_sets
    procedure :: process_imports
    procedure :: idle
    procedure :: apply_optics_map
    procedure :: write_mic_star
    procedure :: write_project
    procedure :: send_status
end type stream_stage_refpick

contains

    !---------------- lifecycle ----------------

    !> Builds the parameters from the normalised @p cline, re-imports a previous run when
    !! restarting, and prepares job submission and the GUI.
    subroutine new( self, cline )
        class(stream_stage_refpick), intent(inout) :: self
        class(cmdline),              intent(inout) :: cline
        call self%kill()
        call self%init_params(cline)
        self%l_exists = .true. ! from here kill() releases what has been built
        call self%init_job_dirs()
        if( self%l_restart ) call self%resume_previous_run()
        call self%init_queue()
        call self%build_worker_cline(cline)
        call self%init_gui(-1, ipc_pipe_refpick_in(2))
        self%last_injection = simple_gettime()
        call self%send_status(string('waiting for picking references'))
    end subroutine new

    !> The stage's project file (made and given a computing environment) and its parameters; the
    !! project must start without micrographs.
    subroutine init_params( self, cline )
        class(stream_stage_refpick), intent(inout) :: self
        class(cmdline),              intent(inout) :: cline
        type(string) :: dir_exec
        ! a restart is recognised by its output directory (or dir_exec), before params%new creates one
        self%l_restart = .false.
        if( cline%defined('outdir') ) self%l_restart = dir_exists(cline%get_carg('outdir'))
        if( cline%defined('dir_exec') )then
            dir_exec = cline%get_carg('dir_exec')
            if( .not. file_exists(dir_exec) ) THROW_HARD('Previous directory does not exist: '//dir_exec%to_char())
            self%l_restart = .true.
        endif
        call create_stream_project(self%spproj, cline, string('reference_picking'))
        call self%params%new(cline)
        ! one queue partition per computing unit; set here because parameter derivation
        ! (derive_parallel_settings) resets split_mode to 'even' whatever the command line says
        self%params%split_mode = 'stream'
        call simple_getcwd(self%cwd)
        call self%spproj%read(self%params%projfile)
        if( self%spproj%os_mic%get_noris() /= 0 )then
            THROW_HARD('stream_pick_extract must start from an empty project (eg from root project folder)')
        endif
        allocate(self%set_projects(0))
    end subroutine init_params

    !> The job sets (job and completed folders) and the folders the workers write to.
    subroutine init_job_dirs( self )
        class(stream_stage_refpick), intent(inout) :: self
        call self%sets%new(string(DIR_STREAM), string(DIR_STREAM_COMPLETED), self%params%numlen)
        call simple_mkdir(PATH_HERE//DIR_PICKER)
        call simple_mkdir(PATH_HERE//DIR_EXTRACT)
    end subroutine init_job_dirs

    !> Restart: re-imports the completed sets, or with clear=yes discards them and the picking and
    !! extraction outputs.
    subroutine resume_previous_run( self )
        class(stream_stage_refpick), intent(inout) :: self
        write(logfhandle,'(A)') '>>> RESTARTING EXISTING JOB'
        call del_file(TERM_STREAM)
        if( self%params%clear .eq. 'yes' )then
            call reset_dir(self%sets%get_completed_dir())
            call reset_dir(string(PATH_HERE//DIR_PICKER))
            call reset_dir(string(PATH_HERE//DIR_EXTRACT))
        endif
        call self%import_previous_sets()
    end subroutine resume_previous_run

    !> The queue environment the pick_extract jobs are submitted through.
    subroutine init_queue( self )
        class(stream_stage_refpick), intent(inout) :: self
        character(len=STDLEN) :: partition_env
        integer               :: envlen
        call get_environment_variable(SIMPLE_STREAM_PICK_PARTITION, partition_env, envlen)
        if( envlen > 0 )then
            call self%qenv%new(self%params, 1, stream=.true., qsys_partition=string(trim(partition_env)))
        else
            call self%qenv%new(self%params, 1, stream=.true.)
        endif
    end subroutine init_queue

    !> The worker command line (the stage's own, as pick_extract run in the job folder without a
    !! new output folder) and the make_pickrefs command line (its pixel size is set once known).
    subroutine build_worker_cline( self, cline )
        class(stream_stage_refpick), intent(inout) :: self
        class(cmdline),              intent(in)    :: cline
        self%cline_exec = cline
        call self%cline_exec%set('prg',     'pick_extract')
        call self%cline_exec%set('mkdir',   'no')
        call self%cline_exec%set('dir',     PATH_PARENT)
        call self%cline_exec%set('extract', 'yes')
        if( cline%defined('box_extract') )then
            call self%cline_exec%set('box_extract', max(self%params%box_extract, DEFAULT_EXTRACT_BOX))
        endif
        if( self%cline_exec%defined('dir_exec') ) call self%cline_exec%delete('dir_exec')
        self%cline_pickrefs = cline
        call self%cline_pickrefs%set('prg',          'make_pickrefs')
        call self%cline_pickrefs%set('mkdir',        'no')
        call self%cline_pickrefs%set('stream',       'no')
        call self%cline_pickrefs%set('nrots',        12)
        call self%cline_pickrefs%set('mirr',         'yes')
        call self%cline_pickrefs%set('trust_header', 'yes')
        call self%cline_pickrefs%set('ncls',         10)
        if( self%cline_pickrefs%defined('dir_exec') ) call self%cline_pickrefs%delete('dir_exec')
    end subroutine build_worker_cline

    !> The GUI metadata objects and the pipe ends to the master (-1: none).
    subroutine init_gui( self, fd_read, fd_write )
        class(stream_stage_refpick), intent(inout) :: self
        integer,                     intent(in)    :: fd_read, fd_write
        call self%meta_status%new(GUI_METADATA_STREAM_REFERENCE_PICKING_TYPE)
        call self%meta_micrograph%new(GUI_METADATA_STREAM_REFERENCE_PICKING_MICROGRAPH_TYPE)
        call self%meta_pickrefs%new(GUI_METADATA_STREAM_REFERENCE_PICKING_CLS2D_TYPE)
        call self%pipe%new(fd_read, fd_write, max_metadata_size(), 'reference_picking')
    end subroutine init_gui

    !> One pass: wait for the upstream folder and the picking references, submit the new upstream
    !! projects, import finished sets, report.
    subroutine iterate( self )
        class(stream_stage_refpick), intent(inout) :: self
        integer :: n_imported
        if( .not. self%l_attached )then
            call self%attach_upstream()
            if( .not. self%l_attached )then
                call sleep(self%wait_s)
                return
            endif
        endif
        if( .not. self%pickrefs_available() )then
            call sleep(self%wait_s)
            return
        endif
        call self%submit_new_projects()
        call self%send_status(string('picking and extracting micrographs'))
        call self%schedule_jobs()
        call self%collect_jobs(n_imported)
        if( n_imported > 0 )then
            call self%process_imports()
        else if( .not. self%l_projects_left )then
            call self%idle()
        endif
        call sleep(self%wait_s)
    end subroutine iterate

    !> .true. once the stream is told to stop.
    logical function finished( self )
        class(stream_stage_refpick), intent(in) :: self
        finished = file_exists(TERM_STREAM)
    end function finished

    !> The project with every set's stacks and particles, and removal of the job scripts.
    subroutine finalize( self )
        class(stream_stage_refpick), intent(inout) :: self
        if( self%spproj%os_mic%get_noris() > 0 ) call self%write_project()
        call qsys_cleanup(self%params)
    end subroutine finalize

    subroutine kill( self )
        class(stream_stage_refpick), intent(inout) :: self
        if( .not. self%l_exists ) return
        call self%spproj%kill
        call self%sets%kill
        call self%qenv%kill
        call self%project_buff%kill
        call self%pipe%kill
        call self%cline_exec%kill
        call self%cline_pickrefs%kill
        call self%meta_status%kill
        call self%meta_micrograph%kill
        call self%meta_pickrefs%kill
        if( allocated(self%set_projects)     ) deallocate(self%set_projects)
        if( allocated(self%restored_sources) ) deallocate(self%restored_sources)
        call self%cwd%kill
        self%n_mics_submitted = 0
        self%nptcls_glob      = 0
        self%n_failed_jobs    = 0
        self%nmic_star        = 0
        self%prev_stacksz     = 0
        self%l_attached       = .false.
        self%l_waiting_logged = .false.
        self%l_pickrefs_found = .false.
        self%l_pickrefs_ready = .false.
        self%l_projects_left  = .false.
        self%l_haschanged     = .false.
        self%l_restart        = .false.
        self%l_exists         = .false.
    end subroutine kill

    !---------------- waiting ----------------

    ! Starts watching the upstream completed-projects folder once preprocessing has created it;
    ! a restart's upstream projects go into the watcher history.
    subroutine attach_upstream( self )
        class(stream_stage_refpick), intent(inout) :: self
        type(string) :: completed
        logical      :: l_ready
        integer      :: i
        completed = self%params%dir_target//'/'//DIR_STREAM_COMPLETED
        l_ready   = dir_exists(self%params%dir_target)
        if( l_ready ) l_ready = dir_exists(completed)
        if( .not. l_ready )then
            if( .not. self%l_waiting_logged )then
                write(logfhandle,'(A)') '>>> WAITING FOR '//completed%to_char()//' TO BE GENERATED'
                self%l_waiting_logged = .true.
            endif
            return
        endif
        write(logfhandle,'(A)') '>>> '//completed%to_char()//' FOUND'
        self%project_buff = stream_watcher(self%settle_s, simple_abspath(completed), spproj=.true., nretries=10)
        if( allocated(self%restored_sources) )then
            do i = 1,size(self%restored_sources)
                if( self%restored_sources(i)%strlen() > 0 ) call self%project_buff%add2history(self%restored_sources(i))
            enddo
            deallocate(self%restored_sources)
        endif
        self%l_attached       = .true.
        self%l_waiting_logged = .false. ! for the wait for the picking references
    end subroutine attach_upstream

    ! .true. once the input picking references exist. The initial analysis publishes them with a
    ! rename, so a file that exists is complete.
    logical function pickrefs_available( self )
        class(stream_stage_refpick), intent(inout) :: self
        pickrefs_available = self%l_pickrefs_found
        if( pickrefs_available ) return
        if( .not. file_exists(self%params%pickrefs) )then
            if( .not. self%l_waiting_logged )then
                write(logfhandle,'(A)') '>>> WAITING FOR '//self%params%pickrefs%to_char()
                self%l_waiting_logged = .true.
            endif
            return
        endif
        write(logfhandle,'(A)') '>>> PERFORMING REFERENCE-BASED PICKING'
        self%l_pickrefs_found = .true.
        self%l_waiting_logged = .false.
        pickrefs_available    = .true.
    end function pickrefs_available

    !---------------- submission ----------------

    ! One job set per newly completed upstream project with an accepted micrograph; the picking
    ! references are made with the first one's pixel size.
    subroutine submit_new_projects( self )
        class(stream_stage_refpick), intent(inout) :: self
        type(string), allocatable :: projects(:)
        integer :: nprojects, iproj, nselected, nmics
        self%l_projects_left = .false.
        call self%project_buff%watch(nprojects, projects, max_nmovies=self%params%nparts)
        if( nprojects == 0 ) return
        nmics = 0
        do iproj = 1,nprojects
            call self%create_set_project(projects(iproj), nselected)
            if( nselected > 0 )then
                if( .not. self%l_pickrefs_ready ) call self%prepare_pickrefs()
                call self%sets%submit(self%qenv, self%cline_exec)
                call self%sets%schedule(self%qenv)
                nmics = nmics + nselected
            endif
            call self%project_buff%add2history(projects(iproj))
        enddo
        ! a full watch may have left upstream projects for the next pass
        self%l_projects_left = nprojects == self%params%nparts
        if( nmics > 0 ) write(logfhandle,'(A,I4,A,A)') '>>> ', nmics, ' NEW MICROGRAPHS ADDED; ', cast_time_char(simple_gettime())
    end subroutine submit_new_projects

    ! Writes the accepted micrographs of the upstream project @p upstream_fname as the next job
    ! set, recording the upstream project as its origin; @p nselected is 0 when none is accepted
    ! (and no set is written).
    subroutine create_set_project( self, upstream_fname, nselected )
        class(stream_stage_refpick), intent(inout) :: self
        class(string),               intent(in)    :: upstream_fname
        integer,                     intent(out)   :: nselected
        character(len=*), parameter :: PICK_PREP_KEYS(4) = [character(len=8) :: 'mic_den', 'mic_topo', 'mic_bin', 'mic_diam']
        type(sp_project) :: upstream, set_proj
        integer          :: imic, cnt, ikey
        call upstream%read_segment('mic', upstream_fname)
        nselected = upstream%os_mic%count_state_gt_zero()
        if( nselected == 0 )then
            call upstream%kill
            return
        endif
        set_proj%compenv = self%spproj%compenv
        set_proj%jobproc = self%spproj%jobproc
        call set_proj%os_mic%new(nselected, is_ptcl=.false.)
        cnt = 0
        do imic = 1,upstream%os_mic%get_noris()
            if( upstream%os_mic%get_state(imic) <= 0 ) cycle
            cnt = cnt + 1
            call set_proj%os_mic%transfer_ori(cnt, upstream%os_mic, imic)
            ! outputs of the picking preprocessing that the picker does not use
            do ikey = 1,size(PICK_PREP_KEYS)
                if( set_proj%os_mic%isthere(cnt, trim(PICK_PREP_KEYS(ikey))) )then
                    call set_proj%os_mic%delete_entry(cnt, trim(PICK_PREP_KEYS(ikey)))
                endif
            enddo
        enddo
        ! the pixel size comes from the data, before the picking references are made
        if( .not. self%l_pickrefs_ready ) self%params%smpd = set_proj%os_mic%get(1, 'smpd')
        self%n_mics_submitted = self%n_mics_submitted + nselected
        call self%sets%write_set(set_proj, self%cline_exec, nselected, source=simple_abspath(upstream_fname))
        call set_proj%kill
        call upstream%kill
    end subroutine create_set_project

    ! Makes the picking references at the micrographs' pixel size (make_pickrefs, in this process),
    ! points the workers at them, and sends them to the GUI. The box decided with them (in
    ! moldiam.txt) is the extraction box.
    subroutine prepare_pickrefs( self )
        class(stream_stage_refpick), intent(inout) :: self
        type(commander_make_pickrefs) :: xmake_pickrefs
        type(cmdline)                 :: cline_make_pickrefs
        type(oris)                    :: moldiam
        type(string)                  :: refs_stk, refs_jpg
        integer, allocatable          :: ref_inds(:)
        integer                       :: nrefs, xtiles, ytiles, i
        call self%send_status(string('preparing picking references'))
        cline_make_pickrefs = self%cline_pickrefs
        call cline_make_pickrefs%set('smpd', self%params%smpd)
        call xmake_pickrefs%execute(cline_make_pickrefs)
        call cline_make_pickrefs%kill
        call self%cline_exec%set('pickrefs', PATH_PARENT//PICKREFS_FBODY//self%params%ext%to_char())
        refs_stk = self%cwd//'/'//PICKREFS_FBODY//self%params%ext%to_char()
        refs_jpg = self%cwd//'/'//PICKREFS_FBODY//'.jpeg'
        call moldiam%new(1, is_ptcl=.false.)
        call moldiam%read(string(STREAM_MOLDIAM))
        self%params%box = moldiam%get_int(1, 'box_for_extract')
        call moldiam%kill
        call mrc2jpeg_tiled(refs_stk, refs_jpg, ntiles=nrefs, n_xtiles=xtiles, n_ytiles=ytiles)
        ref_inds = [(i, i=1,nrefs)]
        call send_cavgs(self%pipe, self%meta_pickrefs, refs_jpg, ref_inds, refs_stk, xtiles, ytiles)
        write(logfhandle,'(A)') '>>> PREPARED PICKING TEMPLATES'
        write(logfhandle,'(A)') '>>> JPEG '//refs_jpg%to_char()
        call qsys_cleanup(self%params)
        self%l_pickrefs_ready = .true.
    end subroutine prepare_pickrefs

    subroutine schedule_jobs( self )
        class(stream_stage_refpick), intent(inout) :: self
        integer :: stacksz
        call self%sets%schedule(self%qenv)
        stacksz = self%qenv%qscripts%get_stacksz()
        if( stacksz /= self%prev_stacksz )then
            self%prev_stacksz = stacksz
            write(logfhandle,'(A,I6)') '>>> MICROGRAPHS TO PROCESS:                 ', self%qenv%qscripts%get_stack_range()
        endif
    end subroutine schedule_jobs

    !---------------- import ----------------

    ! Imports the finished sets and counts the failed jobs; @p n_imported is the number of
    ! micrographs imported now.
    subroutine collect_jobs( self, n_imported )
        class(stream_stage_refpick), intent(inout) :: self
        integer,                     intent(out)   :: n_imported
        type(string), allocatable :: done(:)
        integer :: n_failed
        call self%sets%collect(self%qenv, done, n_failed)
        self%n_failed_jobs = self%n_failed_jobs + n_failed
        call self%import_finished_sets(done, n_imported)
    end subroutine collect_jobs

    ! The finished sets @p job_fnames (absolute paths in the job folder): those with micrographs
    ! (pick_extract keeps only micrographs with particles) move to the completed folder and are
    ! imported.
    subroutine import_finished_sets( self, job_fnames, n_imported )
        class(stream_stage_refpick), intent(inout) :: self
        type(string),                intent(in)    :: job_fnames(:)
        integer,                     intent(out)   :: n_imported
        type(sp_project)          :: job
        type(string), allocatable :: completed(:)
        type(string)              :: completed_fname
        integer :: iset
        n_imported = 0
        allocate(completed(0))
        do iset = 1,size(job_fnames)
            call job%read_segment('mic', job_fnames(iset))
            if( job%os_mic%get_noris() > 0 )then
                call self%sets%complete(job_fnames(iset), completed_fname)
                completed = [completed, completed_fname]
            endif
            call job%kill
        enddo
        if( size(completed) > 0 ) call self%import_sets(completed, n_imported)
    end subroutine import_finished_sets

    ! Appends every micrograph of the completed sets @p fnames to the stage's project; the sets
    ! are remembered, in order, for the project write.
    subroutine import_sets( self, fnames, n_imported )
        class(stream_stage_refpick), intent(inout) :: self
        type(string),                intent(in)    :: fnames(:)
        integer,                     intent(out)   :: n_imported
        integer :: n_old, imic
        n_old = self%spproj%os_mic%get_noris()
        call append_mics_from_projects(self%spproj%os_mic, fnames, .false., n_imported)
        do imic = n_old + 1,n_old + n_imported
            self%nptcls_glob = self%nptcls_glob + self%spproj%os_mic%get_int(imic, 'nptcls')
        enddo
        self%set_projects = [self%set_projects, fnames]
    end subroutine import_sets

    ! Restart: the completed sets of the previous run are imported again, and their upstream
    ! projects are kept for the watcher history (attach_upstream).
    subroutine import_previous_sets( self )
        class(stream_stage_refpick), intent(inout) :: self
        type(string), allocatable :: completed(:), sources(:)
        integer :: nmics
        call self%sets%restore(completed, sources)
        if( size(completed) == 0 ) return
        call self%import_sets(completed, nmics)
        self%nmic_star = self%spproj%os_mic%get_noris()
        if( any(sources%strlen() == 0) )then
            write(logfhandle,'(A)') '>>> WARNING: completed sets without an upstream origin (written before job sets'//&
                &' recorded one); their upstream projects will be picked again'
        endif
        call move_alloc(sources, self%restored_sources)
        write(logfhandle,'(A,I6,A)') '>>> IMPORTED ', nmics, ' PREVIOUSLY PROCESSED MICROGRAPHS'
    end subroutine import_previous_sets

    ! After an import: log, GUI status and thumbnails, STAR snapshot.
    subroutine process_imports( self )
        class(stream_stage_refpick), intent(inout) :: self
        integer :: nmics
        nmics = self%spproj%os_mic%get_noris()
        write(logfhandle,'(A,I8)')       '>>> # MICROGRAPHS PROCESSED & IMPORTED  : ', nmics
        write(logfhandle,'(A,I8)')       '>>> # PARTICLES EXTRACTED               : ', self%nptcls_glob
        write(logfhandle,'(A,I3,A2,I3)') '>>> # OF COMPUTING UNITS IN USE/TOTAL   : ', self%qenv%get_navail_computing_units(),&
                                         &'/ ', self%params%nparts
        if( self%n_failed_jobs > 0 ) write(logfhandle,'(A,I8)') '>>> # DESELECTED MICROGRAPHS/FAILED JOBS: ', self%n_failed_jobs
        call self%send_status(string('finding, picking and extracting micrographs'))
        call send_recent_micrographs(self%pipe, self%meta_micrograph, self%spproj%os_mic, NTHUMB_MAX)
        self%last_injection = simple_gettime()
        self%l_haschanged   = .true.
        if( nmics < STAR_EVERY_NMICS )then
            call self%write_mic_star()
        else if( nmics > self%nmic_star + STAR_STEP_NMICS )then
            call self%write_mic_star()
            self%nmic_star = nmics
        endif
    end subroutine process_imports

    ! Nothing imported and nothing waiting: the project and the STAR file after a long inactivity.
    subroutine idle( self )
        class(stream_stage_refpick), intent(inout) :: self
        if( (simple_gettime() - self%last_injection > INACTIVE_TIME) .and. self%l_haschanged )then
            call self%write_project()
            call self%write_mic_star()
            self%l_haschanged = .false.
        endif
    end subroutine idle

    ! The optics groups of the newest map optics assignment has published in optics_dir, applied
    ! by import index to the micrographs, stacks and particles, with its optics segment; .false.
    ! when there is no map yet. Micrographs newer than the map stay in group 1 until a later map
    ! covers them.
    logical function apply_optics_map( self )
        class(stream_stage_refpick), intent(inout) :: self
        apply_optics_map = .false.
        if( self%params%optics_dir%strlen() == 0 ) return
        if( import_latest_optics_map(self%spproj, self%params%optics_dir) == 0 ) return
        apply_optics_map = self%spproj%os_optics%get_noris() > 0
    end function apply_optics_map

    ! The micrographs STAR file, with the optics groups of the newest map, or one group without.
    subroutine write_mic_star( self )
        class(stream_stage_refpick), intent(inout) :: self
        logical :: l_optics
        if( self%spproj%os_mic%get_noris() == 0 ) return
        l_optics = self%apply_optics_map()
        call self%starproj_stream%stream_export_micrographs(self%params, self%spproj, self%params%outdir, optics_set=l_optics)
    end subroutine write_mic_star

    ! Writes the micrographs, and the stacks and particles of every imported set in import order:
    ! one stack per micrograph, particle ranges renumbered, particles pointing at their stack; all
    ! with the optics groups of the newest map, when there is one.
    subroutine write_project( self )
        class(stream_stage_refpick), intent(inout) :: self
        type(sp_project)     :: set_proj
        integer, allocatable :: fromps(:), nptcls_mic(:)
        integer              :: nmics, nptcls, iset, istk, imic, iptcl, i
        logical              :: l_optics
        write(logfhandle,'(A)') '>>> PROJECT UPDATE'
        nmics = self%spproj%os_mic%get_noris()
        ! stacks
        allocate(fromps(nmics), nptcls_mic(nmics), source=0)
        call self%spproj%os_stk%new(nmics, is_ptcl=.false.)
        nptcls = 0
        imic   = 0
        do iset = 1,size(self%set_projects)
            call set_proj%read_segment('stk', self%set_projects(iset))
            do istk = 1,set_proj%os_stk%get_noris()
                imic = imic + 1
                if( imic > nmics ) THROW_HARD('more stacks than micrographs in the imported sets; write_project')
                fromps(imic)     = set_proj%os_stk%get_fromp(istk)
                nptcls_mic(imic) = set_proj%os_stk%get_top(istk) - fromps(imic) + 1
                call self%spproj%os_stk%transfer_ori(imic, set_proj%os_stk, istk)
                call self%spproj%os_stk%set(imic, 'fromp', nptcls + 1)
                call self%spproj%os_stk%set(imic, 'top',   nptcls + nptcls_mic(imic))
                nptcls = nptcls + nptcls_mic(imic)
            enddo
            call set_proj%kill
        enddo
        if( imic /= nmics ) THROW_HARD('fewer stacks than micrographs in the imported sets; write_project')
        ! particles
        call self%spproj%os_ptcl2D%new(nptcls, is_ptcl=.true.)
        iptcl = 0
        imic  = 0
        do iset = 1,size(self%set_projects)
            call set_proj%read_segment('stk',    self%set_projects(iset))
            call set_proj%read_segment('ptcl2D', self%set_projects(iset))
            do istk = 1,set_proj%os_stk%get_noris()
                imic = imic + 1
                do i = fromps(imic),fromps(imic) + nptcls_mic(imic) - 1
                    iptcl = iptcl + 1
                    call self%spproj%os_ptcl2D%transfer_ori(iptcl, set_proj%os_ptcl2D, i)
                    call self%spproj%os_ptcl2D%set_stkind(iptcl, imic)
                enddo
            enddo
            call set_proj%kill
        enddo
        write(logfhandle,'(A,I8)') '>>> # PARTICLES EXTRACTED:          ', self%spproj%os_ptcl2D%get_noris()
        self%spproj%os_ptcl3D = self%spproj%os_ptcl2D
        call self%spproj%os_ptcl3D%delete_2Dclustering
        ! optics groups on every segment, then the segments
        l_optics = self%apply_optics_map()
        call self%spproj%write_segment_inside('mic',    self%params%projfile)
        call self%spproj%write_segment_inside('stk',    self%params%projfile)
        call self%spproj%write_segment_inside('ptcl2D', self%params%projfile)
        call self%spproj%write_segment_inside('ptcl3D', self%params%projfile)
        if( l_optics ) call self%spproj%write_segment_inside('optics', self%params%projfile)
        call self%spproj%os_ptcl3D%kill
        call self%spproj%write_non_data_segments(self%params%projfile)
    end subroutine write_project

    !---------------- GUI ----------------

    subroutine send_status( self, stage )
        class(stream_stage_refpick), intent(inout) :: self
        type(string),                intent(in)    :: stage
        call self%meta_status%set(stage=stage, micrographs_imported=self%n_mics_submitted,&
            &micrographs_accepted=self%spproj%os_mic%get_noris(), particles_extracted=self%nptcls_glob,&
            &box_size=self%params%box)
        call self%pipe%send_meta(self%meta_status)
    end subroutine send_status

    !---------------- helpers ----------------

    ! an empty folder in place of @p dir
    subroutine reset_dir( dir )
        class(string), intent(in) :: dir
        call simple_rmdir(dir)
        call simple_mkdir(dir)
    end subroutine reset_dir

end module simple_stream_stage_refpick
