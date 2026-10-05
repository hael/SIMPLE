!@descr: state and steps of stream task 6 (pool 2D): import sieved particle sets into a growing pool, drive its 2D classification, report to the GUI
!==============================================================================
! MODULE: simple_stream_stage_pool2D
!
! PURPOSE:
!   The body of stream p06 as a type. Each pass takes the particle sets the
!   sieve has handed off since the last pass, adds them to the pool when it is
!   free, lets the pool run its next 2D iteration (or pauses it while too few
!   new particles arrive), answers the GUI (mask diameter, snapshots) and
!   publishes the pool's classified state for 3D after each completed
!   iteration (doc/policies/stream/stream_3D_ingestion_policy.md). The commander
!   (simple_commanders_stream_p06_pool2D) only normalises the command line and loops over
!   iterate() until finished().
!
!   What happens is delegated:
!     - pool iterations, dimensions, stats -> simple_stream_pool2D_utils
!     - snapshots, publications for 3D,    -> simple_stream_refine2D_utils
!       final project
!     - GUI                                -> simple_stream_pipe, simple_stream_gui_senders
!
!   The pool's state is still the module state of simple_stream_pool2D_utils
!   and simple_stream2D_state, so one process makes one stage; the pool reads
!   the stage's command line through master_cline, which kill() releases.
!
! LIFECYCLE:
!   new(cline) -> { iterate() } until finished() -> finalize() -> kill()
!
! RESTART:
!   A restart (the stage's folder exists, or dir_exec is given) cancels a pool
!   iteration a crashed stage left running and removes the previous pool's
!   files from the stage's folder, REFINE2D_FINISHED and the iteration's exit
!   status with them; the pool starts again from
!   every set the sieve has handed off. Snapshots are kept, and the
!   publications for 3D continue their numbering (3D matches the particles of a
!   publication to its own rows by stack and image, so the restarted pool's
!   order does not matter).
!
! JOBS:
!   Each pool iteration is a queued job with an exit status. A job that exits
!   without finishing is retried once; a second failure stops the stage, which
!   then writes its final project from the last complete iteration. On stop the
!   running iteration is cancelled.
!==============================================================================
module simple_stream_stage_pool2D
use simple_defs,                                only: logfhandle, PATH_HERE, COSMSKHALFWIDTH
use simple_defs_fname,                          only: TERM_STREAM, METADATA_EXT, DIR_SNAPSHOT, REFINE2D_FINISHED, JOB_INFO_EXT
use simple_defs_stream,                         only: DIR_STREAM_COMPLETED, SHORTWAIT, WAITTIME, POOL_EXIT_CODE, POOL_INPUT_PROJFILE
use simple_error,                               only: simple_exception
use simple_string,                              only: string
use simple_string_utils,                        only: int2str, int2str_pad, str2int
use simple_fileio,                              only: basename, del_file, file_exists, get_fbody, simple_abspath, simple_getcwd,&
                                                     &simple_list_files_regexp, swap_suffix
use simple_syslib,                              only: dir_exists, simple_mkdir
use simple_cmdline,                             only: cmdline
use simple_parameters,                          only: parameters
use simple_sp_project,                          only: sp_project
use simple_rec_list,                            only: rec_list, rec_iterator, chunk_rec
use simple_stream_watcher,                      only: stream_watcher
use simple_stream_state,                        only: ipc_pipe_pool2D_in, ipc_pipe_pool2D_out
use simple_stream_utils,                        only: create_stream_project
use simple_stream2D_state,                      only: master_cline, last_complete_iter, pool_jpeg_map, pool_jpeg_pop,&
                                                     &pool_jpeg_res
use simple_stream_refine2D_utils,              only: cleanup_root_folder, terminate_stream2D, write_pool_snapshot,&
                                                     &publish_pool_state, delete_pool_publication
use simple_stream_pool2D_utils,                 only: init_pool_clustering, iterate_pool, generate_pool_stats, get_pool_assigned,&
                                                     &get_pool_cavgs_jpeg, get_pool_cavgs_jpeg_ntilesx, get_pool_cavgs_jpeg_ntilesy,&
                                                     &get_pool_cavgs_mrc, get_pool_iter, get_pool_ptr, get_pool_rejected,&
                                                     &get_pool_resolution, is_pool_available, update_mskdiam, update_pool,&
                                                     &update_pool_aln_params, update_pool_status, is_pool_failed, cancel_pool_job
use simple_gui_metadata_utils,                  only: max_metadata_size
use simple_gui_metadata_types,                  only: GUI_METADATA_STREAM_POOL2D_TYPE, GUI_METADATA_STREAM_POOL2D_CLS2D_TYPE,&
                                                     &GUI_METADATA_STREAM_POOL2D_SNAPSHOT_TYPE,&
                                                     &GUI_METADATA_STREAM_POOL2D_SNAPSHOT_CLS2D_TYPE
use simple_gui_metadata_cavg2D,                 only: gui_metadata_cavg2D
use simple_gui_metadata_stream_pool2D,          only: gui_metadata_stream_pool2D
use simple_gui_metadata_stream_pool2D_snapshot, only: gui_metadata_stream_pool2D_snapshot
use simple_gui_metadata_stream_update,          only: gui_metadata_stream_update
use simple_stream_pipe,                         only: stream_pipe
use simple_stream_gui_senders,                  only: send_cavgs
implicit none

public :: stream_stage_pool2D
private
#include "simple_local_flags.inc"

integer, parameter :: OPTICS_ID_DELTA       = 500 ! optics group ids of the STAR files, per GUI display
integer, parameter :: NPTCLS_PER_CLS_MIN    = 20  ! particles per class before the pool runs or resumes
integer, parameter :: EARLY_RATE_FACTOR     = 50  ! micrographs' worth of particles a pause waits for, iterations 2-20
integer, parameter :: LATE_RATE_FACTOR      = 500 ! the same after iteration 20, and before the first iteration
integer, parameter :: LATE_ITER             = 20  ! last iteration of the early pause rule
integer, parameter :: FINAL_ITER            = 25  ! the final sieve set runs the pool uninterrupted to here
integer, parameter :: MSKDIAM_SWITCH_ITER   = 10  ! iteration from which the sieve's mask diameter applies
integer, parameter :: EXPORT_START_ITER     = 25  ! the pool is published for 3D after each iteration from this one on
                                                   ! (the final run's last iteration, FINAL_ITER, among them)
integer, parameter :: NPUBLICATIONS_KEPT    = 2   ! the newest publications kept on disk; older ones are removed

! Components and steps are public so simple_stream_stage_pool2D_tester can assemble a stage and
! run one step at a time; production code uses new/iterate/finished/finalize/kill.
type :: stream_stage_pool2D
    type(parameters), allocatable             :: params
    type(cmdline), pointer                    :: cline => null() ! the stage's command line, read by the pool (master_cline)
    type(sp_project)                          :: spproj          ! the stage's project
    type(stream_watcher)                      :: project_buff    ! the sets the sieve hands off
    type(rec_list)                            :: setslist        ! one record per set; included once in the pool
    type(stream_pipe)                         :: pipe            ! to and from the master
    type(gui_metadata_stream_pool2D)          :: meta_status
    type(gui_metadata_cavg2D)                 :: meta_cavgs
    type(gui_metadata_stream_pool2D_snapshot) :: meta_snapshot
    type(gui_metadata_cavg2D)                 :: meta_snapshot_cavgs
    ! the latest snapshot, as write_pool_snapshot reports it
    type(string)              :: snapshot_dir, snapshot_filename, snapshot_jpeg, snapshot_mrc
    integer,      allocatable :: snapshot_idx(:), snapshot_pop(:)
    real,         allocatable :: snapshot_res(:)
    integer :: snapshot_ntilesx          = 0
    integer :: snapshot_ntilesy          = 0
    integer :: snapshot_nptcls           = 0 ! selected particles of the last snapshot; 0 when it could not be written
    integer :: last_snapshot_id          = 0
    ! particle counts and the pause rule
    integer :: nptcls_glob               = 0  ! selected particles imported
    integer :: nptcls_glob_state_1       = 0  ! particles of the pool with state > 0
    integer :: nmics                     = 0  ! micrographs of the pool
    integer :: state_1_particle_rate     = 0  ! selected particles per micrograph
    integer :: nptcls_threshold          = 0  ! particles before the first iteration runs
    integer :: nptcls_dynamic_threshold  = 0  ! particles a paused pool waits for
    integer :: iter_last_import          = -1 ! pool iteration of the last import
    integer :: last_sent_iter            = 0  ! iteration whose class averages the GUI has
    integer :: last_export_iteration     = EXPORT_START_ITER - 1
    integer :: last_export_id            = 1
    integer :: optics_id_offset          = 0
    real    :: final_mskdiam             = 0. ! the sieve's mask diameter (A), applied from MSKDIAM_SWITCH_ITER
    real    :: mskdiam                   = 0. ! the pool's mask diameter (A): given, the default, or updated
    real    :: smpd                      = 0. ! the pool's native pixel size (A), from its first import
    integer :: box                       = 0  ! the pool's native box (px), from its first import
    logical :: l_pause                   = .false.
    logical :: l_sieve_final             = .false. ! the sieve's final set is in the pool
    logical :: l_stepwise                = .false. ! import only enough sets to reach the threshold
    logical :: l_mskdiam_read            = .false. ! final_mskdiam has been read from the first set
    logical :: l_restart                 = .false.
    logical :: l_pool_started            = .false. ! init_pool_clustering has run (on the first import)
    logical :: l_attached                = .false. ! the sieve's completed folder exists and is watched
    logical :: l_waiting_logged          = .false.
    logical :: l_exists                  = .false.
    ! waits (s); tests set them to 0, and settle_s to -1 to take files written in the same second
    integer :: settle_s                  = SHORTWAIT ! a set is taken once untouched longer than this;
                                                     ! the sieve hands sets off with a rename
    integer :: wait_s                    = WAITTIME  ! pause at the end of a pass
contains
    procedure :: new
    procedure :: iterate
    procedure :: finished
    procedure :: finalize
    procedure :: kill
    ! the steps of new(), in order
    procedure :: init_params
    procedure :: clean_previous_run
    procedure :: restore_export_id
    procedure :: init_gui
    ! the steps of iterate() and their helpers
    procedure :: attach_upstream
    procedure :: watch_sets
    procedure :: read_final_mskdiam
    procedure :: update_pool_progress
    procedure :: import_sets
    procedure :: transfer_sets
    procedure :: start_pool
    procedure :: apply_pause_policy
    procedure :: unpause
    procedure :: run_iteration
    procedure :: apply_final_mskdiam
    procedure :: set_mskdiam
    procedure :: apply_gui_updates
    procedure :: write_snapshot
    procedure :: export_pool_state
    ! GUI
    procedure :: send_status
    procedure :: send_pool_cavgs
    procedure :: send_snapshot
    ! rules that use no stage state
    procedure, nopass :: pause_rate_factor
    procedure, nopass :: target_nptcls
    procedure, nopass :: runs_to_final
    procedure, nopass :: default_mskdiam
end type stream_stage_pool2D

contains

    !---------------- lifecycle ----------------

    !> Builds the parameters from the normalised @p cline, clears a previous run and prepares the GUI.
    subroutine new( self, cline )
        class(stream_stage_pool2D), intent(inout) :: self
        class(cmdline),             intent(inout) :: cline
        call self%kill()
        call self%init_params(cline)
        self%l_exists = .true. ! from here kill() releases what has been built
        if( self%l_restart ) call self%clean_previous_run()
        call self%restore_export_id()
        call self%init_gui(ipc_pipe_pool2D_out(1), ipc_pipe_pool2D_in(2))
        call self%send_status(string('initialising'))
    end subroutine new

    !> The stage's project file (made and given a computing environment), its parameters and its
    !! copy of the command line; the project must start without micrographs.
    subroutine init_params( self, cline )
        class(stream_stage_pool2D), intent(inout) :: self
        class(cmdline),             intent(inout) :: cline
        type(string) :: outdir, dir_exec
        ! a restart is recognised by its output directory (or dir_exec), before params%new makes one
        self%l_restart = .false.
        if( cline%defined('outdir') )then
            outdir = cline%get_carg('outdir')
            if( outdir%strlen() > 0 ) self%l_restart = dir_exists(outdir)
        endif
        if( cline%defined('dir_exec') )then
            dir_exec = cline%get_carg('dir_exec')
            if( .not. file_exists(dir_exec) ) THROW_HARD('Previous directory does not exist: '//dir_exec%to_char())
            self%l_restart = .true.
        endif
        call create_stream_project(self%spproj, cline, string('pool2D'))
        if( .not. allocated(self%params) ) allocate(self%params)
        call self%params%new(cline)
        self%l_stepwise = self%params%stepwise == 'yes'
        self%mskdiam    = self%params%mskdiam
        allocate(self%cline)
        self%cline = cline
        call self%cline%set('mkdir', 'no')
        if( self%cline%defined('dir_exec') ) call self%cline%delete('dir_exec')
        call self%spproj%read(self%params%projfile)
        if( self%spproj%os_mic%get_noris() /= 0 )then
            THROW_HARD('commander_stream_p06_pool2D must start from an empty project (e.g. from root project folder)')
        endif
        call simple_mkdir(PATH_HERE//DIR_STREAM_COMPLETED)
        self%optics_id_offset = max(self%params%nicedispid - 1, 0) * OPTICS_ID_DELTA
        write(logfhandle,'(A,I8)') '>>> OPTICS ID OFFSET', self%optics_id_offset
    end subroutine init_params

    !> Restart: the previous pool's files, in the stage's folder (params%new has moved there; for
    !! dir_exec too, where the stage it replaces cleaned the folder it was launched from).
    subroutine clean_previous_run( self )
        class(stream_stage_pool2D), intent(inout) :: self
        type(string) :: cwd
        call simple_getcwd(cwd)
        write(logfhandle,'(A,A)') '>>> RESTARTING EXISTING JOB ', cwd%to_char()
        ! an iteration a crashed stage left running would write into the new pool's files
        call cancel_pool_job()
        call cleanup_root_folder() ! with TERM_STREAM
        ! the previous iteration's completion and status would pass for the new pool's first
        call del_file(REFINE2D_FINISHED)
        call del_file(POOL_EXIT_CODE)
        call del_file(POOL_EXIT_CODE//JOB_INFO_EXT)
        call del_file(POOL_INPUT_PROJFILE)
    end subroutine clean_previous_run

    !> The next publication number, one past the highest in the completed folder: after a restart
    !! the sequence continues instead of overwriting a publication 3D may be reading.
    subroutine restore_export_id( self )
        class(stream_stage_pool2D), intent(inout) :: self
        type(string), allocatable :: exports(:)
        type(string) :: fbody
        integer      :: i, id, iostat
        self%last_export_id = 1
        call simple_list_files_regexp(string(PATH_HERE//DIR_STREAM_COMPLETED), '\.simple$', exports)
        if( .not. allocated(exports) ) return
        do i = 1,size(exports)
            fbody = basename(exports(i))
            fbody = get_fbody(fbody, METADATA_EXT, separator=.false.)
            id    = str2int(fbody, iostat)
            if( iostat == 0 ) self%last_export_id = max(self%last_export_id, id + 1)
        enddo
        if( self%last_export_id > 1 ) write(logfhandle,'(A,I6)') '>>> PUBLICATIONS FOR 3D CONTINUE AT ', self%last_export_id
    end subroutine restore_export_id

    !> The GUI metadata objects and the pipe ends to and from the master (-1: none).
    subroutine init_gui( self, fd_read, fd_write )
        class(stream_stage_pool2D), intent(inout) :: self
        integer,                    intent(in)    :: fd_read, fd_write
        call self%meta_status%new(GUI_METADATA_STREAM_POOL2D_TYPE)
        call self%meta_cavgs%new(GUI_METADATA_STREAM_POOL2D_CLS2D_TYPE)
        call self%meta_snapshot%new(GUI_METADATA_STREAM_POOL2D_SNAPSHOT_TYPE)
        call self%meta_snapshot_cavgs%new(GUI_METADATA_STREAM_POOL2D_SNAPSHOT_CLS2D_TYPE)
        call self%pipe%new(fd_read, fd_write, max_metadata_size(), 'pool2D')
    end subroutine init_gui

    !> One pass: wait for the sieve's folder; take new sets; advance, feed or pause the pool; answer
    !! the GUI; export for 3D.
    subroutine iterate( self )
        class(stream_stage_pool2D), intent(inout) :: self
        if( .not. self%l_attached )then
            call self%attach_upstream()
            if( .not. self%l_attached )then
                call self%send_status(string('waiting on particle sieving'))
                call self%apply_gui_updates()
                call sleep(self%wait_s)
                return
            endif
        endif
        call self%watch_sets()
        call self%update_pool_progress()
        ! an iteration's results are back and the next is not dispatched: one consistent state
        call self%export_pool_state()
        call self%import_sets()
        call self%apply_pause_policy()
        call self%run_iteration()
        call self%apply_final_mskdiam()
        call self%send_pool_cavgs()
        call self%apply_gui_updates()
        call sleep(self%wait_s)
    end subroutine iterate

    !> .true. once the stream is told to stop, or once a pool iteration has failed twice.
    logical function finished( self )
        class(stream_stage_pool2D), intent(in) :: self
        finished = file_exists(TERM_STREAM) .or. is_pool_failed()
    end function finished

    !> The last status, the running iteration cancelled, then the final project from the last
    !! complete iteration.
    subroutine finalize( self )
        class(stream_stage_pool2D), intent(inout) :: self
        call self%meta_status%set_user_input(.false.)
        if( is_pool_failed() )then
            call self%send_status(string('stopped: pool iteration '//int2str(get_pool_iter())//' failed twice (see its log)'))
        else
            call self%send_status(string('terminating'))
        endif
        if( self%l_pool_started )then
            call cancel_pool_job()
            call terminate_stream2D(self%params, optics_dir=self%params%optics_dir)
        endif
    end subroutine finalize

    subroutine kill( self )
        class(stream_stage_pool2D), intent(inout) :: self
        if( .not. self%l_exists )then
            if( allocated(self%params) ) deallocate(self%params)
            return
        endif
        if( associated(self%cline) )then
            if( associated(master_cline, self%cline) ) nullify(master_cline)
            call self%cline%kill
            deallocate(self%cline)
        endif
        nullify(self%cline)
        call self%spproj%kill
        call self%project_buff%kill
        call self%setslist%kill
        call self%pipe%kill
        call self%meta_status%kill
        call self%meta_cavgs%kill
        call self%meta_snapshot%kill
        call self%meta_snapshot_cavgs%kill
        call self%snapshot_dir%kill
        call self%snapshot_filename%kill
        call self%snapshot_jpeg%kill
        call self%snapshot_mrc%kill
        if( allocated(self%snapshot_idx) ) deallocate(self%snapshot_idx)
        if( allocated(self%snapshot_pop) ) deallocate(self%snapshot_pop)
        if( allocated(self%snapshot_res) ) deallocate(self%snapshot_res)
        if( allocated(self%params) ) deallocate(self%params)
        self%snapshot_ntilesx         = 0
        self%snapshot_ntilesy         = 0
        self%snapshot_nptcls          = 0
        self%last_snapshot_id         = 0
        self%nptcls_glob              = 0
        self%nptcls_glob_state_1      = 0
        self%nmics                    = 0
        self%state_1_particle_rate    = 0
        self%nptcls_threshold         = 0
        self%nptcls_dynamic_threshold = 0
        self%iter_last_import         = -1
        self%last_sent_iter           = 0
        self%last_export_iteration    = EXPORT_START_ITER - 1
        self%last_export_id           = 1
        self%optics_id_offset         = 0
        self%final_mskdiam            = 0.
        self%mskdiam                  = 0.
        self%smpd                     = 0.
        self%box                      = 0
        self%l_pause                  = .false.
        self%l_sieve_final            = .false.
        self%l_stepwise               = .false.
        self%l_mskdiam_read           = .false.
        self%l_restart                = .false.
        self%l_pool_started           = .false.
        self%l_attached               = .false.
        self%l_waiting_logged         = .false.
        self%l_exists                 = .false.
    end subroutine kill

    !---------------- steps ----------------

    ! Starts watching the sieve's completed folder once the sieve has created it.
    subroutine attach_upstream( self )
        class(stream_stage_pool2D), intent(inout) :: self
        type(string) :: completed
        logical      :: l_ready
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
        self%project_buff     = stream_watcher(self%settle_s, simple_abspath(completed), spproj=.true., nretries=10)
        self%l_attached       = .true.
        self%l_waiting_logged = .false.
    end subroutine attach_upstream

    ! One record per newly handed-off set; the first set gives the sieve's mask diameter.
    subroutine watch_sets( self )
        class(stream_stage_pool2D), intent(inout) :: self
        type(string), allocatable :: projects(:)
        integer :: nprojects, i
        call self%project_buff%watch(nprojects, projects)
        if( nprojects == 0 ) return
        call self%project_buff%add2history(projects)
        do i = 1,nprojects
            call self%setslist%push2chunk_list(projects(i), self%setslist%size() + 1, .true.)
        enddo
        if( .not. self%l_mskdiam_read ) call self%read_final_mskdiam(projects(1))
    end subroutine watch_sets

    ! The mask diameter the sieve's 2D used, from the class averages of set @p projfile; a set
    ! without class averages (the sieve's empty final set) leaves it for a later set.
    subroutine read_final_mskdiam( self, projfile )
        class(stream_stage_pool2D), intent(inout) :: self
        class(string),              intent(in)    :: projfile
        type(sp_project) :: spproj
        type(string)     :: stk
        real             :: smpd
        integer          :: ncls
        call spproj%read_segment('out', projfile)
        call spproj%get_cavgs_stk(stk, ncls, smpd, fail=.false.)
        if( ncls > 0 ) call spproj%get_mskdiam('cavg', self%final_mskdiam)
        call spproj%kill
        if( ncls <= 0 ) return
        self%l_mskdiam_read = .true.
        write(logfhandle,'(A,F8.2)') '>>> FINAL MASK DIAMETER SET TO : ', self%final_mskdiam
    end subroutine read_final_mskdiam

    ! A paused pool refreshes its statistics; a running one is checked for a completed iteration,
    ! whose particle parameters, classes and dimensions come back to the pool.
    subroutine update_pool_progress( self )
        class(stream_stage_pool2D), intent(inout) :: self
        if( .not. self%l_pool_started ) return
        if( self%l_pause )then
            call generate_pool_stats(self%params)
        else
            call update_pool_status(self%params)
            call update_pool(self%params)
        endif
    end subroutine update_pool_progress

    ! The new sets go into the pool when it is free (or not started); the first import starts it.
    ! An import resumes a paused pool once it brings the particles the pause waits for.
    subroutine import_sets( self )
        class(stream_stage_pool2D), intent(inout) :: self
        class(sp_project), pointer :: pool
        integer :: nimported
        if( self%setslist%size() == 0 ) return
        if( self%l_pool_started .and. .not. is_pool_available() ) return
        call get_pool_ptr(pool)
        call self%transfer_sets(pool, nimported)
        if( nimported > 0 )then
            if( .not. self%l_pool_started ) call self%start_pool(pool)
            self%iter_last_import = get_pool_iter()
            if( self%nptcls_glob_state_1 > self%nptcls_dynamic_threshold .or. self%l_sieve_final ) call self%unpause()
        endif
        nullify(pool)
    end subroutine import_sets

    ! Appends the sets not yet included to @p pool: their micrographs and stacks, and their
    ! particles as new ones (no 2D parameters but their shifts). With stepwise=yes only enough sets
    ! are taken for this import's particles to reach the particle threshold; the others wait for a
    ! later import. The sieve's
    ! final set (sieve_final=yes) is noted, and a later set with particles that is not final takes
    ! the note back (the sieve had more particles after all). The sieve's final set may hold no
    ! particles: it only ends the intake.
    subroutine transfer_sets( self, pool, nimported )
        class(stream_stage_pool2D), intent(inout) :: self
        class(sp_project),          intent(inout) :: pool
        integer,                    intent(out)   :: nimported
        type(sp_project), allocatable :: sets(:)
        logical,          allocatable :: l_included(:)
        type(rec_iterator) :: it
        type(chunk_rec)    :: crec
        integer :: nsets, iset, irec, nmics_new, nptcls_new, nsel, nsel_tot, nsel_cum, target_sel
        integer :: pool_nmics, pool_nptcls, imic, jmic, fromp, ind, nptcls, i, iptcl, jptcl
        nimported  = 0
        l_included = self%setslist%get_included_flags()
        if( count(.not. l_included) == 0 ) return
        target_sel = huge(target_sel)
        if( self%l_stepwise )then
            target_sel = max(self%params%ncls * NPTCLS_PER_CLS_MIN, 1)
            if( self%nptcls_threshold > 0 ) target_sel = self%nptcls_threshold
        endif
        ! read the sets in list order
        allocate(sets(count(.not. l_included)))
        nsets      = 0
        nmics_new  = 0
        nptcls_new = 0
        nsel_cum   = 0 ! this import's particles only: one set per import otherwise, once past the threshold
        it         = self%setslist%begin()
        do irec = 1,self%setslist%size()
            call it%get(crec)
            call it%next()
            if( crec%included ) cycle
            nsets = nsets + 1
            call sets(nsets)%read_segment('mic',    crec%projfile)
            call sets(nsets)%read_segment('stk',    crec%projfile)
            call sets(nsets)%read_segment('ptcl2D', crec%projfile)
            call sets(nsets)%read_segment('out',    crec%projfile)
            if( is_final_set(sets(nsets)) )then
                self%l_sieve_final = .true.
                write(logfhandle,'(A,I3)') '>>> FINAL SIEVE SET DETECTED - RUNNING UNINTERRUPTED TO ITERATION ', FINAL_ITER
            else if( self%l_sieve_final .and. sets(nsets)%os_ptcl2D%get_noris() > 0 )then
                self%l_sieve_final = .false.
                write(logfhandle,'(A)') '>>> A SIEVE SET AFTER THE FINAL ONE: THE POOL IS NO LONGER FINAL'
            endif
            nmics_new  = nmics_new  + sets(nsets)%os_mic%get_noris()
            nptcls_new = nptcls_new + sets(nsets)%os_ptcl2D%get_noris()
            nsel_cum   = nsel_cum   + sets(nsets)%os_ptcl2D%get_noris(consider_state=.true.)
            if( self%l_stepwise .and. nsel_cum >= target_sel )then
                write(logfhandle,'(A,I8)') '>>> STEPWISE IMPORT: DEFERRING REMAINING SETS UNTIL POOL PAUSES AGAIN, REACHED ', nsel_cum
                exit
            endif
        enddo
        ! only the sieve's empty final set: noted, nothing to transfer
        if( nmics_new == 0 .and. nptcls_new == 0 )then
            ! the sets read are the first nsets records not yet included
            iset = 0
            it   = self%setslist%begin()
            do irec = 1,self%setslist%size()
                if( iset == nsets ) exit
                call it%get(crec)
                if( .not. crec%included )then
                    iset = iset + 1
                    crec%included = .true.
                    call self%setslist%replace_iterator(it, crec)
                endif
                call it%next()
            enddo
            do iset = 1,nsets
                call sets(iset)%kill
            enddo
            deallocate(sets, l_included)
            return
        endif
        ! room in the pool
        pool_nmics  = pool%os_mic%get_noris()
        pool_nptcls = pool%os_ptcl2D%get_noris()
        if( pool_nmics == 0 )then
            call pool%os_mic%new(nmics_new,     is_ptcl=.false.)
            call pool%os_stk%new(nmics_new,     is_ptcl=.false.)
            call pool%os_ptcl2D%new(nptcls_new, is_ptcl=.true.)
            fromp = 1
        else
            call pool%os_mic%reallocate(pool_nmics + nmics_new)
            call pool%os_stk%reallocate(pool_nmics + nmics_new)
            call pool%os_ptcl2D%reallocate(pool_nptcls + nptcls_new)
            fromp = pool%os_stk%get_top(pool_nmics) + 1
        endif
        ! transfer, the k-th set read being the k-th record not yet included
        imic     = pool_nmics
        nsel_tot = 0
        iset     = 0
        it       = self%setslist%begin()
        do irec = 1,self%setslist%size()
            if( iset == nsets ) exit
            call it%get(crec)
            if( crec%included )then
                call it%next()
                cycle
            endif
            iset = iset + 1
            ind  = 1
            do jmic = 1,sets(iset)%os_mic%get_noris()
                imic = imic + 1
                call pool%os_mic%transfer_ori(imic, sets(iset)%os_mic, jmic)
                call pool%os_stk%transfer_ori(imic, sets(iset)%os_stk, jmic)
                nptcls = sets(iset)%os_stk%get_int(jmic, 'nptcls')
                call pool%os_stk%set(imic, 'fromp', fromp)
                call pool%os_stk%set(imic, 'top',   fromp + nptcls - 1)
                !$omp parallel do private(i,iptcl,jptcl) default(shared) proc_bind(close)
                do i = 1,nptcls
                    iptcl = fromp + i - 1
                    jptcl = ind   + i - 1
                    call pool%os_ptcl2D%transfer_ori(iptcl, sets(iset)%os_ptcl2D, jptcl)
                    call pool%os_ptcl2D%set_stkind(iptcl, imic)
                    call pool%os_ptcl2D%set(iptcl, 'updatecnt', 0)
                    call pool%os_ptcl2D%set(iptcl, 'frac',      0.)
                    call pool%os_ptcl2D%delete_2Dclustering(iptcl, keepshifts=.true.)
                enddo
                !$omp end parallel do
                ind   = ind   + nptcls
                fromp = fromp + nptcls
            enddo
            nsel     = sets(iset)%os_ptcl2D%get_noris(consider_state=.true.)
            nsel_tot = nsel_tot + nsel
            write(logfhandle,'(A,I6,A,I6)') '>>> TRANSFERRED ', nsel, ' PARTICLES FROM SET ', crec%id
            crec%included = .true.
            call self%setslist%replace_iterator(it, crec)
            call it%next()
        enddo
        nimported = nsets
        ! counts
        self%nptcls_glob           = self%nptcls_glob + nsel_tot
        self%nptcls_glob_state_1   = pool%os_ptcl2D%count_state_gt_zero()
        self%nmics                 = pool%os_mic%get_noris()
        self%state_1_particle_rate = ceiling(real(self%nptcls_glob_state_1) / real(max(1, self%nmics)))
        do iset = 1,size(sets)
            call sets(iset)%kill
        enddo
        deallocate(sets, l_included)
    end subroutine transfer_sets

    ! .true. for the sieve's final set (sieve_final=yes in its out segment)
    logical function is_final_set( set )
        class(sp_project), intent(inout) :: set
        is_final_set = .false.
        if( set%os_out%get_noris() < 1 ) return
        if( .not. set%os_out%isthere(1, 'sieve_final') ) return
        is_final_set = set%os_out%get_str(1, 'sieve_final') == 'yes'
    end function is_final_set

    ! The first import: the pool's sampling and box (from the data), a mask diameter when none was
    ! given, and the pool module, which keeps them and a pointer to the stage's command line.
    subroutine start_pool( self, pool )
        class(stream_stage_pool2D), intent(inout) :: self
        class(sp_project),          intent(inout) :: pool
        self%smpd = pool%get_smpd()
        self%box  = pool%get_box()
        if( self%mskdiam <= 0. )then
            self%mskdiam = default_mskdiam(self%box, self%smpd)
            write(logfhandle,'(A,F8.2)') '>>> INITIAL MASK DIAMETER SET TO', self%mskdiam
        endif
        call init_pool_clustering(self%params, self%cline, self%spproj, self%box, self%smpd, self%mskdiam)
        self%l_pool_started = .true.
    end subroutine start_pool

    ! From the sieve's final set to FINAL_ITER the pool runs uninterrupted; otherwise a running pool
    ! pauses once imports stop (pause_rate_factor), until it has the particles target_nptcls gives.
    subroutine apply_pause_policy( self )
        class(stream_stage_pool2D), intent(inout) :: self
        integer :: iter, factor
        iter = get_pool_iter()
        if( runs_to_final(iter, self%l_sieve_final) )then
            call self%unpause()
            return
        endif
        if( self%l_pause ) return
        factor = pause_rate_factor(iter, self%iter_last_import)
        if( factor == 0 ) return
        self%l_pause = is_pool_available()
        if( self%l_pause )then
            self%nptcls_dynamic_threshold = target_nptcls(self%nptcls_glob_state_1, self%params%ncls,&
                &self%state_1_particle_rate, factor)
            write(logfhandle,'(A,I8)') '>>> PAUSING 2D ANALYSIS UNTIL #PTCLS IS ', self%nptcls_dynamic_threshold
        endif
    end subroutine apply_pause_policy

    ! Resumes with new sets or updated parameters.
    subroutine unpause( self )
        class(stream_stage_pool2D), intent(inout) :: self
        if( self%l_pause ) write(logfhandle,'(A)') '>>> RESUMING 2D ANALYSIS'
        self%l_pause = .false.
    end subroutine unpause

    ! Starts the next pool iteration unless the pool is paused or, before the first, has too few
    ! particles (the threshold follows the particle rate until then).
    subroutine run_iteration( self )
        class(stream_stage_pool2D), intent(inout) :: self
        integer :: threshold
        if( get_pool_iter() == 0 )then
            threshold = target_nptcls(0, self%params%ncls, self%state_1_particle_rate, LATE_RATE_FACTOR)
            if( threshold /= self%nptcls_threshold ) write(logfhandle,'(A,I8)') '>>> INITIAL PARTICLE THRESHOLD: ', threshold
            self%nptcls_threshold = threshold
        endif
        if( self%l_pause )then
            call self%send_status(string('paused whilst awaiting new particles'))
        else if( self%nptcls_glob == 0 .or. self%nptcls_threshold == 0 )then
            call self%send_status(string('waiting for sieved particles'))
        else if( self%nptcls_glob_state_1 < self%nptcls_threshold .and. .not. self%l_sieve_final )then
            call self%send_status(string('waiting for minimum number sieved particles ... '//&
                &int2str(ceiling(100. * real(self%nptcls_glob_state_1) / real(self%nptcls_threshold)))//'%'))
        else
            call update_pool_aln_params()
            call iterate_pool(self%params)
            call self%send_status(string('finding and classifying particles'))
        endif
    end subroutine run_iteration

    ! A new mask diameter (A): the stage's, and the running pool's (update_mskdiam).
    subroutine set_mskdiam( self, mskdiam )
        class(stream_stage_pool2D), intent(inout) :: self
        integer,                    intent(in)    :: mskdiam
        self%mskdiam = real(mskdiam)
        if( self%l_pool_started )then
            call update_mskdiam(mskdiam)
        else
            write(logfhandle,'(A,I4,A)') '>>> MASK DIAMETER SET TO', mskdiam, ' A FOR THE POOL TO COME'
        endif
    end subroutine set_mskdiam

    ! From MSKDIAM_SWITCH_ITER the pool uses the mask diameter of the sieve's 2D, once.
    subroutine apply_final_mskdiam( self )
        class(stream_stage_pool2D), intent(inout) :: self
        if( get_pool_iter() < MSKDIAM_SWITCH_ITER .or. self%final_mskdiam <= 0. ) return
        call self%set_mskdiam(nint(self%final_mskdiam))
        self%final_mskdiam = 0.
    end subroutine apply_final_mskdiam

    ! Drains the GUI updates. A new mask diameter applies from the next iteration, resumes a paused
    ! pool and replaces the sieve's pending one; a snapshot request (once the pool runs) is written
    ! and sent back.
    subroutine apply_gui_updates( self )
        class(stream_stage_pool2D), intent(inout) :: self
        type(gui_metadata_stream_update) :: update
        character(len=:), allocatable    :: buffer
        integer :: mskdiam
        do while( self%pipe%receive(buffer) )
            update  = transfer(buffer, update)
            mskdiam = nint(update%get_mskdiam2D_update())
            if( mskdiam > 0 .and. mskdiam /= nint(self%mskdiam) )then
                call self%set_mskdiam(mskdiam)
                self%final_mskdiam = 0.
                call self%unpause()
            endif
            if( update%has_snapshot2D_update() .and. self%l_pool_started ) call self%write_snapshot(update)
        enddo
    end subroutine apply_gui_updates

    ! A snapshot the GUI asked for: the selected classes of an iteration as a project with STAR
    ! files and the newest optics map, in snapshots/<name>/; each request is written once. An
    ! iteration the pool no longer keeps (it keeps the last POOL_NHISTORY) is not written, and is
    ! reported with no particles and no file.
    subroutine write_snapshot( self, update )
        class(stream_stage_pool2D),       intent(inout) :: self
        type(gui_metadata_stream_update), intent(in)    :: update
        integer, allocatable :: selection(:)
        type(string) :: cwd, stem
        integer      :: snapshot_id, iteration
        call update%get_snapshot2D_update(snapshot_id, iteration, selection, self%snapshot_filename)
        if( snapshot_id <= self%last_snapshot_id ) return
        call simple_getcwd(cwd)
        stem              = swap_suffix(self%snapshot_filename, '', METADATA_EXT)
        self%snapshot_dir = cwd//'/'//DIR_SNAPSHOT//stem
        call write_pool_snapshot(iteration, selection, self%snapshot_dir//'/'//self%snapshot_filename,&
            &self%snapshot_dir//'/'//stem, self%params%optics_dir, self%optics_id_offset, self%snapshot_nptcls,&
            &self%snapshot_jpeg, self%snapshot_mrc, self%snapshot_ntilesx, self%snapshot_ntilesy,&
            &self%snapshot_idx, self%snapshot_pop, self%snapshot_res)
        self%last_snapshot_id = snapshot_id
        call self%send_snapshot()
    end subroutine write_snapshot

    ! Once iteration EXPORT_START_ITER or a later one has come back, and before the next is dispatched, the
    ! pool's classified state is published as the next project of the completed folder (stacks
    ! whose particles have been through an iteration; never particles just imported). The newest
    ! NPUBLICATIONS_KEPT publications are kept. One publication per iteration.
    subroutine export_pool_state( self )
        class(stream_stage_pool2D), intent(inout) :: self
        type(string) :: cwd
        integer      :: nstks
        if( .not. self%l_pool_started ) return
        if( .not. is_pool_available() ) return
        if( get_pool_iter() <= self%last_export_iteration ) return
        call simple_getcwd(cwd)
        call publish_pool_state(publication_fname(cwd, self%last_export_id), nstks)
        if( nstks > 0 )then
            if( self%last_export_id > NPUBLICATIONS_KEPT )then
                call delete_pool_publication(publication_fname(cwd, self%last_export_id - NPUBLICATIONS_KEPT))
            endif
            self%last_export_id = self%last_export_id + 1
        endif
        self%last_export_iteration = get_pool_iter()
    end subroutine export_pool_state

    !---------------- GUI ----------------

    ! Iteration, particle counts, mask and resolution.
    subroutine send_status( self, stage )
        class(stream_stage_pool2D), intent(inout) :: self
        type(string),               intent(in)    :: stage
        call self%meta_status%set(stage=stage, iteration=last_complete_iter, particles_imported=self%nptcls_glob,&
            &particles_accepted=get_pool_assigned(), particles_rejected=get_pool_rejected(),&
            &mskdiam=nint(self%mskdiam), mskscale=real(self%box)*self%smpd,&
            &resolution=get_pool_resolution())
        call self%pipe%send_meta(self%meta_status)
    end subroutine send_status

    ! The class averages of each newly completed iteration (from iteration 2), as sprite-sheet
    ! tiles with their resolutions and populations; the GUI may give input from then on.
    subroutine send_pool_cavgs( self )
        class(stream_stage_pool2D), intent(inout) :: self
        type(string) :: jpg
        if( get_pool_iter() <= 1 ) return
        call self%meta_status%set_user_input(.true.)
        if( self%last_sent_iter == last_complete_iter ) return
        if( allocated(pool_jpeg_map) )then
            jpg = get_pool_cavgs_jpeg()
            if( jpg%strlen() > 0 ) call send_cavgs(self%pipe, self%meta_cavgs, jpg, pool_jpeg_map, get_pool_cavgs_mrc(),&
                &get_pool_cavgs_jpeg_ntilesx(), get_pool_cavgs_jpeg_ntilesy(), res=pool_jpeg_res, pop=pool_jpeg_pop)
        endif
        self%last_sent_iter = last_complete_iter
    end subroutine send_pool_cavgs

    ! The latest snapshot: its project and particle count, then its selected classes as tiles.
    subroutine send_snapshot( self )
        class(stream_stage_pool2D), intent(inout) :: self
        ! a snapshot that could not be written goes with no particles and no file
        if( self%snapshot_nptcls > 0 )then
            call self%meta_snapshot%set(id=self%last_snapshot_id, snapshot_filename=self%snapshot_dir//'/'//self%snapshot_filename,&
                &snapshot_nptcls=self%snapshot_nptcls)
        else
            call self%meta_snapshot%set(id=self%last_snapshot_id, snapshot_filename=string(''), snapshot_nptcls=0)
        endif
        call self%pipe%send_meta(self%meta_snapshot)
        if( .not. allocated(self%snapshot_idx) ) return
        if( size(self%snapshot_idx) == 0 .or. self%snapshot_ntilesx <= 0 ) return
        call send_cavgs(self%pipe, self%meta_snapshot_cavgs, self%snapshot_jpeg, self%snapshot_idx, self%snapshot_mrc,&
            &self%snapshot_ntilesx, self%snapshot_ntilesy, res=self%snapshot_res, pop=self%snapshot_pop)
    end subroutine send_snapshot

    !---------------- rules ----------------

    !> In iterations 2 to LATE_ITER a pool pauses once more than one iteration has passed since the
    !! last import, later once one has; returns the rate factor of the particles it then waits for
    !! (target_nptcls), 0 when no pause is due.
    pure integer function pause_rate_factor( iter, iter_last_import )
        integer, intent(in) :: iter, iter_last_import
        pause_rate_factor = 0
        if( iter > 1 .and. iter <= LATE_ITER )then
            if( iter > iter_last_import + 1 ) pause_rate_factor = EARLY_RATE_FACTOR
        else if( iter > LATE_ITER )then
            if( iter > iter_last_import ) pause_rate_factor = LATE_RATE_FACTOR
        endif
    end function pause_rate_factor

    !> The particles to wait for: @p nptcls_now plus the larger of NPTCLS_PER_CLS_MIN per class and
    !! @p factor micrographs' worth at @p rate particles per micrograph.
    pure integer function target_nptcls( nptcls_now, ncls, rate, factor )
        integer, intent(in) :: nptcls_now, ncls, rate, factor
        target_nptcls = nptcls_now + max(ncls * NPTCLS_PER_CLS_MIN, rate * factor)
    end function target_nptcls

    !> Whether the pool runs without pausing: from the sieve's final set until FINAL_ITER.
    pure logical function runs_to_final( iter, l_sieve_final )
        integer, intent(in) :: iter
        logical, intent(in) :: l_sieve_final
        runs_to_final = l_sieve_final .and. iter < FINAL_ITER
    end function runs_to_final

    !> The mask diameter (A) when none is given: the box less the soft edge and a pixel. The stage
    !! this replaces used the pixel count as Angstroms.
    pure real function default_mskdiam( box, smpd )
        integer, intent(in) :: box
        real,    intent(in) :: smpd
        default_mskdiam = smpd * real(2 * nint(real(box/2) - COSMSKHALFWIDTH - 1.))
    end function default_mskdiam

    !> The publication number @p id in the completed folder of the stage directory @p cwd.
    function publication_fname( cwd, id ) result( fname )
        class(string), intent(in) :: cwd
        integer,       intent(in) :: id
        type(string) :: fname
        fname = cwd//'/'//DIR_STREAM_COMPLETED//int2str_pad(id, 5)//METADATA_EXT
    end function publication_fname

end module simple_stream_stage_pool2D
