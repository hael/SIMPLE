!@descr: state and steps of stream task 7 (multistate 3D): import the particles pool 2D publishes, make reference volumes from class averages once, then assign each publication's (sub-)classes to them, refine the particles' poses locally with the states fixed and report to the GUI
!==============================================================================
! MODULE: simple_stream_stage_solve3D
!
! PURPOSE:
!   The body of stream p07 as a type (doc/policies/stream/stream_3D_ingestion_policy.md;
!   doc/refactoring_notes/planned/stream_p07_state_assignment_plan_2026-10-13.md). Pool 2D
!   publishes its classified state; each pass in which no job runs takes the newest publication.
!   The first (in a fresh pool, the sieve's own 2D of the sets handed off) is taken as published,
!   without the pool model: its rows are the first set, whose selection stays for the session.
!
!   The first set makes the states' reference volumes once, from class averages (the route):
!   solve2D classifies it into 150 classes, the pool model selects them (the particles of a
!   rejected class are deselected for the session), cls_expansion splits the kept ones into
!   about 500 in all, and solve3D_cavgs reconstructs the volumes from those. Its result is
!   checked (refs_usable: the even/odd angular distance and state agreement of the class
!   averages, and at least two states each holding 5% of them); passed, the volumes of the
!   states that hold class averages are copied into state_refs/, renumbered from 1, and are the
!   references for the session. Failed, the route runs again (ROUTE_TRIES in all), then
!   solve3D_cavgs runs once on the next publication's pool class averages (the fallback); when
!   every attempt fails, the attempt with the smallest even/odd angular distance is taken, with
!   a warning.
!
!   Every later publication is merged into the stage's rows; its class averages are scored with
!   the pool model, which is logged but not applied. Once the references exist, each is
!   assigned in assign/it_<n>/: cls_expansion splits its classes into about 500 sub-classes
!   (expand/; when it fails, the pool's class averages are matched instead), assign_cavgs_states
!   (simple_private_exec) matches every sub-class average against the references at 15 A, a
!   (sub-)class whose best correlation is an outlier of the publication is rejected
!   (reject_by_fit) and one with a small margin flagged ambiguous (flag_ambiguous), each kept
!   one's pose and state are mapped onto its particles (map2ptcls), and refine3D refine=cont
!   refines each particle's pose locally from there, against its own state, which stays, until
!   it converges, reconstructing the states from the particles (refine3D/). The GUI gets the
!   references, then each assignment's volumes, populations, resolutions and orientations.
!   The commander (simple_commanders_stream_p07_solve3D_multistate) only normalises the
!   command line and loops over iterate() until finished().
!
!   What happens is delegated:
!     - class-average scoring -> simple_cavg_quality_selection
!     - jobs                  -> simple_qsys_async_job (solve2D, cls_expansion, solve3D_cavgs,
!                                assign_cavgs_states, refine3D)
!     - GUI                   -> simple_stream_pipe, simple_stream_gui_senders,
!                                simple_oris_utils (the orientation histograms)
!
!   The stage's rows only ever grow (rows are read by index, never renumbered): a publication's
!   stacks are matched to the rows by stack name, so the order of a restarted pool does not
!   matter.
!
! LIFECYCLE:
!   new(cline) -> { iterate() } until finished() -> finalize() -> kill()
!
!   A publication the stage cannot use (its stacks disagree with the rows, its class averages
!   or FRCs are missing) is passed over with a warning and listed in REJECTED_PUBLICATIONS, and
!   the stage waits for the next one.
!
! RESTART:
!   A stop cancels the running job (finalize). A restart removes a leftover TERM_STREAM, takes
!   the references in state_refs/ when there are (the route is skipped), and starts again from
!   the newest publication not rejected, which is assigned (or makes the references). A job
!   folder left unfinished is first moved aside (fresh_job_dir).
!==============================================================================
module simple_stream_stage_solve3D
use simple_defs,                                      only: logfhandle
use simple_defs_fname,                                only: TERM_STREAM, METADATA_EXT, MRC_EXT, JPG_EXT, PPROC_SUFFIX,&
                                                           &LP_SUFFIX, MIRR_SUFFIX, DIR_SNAPSHOT
use simple_defs_stream,                               only: DIR_STREAM_COMPLETED, SHORTWAIT, WAITTIME, OPTICS_ID_DELTA
use simple_defs_environment,                          only: SIMPLE_STREAM_SOLVE3D_PARTITION
use simple_error,                                     only: simple_exception
use simple_string,                                    only: string
use simple_string_utils,                              only: int2str, int2str_pad, lex_sort
use simple_fileio,                                    only: add2fbody, basename, del_file, file2rarr, file_exists, get_fbody, fname2ext,&
                                                           &get_fpath, simple_abspath, simple_getcwd, swap_suffix, simple_copy_file
use simple_syslib,                                    only: dir_exists, simple_mkdir, simple_list_dirs, simple_rmdir
use simple_math_ft,                                   only: get_resarr
use simple_estimate_ssnr,                             only: get_resolution
use simple_imghead,                                   only: get_mrc_minmax, find_ldim_nptcls
use simple_refine3D_fnames,                           only: refine3D_oris_heatmap_fname, refine3D_reprojs_fname
use simple_oris_utils,                                only: oridist_from_oris
use simple_cmdline,                                   only: cmdline
use simple_parameters,                                only: parameters
use simple_sp_project,                                only: sp_project
use simple_qsys_env,                                  only: qsys_env
use simple_qsys_funs,                                 only: qsys_cleanup
use simple_qsys_async_job,                            only: qsys_async_job, ASYNC_JOB_RUNNING, ASYNC_JOB_DONE, ASYNC_JOB_FAILED
use simple_qsys_job_record,                           only: fresh_job_dir
use simple_rec_list,                                  only: rec_list, rec_iterator, chunk_rec
use simple_stream_watcher,                            only: stream_watcher
use simple_stream_state,                              only: ipc_pipe_solve3D_multistate_in, ipc_pipe_solve3D_multistate_out
use simple_stream_utils,                              only: create_stream_project, init_stream_qenv
use simple_cavg_quality_model,                        only: CAVG_QUALITY_MODEL_POOL_DEFAULT
use simple_state_assignment,                          only: MIN_PTCLS_PER_STATE, SOLVE2D_DIR, CLSEXP_DIR, CAVGS3D_DIR,&
                                                           &ASSIGN_LABEL, REFINE_DIR, REFS_DIR, REFS_RECORD, QUALITY_DIR,&
                                                           &refs_check, assign_stats, expansion_budget, count_active_classes,&
                                                           &fit_mskdiam, usable_job_2D, score_cavgs_with_model,&
                                                           &select_route_classes, model_rejected_rows, check_state_refs,&
                                                           &read_refs_record, map_assignment, write_matched_cavgs,&
                                                           &set_refine_seeds, write_refine_project, map_refine_result,&
                                                           &cavgs_fbody, log_jpeg, set_solve2D_route_cline, set_clsexp_cline,&
                                                           &set_cavgs3D_cline, set_assign_cline, set_refine_cont_cline
use simple_gui_metadata_utils,                        only: max_metadata_size
use simple_gui_metadata_types,                        only: GUI_METADATA_STREAM_SOLVE3D_MULTISTATE_TYPE, GUI_METADATA_VOL3D_TYPE,&
                                                           &GUI_METADATA_STREAM_SOLVE3D_MULTISTATE_REPROJ_TYPE,&
                                                           &GUI_METADATA_STREAM_SOLVE3D_SNAPSHOT_TYPE
use simple_gui_metadata_stream_snapshot,              only: gui_metadata_stream_snapshot
use simple_gui_metadata_stream_update,                only: gui_metadata_stream_update
use simple_gui_metadata_cavg2D,                       only: gui_metadata_cavg2D
use simple_gui_metadata_vol3D,                        only: gui_metadata_vol3D, MAX_FSC_VOL3D, ORIDIST_NBINS_X, ORIDIST_NBINS_Y
use simple_gui_metadata_stream_solve3D_multistate, only: gui_metadata_stream_solve3D_multistate
use simple_stream_pipe,                               only: stream_pipe
use simple_stream_gui_senders,                        only: send_reproj_tiles
implicit none

public :: stream_stage_solve3D
public :: PHASE_IMPORTING, PHASE_SOLVE2D, PHASE_CLSEXP, PHASE_CAVGS3D, PHASE_IDLE, PHASE_EXPAND, PHASE_ASSIGN, PHASE_REFINE
public :: JOB_NONE, JOB_ROUTE, JOB_FALLBACK, JOB_ASSIGN, ROUTE_TRIES
private
#include "simple_local_flags.inc"

! phases of the stage
integer, parameter :: PHASE_IMPORTING = 0 ! no reference volumes and no job
integer, parameter :: PHASE_SOLVE2D   = 1 ! the route's solve2D on the first set
integer, parameter :: PHASE_CLSEXP    = 2 ! the route's cls_expansion of the first set's selected classes
integer, parameter :: PHASE_CAVGS3D   = 3 ! solve3D_cavgs: the route's, or the fallback's on a publication's classes
integer, parameter :: PHASE_IDLE      = 4 ! the reference volumes exist; waiting for a publication
integer, parameter :: PHASE_EXPAND    = 5 ! cls_expansion of a publication's classes
integer, parameter :: PHASE_ASSIGN    = 6 ! assign_cavgs_states: the (sub-)class averages matched against the references
integer, parameter :: PHASE_REFINE    = 7 ! refine3D refine=cont: the particles' poses refined locally, the states fixed
! the job to start next (next_job)
integer, parameter :: JOB_NONE        = 0
integer, parameter :: JOB_ROUTE       = 1 ! the route to the reference volumes, on the first set
integer, parameter :: JOB_FALLBACK    = 2 ! solve3D_cavgs on a later publication's pool class averages
integer, parameter :: JOB_ASSIGN      = 3 ! a publication's assignment (expansion, matching, local refinement)

! the route's attempts, before the fallback (the method's constants: simple_state_assignment)
integer, parameter :: ROUTE_TRIES         = 2
logical, parameter :: L_MODEL_LOGGED      = .true. ! the pool model scores each publication, logged beside the fit

character(len=*), parameter :: ASSIGN_DIR       = 'assign'
character(len=*), parameter :: EXPAND_DIR       = 'expand'
character(len=*), parameter :: RESULT_PROJ      = 'assignment' ! the latest assignment's rows, refined: 3D snapshots' source
! the quality folders kept: the newest NQUALITY_KEPT, and those of the publications the route or
! the fallback started from, listed in QUALITY_RUNS
integer,          parameter :: NQUALITY_KEPT    = 3
character(len=*), parameter :: QUALITY_RUNS     = 'quality_selection/runs.txt'
! the publications passed over as unusable, one file name a line (a restart passes them over too)
character(len=*), parameter :: REJECTED_PUBLICATIONS = 'rejected_publications.txt'

! Components and steps are public so simple_stream_stage_solve3D_tester can assemble a stage
! and run one step at a time; production code uses new/iterate/finished/finalize/kill.
type :: stream_stage_solve3D
    type(parameters), allocatable                   :: params
    ! allocatable (compile-time policy), allocated in init_params and released in kill
    type(sp_project), allocatable                   :: spproj        ! the rows: imported particles, their states and poses
    type(qsys_env),   allocatable                   :: qenv
    type(qsys_async_job)                            :: job           ! the running job
    type(stream_watcher)                            :: project_buff  ! the publications of pool 2D
    type(rec_list)                                  :: setslist      ! one record per publication; included once taken or passed over
    type(stream_pipe)                               :: pipe          ! to the master
    type(gui_metadata_stream_solve3D_multistate)    :: meta_status
    type(gui_metadata_cavg2D)                       :: meta_reproj
    type(gui_metadata_stream_snapshot)              :: meta_snapshot
    type(string),   allocatable :: stk_names(:)    ! the stacks in the pool
    real,           allocatable :: state_res(:)    ! FSC=0.143 resolution per state of the latest refinement (0: none)
    integer,        allocatable :: state_pop(:)    ! population per state of the latest refinement (0: none)
    type(string),   allocatable :: ref_vols(:)     ! the reference volumes, one per state (state_refs/); none before
    real,           allocatable :: try_angdist(:)  ! each attempt's even/odd angular distance (-1: no volumes)
    logical,        allocatable :: first_set(:)    ! the rows the first publication left selected: their selection stays
    logical,        allocatable :: model_rej(:)    ! the selected rows the pool model rejects in the newest publication (logged)
    logical,        allocatable :: fit_desel(:)    ! the rows the latest assignment deselected: a first-set row among them is selected again by the next merge
    type(string)                :: result_projfile ! the latest local refinement's project: 3D snapshots' source
    type(string)                :: assign_dir      ! the running assignment's folder (assign/it_<n>)
    type(string)                :: last_stem       ! the latest publication taken (its quality folder's name)
    integer :: phase              = PHASE_IMPORTING
    integer :: nassign            = 0  ! assignments started
    integer :: nroute_tries       = 0  ! the route's attempts started
    integer :: nptcls_at_last_run = 0  ! particles the latest assignment gave a state
    integer :: nptcls_selected    = 0  ! selected particles in the rows
    integer :: ref_box            = 0  ! the reference volumes' box and pixel size (A)
    real    :: ref_smpd           = 0.
    integer :: last_snapshot_id   = 0  ! the latest 3D snapshot request answered; each is written once
    integer :: optics_id_offset   = 0  ! optics group ids of the snapshots' STAR files, per GUI display
    logical :: l_route_due        = .false. ! the route is to run on the first set
    logical :: l_fallback_pending = .false. ! the route failed: the fallback waits for a later publication
    logical :: l_fallback_due     = .false. ! that publication is in: the fallback is to run
    logical :: l_fallback_done    = .false. ! the fallback has run (once a session)
    logical :: l_assign_due       = .false. ! a publication waits for its assignment
    logical :: l_expanded         = .false. ! the running assignment matches cls_expansion's sub-classes
    logical :: l_final_seen       = .false. ! the newest publication taken is the pool's final one
    logical :: l_dock_refs        = .true.  ! the references are compared after docking (a test without the docking's
                                            ! time sets it .false.: no comparison)
    real    :: mskdiam            = 0. ! pool 2D's mask diameter (A), from the latest publication taken
    logical :: l_restart          = .false.
    logical :: l_attached         = .false. ! pool 2D's completed folder exists and is watched
    logical :: l_waiting_logged   = .false.
    logical :: l_exists           = .false.
    ! waits (s); tests set them to 0, and settle_s to -1 to take files written in the same second
    integer :: settle_s           = SHORTWAIT ! an export is taken once untouched longer than this;
                                              ! pool 2D writes exports with a rename
    integer :: wait_s             = WAITTIME  ! pause at the end of a pass
contains
    procedure :: new
    procedure :: iterate
    procedure :: finished
    procedure :: finalize
    procedure :: kill
    ! the steps of new(), in order
    procedure :: init_params
    procedure :: init_queue
    procedure :: init_gui
    ! the steps of iterate() and their helpers
    procedure :: attach_upstream
    procedure :: watch_sets
    procedure :: take_mskdiam
    procedure :: import_sets
    procedure :: publication_problem
    procedure :: rows_problem
    procedure :: flag_model_rows
    procedure :: take_first_publication
    procedure :: in_first_set
    procedure :: fit_deselected
    procedure :: merge_publication
    procedure :: take_cavgs
    procedure :: stack_index
    procedure :: advance_jobs
    procedure :: job_failed
    procedure :: job_done
    procedure :: set_job_resources
    ! the route to the reference volumes
    procedure :: start_solve2D_first
    procedure :: finish_solve2D_first
    procedure :: select_first_set_classes
    procedure :: start_clsexp
    procedure :: finish_clsexp
    procedure :: adopt_job_2D
    procedure :: start_fallback
    procedure :: start_cavgs3D
    procedure :: finish_cavgs3D
    procedure :: register_attempt
    procedure :: next_attempt
    procedure :: accept_best
    procedure :: accept_refs
    procedure :: restore_refs
    procedure :: n_refs
    ! the assignment of a publication
    procedure :: start_expand
    procedure :: finish_expand
    procedure :: start_assign
    procedure :: finish_assign
    procedure :: take_assignment
    procedure :: write_assignment_cavgs
    procedure :: prepare_refine_seeds
    procedure :: start_refine
    procedure :: finish_refine
    procedure :: take_refine_result
    ! records
    procedure :: record_run_publication
    procedure :: prune_quality_dirs
    procedure :: prune_run_dirs
    procedure :: write_stage_project
    ! GUI
    procedure :: apply_gui_updates
    procedure :: write_snapshot
    procedure :: send_status
    procedure :: send_volumes
    procedure :: send_refs
    procedure :: send_state_volume
    procedure :: read_state_fsc
    ! rules that use no stage state
    procedure, nopass :: next_job
    procedure, nopass :: best_attempt
    procedure, nopass :: is_final_publication
end type stream_stage_solve3D

contains

    !---------------- lifecycle ----------------

    !> Builds the parameters from the normalised @p cline, the queue environment of the jobs and
    !! the GUI.
    subroutine new( self, cline )
        class(stream_stage_solve3D), intent(inout) :: self
        class(cmdline),              intent(inout) :: cline
        call self%kill()
        call self%init_params(cline)
        self%l_exists = .true. ! from here kill() releases what has been built
        call self%init_queue()
        call self%init_gui(ipc_pipe_solve3D_multistate_out(1), ipc_pipe_solve3D_multistate_in(2))
        call self%send_status()
        if( self%n_refs() > 0 ) call self%send_refs()
    end subroutine new

    !> The stage's project file (made and given a computing environment) and its parameters; the
    !! project must start without micrographs. A restart's leftover TERM_STREAM is removed in the
    !! stage's folder, where the loop looks for it. Reference volumes a previous run left in
    !! state_refs/ are taken (restore_refs).
    subroutine init_params( self, cline )
        class(stream_stage_solve3D), intent(inout) :: self
        class(cmdline),              intent(inout) :: cline
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
        if( .not. allocated(self%spproj) ) allocate(self%spproj)
        if( .not. allocated(self%qenv)   ) allocate(self%qenv)
        call create_stream_project(self%spproj, cline, string('3Dmultistate'))
        if( .not. allocated(self%params) ) allocate(self%params)
        call self%params%new(cline)
        if( self%l_restart )then
            write(logfhandle,'(A)') '>>> RESTARTING EXISTING JOB'
            call del_file(TERM_STREAM)
        endif
        call self%spproj%read(self%params%projfile)
        if( self%spproj%os_mic%get_noris() /= 0 )then
            THROW_HARD('commander_stream_p07_solve3D_multistate must start from an empty project (e.g. from root project folder)')
        endif
        allocate(self%stk_names(0), self%try_angdist(0))
        allocate(self%state_res(self%params%nstates), source=0.)
        allocate(self%state_pop(self%params%nstates), source=0)
        self%optics_id_offset = max(self%params%nicedispid - 1, 0) * OPTICS_ID_DELTA
        ! the references of a previous run: the route is skipped
        if( self%restore_refs(string(REFS_DIR)) )then
            self%params%nstates = self%n_refs()
            self%phase          = PHASE_IDLE
            write(logfhandle,'(A,I3,A)') '>>> ', self%n_refs(), ' REFERENCE VOLUMES TAKEN FROM '//REFS_DIR//'; THE ROUTE IS SKIPPED'
        endif
    end subroutine init_params

    !> The queue environment the jobs are submitted through; on the preprocessing partition, as
    !! before.
    subroutine init_queue( self )
        class(stream_stage_solve3D), intent(inout) :: self
        call init_stream_qenv(self%params, self%qenv, string(SIMPLE_STREAM_SOLVE3D_PARTITION))
    end subroutine init_queue

    !> The GUI metadata objects and the pipe ends from and to the master (-1: none).
    subroutine init_gui( self, fd_read, fd_write )
        class(stream_stage_solve3D), intent(inout) :: self
        integer,                     intent(in)    :: fd_read, fd_write
        call self%meta_status%new(GUI_METADATA_STREAM_SOLVE3D_MULTISTATE_TYPE)
        call self%meta_reproj%new(GUI_METADATA_STREAM_SOLVE3D_MULTISTATE_REPROJ_TYPE)
        call self%meta_snapshot%new(GUI_METADATA_STREAM_SOLVE3D_SNAPSHOT_TYPE)
        call self%pipe%new(fd_read, fd_write, max_metadata_size(), 'solve3D_multistate')
    end subroutine init_gui

    !> One pass: wait for pool 2D's folder; take new exports (unless a job runs); advance the
    !! jobs; answer the GUI; report.
    subroutine iterate( self )
        class(stream_stage_solve3D), intent(inout) :: self
        if( .not. self%l_attached )then
            call self%attach_upstream()
            if( .not. self%l_attached )then
                call self%send_status()
                call self%apply_gui_updates()
                call sleep(self%wait_s)
                return
            endif
        endif
        call self%watch_sets()
        ! no import while a job runs on the rows
        if( self%phase == PHASE_IMPORTING .or. self%phase == PHASE_IDLE ) call self%import_sets()
        call self%advance_jobs()
        call self%apply_gui_updates()
        call self%send_status()
        call sleep(self%wait_s)
    end subroutine iterate

    !> .true. once the stream is told to stop.
    logical function finished( self )
        class(stream_stage_solve3D), intent(in) :: self
        finished = file_exists(TERM_STREAM)
    end function finished

    !> The last status, and the stage's project with the latest assignment and the particles
    !! imported since. A running job is cancelled: a restart starts it again in a fresh directory.
    subroutine finalize( self )
        class(stream_stage_solve3D), intent(inout) :: self
        call self%job%cancel()
        call self%send_status(string('terminating'))
        if( self%spproj%os_ptcl3D%get_noris() > 0 ) call self%write_stage_project()
        call qsys_cleanup(self%params)
    end subroutine finalize

    subroutine kill( self )
        class(stream_stage_solve3D), intent(inout) :: self
        if( allocated(self%spproj) )then
            call self%spproj%kill
            deallocate(self%spproj)
        endif
        if( allocated(self%qenv) )then
            call self%qenv%kill
            deallocate(self%qenv)
        endif
        if( .not. self%l_exists )then
            if( allocated(self%params) ) deallocate(self%params)
            return
        endif
        call self%job%kill
        call self%project_buff%kill
        call self%setslist%kill
        call self%pipe%kill
        call self%meta_status%kill
        call self%meta_reproj%kill
        call self%meta_snapshot%kill
        call self%result_projfile%kill
        call self%assign_dir%kill
        call self%last_stem%kill
        if( allocated(self%ref_vols)    ) deallocate(self%ref_vols)
        if( allocated(self%try_angdist) ) deallocate(self%try_angdist)
        if( allocated(self%first_set)   ) deallocate(self%first_set)
        if( allocated(self%model_rej)   ) deallocate(self%model_rej)
        if( allocated(self%fit_desel)   ) deallocate(self%fit_desel)
        if( allocated(self%stk_names)   ) deallocate(self%stk_names)
        if( allocated(self%state_res)   ) deallocate(self%state_res)
        if( allocated(self%state_pop)   ) deallocate(self%state_pop)
        if( allocated(self%params)      ) deallocate(self%params)
        self%phase              = PHASE_IMPORTING
        self%nassign            = 0
        self%nroute_tries       = 0
        self%nptcls_at_last_run = 0
        self%nptcls_selected    = 0
        self%ref_box            = 0
        self%ref_smpd           = 0.
        self%last_snapshot_id   = 0
        self%optics_id_offset   = 0
        self%l_route_due        = .false.
        self%l_fallback_pending = .false.
        self%l_fallback_due     = .false.
        self%l_fallback_done    = .false.
        self%l_assign_due       = .false.
        self%l_expanded         = .false.
        self%l_final_seen       = .false.
        self%l_dock_refs        = .true.
        self%mskdiam            = 0.
        self%l_restart          = .false.
        self%l_attached         = .false.
        self%l_waiting_logged   = .false.
        self%l_exists           = .false.
    end subroutine kill

    !---------------- steps ----------------

    ! Starts watching pool 2D's completed folder once it exists.
    subroutine attach_upstream( self )
        class(stream_stage_solve3D), intent(inout) :: self
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

    ! One record per new export, in export order (the zero-padded names sort that way).
    subroutine watch_sets( self )
        class(stream_stage_solve3D), intent(inout) :: self
        type(string), allocatable :: projects(:)
        integer :: nprojects, i
        call self%project_buff%watch(nprojects, projects)
        if( nprojects == 0 ) return
        call lex_sort(projects)
        call self%project_buff%add2history(projects)
        do i = 1,nprojects
            call self%setslist%push2chunk_list(projects(i), self%setslist%size() + 1, .true.)
        enddo
    end subroutine watch_sets

    ! The mask diameter pool 2D records with the class averages of publication @p set, taken from
    ! every publication: a change applies from the next job (a running one keeps its own) and is
    ! logged. Without one, the latest stays (the box default before any).
    subroutine take_mskdiam( self, set )
        class(stream_stage_solve3D), intent(inout) :: self
        class(sp_project),           intent(inout) :: set
        type(string) :: stk
        real         :: smpd, mskdiam
        integer      :: ncls
        call set%get_cavgs_stk(stk, ncls, smpd, fail=.false.)
        if( ncls <= 0 ) return
        mskdiam = 0.
        call set%get_mskdiam('cavg', mskdiam)
        if( mskdiam <= 0. ) return
        if( abs(mskdiam - self%mskdiam) < 0.5 ) return
        if( self%mskdiam > 0. )then
            write(logfhandle,'(A,F8.2,A,F8.2,A)') '>>> MASK DIAMETER CHANGED FROM ', self%mskdiam, ' TO ', mskdiam,&
                &' A; IT APPLIES FROM THE NEXT JOB'
        else
            write(logfhandle,'(A,F8.2)') '>>> MASK DIAMETER SET TO : ', mskdiam
        endif
        self%mskdiam = mskdiam
    end subroutine take_mskdiam

    ! The newest publication not yet taken. Into a stage without rows it is the first set, taken
    ! as the pool selects it (take_first_publication): without references the route is due on
    ! it, with references (a restart) its assignment. Every later one is scored with the pool
    ! model (logged, not applied) and merged into the rows; once the references exist its
    ! assignment is due, and while the fallback waits for a publication, the fallback is.
    ! Its class averages and FRCs come with its classes (take_cavgs). Older publications not
    ! taken are passed over: the newest holds what they held. A publication the stage cannot use
    ! (publication_problem) is passed over with a warning and listed in REJECTED_PUBLICATIONS, so
    ! a restart passes it over too; the next one is waited for.
    subroutine import_sets( self )
        class(stream_stage_solve3D), intent(inout) :: self
        type(sp_project)     :: set
        type(rec_iterator)   :: it
        type(chunk_rec)      :: crec, newest
        type(string)         :: stem, problem
        integer, allocatable :: model_states(:)
        integer              :: irec, inewest
        logical              :: l_scored
        if( self%setslist%size() == 0 ) return
        ! the records are in publication order: the last one not taken is the newest
        inewest = 0
        it = self%setslist%begin()
        do irec = 1,self%setslist%size()
            call it%get(crec)
            if( .not. crec%included )then
                inewest = irec
                newest  = crec
            endif
            call it%next()
        enddo
        if( inewest == 0 ) return
        stem = get_fbody(basename(newest%projfile), METADATA_EXT, separator=.false.)
        if( is_rejected(basename(newest%projfile)) )then
            write(logfhandle,'(A,A,A)') '>>> PUBLICATION ', stem%to_char(), ' WAS REJECTED BEFORE; PASSED OVER'
        else
            call set%read(newest%projfile)
            problem = self%publication_problem(set)
            if( problem%strlen() > 0 )then
                THROW_WARN('publication '//stem%to_char()//' passed over: '//problem%to_char())
                call add_rejected(basename(newest%projfile))
            else
                call self%take_mskdiam(set)
                if( self%spproj%os_ptcl3D%get_noris() == 0 )then
                    call self%take_first_publication(set, newest%id)
                    if( self%n_refs() > 0 )then
                        self%l_assign_due = .true.
                    else
                        self%l_route_due  = .true.
                        write(logfhandle,'(A)') '>>> THE FIRST SET''S REFERENCE VOLUMES: SOLVE2D, THE POOL MODEL''S SELECTION, '//&
                            &'CLS_EXPANSION AND SOLVE3D_CAVGS'
                    endif
                else
                    ! the pool model's verdict, logged beside the fit's, not applied
                    l_scored = .false.
                    if( L_MODEL_LOGGED ) l_scored = score_cavgs_with_model(set, CAVG_QUALITY_MODEL_POOL_DEFAULT,&
                        &self%mskdiam, string(QUALITY_DIR//'/')//stem, model_states)
                    call self%merge_publication(set, newest%id)
                    if( l_scored )then
                        call self%flag_model_rows(model_states)
                    else if( allocated(self%model_rej) )then
                        deallocate(self%model_rej)
                    endif
                    if( self%n_refs() > 0 )then
                        self%l_assign_due = .true.
                    else if( self%l_fallback_pending )then
                        self%l_fallback_due = .true.
                    endif
                endif
                call self%take_cavgs(set, stem)
                ! the pool's final publication; a later one that is not final (the sieve took
                ! finality back) withdraws it
                self%l_final_seen = is_final_publication(set)
                if( self%l_final_seen ) write(logfhandle,'(A,A)') '>>> THE POOL''S FINAL PUBLICATION: ', stem%to_char()
                self%last_stem = stem
                call self%prune_quality_dirs()
            endif
            call set%kill
        endif
        ! the newest and every older publication are taken
        it = self%setslist%begin()
        do irec = 1,inewest
            call it%get(crec)
            if( .not. crec%included )then
                crec%included = .true.
                call self%setslist%replace_iterator(it, crec)
            endif
            call it%next()
        enddo

    contains

        logical function is_rejected( fname )
            class(string), intent(in) :: fname
            character(len=4096) :: line
            integer :: funit, ios
            is_rejected = .false.
            if( .not. file_exists(string(REJECTED_PUBLICATIONS)) ) return
            open(newunit=funit, file=REJECTED_PUBLICATIONS, status='old', action='read', iostat=ios)
            if( ios /= 0 ) return
            do
                read(funit,'(A)',iostat=ios) line
                if( ios /= 0 ) exit
                if( trim(line) == fname%to_char() )then
                    is_rejected = .true.
                    exit
                endif
            enddo
            close(funit)
        end function is_rejected

        subroutine add_rejected( fname )
            class(string), intent(in) :: fname
            integer :: funit, ios
            open(newunit=funit, file=REJECTED_PUBLICATIONS, status='unknown', position='append', action='write', iostat=ios)
            if( ios /= 0 ) return
            write(funit,'(A)') fname%to_char()
            close(funit)
        end subroutine add_rejected

    end subroutine import_sets

    ! Why publication @p set cannot be used, '' when it can: what merge_publication and take_cavgs
    ! require of it, checked before anything changes (class averages, FRCs, rows_problem).
    function publication_problem( self, set ) result( problem )
        class(stream_stage_solve3D), intent(inout) :: self
        class(sp_project),           intent(inout) :: set
        type(string) :: problem, stk, frcs
        real    :: smpd
        integer :: ncls, ldim(3), nimgs
        problem = ''
        if( set%os_cls2D%get_noris() == 0 )then
            problem = 'no cls2D entries'
            return
        endif
        call set%get_cavgs_stk(stk, ncls, smpd, fail=.false.)
        if( ncls <= 0 )then
            problem = 'no class averages'
            return
        endif
        if( .not. file_exists(stk) )then
            problem = 'its class averages are missing: '//stk%to_char()
            return
        endif
        call find_ldim_nptcls(stk, ldim, nimgs)
        if( nimgs /= set%os_cls2D%get_noris() )then
            problem = '# class averages /= # cls2D entries'
            return
        endif
        call set%get_frcs(frcs, 'frc2D', fail=.false.)
        if( .not. file_exists(frcs) )then
            problem = 'its FRCs are missing'
            return
        endif
        problem = self%rows_problem(set)
    end function publication_problem

    ! Why the stacks of publication @p set disagree with the stage's rows, '' when they agree: a
    ! micrograph per stack, and for each stack the stage holds the same size and image indices
    ! within it that match the rows' (merge_publication).
    function rows_problem( self, set ) result( problem )
        class(stream_stage_solve3D), intent(inout) :: self
        class(sp_project),           intent(inout) :: set
        type(string) :: problem, name
        integer :: nstks, jstk, k, hint, fromp, fromp_set, nptcls, i, jptcl, iimg, iptcl
        problem = ''
        nstks   = set%os_stk%get_noris()
        if( set%os_mic%get_noris() /= nstks )then
            problem = '# micrographs /= # stacks'
            return
        endif
        hint = 0
        do jstk = 1,nstks
            name = set%os_stk%get_str(jstk, 'stk')
            k    = self%stack_index(name, hint + 1)
            if( k == 0 ) cycle
            hint      = k
            fromp_set = set%os_stk%get_fromp(jstk)
            nptcls    = set%os_stk%get_top(jstk) - fromp_set + 1
            fromp     = self%spproj%os_stk%get_fromp(k)
            if( self%spproj%os_stk%get_top(k) - fromp + 1 /= nptcls )then
                problem = 'a stack changed size: '//name%to_char()
                return
            endif
            do i = 0,nptcls - 1
                jptcl = fromp_set + i
                iimg  = i + 1
                if( set%os_ptcl2D%isthere(jptcl, 'indstk') ) iimg = set%os_ptcl2D%get_int(jptcl, 'indstk')
                if( iimg < 1 .or. iimg > nptcls )then
                    problem = 'an image index outside its stack: '//name%to_char()
                    return
                endif
                iptcl = fromp + iimg - 1
                if( self%spproj%os_ptcl2D%isthere(iptcl, 'indstk') )then
                    if( self%spproj%os_ptcl2D%get_int(iptcl, 'indstk') /= iimg )then
                        problem = 'rows and images disagree: '//name%to_char()
                        return
                    endif
                endif
            enddo
        enddo
    end function rows_problem

    ! The selected rows whose class (the newest publication's, merged) the pool model rejects,
    ! @p cls_states its verdict per class: logged beside the fit's rejection (take_assignment).
    subroutine flag_model_rows( self, cls_states )
        class(stream_stage_solve3D), intent(inout) :: self
        integer,                     intent(in)    :: cls_states(:)
        if( allocated(self%model_rej) ) deallocate(self%model_rej)
        self%model_rej = model_rejected_rows(self%spproj, cls_states)
    end subroutine flag_model_rows

    ! The first publication @p set (number @p id), into a stage without rows: merged as published
    ! (the sieve's selection, or the pool's when p06 published its own), without the pool model's
    ! decision. Its selected rows are the first set: their selection stays for the session, and
    ! later publications change only their 2D parameters (merge_publication).
    subroutine take_first_publication( self, set, id )
        class(stream_stage_solve3D), intent(inout) :: self
        class(sp_project),           intent(inout) :: set
        integer,                     intent(in)    :: id
        integer :: i, n
        if( self%spproj%os_ptcl3D%get_noris() > 0 ) THROW_HARD('the first publication needs a stage without rows')
        if( allocated(self%first_set) ) deallocate(self%first_set)
        call self%merge_publication(set, id)
        n = self%spproj%os_ptcl3D%get_noris()
        allocate(self%first_set(n))
        do i = 1,n
            self%first_set(i) = self%spproj%os_ptcl3D%get_state(i) > 0
        enddo
        self%nptcls_selected = count(self%first_set)
        write(logfhandle,'(A,I6,A,I8,A)') '>>> PUBLICATION ', id, ' IS THE FIRST SET: ', self%nptcls_selected,&
            &' PARTICLES AS PUBLISHED'
    end subroutine take_first_publication

    ! .true. for row @p iptcl of the first set, whose selection stays for the session.
    logical function in_first_set( self, iptcl )
        class(stream_stage_solve3D), intent(in) :: self
        integer,                     intent(in) :: iptcl
        in_first_set = .false.
        if( .not. allocated(self%first_set) ) return
        if( iptcl < 1 .or. iptcl > size(self%first_set) ) return
        in_first_set = self%first_set(iptcl)
    end function in_first_set

    ! .true. for row @p iptcl when the latest assignment deselected it (take_assignment): its
    ! sub-class was rejected by fit or it had no class of the matching.
    logical function fit_deselected( self, iptcl )
        class(stream_stage_solve3D), intent(in) :: self
        integer,                     intent(in) :: iptcl
        fit_deselected = .false.
        if( .not. allocated(self%fit_desel) ) return
        if( iptcl < 1 .or. iptcl > size(self%fit_desel) ) return
        fit_deselected = self%fit_desel(iptcl)
    end function fit_deselected

    ! Merges the publication @p set (number @p id) into the stage's rows, which only ever grow. A
    ! stack the stage holds is matched by name and its particles by image index in the stack:
    ! they take the publication's 2D parameters and selection in place, and keep their 3D
    ! parameters, state, CTF and optics group (the next assignment gives the state). A row of the
    ! first set takes the 2D parameters only: its selection stays, but a deselection by the latest
    ! assignment's fit lasts only until this publication (the row is selected again, as state 1,
    ! and judged again by the next assignment). A new stack is appended with
    ! its micrograph and particles (2D and 3D). A stack the publication lacks keeps its rows,
    ! deselected (the first set's keep their selection). The classes are the publication's;
    ! take_cavgs then registers its class averages and FRCs with them.
    subroutine merge_publication( self, set, id )
        class(stream_stage_solve3D), intent(inout) :: self
        class(sp_project),           intent(inout) :: set
        integer,                     intent(in)    :: id
        integer,      allocatable :: kstk(:)   ! the stage's stack of each publication stack, 0 when new
        logical,      allocatable :: l_seen(:) ! the stage's stacks the publication holds
        type(string), allocatable :: new_names(:)
        type(string) :: name
        integer :: nstks, nstks_pool, jstk, k, hint, nnew, nptcls_new, nptcls, i, iptcl, jptcl
        integer :: fromp, fromp_set, pool_nmics, pool_nptcls, imic, nupdated, ndeselected, s, iimg
        nstks = set%os_stk%get_noris()
        if( set%os_mic%get_noris() /= nstks ) THROW_HARD('# micrographs /= # stacks in publication '//int2str(id))
        nstks_pool = size(self%stk_names)
        allocate(kstk(nstks), source=0)
        allocate(l_seen(nstks_pool), source=.false.)
        ! match the stacks; a restarted pool may list them in another order
        hint       = 0
        nnew       = 0
        nptcls_new = 0
        do jstk = 1,nstks
            name = set%os_stk%get_str(jstk, 'stk')
            k    = self%stack_index(name, hint + 1)
            kstk(jstk) = k
            if( k > 0 )then
                l_seen(k) = .true.
                hint      = k
            else
                nnew       = nnew + 1
                nptcls_new = nptcls_new + set%os_stk%get_top(jstk) - set%os_stk%get_fromp(jstk) + 1
            endif
        enddo
        ! the stacks the stage holds: 2D parameters and selection in place
        nupdated = 0
        do jstk = 1,nstks
            k = kstk(jstk)
            if( k == 0 ) cycle
            fromp_set = set%os_stk%get_fromp(jstk)
            nptcls    = set%os_stk%get_top(jstk) - fromp_set + 1
            fromp     = self%spproj%os_stk%get_fromp(k)
            if( self%spproj%os_stk%get_top(k) - fromp + 1 /= nptcls ) THROW_HARD('a stack changed size: '//self%stk_names(k)%to_char())
            do i = 0,nptcls - 1
                jptcl = fromp_set + i
                ! the particle's image in its stack, as the publication records it (indstk); the
                ! stage's rows of a stack hold its images in order, and must record the same one
                iimg = i + 1
                if( set%os_ptcl2D%isthere(jptcl, 'indstk') ) iimg = set%os_ptcl2D%get_int(jptcl, 'indstk')
                if( iimg < 1 .or. iimg > nptcls ) THROW_HARD('an image index outside its stack in publication '//int2str(id))
                iptcl = fromp + iimg - 1
                if( self%spproj%os_ptcl2D%isthere(iptcl, 'indstk') )then
                    if( self%spproj%os_ptcl2D%get_int(iptcl, 'indstk') /= iimg ) THROW_HARD('rows and images disagree: '//self%stk_names(k)%to_char())
                endif
                if( self%in_first_set(iptcl) )then
                    ! the selection is the first set's own; only the 2D parameters change
                    s = self%spproj%os_ptcl2D%get_state(iptcl)
                    call self%spproj%os_ptcl2D%transfer_2Dparams(iptcl, set%os_ptcl2D, jptcl)
                    ! the fit's deselection lasts until the next publication: judged again
                    if( self%fit_deselected(iptcl) )then
                        s = 1
                        call self%spproj%os_ptcl3D%set_state(iptcl, 1)
                        self%fit_desel(iptcl) = .false.
                    endif
                    call self%spproj%os_ptcl2D%set_state(iptcl, s)
                    cycle
                endif
                s     = set%os_ptcl2D%get_state(jptcl)
                call self%spproj%os_ptcl2D%transfer_2Dparams(iptcl, set%os_ptcl2D, jptcl)
                call self%spproj%os_ptcl2D%set_state(iptcl, s)
                ! the 3D state is an assignment's: a deselected particle loses it, a selected one
                ! keeps it, and one without (state 0) is selected as state 1 until it is assigned
                if( s == 0 )then
                    call self%spproj%os_ptcl3D%set_state(iptcl, 0)
                else if( self%spproj%os_ptcl3D%get_state(iptcl) == 0 )then
                    call self%spproj%os_ptcl3D%set_state(iptcl, 1)
                endif
            enddo
            nupdated = nupdated + nptcls
        enddo
        ! the stacks the publication lacks: kept, deselected
        ndeselected = 0
        do k = 1,nstks_pool
            if( l_seen(k) ) cycle
            fromp = self%spproj%os_stk%get_fromp(k)
            do iptcl = fromp,self%spproj%os_stk%get_top(k)
                if( self%in_first_set(iptcl) ) cycle
                call self%spproj%os_ptcl2D%set_state(iptcl, 0)
                call self%spproj%os_ptcl3D%set_state(iptcl, 0)
                ndeselected = ndeselected + 1
            enddo
        enddo
        ! the new stacks: appended, with their micrographs and particles
        if( nnew > 0 )then
            pool_nmics  = self%spproj%os_mic%get_noris()
            pool_nptcls = self%spproj%os_ptcl3D%get_noris()
            if( pool_nmics == 0 )then
                call self%spproj%os_mic%new(nnew,          is_ptcl=.false.)
                call self%spproj%os_stk%new(nnew,          is_ptcl=.false.)
                call self%spproj%os_ptcl2D%new(nptcls_new, is_ptcl=.true.)
                call self%spproj%os_ptcl3D%new(nptcls_new, is_ptcl=.true.)
                fromp = 1
            else
                call self%spproj%os_mic%reallocate(pool_nmics + nnew)
                call self%spproj%os_stk%reallocate(pool_nmics + nnew)
                call self%spproj%os_ptcl2D%reallocate(pool_nptcls + nptcls_new)
                call self%spproj%os_ptcl3D%reallocate(pool_nptcls + nptcls_new)
                fromp = self%spproj%os_stk%get_top(pool_nmics) + 1
            endif
            allocate(new_names(nnew))
            imic = pool_nmics
            do jstk = 1,nstks
                if( kstk(jstk) /= 0 ) cycle
                fromp_set = set%os_stk%get_fromp(jstk)
                nptcls    = set%os_stk%get_top(jstk) - fromp_set + 1
                imic      = imic + 1
                new_names(imic - pool_nmics) = set%os_stk%get_str(jstk, 'stk')
                call self%spproj%os_mic%transfer_ori(imic, set%os_mic, jstk)
                call self%spproj%os_stk%transfer_ori(imic, set%os_stk, jstk)
                call self%spproj%os_stk%set(imic, 'fromp', fromp)
                call self%spproj%os_stk%set(imic, 'top',   fromp + nptcls - 1)
                do i = 0,nptcls - 1
                    iptcl = fromp     + i
                    jptcl = fromp_set + i
                    call self%spproj%os_ptcl2D%transfer_ori(iptcl, set%os_ptcl2D, jptcl)
                    call self%spproj%os_ptcl2D%set_stkind(iptcl, imic)
                    call self%spproj%os_ptcl3D%transfer_ori(iptcl, set%os_ptcl3D, jptcl)
                    call self%spproj%os_ptcl3D%set_stkind(iptcl, imic)
                enddo
                fromp = fromp + nptcls
            enddo
            self%stk_names = [self%stk_names, new_names]
            deallocate(new_names)
        endif
        ! the publication's classes, and its optics table (the newest map's, whose group ids
        ! stay those of earlier maps)
        self%spproj%os_cls2D = set%os_cls2D
        if( set%os_optics%get_noris() > 0 ) self%spproj%os_optics = set%os_optics
        self%nptcls_selected = self%spproj%os_ptcl3D%count_state_gt_zero()
        write(logfhandle,'(A,I6,A,I8,A,I8,A,I8,A,I8,A)') '>>> PUBLICATION ', id, ': ', nupdated, ' PARTICLES UPDATED, ',&
            &nptcls_new, ' APPENDED, ', ndeselected, ' DESELECTED (STACK NOT PUBLISHED); ', self%nptcls_selected, ' SELECTED'
        deallocate(kstk, l_seen)
    end subroutine merge_publication

    ! The stage's stack named @p name, 0 when it has none. The search starts at @p hint and wraps
    ! around: a publication lists the stacks in the order the stage holds them, so the next match
    ! is usually the next stack.
    integer function stack_index( self, name, hint )
        class(stream_stage_solve3D), intent(in) :: self
        class(string),               intent(in) :: name
        integer,                     intent(in) :: hint
        integer :: n, i, k
        stack_index = 0
        n = size(self%stk_names)
        if( n == 0 ) return
        do i = 0,n - 1
            k = mod(max(hint, 1) - 1 + i, n) + 1
            if( self%stk_names(k) == name )then
                stack_index = k
                return
            endif
        enddo
    end function stack_index

    ! The class averages and FRCs of publication @p set, the stage's classes since merge_publication,
    ! copied into its quality folder (QUALITY_DIR/<stem>) with their halves, and registered in the
    ! stage's out segment in place of the earlier ones; the volumes and FSCs stay. The route's
    ! solve2D and the fallback's solve3D_cavgs read them, and an assignment matches them when the
    ! expansion fails. The copies outlive pool 2D's retention of its publications, and are pruned
    ! with the quality folder (prune_quality_dirs).
    subroutine take_cavgs( self, set, stem )
        class(stream_stage_solve3D), intent(inout) :: self
        class(sp_project),           intent(inout) :: set
        class(string),               intent(in)    :: stem
        type(string) :: stk, frcs, dir, stk_copy, frcs_copy, ext
        real         :: smpd, mskdiam
        integer      :: ncls
        call set%get_cavgs_stk(stk, ncls, smpd, fail=.false.)
        if( ncls <= 0 .or. .not. file_exists(stk) ) THROW_HARD('no class averages in the export '//stem%to_char())
        call set%get_frcs(frcs, 'frc2D', fail=.false.)
        if( .not. file_exists(frcs) ) THROW_HARD('no FRCs in the export '//stem%to_char())
        dir = string(QUALITY_DIR//'/')//stem
        call simple_mkdir(QUALITY_DIR)
        call simple_mkdir(dir)
        stk_copy  = dir//'/'//basename(stk)
        frcs_copy = dir//'/'//basename(frcs)
        call simple_copy_file(stk,  stk_copy)
        call simple_copy_file(frcs, frcs_copy)
        ! the even and odd halves beside it, which solve3D_cavgs reads
        ext = string('.')//fname2ext(stk)
        if( file_exists(add2fbody(stk, ext%to_char(), '_even')) ) call simple_copy_file(add2fbody(stk, ext%to_char(), '_even'), add2fbody(stk_copy, ext%to_char(), '_even'))
        if( file_exists(add2fbody(stk, ext%to_char(), '_odd'))  ) call simple_copy_file(add2fbody(stk, ext%to_char(), '_odd'),  add2fbody(stk_copy, ext%to_char(), '_odd'))
        mskdiam = 0.
        call set%get_mskdiam('cavg', mskdiam)
        if( mskdiam > 0. )then
            call self%spproj%add_cavgs2os_out(stk_copy, smpd, 'cavg', mskdiam=mskdiam)
        else
            call self%spproj%add_cavgs2os_out(stk_copy, smpd, 'cavg')
        endif
        call self%spproj%add_frcs2os_out(frcs_copy, 'frc2D')
    end subroutine take_cavgs

    ! A running job is checked: failed or done, the next step follows in the same pass
    ! (job_failed, job_done). Without a job, next_job decides whether the route, the fallback or
    ! an assignment starts; it starts in the pass that takes the publication it needs.
    subroutine advance_jobs( self )
        class(stream_stage_solve3D), intent(inout) :: self
        if( self%phase /= PHASE_IMPORTING .and. self%phase /= PHASE_IDLE )then
            select case(self%job%status())
                case(ASYNC_JOB_RUNNING)
                    return
                case(ASYNC_JOB_FAILED)
                    call self%job_failed()
                case(ASYNC_JOB_DONE)
                    call self%job_done()
            end select
            return
        endif
        select case(next_job(self%phase, self%n_refs(), self%l_route_due, self%l_fallback_due, self%l_assign_due,&
            &self%nptcls_selected, self%params%nstates))
            case(JOB_ROUTE)
                call self%start_solve2D_first()
            case(JOB_FALLBACK)
                call self%start_fallback()
            case(JOB_ASSIGN)
                call self%start_expand()
        end select
    end subroutine advance_jobs

    ! A failed job: the route goes on without the step (solve3D_cavgs on the first set's classes
    ! without solve2D, on the selected classes without cls_expansion); a failed solve3D_cavgs is
    ! an attempt without volumes (next_attempt); a publication whose expansion failed has its
    ! pool class averages matched; a failed matching leaves the states as they are until the next
    ! publication, and a failed refinement the assignment at the classes' poses.
    subroutine job_failed( self )
        class(stream_stage_solve3D), intent(inout) :: self
        type(string) :: logfile
        integer      :: iphase
        logfile = self%job%get_log()
        iphase  = self%phase
        call self%job%kill()
        select case(iphase)
            case(PHASE_SOLVE2D)
                write(logfhandle,'(A,A,A)') '>>> WARNING: THE FIRST SET''S SOLVE2D FAILED (SEE ', logfile%to_char(),&
                    &'); SOLVE3D_CAVGS RUNS ON ITS OWN CLASS AVERAGES'
                self%phase = PHASE_IMPORTING
                call self%start_cavgs3D()
            case(PHASE_CLSEXP)
                write(logfhandle,'(A,A,A)') '>>> WARNING: CLS_EXPANSION FAILED (SEE ', logfile%to_char(),&
                    &'); SOLVE3D_CAVGS RUNS ON THE SELECTED CLASSES'
                self%phase = PHASE_IMPORTING
                call self%start_cavgs3D()
            case(PHASE_CAVGS3D)
                write(logfhandle,'(A,A,A)') '>>> WARNING: SOLVE3D_CAVGS FAILED (SEE ', logfile%to_char(), ')'
                self%phase       = PHASE_IMPORTING
                self%try_angdist = [self%try_angdist, -1.]
                call self%next_attempt()
            case(PHASE_EXPAND)
                write(logfhandle,'(A,A,A)') '>>> WARNING: THE PUBLICATION''S CLS_EXPANSION FAILED (SEE ', logfile%to_char(),&
                    &'); ITS POOL CLASS AVERAGES ARE MATCHED'
                self%phase      = PHASE_IDLE
                self%l_expanded = .false.
                call self%start_assign(self%spproj)
            case(PHASE_ASSIGN)
                write(logfhandle,'(A,A,A)') '>>> WARNING: THE MATCHING FAILED (SEE ', logfile%to_char(),&
                    &'); THE STATES STAY AS THEY WERE UNTIL THE NEXT PUBLICATION'
                self%phase = PHASE_IDLE
            case(PHASE_REFINE)
                write(logfhandle,'(A,A,A)') '>>> WARNING: THE LOCAL REFINEMENT FAILED (SEE ', logfile%to_char(),&
                    &'); THE ASSIGNMENT STAYS AT THE CLASSES'' POSES, THE VOLUMES ARE THE PREVIOUS ONES'
                self%phase = PHASE_IDLE
        end select
    end subroutine job_failed

    ! A job done: its result is taken and the next step starts in the same pass.
    subroutine job_done( self )
        class(stream_stage_solve3D), intent(inout) :: self
        select case(self%phase)
            case(PHASE_SOLVE2D)
                call self%finish_solve2D_first()
            case(PHASE_CLSEXP)
                call self%finish_clsexp()
            case(PHASE_CAVGS3D)
                call self%finish_cavgs3D()
            case(PHASE_EXPAND)
                call self%finish_expand()
            case(PHASE_ASSIGN)
                call self%finish_assign()
            case(PHASE_REFINE)
                call self%finish_refine()
        end select
    end subroutine job_done

    ! The 3D jobs' resources on @p cline_job: their threads, their parts when @p l_parts and more
    ! than one (nparts makes a job distributed even at 1; one part runs in shared memory), and the
    ! persistent worker server.
    subroutine set_job_resources( self, cline_job, l_parts )
        class(stream_stage_solve3D), intent(inout) :: self
        class(cmdline),              intent(inout) :: cline_job
        logical,                     intent(in)    :: l_parts
        type(string) :: server_address
        call cline_job%set('nthr', self%params%nthr3D)
        if( l_parts .and. self%params%nparts3D > 1 ) call cline_job%set('nparts', self%params%nparts3D)
        call cline_job%set('worker_priority', 'high')
        server_address = self%qenv%get_persistent_worker_server_address()
        if( server_address%strlen() > 0 ) call cline_job%set('worker_server', server_address)
        if( self%qenv%get_persistent_worker_nthr() > 0 ) call cline_job%set('worker_server_nthr', self%qenv%get_persistent_worker_nthr())
    end subroutine set_job_resources

    !---------------- the route to the reference volumes ----------------

    ! The route's solve2D, in solve2D/: the first set's selected particles classified into
    ! NCLS2D_FIRST classes with the pool's mask diameter, on the 3D jobs' resources, otherwise with
    ! solve2D's defaults. Each attempt starts afresh from here.
    subroutine start_solve2D_first( self )
        class(stream_stage_solve3D), intent(inout) :: self
        type(cmdline) :: cline_job
        type(string)  :: cwd, dir
        integer       :: i, nfirst, nfirst_all
        self%nroute_tries = self%nroute_tries + 1
        self%l_route_due  = .false.
        write(logfhandle,'(A,I2,A,I2)') '>>> THE ROUTE TO THE REFERENCE VOLUMES, ATTEMPT ', self%nroute_tries, ' OF ', ROUTE_TRIES
        ! an earlier attempt's model selection deselected first-set rows for the session
        nfirst = 0
        do i = 1,self%spproj%os_ptcl3D%get_noris()
            if( self%in_first_set(i) .and. self%spproj%os_ptcl3D%get_state(i) > 0 ) nfirst = nfirst + 1
        enddo
        nfirst_all = 0
        if( allocated(self%first_set) ) nfirst_all = count(self%first_set)
        write(logfhandle,'(A,I8,A,I8,A,I8,A)') '>>> THE FIRST SET: ', nfirst, ' OF ITS ', nfirst_all,&
            &' PARTICLES SELECTED (', self%nptcls_selected, ' SELECTED IN ALL)'
        call simple_getcwd(cwd)
        dir = cwd//'/'//SOLVE2D_DIR
        call fresh_job_dir(dir, SOLVE2D_DIR) ! never the directory of a job left unfinished
        call self%spproj%write(dir//'/'//SOLVE2D_DIR//METADATA_EXT)
        call set_solve2D_route_cline(cline_job, SOLVE2D_DIR//METADATA_EXT, self%mskdiam)
        call self%set_job_resources(cline_job, .true.)
        call cline_job%printline()
        call self%record_run_publication()
        call self%job%start(self%qenv, cline_job, dir, SOLVE2D_DIR, exec_bin=string('simple_exec'))
        self%phase = PHASE_SOLVE2D
        call cline_job%kill
    end subroutine start_solve2D_first

    ! The route's solve2D is done: its classes, class averages and the particles' 2D parameters
    ! become the rows' (adopt_job_2D), the classes are selected with the pool model
    ! (select_first_set_classes), and cls_expansion starts on the selected ones. Without a usable
    ! result, solve3D_cavgs runs on the first set's own class averages.
    subroutine finish_solve2D_first( self )
        class(stream_stage_solve3D), intent(inout) :: self
        type(string) :: projfile
        projfile = self%job%get_dir()//'/'//SOLVE2D_DIR//METADATA_EXT
        call self%job%kill
        self%phase = PHASE_IMPORTING
        if( .not. self%adopt_job_2D(projfile) )then
            write(logfhandle,'(A,A,A)') '>>> WARNING: NO USABLE RESULT FROM THE FIRST SET''S SOLVE2D (', projfile%to_char(),&
                &'); SOLVE3D_CAVGS RUNS ON ITS OWN CLASS AVERAGES'
            call self%start_cavgs3D()
            return
        endif
        call self%select_first_set_classes()
        call self%write_stage_project()
        call self%start_clsexp()
    end subroutine finish_solve2D_first

    ! The first set's classes selected with the pool model, and the selection mapped to their
    ! particles (select_route_classes); the stacks of selected and rejected class averages go to
    ! solve2D/quality_selection. The particles of a rejected class are deselected for the session:
    ! they stay first-set rows, whose selection later publications do not change. When the model
    ! cannot score the classes or keeps none, every populated class is selected.
    subroutine select_first_set_classes( self )
        class(stream_stage_solve3D), intent(inout) :: self
        call select_route_classes(self%spproj, .true., self%mskdiam, string(SOLVE2D_DIR//'/'//QUALITY_DIR))
        self%nptcls_selected = self%spproj%os_ptcl3D%count_state_gt_zero()
    end subroutine select_first_set_classes

    ! The route's cls_expansion of the first set's selected classes, in cls_expansion/: the classes
    ! of the selected particles split into expansion_budget sub-classes in all (cls_expansion
    ! spreads them over groups of similar class averages by population and keeps small classes
    ! whole), their averages kernel pre-images of members (its settings are its defaults), with the mask diameter fitted to the particles' box, on the 3D jobs' resources.
    subroutine start_clsexp( self )
        class(stream_stage_solve3D), intent(inout) :: self
        type(cmdline) :: cline_job
        type(string)  :: cwd, dir
        integer       :: nsplit, nsub
        call simple_getcwd(cwd)
        dir = cwd//'/'//CLSEXP_DIR
        call fresh_job_dir(dir, CLSEXP_DIR) ! never the directory of a job left unfinished
        call self%spproj%write(dir//'/'//CLSEXP_DIR//METADATA_EXT)
        nsplit = count_active_classes(self%spproj%os_ptcl2D%get_all_asint('state'), self%spproj%os_ptcl2D%get_all_asint('class'))
        nsub   = expansion_budget(nsplit)
        write(logfhandle,'(A,I6,A,I6,A)') '>>> CLS_EXPANSION OF THE FIRST SET: ', nsplit, ' CLASSES INTO ', nsub,&
            &' SUBCLASSES IN ALL'
        call set_clsexp_cline(cline_job, CLSEXP_DIR//METADATA_EXT, nsub,&
            &fit_mskdiam(self%mskdiam, self%spproj%get_box(), self%spproj%get_smpd()))
        call self%set_job_resources(cline_job, .true.)
        call cline_job%printline()
        call self%job%start(self%qenv, cline_job, dir, CLSEXP_DIR, exec_bin=string('simple_exec'))
        self%phase = PHASE_CLSEXP
        call cline_job%kill
    end subroutine start_clsexp

    ! The route's cls_expansion is done: its subclasses, their class averages and the particles'
    ! labels become the rows' (adopt_job_2D), and solve3D_cavgs starts on them. Without a usable
    ! result, solve3D_cavgs runs on the selected classes of solve2D.
    subroutine finish_clsexp( self )
        class(stream_stage_solve3D), intent(inout) :: self
        type(string) :: projfile
        projfile = self%job%get_dir()//'/'//CLSEXP_DIR//METADATA_EXT
        call self%job%kill
        self%phase = PHASE_IMPORTING
        if( self%adopt_job_2D(projfile) )then
            write(logfhandle,'(A,I6,A)') '>>> THE FIRST SET EXPANDED INTO ', self%spproj%os_cls2D%get_noris(), ' CLASSES'
        else
            write(logfhandle,'(A,A,A)') '>>> WARNING: NO USABLE RESULT FROM CLS_EXPANSION (', projfile%to_char(),&
                &'); SOLVE3D_CAVGS RUNS ON THE SELECTED CLASSES'
        endif
        call self%write_stage_project()
        call self%start_cavgs3D()
    end subroutine finish_clsexp

    ! The 2D of a route job's project @p projfile (solve2D, cls_expansion) into the rows: each
    ! selected row's 2D parameters and selection (in both segments), the classes, and the class
    ! averages and FRCs, which solve3D_cavgs reads with the halves beside them. .false., with the
    ! rows unchanged, when the project is missing or unusable (usable_job_2D).
    logical function adopt_job_2D( self, projfile ) result( l_ok )
        class(stream_stage_solve3D), intent(inout) :: self
        class(string),               intent(in)    :: projfile
        type(sp_project) :: res
        type(string)     :: stk, frcs
        real             :: smpd, mskdiam
        integer          :: ncls, n, i, s
        l_ok = .false.
        if( .not. file_exists(projfile) ) return
        call res%read_segment('ptcl2D', projfile)
        call res%read_segment('cls2D',  projfile)
        call res%read_segment('out',    projfile)
        n = self%spproj%os_ptcl2D%get_noris()
        if( .not. usable_job_2D(res, n) )then
            call res%kill
            return
        endif
        call res%get_cavgs_stk(stk, ncls, smpd, fail=.false.)
        call res%get_frcs(frcs, 'frc2D', fail=.false.)
        do i = 1,n
            if( self%spproj%os_ptcl2D%get_state(i) <= 0 ) cycle
            s = res%os_ptcl2D%get_state(i)
            call self%spproj%os_ptcl2D%transfer_2Dparams(i, res%os_ptcl2D, i)
            call self%spproj%os_ptcl2D%set_state(i, s)
            call self%spproj%os_ptcl3D%set_state(i, s)
        enddo
        self%spproj%os_cls2D = res%os_cls2D
        mskdiam = 0.
        call res%get_mskdiam('cavg', mskdiam)
        if( mskdiam > 0. )then
            call self%spproj%add_cavgs2os_out(stk, smpd, 'cavg', mskdiam=mskdiam)
        else
            call self%spproj%add_cavgs2os_out(stk, smpd, 'cavg')
        endif
        call self%spproj%add_frcs2os_out(frcs, 'frc2D')
        self%nptcls_selected = self%spproj%os_ptcl3D%count_state_gt_zero()
        l_ok = .true.
        call res%kill
    end function adopt_job_2D

    ! The fallback: solve3D_cavgs once on the pool class averages of the publication just merged,
    ! after the route failed on every attempt.
    subroutine start_fallback( self )
        class(stream_stage_solve3D), intent(inout) :: self
        self%l_fallback_due     = .false.
        self%l_fallback_pending = .false.
        self%l_fallback_done    = .true.
        write(logfhandle,'(A)') '>>> THE FALLBACK: SOLVE3D_CAVGS ON THE PUBLICATION''S POOL CLASS AVERAGES'
        call self%start_cavgs3D()
    end subroutine start_fallback

    ! solve3D_cavgs on the stage's class averages, in solve3D_cavgs/: the route's subclasses (the
    ! kept classes of solve2D when cls_expansion gave nothing, the first set's own classes
    ! without solve2D), or the fallback's publication's pool classes, with their halves and FRCs,
    ! for nstates volumes. The quality folder of the publication is kept (record_run_publication).
    subroutine start_cavgs3D( self )
        class(stream_stage_solve3D), intent(inout) :: self
        type(cmdline) :: cline_job
        type(string)  :: cwd, dir
        call simple_getcwd(cwd)
        dir = cwd//'/'//CAVGS3D_DIR
        call fresh_job_dir(dir, CAVGS3D_DIR) ! never the directory of a job left unfinished
        call self%spproj%write(dir//'/'//CAVGS3D_DIR//METADATA_EXT)
        call set_cavgs3D_cline(cline_job, CAVGS3D_DIR//METADATA_EXT, self%params%nstates, self%mskdiam)
        call self%set_job_resources(cline_job, .true.)
        call cline_job%printline()
        call self%record_run_publication()
        call self%job%start(self%qenv, cline_job, dir, CAVGS3D_DIR, exec_bin=string('simple_exec'))
        self%phase = PHASE_CAVGS3D
        call cline_job%kill
    end subroutine start_cavgs3D

    ! solve3D_cavgs is done: its volumes are an attempt (register_attempt), checked; passed, they
    ! are the references (accept_refs), failed, the next attempt follows (next_attempt).
    subroutine finish_cavgs3D( self )
        class(stream_stage_solve3D), intent(inout) :: self
        type(sp_project) :: res
        type(string)     :: projfile
        integer          :: itry
        logical          :: l_pass
        projfile = self%job%get_dir()//'/'//CAVGS3D_DIR//METADATA_EXT
        call self%job%kill
        self%phase = PHASE_IMPORTING
        call self%send_status(string('checking the initial volumes'))
        itry   = size(self%try_angdist) + 1
        l_pass = .false.
        if( file_exists(projfile) )then
            call res%read_segment('cls3D', projfile)
            call res%read_segment('out',   projfile)
            l_pass = self%register_attempt(res, itry)
            call res%kill
        else
            write(logfhandle,'(A,A,A)') '>>> WARNING: NO PROJECT FROM SOLVE3D_CAVGS (', projfile%to_char(), ')'
            self%try_angdist = [self%try_angdist, -1.]
        endif
        if( l_pass )then
            call self%accept_refs(itry, .false.)
        else
            call self%next_attempt()
        endif
    end subroutine finish_cavgs3D

    ! The volumes of solve3D_cavgs' result @p res as attempt @p itry, checked (check_state_refs):
    ! the states that hold class averages renumbered from 1, their volumes, the JPEGs beside them,
    ! the class averages per state and the record written into state_refs/try_<itry>/; every pair
    ! of states compared after docking in both hands unless l_dock_refs is off. The attempt's
    ! checked even/odd angular distance is appended to try_angdist (-1 without volumes; huge when
    ! not measured). .true. when the attempt passes the check (refs_usable).
    logical function register_attempt( self, res, itry ) result( l_pass )
        class(stream_stage_solve3D), intent(inout) :: self
        class(sp_project),           intent(inout) :: res
        integer,                     intent(in)    :: itry
        type(refs_check) :: chk
        call check_state_refs(res, size(self%state_res), self%mskdiam, string(REFS_DIR//'/try_')//int2str(itry),&
            &'SOLVE3D_CAVGS, ATTEMPT '//int2str(itry), self%l_dock_refs, chk)
        if( .not. chk%l_vols )then
            self%try_angdist = [self%try_angdist, -1.]
        else if( chk%angdist >= 0. )then
            self%try_angdist = [self%try_angdist, chk%angdist]
        else
            self%try_angdist = [self%try_angdist, huge(chk%angdist)]
        endif
        l_pass = chk%l_pass
    end function register_attempt

    ! After an attempt that failed: the route again while it has attempts left; then the fallback
    ! on the next publication; when the pool has published its last, or after the fallback, the
    ! best attempt (accept_best).
    subroutine next_attempt( self )
        class(stream_stage_solve3D), intent(inout) :: self
        if( self%l_fallback_done )then
            call self%accept_best()
        else if( self%nroute_tries < ROUTE_TRIES )then
            self%l_route_due = .true.
            write(logfhandle,'(A)') '>>> THE ROUTE RUNS AGAIN, FROM ITS SOLVE2D'
        else if( self%l_final_seen )then
            write(logfhandle,'(A)') '>>> THE POOL PUBLISHES NO MORE: NO FALLBACK'
            call self%accept_best()
        else
            self%l_fallback_pending = .true.
            write(logfhandle,'(A)') '>>> SOLVE3D_CAVGS RUNS ON THE NEXT PUBLICATION''S POOL CLASS AVERAGES'
        endif
    end subroutine next_attempt

    ! No attempt passed the check: the one with the smallest even/odd angular distance is taken,
    ! with a warning, rather than run with no states. Without any volumes the stage assigns none.
    subroutine accept_best( self )
        class(stream_stage_solve3D), intent(inout) :: self
        integer :: ibest
        ibest = best_attempt(self%try_angdist)
        if( ibest == 0 )then
            THROW_WARN('no attempt made reference volumes; the stage assigns no states')
            return
        endif
        write(logfhandle,'(A,I3,A)') '>>> WARNING: NO ATTEMPT PASSED THE CHECK; ATTEMPT ', ibest,&
            &', WITH THE SMALLEST EVEN/ODD ANGULAR DISTANCE, GIVES THE REFERENCES'
        call self%accept_refs(ibest, .true.)
    end subroutine accept_best

    ! Attempt @p itry's volumes, the JPEGs beside them and its record, copied from
    ! state_refs/try_<itry>/ into state_refs/, are the references from here on (restore_refs):
    ! the stage's states are theirs, and the GUI gets them (send_refs). @p l_forced when they did
    ! not pass the check. When the pool's final publication is in, it is assigned at once.
    subroutine accept_refs( self, itry, l_forced )
        class(stream_stage_solve3D), intent(inout) :: self
        integer,                     intent(in)    :: itry
        logical,                     intent(in)    :: l_forced
        type(string) :: dir, src
        integer      :: s
        dir = string(REFS_DIR//'/try_')//int2str(itry)
        if( .not. self%restore_refs(dir) ) THROW_HARD('the volumes of attempt '//int2str(itry)//' are missing')
        do s = 1,self%n_refs()
            call simple_copy_file(self%ref_vols(s), string(REFS_DIR//'/')//basename(self%ref_vols(s)))
            src = dir//'/'//refine3D_reprojs_fname(s)
            if( file_exists(src) ) call simple_copy_file(src, string(REFS_DIR//'/')//refine3D_reprojs_fname(s))
            src = dir//'/'//refine3D_oris_heatmap_fname(s)
            if( file_exists(src) ) call simple_copy_file(src, string(REFS_DIR//'/')//refine3D_oris_heatmap_fname(s))
            ! the class averages the state was made from
            src = dir//'/'//cavgs_fbody(s)//MRC_EXT
            if( file_exists(src) ) call simple_copy_file(src, string(REFS_DIR//'/'//cavgs_fbody(s)//MRC_EXT))
            src = dir//'/'//cavgs_fbody(s)//JPG_EXT
            if( file_exists(src) )then
                call simple_copy_file(src, string(REFS_DIR//'/'//cavgs_fbody(s)//JPG_EXT))
                call log_jpeg('REFERENCE STATE '//int2str_pad(s, 2)//' CLASS AVERAGES', string(REFS_DIR//'/'//cavgs_fbody(s)//JPG_EXT))
            endif
        enddo
        call simple_copy_file(dir//'/'//REFS_RECORD, string(REFS_DIR//'/'//REFS_RECORD))
        if( .not. self%restore_refs(string(REFS_DIR)) ) THROW_HARD('the references could not be copied to '//REFS_DIR)
        self%params%nstates = self%n_refs()
        self%phase          = PHASE_IDLE
        if( l_forced )then
            write(logfhandle,'(A,I3,A,I3,A)') '>>> THE REFERENCES: ', self%n_refs(), ' STATES (ATTEMPT ', itry, ', NOT PASSED)'
        else
            write(logfhandle,'(A,I3,A,I3,A)') '>>> THE REFERENCES: ', self%n_refs(), ' STATES (ATTEMPT ', itry, ')'
        endif
        call self%send_refs()
        if( self%l_final_seen ) self%l_assign_due = .true.
    end subroutine accept_refs

    ! The reference volumes recorded in @p dir (REFS_RECORD), taken when every one is on disk
    ! (read_refs_record): ref_vols, ref_box and ref_smpd. .false., with none taken, otherwise.
    logical function restore_refs( self, dir ) result( l_ok )
        class(stream_stage_solve3D), intent(inout) :: self
        class(string),               intent(in)    :: dir
        type(string), allocatable :: vols(:)
        real    :: smpd
        integer :: box
        l_ok = read_refs_record(dir, size(self%state_res), vols, box, smpd)
        if( .not. l_ok ) return
        call move_alloc(vols, self%ref_vols)
        self%ref_box  = box
        self%ref_smpd = smpd
    end function restore_refs

    ! The number of reference volumes (0 before they exist).
    integer function n_refs( self )
        class(stream_stage_solve3D), intent(in) :: self
        n_refs = 0
        if( allocated(self%ref_vols) ) n_refs = size(self%ref_vols)
    end function n_refs

    !---------------- the assignment of a publication ----------------

    ! A publication's assignment starts in assign/it_<n>/: cls_expansion of its classes, in
    ! expand/, the classes of the selected particles split into expansion_budget sub-classes in all
    ! (as the route's, start_clsexp), with the mask diameter fitted to the particles' box, on the 3D
    ! jobs' resources.
    subroutine start_expand( self )
        class(stream_stage_solve3D), intent(inout) :: self
        type(cmdline) :: cline_job
        type(string)  :: cwd, dir
        integer       :: nsplit, nsub
        self%nassign      = self%nassign + 1
        self%l_assign_due = .false.
        call simple_getcwd(cwd)
        call simple_mkdir(cwd//'/'//ASSIGN_DIR)
        self%assign_dir = cwd//'/'//ASSIGN_DIR//'/it_'//int2str(self%nassign)
        call fresh_job_dir(self%assign_dir, ASSIGN_LABEL) ! never the directory of a job left unfinished
        dir = self%assign_dir//'/'//EXPAND_DIR
        call fresh_job_dir(dir, CLSEXP_DIR)
        call self%spproj%write(dir//'/'//CLSEXP_DIR//METADATA_EXT)
        nsplit = count_active_classes(self%spproj%os_ptcl2D%get_all_asint('state'), self%spproj%os_ptcl2D%get_all_asint('class'))
        nsub   = expansion_budget(nsplit)
        write(logfhandle,'(A,I6,A,I6,A,I6,A)') '>>> ASSIGNMENT ', self%nassign, ': CLS_EXPANSION OF ', nsplit, ' CLASSES INTO ',&
            &nsub, ' SUBCLASSES IN ALL'
        call set_clsexp_cline(cline_job, CLSEXP_DIR//METADATA_EXT, nsub,&
            &fit_mskdiam(self%mskdiam, self%spproj%get_box(), self%spproj%get_smpd()))
        call self%set_job_resources(cline_job, .true.)
        call cline_job%printline()
        call self%job%start(self%qenv, cline_job, dir, CLSEXP_DIR, exec_bin=string('simple_exec'))
        self%phase = PHASE_EXPAND
        call cline_job%kill
    end subroutine start_expand

    ! The publication's expansion is done: its sub-classes are matched (start_assign on its
    ! project); without a usable result, the pool's class averages are.
    subroutine finish_expand( self )
        class(stream_stage_solve3D), intent(inout) :: self
        type(sp_project) :: res
        type(string)     :: projfile
        projfile = self%job%get_dir()//'/'//CLSEXP_DIR//METADATA_EXT
        call self%job%kill
        self%phase      = PHASE_IDLE
        self%l_expanded = .false.
        if( file_exists(projfile) )then
            call res%read(projfile)
            self%l_expanded = usable_job_2D(res, self%spproj%os_ptcl2D%get_noris())
        endif
        if( self%l_expanded )then
            write(logfhandle,'(A,I6,A)') '>>> THE PUBLICATION EXPANDED INTO ', res%os_cls2D%get_noris(), ' SUBCLASSES'
            call self%start_assign(res)
        else
            write(logfhandle,'(A,A,A)') '>>> WARNING: NO USABLE RESULT FROM THE PUBLICATION''S CLS_EXPANSION (', projfile%to_char(),&
                &'); ITS POOL CLASS AVERAGES ARE MATCHED'
            call self%start_assign(self%spproj)
        endif
        call res%kill
    end subroutine finish_expand

    ! The matching, in the assignment's folder: assign_cavgs_states (simple_private_exec) on a copy
    ! of @p proj (the expansion's project, or the stage's), every class average of it against the
    ! references at LP_ASSIGN over NSPACE_ASSIGN directions per state, with the mask diameter
    ! fitted to the class averages' box, on the 3D jobs' threads.
    subroutine start_assign( self, proj )
        class(stream_stage_solve3D), intent(inout) :: self
        class(sp_project),           intent(inout) :: proj
        type(cmdline) :: cline_job
        type(string)  :: stk
        real          :: smpd, box_stk
        integer       :: ncls
        call proj%write(self%assign_dir//'/'//ASSIGN_LABEL//METADATA_EXT)
        box_stk = 0.
        call proj%get_cavgs_stk(stk, ncls, smpd, fail=.false., box=box_stk)
        write(logfhandle,'(A,I6,A,I3,A)') '>>> MATCHING ', ncls, ' CLASS AVERAGES AGAINST THE ', self%n_refs(), ' REFERENCES'
        call set_assign_cline(cline_job, ASSIGN_LABEL//METADATA_EXT, self%ref_vols, fit_mskdiam(self%mskdiam, nint(box_stk), smpd))
        call cline_job%set('nthr', self%params%nthr3D)
        call cline_job%printline()
        call self%job%start(self%qenv, cline_job, self%assign_dir, ASSIGN_LABEL, exec_bin=string('simple_private_exec'))
        self%phase = PHASE_ASSIGN
        call cline_job%kill
    end subroutine start_assign

    ! The matching is done: the rejection, the flags and the mapping onto the rows
    ! (take_assignment), recorded in the matched project, the (sub-)class averages grouped by state
    ! (write_assignment_cavgs), then the local refinement. Without a usable result the states stay
    ! as they were.
    subroutine finish_assign( self )
        class(stream_stage_solve3D), intent(inout) :: self
        type(sp_project)     :: res
        type(string)         :: projfile
        integer, allocatable :: states_matched(:)
        projfile = self%assign_dir//'/'//ASSIGN_LABEL//METADATA_EXT
        call self%job%kill
        self%phase = PHASE_IDLE
        if( .not. file_exists(projfile) )then
            write(logfhandle,'(A,A,A)') '>>> WARNING: NO PROJECT FROM THE MATCHING (', projfile%to_char(), '); THE STATES STAY'
            return
        endif
        call res%read(projfile)
        ! the states the matching gave, before the rejection
        if( res%os_cls3D%get_noris() > 0 )then
            states_matched = res%os_cls3D%get_all_asint('state')
        else
            allocate(states_matched(0))
        endif
        if( .not. self%take_assignment(res) )then
            write(logfhandle,'(A,A,A)') '>>> WARNING: NO USABLE MATCHING IN ', projfile%to_char(), '; THE STATES STAY'
            call res%kill
            return
        endif
        call res%write(projfile)
        call self%write_assignment_cavgs(res, states_matched)
        call res%kill
        call self%write_stage_project()
        call self%start_refine()
    end subroutine finish_assign

    ! The matching in @p res (a copy of the rows, with the matched classes' cls2D, cls3D and the
    ! particles' class labels) onto the rows (map_assignment): the outliers of the fit rejected,
    ! the small margins flagged, each particle given its (sub-)class's state and pose, in both
    ! particle segments. A row deselected here is recorded (fit_desel): one of the first set is
    ! selected again by the next merge. The counts are logged beside the pool model's verdict
    ! (model_rej). .false., with the rows unchanged, when @p res holds other rows or no matching.
    logical function take_assignment( self, res ) result( l_ok )
        class(stream_stage_solve3D), intent(inout) :: self
        class(sp_project),           intent(inout) :: res
        type(assign_stats) :: stats
        if( self%l_expanded )then
            l_ok = map_assignment(self%spproj, res, self%n_refs(), 'SUBCLASSES', stats, self%fit_desel, self%model_rej)
        else
            l_ok = map_assignment(self%spproj, res, self%n_refs(), 'POOL CLASSES', stats, self%fit_desel, self%model_rej)
        endif
        if( .not. l_ok ) return
        self%nptcls_selected    = self%spproj%os_ptcl3D%count_state_gt_zero()
        self%nptcls_at_last_run = stats%nassigned
    end function take_assignment

    ! The matched (sub-)class averages of @p res grouped into the assignment's folder
    ! (write_matched_cavgs): a stack and a JPEG per state (the ambiguous ones included), one of
    ! those the fit rejected (a state given by the matching, @p states_matched, and none after the
    ! rejection) and one of the ambiguous ones. None when the class averages are not on disk.
    subroutine write_assignment_cavgs( self, res, states_matched )
        class(stream_stage_solve3D), intent(inout) :: self
        class(sp_project),           intent(inout) :: res
        integer,                     intent(in)    :: states_matched(:)
        call write_matched_cavgs(res, states_matched, self%n_refs(), self%assign_dir, 'ASSIGNMENT '//int2str(self%nassign))
    end subroutine write_assignment_cavgs

    ! The seeds of the local refinement on the rows (set_refine_seeds; refine=cont continues from a
    ! discrete pose): the half-sets of the project's partition (eo_transaction: one half per
    ! micrograph, never changed, so every assignment has the same halves) and the projection
    ! direction closest to each mapped pose in the job's grid (nspace directions, c1).
    subroutine prepare_refine_seeds( self )
        class(stream_stage_solve3D), intent(inout) :: self
        call set_refine_seeds(self%spproj, self%params%nspace)
    end subroutine prepare_refine_seeds

    ! The local refinement, in the assignment's refine3D/: one iteration (maxits=1) of refine3D
    ! refine=cont over every selected row from its mapped pose (prepare_refine_seeds), each
    ! particle's pose refined continuously against its own state, which stays, within
    ! refine=cont's bounds (athres_cont, trs); the iteration ends with the states reconstructed
    ! from the refined poses. The starting volumes are reconstructed from the mapped poses; the
    ! low-pass limits come from the FSC; no sigma2 (objfun=cc). Distributed
    ! even with one part: the distributed refine3D reconstructs its starting volumes, the
    ! shared-memory one needs them given. The job holds the states the particles hold,
    ! renumbered from 1 in its project only (write_refine_project), since refine3D reconstructs no
    ! volume for an empty state; take_refine_result maps them back.
    subroutine start_refine( self )
        class(stream_stage_solve3D), intent(inout) :: self
        type(cmdline) :: cline_job
        type(string)  :: dir
        integer       :: nheld
        call self%prepare_refine_seeds()
        if( self%spproj%os_ptcl3D%count_state_gt_zero() < 1 )then
            write(logfhandle,'(A)') '>>> WARNING: NO PARTICLE HAS A STATE; NO LOCAL REFINEMENT'
            return
        endif
        dir = self%assign_dir//'/'//REFINE_DIR
        call fresh_job_dir(dir, REFINE_DIR) ! never the directory of a job left unfinished
        call write_refine_project(self%spproj, dir//'/'//REFINE_DIR//METADATA_EXT, nheld)
        write(logfhandle,'(A,I8,A,I3,A)') '>>> REFINING ', self%nptcls_selected, ' PARTICLES LOCALLY IN ', nheld,&
            &' FIXED STATES'
        call set_refine_cont_cline(cline_job, REFINE_DIR//METADATA_EXT, nheld, self%params%nspace, self%mskdiam)
        call self%set_job_resources(cline_job, .false.)
        call cline_job%set('nparts',   self%params%nparts3D)
        call cline_job%printline()
        call self%job%start(self%qenv, cline_job, dir, REFINE_DIR, exec_bin=string('simple_exec'))
        self%phase = PHASE_REFINE
        call cline_job%kill
    end subroutine start_refine

    ! The local refinement is done: its volumes, FSCs and the refined poses are the stage's
    ! (take_refine_result), the stage's project is written beside it as the snapshots' source
    ! (result_projfile, RESULT_PROJ in the assignment's folder; refine3D's own project is kept as
    ! the job left it), and the GUI gets the states (send_volumes). Older assignment folders go
    ! (prune_run_dirs).
    subroutine finish_refine( self )
        class(stream_stage_solve3D), intent(inout) :: self
        type(sp_project) :: res
        type(string)     :: dir, projfile
        integer          :: nreg
        dir      = self%job%get_dir()
        projfile = dir//'/'//REFINE_DIR//METADATA_EXT
        call self%job%kill
        self%phase = PHASE_IDLE
        nreg = 0
        if( file_exists(projfile) )then
            call res%read_segment('ptcl3D', projfile)
            call res%read_segment('out',    projfile)
            nreg = self%take_refine_result(res)
            call res%kill
        endif
        if( nreg == 0 )then
            write(logfhandle,'(A,A,A)') '>>> WARNING: THE LOCAL REFINEMENT IN ', dir%to_char(), ' GAVE NO USABLE RESULT'
            return
        endif
        self%result_projfile = self%assign_dir//'/'//RESULT_PROJ//METADATA_EXT
        call self%spproj%write(self%result_projfile)
        call self%send_volumes()
        call self%write_stage_project()
        call self%prune_run_dirs()
    end subroutine finish_refine

    ! The local refinement's project @p res onto the stage, its states mapped back to the stage's
    ! (map_refine_result: start_refine renumbered the states the rows hold): each state's volume
    ! and FSC it registered, in place of every earlier state's, and each selected row's refined
    ! pose where its state is the row's. The number of states registered; 0, with the stage
    ! unchanged, when @p res holds other rows or no volume.
    integer function take_refine_result( self, res ) result( nreg )
        class(stream_stage_solve3D), intent(inout) :: self
        class(sp_project),           intent(inout) :: res
        nreg = map_refine_result(self%spproj, res, size(self%state_res))
    end function take_refine_result

    !---------------- records ----------------

    ! The publication the route or the fallback starts from keeps its quality folder (QUALITY_RUNS).
    subroutine record_run_publication( self )
        class(stream_stage_solve3D), intent(inout) :: self
        integer :: funit, ios
        if( self%last_stem%strlen() == 0 ) return
        call simple_mkdir(QUALITY_DIR)
        open(newunit=funit, file=QUALITY_RUNS, position='append', action='write', iostat=ios)
        if( ios /= 0 ) return
        write(funit,'(A)') self%last_stem%to_char()
        close(funit)
    end subroutine record_run_publication

    ! The quality folders kept are the newest NQUALITY_KEPT and those of the publications the
    ! route or the fallback started from; the others go.
    subroutine prune_quality_dirs( self )
        class(stream_stage_solve3D), intent(inout) :: self
        type(string), allocatable :: dirs(:), runs(:)
        character(len=256) :: line
        integer :: i, j, funit, ios, nruns
        logical :: l_run
        if( .not. dir_exists(QUALITY_DIR) ) return
        dirs = simple_list_dirs(QUALITY_DIR)
        if( size(dirs) <= NQUALITY_KEPT ) return
        call lex_sort(dirs) ! 5-digit publication ids: the newest last
        ! the publications that started a run, counted then read: growing an array of strings by
        ! a constructor that holds the array (runs = [runs, ...]) crashed some gfortran builds
        nruns = 0
        open(newunit=funit, file=QUALITY_RUNS, status='old', action='read', iostat=ios)
        if( ios == 0 )then
            do
                read(funit,'(A)',iostat=ios) line
                if( ios /= 0 ) exit
                nruns = nruns + 1
            enddo
            allocate(runs(nruns))
            rewind(funit)
            do j = 1,nruns
                read(funit,'(A)',iostat=ios) line
                if( ios /= 0 )then
                    nruns = j - 1
                    exit
                endif
                runs(j) = trim(line)
            enddo
            close(funit)
        else
            allocate(runs(0))
        endif
        do i = 1,size(dirs) - NQUALITY_KEPT
            l_run = .false.
            do j = 1,nruns
                if( runs(j) == dirs(i)%to_char() ) l_run = .true.
            enddo
            if( l_run ) cycle
            call simple_rmdir(QUALITY_DIR//'/'//dirs(i)%to_char())
        enddo
    end subroutine prune_quality_dirs

    ! Once an assignment is refined, the assignment folders before it and the folders set aside
    ! for jobs left unfinished go.
    subroutine prune_run_dirs( self )
        class(stream_stage_solve3D), intent(inout) :: self
        type(string), allocatable :: dirs(:)
        integer :: i, k
        do k = 1,self%nassign - 1
            if( dir_exists(ASSIGN_DIR//'/it_'//int2str(k)) ) call simple_rmdir(ASSIGN_DIR//'/it_'//int2str(k))
        enddo
        dirs = simple_list_dirs('.')
        do i = 1,size(dirs)
            if( dirs(i)%has_substr(CAVGS3D_DIR//'_unfinished') ) call simple_rmdir(dirs(i)%to_char())
            if( dirs(i)%has_substr(SOLVE2D_DIR//'_unfinished') ) call simple_rmdir(dirs(i)%to_char())
            if( dirs(i)%has_substr(CLSEXP_DIR//'_unfinished')  ) call simple_rmdir(dirs(i)%to_char())
        enddo
        if( dir_exists(ASSIGN_DIR) )then
            dirs = simple_list_dirs(ASSIGN_DIR)
            do i = 1,size(dirs)
                if( dirs(i)%has_substr('_unfinished') ) call simple_rmdir(ASSIGN_DIR//'/'//dirs(i)%to_char())
            enddo
        endif
    end subroutine prune_run_dirs

    ! The rows as the stage's project, written as a temporary file and renamed.
    subroutine write_stage_project( self )
        class(stream_stage_solve3D), intent(inout) :: self
        call self%spproj%write(self%params%projfile, tempfile=.true.)
    end subroutine write_stage_project

    !---------------- GUI ----------------

    ! The GUI's updates, drained once a pass: a 3D snapshot request is written (write_snapshot);
    ! the fields meant for other stages are ignored.
    subroutine apply_gui_updates( self )
        class(stream_stage_solve3D), intent(inout) :: self
        type(gui_metadata_stream_update) :: update
        character(len=:), allocatable    :: buffer
        do while( self%pipe%receive(buffer) )
            update = transfer(buffer, update)
            if( update%has_snapshot3D_update() ) call self%write_snapshot(update)
        enddo
    end subroutine apply_gui_updates

    ! A 3D snapshot the GUI asked for, from the latest local refinement's project: the particles of
    ! the selected states merged into state 1 (the others deselected, and the state volumes and
    ! FSCs out of the project), written with micrograph and particle STAR files and the optics
    ! table to snapshots/<name>/, beside copies of the selected states' volumes
    ! (vol_state<NN>.mrc). Each request is written once. Before a refinement, or with no
    ! particle in the selected states, it is reported with no particles and no file.
    subroutine write_snapshot( self, update )
        class(stream_stage_solve3D),      intent(inout) :: self
        type(gui_metadata_stream_update), intent(in)    :: update
        type(sp_project)     :: snap
        integer, allocatable :: selection(:), states(:)
        type(string)         :: filename, stem, cwd, dir, projfile, volpath
        real                 :: smpd
        integer              :: snapshot_id, nptcls, iptcl, istate, box
        call update%get_snapshot3D_update(snapshot_id, selection, filename)
        if( snapshot_id <= self%last_snapshot_id ) return
        self%last_snapshot_id = snapshot_id
        nptcls   = 0
        projfile = string('')
        if( self%result_projfile%strlen() == 0 .or. .not. file_exists(self%result_projfile) )then
            write(logfhandle,'(A,I6,A)') '>>> 3D SNAPSHOT ', snapshot_id, ' REQUESTED BEFORE A 3D RESULT; NOT WRITTEN'
        else
            call snap%read(self%result_projfile)
            ! the selected states' particles, merged into state 1
            states = snap%os_ptcl3D%get_all_asint('state')
            do iptcl = 1,size(states)
                if( any(selection == states(iptcl)) )then
                    states(iptcl) = 1
                else
                    states(iptcl) = 0
                endif
                call snap%os_ptcl3D%set_state(iptcl, states(iptcl))
            enddo
            if( snap%os_ptcl2D%get_noris() == size(states) )then
                do iptcl = 1,size(states)
                    call snap%os_ptcl2D%set_state(iptcl, states(iptcl))
                enddo
            endif
            nptcls = count(states == 1)
            if( nptcls == 0 )then
                write(logfhandle,'(A,I6,A)') '>>> 3D SNAPSHOT ', snapshot_id, ': NO PARTICLE IN THE SELECTED STATES; NOT WRITTEN'
            else
                call simple_getcwd(cwd)
                stem = swap_suffix(filename, '', METADATA_EXT)
                dir  = cwd//'/'//DIR_SNAPSHOT//stem
                if( .not. dir_exists(cwd//'/'//DIR_SNAPSHOT) ) call simple_mkdir(cwd//'/'//DIR_SNAPSHOT)
                if( .not. dir_exists(dir) )                    call simple_mkdir(dir)
                ! the selected states' volumes, copied beside the project, which has one merged state
                do istate = 1,self%params%nstates
                    if( .not. snap%isthere_in_osout('vol', istate) ) cycle
                    if( any(selection == istate) )then
                        call snap%get_vol('vol', istate, volpath, smpd, box)
                        if( file_exists(volpath) ) call simple_copy_file(volpath, dir//'/vol_state'//int2str_pad(istate,2)//MRC_EXT)
                    endif
                enddo
                do istate = 1,self%params%nstates
                    call snap%remove_state_artifacts_from_osout(istate)
                enddo
                projfile = dir//'/'//filename
                write(logfhandle,'(A,I6,A,I8,A,A)') '>>> WRITING 3D SNAPSHOT ', snapshot_id, ' OF ', nptcls, ' PARTICLES: ',&
                    &projfile%to_char()
                call snap%write(projfile)
                call snap%write_mics_star(dir//'/'//stem//'_micrographs.star', optics_offset=self%optics_id_offset)
                call snap%write_ptcl2D_star(dir//'/'//stem//'_particles.star', optics_offset=self%optics_id_offset)
            endif
            call snap%kill
        endif
        ! a snapshot that could not be written goes with no particles and no file
        if( nptcls > 0 )then
            call self%meta_snapshot%set(id=snapshot_id, snapshot_filename=projfile, snapshot_nptcls=nptcls, states=selection)
        else
            call self%meta_snapshot%set(id=snapshot_id, snapshot_filename=string(''), snapshot_nptcls=0, states=selection)
        endif
        call self%pipe%send_meta(self%meta_snapshot)
    end subroutine write_snapshot

    ! The phase, the assignments, the particle counts and, once the references exist, each
    ! state's population and resolution (of the latest refinement; none for the references).
    subroutine send_status( self, stage )
        class(stream_stage_solve3D), intent(inout) :: self
        type(string), optional,      intent(in)    :: stage
        type(string) :: stage_here
        integer      :: istate, progress
        select case(self%phase)
            case(PHASE_IMPORTING)
                stage_here = 'waiting for classified particles'
                if( .not. self%l_attached ) stage_here = 'waiting on pool 2D'
            case(PHASE_SOLVE2D)
                stage_here = 'running solve2D on the first set'
            case(PHASE_CLSEXP)
                stage_here = 'running cls_expansion on the first set'
            case(PHASE_CAVGS3D)
                stage_here = 'running solve3D_cavgs'
            case(PHASE_EXPAND)
                stage_here = 'expanding the classes'
            case(PHASE_ASSIGN)
                stage_here = 'assigning class averages'
            case(PHASE_REFINE)
                stage_here = 'refining the poses locally'
            case default
                stage_here = 'idle'
        end select
        if( present(stage) ) stage_here = stage
        ! the GUI's 3D progress: 0 none yet, 1 the reference volumes on their way, 2 they exist
        progress = 0
        if( self%n_refs() > 0 )then
            progress = 2
        else if( self%phase == PHASE_SOLVE2D .or. self%phase == PHASE_CLSEXP .or. self%phase == PHASE_CAVGS3D )then
            progress = 1
        endif
        ! particles_imported is the particles selected now; particles_at_last_refine those the
        ! latest assignment gave a state
        call self%meta_status%set(stage=stage_here, solve3D_stage=progress,&
            &refine_iteration=self%nassign, nstates=self%params%nstates,&
            &particles_imported=self%spproj%os_ptcl3D%count_state_gt_zero(), particles_at_last_refine=self%nptcls_at_last_run,&
            &resolution=0.)
        ! each state's population and resolution as sent with its volume (send_volumes, send_refs)
        if( self%n_refs() > 0 .and. allocated(self%state_res) .and. allocated(self%state_pop) )then
            do istate = 1,min(self%params%nstates, size(self%state_res))
                call self%meta_status%set_state_stats(istate, self%state_pop(istate), self%state_res(istate))
            enddo
        endif
        call self%pipe%send_meta(self%meta_status)
    end subroutine send_status

    ! Each state with a volume in the stage's project (send_state_volume, with its FSC). The
    ! resolutions and populations are kept for the status. Called on a local refinement's result,
    ! before any particle merged since joins a state.
    subroutine send_volumes( self )
        class(stream_stage_solve3D), intent(inout) :: self
        type(string) :: volpath
        real         :: smpd
        integer      :: istate, box
        self%state_res = 0.
        self%state_pop = 0
        do istate = 1,self%params%nstates
            if( .not. self%spproj%isthere_in_osout('vol', istate) ) cycle
            call self%spproj%get_vol('vol', istate, volpath, smpd, box)
            if( volpath%strlen() == 0 ) cycle
            call self%send_state_volume(istate, volpath, smpd, box, istate, .true.)
        enddo
    end subroutine send_volumes

    ! The reference volumes (state_refs/) to the GUI, without FSCs, populations or orientations
    ! (send_state_volume): they come from class averages, not from the particles.
    subroutine send_refs( self )
        class(stream_stage_solve3D), intent(inout) :: self
        integer :: s, nsent
        self%state_res = 0.
        self%state_pop = 0
        nsent = 0
        do s = 1,self%n_refs()
            if( .not. file_exists(self%ref_vols(s)) ) cycle
            call self%send_state_volume(s, self%ref_vols(s), self%ref_smpd, self%ref_box, s, .false., l_rows=.false.)
            nsent = nsent + 1
        enddo
        write(logfhandle,'(A,I3,A)') '>>> ', nsent, ' REFERENCE VOLUME(S) SENT TO THE GUI'
    end subroutine send_refs

    ! State @p istate's volume @p volpath (pixel size @p smpd, box @p box) to the GUI: its products,
    ! the minimum and maximum of each volume, its three reprojections as tiles, with @p l_with_fsc
    ! its resolution and FSC curve, and unless @p l_rows is .false. its population and the
    ! orientation distribution of its particles. The reprojection and heatmap JPEGs beside the
    ! volume are named after the state @p ifile, which may differ from the stage's.
    subroutine send_state_volume( self, istate, volpath, smpd, box, ifile, l_with_fsc, l_rows )
        class(stream_stage_solve3D), intent(inout) :: self
        integer,                     intent(in)    :: istate, box, ifile
        type(string),                intent(in)    :: volpath
        real,                        intent(in)    :: smpd
        logical,                     intent(in)    :: l_with_fsc
        logical, optional,           intent(in)    :: l_rows
        type(gui_metadata_vol3D) :: meta, fresh_meta
        type(string)             :: pprocpath, lppath, pprocmirrpath, reprojpath, oridistpath
        real, allocatable        :: fsc(:), res(:)
        real                     :: res05, res0143, vmin, vmax
        integer                  :: pop, n
        integer                  :: hist(ORIDIST_NBINS_X, ORIDIST_NBINS_Y)
        logical                  :: l_fsc, l_rows_here
        l_rows_here = .true.
        if( present(l_rows) ) l_rows_here = l_rows
        pop = 0
        if( l_rows_here ) pop = self%spproj%os_ptcl3D%get_pop(istate, 'state')
        self%state_pop(istate) = pop
        ! products of a postprocessing, when on disk (none runs in this stage)
        pprocpath = add2fbody(volpath, MRC_EXT, PPROC_SUFFIX)
        if( .not. file_exists(pprocpath) ) pprocpath = ''
        lppath = add2fbody(volpath, MRC_EXT, LP_SUFFIX)
        if( .not. file_exists(lppath) ) lppath = ''
        pprocmirrpath = ''
        if( pprocpath%strlen() > 0 )then
            pprocmirrpath = add2fbody(pprocpath, MRC_EXT, MIRR_SUFFIX)
            if( .not. file_exists(pprocmirrpath) ) pprocmirrpath = ''
        endif
        ! JPEGs beside the volume
        reprojpath = get_fpath(volpath)//refine3D_reprojs_fname(ifile)
        if( .not. file_exists(reprojpath) ) reprojpath = ''
        oridistpath = get_fpath(volpath)//refine3D_oris_heatmap_fname(ifile)
        if( .not. file_exists(oridistpath) ) oridistpath = ''
        if( reprojpath%strlen() > 0 ) call send_reproj_tiles(self%pipe, self%meta_reproj, reprojpath, volpath,&
            &istate, self%params%nstates, pop)
        l_fsc = .false.
        if( l_with_fsc ) l_fsc = self%read_state_fsc(istate, smpd, fsc, res, res05, res0143)
        ! every field back to its default: a state without an FSC curve or a product must not
        ! carry the previous state's
        meta = fresh_meta
        call meta%new(GUI_METADATA_VOL3D_TYPE)
        if( l_fsc )then
            self%state_res(istate) = res0143
            call meta%set(reprojpath, volpath, pprocpath, lppath, pprocmirrpath, istate, box, smpd, istate,&
                &self%params%nstates, res0143=res0143, res05=res05, pop=pop, oridistpath=oridistpath)
            n = min(size(fsc), size(res), MAX_FSC_VOL3D)
            call meta%set_fsc(1. / res(:n), fsc(:n))
        else
            call meta%set(reprojpath, volpath, pprocpath, lppath, pprocmirrpath, istate, box, smpd, istate,&
                &self%params%nstates, pop=pop, oridistpath=oridistpath)
        endif
        ! header minimum and maximum, so the GUI need not open the volumes
        call get_mrc_minmax(volpath, vmin, vmax)
        call meta%set_minmax('volpath', vmin, vmax)
        if( pprocpath%strlen() > 0 )then
            call get_mrc_minmax(pprocpath, vmin, vmax)
            call meta%set_minmax('pprocpath', vmin, vmax)
        endif
        if( lppath%strlen() > 0 )then
            call get_mrc_minmax(lppath, vmin, vmax)
            call meta%set_minmax('lppath', vmin, vmax)
        endif
        if( pprocmirrpath%strlen() > 0 )then
            call get_mrc_minmax(pprocmirrpath, vmin, vmax)
            call meta%set_minmax('pprocmirrpath', vmin, vmax)
        endif
        if( l_rows_here )then
            call oridist_from_oris(self%spproj%os_ptcl3D, istate, hist)
            call meta%set_oridist(hist)
        endif
        call self%pipe%send_meta(meta)
        call meta%kill
    end subroutine send_state_volume

    ! The FSC of state @p istate in the stage's project, its resolutions at 1/(Fourier index) for
    ! pixel size @p smpd, and the FSC=0.5 and FSC=0.143 resolutions; .false. without an FSC file.
    logical function read_state_fsc( self, istate, smpd, fsc, res, res05, res0143 ) result( l_fsc )
        class(stream_stage_solve3D), intent(in)    :: self
        integer,                     intent(in)    :: istate
        real,                        intent(in)    :: smpd
        real, allocatable,           intent(inout) :: fsc(:), res(:)
        real,                        intent(out)   :: res05, res0143
        type(string) :: fsc_fname
        integer      :: fsc_box, n
        l_fsc   = .false.
        res05   = 0.
        res0143 = 0.
        if( .not. self%spproj%isthere_in_osout('fsc', istate) ) return
        call self%spproj%get_fsc(istate, fsc_fname, fsc_box)
        if( fsc_fname%strlen() == 0 ) return
        if( .not. file_exists(fsc_fname) ) return
        fsc = file2rarr(fsc_fname)
        res = get_resarr(fsc_box, smpd)
        n   = min(size(fsc), size(res))
        if( n == 0 ) return
        fsc = fsc(:n)
        res = res(:n)
        call get_resolution(fsc, res, res05, res0143)
        l_fsc = .true.
    end function read_state_fsc

    !---------------- rules ----------------

    !> The job to start in @p phase (none while a job runs). Without references (@p nrefs 0), once
    !! @p nselected particles are selected, MIN_PTCLS_PER_STATE for each of @p nstates: the route
    !! when @p l_route_due, else the fallback when @p l_fallback_due. With references, a
    !! publication's assignment when @p l_assign_due and MIN_PTCLS_PER_STATE particles per
    !! reference are selected.
    pure integer function next_job( phase, nrefs, l_route_due, l_fallback_due, l_assign_due, nselected, nstates )
        integer, intent(in) :: phase, nrefs, nselected, nstates
        logical, intent(in) :: l_route_due, l_fallback_due, l_assign_due
        next_job = JOB_NONE
        if( phase /= PHASE_IMPORTING .and. phase /= PHASE_IDLE ) return
        if( nrefs < 1 )then
            if( nselected < MIN_PTCLS_PER_STATE * nstates ) return
            if( l_route_due )then
                next_job = JOB_ROUTE
            else if( l_fallback_due )then
                next_job = JOB_FALLBACK
            endif
        else if( l_assign_due .and. nselected >= MIN_PTCLS_PER_STATE * nrefs )then
            next_job = JOB_ASSIGN
        endif
    end function next_job

    !> The attempt with the smallest even/odd angular distance @p angdists among those with
    !! volumes (not negative), 0 when there is none; the first of equals.
    pure integer function best_attempt( angdists )
        real, intent(in) :: angdists(:)
        integer :: i
        best_attempt = 0
        do i = 1,size(angdists)
            if( angdists(i) < 0. ) cycle
            if( best_attempt == 0 )then
                best_attempt = i
            else if( angdists(i) < angdists(best_attempt) )then
                best_attempt = i
            endif
        enddo
    end function best_attempt

    !> .true. for the pool's final publication (pool_final=yes in an entry of its out segment).
    logical function is_final_publication( set )
        class(sp_project), intent(inout) :: set
        type(string) :: val
        integer      :: i
        is_final_publication = .false.
        do i = 1,set%os_out%get_noris()
            if( .not. set%os_out%isthere(i, 'pool_final') ) cycle
            val = set%os_out%get_str(i, 'pool_final')
            if( val == 'yes' )then
                is_final_publication = .true.
                return
            endif
        enddo
    end function is_final_publication

end module simple_stream_stage_solve3D
