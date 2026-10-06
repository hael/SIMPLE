!@descr: state and steps of stream task 7 (multistate 3D): import the particles pool 2D exports, run solve3D and then solve3D_addon as they grow, report to the GUI
!==============================================================================
! MODULE: simple_stream_stage_solve3D
!
! PURPOSE:
!   The body of stream p07 as a type. Pool 2D publishes its classified state
!   after each completed iteration (doc/policies/stream/stream_3D_ingestion_policy.md).
!   Each pass in which no job runs takes the newest publication. The first one
!   (in a fresh pool, the one after iteration 10, with the sieve's mask
!   diameter) is the first set: its particles are classified again by solve2D
!   and selected with the model and compatibility filter p03 uses, and that
!   selection stays theirs for the session. Every later publication's class
!   averages are selected once with the pool model, and the selection and 2D
!   parameters are merged into the stage's rows (the first set's rows take only
!   the 2D parameters). Its class averages and FRCs, copied into the stage's
!   folder, come with its classes. The selected particles start solve3D; once it is
!   done, every growth of the rows starts an solve3D_addon run from the
!   latest result. An addon run whose verdict has a REGRESSED state is rolled
!   back: the previous result stays. Once the pool's final publication is in
!   (pool_final), one multistate refine3D realigns every selected particle from
!   the current state volumes. After each run the stage's project holds the
!   result, and the GUI gets each state's volume, resolution, reprojections and
!   orientation distribution. The commander (simple_commanders_stream_p07_solve3D_multistate)
!   only normalises the command line and loops over iterate() until finished().
!
!   What happens is delegated:
!     - class-average selection -> simple_cavg_quality_selection, simple_class_compatibility
!     - jobs                    -> simple_qsys_async_job (solve2D, solve3D, solve3D_addon)
!     - GUI                     -> simple_stream_pipe, simple_stream_gui_senders,
!                                  simple_oris_utils (the orientation histograms)
!
!   The stage's rows only ever grow, as solve3D_addon requires (rows are
!   read by index, never renumbered): a publication's stacks are matched to
!   the rows by stack name, so the order of a restarted pool does not matter.
!
! LIFECYCLE:
!   new(cline) -> { iterate() } until finished() -> finalize() -> kill()
!
!   A publication the stage cannot use (its stacks disagree with the rows, its
!   class averages or FRCs are missing) is passed over with a warning and listed in
!   REJECTED_PUBLICATIONS, and the stage waits for the next one.
!
! RESTART:
!   A stop cancels the running job (finalize). A restart removes a leftover
!   TERM_STREAM and starts again from the newest publication not rejected,
!   which is then the first set; solve2D and solve3D run again in their
!   folders, which are first moved aside when they hold a job left unfinished
!   (fresh_job_dir), as is an addon iteration folder.
!==============================================================================
module simple_stream_stage_solve3D
use simple_defs,                                      only: logfhandle, COSMSKHALFWIDTH
use simple_defs_fname,                                only: TERM_STREAM, METADATA_EXT, MRC_EXT, JPG_EXT, PPROC_SUFFIX,&
                                                           &LP_SUFFIX, MIRR_SUFFIX, DIR_SNAPSHOT
use simple_defs_stream,                               only: DIR_STREAM_COMPLETED, SHORTWAIT, WAITTIME, OPTICS_ID_DELTA
use simple_defs_environment,                          only: SIMPLE_STREAM_SOLVE3D_PARTITION
use simple_error,                                     only: simple_exception
use simple_string,                                    only: string
use simple_string_utils,                              only: int2str, int2str_pad, lex_sort
use simple_fileio,                                    only: add2fbody, basename, del_file, file2rarr, file_exists, get_fbody,&
                                                           &get_fpath, simple_abspath, simple_getcwd, swap_suffix, simple_copy_file
use simple_syslib,                                    only: dir_exists, simple_mkdir, simple_list_dirs, simple_rmdir
use simple_math,                                      only: round2even
use simple_math_ft,                                   only: get_resarr
use simple_estimate_ssnr,                             only: get_resolution
use simple_imghead,                                   only: get_mrc_minmax, find_ldim_nptcls
use simple_refine3D_fnames,                           only: refine3D_oris_heatmap_fname, refine3D_reprojs_fname
use simple_oris_utils,                                only: oridist_from_oris
use simple_cmdline,                                   only: cmdline
use simple_parameters,                                only: parameters
use simple_sp_project,                                only: sp_project
use simple_image,                                     only: image
use simple_imgarr_utils,                              only: dealloc_imgarr
use simple_qsys_env,                                  only: qsys_env
use simple_qsys_funs,                                 only: qsys_cleanup
use simple_qsys_async_job,                            only: qsys_async_job, ASYNC_JOB_RUNNING, ASYNC_JOB_DONE, ASYNC_JOB_FAILED
use simple_qsys_job_record,                           only: fresh_job_dir
use simple_rec_list,                                  only: rec_list, rec_iterator, chunk_rec
use simple_stream_watcher,                            only: stream_watcher
use simple_stream_state,                              only: ipc_pipe_solve3D_multistate_in, ipc_pipe_solve3D_multistate_out
use simple_stream_utils,                              only: create_stream_project, init_stream_qenv
use simple_cavg_quality_model,                        only: cavg_quality_model, CAVG_QUALITY_MODEL_POOL_DEFAULT,&
                                                           &CAVG_QUALITY_MODEL_CHUNK_DEFAULT
use simple_class_compatibility,                       only: class_compatibility
use simple_cavg_quality_types,                        only: cavg_quality_result
use simple_cavg_quality_selection,                    only: score_project_cavgs, write_cavg_selection_stacks
use simple_gui_utils,                                 only: mrc2jpeg_tiled
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
public :: PHASE_IMPORTING, PHASE_SOLVE3D, PHASE_IDLE, PHASE_ADDON, PHASE_FINAL, PHASE_SOLVE2D
public :: JOB_NONE, JOB_SOLVE3D, JOB_ADDON
private
#include "simple_local_flags.inc"

! phases of the stage
integer, parameter :: PHASE_IMPORTING = 0 ! no 3D yet
integer, parameter :: PHASE_SOLVE3D   = 1 ! solve3D runs
integer, parameter :: PHASE_IDLE      = 2 ! a result exists; waiting for more particles
integer, parameter :: PHASE_ADDON     = 3 ! solve3D_addon runs
integer, parameter :: PHASE_FINAL     = 4 ! the final refine3D runs
integer, parameter :: PHASE_SOLVE2D   = 5 ! solve2D runs on the first set (no 3D yet)
! the job to start next (next_job)
integer, parameter :: JOB_NONE        = 0
integer, parameter :: JOB_SOLVE3D     = 1
integer, parameter :: JOB_ADDON       = 2

! the first run waits for this many selected particles per state, and an addon run for a cohort of
! as many, solve3D_addon's own floor, and at least ADDON_COHORT_FRAC of the frozen particles
integer, parameter :: MIN_PTCLS_PER_STATE = 5
real,    parameter :: ADDON_COHORT_FRAC   = 0.10
! the first set's solve2D, sized as p03 sizes its own: one class per NPTCLS_PER_CLS2D selected
! particles within NCLS2D_MIN..NCLS2D_MAX, a sample of at least NSAMPLE2D, LPSTOP2D (A)
integer, parameter :: NPTCLS_PER_CLS2D    = 100
integer, parameter :: NCLS2D_MIN          = 10
integer, parameter :: NCLS2D_MAX          = 100
integer, parameter :: NSAMPLE2D           = 2000
real,    parameter :: LPSTOP2D            = 8.

character(len=*), parameter :: SOLVE2D_DIR      = 'solve2D'
character(len=*), parameter :: SOLVE3D_DIR      = 'solve3D'
character(len=*), parameter :: ADDON_DIR        = 'solve3D_addon'
character(len=*), parameter :: FINAL_DIR        = 'refine3D_final'
character(len=*), parameter :: QUALITY_DIR      = 'quality_selection'
character(len=*), parameter :: FIRST_SET_DIR    = 'first_set' ! the first set's selection, in QUALITY_DIR
character(len=*), parameter :: SELECTED_CAVGS   = 'quality_selected_cavgs'
character(len=*), parameter :: REJECTED_CAVGS   = 'quality_rejected_cavgs'
! the quality folders kept: the newest NQUALITY_KEPT, and those of the publications a run started
! from, listed in QUALITY_RUNS (stream fix plan, decision 31)
integer,          parameter :: NQUALITY_KEPT    = 3
character(len=*), parameter :: QUALITY_RUNS     = 'quality_selection/runs.txt'
! the publications passed over as unusable, one file name a line (a restart passes them over too)
character(len=*), parameter :: REJECTED_PUBLICATIONS = 'rejected_publications.txt'

! Components and steps are public so simple_stream_stage_solve3D_tester can assemble a stage
! and run one step at a time; production code uses new/iterate/finished/finalize/kill.
type :: stream_stage_solve3D
    type(parameters), allocatable                   :: params
    ! allocatable (compile-time policy), allocated in init_params and released in kill
    type(sp_project), allocatable                   :: spproj        ! the pool: imported particles, then each run's result
    type(qsys_env),   allocatable                   :: qenv
    type(qsys_async_job)                            :: job           ! the running 3D job
    type(stream_watcher)                            :: project_buff  ! the publications of pool 2D
    type(rec_list)                                  :: setslist      ! one record per publication; included once taken or passed over
    type(stream_pipe)                               :: pipe          ! to the master
    type(gui_metadata_stream_solve3D_multistate) :: meta_status
    type(gui_metadata_cavg2D)                       :: meta_reproj
    type(gui_metadata_stream_snapshot)              :: meta_snapshot
    type(string),   allocatable :: stk_names(:)    ! the stacks in the pool
    real,           allocatable :: state_res(:)    ! FSC=0.143 resolution per state of the latest run (0: none)
    type(string)                :: frozen_projfile ! the latest run's project, which the next addon run builds on
    type(string)                :: result_projfile ! the latest finished run's project (solve3D, addon or final): 3D snapshots' source
    logical,        allocatable :: frozen_active(:) ! the rows active in that project: the frozen particles
    ! the first set (follow-up plan, decisions 31-38): its rows, whose selection the solve2D
    ! selection decides for the session, and the pool model's selection of them, for a failed solve2D
    logical,        allocatable :: first_set(:)
    logical,        allocatable :: first_fallback(:)
    type(string)                :: addon_verdict   ! the latest addon run's verdict per state, for the status
    type(string)                :: last_stem       ! the latest publication taken (its quality folder's name)
    integer :: phase              = PHASE_IMPORTING
    integer :: naddon_runs        = 0
    integer :: nptcls_at_last_run = 0 ! particles in the pool when the latest run started
    integer :: nptcls_selected    = 0 ! selected particles in the rows
    integer :: ncohort_refused    = -1 ! the next addon run needs a larger cohort (after a failed or rolled-back one)
    integer :: n_rollbacks        = 0  ! addon runs rolled back in a row
    integer :: nptcls_at_full     = -1 ! selected particles at the last run that aligned all (solve3D, the final run)
    integer :: last_snapshot_id   = 0  ! the latest 3D snapshot request answered; each is written once
    integer :: optics_id_offset   = 0  ! optics group ids of the snapshots' STAR files, per GUI display
    logical :: l_final_pending    = .false. ! the pool's final publication is in; the final run is due
    logical :: l_solve2D_due      = .false. ! the first set is in; its solve2D is to start
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
    procedure :: select_cavgs
    procedure :: select_with_model
    procedure :: take_first_set
    procedure :: merge_first_set
    procedure :: in_first_set
    procedure :: merge_publication
    procedure :: take_cavgs
    procedure :: stack_index
    procedure :: advance_jobs
    procedure :: start_solve2D
    procedure :: finish_solve2D
    procedure :: fallback_first_set
    procedure :: start_solve3D
    procedure :: start_addon
    procedure :: finish_run
    procedure :: start_final_run
    procedure :: roll_back_addon
    procedure :: read_addon_verdict
    procedure :: count_cohort
    procedure :: count_frozen
    procedure :: record_run_publication
    procedure :: prune_quality_dirs
    procedure :: prune_run_dirs
    procedure :: write_stage_project
    ! GUI
    procedure :: apply_gui_updates
    procedure :: write_snapshot
    procedure :: send_status
    procedure :: send_volumes
    procedure :: read_state_fsc
    ! rules that use no stage state
    procedure, nopass :: next_job
    procedure, nopass :: retry_cohort
    procedure, nopass :: is_final_publication
    procedure, nopass :: fit_mskdiam
end type stream_stage_solve3D

contains

    !---------------- lifecycle ----------------

    !> Builds the parameters from the normalised @p cline, the queue environment of the 3D jobs
    !! and the GUI.
    subroutine new( self, cline )
        class(stream_stage_solve3D), intent(inout) :: self
        class(cmdline),                 intent(inout) :: cline
        call self%kill()
        call self%init_params(cline)
        self%l_exists = .true. ! from here kill() releases what has been built
        call self%init_queue()
        call self%init_gui(ipc_pipe_solve3D_multistate_out(1), ipc_pipe_solve3D_multistate_in(2))
        call self%send_status()
    end subroutine new

    !> The stage's project file (made and given a computing environment) and its parameters; the
    !! project must start without micrographs. A restart's leftover TERM_STREAM is removed in the
    !! stage's folder, where the loop looks for it.
    subroutine init_params( self, cline )
        class(stream_stage_solve3D), intent(inout) :: self
        class(cmdline),                 intent(inout) :: cline
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
        allocate(self%stk_names(0))
        allocate(self%state_res(self%params%nstates), source=0.)
        self%optics_id_offset = max(self%params%nicedispid - 1, 0) * OPTICS_ID_DELTA
    end subroutine init_params

    !> The queue environment the 3D jobs are submitted through; on the preprocessing partition,
    !! as before.
    subroutine init_queue( self )
        class(stream_stage_solve3D), intent(inout) :: self
        call init_stream_qenv(self%params, self%qenv, string(SIMPLE_STREAM_SOLVE3D_PARTITION))
    end subroutine init_queue

    !> The GUI metadata objects and the pipe ends from and to the master (-1: none).
    subroutine init_gui( self, fd_read, fd_write )
        class(stream_stage_solve3D), intent(inout) :: self
        integer,                        intent(in)    :: fd_read, fd_write
        call self%meta_status%new(GUI_METADATA_STREAM_SOLVE3D_MULTISTATE_TYPE)
        call self%meta_reproj%new(GUI_METADATA_STREAM_SOLVE3D_MULTISTATE_REPROJ_TYPE)
        call self%meta_snapshot%new(GUI_METADATA_STREAM_SOLVE3D_SNAPSHOT_TYPE)
        call self%pipe%new(fd_read, fd_write, max_metadata_size(), 'solve3D_multistate')
    end subroutine init_gui

    !> One pass: wait for pool 2D's folder; take new exports (unless a job runs); advance the 3D
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
        ! no import while a job runs on the rows, nor between the first set and its solve2D
        if( self%phase /= PHASE_SOLVE3D .and. self%phase /= PHASE_ADDON .and. self%phase /= PHASE_SOLVE2D&
            &.and. .not. self%l_solve2D_due ) call self%import_sets()
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

    !> The last status, and the stage's project with the latest result and the particles imported
    !! since. A running job is cancelled: a restart starts a new run in a fresh directory.
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
        call self%frozen_projfile%kill
        call self%result_projfile%kill
        call self%addon_verdict%kill
        call self%last_stem%kill
        if( allocated(self%frozen_active) ) deallocate(self%frozen_active)
        if( allocated(self%first_set) ) deallocate(self%first_set)
        if( allocated(self%first_fallback) ) deallocate(self%first_fallback)
        if( allocated(self%stk_names) ) deallocate(self%stk_names)
        if( allocated(self%state_res) ) deallocate(self%state_res)
        if( allocated(self%params) ) deallocate(self%params)
        self%phase              = PHASE_IMPORTING
        self%naddon_runs        = 0
        self%nptcls_at_last_run = 0
        self%nptcls_selected    = 0
        self%ncohort_refused    = -1
        self%n_rollbacks        = 0
        self%nptcls_at_full     = -1
        self%last_snapshot_id   = 0
        self%optics_id_offset   = 0
        self%l_final_pending    = .false.
        self%l_solve2D_due      = .false.
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
    ! every publication: a change applies from the next run (a running one keeps its own) and is
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
                &' A; IT APPLIES FROM THE NEXT RUN'
        else
            write(logfhandle,'(A,F8.2)') '>>> MASK DIAMETER SET TO : ', mskdiam
        endif
        self%mskdiam = mskdiam
    end subroutine take_mskdiam

    ! The newest publication not yet taken: its class averages are selected and it is merged into
    ! the rows; into a stage without rows it is the first set (take_first_set). Its class averages
    ! and FRCs come with its classes (take_cavgs). Older publications
    ! not taken are passed over: the newest holds what they held. A publication the stage cannot
    ! use (publication_problem) is passed over with a warning and listed in REJECTED_PUBLICATIONS,
    ! so a restart passes it over too; the next one is waited for.
    subroutine import_sets( self )
        class(stream_stage_solve3D), intent(inout) :: self
        type(sp_project)   :: set
        type(rec_iterator) :: it
        type(chunk_rec)    :: crec, newest
        type(string)       :: stem, problem
        integer            :: irec, inewest
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
                    call self%take_first_set(set, stem, newest%id)
                else
                    call self%select_cavgs(set, stem)
                    call self%merge_publication(set, newest%id)
                endif
                call self%take_cavgs(set, stem)
                ! the pool's final publication makes the final run due; a later one that is not
                ! final (the sieve took finality back) withdraws it
                self%l_final_pending = is_final_publication(set)
                if( self%l_final_pending ) write(logfhandle,'(A,A)') '>>> THE POOL''S FINAL PUBLICATION: ', stem%to_char()
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

    ! Why publication @p set cannot be used, '' when it can: what select_cavgs, merge_publication and
    ! take_cavgs require of it, checked before any changes anything (class averages, FRCs,
    ! rows_problem).
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

    ! The class averages of publication @p set are scored once with the pool quality model, and
    ! the selection is mapped to its particles by class; the selected and rejected class averages
    ! are written with JPEGs in QUALITY_DIR/<stem>.
    subroutine select_cavgs( self, set, stem )
        class(stream_stage_solve3D), intent(inout) :: self
        class(sp_project),              intent(inout) :: set
        class(string),                  intent(in)    :: stem
        if( set%os_cls2D%get_noris() == 0 ) THROW_HARD('no cls2D entries in the export '//stem%to_char())
        if( .not. self%select_with_model(set, CAVG_QUALITY_MODEL_POOL_DEFAULT, string(QUALITY_DIR//'/')//stem, .false.) )&
            &THROW_HARD('no class averages in the export '//stem%to_char())
    end subroutine select_cavgs

    ! The class averages of @p proj scored with the quality model @p preset, with a mask diameter
    ! fitted to their box; the selection mapped to its selected particles by class; with
    ! @p l_compat, the class-compatibility filter after it, as p03 applies both after its solve2D.
    ! The selected and rejected class averages are written with JPEGs in @p dir. .false. when
    ! @p proj has no class averages.
    logical function select_with_model( self, proj, preset, dir, l_compat ) result( l_ok )
        class(stream_stage_solve3D), intent(inout) :: self
        class(sp_project),           intent(inout) :: proj
        character(len=*),            intent(in)    :: preset
        class(string),               intent(in)    :: dir
        logical,                     intent(in)    :: l_compat
        type(cavg_quality_model)  :: model
        type(cavg_quality_result) :: quality
        type(class_compatibility) :: compatibility
        type(image), allocatable  :: cavg_imgs(:)
        integer,     allocatable  :: states(:)
        type(string)              :: stk, fname
        real                      :: smpd, mskdiam, box_stk
        integer                   :: ncls
        l_ok    = .false.
        box_stk = 0.
        call proj%get_cavgs_stk(stk, ncls, smpd, fail=.false., box=box_stk) ! the os_out entry's box, as a real
        if( ncls <= 0 .or. proj%os_cls2D%get_noris() == 0 ) return
        mskdiam = fit_mskdiam(self%mskdiam, nint(box_stk), smpd)
        call model%init_preset(preset)
        call score_project_cavgs(proj, model, mskdiam, cavg_imgs, quality, smpd=smpd)
        call model%kill
        if( size(cavg_imgs) /= proj%os_cls2D%get_noris() ) THROW_HARD('# class averages /= # cls2D entries in '//dir%to_char())
        call proj%map_cavgs_selection(quality%states)
        call quality%kill
        if( l_compat )then
            call compatibility%new()
            call compatibility%train(proj)
            call compatibility%infer(proj)
            call compatibility%kill()
        endif
        states = proj%os_cls2D%get_all_asint('state')
        call simple_mkdir(QUALITY_DIR)
        call simple_mkdir(dir)
        call write_cavg_selection_stacks(cavg_imgs, states, dir//'/'//SELECTED_CAVGS//MRC_EXT,&
            &dir//'/'//REJECTED_CAVGS//MRC_EXT)
        call mrc2jpeg_tiled(dir//'/'//SELECTED_CAVGS//MRC_EXT, dir//'/'//SELECTED_CAVGS//JPG_EXT)
        call mrc2jpeg_tiled(dir//'/'//REJECTED_CAVGS//MRC_EXT, dir//'/'//REJECTED_CAVGS//JPG_EXT)
        fname = simple_abspath(dir//'/'//SELECTED_CAVGS//JPG_EXT, check_exists=.false.)
        if( file_exists(fname) ) write(logfhandle,'(A,A)') '>>> QUALITY SELECTED CLASS AVERAGES JPEG ', fname%to_char()
        fname = simple_abspath(dir//'/'//REJECTED_CAVGS//JPG_EXT, check_exists=.false.)
        if( file_exists(fname) ) write(logfhandle,'(A,A)') '>>> QUALITY REJECTED CLASS AVERAGES JPEG ', fname%to_char()
        write(logfhandle,'(A,A,A,I6,A,I6)') '>>> CAVG QUALITY (', preset, ') SELECTED / REJECTED : ', count(states > 0), ' / ',&
            &count(states <= 0)
        call dealloc_imgarr(cavg_imgs)
        l_ok = .true.
    end function select_with_model

    ! The first publication, into a stage without rows: every particle it selects (those an
    ! iteration has updated) is the first set, which solve2D classifies again before solve3D
    ! (follow-up plan, decisions 31-38). The pool model's selection of the publication is made
    ! and kept for a failed solve2D (fallback_first_set), but the rows take the publication's own.
    subroutine take_first_set( self, set, stem, id )
        class(stream_stage_solve3D), intent(inout) :: self
        class(sp_project),           intent(inout) :: set
        class(string),               intent(in)    :: stem
        integer,                     intent(in)    :: id
        logical, allocatable :: l_updated(:), l_pool(:)
        integer :: i, n, s
        n = set%os_ptcl2D%get_noris()
        allocate(l_updated(n), l_pool(n))
        do i = 1,n
            l_updated(i) = set%os_ptcl2D%get_state(i) > 0
        enddo
        call self%select_cavgs(set, stem)
        do i = 1,n
            l_pool(i) = set%os_ptcl2D%get_state(i) > 0
            ! back to the publication's own selection
            s = merge(1, 0, l_updated(i))
            call set%os_ptcl2D%set_state(i, s)
            call set%os_ptcl3D%set_state(i, s)
        enddo
        call self%merge_first_set(set, id, l_pool)
        deallocate(l_updated, l_pool)
    end subroutine take_first_set

    ! Merges the first publication @p set (number @p id) into the stage, which has no rows yet, as
    ! it selects its particles. The rows it selects are the first set, due for solve2D;
    ! @p l_pool (per particle of @p set) is the pool model's selection, kept for a failed solve2D.
    subroutine merge_first_set( self, set, id, l_pool )
        class(stream_stage_solve3D), intent(inout) :: self
        class(sp_project),           intent(inout) :: set
        integer,                     intent(in)    :: id
        logical,                     intent(in)    :: l_pool(:)
        integer :: i, n
        if( self%spproj%os_ptcl3D%get_noris() > 0 ) THROW_HARD('the first set needs a stage without rows')
        call self%merge_publication(set, id)
        ! into a stage without rows, the stacks are appended in the publication's order: row i is
        ! its particle i
        n = self%spproj%os_ptcl3D%get_noris()
        if( n /= size(l_pool) ) THROW_HARD('the first set''s rows are not the publication''s particles')
        if( allocated(self%first_set) ) deallocate(self%first_set)
        if( allocated(self%first_fallback) ) deallocate(self%first_fallback)
        allocate(self%first_set(n), self%first_fallback(n))
        do i = 1,n
            self%first_set(i) = self%spproj%os_ptcl2D%get_state(i) > 0
        enddo
        self%first_fallback = l_pool .and. self%first_set
        self%l_solve2D_due  = any(self%first_set)
        write(logfhandle,'(A,I6,A,I8,A)') '>>> PUBLICATION ', id, ' IS THE FIRST SET: SOLVE2D ON ', count(self%first_set),&
            &' PARTICLES BEFORE SOLVE3D'
    end subroutine merge_first_set

    ! .true. for row @p iptcl of the first set, whose selection solve2D decided for the session.
    logical function in_first_set( self, iptcl )
        class(stream_stage_solve3D), intent(in) :: self
        integer,                     intent(in) :: iptcl
        in_first_set = .false.
        if( .not. allocated(self%first_set) ) return
        if( iptcl < 1 .or. iptcl > size(self%first_set) ) return
        in_first_set = self%first_set(iptcl)
    end function in_first_set

    ! Merges the publication @p set (number @p id) into the stage's rows, which only ever grow
    ! (solve3D_addon reads them by index). A stack the stage holds is matched by name and its
    ! particles by image index in the stack: they take the publication's 2D parameters and
    ! selection in place, and keep their 3D parameters, multistate label, CTF and optics group.
    ! A row of the first set takes the 2D parameters only: its selection stays solve2D's. A new
    ! stack is appended with its micrograph and particles (2D and 3D). A stack the publication
    ! lacks keeps its rows, deselected (the first set's keep their selection). The classes are
    ! the publication's; take_cavgs then registers its class averages and FRCs with them.
    subroutine merge_publication( self, set, id )
        class(stream_stage_solve3D), intent(inout) :: self
        class(sp_project),              intent(inout) :: set
        integer,                        intent(in)    :: id
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
                    call self%spproj%os_ptcl2D%set_state(iptcl, s)
                    cycle
                endif
                s     = set%os_ptcl2D%get_state(jptcl)
                call self%spproj%os_ptcl2D%transfer_2Dparams(iptcl, set%os_ptcl2D, jptcl)
                call self%spproj%os_ptcl2D%set_state(iptcl, s)
                ! the 3D state is a run's multistate label: a deselected particle loses it, a
                ! selected one keeps it, and one without a label (state 0) is selected as state 1
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
        class(string),                  intent(in) :: name
        integer,                        intent(in) :: hint
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
    ! copied into its quality folder (QUALITY_DIR/<stem>) and registered in the stage's out segment
    ! in place of the earlier ones; the volumes and FSCs stay. A run's class-average balancing
    ! (balance=cavg) reads the class averages and FRCs of its classes. The copies outlive pool 2D's
    ! retention of its publications, and are pruned with the quality folder, which stays as long
    ! as a run started from the publication (prune_quality_dirs).
    subroutine take_cavgs( self, set, stem )
        class(stream_stage_solve3D), intent(inout) :: self
        class(sp_project),           intent(inout) :: set
        class(string),               intent(in)    :: stem
        type(string) :: stk, frcs, dir, stk_copy, frcs_copy
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
        mskdiam = 0.
        call set%get_mskdiam('cavg', mskdiam)
        if( mskdiam > 0. )then
            call self%spproj%add_cavgs2os_out(stk_copy, smpd, 'cavg', mskdiam=mskdiam)
        else
            call self%spproj%add_cavgs2os_out(stk_copy, smpd, 'cavg')
        endif
        call self%spproj%add_frcs2os_out(frcs_copy, 'frc2D')
    end subroutine take_cavgs

    ! A running job is checked: done, its result becomes the pool (the first set's solve2D: its
    ! selection); a failed solve2D falls back to the pool model's selection, a failed solve3D
    ! stops the stage, a failed addon run leaves the latest result the base. The first set's
    ! selection (or its fallback) goes on to next_job in the same pass, so solve3D starts on it
    ! before a later publication is taken. Without a job, the first set's solve2D starts when due;
    ! otherwise next_job decides whether to start solve3D or an addon run.
    subroutine advance_jobs( self )
        class(stream_stage_solve3D), intent(inout) :: self
        type(string) :: logfile
        if( self%phase == PHASE_SOLVE3D .or. self%phase == PHASE_ADDON .or. self%phase == PHASE_FINAL&
            &.or. self%phase == PHASE_SOLVE2D )then
            select case(self%job%status())
                case(ASYNC_JOB_RUNNING)
                    return
                case(ASYNC_JOB_FAILED)
                    logfile = self%job%get_log()
                    if( self%phase == PHASE_SOLVE2D )then
                        write(logfhandle,'(A,A,A)') '>>> WARNING: SOLVE2D OF THE FIRST SET FAILED (SEE ', logfile%to_char(),&
                            &'); THE POOL MODEL''S SELECTION APPLIES'
                        call self%job%kill()
                        call self%fallback_first_set()
                    else
                        if( self%phase == PHASE_SOLVE3D ) THROW_HARD('solve3D failed; see '//logfile%to_char())
                        if( self%phase == PHASE_FINAL )then
                            ! the latest result stays; the next final publication may try again
                            write(logfhandle,'(A,A,A)') '>>> WARNING: THE FINAL REFINE3D FAILED (SEE ', logfile%to_char(),&
                                &'); THE LATEST RESULT STAYS'
                            call self%job%kill()
                            self%phase = PHASE_IDLE
                            return
                        endif
                        ! a refused or failed addon run: the latest result stays the base, and the
                        ! next run waits for a larger cohort
                        self%ncohort_refused = self%count_cohort()
                        write(logfhandle,'(A,A,A,I8,A)') '>>> WARNING: THE ADDON RUN FAILED (SEE ', logfile%to_char(),&
                            &'); THE NEXT WAITS FOR A COHORT LARGER THAN ', self%ncohort_refused, ' PARTICLES'
                        call self%job%kill()
                        self%phase = PHASE_IDLE
                    endif
                case(ASYNC_JOB_DONE)
                    if( self%phase == PHASE_SOLVE2D )then
                        call self%finish_solve2D()
                    else
                        call self%finish_run()
                    endif
            end select
            ! only the first set's solve2D leaves the stage importing: next_job may start solve3D
            ! on its selection now; with too few particles for it, imports go on in the next pass
            if( self%phase /= PHASE_IMPORTING ) return
        endif
        if( self%l_solve2D_due )then
            call self%start_solve2D()
            return
        endif
        ! the final run, once a result exists, unless the last full alignment had these particles
        if( self%l_final_pending .and. self%phase == PHASE_IDLE )then
            if( self%nptcls_selected /= self%nptcls_at_full )then
                call self%start_final_run()
                return
            endif
            write(logfhandle,'(A)') '>>> THE LAST FULL ALIGNMENT HAD THESE PARTICLES; NO FINAL RUN'
            self%l_final_pending = .false.
        endif
        select case(next_job(self%phase, self%nptcls_selected, self%count_cohort(), self%count_frozen(),&
            &self%params%nstates, self%ncohort_refused))
            case(JOB_SOLVE3D)
                call self%start_solve3D()
            case(JOB_ADDON)
                call self%start_addon()
        end select
    end subroutine advance_jobs

    ! The first set's solve2D, in solve2D/: the selected particles classified again with the
    ! pool's mask diameter, sized as p03 sizes its own (decision 35), on the 3D jobs' resources.
    subroutine start_solve2D( self )
        class(stream_stage_solve3D), intent(inout) :: self
        type(cmdline) :: cline_job
        type(string)  :: cwd, dir, server_address
        integer       :: nptcls, ncls_job, nsample_job
        call simple_getcwd(cwd)
        dir = cwd//'/'//SOLVE2D_DIR
        call fresh_job_dir(dir, SOLVE2D_DIR) ! never the directory of a job left unfinished
        call self%spproj%write(dir//'/'//SOLVE2D_DIR//METADATA_EXT)
        nptcls      = self%spproj%os_ptcl2D%count_state_gt_zero()
        ncls_job    = min(NCLS2D_MAX, max(NCLS2D_MIN, nptcls / NPTCLS_PER_CLS2D))
        nsample_job = ((nptcls / 5 + 999) / 1000) * 1000 ! round up to the nearest 1000
        call cline_job%set('prg',             'solve2D')
        call cline_job%set('mkdir',           'no')
        call cline_job%set('projfile',        SOLVE2D_DIR//METADATA_EXT)
        call cline_job%set('ncls',            ncls_job)
        call cline_job%set('sigma_est',       'global')
        call cline_job%set('center',          'yes')
        call cline_job%set('autoscale',       'yes')
        call cline_job%set('nsample',         max(NSAMPLE2D, nsample_job))
        call cline_job%set('lpstop',          LPSTOP2D)
        if( self%mskdiam > 0. ) call cline_job%set('mskdiam', self%mskdiam)
        call cline_job%set('nparts',          self%params%nparts3D)
        call cline_job%set('nthr',            self%params%nthr3D)
        call cline_job%set('worker_priority', 'high')
        call cline_job%set('cache',           'yes')
        server_address = self%qenv%get_persistent_worker_server_address()
        if( server_address%strlen() > 0 ) call cline_job%set('worker_server', server_address)
        if( self%qenv%get_persistent_worker_nthr() > 0 ) call cline_job%set('worker_server_nthr', self%qenv%get_persistent_worker_nthr())
        call cline_job%printline()
        call self%record_run_publication()
        call self%job%start(self%qenv, cline_job, dir, SOLVE2D_DIR, exec_bin=string('simple_exec'))
        self%l_solve2D_due = .false.
        self%phase         = PHASE_SOLVE2D
        call cline_job%kill
    end subroutine start_solve2D

    ! The first set's solve2D is done: its classes, class averages and 2D parameters become the
    ! rows', its class averages are selected with the chunk model and the compatibility filter
    ! (as p03 selects after its solve2D), and its shifts are the 3D start. Each first-set row's
    ! selection is then its own for the session. Without class averages, the pool model's
    ! selection applies (fallback_first_set).
    subroutine finish_solve2D( self )
        class(stream_stage_solve3D), intent(inout) :: self
        type(string) :: projfile
        integer      :: i, s, n
        projfile = self%job%get_dir()//'/'//SOLVE2D_DIR//METADATA_EXT
        call self%job%kill
        if( .not. file_exists(projfile) ) THROW_HARD('no project from solve2D: '//projfile%to_char())
        n = self%spproj%os_ptcl2D%get_noris()
        call self%spproj%read_segment('ptcl2D', projfile)
        call self%spproj%read_segment('cls2D',  projfile)
        call self%spproj%read_segment('out',    projfile)
        if( self%spproj%os_ptcl2D%get_noris() /= n ) THROW_HARD('solve2D changed the rows: '//projfile%to_char())
        if( .not. self%select_with_model(self%spproj, CAVG_QUALITY_MODEL_CHUNK_DEFAULT,&
            &string(QUALITY_DIR//'/'//FIRST_SET_DIR), .true.) )then
            write(logfhandle,'(A)') '>>> WARNING: SOLVE2D OF THE FIRST SET LEFT NO CLASS AVERAGES; THE POOL MODEL''S SELECTION APPLIES'
            call self%fallback_first_set()
            return
        endif
        ! solve2D's shifts start the 3D, as a publication's ptcl3D carries the pool's
        call self%spproj%map2Dshifts23D()
        ! a first-set row is selected as solve2D and the selection leave it; no other row is
        ! selected before a later publication (the pool had not updated them)
        do i = 1,n
            s = 0
            if( self%in_first_set(i) ) s = merge(1, 0, self%spproj%os_ptcl2D%get_state(i) > 0)
            call self%spproj%os_ptcl2D%set_state(i, s)
            call self%spproj%os_ptcl3D%set_state(i, s)
        enddo
        if( allocated(self%first_fallback) ) deallocate(self%first_fallback)
        self%nptcls_selected = self%spproj%os_ptcl3D%count_state_gt_zero()
        write(logfhandle,'(A,I8,A,I8,A)') '>>> FIRST SET SELECTED: ', self%nptcls_selected, ' OF ', count(self%first_set),&
            &' PARTICLES'
        self%phase = PHASE_IMPORTING
        call self%write_stage_project()
    end subroutine finish_solve2D

    ! The first set without a solve2D selection (it failed, or left no class averages): its rows
    ! take the pool model's selection of the first publication and are no longer the first set,
    ! so the pool governs them like every other row (decision 36).
    subroutine fallback_first_set( self )
        class(stream_stage_solve3D), intent(inout) :: self
        integer :: i, s
        if( allocated(self%first_fallback) )then
            do i = 1,min(size(self%first_fallback), self%spproj%os_ptcl3D%get_noris())
                s = merge(1, 0, self%first_fallback(i))
                call self%spproj%os_ptcl2D%set_state(i, s)
                call self%spproj%os_ptcl3D%set_state(i, s)
            enddo
            deallocate(self%first_fallback)
        endif
        if( allocated(self%first_set) ) deallocate(self%first_set)
        self%l_solve2D_due   = .false.
        self%nptcls_selected = self%spproj%os_ptcl3D%count_state_gt_zero()
        self%phase           = PHASE_IMPORTING
        write(logfhandle,'(A,I8,A)') '>>> THE POOL MODEL SELECTS ', self%nptcls_selected, ' PARTICLES OF THE FIRST PUBLICATION'
    end subroutine fallback_first_set

    ! solve3D on the pool, in solve3D/.
    subroutine start_solve3D( self )
        class(stream_stage_solve3D), intent(inout) :: self
        type(cmdline) :: cline_job
        type(string)  :: cwd, dir, server_address
        call simple_getcwd(cwd)
        dir = cwd//'/'//SOLVE3D_DIR
        call fresh_job_dir(dir, SOLVE3D_DIR) ! never the directory of a job left unfinished
        call self%spproj%write(dir//'/'//SOLVE3D_DIR//METADATA_EXT)
        self%nptcls_at_full = self%nptcls_selected
        call cline_job%set('prg',             'solve3D')
        call cline_job%set('mkdir',           'no')
        call cline_job%set('projfile',        SOLVE3D_DIR//METADATA_EXT)
        call cline_job%set('pgrp',            'c1')
        call cline_job%set('nstates',         self%params%nstates)
        call cline_job%set('nstages',         self%params%nstages)
        call cline_job%set('lpstart',         self%params%lpstart)
        call cline_job%set('lpstop',          self%params%lpstop)
        call cline_job%set('force_lp_range',  'yes')
        if( self%mskdiam > 0. ) call cline_job%set('mskdiam', self%mskdiam)
        call cline_job%set('nparts',          self%params%nparts3D)
        call cline_job%set('nthr',            self%params%nthr3D)
        call cline_job%set('worker_priority', 'high')
        server_address = self%qenv%get_persistent_worker_server_address()
        if( server_address%strlen() > 0 ) call cline_job%set('worker_server', server_address)
        if( self%qenv%get_persistent_worker_nthr() > 0 ) call cline_job%set('worker_server_nthr', self%qenv%get_persistent_worker_nthr())
        call cline_job%printline()
        self%nptcls_at_last_run = self%spproj%os_ptcl3D%get_noris()
        call self%record_run_publication()
        call self%job%start(self%qenv, cline_job, dir, SOLVE3D_DIR, exec_bin=string('simple_exec'))
        self%phase = PHASE_SOLVE3D
        call cline_job%kill
    end subroutine start_solve3D

    ! solve3D_addon on the grown pool, from the latest result, in solve3D_addon/it_<n>/.
    subroutine start_addon( self )
        class(stream_stage_solve3D), intent(inout) :: self
        type(cmdline) :: cline_job
        type(string)  :: cwd, dir, server_address
        self%naddon_runs = self%naddon_runs + 1
        write(logfhandle,'(A,I0)') '>>> ENTERING ADDON STAGE ', self%naddon_runs
        call simple_getcwd(cwd)
        call simple_mkdir(cwd//'/'//ADDON_DIR)
        dir = cwd//'/'//ADDON_DIR//'/it_'//int2str(self%naddon_runs)
        call fresh_job_dir(dir, ADDON_DIR) ! never the directory of a job left unfinished
        call self%spproj%write(dir//'/'//ADDON_DIR//METADATA_EXT)
        call cline_job%set('prg',             'solve3D_addon')
        call cline_job%set('mkdir',           'no')
        call cline_job%set('projfile',        ADDON_DIR//METADATA_EXT)
        call cline_job%set('projfile_frozen', self%frozen_projfile)
        call cline_job%set('nparts',          self%params%nparts3D)
        call cline_job%set('nthr',            self%params%nthr3D)
        call cline_job%set('worker_priority', 'high')
        server_address = self%qenv%get_persistent_worker_server_address()
        if( server_address%strlen() > 0 ) call cline_job%set('worker_server', server_address)
        if( self%qenv%get_persistent_worker_nthr() > 0 ) call cline_job%set('worker_server_nthr', self%qenv%get_persistent_worker_nthr())
        call cline_job%printline()
        self%nptcls_at_last_run = self%spproj%os_ptcl3D%get_noris()
        call self%record_run_publication()
        call self%job%start(self%qenv, cline_job, dir, ADDON_DIR, exec_bin=string('simple_exec'))
        self%phase = PHASE_ADDON
        call cline_job%kill
    end subroutine start_addon

    ! The finished run's result. An addon run's verdict is read first, and a REGRESSED state rolls
    ! the run back (roll_back_addon). Otherwise the job's data become the stage's (the stage keeps
    ! its own project information and computing environment, which the job's project names after
    ! the job, and its optics table) and go to the GUI state by state. A solve3D or addon result is
    ! also the base of the next addon run; the final refine3D's is not, as solve3D_addon needs a
    ! solve3D base: when the sieve takes finality back, the addon runs go on from theirs.
    subroutine finish_run( self )
        class(stream_stage_solve3D), intent(inout) :: self
        type(string)      :: dir, projfile
        real, allocatable :: res_before(:)
        integer           :: i, n
        logical           :: l_addon, l_final, l_regressed
        dir     = self%job%get_dir()
        l_addon = self%phase == PHASE_ADDON
        l_final = self%phase == PHASE_FINAL
        if( l_addon )then
            projfile = dir//'/'//ADDON_DIR//METADATA_EXT
        else if( l_final )then
            projfile = dir//'/'//FINAL_DIR//METADATA_EXT
        else
            projfile = dir//'/'//SOLVE3D_DIR//METADATA_EXT
        endif
        if( .not. file_exists(projfile) ) THROW_HARD('no project from the 3D job: '//projfile%to_char())
        if( l_addon )then
            call self%read_addon_verdict(dir, l_regressed)
            if( l_regressed )then
                call self%roll_back_addon()
                return
            endif
            self%n_rollbacks = 0
        endif
        call self%spproj%read_segment('mic',    projfile)
        call self%spproj%read_segment('stk',    projfile)
        call self%spproj%read_segment('ptcl2D', projfile)
        call self%spproj%read_segment('ptcl3D', projfile)
        call self%spproj%read_segment('cls2D',  projfile)
        call self%spproj%read_segment('cls3D',  projfile)
        call self%spproj%read_segment('out',    projfile)
        self%result_projfile = projfile
        if( .not. l_final )then
            self%frozen_projfile = projfile
            ! the rows active in this result are the next addon run's frozen particles
            n = self%spproj%os_ptcl3D%get_noris()
            if( allocated(self%frozen_active) ) deallocate(self%frozen_active)
            allocate(self%frozen_active(n))
            do i = 1,n
                self%frozen_active(i) = self%spproj%os_ptcl3D%get_state(i) > 0
            enddo
            self%ncohort_refused = -1
            call self%prune_run_dirs()
        endif
        call self%job%kill
        self%phase = PHASE_IDLE
        call self%write_stage_project()
        if( allocated(self%state_res) ) res_before = self%state_res
        call self%send_volumes()
        if( l_final .and. allocated(res_before) )then
            do i = 1,min(size(res_before), size(self%state_res))
                write(logfhandle,'(A,I3,A,F8.2,A,F8.2,A)') '>>> FINAL REFINE3D STATE ', i, ': RESOLUTION ', res_before(i),&
                    &' A BEFORE, ', self%state_res(i), ' A AFTER'
            enddo
        endif
    end subroutine finish_run

    ! The final run (decisions 23 and 24 of the follow-up plan): a multistate refine3D of every
    ! selected particle, started from the current state volumes, in refine3D_final/, once the
    ! pool's final publication is in. refine3D sets its low-pass limits from the FSC, capped by
    ! the stage's lpstop.
    subroutine start_final_run( self )
        class(stream_stage_solve3D), intent(inout) :: self
        type(cmdline) :: cline_job
        type(string)  :: cwd, dir, server_address, volpath
        real          :: smpd
        integer       :: box, istate
        write(logfhandle,'(A,I8,A)') '>>> STARTING THE FINAL REFINE3D OF ', self%nptcls_selected, ' SELECTED PARTICLES'
        call simple_getcwd(cwd)
        dir = cwd//'/'//FINAL_DIR
        call fresh_job_dir(dir, FINAL_DIR) ! never the directory of a job left unfinished
        call self%spproj%write(dir//'/'//FINAL_DIR//METADATA_EXT)
        call cline_job%set('prg',             'refine3D')
        call cline_job%set('mkdir',           'no')
        call cline_job%set('projfile',        FINAL_DIR//METADATA_EXT)
        call cline_job%set('pgrp',            'c1')
        call cline_job%set('nstates',         self%params%nstates)
        do istate = 1,self%params%nstates
            if( .not. self%spproj%isthere_in_osout('vol', istate) ) cycle
            call self%spproj%get_vol('vol', istate, volpath, smpd, box)
            if( volpath%strlen() > 0 ) call cline_job%set('vol'//int2str(istate), simple_abspath(volpath))
        enddo
        call cline_job%set('lpstop',          self%params%lpstop)
        if( self%mskdiam > 0. ) call cline_job%set('mskdiam', self%mskdiam)
        call cline_job%set('nparts',          self%params%nparts3D)
        call cline_job%set('nthr',            self%params%nthr3D)
        call cline_job%set('worker_priority', 'high')
        server_address = self%qenv%get_persistent_worker_server_address()
        if( server_address%strlen() > 0 ) call cline_job%set('worker_server', server_address)
        if( self%qenv%get_persistent_worker_nthr() > 0 ) call cline_job%set('worker_server_nthr', self%qenv%get_persistent_worker_nthr())
        call cline_job%printline()
        self%nptcls_at_last_run = self%spproj%os_ptcl3D%get_noris()
        self%nptcls_at_full     = self%nptcls_selected
        self%l_final_pending    = .false.
        call self%record_run_publication()
        call self%job%start(self%qenv, cline_job, dir, FINAL_DIR, exec_bin=string('simple_exec'))
        self%phase = PHASE_FINAL
        call cline_job%kill
    end subroutine start_final_run

    ! A REGRESSED addon run (decision 22 of the follow-up plan) is not adopted: the frozen base, its
    ! active rows, the stage's project and the GUI's volumes stay the previous run's (the run's
    ! folder stays until a later run is adopted). The next attempt waits until the cohort has grown
    ! by the cadence step again (retry_cohort).
    subroutine roll_back_addon( self )
        class(stream_stage_solve3D), intent(inout) :: self
        self%n_rollbacks     = self%n_rollbacks + 1
        self%ncohort_refused = retry_cohort(self%count_cohort(), self%count_frozen(), self%params%nstates)
        self%addon_verdict   = self%addon_verdict//' - rolled back ('//int2str(self%n_rollbacks)//' in a row)'
        write(logfhandle,'(A,I3,A,I8,A)') '>>> THE ADDON RUN REGRESSED A STATE AND IS ROLLED BACK (', self%n_rollbacks,&
            &' IN A ROW); THE NEXT WAITS FOR A COHORT LARGER THAN ', self%ncohort_refused, ' PARTICLES'
        call self%job%kill
        self%phase = PHASE_IDLE
    end subroutine roll_back_addon

    ! The verdict of the addon run in @p dir (solve3D_addon_report.txt), per state: logged and kept
    ! for the status; @p l_regressed when a state REGRESSED, which rolls the run back (finish_run).
    ! A run that wrote no report is not rolled back.
    subroutine read_addon_verdict( self, dir, l_regressed )
        use simple_solve3D_addon_report, only: solve3D_addon_report, ADDON_REPORT_FNAME
        class(stream_stage_solve3D), intent(inout) :: self
        class(string),               intent(in)    :: dir
        logical,                     intent(out)   :: l_regressed
        type(solve3D_addon_report) :: report
        type(string)               :: fname
        integer                    :: istate
        l_regressed = .false.
        call self%addon_verdict%kill()
        fname = dir//'/'//ADDON_REPORT_FNAME
        if( .not. file_exists(fname) )then
            THROW_WARN('the addon run wrote no report: '//fname%to_char())
            return
        endif
        call report%read(fname)
        self%addon_verdict = 'addon run '//int2str(self%naddon_runs)//':'
        do istate = 1,report%get_nstates()
            self%addon_verdict = self%addon_verdict//' state '//int2str(istate)//' '//trim(report%get_verdict(istate))
        enddo
        write(logfhandle,'(A,A)') '>>> VERDICT OF THE ', self%addon_verdict%to_char()
        l_regressed = report%any_regressed()
        call report%kill()
    end subroutine read_addon_verdict

    ! The addon's cohort: the selected rows (3D state > 0) that were not active in the frozen
    ! solution, appended since or selected again.
    integer function count_cohort( self )
        class(stream_stage_solve3D), intent(in) :: self
        integer :: i, nfrozen_rows
        count_cohort = 0
        nfrozen_rows = 0
        if( allocated(self%frozen_active) ) nfrozen_rows = size(self%frozen_active)
        do i = 1,self%spproj%os_ptcl3D%get_noris()
            if( self%spproj%os_ptcl3D%get_state(i) == 0 ) cycle
            if( i <= nfrozen_rows )then
                if( self%frozen_active(i) ) cycle
            endif
            count_cohort = count_cohort + 1
        enddo
    end function count_cohort

    ! The frozen particles: the rows active in the latest result.
    integer function count_frozen( self )
        class(stream_stage_solve3D), intent(in) :: self
        count_frozen = 0
        if( allocated(self%frozen_active) ) count_frozen = count(self%frozen_active)
    end function count_frozen

    ! The publication a run starts from keeps its quality folder (QUALITY_RUNS).
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

    ! The quality folders kept are the newest NQUALITY_KEPT and those of the publications a run
    ! started from; the others go.
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

    ! Once a run has completed and is the frozen base, the addon iteration folders before it and
    ! the folders set aside for jobs left unfinished go (stream fix plan, decision 31).
    subroutine prune_run_dirs( self )
        class(stream_stage_solve3D), intent(inout) :: self
        type(string), allocatable :: dirs(:)
        integer :: i, k
        do k = 1,self%naddon_runs - 1
            if( dir_exists(ADDON_DIR//'/it_'//int2str(k)) ) call simple_rmdir(ADDON_DIR//'/it_'//int2str(k))
        enddo
        dirs = simple_list_dirs('.')
        do i = 1,size(dirs)
            if( dirs(i)%has_substr(SOLVE3D_DIR//'_unfinished') ) call simple_rmdir(dirs(i)%to_char())
            if( dirs(i)%has_substr(SOLVE2D_DIR//'_unfinished') ) call simple_rmdir(dirs(i)%to_char())
        enddo
        if( dir_exists(ADDON_DIR) )then
            dirs = simple_list_dirs(ADDON_DIR)
            do i = 1,size(dirs)
                if( dirs(i)%has_substr('_unfinished') ) call simple_rmdir(ADDON_DIR//'/'//dirs(i)%to_char())
            enddo
        endif
    end subroutine prune_run_dirs

    ! The pool as the stage's project, written as a temporary file and renamed.
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

    ! A 3D snapshot the GUI asked for, from the latest finished run's result: the particles of the
    ! selected states merged into state 1 (the others deselected, and the state volumes and FSCs
    ! out of the project), written with micrograph and particle STAR files and the optics table to
    ! snapshots/<name>/, beside copies of the selected states' volumes (vol_state<NN>.mrc). Each
    ! request is written once. Before a result, or with no particle in the selected states, it is
    ! reported with no particles and no file.
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

    ! The phase, the run count, the particle counts and, once a run is done, each state's
    ! population and resolution.
    subroutine send_status( self, stage )
        class(stream_stage_solve3D), intent(inout) :: self
        type(string), optional,         intent(in)    :: stage
        type(string) :: stage_here
        integer      :: istate, progress
        select case(self%phase)
            case(PHASE_IMPORTING)
                stage_here = 'waiting for classified particles'
                if( .not. self%l_attached ) stage_here = 'waiting on pool 2D'
            case(PHASE_SOLVE2D)
                stage_here = 'running solve2D on the first set'
            case(PHASE_SOLVE3D)
                stage_here = 'running solve3D'
            case(PHASE_ADDON)
                stage_here = 'running solve3D_addon'
            case(PHASE_FINAL)
                stage_here = 'running final refine3D'
            case default
                stage_here = 'idle'
        end select
        ! the latest addon run's verdict, per state
        if( self%addon_verdict%strlen() > 0 ) stage_here = stage_here//'; '//self%addon_verdict%to_char()
        if( present(stage) ) stage_here = stage
        ! the GUI's 3D progress: 0 none yet (the first set's solve2D among it), 1 solve3D, 2 a result
        progress = min(self%phase, PHASE_IDLE)
        if( self%phase == PHASE_SOLVE2D ) progress = PHASE_IMPORTING
        call self%meta_status%set(stage=stage_here, solve3D_stage=progress,&
            &refine_iteration=self%naddon_runs, nstates=self%params%nstates,&
            &particles_imported=self%spproj%os_ptcl3D%count_state_gt_zero(), particles_at_last_refine=self%nptcls_at_last_run,&
            &resolution=0.)
        if( self%frozen_projfile%strlen() > 0 .and. allocated(self%state_res) )then
            do istate = 1,self%params%nstates
                call self%meta_status%set_state_stats(istate, self%spproj%os_ptcl3D%get_pop(istate, 'state'),&
                    &self%state_res(istate))
            enddo
        endif
        call self%pipe%send_meta(self%meta_status)
    end subroutine send_status

    ! Each state with a volume in the pool's project: its products, population, resolution and FSC
    ! curve, the minimum and maximum of each volume, the orientation distribution, and its three
    ! reprojections as tiles. The resolutions are kept for the status.
    subroutine send_volumes( self )
        class(stream_stage_solve3D), intent(inout) :: self
        type(gui_metadata_vol3D) :: meta, fresh_meta
        type(string)             :: volpath, pprocpath, lppath, pprocmirrpath, reprojpath, oridistpath
        real, allocatable        :: fsc(:), res(:)
        real                     :: smpd, res05, res0143, vmin, vmax
        integer                  :: istate, box, pop, n
        integer                  :: hist(ORIDIST_NBINS_X, ORIDIST_NBINS_Y)
        logical                  :: l_fsc
        self%state_res = 0.
        do istate = 1,self%params%nstates
            if( .not. self%spproj%isthere_in_osout('vol', istate) ) cycle
            call self%spproj%get_vol('vol', istate, volpath, smpd, box)
            if( volpath%strlen() == 0 ) cycle
            pop = self%spproj%os_ptcl3D%get_pop(istate, 'state')
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
            ! JPEGs the 3D writes beside the volume
            reprojpath = get_fpath(volpath)//refine3D_reprojs_fname(istate)
            if( .not. file_exists(reprojpath) ) reprojpath = ''
            oridistpath = get_fpath(volpath)//refine3D_oris_heatmap_fname(istate)
            if( .not. file_exists(oridistpath) ) oridistpath = ''
            if( reprojpath%strlen() > 0 ) call send_reproj_tiles(self%pipe, self%meta_reproj, reprojpath, volpath,&
                &istate, self%params%nstates, pop)
            l_fsc = self%read_state_fsc(istate, smpd, fsc, res, res05, res0143)
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
            call oridist_from_oris(self%spproj%os_ptcl3D, istate, hist)
            call meta%set_oridist(hist)
            call self%pipe%send_meta(meta)
            call meta%kill
        enddo
    end subroutine send_volumes

    ! The FSC of state @p istate in the pool's project, its resolutions at 1/(Fourier index) for
    ! pixel size @p smpd, and the FSC=0.5 and FSC=0.143 resolutions; .false. without an FSC file.
    logical function read_state_fsc( self, istate, smpd, fsc, res, res05, res0143 ) result( l_fsc )
        class(stream_stage_solve3D), intent(in)    :: self
        integer,                        intent(in)    :: istate
        real,                           intent(in)    :: smpd
        real, allocatable,              intent(inout) :: fsc(:), res(:)
        real,                           intent(out)   :: res05, res0143
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

    !> The job to start in @p phase: solve3D once @p nselected particles are selected, at least
    !! MIN_PTCLS_PER_STATE per state; an addon run once the cohort (@p ncohort, selected particles
    !! not frozen) reaches max(MIN_PTCLS_PER_STATE * @p nstates, ADDON_COHORT_FRAC of the
    !! @p nfrozen frozen particles) and exceeds the cohort a failed addon run had
    !! (@p ncohort_refused, -1 for none); none while a job runs.
    pure integer function next_job( phase, nselected, ncohort, nfrozen, nstates, ncohort_refused )
        integer, intent(in) :: phase, nselected, ncohort, nfrozen, nstates, ncohort_refused
        next_job = JOB_NONE
        select case(phase)
            case(PHASE_IMPORTING)
                if( nselected >= MIN_PTCLS_PER_STATE * nstates ) next_job = JOB_SOLVE3D
            case(PHASE_IDLE)
                if( ncohort <= ncohort_refused ) return
                if( ncohort >= max(MIN_PTCLS_PER_STATE * nstates, ceiling(ADDON_COHORT_FRAC * real(nfrozen))) )&
                    &next_job = JOB_ADDON
        end select
    end function next_job

    !> After a rolled-back addon run of cohort @p ncohort: the cohort the next attempt must exceed,
    !! one cadence step (next_job's threshold, with @p nfrozen frozen particles) beyond it less one.
    pure integer function retry_cohort( ncohort, nfrozen, nstates )
        integer, intent(in) :: ncohort, nfrozen, nstates
        retry_cohort = ncohort + max(MIN_PTCLS_PER_STATE * nstates, ceiling(ADDON_COHORT_FRAC * real(nfrozen))) - 1
    end function retry_cohort

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

    !> The mask diameter (A) for class averages of @p box pixels at @p smpd: @p mskdiam, or the
    !! box default when it is not positive or too large for the box, as parameters falls back.
    pure real function fit_mskdiam( mskdiam, box, smpd )
        real,    intent(in) :: mskdiam, smpd
        integer, intent(in) :: box
        real :: msk, msk_default
        fit_mskdiam = mskdiam
        msk_default = round2even((real(box) - COSMSKHALFWIDTH) / 2.)
        if( msk_default <= 0. ) return
        msk = 0.
        if( smpd > 0. ) msk = round2even((mskdiam / smpd) / 2.)
        if( msk < 0.1 .or. msk > msk_default ) fit_mskdiam = (real(box) - COSMSKHALFWIDTH) * smpd
    end function fit_mskdiam

end module simple_stream_stage_solve3D
