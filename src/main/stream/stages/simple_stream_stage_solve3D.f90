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
use simple_defs,                                      only: logfhandle, COSMSKHALFWIDTH
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
use simple_math,                                      only: round2even
use, intrinsic :: ieee_arithmetic,                    only: ieee_is_finite
use simple_math_ft,                                   only: get_resarr
use simple_estimate_ssnr,                             only: get_resolution
use simple_imghead,                                   only: get_mrc_minmax, find_ldim_nptcls
use simple_refine3D_fnames,                           only: refine3D_oris_heatmap_fname, refine3D_reprojs_fname
use simple_oris_utils,                                only: oridist_from_oris
use simple_oris,                                      only: oris
use simple_sym,                                       only: sym
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
use simple_cavg_quality_model,                        only: cavg_quality_model, CAVG_QUALITY_MODEL_POOL_DEFAULT
use simple_cavg_quality_types,                        only: cavg_quality_result
use simple_cavg_quality_selection,                    only: score_project_cavgs, write_cavg_selection_stacks
use simple_dock_vols,                                 only: dock_vols
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
public :: PHASE_IMPORTING, PHASE_SOLVE2D, PHASE_CLSEXP, PHASE_CAVGS3D, PHASE_IDLE, PHASE_EXPAND, PHASE_ASSIGN, PHASE_REFINE
public :: JOB_NONE, JOB_ROUTE, JOB_FALLBACK, JOB_ASSIGN, MIN_PTCLS_PER_STATE, ROUTE_TRIES
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

! a job needs this many selected particles per state
integer, parameter :: MIN_PTCLS_PER_STATE = 5
! the route: solve2D into NCLS2D_FIRST classes, then cls_expansion of the classes the pool model
! keeps into NCLS_EXPANDED sub-classes in all (each publication is expanded into as many); never
! fewer than the classes split, which cls_expansion refuses (expansion_budget)
integer, parameter :: NCLS2D_FIRST        = 150
integer, parameter :: NCLS_EXPANDED       = 500
! cls_expansion's sub-class averages, for the route and every publication alike (the references are
! made from sub-class averages and matched against sub-class averages): kernel pre-images of members
character(len=*), parameter :: SUBCLS_AVG = 'kernel'
! solve3D_cavgs' low-pass ramp (A), fixed (force_lp_range) rather than derived from the class FRCs
real,    parameter :: LPSTART_CAVGS3D     = 20.
real,    parameter :: LPSTOP_CAVGS3D      = 12.
! the check of the reference volumes (provisional, until a validation run sets them)
real,    parameter :: EO_ANGDIST_MAX      = 50.  ! the class averages' even/odd angular distance at most (degrees)
real,    parameter :: EO_OVERLAP_MIN      = 70.  ! their even/odd state agreement at least (%)
real,    parameter :: STATE_FRAC_MIN      = 0.05 ! a state counts when it holds this fraction of the class averages
integer, parameter :: NSTATES_MIN         = 2    ! states that count, at least
real,    parameter :: REFS_PAIR_CC_MAX    = 0.97 ! two states' correlation after docking, at most (distinct states)
real,    parameter :: DOCK_HP_REFS        = 100. ! the docking band: from this high-pass to the matching's limit (A)
integer, parameter :: ROUTE_TRIES         = 2    ! the route's attempts, before the fallback
! the assignment (provisional)
real,    parameter :: LP_ASSIGN           = 15.  ! the matching's low-pass limit (A), the same for every state
integer, parameter :: NSPACE_ASSIGN       = 1000 ! the matching's directions per state, about 6 degrees apart
real,    parameter :: FIT_NMAD            = 3.   ! a best correlation below the median less this many MADs is rejected
real,    parameter :: FIT_REJECT_MAX      = 0.30 ! of a publication's (sub-)classes rejected by fit, at most
real,    parameter :: MARGIN_AMBIGUOUS    = 0.02 ! a margin below this flags a (sub-)class ambiguous
logical, parameter :: L_MODEL_LOGGED      = .true. ! the pool model scores each publication, logged beside the fit

character(len=*), parameter :: SOLVE2D_DIR      = 'solve2D'
character(len=*), parameter :: CLSEXP_DIR       = 'cls_expansion'
character(len=*), parameter :: CAVGS3D_DIR      = 'solve3D_cavgs'
character(len=*), parameter :: ASSIGN_DIR       = 'assign'
character(len=*), parameter :: ASSIGN_LABEL     = 'assign_cavgs_states'
character(len=*), parameter :: EXPAND_DIR       = 'expand'
character(len=*), parameter :: REFINE_DIR       = 'refine3D'
character(len=*), parameter :: RESULT_PROJ      = 'assignment' ! the latest assignment's rows, refined: 3D snapshots' source
character(len=*), parameter :: REFS_DIR         = 'state_refs'
character(len=*), parameter :: REFS_RECORD      = 'state_refs.txt'
character(len=*), parameter :: QUALITY_DIR      = 'quality_selection'
character(len=*), parameter :: SELECTED_CAVGS   = 'quality_selected_cavgs'
character(len=*), parameter :: REJECTED_CAVGS   = 'quality_rejected_cavgs'
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
    procedure :: score_with_model
    procedure :: select_with_model
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
    procedure :: docked_pair_cc
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
    procedure, nopass :: refs_usable
    procedure, nopass :: reject_by_fit
    procedure, nopass :: flag_ambiguous
    procedure, nopass :: best_attempt
    procedure, nopass :: is_final_publication
    procedure, nopass :: expansion_budget
    procedure, nopass :: count_active_classes
    procedure, nopass :: compact_states
    procedure, nopass :: fit_mskdiam
    procedure, nopass :: usable_job_2D
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
                    if( L_MODEL_LOGGED ) l_scored = self%score_with_model(set, CAVG_QUALITY_MODEL_POOL_DEFAULT,&
                        &string(QUALITY_DIR//'/')//stem, model_states)
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

    ! The class averages of @p proj scored with the quality model @p preset, with a mask diameter
    ! fitted to their box: @p cls_states is the model's verdict per class (0 rejected). The
    ! selected and rejected class averages are written with JPEGs in @p dir. Nothing in @p proj
    ! changes. .false. when @p proj has no class averages.
    logical function score_with_model( self, proj, preset, dir, cls_states ) result( l_ok )
        class(stream_stage_solve3D), intent(inout) :: self
        class(sp_project),           intent(inout) :: proj
        character(len=*),            intent(in)    :: preset
        class(string),               intent(in)    :: dir
        integer, allocatable,        intent(inout) :: cls_states(:)
        type(cavg_quality_model)  :: model
        type(cavg_quality_result) :: quality
        type(image), allocatable  :: cavg_imgs(:)
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
        cls_states = quality%states
        call quality%kill
        call simple_mkdir(QUALITY_DIR)
        call simple_mkdir(dir)
        call write_cavg_selection_stacks(cavg_imgs, cls_states, dir//'/'//SELECTED_CAVGS//MRC_EXT,&
            &dir//'/'//REJECTED_CAVGS//MRC_EXT)
        call mrc2jpeg_tiled(dir//'/'//SELECTED_CAVGS//MRC_EXT, dir//'/'//SELECTED_CAVGS//JPG_EXT)
        call mrc2jpeg_tiled(dir//'/'//REJECTED_CAVGS//MRC_EXT, dir//'/'//REJECTED_CAVGS//JPG_EXT)
        fname = simple_abspath(dir//'/'//SELECTED_CAVGS//JPG_EXT, check_exists=.false.)
        if( file_exists(fname) ) write(logfhandle,'(A,A)') '>>> QUALITY SELECTED CLASS AVERAGES JPEG ', fname%to_char()
        fname = simple_abspath(dir//'/'//REJECTED_CAVGS//JPG_EXT, check_exists=.false.)
        if( file_exists(fname) ) write(logfhandle,'(A,A)') '>>> QUALITY REJECTED CLASS AVERAGES JPEG ', fname%to_char()
        write(logfhandle,'(A,A,A,I6,A,I6)') '>>> CAVG QUALITY (', preset, ') SELECTED / REJECTED : ', count(cls_states > 0), ' / ',&
            &count(cls_states <= 0)
        call dealloc_imgarr(cavg_imgs)
        l_ok = .true.
    end function score_with_model

    ! The class averages of @p proj scored with the quality model @p preset (score_with_model) and
    ! the selection mapped to its selected particles by class. .false. when @p proj has no class
    ! averages.
    logical function select_with_model( self, proj, preset, dir ) result( l_ok )
        class(stream_stage_solve3D), intent(inout) :: self
        class(sp_project),           intent(inout) :: proj
        character(len=*),            intent(in)    :: preset
        class(string),               intent(in)    :: dir
        integer, allocatable :: cls_states(:)
        l_ok = self%score_with_model(proj, preset, dir, cls_states)
        if( l_ok ) call proj%map_cavgs_selection(cls_states)
    end function select_with_model

    ! The selected rows whose class (the newest publication's, merged) the pool model rejects,
    ! @p cls_states its verdict per class: logged beside the fit's rejection (take_assignment).
    subroutine flag_model_rows( self, cls_states )
        class(stream_stage_solve3D), intent(inout) :: self
        integer,                     intent(in)    :: cls_states(:)
        integer :: n, i, icls
        n = self%spproj%os_ptcl2D%get_noris()
        if( allocated(self%model_rej) ) deallocate(self%model_rej)
        allocate(self%model_rej(n), source=.false.)
        do i = 1,n
            if( self%spproj%os_ptcl3D%get_state(i) <= 0 ) cycle
            icls = self%spproj%os_ptcl2D%get_class(i)
            if( icls < 1 .or. icls > size(cls_states) ) cycle
            self%model_rej(i) = cls_states(icls) <= 0
        enddo
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
        call cline_job%set('prg',      'solve2D')
        call cline_job%set('mkdir',    'no')
        call cline_job%set('projfile', SOLVE2D_DIR//METADATA_EXT)
        call cline_job%set('ncls',     NCLS2D_FIRST)
        if( self%mskdiam > 0. ) call cline_job%set('mskdiam', self%mskdiam)
        call cline_job%set('cache',    'yes')
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

    ! The first set's classes selected with the pool model (select_with_model), and the selection
    ! mapped to their particles; the stacks of selected and rejected class averages go to
    ! solve2D/quality_selection. The particles of a rejected class are deselected for the session:
    ! they stay first-set rows, whose selection later publications do not change. When the model
    ! cannot score the classes or keeps none, every populated class is selected.
    subroutine select_first_set_classes( self )
        class(stream_stage_solve3D), intent(inout) :: self
        integer, allocatable :: states2D(:), states3D(:), cls_states(:)
        logical :: l_kept
        states2D = self%spproj%os_ptcl2D%get_all_asint('state')
        states3D = self%spproj%os_ptcl3D%get_all_asint('state')
        l_kept   = self%select_with_model(self%spproj, CAVG_QUALITY_MODEL_POOL_DEFAULT, string(SOLVE2D_DIR//'/'//QUALITY_DIR))
        if( l_kept ) l_kept = any(self%spproj%os_cls2D%get_all_asint('state') > 0)
        if( .not. l_kept )then
            write(logfhandle,'(A)') '>>> WARNING: THE POOL MODEL KEPT NO CLASS OF THE FIRST SET; EVERY POPULATED CLASS IS USED'
            ! map_cavgs_selection only deselects: the particles' states from before it first
            call self%spproj%os_ptcl2D%set_all('state', states2D)
            call self%spproj%os_ptcl3D%set_all('state', states3D)
            cls_states = merge(1, 0, self%spproj%os_cls2D%get_all_asint('pop') > 0)
            if( .not. any(cls_states > 0) ) cls_states = 1
            call self%spproj%map_cavgs_selection(cls_states)
        endif
        self%nptcls_selected = self%spproj%os_ptcl3D%count_state_gt_zero()
        write(logfhandle,'(A,I6,A,I6,A,I8,A)') '>>> THE FIRST SET KEEPS ', count(self%spproj%os_cls2D%get_all_asint('state') > 0),&
            &' OF ', self%spproj%os_cls2D%get_noris(), ' CLASSES AND ', self%nptcls_selected, ' PARTICLES'
    end subroutine select_first_set_classes

    ! The route's cls_expansion of the first set's selected classes, in cls_expansion/: the classes
    ! of the selected particles split into expansion_budget sub-classes in all (cls_expansion
    ! spreads them over groups of similar class averages by population and keeps small classes
    ! whole), their averages kernel pre-images of members (SUBCLS_AVG; its other settings are its
    ! defaults), with the mask diameter fitted to the particles' box, on the 3D jobs' resources.
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
        call cline_job%set('prg',      'cls_expansion')
        call cline_job%set('mkdir',    'no')
        call cline_job%set('projfile', CLSEXP_DIR//METADATA_EXT)
        call cline_job%set('ncls',     nsub)
        call cline_job%set('subcls_avg', SUBCLS_AVG)
        call cline_job%set('mskdiam',  fit_mskdiam(self%mskdiam, self%spproj%get_box(), self%spproj%get_smpd()))
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
        call cline_job%set('prg',               'solve3D_cavgs')
        call cline_job%set('mkdir',             'no')
        call cline_job%set('projfile',          CAVGS3D_DIR//METADATA_EXT)
        call cline_job%set('pgrp',              'c1')
        call cline_job%set('nstates',           self%params%nstates)
        call cline_job%set('prune',             'no')
        ! a fixed low-pass ramp, not one derived from the class FRCs
        call cline_job%set('lpstart',           LPSTART_CAVGS3D)
        call cline_job%set('lpstop',            LPSTOP_CAVGS3D)
        call cline_job%set('force_lp_range',    'yes')
        ! the angular sampling window narrowed after the symmetry search, so the poses settle
        call cline_job%set('prob_athres_taper', 'yes')
        if( self%mskdiam > 0. ) call cline_job%set('mskdiam', self%mskdiam)
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

    ! The volumes of solve3D_cavgs' result @p res as attempt @p itry: the states that hold class
    ! averages (cls3D) renumbered from 1 (compact_states), their volumes (vol_cavg entries) and the
    ! JPEGs beside them copied into state_refs/try_<itry>/ with a record (REFS_RECORD: the count,
    ! box and pixel size; the angular distance checked, the state agreement, the largest correlation
    ! between two docked states and the mean angular distance over every class; the file names),
    ! and the attempt's checked even/odd angular distance appended to try_angdist (-1
    ! without volumes; huge when not measured). The states are refined independently, so each is
    ! in a frame and hand of its own: every pair is compared after docking, in both hands
    ! (docked_pair_cc); a pair whose correlation is not finite counts as identical. .true. when
    ! the attempt passes the check (refs_usable).
    logical function register_attempt( self, res, itry ) result( l_pass )
        class(stream_stage_solve3D), intent(inout) :: self
        class(sp_project),           intent(inout) :: res
        integer,                     intent(in)    :: itry
        integer, allocatable :: cls_states(:), newind(:), counts(:)
        type(string), allocatable :: names(:)
        type(string)         :: dir, volpath, src
        real                 :: angdist, angdist_all, overlap, smpd, smpd_try, cc, cc_direct, cc_mirror, cc_max, mskdiam
        integer              :: nheld, s, s2, i, box, box_try, ncopied, funit, ios
        l_pass      = .false.
        angdist     = -1.
        angdist_all = -1.
        overlap     = -1.
        allocate(names(size(self%state_res)))
        if( res%os_cls3D%get_noris() > 0 ) cls_states = res%os_cls3D%get_all_asint('state')
        if( .not. allocated(cls_states) ) allocate(cls_states(0))
        call compact_states(cls_states, newind, nheld)
        if( nheld < 1 .or. nheld > size(names) )then
            write(logfhandle,'(A,I3,A)') '>>> WARNING: SOLVE3D_CAVGS (ATTEMPT ', itry, ') LEFT NO CLASS AVERAGE WITH A STATE'
            self%try_angdist = [self%try_angdist, -1.]
            return
        endif
        allocate(counts(nheld), source=0)
        do i = 1,size(cls_states)
            if( cls_states(i) > 0 ) counts(cls_states(i)) = counts(cls_states(i)) + 1
        enddo
        ! the even/odd agreement solve3D_cavgs records with its volumes: the same values on every
        ! vol_cavg entry (solve3D_cavgs_policy.md, section 7), so the first is read. The angular
        ! distance checked is the median over the class averages whose halves are in one state
        ! (eo_angdist_median); the mean over every class (eo_angdist), which compares poses across
        ! states' frames too, is logged, and checked only when the median is missing
        do i = 1,res%os_out%get_noris()
            if( .not. res%os_out%isthere(i, 'imgkind') ) cycle
            if( res%os_out%get_str(i, 'imgkind') /= 'vol_cavg' ) cycle
            if( res%os_out%isthere(i, 'eo_angdist') )        angdist_all = res%os_out%get(i, 'eo_angdist')
            if( res%os_out%isthere(i, 'eo_angdist_median') ) angdist     = res%os_out%get(i, 'eo_angdist_median')
            if( res%os_out%isthere(i, 'eo_state_overlap') )  overlap     = res%os_out%get(i, 'eo_state_overlap')
            exit
        enddo
        if( angdist < 0. ) angdist = angdist_all
        ! the volumes, renumbered
        dir = string(REFS_DIR//'/try_')//int2str(itry)
        call simple_mkdir(REFS_DIR)
        call simple_mkdir(dir)
        ncopied  = 0
        box_try  = 0
        smpd_try = 0.
        do s = 1,size(newind)
            if( newind(s) < 1 ) cycle
            if( .not. res%isthere_in_osout('vol_cavg', s) ) cycle
            call res%get_vol('vol_cavg', s, volpath, smpd, box)
            if( .not. file_exists(volpath) ) cycle
            names(newind(s)) = 'vol_state'//int2str_pad(newind(s), 2)//MRC_EXT
            call simple_copy_file(volpath, dir//'/'//names(newind(s))%to_char())
            src = get_fpath(volpath)//refine3D_reprojs_fname(s)
            if( file_exists(src) ) call simple_copy_file(src, dir//'/'//refine3D_reprojs_fname(newind(s)))
            src = get_fpath(volpath)//refine3D_oris_heatmap_fname(s)
            if( file_exists(src) ) call simple_copy_file(src, dir//'/'//refine3D_oris_heatmap_fname(newind(s)))
            ncopied  = ncopied + 1
            box_try  = box
            smpd_try = smpd
        enddo
        if( ncopied /= nheld )then
            write(logfhandle,'(A,I3,A,I3,A,I3,A)') '>>> WARNING: SOLVE3D_CAVGS (ATTEMPT ', itry, ') HAS VOLUMES FOR ', ncopied,&
                &' OF THE ', nheld, ' STATES HOLDING CLASS AVERAGES'
            self%try_angdist = [self%try_angdist, -1.]
            return
        endif
        ! every pair of states compared after docking: near-duplicates are no states
        mskdiam = fit_mskdiam(self%mskdiam, box_try, smpd_try)
        cc_max  = -1.
        if( self%l_dock_refs )then
            do s = 1,nheld - 1
                do s2 = s + 1,nheld
                    cc = self%docked_pair_cc(dir//'/'//names(s)%to_char(), dir//'/'//names(s2)%to_char(), smpd_try, mskdiam,&
                        &cc_direct, cc_mirror)
                    write(logfhandle,'(A,I3,A,I3,A,F7.4,A,F7.4,A)') '>>>   STATES ', s, ' AND ', s2, ': CORRELATION AFTER DOCKING ',&
                        &cc_direct, ' (', cc_mirror, ' MIRRORED)'
                    cc_max = max(cc_max, cc)
                enddo
            enddo
        else
            write(logfhandle,'(A)') '>>> THE STATES ARE NOT COMPARED (NO DOCKING)'
        endif
        open(newunit=funit, file=dir%to_char()//'/'//REFS_RECORD, status='replace', action='write', iostat=ios)
        if( ios /= 0 ) THROW_HARD('cannot write the record of '//dir%to_char())
        write(funit,'(I0,1X,I0,1X,F12.6)') nheld, box_try, smpd_try
        write(funit,'(F12.4,1X,F12.4,1X,F12.4,1X,F12.4)') angdist, overlap, cc_max, angdist_all
        do s = 1,nheld
            write(funit,'(A)') names(s)%to_char()
        enddo
        close(funit)
        if( angdist >= 0. )then
            self%try_angdist = [self%try_angdist, angdist]
        else
            self%try_angdist = [self%try_angdist, huge(angdist)]
        endif
        l_pass = refs_usable(angdist, overlap, counts, cc_max)
        write(logfhandle,'(A,I3,A,I3,A,F7.2,A,F7.2,A,F7.2,A,F7.4)') '>>> SOLVE3D_CAVGS, ATTEMPT ', itry, ': ', nheld,&
            &' STATES HOLD CLASS AVERAGES; EVEN/ODD ANGULAR DISTANCE ', angdist, ' DEGREES (MEDIAN, HALVES IN ONE STATE; ',&
            &angdist_all, ' MEAN OVER ALL), STATE AGREEMENT ', overlap, '%, LARGEST CORRELATION BETWEEN TWO DOCKED STATES ', cc_max
        do s = 1,nheld
            write(logfhandle,'(A,I3,A,I6,A)') '>>>   STATE ', s, ': ', counts(s), ' CLASS AVERAGES'
        enddo
        if( l_pass )then
            write(logfhandle,'(A)') '>>> THE VOLUMES PASS THE CHECK'
        else
            write(logfhandle,'(A)') '>>> THE VOLUMES FAIL THE CHECK'
        endif
    end function register_attempt

    ! The correlation of the volumes @p vol1 and @p vol2 (pixel size @p smpd, mask diameter
    ! @p mskdiam in A) once @p vol2 is docked onto @p vol1 over rotations and shifts, from
    ! DOCK_HP_REFS to the matching's limit (LP_ASSIGN): the better of @p vol2 as it is
    ! (@p cc_direct) and mirrored (@p cc_mirror), since independently refined states may differ in
    ! hand. A correlation that is not finite counts as identical (1), with a warning.
    real function docked_pair_cc( self, vol1, vol2, smpd, mskdiam, cc_direct, cc_mirror ) result( cc )
        class(stream_stage_solve3D), intent(inout) :: self
        class(string),               intent(in)    :: vol1, vol2
        real,                        intent(in)    :: smpd, mskdiam
        real,                        intent(out)   :: cc_direct, cc_mirror
        type(dock_vols) :: docker
        type(image)     :: vol
        type(string)    :: mirror
        real    :: eul(3), shift(3)
        integer :: ldim(3), nimgs
        ! vol2 mirrored, beside it
        mirror = add2fbody(vol2, MRC_EXT, '_mirror')
        call find_ldim_nptcls(vol2, ldim, nimgs)
        call vol%new(ldim, smpd)
        call vol%read(vol2)
        call vol%mirror('x')
        call vol%write(mirror)
        call vol%kill
        call docker%new(vol1, vol2, smpd, DOCK_HP_REFS, LP_ASSIGN, mskdiam)
        call docker%srch()
        call docker%get_dock_info(eul, shift, cc_direct)
        call docker%kill()
        call docker%new(vol1, mirror, smpd, DOCK_HP_REFS, LP_ASSIGN, mskdiam)
        call docker%srch()
        call docker%get_dock_info(eul, shift, cc_mirror)
        call docker%kill()
        call del_file(mirror)
        if( .not. ieee_is_finite(cc_direct) .or. .not. ieee_is_finite(cc_mirror) )then
            THROW_WARN('a docked correlation is not finite; the states count as identical')
            if( .not. ieee_is_finite(cc_direct) ) cc_direct = 1.
            if( .not. ieee_is_finite(cc_mirror) ) cc_mirror = 1.
        endif
        cc = max(cc_direct, cc_mirror)
    end function docked_pair_cc

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

    ! The reference volumes recorded in @p dir (REFS_RECORD), taken when every one is on disk:
    ! ref_vols, ref_box and ref_smpd. .false., with none taken, otherwise.
    logical function restore_refs( self, dir ) result( l_ok )
        class(stream_stage_solve3D), intent(inout) :: self
        class(string),               intent(in)    :: dir
        type(string), allocatable :: vols(:)
        character(len=1024) :: line
        real    :: smpd
        integer :: funit, ios, n, box, s
        l_ok = .false.
        if( .not. file_exists(dir//'/'//REFS_RECORD) ) return
        open(newunit=funit, file=dir%to_char()//'/'//REFS_RECORD, status='old', action='read', iostat=ios)
        if( ios /= 0 ) return
        read(funit,*,iostat=ios) n, box, smpd
        if( ios /= 0 .or. n < 1 .or. n > size(self%state_res) )then
            close(funit)
            return
        endif
        read(funit,'(A)',iostat=ios) line ! the even/odd values
        allocate(vols(n))
        do s = 1,n
            read(funit,'(A)',iostat=ios) line
            if( ios /= 0 ) exit
            vols(s) = dir//'/'//trim(adjustl(line))
            if( .not. file_exists(vols(s)) ) ios = 1
            if( ios /= 0 ) exit
            vols(s) = simple_abspath(vols(s))
        enddo
        close(funit)
        if( ios /= 0 ) return
        if( allocated(self%ref_vols) ) deallocate(self%ref_vols)
        allocate(self%ref_vols(n))
        do s = 1,n
            self%ref_vols(s) = vols(s)
        enddo
        self%ref_box  = box
        self%ref_smpd = smpd
        l_ok = .true.
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
        call cline_job%set('prg',      'cls_expansion')
        call cline_job%set('mkdir',    'no')
        call cline_job%set('projfile', CLSEXP_DIR//METADATA_EXT)
        call cline_job%set('ncls',     nsub)
        call cline_job%set('subcls_avg', SUBCLS_AVG)
        call cline_job%set('mskdiam',  fit_mskdiam(self%mskdiam, self%spproj%get_box(), self%spproj%get_smpd()))
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
        integer       :: ncls, s
        call proj%write(self%assign_dir//'/'//ASSIGN_LABEL//METADATA_EXT)
        box_stk = 0.
        call proj%get_cavgs_stk(stk, ncls, smpd, fail=.false., box=box_stk)
        write(logfhandle,'(A,I6,A,I3,A)') '>>> MATCHING ', ncls, ' CLASS AVERAGES AGAINST THE ', self%n_refs(), ' REFERENCES'
        call cline_job%set('prg',      ASSIGN_LABEL)
        call cline_job%set('mkdir',    'no')
        call cline_job%set('projfile', ASSIGN_LABEL//METADATA_EXT)
        call cline_job%set('nstates',  self%n_refs())
        do s = 1,self%n_refs()
            call cline_job%set('vol'//int2str(s), self%ref_vols(s))
        enddo
        call cline_job%set('mskdiam',  fit_mskdiam(self%mskdiam, nint(box_stk), smpd))
        call cline_job%set('lp',       LP_ASSIGN)
        call cline_job%set('nspace',   NSPACE_ASSIGN)
        call cline_job%set('nthr',     self%params%nthr3D)
        call cline_job%printline()
        call self%job%start(self%qenv, cline_job, self%assign_dir, ASSIGN_LABEL, exec_bin=string('simple_private_exec'))
        self%phase = PHASE_ASSIGN
        call cline_job%kill
    end subroutine start_assign

    ! The matching is done: the rejection, the flags and the mapping onto the rows
    ! (take_assignment), recorded in the matched project, then the local refinement. Without a
    ! usable result the states stay as they were.
    subroutine finish_assign( self )
        class(stream_stage_solve3D), intent(inout) :: self
        type(sp_project) :: res
        type(string)     :: projfile
        projfile = self%assign_dir//'/'//ASSIGN_LABEL//METADATA_EXT
        call self%job%kill
        self%phase = PHASE_IDLE
        if( .not. file_exists(projfile) )then
            write(logfhandle,'(A,A,A)') '>>> WARNING: NO PROJECT FROM THE MATCHING (', projfile%to_char(), '); THE STATES STAY'
            return
        endif
        call res%read(projfile)
        if( .not. self%take_assignment(res) )then
            write(logfhandle,'(A,A,A)') '>>> WARNING: NO USABLE MATCHING IN ', projfile%to_char(), '; THE STATES STAY'
            call res%kill
            return
        endif
        call res%write(projfile)
        call res%kill
        call self%write_stage_project()
        call self%start_refine()
    end subroutine finish_assign

    ! The matching in @p res (a copy of the rows, with the matched classes' cls2D, cls3D and the
    ! particles' class labels) onto the rows. Over the classes with a state and particles, those
    ! whose best correlation is an outlier are rejected (reject_by_fit: state 0) and the kept ones
    ! with a small margin flagged (flag_ambiguous: cls2D key 'ambiguous'); each kept class's pose
    ! and state are mapped onto its selected particles (map2ptcls), and each row selected before
    ! takes its particle's state, pose and correlation (state 0 when its class is rejected or it
    ! has none), in both particle segments: a particle follows its (sub-)class. A row deselected
    ! here is recorded (fit_desel): one of the first set is selected again by the next merge. The
    ! fit's statistics, each rejected class's correlation and population, and how the correlation
    ! follows the population are logged, and the counts beside the pool model's verdict
    ! (model_rej). .false., with the rows unchanged, when @p res holds other rows or no matching.
    logical function take_assignment( self, res ) result( l_ok )
        class(stream_stage_solve3D), intent(inout) :: self
        class(sp_project),           intent(inout) :: res
        integer, allocatable :: cls_states(:), counts(:)
        real,    allocatable :: corrs(:), margins(:), pops(:)
        logical, allocatable :: l_cand(:), l_rej(:), l_amb(:), l_sel(:), l_mapped(:)
        real    :: med, mad, thres
        integer :: n, ncls, icls, i, s, nboth, nmodel, nfit, nrej_ptcls, nmax
        l_ok = .false.
        n    = self%spproj%os_ptcl3D%get_noris()
        ncls = res%os_cls2D%get_noris()
        if( n == 0 .or. res%os_ptcl2D%get_noris() /= n ) return
        if( ncls == 0 .or. res%os_cls3D%get_noris() /= ncls ) return
        cls_states = res%os_cls3D%get_all_asint('state')
        corrs      = res%os_cls3D%get_all('corr')
        margins    = res%os_cls3D%get_all('margin')
        ! the candidates: classes with a state and particles
        allocate(l_cand(ncls))
        do icls = 1,ncls
            l_cand(icls) = cls_states(icls) > 0 .and. cls_states(icls) <= self%n_refs()
            if( res%os_cls2D%isthere(icls, 'pop') ) l_cand(icls) = l_cand(icls) .and. res%os_cls2D%get_int(icls, 'pop') > 0
        enddo
        if( .not. any(l_cand) ) return
        call reject_by_fit(corrs, l_cand, l_rej, med, mad, thres, nmax)
        l_amb = flag_ambiguous(margins, l_cand .and. .not. l_rej)
        ! the populations, for the log: whether the rejection follows the classes' sizes
        allocate(pops(ncls), source=0.)
        do icls = 1,ncls
            if( res%os_cls2D%isthere(icls, 'pop') ) pops(icls) = real(res%os_cls2D%get_int(icls, 'pop'))
        enddo
        do icls = 1,ncls
            if( .not. l_cand(icls) .or. l_rej(icls) ) call res%os_cls3D%set_state(icls, 0)
            call res%os_cls2D%set(icls, 'ambiguous', merge(1, 0, l_amb(icls)))
        enddo
        ! the particles each kept class maps onto: selected, with a class of the matching
        allocate(l_sel(n), l_mapped(n))
        do i = 1,n
            l_sel(i)    = self%spproj%os_ptcl3D%get_state(i) > 0
            icls        = res%os_ptcl2D%get_class(i)
            l_mapped(i) = res%os_ptcl2D%get_state(i) > 0 .and. icls >= 1 .and. icls <= ncls
        enddo
        call res%map2ptcls()
        ! map2ptcls maps nothing onto particles without a 2D alignment (a virgin ptcl2D)
        if( res%is_virgin_field('ptcl3D') ) return
        allocate(counts(self%n_refs()), source=0)
        if( allocated(self%fit_desel) ) deallocate(self%fit_desel)
        allocate(self%fit_desel(n), source=.false.)
        nrej_ptcls = 0
        nboth      = 0
        nmodel     = 0
        nfit       = 0
        do i = 1,n
            if( .not. l_sel(i) ) cycle
            s = 0
            if( l_mapped(i) ) s = res%os_ptcl3D%get_state(i)
            if( s > self%n_refs() ) s = 0
            if( s > 0 )then
                call self%spproj%os_ptcl3D%set_euler(i, res%os_ptcl3D%get_euler(i))
                call self%spproj%os_ptcl3D%set_shift(i, res%os_ptcl3D%get_2Dshift(i))
                call self%spproj%os_ptcl3D%set(i, 'corr', res%os_ptcl3D%get(i, 'corr'))
                counts(s) = counts(s) + 1
            else
                nrej_ptcls = nrej_ptcls + 1
                self%fit_desel(i) = .true.
            endif
            call self%spproj%os_ptcl3D%set_state(i, s)
            call self%spproj%os_ptcl2D%set_state(i, s)
            ! the fit's rejection beside the pool model's
            if( allocated(self%model_rej) )then
                if( i <= size(self%model_rej) )then
                    if( s == 0 .and. self%model_rej(i) ) nboth  = nboth  + 1
                    if( s >  0 .and. self%model_rej(i) ) nmodel = nmodel + 1
                    if( s == 0 .and. .not. self%model_rej(i) ) nfit = nfit + 1
                endif
            endif
        enddo
        self%nptcls_selected    = self%spproj%os_ptcl3D%count_state_gt_zero()
        self%nptcls_at_last_run = sum(counts)
        if( self%l_expanded )then
            write(logfhandle,'(A,I6,A)') '>>> ASSIGNMENT OF ', ncls, ' SUBCLASSES'
        else
            write(logfhandle,'(A,I6,A)') '>>> ASSIGNMENT OF ', ncls, ' POOL CLASSES'
        endif
        write(logfhandle,'(A,I6,A,I6,A,I6,A)') '>>> ', count(l_cand), ' CANDIDATES: ', count(l_rej), ' REJECTED BY FIT, ',&
            &count(l_amb), ' AMBIGUOUS (ASSIGNED)'
        call log_fit(corrs, pops, l_cand, l_rej, med, mad, thres, nmax)
        do s = 1,self%n_refs()
            write(logfhandle,'(A,I3,A,I8,A)') '>>>   STATE ', s, ': ', counts(s), ' PARTICLES'
        enddo
        write(logfhandle,'(A,I8,A)') '>>> ', nrej_ptcls, ' PARTICLES DESELECTED'
        if( allocated(self%model_rej) )then
            write(logfhandle,'(A,I8,A,I8,A,I8)') '>>> PARTICLES REJECTED BY THE FIT AND THE POOL MODEL: ', nboth,&
                &'; THE MODEL ALONE (NOT APPLIED): ', nmodel, '; THE FIT ALONE: ', nfit
        endif
        l_ok = .true.

    contains

        ! the fit's statistics; each rejected class's correlation, population and state; the
        ! median population of the kept and rejected candidates; and the correlation of the best
        ! correlation with the log of the population over the candidates (a strong positive one:
        ! the rejection follows the classes' sizes more than their fit)
        subroutine log_fit( corrs, pops, l_cand, l_rej, med, mad, thres, nmax )
            real,    intent(in) :: corrs(:), pops(:), med, mad, thres
            logical, intent(in) :: l_cand(:), l_rej(:)
            integer, intent(in) :: nmax
            real, allocatable :: x(:), y(:)
            real    :: r, sx, sy
            integer :: j
            write(logfhandle,'(A,F7.4,A,F7.4,A,F7.4,A,I6,A,I6,A)') '>>> FIT: MEDIAN BEST CORRELATION ', med, ', MAD ', mad,&
                &', THRESHOLD ', thres, ' (', count(l_cand .and. corrs < thres), ' CANDIDATES BELOW IT, AT MOST ', nmax, ' REJECTED)'
            do j = 1,size(corrs)
                if( .not. l_rej(j) ) cycle
                write(logfhandle,'(A,I6,A,F7.4,A,I8,A,I3)') '>>>   REJECTED CLASS ', j, ': CORRELATION ', corrs(j),&
                    &', POPULATION ', nint(pops(j)), ', STATE ', cls_states(j)
            enddo
            if( .not. any(pops > 0.) ) return
            if( any(l_cand .and. .not. l_rej) .and. any(l_rej) )then
                write(logfhandle,'(A,F10.1,A,F10.1)') '>>> MEDIAN POPULATION OF THE KEPT / REJECTED CANDIDATES: ',&
                    &median_of(pack(pops, l_cand .and. .not. l_rej)), ' / ', median_of(pack(pops, l_rej))
            endif
            if( count(l_cand) < 3 ) return
            x  = pack(corrs, l_cand)
            y  = log(max(pack(pops, l_cand), 1.))
            x  = x - sum(x) / real(size(x))
            y  = y - sum(y) / real(size(y))
            sx = sqrt(sum(x**2))
            sy = sqrt(sum(y**2))
            r  = 0.
            if( sx > 0. .and. sy > 0. ) r = sum(x * y) / (sx * sy)
            write(logfhandle,'(A,F7.4)') '>>> CORRELATION OF THE BEST CORRELATION WITH LOG POPULATION OVER THE CANDIDATES: ', r
        end subroutine log_fit

    end function take_assignment

    ! The seeds of the local refinement on the rows (refine=cont continues from a discrete pose:
    ! state, half-set, projection direction): the half-sets by row parity (partition_eo), the
    ! same for every assignment since rows are never renumbered, and the projection direction
    ! closest to each mapped pose in the job's grid (nspace directions, c1). refine3D clears the
    ! update counts at its first iteration, so its starting reconstruction takes every row with a
    ! state.
    subroutine prepare_refine_seeds( self )
        class(stream_stage_solve3D), intent(inout) :: self
        type(oris) :: eulspace
        type(sym)  :: c1
        call self%spproj%os_ptcl3D%partition_eo()
        call c1%new('c1')
        call eulspace%new(self%params%nspace, is_ptcl=.false.)
        call c1%build_refspiral(eulspace)
        call self%spproj%os_ptcl3D%set_projs(eulspace)
        call eulspace%kill
        call c1%kill
    end subroutine prepare_refine_seeds

    ! The local refinement, in the assignment's refine3D/: one iteration (maxits=1) of refine3D
    ! refine=cont over every selected row from its mapped pose (prepare_refine_seeds), each
    ! particle's pose refined continuously against its own state, which stays, within
    ! refine=cont's bounds (athres_cont, trs); the iteration ends with the states reconstructed
    ! from the refined poses. The starting volumes are reconstructed from the mapped poses; the
    ! low-pass limits come from the FSC; no sigma2 (objfun=cc). Distributed
    ! even with one part: the distributed refine3D reconstructs its starting volumes, the
    ! shared-memory one needs them given. The job holds the states the particles hold,
    ! renumbered from 1 (compact_states) in its project only, since refine3D reconstructs no
    ! volume for an empty state; take_refine_result maps them back.
    subroutine start_refine( self )
        class(stream_stage_solve3D), intent(inout) :: self
        type(cmdline)        :: cline_job
        type(string)         :: dir
        integer, allocatable :: states(:), states_job(:), newind(:)
        integer              :: nheld
        call self%prepare_refine_seeds()
        states     = self%spproj%os_ptcl3D%get_all_asint('state')
        states_job = states
        call compact_states(states_job, newind, nheld)
        if( nheld < 1 )then
            write(logfhandle,'(A)') '>>> WARNING: NO PARTICLE HAS A STATE; NO LOCAL REFINEMENT'
            return
        endif
        dir = self%assign_dir//'/'//REFINE_DIR
        call fresh_job_dir(dir, REFINE_DIR) ! never the directory of a job left unfinished
        call self%spproj%os_ptcl3D%set_all('state', states_job)
        call self%spproj%write(dir//'/'//REFINE_DIR//METADATA_EXT)
        call self%spproj%os_ptcl3D%set_all('state', states)
        write(logfhandle,'(A,I8,A,I3,A)') '>>> REFINING ', self%nptcls_selected, ' PARTICLES LOCALLY IN ', nheld,&
            &' FIXED STATES'
        call cline_job%set('prg',      'refine3D')
        call cline_job%set('mkdir',    'no')
        call cline_job%set('projfile', REFINE_DIR//METADATA_EXT)
        call cline_job%set('refine',   'cont')
        call cline_job%set('maxits',   1)
        call cline_job%set('pgrp',     'c1')
        call cline_job%set('nstates',  nheld)
        call cline_job%set('nspace',   self%params%nspace)
        call cline_job%set('objfun',   'cc')
        if( self%mskdiam > 0. ) call cline_job%set('mskdiam', self%mskdiam)
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
    ! (start_refine renumbered the states the rows hold, compact_states): each state's volume and
    ! FSC it registered, in place of every earlier state's, and each selected row's refined pose
    ! (Euler angles and shift) where its state is the row's. The number of states registered; 0,
    ! with the stage unchanged, when @p res holds other rows or no volume.
    integer function take_refine_result( self, res ) result( nreg )
        class(stream_stage_solve3D), intent(inout) :: self
        class(sp_project),           intent(inout) :: res
        integer, allocatable :: states_job(:), newind(:)
        logical, allocatable :: l_vol(:)
        type(string) :: vol, fsc
        real         :: smpd
        integer      :: n, s, i, box, fsc_box, nposes, nheld
        nreg = 0
        n    = self%spproj%os_ptcl3D%get_noris()
        if( n == 0 .or. res%os_ptcl3D%get_noris() /= n ) return
        ! the job's state of each stage state (0: none)
        states_job = self%spproj%os_ptcl3D%get_all_asint('state')
        call compact_states(states_job, newind, nheld)
        allocate(l_vol(size(newind)), source=.false.)
        do s = 1,size(newind)
            if( newind(s) < 1 ) cycle
            if( .not. res%isthere_in_osout('vol', newind(s)) .or. .not. res%isthere_in_osout('fsc', newind(s)) ) cycle
            call res%get_vol('vol', newind(s), vol, smpd, box)
            call res%get_fsc(newind(s), fsc, fsc_box)
            l_vol(s) = file_exists(vol) .and. file_exists(fsc)
        enddo
        if( .not. any(l_vol) ) return
        do s = 1,size(self%state_res)
            call self%spproj%remove_state_artifacts_from_osout(s)
        enddo
        do s = 1,size(l_vol)
            if( .not. l_vol(s) ) cycle
            call res%get_vol('vol', newind(s), vol, smpd, box)
            call res%get_fsc(newind(s), fsc, fsc_box)
            call self%spproj%add_vol2os_out(vol, smpd, s, 'vol')
            call self%spproj%add_fsc2os_out(fsc, s, fsc_box)
            nreg = nreg + 1
        enddo
        nposes = 0
        do i = 1,n
            s = self%spproj%os_ptcl3D%get_state(i)
            if( s <= 0 .or. s > size(newind) ) cycle
            if( res%os_ptcl3D%get_state(i) /= newind(s) ) cycle
            call self%spproj%os_ptcl3D%set_euler(i, res%os_ptcl3D%get_euler(i))
            call self%spproj%os_ptcl3D%set_shift(i, res%os_ptcl3D%get_2Dshift(i))
            nposes = nposes + 1
        enddo
        write(logfhandle,'(A,I3,A,I8,A)') '>>> LOCAL REFINEMENT: ', nreg, ' STATE VOLUMES, ', nposes, ' REFINED POSES'
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

    !> .true. when reference volumes may be used: the class averages' even/odd angular distance
    !! @p angdist (degrees) at most EO_ANGDIST_MAX (not negative: measured), their even/odd state
    !! agreement @p overlap (%) at least EO_OVERLAP_MIN, the largest correlation between two docked
    !! states @p cc_max at most REFS_PAIR_CC_MAX (distinct states; negative: not compared), and at
    !! least NSTATES_MIN states each holding STATE_FRAC_MIN of the class averages (@p counts per
    !! state).
    pure logical function refs_usable( angdist, overlap, counts, cc_max )
        real,    intent(in) :: angdist, overlap, cc_max
        integer, intent(in) :: counts(:)
        integer :: ntot
        refs_usable = .false.
        if( angdist < 0. .or. angdist > EO_ANGDIST_MAX ) return
        if( overlap < EO_OVERLAP_MIN ) return
        if( cc_max > REFS_PAIR_CC_MAX ) return
        ntot = sum(counts)
        if( ntot < 1 ) return
        refs_usable = count(real(counts) >= STATE_FRAC_MIN * real(ntot)) >= NSTATES_MIN
    end function refs_usable

    !> Of the candidate classes @p l_cand, those whose best correlation @p corrs falls below the
    !! candidates' median less FIT_NMAD times their median absolute deviation are rejected
    !! (@p l_rej), the worst first, at most FIT_REJECT_MAX of the candidates (@p nmax_out). The
    !! median, the deviation and the threshold are returned for the log (0 without candidates).
    pure subroutine reject_by_fit( corrs, l_cand, l_rej, med_out, mad_out, thres_out, nmax_out )
        real,                 intent(in)    :: corrs(:)
        logical,              intent(in)    :: l_cand(:)
        logical, allocatable, intent(inout) :: l_rej(:)
        real,    optional,    intent(out)   :: med_out, mad_out, thres_out
        integer, optional,    intent(out)   :: nmax_out
        real, allocatable :: vals(:)
        real    :: med, mad, thres, worst
        integer :: ncand, nmax, i, iworst
        if( allocated(l_rej) ) deallocate(l_rej)
        allocate(l_rej(size(corrs)), source=.false.)
        if( present(med_out)   ) med_out   = 0.
        if( present(mad_out)   ) mad_out   = 0.
        if( present(thres_out) ) thres_out = 0.
        if( present(nmax_out)  ) nmax_out  = 0
        ncand = count(l_cand)
        if( ncand < 1 ) return
        vals  = pack(corrs, l_cand)
        med   = median_of(vals)
        vals  = abs(vals - med)
        mad   = median_of(vals)
        thres = med - FIT_NMAD * mad
        nmax  = int(FIT_REJECT_MAX * real(ncand))
        if( present(med_out)   ) med_out   = med
        if( present(mad_out)   ) mad_out   = mad
        if( present(thres_out) ) thres_out = thres
        if( present(nmax_out)  ) nmax_out  = nmax
        ! the worst first, while below the threshold and under the cap
        do while( count(l_rej) < nmax )
            iworst = 0
            worst  = thres
            do i = 1,size(corrs)
                if( .not. l_cand(i) .or. l_rej(i) ) cycle
                if( corrs(i) < worst )then
                    worst  = corrs(i)
                    iworst = i
                endif
            enddo
            if( iworst == 0 ) exit
            l_rej(iworst) = .true.
        enddo
    end subroutine reject_by_fit

    !> The median of @p x (none empty).
    pure real function median_of( x )
        real, intent(in) :: x(:)
        real    :: y(size(x)), t
        integer :: n, i, j
        n = size(x)
        y = x
        do i = 2,n
            t = y(i)
            j = i - 1
            do while( j >= 1 )
                if( y(j) <= t ) exit
                y(j+1) = y(j)
                j = j - 1
            enddo
            y(j+1) = t
        enddo
        if( mod(n, 2) == 1 )then
            median_of = y((n + 1) / 2)
        else
            median_of = 0.5 * (y(n / 2) + y(n / 2 + 1))
        endif
    end function median_of

    !> The classes @p l_kept whose @p margins (best less second-best state) fall below
    !! MARGIN_AMBIGUOUS: ambiguous, assigned as any other.
    pure function flag_ambiguous( margins, l_kept ) result( l_amb )
        real,    intent(in) :: margins(:)
        logical, intent(in) :: l_kept(:)
        logical :: l_amb(size(margins))
        l_amb = l_kept .and. margins < MARGIN_AMBIGUOUS
    end function flag_ambiguous

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

    !> The total number of sub-classes cls_expansion splits @p nsplit classes into (its ncls):
    !! NCLS_EXPANDED, and never fewer than the classes, which cls_expansion refuses.
    pure integer function expansion_budget( nsplit ) result( nsub )
        integer, intent(in) :: nsplit
        nsub = max(NCLS_EXPANDED, nsplit)
    end function expansion_budget

    !> The classes cls_expansion splits: the distinct positive class labels @p classes of the rows
    !! with a positive state @p states (its active particles).
    pure integer function count_active_classes( states, classes ) result( n )
        integer, intent(in) :: states(:), classes(:)
        logical, allocatable :: l_seen(:)
        integer :: i, cmax
        n    = 0
        cmax = 0
        do i = 1,min(size(states), size(classes))
            if( states(i) > 0 ) cmax = max(cmax, classes(i))
        enddo
        if( cmax < 1 ) return
        allocate(l_seen(cmax), source=.false.)
        do i = 1,min(size(states), size(classes))
            if( states(i) > 0 .and. classes(i) > 0 ) l_seen(classes(i)) = .true.
        enddo
        n = count(l_seen)
    end function count_active_classes

    !> The positive states of @p states renumbered from 1 in ascending order, in place; @p newind(s)
    !! is state s's new number (0 for a state none holds) and @p nstates its count.
    pure subroutine compact_states( states, newind, nstates )
        integer,              intent(inout) :: states(:)
        integer, allocatable, intent(inout) :: newind(:)
        integer,              intent(out)   :: nstates
        integer :: s, smax, i
        nstates = 0
        smax    = 0
        do i = 1,size(states)
            smax = max(smax, states(i))
        enddo
        if( allocated(newind) ) deallocate(newind)
        allocate(newind(max(smax, 0)), source=0)
        do s = 1,smax
            if( any(states == s) )then
                nstates   = nstates + 1
                newind(s) = nstates
            endif
        enddo
        do i = 1,size(states)
            if( states(i) > 0 ) states(i) = newind(states(i))
        enddo
    end subroutine compact_states

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

    !> .true. when the 2D job's project @p res (solve2D, cls_expansion) can be used: @p n particles
    !! (the stage's rows), one class average per class, and the class averages, their halves and
    !! the FRCs on disk.
    logical function usable_job_2D( res, n ) result( l_ok )
        class(sp_project), intent(inout) :: res
        integer,           intent(in)    :: n
        type(string) :: stk, frcs, ext
        real         :: smpd
        integer      :: ncls
        l_ok = .false.
        ncls = 0
        call res%get_cavgs_stk(stk, ncls, smpd, fail=.false.)
        if( res%os_ptcl2D%get_noris() /= n .or. ncls <= 0 .or. ncls /= res%os_cls2D%get_noris() ) return
        call res%get_frcs(frcs, 'frc2D', fail=.false.)
        ext  = string('.')//fname2ext(stk)
        l_ok = file_exists(stk) .and. file_exists(frcs)
        if( l_ok ) l_ok = file_exists(add2fbody(stk, ext%to_char(), '_even')) .and.&
            &file_exists(add2fbody(stk, ext%to_char(), '_odd'))
    end function usable_job_2D

end module simple_stream_stage_solve3D
