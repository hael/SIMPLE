!@descr: state and steps of stream task 3 (initial analysis): pick, extract, classify and sieve the first micrographs, and make the picking references
!==============================================================================
! MODULE: simple_stream_stage_initial_analysis
!
! PURPOSE:
!   The body of stream p03 as a type. It runs a fixed plan of two cycles over
!   the first preprocessed micrographs and ends by writing the picking
!   references the reference-based picking stage waits for:
!
!   cycle 1 ("init"), on the first NMICS_PLAN(1) accepted micrographs:
!     set up -> wait for the micrographs and pick them (which decides the
!     diameter bins and the box) -> extract -> solve2D -> select classes
!   meanwhile, once the bins are known ("all"): every imported project up to
!   NMICS_PLAN(2) micrographs is picked and extracted with the same bins and
!   box, and its particles are fed to a particle sieve
!   cycle 2 ("all"), once the sieve has every particle:
!     sieve finished -> solve2D -> select classes -> balance classes
!     -> solve3D_cavgs -> reprojections of the best state as references
!
!   The picking references are published once per run, as OPENING2D_PICKREFS
!   in the stage directory, the file the master points reference picking at:
!   written in full under another name and renamed into place, so reference
!   picking never reads a partial stack. Once published they are final.
!
!   A selection of references made in the GUI (for cycle 1 or 2) takes
!   precedence over the 3D route. Each pass reads the GUI's updates before
!   the cycle steps, so a selection pre-empts a 3D result collected in the
!   same pass. A selection that publishes references ends the stage, at once
!   and wherever the plan is; one that selects nothing is logged and ignored.
!   The rest of the plan is skipped, nothing more is imported, and jobs
!   already submitted (the "all" extractions, the sieve's chunks, solve2D/3D)
!   are not stopped: they run on unattended, and no result of theirs is
!   published.
!   The commander (simple_commanders_stream_p03_initial_analysis) only normalises the
!   command line and loops over iterate() until finished().
!
!   What happens is delegated:
!     - picking                 -> segdiam_bin_picker
!     - extraction, 2D, 3D jobs -> qsys_async_job (done or failed, from the exit status)
!     - particle sieve          -> ptcl_sieve
!     - class-average selection -> simple_cavg_quality_selection, class_compatibility
!     - micrograph import       -> import_new_projects, simple_mic_import
!     - micrograph rejection    -> simple_mic_selection
!     - GUI                     -> simple_stream_pipe, simple_stream_gui_senders
!
! LIFECYCLE:
!   new(cline) = init_params -> init_queue -> init_gui(fd_read, fd_write) -> restore_pickrefs;
!   { iterate() } until finished() -> finalize() -> kill()
!   A pass checks for SIGTERM (simple_stream_sigterm) after the import, between
!   the projects of the in-process picking, and before the cycle steps, so it
!   starts no job once a stop is requested.
!
! RESTART:
!   Recognised by the output directory and logged. When the picking references
!   are published already (a user's selection or the 3D route's), they are
!   sent to the GUI again and the stage is finished at once: reference picking
!   may be using them, and a user's selection is never replaced. Otherwise
!   nothing is restored: the plan starts again from cycle 1 and every completed
!   upstream project is imported again. The sieve takes up the chunks left in
!   its folders, and the last numbered solve3D_cavgs directory holds the 3D
!   result.
!
! TESTS:
!   simple_stream_stage_initial_analysis_tester (unit_stream, "initial
!   analysis"). Components and steps are public for it; the waits (settle_s,
!   wait_s) are components it sets, and balance_classes and
!   find_final_solve3D_cavgs_dir, which use no stage state, are bound with
!   nopass so it can call them.
!
! METHOD:
!   Unchanged from the stage this replaces. balance_classes replicates class
!   averages to TARGET_NCLS rows, and the references are reprojections of a
!   3-state ab initio volume; see stream_area_review_2026-09-30.md, B1/B2 and
!   E2/E3, for the proposed changes.
!==============================================================================
module simple_stream_stage_initial_analysis
use simple_core_module_api
use simple_defs_environment,              only: SIMPLE_STREAM_PREPROC_PARTITION
use simple_cmdline,                       only: cmdline
use simple_parameters,                    only: parameters
use simple_sp_project,                    only: sp_project
use simple_image,                         only: image
use simple_image_bin,                     only: image_bin
use simple_imghead,                       only: get_mrc_minmax
use simple_procimgstk,                    only: scale_imgfile
use simple_gui_utils,                     only: mrc2jpeg_tiled
use simple_imgarr_utils,                  only: read_cavgs_into_imgarr, read_stk_into_imgarr, dealloc_imgarr
use simple_qsys_env,                      only: qsys_env
use simple_qsys_async_job,                only: qsys_async_job, ASYNC_JOB_IDLE, ASYNC_JOB_DONE, ASYNC_JOB_FAILED
use simple_rec_list,                      only: rec_list, project_rec
use simple_stream_watcher,                only: stream_watcher
use simple_stream_utils,                  only: create_stream_project, import_new_projects
use simple_stream_state,                  only: ipc_pipe_initial_analysis_in, ipc_pipe_initial_analysis_out
use simple_ptcl_sieve,                    only: ptcl_sieve
use simple_class_compatibility,           only: class_compatibility, support_model_metrics
use simple_cavg_quality_model,            only: cavg_quality_model, CAVG_QUALITY_MODEL_CHUNK_DEFAULT
use simple_cavg_quality_types,            only: cavg_quality_result
use simple_cavg_quality_selection,        only: score_project_cavgs, write_cavg_selection_stacks, write_cavg_stack
use simple_commanders_reproject,          only: commander_reproject
use simple_segdiam_bin_picker,            only: segdiam_bin_picker
use simple_stream_sigterm,                only: sigterm_received
use simple_gui_metadata_utils,            only: max_metadata_size
use simple_gui_metadata_types,            only: GUI_METADATA_STREAM_INITIAL_PICKING_TYPE,&
                                               &GUI_METADATA_STREAM_OPENING2D_TYPE,&
                                               &GUI_METADATA_STREAM_INITIAL_PICKING_MICROGRAPH_TYPE,&
                                               &GUI_METADATA_STREAM_OPENING2D_CLS2D_TYPE,&
                                               &GUI_METADATA_STREAM_OPENING2D_CLS2D_FINAL_TYPE,&
                                               &GUI_METADATA_STREAM_OPENING2D_VOL3D_TYPE
use simple_gui_metadata_cavg2D,           only: gui_metadata_cavg2D
use simple_gui_metadata_micrograph,       only: gui_metadata_micrograph
use simple_gui_metadata_vol3D,            only: gui_metadata_vol3D
use simple_gui_metadata_stream_update,    only: gui_metadata_stream_update
use simple_gui_metadata_stream_picking,   only: gui_metadata_stream_picking
use simple_gui_metadata_stream_opening2D, only: gui_metadata_stream_opening2D
use simple_stream_pipe,                   only: stream_pipe
use simple_stream_gui_senders,            only: send_cavgs, send_recent_micrographs
use simple_mic_import,                    only: append_mics_from_projects
use simple_mic_selection,                 only: reject_mics_without_particles
implicit none

public :: stream_stage_initial_analysis
private
#include "simple_local_flags.inc"

integer, parameter :: NMICS_PLAN(2)       = [100, 500] ! micrographs of the init cycle, and of the "all" set
integer, parameter :: MAX_PROJECTS_IMPORT = 50         ! completed upstream projects taken per pass
integer, parameter :: NTHUMB_MAX          = 10         ! most recent micrograph thumbnails sent to the GUI
integer, parameter :: NCLS_MIN = 10, NCLS_MAX = 100    ! solve2D class-count bounds
integer, parameter :: NSAMPLE2D           = 2000       ! minimum solve2D sample
real,    parameter :: LPSTOP2D            = 8.         ! solve2D low-pass stop (A)
integer, parameter :: NSTATES3D           = 3          ! solve3D_cavgs states
integer, parameter :: TARGET_NCLS         = 501        ! rows after class balancing
character(len=*), parameter :: PICKREFS_SELECTION = 'pickrefs_selection.mrcs' ! a GUI selection, written in full before it is published

! cycle 1 steps
integer, parameter :: INIT_SETUP = 0, INIT_PICK = 1, INIT_EXTRACT = 2, INIT_CLASSIFY = 3, INIT_SELECT = 4
! cycle 2 steps; ALL_COLLECT lasts until the sieve has every particle of the "all" set
integer, parameter :: ALL_COLLECT = 0, ALL_SIEVE = 1, ALL_CLASSIFY = 2, ALL_SELECT = 3, ALL_BALANCE = 4,&
                     &ALL_SOLVE3D = 5, ALL_DONE = 6

! Components and steps are public so simple_stream_stage_initial_analysis_tester can assemble a
! stage and run one step at a time; production code uses new/iterate/finished/finalize/kill.
type :: stream_stage_initial_analysis
    type(parameters)                    :: params
    type(qsys_env)                      :: qenv
    type(sp_project)                    :: spproj                  ! cycle 1 project
    type(sp_project)                    :: spproj_part             ! one project of the "all" set
    type(sp_project)                    :: spproj_all              ! cycle 2 project
    type(stream_watcher)                :: project_buff
    type(rec_list)                      :: project_list            ! one record per accepted imported micrograph
    type(rec_list)                      :: extracted_project_list  ! one record per extracted micrograph of the "all" set
    type(segdiam_bin_picker)            :: picker
    type(ptcl_sieve)                    :: sieve
    type(qsys_async_job)                :: job                     ! the current cycle step's extraction, 2D or 3D run
    type(qsys_async_job),   allocatable :: extract_jobs(:)         ! the "all" extractions, one per project
    logical,                allocatable :: extract_collected(:)
    type(stream_pipe)                   :: pipe
    type(gui_metadata_stream_picking)   :: meta_picking
    type(gui_metadata_stream_opening2D) :: meta_opening2D
    type(gui_metadata_micrograph)       :: meta_micrograph
    type(gui_metadata_cavg2D)           :: meta_cavg2D, meta_pickrefs
    type(string)                        :: cwd                     ! absolute stage directory
    integer :: icycle            = 1
    integer :: step1             = INIT_SETUP
    integer :: step2             = ALL_COLLECT
    integer :: nmics_target      = NMICS_PLAN(1)
    integer :: n_mics_imported   = 0
    integer :: n_ptcls_imported  = 0
    integer :: n_extract_started = 0
    integer :: n_extract_done    = 0
    integer :: box               = 0       ! picking and extraction box (px)
    real    :: mskdiam           = 0.      ! mask diameter (A) decided with the box
    integer :: vis_cycle         = 0
    logical :: l_attached        = .false.
    logical :: l_waiting_logged  = .false.
    logical :: l_sieve_active    = .false.
    logical :: l_done            = .false.
    logical :: l_exists          = .false.
    ! waits (s); tests set them to 0, and settle_s to -1 to take files written in the same second
    integer :: settle_s          = SHORTWAIT ! an upstream project is taken once untouched longer than this;
                                             ! preprocessing moves finished projects in with a rename
    integer :: wait_s            = WAITTIME ! pause at the end of a pass
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
    procedure :: restore_pickrefs
    ! the steps of iterate() and their helpers
    procedure :: attach_upstream
    procedure :: import_projects
    procedure :: pick_extract_all
    procedure :: collect_extractions
    procedure :: start_sieve
    procedure :: run_cycle1
    procedure :: run_cycle2
    procedure :: rebuild_init_mics
    procedure :: select_and_send
    procedure :: finish_solve3D
    procedure :: apply_gui_updates
    procedure :: save_pickrefs_selection
    procedure :: publish_pickrefs
    procedure :: send_pickrefs
    procedure :: send_picking_status
    procedure :: send_opening2D_status
    procedure :: cycle_projfile
    procedure :: all_projfile
    ! steps that use no stage state, bound so the tester can reach them
    procedure, nopass :: balance_classes
    procedure, nopass :: find_final_solve3D_cavgs_dir
end type stream_stage_initial_analysis

contains

    !---------------- lifecycle ----------------

    !> Builds the parameters from the normalised @p cline and prepares job submission and the GUI.
    subroutine new( self, cline )
        class(stream_stage_initial_analysis), intent(inout) :: self
        class(cmdline),                       intent(inout) :: cline
        call self%kill()
        call self%init_params(cline)
        call self%init_queue()
        call self%init_gui(ipc_pipe_initial_analysis_out(1), ipc_pipe_initial_analysis_in(2))
        allocate(self%extract_jobs(0), self%extract_collected(0))
        call self%restore_pickrefs()
        self%l_exists = .true.
    end subroutine new

    !> The stage's project file (made and given a computing environment) and its parameters.
    subroutine init_params( self, cline )
        class(stream_stage_initial_analysis), intent(inout) :: self
        class(cmdline),                       intent(inout) :: cline
        type(string) :: outdir
        ! a restart is recognised by its output directory, before params%new creates one
        ! (this stage is stopped by SIGTERM only and does not use TERM_STREAM)
        outdir = cline%get_carg('outdir')
        if( .not. (outdir == '') )then
            if( dir_exists(outdir) ) write(logfhandle,'(A)') '>>> RESTARTING EXISTING JOB'
        endif
        call create_stream_project(self%spproj, cline, string('opening_2D'))
        call self%params%new(cline)
        call simple_getcwd(self%cwd)
    end subroutine init_params

    !> The queue environment the extraction, 2D and 3D jobs are submitted through.
    subroutine init_queue( self )
        class(stream_stage_initial_analysis), intent(inout) :: self
        character(len=STDLEN) :: partition_env
        integer               :: envlen
        call get_environment_variable(SIMPLE_STREAM_PREPROC_PARTITION, partition_env, envlen)
        if( envlen > 0 )then
            call self%qenv%new(self%params, 1, stream=.true., qsys_partition=string(trim(partition_env)))
        else
            call self%qenv%new(self%params, 1, stream=.true.)
        endif
    end subroutine init_queue

    !> The GUI metadata objects and the pipe ends to the master (-1: none).
    subroutine init_gui( self, fd_read, fd_write )
        class(stream_stage_initial_analysis), intent(inout) :: self
        integer,                              intent(in)    :: fd_read, fd_write
        call self%meta_picking%new(GUI_METADATA_STREAM_INITIAL_PICKING_TYPE)
        call self%meta_opening2D%new(GUI_METADATA_STREAM_OPENING2D_TYPE)
        call self%meta_micrograph%new(GUI_METADATA_STREAM_INITIAL_PICKING_MICROGRAPH_TYPE)
        call self%meta_cavg2D%new(GUI_METADATA_STREAM_OPENING2D_CLS2D_TYPE)
        call self%meta_pickrefs%new(GUI_METADATA_STREAM_OPENING2D_CLS2D_FINAL_TYPE)
        call self%pipe%new(fd_read, fd_write, max_metadata_size(), 'initial_analysis')
    end subroutine init_gui

    !> On a restart: picking references an earlier run published are final. They go to the GUI
    !! again and the stage is finished, so the plan does not run again and cannot replace them.
    subroutine restore_pickrefs( self )
        class(stream_stage_initial_analysis), intent(inout) :: self
        if( .not. file_exists(self%cwd//'/'//OPENING2D_PICKREFS) ) return
        write(logfhandle,'(A)') '>>> PICKING REFERENCES ALREADY PUBLISHED: '//OPENING2D_PICKREFS//'; NOTHING TO DO'
        call self%send_pickrefs()
        self%l_done = .true.
    end subroutine restore_pickrefs

    !> One pass: import, pick and extract the "all" set, advance the cycle, take GUI updates.
    subroutine iterate( self )
        class(stream_stage_initial_analysis), intent(inout) :: self
        if( .not. self%l_attached )then
            call self%attach_upstream()
            if( .not. self%l_attached )then
                call sleep(self%wait_s)
                return
            endif
        endif
        call self%import_projects()
        if( sigterm_received() ) return
        if( self%picker%bins_set() .and. self%project_list%size() > 0 ) call self%pick_extract_all()
        ! a selection made in the GUI comes first: it pre-empts the 3D route, even one finishing now
        call self%apply_gui_updates()
        ! no new job once a stop is requested
        if( self%l_done .or. sigterm_received() ) return
        if( self%icycle == 1 )then
            call self%run_cycle1()
        else
            call self%run_cycle2()
        endif
        if( self%l_done .or. sigterm_received() ) return
        call sleep(self%wait_s)
    end subroutine iterate

    !> .true. once the picking references are published, by the 3D route, from a GUI selection
    !! or by an earlier run.
    logical function finished( self )
        class(stream_stage_initial_analysis), intent(in) :: self
        finished = self%l_done
    end function finished

    subroutine finalize( self )
        class(stream_stage_initial_analysis), intent(inout) :: self
        call self%send_opening2D_status(string('terminating'), self%box, self%vis_cycle)
    end subroutine finalize

    subroutine kill( self )
        class(stream_stage_initial_analysis), intent(inout) :: self
        integer :: i
        if( .not. self%l_exists ) return
        call self%spproj%kill
        call self%spproj_part%kill
        call self%spproj_all%kill
        call self%project_buff%kill
        call self%picker%kill
        if( self%l_sieve_active ) call self%sieve%kill
        call self%job%kill
        if( allocated(self%extract_jobs) )then
            do i = 1,size(self%extract_jobs)
                call self%extract_jobs(i)%kill
            enddo
            deallocate(self%extract_jobs)
        endif
        if( allocated(self%extract_collected) ) deallocate(self%extract_collected)
        call self%qenv%kill
        call self%pipe%kill
        call self%meta_picking%kill
        call self%meta_opening2D%kill
        call self%meta_micrograph%kill
        call self%meta_cavg2D%kill
        call self%meta_pickrefs%kill
        call self%cwd%kill
        self%icycle            = 1
        self%step1             = INIT_SETUP
        self%step2             = ALL_COLLECT
        self%nmics_target      = NMICS_PLAN(1)
        self%n_mics_imported   = 0
        self%n_ptcls_imported  = 0
        self%n_extract_started = 0
        self%n_extract_done    = 0
        self%box               = 0
        self%mskdiam           = 0.
        self%vis_cycle         = 0
        self%l_attached        = .false.
        self%l_waiting_logged  = .false.
        self%l_sieve_active    = .false.
        self%l_done            = .false.
        self%l_exists          = .false.
    end subroutine kill

    !---------------- import, and the "all" set ----------------

    ! Starts watching the upstream completed-projects folder once preprocessing has created it.
    subroutine attach_upstream( self )
        class(stream_stage_initial_analysis), intent(inout) :: self
        type(string) :: completed
        logical      :: l_ready
        completed = self%params%dir_target//'/'//DIR_STREAM_COMPLETED
        l_ready   = dir_exists(self%params%dir_target)
        if( l_ready ) l_ready = dir_exists(self%params%dir_target//'/'//DIR_STREAM)
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
        self%l_attached   = .true.
    end subroutine attach_upstream

    ! One record per accepted micrograph of the newly completed upstream projects.
    subroutine import_projects( self )
        class(stream_stage_initial_analysis), intent(inout) :: self
        type(string), allocatable :: projects(:)
        integer :: nprojects
        call self%project_buff%watch(nprojects, projects, max_nmovies=MAX_PROJECTS_IMPORT)
        if( nprojects == 0 ) return
        call import_new_projects(self%project_list, projects, self%n_mics_imported, self%n_ptcls_imported,&
            &ignore_ptcls=.true., check_state=.true.)
        call self%project_buff%add2history(projects)
        write(logfhandle,'(A,I6)') '>>> # MICROGRAPHS IMPORTED : ', self%n_mics_imported
        write(logfhandle,'(A,A)')  '>>> LAST IMPORT AT         : ', cast_time_char(simple_gettime())
    end subroutine import_projects

    ! Picks and extracts every project not yet picked, up to NMICS_PLAN(2) micrographs, with the
    ! bins and box of cycle 1; collects finished extractions into the sieve.
    subroutine pick_extract_all( self )
        class(stream_stage_initial_analysis), intent(inout) :: self
        type(qsys_async_job) :: new_job
        type(project_rec)    :: prec
        type(string)         :: cur_projname, proj_local
        integer              :: iproj, i, nmics
        call simple_mkdir(DIR_STREAM)
        call simple_mkdir(DIR_STREAM//'all')
        write(logfhandle,'(A,I0)') '>>> PICKING AND EXTRACTING ALL UNPROCESSED PROJECTS IN project_list; # INCLUDED : ',&
            &count(self%project_list%get_included_flags())
        if( count(self%project_list%get_included_flags()) < NMICS_PLAN(2) )then
            do iproj = 1, self%project_list%size()
                ! picking runs in-process; stop between projects, each of which is complete
                if( sigterm_received() ) return
                call self%project_list%at(iproj, prec)
                if( prec%included ) cycle
                self%n_extract_started = self%n_extract_started + 1
                cur_projname = prec%projname
                proj_local   = self%all_projfile(self%n_extract_started)
                write(logfhandle,'(A,A,A)') '>>> PICKING PROJECT : ', cur_projname%to_char(), proj_local%to_char()
                call simple_copy_file(cur_projname, proj_local)
                call self%spproj_part%read(proj_local)
                nmics = self%spproj_part%os_mic%get_noris()
                if( nmics == 0 )then
                    ! nothing to extract: counted as collected so the sieve's final ingestion is not held up
                    call self%spproj_part%kill()
                    self%extract_jobs      = [self%extract_jobs, new_job]
                    self%extract_collected = [self%extract_collected, .true.]
                    self%n_extract_done    = self%n_extract_done + 1
                    cycle
                endif
                call self%send_picking_status(string('picking particles'))
                call simple_mkdir(DIR_PICKER//'all')
                call simple_chdir(DIR_PICKER//'all')
                call self%picker%pick(self%spproj_part, self%params%pcontrast, nmics)
                call simple_chdir(self%cwd)
                call self%send_picking_status(string('extracting particles'))
                call start_extract(self%qenv, new_job, proj_local, string(DIR_EXTRACT//'all'), self%box,&
                    &self%n_extract_started, self%spproj_part%os_mic%get_noris())
                self%extract_jobs      = [self%extract_jobs, new_job]
                self%extract_collected = [self%extract_collected, .false.]
                call new_job%kill()
                call self%send_picking_status(string('complete'))
                call self%spproj_part%kill()
                ! flag every record of this project as picked, so it is never picked again
                do i = 1, self%project_list%size()
                    call self%project_list%at(i, prec)
                    if( prec%included ) cycle
                    if( prec%projname%to_char() /= cur_projname%to_char() ) cycle
                    prec%included = .true.
                    prec%projname = proj_local
                    call self%project_list%replace_at(i, prec)
                end do
            end do
        end if
        call self%collect_extractions()
        if( self%l_sieve_active )then
            call self%sieve%cycle(self%extracted_project_list)
            if( count(self%project_list%get_included_flags()) >= NMICS_PLAN(2) .and.&
                &self%n_extract_done == self%n_extract_started .and. self%step2 == ALL_COLLECT )then
                call self%sieve%set_final_ingestion()
                self%step2 = ALL_SIEVE
                write(logfhandle,'(A)') '>>> ALL PROJECTS PICKED AND EXTRACTED, SIEVE FINAL INGESTION SET'
            end if
        end if
    end subroutine pick_extract_all

    ! Imports each finished "all" extraction once: one sieve record per extracted micrograph.
    ! A failed extraction is logged and left out.
    subroutine collect_extractions( self )
        class(stream_stage_initial_analysis), intent(inout) :: self
        type(project_rec) :: extracted_prec
        type(string)      :: proj_local, job_log
        integer           :: iextract, imic, job_status
        do iextract = 1, self%n_extract_started
            if( self%extract_collected(iextract) ) cycle
            job_status = self%extract_jobs(iextract)%status()
            if( job_status == ASYNC_JOB_FAILED )then
                job_log = self%extract_jobs(iextract)%get_log()
                write(logfhandle,'(A,I0,A,A)') '>>> EXTRACTION ', iextract, ' FAILED; ITS MICROGRAPHS ARE LEFT OUT. LOG: ',&
                    &job_log%to_char()
                self%extract_collected(iextract) = .true.
                self%n_extract_done = self%n_extract_done + 1
                cycle
            endif
            if( job_status /= ASYNC_JOB_DONE ) cycle
            proj_local = self%all_projfile(iextract)
            call finish_extract(self%spproj_part, proj_local, string(DIR_EXTRACT//'all'))
            do imic = 1, self%spproj_part%os_mic%get_noris()
                if( self%spproj_part%os_mic%get_int(imic, 'state') <= 0 ) cycle
                if( .not. self%spproj_part%os_mic%isthere(imic, 'nptcls') ) cycle
                if( self%spproj_part%os_mic%get_int(imic, 'nptcls') <= 0 ) cycle
                extracted_prec%id         = self%extracted_project_list%size() + 1
                extracted_prec%projname   = proj_local
                extracted_prec%micind     = imic
                extracted_prec%nptcls     = self%spproj_part%os_mic%get_int(imic, 'nptcls')
                extracted_prec%nptcls_sel = extracted_prec%nptcls
                extracted_prec%included   = .false.
                call self%extracted_project_list%push_back(extracted_prec)
            end do
            self%extract_collected(iextract) = .true.
            self%n_extract_done = self%n_extract_done + 1
            ! the sieve starts once, while the "all" set is still being collected; never again after it is killed
            if( .not. self%l_sieve_active .and. self%step2 == ALL_COLLECT ) call self%start_sieve()
            call self%send_picking_status(string('complete'))
            call send_recent_micrographs(self%pipe, self%meta_micrograph, self%spproj_part%os_mic, NTHUMB_MAX)
        end do
    end subroutine collect_extractions

    ! The particle sieve of the "all" set, started with the first finished extraction.
    subroutine start_sieve( self )
        class(stream_stage_initial_analysis), intent(inout) :: self
        type(parameters) :: params_sieve
        ! the sieve takes a parameters object; this stage's settings for it
        params_sieve             = self%params
        params_sieve%mskdiam     = 0.
        params_sieve%lpstart     = 0.
        params_sieve%single_pass = 'yes'
        params_sieve%nmics       = 100
        params_sieve%nchunks     = 4
        params_sieve%workers     = params_sieve%nchunks
        params_sieve%nthr        = 16
        params_sieve%worker_nthr = params_sieve%nthr
        call simple_mkdir(self%cwd//'/spprojs_sieved')
        call self%sieve%new(params_sieve, self%cwd//'/spprojs_sieved')
        call self%sieve%cycle(self%extracted_project_list)
        call self%sieve%cycle(self%extracted_project_list)
        self%l_sieve_active = .true.
    end subroutine start_sieve

    !---------------- the two cycles ----------------

    ! Cycle 1 on the first NMICS_PLAN(1) micrographs. Several steps can follow one another in one pass.
    subroutine run_cycle1( self )
        class(stream_stage_initial_analysis), intent(inout) :: self
        type(string) :: projfile, job_log
        projfile = self%cycle_projfile(1)
        if( self%step1 == INIT_SETUP )then
            self%nmics_target = NMICS_PLAN(1)
            call simple_chdir(self%cwd)
            call simple_mkdir(self%cwd//'/'//DIR_STREAM)
            call simple_mkdir(self%cwd//'/'//DIR_STREAM//'init')
            call self%spproj%read(self%params%projfile)
            write(logfhandle, '(A,I6,A,I6,A)') '>>> INITIATED INITIAL ANALYSIS CYCLE ', 1, ' WITH AT LEAST ',&
                &self%nmics_target, ' MICROGRAPHS'
            call self%spproj%update_projinfo(projfile)
            call self%spproj%write()
            call self%send_picking_status(string('waiting for micrographs'))
            call self%send_opening2D_status(string('waiting for particles'), self%box, self%vis_cycle)
            self%step1 = INIT_PICK
        endif
        if( self%step1 == INIT_PICK )then
            if( self%n_mics_imported > self%nmics_target )then
                write(logfhandle,'(A,I6,A,I6)') '>>> IMPORTED SUFFICIENT MICROGRAPHS: ', self%n_mics_imported, ' >= ',&
                    &self%nmics_target
                call self%rebuild_init_mics()
                call self%send_picking_status(string('picking particles'))
                call simple_mkdir(DIR_PICKER//'init')
                call simple_chdir(DIR_PICKER//'init')
                ! the first pick decides the diameter bins and the box for every later one
                call self%picker%new()
                call self%picker%pick(self%spproj, self%params%pcontrast, self%nmics_target)
                if( self%picker%get_box() > 0 )then
                    self%box     = self%picker%get_box()
                    self%mskdiam = self%picker%get_mskdiam()
                endif
                call send_recent_micrographs(self%pipe, self%meta_micrograph, self%spproj%os_mic, NTHUMB_MAX)
                call simple_chdir(self%cwd)
                self%step1 = INIT_EXTRACT
            endif
        endif
        if( self%step1 == INIT_EXTRACT )then
            select case( self%job%status() )
                case( ASYNC_JOB_IDLE )
                    call self%send_picking_status(string('extracting particles'))
                    call start_extract(self%qenv, self%job, projfile, string(DIR_EXTRACT//'init'), self%box, 1,&
                        &self%spproj%os_mic%get_noris())
                case( ASYNC_JOB_DONE )
                    write(logfhandle,'(A)') '>>> PARTICLE EXTRACTION COMPLETE'
                    call finish_extract(self%spproj, projfile, string(DIR_EXTRACT//'init'))
                    call self%send_picking_status(string('complete'))
                    call self%job%kill()
                    self%step1 = INIT_CLASSIFY
                case( ASYNC_JOB_FAILED )
                    job_log = self%job%get_log()
                    THROW_HARD('initial extraction failed; see '//job_log%to_char())
            end select
        endif
        if( self%step1 == INIT_CLASSIFY )then
            select case( self%job%status() )
                case( ASYNC_JOB_IDLE )
                    call self%send_opening2D_status(string('classifying particles'), self%box, self%vis_cycle)
                    call start_solve2D(self%qenv, self%job, projfile, string('solve2D/init'),&
                        &self%spproj%os_ptcl2D%count_state_gt_zero(), self%params%nptcls_per_cls)
                case( ASYNC_JOB_DONE )
                    call self%spproj%kill()
                    call self%spproj%read(projfile)
                    call self%send_picking_status(string('complete'))
                    call self%job%kill()
                    self%step1 = INIT_SELECT
                case( ASYNC_JOB_FAILED )
                    job_log = self%job%get_log()
                    THROW_HARD('initial solve2D failed; see '//job_log%to_char())
            end select
        endif
        if( self%step1 == INIT_SELECT )then
            self%vis_cycle = 1
            call self%send_opening2D_status(string('evaluating class average quality'), self%box, self%vis_cycle)
            call self%select_and_send(1, string('quality_selection/init'))
            self%icycle = 2
        endif
    end subroutine run_cycle1

    ! Cycle 2 on the "all" set, once the sieve has every particle.
    subroutine run_cycle2( self )
        class(stream_stage_initial_analysis), intent(inout) :: self
        type(string) :: projfile, job_log
        projfile = self%cycle_projfile(2)
        if( self%step2 == ALL_SIEVE )then
            if( self%sieve%get_finished() )then
                write(logfhandle,'(A)') '>>> ALL SIEVE CHUNKS PROCESSED, COMBINING RESULTS...'
                call self%sieve%combine_completed_chunks(projfile, with_sigma2=.false.)
                call self%sieve%kill()
                self%l_sieve_active = .false.
                call self%spproj_all%kill()
                call self%spproj_all%read(projfile)
                self%step2 = ALL_CLASSIFY
            endif
        endif
        if( self%step2 == ALL_CLASSIFY )then
            select case( self%job%status() )
                case( ASYNC_JOB_IDLE )
                    call self%send_opening2D_status(string('classifying particles'), self%box, self%vis_cycle)
                    call start_solve2D(self%qenv, self%job, projfile, string('solve2D/all'),&
                        &self%spproj_all%os_ptcl2D%count_state_gt_zero(), self%params%nptcls_per_cls)
                case( ASYNC_JOB_DONE )
                    call self%spproj_all%kill()
                    call self%spproj_all%read(projfile)
                    call self%send_picking_status(string('complete'))
                    call self%job%kill()
                    self%step2 = ALL_SELECT
                case( ASYNC_JOB_FAILED )
                    job_log = self%job%get_log()
                    THROW_HARD('solve2D of the full set failed; see '//job_log%to_char())
            end select
        endif
        if( self%step2 == ALL_SELECT )then
            self%vis_cycle = 2
            call self%send_opening2D_status(string('evaluating class average quality'), self%box, self%vis_cycle)
            call self%select_and_send(2, string('quality_selection'))
            self%step2 = ALL_BALANCE
        endif
        if( self%step2 == ALL_BALANCE )then
            call self%send_opening2D_status(string('balancing classes'), self%box, self%vis_cycle)
            call balance_classes(self%spproj_all, projfile, string('balance_classes/all'))
            self%step2 = ALL_SOLVE3D
        endif
        if( self%step2 == ALL_SOLVE3D )then
            select case( self%job%status() )
                case( ASYNC_JOB_IDLE )
                    call self%send_opening2D_status(string('solve3D and reproject'), self%box, self%vis_cycle)
                    call start_solve3D(self%qenv, self%job, projfile, string('solve3D/all'), nint(self%mskdiam))
                case( ASYNC_JOB_DONE )
                    call self%finish_solve3D(projfile, string('solve3D/all'))
                    call self%send_picking_status(string('complete'))
                    call self%job%kill()
                    self%step2  = ALL_DONE
                    self%l_done = .true.
                case( ASYNC_JOB_FAILED )
                    job_log = self%job%get_log()
                    THROW_HARD('solve3D_cavgs failed; see '//job_log%to_char())
            end select
        endif
    end subroutine run_cycle2

    ! The cycle 1 micrograph segment: the accepted micrographs of every imported project, in order.
    subroutine rebuild_init_mics( self )
        class(stream_stage_initial_analysis), intent(inout) :: self
        type(string), allocatable :: fnames(:)
        integer :: nappended
        fnames = unique_projnames(self%project_list)
        call self%spproj%os_mic%kill()
        call append_mics_from_projects(self%spproj%os_mic, fnames, .true., nappended)
    end subroutine rebuild_init_mics

    ! Class-average selection on the project of cycle @p icycle, and the selected classes to the GUI.
    subroutine select_and_send( self, icycle, outdir )
        class(stream_stage_initial_analysis), intent(inout) :: self
        integer,                              intent(in)    :: icycle
        class(string),                        intent(in)    :: outdir
        integer, allocatable :: inds(:)
        type(string)         :: projfile, stk, jpg
        integer              :: n_selected, xtiles, ytiles
        projfile = self%cycle_projfile(icycle)
        jpg      = self%cwd//'/'//outdir//'/quality_cavgs'//JPG_EXT
        if( icycle == 1 )then
            call select_project_cavgs(self%spproj, projfile, outdir, self%mskdiam, n_selected, inds, stk, xtiles, ytiles)
            if( .not. allocated(inds) ) return
            if( size(inds) > 0 ) call send_cavgs(self%pipe, self%meta_cavg2D, jpg, inds, stk, xtiles, ytiles,&
                &os_cls2D=self%spproj%os_cls2D)
        else
            call select_project_cavgs(self%spproj_all, projfile, outdir, self%mskdiam, n_selected, inds, stk, xtiles, ytiles)
            if( .not. allocated(inds) ) return
            if( size(inds) > 0 ) call send_cavgs(self%pipe, self%meta_cavg2D, jpg, inds, stk, xtiles, ytiles,&
                &os_cls2D=self%spproj_all%os_cls2D)
        endif
    end subroutine select_and_send

    ! Picks the solve3D_cavgs state with the widest view coverage (else the most populated
    ! one with a volume), reprojects it, and publishes the reprojections, rescaled to the
    ! particle sampling, as the picking references; the volume and the references go to the GUI.
    subroutine finish_solve3D( self, projfile, outdir )
        class(stream_stage_initial_analysis), intent(inout) :: self
        class(string),                        intent(in)    :: projfile, outdir
        integer, allocatable      :: states(:), projs(:), nunique_proj(:), uniqbuf(:)
        type(commander_reproject) :: xreproject
        type(cmdline)             :: cline_reproject
        type(string)              :: cwd, volpath, final_dir, reprojdir, volpath_abs, empty_path
        type(image)               :: vol_shape
        type(image_bin)           :: mskvol_shape
        type(gui_metadata_vol3D)  :: meta_vol3D
        integer :: ldim(3), ldim_new(3), nvols, nuniq, i_cls3d, proj_here, xtiles, ytiles
        integer :: ivol, bestvol, pop, bestpop
        real    :: vol_smpd, minval3D, maxval3D, smpd_part
        logical :: l_published
        call simple_getcwd(cwd)
        call simple_chdir(outdir)
        call find_final_solve3D_cavgs_dir(final_dir)
        call self%spproj_all%kill()
        if( final_dir%strlen() > 0 )then
            ! mkdir=yes makes a numbered '<n>_solve3D_cavgs' directory per restart; the last one holds the result
            write(logfhandle,'(A,A)') '>>> SOLVE3D_CAVGS RESTART OUTPUT DIRECTORY: ', final_dir%to_char()
            call simple_chdir(final_dir)
            call self%spproj_all%read(basename(projfile))
            reprojdir = cwd//'/'//outdir//'/'//final_dir
        else
            call self%spproj_all%read(projfile)
            reprojdir = cwd//'/'//outdir
        endif
        ! shape descriptors of every populated state volume
        do ivol = 1, NSTATES3D
            if( self%spproj_all%os_cls3D%get_pop(ivol, 'state') == 0 ) cycle
            volpath = string('recvol_state'//int2str_pad(ivol,2)//MRC_EXT)
            if( .not. file_exists(volpath) ) cycle
            call find_ldim_nptcls(volpath, ldim, nuniq)
            call vol_shape%new(ldim, find_img_smpd(volpath))
            call vol_shape%read(volpath)
            write(logfhandle,'(A,I0)') '>>> VOLUME SHAPE DESCRIPTORS FOR STATE=', ivol
            call mskvol_shape%vol_shape_descr(vol_shape, 20.0, real(nint(self%mskdiam)))
            call mskvol_shape%kill_bimg
            call vol_shape%kill
        enddo
        ! fallback when the cls3D arrays cannot rank the states: the most populated state with a volume
        bestvol = 0
        bestpop = 0
        do ivol = 1, NSTATES3D
            pop = self%spproj_all%os_cls3D%get_pop(ivol, 'state')
            if( pop <= bestpop ) cycle
            volpath = string('recvol_state'//int2str_pad(ivol,2)//MRC_EXT)
            if( .not. file_exists(volpath) ) cycle
            bestvol = ivol
            bestpop = pop
        enddo
        ! the state whose classes cover the most distinct projection directions
        if( self%spproj_all%os_cls3D%isthere('state') .and. self%spproj_all%os_cls3D%isthere('proj') )then
            states = self%spproj_all%os_cls3D%get_all_asint('state')
            projs  = self%spproj_all%os_cls3D%get_all_asint('proj')
            nvols  = 0
            if( size(states) > 0 ) nvols = maxval(states)
            if( size(states) == size(projs) .and. nvols >= 1 )then
                allocate(nunique_proj(nvols), source=0)
                allocate(uniqbuf(max(1, size(projs))), source=0)
                do ivol = 1, nvols
                    nuniq = 0
                    do i_cls3d = 1, size(states)
                        if( states(i_cls3d) /= ivol ) cycle
                        proj_here = projs(i_cls3d)
                        if( nuniq == 0 )then
                            nuniq = 1
                            uniqbuf(1) = proj_here
                        else if( .not. any(uniqbuf(1:nuniq) == proj_here) )then
                            nuniq = nuniq + 1
                            uniqbuf(nuniq) = proj_here
                        endif
                    enddo
                    nunique_proj(ivol) = nuniq
                enddo
                bestvol = maxloc(nunique_proj, 1)
                write(logfhandle,'(A,I0,A,I0,A)') '>>> BEST VOLUME BY UNIQUE PROJ IN CLS3D: STATE=', bestvol, &
                    ' (NUNIQUE_PROJ=', nunique_proj(bestvol), ')'
            else
                write(logfhandle,'(A)') '>>> WARNING: cannot rank volumes by unique proj, nonconforming or empty cls3D arrays'
                write(logfhandle,'(A,I0)') '>>> FALLING BACK TO THE MOST POPULATED VOLUME, STATE=', bestvol
            endif
        else
            write(logfhandle,'(A)') '>>> WARNING: missing cls3D state/proj or empty volume set; unique-proj ranking skipped'
            write(logfhandle,'(A,I0)') '>>> FALLING BACK TO THE MOST POPULATED VOLUME, STATE=', bestvol
        endif
        if( bestvol == 0 ) THROW_HARD('No populated solve3D state with a reconstructed volume')
        volpath = string('recvol_state'//int2str_pad(bestvol,2)//MRC_EXT)
        if( .not. file_exists(volpath) ) THROW_HARD('Expected solve3D output volume not found: '//volpath%to_char())
        call find_ldim_nptcls(volpath, ldim, nuniq)
        vol_smpd = find_img_smpd(volpath)
        ! reprojections of the chosen state
        call cline_reproject%set('prg',     'reproject')
        call cline_reproject%set('vol1',    volpath)
        call cline_reproject%set('smpd',    vol_smpd)
        call cline_reproject%set('nspace',  50)
        call cline_reproject%set('pgrp',    'c1')
        call cline_reproject%set('nthr',    8)
        call cline_reproject%set('mskdiam', nint(self%mskdiam))
        call cline_reproject%printline()
        call xreproject%execute(cline_reproject)
        call cline_reproject%kill()
        call mrc2jpeg_tiled(string('reprojs.mrcs'), string('reprojs'//JPG_EXT), n_xtiles=xtiles, n_ytiles=ytiles)
        ! the chosen volume to the GUI (no FSC or postprocessed products at this stage)
        call meta_vol3D%new(GUI_METADATA_STREAM_OPENING2D_VOL3D_TYPE)
        volpath_abs = simple_abspath(volpath)
        empty_path  = string('')
        call meta_vol3D%set(reprojdir//'/reprojs'//JPG_EXT, volpath_abs, empty_path, empty_path, empty_path, &
            &bestvol, ldim(1), vol_smpd, 1, 1)
        call get_mrc_minmax(volpath, minval3D, maxval3D)
        call meta_vol3D%set_minmax('volpath', minval3D, maxval3D)
        call self%pipe%send_meta(meta_vol3D)
        ! reprojections rescaled to the particle sampling are the picking references
        smpd_part   = self%spproj_all%os_stk%get(1, 'smpd')
        ldim_new(1) = round2even(real(ldim(1)) * vol_smpd / smpd_part)
        ldim_new(2) = ldim_new(1)
        ldim_new(3) = 1
        write(logfhandle,'(A,I0,A,I0,A)') '>>> RESCALING AND CLIPPING REPROJECTIONS TO ', ldim_new(1), ' PIXEL BOX (',&
            &self%spproj_all%os_stk%get_int(1, 'box'), ' A) FOR PICKING REFERENCES'
        call scale_imgfile(string('reprojs.mrcs'), string('reprojs_rescaled.mrcs'), vol_smpd, ldim_new, smpd_part)
        call self%publish_pickrefs(simple_abspath(string('reprojs_rescaled.mrcs')), 'SOLVE3D', l_published)
        call meta_vol3D%kill()
        call simple_chdir(cwd)
        call self%spproj_all%kill()
        call self%spproj_all%read(projfile)
    end subroutine finish_solve3D

    !---------------- GUI ----------------

    ! Drains the GUI updates; a selection of picking references that publishes them ends the stage.
    subroutine apply_gui_updates( self )
        class(stream_stage_initial_analysis), intent(inout) :: self
        type(gui_metadata_stream_update) :: update
        character(len=:), allocatable    :: buffer
        logical :: l_published
        do while( self%pipe%receive(buffer) )
            update = transfer(buffer, update)
            if( update%get_pickrefs_cycle() > 0 .and. update%get_pickrefs_selection_length() > 0 )then
                call self%save_pickrefs_selection(update%get_pickrefs_selection(), update%get_pickrefs_cycle(), l_published)
                if( l_published )then
                    self%l_done = .true.
                    return
                endif
            endif
        enddo
    end subroutine apply_gui_updates

    ! Publishes the GUI-selected classes of cycle @p icycle (1-based class indices) as the picking
    ! references. @p l_published is .false. when nothing is published: the cycle has no class
    ! averages yet, or no index is one of its classes.
    subroutine save_pickrefs_selection( self, selection, icycle, l_published )
        class(stream_stage_initial_analysis), intent(inout) :: self
        integer,                              intent(in)    :: selection(:)
        integer,                              intent(in)    :: icycle
        logical,                              intent(out)   :: l_published
        type(image), allocatable :: cavg_imgs(:)
        integer,     allocatable :: states(:)
        type(sp_project)         :: spproj_sel
        type(string)             :: projfile, stk, stk_dummy
        real                     :: smpd_dummy
        integer                  :: ncls, icls
        l_published = .false.
        ncls        = 0
        projfile    = self%cycle_projfile(icycle)
        if( file_exists(projfile) )then
            call spproj_sel%read(projfile)
            call spproj_sel%get_cavgs_stk(stk_dummy, ncls, smpd_dummy, fail=.false.)
        endif
        if( ncls <= 0 )then
            write(logfhandle,'(A,I0)') '>>> WARNING: no class averages to select picking references from, cycle=', icycle
            call spproj_sel%kill()
            return
        endif
        allocate(states(ncls), source=0)
        do icls = 1, size(selection)
            if( selection(icls) >= 1 .and. selection(icls) <= ncls ) states(selection(icls)) = 1
        enddo
        if( count(states > 0) == 0 )then
            write(logfhandle,'(A,I0)') '>>> WARNING: the picking reference selection names no class of cycle ', icycle
            call spproj_sel%kill()
            return
        endif
        cavg_imgs = read_cavgs_into_imgarr(spproj_sel)
        stk       = self%cwd//'/'//PICKREFS_SELECTION
        call write_cavg_stack(cavg_imgs, states > 0, stk)
        call dealloc_imgarr(cavg_imgs)
        call spproj_sel%kill()
        call self%publish_pickrefs(stk, 'THE USER SELECTION', l_published)
        if( .not. l_published ) call del_file(stk)
    end subroutine save_pickrefs_selection

    ! Publishes the complete stack @p stk (made by @p source) as the picking references: it is
    ! renamed into place, so reference picking never reads a partial stack, and sent to the GUI.
    ! References are published once per run: @p l_published is .false., and @p stk left where it
    ! is, when they are published already.
    subroutine publish_pickrefs( self, stk, source, l_published )
        class(stream_stage_initial_analysis), intent(inout) :: self
        class(string),                        intent(in)    :: stk
        character(len=*),                     intent(in)    :: source
        logical,                              intent(out)   :: l_published
        type(string) :: pickrefs
        pickrefs    = self%cwd//'/'//OPENING2D_PICKREFS
        l_published = .not. file_exists(pickrefs)
        if( .not. l_published )then
            write(logfhandle,'(A)') '>>> PICKING REFERENCES ALREADY PUBLISHED; THOSE FROM '//source//' ARE NOT USED'
            return
        endif
        call simple_rename(stk, pickrefs)
        write(logfhandle,'(A)') '>>> PICKING REFERENCES PUBLISHED FROM '//source//': '//pickrefs%to_char()
        call self%send_pickrefs()
    end subroutine publish_pickrefs

    ! The published picking references and their sprite sheet to the GUI.
    subroutine send_pickrefs( self )
        class(stream_stage_initial_analysis), intent(inout) :: self
        integer, allocatable :: ref_inds(:)
        type(string)         :: pickrefs, jpg
        integer              :: nrefs, xtiles, ytiles, i
        pickrefs = self%cwd//'/'//OPENING2D_PICKREFS
        jpg      = self%cwd//'/'//swap_suffix(OPENING2D_PICKREFS, JPG_EXT, STK_EXT)
        call mrc2jpeg_tiled(pickrefs, jpg, ntiles=nrefs, n_xtiles=xtiles, n_ytiles=ytiles)
        ref_inds = [(i, i=1,nrefs)]
        call send_cavgs(self%pipe, self%meta_pickrefs, jpg, ref_inds, pickrefs, xtiles, ytiles)
    end subroutine send_pickrefs

    ! Picking progress, counted over the whole run (both the cycle 1 project and the "all" set).
    subroutine send_picking_status( self, stage )
        class(stream_stage_initial_analysis), intent(inout) :: self
        type(string),                         intent(in)    :: stage
        integer :: nmics, ntarget, nptcls
        nmics   = max(self%project_list%size(), self%spproj%os_mic%count_state_gt_zero())
        ntarget = max(self%nmics_target, nmics)
        nptcls  = max(self%extracted_project_list%get_nptcls_tot(), self%spproj%os_ptcl2D%get_noris())
        call self%meta_picking%set(stage=stage, micrographs_imported=ntarget, micrographs_accepted=nmics,&
            &particles_extracted=nptcls, box_size=self%box)
        call self%pipe%send_meta(self%meta_picking)
    end subroutine send_picking_status

    ! 2D progress; sent even before any particle exists (smpd is then 0).
    subroutine send_opening2D_status( self, stage, box_size, icycle )
        class(stream_stage_initial_analysis), intent(inout) :: self
        type(string),                         intent(in)    :: stage
        integer,                              intent(in)    :: box_size, icycle
        real    :: smpd
        integer :: nptcls
        nptcls = self%spproj%os_ptcl2D%get_noris()
        smpd   = 0.0
        if( nptcls > 0 ) smpd = self%spproj%get_smpd()
        call self%meta_opening2D%set(stage=stage, particles_imported=nptcls,&
            &particles_accepted=self%spproj%os_ptcl2D%count_state_gt_zero(), mask_diam=nint(self%mskdiam),&
            &mask_scale=box_size*smpd, box_size=box_size, cycle=icycle)
        call self%pipe%send_meta(self%meta_opening2D)
    end subroutine send_opening2D_status

    !---------------- paths ----------------

    function cycle_projfile( self, icycle ) result( fname )
        class(stream_stage_initial_analysis), intent(in) :: self
        integer,                              intent(in) :: icycle
        type(string) :: fname
        if( icycle == 1 )then
            fname = self%cwd//'/'//DIR_STREAM//'init/init'//METADATA_EXT
        else
            fname = self%cwd//'/'//DIR_STREAM//'all/all'//METADATA_EXT
        endif
    end function cycle_projfile

    ! the local copy of the n-th project of the "all" set
    function all_projfile( self, n ) result( fname )
        class(stream_stage_initial_analysis), intent(in) :: self
        integer,                              intent(in) :: n
        type(string) :: fname
        fname = self%cwd//'/'//DIR_STREAM//'all/'//int2str_pad(n, 5)//METADATA_EXT
    end function all_projfile

    !---------------- jobs (no stage state: called with the stage's components) ----------------

    subroutine start_extract( qenv, job, projfile, outdir, box, part, nmics )
        class(qsys_env),      intent(inout) :: qenv
        type(qsys_async_job), intent(inout) :: job
        class(string),        intent(in)    :: projfile, outdir
        integer,              intent(in)    :: box, part, nmics
        type(cmdline) :: cline
        type(string)  :: server_address
        server_address = qenv%get_persistent_worker_server_address()
        call cline%set('prg',             'extract')
        call cline%set('box',             box)
        call cline%set('nthr',            4)
        call cline%set('nparts',          1)
        call cline%set('part',            part)
        call cline%set('mkdir',           'no')
        call cline%set('dir',             '.')
        call cline%set('stream',          'yes')
        call cline%set('fromp',           1)
        call cline%set('top',             nmics)
        call cline%set('projfile',        projfile)
        call cline%set('worker_priority', 'high')
        if( server_address%strlen() > 0 ) call cline%set('worker_server', server_address)
        call cline%printline()
        call job%start(qenv, cline, outdir, 'extract_'//int2str(part))
        call cline%kill
    end subroutine start_extract

    ! Re-reads the extracted project and rejects micrographs left without particles or box file.
    subroutine finish_extract( spproj, projfile, outdir )
        type(sp_project), intent(inout) :: spproj
        class(string),    intent(in)    :: projfile, outdir
        type(string) :: cwd
        integer      :: nrejected
        call simple_getcwd(cwd)
        call simple_chdir(outdir)
        call spproj%kill()
        call spproj%read(projfile)
        call reject_mics_without_particles(spproj%os_mic, nrejected)
        call spproj%write(projfile)
        call simple_chdir(cwd)
    end subroutine finish_extract

    subroutine start_solve2D( qenv, job, projfile, outdir, nptcls, nptcls_per_cls )
        class(qsys_env),      intent(inout) :: qenv
        type(qsys_async_job), intent(inout) :: job
        class(string),        intent(in)    :: projfile, outdir
        integer,              intent(in)    :: nptcls, nptcls_per_cls
        type(cmdline) :: cline
        type(string)  :: server_address
        integer       :: ncls_job, nsample_job
        ncls_job    = min(NCLS_MAX, max(NCLS_MIN, nptcls/nptcls_per_cls))
        nsample_job = ((nptcls/5 + 999) / 1000) * 1000 ! round up to the nearest 1000
        server_address = qenv%get_persistent_worker_server_address()
        call simple_mkdir('solve2D')
        call cline%set('prg',             'solve2D')
        call cline%set('mkdir',           'no')
        call cline%set('ncls',            ncls_job)
        call cline%set('sigma_est',       'global')
        call cline%set('center',          'yes')
        call cline%set('autoscale',       'yes')
        call cline%set('nsample',         max(NSAMPLE2D, nsample_job))
        call cline%set('lpstop',          LPSTOP2D)
        call cline%set('mskdiam',         999.)
        call cline%set('nthr',            16)
        call cline%set('nparts',          1)
        call cline%set('projfile',        projfile)
        call cline%set('worker_priority', 'high')
        call cline%set('cache',           'yes')
        if( server_address%strlen() > 0 ) call cline%set('worker_server', server_address)
        call cline%printline()
        call job%start(qenv, cline, outdir, 'solve2D', exec_bin=string('simple_exec'))
        call cline%kill
    end subroutine start_solve2D

    subroutine start_solve3D( qenv, job, projfile, outdir, mskdiam )
        class(qsys_env),      intent(inout) :: qenv
        type(qsys_async_job), intent(inout) :: job
        class(string),        intent(in)    :: projfile, outdir
        integer,              intent(in)    :: mskdiam
        type(cmdline) :: cline
        call simple_mkdir('solve3D')
        call cline%set('prg',                'solve3D_cavgs')
        call cline%set('pgrp',               'c1')
        call cline%set('nstates',            NSTATES3D)
        call cline%set('lpstop',             8)
        call cline%set('mskdiam',            mskdiam)
        call cline%set('lpstart_ini3D',      100)
        call cline%set('lpstop_ini3D',       20)
        call cline%set('prune',              'no')
        call cline%set('nthr',               16)
        call cline%set('nstages',            3)
        call cline%set('nrestarts_collapse', 3)
        call cline%set('projfile',           projfile)
        call cline%printline()
        call job%start(qenv, cline, outdir, 'solve3D', exec_bin=string('simple_exec'))
        call cline%kill
    end subroutine start_solve3D

    !---------------- class averages ----------------

    ! Scores the class averages of @p spproj (chunk model, with their pixel size for the mask
    ! radius of the relational feature), maps the selection to the particles, filters by class
    ! compatibility, and writes the project; returns the selected classes' indices, stack and
    ! sprite-sheet layout for the GUI.
    subroutine select_project_cavgs( spproj, projfile, outdir, mskdiam, n_selected, cavg_inds, cavgs_stk, xtiles, ytiles )
        type(sp_project),     intent(inout) :: spproj
        class(string),        intent(in)    :: projfile, outdir
        real,                 intent(in)    :: mskdiam
        integer,              intent(out)   :: n_selected, xtiles, ytiles
        integer, allocatable, intent(inout) :: cavg_inds(:)
        type(string),         intent(inout) :: cavgs_stk
        type(image), allocatable    :: cavg_imgs(:)
        type(cavg_quality_model)    :: model
        type(cavg_quality_result)   :: quality
        type(class_compatibility)   :: compatibility
        type(support_model_metrics) :: metrics
        type(string)                :: cwd
        integer :: ncls, nrejected
        real    :: smpd
        n_selected = 0
        xtiles     = 0
        ytiles     = 0
        if( allocated(cavg_inds) ) deallocate(cavg_inds)
        call simple_getcwd(cwd)
        call simple_mkdir('quality_selection')
        call simple_mkdir(outdir)
        call simple_chdir(outdir)
        call spproj%get_cavgs_stk(cavgs_stk, ncls, smpd, fail=.false.)
        if( ncls <= 0 )then
            write(logfhandle,'(A)') '>>> WARNING: no class averages available for quality selection; skipping'
            call simple_chdir(cwd)
            return
        endif
        call model%init_preset(CAVG_QUALITY_MODEL_CHUNK_DEFAULT)
        call score_project_cavgs(spproj, model, mskdiam, cavg_imgs, quality, smpd=smpd)
        call model%kill()
        n_selected = count(quality%states > 0)
        call write_cavg_selection_stacks(cavg_imgs, quality%states, string('quality_selected_cavgs'//MRC_EXT),&
            &string('quality_rejected_cavgs'//MRC_EXT))
        call spproj%map_cavgs_selection(quality%states)
        call compatibility%new()
        call compatibility%train(spproj)
        call compatibility%get_support_model_metrics(metrics)
        write(logfhandle,'(A,3(F10.4,1X),A,3(F10.4,1X),A,L1,A,L1,A,L1)') '>>> COARSE COMPAT METRICS a/b/c=',&
            &metrics%axis_a, metrics%axis_b, metrics%axis_c, ' da/db/dc=', metrics%delta_a, metrics%delta_b, metrics%delta_c,&
            &' valid=', metrics%valid, ' delta_valid=', metrics%delta_valid, ' converged=', metrics%converged
        call compatibility%infer(spproj)
        call compatibility%kill()
        call reject_mics_without_particles(spproj%os_mic, nrejected)
        call spproj%cavgs2jpg(cavg_inds, string('quality_cavgs')//JPG_EXT, xtiles, ytiles, ignore_states=.false.)
        if( allocated(cavg_inds) ) cavg_inds = pack(cavg_inds, cavg_inds > 0) ! unselected classes are 0
        call dealloc_imgarr(cavg_imgs)
        call spproj%write(projfile)
        call simple_chdir(cwd)
    end subroutine select_project_cavgs

    ! Replicates the selected class averages (and their even/odd stacks) in proportion to
    ! their populations, up to TARGET_NCLS rows; see the module header on this method.
    subroutine balance_classes( spproj, projfile, outdir )
        type(sp_project), intent(inout) :: spproj
        class(string),    intent(in)    :: projfile, outdir
        type(string)              :: cavgsstk, balanced_stk, odd_stk, even_stk, sigma2_stk
        type(string)              :: odd_balanced_stk, even_balanced_stk, sigma2_balanced_stk, cwd
        type(image), allocatable  :: cavg_imgs(:)
        type(oris)                :: os_cls2D_src
        integer,     allocatable  :: src_inds(:), pops(:), reps(:), extra(:)
        real,        allocatable  :: frac(:)
        logical,     allocatable  :: picked(:)
        integer :: ncls_all, nsrc, total_pop, icls, j, out_ind, istk, rem, imax, iout, n_balanced
        real    :: smpd_dummy
        call simple_getcwd(cwd)
        call simple_mkdir('balance_classes')
        call simple_mkdir(outdir)
        call simple_chdir(outdir)
        call spproj%get_cavgs_stk(cavgsstk, ncls_all, smpd_dummy, out_ind=out_ind, fail=.false.)
        if( ncls_all <= 0 )then
            write(logfhandle,'(A)') '>>> WARNING: no class averages available for balancing; skipping'
            call simple_chdir(cwd)
            return
        endif
        allocate(src_inds(spproj%os_cls2D%count_state_gt_zero()))
        nsrc = 0
        do icls = 1, spproj%os_cls2D%get_noris()
            if( nint(spproj%os_cls2D%get(icls, 'state')) > 0 )then
                nsrc = nsrc + 1
                src_inds(nsrc) = icls
            endif
        enddo
        if( nsrc <= 0 .or. nsrc >= TARGET_NCLS )then
            call simple_chdir(cwd)
            return
        endif
        cavg_imgs = read_cavgs_into_imgarr(spproj)
        call os_cls2D_src%new(ncls_all, is_ptcl=.false.)
        do icls = 1, ncls_all
            call os_cls2D_src%transfer_ori(icls, spproj%os_cls2D, icls)
        enddo
        ! replication counts proportional to population, the remainder to the largest fractions
        allocate(pops(nsrc), reps(nsrc), extra(nsrc), frac(nsrc), picked(nsrc))
        pops = 1
        do icls = 1, nsrc
            if( spproj%os_cls2D%isthere(src_inds(icls), 'pop') )then
                pops(icls) = max(1, nint(spproj%os_cls2D%get(src_inds(icls), 'pop')))
            endif
        enddo
        total_pop = sum(pops)
        rem   = TARGET_NCLS - nsrc
        frac  = real(rem) * real(pops) / real(total_pop)
        extra = int(frac)
        reps  = 1 + extra
        rem   = rem - sum(extra)
        if( rem > 0 )then
            frac   = frac - real(extra)
            picked = .false.
            do j = 1, rem
                imax = maxloc(frac, dim=1, mask=.not. picked)
                reps(imax)   = reps(imax) + 1
                picked(imax) = .true.
            enddo
        endif
        n_balanced   = sum(reps)
        balanced_stk = 'cavgs_balanced'//MRC_EXT
        if( file_exists(balanced_stk) ) call del_file(balanced_stk)
        call spproj%os_cls2D%new(n_balanced, is_ptcl=.false.)
        istk = 0
        do icls = 1, nsrc
            do j = 1, reps(icls)
                istk = istk + 1
                call cavg_imgs(src_inds(icls))%write(balanced_stk, istk)
                call spproj%os_cls2D%transfer_ori(istk, os_cls2D_src, src_inds(icls))
                call spproj%os_cls2D%set(istk, 'indstk', istk)
                call spproj%os_cls2D%set(istk, 'state', 1.)
            enddo
        enddo
        if( spproj%os_cls3D%get_noris() /= n_balanced )then
            call spproj%os_cls3D%new(n_balanced, is_ptcl=.false.)
            do icls = 1, n_balanced
                call spproj%os_cls3D%transfer_ori(icls, spproj%os_cls2D, icls)
            enddo
        endif
        call spproj%os_out%set(out_ind, 'stk',        simple_abspath(balanced_stk))
        call spproj%os_out%set(out_ind, 'nptcls',     n_balanced)
        call spproj%os_out%set(out_ind, 'nptcls_stk', n_balanced)
        call spproj%os_out%set(out_ind, 'fromp',      1)
        call spproj%os_out%set(out_ind, 'top',        n_balanced)
        ! even/odd class-average stacks, when present, follow the same replication
        odd_stk = swap_suffix(cavgsstk, string('_odd'//MRC_EXT), string(MRC_EXT))
        if( file_exists(odd_stk) )then
            odd_balanced_stk = 'cavgs_balanced_odd'//MRC_EXT
            call duplicate_balanced_stack(odd_stk, odd_balanced_stk, ncls_all, nsrc, src_inds, reps)
            call update_os_out_stk(spproj, odd_stk, odd_balanced_stk, n_balanced)
        endif
        even_stk = swap_suffix(cavgsstk, string('_even'//MRC_EXT), string(MRC_EXT))
        if( file_exists(even_stk) )then
            even_balanced_stk = 'cavgs_balanced_even'//MRC_EXT
            call duplicate_balanced_stack(even_stk, even_balanced_stk, ncls_all, nsrc, src_inds, reps)
            call update_os_out_stk(spproj, even_stk, even_balanced_stk, n_balanced)
        endif
        ! the sigma2 output is carried forward under the balanced name
        do iout = 1, spproj%os_out%get_noris()
            if( .not. spproj%os_out%isthere(iout, 'imgkind') ) cycle
            if( spproj%os_out%get_str(iout, 'imgkind') /= 'sigma2' ) cycle
            if( .not. spproj%os_out%isthere(iout, 'sigma2') ) cycle
            sigma2_stk = spproj%os_out%get_str(iout, 'sigma2')
            if( .not. file_exists(sigma2_stk) ) cycle
            sigma2_balanced_stk = 'cavgs_balanced'//STAR_EXT
            call simple_copy_file(sigma2_stk, sigma2_balanced_stk)
            call spproj%os_out%set(iout, 'sigma2', simple_abspath(sigma2_balanced_stk))
        enddo
        call spproj%write(projfile)
        call simple_chdir(cwd)
        write(logfhandle, '(A,I0,A)') '>>> BALANCED CLASS AVERAGES TO ', n_balanced, ' ENTRIES'
        call os_cls2D_src%kill()
        call dealloc_imgarr(cavg_imgs)
    end subroutine balance_classes

    subroutine duplicate_balanced_stack( stk_in, stk_out, ncls_all, nsrc, src_inds, reps )
        type(string), intent(in) :: stk_in, stk_out
        integer,      intent(in) :: ncls_all, nsrc
        integer,      intent(in) :: src_inds(:), reps(:)
        type(image), allocatable :: imgs_in(:)
        integer :: ii, jj, kk
        imgs_in = read_stk_into_imgarr(stk_in)
        if( size(imgs_in) /= ncls_all )then
            call dealloc_imgarr(imgs_in)
            return
        endif
        if( file_exists(stk_out) ) call del_file(stk_out)
        kk = 0
        do ii = 1, nsrc
            do jj = 1, reps(ii)
                kk = kk + 1
                call imgs_in(src_inds(ii))%write(stk_out, kk)
            enddo
        enddo
        call dealloc_imgarr(imgs_in)
    end subroutine duplicate_balanced_stack

    subroutine update_os_out_stk( spproj, old_stk, new_stk, nstk )
        type(sp_project), intent(inout) :: spproj
        type(string),     intent(in)    :: old_stk, new_stk
        integer,          intent(in)    :: nstk
        type(string) :: stk_here
        integer      :: io
        do io = 1, spproj%os_out%get_noris()
            if( .not. spproj%os_out%isthere(io, 'stk') ) cycle
            stk_here = spproj%os_out%get_str(io, 'stk')
            if( simple_abspath(stk_here) /= simple_abspath(old_stk) ) cycle
            call spproj%os_out%set(io, 'stk',        simple_abspath(new_stk))
            call spproj%os_out%set(io, 'nptcls',     nstk)
            call spproj%os_out%set(io, 'nptcls_stk', nstk)
            call spproj%os_out%set(io, 'fromp',      1)
            call spproj%os_out%set(io, 'top',        nstk)
        enddo
    end subroutine update_os_out_stk

    ! The highest-numbered '<n>_solve3D_cavgs' restart directory under the working
    ! directory, or an empty string when there is none.
    subroutine find_final_solve3D_cavgs_dir( final_dir )
        type(string), intent(out) :: final_dir
        character(len=*), parameter   :: SUFFIX = '_solve3D_cavgs'
        type(string),     allocatable :: dirs(:)
        character(len=:), allocatable :: dname
        integer :: idir, us, n, best_n, io_stat
        final_dir = ''
        dirs = simple_list_dirs('.')
        if( .not. allocated(dirs) ) return
        best_n = -1
        do idir = 1, size(dirs)
            dname = trim(dirs(idir)%to_char())
            us = index(dname, SUFFIX)
            if( us < 2 ) cycle                 ! at least one digit before the suffix
            if( dname(us:) /= SUFFIX ) cycle    ! the suffix ends the name
            n = str2int(dname(1:us-1), io_stat)
            if( io_stat /= 0 ) cycle
            if( n > best_n )then
                best_n    = n
                final_dir = dirs(idir)
            endif
        enddo
    end subroutine find_final_solve3D_cavgs_dir

    ! Project files of a record list, once each, in order; a project's records are consecutive.
    function unique_projnames( list ) result( fnames )
        type(rec_list), intent(inout) :: list
        type(string), allocatable :: fnames(:)
        type(project_rec) :: prec
        integer :: i
        allocate(fnames(0))
        do i = 1, list%size()
            call list%at(i, prec)
            if( size(fnames) > 0 )then
                if( fnames(size(fnames))%to_char() == prec%projname%to_char() ) cycle
            endif
            fnames = [fnames, prec%projname]
        enddo
    end function unique_projnames

end module simple_stream_stage_initial_analysis
