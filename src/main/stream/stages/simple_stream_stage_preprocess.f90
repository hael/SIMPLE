!@descr: state and steps of stream task 1 (preprocessing): watch movies, schedule preprocess jobs, import results, report to the GUI
!==============================================================================
! MODULE: simple_stream_stage_preprocess
!
! PURPOSE:
!   The body of stream p01 as a type: the state the stage carries between
!   iterations is explicit, and each step is a type-bound procedure that can
!   be called, and tested, on its own. The commander
!   (simple_commanders_stream_p01_preprocess) only normalises the command line and loops
!   over iterate() until finished().
!
!   The stage decides when things happen (watching, scheduling, importing,
!   reporting). What happens is delegated:
!     - micrograph rejection    -> simple_mic_selection
!     - GUI plots               -> simple_stream_meta_plots
!     - gain flip / generation  -> gain_flip_analyzer, simple_motion_gain_helpers, flip_gain
!       (the stage only supplies batches of movies)
!     - pipe framing            -> simple_stream_pipe
!     - movie sets and their jobs (naming, submission, collection, restart)
!                               -> simple_stream_job_sets
!
! LIFECYCLE:
!   new(cline) -> { iterate() } until finished() -> finalize() -> kill()
!
! PARAMETERS:
!   params is built once, in new(). After that only values known at run time
!   are changed (the resolved gain reference, GUI threshold updates), each in
!   one place that also updates the command line the worker jobs receive.
!
! RESTART:
!   Recognised by the execution folder existing (outdir, or dir_exec). Every
!   movie of the completed sets, accepted or rejected, goes into the watcher's
!   history, and import indices continue after the highest given; the accepted
!   micrographs are imported again. The thresholds the GUI set come back from
!   gui_thresholds.txt, and a generated gain reference is reused. Set
!   numbering continues (simple_stream_job_sets); unfinished sets are dropped
!   and their movies submitted again.
!==============================================================================
module simple_stream_stage_preprocess
use simple_defs,                           only: logfhandle, STDLEN, PATH_HERE
use simple_defs_fname,                     only: DIR_CTF_ESTIMATE, DIR_MOTION_CORRECT, TERM_STREAM, GAIN_THUMBNAIL
use simple_defs_stream,                    only: DIR_STREAM, DIR_STREAM_COMPLETED, STREAM_NMOVS_SET, LONGTIME,&
                                                &SHORTWAIT, WAITTIME, INACTIVE_TIME, CTFRES_BINS, ICESCORE_BINS, ASTIG_BINS,&
                                                &STREAM_IDLE_MARKER, STREAM_FINISHED_MARKER, MOVIES_IDLE_TIME_S
use simple_defs_environment,               only: SIMPLE_STREAM_PREPROC_PARTITION
use simple_type_defs,                      only: ctfparams, CTFFLAG_YES
use simple_error,                          only: simple_exception
use simple_string,                         only: string
use simple_fileio,                         only: basename, del_file, file_exists, filepath, fname2format, simple_touch,&
                                                &simple_abspath, simple_getcwd, stemname
use simple_syslib,                         only: simple_mkdir, dir_exists, simple_rename
use simple_timer,                          only: simple_gettime, cast_time_char
use simple_cmdline,                        only: cmdline
use simple_image,                          only: image
use simple_oris,                           only: oris
use simple_parameters,                     only: parameters
use simple_sp_project,                     only: sp_project
use simple_qsys_env,                       only: qsys_env
use simple_qsys_funs,                      only: qsys_cleanup
use simple_starproject_stream,             only: starproject_stream
use simple_stream_watcher,                 only: stream_watcher, sniff_folders_SJ, workout_directory_structure
use simple_stream_state,                   only: ipc_pipe_preprocess_in, ipc_pipe_preprocess_out
use simple_motion_correct_utils,           only: flip_gain
use simple_motion_gain_analysis,           only: gain_flip_analyzer
use simple_motion_gain_helpers,            only: gainref_to_jpg, read_movies_and_sum_frames, add_movies_to_gain_sum,&
                                                &write_gain_from_sum
use simple_gui_metadata_utils,             only: max_metadata_size
use simple_gui_metadata_types,             only: GUI_METADATA_STREAM_PREPROCESS_TYPE,&
                                                &GUI_METADATA_STREAM_PREPROCESS_MICROGRAPH_TYPE,&
                                                &GUI_METADATA_STREAM_PREPROCESS_HISTOGRAM_CTFRES_TYPE,&
                                                &GUI_METADATA_STREAM_PREPROCESS_HISTOGRAM_ICEFRAC_TYPE,&
                                                &GUI_METADATA_STREAM_PREPROCESS_HISTOGRAM_ASTIG_TYPE,&
                                                &GUI_METADATA_STREAM_PREPROCESS_TIMEPLOT_CTFRES_TYPE,&
                                                &GUI_METADATA_STREAM_PREPROCESS_TIMEPLOT_ASTIG_TYPE,&
                                                &GUI_METADATA_STREAM_PREPROCESS_TIMEPLOT_DF_TYPE,&
                                                &GUI_METADATA_STREAM_PREPROCESS_TIMEPLOT_RATE_TYPE
use simple_gui_metadata_micrograph,        only: gui_metadata_micrograph
use simple_gui_metadata_histogram,         only: gui_metadata_histogram
use simple_gui_metadata_timeplot,          only: gui_metadata_timeplot
use simple_gui_metadata_stream_update,     only: gui_metadata_stream_update
use simple_gui_metadata_stream_preprocess, only: gui_metadata_stream_preprocess
use simple_stream_pipe,                    only: stream_pipe
use simple_mic_selection,                  only: reject_mics_by_thresholds
use simple_mic_import,                     only: append_mics_from_projects
use simple_stream_job_sets,                only: stream_job_sets
use simple_stream_meta_plots,              only: set_histogram_from_oris, set_timeplot_from_oris, set_rate_timeplot
use simple_stream_sigterm,                 only: sigterm_received
implicit none

public :: stream_stage_preprocess
private
#include "simple_local_flags.inc"

! the micrograph rejection thresholds last set from the GUI, kept in the stage's folder for a restart
character(len=*), parameter :: GUI_THRESHOLDS = 'gui_thresholds.txt'

integer, parameter :: PLOT_WINDOW      = 500   ! micrographs per point of the windowed time plots
integer, parameter :: NTHUMBS          = 10    ! most recent micrograph thumbnails sent to the GUI
integer, parameter :: STAR_EVERY_NMICS = 1000  ! below this many micrographs the STAR file is rewritten on every import...
integer, parameter :: STAR_STEP_NMICS  = 100   ! ...above it, every this many new micrographs
integer, parameter :: GAIN_BATCH_NMOVIES = 10  ! movies per batch handed to the gain analysis/estimation
character(len=*), parameter :: GENERATED_GAINREF = 'gainref_generated.mrc'

! Components and steps are public so simple_stream_stage_preprocess_tester can assemble a
! stage and run one step at a time; production code uses new/iterate/finished/finalize/kill.
type :: stream_stage_preprocess
    type(parameters), allocatable        :: params
    type(cmdline)                        :: cline_exec            ! command line of the preprocess worker jobs
    type(qsys_env)                       :: qenv
    type(sp_project)                     :: spproj_glob           ! every imported micrograph
    type(stream_watcher)                 :: movie_buff
    type(starproject_stream)             :: starproj_stream
    type(stream_pipe)                    :: pipe                  ! to and from the master
    type(gui_metadata_stream_preprocess) :: meta_status
    type(gui_metadata_micrograph)        :: meta_micrograph
    type(gui_metadata_histogram)         :: meta_hist_ctfres, meta_hist_icefrac, meta_hist_astig
    type(gui_metadata_timeplot)          :: meta_plot_ctfres, meta_plot_astig, meta_plot_df, meta_plot_rate
    type(stream_job_sets)                :: sets                  ! the movie sets the worker jobs run on
    ! settled after params%new: the gain reference and its flip once resolved (resolve_gain), and
    ! the micrograph rejection thresholds the GUI may change (set_threshold)
    type(string)          :: gainref
    character(len=STDLEN) :: flipgain     = 'no'
    real                  :: ctfres_thres  = 0.
    real                  :: astig_thres   = 0.
    real                  :: icefrac_thres = 0.
    integer :: import_counter      = 0       ! last import index given to a movie
    integer :: last_movie_time     = 0       ! when a new movie was last seen
    integer :: nmic_star           = 0       ! micrographs in the last STAR snapshot (above STAR_EVERY_NMICS)
    integer :: n_failed_jobs       = 0
    integer :: prev_stacksz        = 0
    integer :: last_injection      = 0       ! time of the last import
    integer :: nmovs2importperiter = 0
    logical :: l_sj_dirs           = .false. ! movies arrive in <dir_movies>/<xx>/Data/ folders
    logical :: l_xml_meta          = .false. ! per-movie XML metadata in dir_meta
    logical :: l_movies_left       = .false. ! the watcher returned more movies than were submitted
    logical :: l_haschanged        = .false. ! imports since the last idle STAR snapshot
    logical :: l_nmics_reached     = .false.
    logical :: l_restart           = .false. ! the output directory existed before params%new
    logical :: l_exists            = .false.
    ! waits (s); tests set them to 0
    integer :: settle_s            = LONGTIME  ! a file is taken once untouched this long
    integer :: wait_s              = WAITTIME  ! idle pause between iterations
    integer :: sniff_wait_s        = SHORTWAIT ! pause while waiting for the first movie or a gain batch
contains
    procedure :: new
    procedure :: iterate
    procedure :: finished
    procedure :: finalize
    procedure :: kill
    ! the steps of new(), in order
    procedure :: init_params
    procedure :: init_movie_watcher
    procedure :: resume_previous_run
    procedure :: init_job_dirs
    procedure :: init_queue
    procedure :: resolve_gain
    procedure :: build_worker_cline
    procedure :: init_gui
    ! the steps of iterate() and their helpers
    procedure :: import_previous_projects
    procedure :: detect_gain_flip
    procedure :: generate_gain_from_movies
    procedure :: next_movie_batch
    procedure :: submit_new_movies
    procedure :: create_movies_set_project
    procedure :: schedule_jobs
    procedure :: collect_jobs
    procedure :: import_completed
    procedure :: process_imports
    procedure :: send_plots
    procedure :: idle
    procedure :: send_status
    procedure :: apply_gui_updates
    procedure :: set_threshold
    procedure :: apply_thresholds
    procedure :: update_idle_marker
    procedure :: save_thresholds
    procedure :: restore_thresholds
    procedure :: write_mic_star_and_field
end type stream_stage_preprocess

contains

    !---------------- lifecycle ----------------

    !> Builds the parameters from the normalised @p cline, attaches to the movie directory,
    !! re-imports a previous run when restarting, and prepares job submission.
    subroutine new( self, cline )
        class(stream_stage_preprocess), intent(inout) :: self
        class(cmdline),                 intent(inout) :: cline
        call self%kill()
        call self%init_params(cline)
        self%l_exists = .true. ! from here kill() releases what has been built
        call self%init_movie_watcher()
        if( sigterm_received() ) return ! stopped while waiting for the first movie
        call self%init_job_dirs()
        if( self%l_restart ) call self%resume_previous_run(cline)
        call self%init_queue()
        ! the GUI hears from the stage while it waits for the movies of the gain step
        call self%init_gui(ipc_pipe_preprocess_out(1), ipc_pipe_preprocess_in(2))
        call self%resolve_gain(cline)
        if( sigterm_received() ) return ! stopped during the gain step
        call self%build_worker_cline(cline)
        if( self%l_restart ) call self%restore_thresholds()
        self%last_injection      = simple_gettime()
        self%last_movie_time     = simple_gettime()
        self%nmovs2importperiter = 2 * self%params%nparts * STREAM_NMOVS_SET
        ! the markers downstream ends its intake on are those of this run
        call del_file(STREAM_IDLE_MARKER)
        call del_file(STREAM_FINISHED_MARKER)
    end subroutine new

    !> The stage's project file (created when missing) and its parameters; the global
    !! project must start empty.
    subroutine init_params( self, cline )
        class(stream_stage_preprocess), intent(inout) :: self
        class(cmdline),                 intent(inout) :: cline
        type(string) :: projfile
        ! a restart is recognised by its execution folder existing before params%new creates or
        ! enters it: outdir, or dir_exec (which simple_stream takes as the execution folder)
        self%l_restart = .false.
        if( cline%defined('outdir') )then
            self%l_restart = dir_exists(cline%get_carg('outdir'))
        endif
        if( cline%defined('dir_exec') )then
            self%l_restart = self%l_restart .or. dir_exists(cline%get_carg('dir_exec'))
        endif
        projfile = cline%get_carg('projfile')
        if( .not. file_exists(projfile) )then
            call self%spproj_glob%update_projinfo(cline)
            call self%spproj_glob%update_compenv(cline)
            call self%spproj_glob%write
        endif
        if( .not. allocated(self%params) ) allocate(self%params)
        ! one queue partition per computing unit; not passed on to the workers' command lines
        call cline%set('split_mode', 'stream')
        call self%params%new(cline)
        call cline%delete('split_mode')
        self%gainref       = self%params%gainref
        self%flipgain      = self%params%flipgain
        self%ctfres_thres  = self%params%ctfresthreshold
        self%astig_thres   = self%params%astigthreshold
        self%icefrac_thres = self%params%icefracthreshold
        self%l_xml_meta = cline%defined('dir_meta')
        call self%spproj_glob%read(self%params%projfile)
        if( self%spproj_glob%os_mic%get_noris() /= 0 )then
            THROW_HARD('PREPROCESS_STREAM must start from an empty project (eg from root project folder)')
        endif
    end subroutine init_params

    !> Restart: re-imports the accepted micrographs of the previous run and rewrites the STAR file.
    subroutine resume_previous_run( self, cline )
        class(stream_stage_preprocess), intent(inout) :: self
        class(cmdline),                 intent(inout) :: cline
        write(logfhandle,'(A)') '>>> RESTARTING EXISTING JOB'
        call del_file(TERM_STREAM)
        if( cline%defined('dir_exec') ) call cline%delete('dir_exec')
        call self%import_previous_projects()
        self%nmic_star = self%spproj_glob%os_mic%get_noris()
        call self%write_mic_star_and_field(write_field=.true.)
    end subroutine resume_previous_run

    !> The job sets (job and completed folders) and the directories the worker outputs use.
    subroutine init_job_dirs( self )
        class(stream_stage_preprocess), intent(inout) :: self
        call self%sets%new(string(PATH_HERE//DIR_STREAM), string(PATH_HERE//DIR_STREAM_COMPLETED), self%params%numlen)
        call simple_mkdir(filepath(PATH_HERE, DIR_CTF_ESTIMATE))
        call simple_mkdir(filepath(PATH_HERE, DIR_MOTION_CORRECT))
    end subroutine init_job_dirs

    !> The queue environment the preprocess jobs are submitted through.
    subroutine init_queue( self )
        class(stream_stage_preprocess), intent(inout) :: self
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
        class(stream_stage_preprocess), intent(inout) :: self
        integer,                        intent(in)    :: fd_read, fd_write
        call self%meta_status%new(GUI_METADATA_STREAM_PREPROCESS_TYPE)
        call self%meta_micrograph%new(GUI_METADATA_STREAM_PREPROCESS_MICROGRAPH_TYPE)
        call self%meta_hist_ctfres%new(GUI_METADATA_STREAM_PREPROCESS_HISTOGRAM_CTFRES_TYPE)
        call self%meta_hist_icefrac%new(GUI_METADATA_STREAM_PREPROCESS_HISTOGRAM_ICEFRAC_TYPE)
        call self%meta_hist_astig%new(GUI_METADATA_STREAM_PREPROCESS_HISTOGRAM_ASTIG_TYPE)
        call self%meta_plot_ctfres%new(GUI_METADATA_STREAM_PREPROCESS_TIMEPLOT_CTFRES_TYPE)
        call self%meta_plot_astig%new(GUI_METADATA_STREAM_PREPROCESS_TIMEPLOT_ASTIG_TYPE)
        call self%meta_plot_df%new(GUI_METADATA_STREAM_PREPROCESS_TIMEPLOT_DF_TYPE)
        call self%meta_plot_rate%new(GUI_METADATA_STREAM_PREPROCESS_TIMEPLOT_RATE_TYPE)
        call self%pipe%new(fd_read, fd_write, max_metadata_size(), 'preprocess')
    end subroutine init_gui

    !> One pass: submit new movie sets, import finished jobs, report, take GUI updates.
    subroutine iterate( self )
        class(stream_stage_preprocess), intent(inout) :: self
        integer :: n_imported
        call self%submit_new_movies()
        call self%schedule_jobs()
        call self%collect_jobs(n_imported)
        if( n_imported > 0 )then
            call self%process_imports()
        else if( .not. self%l_movies_left )then
            call self%idle()
        endif
        call self%update_idle_marker()
        call self%send_status()
        call self%apply_gui_updates()
        if( self%params%nmics > 0 )then
            if( self%spproj_glob%os_mic%get_noris() >= self%params%nmics .and. .not. self%l_nmics_reached )then
                write(logfhandle,'(A,I8)') '>>> TERMINATING PROCESS: requested number of micrographs reached: ', self%params%nmics
                self%l_nmics_reached = .true.
            endif
        endif
        call flush(logfhandle)
    end subroutine iterate

    !> .true. once the stream is told to stop or the requested number of micrographs is imported.
    logical function finished( self )
        class(stream_stage_preprocess), intent(in) :: self
        finished = self%l_nmics_reached .or. file_exists(TERM_STREAM)
    end function finished

    !> Cancels the jobs in flight, writes the final STAR file and removes the job scripts.
    subroutine finalize( self )
        class(stream_stage_preprocess), intent(inout) :: self
        call self%sets%cancel(self%qenv) ! a restart sets the folder of unfinished sets aside
        if( self%spproj_glob%os_mic%get_noris() > 0 )then
            call self%write_mic_star_and_field(write_field=.true., copy_optics=.true.)
        endif
        call qsys_cleanup(self%params)
        call simple_touch(STREAM_FINISHED_MARKER) ! downstream ends its intake
    end subroutine finalize

    subroutine kill( self )
        class(stream_stage_preprocess), intent(inout) :: self
        if( .not. self%l_exists )then
            if( allocated(self%params) ) deallocate(self%params)
            return
        endif
        call self%spproj_glob%kill
        call self%qenv%kill
        call self%movie_buff%kill
        call self%pipe%kill
        call self%cline_exec%kill
        call self%meta_status%kill
        call self%meta_micrograph%kill
        call self%gainref%kill
        self%flipgain      = 'no'
        self%ctfres_thres  = 0.
        self%astig_thres   = 0.
        self%icefrac_thres = 0.
        call self%meta_hist_ctfres%kill
        call self%meta_hist_icefrac%kill
        call self%meta_hist_astig%kill
        call self%meta_plot_ctfres%kill
        call self%meta_plot_astig%kill
        call self%meta_plot_df%kill
        call self%meta_plot_rate%kill
        call self%sets%kill
        if( allocated(self%params) ) deallocate(self%params)
        self%import_counter     = 0
        self%nmic_star          = 0
        self%n_failed_jobs      = 0
        self%prev_stacksz       = 0
        self%l_movies_left      = .false.
        self%l_haschanged       = .false.
        self%l_nmics_reached    = .false.
        self%l_restart          = .false.
        self%l_exists           = .false.
    end subroutine kill

    !---------------- set-up ----------------

    ! Waits for the first movie or movie folder, which tells the directory layout, and attaches the watcher.
    ! Returns without a watcher when SIGTERM arrives first.
    subroutine init_movie_watcher( self )
        class(stream_stage_preprocess), intent(inout) :: self
        logical :: l_dir_found
        integer :: nwaits
        l_dir_found = .false.
        nwaits      = 0
        do while( .not. l_dir_found )
            call workout_directory_structure(self%params%dir_movies, l_dir_found, self%l_sj_dirs)
            if( l_dir_found ) exit
            if( sigterm_received() ) return
            call sleep(self%sniff_wait_s)
            nwaits = nwaits + 1
            if( mod(nwaits*max(1,self%sniff_wait_s), 60) == 0 )then
                write(logfhandle,'(A,I3,A)') '>>> NO MOVIE HAS BEEN DETECTED FOR ', nint(real(nwaits*max(1,self%sniff_wait_s))/60.), ' MINS'
            endif
        enddo
        call new_movie_watcher(self%params%dir_movies, self%l_sj_dirs, self%settle_s, self%movie_buff, l_dir_found)
        if( .not. l_dir_found ) THROW_HARD('Fatal error directory structure')
    end subroutine init_movie_watcher

    ! Re-imports the accepted micrographs of a previous run and puts their movies in the watcher
    ! history; the set numbering continues after every completed set (job sets restore).
    subroutine import_previous_projects( self )
        class(stream_stage_preprocess), intent(inout) :: self
        type(string), allocatable :: completed_fnames(:)
        type(oris) :: os_done
        integer    :: nmics, ndone, imic
        call self%sets%restore(completed_fnames)
        if( size(completed_fnames) == 0 ) return
        ! every movie of the completed sets is done, accepted or rejected: none is processed again,
        ! and no import index already handed downstream is given again
        call append_mics_from_projects(os_done, completed_fnames, .false., ndone)
        do imic = 1,os_done%get_noris()
            call self%movie_buff%add2history(os_done%get_str(imic, 'movie'))
            if( os_done%isthere(imic, 'importind') )then
                self%import_counter = max(self%import_counter, os_done%get_int(imic, 'importind'))
            endif
        enddo
        call os_done%kill()
        ! the accepted micrographs make the stage's project again
        call append_mics_from_projects(self%spproj_glob%os_mic, completed_fnames, .true., nmics)
        if( nmics == 0 ) write(logfhandle,'(A)') '>>> NO ACCEPTED MICROGRAPHS IN THE PREVIOUS RUN'
        write(logfhandle,'(A,I6,A,I6,A)') '>>> IMPORTED ', nmics, ' PREVIOUSLY ACCEPTED MICROGRAPHS OF ', ndone,&
            &' MOVIES PROCESSED'
    end subroutine import_previous_projects

    ! Resolves flip_auto and generate from the first movies, then flips the gain reference once
    ! here; parameters and command line are updated together.
    subroutine resolve_gain( self, cline )
        class(stream_stage_preprocess), intent(inout) :: self
        class(cmdline),                 intent(inout) :: cline
        type(string) :: cwd, thumb
        select case( trim(self%flipgain) )
            case( 'flip_auto' )
                self%flipgain = self%detect_gain_flip()
                call cline%set('flipgain', trim(self%flipgain))
            case( 'generate' )
                if( self%l_restart .and. file_exists(GENERATED_GAINREF) )then
                    ! a restart keeps the gain reference generated before
                    self%gainref = simple_abspath(GENERATED_GAINREF)
                    write(logfhandle,'(A,A)') '>>> GAIN GENERATION: reusing ', self%gainref%to_char()
                else
                    self%gainref = self%generate_gain_from_movies()
                endif
                self%flipgain = 'no'
                call cline%set('gainref',  self%gainref)
                call cline%set('flipgain', 'no')
        end select
        if( sigterm_received() ) return ! stopped while waiting for movies: no gain to resolve
        call flip_gain(cline, self%gainref, trim(self%flipgain))
        if( cline%defined('gainref') )then
            if( .not. file_exists(GAIN_THUMBNAIL) )then
                call simple_getcwd(cwd)
                thumb = cwd//'/'//GAIN_THUMBNAIL
                call gainref_to_jpg(self%gainref, thumb)
                write(logfhandle,'(A)') '>>> GAIN REFERENCE'
                write(logfhandle,'(A)') '>>> JPEG '//thumb%to_char()
            endif
        endif
    end subroutine resolve_gain

    ! Flip ('no', 'x', 'y' or 'xy') that best matches the gain reference to batches of the first
    ! movies; 'no' when there is no gain reference or too few movies arrive. Uses its own watcher,
    ! so the movies looked at here are still preprocessed.
    function detect_gain_flip( self ) result( flipgain )
        class(stream_stage_preprocess), intent(inout) :: self
        character(len=:), allocatable :: flipgain
        integer, parameter :: MAX_BATCHES = 8
        integer, parameter :: MAX_WAITS   = 180
        type(stream_watcher)      :: watcher
        type(gain_flip_analyzer)  :: analyzer
        type(image)               :: batch_sum
        type(string), allocatable :: batch(:)
        logical :: l_ok, ran_analysis
        integer :: nwaits, nbatches, part_movies, part_frames
        flipgain = 'no'
        if( self%gainref == '' )then
            write(logfhandle,'(A)') '>>> GAIN AUTO: gainref is empty; defaulting to no flip'
            return
        endif
        if( .not. file_exists(self%gainref) )then
            write(logfhandle,'(A)') '>>> GAIN AUTO: gainref not found; defaulting to no flip'
            return
        endif
        write(logfhandle,'(A,I0,A)') '>>> GAIN AUTO: waiting for movie batches of ', GAIN_BATCH_NMOVIES, ' for orientation analysis'
        call new_movie_watcher(self%params%dir_movies, self%l_sj_dirs, self%settle_s, watcher, l_ok)
        if( .not. l_ok )then
            write(logfhandle,'(A)') '>>> GAIN AUTO: no movie folders detected; defaulting to no flip'
            return
        endif
        call analyzer%new(self%gainref, self%params%smpd)
        nwaits   = 0
        nbatches = 0
        do while( nbatches < MAX_BATCHES )
            if( analyzer%get_converged() ) exit
            if( sigterm_received() ) exit
            call self%next_movie_batch(watcher, batch, l_ok)
            if( .not. l_ok )then
                nwaits = nwaits + 1
                if( nwaits > MAX_WAITS ) exit
                call self%send_status(string('waiting for movies for the gain orientation'))
                call sleep(self%sniff_wait_s)
                cycle
            endif
            call read_movies_and_sum_frames(batch, self%params%smpd, batch_sum, part_movies, part_frames)
            call analyzer%analyze_if_due(batch_sum, part_frames, part_movies, ran_analysis)
            call watcher%add2history(batch)
            call batch_sum%kill()
            deallocate(batch)
            nwaits   = 0
            nbatches = nbatches + 1
        enddo
        flipgain = analyzer%get_flip_mode()
        write(logfhandle,'(A,A)') '>>> GAIN AUTO: selected flip option = ', flipgain
        call analyzer%kill()
        call watcher%kill()
    end function detect_gain_flip

    ! Writes a gain reference estimated from the first movies to the working directory and
    ! returns its absolute path. Uses its own watcher, so these movies are still preprocessed.
    function generate_gain_from_movies( self ) result( gainref )
        class(stream_stage_preprocess), intent(inout) :: self
        type(string) :: gainref
        integer, parameter :: N_MOVIES_TARGET = 1000
        type(stream_watcher)      :: watcher
        type(image)               :: gain_sum
        type(string), allocatable :: batch(:)
        type(string)              :: fname
        logical :: l_ok
        integer :: movies_used, frames_used
        gainref = ''
        write(logfhandle,'(A,I0,A)') '>>> GAIN GENERATION: waiting for ', N_MOVIES_TARGET, ' movies'
        call new_movie_watcher(self%params%dir_movies, self%l_sj_dirs, self%settle_s, watcher, l_ok)
        if( .not. l_ok ) THROW_HARD('GAIN GENERATION: could not detect movie folders')
        movies_used = 0
        frames_used = 0
        ! a session may pause for a grid exchange: the stage waits as long as it takes, and stops
        ! when asked to
        do while( movies_used < N_MOVIES_TARGET )
            if( sigterm_received() )then
                call gain_sum%kill()
                call watcher%kill()
                return
            endif
            call self%next_movie_batch(watcher, batch, l_ok)
            if( .not. l_ok )then
                call self%send_status(string('waiting for movies for the gain reference'))
                call sleep(self%sniff_wait_s)
                cycle
            endif
            call add_movies_to_gain_sum(batch, self%params%smpd, gain_sum, movies_used, frames_used)
            write(logfhandle,'(A,I0,A,I0,A,I0)') '>>> GAIN GENERATION: movies collected ', movies_used, '/', N_MOVIES_TARGET,&
                &'; frames accumulated ', frames_used
            call watcher%add2history(batch)
            deallocate(batch)
        enddo
        fname = GENERATED_GAINREF
        call write_gain_from_sum(gain_sum, frames_used, fname)
        gainref = simple_abspath(GENERATED_GAINREF)
        write(logfhandle,'(A,A)') '>>> GAIN GENERATION: wrote generated gainref ', gainref%to_char()
        call gain_sum%kill()
        call watcher%kill()
    end function generate_gain_from_movies

    ! The next GAIN_BATCH_NMOVIES non-EER movies from watcher, or l_ok=.false. when not enough have arrived.
    subroutine next_movie_batch( self, watcher, batch, l_ok )
        class(stream_stage_preprocess), intent(in)    :: self
        type(stream_watcher),           intent(inout) :: watcher
        type(string), allocatable,      intent(inout) :: batch(:)
        logical,                        intent(out)   :: l_ok
        type(string), allocatable :: movies(:)
        integer :: nmovies, imov, nvalid
        l_ok = .false.
        if( allocated(batch) ) deallocate(batch)
        call watcher%detect_and_add_dirs(self%params%dir_movies, self%l_sj_dirs)
        call watcher%watch(nmovies, movies, max_nmovies=5*GAIN_BATCH_NMOVIES)
        if( nmovies < GAIN_BATCH_NMOVIES ) return
        allocate(batch(GAIN_BATCH_NMOVIES))
        nvalid = 0
        do imov = 1,nmovies
            if( fname2format(movies(imov)) == 'K' ) cycle ! EER
            nvalid        = nvalid + 1
            batch(nvalid) = movies(imov)
            if( nvalid == GAIN_BATCH_NMOVIES ) exit
        enddo
        if( nvalid < GAIN_BATCH_NMOVIES )then
            deallocate(batch)
            return
        endif
        l_ok = .true.
    end subroutine next_movie_batch

    ! The worker command line: the stage's own, run in the job directory without a new
    ! output folder, on one movie set, with the gain reference already resolved.
    subroutine build_worker_cline( self, cline )
        class(stream_stage_preprocess), intent(inout) :: self
        class(cmdline),                 intent(in)    :: cline
        self%cline_exec = cline
        call self%cline_exec%set('prg',   'preprocess')
        call self%cline_exec%set('mkdir', 'no')
        call self%cline_exec%set('dir',   '../')
        call self%cline_exec%set('fromp', 1)
        call self%cline_exec%set('top',   STREAM_NMOVS_SET)
        ! gainref now names the reference as it must be used: workers must not flip it again
        if( self%cline_exec%defined('flipgain') ) call self%cline_exec%set('flipgain', 'no')
    end subroutine build_worker_cline


    !---------------- submission ----------------

    ! Groups newly detected movies into sets of STREAM_NMOVS_SET and queues one preprocess job per set.
    subroutine submit_new_movies( self )
        class(stream_stage_preprocess), intent(inout) :: self
        type(string), allocatable :: movies(:)
        integer :: nmovies, nsets, iset, ifirst, ilast, imovie, cnt
        self%l_movies_left = .false.
        call self%movie_buff%detect_and_add_dirs(self%params%dir_movies, self%l_sj_dirs)
        call self%movie_buff%watch(nmovies, movies, max_nmovies=self%nmovs2importperiter)
        if( nmovies > 0 )then
            ! a movie arrived: preprocessing is not idle
            self%last_movie_time = simple_gettime()
            if( file_exists(STREAM_IDLE_MARKER) ) call del_file(STREAM_IDLE_MARKER)
        endif
        if( nmovies < STREAM_NMOVS_SET ) return
        nsets = nmovies / STREAM_NMOVS_SET
        cnt   = 0
        do iset = 1,nsets
            ifirst = (iset - 1) * STREAM_NMOVS_SET + 1
            ilast  = iset * STREAM_NMOVS_SET
            call self%create_movies_set_project(movies(ifirst:ilast))
            call self%sets%submit(self%qenv, self%cline_exec)
            do imovie = ifirst,ilast
                call self%movie_buff%add2history(movies(imovie))
                cnt = cnt + 1
            enddo
            if( cnt == min(self%nmovs2importperiter, nmovies) ) exit
        enddo
        write(logfhandle,'(A,I4,A,A)') '>>> ', cnt, ' NEW MOVIES ADDED; ', cast_time_char(simple_gettime())
        self%l_movies_left = cnt /= nmovies
    end subroutine submit_new_movies

    ! Writes the project of one movie set (absolute movie paths) as the next job set and points
    ! the worker command line at it.
    subroutine create_movies_set_project( self, movie_names )
        class(stream_stage_preprocess), intent(inout) :: self
        type(string),                   intent(in)    :: movie_names(STREAM_NMOVS_SET)
        type(sp_project) :: spproj_here
        type(ctfparams)  :: ctfvars
        type(string)     :: xmlfile, xmldir
        integer          :: imov
        spproj_here%compenv = self%spproj_glob%compenv
        spproj_here%jobproc = self%spproj_glob%jobproc
        ctfvars%ctfflag = CTFFLAG_YES
        ctfvars%smpd    = self%params%smpd
        ctfvars%cs      = self%params%cs
        ctfvars%kv      = self%params%kv
        ctfvars%fraca   = self%params%fraca
        call spproj_here%add_movies(movie_names(1:STREAM_NMOVS_SET), ctfvars, verbose=.false.)
        do imov = 1,STREAM_NMOVS_SET
            self%import_counter = self%import_counter + 1
            call spproj_here%os_mic%set(imov, 'importind', real(self%import_counter))
            call spproj_here%os_mic%set(imov, 'tiltgrp',   0.0)
            call spproj_here%os_mic%set(imov, 'shiftx',    0.0)
            call spproj_here%os_mic%set(imov, 'shifty',    0.0)
            call spproj_here%os_mic%set(imov, 'flsht',     0.0)
            if( self%l_xml_meta )then
                if( self%l_sj_dirs )then
                    xmldir = stemname(movie_names(imov))
                else
                    xmldir = self%params%dir_meta
                endif
                xmlfile = basename(movie_names(imov))
                if( xmlfile%substr_ind('_fractions') > 0 ) xmlfile = xmlfile%to_char([1, xmlfile%substr_ind('_fractions') - 1])
                if( xmlfile%substr_ind('_EER')       > 0 ) xmlfile = xmlfile%to_char([1, xmlfile%substr_ind('_EER')       - 1])
                xmlfile = xmldir//'/'//xmlfile//'.xml'
                call spproj_here%os_mic%set(imov, 'meta', xmlfile)
            endif
        enddo
        call self%sets%write_set(spproj_here, self%cline_exec, STREAM_NMOVS_SET)
        call spproj_here%kill
    end subroutine create_movies_set_project

    subroutine schedule_jobs( self )
        class(stream_stage_preprocess), intent(inout) :: self
        integer :: stacksz
        call self%sets%schedule(self%qenv)
        stacksz = self%qenv%qscripts%get_stacksz()
        if( stacksz /= self%prev_stacksz )then
            self%prev_stacksz = stacksz
            write(logfhandle,'(A,I6)') '>>> MOVIES TO PROCESS:                ', stacksz * STREAM_NMOVS_SET
        endif
    end subroutine schedule_jobs

    !---------------- import ----------------

    ! Imports finished jobs and counts failed ones; @p n_imported is the number of micrographs imported now.
    subroutine collect_jobs( self, n_imported )
        class(stream_stage_preprocess), intent(inout) :: self
        integer,                        intent(out)   :: n_imported
        type(string), allocatable :: done(:)
        integer :: n_failed
        n_imported = 0
        call self%sets%collect(self%qenv, done, n_failed)
        if( size(done) > 0 ) call self%import_completed(done, n_imported)
        self%n_failed_jobs = self%n_failed_jobs + n_failed
    end subroutine collect_jobs

    ! Appends the accepted micrographs of the finished sets @p job_fnames (absolute paths) to the
    ! global project, applies the thresholds to each set's own project, and moves the sets with
    ! an accepted micrograph to the completed folder.
    subroutine import_completed( self, job_fnames, n_imported )
        class(stream_stage_preprocess), intent(inout) :: self
        type(string),                   intent(in)    :: job_fnames(:)
        integer,                        intent(out)   :: n_imported
        type(sp_project), allocatable :: job_projs(:)
        logical,          allocatable :: mics_mask(:)
        type(string) :: fname
        integer      :: n_jobs, n_old, nmics, iproj, i, imic, j, nrejected
        n_imported = 0
        n_jobs     = size(job_fnames)
        if( n_jobs == 0 ) return
        n_old = self%spproj_glob%os_mic%get_noris()
        nmics = STREAM_NMOVS_SET * n_jobs
        allocate(job_projs(n_jobs), mics_mask(nmics))
        do iproj = 1,n_jobs
            call job_projs(iproj)%read_segment('mic', job_fnames(iproj))
            do i = 1,STREAM_NMOVS_SET
                mics_mask((iproj - 1) * STREAM_NMOVS_SET + i) = job_projs(iproj)%os_mic%get_state(i) == 1
            enddo
        enddo
        n_imported         = count(mics_mask)
        self%n_failed_jobs = self%n_failed_jobs + (nmics - n_imported)
        if( n_imported > 0 )then
            if( n_old == 0 )then
                call self%spproj_glob%os_mic%new(n_imported, is_ptcl=.false.)
            else
                call self%spproj_glob%os_mic%reallocate(n_old + n_imported)
            endif
            imic = 0
            j    = n_old
            do iproj = 1,n_jobs
                do i = 1,STREAM_NMOVS_SET
                    imic = imic + 1
                    if( .not. mics_mask(imic) ) cycle
                    j = j + 1
                    call self%spproj_glob%os_mic%transfer_ori(j, job_projs(iproj)%os_mic, i)
                enddo
                call self%apply_thresholds(job_projs(iproj)%os_mic, nrejected)
                if( nrejected > 0 ) call job_projs(iproj)%write_segment_inside('mic', job_fnames(iproj))
            enddo
        endif
        do iproj = 1,n_jobs
            imic = (iproj - 1) * STREAM_NMOVS_SET + 1
            if( any(mics_mask(imic:imic+STREAM_NMOVS_SET-1)) ) call self%sets%complete(job_fnames(iproj), fname)
            call job_projs(iproj)%kill
        enddo
        deallocate(job_projs, mics_mask)
    end subroutine import_completed

    ! After an import: thresholds, log, GUI plots and thumbnails, STAR snapshot.
    subroutine process_imports( self )
        class(stream_stage_preprocess), intent(inout) :: self
        integer :: nmics, nrejected
        nmics = self%spproj_glob%os_mic%get_noris()
        call self%apply_thresholds(self%spproj_glob%os_mic, nrejected)
        write(logfhandle,'(A,I8)')       '>>> # MOVIES PROCESSED & IMPORTED       : ', nmics
        write(logfhandle,'(A,I3,A2,I3)') '>>> # OF COMPUTING UNITS IN USE/TOTAL   : ', self%qenv%get_navail_computing_units(),&
                                         &'/ ', self%params%nparts
        if( self%n_failed_jobs > 0 ) write(logfhandle,'(A,I8)') '>>> # DESELECTED MICROGRAPHS/FAILED JOBS: ', self%n_failed_jobs
        ! plots follow the first status message, as the GUI expects
        if( self%meta_status%assigned() ) call self%send_plots()
        self%last_injection = simple_gettime()
        self%l_haschanged   = .true.
        if( nmics < STAR_EVERY_NMICS )then
            call self%write_mic_star_and_field()
        else if( nmics > self%nmic_star + STAR_STEP_NMICS )then
            call self%write_mic_star_and_field()
            self%nmic_star = nmics
        endif
    end subroutine process_imports

    subroutine send_plots( self )
        class(stream_stage_preprocess), intent(inout) :: self
        integer :: nmics, i_max, iori, imic
        associate( os_mic => self%spproj_glob%os_mic )
            if( os_mic%isthere('ctfres') )then
                call set_histogram_from_oris(self%meta_hist_ctfres, os_mic, 'ctfres', CTFRES_BINS)
                call self%pipe%send_meta(self%meta_hist_ctfres)
                call set_timeplot_from_oris(self%meta_plot_ctfres, os_mic, 'ctfres', PLOT_WINDOW)
                call self%pipe%send_meta(self%meta_plot_ctfres)
            endif
            if( os_mic%isthere('icefrac') )then
                call set_histogram_from_oris(self%meta_hist_icefrac, os_mic, 'icefrac', ICESCORE_BINS)
                call self%pipe%send_meta(self%meta_hist_icefrac)
            endif
            if( os_mic%isthere('astig') )then
                call set_histogram_from_oris(self%meta_hist_astig, os_mic, 'astig', ASTIG_BINS)
                call self%pipe%send_meta(self%meta_hist_astig)
                call set_timeplot_from_oris(self%meta_plot_astig, os_mic, 'astig', PLOT_WINDOW)
                call self%pipe%send_meta(self%meta_plot_astig)
            endif
            call set_timeplot_from_oris(self%meta_plot_df, os_mic, 'df', PLOT_WINDOW)
            call self%pipe%send_meta(self%meta_plot_df)
            if( allocated(self%movie_buff%ratehistory) )then
                call set_rate_timeplot(self%meta_plot_rate, self%movie_buff%ratehistory)
                call self%pipe%send_meta(self%meta_plot_rate)
            endif
            if( os_mic%isthere('thumb') )then
                nmics = os_mic%get_noris()
                i_max = min(nmics, NTHUMBS)
                do iori = 1,i_max
                    imic = nmics - i_max + iori
                    call self%meta_micrograph%set(path=os_mic%get_str(imic, 'thumb'), dfx=os_mic%get(imic, 'dfx'),&
                        &dfy=os_mic%get(imic, 'dfy'), ctfres=os_mic%get(imic, 'ctfres'), i_max=i_max, i=iori)
                    call self%pipe%send_meta(self%meta_micrograph)
                enddo
            endif
        end associate
    end subroutine send_plots

    ! STREAM_IDLE in the stage's folder once no new movie has been seen for MOVIES_IDLE_TIME_S and
    ! every submitted set has been collected (none queued or running); removed as soon as a movie
    ! arrives (submit_new_movies). Downstream stages end their intake on it.
    subroutine update_idle_marker( self )
        class(stream_stage_preprocess), intent(inout) :: self
        logical :: l_drained
        if( file_exists(STREAM_IDLE_MARKER) ) return
        if( simple_gettime() - self%last_movie_time < MOVIES_IDLE_TIME_S ) return
        l_drained = self%qenv%qscripts%get_stacksz() == 0
        if( l_drained ) l_drained = self%qenv%get_navail_computing_units() >= self%params%nparts
        if( .not. l_drained ) return
        call simple_touch(STREAM_IDLE_MARKER)
        write(logfhandle,'(A,I0,A)') '>>> NO NEW MOVIE FOR ', MOVIES_IDLE_TIME_S, ' S AND NOTHING LEFT TO PROCESS: IDLE'
    end subroutine update_idle_marker

    ! Nothing imported and nothing pending: snapshot after a long inactivity, otherwise wait.
    subroutine idle( self )
        class(stream_stage_preprocess), intent(inout) :: self
        if( (simple_gettime() - self%last_injection > INACTIVE_TIME) .and. self%l_haschanged )then
            call self%write_mic_star_and_field()
            self%l_haschanged = .false.
        else
            call sleep(self%wait_s)
        endif
    end subroutine idle

    !---------------- GUI ----------------

    subroutine send_status( self, stage )
        class(stream_stage_preprocess), intent(inout) :: self
        type(string), optional,         intent(in)    :: stage
        type(string) :: stage_here
        real    :: average_ctfres, average_astig, average_icefrac
        integer :: nmics
        stage_here = 'finding and processing new movies'
        if( present(stage) ) stage_here = stage
        associate( os_mic => self%spproj_glob%os_mic )
            nmics           = os_mic%get_noris()
            average_ctfres  = 0.
            average_astig   = 0.
            average_icefrac = 0.
            if( os_mic%isthere('ctfres')  ) average_ctfres  = os_mic%get_avg('ctfres')
            if( os_mic%isthere('icefrac') ) average_icefrac = os_mic%get_avg('icefrac')
            if( os_mic%isthere('astig')   ) average_astig   = os_mic%get_avg('astig')
            call self%meta_status%set(stage=stage_here,                                         &
                movies_imported     = self%movie_buff%n_history,                                &
                movies_processed    = nmics + self%n_failed_jobs,                               &
                movies_rejected     = self%n_failed_jobs + nmics - os_mic%count_state_gt_zero(),&
                movies_rate         = self%movie_buff%rate,                                     &
                average_ctf_res     = average_ctfres,                                           &
                average_ice_score   = average_icefrac,                                          &
                average_astigmatism = average_astig,                                            &
                cutoff_ctf_res      = self%ctfres_thres,                                        &
                cutoff_ice_score    = self%icefrac_thres,                                       &
                cutoff_astigmatism  = self%astig_thres)
        end associate
        call self%pipe%send_meta(self%meta_status)
    end subroutine send_status

    ! Drains every queued GUI update; threshold changes apply to later imports and to the global project.
    subroutine apply_gui_updates( self )
        class(stream_stage_preprocess), intent(inout) :: self
        type(gui_metadata_stream_update) :: update
        character(len=:), allocatable    :: buffer
        real :: val
        do while( self%pipe%receive(buffer) )
            update = transfer(buffer, update)
            val    = update%get_ctfres_update()
            if( abs(val) > 0.001 .and. val /= self%ctfres_thres  ) call self%set_threshold('ctfresthreshold',  val)
            val    = update%get_astigmatism_update()
            if( abs(val) > 0.001 .and. val /= self%astig_thres   ) call self%set_threshold('astigthreshold',   val)
            val    = update%get_icescore_update()
            if( abs(val) > 0.001 .and. val /= self%icefrac_thres ) call self%set_threshold('icefracthreshold', val)
        enddo
    end subroutine apply_gui_updates

    ! The single place a rejection threshold changes after new(): the stage's threshold and the
    ! worker command line together.
    subroutine set_threshold( self, key, val )
        class(stream_stage_preprocess), intent(inout) :: self
        character(len=*),               intent(in)    :: key
        real,                           intent(in)    :: val
        select case( key )
            case( 'ctfresthreshold' )
                self%ctfres_thres  = val
                write(logfhandle,'(A,F8.2)') '>>> CTF RESOLUTION THRESHOLD UPDATED TO: ', val
            case( 'astigthreshold' )
                self%astig_thres   = val
                write(logfhandle,'(A,F8.2)') '>>> ASTIGMATISM THRESHOLD UPDATED TO: ', val
            case( 'icefracthreshold' )
                self%icefrac_thres = val
                write(logfhandle,'(A,F8.2)') '>>> ICE SCORE THRESHOLD UPDATED TO: ', val
            case DEFAULT
                THROW_HARD('not a micrograph rejection threshold: '//key)
        end select
        call self%cline_exec%set(key, val)
        call self%save_thresholds()
    end subroutine set_threshold

    !> The thresholds, kept for a restart (written as a temporary file and renamed).
    subroutine save_thresholds( self )
        class(stream_stage_preprocess), intent(inout) :: self
        integer :: funit, ios
        open(newunit=funit, file=GUI_THRESHOLDS//'.tmp', status='replace', action='write', iostat=ios)
        if( ios /= 0 ) return
        write(funit,*) self%ctfres_thres, self%astig_thres, self%icefrac_thres
        close(funit)
        call simple_rename(GUI_THRESHOLDS//'.tmp', GUI_THRESHOLDS)
    end subroutine save_thresholds

    !> Restart: the thresholds the GUI set before hold from the first pass, before any answer of
    !! the GUI reaches the restarted stage.
    subroutine restore_thresholds( self )
        class(stream_stage_preprocess), intent(inout) :: self
        real    :: ctfres, astig, icefrac
        integer :: funit, ios
        if( .not. file_exists(GUI_THRESHOLDS) ) return
        open(newunit=funit, file=GUI_THRESHOLDS, status='old', action='read', iostat=ios)
        if( ios /= 0 ) return
        read(funit,*,iostat=ios) ctfres, astig, icefrac
        close(funit)
        if( ios /= 0 ) return
        if( ctfres  /= self%ctfres_thres  ) call self%set_threshold('ctfresthreshold',  ctfres)
        if( astig   /= self%astig_thres   ) call self%set_threshold('astigthreshold',   astig)
        if( icefrac /= self%icefrac_thres ) call self%set_threshold('icefracthreshold', icefrac)
    end subroutine restore_thresholds

    !---------------- helpers ----------------

    subroutine write_mic_star_and_field( self, write_field, copy_optics )
        class(stream_stage_preprocess), intent(inout) :: self
        logical, optional,              intent(in)    :: write_field, copy_optics
        logical :: l_write_field, l_copy_optics
        l_write_field = .false.
        l_copy_optics = .false.
        if( present(write_field) ) l_write_field = write_field
        if( present(copy_optics) ) l_copy_optics = copy_optics
        if( l_copy_optics )then
            call self%starproj_stream%copy_micrographs_optics(self%spproj_glob, verbose=.false.)
            call self%starproj_stream%stream_export_micrographs(self%params, self%spproj_glob, self%params%cwd, optics_set=.true.)
        else
            call self%starproj_stream%stream_export_micrographs(self%params, self%spproj_glob, self%params%cwd)
        endif
        if( l_write_field )then
            call self%spproj_glob%write_segment_inside('mic', self%params%projfile)
            call self%spproj_glob%write_non_data_segments(self%params%projfile)
        endif
    end subroutine write_mic_star_and_field

    ! A watcher on the movie directory, or on its first sub-folder in the multi-folder layout;
    ! l_ok is .false. when that layout has no sub-folder yet. Shared by the main loop's watcher
    ! and the gain step's own one; a free procedure so the stage can pass its own movie_buff.
    subroutine new_movie_watcher( dir_movies, l_sj_dirs, settle_s, watcher, l_ok )
        type(string),         intent(in)    :: dir_movies
        logical,              intent(in)    :: l_sj_dirs
        integer,              intent(in)    :: settle_s
        type(stream_watcher), intent(inout) :: watcher
        logical,              intent(out)   :: l_ok
        type(string), allocatable :: dirs(:)
        l_ok = .true.
        if( l_sj_dirs )then
            call sniff_folders_SJ(dir_movies, l_ok, dirs)
            if( .not. l_ok ) return
            watcher = stream_watcher(settle_s, dirs(1), suffix_filter=string('_fractions'))
            call watcher%detect_and_add_dirs(dir_movies, l_sj_dirs)
            deallocate(dirs)
        else
            watcher = stream_watcher(settle_s, dir_movies)
        endif
    end subroutine new_movie_watcher

    ! The current rejection thresholds applied to a micrograph segment.
    subroutine apply_thresholds( self, os_mic, nrejected )
        class(stream_stage_preprocess), intent(in)    :: self
        class(oris),                    intent(inout) :: os_mic
        integer,                        intent(out)   :: nrejected
        call reject_mics_by_thresholds(os_mic, nrejected, ctfres=self%ctfres_thres,&
            &icefrac=self%icefrac_thres, astig=self%astig_thres)
    end subroutine apply_thresholds

end module simple_stream_stage_preprocess
