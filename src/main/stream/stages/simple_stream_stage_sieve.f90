!@descr: state and steps of stream task 5 (particle sieving): import extracted particles, run the particle sieve, report to the GUI
!==============================================================================
! MODULE: simple_stream_stage_sieve
!
! PURPOSE:
!   The body of stream p05 as a type. Each pass imports the particles of the
!   reference-picking sets completed since the last pass and drives one cycle
!   of the particle sieve (collect and reject, coarse and fine chunks,
!   submission). The sieve hands its finished chunks to the completed folder,
!   with the groups of the newest optics map, where pool 2D picks them up. The
!   commander (simple_commanders_stream_p05_sieve_cavgs) only normalises the command line
!   and loops over iterate() until finished().
!
!   The mask diameter is the one make_pickrefs decided, read from moldiam.txt
!   in the reference-picking stage's directory (dir_target).
!
!   What happens is delegated:
!     - chunking, 2D, rejection, hand-off -> ptcl_sieve
!     - project import                    -> import_new_projects
!     - GUI                               -> simple_stream_pipe, simple_stream_gui_senders
!
! LIFECYCLE:
!   new(cline) -> { iterate() } until finished() -> finalize() -> kill()
!
! RESTART:
!   Recognised by the output directory and logged; a leftover TERM_STREAM is
!   removed, so the restarted stage runs. The sieve restores its chunks from
!   its folders. The projects it has already chunked (its imported_projects.txt)
!   go into the watcher history; projects imported but not yet chunked are
!   imported again.
!==============================================================================
module simple_stream_stage_sieve
use simple_defs,                                 only: logfhandle, PATH_HERE
use simple_defs_fname,                           only: TERM_STREAM, STREAM_MOLDIAM
use simple_defs_stream,                          only: DIR_STREAM_COMPLETED, SHORTWAIT, WAITTIME
use simple_defs_environment,                     only: SIMPLE_STREAM_CHUNK_PARTITION
use simple_error,                                only: simple_exception
use simple_string,                               only: string
use simple_fileio,                               only: del_file, file_exists, read_filetable, simple_abspath
use simple_syslib,                               only: dir_exists, simple_mkdir
use simple_timer,                                only: simple_gettime, cast_time_char
use simple_cmdline,                              only: cmdline
use simple_oris,                                 only: oris
use simple_parameters,                           only: parameters
use simple_sp_project,                           only: sp_project
use simple_qsys_env,                             only: qsys_env
use simple_rec_list,                             only: rec_list
use simple_stream_watcher,                       only: stream_watcher
use simple_stream_state,                         only: ipc_pipe_sieve_cavgs_in
use simple_stream_utils,                         only: create_stream_project, init_stream_qenv, import_new_projects
use simple_ptcl_sieve,                           only: ptcl_sieve
use simple_gui_metadata_utils,                   only: max_metadata_size
use simple_gui_metadata_types,                   only: GUI_METADATA_STREAM_PARTICLE_SIEVING_TYPE,&
                                                      &GUI_METADATA_STREAM_PARTICLE_SIEVING_CLS2D_TYPE
use simple_gui_metadata_cavg2D,                  only: gui_metadata_cavg2D
use simple_gui_metadata_stream_particle_sieving, only: gui_metadata_stream_particle_sieving
use simple_stream_pipe,                          only: stream_pipe
use simple_stream_gui_senders,                   only: send_cavgs
implicit none

public :: stream_stage_sieve
private
#include "simple_local_flags.inc"

integer,          parameter :: MAX_PROJECTS_IMPORT       = 20      ! completed upstream sets taken per pass
integer,          parameter :: FINAL_INGESTION_IDLE_TIME = 10 * 60 ! idle time (s) after the last import before final ingestion
character(len=*), parameter :: IMPORTED_PROJECTS         = 'imported_projects.txt' ! written by the sieve

! Components and steps are public so simple_stream_stage_sieve_tester can assemble a stage and
! run one step at a time; production code uses new/iterate/finished/finalize/kill.
type :: stream_stage_sieve
    type(parameters)                           :: params
    type(qsys_env)                             :: qenv          ! starts the persistent workers the sieve's chunk jobs run on
    type(sp_project)                           :: spproj        ! the stage's project
    type(stream_watcher)                       :: project_buff  ! completed reference-picking sets
    type(rec_list)                             :: project_list  ! one record per imported micrograph
    type(ptcl_sieve), allocatable              :: sieve
    type(stream_pipe)                          :: pipe          ! to the master
    type(gui_metadata_stream_particle_sieving) :: meta_status
    type(gui_metadata_cavg2D)                  :: meta_cavgs
    type(string), allocatable :: restored_imports(:)            ! restart: projects the sieve has already chunked
    integer,      allocatable :: latest_inds(:), latest_pops(:), latest_selection(:)
    real,         allocatable :: latest_res(:)
    type(string)              :: latest_jpeg, latest_stk        ! the sieve's latest class averages
    integer :: latest_xtiles    = 0
    integer :: latest_ytiles    = 0
    integer :: n_mics_imported  = 0
    integer :: n_ptcls_imported = 0
    integer :: last_import_time = 0
    logical :: l_attached       = .false. ! the upstream completed-sets folder exists and is watched
    logical :: l_waiting_logged = .false.
    logical :: l_sieve_active   = .false. ! the sieve is made (on the first import)
    logical :: l_restart        = .false. ! the output directory existed before params%new
    logical :: l_exists         = .false.
    ! waits (s); tests set them to 0, and settle_s to -1 to take files written in the same second
    integer :: settle_s         = SHORTWAIT ! a completed set is taken once untouched longer than this;
                                            ! reference picking moves finished sets in with a rename
    integer :: wait_s           = WAITTIME  ! pause at the end of a pass
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
    procedure :: restore_imports
    ! the steps of iterate() and their helpers
    procedure :: attach_upstream
    procedure :: import_projects
    procedure :: start_sieve
    procedure :: read_mask_diameter
    procedure :: send_status
    procedure :: send_latest_cavgs
end type stream_stage_sieve

contains

    !---------------- lifecycle ----------------

    !> Builds the parameters from the normalised @p cline, starts the chunk workers, and prepares
    !! the GUI and a restart's history.
    subroutine new( self, cline )
        class(stream_stage_sieve), intent(inout) :: self
        class(cmdline),            intent(inout) :: cline
        call self%kill()
        call self%init_params(cline)
        self%l_exists = .true. ! from here kill() releases what has been built
        call self%init_queue()
        call self%init_gui(-1, ipc_pipe_sieve_cavgs_in(2))
        call self%restore_imports()
        self%last_import_time = simple_gettime()
        call self%send_status(string('initialising'))
    end subroutine new

    !> The stage's project file (made and given a computing environment) and its parameters; the
    !! project must start without micrographs. A restart's leftover TERM_STREAM is removed in the
    !! stage's folder, where the loop looks for it.
    subroutine init_params( self, cline )
        class(stream_stage_sieve), intent(inout) :: self
        class(cmdline),            intent(inout) :: cline
        type(string) :: outdir
        ! a restart is recognised by its output directory, before params%new makes one
        self%l_restart = .false.
        if( cline%defined('outdir') )then
            outdir = cline%get_carg('outdir')
            if( outdir%strlen() > 0 ) self%l_restart = dir_exists(outdir)
        endif
        call create_stream_project(self%spproj, cline, string('sieve_cavgs'))
        call self%params%new(cline)
        if( self%l_restart )then
            write(logfhandle,'(A)') '>>> RESTARTING EXISTING JOB'
            call del_file(TERM_STREAM)
        endif
        call self%spproj%read(self%params%projfile)
        if( self%spproj%os_mic%get_noris() /= 0 )then
            THROW_HARD('commander_stream_p05_sieve_cavgs must start from an empty project (e.g. from root project folder)')
        endif
        call simple_mkdir(PATH_HERE//DIR_STREAM_COMPLETED)
    end subroutine init_params

    !> The queue environment on the chunk partition. The sieve submits through its own; this one
    !! starts the persistent workers (workers = nchunks) that the sieve's then reuses.
    subroutine init_queue( self )
        class(stream_stage_sieve), intent(inout) :: self
        call init_stream_qenv(self%params, self%qenv, string(SIMPLE_STREAM_CHUNK_PARTITION))
    end subroutine init_queue

    !> The GUI metadata objects and the pipe ends to the master (-1: none).
    subroutine init_gui( self, fd_read, fd_write )
        class(stream_stage_sieve), intent(inout) :: self
        integer,                   intent(in)    :: fd_read, fd_write
        call self%meta_status%new(GUI_METADATA_STREAM_PARTICLE_SIEVING_TYPE)
        call self%meta_cavgs%new(GUI_METADATA_STREAM_PARTICLE_SIEVING_CLS2D_TYPE)
        call self%pipe%new(fd_read, fd_write, max_metadata_size(), 'particle_sieving')
    end subroutine init_gui

    !> Restart: the projects the sieve has already chunked, for the watcher history (attach_upstream).
    subroutine restore_imports( self )
        class(stream_stage_sieve), intent(inout) :: self
        if( .not. file_exists(IMPORTED_PROJECTS) ) return
        call self%send_status(string('importing previous run'))
        call read_filetable(string(IMPORTED_PROJECTS), self%restored_imports)
    end subroutine restore_imports

    !> One pass: wait for the upstream folder, import new sets, run the sieve, report.
    subroutine iterate( self )
        class(stream_stage_sieve), intent(inout) :: self
        if( .not. self%l_attached )then
            call self%attach_upstream()
            if( .not. self%l_attached )then
                call self%send_status(string('waiting on reference picking'))
                call sleep(self%wait_s)
                return
            endif
        endif
        call self%import_projects()
        if( self%l_sieve_active )then
            ! final ingestion once no set has arrived for a while; the next import undoes it
            if( simple_gettime() - self%last_import_time >= FINAL_INGESTION_IDLE_TIME ) call self%sieve%set_final_ingestion()
            call self%sieve%cycle(self%project_list)
        else if( self%project_list%size() > 0 )then
            call self%start_sieve()
        endif
        if( self%n_ptcls_imported > 0 )then
            call self%send_status(string('importing and sieving particles'))
        else
            call self%send_status(string('waiting on reference picking'))
        endif
        call self%send_latest_cavgs()
        call sleep(self%wait_s)
    end subroutine iterate

    !> .true. once the stream is told to stop.
    logical function finished( self )
        class(stream_stage_sieve), intent(in) :: self
        finished = file_exists(TERM_STREAM)
    end function finished

    !> The last status: no more user input.
    subroutine finalize( self )
        class(stream_stage_sieve), intent(inout) :: self
        call self%meta_status%set_user_input(.false.)
        call self%send_status(string('terminating'))
    end subroutine finalize

    subroutine kill( self )
        class(stream_stage_sieve), intent(inout) :: self
        if( .not. self%l_exists ) return
        if( allocated(self%sieve) )then
            if( self%l_sieve_active ) call self%sieve%kill
            deallocate(self%sieve)
        endif
        call self%project_buff%kill
        call self%project_list%kill
        call self%qenv%kill
        call self%spproj%kill
        call self%pipe%kill
        call self%meta_status%kill
        call self%meta_cavgs%kill
        if( allocated(self%restored_imports) ) deallocate(self%restored_imports)
        if( allocated(self%latest_inds)      ) deallocate(self%latest_inds)
        if( allocated(self%latest_pops)      ) deallocate(self%latest_pops)
        if( allocated(self%latest_selection) ) deallocate(self%latest_selection)
        if( allocated(self%latest_res)       ) deallocate(self%latest_res)
        call self%latest_jpeg%kill
        call self%latest_stk%kill
        self%latest_xtiles    = 0
        self%latest_ytiles    = 0
        self%n_mics_imported  = 0
        self%n_ptcls_imported = 0
        self%l_attached       = .false.
        self%l_waiting_logged = .false.
        self%l_sieve_active   = .false.
        self%l_restart        = .false.
        self%l_exists         = .false.
    end subroutine kill

    !---------------- steps ----------------

    ! Starts watching the upstream completed-sets folder once reference picking has created it;
    ! a restart's already chunked projects go into the watcher history.
    subroutine attach_upstream( self )
        class(stream_stage_sieve), intent(inout) :: self
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
        if( allocated(self%restored_imports) )then
            do i = 1,size(self%restored_imports)
                call self%project_buff%add2history(self%restored_imports(i))
            enddo
            deallocate(self%restored_imports)
        endif
        self%l_attached       = .true.
        self%l_waiting_logged = .false.
    end subroutine attach_upstream

    ! One record per micrograph of the newly completed upstream sets.
    subroutine import_projects( self )
        class(stream_stage_sieve), intent(inout) :: self
        type(string), allocatable :: projects(:)
        integer :: nprojects
        call self%project_buff%watch(nprojects, projects, max_nmovies=MAX_PROJECTS_IMPORT)
        if( nprojects == 0 ) return
        call import_new_projects(self%project_list, projects, self%n_mics_imported, self%n_ptcls_imported)
        call self%project_buff%add2history(projects)
        if( self%l_sieve_active ) call self%sieve%unset_final_ingestion()
        self%last_import_time = simple_gettime()
        write(logfhandle,'(A,I6,I9)') '>>> # MICROGRAPHS / PARTICLES IMPORTED : ', self%n_mics_imported, self%n_ptcls_imported
        write(logfhandle,'(A,A)')     '>>> LAST IMPORT AT                     : ', cast_time_char(self%last_import_time)
    end subroutine import_projects

    ! The sieve, made on the first import with the mask diameter of the picking references and the
    ! optics directory for its hand-offs, then two warm-up cycles.
    subroutine start_sieve( self )
        class(stream_stage_sieve), intent(inout) :: self
        call self%read_mask_diameter()
        if( .not. allocated(self%sieve) ) allocate(self%sieve)
        if( self%params%optics_dir%strlen() > 0 )then
            call self%sieve%new(self%params, string(PATH_HERE//DIR_STREAM_COMPLETED), optics_dir=self%params%optics_dir)
        else
            call self%sieve%new(self%params, string(PATH_HERE//DIR_STREAM_COMPLETED))
        endif
        self%l_sieve_active = .true.
        call self%sieve%cycle(self%project_list)
        call self%sieve%cycle(self%project_list)
    end subroutine start_sieve

    ! The mask diameter make_pickrefs decided, from moldiam.txt in the reference-picking
    ! stage's directory; it is written before reference picking completes any set.
    subroutine read_mask_diameter( self )
        class(stream_stage_sieve), intent(inout) :: self
        type(oris)   :: moldiam
        type(string) :: fname
        fname = self%params%dir_target//'/'//STREAM_MOLDIAM
        if( .not. file_exists(fname) ) THROW_HARD('no mask diameter from reference picking: '//fname%to_char())
        call moldiam%new(1, is_ptcl=.false.)
        call moldiam%read(fname)
        self%params%mskdiam = moldiam%get(1, 'mskdiam')
        call moldiam%kill
        write(logfhandle,'(A,F8.2)') '>>> MASK DIAMETER SET TO : ', self%params%mskdiam
    end subroutine read_mask_diameter

    !---------------- GUI ----------------

    ! Particle counts and the classes of the latest product the sieve selected.
    subroutine send_status( self, stage )
        class(stream_stage_sieve), intent(inout) :: self
        type(string),              intent(in)    :: stage
        integer :: i, naccepted, nrejected
        naccepted = 0
        nrejected = 0
        if( allocated(self%sieve) )then
            naccepted = self%sieve%get_n_accepted_ptcls()
            nrejected = self%sieve%get_n_rejected_ptcls()
        endif
        call self%meta_status%set(stage=stage, particles_imported=self%n_ptcls_imported,&
            &particles_accepted=naccepted, particles_rejected=nrejected)
        call self%meta_status%clear_selection()
        if( allocated(self%latest_inds) .and. allocated(self%latest_selection) )then
            do i = 1,size(self%latest_inds)
                if( self%latest_selection(i) /= 0 ) call self%meta_status%set_selection(self%latest_inds(i))
            enddo
        endif
        call self%pipe%send_meta(self%meta_status)
    end subroutine send_status

    ! The sieve's latest class averages as sprite-sheet tiles, with their resolutions and populations.
    subroutine send_latest_cavgs( self )
        class(stream_stage_sieve), intent(inout) :: self
        if( .not. self%l_sieve_active ) return
        if( .not. self%sieve%get_latest(self%latest_inds, self%latest_pops, self%latest_res, self%latest_jpeg,&
            &self%latest_stk, self%latest_xtiles, self%latest_ytiles, self%latest_selection) ) return
        call send_cavgs(self%pipe, self%meta_cavgs, self%latest_jpeg, self%latest_inds, self%latest_stk,&
            &self%latest_xtiles, self%latest_ytiles, res=self%latest_res, pop=self%latest_pops)
    end subroutine send_latest_cavgs

end module simple_stream_stage_sieve
