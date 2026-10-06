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
!   Final ingestion (the sieve stages its leftover particles and ends with a
!   final set) is set while reference picking is idle or stopped (its
!   STREAM_IDLE or STREAM_FINISHED marker in dir_target), once a watch a settle
!   time later has found nothing new, and withdrawn when the marker goes or a
!   new set arrives.
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
!   its folders, and is made as soon as the upstream folder is attached, even
!   with nothing new to import. Every set it has chunked from (its
!   chunked_mics.txt) is imported again with the micrographs it chunked marked,
!   so the rest of a partly chunked set is still sieved, and goes into the
!   watcher history; sets imported but not chunked from are imported again.
!   On stop, the 2D jobs of the running chunks are cancelled.
!==============================================================================
module simple_stream_stage_sieve
use simple_defs,                                 only: logfhandle, PATH_HERE, STDLEN
use simple_defs_fname,                           only: TERM_STREAM, STREAM_MOLDIAM
use simple_defs_stream,                          only: DIR_STREAM_COMPLETED, SHORTWAIT, WAITTIME
use simple_defs_environment,                     only: SIMPLE_STREAM_CHUNK_PARTITION
use simple_error,                                only: simple_exception
use simple_string,                               only: string
use simple_string_utils,                         only: int2str
use simple_fileio,                               only: del_file, file_exists, simple_abspath
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
use simple_stream_utils,                         only: create_stream_project, init_stream_qenv, import_new_projects, upstream_done
use simple_ptcl_sieve,                           only: ptcl_sieve, ptcl_sieve_settings, sieve_settings, CHUNKED_MICS,&
                                                      &read_chunked_mics, rebuild_chunked_mics
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

! Components and steps are public so simple_stream_stage_sieve_tester can assemble a stage and
! run one step at a time; production code uses new/iterate/finished/finalize/kill.
type :: stream_stage_sieve
    type(parameters), allocatable              :: params
    ! allocatable (compile-time policy), allocated in init_params and released in kill
    type(qsys_env),   allocatable              :: qenv          ! starts the persistent workers the sieve's chunk jobs run on
    type(sp_project), allocatable              :: spproj        ! the stage's project
    type(stream_watcher)                       :: project_buff  ! completed reference-picking sets
    type(rec_list)                             :: project_list  ! one record per imported micrograph
    type(ptcl_sieve), allocatable              :: sieve
    type(stream_pipe)                          :: pipe          ! to the master
    type(gui_metadata_stream_particle_sieving) :: meta_status
    type(gui_metadata_cavg2D)                  :: meta_cavgs
    type(string), allocatable :: restored_imports(:)            ! restart: the sets the sieve has chunked from
    integer,      allocatable :: latest_inds(:), latest_pops(:), latest_selection(:)
    real,         allocatable :: latest_res(:)
    type(string)              :: latest_jpeg, latest_stk        ! the sieve's latest class averages
    real    :: mskdiam          = 0.      ! the picking references' mask diameter (moldiam.txt), the sieve's
    integer :: latest_xtiles    = 0
    integer :: latest_ytiles    = 0
    integer :: n_mics_imported  = 0
    integer :: n_ptcls_imported = 0
    integer :: last_import_time = 0
    integer :: last_watch       = 0       ! time of the last watch of the upstream folder
    integer :: upstream_done_since = 0    ! when reference picking was first seen idle or stopped; 0: it is not
    logical :: l_final          = .false. ! final ingestion is set
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
    procedure :: update_final_ingestion
    procedure :: resumable
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
        if( .not. allocated(self%spproj) ) allocate(self%spproj)
        if( .not. allocated(self%qenv)   ) allocate(self%qenv)
        call create_stream_project(self%spproj, cline, string('sieve_cavgs'))
        if( .not. allocated(self%params) ) allocate(self%params)
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

    !> Restart: every set the sieve has chunked from is imported again, in the order it was
    !! chunked, and the micrographs it put in a chunk are marked chunked; the rest of a partly
    !! chunked set is left for the sieve. The sets go into the watcher history (attach_upstream).
    !! The records of a set are pushed in micrograph order, so micrograph m of a set whose first
    !! record is at n+1 is record n+m; the marked records come first, as the sieve's slicing needs.
    subroutine restore_imports( self )
        class(stream_stage_sieve), intent(inout) :: self
        type(string), allocatable :: projnames(:), sets(:), one_set(:)
        integer,      allocatable :: micinds(:)
        integer :: i, j, nsets, n_before, irec
        ! from the chunks that exist: a crash may have come between a chunk and the list
        call rebuild_chunked_mics()
        call read_chunked_mics(string(CHUNKED_MICS), projnames, micinds)
        if( .not. allocated(projnames) ) return
        if( size(projnames) == 0 ) return
        call self%send_status(string('importing previous run'))
        allocate(sets(size(projnames)))
        nsets = 0
        i     = 1
        do while( i <= size(projnames) )
            ! one set's micrographs are consecutive
            j = i
            do while( j < size(projnames) )
                if( projnames(j+1) /= projnames(i) ) exit
                j = j + 1
            end do
            if( any_set(projnames(i)) )then
                THROW_WARN('a set is listed twice in '//CHUNKED_MICS//'; its later lines are ignored')
            else if( file_exists(projnames(i)) )then
                nsets       = nsets + 1
                sets(nsets) = projnames(i)
                n_before    = self%project_list%size()
                one_set     = [projnames(i)]
                call import_new_projects(self%project_list, one_set, self%n_mics_imported, self%n_ptcls_imported)
                do irec = i,j
                    if( micinds(irec) < 1 ) cycle
                    if( n_before + micinds(irec) > self%project_list%size() ) cycle
                    call self%project_list%set_included_flags([n_before + micinds(irec), n_before + micinds(irec)])
                end do
            else
                THROW_WARN('a set the sieve chunked from is gone: '//projnames(i)%to_char())
            endif
            i = j + 1
        end do
        if( nsets > 0 ) self%restored_imports = sets(:nsets)
        write(logfhandle,'(A,I6,A,I8,A)') '>>> RESTORED ', nsets, ' SETS WITH ', size(projnames), ' CHUNKED MICROGRAPHS'

    contains

        logical function any_set( projname )
            type(string), intent(in) :: projname
            integer :: k
            any_set = .false.
            do k = 1,nsets
                if( sets(k) == projname )then
                    any_set = .true.
                    return
                endif
            end do
        end function any_set

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
            call self%update_final_ingestion()
            call self%sieve%cycle(self%project_list)
        else if( self%project_list%size() > 0 .or. self%resumable() )then
            ! a restart resumes the sieve's chunks with nothing new to import
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

    !> Cancels the running chunk jobs; the last status: no more user input.
    subroutine finalize( self )
        class(stream_stage_sieve), intent(inout) :: self
        if( self%l_sieve_active ) call self%sieve%cancel()
        call self%meta_status%set_user_input(.false.)
        call self%send_status(string('terminating'))
    end subroutine finalize

    subroutine kill( self )
        class(stream_stage_sieve), intent(inout) :: self
        if( .not. self%l_exists )then
            if( allocated(self%params) ) deallocate(self%params)
            call release_heavy
            return
        endif
        if( allocated(self%sieve) )then
            if( self%l_sieve_active ) call self%sieve%kill
            deallocate(self%sieve)
        endif
        call self%project_buff%kill
        call self%project_list%kill
        call release_heavy
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
        if( allocated(self%params) ) deallocate(self%params)
        self%latest_xtiles    = 0
        self%mskdiam          = 0.
        self%latest_ytiles    = 0
        self%n_mics_imported  = 0
        self%n_ptcls_imported = 0
        self%l_attached       = .false.
        self%l_waiting_logged = .false.
        self%l_sieve_active   = .false.
        self%l_restart        = .false.
        self%l_final          = .false.
        self%last_watch       = 0
        self%upstream_done_since = 0
        self%l_exists         = .false.

    contains

        subroutine release_heavy
            if( allocated(self%spproj) )then
                call self%spproj%kill
                deallocate(self%spproj)
            endif
            if( allocated(self%qenv) )then
                call self%qenv%kill
                deallocate(self%qenv)
            endif
        end subroutine release_heavy

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
        self%last_watch = simple_gettime()
        call self%project_buff%watch(nprojects, projects, max_nmovies=MAX_PROJECTS_IMPORT)
        ! a capped watch may have left sets for the next pass
        if( nprojects == MAX_PROJECTS_IMPORT ) self%upstream_done_since = 0
        if( nprojects == 0 ) return
        call import_new_projects(self%project_list, projects, self%n_mics_imported, self%n_ptcls_imported)
        call self%project_buff%add2history(projects)
        ! new sets: final ingestion waits for another quiet watch
        self%upstream_done_since = 0
        if( self%l_final )then
            call self%sieve%unset_final_ingestion()
            self%l_final = .false.
            write(logfhandle,'(A)') '>>> NEW SETS: FINAL INGESTION WITHDRAWN'
        endif
        self%last_import_time = simple_gettime()
        write(logfhandle,'(A,I6,I9)') '>>> # MICROGRAPHS / PARTICLES IMPORTED : ', self%n_mics_imported, self%n_ptcls_imported
        write(logfhandle,'(A,A)')     '>>> LAST IMPORT AT                     : ', cast_time_char(self%last_import_time)
    end subroutine import_projects

    ! The sieve, made on the first import with the mask diameter of the picking references, the
    ! chunk partition and the optics directory for its hand-offs, then two warm-up cycles.
    subroutine start_sieve( self )
        class(stream_stage_sieve), intent(inout) :: self
        type(ptcl_sieve_settings) :: settings
        character(len=STDLEN)     :: partition_env
        integer                   :: envlen
        call self%read_mask_diameter()
        settings         = sieve_settings(self%params)
        settings%mskdiam = self%mskdiam
        call get_environment_variable(SIMPLE_STREAM_CHUNK_PARTITION, partition_env, envlen)
        if( envlen > 0 ) settings%partition = trim(partition_env)
        if( .not. allocated(self%sieve) ) allocate(self%sieve)
        if( self%params%optics_dir%strlen() > 0 )then
            call self%sieve%new(self%params, settings, string(PATH_HERE//DIR_STREAM_COMPLETED), optics_dir=self%params%optics_dir)
        else
            call self%sieve%new(self%params, settings, string(PATH_HERE//DIR_STREAM_COMPLETED))
        endif
        self%l_sieve_active = .true.
        call self%sieve%cycle(self%project_list)
        call self%sieve%cycle(self%project_list)
    end subroutine start_sieve

    ! Final ingestion while reference picking is idle or stopped and a watch made a settle time
    ! after that found nothing (every set it handed on before its marker has settled and been
    ! taken); withdrawn when it is neither.
    subroutine update_final_ingestion( self )
        class(stream_stage_sieve), intent(inout) :: self
        if( .not. upstream_done(self%params%dir_target) )then
            self%upstream_done_since = 0
            if( self%l_final )then
                call self%sieve%unset_final_ingestion()
                self%l_final = .false.
                write(logfhandle,'(A)') '>>> REFERENCE PICKING IS ACTIVE AGAIN: FINAL INGESTION WITHDRAWN'
            endif
            return
        endif
        if( self%l_final ) return
        if( self%upstream_done_since == 0 )then
            self%upstream_done_since = simple_gettime()
            return
        endif
        if( self%last_watch - self%upstream_done_since <= max(self%settle_s, 0) ) return
        call self%sieve%set_final_ingestion()
        self%l_final = .true.
        write(logfhandle,'(A)') '>>> REFERENCE PICKING IS IDLE OR STOPPED AND EVERY SET IS TAKEN: FINAL INGESTION'
    end subroutine update_final_ingestion

    ! A restart whose sieve has chunks to take up: the mask diameter of the picking references
    ! exists (written before reference picking completes any set).
    logical function resumable( self )
        class(stream_stage_sieve), intent(in) :: self
        resumable = .false.
        if( .not. self%l_restart ) return
        resumable = file_exists(self%params%dir_target//'/'//STREAM_MOLDIAM)
    end function resumable

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
        self%mskdiam = moldiam%get(1, 'mskdiam')
        call moldiam%kill
        write(logfhandle,'(A,F8.2)') '>>> MASK DIAMETER SET TO : ', self%mskdiam
    end subroutine read_mask_diameter

    !---------------- GUI ----------------

    ! Particle counts and the classes of the latest product the sieve selected.
    subroutine send_status( self, stage )
        class(stream_stage_sieve), intent(inout) :: self
        type(string),              intent(in)    :: stage
        type(string) :: stage_text
        integer :: i, naccepted, nrejected, nfailed
        naccepted  = 0
        nrejected  = 0
        nfailed    = 0
        stage_text = stage
        if( allocated(self%sieve) )then
            naccepted = self%sieve%get_n_accepted_ptcls()
            nrejected = self%sieve%get_n_rejected_ptcls()
            nfailed   = self%sieve%get_n_failed_chunks()
        endif
        ! chunks whose 2D job failed twice are dropped with their particles (ptcl_sieve_policy.md)
        if( nfailed > 0 ) stage_text = stage//'; '//int2str(nfailed)//' chunk(s) failed, '//&
            &int2str(self%sieve%get_n_failed_ptcls())//' particles dropped'
        call self%meta_status%set(stage=stage_text, particles_imported=self%n_ptcls_imported,&
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
