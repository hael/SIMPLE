!@descr: task 00 in the stream pipeline: the master that runs the stages for the GUI
!==============================================================================
! MODULE: simple_commanders_stream_p00_master
!
! PURPOSE:
!   Runs the stream for the GUI: forks the seven stages, keeps the latest
!   metadata each one sends, posts a heartbeat to the GUI every few seconds,
!   and applies what the GUI answers (stop the stream, stop or restart a
!   stage, updates for the stages). What it uses:
!     - the stages            -> simple_stream_master_stage (one per id of simple_stream_master_stage_ids),
!                                each forked with its commander (stage_commander)
!     - their messages        -> simple_stream_master_meta_store, filled by a listener thread
!     - the GUI's answer      -> simple_stream_master_gui_commands
!     - the heartbeat         -> simple_gui_assembler, simple_http_post
!     - SIGTERM and SIGINT    -> simple_stream_sigterm (a flag, polled by the loop)
!
! STOPPING:
!   The GUI's 'terminate', SIGTERM or Ctrl-C stop the stream in order: the
!   optics assignment first (waited for up to a minute, so the other stages
!   see its last map), then every other stage but the large ones together,
!   then the large ones (STOP_ONE_BY_ONE: multistate 3D, pool 2D, reference
!   picking) one after another, each asked once the one before has stopped, so
!   their final writes do not hold their memory at the same time. A stage still
!   running STOP_TIMEOUT_S after it was asked is killed. The master exits once
!   all have stopped, after a last heartbeat.
!==============================================================================
module simple_commanders_stream_p00_master
use, intrinsic :: iso_c_binding,                        only: c_ptr, c_null_ptr, c_loc, c_f_pointer, c_funloc
use unix,                                               only: c_pthread_t, c_pthread_mutex_t, c_pthread_create, c_pthread_join,&
                                                             &c_pthread_mutex_init, c_pthread_mutex_destroy, c_pthread_mutex_lock,&
                                                             &c_pthread_mutex_unlock, c_usleep
use simple_defs,                                        only: logfhandle
use simple_defs_fname,                                  only: METADATA_EXT
use simple_defs_stream,                                 only: PREPROC_JOB_NAME, OPTICS_JOB_NAME, INITIAL_ANALYSIS_JOB_NAME, REFPICK_JOB_NAME,&
                                                             &INITIAL_ANALYSIS_PICKREFS, SIEVING_JOB_NAME, CLASS2D_JOB_NAME, MULTISTATE3D_JOB_NAME, PREPROC_NINIPICK,&
                                                             &STREAM_IDLE_MARKER, STREAM_FINISHED_MARKER
use simple_error,                                       only: simple_exception
use simple_string,                                      only: string
use simple_fileio,                                      only: simple_getcwd, file_exists, simple_touch
use simple_syslib,                                      only: dir_exists, symlink, simple_abspath
use simple_timer,                                       only: simple_gettime
use simple_cmdline,                                     only: cmdline
use simple_parameters,                                  only: parameters
use simple_commander_base,                              only: commander_base
use simple_qsys_env,                                    only: qsys_env
use simple_http_post,                                   only: http_post, http_response
use simple_gui_assembler,                               only: gui_assembler
use simple_gui_metadata_utils,                          only: max_metadata_size
use simple_memory_monitor,                              only: mem_monitor_init, mem_monitor_finish
use simple_stream_master_stage_ids,                     only: NSTAGES, STAGE_PREPROCESS, STAGE_ASSIGN_OPTICS, STAGE_INITIAL_ANALYSIS,&
                                                             &STAGE_REFERENCE_PICKING, STAGE_PARTICLE_SIEVING, STAGE_POOL2D,&
                                                             &STAGE_SOLVE3D
use simple_stream_master_stage,                         only: stream_master_stage, fork_gui_status
use simple_stream_master_resources,                     only: stream_resources, stream_resources_from_env
use simple_stream_master_meta_store,                    only: stream_master_meta_store
use simple_stream_master_gui_commands,                  only: stream_master_gui_commands
use simple_stream_sigterm,                              only: install_sigterm_handler, restore_sigterm_handler, sigterm_received
implicit none

public :: commander_stream_p00_master
private
#include "simple_local_flags.inc"

! what the master gives the stages (its choices for a stream run from the GUI); their threads
! and parts are in simple_stream_master_resources
integer, parameter :: POOL2D_NCLS           = 150
! the initial analysis' 3D route settings (decision 20): forwarded to it only when given, kept off
! preprocessing's command line (it is the master's own); nthr3D_pickrefs also has a table default
character(len=*), parameter :: PICKREFS_3D_KEYS(7) = [character(len=18) :: 'nstates_pickrefs', 'nstages_pickrefs',&
    &'lpstop_pickrefs', 'nspace_pickrefs', 'nrestarts_collapse', 'lpstart_ini3D', 'lpstop_ini3D']
! the master itself
integer, parameter :: MASTER_NCUNITS        = 16    ! computing units of the master's queue environment
integer, parameter :: MASTER_QSYS_NTHR      = 16
integer, parameter :: HEARTBEAT_S           = 5     ! between heartbeats
integer, parameter :: MEMLOG_EVERY          = 12    ! heartbeats between memory logs
integer, parameter :: OPTICS_STOP_TIMEOUT_S = 60    ! optics assignment stops first, waited for up to this
integer, parameter :: STOP_TIMEOUT_S        = 600   ! a stage still running this long after it was asked to stop is killed
! the large stages, stopped one after another (downstream first), each when the one before has stopped
integer, parameter :: STOP_ONE_BY_ONE(3)    = [STAGE_SOLVE3D, STAGE_POOL2D, STAGE_REFERENCE_PICKING]
integer, parameter :: STARTUP_STOP_TIMEOUT_S = 60   ! after a start-up failure, the stages are killed after this
integer, parameter :: LISTENER_SLEEP_US     = 10000 ! the listener's pause between drains

! What the listener thread shares with the main loop; it gets the address.
type :: master_shared
    type(stream_master_meta_store) :: store
    type(stream_master_stage) :: stages(NSTAGES)
    type(c_pthread_mutex_t) :: meta_mutex      ! the store, and the stages' readers
    type(c_pthread_mutex_t) :: stop_mutex      ! l_listener_stop
    logical                 :: l_listener_stop = .false.
end type master_shared

type, extends(commander_base) :: commander_stream_p00_master
  contains
    procedure :: execute => exec_stream_p00_master
end type commander_stream_p00_master

contains

    subroutine exec_stream_p00_master( self, cline )
        class(commander_stream_p00_master), intent(inout) :: self
        class(cmdline),           intent(inout) :: cline
        type(master_shared), target   :: shared
        type(parameters)              :: params
        type(cmdline)                 :: clines(NSTAGES)
        type(http_post)               :: post
        type(http_response)           :: response
        type(gui_assembler)           :: assembler
        type(qsys_env)                :: qsys
        type(stream_master_gui_commands)     :: commands
        type(c_pthread_t)             :: listener
        type(c_ptr)                   :: thread_ret
        type(string)                  :: request, cwd
        character(len=:), allocatable :: update_buffer
        logical :: l_existing_pickrefs, l_existing_box, l_existing_preprocess, l_updates, l_linked
        logical :: l_stop, l_stopping, l_last_loop
        logical :: l_restart_seen(NSTAGES) ! a restart request in the last answer, acted on once
        integer :: id, nmics_stop, loop_counter, rc, max_frame_bytes
        integer :: stop_times(NSTAGES) ! when each stage was asked to stop (0: not yet)
        ! the command line
        l_existing_pickrefs = cline%defined('pickrefs')
        l_existing_box      = cline%defined('box_extract')
        nmics_stop          = 0
        if( cline%defined('nmics') )then
            nmics_stop = cline%get_iarg('nmics')
            call cline%delete('nmics')
        endif
        if( .not. cline%defined('smpd_downscale') ) call cline%set('smpd_downscale', 1.3)
        if( .not. cline%defined('memreport')      ) call cline%set('memreport',      'yes')
        call cline%printline()
        call params%new(cline)
        call simple_getcwd(cwd)
        l_existing_preprocess = params%dir_preprocess%strlen() > 0
        if( l_existing_preprocess )then
            if( .not. dir_exists(params%dir_preprocess) ) THROW_HARD('Preprocessing directory '//params%dir_preprocess%to_char()//' does not exist.')
        endif
        ! the GUI link
        call post%new(params%niceserver)
        call assembler%new(params%niceprocid)
        ! locks shared with the listener
        if( c_pthread_mutex_init(shared%meta_mutex, c_null_ptr) /= 0 ) THROW_HARD('failed to initialise metadata mutex')
        if( c_pthread_mutex_init(shared%stop_mutex, c_null_ptr) /= 0 ) THROW_HARD('failed to initialise stop mutex')
        ! the queue environment, with the persistent worker server when the environment asks for it
        params%qsys_name = '' ! read from the environment variables
        params%ncunits   = MASTER_NCUNITS
        call qsys%new(params, 1, qsys_nthr=MASTER_QSYS_NTHR, stream=.true.)
        ! the stages, and every pipe before any stage is forked
        max_frame_bytes = max_metadata_size()
        write(logfhandle,*) 'Max metadata size: ', max_frame_bytes
        call shared%store%new()
        call make_stage_clines()
        do id = 1,NSTAGES
            l_updates = id == STAGE_PREPROCESS .or. id == STAGE_INITIAL_ANALYSIS .or. id == STAGE_POOL2D .or.&
                &id == STAGE_SOLVE3D
            call shared%stages(id)%new(id, stage_commander(id), clines(id), l_updates, max_frame_bytes)
        enddo
        ! fork the stages; preprocessing and the initial analysis are skipped when given their outputs.
        ! The listener thread is started after them: a child gets a copy of the parent's memory but
        ! not of its threads, so the first forks are made before the listener can be caught
        ! mid-write. The persistent-worker server's thread (qsys%new) does run already: each stage
        ! forgets the server it inherits before its commander runs (stream_master_stage_fork) and
        ! reaches it as a client through worker_server
        if( l_existing_preprocess )then
            rc = symlink(params%dir_preprocess%to_char()//achar(0), PREPROC_JOB_NAME//achar(0))
            if( rc /= 0 )then
                ! a master restarted in the same folder finds its own link
                l_linked = dir_exists(string(PREPROC_JOB_NAME))
                if( l_linked ) l_linked = simple_abspath(PREPROC_JOB_NAME) == simple_abspath(params%dir_preprocess)
                if( .not. l_linked ) THROW_HARD('failed to create symlink for existing preprocessing directory')
            endif
            ! nothing more comes from the earlier preprocessing: the stages that end their intake
            ! when preprocessing goes idle or stops see it stopped
            if( .not. file_exists(PREPROC_JOB_NAME//'/'//STREAM_IDLE_MARKER) .and.&
                &.not. file_exists(PREPROC_JOB_NAME//'/'//STREAM_FINISHED_MARKER) )then
                call simple_touch(PREPROC_JOB_NAME//'/'//STREAM_FINISHED_MARKER)
                write(logfhandle,'(A)') '>>> THE EXISTING PREPROCESSING IS MARKED FINISHED'
            endif
            call shared%stages(STAGE_PREPROCESS)%skip()
        else
            call shared%stages(STAGE_PREPROCESS)%start()
        endif
        call shared%stages(STAGE_ASSIGN_OPTICS)%start()
        call shared%stages(STAGE_REFERENCE_PICKING)%start()
        call shared%stages(STAGE_PARTICLE_SIEVING)%start()
        call shared%stages(STAGE_POOL2D)%start()
        call shared%stages(STAGE_SOLVE3D)%start()
        if( l_existing_pickrefs )then
            call shared%stages(STAGE_INITIAL_ANALYSIS)%skip()
        else
            call shared%stages(STAGE_INITIAL_ANALYSIS)%start()
        endif
        do id = 1,NSTAGES
            if( id == STAGE_PREPROCESS       .and. l_existing_preprocess ) cycle
            if( id == STAGE_INITIAL_ANALYSIS .and. l_existing_pickrefs   ) cycle
            if( .not. shared%stages(id)%is_running() ) call abort_startup('failed to fork '//shared%stages(id)%get_label())
        enddo
        ! the listener thread; the stages' first messages wait in their pipes until it reads them
        rc = c_pthread_create(listener, c_null_ptr, c_funloc(metadata_listener), c_loc(shared))
        if( rc /= 0 ) call abort_startup('failed to create metadata listener thread')
        ! the handlers after the forks (a forked stage resets them anyway)
        call install_sigterm_handler(also_sigint=.true.)
        call mem_monitor_init(cline, 'simple_stream: master')
        ! heartbeats until the stream has stopped
        loop_counter = 0
        l_stop       = .false.
        l_stopping   = .false.
        l_last_loop  = .false.
        l_restart_seen = .false.
        stop_times   = 0
        do
            loop_counter = loop_counter + 1
            call send_heartbeat()
            if( mod(loop_counter, MEMLOG_EVERY) == 0 ) call log_memory()
            if( l_last_loop ) exit
            if( l_stop .or. sigterm_received() )then
                call stop_stages()
                call assembler%set_stoptime()
                if( .not. l_last_loop ) call sleep(HEARTBEAT_S)
            else
                call sleep(HEARTBEAT_S)
            endif
            call flush(logfhandle)
        enddo
        ! stop the listener, then release everything
        call lock(shared%stop_mutex)
        shared%l_listener_stop = .true.
        call unlock(shared%stop_mutex)
        if( c_pthread_join(listener, thread_ret) /= 0 ) THROW_WARN('failed to join metadata listener thread')
        call restore_sigterm_handler()
        call commands%kill()
        call qsys%kill()
        call post%kill()
        call assembler%kill()
        call mem_monitor_finish()
        do id = 1,NSTAGES
            call shared%stages(id)%kill()
            call clines(id)%kill()
        enddo
        call shared%store%kill()
        if( c_pthread_mutex_destroy(shared%meta_mutex) /= 0 ) THROW_WARN('failed to destroy metadata mutex')
        if( c_pthread_mutex_destroy(shared%stop_mutex) /= 0 ) THROW_WARN('failed to destroy stop mutex')
        call flush(logfhandle)
        ! back to the entry point, which stops the persistent worker server this process started
        ! and closes the log (a forked stage exits directly instead: the server is not its own)

    contains

        ! One heartbeat: the stages' state to the GUI, and its answer applied. A heartbeat the GUI
        ! did not accept (no answer, or a status other than 200) may not have been read: the next
        ! one sends everything again, not only what changed since.
        subroutine send_heartbeat()
            call assembler%assemble_stream_heartbeat(fork_gui_status(shared%stages(STAGE_PREPROCESS)%fork),&
                &fork_gui_status(shared%stages(STAGE_ASSIGN_OPTICS)%fork),&
                &fork_gui_status(shared%stages(STAGE_INITIAL_ANALYSIS)%fork),&
                &fork_gui_status(shared%stages(STAGE_REFERENCE_PICKING)%fork),&
                &fork_gui_status(shared%stages(STAGE_PARTICLE_SIEVING)%fork),&
                &fork_gui_status(shared%stages(STAGE_POOL2D)%fork), fork_gui_status(shared%stages(STAGE_SOLVE3D)%fork),&
                &n_active_persistent_workers=qsys%get_n_active_persistent_workers())
            call lock(shared%meta_mutex)
            call shared%store%assemble(assembler)
            call unlock(shared%meta_mutex)
            request = assembler%to_string()
            if( post%request(response, request) )then
                if( response%code == 200 )then
                    if( commands%parse(response%content%to_char()) ) call apply_commands()
                else
                    call assembler%clear_hashes()
                endif
            else
                call assembler%clear_hashes()
            endif
            call request%kill()
            call response%content%kill()
            call response%content_type%kill()
            call qsys%service_persistent_worker_warmup()
        end subroutine send_heartbeat

        ! A start-up failure once stages are forked: they are asked to stop, and killed (their
        ! jobs cancelled) when still running after STARTUP_STOP_TIMEOUT_S, before the master
        ! stops. The persistent-worker server ends with this process, and its workers, which
        ! leave when they lose it, with it.
        subroutine abort_startup( msg )
            character(len=*), intent(in) :: msg
            integer, parameter :: POLL_US = 200000
            integer :: jd, ipoll, rc_wait
            logical :: l_any
            write(logfhandle,'(A)') '>>> START-UP FAILED: '//msg//'; STOPPING THE STAGES'
            do jd = 1,NSTAGES
                if( shared%stages(jd)%is_running() ) call shared%stages(jd)%request_stop()
            enddo
            do ipoll = 1,(STARTUP_STOP_TIMEOUT_S * 1000000) / POLL_US
                l_any = .false.
                do jd = 1,NSTAGES
                    if( shared%stages(jd)%is_running() ) l_any = .true.
                enddo
                if( .not. l_any ) exit
                rc_wait = c_usleep(POLL_US)
            enddo
            do jd = 1,NSTAGES
                call shared%stages(jd)%force_stop()
            enddo
            call flush(logfhandle)
            THROW_HARD(msg)
        end subroutine abort_startup

        ! What the GUI asked in its answer to the last heartbeat. A restarted stage starts on clean
        ! pipes, and is forked, under the listener's lock: the listener reads the pipes and logs
        ! only while it holds it, so the child is never forked with a log write half done. NICE
        ! keeps a restart key in its answers until it sees the stage running, so a request is acted
        ! on once, and again only after the key has left an answer; and never once the stream is
        ! stopping. The updates of this answer go to the running stages that read them.
        subroutine apply_commands()
            logical :: l_new_request
            if( commands%l_terminate_all ) l_stop = .true.
            do id = 1,NSTAGES
                if( commands%l_terminate(id) ) call shared%stages(id)%request_stop()
            enddo
            do id = 1,NSTAGES
                l_new_request      = commands%l_restart(id) .and. .not. l_restart_seen(id)
                l_restart_seen(id) = commands%l_restart(id)
                if( .not. l_new_request ) cycle
                if( l_stop .or. l_stopping .or. sigterm_received() ) cycle
                if( shared%stages(id)%is_running() ) cycle
                ! a skipped stage's output is the user's earlier run: it is never run here
                if( shared%stages(id)%is_skipped() )then
                    write(logfhandle,'(A)') '>>> RESTART OF SKIPPED '//shared%stages(id)%get_label()//' IGNORED'
                    cycle
                endif
                call lock(shared%meta_mutex)
                call shared%stages(id)%discard_pipes(max_frame_bytes)
                call shared%store%clear_stage(id)
                call shared%stages(id)%start()
                call unlock(shared%meta_mutex)
            enddo
            if( commands%update%assigned() )then
                call commands%update%serialise(update_buffer)
                do id = 1,NSTAGES
                    call shared%stages(id)%send_update(update_buffer)
                enddo
            endif
        end subroutine apply_commands

        ! One step of the stop: the first time, the optics assignment and then every other stage but
        ! the large ones (STOP_ONE_BY_ONE) are asked to stop; every time, the first large stage still
        ! running is asked once those before it have stopped, the stages asked and still running are
        ! asked again and listed, and each is killed once STOP_TIMEOUT_S has passed since it was
        ! asked. A large stage waiting its turn runs on. l_last_loop once none runs.
        subroutine stop_stages()
            integer :: k
            if( .not. l_stopping )then
                write(logfhandle,'(A)') 'TERMINATE '
                ! stopping from here on: the heartbeats of the wait act on no restart request
                l_stopping = .true.
                if( shared%stages(STAGE_ASSIGN_OPTICS)%is_running() )then
                    call shared%stages(STAGE_ASSIGN_OPTICS)%request_stop()
                    stop_times(STAGE_ASSIGN_OPTICS) = simple_gettime()
                    call wait_for_stop(shared%stages(STAGE_ASSIGN_OPTICS), OPTICS_STOP_TIMEOUT_S)
                endif
                do id = 1,NSTAGES
                    if( any(STOP_ONE_BY_ONE == id) ) cycle
                    call ask_to_stop(id)
                enddo
            endif
            ! the large stages one after another: the first still running is asked, the later wait
            do k = 1,size(STOP_ONE_BY_ONE)
                if( .not. shared%stages(STOP_ONE_BY_ONE(k))%is_running() ) cycle
                call ask_to_stop(STOP_ONE_BY_ONE(k))
                exit
            enddo
            l_last_loop = .true.
            do id = 1,NSTAGES
                if( .not. shared%stages(id)%is_running() ) cycle
                l_last_loop = .false.
                if( stop_times(id) == 0 ) cycle ! a large stage waiting its turn
                call shared%stages(id)%request_stop()
                write(logfhandle,'(A)') shared%stages(id)%get_label()//' STILL RUNNING. WAITING FOR TERMINATION'
                if( simple_gettime() - stop_times(id) > STOP_TIMEOUT_S ) call shared%stages(id)%force_stop()
            enddo
        end subroutine stop_stages

        ! Asks stage @p id to stop, once, when it runs, and records when.
        subroutine ask_to_stop( id )
            integer, intent(in) :: id
            if( stop_times(id) > 0 ) return
            if( .not. shared%stages(id)%is_running() ) return
            call shared%stages(id)%request_stop()
            stop_times(id) = simple_gettime()
            write(logfhandle,'(A)') shared%stages(id)%get_label()//' ASKED TO STOP'
        end subroutine ask_to_stop

        ! Waits for @p stage to stop, for up to @p timeout_s, with a heartbeat every HEARTBEAT_S so
        ! the GUI does not lose the master meanwhile.
        subroutine wait_for_stop( stage, timeout_s )
            type(stream_master_stage), intent(inout) :: stage
            integer,                 intent(in)    :: timeout_s
            integer, parameter :: POLL_US = 200000
            integer :: ipoll, rc_wait
            do ipoll = 1,(timeout_s * 1000000) / POLL_US
                if( .not. stage%is_running() ) return
                if( mod(ipoll, (HEARTBEAT_S * 1000000) / POLL_US) == 0 ) call send_heartbeat()
                rc_wait = c_usleep(POLL_US)
            enddo
            write(logfhandle,'(A)') stage%get_label()//' DID NOT TERMINATE WITHIN TIMEOUT'
        end subroutine wait_for_stop

        subroutine log_memory()
            integer :: pending, jd
            call lock(shared%meta_mutex)
            pending = 0
            do jd = 1,NSTAGES
                pending = pending + shared%stages(jd)%reader%get_pending_bytes()
            enddo
            call shared%store%log_state(pending)
            call unlock(shared%meta_mutex)
        end subroutine log_memory

        ! The stages' command lines: their programs, folders and links, and the master's settings.
        ! Preprocessing gets the master's own command line with the user's preprocessing options.
        subroutine make_stage_clines()
            type(stream_resources) :: res
            type(string)           :: server_address
            integer                :: ikey
            server_address = qsys%get_persistent_worker_server_address()
            ! the threads and parts of every stage and its jobs: defaults, the stages' environment
            ! variables over them
            res = stream_resources_from_env()
            call res%log()
            ! preprocessing
            clines(STAGE_PREPROCESS) = cline
            associate( c => clines(STAGE_PREPROCESS) )
                call c%set('prg',      'preproc')
                call c%set('projfile', PREPROC_JOB_NAME//METADATA_EXT)
                call c%set('outdir',   PREPROC_JOB_NAME)
                call c%set('ninipick', PREPROC_NINIPICK)
                call c%set('nparts',   res%preprocess_nparts)
                call c%set('nthr',     res%preprocess_nthr)
                call c%set('mkdir',    'yes')
                call c%delete('niceserver')
                call c%delete('niceprocid')
                call c%delete('box_extract')
                call c%delete('pickrefs')
                call c%delete('nthr3D_pickrefs')
                call c%delete('sieve_ini3D')
                do ikey = 1,size(PICKREFS_3D_KEYS)
                    call c%delete(trim(PICKREFS_3D_KEYS(ikey)))
                enddo
                if( nmics_stop > 0 ) call c%set('nmics', nmics_stop)
                if( server_address%strlen() > 0 ) call c%set('worker_server', server_address)
            end associate
            ! optics assignment, with the beam-tilt options entered with the preprocessing ones
            associate( c => clines(STAGE_ASSIGN_OPTICS) )
                call c%set('prg',        'assign_optics')
                call c%set('projfile',   OPTICS_JOB_NAME//METADATA_EXT)
                call c%set('outdir',     OPTICS_JOB_NAME)
                call c%set('dir_target', PREPROC_JOB_NAME)
                call c%set('nthr',       res%optics_nthr)
                call c%set('mkdir',      'yes')
                if( cline%defined('beamtilt')   ) call c%set('beamtilt',   trim(params%beamtilt))
                if( cline%defined('tilt_thres') ) call c%set('tilt_thres', params%tilt_thres)
            end associate
            ! initial analysis
            associate( c => clines(STAGE_INITIAL_ANALYSIS) )
                call c%set('prg',             'gen_pickrefs')
                call c%set('projfile',        INITIAL_ANALYSIS_JOB_NAME//METADATA_EXT)
                call c%set('outdir',          INITIAL_ANALYSIS_JOB_NAME)
                call c%set('dir_target',      PREPROC_JOB_NAME)
                call c%set('optics_dir',      cwd//'/'//OPTICS_JOB_NAME)
                call c%set('nthr',            res%initial_analysis_nthr)
                call c%set('nthr2D',          res%initial_analysis_nthr2D)
                call c%set('nparts',          res%initial_analysis_nparts)
                call c%set('nchunks',         res%initial_analysis_nchunks)
                ! the user's 3D threads win over the table's
                if( cline%defined('nthr3D_pickrefs') )then
                    call c%set('nthr3D_pickrefs', params%nthr3D_pickrefs)
                else
                    call c%set('nthr3D_pickrefs', res%initial_analysis_nthr3D)
                endif
                ! the 3D route's settings the user gave; the stage's commander defaults the others
                do ikey = 1,size(PICKREFS_3D_KEYS)
                    call c%copy_arg(cline, trim(PICKREFS_3D_KEYS(ikey)))
                enddo
                call c%set('mkdir',           'yes')
                call c%set('worker_priority', 'high')
                if( server_address%strlen() > 0 ) call c%set('worker_server', server_address)
            end associate
            ! reference picking, with the given references or those of the initial analysis
            associate( c => clines(STAGE_REFERENCE_PICKING) )
                call c%set('prg',        'pick_extract')
                call c%set('projfile',   REFPICK_JOB_NAME//METADATA_EXT)
                call c%set('outdir',     REFPICK_JOB_NAME)
                call c%set('dir_target', PREPROC_JOB_NAME)
                call c%set('optics_dir', cwd//'/'//OPTICS_JOB_NAME)
                call c%set('nthr',       res%refpick_nthr)
                call c%set('nparts',     res%refpick_nparts)
                call c%set('mkdir',      'yes')
                if( l_existing_pickrefs )then
                    call c%set('pickrefs', params%pickrefs)
                else
                    call c%set('pickrefs', '../'//INITIAL_ANALYSIS_JOB_NAME//'/'//INITIAL_ANALYSIS_PICKREFS)
                endif
                if( l_existing_box       ) call c%set('box_extract', params%box_extract)
                if( params%thres > 0.0   ) call c%set('thres',       params%thres)
                if( server_address%strlen() > 0 ) call c%set('worker_server', server_address)
            end associate
            ! particle sieving
            associate( c => clines(STAGE_PARTICLE_SIEVING) )
                call c%set('prg',             'sieve_cavgs')
                call c%set('projfile',        SIEVING_JOB_NAME//METADATA_EXT)
                call c%set('outdir',          SIEVING_JOB_NAME)
                call c%set('dir_target',      REFPICK_JOB_NAME)
                call c%set('optics_dir',      cwd//'/'//OPTICS_JOB_NAME)
                call c%set('nthr',            res%sieve_nthr)
                call c%set('nchunks',         res%sieve_nchunks)
                call c%set('mkdir',           'yes')
                call c%set('worker_priority', 'high')
                if( server_address%strlen() > 0 ) call c%set('worker_server', server_address)
            end associate
            ! pool 2D
            associate( c => clines(STAGE_POOL2D) )
                call c%set('prg',             'pool2D')
                call c%set('projfile',        CLASS2D_JOB_NAME//METADATA_EXT)
                call c%set('outdir',          CLASS2D_JOB_NAME)
                call c%set('dir_target',      SIEVING_JOB_NAME)
                call c%set('optics_dir',      cwd//'/'//OPTICS_JOB_NAME)
                call c%set('nthr',            res%pool2D_nthr)
                call c%set('nparts',          res%pool2D_nparts)
                call c%set('ncls',            POOL2D_NCLS)
                call c%set('mkdir',           'yes')
                call c%set('nicedispid',      params%nicedispid)
                ! the first publication from the sieve's class averages (multistate 3D follows its marker)
                if( params%sieve_ini3D == 'yes' ) call c%set('sieve_ini3D', 'yes')
                call c%set('worker_priority', 'high')
                if( server_address%strlen() > 0 ) call c%set('worker_server', server_address)
            end associate
            ! multistate 3D
            associate( c => clines(STAGE_SOLVE3D) )
                call c%set('prg',             'solve3D_stream')
                call c%set('projfile',        MULTISTATE3D_JOB_NAME//METADATA_EXT)
                call c%set('outdir',          MULTISTATE3D_JOB_NAME)
                call c%set('dir_target',      CLASS2D_JOB_NAME)
                call c%set('nthr',            res%solve3D_nthr)
                call c%set('nthr3D',          res%solve3D_nthr3D)
                call c%set('nparts3D',        res%solve3D_nparts3D)
                call c%set('mkdir',           'yes')
                call c%set('nicedispid',      params%nicedispid)
                call c%set('worker_priority', 'high')
                if( server_address%strlen() > 0 ) call c%set('worker_server', server_address)
            end associate
            ! the worker server's threads per worker beside its address: a stage's queues check their
            ! jobs' claims against it (qsys_env)
            if( server_address%strlen() > 0 .and. qsys%get_persistent_worker_nthr() > 0 )then
                do id = 1,NSTAGES
                    if( clines(id)%defined('worker_server') )&
                        &call clines(id)%set('worker_server_nthr', qsys%get_persistent_worker_nthr())
                enddo
            endif
            ! every stage reports its memory when the master does
            if( params%memreport == 'yes' )then
                do id = 1,NSTAGES
                    call clines(id)%set('memreport',          'yes')
                    call clines(id)%set('memreport_interval', params%memreport_interval)
                enddo
            endif
        end subroutine make_stage_clines

    end subroutine exec_stream_p00_master

    ! The commander that runs stage @p id in its forked process. The stage commanders are imported
    ! here only, the one place that needs them (compile-time policy).
    function stage_commander( id ) result( commander )
        use simple_commanders_stream_p01_preprocess,         only: commander_stream_p01_preprocess
        use simple_commanders_stream_p02_assign_optics,      only: commander_stream_p02_assign_optics
        use simple_commanders_stream_p03_initial_analysis,   only: commander_stream_p03_initial_analysis
        use simple_commanders_stream_p04_refpick_extract,    only: commander_stream_p04_refpick_extract
        use simple_commanders_stream_p05_sieve_cavgs,        only: commander_stream_p05_sieve_cavgs
        use simple_commanders_stream_p06_pool2D,             only: commander_stream_p06_pool2D
        use simple_commanders_stream_p07_solve3D_multistate, only: commander_stream_p07_solve3D_multistate
        integer, intent(in) :: id
        class(commander_base), allocatable :: commander
        select case(id)
            case(STAGE_PREPROCESS);        allocate(commander_stream_p01_preprocess            :: commander)
            case(STAGE_ASSIGN_OPTICS);     allocate(commander_stream_p02_assign_optics         :: commander)
            case(STAGE_INITIAL_ANALYSIS);  allocate(commander_stream_p03_initial_analysis      :: commander)
            case(STAGE_REFERENCE_PICKING); allocate(commander_stream_p04_refpick_extract       :: commander)
            case(STAGE_PARTICLE_SIEVING);  allocate(commander_stream_p05_sieve_cavgs           :: commander)
            case(STAGE_POOL2D);            allocate(commander_stream_p06_pool2D                :: commander)
            case(STAGE_SOLVE3D);           allocate(commander_stream_p07_solve3D_multistate    :: commander)
            case default;                  THROW_HARD('unknown stream stage id')
        end select
    end function stage_commander

    ! The listener thread: moves every message the stages have sent into the store, under the
    ! metadata lock, until the main loop asks it to stop.
    subroutine metadata_listener( arg ) bind(c)
        type(c_ptr), value, intent(in) :: arg
        type(master_shared), pointer  :: shared
        character(len=:), allocatable :: buffer
        logical :: l_continue
        integer :: id, rc
        call c_f_pointer(arg, shared)
        l_continue = .true.
        do while( l_continue )
            call lock(shared%meta_mutex)
            do id = 1,NSTAGES
                do while( shared%stages(id)%reader%receive(buffer) )
                    call shared%store%store(buffer)
                enddo
            enddo
            call unlock(shared%meta_mutex)
            call lock(shared%stop_mutex)
            l_continue = .not. shared%l_listener_stop
            call unlock(shared%stop_mutex)
            rc = c_usleep(LISTENER_SLEEP_US)
        enddo
    end subroutine metadata_listener

    subroutine lock( mutex )
        type(c_pthread_mutex_t), intent(inout) :: mutex
        if( c_pthread_mutex_lock(mutex) /= 0 ) THROW_HARD('failed to lock a master mutex')
    end subroutine lock

    subroutine unlock( mutex )
        type(c_pthread_mutex_t), intent(inout) :: mutex
        if( c_pthread_mutex_unlock(mutex) /= 0 ) THROW_HARD('failed to unlock a master mutex')
    end subroutine unlock

end module simple_commanders_stream_p00_master
