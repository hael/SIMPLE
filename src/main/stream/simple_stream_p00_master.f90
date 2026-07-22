!@descr: task 00 in the stream pipeline: master controller used when running from GUI
!==============================================================================
! MODULE: simple_stream_p00_master
!
! PURPOSE:
!   Orchestrates the full stream pipeline from the GUI side. The master process
!   launches stage workers, aggregates metadata from stage->master pipes,
!   assembles heartbeat/progress payloads for NICE, and routes GUI updates back
!   to workers via master->stage pipes.
!
! RESPONSIBILITIES:
!   - Initialise stream stage command lines and metadata containers.
!   - Create and manage IPC pipes for all stages.
!   - Spawn/terminate/restart forked stage processes.
!   - Run a metadata listener thread and merge framed pipe messages.
!   - Poll NICE for control updates and broadcast framed update messages.
!   - Perform orderly shutdown (thread join, pipe close, mutex cleanup).
!==============================================================================
module simple_stream_p00_master
use unix
use simple_syslib,                         only: symlink
use simple_stream_api
use simple_stream_state             
use simple_stream_p01_preprocess_new,      only: stream_p01_preprocess
use simple_stream_p02_assign_optics_new,   only: stream_p02_assign_optics
use simple_stream_p03_initial_analysis,    only: stream_p03_initial_analysis
use simple_stream_p04_refpick_extract_new, only: stream_p04_refpick_extract
use simple_stream_p05_sieve_cavgs_new,     only: stream_p05_sieve_cavgs
use simple_stream_p06_pool2D_new,          only: stream_p06_pool2D
use simple_http_post,                      only: http_post, http_response
use simple_forked_process,                 only: forked_process, FORK_STATUS_RUNNING
use simple_gui_metadata_api
use simple_gui_assembler,                  only: gui_assembler
use simple_gui_metadata_utils,             only: max_metadata_size

implicit none

type(c_pthread_mutex_t), save :: terminate_mutex
type(c_pthread_mutex_t), save :: meta_mutex

public :: stream_p00_master
private
#include "simple_local_flags.inc"

!================ FORKED PROCESSES =================

type, extends(forked_process) :: preprocess_fork
    contains
    procedure :: execute => xpreprocess
end type preprocess_fork

type, extends(forked_process) :: assign_optics_fork
    contains
    procedure :: execute => xassign_optics
end type assign_optics_fork

type, extends(forked_process) :: initial_analysis_fork
    contains
    procedure :: execute => xinitial_analysis
end type initial_analysis_fork

type, extends(forked_process) :: reference_picking_fork
    contains
    procedure :: execute => xreference_picking
end type reference_picking_fork

type, extends(forked_process) :: particle_sieving_fork
    contains
    procedure :: execute => xparticle_sieving
end type particle_sieving_fork

type, extends(forked_process) :: pool2D_fork
    contains
    procedure :: execute => xpool2D
end type pool2D_fork

!================ STATE TYPES =================

type pipe_rx_state
    character(len=:), allocatable :: pending
    integer                       :: expected_len = -1
end type pipe_rx_state

!=========================================================

type, extends(commander_base) :: stream_p00_master

contains

    procedure :: execute  => exec_stream_p00_master

end type stream_p00_master

contains

    subroutine xpreprocess( self, cline )
        class(preprocess_fork), intent(inout) :: self
        class(cmdline),         intent(inout) :: cline
        type(stream_p01_preprocess)           :: commander
        call close_child_pipe_fds(ipc_pipe_preprocess_in(2), ipc_pipe_preprocess_out(1))
        call commander%execute(cline)
    end subroutine xpreprocess

    subroutine xassign_optics( self, cline )
        class(assign_optics_fork), intent(inout) :: self
        class(cmdline),            intent(inout) :: cline
        type(stream_p02_assign_optics)           :: commander
        call close_child_pipe_fds(ipc_pipe_assign_optics_in(2), ipc_pipe_assign_optics_out(1))
        call commander%execute(cline)
    end subroutine xassign_optics

    subroutine xinitial_analysis( self, cline )
        class(initial_analysis_fork), intent(inout) :: self
        class(cmdline),        intent(inout) :: cline
        type(stream_p03_initial_analysis)           :: commander
        call close_child_pipe_fds(ipc_pipe_initial_analysis_in(2), ipc_pipe_initial_analysis_out(1))
        call commander%execute(cline)
    end subroutine xinitial_analysis

    subroutine xreference_picking( self, cline )
        class(reference_picking_fork), intent(inout) :: self
        class(cmdline),                intent(inout) :: cline
        type(stream_p04_refpick_extract)             :: commander
        call close_child_pipe_fds(ipc_pipe_refpick_in(2), ipc_pipe_refpick_out(1))
        call commander%execute(cline)
    end subroutine xreference_picking
 
    subroutine xparticle_sieving( self, cline )
        class(particle_sieving_fork), intent(inout) :: self
        class(cmdline),               intent(inout) :: cline
        type(stream_p05_sieve_cavgs)                :: commander
        call close_child_pipe_fds(ipc_pipe_sieve_cavgs_in(2), ipc_pipe_sieve_cavgs_out(1))
        call commander%execute(cline)
    end subroutine xparticle_sieving  

    subroutine xpool2D( self, cline )
        class(pool2D_fork),   intent(inout) :: self
        class(cmdline),       intent(inout) :: cline
        type(stream_p06_pool2D)             :: commander
        call close_child_pipe_fds(ipc_pipe_pool2D_in(2), ipc_pipe_pool2D_out(1))
        call commander%execute(cline)
    end subroutine xpool2D   

    subroutine close_child_pipe_fds( keep_fd, keep_fd2 )
        integer, intent(in)           :: keep_fd
        integer, intent(in), optional :: keep_fd2
        call close_pipe_except_fd(ipc_pipe_preprocess_in, keep_fd, keep_fd2)
        call close_pipe_except_fd(ipc_pipe_preprocess_out, keep_fd, keep_fd2)
        call close_pipe_except_fd(ipc_pipe_assign_optics_in, keep_fd, keep_fd2)
        call close_pipe_except_fd(ipc_pipe_assign_optics_out, keep_fd, keep_fd2)
        call close_pipe_except_fd(ipc_pipe_initial_analysis_in, keep_fd, keep_fd2)
        call close_pipe_except_fd(ipc_pipe_initial_analysis_out, keep_fd, keep_fd2)
        call close_pipe_except_fd(ipc_pipe_refpick_in, keep_fd, keep_fd2)
        call close_pipe_except_fd(ipc_pipe_refpick_out, keep_fd, keep_fd2)
        call close_pipe_except_fd(ipc_pipe_sieve_cavgs_in, keep_fd, keep_fd2)
        call close_pipe_except_fd(ipc_pipe_sieve_cavgs_out, keep_fd, keep_fd2)
        call close_pipe_except_fd(ipc_pipe_pool2D_in, keep_fd, keep_fd2)
        call close_pipe_except_fd(ipc_pipe_pool2D_out, keep_fd, keep_fd2)
    end subroutine close_child_pipe_fds

    subroutine close_pipe_except_fd(pipe, keep_fd, keep_fd2)
        integer, intent(inout) :: pipe(2)
        integer, intent(in)    :: keep_fd
        integer, intent(in), optional :: keep_fd2

        if( pipe(1) >= 0 .and. pipe(1) /= keep_fd ) then
            if( .not.(present(keep_fd2) .and. pipe(1) == keep_fd2) ) call close_fd_silent(pipe(1))
        endif
        if( pipe(2) >= 0 .and. pipe(2) /= keep_fd ) then
            if( .not.(present(keep_fd2) .and. pipe(2) == keep_fd2) ) call close_fd_silent(pipe(2))
        endif
    end subroutine close_pipe_except_fd

    subroutine close_fd_silent(fd)
        integer, intent(inout) :: fd
            integer                :: rc_close

        if( fd < 0 ) return
            rc_close = c_close(fd)
        fd = -1
    end subroutine close_fd_silent

    subroutine exec_stream_p00_master( self, cline )
        class(stream_p00_master), intent(inout)    :: self
        class(cmdline),           intent(inout)    :: cline
        integer, parameter                        :: N_STREAM_PIPES = 6
        type(parameters)                           :: params
        type(cmdline)                              :: cline_preprocess, cline_assign_optics
        type(cmdline)                              :: cline_opening2D, cline_reference_picking
        type(cmdline)                              :: cline_particle_sieving, cline_pool2D
        type(http_post)                            :: post
        type(http_response)                        :: response
        type(string)                               :: request, cwd
        type(json_core)                            :: json
        type(json_value),              pointer     :: json_response_ptr => null()
        type(gui_assembler)                        :: assembler
        type(qsys_env)                             :: qsys
        ! gui metadata
        type(gui_metadata_stream_update)             :: meta_update
        type(gui_metadata_stream_preprocess)         :: meta_preprocess
        type(gui_metadata_stream_optics_assignment)  :: meta_optics_assignment
        type(gui_metadata_stream_picking)            :: meta_initial_picking
        type(gui_metadata_stream_opening2D)          :: meta_opening2D
        type(gui_metadata_stream_picking)            :: meta_reference_picking
        type(gui_metadata_stream_particle_sieving)   :: meta_particle_sieving
        type(gui_metadata_stream_pool2D)             :: meta_pool2D
        type(gui_metadata_stream_pool2D_snapshot)    :: meta_pool2D_snapshot
        type(gui_metadata_micrograph),   allocatable :: meta_preprocess_micrographs(:)
        type(gui_metadata_histogram),    allocatable :: meta_preprocess_histograms(:)
        type(gui_metadata_timeplot),     allocatable :: meta_preprocess_timeplots(:)
        type(gui_metadata_optics_group), allocatable :: meta_optics_assignment_optics_groups(:)
        type(gui_metadata_micrograph),   allocatable :: meta_initial_picking_micrographs(:)
        type(gui_metadata_micrograph),   allocatable :: meta_reference_picking_micrographs(:)
        type(gui_metadata_cavg2D),       allocatable :: meta_opening2D_cavgs2D(:), meta_opening2D_final_cavgs2D(:)
        type(gui_metadata_cavg2D),       allocatable :: meta_reference_picking_cavgs2D(:), meta_pool2D_cavgs2D(:)
        type(gui_metadata_cavg2D),       allocatable :: meta_particle_sieving_cavgs2D(:), meta_particle_sieving_ref_cavgs2D(:)
        ! forked processes
        type(preprocess_fork)                      :: fork_preprocess
        type(assign_optics_fork)                   :: fork_assign_optics
        type(initial_analysis_fork)                :: fork_initial_analysis
        type(reference_picking_fork)               :: fork_reference_picking
        type(particle_sieving_fork)                :: fork_particle_sieving
        type(pool2D_fork)                          :: fork_pool2D
        type(c_pthread_t)                          :: meta_listener_thread
        type(c_ptr)                                :: ptr
        character(len=:),              allocatable :: meta_buffer
        character(kind=CK, len=:),     allocatable :: str_val
        integer,                       allocatable :: i_arr(:)
        type(json_value),              pointer     :: json_child_ptr
        logical                                    :: l_terminate=.false., l_last_loop=.false., l_found, l_test=.false., l_terminate_loop=.false.
        logical                                    :: got_snapshot_id, got_snapshot_iter, got_snapshot_sel, got_snapshot_file
        logical                                    :: l_existing_pickrefs, l_existing_box, l_existing_preprocess
        integer                                    :: stat, rc, max_msgsize, i_val, snapshot_id
        real(kind=dp)                              :: r_val
        type(pipe_rx_state)                        :: rx_state(N_STREAM_PIPES)
        ! check cline arguments 
        l_existing_pickrefs   = .false.
        l_existing_box        = .false.
        l_existing_preprocess = .false.
        if( cline%defined('pickrefs') )    l_existing_pickrefs = .true.
        if( cline%defined('box_extract') ) l_existing_box      = .true.
        ! init params
        call cline%printline()
        call params%new(cline)
        call simple_getcwd(cwd)
        ! check for existing preprocessing directory
        if( params%dir_preprocess%strlen() > 0 ) then
            if( .not.dir_exists(params%dir_preprocess)) THROW_HARD('Preprocessing directory '//params%dir_preprocess%to_char()//' does not exist.')
            l_existing_preprocess = .true.
        endif
        ! init assembler and http handler
        call post%new(params%niceserver)
        call assembler%new(params%niceprocid)
        ! init mutexes
        if( c_pthread_mutex_init(meta_mutex,      c_null_ptr) /= 0 ) THROW_HARD('failed to initialise metadata mutex' )
        if( c_pthread_mutex_init(terminate_mutex, c_null_ptr) /= 0 ) THROW_HARD('failed to initialise terminate mutex')
        ! start persistent worker server if requested by params
        params%qsys_name = '' ! force qsys_env to read from env vars so we can control with params
        params%ncunits   = 8  ! set to 8 for now to ensure enough threads for stream processes; can be overridden by env var or compenv
        call qsys%new(params, 1, qsys_nthr=16, stream=.true.)
        ! init update metadata
        call meta_update%new(GUI_METADATA_STREAM_UPDATE_TYPE)
        ! init metadata 
        call init_metadata_preprocess()
        call init_metadata_assign_optics()
        call init_metadata_opening2D()
        call init_metadata_reference_picking()
        call init_metadata_particle_sieving()
        call init_metadata_pool2D()
        ! init cmdlines
        call init_cline_preprocess()
        call init_cline_assign_optics()
        call init_cline_opening2D()
        call init_cline_reference_picking()
        call init_cline_particle_sieving()
        call init_cline_pool2D()
        ! create ipc pipes
        max_msgsize = max_metadata_size()
        write(logfhandle, *)"Max metadata size: ", max_msgsize
        call init_ipc_pipe(ipc_pipe_preprocess_in)
        call init_ipc_pipe(ipc_pipe_preprocess_out)
        call init_ipc_pipe(ipc_pipe_assign_optics_in)
        call init_ipc_pipe(ipc_pipe_assign_optics_out)
        call init_ipc_pipe(ipc_pipe_initial_analysis_in)
        call init_ipc_pipe(ipc_pipe_initial_analysis_out)
        call init_ipc_pipe(ipc_pipe_refpick_in)
        call init_ipc_pipe(ipc_pipe_refpick_out)
        call init_ipc_pipe(ipc_pipe_sieve_cavgs_in)
        call init_ipc_pipe(ipc_pipe_sieve_cavgs_out)
        call init_ipc_pipe(ipc_pipe_pool2D_in)
        call init_ipc_pipe(ipc_pipe_pool2D_out)
        ! spawn metadata listener thread
        stat = c_pthread_create(thread        = meta_listener_thread, &
                                attr          = c_null_ptr, &
                                start_routine = c_funloc(metadata_listener), &
                                arg           = c_null_ptr)
        if( stat /= 0 ) THROW_HARD('failed to create metadata listener thread')
        ! fork stream processes
        if( l_existing_preprocess ) then
            rc = symlink(params%dir_preprocess%to_char()//achar(0), PREPROC_JOB_NAME)
            if( rc /= 0 ) THROW_HARD('failed to create symlink for existing preprocessing directory')
            call fork_preprocess%skip()
        else
            call fork_preprocess%start(name=string(PREPROC_JOB_NAME), logfile=string(PREPROC_JOB_NAME//'.log'), cline=cline_preprocess, restart=.false.)
        endif
        call fork_assign_optics%start(    name=string(OPTICS_JOB_NAME),    logfile=string(OPTICS_JOB_NAME//'.log'),    cline=cline_assign_optics,    restart=.false.)
        call fork_reference_picking%start(name=string(REFPICK_JOB_NAME),   logfile=string(REFPICK_JOB_NAME//'.log'),   cline=cline_reference_picking,restart=.false.)
        call fork_particle_sieving%start( name=string(SIEVING_JOB_NAME),   logfile=string(SIEVING_JOB_NAME//'.log'),   cline=cline_particle_sieving, restart=.false.)
        call fork_pool2D%start(           name=string(CLASS2D_JOB_NAME),   logfile=string(CLASS2D_JOB_NAME//'.log'),   cline=cline_pool2D,           restart=.false.)
        if( l_existing_pickrefs ) then
            call fork_initial_analysis%skip()
        else
            call fork_initial_analysis%start(name=string(OPENING2D_JOB_NAME), logfile=string(OPENING2D_JOB_NAME//'.log'), cline=cline_opening2D, restart=.false.)
        endif
        ! test stream processes started successfully
        if( .not. l_existing_preprocess ) then
            if( fork_preprocess%status() /= FORK_STATUS_RUNNING ) THROW_HARD('failed to fork preprocessing'    )
        endif
        if( fork_assign_optics%status()     /= FORK_STATUS_RUNNING ) THROW_HARD('failed to fork assign optics'    )
        if( fork_reference_picking%status() /= FORK_STATUS_RUNNING ) THROW_HARD('failed to fork reference picking')
        if( fork_particle_sieving%status()  /= FORK_STATUS_RUNNING ) THROW_HARD('failed to fork particle sieving' )
        if( fork_pool2D%status()            /= FORK_STATUS_RUNNING ) THROW_HARD('failed to fork pool2D'           )
        if( .not. l_existing_pickrefs ) then
           if( fork_initial_analysis%status() /= FORK_STATUS_RUNNING ) THROW_HARD('failed to fork opening2D')
        endif
        ! attach signal handlers after fork else propagated to processes
        call signal(SIGTERM, sigterm_handler)
        call signal(SIGINT,   sigint_handler)
        ! main loop
        do while( .true. )
            ! heartbeat
            call assembler%assemble_stream_heartbeat(fork_preprocess, fork_assign_optics, fork_initial_analysis, fork_reference_picking, fork_particle_sieving, fork_pool2D)
            ! processes
            if( c_pthread_mutex_lock(meta_mutex) /= 0 ) THROW_HARD('failed to lock meta mutex')
            call assembler%assemble_stream_preprocess(meta_preprocess, meta_preprocess_micrographs, meta_preprocess_histograms, meta_preprocess_timeplots)
            call assembler%assemble_stream_optics_assignment(meta_optics_assignment, meta_optics_assignment_optics_groups)
            call assembler%assemble_stream_initial_picking(meta_initial_picking, meta_initial_picking_micrographs)
            call assembler%assemble_stream_opening2D(meta_opening2D, meta_opening2D_cavgs2D, meta_opening2D_final_cavgs2D)
            call assembler%assemble_stream_reference_picking(meta_reference_picking, meta_reference_picking_micrographs, meta_reference_picking_cavgs2D)
            call assembler%assemble_stream_particle_sieving(meta_particle_sieving, meta_particle_sieving_cavgs2D, meta_particle_sieving_ref_cavgs2D)
            call assembler%assemble_stream_pool2D(meta_pool2D, meta_pool2D_cavgs2D, meta_pool2D_snapshot)
            if( c_pthread_mutex_unlock(meta_mutex) /= 0 ) THROW_HARD('failed to unlock meta mutex')
            ! stringify assembled json
            request = assembler%to_string()
            write(*,*) request%to_char()
            ! send
            if( post%request(response, request) ) then
                write(*, *) "POST", request%to_char()
                if( response%code == 200) then
                    write(*, *) "RESPONSE", response%content%to_char()
                    ! parse response JSON
                    call json%parse(json_response_ptr, response%content%to_char())
                    if( json%failed()) then
                        write(logfhandle, '(A,A)') "FAILED TO PARSE JSON RESPONSE ", response%content%to_char()
                        call json%clear_exceptions()
                        nullify(json_response_ptr)
                    else
                        ! check for master process termination
                        call json%get(json_response_ptr, 'terminate', l_test, l_found)
                        if( l_found .and. l_test ) then
                            l_terminate_loop = .true.
                        endif
                        ! check for forked process termination
                        call json%get(json_response_ptr, 'terminate_preprocess', l_test, l_found)
                        if( l_found .and. l_test ) then
                            if( fork_preprocess%status() == FORK_STATUS_RUNNING ) call fork_preprocess%terminate()
                        endif
                        call json%get(json_response_ptr, 'terminate_optics_assignment', l_test, l_found)
                        if( l_found .and. l_test ) then
                            if( fork_assign_optics%status() == FORK_STATUS_RUNNING ) call fork_assign_optics%terminate()
                        endif
                        call json%get(json_response_ptr, 'terminate_opening2D', l_test, l_found)
                        if( l_found .and. l_test ) then
                            if( fork_initial_analysis%status() == FORK_STATUS_RUNNING ) call fork_initial_analysis%terminate()
                        endif
                        call json%get(json_response_ptr, 'terminate_reference_picking', l_test, l_found)
                        if( l_found .and. l_test ) then
                            if( fork_reference_picking%status() == FORK_STATUS_RUNNING ) call fork_reference_picking%terminate()
                        endif
                        call json%get(json_response_ptr, 'terminate_particle_sieving', l_test, l_found)
                        if( l_found .and. l_test ) then
                            if( fork_particle_sieving%status() == FORK_STATUS_RUNNING ) call fork_particle_sieving%terminate()
                        endif
                        call json%get(json_response_ptr, 'terminate_pool2D', l_test, l_found)
                        if( l_found .and. l_test ) then
                            if( fork_pool2D%status() == FORK_STATUS_RUNNING ) call fork_pool2D%terminate()
                        endif
                        ! check for forked process restart
                        call json%get(json_response_ptr, 'restart_preprocess', l_test, l_found)
                        if( l_found .and. l_test ) then
                            if( fork_preprocess%status() /= FORK_STATUS_RUNNING ) then
                                call fork_preprocess%start(name=string(PREPROC_JOB_NAME), logfile=string(PREPROC_JOB_NAME//'.log'), cline=cline_preprocess, restart=.true.)
                            endif
                        endif
                        call json%get(json_response_ptr, 'restart_optics_assignment', l_test, l_found)
                        if( l_found .and. l_test ) then
                            if( fork_assign_optics%status() /= FORK_STATUS_RUNNING ) then
                                call fork_assign_optics%start(name=string(OPTICS_JOB_NAME), logfile=string(OPTICS_JOB_NAME//'.log'),  cline=cline_assign_optics, restart=.true.)
                            endif
                        endif
                        call json%get(json_response_ptr, 'restart_opening2D', l_test, l_found)
                        if( l_found .and. l_test ) then
                            if( fork_initial_analysis%status() /= FORK_STATUS_RUNNING ) then
                                call fork_initial_analysis%start(name=string(OPENING2D_JOB_NAME), logfile=string(OPENING2D_JOB_NAME//'.log'),  cline=cline_opening2D, restart=.true.)
                            endif
                        endif
                        call json%get(json_response_ptr, 'restart_reference_picking', l_test, l_found)
                        if( l_found .and. l_test ) then
                            if( fork_reference_picking%status() /= FORK_STATUS_RUNNING ) then
                                call fork_reference_picking%start(name=string(REFPICK_JOB_NAME), logfile=string(REFPICK_JOB_NAME//'.log'),  cline=cline_reference_picking, restart=.true.)
                            endif
                        endif
                        call json%get(json_response_ptr, 'restart_particle_sieving', l_test, l_found)
                        if( l_found .and. l_test ) then
                            if( fork_particle_sieving%status() /= FORK_STATUS_RUNNING ) then
                                call fork_particle_sieving%start(name=string(SIEVING_JOB_NAME), logfile=string(SIEVING_JOB_NAME//'.log'),  cline=cline_particle_sieving, restart=.true.)
                            endif
                        endif
                        call json%get(json_response_ptr, 'restart_pool2D', l_test, l_found)
                        if( l_found .and. l_test ) then
                            if( fork_pool2D%status() /= FORK_STATUS_RUNNING ) then
                                call fork_pool2D%start(name=string(CLASS2D_JOB_NAME), logfile=string(CLASS2D_JOB_NAME//'.log'),  cline=cline_pool2D, restart=.true.)
                            endif
                        endif
                        ! gather update payload from HTTP response
                        call json%get(json_response_ptr, 'ctfresthreshold', r_val, l_found)
                        if(l_found) call meta_update%set_ctfres_update(real(r_val))
                        call json%get(json_response_ptr, 'astigthreshold', r_val, l_found)
                        if(l_found) call meta_update%set_astigmatism_update(real(r_val))
                        call json%get(json_response_ptr, 'icefracthreshold', r_val, l_found)
                        if(l_found) call meta_update%set_icescore_update(real(r_val))
                        call json%get(json_response_ptr, 'increase_nmics', i_val, l_found)
                        if(l_found) call meta_update%set_increase_nmics(i_val)
                        call json%get(json_response_ptr, 'pickrefs_selection', i_arr, l_found)
                        if(l_found) call meta_update%set_pickrefs_selection(i_arr)
                        call json%get(json_response_ptr, 'pickrefs_clusters', i_arr, l_found)
                        if(l_found) call meta_update%set_pickrefs_clusters(i_arr)
                        call json%get(json_response_ptr, 'ref_selection', i_arr, l_found)
                        if(l_found) call meta_update%set_sieverefs_selection(i_arr)
                        call json%get(json_response_ptr, 'mskdiam2D', r_val, l_found)
                        if(l_found) call meta_update%set_mskdiam2D_update(real(r_val))
                        nullify(json_child_ptr)
                        call json%get(json_response_ptr, 'snapshot2D', json_child_ptr, l_found)
                        if( l_found .and. associated(json_child_ptr) ) then
                            if( allocated(i_arr) )   deallocate(i_arr)
                            if( allocated(str_val) ) deallocate(str_val)
                            call json%get(json_child_ptr, 'id',        snapshot_id, got_snapshot_id)
                            call json%get(json_child_ptr, 'iteration', i_val,       got_snapshot_iter)
                            call json%get(json_child_ptr, 'selection', i_arr,       got_snapshot_sel)
                            call json%get(json_child_ptr, 'filename',  str_val,     got_snapshot_file)
                            if( got_snapshot_id .and. got_snapshot_iter .and. got_snapshot_sel .and. got_snapshot_file .and. &
                                & allocated(i_arr) .and. allocated(str_val) ) &
                                call meta_update%set_snapshot2D_update(snapshot_id, i_val, i_arr, string(str_val))
                            if( allocated(str_val) ) deallocate(str_val)
                            if( allocated(i_arr)   ) deallocate(i_arr)
                            nullify(json_child_ptr)
                        end if
                        if( meta_update%assigned() ) then
                            call meta_update%serialise(meta_buffer)
                            call send_update_to_stage_pipes(meta_buffer)
                        endif
                        if(allocated(i_arr)) deallocate(i_arr)
                    endif
                endif
            else
                ! clearing hashes ensure full json is sent on next attempt
                call assembler%clear_hashes()
            endif
            ! clean up
            if( associated(json_response_ptr) ) call json%destroy(json_response_ptr)
            call qsys%service_persistent_worker_warmup()
            ! exit if l_last_loop
            if( l_last_loop ) exit
            ! exit if terminate received all processes stopped
            if( l_terminate_loop ) then
                write(logfhandle, '(A)') "TERMINATE "
                ! send sigterm to all running forked processes
                if( fork_preprocess%status()        == FORK_STATUS_RUNNING ) call fork_preprocess%terminate()
                if( fork_assign_optics%status()     == FORK_STATUS_RUNNING ) call fork_assign_optics%terminate()
                if( fork_initial_analysis%status()  == FORK_STATUS_RUNNING ) call fork_initial_analysis%terminate()
                if( fork_reference_picking%status() == FORK_STATUS_RUNNING ) call fork_reference_picking%terminate()
                if( fork_particle_sieving%status()  == FORK_STATUS_RUNNING ) call fork_particle_sieving%terminate()
                if( fork_pool2D%status()            == FORK_STATUS_RUNNING ) call fork_pool2D%terminate()
                l_last_loop = .true.
                ! if processes are still running set last_loop back to false
                if( fork_preprocess%status()        == FORK_STATUS_RUNNING ) l_last_loop = .false.
                if( fork_assign_optics%status()     == FORK_STATUS_RUNNING ) l_last_loop = .false.
                if( fork_initial_analysis%status()  == FORK_STATUS_RUNNING ) l_last_loop = .false.
                if( fork_reference_picking%status() == FORK_STATUS_RUNNING ) l_last_loop = .false.
                if( fork_particle_sieving%status()  == FORK_STATUS_RUNNING ) l_last_loop = .false.
                if( fork_pool2D%status()            == FORK_STATUS_RUNNING ) l_last_loop = .false.
                ! set stoptime in assembler
                call assembler%set_stoptime()
            else
                 ! sleep for 5 seconds between requests
                call sleep(5)
            endif
            call flush(logfhandle)
        enddo
        ! set l_terminate to true -> terminates listener
        if( c_pthread_mutex_lock(terminate_mutex) /= 0 ) THROW_HARD('failed to lock terminate mutex')
        l_terminate = .true.
        if( c_pthread_mutex_unlock(terminate_mutex) /= 0 ) THROW_HARD('failed to unlock terminate mutex')
        ! join listener
        if( c_pthread_join(meta_listener_thread, ptr) /= 0) THROW_WARN('failed to join metadata listener thread' )
        ! cleanup
        call qsys%kill()
        call post%kill()
        call assembler%kill()
        ! close IPC pipes
        call kill_ipc_pipe(ipc_pipe_preprocess_in)
        call kill_ipc_pipe(ipc_pipe_preprocess_out)
        call kill_ipc_pipe(ipc_pipe_assign_optics_in)
        call kill_ipc_pipe(ipc_pipe_assign_optics_out)
        call kill_ipc_pipe(ipc_pipe_initial_analysis_in)
        call kill_ipc_pipe(ipc_pipe_initial_analysis_out)
        call kill_ipc_pipe(ipc_pipe_refpick_in)
        call kill_ipc_pipe(ipc_pipe_refpick_out)
        call kill_ipc_pipe(ipc_pipe_sieve_cavgs_in)
        call kill_ipc_pipe(ipc_pipe_sieve_cavgs_out)
        call kill_ipc_pipe(ipc_pipe_pool2D_in)
        call kill_ipc_pipe(ipc_pipe_pool2D_out)
        ! destroy mutexes
        if( c_pthread_mutex_destroy(meta_mutex)      /= 0) THROW_WARN('failed to destroy metadata mutex' )
        if( c_pthread_mutex_destroy(terminate_mutex) /= 0) THROW_WARN('failed to destroy terminate mutex')
        call flush(logfhandle)
        call exit(EXIT_SUCCESS)

    contains

        subroutine sigterm_handler()
            write(logfhandle, '(A)') 'SIGTERM RECEIVED (MASTER)'
            l_terminate_loop = .true.
        end subroutine sigterm_handler

        subroutine sigint_handler()
            integer :: my_rc
            write(logfhandle, '(A)') 'SIGINT RECEIVED (MASTER)'
            my_rc = c_pthread_mutex_lock(terminate_mutex)
            l_terminate = .true.
            my_rc = c_pthread_mutex_unlock(terminate_mutex)
            ! this is needed to prevent complete node death
            if( c_pthread_join(meta_listener_thread, ptr) /= 0) THROW_WARN('failed to join metadata listener thread')
            call exit(1)
            !l_terminate_loop = .true.
        end subroutine sigint_handler

        subroutine metadata_listener()
            character(len=:),   allocatable :: my_buffer
            type(gui_metadata_micrograph)   :: meta_mic_tmp
            type(gui_metadata_optics_group) :: meta_optics_group_tmp
            type(gui_metadata_cavg2D)       :: meta_cavg2D_tmp 
            integer                         :: my_rc, my_buffer_type, my_i
            logical                         :: my_l_continue, my_l_reinit
            my_l_continue = .true.
            do while( my_l_continue )
                my_rc = c_pthread_mutex_lock(meta_mutex)
                do while( read_any_pipe_message(my_buffer) )
                        if( allocated(my_buffer) ) then
                            my_buffer_type = transfer(my_buffer, my_buffer_type)
                            select case(my_buffer_type)
                                case(GUI_METADATA_STREAM_PREPROCESS_TYPE)
                                    meta_preprocess = transfer(my_buffer, meta_preprocess)
                                case( GUI_METADATA_STREAM_PREPROCESS_HISTOGRAM_ASTIG_TYPE )
                                    meta_preprocess_histograms(1) = transfer(my_buffer, meta_preprocess_histograms(1))
                                case( GUI_METADATA_STREAM_PREPROCESS_HISTOGRAM_CTFRES_TYPE )
                                    meta_preprocess_histograms(2) = transfer(my_buffer, meta_preprocess_histograms(2))
                                case( GUI_METADATA_STREAM_PREPROCESS_HISTOGRAM_ICEFRAC_TYPE )
                                    meta_preprocess_histograms(3) = transfer(my_buffer, meta_preprocess_histograms(3))  
                                case( GUI_METADATA_STREAM_PREPROCESS_TIMEPLOT_ASTIG_TYPE )
                                    meta_preprocess_timeplots(1) = transfer(my_buffer, meta_preprocess_timeplots(1)) 
                                case( GUI_METADATA_STREAM_PREPROCESS_TIMEPLOT_CTFRES_TYPE )
                                    meta_preprocess_timeplots(2) = transfer(my_buffer, meta_preprocess_timeplots(2)) 
                                case( GUI_METADATA_STREAM_PREPROCESS_TIMEPLOT_DF_TYPE )
                                    meta_preprocess_timeplots(3) = transfer(my_buffer, meta_preprocess_timeplots(3)) 
                                case( GUI_METADATA_STREAM_PREPROCESS_TIMEPLOT_RATE_TYPE )
                                    meta_preprocess_timeplots(4) = transfer(my_buffer, meta_preprocess_timeplots(4))
                                case( GUI_METADATA_STREAM_OPTICS_ASSIGNMENT_TYPE )
                                    meta_optics_assignment = transfer(my_buffer, meta_optics_assignment)
                                case( GUI_METADATA_STREAM_INITIAL_PICKING_TYPE )
                                    meta_initial_picking = transfer(my_buffer, meta_initial_picking)
                                case( GUI_METADATA_STREAM_OPENING2D_TYPE )
                                    meta_opening2D = transfer(my_buffer, meta_opening2D)
                                case( GUI_METADATA_STREAM_REFERENCE_PICKING_TYPE )
                                    meta_reference_picking = transfer(my_buffer, meta_reference_picking)
                                case( GUI_METADATA_STREAM_PARTICLE_SIEVING_TYPE )
                                    meta_particle_sieving = transfer(my_buffer, meta_particle_sieving)    
                                case( GUI_METADATA_STREAM_POOL2D_TYPE )
                                    meta_pool2D = transfer(my_buffer, meta_pool2D)  
                                case( GUI_METADATA_STREAM_POOL2D_SNAPSHOT_TYPE )
                                    meta_pool2D_snapshot = transfer(my_buffer, meta_pool2D_snapshot)         
                                case( GUI_METADATA_STREAM_PREPROCESS_MICROGRAPH_TYPE )
                                    my_l_reinit = .false.
                                    ! deserialise temporary copy of mic meta data
                                    meta_mic_tmp = transfer(my_buffer, meta_mic_tmp)
                                    ! allocate or resize meta_preprocess_micrographs as neccesary based on i_max
                                    if( .not.allocated(meta_preprocess_micrographs) ) then
                                        my_l_reinit = .true.
                                    else if( size(meta_preprocess_micrographs) /= meta_mic_tmp%get_i_max() ) then
                                        deallocate(meta_preprocess_micrographs)
                                        my_l_reinit = .true.
                                    endif
                                    if( my_l_reinit ) then
                                        ! allocate and initialise each object in meta_preprocess_micrographs
                                        allocate(meta_preprocess_micrographs(meta_mic_tmp%get_i_max()))
                                        do my_i=1, size(meta_preprocess_micrographs)
                                            call meta_preprocess_micrographs(my_i)%new(GUI_METADATA_STREAM_PREPROCESS_MICROGRAPH_TYPE)
                                            if( .not.meta_preprocess_micrographs(my_i)%initialized() ) THROW_HARD('failed to initialise preprocess micrograph metadata')
                                        enddo
                                    endif
                                    ! place the already-deserialised tmp object into the correct slot
                                    meta_preprocess_micrographs(meta_mic_tmp%get_i()) = meta_mic_tmp
                                case( GUI_METADATA_STREAM_INITIAL_PICKING_MICROGRAPH_TYPE )
                                    write(*,*) 'GUI_METADATA_STREAM_PICKING_MICROGRAPH_TYPE'
                                    my_l_reinit = .false.
                                    ! deserialise temporary copy of mic meta data
                                    meta_mic_tmp = transfer(my_buffer, meta_mic_tmp)
                                    ! allocate or resize meta_preprocess_micrographs as neccesary based on i_max
                                    if( .not.allocated(meta_initial_picking_micrographs) ) then
                                        my_l_reinit = .true.
                                    else if( size(meta_initial_picking_micrographs) /= meta_mic_tmp%get_i_max() ) then
                                        deallocate(meta_initial_picking_micrographs)
                                        my_l_reinit = .true.
                                    endif
                                    if( my_l_reinit ) then
                                        ! allocate and initialise each object in meta_initial_picking_micrographs
                                        allocate(meta_initial_picking_micrographs(meta_mic_tmp%get_i_max()))
                                        do my_i=1, size(meta_initial_picking_micrographs)
                                            call meta_initial_picking_micrographs(my_i)%new(GUI_METADATA_STREAM_INITIAL_PICKING_MICROGRAPH_TYPE)
                                            if( .not.meta_initial_picking_micrographs(my_i)%initialized() ) THROW_HARD('failed to initialise initial picking micrograph metadata')
                                        enddo
                                    endif
                                    ! place the already-deserialised tmp object into the correct slot
                                    meta_initial_picking_micrographs(meta_mic_tmp%get_i()) = meta_mic_tmp
                                case( GUI_METADATA_STREAM_REFERENCE_PICKING_MICROGRAPH_TYPE )
                                    write(*,*) 'GUI_METADATA_STREAM_REFERENCE_PICKING_MICROGRAPH_TYPE'
                                    my_l_reinit = .false.
                                    ! deserialise temporary copy of mic meta data
                                    meta_mic_tmp = transfer(my_buffer, meta_mic_tmp)
                                    ! allocate or resize meta_preprocess_micrographs as neccesary based on i_max
                                    if( .not.allocated(meta_reference_picking_micrographs) ) then
                                        my_l_reinit = .true.
                                    else if( size(meta_reference_picking_micrographs) /= meta_mic_tmp%get_i_max() ) then
                                        deallocate(meta_reference_picking_micrographs)
                                        my_l_reinit = .true.
                                    endif
                                    if( my_l_reinit ) then
                                        ! allocate and initialise each object in meta_reference_picking_micrographs
                                        allocate(meta_reference_picking_micrographs(meta_mic_tmp%get_i_max()))
                                        do my_i=1, size(meta_reference_picking_micrographs)
                                            call meta_reference_picking_micrographs(my_i)%new(GUI_METADATA_STREAM_REFERENCE_PICKING_MICROGRAPH_TYPE)
                                            if( .not.meta_reference_picking_micrographs(my_i)%initialized() ) THROW_HARD('failed to initialise reference picking micrograph metadata')
                                        enddo
                                    endif
                                    ! place the already-deserialised tmp object into the correct slot
                                    meta_reference_picking_micrographs(meta_mic_tmp%get_i()) = meta_mic_tmp
                                case( GUI_METADATA_STREAM_OPTICS_ASSIGNMENT_OPTICS_GROUP_TYPE )
                                    my_l_reinit = .false.
                                    ! deserialise temporary copy of optics group
                                    meta_optics_group_tmp = transfer(my_buffer, meta_optics_group_tmp)
                                    ! allocate or resize meta_optics_assignment_optics_groups as neccesary based on i_max
                                    if( .not.allocated(meta_optics_assignment_optics_groups) ) then
                                        my_l_reinit = .true.
                                    else if( size(meta_optics_assignment_optics_groups) /= meta_optics_group_tmp%get_i_max() ) then
                                        deallocate(meta_optics_assignment_optics_groups)
                                        my_l_reinit = .true.
                                    endif
                                    if( my_l_reinit ) then
                                        ! allocate and initialise each object in meta_optics_assignment_optics_groups
                                        allocate(meta_optics_assignment_optics_groups(meta_optics_group_tmp%get_i_max()))
                                        do my_i=1, size(meta_optics_assignment_optics_groups)
                                            call meta_optics_assignment_optics_groups(my_i)%new(GUI_METADATA_STREAM_OPTICS_ASSIGNMENT_OPTICS_GROUP_TYPE)
                                            if( .not.meta_optics_assignment_optics_groups(my_i)%initialized() ) THROW_HARD('failed to initialise optics assignment optics group metadata')
                                        enddo
                                    endif
                                    ! place the already-deserialised tmp object into the correct slot
                                    meta_optics_assignment_optics_groups(meta_optics_group_tmp%get_i()) = meta_optics_group_tmp
                                case( GUI_METADATA_STREAM_OPENING2D_CLS2D_TYPE )
                                    my_l_reinit = .false.
                                    ! deserialise temporary copy of cavg2D metadata to read routing fields
                                    meta_cavg2D_tmp = transfer(my_buffer, meta_cavg2D_tmp)
                                    ! allocate or resize meta_opening2D_cavgs2D as necessary based on i_max
                                    if( .not.allocated(meta_opening2D_cavgs2D) ) then
                                        my_l_reinit = .true.
                                    else if( size(meta_opening2D_cavgs2D) /= meta_cavg2D_tmp%get_i_max() ) then
                                        deallocate(meta_opening2D_cavgs2D)
                                        my_l_reinit = .true.
                                    endif
                                    if( my_l_reinit ) then
                                        ! allocate and initialise each object in meta_opening2D_cavgs2D
                                        allocate(meta_opening2D_cavgs2D(meta_cavg2D_tmp%get_i_max()))
                                        do my_i=1, size(meta_opening2D_cavgs2D)
                                            call meta_opening2D_cavgs2D(my_i)%new(GUI_METADATA_STREAM_OPENING2D_CLS2D_TYPE)
                                            if( .not.meta_opening2D_cavgs2D(my_i)%initialized() ) THROW_HARD('failed to initialise opening2D cavg2D metadata')
                                        enddo
                                    endif
                                    ! place the already-deserialised tmp object into the correct slot
                                    meta_opening2D_cavgs2D(meta_cavg2D_tmp%get_i()) = meta_cavg2D_tmp
                                case( GUI_METADATA_STREAM_OPENING2D_CLS2D_FINAL_TYPE )
                                    my_l_reinit = .false.
                                    ! deserialise temporary copy of cavg2D metadata to read routing fields
                                    meta_cavg2D_tmp = transfer(my_buffer, meta_cavg2D_tmp)
                                    ! allocate or resize meta_opening2D_final_cavgs2D as necessary based on i_max
                                    if( .not.allocated(meta_opening2D_final_cavgs2D) ) then
                                        my_l_reinit = .true.
                                    else if( size(meta_opening2D_final_cavgs2D) /= meta_cavg2D_tmp%get_i_max() ) then
                                        deallocate(meta_opening2D_final_cavgs2D)
                                        my_l_reinit = .true.
                                    endif
                                    if( my_l_reinit ) then
                                        ! allocate and initialise each object in meta_opening2D_final_cavgs2D
                                        allocate(meta_opening2D_final_cavgs2D(meta_cavg2D_tmp%get_i_max()))
                                        do my_i=1, size(meta_opening2D_final_cavgs2D)
                                            call meta_opening2D_final_cavgs2D(my_i)%new(GUI_METADATA_STREAM_OPENING2D_CLS2D_TYPE)
                                            if( .not.meta_opening2D_final_cavgs2D(my_i)%initialized() ) THROW_HARD('failed to initialise opening2D cavg2D metadata')
                                        enddo
                                    endif
                                    ! place the already-deserialised tmp object into the correct slot
                                    meta_opening2D_final_cavgs2D(meta_cavg2D_tmp%get_i()) = meta_cavg2D_tmp  
                                case( GUI_METADATA_STREAM_REFERENCE_PICKING_CLS2D_TYPE )
                                    my_l_reinit = .false.
                                    ! deserialise temporary copy of cavg2D metadata to read routing fields
                                    meta_cavg2D_tmp = transfer(my_buffer, meta_cavg2D_tmp)
                                    ! allocate or resize meta_reference_picking_cavgs2D as necessary based on i_max
                                    if( .not.allocated(meta_reference_picking_cavgs2D) ) then
                                        my_l_reinit = .true.
                                    else if( size(meta_reference_picking_cavgs2D) /= meta_cavg2D_tmp%get_i_max() ) then
                                        deallocate(meta_reference_picking_cavgs2D)
                                        my_l_reinit = .true.
                                    endif
                                    if( my_l_reinit ) then
                                        ! allocate and initialise each object in meta_reference_picking_cavgs2D
                                        allocate(meta_reference_picking_cavgs2D(meta_cavg2D_tmp%get_i_max()))
                                        do my_i=1, size(meta_reference_picking_cavgs2D)
                                            call meta_reference_picking_cavgs2D(my_i)%new(GUI_METADATA_STREAM_REFERENCE_PICKING_CLS2D_TYPE)
                                            if( .not.meta_reference_picking_cavgs2D(my_i)%initialized() ) THROW_HARD('failed to initialise opening2D cavg2D metadata')
                                        enddo
                                    endif
                                    ! place the already-deserialised tmp object into the correct slot
                                    meta_reference_picking_cavgs2D(meta_cavg2D_tmp%get_i()) = meta_cavg2D_tmp

                                case( GUI_METADATA_STREAM_PARTICLE_SIEVING_CLS2D_TYPE )
                                    my_l_reinit = .false.
                                    ! deserialise temporary copy of cavg2D metadata to read routing fields
                                    meta_cavg2D_tmp = transfer(my_buffer, meta_cavg2D_tmp)
                                    ! allocate or resize meta_particle_sieving_cavgs2D as necessary based on i_max
                                    if( .not.allocated(meta_particle_sieving_cavgs2D) ) then
                                        my_l_reinit = .true.
                                    else if( size(meta_particle_sieving_cavgs2D) /= meta_cavg2D_tmp%get_i_max() ) then
                                        deallocate(meta_particle_sieving_cavgs2D)
                                        my_l_reinit = .true.
                                    endif
                                    if( my_l_reinit ) then
                                        ! allocate and initialise each object in meta_particle_sieving_cavgs2D
                                        allocate(meta_particle_sieving_cavgs2D(meta_cavg2D_tmp%get_i_max()))
                                        do my_i=1, size(meta_particle_sieving_cavgs2D)
                                            call meta_particle_sieving_cavgs2D(my_i)%new(GUI_METADATA_STREAM_PARTICLE_SIEVING_CLS2D_TYPE)
                                            if( .not.meta_particle_sieving_cavgs2D(my_i)%initialized() ) THROW_HARD('failed to initialise particle sieving cavg2D metadata')
                                        enddo
                                    endif
                                    ! place the already-deserialised tmp object into the correct slot
                                    meta_particle_sieving_cavgs2D(meta_cavg2D_tmp%get_i()) = meta_cavg2D_tmp
                                case( GUI_METADATA_STREAM_POOL2D_CLS2D_TYPE )
                                    my_l_reinit = .false.
                                    ! deserialise temporary copy of cavg2D metadata to read routing fields
                                    meta_cavg2D_tmp = transfer(my_buffer, meta_cavg2D_tmp)
                                    ! allocate or resize meta_pool2D_cavgs2D as necessary based on i_max
                                    if( .not.allocated(meta_pool2D_cavgs2D) ) then
                                        my_l_reinit = .true.
                                    else if( size(meta_pool2D_cavgs2D) /= meta_cavg2D_tmp%get_i_max() ) then
                                        deallocate(meta_pool2D_cavgs2D)
                                        my_l_reinit = .true.
                                    endif
                                    if( my_l_reinit ) then
                                        ! allocate and initialise each object in meta_pool2D_cavgs2D
                                        allocate(meta_pool2D_cavgs2D(meta_cavg2D_tmp%get_i_max()))
                                        do my_i=1, size(meta_pool2D_cavgs2D)
                                            call meta_pool2D_cavgs2D(my_i)%new(GUI_METADATA_STREAM_POOL2D_CLS2D_TYPE)
                                            if( .not.meta_pool2D_cavgs2D(my_i)%initialized() ) THROW_HARD('failed to initialise particle sieving ref cavg2D metadata')
                                        enddo
                                    endif
                                    ! place the already-deserialised tmp object into the correct slot
                                    meta_pool2D_cavgs2D(meta_cavg2D_tmp%get_i()) = meta_cavg2D_tmp
                            end select
                            deallocate(my_buffer)
                        end if
                    end do
                my_rc = c_pthread_mutex_unlock(meta_mutex)
                if( c_pthread_mutex_lock(terminate_mutex) /= 0 ) THROW_HARD('failed to lock terminate mutex')
                if( l_terminate ) my_l_continue = .false.
                if( c_pthread_mutex_unlock(terminate_mutex) /= 0 ) THROW_HARD('failed to unlock terminate mutex')
                ! sleep briefly to avoid busy-polling when queue is idle
                rc = c_usleep(10000)
            end do
        end subroutine metadata_listener

        logical function read_any_pipe_message(buffer)
            character(len=:), allocatable, intent(inout) :: buffer
            character(kind=c_char), target               :: raw(max_msgsize)
            integer                                      :: ipipe
            integer                                      :: pipe_fds(N_STREAM_PIPES)

            read_any_pipe_message = .false.
            if( allocated(buffer) ) deallocate(buffer)

            pipe_fds = [ipc_pipe_preprocess_in(1),    &
                        ipc_pipe_assign_optics_in(1), &
                        ipc_pipe_initial_analysis_in(1), &
                        ipc_pipe_refpick_in(1),       &
                        ipc_pipe_sieve_cavgs_in(1),   &
                        ipc_pipe_pool2D_in(1)]

            ! First, emit any fully assembled frame already buffered.
            do ipipe = 1, N_STREAM_PIPES
                call try_extract_framed_message(rx_state(ipipe), buffer, read_any_pipe_message)
                if( read_any_pipe_message ) return
            end do

            do ipipe = 1, N_STREAM_PIPES
                call try_read_from_fd(pipe_fds(ipipe), raw, rx_state(ipipe), buffer, read_any_pipe_message)
                if( read_any_pipe_message ) return
            end do

        end function read_any_pipe_message

        subroutine try_read_from_fd(fd, raw, state, buffer, got_message)
            integer, intent(in)                           :: fd
            character(kind=c_char), target, intent(inout) :: raw(:)
            type(pipe_rx_state), intent(inout)            :: state
            character(len=:), allocatable, intent(inout)  :: buffer
            logical, intent(inout)                        :: got_message
            integer                                       :: nread, i
            character(len=:), allocatable                 :: chunk

            if( fd < 0 .or. got_message ) return
            nread = c_read(fd, c_loc(raw(1)), int(size(raw), c_size_t))
            if( nread > 0 ) then
                allocate(character(len=nread) :: chunk)
                do i = 1, nread
                    chunk(i:i) = transfer(raw(i), 'a')
                end do
                call append_pending(state, chunk)
            endif
            call try_extract_framed_message(state, buffer, got_message)
        end subroutine try_read_from_fd

        subroutine append_pending(state, chunk)
            type(pipe_rx_state), intent(inout) :: state
            character(len=*), intent(in)       :: chunk

            if( len(chunk) <= 0 ) return
            if( allocated(state%pending) ) then
                state%pending = state%pending // chunk
            else
                allocate(character(len=len(chunk)) :: state%pending)
                state%pending = chunk
            endif
        end subroutine append_pending

        subroutine try_extract_framed_message(state, buffer, got_message)
            type(pipe_rx_state), intent(inout)           :: state
            character(len=:), allocatable, intent(inout) :: buffer
            logical, intent(inout)                       :: got_message
            integer(c_int), target                       :: msg_len_c
            integer                                      :: header_bytes

            if( got_message ) return
            header_bytes = sizeof(msg_len_c)

            if( state%expected_len < 0 ) then
                if( .not. allocated(state%pending) ) return
                if( len(state%pending) < header_bytes ) return
                msg_len_c = transfer(state%pending(1:header_bytes), msg_len_c)
                state%expected_len = int(msg_len_c)
                if( state%expected_len <= 0 .or. state%expected_len > max_msgsize ) then
                    THROW_HARD('invalid framed metadata length read from stream pipe')
                endif
                if( len(state%pending) == header_bytes ) then
                    deallocate(state%pending)
                else
                    state%pending = state%pending(header_bytes + 1:)
                endif
            endif

            if( .not. allocated(state%pending) ) return
            if( len(state%pending) < state%expected_len ) return

            if( allocated(buffer) ) deallocate(buffer)
            allocate(character(len=state%expected_len) :: buffer)
            buffer = state%pending(1:state%expected_len)

            if( len(state%pending) == state%expected_len ) then
                deallocate(state%pending)
            else
                state%pending = state%pending(state%expected_len + 1:)
            endif
            state%expected_len = -1
            got_message = .true.
        end subroutine try_extract_framed_message

        subroutine init_cline_preprocess()
            cline_preprocess = cline
            call cline_preprocess%set('prg',                           'preproc')
            call cline_preprocess%set('projfile', PREPROC_JOB_NAME//METADATA_EXT)
            call cline_preprocess%set('outdir',                 PREPROC_JOB_NAME)
            call cline_preprocess%set('ninipick',               PREPROC_NINIPICK)
            call cline_preprocess%set('nparts',                                8)
            call cline_preprocess%set('nthr',                                  8)
            call cline_preprocess%set('mkdir',                             'yes')
            call cline_preprocess%delete( 'niceserver')
            call cline_preprocess%delete( 'niceprocid')
            call cline_preprocess%delete('box_extract')
            call cline_preprocess%delete(   'pickrefs')
        end subroutine init_cline_preprocess

        subroutine init_metadata_preprocess()
            ! preprocess
            call meta_preprocess%new(GUI_METADATA_STREAM_PREPROCESS_TYPE)
            if( .not.meta_preprocess%initialized() ) THROW_HARD('failed to initialise preprocess metadata')
            ! preprocess micrographs - allocated on receive
            ! preprocess_histograms
            allocate(meta_preprocess_histograms(3))
            call meta_preprocess_histograms(1)%new(GUI_METADATA_STREAM_PREPROCESS_HISTOGRAM_ASTIG_TYPE)
            if( .not.meta_preprocess_histograms(1)%initialized() ) THROW_HARD('failed to initialise preprocess astigmatism histogram metadata')
            call meta_preprocess_histograms(2)%new(GUI_METADATA_STREAM_PREPROCESS_HISTOGRAM_CTFRES_TYPE)
            if( .not.meta_preprocess_histograms(2)%initialized() ) THROW_HARD('failed to initialise preprocess ctfres histogram metadata')
            call meta_preprocess_histograms(3)%new(GUI_METADATA_STREAM_PREPROCESS_HISTOGRAM_ICEFRAC_TYPE)
            if( .not.meta_preprocess_histograms(3)%initialized() ) THROW_HARD('failed to initialise preprocess icefrac histogram metadata')
            ! preprocess timeplots
            allocate(meta_preprocess_timeplots(4))
            call meta_preprocess_timeplots(1)%new(GUI_METADATA_STREAM_PREPROCESS_TIMEPLOT_ASTIG_TYPE)
            if( .not.meta_preprocess_timeplots(1)%initialized() ) THROW_HARD('failed to initialise preprocess astigmatism timeplot metadata')
            call meta_preprocess_timeplots(2)%new(GUI_METADATA_STREAM_PREPROCESS_TIMEPLOT_CTFRES_TYPE)
            if( .not.meta_preprocess_timeplots(2)%initialized() ) THROW_HARD('failed to initialise preprocess ctfres timeplot metadata')
            call meta_preprocess_timeplots(3)%new(GUI_METADATA_STREAM_PREPROCESS_TIMEPLOT_DF_TYPE)
            if( .not.meta_preprocess_timeplots(3)%initialized() ) THROW_HARD('failed to initialise preprocess df timeplot metadata')
            call meta_preprocess_timeplots(4)%new(GUI_METADATA_STREAM_PREPROCESS_TIMEPLOT_RATE_TYPE)
            if( .not.meta_preprocess_timeplots(4)%initialized() ) THROW_HARD('failed to initialise preprocess rate timeplot metadata')
        end subroutine init_metadata_preprocess

        subroutine init_cline_assign_optics()
            call cline_assign_optics%set('prg',                    'assign_optics')
            call cline_assign_optics%set('projfile', OPTICS_JOB_NAME//METADATA_EXT)
            call cline_assign_optics%set('outdir',                 OPTICS_JOB_NAME)
            call cline_assign_optics%set('dir_target',            PREPROC_JOB_NAME)
            call cline_assign_optics%set('nthr',                                 1)
            call cline_assign_optics%set('mkdir',                            'yes')
        end subroutine init_cline_assign_optics

        subroutine init_metadata_assign_optics()
            ! assign optics
            call meta_optics_assignment%new(GUI_METADATA_STREAM_OPTICS_ASSIGNMENT_TYPE)
            if( .not.meta_optics_assignment%initialized() ) THROW_HARD('failed to initialise optics assignment metadata')
        end subroutine init_metadata_assign_optics

        subroutine init_cline_opening2D()
            call cline_opening2D%set('prg',                        'gen_pickrefs')
            call cline_opening2D%set('projfile', OPENING2D_JOB_NAME//METADATA_EXT)
            call cline_opening2D%set('outdir',                 OPENING2D_JOB_NAME)
            call cline_opening2D%set('dir_target',               PREPROC_JOB_NAME)
            call cline_opening2D%set('optics_dir',  cwd // '/' // OPTICS_JOB_NAME)
            call cline_opening2D%set('nthr',                                   32)
            call cline_opening2D%set('mkdir',                               'yes')
            call cline_opening2D%set('worker_priority',                    'high')
        end subroutine init_cline_opening2D

        subroutine init_metadata_opening2D()
            ! initial picking
            call meta_initial_picking%new(GUI_METADATA_STREAM_INITIAL_PICKING_TYPE)
            if( .not.meta_initial_picking%initialized() ) THROW_HARD('failed to initialise initial picking metadata')
            ! initial picking micrographs - allocated on receive
            ! opening2D
            call meta_opening2D%new(GUI_METADATA_STREAM_OPENING2D_TYPE)
            if( .not.meta_opening2D%initialized() ) THROW_HARD('failed to initialise opening2D metadata')
            ! opening2D cavgs2D - allocated on receive
        end subroutine init_metadata_opening2D

        subroutine init_cline_reference_picking()
            call cline_reference_picking%set('prg',                        'pick_extract')
            call cline_reference_picking%set('projfile',   REFPICK_JOB_NAME//METADATA_EXT)
            call cline_reference_picking%set('outdir',                   REFPICK_JOB_NAME)
            call cline_reference_picking%set('dir_target',               PREPROC_JOB_NAME)
            call cline_reference_picking%set('optics_dir',  cwd // '/' // OPTICS_JOB_NAME)
            call cline_reference_picking%set('nthr',                                    8)
            call cline_reference_picking%set('mkdir',                               'yes')
            call cline_reference_picking%set('nparts',                                  8)
            if( l_existing_pickrefs ) then
                call cline_reference_picking%set('pickrefs',              params%pickrefs)
            else
                call cline_reference_picking%set('pickrefs',    '../'//OPENING2D_JOB_NAME//'/selected_references.mrcs') 
            end if
            if( l_existing_box ) call cline_reference_picking%set('box_extract', params%box_extract)
        end subroutine init_cline_reference_picking

        subroutine init_metadata_reference_picking()
            ! reference picking
            call meta_reference_picking%new(GUI_METADATA_STREAM_REFERENCE_PICKING_TYPE)
            if( .not.meta_reference_picking%initialized() ) THROW_HARD('failed to initialise reference picking metadata')
        end subroutine init_metadata_reference_picking

        subroutine init_cline_particle_sieving()
            type(string) :: server_address
            server_address = qsys%get_persistent_worker_server_address()
            call cline_particle_sieving%set('prg',                         'sieve_cavgs')
            call cline_particle_sieving%set('projfile',   SIEVING_JOB_NAME//METADATA_EXT)
            call cline_particle_sieving%set('outdir',                   SIEVING_JOB_NAME)
            call cline_particle_sieving%set('dir_target',               REFPICK_JOB_NAME)
            call cline_particle_sieving%set('optics_dir',  cwd // '/' // OPTICS_JOB_NAME)
            call cline_particle_sieving%set('nthr',                                   16)
            call cline_particle_sieving%set('mkdir',                               'yes')
          !  call cline_particle_sieving%set('nparts',                                  8)
            call cline_particle_sieving%set('nchunks',                                 4)
            call cline_particle_sieving%set('worker_priority',                    'high')
            if( l_existing_pickrefs ) then
                call cline_particle_sieving%set('pickrefs',              params%pickrefs)
            else
                call cline_particle_sieving%set('pickrefs',    '../'//OPENING2D_JOB_NAME//'/selected_references.mrcs') 
            end if
            if( server_address%strlen() > 0 ) call cline_particle_sieving%set('worker_server', server_address)
        end subroutine init_cline_particle_sieving

        subroutine init_metadata_particle_sieving()
            ! particle sieving
            call meta_particle_sieving%new(GUI_METADATA_STREAM_PARTICLE_SIEVING_TYPE)
            if( .not.meta_particle_sieving%initialized() ) THROW_HARD('failed to initialise particle sieving metadata')
        end subroutine init_metadata_particle_sieving

        subroutine init_cline_pool2D()
            type(string) :: server_address
            server_address = qsys%get_persistent_worker_server_address()
            call cline_pool2D%set('prg',                       'abinitio2D_stream')
            call cline_pool2D%set('projfile',       CLASS2D_JOB_NAME//METADATA_EXT)
            call cline_pool2D%set('outdir',                       CLASS2D_JOB_NAME)
            call cline_pool2D%set('dir_target',                   SIEVING_JOB_NAME)
            call cline_pool2D%set('optics_dir',      cwd // '/' // OPTICS_JOB_NAME)
            call cline_pool2D%set('projfile_optics', OPTICS_JOB_NAME//METADATA_EXT)
            call cline_pool2D%set('nthr',                                        8)
            call cline_pool2D%set('mkdir',                                   'yes')
            call cline_pool2D%set('nparts',                                      6)
            call cline_pool2D%set('ncls',                                      150)
            call cline_pool2D%set('nicedispid',                  params%nicedispid)
            call cline_pool2D%set('worker_priority',                        'high')
            if( server_address%strlen() > 0 ) call cline_pool2D%set('worker_server', server_address)
        end subroutine init_cline_pool2D

        subroutine init_metadata_pool2D()
            ! pool2D
            call meta_pool2D%new(GUI_METADATA_STREAM_POOL2D_TYPE)
            if( .not.meta_pool2D%initialized() ) THROW_HARD('failed to initialise pool2D metadata')
            call meta_pool2D_snapshot%new(GUI_METADATA_STREAM_POOL2D_SNAPSHOT_TYPE)
            if( .not.meta_pool2D_snapshot%initialized() ) THROW_HARD('failed to initialise pool2D snapshot metadata')
        end subroutine init_metadata_pool2D

        subroutine init_ipc_pipe( pipe )
            integer, intent(inout) :: pipe(2)
            integer :: rc_pipe, flags
            rc_pipe = c_pipe(pipe)
            if( rc_pipe /= 0 ) THROW_HARD('failed to create IPC pipe')
            ! Keep both ends open across forks; stage processes write to *_in(2),
            ! and master listener reads from *_in(1).
            flags = c_fcntl(pipe(1), F_GETFL, 0)
            if( flags < 0 ) THROW_HARD('failed to get IPC pipe read-end flags')
            rc_pipe = c_fcntl(pipe(1), F_SETFL, ior(flags, O_NONBLOCK))
            if( rc_pipe < 0 ) THROW_HARD('failed to set IPC pipe read-end non-blocking mode')
            flags = c_fcntl(pipe(2), F_GETFL, 0)
            if( flags < 0 ) THROW_HARD('failed to get IPC pipe write-end flags')
            rc_pipe = c_fcntl(pipe(2), F_SETFL, ior(flags, O_NONBLOCK))
            if( rc_pipe < 0 ) THROW_HARD('failed to set IPC pipe write-end non-blocking mode')
        end subroutine init_ipc_pipe

        subroutine kill_ipc_pipe( pipe )
            integer, intent(inout) :: pipe(2)
            integer :: rc_close

            if( pipe(1) >= 0 ) then
                rc_close = c_close(pipe(1))
                if( rc_close /= 0 ) THROW_WARN('failed to close IPC pipe read-end')
                pipe(1) = -1
            endif
            if( pipe(2) >= 0 ) then
                rc_close = c_close(pipe(2))
                if( rc_close /= 0 ) THROW_WARN('failed to close IPC pipe write-end')
                pipe(2) = -1
            endif
        end subroutine kill_ipc_pipe

        subroutine send_update_to_stage_pipes(buffer)
            character(len=*), intent(in) :: buffer
            call send_framed_to_pipe(ipc_pipe_preprocess_out(2), buffer)
            call send_framed_to_pipe(ipc_pipe_assign_optics_out(2), buffer)
            call send_framed_to_pipe(ipc_pipe_initial_analysis_out(2), buffer)
            call send_framed_to_pipe(ipc_pipe_refpick_out(2), buffer)
            call send_framed_to_pipe(ipc_pipe_sieve_cavgs_out(2), buffer)
            call send_framed_to_pipe(ipc_pipe_pool2D_out(2), buffer)
        end subroutine send_update_to_stage_pipes

        subroutine send_framed_to_pipe(fd, buffer)
            integer, intent(in)                          :: fd
            character(len=*), intent(in)                 :: buffer
            character(len=:), allocatable                :: framed
            character(kind=c_char), allocatable, target  :: cbuf(:)
            integer(c_int), target                       :: msg_len
            integer(c_int)                               :: nwritten
            integer                                      :: nbytes, header_bytes, framed_nbytes
            integer                                      :: sent, retry_count, rc_sleep, err_no, ich
            integer, parameter                           :: MAX_RETRIES = 200
            integer, parameter                           :: RETRY_SLEEP_US = 10000

            if( fd < 0 ) return
            nbytes = len(buffer)
            if( nbytes <= 0 ) return

            msg_len = int(nbytes, c_int)
            header_bytes = sizeof(msg_len)
            framed_nbytes = header_bytes + nbytes
            allocate(character(len=framed_nbytes) :: framed)
            framed(1:header_bytes) = transfer(msg_len, framed(1:header_bytes))
            framed(header_bytes + 1:) = buffer

            allocate(cbuf(framed_nbytes))
            do ich = 1, framed_nbytes
                cbuf(ich) = transfer(framed(ich:ich), cbuf(ich))
            end do

            sent = 0
            retry_count = 0
            do while( sent < framed_nbytes )
                nwritten = c_write(fd, c_loc(cbuf(sent + 1)), int(framed_nbytes - sent, c_size_t))
                if( nwritten > 0 ) then
                    sent = sent + int(nwritten)
                    retry_count = 0
                    cycle
                endif

                err_no = ierrno()
                if( err_no == int(EAGAIN) .or. err_no == int(EWOULDBLOCK) ) then
                    retry_count = retry_count + 1
                    if( retry_count > MAX_RETRIES ) then
                        THROW_WARN('failed to send framed update to stream pipe: retry limit exceeded')
                        exit
                    endif
                    rc_sleep = c_usleep(RETRY_SLEEP_US)
                    cycle
                endif

                if( err_no == int(EPIPE) .or. err_no == int(EBADF) ) then
                    THROW_WARN('failed to send framed update to stream pipe: no active reader')
                else
                    THROW_WARN('failed to send framed update to stream pipe')
                endif
                exit
            end do

            if( allocated(cbuf) ) deallocate(cbuf)
        end subroutine send_framed_to_pipe

    end subroutine exec_stream_p00_master

end module simple_stream_p00_master