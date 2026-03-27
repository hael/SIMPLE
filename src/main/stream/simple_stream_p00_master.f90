!@descr: task 00 in the stream pipeline: master controller used when running from GUI
module simple_stream_p00_master
use unix
use simple_stream_api
use simple_stream_mq_defs             
use simple_stream_p01_preprocess_new,    only: stream_p01_preprocess
use simple_stream_p02_assign_optics_new, only: stream_p02_assign_optics
use simple_stream_p03_opening2D_new,     only: stream_p03_opening2D
use simple_stream_p04_refpick_extract_new, only: stream_p04_refpick_extract
use simple_http_post,                    only: http_post, http_response
use simple_forked_process,               only: forked_process, FORK_STATUS_RUNNING
use simple_gui_metadata_api
use simple_gui_assembler,                only: gui_assembler
use simple_gui_metadata_utils,           only: max_metadata_size

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

type, extends(forked_process) :: opening2D_fork
    contains
    procedure :: execute => xopening2D
end type opening2D_fork

type, extends(forked_process) :: reference_picking_fork
    contains
    procedure :: execute => xreference_picking
end type reference_picking_fork

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
        call commander%execute(cline)
    end subroutine xpreprocess

    subroutine xassign_optics( self, cline )
        class(assign_optics_fork), intent(inout) :: self
        class(cmdline),            intent(inout) :: cline
        type(stream_p02_assign_optics)           :: commander
        call commander%execute(cline)
    end subroutine xassign_optics

    subroutine xopening2D( self, cline )
        class(opening2D_fork), intent(inout) :: self
        class(cmdline),        intent(inout) :: cline
        type(stream_p03_opening2D)           :: commander
        call commander%execute(cline)
    end subroutine xopening2D

    subroutine xreference_picking( self, cline )
        class(reference_picking_fork), intent(inout) :: self
        class(cmdline),                intent(inout) :: cline
        type(stream_p04_refpick_extract)             :: commander
        call commander%execute(cline)
    end subroutine xreference_picking

    subroutine exec_stream_p00_master( self, cline )
        class(stream_p00_master), intent(inout)    :: self
        class(cmdline),           intent(inout)    :: cline
        type(parameters)                           :: params
        type(cmdline)                              :: cline_preprocess, cline_assign_optics
        type(cmdline)                              :: cline_opening2D, cline_reference_picking
        type(http_post)                            :: post
        type(http_response)                        :: response
        type(string)                               :: request
        type(json_core)                            :: json
        type(json_value),              pointer     :: json_response_ptr => null()
        type(gui_assembler)                        :: assembler
        ! gui metadata
        type(gui_metadata_stream_update)             :: meta_update
        type(gui_metadata_stream_preprocess)         :: meta_preprocess
        type(gui_metadata_stream_optics_assignment)  :: meta_optics_assignment
        type(gui_metadata_stream_picking)            :: meta_initial_picking
        type(gui_metadata_stream_opening2D)          :: meta_opening2D
        type(gui_metadata_stream_picking)            :: meta_reference_picking
        type(gui_metadata_micrograph),   allocatable :: meta_preprocess_micrographs(:)
        type(gui_metadata_histogram),    allocatable :: meta_preprocess_histograms(:)
        type(gui_metadata_timeplot),     allocatable :: meta_preprocess_timeplots(:)
        type(gui_metadata_optics_group), allocatable :: meta_optics_assignment_optics_groups(:)
        type(gui_metadata_micrograph),   allocatable :: meta_initial_picking_micrographs(:)
        type(gui_metadata_micrograph),   allocatable :: meta_reference_picking_micrographs(:)
        type(gui_metadata_cavg2D),       allocatable :: meta_opening2D_cavgs2D(:), meta_opening2D_final_cavgs2D(:)
        type(gui_metadata_cavg2D),       allocatable :: meta_reference_picking_cavgs2D(:)
        ! forked processes
        type(preprocess_fork)                      :: fork_preprocess
        type(assign_optics_fork)                   :: fork_assign_optics
        type(opening2D_fork)                       :: fork_opening2D
        type(reference_picking_fork)               :: fork_reference_picking
        type(c_pthread_t)                          :: meta_listener_thread
        type(c_ptr)                                :: ptr
        character(len=:),              allocatable :: meta_buffer
        integer,                       allocatable :: i_arr(:)
        logical                                    :: l_terminate=.false., l_last_loop=.false., l_found, l_test=.false., l_terminate_loop=.false.
        integer                                    :: stat, i, rc, max_msgsize, i_val
        real(kind=dp)                              :: r_val
        ! init params
        call params%new(cline)
        ! init assembler and http handler
        call post%new(params%niceserver)
        call assembler%new(params%niceprocid)
        ! init mutexes
        if( c_pthread_mutex_init(meta_mutex,      c_null_ptr) /= 0 ) THROW_HARD('failed to initialise metadata mutex' )
        if( c_pthread_mutex_init(terminate_mutex, c_null_ptr) /= 0 ) THROW_HARD('failed to initialise terminate mutex')
        ! init update metadata
        call meta_update%new(GUI_METADATA_STREAM_UPDATE_TYPE)
        ! init metadata 
        call init_metadata_preprocess()
        call init_metadata_assign_optics()
        call init_metadata_opening2D()
        call init_metadata_reference_picking()
        ! init cmdlines
        call init_cline_preprocess()
        call init_cline_assign_optics()
        call init_cline_opening2D()
        call init_cline_reference_picking()
        ! create message queues
        max_msgsize = max_metadata_size()
        call mq_stream_master_in%new(name=string('stream_master_in'), max_msgsize=max_msgsize)
        if( .not.mq_stream_master_in%is_active() ) THROW_HARD('failed to create stream_master_in message queue')
        call mq_stream_master_out%new(name=string('stream_master_out'), max_msgsize=max_msgsize)
        if( .not.mq_stream_master_out%is_active() ) THROW_HARD('failed to create stream_master_out message queue')
        ! spawn metadata listener thread
        stat = c_pthread_create(thread        = meta_listener_thread, &
                                attr          = c_null_ptr, &
                                start_routine = c_funloc(metadata_listener), &
                                arg           = c_null_ptr)
        ! fork and test stream processes
        call fork_preprocess%start(name=string(PREPROC_JOB_NAME),   logfile=string(PREPROC_JOB_NAME//'.log'),   cline=cline_preprocess,    restart=.true.)
        call fork_assign_optics%start(name=string(OPTICS_JOB_NAME), logfile=string(OPTICS_JOB_NAME//'.log'),    cline=cline_assign_optics, restart=.true.)
        call fork_opening2D%start(name=string(OPENING2D_JOB_NAME),  logfile=string(OPENING2D_JOB_NAME//'.log'), cline=cline_opening2D,     restart=.true.)
        call fork_reference_picking%start(name=string(REFPICK_JOB_NAME),  logfile=string(REFPICK_JOB_NAME//'.log'), cline=cline_reference_picking, restart=.true.)
        if( fork_preprocess%status()    /= FORK_STATUS_RUNNING ) THROW_HARD('failed to fork preprocessing')
        if( fork_assign_optics%status() /= FORK_STATUS_RUNNING ) THROW_HARD('failed to fork assign optics')
        if( fork_opening2D%status()     /= FORK_STATUS_RUNNING ) THROW_HARD('failed to fork opening2D'    )
        if( fork_reference_picking%status()     /= FORK_STATUS_RUNNING ) THROW_HARD('failed to fork reference picking'    )
        ! attach signal handlers after fork else propagated to processes
        call signal(SIGTERM, sigterm_handler)
        call signal(SIGINT,   sigint_handler)
        ! main loop
        do while( .true. )
            ! heartbeat
            call assembler%assemble_stream_heartbeat(fork_preprocess, fork_assign_optics, fork_opening2D, fork_reference_picking)
            ! processes
            if( c_pthread_mutex_lock(meta_mutex) /= 0 ) THROW_HARD('failed to lock meta mutex')
            call assembler%assemble_stream_preprocess(meta_preprocess, meta_preprocess_micrographs, meta_preprocess_histograms, meta_preprocess_timeplots)
            call assembler%assemble_stream_optics_assignment(meta_optics_assignment, meta_optics_assignment_optics_groups)
            call assembler%assemble_stream_initial_picking(meta_initial_picking, meta_initial_picking_micrographs)
            call assembler%assemble_stream_opening2D(meta_opening2D, meta_opening2D_cavgs2D, meta_opening2D_final_cavgs2D)
            call assembler%assemble_stream_reference_picking(meta_reference_picking, meta_reference_picking_micrographs, meta_reference_picking_cavgs2D)
            if( c_pthread_mutex_unlock(meta_mutex) /= 0 ) THROW_HARD('failed to unlock meta mutex')
            ! stringify assembled json
            request = assembler%to_string()
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
                            if( fork_opening2D%status() == FORK_STATUS_RUNNING ) call fork_opening2D%terminate()
                        endif
                        call json%get(json_response_ptr, 'terminate_reference_picking', l_test, l_found)
                        if( l_found .and. l_test ) then
                            if( fork_reference_picking%status() == FORK_STATUS_RUNNING ) call fork_reference_picking%terminate()
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
                            if( fork_opening2D%status() /= FORK_STATUS_RUNNING ) then
                                call fork_opening2D%start(name=string(OPENING2D_JOB_NAME), logfile=string(OPENING2D_JOB_NAME//'.log'),  cline=cline_opening2D, restart=.true.)
                            endif
                        endif
                        call json%get(json_response_ptr, 'restart_reference_picking', l_test, l_found)
                        if( l_found .and. l_test ) then
                            if( fork_reference_picking%status() /= FORK_STATUS_RUNNING ) then
                                call fork_reference_picking%start(name=string(REFPICK_JOB_NAME), logfile=string(REFPICK_JOB_NAME//'.log'),  cline=cline_reference_picking, restart=.true.)
                            endif
                        endif
                        ! wait to get update message from outbound queue and destroy
                        if( mq_stream_master_out%is_active() ) then
                            if( mq_stream_master_out%receive_timed(meta_buffer, 1) ) then
                                write(logfhandle, *) "GOT EXISTING UPDATE"
                                if( allocated(meta_buffer) ) deallocate(meta_buffer)
                            endif
                        endif
                        ! add update message from outbound message queue 
                        call json%get(json_response_ptr, 'ctfresthreshold', r_val, l_found)
                        if(l_found) call meta_update%set_ctfres_update(real(r_val))
                        call json%get(json_response_ptr, 'astigmatism', r_val, l_found)
                        if(l_found) call meta_update%set_astigmatism_update(real(r_val))
                        call json%get(json_response_ptr, 'icescore', r_val, l_found)
                        if(l_found) call meta_update%set_icescore_update(real(r_val))
                        call json%get(json_response_ptr, 'increase_nmics', i_val, l_found)
                        if(l_found) call meta_update%set_increase_nmics(i_val)
                        call json%get(json_response_ptr, 'pickrefs_selection', i_arr, l_found)
                        if(l_found) call meta_update%set_pickrefs_selection(i_arr)
                        if( meta_update%assigned() .and. mq_stream_master_out%is_active() ) then
                            write(logfhandle, *) "SEND UPDATE"
                            call meta_update%serialise(meta_buffer)
                            call mq_stream_master_out%send(meta_buffer)
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
            ! exit if l_last_loop
            if( l_last_loop ) exit
            ! exit if terminate received all processes stopped
            if( l_terminate_loop ) then
                write(logfhandle, '(A)') "TERMINATE "
                ! send sigterm to all running forked processes
                if( fork_preprocess%status()    == FORK_STATUS_RUNNING ) call fork_preprocess%terminate()
                if( fork_assign_optics%status() == FORK_STATUS_RUNNING ) call fork_assign_optics%terminate()
                if( fork_opening2D%status()     == FORK_STATUS_RUNNING ) call fork_opening2D%terminate()
                if( fork_reference_picking%status() == FORK_STATUS_RUNNING ) call fork_reference_picking%terminate()
                l_last_loop = .true.
                ! if processes are still running set last_loop back to false
                if( fork_preprocess%status()    == FORK_STATUS_RUNNING ) l_last_loop = .false.
                if( fork_assign_optics%status() == FORK_STATUS_RUNNING ) l_last_loop = .false.
                if( fork_opening2D%status()     == FORK_STATUS_RUNNING ) l_last_loop = .false.
                if( fork_reference_picking%status() == FORK_STATUS_RUNNING ) l_last_loop = .false.
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
        call post%kill()
        call assembler%kill()
        ! kill message queue
        if( mq_stream_master_in%is_active() )  call mq_stream_master_in%kill()
        if( mq_stream_master_out%is_active() ) call mq_stream_master_out%kill()
        ! destroy mutexes
        if( c_pthread_mutex_destroy(meta_mutex)      /= 0) THROW_WARN('failed to destroy metadata mutex' )
        if( c_pthread_mutex_destroy(terminate_mutex) /= 0) THROW_WARN('failed to destroy terminate mutex')
        call flush(logfhandle)
        call exit(EXIT_SUCCESS)

    contains

        subroutine sigterm_handler()
            integer :: my_rc
            write(logfhandle, '(A)') 'SIGTERM RECEIVED (MASTER)'
            l_terminate_loop = .true.
        end subroutine sigterm_handler

        subroutine sigint_handler()
            integer :: my_rc
            write(logfhandle, '(A)') 'SIGINT RECEIVED (MASTER)'
            my_rc = c_pthread_mutex_lock(terminate_mutex)
            l_terminate = .true.
            my_rc = c_pthread_mutex_unlock(terminate_mutex)
            ! join listener
            if( c_pthread_join(meta_listener_thread, ptr) /= 0) THROW_WARN('failed to join metadata listener thread' )
            call exit(1)
        end subroutine sigint_handler

        subroutine metadata_listener()
            type(gui_metadata_micrograph)   :: meta_mic_tmp
            type(gui_metadata_optics_group) :: meta_optics_group_tmp
            type(gui_metadata_cavg2D)       :: meta_cavg2D_tmp
            character(len=:), allocatable :: my_buffer
            integer                       :: my_rc, my_buffer_type, my_i
            logical                       :: my_l_continue, my_l_reinit
            my_l_continue = .true.
            do while( my_l_continue )
                my_rc = c_pthread_mutex_lock(meta_mutex)
                if( mq_stream_master_in%is_active() ) then
                    do while( mq_stream_master_in%receive(my_buffer) )
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
                            end select
                            deallocate(my_buffer)
                        end if
                    end do
                end if
                my_rc = c_pthread_mutex_unlock(meta_mutex)
                if( c_pthread_mutex_lock(terminate_mutex) /= 0 ) THROW_HARD('failed to lock terminate mutex')
                if( l_terminate ) my_l_continue = .false.
                if( c_pthread_mutex_unlock(terminate_mutex) /= 0 ) THROW_HARD('failed to unlock terminate mutex')
                ! sleep for 1us to save cpu cycles
                rc = c_usleep(1)
            end do
        end subroutine metadata_listener

        subroutine init_cline_preprocess()
            cline_preprocess = cline
            call cline_preprocess%set('prg',                           'preproc')
            call cline_preprocess%set('projfile', PREPROC_JOB_NAME//METADATA_EXT)
            call cline_preprocess%set('outdir',                 PREPROC_JOB_NAME)
            call cline_preprocess%set('ninipick',               PREPROC_NINIPICK)
            call cline_preprocess%set('nparts',                                2)
            call cline_preprocess%set('nthr',                                 16)
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
            call cline_opening2D%set('optics_dir',                OPTICS_JOB_NAME)
            call cline_opening2D%set('nthr',                                   32)
            call cline_opening2D%set('mkdir',                               'yes')
            call cline_opening2D%set('nmics',                                  20) ! for testing only
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
            call cline_reference_picking%set('optics_dir',                OPTICS_JOB_NAME)
            call cline_reference_picking%set('nthr',                                    1)
            call cline_reference_picking%set('mkdir',                               'yes')
            call cline_reference_picking%set('nparts',                                  5)
            call cline_reference_picking%set('pickrefs',    '../'//OPENING2D_JOB_NAME//'/selected_references.mrcs') 
        end subroutine init_cline_reference_picking

        subroutine init_metadata_reference_picking()
            ! initial picking
            call meta_reference_picking%new(GUI_METADATA_STREAM_REFERENCE_PICKING_TYPE)
            if( .not.meta_reference_picking%initialized() ) THROW_HARD('failed to initialise reference picking metadata')
        end subroutine init_metadata_reference_picking

    end subroutine exec_stream_p00_master

end module simple_stream_p00_master