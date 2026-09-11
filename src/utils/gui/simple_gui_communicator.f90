!@descr: communication between SIMPLE and the GUI
module simple_gui_communicator
use unix
use json_kinds
use json_module
use simple_sp_project,           only: sp_project
use simple_gui_assembler,        only: gui_assembler
use simple_http_post,            only: http_post, http_response
use simple_core_module_api
use simple_parameters,           only: parameters
use simple_gui_metadata_project, only: gui_metadata_project
use simple_gui_metadata_types,   only: GUI_METADATA_PROJECT_TYPE
implicit none

! payload handed to the worker thread through c_pthread_create's arg pointer
type, bind(c) :: gui_comm_args
    integer(c_int)           :: procid
    character(kind=c_char)   :: url(STDLEN)
    type(c_pthread_mutex_t)  :: terminate_mutex
    type(c_pthread_mutex_t)  :: metadata_mutex
    logical(c_bool)          :: l_terminate
    type(c_ptr)              :: metadata_project ! points at gui_project_metadata_inst
end type gui_comm_args

type(gui_comm_args),        target, save :: gui_comm_args_inst
type(gui_metadata_project), target, save :: gui_project_metadata_inst

public :: gui_communicator
private
#include "simple_local_flags.inc"

type gui_communicator
    private
    type(c_pthread_t) :: comm_thread
    logical           :: is_active     = .false.
    logical           :: remote_active = .false.
contains
    procedure :: new
    procedure :: kill => kill_gui_communicator
    procedure :: add_metadata_1
    procedure :: add_metadata_2
    generic   :: add_metadata => add_metadata_1, add_metadata_2
end type gui_communicator

contains

    subroutine new( self, params )
        class(gui_communicator), intent(inout) :: self
        type(parameters),        intent(inout) :: params
        integer(c_int)                         :: stat
        character(len=:), allocatable          :: niceserver_url
        integer                                :: url_len
        if( self%is_active ) return
        self%is_active = .true.
        ! pack url + procid for the worker thread, read back via c_f_pointer on arg
        niceserver_url                      = params%niceserver%to_char()
        url_len                             = min(len(niceserver_url), STDLEN - 1)
        gui_comm_args_inst%procid           = params%niceprocid
        gui_comm_args_inst%url              = c_null_char
        gui_comm_args_inst%l_terminate      = .false.
        gui_comm_args_inst%url(1:url_len)   = transfer(niceserver_url(1:url_len), gui_comm_args_inst%url(1:url_len))
        gui_comm_args_inst%metadata_project = c_loc(gui_project_metadata_inst)
        if( c_pthread_mutex_init(gui_comm_args_inst%terminate_mutex, c_null_ptr) /= 0 ) THROW_HARD('failed to initialise termination mutex')
        if( c_pthread_mutex_init(gui_comm_args_inst%metadata_mutex,  c_null_ptr) /= 0 ) THROW_HARD('failed to initialise metadata mutex')
        ! initialise the metadata
        call gui_project_metadata_inst%new(GUI_METADATA_PROJECT_TYPE)
        if( params%niceserver .ne. "" .and. params%niceprocid > 0 ) then
            ! spawn metadata listener thread
            stat = c_pthread_create(thread        = self%comm_thread, &
                                    attr          = c_null_ptr, &
                                    start_routine = c_funloc(communication_worker), &
                                    arg           = c_loc(gui_comm_args_inst))
            if( stat /= 0 ) THROW_HARD('failed to create metadata listener thread')
            ! clear the niceprocid and niceserver in the parameters to indicate they have been consumed
            params%niceprocid  = -1
            params%niceserver  = ""
            self%remote_active = .true.
        end if
    end subroutine new

    subroutine kill_gui_communicator( self )
        class(gui_communicator), intent(inout) :: self
        type(c_ptr)                            :: ptr
        type(gui_assembler)                    :: assembler
        type(string)                           :: metadata_backup
        if( .not. self%is_active) return
        if( c_pthread_mutex_lock(gui_comm_args_inst%terminate_mutex) /= 0    ) THROW_HARD('failed to lock terminate mutex')
        gui_comm_args_inst%l_terminate = .true.
        if( c_pthread_mutex_unlock(gui_comm_args_inst%terminate_mutex) /= 0  ) THROW_HARD('failed to unlock terminate mutex')
        if( self%remote_active) then
            if( c_pthread_join(self%comm_thread, ptr) /= 0 ) then
                THROW_WARN('failed to join communication worker thread')
                return ! thread may still be running, unsafe to destroy mutexes or kill shared metadata
            end if
        end if
        ! write metadata to file
        call assembler%new(gui_comm_args_inst%procid)
        call assembler%set_stoptime()
        call assembler%assemble_batch_heartbeat()
        if( c_pthread_mutex_lock(gui_comm_args_inst%metadata_mutex) /= 0 ) THROW_HARD('failed to lock metadata mutex')
        call assembler%assemble_batch_metadata(gui_project_metadata_inst)
        if( c_pthread_mutex_unlock(gui_comm_args_inst%metadata_mutex) /= 0 ) THROW_HARD('failed to unlock metadata mutex')
        metadata_backup = assembler%to_string()
        call write_singlelineoftext(string('metadata.json'), metadata_backup)
        if( c_pthread_mutex_destroy(gui_comm_args_inst%terminate_mutex) /= 0 ) THROW_WARN('failed to destroy terminate mutex')
        if( c_pthread_mutex_destroy(gui_comm_args_inst%metadata_mutex) /= 0  ) THROW_WARN('failed to destroy metadata mutex')
        call gui_project_metadata_inst%kill()
        self%is_active     = .false.
        self%remote_active = .false.
    end subroutine kill_gui_communicator

    subroutine add_metadata_1( self, spproj, oritype, stage2D, selection )
        class(gui_communicator),    intent(inout) :: self
        type(sp_project),           intent(inout) :: spproj
        character(len=*), optional, intent(in)    :: oritype
        integer,          optional, intent(in)    :: stage2D
        logical,          optional, intent(in)    :: selection
        character(len=SHORTSTRLEN)                :: md_oritype
        integer                                   :: i_stage2D
        logical                                   :: l_selection
        if( .not. self%is_active) return
        i_stage2D = 1
        md_oritype = 'all'
        l_selection = .false.
        if( present(oritype)   ) md_oritype = oritype
        if( present(stage2D)   ) i_stage2D = stage2D
        if( present(selection) ) l_selection = selection
        if( c_pthread_mutex_lock(gui_comm_args_inst%metadata_mutex) /= 0   ) THROW_HARD('failed to lock metadata mutex')
        call gui_project_metadata_inst%set(spproj, md_oritype, i_stage2D, l_selection)
        if( c_pthread_mutex_unlock(gui_comm_args_inst%metadata_mutex) /= 0 ) THROW_HARD('failed to unlock metadata mutex')
    end subroutine add_metadata_1

    subroutine add_metadata_2( self, projfile, oritype, stage2D, selection )
        class(gui_communicator),    intent(inout) :: self
        type(string),               intent(in)    :: projfile
        character(len=*), optional, intent(in)    :: oritype
        integer,          optional, intent(in)    :: stage2D
        logical,          optional, intent(in)    :: selection
        type(sp_project)                          :: spproj
        character(len=SHORTSTRLEN)                :: md_oritype
        integer                                   :: i_stage2D
        logical                                   :: l_selection
        if( .not. self%is_active        ) return
        if( .not. file_exists(projfile) ) return
        l_selection = .false.
        if( present(selection) ) l_selection = selection
        if( present(oritype) ) then
            if( oritype == 'mic' ) then
                call spproj%read_segment('mic', projfile)
            else if( oritype == 'ptcl' ) then
                call spproj%read_segment('mic',    projfile)
                call spproj%read_segment('stk',    projfile)
                call spproj%read_segment('ptcl2D', projfile)    
            else if( oritype == 'cls2D' ) then
                call spproj%read_segment('cls2D', projfile)
            end if
        else
            call spproj%read_segment('mic',      projfile)
            call spproj%read_segment('cls2D',    projfile)
        end if
        call spproj%read_segment('projinfo', projfile)
        call spproj%read_segment('out',      projfile)
        i_stage2D  = 1
        md_oritype = 'all'
        if( present(stage2D) ) i_stage2D  = stage2D
        if( present(oritype) ) md_oritype = oritype
        if( c_pthread_mutex_lock(gui_comm_args_inst%metadata_mutex) /= 0   ) THROW_HARD('failed to lock metadata mutex')
        call gui_project_metadata_inst%set(spproj, md_oritype, i_stage2D, l_selection)
        if( c_pthread_mutex_unlock(gui_comm_args_inst%metadata_mutex) /= 0 ) THROW_HARD('failed to unlock metadata mutex')
        call spproj%kill()
    end subroutine add_metadata_2

    subroutine communication_worker( carg ) bind(c)
        integer, parameter            :: SEND_INTERVAL_SEC = 10
        type(c_ptr), value            :: carg
        type(gui_comm_args),        pointer  :: comm_args
        type(gui_metadata_project), pointer  :: metadata_project  => null()
        type(json_value),           pointer  :: json_response_ptr => null()
        type(http_post)               :: post
        type(http_response)           :: response
        type(gui_assembler)           :: assembler
        type(string)                  :: my_niceserver_url, my_request
        type(json_core)               :: json
        character(len=STDLEN)         :: url_buffer
        integer                       :: my_rc, last_send_time, now
        logical                       :: my_l_continue, my_l_terminate, my_l_found, my_l_terminate_requested
        call c_f_pointer(carg, comm_args)
        call c_f_pointer(comm_args%metadata_project, metadata_project)
        url_buffer          = transfer(comm_args%url, url_buffer)
        my_niceserver_url   = string(url_buffer(1:index(url_buffer, c_null_char) - 1))
        my_l_continue       = .true.
        my_l_terminate      = .false.
        last_send_time      = 0
        call post%new(my_niceserver_url)
        call assembler%new(comm_args%procid)
        write(logfhandle, '(A,A)') 'Starting communication worker with server URL: ', my_niceserver_url%to_char()
        do while( my_l_continue )
            ! checked every pass, independent of the send throttle, so shutdown latency is bounded by the sleep below
            if( c_pthread_mutex_lock(comm_args%terminate_mutex) /= 0 ) THROW_HARD('failed to lock terminate mutex')
            if( comm_args%l_terminate ) my_l_terminate = .true.
            if( c_pthread_mutex_unlock(comm_args%terminate_mutex) /= 0 ) THROW_HARD('failed to unlock terminate mutex')
            now = int(c_time(0_c_long))
            if( .not. my_l_terminate .and. now - last_send_time < SEND_INTERVAL_SEC ) then
                my_rc = c_usleep(500000)
                cycle
            end if
            if( my_l_terminate ) then
                call assembler%set_stoptime()
                my_l_continue = .false.
            end if
            call assembler%assemble_batch_heartbeat()
            if( c_pthread_mutex_lock(comm_args%metadata_mutex) /= 0 ) THROW_HARD('failed to lock metadata mutex')
            call assembler%assemble_batch_metadata(metadata_project)
            if( c_pthread_mutex_unlock(comm_args%metadata_mutex) /= 0 ) THROW_HARD('failed to unlock metadata mutex')
            last_send_time = now
            my_request     = assembler%to_string()
            if( post%request(response, my_request) ) then
                if( response%code == 200) then
                    ! parse response JSON
                    call json%parse(json_response_ptr, response%content%to_char())
                    if( json%failed()) then
                        write(logfhandle, '(A,A)') "FAILED TO PARSE JSON RESPONSE ", response%content%to_char()
                        call json%clear_exceptions()
                    else
                        ! check for master process termination
                        call json%get(json_response_ptr, 'terminate', my_l_terminate_requested, my_l_found)
                        if( my_l_found .and. my_l_terminate_requested ) then
                            ! threads share the process's signal handlers, so this reaches the master's SIGTERM handler
                            my_rc = c_kill(c_getpid(), SIGTERM)
                        endif
                    end if
                    call safe_destroy_json_ptr(json_response_ptr)
                else
                    write(logfhandle, '(A,I0)') "HTTP RESPONSE CODE: ", response%code
                end if
            else
                write(logfhandle, '(A)') "FAILED TO SEND HTTP REQUEST"
            end if
            ! sleep briefly to avoid busy-polling when queue is idle
            my_rc = c_usleep(10000)
        end do
        call assembler%kill()
        call post%kill()

        contains

            subroutine safe_destroy_json_ptr(ptr_json)
                type(json_value), pointer, intent(inout) :: ptr_json
                if( associated(ptr_json) ) call json%destroy(ptr_json)
                nullify(ptr_json)
            end subroutine safe_destroy_json_ptr

    end subroutine communication_worker

end module simple_gui_communicator
