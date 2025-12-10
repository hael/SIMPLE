module simple_stream_communicator
use, intrinsic :: iso_c_binding 
include 'simple_lib.f08'
use json_kinds
use json_module
use simple_fileio
use unix
implicit none

public :: stream_http_communicator
private

type :: stream_http_communicator
    private
    type(string)              :: url
    type(string)              :: checksum
    type(string)              :: tmp_recv
    type(string)              :: tmp_send
    integer                   :: id        = 0
    integer                   :: bg_pid    = 0
    integer                   :: bg_pipe(2)
    logical                   :: active    = .false.
    logical                   :: exit      = .false.
    logical                   :: stop      = .false.
    type(json_core)           :: json
    type(json_value), pointer :: job_json, heartbeat_json, update_arguments
 contains
    ! LIFECYCLE
    procedure           :: create
    procedure           :: term
    procedure           :: send_heartbeat
    procedure           :: send_jobstats
    procedure           :: background_heartbeat
    procedure           :: join_background_heartbeat
    procedure, private  :: evaluate_checksum
    procedure, private  :: curl_request
    ! STATUS CHECKS AND JSON ARGUMENT MANAGEMENT
    procedure           :: active_status
    procedure           :: exit_status
    procedure           :: stop_status
    procedure           :: arg_associated
    procedure           :: destroy_arg
    ! JSON MODIFIERS
    procedure, private  :: update_json_1, update_json_2, update_json_3, update_json_4
    generic             :: update_json => update_json_1, update_json_2, update_json_3, update_json_4
    procedure, private  :: add_to_json_1, add_to_json_2, add_to_json_3, add_to_json_4, add_to_json_5
    generic             :: add_to_json => add_to_json_1, add_to_json_2, add_to_json_3, add_to_json_4, add_to_json_5
    procedure           :: remove_from_json_if_present
    ! JSON GETTERS
    procedure, private  :: get_json_arg_1, get_json_arg_2, get_json_arg_3, get_json_arg_4, get_json_arg_5
    generic             :: get_json_arg => get_json_arg_1, get_json_arg_2, get_json_arg_3, get_json_arg_4, get_json_arg_5
end type stream_http_communicator

contains

    ! LIFECYCLE

    subroutine create(self, jobid, url, name)
        class(stream_http_communicator), intent(inout) :: self
        character(*),                    intent(in)    :: url
        character(*), optional,          intent(in)    :: name
        integer,                         intent(in)    :: jobid
        character(len=STDLEN)                          :: l_name
        ! initialize dynamicstrings
        self%url       = ""
        self%checksum  = ""
        self%tmp_recv  = "__stream_comm_recv__"
        self%tmp_send  = "__stream_comm_send__"
        ! proceed
        l_name = 'job'
        call self%json%initialize(no_whitespace=.true.)
        if(present(name)) l_name = name
        call self%json%create_object(self%job_json,         trim(l_name))
        call self%json%create_object(self%heartbeat_json,   trim(l_name) // "_heartbeat")
        call self%json%create_object(self%update_arguments, "")
        if(jobid > 0 .and. url .ne. "" .and. trim(l_name) .ne. "") then
            self%id        = jobid
            self%url       = trim(url)
            self%tmp_recv  = "__stream_comm_recv_" // trim(l_name) // "__"
            self%tmp_send  = "__stream_comm_send_" // trim(l_name) // "__"
            self%active    = .true.
        endif
    end subroutine create

    subroutine term(self)
        class(stream_http_communicator), intent(inout) :: self
        type(json_value),                pointer       :: jobstats_json
        integer                                        :: file_unit, stat
        if( .not.self%active ) return
        call self%json%add(self%job_json, "terminate", .true.)
        call self%json%create_object(jobstats_json,'')
        call self%json%add(jobstats_json, "jobid", self%id)
        call self%json%add(jobstats_json, self%job_json)
        ! Open the temporary file.
        if(file_exists(self%tmp_send)) call del_file(self%tmp_send)
        call fopen(action='write', file=self%tmp_send, iostat=stat, funit=file_unit, status='new')
        if (stat /= 0) then
            call simple_exception("Error opening file "//self%tmp_send%to_char(), __FILENAME__ , __LINE__,  l_stop=.false.)
            return
        end if
        call self%json%print(jobstats_json, file_unit)
        ! Close and delete file.
        if( is_open(file_unit) ) close (file_unit)
        if(self%evaluate_checksum()) then
            call self%send_heartbeat()
        else
            call self%curl_request()
        endif
    end subroutine term

    subroutine send_heartbeat(self)
        class(stream_http_communicator), intent(inout) :: self
        type(json_value),                pointer       :: heartbeat_json
        character(kind=CK,len=:), allocatable          :: heartbeat_json_str
        if(self%active) then
            call self%json%create_object(heartbeat_json,'')
            call self%json%add(heartbeat_json, "jobid",     self%id)
            call self%json%add(heartbeat_json, self%heartbeat_json)
            call self%json%print_to_string(heartbeat_json, heartbeat_json_str)
            call self%curl_request(heartbeat_json_str)
        endif
    end subroutine send_heartbeat

    subroutine send_jobstats(self)
        class(stream_http_communicator), intent(inout) :: self
        type(json_value),                pointer       :: jobstats_json
        integer                                        :: file_unit, stat
        if(self%active) then
            call self%json%create_object(jobstats_json,'')
            call self%json%add(jobstats_json, "jobid", self%id)
            call self%json%add(jobstats_json, self%job_json)
            ! Open the temporary file.
            if(file_exists(self%tmp_send)) call del_file(self%tmp_send)
            call fopen(action='readwrite', file=self%tmp_send, iostat=stat, funit=file_unit, status='new')
            if (stat /= 0) then
                call simple_exception("Error opening file "//self%tmp_send%to_char(), __FILENAME__ , __LINE__,  l_stop=.false.)
                return
            end if
            call self%json%print(jobstats_json, file_unit)
            ! Close file.
            if( is_open(file_unit) ) close (file_unit)
            if(self%evaluate_checksum()) then
                call self%send_heartbeat()
            else
                call self%curl_request()
            endif
        endif
    end subroutine send_jobstats

    subroutine background_heartbeat(self)
        class(stream_http_communicator), intent(inout) :: self
        character,                       target        :: rcv, snd
        integer                                        :: rc, i
        if(.not. self%active) return ! this process isn't communicating with nice
        ! create pipe
        rc = c_pipe(self%bg_pipe)
        if (rc < 0) call simple_exception("Creating background heartbeat pipe failed", __FILENAME__ , __LINE__)
        ! set pipe reads to non-blocking
        rc = c_fcntl(self%bg_pipe(1), F_SETFL, O_NONBLOCK); 
        if (rc < 0) call simple_exception("Updating background heartbeat pipe failed", __FILENAME__ , __LINE__)
        ! fork process
        self%bg_pid = c_fork()
        if (self%bg_pid < 0) then
            ! Fork failed.
            call simple_exception("Background heartbeat fork failed", __FILENAME__ , __LINE__)
        else if (self%bg_pid == 0) then
            ! Child process infinite do loop until term signal recieved
            i = 0
            do
                ! Send heartbeat every 10 seconds
                if(i == 10) then
                    call self%send_heartbeat()
                    i = 0
                endif
                ! Read character from pipe. If present and T => terminate
                if(c_read(self%bg_pipe(1), c_loc(rcv), 1_c_size_t) > 0) then
                    if(rcv == 'T') exit
                endif
                ! Sleep 1 second
                call sleep(1)
                i = i + 1
            end do
            ! send parent X if self%exit is false, else Y
            snd = 'Y'
            if(self%exit) snd = 'X'
            rc = c_write(self%bg_pipe(2), c_loc(snd), len(snd, kind=c_size_t))
            if (rc < 0) call simple_exception("Child write to background heartbeat pipe failed", __FILENAME__ , __LINE__)
            ! close pipes and exit
            rc = c_close(self%bg_pipe(1))
            if (rc < 0) call simple_exception("Child close background heartbeat pipe failed 1", __FILENAME__ , __LINE__)
            rc = c_close(self%bg_pipe(2))
            if (rc < 0) call simple_exception("Child close background heartbeat pipe failed 2", __FILENAME__ , __LINE__)
            call c_exit(0)
        else
            ! Parent process. continue as normal
        end if
    end subroutine background_heartbeat

    subroutine join_background_heartbeat(self)
        class(stream_http_communicator), intent(inout) :: self
        character,                       target        :: rcv, snd
        integer                                        :: pid, rc
        if(.not. self%active) return ! this process isn't communicating with nice
        snd = 'T'
        ! send child T to terminate
        rc = c_write(self%bg_pipe(2), c_loc(snd), len(snd, kind=c_size_t))
        if (rc < 0) call simple_exception("Parent write to background heartbeat pipe failed", __FILENAME__ , __LINE__)
        ! wait for child to terminate
        pid = c_waitpid(self%bg_pid, rc, 0)
        if (rc < 0) call simple_exception("Parent wait for background heartbeat termination failed", __FILENAME__ , __LINE__)
        ! Read character from pipe. If present and X => set self%exit to true
        if(c_read(self%bg_pipe(1), c_loc(rcv), 1_c_size_t) > 0) then
            if(rcv == 'X') self%exit = .true.
        endif
        ! close pipes
        rc = c_close(self%bg_pipe(1))
        if (rc < 0) call simple_exception("Parent close background heartbeat pipe failed 1", __FILENAME__ , __LINE__)
        rc = c_close(self%bg_pipe(2))
        if (rc < 0) call simple_exception("Parent close background heartbeat pipe failed 2", __FILENAME__ , __LINE__)
    end subroutine join_background_heartbeat

    subroutine curl_request(self, data)
        class(stream_http_communicator), intent(inout) :: self
        character(*), optional,          intent(in)    :: data
        type(json_value), pointer :: response_json
        type(string) :: cmd, buf
        integer      :: file_unit, stat
        logical      :: found, terminate
        if(present(data)) then
            cmd = 'curl --header "Content-Type: application/json" --silent --request GET --data ' // "'" // data // "'" // ' -o ' // self%tmp_recv%to_char() // ' ' // self%url%to_char()
        else
            cmd = 'curl --header "Content-Type: application/json" --silent --request GET --data ' // "'@" // self%tmp_send%to_char() // "'" // ' -o ' // self%tmp_recv%to_char() // ' ' // self%url%to_char()
        endif
        ! Run cURL from command-line.
        call execute_command_line(cmd%to_char(), exitstat=stat)
        if (stat /= 0) then
            call simple_exception("Error HTTP request failed", __FILENAME__ , __LINE__,  l_stop=.false.)
            return
        end if
        ! Open the temporary file.
        call fopen(action='readwrite', file=self%tmp_recv, iostat=stat, funit=file_unit, status='old')
        if (stat /= 0) then
            call simple_exception("Error opening file "//self%tmp_recv%to_char(), __FILENAME__ , __LINE__,  l_stop=.false.)
            return
        end if
        ! Read and output contents of file.
        do
            call buf%readline(file_unit, stat)
            if (stat /= 0) exit
        !    print '(a)', trim(buf) ! uncomment to view raw json responses
        end do
        ! Close and delete file.
        if( is_open(file_unit) ) close (file_unit, status='delete')
        ! parse response
        call self%json%parse(response_json, buf%to_char())
        ! terminate if requested
        call self%json%get(response_json, 'terminate', terminate, found)
        call self%json%clone(response_json, self%update_arguments)
        if(found) then
            if(terminate) then
                self%exit = .true.
            endif
        endif
        call cmd%kill
        call buf%kill
    end subroutine curl_request

    function evaluate_checksum(self)
        class(stream_http_communicator), intent(inout) :: self
        character(len=512)            :: buf
        character(len=:), allocatable :: cmd
        integer                       :: file_unit, stat
        logical                       :: evaluate_checksum
        cmd = 'md5sum ' // self%tmp_send%to_char() // ' >> ' // self%tmp_recv%to_char()
        call execute_command_line(cmd, exitstat=stat)
        if (stat /= 0) then
            call simple_exception("Error calculating MD5 sum", __FILENAME__ , __LINE__,  l_stop=.false.)
            evaluate_checksum = .false.
            return
        end if
        ! Open the temporary file.
        call fopen(action='readwrite', file=self%tmp_recv, iostat=stat, funit=file_unit, status='old')
        if (stat /= 0) then
            call simple_exception("Error reading file "//self%tmp_recv%to_char(), __FILENAME__ , __LINE__,  l_stop=.false.)
            evaluate_checksum = .false.
        end if
        ! Read and output contents of file.
        do
            read (file_unit, '(a)', iostat=stat) buf
            if (stat /= 0) exit
        end do
        ! Close and delete file.
        if( is_open(file_unit) ) close (file_unit, status='delete')
        if(trim(buf) == self%checksum%to_char() ) then
            evaluate_checksum = .true.
        else
            self%checksum = trim(buf)
            evaluate_checksum = .false.
        endif
    end function evaluate_checksum

    ! STATUS CHECKS AND JSON ARGUMENT MANAGEMENT

    pure logical function active_status( self )
        class(stream_http_communicator), intent(in) :: self
        active_status = self%active
    end function active_status

    pure logical function exit_status( self )
        class(stream_http_communicator), intent(in) :: self
        exit_status = self%exit
    end function exit_status

    pure logical function stop_status( self )
        class(stream_http_communicator), intent(in) :: self
        stop_status = self%stop
    end function stop_status

    pure logical function arg_associated( self )
        class(stream_http_communicator), intent(in) :: self
        arg_associated = associated(self%update_arguments)
    end function arg_associated

    subroutine destroy_arg( self )
        class(stream_http_communicator), intent(inout) :: self
        call self%json%destroy(self%update_arguments)
        nullify(self%update_arguments)
    end subroutine destroy_arg

    ! JSON MODIFIERS

    subroutine update_json_1(self, key, value, found )
        class(stream_http_communicator), intent(inout) :: self
        character(len=*),                intent(in)    :: key, value
        logical,                         intent(out)   :: found
        call self%json%update(self%job_json, trim(key), trim(value), found)
    end subroutine update_json_1

    subroutine update_json_2(self, key, ival, found )
        class(stream_http_communicator), intent(inout) :: self
        character(len=*),                intent(in)    :: key
        integer,                         intent(in)    :: ival
        logical,                         intent(out)   :: found
        call self%json%update(self%job_json, trim(key), ival, found)
    end subroutine update_json_2

    subroutine update_json_3(self, key, rval, found )
        class(stream_http_communicator), intent(inout) :: self
        character(len=*),                intent(in)    :: key
        real(dp),                        intent(in)    :: rval
        logical,                         intent(out)   :: found
        call self%json%update(self%job_json, trim(key), rval, found)
    end subroutine update_json_3
   
    subroutine update_json_4(self, key, lval, found )
        class(stream_http_communicator), intent(inout) :: self
        character(len=*),                intent(in)    :: key
        logical,                         intent(in)    :: lval
        logical,                         intent(out)   :: found
        call self%json%update(self%job_json, trim(key), lval, found)
    end subroutine update_json_4

    subroutine add_to_json_1( self, key, ival )
        class(stream_http_communicator), intent(inout) :: self
        character(len=*),                intent(in)    :: key
        integer,                         intent(in)    :: ival
        call self%json%add(self%job_json, trim(key), ival)
    end subroutine add_to_json_1

    subroutine add_to_json_2( self, key, rval )
        class(stream_http_communicator), intent(inout) :: self
        character(len=*),                intent(in)    :: key
        real(dp),                        intent(in)    :: rval
        call self%json%add(self%job_json, trim(key), rval)
    end subroutine add_to_json_2

    subroutine add_to_json_3( self, key, value )
        class(stream_http_communicator), intent(inout) :: self
        character(len=*),                intent(in)    :: key, value
        call self%json%add(self%job_json, trim(key), trim(value))
    end subroutine add_to_json_3

    subroutine add_to_json_4( self, jval )
        class(stream_http_communicator), intent(inout) :: self
        type(json_value), pointer,       intent(in)    :: jval
        call self%json%add(self%job_json, jval)
    end subroutine add_to_json_4

    subroutine add_to_json_5( self, key, lval )
        class(stream_http_communicator), intent(inout) :: self
        character(len=*),                intent(in)    :: key
        logical,                         intent(in)    :: lval
        call self%json%add(self%job_json, trim(key), lval)
    end subroutine add_to_json_5

    subroutine remove_from_json_if_present( self, str )
        class(stream_http_communicator), intent(inout) :: self
        character(len=*),                intent(in)    :: str
        call self%json%remove_if_present(self%job_json, trim(str))
    end subroutine remove_from_json_if_present

    ! JSON GETTERS

    subroutine get_json_arg_1( self, key, ival, found )
        class(stream_http_communicator), intent(inout) :: self
        character(len=*),                intent(in)    :: key
        integer,                         intent(inout) :: ival
        logical,                         intent(out)   :: found
        call self%json%get(self%update_arguments, trim(key), ival, found)
    end subroutine get_json_arg_1

    subroutine get_json_arg_2( self, key, rval, found )
        class(stream_http_communicator), intent(inout) :: self
        character(len=*),                intent(in)    :: key
        real(dp),                        intent(inout) :: rval
        logical,                         intent(out)   :: found
        call self%json%get(self%update_arguments, trim(key), rval, found)
    end subroutine get_json_arg_2

    subroutine get_json_arg_3( self, key, str, found )
        class(stream_http_communicator),       intent(inout) :: self
        character(len=*),                      intent(in)    :: key
        character(kind=CK,len=:), allocatable, intent(inout) :: str
        logical,                               intent(out)   :: found
        call self%json%get(self%update_arguments, trim(key), str, found)
    end subroutine get_json_arg_3

    subroutine get_json_arg_4( self, key, lval, found )
        class(stream_http_communicator), intent(inout) :: self
        character(len=*),                intent(in)    :: key
        logical,                         intent(inout) :: lval
        logical,                         intent(out)   :: found
        call self%json%get(self%update_arguments, trim(key), lval, found)
    end subroutine get_json_arg_4

    subroutine get_json_arg_5( self, key, iarr, found )
        class(stream_http_communicator), intent(inout) :: self
        character(len=*),                intent(in)    :: key
        integer, allocatable,            intent(inout) :: iarr(:)
        logical,                         intent(out)   :: found
        if( allocated(iarr) ) deallocate(iarr)
        call self%json%get(self%update_arguments, trim(key), iarr, found)
    end subroutine get_json_arg_5

end module simple_stream_communicator