module simple_stream_communicator
    use, intrinsic :: iso_c_binding 
    include 'simple_lib.f08'
    use json_kinds
    use json_module
    use simple_fileio

    implicit none

    type, public :: stream_http_communicator

        character(len=XLONGSTRLEN),          private :: url       = ""
        character(len=XLONGSTRLEN),          private :: checksum  = ""
        character(len=LONGSTRLEN),           private :: tmp_recv  = "__stream_comm_recv__"
        character(len=LONGSTRLEN),           private :: tmp_send  = "__stream_comm_send__"
        integer,                             private :: id        = 0
        logical,                             public  :: exit      = .false.
        logical,                             public  :: stop      = .false.
        type(json_core),                     public  :: json
        type(json_value),           pointer, public  :: job_json, heartbeat_json, update_arguments

        contains
            procedure, public   :: create
            procedure, public   :: term
            procedure, public   :: send_heartbeat
            procedure, public   :: send_jobstats
            procedure, private  :: evaluate_checksum
            procedure, private  :: curl_request

    end type stream_http_communicator

    contains

        subroutine create(self, jobid, url, name)
            class(stream_http_communicator), intent(inout) :: self
            character(*),                    intent(in)    :: url, name
            integer,                         intent(in)    :: jobid
            call self%json%initialize(no_whitespace=.true.)
            call self%json%create_object(self%job_json,         name)
            call self%json%create_object(self%heartbeat_json,   name // "_heartbeat")
            call self%json%create_object(self%update_arguments, "")
            if(jobid > 0 .and. url .ne. "" .and. name .ne. "") then
                self%id   = jobid
                self%url  = url
                self%tmp_recv  = "__stream_comm_recv_" // name // "__"
                self%tmp_send  = "__stream_comm_send_" // name // "__"
            endif
        end subroutine create

        subroutine term(self)
            class(stream_http_communicator), intent(inout) :: self
            type(json_value),                pointer       :: jobstats_json
            integer                                        :: file_unit, stat
            call self%json%add(self%job_json, "terminate", .true.)
            call self%json%create_object(jobstats_json,'')
            call self%json%add(jobstats_json, "jobid", self%id)
            call self%json%add(jobstats_json, self%job_json)
            ! Open the temporary file.
            if(file_exists(trim(self%tmp_send))) call del_file(trim(self%tmp_send))
            call fopen(action='write', file=trim(self%tmp_send), iostat=stat, funit=file_unit, status='new')
            if (stat /= 0) then
                print '("Error: opening file ", a, " failed: ", i0)', trim(self%tmp_send), stat
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
            call self%json%create_object(heartbeat_json,'')
            call self%json%add(heartbeat_json, "jobid",     self%id)
            call self%json%add(heartbeat_json, self%heartbeat_json)
            call self%json%print_to_string(heartbeat_json, heartbeat_json_str)
            call self%curl_request(heartbeat_json_str)
        end subroutine send_heartbeat

        subroutine send_jobstats(self)
            class(stream_http_communicator), intent(inout) :: self
            type(json_value),                pointer       :: jobstats_json
            integer                                        :: file_unit, stat
            call self%json%create_object(jobstats_json,'')
            call self%json%add(jobstats_json, "jobid", self%id)
            call self%json%add(jobstats_json, self%job_json)
            ! Open the temporary file.
            if(file_exists(trim(self%tmp_send))) call del_file(trim(self%tmp_send))
            call fopen(action='readwrite', file=self%tmp_send, iostat=stat, funit=file_unit, status='new')
            if (stat /= 0) then
                print '("Error: opening file ", a, " failed: ", i0)', trim(self%tmp_send), stat
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
        end subroutine send_jobstats

        subroutine curl_request(self, data)
            class(stream_http_communicator),           intent(inout) :: self
            character(*),                    optional, intent(in)    :: data
            type(json_value),                pointer                 :: response_json
            character(len=512)            :: buf
            character(len=:), allocatable :: cmd, tmp_file
            integer                       :: file_unit, stat
            logical                       :: found, terminate
            if(present(data)) then
                cmd = 'curl --header "Content-Type: application/json" --silent --request GET --data ' // "'" // data // "'" // ' -o ' // trim(self%tmp_recv) // ' ' // trim(self%url)
            else
                cmd = 'curl --header "Content-Type: application/json" --silent --request GET --data ' // "'@" // trim(self%tmp_send) // "'" // ' -o ' // trim(self%tmp_recv) // ' ' // trim(self%url)
            endif
            ! Run cURL from command-line.
            call execute_command_line(cmd, exitstat=stat)
            if (stat /= 0) then
                print '("Error: HTTP request failed: ", i0)', stat
                return
            end if
            ! Open the temporary file.
            call fopen(action='readwrite', file=self%tmp_recv, iostat=stat, funit=file_unit, status='old')
            if (stat /= 0) then
                print '("Error: reading file ", a, " failed: ", i0)', tmp_file, stat
                return
            end if
            ! Read and output contents of file.
            do
                read (file_unit, '(a)', iostat=stat) buf
                if (stat /= 0) exit
                print '(a)', trim(buf)
            end do
            ! Close and delete file.
            if( is_open(file_unit) ) close (file_unit, status='delete')
            ! parse response
            call self%json%parse(response_json, trim(buf))
            ! terminate if requested
            call self%json%get(response_json, 'terminate', terminate, found)
            call self%json%clone(response_json, self%update_arguments)
            if(found) then
                if(terminate) then
                    self%exit = .true.
                endif
            endif
        end subroutine curl_request

        function evaluate_checksum(self)
            class(stream_http_communicator), intent(inout) :: self
            character(len=512)            :: buf
            character(len=:), allocatable :: cmd
            integer                       :: file_unit, stat
            logical                       :: evaluate_checksum, found, terminate
            cmd = 'md5sum ' // trim(self%tmp_send) // ' >> ' // trim(self%tmp_recv)
            call execute_command_line(cmd, exitstat=stat)
            if (stat /= 0) then
                print '("Error: MD5 sum failed: ", i0)', stat
                evaluate_checksum = .false.
                return
            end if
            ! Open the temporary file.
            call fopen(action='readwrite', file=trim(self%tmp_recv), iostat=stat, funit=file_unit, status='old')
            if (stat /= 0) then
                print '("Error: reading file ", a, " failed: ", i0)', trim(self%tmp_recv), stat
                evaluate_checksum = .false.
            end if
            ! Read and output contents of file.
            do
                read (file_unit, '(a)', iostat=stat) buf
                if (stat /= 0) exit
            end do
            ! Close and delete file.
            if( is_open(file_unit) ) close (file_unit, status='delete')
            if(trim(buf) == trim(self%checksum)) then
                evaluate_checksum = .true.
            else
                self%checksum = trim(buf)
                evaluate_checksum = .false.
            endif
        end function evaluate_checksum

end module simple_stream_communicator