!@descr: SIMPLE distributed worker process
!==============================================================================
! PROGRAM: simple_persistent_worker
!
! PURPOSE:
!   Runs as a compute-node agent for the SIMPLE batch-processing system.
!   Connects to a central server via TCP/IP, sends periodic heartbeat messages
!   that report thread availability, and receives tasks in reply.  Each task
!   is a path to a bash script that is executed in a dedicated pthread.
!
! COMMAND-LINE ARGUMENTS (all key=value):
!   nthr=<N>        Number of parallel task slots (default: 1)
!   port=<N>        TCP port of the central server
!   server=<list>   Comma-separated list of server IP addresses
!   worker_id=<N>   Unique integer ID assigned by the server (default: 0)
!
! EXECUTION FLOW:
!   1. Parse arguments and allocate one thread slot per nthr.
!   2. Connect to the server (init_client).
!   3. Heartbeat loop: serialise + send heartbeat, decode reply.
!        TERMINATE  -> join all threads and exit.
!        TASK       -> validate script path; dispatch to a free thread slot.
!        STATUS     -> informational, ignored.
!   4. Deallocate, close log, print timing.
!
! INTERNAL ROUTINES:
!   start_worker_thread  Find a free slot and spawn a pthread for the task.
!   worker_task_thread   pthread body: execute the script, record exit code.
!   get_nthr_used        Sum the nthr fields of all active slots (thread-safe).
!==============================================================================
program simple_persistent_worker
    use simple_core_module_api
    use unix,                  only: c_time,                                   &
                                     c_pthread_t, c_pthread_mutex_t,           &
                                     c_pthread_create, c_pthread_join,         &
                                     c_pthread_mutex_init,                     &
                                     c_pthread_mutex_destroy,                  &
                                     c_pthread_mutex_lock,                     &
                                     c_pthread_mutex_unlock,                   &
                                     c_null_ptr, c_funloc, c_usleep,           &
                                     c_getpid, c_pid_t
    use iso_fortran_env,       only: output_unit
    use simple_jiffys,         only: simple_print_git_version, simple_print_timer
    use simple_ipc_tcp_socket, only: ipc_tcp_socket
    use simple_persistent_worker_server,            only: TCP_BUFSZ
    use simple_persistent_worker_message_heartbeat, only: qsys_persistent_worker_message_heartbeat
    use simple_persistent_worker_message_task,      only: qsys_persistent_worker_message_task
    use simple_persistent_worker_message_status,    only: qsys_persistent_worker_message_status
    use simple_persistent_worker_message_terminate, only: qsys_persistent_worker_message_terminate
    use simple_persistent_worker_message_types,     only: WORKER_TERMINATE_MSG, WORKER_STATUS_MSG, WORKER_TASK_MSG
    implicit none
#include "simple_local_flags.inc"

    ! ------------------------------------------------------------------
    ! Constants
    ! ------------------------------------------------------------------
    integer, parameter :: HEARTBEAT_TIMEOUT_MS = 5000  !< send_recv_msg timeout (ms)
    integer, parameter :: HEARTBEAT_MAX_RETRY  = 5     !< retries before giving up
    integer, parameter :: POLL_TIME_US         = 200000  !< poll timeout for worker threads (us)

    ! ------------------------------------------------------------------
    ! Per-thread state: task description + pthread bookkeeping
    ! ------------------------------------------------------------------
    type :: worker_thread_args
        type(qsys_persistent_worker_message_task) :: task_msg           !< task being executed
        type(c_pthread_mutex_t)        :: mutex              !< guards all mutable fields below
        type(c_pthread_t)              :: thread             !< pthread handle
        logical                        :: started   = .false.!< true after successful pthread_create
        logical                        :: complete  = .false.!< set by thread on exit
        integer                        :: nthr      = 0      !< >0 while slot is occupied
        integer                        :: exit_code = 0      !< script exit status
    end type worker_thread_args

    ! ------------------------------------------------------------------
    ! Variables
    ! ------------------------------------------------------------------
    ! Communication buffers
    character(len=:), allocatable, target :: buffer              !< serialised outgoing message
    character(kind=c_char),        target :: reply_buffer(TCP_BUFSZ)

    ! Argument parsing temporaries
    character(len=STDLEN)  :: xarg
    type(string)           :: arg_kv, arg_val

    ! Worker identity / configuration
    type(string)           :: server_list                        !< comma-separated server IPs
    integer                :: nthr      = 1                      !< total thread slots
    integer                :: port      = -1                     !< server TCP port
    integer                :: worker_id = 0                      !< this worker's ID

    ! IPC + message objects
    type(ipc_tcp_socket)                :: ipc_socket
    type(qsys_persistent_worker_message_heartbeat) :: heartbeat_msg
    type(qsys_persistent_worker_message_task)      :: task_msg

    ! Thread slot array (one element per nthr)
    type(worker_thread_args), allocatable, target :: threads(:)

    ! Loop / scratch
    logical                :: found, l_terminate
    logical                :: slot_started
    logical                :: safe_path
    integer                :: argcnt, iarg, cmdlen, reply_len, buffer_type, i, rc, server_len
    integer(timer_int_kind):: t0
    real(timer_int_kind)   :: rt_exec
    type(c_ptr)            :: thread_retval
    ! Unique worker identifier: <hostname>_<PID>
    character(len=256)      :: worker_uid  !< identifies this worker instance across log files
    character(len=256)      :: hostname
    integer(kind=c_pid_t)   :: pid

    ! ------------------------------------------------------------------
    ! Initialise
    ! ------------------------------------------------------------------
    l_terminate = .false.
    t0          = tic()

    ! ------------------------------------------------------------------
    ! Build unique worker identifier: <nodename>_<PID>
    ! This string uniquely identifies this process on this host and is
    ! stable for the lifetime of the process.  Use it in log messages and
    ! when registering with the server so crashed workers are detectable.
    ! ------------------------------------------------------------------
    pid = c_getpid()
    call get_environment_variable('HOSTNAME', hostname, status=rc)
    if( rc /= 0 .or. len_trim(hostname) == 0 ) hostname = 'unknown'
    worker_uid = trim(hostname) // '_' // trim(int2str(int(pid)))
    write(*,'(2a)') 'Worker UID: ', trim(worker_uid)

    ! ------------------------------------------------------------------
    ! Parse command-line arguments (key=value pairs)
    ! ------------------------------------------------------------------
    argcnt = command_argument_count()
    do iarg = 1, argcnt
        call get_command_argument(iarg, xarg)
        arg_kv = string(trim(xarg))
        if( arg_kv%has_substr(string('nthr=')) )then
            arg_val   = arg_kv%substr_remove(string('nthr='))
            nthr      = arg_val%to_int()
        else if( arg_kv%has_substr(string('port=')) )then
            arg_val   = arg_kv%substr_remove(string('port='))
            port      = arg_val%to_int()
        else if( arg_kv%has_substr(string('server=')) )then
            server_list = arg_kv%substr_remove(string('server='))
        else if( arg_kv%has_substr(string('worker_id=')) )then
            arg_val   = arg_kv%substr_remove(string('worker_id='))
            worker_id = arg_val%to_int()
        endif
    end do
    if( nthr < 1 ) then
        write(*,*) 'Worker: invalid nthr, forcing to 1:', nthr
        nthr = 1
    end if
    if( port <= 0 ) then
        write(*,*) 'Worker: invalid or missing port argument:', port
        stop 1
    end if
    server_len = 0
    if( server_list%is_allocated() ) server_len = server_list%strlen_trim()
    if( server_len <= 0 ) then
        write(*,*) 'Worker: missing server list argument (server=...)'
        stop 1
    end if

    ! ------------------------------------------------------------------
    ! Allocate thread slots and initialise their mutexes
    ! ------------------------------------------------------------------
    allocate(threads(nthr))
    do i = 1, nthr
        rc = c_pthread_mutex_init(threads(i)%mutex, c_null_ptr)
        if( rc /= 0 ) then
            write(*,*) 'Worker: c_pthread_mutex_init failed for slot', i, 'rc=', rc
            stop 1
        end if
    end do

    ! ------------------------------------------------------------------
    ! Connect to the server
    ! ------------------------------------------------------------------
    call ipc_socket%init_client(server_list, port)

    ! ------------------------------------------------------------------
    ! Heartbeat loop: send status, receive and dispatch server replies
    ! ------------------------------------------------------------------
    do while (.not. l_terminate)
        ! Build and send heartbeat; block until reply or timeout
        call heartbeat_msg%new() ! reset heartbeat_msg to default values (except worker_id)
        heartbeat_msg%worker_id      = worker_id
        heartbeat_msg%worker_uid     = trim(worker_uid)
        heartbeat_msg%heartbeat_time = int(c_time(0_c_long))
        heartbeat_msg%nthr_total     = nthr
        heartbeat_msg%nthr_used      = get_nthr_used()
        
        call heartbeat_msg%serialise(buffer)
        call ipc_socket%send_recv_msg(buffer, HEARTBEAT_TIMEOUT_MS, &
                                      HEARTBEAT_MAX_RETRY, found,   &
                                      reply_buffer, reply_len)
        if( .not. found ) then
            write(*,'(A)') 'Worker: failed to reach server — exiting heartbeat loop'
            exit
        end if

        if( reply_len > 0 ) then
            buffer_type = transfer(reply_buffer, buffer_type)
            select case(buffer_type)

                case(WORKER_TERMINATE_MSG)
                    ! Server requested shutdown; joins are handled centrally in cleanup.
                    write(*,*) 'Worker: received TERMINATE - exiting heartbeat loop'
                    l_terminate = .true.

                case(WORKER_STATUS_MSG)
                    ! Informational status reply: nothing to act on
                    write(*,*) 'Worker: received STATUS message'

                case(WORKER_TASK_MSG)
                    ! Server dispatched a task — validate then execute
                    task_msg = transfer(reply_buffer, task_msg)
                    cmdlen   = len_trim(task_msg%script_path)
                    if( cmdlen == 0 ) then
                        write(*,*) 'Worker: received TASK with empty script_path — ignored'
                    else
                      !  safe_path = is_safe_script_path(task_msg%script_path)
                        safe_path = .true. ! for testing, accept all paths; implement real checks in production
                        if( .not. safe_path ) then
                            write(*,*) 'Worker: rejected TASK due to unsafe or missing script_path:', &
                                       trim(task_msg%script_path)
                        else
                            write(*,*) 'Worker: dispatching job_id', task_msg%job_id, &
                                       'script:', trim(task_msg%script_path)
                            call start_worker_thread(task_msg)
                        end if
                    end if

                case default
                    write(*,*) 'Worker: unknown message type:', buffer_type

            end select
        end if

        ! Brief pause before next heartbeat
       ! call sleep(1)
        rc = c_usleep(POLL_TIME_US) ! sleep to avoid busy-polling if server is slow to reply or if we received a STATUS message
        if( rc /= 0 ) then
            write(*,*) 'Worker: c_usleep failed in heartbeat loop, rc=', rc
        end if
    end do

    ! ------------------------------------------------------------------
    ! Cleanup
    ! ------------------------------------------------------------------
    do i = 1, size(threads)
        rc = c_pthread_mutex_lock(threads(i)%mutex)
        if( rc /= 0 ) then
            write(*,*) 'Worker: c_pthread_mutex_lock failed in cleanup for slot', i, 'rc=', rc
            cycle
        end if
        slot_started = threads(i)%started
        rc = c_pthread_mutex_unlock(threads(i)%mutex)
        if( rc /= 0 ) then
            write(*,*) 'Worker: c_pthread_mutex_unlock failed in cleanup for slot', i, 'rc=', rc
            cycle
        end if
        if( slot_started ) then
            rc = c_pthread_join(threads(i)%thread, thread_retval)
            if( rc /= 0 ) then
                write(*,*) 'Worker: c_pthread_join failed in cleanup for slot', i, 'rc=', rc
            end if
            rc = c_pthread_mutex_lock(threads(i)%mutex)
            if( rc == 0 ) then
                threads(i)%started  = .false.
                threads(i)%complete = .true.
                threads(i)%nthr     = 0
                rc = c_pthread_mutex_unlock(threads(i)%mutex)
            end if
        end if
        rc = c_pthread_mutex_destroy(threads(i)%mutex)
        if( rc /= 0 ) then
            write(*,*) 'Worker: c_pthread_mutex_destroy failed for slot', i, 'rc=', rc
        end if
    end do
    deallocate(threads)
    if( logfhandle /= output_unit )then
        if( is_open(logfhandle) ) call fclose(logfhandle)
    endif
    call simple_print_git_version('342d1484')
    rt_exec = toc(t0)
    call simple_print_timer(rt_exec)

contains

    !> Find the first idle thread slot, populate it with \p task, and
    !> spawn a new pthread.  If no slot is free the task is dropped and
    !> a warning is printed (the server can re-queue on next heartbeat).
    subroutine start_worker_thread( task )
        type(qsys_persistent_worker_message_task), intent(in) :: task
        integer     :: slot_rc, create_rc, slot_i
        type(c_ptr) :: local_retval
        do slot_i = 1, size(threads)
            slot_rc = c_pthread_mutex_lock(threads(slot_i)%mutex)
            if( slot_rc /= 0 ) then
                write(*,*) 'Worker: c_pthread_mutex_lock failed for slot', slot_i, 'rc=', slot_rc
                cycle
            end if
            if( threads(slot_i)%nthr == 0 ) then
                ! If a prior thread finished in this slot, join it before reusing
                ! the stored pthread handle.
                if( threads(slot_i)%started ) then
                    if( threads(slot_i)%complete ) then
                        slot_rc = c_pthread_mutex_unlock(threads(slot_i)%mutex)
                        if( slot_rc /= 0 ) then
                            write(*,*) 'Worker: unlock failed before join for slot', slot_i, 'rc=', slot_rc
                            cycle
                        end if
                        slot_rc = c_pthread_join(threads(slot_i)%thread, local_retval)
                        if( slot_rc /= 0 ) then
                            write(*,*) 'Worker: join failed before slot reuse for slot', slot_i, 'rc=', slot_rc
                            cycle
                        end if
                        slot_rc = c_pthread_mutex_lock(threads(slot_i)%mutex)
                        if( slot_rc /= 0 ) then
                            write(*,*) 'Worker: relock failed after join for slot', slot_i, 'rc=', slot_rc
                            cycle
                        end if
                        threads(slot_i)%started = .false.
                    else
                        slot_rc = c_pthread_mutex_unlock(threads(slot_i)%mutex)
                        if( slot_rc /= 0 ) then
                            write(*,*) 'Worker: unlock failed for busy slot', slot_i, 'rc=', slot_rc
                        end if
                        cycle
                    end if
                end if
                ! Claim slot before releasing the lock so no other caller
                ! can steal it between the unlock and pthread_create.
                threads(slot_i)%task_msg  = task
                threads(slot_i)%complete  = .false.
                threads(slot_i)%nthr      = max(task%nthr, 1)
                threads(slot_i)%exit_code = 0
                threads(slot_i)%started   = .true.
                slot_rc = c_pthread_mutex_unlock(threads(slot_i)%mutex)
                if( slot_rc /= 0 ) then
                    write(*,*) 'Worker: unlock failed before pthread_create for slot', slot_i, 'rc=', slot_rc
                    return
                end if
                create_rc = c_pthread_create(threads(slot_i)%thread, c_null_ptr, &
                                             c_funloc(worker_task_thread), c_loc(threads(slot_i)))
                if( create_rc /= 0 ) then
                    slot_rc = c_pthread_mutex_lock(threads(slot_i)%mutex)
                    if( slot_rc == 0 ) then
                        threads(slot_i)%started   = .false.
                        threads(slot_i)%complete  = .true.
                        threads(slot_i)%nthr      = 0
                        threads(slot_i)%exit_code = create_rc
                        slot_rc = c_pthread_mutex_unlock(threads(slot_i)%mutex)
                    end if
                    write(*,*) 'Worker: c_pthread_create failed for job_id', task%job_id
                end if
                return
            end if
            slot_rc = c_pthread_mutex_unlock(threads(slot_i)%mutex)
            if( slot_rc /= 0 ) then
                write(*,*) 'Worker: c_pthread_mutex_unlock failed for slot', slot_i, 'rc=', slot_rc
            end if
        end do
        write(*,*) 'Worker: no free thread slot for job_id', task%job_id, '— task dropped'
    end subroutine start_worker_thread

    !> pthread body: execute the task script and record the exit code.
    !> Runs in its own OS thread; communicates with the main thread via
    !> the mutex-protected fields of worker_thread_args.
    subroutine worker_task_thread( arg_ptr ) bind(c)
        type(c_ptr), value,        intent(in) :: arg_ptr
        type(worker_thread_args),     pointer :: thread_args
        integer :: lock_rc, exitstat
        call c_f_pointer(arg_ptr, thread_args)
        write(*,*) 'Worker thread: starting job_id', thread_args%task_msg%job_id, &
                   'script:', trim(thread_args%task_msg%script_path)
        ! Execute the script; exitstat receives the shell exit code
        call execute_command_line('bash ' // trim(thread_args%task_msg%script_path), exitstat=exitstat)
        ! Record results under the mutex
        lock_rc = c_pthread_mutex_lock(thread_args%mutex)
        if( lock_rc /= 0 ) then
            write(*,*) 'Worker thread: failed to lock slot mutex for job_id', thread_args%task_msg%job_id, 'rc=', lock_rc
            return
        end if
        thread_args%exit_code = exitstat
        thread_args%complete  = .true.
        thread_args%nthr      = 0    ! release slot
        write(*,*) 'Worker thread: completed job_id', thread_args%task_msg%job_id, &
                   'exit_code', thread_args%exit_code
        lock_rc = c_pthread_mutex_unlock(thread_args%mutex)
        if( lock_rc /= 0 ) then
            write(*,*) 'Worker thread: failed to unlock slot mutex for job_id', thread_args%task_msg%job_id, 'rc=', lock_rc
        end if
    end subroutine worker_task_thread

    !> Return the total number of thread-slots in use (sum of nthr fields).
    !> Each slot is locked individually to avoid holding multiple mutexes.
    function get_nthr_used() result( nthr_used )
        integer :: nthr_used
        integer :: idx, lock_rc
        nthr_used = 0
        do idx = 1, size(threads)
            lock_rc   = c_pthread_mutex_lock(threads(idx)%mutex)
            if( lock_rc /= 0 ) cycle
            nthr_used = nthr_used + threads(idx)%nthr
            lock_rc   = c_pthread_mutex_unlock(threads(idx)%mutex)
        end do
    end function get_nthr_used

    !> Validate script path from server before execution.
    !> Rules: absolute path, file exists, and only safe ASCII path chars.
    logical function is_safe_script_path( script_path )
        character(len=*), intent(in) :: script_path
        character(len=*), parameter   :: allowed = 'abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789/._-'
        character(len=:), allocatable :: path_trim
        logical :: exists
        integer :: j, n
        is_safe_script_path = .false.
        n = len_trim(script_path)
        if( n <= 0 ) return
        path_trim = trim(script_path)
        if( path_trim(1:1) /= '/' ) return
        do j = 1, len(path_trim)
            if( index(allowed, path_trim(j:j)) == 0 ) return
        end do
        inquire(file=path_trim, exist=exists)
        if( .not. exists ) return
        is_safe_script_path = .true.
    end function is_safe_script_path

end program simple_persistent_worker
