!==============================================================================
! MODULE: simple_persistent_worker_server
!
! PURPOSE:
!   Provides the TCP accept-loop server that coordinates SIMPLE worker
!   processes.  Workers connect periodically with heartbeat messages;
!   the server replies with a task to execute or a TERMINATE command.
!
!   The module exposes:
!     - persistent_worker_server  — the server object (TCP socket + task queues)
!     - TCP_BUFSZ                  — wire buffer size constant used by both sides
!
! DESIGN CONTRACT:
!   1.  Mutex ownership: the mutex inside listener_args is owned by
!       persistent_worker_server (init in new(), destroy in kill()), not by
!       ipc_tcp_socket, which only borrows the pointer.
!   2.  job_count and task%job_id assignment for WORKER_NEW_TASK_MSG happen
!       in worker_listener_thread, where task queues are thread-local.
!   3.  self%port is zeroed only after ipc_socket_server%kill() returns (i.e.
!       after the listener thread has been joined), so is_running() stays
!       true until shutdown is complete.
!   4.  pack_task_queue() is called only when a task has actually been
!       dispatched and shutdown is not in progress, to avoid unnecessary
!       queue scans on idle heartbeats.
!   5.  Network replies (repl_msg) are issued OUTSIDE the mutex to
!       prevent blocking indefinitely while holding the lock.
!   6.  Debug logging (heartbeat lines, per-task skip messages, idle
!       replies) is controlled by the DEBUG compile-time constant.
!       Set DEBUG=.true. in the parameter declaration to enable verbose output.
!
! DEPENDENCIES:
!   simple_core_module_api          — logfhandle, int2str
!   simple_string                   — string type
!   simple_persistent_worker_message_* — all wire message types and type enumerator
!   simple_ipc_tcp_socket_server     — ipc_tcp_socket_server, listener_args, repl_msg
!   unix                            — pthreads, c_usleep, c_null_ptr
!   iso_c_binding                   — C interop types and procedures
!==============================================================================
module simple_persistent_worker_server
    use simple_core_module_api
    use simple_string,                              only: string
    use simple_persistent_worker_message_types,     only: WORKER_STATUS_MSG, WORKER_TERMINATE_MSG, WORKER_HEARTBEAT_MSG, WORKER_NEW_TASK_MSG
    use simple_persistent_worker_message_task,      only: qsys_persistent_worker_message_task
    use simple_persistent_worker_message_status,    only: qsys_persistent_worker_message_status
    use simple_persistent_worker_message_heartbeat, only: qsys_persistent_worker_message_heartbeat
    use simple_persistent_worker_message_terminate, only: qsys_persistent_worker_message_terminate
    use simple_ipc_tcp_socket_server,               only: ipc_tcp_socket_server, listener_args, repl_msg, c_pollfd, poll_fds
    use simple_ipc_tcp_socket_client,               only: ipc_tcp_socket_client
    use unix,                                       only: c_pthread_mutex_lock, c_pthread_mutex_unlock,  &
                                                    c_pthread_mutex_init, c_pthread_mutex_destroy, &
                                                    c_usleep, c_null_ptr, c_accept, c_socklen_t, c_close, c_read
    use simple_ipc_tcp_socket_helpers,              only: fd_is_healthy, POLLIN
    use iso_c_binding,                              only: c_ptr, c_int, c_char, c_short, c_size_t, c_associated, c_f_pointer, c_loc, c_funloc
    implicit none

    public :: persistent_worker
    public :: persistent_worker_server
    public :: TCP_BUFSZ
    private

    ! Named constants for status reply codes sent to workers
    integer, parameter :: WORKER_STATUS_IDLE  = 0       !< server has no pending task
    integer, parameter :: WORKER_STATUS_ERROR = 2       !< message was malformed or invalid
    integer, parameter :: TASK_QUEUE_SIZE     = 1000    !< max tasks per priority level; arbitrary large value
    integer, parameter :: TCP_BUFSZ           = 1460    !< TCP MTU (1500) - IP header (20) - TCP header (20)
    integer, parameter :: KILL_WAIT_TIME_US   = 2000000 !< wait time after kill() before closing socket (us)
    integer, parameter :: POLL_TIMEOUT_MS     = 100     !< timeout for poll() in milliseconds
    integer, parameter :: SWEEP_PERIOD        = 1       !< sweep dead fds every N poll wakeups (1 = every wakeup)
    logical, parameter :: DEBUG               = .false. !< set to .true. to enable verbose debug logging in the listener thread

    !> Module-level singleton tracking the persistent worker pool.
    !! Only accessed from the main (non-listener) thread; no mutex needed.
    type persistent_worker_runtime
        integer                                 :: n_workers       = 0  !< number of worker slots launched
        integer                                 :: nthr_per_worker = 0  !< thread count per worker process
        type(string)                            :: launch_backend       !< base scheduler used to launch workers
        type(persistent_worker_server), pointer :: server => null()     !< owning server object
    end type persistent_worker_runtime

    type(persistent_worker_runtime) :: persistent_worker

    !> Shared state accessed by both the caller thread (via queue_task)
    !> and the listener pthread — always under mutex protection.
    type persistent_worker_data
        integer :: n_workers   = 0       !< highest worker_id registered so far
        logical :: l_terminate = .false. !< set .true. by kill(); listener sends TERMINATE to all workers
    end type persistent_worker_data

    !> TCP accept-loop server that receives worker heartbeats and dispatches
    !> tasks to workers.  Owns the mutex in listener_args; start with new(),
    !> stop with kill().
    type persistent_worker_server
        integer                               :: port          = 0        !< TCP listen port; 0 when not running
        integer                               :: n_workers     = 0 
        integer                               :: nthr_workers  = 0        !< total worker thread slots to support
        integer                               :: job_count     = 0        !< monotonic counter for unique job_id assignment
        type(string)                          :: host_ips                 !< comma-separated local IP addresses
        type(ipc_tcp_socket_server)           :: ipc_socket_server        !< underlying TCP socket and listener thread
        type(ipc_tcp_socket_client)           :: ipc_socket_client
        type(listener_args),          pointer :: listener_args => null()  !< C-struct forwarded to the listener pthread
        type(persistent_worker_data), pointer :: worker_data   => null()  !< shared state (workers + queues); mutex-protected
    contains
        procedure :: new          !< initialise server, bind socket, start listener thread
        procedure :: kill         !< signal workers, join thread, destroy mutex, free memory
        procedure :: queue_task   !< submit a task request to the listener thread
        procedure :: get_port     !< return the TCP port being listened on
        procedure :: get_host_ips !< return the local IP address list
        procedure :: is_running   !< .true. while the listener thread is active
    end type persistent_worker_server

contains

    !> Initialise the server: allocate shared state, initialise the mutex,
    !> bind a TCP socket, and spawn the listener pthread.  No-op if already running.
    !> \param[in] nthr_workers  number of worker thread slots to support (must be > 0)
    subroutine new( self, n_workers, nthr_workers, client_only )
        class(persistent_worker_server), intent(inout) :: self
        integer,                            intent(in) :: n_workers
        integer,                            intent(in) :: nthr_workers
        type(string), optional,             intent(in) :: client_only
        integer(kind=c_int)                            :: rc
        character(len=STDLEN)                          :: addr_char, host_part, port_part
        integer                                        :: sep, ios
        if( n_workers < 1 ) then
            write(logfhandle,'(A,I0)') '>>> PERSISTENT_WORKER_SERVER new: n_workers must be > 0; got ', n_workers
            return
        end if
        if( nthr_workers < 1 ) then
            write(logfhandle,'(A,I0)') '>>> PERSISTENT_WORKER_SERVER new: nthr_workers must be > 0; got ', nthr_workers
            return
        end if
        if( associated(self%worker_data) .or. associated(self%listener_args) ) then
            if( DEBUG ) write(logfhandle,'(A)') '>>> PERSISTENT_WORKER_SERVER new: listener already running; no-op'
            return
        end if
        allocate(self%listener_args)
        allocate(self%worker_data)
        ! Initialise the mutex here — ipc_tcp_socket borrows listener_args by pointer
        ! and must not init/destroy the mutex itself (symmetric with kill()).
        rc = c_pthread_mutex_init(self%listener_args%mutex, c_null_ptr)
        if( rc /= 0 ) then
            write(logfhandle,'(A,I0)') '>>> PERSISTENT_WORKER_SERVER new: c_pthread_mutex_init failed, rc=', rc
            deallocate(self%listener_args)
            deallocate(self%worker_data)
            nullify(self%listener_args)
            nullify(self%worker_data)
            return
        end if
        self%worker_data%n_workers  = n_workers      
        self%n_workers              = n_workers
        self%nthr_workers           = nthr_workers
        self%listener_args%data_ptr = c_loc(self%worker_data)
        if( present(client_only) ) then
            ! client_only is expected as "hostlist:port" from qsys_env%get_persistent_worker_server_address.
            addr_char = trim(client_only%to_char())
            sep       = index(addr_char, ':', back=.true.)
            if( sep <= 1 .or. sep >= len_trim(addr_char) ) then
                write(logfhandle,'(A,A)') '>>> PERSISTENT_WORKER_SERVER new: invalid client_only address: ', trim(addr_char)
                self%port     = 0
                self%host_ips = ''
            else
                host_part = adjustl(addr_char(:sep-1))
                port_part = adjustl(addr_char(sep+1:len_trim(addr_char)))
                read(port_part, *, iostat=ios) self%port
                if( ios /= 0 .or. self%port <= 0 ) then
                    write(logfhandle,'(A,A)') '>>> PERSISTENT_WORKER_SERVER new: invalid client_only port in address: ', trim(addr_char)
                    self%port     = 0
                    self%host_ips = ''
                else
                    self%host_ips = trim(host_part)
                    write(logfhandle,'(A,A,A,I0)') '>>> PERSISTENT_WORKER_SERVER new: client-only mode using server ', &
                        trim(host_part), ':', self%port
                end if
            end if
        else
            call self%ipc_socket_server%new(c_funloc(worker_listener_thread), c_loc(self%listener_args))
            write(logfhandle,'(A)') '>>> PERSISTENT_WORKER: SERVER STARTED'
            self%port     = self%ipc_socket_server%get_port()
            self%host_ips = self%ipc_socket_server%get_server_ips()
        end if
        call self%ipc_socket_client%new(self%host_ips, self%port)
        write(logfhandle,'(A)') '>>> PERSISTENT_WORKER: DISPATCH CLIENT CONNECTED'
    end subroutine new

    !> Gracefully stop the server: signal workers with TERMINATE, wait
    !> KILL_WAIT_TIME_US microseconds, join the listener thread, destroy
    !> the mutex, and release all heap allocations.
    subroutine kill( self )
        class(persistent_worker_server), intent(inout) :: self
        integer(kind=c_int)                      :: rc
        self%nthr_workers = 0
        self%job_count    = 0
        ! host_ips%kill is a no-op when the string was never allocated, so this is safe
        ! to call unconditionally before the null-pointer guard below.
        call self%host_ips%kill()
        if( .not. associated(self%listener_args) .or. .not. associated(self%worker_data) ) then
            ! Listener was never started; clean up what we can.
            call self%ipc_socket_server%kill()
            self%port = 0
            return
        end if
        ! Signal the listener thread; it will forward TERMINATE to each worker.
        rc = c_pthread_mutex_lock(self%listener_args%mutex)
        if( rc /= 0 ) write(logfhandle,'(A,I0)') '>>> PERSISTENT_WORKER_SERVER kill: mutex_lock failed, rc=', rc
        self%worker_data%l_terminate = .true.
        rc = c_pthread_mutex_unlock(self%listener_args%mutex)
        if( rc /= 0 ) write(logfhandle,'(A,I0)') '>>> PERSISTENT_WORKER_SERVER kill: mutex_unlock failed, rc=', rc
        ! Wait to allow workers to receive TERMINATE before we close the socket.
        rc = c_usleep(KILL_WAIT_TIME_US)
        call self%ipc_socket_server%kill()
        self%port = 0  ! zero only after listener thread has been joined
        call self%ipc_socket_client%kill()  ! close the short-lived connection
        ! Destroy the mutex here — ipc_tcp_socket borrows the mutex by pointer;
        ! the owner (this type) must destroy it only after the listener thread has
        ! been joined (which ipc_socket_server%kill() guarantees above).
        rc = c_pthread_mutex_destroy(self%listener_args%mutex)
        if( rc /= 0 ) write(logfhandle,'(A,I0)') '>>> PERSISTENT_WORKER_SERVER kill: mutex_destroy failed, rc=', rc
        deallocate(self%worker_data)
        nullify(self%worker_data)
        deallocate(self%listener_args)
        nullify(self%listener_args)
    end subroutine kill

    !> Submit \p task to the listener thread over a short-lived TCP connection.
    !> Current wire path enqueues into the listener's normal-priority queue.
    !> Returns .false. if no valid status reply is received or if the listener
    !> rejects the task.
    function queue_task( self, task, priority ) result( queued )
        class(persistent_worker_server),      intent(inout) :: self
        type(qsys_persistent_worker_message_task), intent(inout) :: task
        type(string),                   intent(in)    :: priority
        logical                                       :: queued
        type(qsys_persistent_worker_message_status)        :: status_reply
        character(len=:), allocatable                :: snd_buffer
        character(kind=c_char),                      target :: rcv_buffer(TCP_BUFSZ)
        integer             :: nread, expected_status_bytes
        logical             :: sent
        queued        = .false.
        task%msg_type = WORKER_NEW_TASK_MSG  ! set the message type for the listener to identify this as a new task submission from the master
        if( DEBUG ) write(logfhandle,'(A,A)') '>>> PERSISTENT_WORKER_SERVER queue_task: requested priority ', priority%to_char()
        call task%serialise(snd_buffer)  ! serialize the task into a wire buffer for sending to the listener thread
        call self%ipc_socket_client%send_recv_msg(snd_buffer, rcv_buffer, sent, nread)  ! send the task as a message and wait for the reply
        if( .not. sent .or. nread < 1 ) then
            write(logfhandle,'(A)') '>>> PERSISTENT_WORKER_SERVER queue_task: no reply from listener; task dropped'
            return
        end if
        expected_status_bytes = storage_size(status_reply) / 8
        if( nread < expected_status_bytes ) then
            write(logfhandle,'(A,I0,A,I0)') '>>> PERSISTENT_WORKER_SERVER queue_task: short status reply nread=', nread, &
                    ' expected>=', expected_status_bytes
            return
        end if
        ! Listener replies with a status message; validate type before consuming status.
        status_reply = transfer(rcv_buffer, status_reply)
        if( status_reply%msg_type /= WORKER_STATUS_MSG ) then
            write(logfhandle,'(A,I0)') '>>> PERSISTENT_WORKER_SERVER queue_task: unexpected reply msg_type=', status_reply%msg_type
            return
        end if
        queued = ( status_reply%status == WORKER_STATUS_IDLE )
        if( .not. queued ) &
            write(logfhandle,'(A,I0)') '>>> PERSISTENT_WORKER_SERVER queue_task: listener rejected task, status=', status_reply%status
    end function queue_task

    !> Return the TCP port the server is listening on, or 0 if not running.
    function get_port( self ) result ( port )
        class(persistent_worker_server), intent(in)  :: self
        integer                                :: port
        port = self%port
    end function get_port

    !> Return the comma-separated list of local IP addresses for this host.
    function get_host_ips( self ) result ( host_ips )
        class(persistent_worker_server), intent(in)  :: self
        type(string)                           :: host_ips
        host_ips = self%host_ips
    end function get_host_ips

    !> Return .true. while the listener thread is active (port > 0).
    function is_running( self ) result ( running )
        class(persistent_worker_server), intent(in)  :: self
        logical                                :: running
        running = self%port > 0
    end function is_running


    ! ------------------------------------------------------------------
    ! Listener pthread body (bind(c) — spawned by ipc_tcp_socket)
    ! ------------------------------------------------------------------

    !> Accept-loop pthread.  Receives messages from workers, updates shared
    !> persistent_worker_data, and dispatches tasks or TERMINATE replies.
    !>
    !> Critical section strategy:
    !>   - Lock args%mutex to read/update status and select a task.
    !>   - Copy the selected task to local dispatch_task.
    !>   - Unlock BEFORE calling repl_msg so the mutex is never held
    !>     during a blocking network write.
    subroutine worker_listener_thread( arg_ptr ) bind(c)
        type(c_ptr), value,                     intent(in) :: arg_ptr
        type(listener_args),                       pointer :: args
        type(persistent_worker_data),              pointer :: status
        type(qsys_persistent_worker_message_heartbeat)     :: heartbeat_msg
        type(qsys_persistent_worker_message_heartbeat), allocatable     :: workers(:)
        type(qsys_persistent_worker_message_status)        :: status_msg
        type(qsys_persistent_worker_message_task)          :: dispatch_task
        type(qsys_persistent_worker_message_task)          :: tasks_priority_high(TASK_QUEUE_SIZE)
        type(qsys_persistent_worker_message_task)          :: tasks_priority_norm(TASK_QUEUE_SIZE)
        type(qsys_persistent_worker_message_task)          :: tasks_priority_low(TASK_QUEUE_SIZE)
        type(qsys_persistent_worker_message_terminate)     :: terminate_msg
        character(kind=c_char, len=TCP_BUFSZ),      target :: buf
        character(len=:),                      allocatable :: send_buffer
        type(c_pollfd),                        allocatable :: fds(:), fds_tmp(:)
        integer(kind=c_int)                                :: rc, conn_fd
        integer(kind=c_size_t)                             :: nr
        integer                                            :: nread, buffer_type, i, j, worker_id, nready, nfds
        integer                                            :: n_high, n_norm, n_low, nthr_avail, sweep_counter, job_count
        logical                                            :: ok, replied, shutdown, has_task, send_terminate

        ! ------------------------------------------------------------------
        ! Initialise Fortran pointers from C void* arguments
        ! ------------------------------------------------------------------
        if( .not. c_associated(arg_ptr) ) return
        call c_f_pointer(arg_ptr, args)
        if( c_associated(args%data_ptr) ) call c_f_pointer(args%data_ptr, status)
        if( .not. associated(status) ) then
            ! Misconfigured call-site; signal ready so the parent does not hang.
            write(logfhandle,'(A)') '>>> PERSISTENT_WORKER_SERVER worker_listener_thread: status pointer null — aborting thread'
            call signal_ready(args)
            return
        end if
        if( status%n_workers < 1 ) then
            ! Misconfigured call-site; signal ready so the parent does not hang.
            write(logfhandle,'(A,I0)') '>>> PERSISTENT_WORKER_SERVER worker_listener_thread: invalid n_workers = ', status%n_workers
            call signal_ready(args)
            return
        end if
        allocate(fds(status%n_workers + 10))  ! +1 for server_fd, +1 for kill-sentinel connection from ipc_socket_server%kill()
        allocate(workers(status%n_workers))
        fds(1)%fd     = args%fd
        fds(1)%events = POLLIN
        nfds          = 1
        shutdown      = .false.
        n_high        = 0
        n_norm        = 0
        n_low         = 0
        sweep_counter = 0
        job_count     = 0
        allocate( character(len=TCP_BUFSZ) :: send_buffer )  ! pre-alloc; serialise reallocates if a different size is needed
        ! Signal the parent thread that we are ready to accept connections.
        call signal_ready(args)

        ! ------------------------------------------------------------------
        ! Accept loop: one iteration per poll wakeup
        ! ------------------------------------------------------------------
        accept_loop: do
            ! Wait for activity on any fd.
            call poll_fds(fds, nfds, POLL_TIMEOUT_MS, nready)
            if( nready == 0 ) cycle accept_loop  ! timeout — no activity
            if( nready < 0 ) then
                write(logfhandle,'(A,I0)') '>>> PERSISTENT_WORKER_SERVER worker_listener_thread: c_poll failed, rc=', nready
                exit accept_loop
            end if

                ! Rate-limited dead-fd sweep: only every SWEEP_PERIOD wakeups to amortise
                ! the fd_is_healthy syscall cost (~poll+recv per idle fd).
            sweep_counter = sweep_counter + 1

            ! ---- New incoming connection on the listening socket ----------
            if( iand(fds(1)%revents, POLLIN) /= 0_c_short ) then
                if( nfds < size(fds) ) then
                    conn_fd = c_accept(fds(1)%fd, c_null_ptr, int(0, c_socklen_t))
                    if( conn_fd < 0 ) then
                        write(logfhandle,'(A,I0)') '>>> PERSISTENT_WORKER_SERVER worker_listener_thread: c_accept failed, rc=', conn_fd
                    else
                        nfds = nfds + 1
                        fds(nfds)%fd      = conn_fd
                        fds(nfds)%events  = POLLIN
                        fds(nfds)%revents = 0_c_short
                        call clear_stale_registry_entry()  ! evict any stale workers entry using this fd
                        if( DEBUG ) write(logfhandle,'(A,I0)') '>>> PERSISTENT_WORKER_SERVER worker_listener_thread: accepted new connection, fd=', conn_fd
                    end if
                else
                    ! Poll set full — grow the fds array by doubling, then accept.
                    allocate(fds_tmp(size(fds) * 2))
                    fds_tmp(1:size(fds)) = fds
                    call move_alloc(fds_tmp, fds)
                    if( DEBUG ) write(logfhandle,'(A,I0)') '>>> PERSISTENT_WORKER_SERVER worker_listener_thread: grew fds array to ', size(fds)
                    conn_fd = c_accept(fds(1)%fd, c_null_ptr, int(0, c_socklen_t))
                    if( conn_fd < 0 ) then
                        write(logfhandle,'(A,I0)') '>>> PERSISTENT_WORKER_SERVER worker_listener_thread: c_accept failed after grow, rc=', conn_fd
                    else
                        nfds = nfds + 1
                        fds(nfds)%fd      = conn_fd
                        fds(nfds)%events  = POLLIN
                        fds(nfds)%revents = 0_c_short
                        call clear_stale_registry_entry()
                        if( DEBUG ) write(logfhandle,'(A,I0)') '>>> PERSISTENT_WORKER_SERVER worker_listener_thread: accepted new connection, fd=', conn_fd
                    end if
                end if
            end if

            ! ---- Process readable client fds (skip server fd at index 1) --
            i = 2
            do while( i <= nfds )
                if( iand(fds(i)%revents, POLLIN) == 0_c_short ) then
                    i = i + 1
                    cycle
                end if
                conn_fd = fds(i)%fd
                replied = .false.

                ! Read directly from the already-accepted fd.
                nr = c_read(conn_fd, c_loc(buf), int(TCP_BUFSZ, c_size_t))
                if( nr == 0 ) then  ! EOF: peer closed connection gracefully
                    if( DEBUG ) write(logfhandle,'(A,I0)') '>>> PERSISTENT_WORKER_SERVER worker_listener_thread: EOF on fd ', conn_fd
                    rc = c_close(conn_fd)
                    fds(i) = fds(nfds)
                    nfds   = nfds - 1
                    cycle  ! recheck slot i (now holds what was fds(nfds))
                end if
                if( nr > int(TCP_BUFSZ, c_size_t) ) then  ! ssize_t -1 wraps to HUGE (error)
                    write(logfhandle,'(A,I0)') '>>> PERSISTENT_WORKER_SERVER worker_listener_thread: c_read error on fd ', conn_fd
                    rc = c_close(conn_fd)
                    fds(i) = fds(nfds)
                    nfds   = nfds - 1
                    cycle  ! recheck slot i (now holds what was fds(nfds))
                end if
                nread     = int(nr)
                worker_id = -1

                buffer_type = transfer(buf, buffer_type)
                if( DEBUG ) write(logfhandle,'(A,I0)') '>>> PERSISTENT_WORKER_SERVER worker_listener_thread: received message of type ', buffer_type
                select case(buffer_type)

                    ! ---- Server-side TERMINATE (from kill()) ------------------
                    case(WORKER_TERMINATE_MSG)
                        call handle_terminate_msg()
                        exit  ! exit while loop; shutdown flag exits accept_loop below

                    ! ---- Heartbeat from a remote worker -----------------------
                    case(WORKER_HEARTBEAT_MSG)
                        call handle_heartbeat_msg()

                    case(WORKER_NEW_TASK_MSG)
                        call handle_new_task_msg()
                    ! ---- Unknown message type ---------------------------------
                    case default
                        call handle_unknown_msg()
                end select

                ! If nothing was sent yet, reply with an idle status message
                if( .not. replied ) call send_idle_reply()
                i = i + 1
            end do

            if( sweep_counter >= SWEEP_PERIOD ) then
                sweep_counter = 0
                i = 2
                do while( i <= nfds )
                    if( iand(fds(i)%revents, POLLIN) /= 0_c_short ) then
                        i = i + 1   ! has data — skip health check, handle below
                    else if( .not. fd_is_healthy(fds(i)%fd) ) then
                        if( DEBUG ) write(logfhandle,'(A,I0)') '>>> PERSISTENT_WORKER_SERVER worker_listener_thread: closing dead fd ', fds(i)%fd
                        rc = c_close(fds(i)%fd)
                        fds(i) = fds(nfds)   ! compact: overwrite with last entry
                        nfds   = nfds - 1
                    else
                        i = i + 1
                    end if
                end do
            end if

            if( shutdown .and. DEBUG ) write(logfhandle,'(A)') '>>> PERSISTENT_WORKER_SERVER worker_listener_thread: shutdown flag set — exiting accept loop'
            if( shutdown ) exit accept_loop
        end do accept_loop

        ! Close any remaining accepted client sockets before releasing the poll set.
        do i = 2, nfds
            rc = c_close(fds(i)%fd)
        end do

        if( allocated(fds_tmp) ) deallocate(fds_tmp)
        deallocate(fds)
        deallocate(workers)
        deallocate(send_buffer)

    contains
        ! ------------------------------------------------------------------
        ! Internal helpers
        ! ------------------------------------------------------------------

        !> Compact \p queue in-place; only scans the first \p n_active entries.
        !> Updates \p n_active to the new active count on return.
        subroutine pack_task_queue( queue, n_active )
            type(qsys_persistent_worker_message_task), intent(inout) :: queue(:)
            integer,                                   intent(inout) :: n_active
            integer :: ii, jj
            jj = 0
            do ii = 1, n_active
                if( queue(ii)%job_id > 0 .and. queue(ii)%end_time == 0 ) then
                    jj = jj + 1
                    if( jj /= ii ) queue(jj) = queue(ii)
                end if
            end do
            do ii = jj + 1, n_active
                queue(ii) = qsys_persistent_worker_message_task()
            end do
            n_active = jj
        end subroutine pack_task_queue

        !> Return the index of the first pending, unsubmitted task in \p queue
        !> (among the first \p n_active entries) that fits within \p nthr_avail threads.
        !> Returns 0 if no eligible task is found.
        integer function find_first_dispatchable( queue, n_active, nthr_avail_in )
            type(qsys_persistent_worker_message_task), intent(in) :: queue(:)
            integer,                                   intent(in) :: n_active
            integer,                                   intent(in) :: nthr_avail_in
            integer :: k
            find_first_dispatchable = 0
            do k = 1, n_active
                if( queue(k)%job_id > 0 .and. .not. queue(k)%submitted .and. &
                    queue(k)%nthr <= nthr_avail_in ) then
                    find_first_dispatchable = k
                    return
                end if
            end do
        end function find_first_dispatchable

        !> Lock the listener mutex, set args%ready=1, and unlock.
        !> Centralises the three-line lock/set/unlock repeated at each early-exit path.
        subroutine signal_ready( a )
            type(listener_args), intent(inout) :: a
            integer(kind=c_int)               :: r
            r = c_pthread_mutex_lock(a%mutex)
            a%ready = 1
            r = c_pthread_mutex_unlock(a%mutex)
        end subroutine signal_ready

        !> Send a WORKER_STATUS_IDLE reply when no handler set replied=.true.
        !> All variables accessed via host association from worker_listener_thread.
        subroutine send_idle_reply()
            if( DEBUG ) then
                if( worker_id > 0 ) then
                    write(logfhandle,'(A,A)') '>>> PERSISTENT_WORKER_SERVER no tasks available for worker ', int2str(worker_id)
                else
                    write(logfhandle,'(A)') '>>> PERSISTENT_WORKER_SERVER no tasks available (non-heartbeat or unknown worker)'
                end if
            end if
            call status_msg%new()
            status_msg%status = WORKER_STATUS_IDLE
            call status_msg%serialise(send_buffer)
            call repl_msg(conn_fd, send_buffer, nread, ok, close_after_write=.false.)
            call status_msg%kill()
        end subroutine send_idle_reply

        !> Handle an unrecognised message type: log it and reply with WORKER_STATUS_ERROR.
        !> All variables accessed via host association from worker_listener_thread.
        subroutine handle_unknown_msg()
            write(logfhandle,'(A,A)') '>>> PERSISTENT_WORKER_SERVER worker_listener_thread: unknown message type: ', int2str(buffer_type)
            call status_msg%new()
            status_msg%status = WORKER_STATUS_ERROR
            call status_msg%serialise(send_buffer)
            call repl_msg(conn_fd, send_buffer, nread, ok, close_after_write=.false.)
            call status_msg%kill()
            replied = .true.
        end subroutine handle_unknown_msg

        !> Handle a WORKER_NEW_TASK_MSG sent by queue_task via a short-lived TCP client.
        !> Deserialises the task from buf, assigns a job_id, and appends it to the
        !> normal-priority queue.  No mutex is needed: the task queues and job_count are
        !> thread-local to worker_listener_thread and only ever accessed from here.
        !> Replies with WORKER_STATUS_IDLE to acknowledge receipt; the actual dispatch
        !> happens at the next worker heartbeat.
        !> All variables accessed via host association from worker_listener_thread.
        subroutine handle_new_task_msg()
            type(qsys_persistent_worker_message_task) :: new_task
            integer                                   :: slot
            new_task = transfer(buf, new_task)
            slot     = 0
            ! task queues and job_count are thread-local — no mutex needed
            job_count       = job_count + 1
            new_task%job_id = job_count
            if( n_norm < TASK_QUEUE_SIZE ) then
                n_norm = n_norm + 1
                tasks_priority_norm(n_norm) = new_task
                slot = n_norm
            end if
            if( slot > 0 ) then
                if( DEBUG ) write(logfhandle,'(A,A,A,A)') '>>> PERSISTENT_WORKER_SERVER queued job_id ', int2str(new_task%job_id), &
                        ' in normal-priority slot ', int2str(slot)
            else
                if( DEBUG ) write(logfhandle,'(A,A)') '>>> PERSISTENT_WORKER_SERVER normal-priority queue full; dropped job_id ', int2str(new_task%job_id)
            end if
            ! Reply with idle only when enqueued; queue-full returns error so queue_task can fail.
            call status_msg%new()
            if( slot > 0 ) then
                status_msg%status = WORKER_STATUS_IDLE
            else
                status_msg%status = WORKER_STATUS_ERROR
            end if
            call status_msg%serialise(send_buffer)
            call repl_msg(conn_fd, send_buffer, nread, ok, close_after_write=.false.)
            call status_msg%kill()
            replied = .true.
        end subroutine handle_new_task_msg

        !> Handle a WORKER_TERMINATE_MSG: acknowledge the kill signal sent by
        !> ipc_socket_server%kill(), reply to the sentinel connection, and set
        !> shutdown=.true. so the caller exits the do-while loop.
        !> All variables accessed via host association from worker_listener_thread.
        subroutine handle_terminate_msg()
            if( DEBUG ) write(logfhandle,'(A)') '>>> PERSISTENT_WORKER_SERVER worker_listener_thread: received TERMINATE — exiting'
            call terminate_msg%new()
            terminate_msg%msg_type = WORKER_TERMINATE_MSG
            call terminate_msg%serialise(send_buffer)
            call repl_msg(conn_fd, send_buffer, nread, ok)
            if( .not. ok ) write(logfhandle,'(A)') '>>> PERSISTENT_WORKER_SERVER worker_listener_thread: repl_msg failed for TERMINATE'
            call terminate_msg%kill()
            shutdown = .true.
        end subroutine handle_terminate_msg

        !> Scan the workers registry and clear any entry whose fd matches conn_fd.
        !> Called at accept time to evict stale entries when the OS reuses an fd number
        !> after a worker disconnect.  workers is thread-local so no mutex is needed.
        subroutine clear_stale_registry_entry()
            integer :: k
            do k = 1, size(workers)
                if( workers(k)%fd == conn_fd .and. workers(k)%worker_uid /= '' ) then
                    if( DEBUG ) write(logfhandle,'(A,I0,A,I0,A)') '>>> PERSISTENT_WORKER_SERVER: evicting stale registry entry fd=', &
                            conn_fd, ' (worker slot ', k, ')'
                    workers(k) = qsys_persistent_worker_message_heartbeat()
                end if
            end do
        end subroutine clear_stale_registry_entry

        !> Register or validate the heartbeating worker in the local workers array.
        !> Sets send_terminate=.true. if the worker_id is out of range or the
        !> worker_uid does not match an existing registration (stale slot).
        !> All variables accessed via host association from worker_listener_thread.
        subroutine update_worker_registry()
            if( worker_id >= 1 .and. worker_id <= size(workers) ) then
                if( workers(worker_id)%worker_uid /= '' .and. &
                    workers(worker_id)%worker_uid /= heartbeat_msg%worker_uid ) then
                    send_terminate = .true. ! worker identity changed (possible crash/restart at same slot)
                else
                    if( workers(worker_id)%worker_uid == '' ) then
                        write(logfhandle,'(A,I0,A)') '>>> PERSISTENT_WORKER: WORKER ', worker_id, ' CONNECTED'
                    end if
                    workers(worker_id)    = heartbeat_msg
                    workers(worker_id)%fd = conn_fd  ! record server-side fd; client-supplied value is ignored
                end if
            else
                send_terminate = .true. ! worker_id is out of range
            end if
        end subroutine update_worker_registry

        !> Process a WORKER_HEARTBEAT_MSG received on conn_fd.
        !> All variables accessed via host association from worker_listener_thread.
        !> Sets worker_id, replied; may set send_terminate, has_task.
        subroutine handle_heartbeat_msg()
            heartbeat_msg  = transfer(buf, heartbeat_msg)
            worker_id      = heartbeat_msg%worker_id
            replied        = .false.
            send_terminate = .false.
            has_task       = .false.

            ! --- Critical section: read/update status and pick a task ---
            rc = c_pthread_mutex_lock(args%mutex)

            send_terminate = status%l_terminate

            ! Update worker registry (bounds-checked).
            ! Log lines are deferred to outside the mutex (design contract point 5).
            call update_worker_registry()

            ! Select one task (highest priority first) if not terminating.
            ! dispatch_task is re-initialised here so stale values from a prior
            ! heartbeat cannot be accidentally re-sent.
            call dispatch_task%new()
            nthr_avail = heartbeat_msg%nthr_total - heartbeat_msg%nthr_used
            if( .not. send_terminate ) then
                j = find_first_dispatchable(tasks_priority_high, n_high, nthr_avail)
                if( j > 0 ) then
                    tasks_priority_high(j)%submitted = .true.
                    dispatch_task = tasks_priority_high(j)                                      
                    tasks_priority_high(j)%job_id = 0  ! mark slot inactive to avoid re-dispatch
                    has_task = .true.
                end if
                if( .not. has_task ) then
                    j = find_first_dispatchable(tasks_priority_norm, n_norm, nthr_avail)
                    if( j > 0 ) then
                        tasks_priority_norm(j)%submitted = .true.
                        dispatch_task = tasks_priority_norm(j)
                        tasks_priority_norm(j)%job_id = 0  ! mark slot inactive to avoid re-dispatch
                        has_task = .true.
                    end if
                end if
                if( .not. has_task ) then
                    j = find_first_dispatchable(tasks_priority_low, n_low, nthr_avail)
                    if( j > 0 ) then
                        tasks_priority_low(j)%submitted = .true.
                        dispatch_task = tasks_priority_low(j)
                        tasks_priority_low(j)%job_id = 0  ! mark slot inactive to avoid re-dispatch
                        has_task = .true.
                    end if
                end if
            end if

            if( has_task .and. .not. send_terminate ) then
                call pack_task_queue(tasks_priority_high, n_high)
                call pack_task_queue(tasks_priority_norm, n_norm)
                call pack_task_queue(tasks_priority_low, n_low)
            end if

            rc = c_pthread_mutex_unlock(args%mutex)
            ! --- End of critical section ---

            ! All log lines and network replies are outside the mutex
            ! (design contract: never hold lock during blocking I/O or writes).
            if( worker_id >= 1 .and. worker_id <= status%n_workers ) then
                if( DEBUG ) &
                    write(logfhandle,'(A,A,A,A,A,A,A)') '>>> PERSISTENT_WORKER_SERVER heartbeat worker ', int2str(worker_id), &
                            ' time ', int2str(heartbeat_msg%heartbeat_time), &
                            ' nthr ', int2str(heartbeat_msg%nthr_used), '/'//int2str(heartbeat_msg%nthr_total)
            else
                write(logfhandle,'(A,A)') '>>> PERSISTENT_WORKER_SERVER worker_listener_thread: ignoring out-of-range worker_id ', &
                        int2str(worker_id)
            end if

            call terminate_msg%new()
            if( send_terminate ) then
                if( DEBUG ) write(logfhandle,'(A,A)') '>>> PERSISTENT_WORKER_SERVER sending TERMINATE to worker ', int2str(worker_id)
                terminate_msg%msg_type = WORKER_TERMINATE_MSG
                call terminate_msg%serialise(send_buffer)
                call repl_msg(conn_fd, send_buffer, nread, ok)
                if( .not. ok ) write(logfhandle,'(A,A)') '>>> PERSISTENT_WORKER_SERVER repl_msg failed for TERMINATE to worker ', int2str(worker_id)
                call terminate_msg%kill()
                call dispatch_task%kill()
                replied = .true.
            else if( has_task ) then
                if( DEBUG ) write(logfhandle,'(A,A,A,A)') '>>> PERSISTENT_WORKER_SERVER dispatching job_id ', int2str(dispatch_task%job_id), &
                            ' to worker ', int2str(worker_id)
                call dispatch_task%serialise(send_buffer)
                call repl_msg(conn_fd, send_buffer, nread, ok, close_after_write=.false.)
                if( .not. ok ) write(logfhandle,'(A,A,A,A)') '>>> PERSISTENT_WORKER_SERVER repl_msg failed dispatching job_id ', &
                            int2str(dispatch_task%job_id), ' to worker ', int2str(worker_id)
                call terminate_msg%kill()
                call dispatch_task%kill()
                replied = .true.
            else
                call terminate_msg%kill()
                call dispatch_task%kill()
            end if
        end subroutine handle_heartbeat_msg

    end subroutine worker_listener_thread

end module simple_persistent_worker_server