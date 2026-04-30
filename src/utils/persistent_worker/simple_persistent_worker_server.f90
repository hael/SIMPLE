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
!   2.  job_count is incremented and task%job_id is assigned under the
!       mutex in queue_task() to avoid data races with the listener.
!   3.  self%port is zeroed only after ipc_socket%kill() returns (i.e.
!       after the listener thread has been joined), so is_running() stays
!       true until shutdown is complete.
!   4.  pack_task_queue() is called only when a task has actually been
!       dispatched and shutdown is not in progress, to avoid unnecessary
!       queue scans on idle heartbeats.
!   5.  Network replies (repl_msg) are issued OUTSIDE the mutex to
!       prevent blocking indefinitely while holding the lock.
!   6.  Debug logging (heartbeat lines, per-task skip messages, idle
!       replies) is controlled by the l_debug flag.  Set it with
!       call server%set_debug(.true.) before starting the server or at
!       any time afterwards (the setter acquires the mutex).
!
! DEPENDENCIES:
!   simple_core_module_api          — logfhandle, int2str
!   simple_string                   — string type
!   simple_persistent_worker_message_* — all wire message types and type enumerator
!   simple_ipc_tcp_socket           — ipc_tcp_socket, listener_args, recv_msg, repl_msg
!   unix                            — pthreads, c_usleep, c_null_ptr
!   iso_c_binding                   — C interop types and procedures
!==============================================================================
module simple_persistent_worker_server
    use simple_core_module_api
    use simple_string,                              only: string
    use simple_persistent_worker_message_types,     only: WORKER_TASK_MSG, WORKER_STATUS_MSG, WORKER_TERMINATE_MSG, WORKER_HEARTBEAT_MSG
    use simple_persistent_worker_message_task,      only: qsys_persistent_worker_message_task
    use simple_persistent_worker_message_status,    only: qsys_persistent_worker_message_status
    use simple_persistent_worker_message_heartbeat, only: qsys_persistent_worker_message_heartbeat
    use simple_persistent_worker_message_terminate, only: qsys_persistent_worker_message_terminate
    use simple_ipc_tcp_socket,                      only: ipc_tcp_socket, listener_args, recv_msg, repl_msg
    use unix,                                       only: c_pthread_mutex_lock, c_pthread_mutex_unlock,  &
                                                    c_pthread_mutex_init, c_pthread_mutex_destroy, &
                                                    c_usleep, c_null_ptr
    use iso_c_binding,                              only: c_ptr, c_int, c_char, c_associated, c_f_pointer, c_loc, c_funloc
    implicit none

    public :: persistent_worker
    public :: persistent_worker_server
    public :: TCP_BUFSZ
    private

    ! Named constants for status reply codes sent to workers
    integer, parameter :: WORKER_STATUS_IDLE  = 0       !< server has no pending task
    integer, parameter :: WORKER_STATUS_ERROR = 2       !< message was malformed or invalid
    integer, parameter :: MAX_WORKERS         = 100     !< max worker_id supported by the server; arbitrary large value
    integer, parameter :: TASK_QUEUE_SIZE     = 1000    !< max tasks per priority level; arbitrary large value
    integer, parameter :: TCP_BUFSZ           = 1460    !< TCP MTU (1500) - IP header (20) - TCP header (20)
    integer, parameter :: KILL_WAIT_TIME_US   = 2000000 !< wait time after kill() before closing socket (us)

    !> Module-level singleton tracking the persistent worker pool.
    !! Only accessed from the main (non-listener) thread; no mutex needed.
    type persistent_worker_runtime
        integer :: n_workers       = 0  !< number of worker slots launched
        integer :: nthr_per_worker = 0  !< thread count per worker process
        type(string)                            :: launch_backend !< base scheduler used to launch workers
        type(persistent_worker_server), pointer :: server => null() !< owning server object
    end type persistent_worker_runtime

    type(persistent_worker_runtime) :: persistent_worker

    !> Shared state accessed by both the caller thread (via queue_task)
    !> and the listener pthread — always under mutex protection.
    type persistent_worker_data
        integer :: n_workers   = 0       !< highest worker_id registered so far
        logical :: l_terminate = .false. !< set .true. by kill(); listener sends TERMINATE to all workers
        logical :: l_debug     = .false. !< when .true. the listener emits verbose heartbeat / skip / idle log lines
        !> Latest heartbeat from each worker; indexed by worker_id (1..MAX_WORKERS)
        type(qsys_persistent_worker_message_heartbeat) :: workers(MAX_WORKERS)
        !> High-priority task queue; slots are active while job_id > 0 and end_time == 0
        type(qsys_persistent_worker_message_task) :: tasks_priority_high(TASK_QUEUE_SIZE)
        !> Normal-priority task queue
        type(qsys_persistent_worker_message_task) :: tasks_priority_norm(TASK_QUEUE_SIZE)
        !> Low-priority task queue
        type(qsys_persistent_worker_message_task) :: tasks_priority_low(TASK_QUEUE_SIZE)
    end type persistent_worker_data

    !> TCP accept-loop server that receives worker heartbeats and dispatches
    !> tasks to workers.  Owns the mutex in listener_args; start with new(),
    !> stop with kill().
    type persistent_worker_server
        integer                         :: port         = 0        !< TCP listen port; 0 when not running
        integer                         :: nthr_workers = 0        !< total worker thread slots to support
        integer                         :: job_count    = 0        !< monotonic counter for unique job_id assignment
        logical                         :: l_debug      = .false.  !< enable verbose debug logging (propagated to worker_data)
        type(string)                    :: host_ips                !< comma-separated local IP addresses
        type(ipc_tcp_socket)            :: ipc_socket              !< underlying TCP socket and listener thread
        type(listener_args),    pointer :: listener_args => null() !< C-struct forwarded to the listener pthread
        type(persistent_worker_data), pointer :: worker_data   => null() !< shared state (workers + queues); mutex-protected
    contains
        procedure :: new          !< initialise server, bind socket, start listener thread
        procedure :: kill         !< signal workers, join thread, destroy mutex, free memory
        procedure :: queue_task   !< enqueue a task in the selected priority queue
        procedure :: set_debug    !< enable or disable verbose debug logging
        procedure :: get_port     !< return the TCP port being listened on
        procedure :: get_host_ips !< return the local IP address list
        procedure :: is_running   !< .true. while the listener thread is active
    end type persistent_worker_server

contains

    !> Initialise the server: allocate shared state, initialise the mutex,
    !> bind a TCP socket, and spawn the listener pthread.  No-op if already running.
    !> \param[in] nthr_workers  number of worker thread slots to support (must be > 0)
    subroutine new( self, nthr_workers )
        class(persistent_worker_server), intent(inout) :: self
        integer,                   intent(in)    :: nthr_workers
        integer(kind=c_int)                      :: rc
        if( nthr_workers <= 0 ) then
            write(logfhandle,'(A,I0)') '>>> PERSISTENT_WORKER_SERVER new: nthr_workers must be > 0; got ', nthr_workers
            return
        end if
        if( associated(self%worker_data) .or. associated(self%listener_args) ) then
            write(logfhandle,'(A)') '>>> PERSISTENT_WORKER_SERVER new: listener already running; no-op'
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
        self%listener_args%data_ptr = c_loc(self%worker_data)
        self%worker_data%l_debug = self%l_debug  ! propagate flag before spawning the thread
        self%nthr_workers = nthr_workers
        call self%ipc_socket%init_server(c_funloc(worker_listener_thread), c_loc(self%listener_args))
        self%port     = self%ipc_socket%get_port()
        self%host_ips = self%ipc_socket%get_host_ips()
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
            call self%ipc_socket%kill()
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
        call self%ipc_socket%kill()
        self%port = 0  ! zero only after listener thread has been joined
        ! Destroy the mutex here — ipc_tcp_socket borrows the mutex by pointer;
        ! the owner (this type) must destroy it only after the listener thread has
        ! been joined (which ipc_socket%kill() guarantees above).
        rc = c_pthread_mutex_destroy(self%listener_args%mutex)
        if( rc /= 0 ) write(logfhandle,'(A,I0)') '>>> PERSISTENT_WORKER_SERVER kill: mutex_destroy failed, rc=', rc
        deallocate(self%worker_data)
        nullify(self%worker_data)
        deallocate(self%listener_args)
        nullify(self%listener_args)
    end subroutine kill

    !> Enqueue \p task in the priority queue selected by \p priority
    !> ('high', 'low', or anything else for normal).
    !> Assigns a unique job_id under the mutex.
    !> Returns .false. if the server is not initialised or the selected queue is full.
    function queue_task( self, task, priority ) result( queued )
        class(persistent_worker_server),      intent(inout) :: self
        type(qsys_persistent_worker_message_task), intent(inout) :: task
        type(string),                   intent(in)    :: priority
        logical                                       :: queued
        integer(kind=c_int) :: rc
        integer             :: i
        queued = .false.
        ! Guard on both pointers: they are always set together in new(), but
        ! a partially-initialised state must not dereference listener_args.
        if( .not. associated(self%worker_data) .or. .not. associated(self%listener_args) ) then
            write(logfhandle,'(A)') '>>> PERSISTENT_WORKER_SERVER queue_task: server not initialised; dropping task'
            return
        end if
        ! All branches reach the rc = c_pthread_mutex_unlock at the end of the
        ! function — even the queue-full paths — so the lock is always released.
        rc = c_pthread_mutex_lock(self%listener_args%mutex)
        if( rc /= 0 ) then
            write(logfhandle,'(A,I0)') '>>> PERSISTENT_WORKER_SERVER queue_task: mutex_lock failed, rc=', rc
            return
        end if
        self%job_count = self%job_count + 1  ! assign job_id under the lock to avoid races
        task%job_id    = self%job_count
        select case( priority%to_char() )
            case('high')
                do i = 1, size(self%worker_data%tasks_priority_high)
                    if( self%worker_data%tasks_priority_high(i)%job_id == 0 ) then
                        self%worker_data%tasks_priority_high(i) = task
                        queued = .true.
                        write(logfhandle,'(A,A,A,A)') '>>> PERSISTENT_WORKER_SERVER queued job_id ', &
                            int2str(task%job_id), ' in high-priority slot ', int2str(i)
                        exit
                    end if
                end do
                if( .not. queued ) write(logfhandle,'(A,A)') '>>> PERSISTENT_WORKER_SERVER high-priority queue full; dropped job_id ', int2str(task%job_id)
            case('low')
                do i = 1, size(self%worker_data%tasks_priority_low)
                    if( self%worker_data%tasks_priority_low(i)%job_id == 0 ) then
                        self%worker_data%tasks_priority_low(i) = task
                        queued = .true.
                        write(logfhandle,'(A,A,A,A)') '>>> PERSISTENT_WORKER_SERVER queued job_id ', &
                            int2str(task%job_id), ' in low-priority slot ', int2str(i)
                        exit
                    end if
                end do
                if( .not. queued ) write(logfhandle,'(A,A)') '>>> PERSISTENT_WORKER_SERVER low-priority queue full; dropped job_id ', int2str(task%job_id)
            case default  ! 'norm' or anything unrecognised
                do i = 1, size(self%worker_data%tasks_priority_norm)
                    if( self%worker_data%tasks_priority_norm(i)%job_id == 0 ) then
                        self%worker_data%tasks_priority_norm(i) = task
                        queued = .true.
                        write(logfhandle,'(A,A,A,A)') '>>> PERSISTENT_WORKER_SERVER queued job_id ', &
                            int2str(task%job_id), ' in normal-priority slot ', int2str(i)
                        exit
                    end if
                end do
                if( .not. queued ) write(logfhandle,'(A,A)') '>>> PERSISTENT_WORKER_SERVER normal-priority queue full; dropped job_id ', int2str(task%job_id)
        end select
        rc = c_pthread_mutex_unlock(self%listener_args%mutex)
        if( rc /= 0 ) write(logfhandle,'(A,I0)') '>>> PERSISTENT_WORKER_SERVER queue_task: mutex_unlock failed, rc=', rc
    end function queue_task

    !> Enable or disable verbose debug logging.
    !> Propagates the flag to worker_data under the mutex so the listener
    !> thread sees the change immediately.  Safe to call at any time.
    subroutine set_debug( self, l_debug )
        class(persistent_worker_server), intent(inout) :: self
        logical,                   intent(in)    :: l_debug
        integer(kind=c_int) :: rc
        self%l_debug = l_debug
        if( associated(self%listener_args) .and. associated(self%worker_data) ) then
            rc = c_pthread_mutex_lock(self%listener_args%mutex)
            if( rc /= 0 ) write(logfhandle,'(A,I0)') '>>> PERSISTENT_WORKER_SERVER set_debug: mutex_lock failed, rc=', rc
            self%worker_data%l_debug = l_debug
            rc = c_pthread_mutex_unlock(self%listener_args%mutex)
            if( rc /= 0 ) write(logfhandle,'(A,I0)') '>>> PERSISTENT_WORKER_SERVER set_debug: mutex_unlock failed, rc=', rc
        end if
    end subroutine set_debug

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
        type(c_ptr), value,          intent(in) :: arg_ptr
        type(listener_args),            pointer :: args
        type(persistent_worker_data),         pointer :: status
        type(qsys_persistent_worker_message_heartbeat)     :: heartbeat_msg
        type(qsys_persistent_worker_message_status)        :: status_msg
        type(qsys_persistent_worker_message_task)          :: dispatch_task
        type(qsys_persistent_worker_message_terminate)     :: terminate_msg
        character(kind=c_char, len=TCP_BUFSZ), target :: buf
        character(len=:), allocatable      :: send_buffer
        integer(kind=c_int)  :: rc, conn_fd
        integer              :: nread, buffer_type, i, worker_id
        logical              :: ok, replied, send_terminate, has_task

        ! ------------------------------------------------------------------
        ! Initialise Fortran pointers from C void* arguments
        ! ------------------------------------------------------------------
        if( .not. c_associated(arg_ptr) ) return
        call c_f_pointer(arg_ptr, args)
        if( c_associated(args%data_ptr) ) call c_f_pointer(args%data_ptr, status)
        if( .not. associated(status) ) then
            ! Misconfigured call-site; signal ready so the parent does not hang.
            write(logfhandle,'(A)') '>>> PERSISTENT_WORKER_SERVER worker_listener_thread: status pointer null — aborting thread'
            rc = c_pthread_mutex_lock(args%mutex)
            args%ready = 1
            rc = c_pthread_mutex_unlock(args%mutex)
            return
        end if
        ! Signal the parent thread that we are ready to accept connections.
        rc = c_pthread_mutex_lock(args%mutex)
        args%ready = 1
        rc = c_pthread_mutex_unlock(args%mutex)

        ! ------------------------------------------------------------------
        ! Accept loop: one iteration per incoming connection
        ! ------------------------------------------------------------------
        do
            call recv_msg(args%server_fd, conn_fd, buf, nread, ok, close_after_read=.false.)
            if( .not. ok ) exit

            replied   = .false.
            worker_id = -1  ! set to real value only after a valid heartbeat

            ! Guard against zero-length or oversized messages
            if( nread <= 0 .or. nread > TCP_BUFSZ ) then
                call status_msg%new()
                status_msg%status = WORKER_STATUS_ERROR
                call status_msg%serialise(send_buffer)
                call repl_msg(conn_fd, send_buffer, nread, ok)
                call status_msg%kill()
                replied = .true.
                cycle
            end if

            buffer_type = transfer(buf, buffer_type)
            select case(buffer_type)

                ! ---- Server-side TERMINATE (from kill()) ------------------
                case(WORKER_TERMINATE_MSG)
                    write(logfhandle,'(A)') '>>> PERSISTENT_WORKER_SERVER worker_listener_thread: received TERMINATE — exiting'
                    call terminate_msg%new()
                    call terminate_msg%serialise(send_buffer)
                    call repl_msg(conn_fd, send_buffer, nread, ok)
                    call terminate_msg%kill()
                    replied = .true.
                    exit  ! break out of accept loop

                ! ---- Heartbeat from a remote worker -----------------------
                case(WORKER_HEARTBEAT_MSG)
                    heartbeat_msg  = transfer(buf, heartbeat_msg)
                    worker_id      = heartbeat_msg%worker_id
                    replied        = .false.
                    send_terminate = .false.
                    has_task       = .false.

                    ! --- Critical section: read/update status and pick a task ---
                    rc = c_pthread_mutex_lock(args%mutex)

                    send_terminate = status%l_terminate

                    ! Update worker registry (bounds-checked).
                    ! Log lines are deferred to outside the mutex (see design contract point 5).
                    if( worker_id >= 1 .and. worker_id <= size(status%workers) ) then
                        if( status%n_workers < worker_id ) status%n_workers = worker_id
                        if( status%workers(worker_id)%worker_uid /= '' .and. &
                            status%workers(worker_id)%worker_uid /= heartbeat_msg%worker_uid ) then
                            send_terminate = .true. ! worker identity changed (possible crash/restart at same slot)
                        else
                            status%workers(worker_id) = heartbeat_msg
                        end if
                    else
                        send_terminate = .true. ! worker_id is out of range
                    end if

                    ! Select one task (highest priority first) if not terminating.
                    ! dispatch_task is re-initialised here each iteration so stale
                    ! values from a prior heartbeat cannot be accidentally re-sent.
                    call dispatch_task%new()
                    if( .not. send_terminate ) then
                        do i = 1, size(status%tasks_priority_high)
                            if( has_task ) exit
                            if( status%tasks_priority_high(i)%job_id > 0 .and. &
                                .not. status%tasks_priority_high(i)%submitted ) then
                                if( heartbeat_msg%nthr_total - heartbeat_msg%nthr_used < &
                                    status%tasks_priority_high(i)%nthr ) then
                                    if( status%l_debug ) &
                                        write(logfhandle,'(A,A,A,A,A)') '>>> PERSISTENT_WORKER_SERVER high task ', &
                                        int2str(status%tasks_priority_high(i)%job_id), &
                                        ' needs ', int2str(status%tasks_priority_high(i)%nthr), &
                                        ' threads; skipping'
                                    cycle
                                end if
                                status%tasks_priority_high(i)%submitted = .true.
                                dispatch_task = status%tasks_priority_high(i)
                                has_task = .true.
                            end if
                        end do
                        do i = 1, size(status%tasks_priority_norm)
                            if( has_task ) exit
                            if( status%tasks_priority_norm(i)%job_id > 0 .and. &
                                .not. status%tasks_priority_norm(i)%submitted ) then
                                if( heartbeat_msg%nthr_total - heartbeat_msg%nthr_used < &
                                    status%tasks_priority_norm(i)%nthr ) then
                                    if( status%l_debug ) &
                                        write(logfhandle,'(A,A,A,A,A)') '>>> PERSISTENT_WORKER_SERVER norm task ', &
                                        int2str(status%tasks_priority_norm(i)%job_id), &
                                        ' needs ', int2str(status%tasks_priority_norm(i)%nthr), &
                                        ' threads; skipping'
                                    cycle
                                end if
                                status%tasks_priority_norm(i)%submitted = .true.
                                dispatch_task = status%tasks_priority_norm(i)
                                status%tasks_priority_norm(i)%job_id = 0  ! mark the slot as inactive immediately to avoid races
                                has_task = .true.
                            end if
                        end do
                        do i = 1, size(status%tasks_priority_low)
                            if( has_task ) exit
                            if( status%tasks_priority_low(i)%job_id > 0 .and. &
                                .not. status%tasks_priority_low(i)%submitted ) then
                                if( heartbeat_msg%nthr_total - heartbeat_msg%nthr_used < &
                                    status%tasks_priority_low(i)%nthr ) then
                                    if( status%l_debug ) &
                                        write(logfhandle,'(A,A,A,A,A)') '>>> PERSISTENT_WORKER_SERVER low task ', &
                                        int2str(status%tasks_priority_low(i)%job_id), &
                                        ' needs ', int2str(status%tasks_priority_low(i)%nthr), &
                                        ' threads; skipping'
                                    cycle
                                end if
                                status%tasks_priority_low(i)%submitted = .true.
                                dispatch_task = status%tasks_priority_low(i)
                                has_task = .true.
                            end if
                        end do
                    end if

                    ! Compact queues only when a task was actually dispatched to avoid
                    ! O(3*TASK_QUEUE_SIZE) work on every idle heartbeat.
                    if( has_task .and. .not. send_terminate ) then
                        call pack_task_queue(status%tasks_priority_high)
                        call pack_task_queue(status%tasks_priority_norm)
                        call pack_task_queue(status%tasks_priority_low)
                    end if

                    rc = c_pthread_mutex_unlock(args%mutex)
                    ! --- End of critical section ---

                    ! All log lines and network replies are outside the mutex
                    ! (design contract: never hold lock during blocking I/O or writes).
                    if( worker_id >= 1 .and. worker_id <= MAX_WORKERS ) then
                        if( status%l_debug ) &
                            write(logfhandle,'(A,A,A,A,A,A,A)') '>>> PERSISTENT_WORKER_SERVER heartbeat worker ', int2str(worker_id), &
                                   ' time ', int2str(heartbeat_msg%heartbeat_time), &
                                   ' nthr ', int2str(heartbeat_msg%nthr_used), '/'//int2str(heartbeat_msg%nthr_total)
                    else
                        write(logfhandle,'(A,A)') '>>> PERSISTENT_WORKER_SERVER worker_listener_thread: ignoring out-of-range worker_id ', &
                               int2str(worker_id)
                    end if

                    ! Blocking network reply is outside the mutex.
                    call terminate_msg%new()
                    if( send_terminate ) then
                        write(logfhandle,'(A,A)') '>>> PERSISTENT_WORKER_SERVER sending TERMINATE to worker ', int2str(worker_id)
                        ! Use a terminate message (not the received heartbeat) to carry the
                        ! TERMINATE discriminator; avoids sending stale heartbeat payload.
                        terminate_msg%msg_type = WORKER_TERMINATE_MSG
                        call terminate_msg%serialise(send_buffer)
                        call repl_msg(conn_fd, send_buffer, nread, ok)
                        call terminate_msg%kill()
                        call dispatch_task%kill()
                        replied = .true.
                    else if( has_task ) then
                        write(logfhandle,'(A,A,A,A)') '>>> PERSISTENT_WORKER_SERVER dispatching job_id ', int2str(dispatch_task%job_id), &
                                   ' to worker ', int2str(worker_id)
                        call dispatch_task%serialise(send_buffer)
                        call repl_msg(conn_fd, send_buffer, nread, ok)
                        call terminate_msg%kill()
                        call dispatch_task%kill()
                        replied = .true.
                    else
                        call terminate_msg%kill()
                        call dispatch_task%kill()
                    end if

                ! ---- Unknown message type ---------------------------------
                case default
                    write(logfhandle,'(A,A)') '>>> PERSISTENT_WORKER_SERVER worker_listener_thread: unknown message type: ', int2str(buffer_type)
                    ! Reply with ERROR so the sender knows the message was rejected.
                    call status_msg%new()
                    status_msg%status = WORKER_STATUS_ERROR
                    call status_msg%serialise(send_buffer)
                    call repl_msg(conn_fd, send_buffer, nread, ok)
                    call status_msg%kill()
                    replied = .true.

            end select

            ! If nothing was sent yet, reply with an idle status message
            if( .not. replied ) then
                if( status%l_debug ) then
                    if( worker_id > 0 ) then
                        write(logfhandle,'(A,A)') '>>> PERSISTENT_WORKER_SERVER no tasks available for worker ', int2str(worker_id)
                    else
                        write(logfhandle,'(A)') '>>> PERSISTENT_WORKER_SERVER no tasks available (non-heartbeat or unknown worker)'
                    end if
                end if
                call status_msg%new()
                status_msg%status = WORKER_STATUS_IDLE
                call status_msg%serialise(send_buffer)
                call repl_msg(conn_fd, send_buffer, nread, ok)
                call status_msg%kill()
            end if
        end do

    contains
        ! ------------------------------------------------------------------
        ! Internal helpers
        ! ------------------------------------------------------------------

        !> Compact \p queue in-place, shifting active entries (job_id>0 and
        !> end_time==0) to the front and zero-filling the tail.
        !> A task with end_time==0 is either pending or submitted-but-unconfirmed;
        !> retaining unconfirmed tasks makes worker crashes detectable.
        subroutine pack_task_queue( queue )
            type(qsys_persistent_worker_message_task), intent(inout) :: queue(:)
            integer :: ii, jj
            jj = 0
            do ii = 1, size(queue)
                if( queue(ii)%job_id > 0 .and. queue(ii)%end_time == 0 ) then
                    jj = jj + 1
                    if( jj /= ii ) queue(jj) = queue(ii)
                end if
            end do
            do ii = jj + 1, size(queue)
                queue(ii) = qsys_persistent_worker_message_task()  ! reset to default (job_id=0)
            end do
        end subroutine pack_task_queue

    end subroutine worker_listener_thread

end module simple_persistent_worker_server