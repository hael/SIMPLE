module simple_distr_comm
use, intrinsic :: iso_c_binding
use simple_socket_comm
use unix
implicit none

public :: distr_comm
private

type :: thread_comm
    private
    type(c_pthread_mutex_t) :: lock
    integer                 :: port      = 39000
    logical                 :: listening = .false.
    ! add a message queue - use it to communicate with shmem and distr processes
end type thread_comm

type(thread_comm), target :: server_comm

type :: distr_comm
    private
    integer           :: sock_fd, port
    type(c_pthread_t) :: server_thread
contains
    procedure :: init
end type distr_comm

contains
    
    subroutine init(this)
        class(distr_comm), intent(inout) :: this
        integer :: rc, i_trials
        ! start server thread
        rc = c_pthread_create(this%server_thread, c_null_ptr, c_funloc(server), c_loc(server_comm))
        do i_trials=1, 50
            rc = c_pthread_mutex_lock(server_comm%lock)
            if(server_comm%listening) exit
            rc = c_pthread_mutex_unlock(server_comm%lock)
            call sleep(1)
        end do
        if(server_comm%listening) then
            write(*,*) "thread_comm: server listening on port", server_comm%port
        else
            write(*,*) "thread_comm: failed to start server"
        end if
        rc = c_pthread_mutex_unlock(server_comm%lock)
    end subroutine init
   
    subroutine server(arg) bind(c)
        type(c_ptr), intent(in), value  :: arg ! Client data.
        type(simple_socket)             :: socket
        type(thread_comm), pointer :: server_comm ! Fortran pointer to client data.
        integer             :: fd, rc
        if (.not. c_associated(arg)) return
        call c_f_pointer(arg, server_comm)
        rc = c_pthread_mutex_lock(server_comm%lock)
        call socket%open()
        call socket%bind_any()
        call socket%listen()
        server_comm%listening = .true.
        rc = c_pthread_mutex_unlock(server_comm%lock)
        do
            call socket%accept(fd)
            call socket%read(fd)
            call socket%close(fd)
        end do
        call socket%close
    end subroutine server

end module simple_distr_comm