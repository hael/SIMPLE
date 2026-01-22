program simple_test_socket
use, intrinsic :: iso_c_binding
use simple_socket_comm
use json_kinds
use json_module
use unix, only : c_pthread_t, c_pthread_mutex_t 
use unix, only : c_pthread_create, c_pthread_join
use unix, only : c_pthread_mutex_init, c_pthread_mutex_destroy
use unix, only : c_pthread_mutex_lock, c_pthread_mutex_unlock
implicit none
type(simple_socket)                    :: socket
integer                                :: fd
write(*,*) "Socket server test. Waiting for client message"
call socket%open
call socket%set_options
call socket%bind_any
call socket%listen
do
   call socket%accept(fd)
   call socket%read(fd)
   call socket%close(fd)
end do
call socket%close
end program simple_test_socket
