program simple_test_socket_io
use, intrinsic :: iso_c_binding
use simple_socket_comm
use simple_string_utils
use unix, only : c_pthread_t, c_pthread_create 
implicit none
type(c_pthread_t)   :: server_thread
type(simple_socket) :: client_socket
integer             :: server_rc, i
! start server thread
server_rc = c_pthread_create(server_thread, c_null_ptr, c_funloc(socket_server), c_null_ptr)
! send a client message every 5 seconds for 5 iterations
do i=1,5
   call sleep(10)
   call client_socket%open
   call client_socket%send("TEST MESSAGE FROM THE CLIENT. ITERATION : "//int2str(i))
   call client_socket%close
end do

contains

   subroutine socket_server() bind(c)
      type(simple_socket) :: socket
      integer             :: fd
      write(*,*) "Starting socket server thread"
      call socket%open()
      call socket%bind_any()
      call socket%listen()
      write(*,*) "Socket server listening"
      do
         call socket%accept(fd)
         call socket%read(fd)
         call socket%close(fd)
      end do
      call socket%close
   end subroutine socket_server

end program simple_test_socket_io
