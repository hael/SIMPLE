program simple_test_socket
   use simple_socket_comm
   implicit none

   type(simple_socket) :: socket
   integer             :: fd, i
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
