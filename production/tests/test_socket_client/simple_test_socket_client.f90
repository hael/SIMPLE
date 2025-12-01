program simple_test_socket
use simple_socket_comm
implicit none
type(simple_socket) :: socket
integer             :: fd, i
write(*,*) "Socket client test"
call socket%open
call socket%send("TEST MESSAGE FROM THE CLIENT")
call socket%close
write(*,*) "Sent message"
end program simple_test_socket
