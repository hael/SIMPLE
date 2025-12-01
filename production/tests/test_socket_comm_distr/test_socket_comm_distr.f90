program test_socket_comm_distr
use, intrinsic :: iso_c_binding
use simple_distr_comm
implicit none
type(distr_comm) :: comm
write(*,*) "START"
call comm%init()
call sleep(10)
end program test_socket_comm_distr
