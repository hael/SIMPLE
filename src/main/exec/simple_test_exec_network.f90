!@descr: execution of test network processing commanders
module simple_test_exec_network
use simple_cmdline,                 only: cmdline
use simple_commanders_test_network, only: commander_test_socket_client, &
                                          commander_test_socket_comm_distr, commander_test_socket_io, &
                                          commander_test_socket_server
implicit none

public :: exec_test_network_commander
private

type(commander_test_socket_client)     :: xsocket_client
type(commander_test_socket_comm_distr) :: xsocket_comm_distr
type(commander_test_socket_io)         :: xsocket_io
type(commander_test_socket_server)     :: xsocket_server

contains

    subroutine exec_test_network_commander( which, cline, l_silent, l_did_execute )
        character(len=*),    intent(in)    :: which
        class(cmdline),      intent(inout) :: cline
        logical,             intent(inout) :: l_did_execute
        logical,             intent(out)   :: l_silent
        if( l_did_execute )return
        l_silent      = .false.
        l_did_execute = .true.
        select case(trim(which))
            case( 'socket_client' )
                call xsocket_client%execute(cline)
            case( 'socket_comm_distr' )
                call xsocket_comm_distr%execute(cline)
            case( 'socket_io' )
                call xsocket_io%execute(cline)
            case( 'socket_server' )
                call xsocket_server%execute(cline)
            case default
                l_did_execute = .false.
        end select
    end subroutine exec_test_network_commander

end module simple_test_exec_network
