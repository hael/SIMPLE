!@ descr: module defining the user interfaces for network programs in the simple_test_exec suite
module simple_test_ui_network
use simple_ui_modules
implicit none

type(ui_program), target :: simple_test_socket_client
type(ui_program), target :: simple_test_socket_server
type(ui_program), target :: simple_test_socket_io
type(ui_program), target :: simple_test_socket_comm_distr

contains

    subroutine construct_network_programs(prgtab)
        class(ui_hash), intent(inout) :: prgtab
        call new_simple_test_socket_client(prgtab)
        call new_simple_test_socket_server(prgtab)
        call new_simple_test_socket_io(prgtab)
        call new_simple_test_socket_comm_distr(prgtab)
    end subroutine construct_network_programs

    subroutine print_network_programs(logfhandle)
        integer, intent(in) :: logfhandle
        write(logfhandle,'(A)') format_str('NETWORK', C_UNDERLINED)
        write(logfhandle,'(A)') simple_test_socket_client%name%to_char()
        write(logfhandle,'(A)') simple_test_socket_server%name%to_char()
        write(logfhandle,'(A)') simple_test_socket_io%name%to_char()
        write(logfhandle,'(A)') simple_test_socket_comm_distr%name%to_char()
        write(logfhandle,'(A)') ''
    end subroutine print_network_programs

    subroutine new_simple_test_socket_client( prgtab )
        class(ui_hash), intent(inout) :: prgtab
        ! add to ui_hash
        call add_ui_program('simple_test_socket_client', simple_test_socket_client, prgtab)
    end subroutine new_simple_test_socket_client

    subroutine new_simple_test_socket_server( prgtab )
        class(ui_hash), intent(inout) :: prgtab
        ! add to ui_hash
        call add_ui_program('simple_test_socket_server', simple_test_socket_server, prgtab)
    end subroutine new_simple_test_socket_server

    subroutine new_simple_test_socket_io( prgtab )
        class(ui_hash), intent(inout) :: prgtab
        ! add to ui_hash
        call add_ui_program('simple_test_socket_io', simple_test_socket_io, prgtab)
    end subroutine new_simple_test_socket_io

    subroutine new_simple_test_socket_comm_distr( prgtab )
        class(ui_hash), intent(inout) :: prgtab
        ! add to ui_hash
        call add_ui_program('simple_test_socket_comm_distr', simple_test_socket_comm_distr, prgtab)
    end subroutine new_simple_test_socket_comm_distr

end module simple_test_ui_network
