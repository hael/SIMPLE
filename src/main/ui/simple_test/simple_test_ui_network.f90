!@descr: module defining the user interfaces for network test programs in the simple_test_exec suite
module simple_test_ui_network
use simple_ui_modules
implicit none

type(ui_program), target :: socket_client
type(ui_program), target :: socket_comm_distr
type(ui_program), target :: socket_io
type(ui_program), target :: socket_server

contains

    subroutine construct_test_network_programs( tsttab )
        class(ui_hash), intent(inout) :: tsttab
        call new_socket_client(tsttab)
        call new_socket_comm_distr(tsttab)
        call new_socket_io(tsttab)
        call new_socket_server(tsttab)
    end subroutine construct_test_network_programs

    subroutine print_test_network_programs( logfhandle )
        integer, intent(in) :: logfhandle
        write(logfhandle,'(A)') format_str('NETWORK:', C_UNDERLINED)
        write(logfhandle,'(A)') socket_client%name%to_char()
        write(logfhandle,'(A)') socket_comm_distr%name%to_char()
        write(logfhandle,'(A)') socket_io%name%to_char()
        write(logfhandle,'(A)') socket_server%name%to_char()
        write(logfhandle,'(A)') ''
    end subroutine print_test_network_programs

    subroutine new_socket_client( tsttab )
        class(ui_hash), intent(inout) :: tsttab
        ! PROGRAM SPECIFICATION
        call socket_client%new(&
        &'socket_client',&                     ! name
        &'socket_client ',&                    ! descr_short
        &'is a test program for ',&
        &'simple_test_exec',&                  ! executable
        &.false.)                              ! requires sp_project
        ! INPUT PARAMETER SPECIFICATIONS
        ! image input/output
        !call socket_client%add_input(UI_IO, )
        ! parameter input/output
        !call socket_client%add_input(UI_IMG, )
        ! alternative inputs
        !call socket_client%add_input(UI_PARM, )
        ! search controls
        !call socket_client%add_input(UI_SRCH, )
        ! filter controls
        !call socket_client%add_input(UI_FILT, )
        ! mask controls
        !call socket_client%add_input(UI_MASK, )
        ! computer controls
        !call socket_client%add_input(UI_COMP, )
        ! add to ui_hash
        call add_ui_program('socket_client', socket_client, tsttab)
    end subroutine new_socket_client

    subroutine new_socket_comm_distr( tsttab )
        class(ui_hash), intent(inout) :: tsttab
        ! PROGRAM SPECIFICATION
        call socket_comm_distr%new(&
        &'socket_comm_distr',&                 ! name
        &'socket_comm_distr ',&                ! descr_short
        &'is a test program for ',&
        &'simple_test_exec',&                  ! executable
        &.false.)                              ! requires sp_project
        ! INPUT PARAMETER SPECIFICATIONS
        ! image input/output
        !call socket_comm_distr%add_input(UI_IO, )
        ! parameter input/output
        !call socket_comm_distr%add_input(UI_IMG, )
        ! alternative inputs
        !call socket_comm_distr%add_input(UI_PARM, )
        ! search controls
        !call socket_comm_distr%add_input(UI_SRCH, )
        ! filter controls
        !call socket_comm_distr%add_input(UI_FILT, )
        ! mask controls
        !call socket_comm_distr%add_input(UI_MASK, )
        ! computer controls
        !call socket_comm_distr%add_input(UI_COMP, )
        ! add to ui_hash
        call add_ui_program('socket_comm_distr', socket_comm_distr, tsttab)
    end subroutine new_socket_comm_distr

    subroutine new_socket_io( tsttab )
        class(ui_hash), intent(inout) :: tsttab
        ! PROGRAM SPECIFICATION
        call socket_io%new(&
        &'socket_io',&                         ! name
        &'socket_io ',&                        ! descr_short
        &'is a test program for ',&
        &'simple_test_exec',&                  ! executable
        &.false.)                              ! requires sp_project
        ! INPUT PARAMETER SPECIFICATIONS
        ! image input/output
        !call socket_io%add_input(UI_IO, )
        ! parameter input/output
        !call socket_io%add_input(UI_IMG, )
        ! alternative inputs
        !call socket_io%add_input(UI_PARM, )
        ! search controls
        !call socket_io%add_input(UI_SRCH, )
        ! filter controls
        !call socket_io%add_input(UI_FILT, )
        ! mask controls
        !call socket_io%add_input(UI_MASK, )
        ! computer controls
        !call socket_io%add_input(UI_COMP, )
        ! add to ui_hash
        call add_ui_program('socket_io', socket_io, tsttab)
    end subroutine new_socket_io

    subroutine new_socket_server( tsttab )
        class(ui_hash), intent(inout) :: tsttab
        ! PROGRAM SPECIFICATION
        call socket_server%new(&
        &'socket_server',&                     ! name
        &'socket_server ',&                    ! descr_short
        &'is a test program for ',&
        &'simple_test_exec',&                  ! executable
        &.false.)                              ! requires sp_project
        ! INPUT PARAMETER SPECIFICATIONS
        ! image input/output
        !call socket_server%add_input(UI_IO, )
        ! parameter input/output
        !call socket_server%add_input(UI_IMG, )
        ! alternative inputs
        !call socket_server%add_input(UI_PARM, )
        ! search controls
        !call socket_server%add_input(UI_SRCH, )
        ! filter controls
        !call socket_server%add_input(UI_FILT, )
        ! mask controls
        !call socket_server%add_input(UI_MASK, )
        ! computer controls
        !call socket_server%add_input(UI_COMP, )
        ! add to ui_hash
        call add_ui_program('socket_server', socket_server, tsttab)
    end subroutine new_socket_server

end module simple_test_ui_network
