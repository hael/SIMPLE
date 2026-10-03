!@descr: SIGTERM for the stream stages: a handler that only records the request, polled by the commanders and long stage steps
!==============================================================================
! MODULE: simple_stream_sigterm
!
! PURPOSE:
!   The master stops a stage by sending SIGTERM to its process. The handler
!   only sets a flag (no I/O and no exit inside the handler); the commander's
!   loop, and stage steps that can run for minutes, poll sigterm_received()
!   and stop between steps, so no file is left half-written and no new job is
!   started once a stop is requested.
!
! USE:
!   install_sigterm_handler() at the start of a commander's execute and
!   restore_sigterm_handler() before it returns, so an in-process caller (a
!   test) is not left with a handler that nothing polls. The master passes
!   also_sigint=.true.: Ctrl-C then stops the stream in order, as SIGTERM does.
!==============================================================================
module simple_stream_sigterm
use, intrinsic :: iso_c_binding, only: c_funptr, c_null_funptr
use unix,                        only: SIGTERM, SIGINT, c_signal
implicit none

public :: install_sigterm_handler, restore_sigterm_handler, sigterm_received
private

logical, volatile :: l_sigterm = .false. ! set asynchronously by on_sigterm
logical           :: l_sigint  = .false. ! SIGINT is routed to on_sigterm too

contains

    subroutine install_sigterm_handler( also_sigint )
        logical, optional, intent(in) :: also_sigint
        l_sigterm = .false.
        l_sigint  = .false.
        if( present(also_sigint) ) l_sigint = also_sigint
        call signal(SIGTERM, on_sigterm)
        if( l_sigint ) call signal(SIGINT, on_sigterm)
    end subroutine install_sigterm_handler

    !> SIGTERM (and SIGINT, when routed) terminates the process again (SIG_DFL).
    subroutine restore_sigterm_handler()
        type(c_funptr) :: previous
        previous = c_signal(SIGTERM, c_null_funptr)
        if( l_sigint ) previous = c_signal(SIGINT, c_null_funptr)
        l_sigint = .false.
    end subroutine restore_sigterm_handler

    logical function sigterm_received()
        sigterm_received = l_sigterm
    end function sigterm_received

    subroutine on_sigterm()
        l_sigterm = .true.
    end subroutine on_sigterm

end module simple_stream_sigterm
