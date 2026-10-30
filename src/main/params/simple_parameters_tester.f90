!@descr: parameter-constructor bindings used by focused high-level test selection
module simple_parameters_tester
use simple_string,     only: string
use simple_test_utils, only: assert_char, enter_fixture, leave_fixture, tests_failed
implicit none
private
public :: run_all_parameters_tests
real, parameter :: SMPD = 1.5

contains

    subroutine run_all_parameters_tests()
        use simple_defs, only: CWD_GLOB, CWD_GLOB_ORIG
        type(string)                  :: saved, root
        character(len=:), allocatable :: cwd_before, origin_before
        integer                       :: nfail0
        write(*,'(A)') '**** running all parameters tests ****'
        nfail0 = tests_failed
        if( allocated(CWD_GLOB) ) cwd_before = CWD_GLOB
        if( allocated(CWD_GLOB_ORIG) ) origin_before = CWD_GLOB_ORIG
        call enter_fixture('parameters', saved, root)
        call test_suite_binding()
        call leave_fixture(saved, root, nfail0)
        ! The constructor must not leave directory globals pointing into the removed fixture.
        if( allocated(CWD_GLOB) ) deallocate(CWD_GLOB)
        if( allocated(CWD_GLOB_ORIG) ) deallocate(CWD_GLOB_ORIG)
        if( allocated(cwd_before) ) CWD_GLOB = cwd_before
        if( allocated(origin_before) ) CWD_GLOB_ORIG = origin_before
    end subroutine run_all_parameters_tests

    ! The parameter constructor must retain the selector before the expensive simulation.
    subroutine test_suite_binding()
        use simple_parameters, only: parameters
        use simple_cmdline,    only: cmdline
        class(parameters), allocatable :: params
        class(cmdline), allocatable    :: cline
        write(*,'(A)') 'test_suite_binding'
        allocate(params, cline)
        call cline%set('prg', 'cont_refine3D_1jxy')
        call cline%set('mkdir', 'no')
        call cline%set('box', 32)
        call cline%set('smpd', SMPD)
        call cline%set('mskdiam', 12.)
        call cline%set('nthr', 1)
        call params%new(cline, silent=.true.)
        call assert_char('', trim(params%suite), 'suite parser: default selects the full gate')
        call cline%set('suite', 'nu_handoff')
        call params%new(cline, silent=.true.)
        call assert_char('nu_handoff', trim(params%suite), 'suite parser: focused NU handoff survives construction')
        call cline%set('suite', 'all')
        call params%new(cline, silent=.true.)
        call assert_char('all', trim(params%suite), 'suite parser: explicit full gate survives construction')
        call cline%kill()
    end subroutine test_suite_binding
end module simple_parameters_tester
