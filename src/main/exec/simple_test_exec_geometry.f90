!@descr: execution of test geometry processing commanders
module simple_test_exec_geometry
use simple_cmdline,                  only: cmdline
use simple_commanders_test_geometry, only: commander_test_angres, commander_test_ori_test, &
                                           commander_test_oris_test, commander_test_sym_test, &
                                           commander_test_uniform_euler, commander_test_uniform_rot
implicit none

public :: exec_geometry_commander
private

type(commander_test_angres)        :: xangres
type(commander_test_ori_test)      :: xori_test
type(commander_test_oris_test)     :: xoris_test
type(commander_test_sym_test)      :: xsym_test
type(commander_test_uniform_euler) :: xuniform_euler
type(commander_test_uniform_rot)   :: xuniform_rot

contains

    subroutine exec_geometry_commander(which, cline, l_silent, l_did_execute)
        character(len=*),    intent(in)    :: which
        class(cmdline),      intent(inout) :: cline
        logical,             intent(inout) :: l_did_execute
        logical,             intent(out)   :: l_silent
        if( l_did_execute )return
        l_silent      = .false.
        l_did_execute = .true.
        select case(trim(which))
            case( 'angres' )
                call xangres%execute(cline)
            case( 'ori_test' )
                call xori_test%execute(cline)
            case( 'oris_test' )
                call xoris_test%execute(cline)
            case( 'sym_test' )
                call xsym_test%execute(cline)
            case( 'uniform_euler' )
                call xuniform_euler%execute(cline)
            case( 'uniform_rot' )
                call xuniform_rot%execute(cline)
            case default
                l_did_execute = .false.
        end select
    end subroutine exec_geometry_commander

end module simple_test_exec_geometry
