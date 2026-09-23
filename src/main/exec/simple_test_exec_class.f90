!@descr: execution of the unit-test suite commanders (the fast gate and its umbrella)
module simple_test_exec_class
use simple_cmdline,               only: cmdline
use simple_commanders_test_class, only: commander_test_units, &
                                        commander_test_unit_core, commander_test_unit_ori, commander_test_unit_image, &
                                        commander_test_unit_numerics, commander_test_unit_project, commander_test_unit_ui, &
                                        commander_test_unit_ipc, commander_test_forked_process, &
                                        commander_test_unit_reconstruction, commander_test_lib_reconstruction, &
                                        commander_test_unit_pftc_align2D3D, &
                                        commander_test_unit_cart_align3D, commander_test_lib_cart_align3D
implicit none

public :: exec_test_class_commander
private

type(commander_test_units)          :: xunits
type(commander_test_unit_core)      :: xunit_core
type(commander_test_unit_ori)       :: xunit_ori
type(commander_test_unit_image)     :: xunit_image
type(commander_test_unit_numerics)  :: xunit_numerics
type(commander_test_unit_project)   :: xunit_project
type(commander_test_unit_ui)        :: xunit_ui
type(commander_test_unit_ipc)       :: xunit_ipc
type(commander_test_unit_reconstruction) :: xunit_reconstruction
type(commander_test_lib_reconstruction)  :: xlib_reconstruction
type(commander_test_unit_pftc_align2D3D) :: xunit_pftc_align2D3D
type(commander_test_unit_cart_align3D)   :: xunit_cart_align3D
type(commander_test_lib_cart_align3D)    :: xlib_cart_align3D
type(commander_test_forked_process) :: xforked_process

contains

    subroutine exec_test_class_commander( which, cline, l_silent, l_did_execute )
        character(len=*),    intent(in)    :: which
        class(cmdline),      intent(inout) :: cline
        logical,             intent(inout) :: l_did_execute
        logical,             intent(out)   :: l_silent
        if( l_did_execute )return
        l_silent      = .false.
        l_did_execute = .true.
        select case(trim(which))
            case( 'units' )
                call xunits%execute(cline)
            case( 'unit_core' )
                call xunit_core%execute(cline)
            case( 'unit_ori' )
                call xunit_ori%execute(cline)
            case( 'unit_image' )
                call xunit_image%execute(cline)
            case( 'unit_numerics' )
                call xunit_numerics%execute(cline)
            case( 'unit_project' )
                call xunit_project%execute(cline)
            case( 'unit_ui' )
                call xunit_ui%execute(cline)
            case( 'unit_ipc' )
                call xunit_ipc%execute(cline)
            case( 'unit_reconstruction' )
                call xunit_reconstruction%execute(cline)
            case( 'lib_reconstruction' )
                call xlib_reconstruction%execute(cline)
            case( 'unit_pftc_align2D3D' )
                call xunit_pftc_align2D3D%execute(cline)
            case( 'unit_cart_align3D' )
                call xunit_cart_align3D%execute(cline)
            case( 'lib_cart_align3D' )
                call xlib_cart_align3D%execute(cline)
            case( 'forked_process' )
                call xforked_process%execute(cline)
            case default
                l_did_execute = .false.
        end select
    end subroutine exec_test_class_commander

end module simple_test_exec_class
