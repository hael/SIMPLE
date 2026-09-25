!@descr: aggregates SIMPLE TEST ui program constructors
module simple_ui_test_group
use simple_ui_hash,          only: ui_hash
use simple_test_ui_class,    only: construct_test_class_programs
use simple_test_ui_highlevel,only: construct_test_highlevel_programs
implicit none

public :: add_test_programs
private

contains

    subroutine add_test_programs( tsttab )
        class(ui_hash), intent(inout) :: tsttab
        call construct_test_class_programs(tsttab)
        call construct_test_highlevel_programs(tsttab)
    end subroutine add_test_programs

end module simple_ui_test_group
