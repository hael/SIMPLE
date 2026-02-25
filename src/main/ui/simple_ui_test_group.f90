!@descr: aggregates SIMPLE TEST ui program constructors
module simple_ui_test_group
use simple_ui_hash,          only: ui_hash
use simple_test_ui_class,    only: construct_test_class_programs,    print_test_class_programs
use simple_test_ui_fft,      only: construct_test_fft_programs,      print_test_fft_programs
use simple_test_ui_geometry, only: construct_test_geometry_programs, print_test_geometry_programs
use simple_test_ui_highlevel,only: construct_test_highlevel_programs,print_test_highlevel_programs
use simple_test_ui_io,       only: construct_test_io_programs,       print_test_io_programs
use simple_test_ui_masks,    only: construct_test_masks_programs,    print_test_masks_programs
use simple_test_ui_network,  only: construct_test_network_programs,  print_test_network_programs
use simple_test_ui_numerics, only: construct_test_numerics_programs, print_test_numerics_programs
use simple_test_ui_optimize, only: construct_test_optimize_programs, print_test_optimize_programs
use simple_test_ui_parallel, only: construct_test_parallel_programs, print_test_parallel_programs
use simple_test_ui_stats,    only: construct_test_stats_programs,    print_test_stats_programs
use simple_test_ui_utils,    only: construct_test_utils_programs,    print_test_utils_programs
implicit none

public :: add_test_programs, print_test_programs
private

contains

    subroutine add_test_programs( tsttab )
        class(ui_hash), intent(inout) :: tsttab
        call construct_test_class_programs(tsttab)
        call construct_test_fft_programs(tsttab)
        call construct_test_geometry_programs(tsttab)
        call construct_test_highlevel_programs(tsttab)
        call construct_test_io_programs(tsttab)
        call construct_test_masks_programs(tsttab)
        call construct_test_network_programs(tsttab)
        call construct_test_numerics_programs(tsttab)
        call construct_test_optimize_programs(tsttab)
        call construct_test_parallel_programs(tsttab)
        call construct_test_stats_programs(tsttab)
        call construct_test_utils_programs(tsttab)
    end subroutine add_test_programs

    subroutine print_test_programs( logfhandle )
        integer, intent(in) :: logfhandle
        call print_test_class_programs(logfhandle)
        call print_test_fft_programs(logfhandle)
        call print_test_geometry_programs(logfhandle)
        call print_test_highlevel_programs(logfhandle)
        call print_test_io_programs(logfhandle)
        call print_test_masks_programs(logfhandle)
        call print_test_network_programs(logfhandle)
        call print_test_numerics_programs(logfhandle)
        call print_test_optimize_programs(logfhandle)
        call print_test_parallel_programs(logfhandle)
        call print_test_stats_programs(logfhandle)
        call print_test_utils_programs(logfhandle)
    end subroutine print_test_programs

end module simple_ui_test_group
