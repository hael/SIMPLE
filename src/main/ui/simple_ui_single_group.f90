!@descr: aggregates SINGLE ui program constructors
module simple_ui_single_group
use simple_ui_hash,       only: ui_hash
use single_ui_atom,       only: construct_single_atom_programs,       print_single_atom_programs
use single_ui_map,        only: construct_single_map_programs,        print_single_map_programs
use single_ui_nano2D,     only: construct_single_nano2D_programs,     print_single_nano2D_programs
use single_ui_nano3D,     only: construct_single_nano3D_programs,     print_single_nano3D_programs
use single_ui_trajectory, only: construct_single_trajectory_programs, print_single_trajectory_programs
use single_ui_tseries,    only: construct_single_tseries_programs,    print_single_tseries_programs
use single_ui_validate,   only: construct_single_validate_programs,   print_single_validate_programs
implicit none

public :: add_single_programs, print_single_programs
private

contains

    subroutine add_single_programs( prgtab )
        class(ui_hash), intent(inout) :: prgtab
        call construct_single_atom_programs(prgtab)
        call construct_single_map_programs(prgtab)
        call construct_single_nano2D_programs(prgtab)
        call construct_single_nano3D_programs(prgtab)
        call construct_single_trajectory_programs(prgtab)
        call construct_single_tseries_programs(prgtab)
        call construct_single_validate_programs(prgtab)
    end subroutine add_single_programs

    subroutine print_single_programs( logfhandle )
        integer, intent(in) :: logfhandle
        call print_single_tseries_programs(logfhandle)
        call print_single_trajectory_programs(logfhandle)
        call print_single_nano2D_programs(logfhandle)
        call print_single_nano3D_programs(logfhandle)
        call print_single_atom_programs(logfhandle)
        call print_single_map_programs(logfhandle)
        call print_single_validate_programs(logfhandle)
    end subroutine print_single_programs

end module simple_ui_single_group
