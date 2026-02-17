module single_exec_map
use simple_cmdline,                 only: cmdline
use single_commanders_experimental, only: commander_tsegmaps_core_finder
use simple_commanders_atoms,        only: commander_conv_atom_denoise
implicit none

public :: exec_map_commander
private

type(commander_conv_atom_denoise)    :: xconv_atom_denoise
type(commander_tsegmaps_core_finder) :: xtsegmaps_core_finder

contains

    subroutine exec_map_commander(which, cline, l_silent, l_did_execute)
        character(len=*),    intent(in)    :: which
        class(cmdline),      intent(inout) :: cline
        logical,             intent(inout) :: l_did_execute
        logical,             intent(out)   :: l_silent
        if( l_did_execute )return
        l_silent      = .false.
        l_did_execute = .true.
        select case(trim(which))
            case( 'conv_atom_denoise' )
                call xconv_atom_denoise%execute(cline)
            case( 'tsegmaps_core_finder' )
                call xtsegmaps_core_finder%execute(cline)
            case default    
                l_did_execute = .false.
        end select
    end subroutine exec_map_commander

end module single_exec_map
