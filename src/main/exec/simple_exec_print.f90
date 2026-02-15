module simple_exec_print
use simple_cmdline,           only: cmdline
use simple_commanders_misc,   only: commander_print_fsc, commander_print_magic_boxes, commander_print_dose_weights
use simple_commanders_checks, only: commander_info_image, commander_info_stktab
implicit none

type(commander_info_image)         :: xinfo_image
type(commander_info_stktab)        :: xinfo_stktab
type(commander_print_dose_weights) :: xprint_dose_weights
type(commander_print_fsc)          :: xprint_fsc
type(commander_print_magic_boxes)  :: xprint_magic_boxes

contains

    subroutine exec_print_commander(which, cline, l_silent, l_did_execute)
        character(len=*),    intent(in)    :: which
        class(cmdline),      intent(inout) :: cline
        logical,             intent(out)   :: l_silent
        logical,             intent(out)   :: l_did_execute
        l_silent      = .false.
        l_did_execute = .true.
        select case(trim(which))
            case( 'info_image' )
                call xinfo_image%execute(cline)
            case( 'info_stktab' )
                call xinfo_stktab%execute(cline)
            case( 'print_dose_weights' )
                call xprint_dose_weights%execute(cline)
                l_silent = .true.
            case( 'print_fsc' )
                call xprint_fsc%execute(cline)
                l_silent = .true.
            case( 'print_magic_boxes' )
                call xprint_magic_boxes%execute(cline)
                l_silent = .true.
            case default
                l_did_execute = .false.
        end select
    end subroutine exec_print_commander

end module simple_exec_print
