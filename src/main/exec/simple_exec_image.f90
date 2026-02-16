!@descr: execution of image processing commanders
module simple_exec_image
use simple_cmdline,            only: cmdline
use simple_commanders_imgops,  only: commander_binarize, commander_normalize, commander_scale
use simple_commanders_imgproc, only: commander_ctfops, commander_ctf_phaseflip
use simple_commanders_stkops,  only: commander_convert, commander_stack, commander_stackops
implicit none

public :: exec_image_commander
private

type(commander_binarize)      :: xbinarize
type(commander_convert)       :: xconvert
type(commander_ctf_phaseflip) :: xctf_phaseflip
type(commander_ctfops)        :: xctfops
type(commander_normalize)     :: xnormalize
type(commander_scale)         :: xscale
type(commander_stack)         :: xstack
type(commander_stackops)      :: xstackops

contains

    subroutine exec_image_commander(which, cline, l_silent, l_did_execute)
        character(len=*),    intent(in)    :: which
        class(cmdline),      intent(inout) :: cline
        logical,             intent(inout) :: l_did_execute
        logical,             intent(out)   :: l_silent
        if( l_did_execute )return
        l_silent      = .false.
        l_did_execute = .true.
        select case(trim(which))
            case( 'binarize' )
                call xbinarize%execute(cline)
            case( 'convert' )
                call xconvert%execute(cline)
            case( 'ctf_phaseflip' )
                call xctf_phaseflip%execute(cline)
            case( 'ctfops' )
                call xctfops%execute(cline)
            case( 'normalize' )
                call xnormalize%execute(cline)
            case( 'scale' )
                call xscale%execute(cline)
            case( 'stack' )
                call xstack%execute(cline)
            case( 'stackops' )
                call xstackops%execute(cline)
            case default
                l_did_execute = .false.
        end select
    end subroutine exec_image_commander

end module simple_exec_image
