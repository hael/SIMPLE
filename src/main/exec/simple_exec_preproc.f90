!@descr: execution of preprocessing commanders
module simple_exec_preproc
use simple_cmdline,                only: cmdline
use simple_commanders_starproject, only: commander_assign_optics_groups
use simple_commanders_pick,        only: commander_pick, commander_extract, commander_reextract
use simple_commanders_preprocess,  only: commander_preprocess, commander_motion_correct,&
commander_gen_pspecs_and_thumbs, commander_ctf_estimate
implicit none

public :: exec_preproc_commander
private

type(commander_assign_optics_groups)  :: xassign_optics_groups
type(commander_ctf_estimate)          :: xctf_estimate
type(commander_extract)               :: xextract
type(commander_gen_pspecs_and_thumbs) :: xgen_pspecs_and_thumbs
type(commander_motion_correct)        :: xmotion_correct
type(commander_pick)                  :: xpick
type(commander_preprocess)            :: xpreprocess
type(commander_reextract)             :: xreextract

contains

    subroutine exec_preproc_commander(which, cline, l_silent, l_did_execute)
        character(len=*),    intent(in)    :: which
        class(cmdline),      intent(inout) :: cline
        logical,             intent(inout) :: l_did_execute
        logical,             intent(out)   :: l_silent
        if( l_did_execute )return
        l_silent      = .false.
        l_did_execute = .true.
        select case(trim(which))
            case( 'assign_optics_groups' )
                call xassign_optics_groups%execute(cline)
            case( 'ctf_estimate' )
                call xctf_estimate%execute(cline)
            case( 'extract' )
                call xextract%execute(cline)
            case( 'gen_pspecs_and_thumbs' )
                call xgen_pspecs_and_thumbs%execute(cline)
            case( 'motion_correct' )
                call xmotion_correct%execute(cline)
            case( 'pick' )
                call xpick%execute(cline)
            case( 'preprocess' )
                call xpreprocess%execute(cline)
            case( 'reextract' )
                call xreextract%execute(cline)
            case default
                l_did_execute = .false.
        end select
    end subroutine exec_preproc_commander

end module simple_exec_preproc