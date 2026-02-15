module simple_exec_preproc
use simple_cmdline,                only: cmdline
use simple_commanders_starproject, only: commander_assign_optics_groups
use simple_commanders_pick,        only: commander_pick_distr, commander_extract_distr, commander_reextract_distr
use simple_commanders_preprocess,  only: commander_preprocess_distr, commander_motion_correct_distr,&
commander_gen_pspecs_and_thumbs_distr, commander_ctf_estimate_distr
implicit none

public :: exec_preproc_commander
private

type(commander_assign_optics_groups)        :: xassign_optics_groups
type(commander_ctf_estimate_distr)          :: xctf_estimate_distr
type(commander_extract_distr)               :: xextract_distr
type(commander_gen_pspecs_and_thumbs_distr) :: xgen_pspecs_and_thumbs
type(commander_motion_correct_distr)        :: xmotion_correct_distr
type(commander_pick_distr)                  :: xpick_distr
type(commander_preprocess_distr)            :: xpreprocess
type(commander_reextract_distr)             :: xreextract_distr

contains

    subroutine exec_preproc_commander(which, cline, l_silent, l_did_execute)
        character(len=*),    intent(in)    :: which
        class(cmdline),      intent(inout) :: cline
        logical,             intent(out)   :: l_silent
        logical,             intent(out)   :: l_did_execute
        l_silent      = .false.
        l_did_execute = .true.
        select case(trim(which))
            case( 'assign_optics_groups' )
                call xassign_optics_groups%execute(cline)
            case( 'ctf_estimate' )
                call xctf_estimate_distr%execute(cline)
            case( 'extract' )
                call xextract_distr%execute(cline)
            case( 'gen_pspecs_and_thumbs' )
                call xgen_pspecs_and_thumbs%execute(cline)
            case( 'motion_correct' )
                call xmotion_correct_distr%execute(cline)
            case( 'pick' )
                call xpick_distr%execute(cline)
            case( 'preprocess' )
                call xpreprocess%execute(cline)
            case( 'reextract' )
                call xreextract_distr%execute(cline)
            case default
                l_did_execute = .false.
        end select
    end subroutine exec_preproc_commander

end module simple_exec_preproc