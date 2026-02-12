!@descr: "preproc" UI api (concrete implementation)
module simple_ui_api_preproc
use simple_ui_program, only: ui_program
use simple_ui_hash,    only: ui_hash
use simple_ui_utils,   only: add_ui_program
implicit none

type(ui_program), target :: assign_optics_groups
type(ui_program), target :: ctf_estimate
type(ui_program), target :: extract
type(ui_program), target :: gen_pspecs_and_thumbs
type(ui_program), target :: motion_correct
type(ui_program), target :: pick
type(ui_program), target :: preprocess
type(ui_program), target :: reextract

contains

    subroutine register_simple_ui_preproc(prgtab)
        class(ui_hash), intent(inout) :: prgtab
        call add_ui_program('assign_optics_groups',  assign_optics_groups,  prgtab)
        call add_ui_program('ctf_estimate',          ctf_estimate,          prgtab)
        call add_ui_program('extract',               extract,               prgtab)
        call add_ui_program('gen_pspecs_and_thumbs', gen_pspecs_and_thumbs, prgtab)
        call add_ui_program('motion_correct',        motion_correct,        prgtab)
        call add_ui_program('pick',                  pick,                  prgtab)
        call add_ui_program('preprocess',            preprocess,            prgtab)
        call add_ui_program('reextract',             reextract,             prgtab)
    end subroutine register_simple_ui_preproc

end module simple_ui_api_preproc
