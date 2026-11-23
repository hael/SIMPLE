!@descr: execution of denoising commanders
module simple_exec_denoise
use simple_cmdline,              only: cmdline
use simple_commanders_denoise,   only: commander_cls_expansion
use simple_commanders_flex_pca,  only: commander_flex_pca
use simple_commanders_flex_tree, only: commander_flex_tree
use simple_string,               only: string
implicit none

public :: exec_denoise_commander
private

type(commander_cls_expansion) :: xcls_expansion
type(commander_flex_pca)      :: xflex_pca
type(commander_flex_tree)     :: xflex_tree

contains

     subroutine exec_denoise_commander(which, cline, l_silent, l_did_execute)
        character(len=*),    intent(in)    :: which
        class(cmdline),      intent(inout) :: cline
        logical,             intent(inout) :: l_did_execute
        logical,             intent(out)   :: l_silent
        type(string) :: mode
        if( l_did_execute )return
        l_silent      = .false.
        l_did_execute = .true.
        select case(trim(which))
            case( 'cls_expansion' )
                call xcls_expansion%execute(cline)
            case( 'flex_pca' )
                mode = 'single'
                if( cline%defined('mode') ) mode = cline%get_carg('mode')
                if( mode%to_char() == 'tree' )then
                    call xflex_tree%execute(cline)
                else
                    call xflex_pca%execute(cline)
                endif
            case default
                l_did_execute = .false.
        end select
    end subroutine exec_denoise_commander

end module simple_exec_denoise
