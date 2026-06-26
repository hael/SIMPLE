!@descr: class split commander
module simple_commanders_cls_split
use simple_commanders_api
implicit none

#include "simple_local_flags.inc"

type, extends(commander_base) :: commander_cls_split
  contains
    procedure :: execute => exec_cls_split
end type commander_cls_split

type, extends(commander_base) :: commander_denoise_project
  contains
        procedure :: execute => exec_denoise_project
end type commander_denoise_project

contains

    subroutine exec_cls_split( self, cline )
        use simple_cls_split_strategy
        use simple_commanders_mkcavgs, only: run_make_cavgs_workflow
        class(commander_cls_split), intent(inout) :: self
        class(cmdline),             intent(inout) :: cline
        class(cls_split_strategy), allocatable :: strategy
        type(parameters) :: params
        type(builder)    :: build
        type(cmdline)    :: cline_make
        strategy = create_cls_split_strategy(cline)
        call strategy%initialize(params, build, cline)
        call strategy%execute(params, build, cline)
        call strategy%finalize_run(params, build, cline)
        call strategy%cleanup(params)
        if( allocated(strategy) ) deallocate(strategy)
        if( .not. cline%defined('part') )then
            if( trim(params%oritype) == 'ptcl2D' )then
                cline_make = cline
                call cline_make%set('prg',     'make_cavgs')
                call cline_make%set('projfile', params%projfile%to_char())
                if( .not. cline_make%defined('refs') ) call cline_make%set('refs', 'cls_split_cavgs.mrc')
                call cline_make%set('mkdir',   'no')
                call cline_make%set('oritype', 'ptcl2D')
                call cline_make%delete('class')
                call cline_make%delete('ncls')
                call cline_make%delete('gen_model')
                call run_make_cavgs_workflow(cline_make, from_distr_cmd=.true.)
                call cline_make%kill
            else
                THROW_WARN('cls_split class-average regeneration is only available for oritype=ptcl2D')
            endif
        endif
        call build%kill_general_tbox
        call simple_end('**** SIMPLE_CLS_SPLIT NORMAL STOP ****')
    end subroutine exec_cls_split

    subroutine exec_denoise_project( self, cline )
        use simple_denoise_project_strategy
        class(commander_denoise_project), intent(inout) :: self
        class(cmdline),                    intent(inout) :: cline
        class(denoise_project_strategy), allocatable :: strategy
        type(parameters) :: params
        type(builder)    :: build
        strategy = create_denoise_project_strategy(cline)
        call strategy%initialize(params, build, cline)
        call strategy%execute(params, build, cline)
        call strategy%finalize_run(params, build, cline)
        call strategy%cleanup(params)
        if( allocated(strategy) ) deallocate(strategy)
        call build%kill_general_tbox
        call simple_end('**** SIMPLE_DENOISE_PROJECT NORMAL STOP ****')
    end subroutine exec_denoise_project

end module simple_commanders_cls_split
