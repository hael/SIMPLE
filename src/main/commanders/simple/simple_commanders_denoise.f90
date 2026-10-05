!@descr: class expansion commander (cls_expansion)
module simple_commanders_denoise
use simple_commanders_api
implicit none

#include "simple_local_flags.inc"

type, extends(commander_base) :: commander_cls_expansion
  contains
    procedure :: execute => exec_cls_expansion
end type commander_cls_expansion

contains

    subroutine exec_cls_expansion( self, cline )
        use simple_cls_expansion_strategy
        use simple_commanders_refine2D,  only: commander_refine2D
        use simple_commanders_euclid,    only: commander_calc_pspec
        use simple_sp_project,           only: sp_project
        class(commander_cls_expansion), intent(inout) :: self
        class(cmdline),             intent(inout) :: cline
        class(cls_expansion_strategy), allocatable :: strategy
        type(parameters)         :: params
        type(builder)            :: build
        type(commander_refine2D)   :: xrefine2D
        type(commander_calc_pspec) :: xcalc_pspec
        type(cmdline)            :: cline_c2d, cline_pspec, cline_restore
        type(sp_project)         :: spproj
        type(string)             :: sigma2_path
        integer :: nsub, funit
        logical :: l_sigma2
        ! 1. the split: per-class covariance model, placement, weighted sub-averages
        strategy = create_cls_expansion_strategy(cline)
        call strategy%initialize(params, build, cline)
        call strategy%execute(params, build, cline)
        call strategy%finalize_run(params, build, cline)
        call strategy%cleanup(params)
        if( allocated(strategy) ) deallocate(strategy)
        if( cline%defined('part') )then
            call build%kill_general_tbox
            call simple_end('**** SIMPLE_CLS_EXPANSION NORMAL STOP ****')
            return
        endif
        call build%kill_general_tbox
        ! 2. one greedy in-plane round of refine2D against the sub-averages, classes fixed; the
        !    sigma2 state is estimated first when the project has none that resolves
        call spproj%read_segment('projinfo', params%projfile)
        call spproj%get_sigma2_state_path(sigma2_path, l_sigma2)
        if( l_sigma2 ) l_sigma2 = file_exists(sigma2_path%to_char())
        call spproj%read_segment('cls2D', params%projfile)
        nsub = spproj%os_cls2D%get_noris()
        call spproj%kill
        if( .not. l_sigma2 )then
            write(logfhandle,'(A)') '>>> CLS_EXPANSION: no usable sigma2 state in the project; estimating it with calc_pspec'
            cline_pspec = cline
            call cline_pspec%set('prg',   'calc_pspec')
            call cline_pspec%set('mkdir', 'no')
            call cline_pspec%delete('class')
            call cline_pspec%delete('neigs')
            call cline_pspec%delete('lp')
            call xcalc_pspec%execute(cline_pspec)
            call cline_pspec%kill
        endif
        cline_c2d = cline
        call cline_c2d%set('prg',      'refine2D')
        call cline_c2d%set('projfile', params%projfile%to_char())
        call cline_c2d%set('refs',     FLEX_CLS_CAVGS_FILE)
        call cline_c2d%set('ncls',     nsub)
        call cline_c2d%set('refine',   'inpl')
        call cline_c2d%set('maxits',   1)
        call cline_c2d%set('startit',  3)
        call cline_c2d%set('trs',      3.0)
        call cline_c2d%set('mkdir',    'no')
        call cline_c2d%delete('class')
        call cline_c2d%delete('neigs')
        call cline_c2d%delete('lp')
        write(logfhandle,'(A)') '>>> CLS_EXPANSION: one greedy in-plane round of refine2D against the sub-averages'
        call xrefine2D%execute(cline_c2d)
        call cline_c2d%kill
        call sigma2_path%kill
        ! 3. the restoration from the refined alignment (labels from the project, via the mode file)
        open(newunit=funit, file=FLEX_CLS_LABELS_MODE, status='replace', action='write')
        write(funit,'(A)') 'project'
        close(funit)
        cline_restore = cline
        call cline_restore%set('mkdir', 'no')
        write(logfhandle,'(A)') '>>> CLS_EXPANSION: restoring the sub-averages from the refined alignment'
        strategy = create_cls_expansion_strategy(cline_restore)
        call strategy%initialize(params, build, cline_restore)
        call strategy%execute(params, build, cline_restore)
        call strategy%finalize_run(params, build, cline_restore)
        call strategy%cleanup(params)
        if( allocated(strategy) ) deallocate(strategy)
        call cline_restore%kill
        write(logfhandle,'(A)') '>>> CLS_EXPANSION: weighted sub-class averages in cls_expansion_cavgs.mrc (plus even/odd), &
            &responsibilities in cls_expansion_weights.txt'
        call build%kill_general_tbox
        call simple_end('**** SIMPLE_CLS_EXPANSION NORMAL STOP ****')
    end subroutine exec_cls_expansion

end module simple_commanders_denoise
