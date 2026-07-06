!@descr: linear 3D variability commanders
module simple_commanders_flex
use simple_commanders_api
implicit none

#include "simple_local_flags.inc"

type, extends(commander_base) :: commander_flex_eigenvol
  contains
    procedure :: execute => exec_flex_eigenvol
end type commander_flex_eigenvol

contains

    subroutine exec_flex_eigenvol( self, cline )
        use simple_flex_eigenvol_strategy, only: run_flex_eigenvol_linear
        class(commander_flex_eigenvol), intent(inout) :: self
        class(cmdline),                 intent(inout) :: cline
        type(parameters) :: params
        type(builder)    :: build
        if( .not. cline%defined('mkdir')    ) call cline%set('mkdir',   'yes')
        if( .not. cline%defined('oritype')  ) call cline%set('oritype', 'ptcl3D')
        if( .not. cline%defined('nstates')  ) call cline%set('nstates', 1)
        if( .not. cline%defined('neigs')    ) call cline%set('neigs',   3)
        if( .not. cline%defined('maxits')   ) call cline%set('maxits',  3)
        if( .not. cline%defined('lp')       ) call cline%set('lp',      8.0)
        if( .not. cline%defined('outvol')   ) call cline%set('outvol',  'flex_eigvol_001.mrc')
        call cline%set('ml_reg', 'no')
        call build%init_params_and_build_general_tbox(cline, params, do3d=.true.)
        call run_flex_eigenvol_linear(params, build, cline)
        call build%kill_general_tbox
        call simple_end('**** SIMPLE_FLEX_EIGENVOL NORMAL STOP ****')
    end subroutine exec_flex_eigenvol

end module simple_commanders_flex
