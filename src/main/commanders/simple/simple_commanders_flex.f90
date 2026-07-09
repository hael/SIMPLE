!@descr: linear 3D variability commanders
module simple_commanders_flex
use simple_commanders_api
implicit none

#include "simple_local_flags.inc"

type, extends(commander_base) :: commander_flex_eigenvol
  contains
    procedure :: execute => exec_flex_eigenvol
end type commander_flex_eigenvol

type, extends(commander_base) :: commander_flex_eigenvol_mstep
  contains
    procedure :: execute => exec_flex_eigenvol_mstep
end type commander_flex_eigenvol_mstep

type, extends(commander_base) :: commander_flex_eigenvol_estep
  contains
    procedure :: execute => exec_flex_eigenvol_estep
end type commander_flex_eigenvol_estep

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

    subroutine exec_flex_eigenvol_mstep( self, cline )
        use simple_flex_eigenvol_strategy, only: run_flex_eigenvol_mstep_worker
        class(commander_flex_eigenvol_mstep), intent(inout) :: self
        class(cmdline),                       intent(inout) :: cline
        type(parameters) :: params
        type(builder)    :: build
        if( .not. cline%defined('mkdir')    ) call cline%set('mkdir',   'no')
        if( .not. cline%defined('oritype')  ) call cline%set('oritype', 'ptcl3D')
        if( .not. cline%defined('nstates')  ) call cline%set('nstates', 1)
        call cline%set('ml_reg', 'no')
        call build%init_params_and_build_general_tbox(cline, params, do3d=.true.)
        call run_flex_eigenvol_mstep_worker(params, build, cline)
        call qsys_declare_part_finished(params, string('simple_commanders_flex :: exec_flex_eigenvol_mstep'))
        call build%kill_general_tbox
        call simple_end('**** SIMPLE_FLEX_EIGENVOL_MSTEP NORMAL STOP ****', print_simple=.false.)
    end subroutine exec_flex_eigenvol_mstep

    subroutine exec_flex_eigenvol_estep( self, cline )
        use simple_flex_eigenvol_strategy, only: run_flex_eigenvol_estep_worker
        class(commander_flex_eigenvol_estep), intent(inout) :: self
        class(cmdline),                       intent(inout) :: cline
        type(parameters) :: params
        type(builder)    :: build
        if( .not. cline%defined('mkdir')    ) call cline%set('mkdir',   'no')
        if( .not. cline%defined('oritype')  ) call cline%set('oritype', 'ptcl3D')
        if( .not. cline%defined('nstates')  ) call cline%set('nstates', 1)
        call cline%set('ml_reg', 'no')
        call build%init_params_and_build_general_tbox(cline, params, do3d=.true.)
        call run_flex_eigenvol_estep_worker(params, build, cline)
        call qsys_declare_part_finished(params, string('simple_commanders_flex :: exec_flex_eigenvol_estep'))
        call build%kill_general_tbox
        call simple_end('**** SIMPLE_FLEX_EIGENVOL_ESTEP NORMAL STOP ****', print_simple=.false.)
    end subroutine exec_flex_eigenvol_estep

end module simple_commanders_flex
