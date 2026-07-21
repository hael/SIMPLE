!@descr: diffusion-map 3D variability commanders
module simple_commanders_flex
use simple_commanders_api
implicit none

#include "simple_local_flags.inc"

type, extends(commander_base) :: commander_flex_eigenvol
  contains
    procedure :: execute => exec_flex_eigenvol
end type commander_flex_eigenvol

type, extends(commander_base) :: commander_flex_eigenvol_reconstruct
  contains
    procedure :: execute => exec_flex_eigenvol_reconstruct
end type commander_flex_eigenvol_reconstruct

contains

    subroutine exec_flex_eigenvol( self, cline )
        use simple_flex_eigenvol_strategy, only: run_flex_eigenvol_diffmap
        class(commander_flex_eigenvol), intent(inout) :: self
        class(cmdline),                 intent(inout) :: cline
        type(parameters) :: params
        type(builder)    :: build
        if( .not. cline%defined('mkdir')    ) call cline%set('mkdir',   'yes')
        if( .not. cline%defined('oritype')  ) call cline%set('oritype', 'ptcl3D')
        if( .not. cline%defined('nstates')  ) call cline%set('nstates', 1)
        if( .not. cline%defined('neigs')    ) call cline%set('neigs',   20)
        if( .not. cline%defined('k_nn')     ) call cline%set('k_nn',    10)
        if( .not. cline%defined('nang_nbrs')) call cline%set('nang_nbrs', 1000)
        if( .not. cline%defined('lp')       ) call cline%set('lp',      8.0)
        if( .not. cline%defined('outvol')   ) call cline%set('outvol',  'flex_eigvol_001.mrc')
        call cline%set('ml_reg', 'no')
        call build%init_params_and_build_general_tbox(cline, params, do3d=.true.)
        call run_flex_eigenvol_diffmap(params, build, cline)
        call build%kill_general_tbox
        call simple_end('**** SIMPLE_FLEX_EIGENVOL NORMAL STOP ****')
    end subroutine exec_flex_eigenvol

    subroutine exec_flex_eigenvol_reconstruct( self, cline )
        use simple_flex_eigenvol_strategy, only: run_flex_eigenvol_reconstruct_worker
        class(commander_flex_eigenvol_reconstruct), intent(inout) :: self
        class(cmdline), intent(inout) :: cline
        type(parameters) :: params
        type(builder) :: build
        if( .not.cline%defined('mkdir') ) call cline%set('mkdir','no')
        if( .not.cline%defined('oritype') ) call cline%set('oritype','ptcl3D')
        if( .not.cline%defined('nstates') ) call cline%set('nstates',1)
        call cline%set('ml_reg','no')
        call build%init_params_and_build_general_tbox(cline,params,do3d=.true.)
        call run_flex_eigenvol_reconstruct_worker(params,build,cline)
        call qsys_declare_part_finished(params,string('simple_commanders_flex :: exec_flex_eigenvol_reconstruct'))
        call build%kill_general_tbox
        call simple_end('**** SIMPLE_FLEX_EIGENVOL_RECONSTRUCT NORMAL STOP ****',print_simple=.false.)
    end subroutine exec_flex_eigenvol_reconstruct

end module simple_commanders_flex
