!@descr: diffusion-map 3D variability commanders
module simple_commanders_flex_analysis
use simple_commanders_api
implicit none

#include "simple_local_flags.inc"

type, extends(commander_base) :: commander_flex_analysis
  contains
    procedure :: execute => exec_flex_analysis
end type commander_flex_analysis

type, extends(commander_base) :: commander_flex_pca
  contains
    procedure :: execute => exec_flex_pca
end type commander_flex_pca

contains

    subroutine exec_flex_analysis( self, cline )
        use simple_flex_analysis_strategy, only: flex_analysis_strategy, create_flex_analysis_strategy
        class(commander_flex_analysis), intent(inout) :: self
        class(cmdline),                 intent(inout) :: cline
        type(parameters) :: params
        type(builder)    :: build
        class(flex_analysis_strategy), allocatable :: strategy
        strategy=create_flex_analysis_strategy(cline)
        call strategy%initialize(params,build,cline)
        call strategy%execute(params,build,cline)
        call strategy%finalize_run(params,build,cline)
        call strategy%cleanup(params)
        deallocate(strategy)
        call build%kill_general_tbox
        if( cline%defined('part') )then
            call simple_end('**** SIMPLE_FLEX_ANALYSIS WORKER NORMAL STOP ****',print_simple=.false.)
        else
            call simple_end('**** SIMPLE_FLEX_ANALYSIS NORMAL STOP ****')
        endif
    end subroutine exec_flex_analysis

    subroutine exec_flex_pca( self, cline )
        use simple_flex_pca_model, only: run_flex_pca
        use simple_flex_pca_distr, only: flex_pca_distr_init, flex_pca_distr_kill
        class(commander_flex_pca), intent(inout) :: self
        class(cmdline),                   intent(inout) :: cline
        type(parameters) :: params
        type(builder)    :: build
        ! a worker runs in the master's directory and must not descend into one of its own.
        ! NOT merge('no ','yes',..): merge pads the shorter branch, yielding 'no ' with a trailing blank.
        if( .not.cline%defined('mkdir') )then
            if( cline%defined('part') )then
                call cline%set('mkdir','no')
            else
                call cline%set('mkdir','yes')
            endif
        endif
        if( .not.cline%defined('oritype') )     call cline%set('oritype','ptcl3D')
        if( .not.cline%defined('nstates') )     call cline%set('nstates',1)
        if( .not.cline%defined('npreimages') )  call cline%set('npreimages',7)
        if( .not.cline%defined('neigs') )       call cline%set('neigs',6)
        if( .not.cline%defined('lp') )          call cline%set('lp',12.0)
        if( .not.cline%defined('box_crop') )    call cline%set('box_crop',64)
        if( .not.cline%defined('ptcl_src') )    call cline%set('ptcl_src','raw')
        if( .not.cline%defined('objfun') )      call cline%set('objfun','euclid')
        if( .not.cline%defined('outvol') )      call cline%set('outvol','flex_pca_state_001.mrc')
        call build%init_params_and_build_general_tbox(cline,params,do3d=.true.)
        call flex_pca_distr_init(params,cline)
        call run_flex_pca(params,build,cline)
        call flex_pca_distr_kill(params)
        call build%kill_general_tbox
        if( cline%defined('part') )then
            call simple_end('**** SIMPLE_FLEX_PCA WORKER NORMAL STOP ****',print_simple=.false.)
        else
            call simple_end('**** SIMPLE_FLEX_PCA NORMAL STOP ****')
        endif
    end subroutine exec_flex_pca

end module simple_commanders_flex_analysis
