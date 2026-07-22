!@descr: diffusion-map 3D variability commanders
module simple_commanders_flex_analysis
use simple_commanders_api
implicit none

#include "simple_local_flags.inc"

type, extends(commander_base) :: commander_flex_analysis
  contains
    procedure :: execute => exec_flex_analysis
end type commander_flex_analysis

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

end module simple_commanders_flex_analysis
