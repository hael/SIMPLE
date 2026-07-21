!@descr: diffusion-map 3D variability commanders
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
        use simple_flex_eigenvol_strategy, only: flex_eigenvol_strategy, create_flex_eigenvol_strategy
        class(commander_flex_eigenvol), intent(inout) :: self
        class(cmdline),                 intent(inout) :: cline
        type(parameters) :: params
        type(builder)    :: build
        class(flex_eigenvol_strategy), allocatable :: strategy
        strategy=create_flex_eigenvol_strategy(cline)
        call strategy%initialize(params,build,cline)
        call strategy%execute(params,build,cline)
        call strategy%finalize_run(params,build,cline)
        call strategy%cleanup(params)
        deallocate(strategy)
        call build%kill_general_tbox
        if( cline%defined('part') )then
            call simple_end('**** SIMPLE_FLEX_EIGENVOL WORKER NORMAL STOP ****',print_simple=.false.)
        else
            call simple_end('**** SIMPLE_FLEX_EIGENVOL NORMAL STOP ****')
        endif
    end subroutine exec_flex_eigenvol

end module simple_commanders_flex
