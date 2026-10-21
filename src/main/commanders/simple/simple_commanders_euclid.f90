!@descr: for sigma2 calculations in objfun=euclid 2D and 3D refinement
module simple_commanders_euclid
use simple_commanders_api
use simple_strategy_sigma2, only: run_commit_sigma2
implicit none
private :: run_commit_sigma2
#include "simple_local_flags.inc"

type, extends(commander_base) :: commander_calc_pspec
  contains
    procedure :: execute      => exec_calc_pspec
end type commander_calc_pspec

type, extends(commander_base) :: commander_commit_sigma2
  contains
    procedure :: execute      => exec_commit_sigma2
end type commander_commit_sigma2

contains

    subroutine exec_calc_pspec( self,cline )
        use simple_calc_pspec_strategy, only: calc_pspec_strategy, create_calc_pspec_strategy
        class(commander_calc_pspec), intent(inout) :: self
        class(cmdline), intent(inout) :: cline
        type(parameters) :: params
        type(builder)    :: build
        class(calc_pspec_strategy), allocatable :: strategy
        call cline%set('stream','no')
        if( .not. cline%defined('mkdir')   ) call cline%set('mkdir',  'yes')
        if( .not. cline%defined('oritype') ) call cline%set('oritype','ptcl3D')
        ! Create and run strategy (shared-memory/worker vs distributed master)
        strategy = create_calc_pspec_strategy(cline)
        call strategy%initialize(params, build, cline)
        call strategy%execute(params, build, cline)
        call strategy%finalize_run(params, build, cline)
        call strategy%cleanup(params, build, cline)
        ! Unified termination point
        call simple_end(strategy%end_msg, print_simple=strategy%end_print_simple)
        if( allocated(strategy) ) deallocate(strategy)
    end subroutine exec_calc_pspec

    subroutine exec_commit_sigma2( self, cline )
        class(commander_commit_sigma2), intent(inout) :: self
        class(cmdline),                 intent(inout) :: cline
        call run_commit_sigma2(cline)
    end subroutine exec_commit_sigma2

end module simple_commanders_euclid
