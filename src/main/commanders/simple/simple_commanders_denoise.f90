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

    !> one pass: the subclass budget over the classes (ncls in total), per-class covariance model, placement (hard
    !! labels) and the conditional sub-class averages with their even/odd halves and FRCs. The poses of the averages
    !! are refined downstream (refine3D_states_cavgs), the particles' by particle-level refine3D_states
    subroutine exec_cls_expansion( self, cline )
        use simple_cls_expansion_strategy
        class(commander_cls_expansion), intent(inout) :: self
        class(cmdline),             intent(inout) :: cline
        class(cls_expansion_strategy), allocatable :: strategy
        type(parameters) :: params
        type(builder)    :: build
        strategy = create_cls_expansion_strategy(cline)
        call strategy%initialize(params, build, cline)
        call strategy%execute(params, build, cline)
        call strategy%finalize_run(params, build, cline)
        call strategy%cleanup(params)
        if( allocated(strategy) ) deallocate(strategy)
        call build%kill_general_tbox
        if( .not. cline%defined('part') )then
            write(logfhandle,'(A)') '>>> CLS_EXPANSION: sub-class averages (subcls_avg='//trim(params%subcls_avg)// &
                &') in cls_expansion_cavgs.mrc (plus even/odd), their FRCs in cls_expansion_frcs.bin, the particles'' &
                &subclasses in the project (class; parent in cluster)'
        endif
        call simple_end('**** SIMPLE_CLS_EXPANSION NORMAL STOP ****')
    end subroutine exec_cls_expansion

end module simple_commanders_denoise
