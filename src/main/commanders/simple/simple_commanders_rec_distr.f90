!@descr: volume assembly commander: restores the even, odd and merged state volumes from partial reconstructions
module simple_commanders_rec_distr
use simple_commanders_api
use simple_volassemble_strategy, only: run_volassemble
implicit none
private
public :: commander_volassemble
#include "simple_local_flags.inc"

type, extends(commander_base) :: commander_volassemble
  contains
    procedure :: execute => exec_volassemble
end type commander_volassemble

contains

    subroutine exec_volassemble( self, cline )
        class(commander_volassemble), intent(inout) :: self
        class(cmdline),               intent(inout) :: cline
        call run_volassemble(cline)
    end subroutine exec_volassemble

end module simple_commanders_rec_distr
