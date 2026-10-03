!@descr: task 2 in the stream pipeline: assign optics groups to streamed micrographs
!==============================================================================
! MODULE: simple_commanders_stream_p02_assign_optics
!
! PURPOSE:
!   The commander of stream stage 2. It normalises the command line, then
!   drives a stream_stage_optics until the stream stops; everything the stage
!   does lives in simple_stream_stage_optics.
!
! ENTRY POINT:
!   commander_stream_p02_assign_optics%execute(cline)
!
! SIGNALS:
!   SIGTERM only sets a flag (simple_stream_sigterm); the loop sees it before
!   the next iteration and the stage finalises normally. The default action is
!   restored before execute returns.
!==============================================================================
module simple_commanders_stream_p02_assign_optics
use simple_defs,                only: logfhandle
use simple_error,               only: simple_exception
use simple_jiffys,              only: simple_end
use simple_cmdline,             only: cmdline
use simple_commander_base,      only: commander_base
use simple_stream_stage_optics, only: stream_stage_optics
use simple_stream_sigterm,      only: install_sigterm_handler, restore_sigterm_handler, sigterm_received
implicit none

public :: commander_stream_p02_assign_optics
private
#include "simple_local_flags.inc"

type, extends(commander_base) :: commander_stream_p02_assign_optics
  contains
    procedure :: execute => exec_stream_p02_assign_optics
end type commander_stream_p02_assign_optics

contains

    subroutine exec_stream_p02_assign_optics( self, cline )
        class(commander_stream_p02_assign_optics), intent(inout) :: self
        class(cmdline),                  intent(inout) :: cline
        type(stream_stage_optics) :: stage
        call install_sigterm_handler()
        call cline%printline()
        call flush(logfhandle)
        call cline%set('mkdir', 'yes')
        if( .not. cline%defined('dir_target') ) THROW_HARD('DIR_TARGET must be defined!')
        if( .not. cline%defined('outdir')     ) call cline%set('outdir', '')
        call stage%new(cline)
        do
            if( sigterm_received() )then
                write(logfhandle,'(A)') 'SIGTERM RECEIVED (OPTICS ASSIGNMENT)'
                exit
            endif
            if( stage%finished() ) exit
            call stage%iterate()
        enddo
        write(logfhandle,'(A)') '>>> TERMINATING PROCESS'
        call stage%finalize()
        call stage%kill()
        call restore_sigterm_handler()
        call simple_end('**** SIMPLE_STREAM_ASSIGN_OPTICS NORMAL STOP ****')
    end subroutine exec_stream_p02_assign_optics

end module simple_commanders_stream_p02_assign_optics
