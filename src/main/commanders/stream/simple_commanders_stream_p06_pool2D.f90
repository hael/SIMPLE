!@descr: task 6 in the stream pipeline: global 2D classification of the pooled particles from sieving
!==============================================================================
! MODULE: simple_commanders_stream_p06_pool2D
!
! PURPOSE:
!   The commander of stream stage 6 (the master launches it as
!   prg=pool2D). It normalises the command line, then drives a
!   stream_stage_pool2D until the stream stops; everything the stage does
!   lives in simple_stream_stage_pool2D.
!
! ENTRY POINT:
!   commander_stream_p06_pool2D%execute(cline)
!
! SIGNALS:
!   SIGTERM only sets a flag (simple_stream_sigterm); the loop sees it before
!   the next iteration and the stage finalises normally (final project from
!   the last complete iteration). The default action is restored before
!   execute returns.
!==============================================================================
module simple_commanders_stream_p06_pool2D
use simple_defs,                only: logfhandle
use simple_defs_fname,          only: METADATA_EXT
use simple_jiffys,              only: simple_end
use simple_cmdline,             only: cmdline
use simple_commander_base,      only: commander_base
use simple_stream_stage_pool2D, only: stream_stage_pool2D
use simple_stream_sigterm,      only: install_sigterm_handler, restore_sigterm_handler, sigterm_received
implicit none

public :: commander_stream_p06_pool2D
public :: set_pool2D_cline ! the stage's command-line defaults, for the chained stream tests
private

type, extends(commander_base) :: commander_stream_p06_pool2D
  contains
    procedure :: execute => exec_stream_p06_pool2D
end type commander_stream_p06_pool2D

contains

    subroutine exec_stream_p06_pool2D( self, cline )
        class(commander_stream_p06_pool2D), intent(inout) :: self
        class(cmdline),           intent(inout) :: cline
        type(stream_stage_pool2D) :: stage
        call install_sigterm_handler()
        call set_pool2D_cline(cline)
        call stage%new(cline)
        do
            if( sigterm_received() )then
                write(logfhandle,'(A)') 'SIGTERM RECEIVED (POOL 2D)'
                exit
            endif
            if( stage%finished() ) exit
            call stage%iterate()
        enddo
        write(logfhandle,'(A)') '>>> TERMINATING PROCESS'
        call stage%finalize()
        call stage%kill()
        call restore_sigterm_handler()
        call simple_end('**** SIMPLE_STREAM_SOLVE2D NORMAL STOP ****')
    end subroutine exec_stream_p06_pool2D

    ! Everything the stage needs on its command line before params%new: the pool's fixed 2D
    ! settings, its defaults, and the project the GUI reads when none is named.
    subroutine set_pool2D_cline( cline )
        class(cmdline), intent(inout) :: cline
        call cline%set('oritype',     'mic')
        call cline%set('mkdir',       'yes')
        call cline%set('autoscale',   'yes')
        call cline%set('reject_mics', 'no')
        call cline%set('refine',      'snhc_smpl')
        call cline%set('ml_reg',      'no')
        call cline%set('objfun',      'euclid')
        call cline%set('sigma_est',   'global')
        call cline%set('cls_init',    'rand')
        call cline%set('numlen',      5)
        if( .not. cline%defined('dynreslim')       ) call cline%set('dynreslim',       'yes')
        if( .not. cline%defined('stepwise')        ) call cline%set('stepwise',        'yes')
        if( .not. cline%defined('center')          ) call cline%set('center',          'yes')
        if( .not. cline%defined('ncls')            ) call cline%set('ncls',            200)
        if( .not. cline%defined('projfile') )then
            call cline%set('projname', 'stream_solve2D')
            call cline%set('projfile', 'stream_solve2D'//METADATA_EXT)
        endif
    end subroutine set_pool2D_cline

end module simple_commanders_stream_p06_pool2D
