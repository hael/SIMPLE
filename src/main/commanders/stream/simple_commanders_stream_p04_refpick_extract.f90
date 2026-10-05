!@descr: task 4 in the stream pipeline: reference-based picking and extraction
!==============================================================================
! MODULE: simple_commanders_stream_p04_refpick_extract
!
! PURPOSE:
!   The commander of stream stage 4. It normalises the command line, then
!   drives a stream_stage_refpick until the stream stops; everything the stage
!   does lives in simple_stream_stage_refpick.
!
! ENTRY POINT:
!   commander_stream_p04_refpick_extract%execute(cline)
!
! SIGNALS:
!   SIGTERM only sets a flag (simple_stream_sigterm); the loop sees it before
!   the next iteration and the stage finalises normally. The default action is
!   restored before execute returns.
!==============================================================================
module simple_commanders_stream_p04_refpick_extract
use simple_defs,                 only: logfhandle, PICK_LP_DEFAULT
use simple_error,                only: simple_exception
use simple_jiffys,               only: simple_end
use simple_cmdline,              only: cmdline
use simple_commander_base,       only: commander_base
use simple_stream_stage_refpick, only: stream_stage_refpick
use simple_stream_sigterm,       only: install_sigterm_handler, restore_sigterm_handler, sigterm_received
implicit none

public :: commander_stream_p04_refpick_extract
public :: set_refpick_cline ! the stage's command-line defaults, for the chained stream tests
private
#include "simple_local_flags.inc"

type, extends(commander_base) :: commander_stream_p04_refpick_extract
  contains
    procedure :: execute => exec_stream_p04_refpick_extract
end type commander_stream_p04_refpick_extract

contains

    subroutine exec_stream_p04_refpick_extract( self, cline )
        class(commander_stream_p04_refpick_extract), intent(inout) :: self
        class(cmdline),                    intent(inout) :: cline
        type(stream_stage_refpick) :: stage
        call install_sigterm_handler()
        call set_refpick_cline(cline)
        call stage%new(cline)
        do
            if( sigterm_received() )then
                write(logfhandle,'(A)') 'SIGTERM RECEIVED (REFERENCE PICKING)'
                exit
            endif
            if( stage%finished() ) exit
            call stage%iterate()
        enddo
        write(logfhandle,'(A)') '>>> TERMINATING PROCESS'
        call stage%finalize()
        call stage%kill()
        call restore_sigterm_handler()
        call simple_end('**** SIMPLE_STREAM_PICK_EXTRACT NORMAL STOP ****')
    end subroutine exec_stream_p04_refpick_extract

    ! Everything the stage needs on its command line before params%new: fixed settings,
    ! picking and extraction defaults (threads and parts: the master).
    subroutine set_refpick_cline( cline )
        class(cmdline), intent(inout) :: cline
        if( .not. cline%defined('pickrefs') ) THROW_HARD('pickrefs must be defined')
        ! fixed for this stage
        call cline%set('oritype', 'mic')
        call cline%set('mkdir',   'yes')
        call cline%set('picker',  'new')
        call cline%set('numlen',  5)
        call cline%set('stream',  'yes')
        if( .not. cline%defined('outdir')         ) call cline%set('outdir',         '')
        if( .not. cline%defined('walltime')       ) call cline%set('walltime',       29*60) ! 29 minutes
        ! picking; the molecular diameter comes from the picking references
        if( .not. cline%defined('lp_pick')        ) call cline%set('lp_pick',        PICK_LP_DEFAULT)
        if( .not. cline%defined('pick_roi')       ) call cline%set('pick_roi',       'yes')
        if( .not. cline%defined('backgr_subtr')   ) call cline%set('backgr_subtr',   'no')
        if( .not. cline%defined('thres')          ) call cline%set('thres',          0.0)
        if( cline%defined('moldiam') )then
            call cline%delete('moldiam')
            write(logfhandle,'(A)') '>>> MOLDIAM IGNORED'
        endif
        ! extraction
        if( .not. cline%defined('pcontrast')      ) call cline%set('pcontrast',      'black')
        if( .not. cline%defined('extractfrommov') ) call cline%set('extractfrommov', 'no')
    end subroutine set_refpick_cline

end module simple_commanders_stream_p04_refpick_extract
