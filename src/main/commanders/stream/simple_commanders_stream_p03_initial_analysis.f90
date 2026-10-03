!@descr: task 3 in the stream pipeline: the first 2D analysis of segmentation-picked particles, ending in the picking references
!==============================================================================
! MODULE: simple_commanders_stream_p03_initial_analysis
!
! PURPOSE:
!   The commander of stream stage 3, initial analysis (the master launches it
!   as prg=gen_pickrefs, an alias of initial_analysis in simple_stream). It
!   normalises the command line, then drives a stream_stage_initial_analysis
!   until the picking references are written, by the 3D route or from a GUI
!   selection, which stops the process early (jobs already submitted keep
!   running); everything the stage does lives in
!   simple_stream_stage_initial_analysis.
!
! ENTRY POINT:
!   commander_stream_p03_initial_analysis%execute(cline)
!
! SIGNALS:
!   SIGTERM only sets a flag (simple_stream_sigterm); the loop sees it before
!   the next iteration, and the stage's long steps stop between projects and
!   start no new job once it is set. The stage it replaces called exit(0)
!   inside the handler, which could truncate a project file being written.
!   The default action is restored before execute returns.
!==============================================================================
module simple_commanders_stream_p03_initial_analysis
use simple_defs,                          only: logfhandle
use simple_error,                         only: simple_exception
use simple_jiffys,                        only: simple_end
use simple_cmdline,                       only: cmdline
use simple_commander_base,                only: commander_base
use simple_stream_stage_initial_analysis, only: stream_stage_initial_analysis
use simple_stream_sigterm,                only: install_sigterm_handler, restore_sigterm_handler, sigterm_received
implicit none

public :: commander_stream_p03_initial_analysis
private
#include "simple_local_flags.inc"

integer, parameter :: NPARTS2D = 8 ! persistent workers for the stage's jobs

type, extends(commander_base) :: commander_stream_p03_initial_analysis
  contains
    procedure :: execute => exec_stream_p03_initial_analysis
end type commander_stream_p03_initial_analysis

contains

    subroutine exec_stream_p03_initial_analysis( self, cline )
        class(commander_stream_p03_initial_analysis), intent(inout) :: self
        class(cmdline),                     intent(inout) :: cline
        type(stream_stage_initial_analysis) :: stage
        call install_sigterm_handler()
        call set_initial_analysis_cline(cline)
        call stage%new(cline)
        do
            if( sigterm_received() )then
                write(logfhandle,'(A)') 'SIGTERM RECEIVED (INITIAL ANALYSIS)'
                exit
            endif
            if( stage%finished() ) exit
            call stage%iterate()
        enddo
        write(logfhandle,'(A)') '>>> TERMINATING PROCESS'
        call stage%finalize()
        call stage%kill()
        call restore_sigterm_handler()
        call simple_end('**** SIMPLE_STREAM_INITIAL_ANALYSIS NORMAL STOP ****')
    end subroutine exec_stream_p03_initial_analysis

    ! Everything the stage needs on its command line before params%new.
    subroutine set_initial_analysis_cline( cline )
        class(cmdline), intent(inout) :: cline
        if( .not. cline%defined('dir_target')     ) THROW_HARD('DIR_TARGET must be defined!')
        if( .not. cline%defined('mkdir')          ) call cline%set('mkdir',          'yes')
        if( .not. cline%defined('nptcls_per_cls') ) call cline%set('nptcls_per_cls', 100)
        if( .not. cline%defined('pick_roi')       ) call cline%set('pick_roi',       'yes')
        if( .not. cline%defined('outdir')         ) call cline%set('outdir',         '')
        ! automask2D
        if( .not. cline%defined('ngrow')          ) call cline%set('ngrow',          3)
        if( .not. cline%defined('winsz')          ) call cline%set('winsz',          5.)
        if( .not. cline%defined('amsklp')         ) call cline%set('amsklp',         20.)
        if( .not. cline%defined('edge')           ) call cline%set('edge',           6)
        ! workers for the stage's jobs, a quarter of the stage's threads each
        call cline%set('workers', NPARTS2D)
        if( cline%defined('nthr') ) call cline%set('worker_nthr', max(1, floor(real(cline%get_iarg('nthr'))/4.)))
    end subroutine set_initial_analysis_cline

end module simple_commanders_stream_p03_initial_analysis
