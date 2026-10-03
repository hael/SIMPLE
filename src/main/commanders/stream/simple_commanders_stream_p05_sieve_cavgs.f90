!@descr: task 5 in the stream pipeline: continuous particle sieving with staged chunk generation and class-average rejection
!==============================================================================
! MODULE: simple_commanders_stream_p05_sieve_cavgs
!
! PURPOSE:
!   The commander of stream stage 5. It normalises the command line, then
!   drives a stream_stage_sieve until the stream stops; everything the stage
!   does lives in simple_stream_stage_sieve.
!
! ENTRY POINT:
!   commander_stream_p05_sieve_cavgs%execute(cline)
!
! SIGNALS:
!   SIGTERM only sets a flag (simple_stream_sigterm); the loop sees it before
!   the next iteration and the stage finalises normally. The default action is
!   restored before execute returns.
!==============================================================================
module simple_commanders_stream_p05_sieve_cavgs
use simple_defs,               only: logfhandle
use simple_jiffys,             only: simple_end
use simple_cmdline,            only: cmdline
use simple_commander_base,     only: commander_base
use simple_ptcl_sieve,         only: DEFAULT_COARSE_POP_THRESHOLD, DEFAULT_FINE_POP_THRESHOLD, DEFAULT_COARSE_BOX,&
                                    &DEFAULT_FINE_BOX, DEFAULT_COARSE_NSAMPLE, DEFAULT_FINE_NSAMPLE, DEFAULT_LPSTART,&
                                    &DEFAULT_COARSE_LP, DEFAULT_FINE_LP, DEFAULT_NCLS
use simple_stream_stage_sieve, only: stream_stage_sieve
use simple_stream_sigterm,     only: install_sigterm_handler, restore_sigterm_handler, sigterm_received
implicit none

public :: commander_stream_p05_sieve_cavgs
private

type, extends(commander_base) :: commander_stream_p05_sieve_cavgs
  contains
    procedure :: execute => exec_stream_p05_sieve_cavgs
end type commander_stream_p05_sieve_cavgs

contains

    subroutine exec_stream_p05_sieve_cavgs( self, cline )
        class(commander_stream_p05_sieve_cavgs), intent(inout) :: self
        class(cmdline),                intent(inout) :: cline
        type(stream_stage_sieve) :: stage
        call install_sigterm_handler()
        call set_sieve_cline(cline)
        call stage%new(cline)
        do
            if( sigterm_received() )then
                write(logfhandle,'(A)') 'SIGTERM RECEIVED (PARTICLE SIEVING)'
                exit
            endif
            if( stage%finished() ) exit
            call stage%iterate()
        enddo
        write(logfhandle,'(A)') '>>> TERMINATING PROCESS'
        call stage%finalize()
        call stage%kill()
        call restore_sigterm_handler()
        call simple_end('**** SIMPLE_STREAM_SIEVE_CAVGS NORMAL STOP ****')
    end subroutine exec_stream_p05_sieve_cavgs

    ! Everything the stage needs on its command line before params%new: the sieve's tuning
    ! defaults, and one persistent worker per chunk with the stage's thread count.
    subroutine set_sieve_cline( cline )
        class(cmdline), intent(inout) :: cline
        call cline%set('mkdir', 'yes')
        if( .not. cline%defined('walltime')       ) call cline%set('walltime',       29 * 60)
        if( .not. cline%defined('outdir')         ) call cline%set('outdir',         '')
        if( .not. cline%defined('nmics')          ) call cline%set('nmics',          100)
        if( .not. cline%defined('nptcls_coarse')  ) call cline%set('nptcls_coarse',  DEFAULT_COARSE_POP_THRESHOLD)
        if( .not. cline%defined('nptcls_fine')    ) call cline%set('nptcls_fine',    DEFAULT_FINE_POP_THRESHOLD)
        if( .not. cline%defined('maxnchunks')     ) call cline%set('maxnchunks',     0)
        if( .not. cline%defined('lpstart')        ) call cline%set('lpstart',        DEFAULT_LPSTART)
        if( .not. cline%defined('lpstop_coarse')  ) call cline%set('lpstop_coarse',  DEFAULT_COARSE_LP)
        if( .not. cline%defined('lpstop_fine')    ) call cline%set('lpstop_fine',    DEFAULT_FINE_LP)
        if( .not. cline%defined('box_coarse')     ) call cline%set('box_coarse',     DEFAULT_COARSE_BOX)
        if( .not. cline%defined('box_fine')       ) call cline%set('box_fine',       DEFAULT_FINE_BOX)
        if( .not. cline%defined('nsample_coarse') ) call cline%set('nsample_coarse', DEFAULT_COARSE_NSAMPLE)
        if( .not. cline%defined('nsample_fine')   ) call cline%set('nsample_fine',   DEFAULT_FINE_NSAMPLE)
        if( .not. cline%defined('ncls_coarse')    ) call cline%set('ncls_coarse',    DEFAULT_NCLS)
        if( .not. cline%defined('ncls_fine')      ) call cline%set('ncls_fine',      DEFAULT_NCLS)
        if( .not. cline%defined('use_model')      ) call cline%set('use_model',      'yes')
        if( .not. cline%defined('single_pass')    ) call cline%set('single_pass',    'no')
        ! one persistent worker per chunk, each with the stage's thread count
        if( cline%defined('nchunks') ) call cline%set('workers',     cline%get_iarg('nchunks'))
        if( cline%defined('nthr')    ) call cline%set('worker_nthr', cline%get_iarg('nthr'))
    end subroutine set_sieve_cline

end module simple_commanders_stream_p05_sieve_cavgs
