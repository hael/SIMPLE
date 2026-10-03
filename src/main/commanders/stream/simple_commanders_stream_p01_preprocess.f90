!@descr: task 1 in the stream pipeline: pre-processing (movie registration, CTF estimation, segmentation-based picking)
!==============================================================================
! MODULE: simple_commanders_stream_p01_preprocess
!
! PURPOSE:
!   The commander of stream stage 1. It normalises the command line, then
!   drives a stream_stage_preprocess until the stream stops; everything the
!   stage does lives in simple_stream_stage_preprocess.
!
! ENTRY POINT:
!   commander_stream_p01_preprocess%execute(cline)
!
! SIGNALS:
!   SIGTERM only sets a flag (simple_stream_sigterm); the loop sees it before
!   the next iteration and the stage finalises normally. The wait for the
!   first movie in new() also stops on it. The default action is restored
!   before execute returns.
!==============================================================================
module simple_commanders_stream_p01_preprocess
use simple_defs,                    only: logfhandle, STDLEN, HP_CTF_ESTIMATE, LP_CTF_ESTIMATE, DFMIN_DEFAULT, DFMAX_DEFAULT
use simple_defs_stream,             only: STREAM_CTFRES_THRESHOLD
use simple_defs_environment,        only: SIMPLE_STREAM_PREPROC_NTHR, SIMPLE_STREAM_PREPROC_NPARTS
use simple_error,                   only: simple_exception
use simple_string,                  only: string
use simple_string_utils,            only: str2int
use simple_jiffys,                  only: simple_end
use simple_cmdline,                 only: cmdline
use simple_commander_base,          only: commander_base
use simple_stream_stage_preprocess, only: stream_stage_preprocess
use simple_stream_sigterm,          only: install_sigterm_handler, restore_sigterm_handler, sigterm_received
implicit none

public :: commander_stream_p01_preprocess
private
#include "simple_local_flags.inc"

type, extends(commander_base) :: commander_stream_p01_preprocess
  contains
    procedure :: execute => exec_stream_p01_preprocess
end type commander_stream_p01_preprocess

contains

    subroutine exec_stream_p01_preprocess( self, cline )
        class(commander_stream_p01_preprocess), intent(inout) :: self
        class(cmdline),               intent(inout) :: cline
        type(stream_stage_preprocess) :: stage
        call install_sigterm_handler()
        call set_preprocess_stream_cline(cline)
        call stage%new(cline)
        do
            if( sigterm_received() )then
                write(logfhandle,'(A)') 'SIGTERM RECEIVED (PREPROCESSING)'
                exit
            endif
            if( stage%finished() ) exit
            call stage%iterate()
        enddo
        write(logfhandle,'(A)') '>>> TERMINATING PROCESS'
        call stage%finalize()
        call stage%kill()
        call restore_sigterm_handler()
        call simple_end('**** SIMPLE_STREAM_PREPROC NORMAL STOP ****')
    end subroutine exec_stream_p01_preprocess

    ! Everything the stage needs on its command line before params%new: fixed settings,
    ! defaults, environment overrides, and the gain options that need no movies.
    subroutine set_preprocess_stream_cline( cline )
        class(cmdline), intent(inout) :: cline
        character(len=STDLEN) :: env_val
        integer               :: envlen
        ! fixed for this stage
        call cline%set('oritype',     'mic')
        call cline%set('mkdir',       'yes')
        call cline%set('reject_mics', 'no')
        call cline%set('numlen',      5)
        call cline%set('stream',      'yes')
        if( .not. cline%defined('walltime')        ) call cline%set('walltime',        29.0*60.0) ! 29 minutes
        if( .not. cline%defined('nmics')           ) call cline%set('nmics',           0)
        ! motion correction
        if( .not. cline%defined('trs')             ) call cline%set('trs',             20.)
        if( .not. cline%defined('lpstart')         ) call cline%set('lpstart',         8.)
        if( .not. cline%defined('lpstop')          ) call cline%set('lpstop',          5.)
        if( .not. cline%defined('bfac')            ) call cline%set('bfac',            50.)
        if( .not. cline%defined('mcconvention')    ) call cline%set('mcconvention',    'simple')
        if( .not. cline%defined('eer_upsampling')  ) call cline%set('eer_upsampling',  1)
        if( .not. cline%defined('algorithm')       ) call cline%set('algorithm',       'patch')
        if( .not. cline%defined('mcpatch')         ) call cline%set('mcpatch',         'yes')
        if( .not. cline%defined('mcpatch_thres')   ) call cline%set('mcpatch_thres',   'yes')
        if( .not. cline%defined('tilt_thres')      ) call cline%set('tilt_thres',      0.05)
        if( .not. cline%defined('beamtilt')        ) call cline%set('beamtilt',        'no')
        ! CTF estimation
        if( .not. cline%defined('pspecsz')         ) call cline%set('pspecsz',         512)
        if( .not. cline%defined('hp_ctf_estimate') ) call cline%set('hp_ctf_estimate', HP_CTF_ESTIMATE)
        if( .not. cline%defined('lp_ctf_estimate') ) call cline%set('lp_ctf_estimate', LP_CTF_ESTIMATE)
        if( .not. cline%defined('dfmin')           ) call cline%set('dfmin',           DFMIN_DEFAULT)
        if( .not. cline%defined('dfmax')           ) call cline%set('dfmax',           DFMAX_DEFAULT)
        if( .not. cline%defined('ctfpatch')        ) call cline%set('ctfpatch',        'yes')
        if( .not. cline%defined('ctfresthreshold') ) call cline%set('ctfresthreshold', STREAM_CTFRES_THRESHOLD)
        ! environment overrides
        call get_environment_variable(SIMPLE_STREAM_PREPROC_NTHR, env_val, envlen)
        if( envlen > 0 ) call cline%set('nthr', str2int(env_val))
        call get_environment_variable(SIMPLE_STREAM_PREPROC_NPARTS, env_val, envlen)
        if( envlen > 0 ) call cline%set('nparts', str2int(env_val))
        call set_flipgain(cline)
    end subroutine set_preprocess_stream_cline

    ! Maps the GUI gain options onto flip_gain's modes. flip_auto and generate need
    ! movies and are resolved by the stage.
    subroutine set_flipgain( cline )
        class(cmdline), intent(inout) :: cline
        type(string)                  :: flipgain
        character(len=:), allocatable :: error_message
        if( .not. cline%defined('flipgain') ) return
        flipgain = cline%get_carg('flipgain')
        select case( flipgain%to_char() )
            case( 'none', 'no' )
                call cline%set('flipgain', 'no')
            case( 'flip_x', 'x' )
                call cline%set('flipgain', 'x')
            case( 'flip_y', 'y' )
                call cline%set('flipgain', 'y')
            case( 'flip_xy', 'xy' )
                call cline%set('flipgain', 'xy')
            case( 'yx' )
                ! flip_gain mode, kept as given
            case( 'flip_auto', 'generate' )
                ! resolved by the stage once movies arrive
            case DEFAULT
                error_message = 'Unknown gain processing option: '//flipgain%to_char()
                THROW_HARD(error_message)
        end select
    end subroutine set_flipgain

end module simple_commanders_stream_p01_preprocess
