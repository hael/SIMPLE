!@descr: task 3 in the stream pipeline: the first 2D analysis of segmentation-picked particles, ending in the picking references
!==============================================================================
! MODULE: simple_commanders_stream_p03_initial_analysis
!
! PURPOSE:
!   The commander of stream stage 3, initial analysis (the master launches it
!   as prg=gen_pickrefs, an alias of initial_analysis in simple_stream). It
!   normalises the command line, then drives a stream_stage_initial_analysis
!   until the picking references are written, by the 3D route or from a GUI
!   selection, which stops the process early (the jobs already submitted are
!   cancelled); everything the stage does lives in
!   simple_stream_stage_initial_analysis.
!
!   The 3D route's settings (nstates_pickrefs, nstages_pickrefs,
!   lpstop_pickrefs, nspace_pickrefs, nthr3D_pickrefs, and nrestarts_collapse,
!   lpstart_ini3D, lpstop_ini3D) and the jobs' resources (nthr2D, nparts,
!   nchunks) come from the master, or the defaults below when run on its own.
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
public :: set_initial_analysis_cline ! for simple_stream_stage_initial_analysis_tester
private
#include "simple_local_flags.inc"

integer, parameter :: NPARTS2D = 8 ! persistent workers for the stage's jobs
! the defaults when the master gives none (the stage run on its own; simple_stream_master_resources)
integer, parameter :: NTHR2D_JOBS  = 16 ! threads of the 2D jobs and of the sieve's chunks
integer, parameter :: NCHUNKS      = 4  ! the sieve's chunks at once
! the 3D route to the picking references (stream fix plan, decision 20)
integer, parameter :: NSTATES_PICKREFS   = 3
integer, parameter :: NSTAGES_PICKREFS   = 3
real,    parameter :: LPSTOP_PICKREFS    = 8.   ! A
integer, parameter :: NSPACE_PICKREFS    = 50   ! reprojections
integer, parameter :: NTHR3D_PICKREFS    = 16   ! threads of the 3D job and the reprojection
integer, parameter :: NRESTARTS_COLLAPSE = 3
real,    parameter :: LPSTART_INI3D      = 100. ! A
real,    parameter :: LPSTOP_INI3D       = 20.  ! A

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

    ! Everything the stage needs on its command line before params%new: defaults, the 3D route's
    ! settings and its jobs' resources unless given, and their range checks.
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
        ! the jobs' resources
        if( .not. cline%defined('nthr2D')          ) call cline%set('nthr2D',          NTHR2D_JOBS)
        if( .not. cline%defined('nchunks')         ) call cline%set('nchunks',         NCHUNKS)
        ! the 3D route to the picking references; 0 means the default
        call default_unless_positive('nstates_pickrefs',   real(NSTATES_PICKREFS))
        call default_unless_positive('nstages_pickrefs',   real(NSTAGES_PICKREFS))
        call default_unless_positive('lpstop_pickrefs',    LPSTOP_PICKREFS)
        call default_unless_positive('nspace_pickrefs',    real(NSPACE_PICKREFS))
        call default_unless_positive('nthr3D_pickrefs',    real(NTHR3D_PICKREFS))
        call default_unless_positive('nrestarts_collapse', real(NRESTARTS_COLLAPSE))
        call default_unless_positive('lpstart_ini3D',      LPSTART_INI3D)
        call default_unless_positive('lpstop_ini3D',       LPSTOP_INI3D)
        call check_pickrefs_settings(cline)
        ! the threads each of the stage's jobs claims on a persistent worker (the master's server,
        ! or the stage's own when it runs alone): as many as the largest of them uses, the 2D jobs
        ! and the sieve's chunks (nthr2D) or the 3D job (nthr3D_pickrefs); the extractions use
        ! fewer. The NPARTS2D workers start only when the stage runs alone
        call cline%set('workers', NPARTS2D)
        call cline%set('worker_nthr', max(cline%get_iarg('nthr2D'), cline%get_iarg('nthr3D_pickrefs')))

    contains

        subroutine default_unless_positive( key, val )
            character(len=*), intent(in) :: key
            real,             intent(in) :: val
            if( cline%defined(key) )then
                if( cline%get_rarg(key) > 0. ) return
            endif
            call cline%set(key, val)
        end subroutine default_unless_positive

    end subroutine set_initial_analysis_cline

    ! The 3D route's settings, once defaulted: two states at least (solve3D_cavgs' conditional
    ! restarts need them, and the state choice a choice); lpstart_ini3D > lpstop_ini3D >=
    ! lpstop_pickrefs > 0; positive counts.
    subroutine check_pickrefs_settings( cline )
        class(cmdline), intent(in) :: cline
        if( cline%get_iarg('nstates_pickrefs') < 2 ) THROW_HARD('nstates_pickrefs must be at least 2')
        if( cline%get_iarg('nstages_pickrefs') < 1 ) THROW_HARD('nstages_pickrefs must be at least 1')
        if( cline%get_iarg('nspace_pickrefs')  < 1 ) THROW_HARD('nspace_pickrefs must be at least 1')
        if( cline%get_iarg('nthr3D_pickrefs')  < 1 ) THROW_HARD('nthr3D_pickrefs must be at least 1')
        if( cline%get_iarg('nrestarts_collapse') < 1 ) THROW_HARD('nrestarts_collapse must be at least 1')
        if( cline%get_rarg('lpstop_ini3D') >= cline%get_rarg('lpstart_ini3D') ) THROW_HARD('lpstart_ini3D must exceed lpstop_ini3D')
        if( cline%get_rarg('lpstop_pickrefs') > cline%get_rarg('lpstop_ini3D') ) THROW_HARD('lpstop_pickrefs must not exceed lpstop_ini3D')
        if( cline%get_iarg('nthr2D')  < 1 ) THROW_HARD('nthr2D must be at least 1')
        if( cline%get_iarg('nchunks') < 1 ) THROW_HARD('nchunks must be at least 1')
    end subroutine check_pickrefs_settings

end module simple_commanders_stream_p03_initial_analysis
