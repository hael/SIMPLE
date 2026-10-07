!@descr: task 7 in the stream pipeline: multistate 3D reconstruction of the particles pool 2D exports
!==============================================================================
! MODULE: simple_commanders_stream_p07_solve3D_multistate
!
! PURPOSE:
!   The commander of stream stage 7 (the master launches it as
!   prg=solve3D_stream). It normalises the command line, then drives a
!   stream_stage_solve3D until the stream stops; everything the stage does
!   lives in simple_stream_stage_solve3D.
!
! 3D JOBS:
!   The settings of the solve3D and solve3D_addon runs are the named
!   constants below; a command line overrides them with nstates, nstages,
!   nptcls3D_max (the first set's particles at most, whole stacks in order),
!   and nparts3D and nthr3D (the jobs' parts and threads; the stage's own
!   nparts and nthr are the master's settings for the stage). The 3D jobs set
!   their own low-pass limits.
!
! ENTRY POINT:
!   commander_stream_p07_solve3D_multistate%execute(cline)
!
! SIGNALS:
!   SIGTERM only sets a flag (simple_stream_sigterm); the loop sees it before
!   the next iteration and the stage finalises normally. The default action is
!   restored before execute returns.
!==============================================================================
module simple_commanders_stream_p07_solve3D_multistate
use simple_defs,                    only: logfhandle
use simple_error,                   only: simple_exception
use simple_defs_fname,              only: METADATA_EXT
use simple_defs_stream,             only: NPTCLS_FIRST3D
use simple_jiffys,                  only: simple_end
use simple_cmdline,                 only: cmdline
use simple_commander_base,          only: commander_base
use simple_stream_stage_solve3D, only: stream_stage_solve3D, MIN_PTCLS_PER_STATE
use simple_stream_sigterm,          only: install_sigterm_handler, restore_sigterm_handler, sigterm_received
use simple_gui_metadata_stream_solve3D_multistate, only: MAX_STATES_SOLVE3D_MULTISTATE
implicit none
#include "simple_local_flags.inc"

public :: commander_stream_p07_solve3D_multistate
public :: set_solve3D_cline ! the stage's command-line defaults, for the chained stream tests
private

integer, parameter :: NSTATES3D    = 3      ! states of the 3D
integer, parameter :: NSTAGES3D    = 5      ! solve3D stages
integer, parameter :: NPARTS3D     = 8      ! parts and threads of each 3D job
integer, parameter :: NTHR3D       = 8
integer, parameter :: NPTCLS3D_MAX = NPTCLS_FIRST3D ! the first solve3D's particles at most, as many as p06's first publication waits for

type, extends(commander_base) :: commander_stream_p07_solve3D_multistate
  contains
    procedure :: execute => exec_stream_p07_solve3D_multistate
end type commander_stream_p07_solve3D_multistate

contains

    subroutine exec_stream_p07_solve3D_multistate( self, cline )
        class(commander_stream_p07_solve3D_multistate), intent(inout) :: self
        class(cmdline),                          intent(inout) :: cline
        type(stream_stage_solve3D) :: stage
        call install_sigterm_handler()
        call set_solve3D_cline(cline)
        call stage%new(cline)
        do
            if( sigterm_received() )then
                write(logfhandle,'(A)') 'SIGTERM RECEIVED (SOLVE3D)'
                exit
            endif
            if( stage%finished() ) exit
            call stage%iterate()
        enddo
        write(logfhandle,'(A)') '>>> TERMINATING PROCESS'
        call stage%finalize()
        call stage%kill()
        call restore_sigterm_handler()
        call simple_end('**** SIMPLE_STREAM_SOLVE3D_MULTISTATE NORMAL STOP ****')
    end subroutine exec_stream_p07_solve3D_multistate

    ! Everything the stage needs on its command line before params%new: the 3D job settings
    ! unless given, and the project name when none is. The states are checked: a multistate run
    ! needs two at least, and the GUI's status holds MAX_STATES_SOLVE3D_MULTISTATE (20). So is
    ! the first solve3D's cap: it must allow the run's minimum, MIN_PTCLS_PER_STATE per state.
    subroutine set_solve3D_cline( cline )
        class(cmdline), intent(inout) :: cline
        call cline%set('oritype', 'mic')
        call cline%set('mkdir',   'yes')
        if( .not. cline%defined('nstates')  ) call cline%set('nstates',  NSTATES3D)
        if( cline%get_iarg('nstates') < 2 ) THROW_HARD('multistate 3D needs nstates of 2 at least')
        if( cline%get_iarg('nstates') > MAX_STATES_SOLVE3D_MULTISTATE ) THROW_HARD('multistate 3D takes nstates of 20 at most')
        if( .not. cline%defined('nstages')  ) call cline%set('nstages',  NSTAGES3D)
        if( .not. cline%defined('nparts3D') ) call cline%set('nparts3D', NPARTS3D)
        if( .not. cline%defined('nthr3D')   ) call cline%set('nthr3D',   NTHR3D)
        if( .not. cline%defined('nptcls3D_max') ) call cline%set('nptcls3D_max', NPTCLS3D_MAX)
        if( cline%get_iarg('nptcls3D_max') < MIN_PTCLS_PER_STATE * cline%get_iarg('nstates') ) THROW_HARD('nptcls3D_max is below the first solve3D''s minimum of 5 particles per state')
        if( .not. cline%defined('projfile') )then
            call cline%set('projname', 'stream_3Dmultistate')
            call cline%set('projfile', 'stream_3Dmultistate'//METADATA_EXT)
        endif
    end subroutine set_solve3D_cline

end module simple_commanders_stream_p07_solve3D_multistate
