!@descr: developer wrapper for table-free streaming-SGD ab initio 2D classification
module simple_commanders_abinitio2d_sgd
use simple_commanders_api
use simple_commanders_abinitio2D, only: execute_abinitio2D_staged
implicit none
#include "simple_local_flags.inc"

public :: commander_abinitio2d_sgd, prepare_abinitio2d_sgd
private

! The public abinitio2D_sgd command deliberately has a small surface: it
! selects the development stream, validates its compatible objective, and
! delegates the staged workflow to the established abinitio2D commander.
! This keeps standard abinitio2D free of experimental user-facing controls.
type, extends(commander_base) :: commander_abinitio2d_sgd
    contains
    procedure :: execute => exec_abinitio2d_sgd
end type commander_abinitio2d_sgd

contains

    subroutine exec_abinitio2d_sgd( self, cline )
        class(commander_abinitio2d_sgd), intent(inout) :: self
        class(cmdline),                  intent(inout) :: cline
        integer :: start_stage, stop_stage, checkpoint_last_iter
        logical :: l_checkpoint
        call prepare_abinitio2d_sgd(cline)
        call get_checkpoint_request(start_stage, stop_stage, checkpoint_last_iter, l_checkpoint)
        call execute_abinitio2D_staged(cline, start_stage, stop_stage, checkpoint_last_iter, &
            &l_checkpoint, prepare_terminal_probability_command)
    end subroutine exec_abinitio2d_sgd

    subroutine prepare_abinitio2d_sgd( cline )
        class(cmdline), intent(inout) :: cline

        ! The streaming path uses a raw, noise-normalised Euclidean loss.
        ! Its analytical shift derivative is not defined for CC or for the
        ! separate continuous in-plane optimisation path.
        if( .not. cline%defined('objfun') ) call cline%set('objfun', 'euclid')
        if( cline%get_carg('objfun') .ne. 'euclid' )then
            THROW_HARD('abinitio2D_sgd requires objfun=euclid')
        endif
        if( .not. cline%defined('inpl_cont') ) call cline%set('inpl_cont', 'no')
        if( cline%get_carg('inpl_cont') .ne. 'no' )then
            THROW_HARD('abinitio2D_sgd does not support inpl_cont=yes')
        endif

        ! off is the conventional abinitio2D mode.  This development command
        ! requires an active stream and permits either fully streamed stage 4
        ! or the legacy/stream alternating schedule.
        if( .not. cline%defined('sgd_stage4_mode') ) call cline%set('sgd_stage4_mode', 'on')
        if( cline%get_carg('sgd_stage4_mode') .ne. 'on' .and. &
            &cline%get_carg('sgd_stage4_mode') .ne. 'alternate' )then
            THROW_HARD('abinitio2D_sgd requires sgd_stage4_mode=on or alternate')
        endif

        ! These fields are internal compatibility inputs consumed by the
        ! controller and strategy layers; table is intentionally unsupported.
        call cline%set('sgd',      'yes')
        call cline%set('sgd_path', 'stream')

        ! Delegate to the standard orchestration label so parameters, project
        ! handling, stage directories, and terminal cleanup retain their
        ! established abinitio2D lifecycle.
        call cline%set('prg', 'abinitio2D')
    end subroutine prepare_abinitio2d_sgd

    subroutine get_checkpoint_request( start_stage, stop_stage, checkpoint_last_iter, l_checkpoint )
        integer, intent(out) :: start_stage, stop_stage, checkpoint_last_iter
        logical, intent(out) :: l_checkpoint
        character(len=32) :: env_value
        integer :: env_status

        start_stage = 1
        stop_stage = 0
        checkpoint_last_iter = 0
        l_checkpoint = .false.
        call get_environment_variable('JOINT2D_SGD_CHECKPOINT_START_STAGE', env_value, status=env_status)
        if( env_status == 0 .and. len_trim(env_value) > 0 )then
            read(env_value,*) start_stage
            l_checkpoint = .true.
        endif
        call get_environment_variable('JOINT2D_SGD_CHECKPOINT_STOP_STAGE', env_value, status=env_status)
        if( env_status == 0 .and. len_trim(env_value) > 0 )then
            read(env_value,*) stop_stage
            l_checkpoint = .true.
        endif
        call get_environment_variable('JOINT2D_SGD_CHECKPOINT_LAST_ITER', env_value, status=env_status)
        if( env_status == 0 .and. len_trim(env_value) > 0 )then
            read(env_value,*) checkpoint_last_iter
            l_checkpoint = .true.
        endif
    end subroutine get_checkpoint_request

    subroutine prepare_terminal_probability_command( cline )
        class(cmdline), intent(inout) :: cline
        ! The terminal pass is conventional probability refinement.  Remove
        ! development-only stage settings so parameter defaults select it.
        call cline%delete('sgd')
        call cline%delete('sgd_path')
        call cline%delete('sgd_stage4_mode')
    end subroutine prepare_terminal_probability_command

end module simple_commanders_abinitio2d_sgd
