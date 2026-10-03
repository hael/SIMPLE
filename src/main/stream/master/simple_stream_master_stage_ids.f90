!@descr: the stream stages the master runs: their ids, job names, GUI keys, and their pipes in simple_stream_state
!==============================================================================
! MODULE: simple_stream_master_stage_ids
!
! PURPOSE:
!   Gives each stage of the stream an id (1..NSTAGES, in the order the GUI
!   heartbeat lists them) and maps the id to its job name (folder and log),
!   the stem of its GUI keys (terminate_<key>, restart_<key>), and its two
!   pipes, which stay the named arrays of simple_stream_state that the stages
!   read. A stage writes *_in(2) and reads *_out(1); the master reads *_in(1)
!   and writes *_out(2), and keeps every end open across forks.
!
!   Stateless: the pipes are simple_stream_state's.
!==============================================================================
module simple_stream_master_stage_ids
use unix,                only: c_pipe, c_close, c_fcntl, F_GETFL, F_SETFL, O_NONBLOCK
use simple_error,        only: simple_exception
use simple_defs_stream,  only: PREPROC_JOB_NAME, OPTICS_JOB_NAME, OPENING2D_JOB_NAME, REFPICK_JOB_NAME,&
                               &SIEVING_JOB_NAME, CLASS2D_JOB_NAME, MULTISTATE3D_JOB_NAME
use simple_stream_state, only: ipc_pipe_preprocess_in, ipc_pipe_preprocess_out, ipc_pipe_assign_optics_in,&
                               &ipc_pipe_assign_optics_out, ipc_pipe_initial_analysis_in, ipc_pipe_initial_analysis_out,&
                               &ipc_pipe_refpick_in, ipc_pipe_refpick_out, ipc_pipe_sieve_cavgs_in, ipc_pipe_sieve_cavgs_out,&
                               &ipc_pipe_pool2D_in, ipc_pipe_pool2D_out, ipc_pipe_solve3D_multistate_in,&
                               &ipc_pipe_solve3D_multistate_out
implicit none

public :: NSTAGES, STAGE_PREPROCESS, STAGE_ASSIGN_OPTICS, STAGE_INITIAL_ANALYSIS, STAGE_REFERENCE_PICKING
public :: STAGE_PARTICLE_SIEVING, STAGE_POOL2D, STAGE_SOLVE3D
public :: stage_job_name, stage_gui_key, stage_label
public :: open_stage_pipes, close_stage_pipes, close_other_pipe_ends, master_fds, stage_fds
private
#include "simple_local_flags.inc"

integer, parameter :: STAGE_PREPROCESS        = 1
integer, parameter :: STAGE_ASSIGN_OPTICS     = 2
integer, parameter :: STAGE_INITIAL_ANALYSIS  = 3
integer, parameter :: STAGE_REFERENCE_PICKING = 4
integer, parameter :: STAGE_PARTICLE_SIEVING  = 5
integer, parameter :: STAGE_POOL2D            = 6
integer, parameter :: STAGE_SOLVE3D        = 7
integer, parameter :: NSTAGES                 = 7

! what pipe_op does to a pipe
integer, parameter :: OP_OPEN        = 1
integer, parameter :: OP_CLOSE       = 2
integer, parameter :: OP_CLOSE_OTHER = 3 ! close the ends that are not kept

contains

    !> The job name: the stage's folder, and its log is <name>.log.
    function stage_job_name( id ) result( name )
        integer, intent(in) :: id
        character(len=:), allocatable :: name
        select case(id)
            case(STAGE_PREPROCESS);        name = PREPROC_JOB_NAME
            case(STAGE_ASSIGN_OPTICS);     name = OPTICS_JOB_NAME
            case(STAGE_INITIAL_ANALYSIS);  name = OPENING2D_JOB_NAME
            case(STAGE_REFERENCE_PICKING); name = REFPICK_JOB_NAME
            case(STAGE_PARTICLE_SIEVING);  name = SIEVING_JOB_NAME
            case(STAGE_POOL2D);            name = CLASS2D_JOB_NAME
            case(STAGE_SOLVE3D);        name = MULTISTATE3D_JOB_NAME
            case default;                  THROW_HARD('unknown stream stage id')
        end select
    end function stage_job_name

    !> The stem of the stage's GUI keys: terminate_<key>, restart_<key>.
    function stage_gui_key( id ) result( key )
        integer, intent(in) :: id
        character(len=:), allocatable :: key
        select case(id)
            case(STAGE_PREPROCESS);        key = 'preprocess'
            case(STAGE_ASSIGN_OPTICS);     key = 'optics_assignment'
            case(STAGE_INITIAL_ANALYSIS);  key = 'opening2D'
            case(STAGE_REFERENCE_PICKING); key = 'reference_picking'
            case(STAGE_PARTICLE_SIEVING);  key = 'particle_sieving'
            case(STAGE_POOL2D);            key = 'pool2D'
            case(STAGE_SOLVE3D);        key = 'solve3D_multistate'
            case default;                  THROW_HARD('unknown stream stage id')
        end select
    end function stage_gui_key

    !> The stage's name in the master's log.
    function stage_label( id ) result( label )
        integer, intent(in) :: id
        character(len=:), allocatable :: label
        select case(id)
            case(STAGE_PREPROCESS);        label = 'PREPROCESS'
            case(STAGE_ASSIGN_OPTICS);     label = 'ASSIGN OPTICS'
            case(STAGE_INITIAL_ANALYSIS);  label = 'INITIAL ANALYSIS'
            case(STAGE_REFERENCE_PICKING); label = 'REFERENCE PICKING'
            case(STAGE_PARTICLE_SIEVING);  label = 'PARTICLE SIEVING'
            case(STAGE_POOL2D);            label = 'POOL2D'
            case(STAGE_SOLVE3D);        label = 'SOLVE3D MULTISTATE'
            case default;                  THROW_HARD('unknown stream stage id')
        end select
    end function stage_label

    !> Creates the stage's two pipes with every end non-blocking.
    subroutine open_stage_pipes( id )
        integer, intent(in) :: id
        call pipe_op(id, OP_OPEN, [-1, -1])
    end subroutine open_stage_pipes

    !> Closes the stage's pipe ends this process holds.
    subroutine close_stage_pipes( id )
        integer, intent(in) :: id
        call pipe_op(id, OP_CLOSE, [-1, -1])
    end subroutine close_stage_pipes

    !> In the forked process of stage @p id: closes every pipe end of every stage but the two
    !! this stage uses, so a stage that exits closes its pipes for good.
    subroutine close_other_pipe_ends( id )
        integer, intent(in) :: id
        integer :: keep(2), jd
        call stage_fds(id, keep(1), keep(2))
        do jd = 1,NSTAGES
            call pipe_op(jd, OP_CLOSE_OTHER, keep)
        enddo
    end subroutine close_other_pipe_ends

    !> The ends the master uses: it reads what the stage sends and writes the GUI updates.
    subroutine master_fds( id, fd_read, fd_write )
        integer, intent(in)  :: id
        integer, intent(out) :: fd_read, fd_write
        integer :: fds_in(2), fds_out(2)
        call stage_pipes(id, fds_in, fds_out)
        fd_read  = fds_in(1)
        fd_write = fds_out(2)
    end subroutine master_fds

    !> The ends the stage uses: it writes its messages and reads the GUI updates.
    subroutine stage_fds( id, fd_write, fd_read )
        integer, intent(in)  :: id
        integer, intent(out) :: fd_write, fd_read
        integer :: fds_in(2), fds_out(2)
        call stage_pipes(id, fds_in, fds_out)
        fd_write = fds_in(2)
        fd_read  = fds_out(1)
    end subroutine stage_fds

    ! A copy of the stage's two pipes.
    subroutine stage_pipes( id, fds_in, fds_out )
        integer, intent(in)  :: id
        integer, intent(out) :: fds_in(2), fds_out(2)
        select case(id)
            case(STAGE_PREPROCESS)
                fds_in = ipc_pipe_preprocess_in;           fds_out = ipc_pipe_preprocess_out
            case(STAGE_ASSIGN_OPTICS)
                fds_in = ipc_pipe_assign_optics_in;        fds_out = ipc_pipe_assign_optics_out
            case(STAGE_INITIAL_ANALYSIS)
                fds_in = ipc_pipe_initial_analysis_in;     fds_out = ipc_pipe_initial_analysis_out
            case(STAGE_REFERENCE_PICKING)
                fds_in = ipc_pipe_refpick_in;              fds_out = ipc_pipe_refpick_out
            case(STAGE_PARTICLE_SIEVING)
                fds_in = ipc_pipe_sieve_cavgs_in;          fds_out = ipc_pipe_sieve_cavgs_out
            case(STAGE_POOL2D)
                fds_in = ipc_pipe_pool2D_in;               fds_out = ipc_pipe_pool2D_out
            case(STAGE_SOLVE3D)
                fds_in = ipc_pipe_solve3D_multistate_in; fds_out = ipc_pipe_solve3D_multistate_out
            case default
                THROW_HARD('unknown stream stage id')
        end select
    end subroutine stage_pipes

    ! Applies op to the stage's two named pipe arrays in place.
    subroutine pipe_op( id, op, keep )
        integer, intent(in) :: id, op, keep(2)
        select case(id)
            case(STAGE_PREPROCESS)
                call apply(ipc_pipe_preprocess_in,           op, keep)
                call apply(ipc_pipe_preprocess_out,          op, keep)
            case(STAGE_ASSIGN_OPTICS)
                call apply(ipc_pipe_assign_optics_in,        op, keep)
                call apply(ipc_pipe_assign_optics_out,       op, keep)
            case(STAGE_INITIAL_ANALYSIS)
                call apply(ipc_pipe_initial_analysis_in,     op, keep)
                call apply(ipc_pipe_initial_analysis_out,    op, keep)
            case(STAGE_REFERENCE_PICKING)
                call apply(ipc_pipe_refpick_in,              op, keep)
                call apply(ipc_pipe_refpick_out,             op, keep)
            case(STAGE_PARTICLE_SIEVING)
                call apply(ipc_pipe_sieve_cavgs_in,          op, keep)
                call apply(ipc_pipe_sieve_cavgs_out,         op, keep)
            case(STAGE_POOL2D)
                call apply(ipc_pipe_pool2D_in,               op, keep)
                call apply(ipc_pipe_pool2D_out,              op, keep)
            case(STAGE_SOLVE3D)
                call apply(ipc_pipe_solve3D_multistate_in,  op, keep)
                call apply(ipc_pipe_solve3D_multistate_out, op, keep)
            case default
                THROW_HARD('unknown stream stage id')
        end select
    end subroutine pipe_op

    subroutine apply( pipe, op, keep )
        integer, intent(inout) :: pipe(2)
        integer, intent(in)    :: op, keep(2)
        integer :: rc, flags, iend
        select case(op)
            case(OP_OPEN)
                rc = c_pipe(pipe)
                if( rc /= 0 ) THROW_HARD('failed to create IPC pipe')
                do iend = 1,2
                    flags = c_fcntl(pipe(iend), F_GETFL, 0)
                    if( flags < 0 ) THROW_HARD('failed to get IPC pipe flags')
                    rc = c_fcntl(pipe(iend), F_SETFL, ior(flags, O_NONBLOCK))
                    if( rc < 0 ) THROW_HARD('failed to make an IPC pipe end non-blocking')
                enddo
            case(OP_CLOSE)
                do iend = 1,2
                    if( pipe(iend) < 0 ) cycle
                    rc = c_close(pipe(iend))
                    if( rc /= 0 ) THROW_WARN('failed to close an IPC pipe end')
                    pipe(iend) = -1
                enddo
            case(OP_CLOSE_OTHER)
                do iend = 1,2
                    if( pipe(iend) < 0 ) cycle
                    if( any(keep == pipe(iend)) ) cycle
                    rc = c_close(pipe(iend))
                    pipe(iend) = -1
                enddo
        end select
    end subroutine apply

end module simple_stream_master_stage_ids
