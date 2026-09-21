!@descr: thread-local accumulation and iteration reporting for pose_cont refinement
module simple_pose_cont_run_stats
use simple_defs, only: dp, logfhandle
use simple_error, only: simple_exception
use simple_pose_cont_refine3D_adapter, only: pose_cont_config, pose_cont_transaction_result, &
    &LM_ACCEPTED_IMPROVEMENT, LM_FINITE_NO_IMPROVEMENT, LM_NO_RELIABLE_UPDATE, &
    &LM_STEP_BOUND_REJECTED, LM_INVALID_NUMERICS, LM_ITERATION_LIMIT, &
    &POSE_CONT_INVALID_PREPARATION, &
    &POSE_CONT_ROUTE_SHIFT_THEN_JOINT, POSE_CONT_ROUTE_JOINT, &
    &POSE_CONT_OBJECTIVE_CART_EUCLID, POSE_CONT_OBJECTIVE_CART_NCC
implicit none
private

#include "simple_local_flags.inc"

public :: pose_cont_run_stats

!> Thread-local pose_cont accounting reduced after the particle OpenMP loop.
type :: pose_cont_run_stats
    private
    integer :: attempted = 0
    integer :: improved = 0
    integer :: finite_no_improvement = 0
    integer :: bound_rejected = 0
    integer :: no_reliable_update = 0
    integer :: invalid_numerics = 0
    integer :: iteration_limit = 0
    integer :: invalid_preparation = 0
    integer :: unknown_status = 0
    integer :: iterations = 0
    integer :: proposals = 0
    integer :: accepted_proposals = 0
    integer :: bound_hits = 0
    integer :: objective_pairs = 0
    real(dp) :: rotation_motion_sum = 0._dp
    real(dp) :: shift_motion_sum = 0._dp
    real(dp) :: rotation_motion_max = 0._dp
    real(dp) :: shift_motion_max = 0._dp
    real(dp) :: objective_before_sum = 0._dp
    real(dp) :: objective_after_sum = 0._dp
    real(dp) :: runtime_sum = 0._dp
    real(dp) :: runtime_max = 0._dp
contains
    procedure, public :: record => record_pose_cont_stats
    procedure, public :: merge => merge_pose_cont_stats
    procedure, public :: report => report_pose_cont_stats
end type pose_cont_run_stats

contains

    !> Add one completed transaction to a thread-private accumulator.
    subroutine record_pose_cont_stats(self, result, elapsed_seconds)
        class(pose_cont_run_stats), intent(inout) :: self
        type(pose_cont_transaction_result), intent(in) :: result
        real(dp), intent(in) :: elapsed_seconds

        self%attempted = self%attempted + 1
        select case (result%status)
        case (LM_ACCEPTED_IMPROVEMENT)
            self%improved = self%improved + 1
        case (LM_FINITE_NO_IMPROVEMENT)
            self%finite_no_improvement = self%finite_no_improvement + 1
        case (LM_STEP_BOUND_REJECTED)
            self%bound_rejected = self%bound_rejected + 1
        case (LM_NO_RELIABLE_UPDATE)
            self%no_reliable_update = self%no_reliable_update + 1
        case (LM_INVALID_NUMERICS)
            self%invalid_numerics = self%invalid_numerics + 1
        case (LM_ITERATION_LIMIT)
            self%iteration_limit = self%iteration_limit + 1
        case (POSE_CONT_INVALID_PREPARATION)
            self%invalid_preparation = self%invalid_preparation + 1
        case default
            self%unknown_status = self%unknown_status + 1
        end select
        self%iterations = self%iterations + result%iterations
        self%proposals = self%proposals + result%attempts
        self%accepted_proposals = self%accepted_proposals + result%accepts
        self%bound_hits = self%bound_hits + result%bound_hits
        self%rotation_motion_sum = self%rotation_motion_sum + result%cumulative_rotation
        self%shift_motion_sum = self%shift_motion_sum + result%cumulative_shift
        self%rotation_motion_max = max(self%rotation_motion_max, result%cumulative_rotation)
        self%shift_motion_max = max(self%shift_motion_max, result%cumulative_shift)
        if (result%objective_before >= 0._dp .and. result%objective_after >= 0._dp) then
            self%objective_pairs = self%objective_pairs + 1
            self%objective_before_sum = self%objective_before_sum + result%objective_before
            self%objective_after_sum = self%objective_after_sum + result%objective_after
        end if
        self%runtime_sum = self%runtime_sum + elapsed_seconds
        self%runtime_max = max(self%runtime_max, elapsed_seconds)
    end subroutine record_pose_cont_stats

    !> Merge one thread-private accumulator into an iteration total.
    pure subroutine merge_pose_cont_stats(self, other)
        class(pose_cont_run_stats), intent(inout) :: self
        class(pose_cont_run_stats), intent(in) :: other

        self%attempted = self%attempted + other%attempted
        self%improved = self%improved + other%improved
        self%finite_no_improvement = self%finite_no_improvement + other%finite_no_improvement
        self%bound_rejected = self%bound_rejected + other%bound_rejected
        self%no_reliable_update = self%no_reliable_update + other%no_reliable_update
        self%invalid_numerics = self%invalid_numerics + other%invalid_numerics
        self%iteration_limit = self%iteration_limit + other%iteration_limit
        self%invalid_preparation = self%invalid_preparation + other%invalid_preparation
        self%unknown_status = self%unknown_status + other%unknown_status
        self%iterations = self%iterations + other%iterations
        self%proposals = self%proposals + other%proposals
        self%accepted_proposals = self%accepted_proposals + other%accepted_proposals
        self%bound_hits = self%bound_hits + other%bound_hits
        self%objective_pairs = self%objective_pairs + other%objective_pairs
        self%rotation_motion_sum = self%rotation_motion_sum + other%rotation_motion_sum
        self%shift_motion_sum = self%shift_motion_sum + other%shift_motion_sum
        self%rotation_motion_max = max(self%rotation_motion_max, other%rotation_motion_max)
        self%shift_motion_max = max(self%shift_motion_max, other%shift_motion_max)
        self%objective_before_sum = self%objective_before_sum + other%objective_before_sum
        self%objective_after_sum = self%objective_after_sum + other%objective_after_sum
        self%runtime_sum = self%runtime_sum + other%runtime_sum
        self%runtime_max = max(self%runtime_max, other%runtime_max)
    end subroutine merge_pose_cont_stats

    !> Print and persist one compact iteration-level pose_cont summary.
    subroutine report_pose_cont_stats(self, config, iteration)
        class(pose_cont_run_stats), intent(in) :: self
        type(pose_cont_config), intent(in) :: config
        integer, intent(in) :: iteration
        character(len=32) :: route, objective_name
        character(len=64) :: filename
        real(dp) :: denominator, objective_denominator
        integer :: unit, terminal_count, invalid_unreliable

        invalid_unreliable = self%no_reliable_update + self%invalid_numerics + &
            &self%iteration_limit + self%invalid_preparation + self%unknown_status
        terminal_count = self%improved + self%finite_no_improvement + &
            &self%bound_rejected + invalid_unreliable
        if (terminal_count /= self%attempted) &
            &THROW_HARD('pose_cont terminal accounting does not balance')
        select case (config%route)
        case (POSE_CONT_ROUTE_SHIFT_THEN_JOINT)
            route = 'shift_then_joint'
        case (POSE_CONT_ROUTE_JOINT)
            route = 'joint'
        case default
            route = 'invalid'
        end select
        select case (config%objective)
        case (POSE_CONT_OBJECTIVE_CART_EUCLID)
            objective_name = 'cart_euclid'
        case (POSE_CONT_OBJECTIVE_CART_NCC)
            objective_name = 'cart_ncc'
        case default
            objective_name = 'invalid'
        end select
        denominator = real(max(self%attempted, 1), dp)
        objective_denominator = real(max(self%objective_pairs, 1), dp)
        write (logfhandle, '(A,1X,A,1X,A,1X,A,I0,A,I0,A,I0,A,I0,A,I0)') '>>> POSE_CONT', &
            &trim(route), trim(objective_name), 'attempted=', self%attempted, &
            &' improved=', self%improved, ' finite_no_improvement=', &
            &self%finite_no_improvement, ' bound_rejected=', self%bound_rejected, &
            &' invalid_unreliable=', invalid_unreliable
        write (logfhandle, '(A,5(A,I0))') '>>> POSE_CONT invalid detail', &
            &' no_reliable_update=', self%no_reliable_update, &
            &' invalid_numerics=', self%invalid_numerics, &
            &' iteration_limit=', self%iteration_limit, &
            &' invalid_preparation=', self%invalid_preparation, &
            &' unknown_status=', self%unknown_status
        write (logfhandle, '(A,I0,A,I0,A,I0,A,I0,A,4(ES12.4,1X))') '>>> POSE_CONT iterations=', &
            &self%iterations, ' proposals=', self%proposals, &
            &' accepted_proposals=', self%accepted_proposals, ' bound_hits=', self%bound_hits, &
            &' motion mean/max rotation shift=', self%rotation_motion_sum/denominator, &
            &self%rotation_motion_max, self%shift_motion_sum/denominator, self%shift_motion_max
        write (logfhandle, '(A,2(ES12.4,1X),A,3(ES12.4,1X))') '>>> POSE_CONT objective mean before/after=', &
            &self%objective_before_sum/objective_denominator, &
            &self%objective_after_sum/objective_denominator, &
            &'transaction seconds sum/mean/max=', self%runtime_sum, &
            &self%runtime_sum/denominator, self%runtime_max
        write (filename, '(A,I3.3,A)') 'POSE_CONT_STATS_ITER', iteration, '.txt'
        open (newunit=unit, file=trim(filename), status='replace', action='write')
        write (unit, '(A)') 'route='//trim(route)
        write (unit, '(A)') 'objective='//trim(objective_name)
        write (unit, '(A,I0)') 'attempted=', self%attempted
        write (unit, '(A,I0)') 'improved=', self%improved
        write (unit, '(A,I0)') 'finite_no_improvement=', self%finite_no_improvement
        write (unit, '(A,I0)') 'bound_rejected=', self%bound_rejected
        write (unit, '(A,I0)') 'invalid_unreliable=', invalid_unreliable
        write (unit, '(A,I0)') 'no_reliable_update=', self%no_reliable_update
        write (unit, '(A,I0)') 'invalid_numerics=', self%invalid_numerics
        write (unit, '(A,I0)') 'iteration_limit=', self%iteration_limit
        write (unit, '(A,I0)') 'invalid_preparation=', self%invalid_preparation
        write (unit, '(A,I0)') 'unknown_status=', self%unknown_status
        write (unit, '(A,I0)') 'iterations=', self%iterations
        write (unit, '(A,I0)') 'proposals=', self%proposals
        write (unit, '(A,I0)') 'accepted_proposals=', self%accepted_proposals
        write (unit, '(A,I0)') 'bound_hits=', self%bound_hits
        write (unit, '(A,I0)') 'objective_pairs=', self%objective_pairs
        write (unit, '(A,ES16.8)') 'mean_objective_before=', &
            &self%objective_before_sum/objective_denominator
        write (unit, '(A,ES16.8)') 'mean_objective_after=', &
            &self%objective_after_sum/objective_denominator
        write (unit, '(A,ES16.8)') 'mean_rotation_motion=', self%rotation_motion_sum/denominator
        write (unit, '(A,ES16.8)') 'max_rotation_motion=', self%rotation_motion_max
        write (unit, '(A,ES16.8)') 'mean_shift_motion=', self%shift_motion_sum/denominator
        write (unit, '(A,ES16.8)') 'max_shift_motion=', self%shift_motion_max
        ! Sum is aggregate worker time, not parallel wall-clock elapsed time.
        write (unit, '(A,ES16.8)') 'transaction_seconds_sum=', self%runtime_sum
        write (unit, '(A,ES16.8)') 'transaction_seconds_mean=', self%runtime_sum/denominator
        write (unit, '(A,ES16.8)') 'transaction_seconds_max=', self%runtime_max
        close (unit)
        write (logfhandle, '(A)') '>>> POSE_CONT detailed statistics written to '//trim(filename)
    end subroutine report_pose_cont_stats

end module simple_pose_cont_run_stats
