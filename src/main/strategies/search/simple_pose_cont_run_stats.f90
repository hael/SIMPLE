!@descr: thread-local accumulation and iteration reporting for pose_cont refinement
module simple_pose_cont_run_stats
use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
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

integer, parameter :: POSE_CONT_STATS_FORMAT_VERSION = 1

public :: pose_cont_run_stats, aggregate_pose_cont_stats_files

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
    integer :: objective_increased = 0
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
        if (ieee_is_finite(result%cumulative_rotation) .and. result%cumulative_rotation >= 0._dp) then
            self%rotation_motion_sum = self%rotation_motion_sum + result%cumulative_rotation
            self%rotation_motion_max = max(self%rotation_motion_max, result%cumulative_rotation)
        end if
        if (ieee_is_finite(result%cumulative_shift) .and. result%cumulative_shift >= 0._dp) then
            self%shift_motion_sum = self%shift_motion_sum + result%cumulative_shift
            self%shift_motion_max = max(self%shift_motion_max, result%cumulative_shift)
        end if
        if (ieee_is_finite(result%objective_before) .and. &
            &ieee_is_finite(result%objective_after) .and. &
            &result%objective_before >= 0._dp .and. result%objective_after >= 0._dp) then
            self%objective_pairs = self%objective_pairs + 1
            self%objective_before_sum = self%objective_before_sum + result%objective_before
            self%objective_after_sum = self%objective_after_sum + result%objective_after
            if (materially_increased(result%objective_before, result%objective_after)) &
                &self%objective_increased = self%objective_increased + 1
        end if
        if (ieee_is_finite(elapsed_seconds) .and. elapsed_seconds >= 0._dp) then
            self%runtime_sum = self%runtime_sum + elapsed_seconds
            self%runtime_max = max(self%runtime_max, elapsed_seconds)
        end if
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
        self%objective_increased = self%objective_increased + other%objective_increased
        self%rotation_motion_sum = self%rotation_motion_sum + other%rotation_motion_sum
        self%shift_motion_sum = self%shift_motion_sum + other%shift_motion_sum
        self%rotation_motion_max = max(self%rotation_motion_max, other%rotation_motion_max)
        self%shift_motion_max = max(self%shift_motion_max, other%shift_motion_max)
        self%objective_before_sum = self%objective_before_sum + other%objective_before_sum
        self%objective_after_sum = self%objective_after_sum + other%objective_after_sum
        self%runtime_sum = self%runtime_sum + other%runtime_sum
        self%runtime_max = max(self%runtime_max, other%runtime_max)
    end subroutine merge_pose_cont_stats

    !> Persist one process-private report, or the final report for a one-part run.
    subroutine report_pose_cont_stats(self, config, iteration, part, nparts)
        class(pose_cont_run_stats), intent(in) :: self
        type(pose_cont_config), intent(in) :: config
        integer, intent(in) :: iteration, part, nparts
        character(len=64) :: filename

        if (nparts == 1) then
            if (part /= 1) THROW_HARD('invalid pose_cont statistics part index')
            write (filename, '(A,I3.3,A)') 'POSE_CONT_STATS_ITER', iteration, '.txt'
            call write_pose_cont_summary(self, config, filename, iteration, nparts)
        else
            if (part < 1 .or. part > nparts) &
                &THROW_HARD('invalid pose_cont statistics part index')
            write (filename, '(A,I3.3,A,I3.3,A)') &
                &'POSE_CONT_STATS_ITER', iteration, '_PART', part, '.txt'
            call write_pose_cont_partition(self, config, filename, iteration, part, nparts)
        end if
    end subroutine report_pose_cont_stats

    !> Merge process-private reports after all distributed workers finish.
    subroutine aggregate_pose_cont_stats_files(config, iteration, nparts)
        type(pose_cont_config), intent(in) :: config
        integer, intent(in) :: iteration, nparts
        type(pose_cont_run_stats) :: total, part_stats
        character(len=64) :: filename
        integer :: part

        if (nparts == 1) return
        do part = 1, nparts
            write (filename, '(A,I3.3,A,I3.3,A)') &
                &'POSE_CONT_STATS_ITER', iteration, '_PART', part, '.txt'
            call read_pose_cont_partition(filename, config, iteration, part, nparts, part_stats)
            call total%merge(part_stats)
        end do
        write (filename, '(A,I3.3,A)') 'POSE_CONT_STATS_ITER', iteration, '.txt'
        call write_pose_cont_summary(total, config, filename, iteration, nparts)
        do part = 1, nparts
            write (filename, '(A,I3.3,A,I3.3,A)') &
                &'POSE_CONT_STATS_ITER', iteration, '_PART', part, '.txt'
            call delete_pose_cont_partition(filename)
        end do
    end subroutine aggregate_pose_cont_stats_files

    subroutine write_pose_cont_summary(self, config, filename, iteration, nparts)
        class(pose_cont_run_stats), intent(in) :: self
        type(pose_cont_config), intent(in) :: config
        character(len=*), intent(in) :: filename
        integer, intent(in) :: iteration, nparts
        character(len=32) :: route, objective_name
        real(dp) :: denominator, objective_denominator
        integer :: unit, invalid_unreliable

        call validate_pose_cont_stats(self)
        invalid_unreliable = pose_cont_invalid_count(self)
        call pose_cont_stats_labels(config, route, objective_name)
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
        write (logfhandle, '(A,2(ES12.4,1X),A,3(ES12.4,1X))') &
            &'>>> POSE_CONT objective mean before/after=', &
            &self%objective_before_sum/objective_denominator, &
            &self%objective_after_sum/objective_denominator, &
            &'transaction seconds sum/mean/max=', self%runtime_sum, &
            &self%runtime_sum/denominator, self%runtime_max
        if (self%objective_increased > 0) write (logfhandle, '(A,I0)') &
            &'>>> POSE_CONT WARNING: materially increased terminal objectives=', &
            &self%objective_increased
        open (newunit=unit, file=trim(filename), status='replace', action='write')
        write (unit, '(A,I0)') 'format_version=', POSE_CONT_STATS_FORMAT_VERSION
        write (unit, '(A,I0)') 'iteration=', iteration
        write (unit, '(A,I0)') 'nparts=', nparts
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
        write (unit, '(A,I0)') 'objective_increased=', self%objective_increased
        write (unit, '(A,ES24.16)') 'mean_objective_before=', &
            &self%objective_before_sum/objective_denominator
        write (unit, '(A,ES24.16)') 'mean_objective_after=', &
            &self%objective_after_sum/objective_denominator
        write (unit, '(A,ES24.16)') 'mean_rotation_motion=', self%rotation_motion_sum/denominator
        write (unit, '(A,ES24.16)') 'max_rotation_motion=', self%rotation_motion_max
        write (unit, '(A,ES24.16)') 'mean_shift_motion=', self%shift_motion_sum/denominator
        write (unit, '(A,ES24.16)') 'max_shift_motion=', self%shift_motion_max
        ! Sum is aggregate worker time, not parallel wall-clock elapsed time.
        write (unit, '(A,ES24.16)') 'transaction_seconds_sum=', self%runtime_sum
        write (unit, '(A,ES24.16)') 'transaction_seconds_mean=', self%runtime_sum/denominator
        write (unit, '(A,ES24.16)') 'transaction_seconds_max=', self%runtime_max
        close (unit)
        write (logfhandle, '(A)') '>>> POSE_CONT detailed statistics written to '//trim(filename)
    end subroutine write_pose_cont_summary

    subroutine write_pose_cont_partition(self, config, filename, iteration, part, nparts)
        class(pose_cont_run_stats), intent(in) :: self
        type(pose_cont_config), intent(in) :: config
        character(len=*), intent(in) :: filename
        integer, intent(in) :: iteration, part, nparts
        character(len=32) :: route, objective_name
        integer :: unit

        call validate_pose_cont_stats(self)
        call pose_cont_stats_labels(config, route, objective_name)
        open (newunit=unit, file=trim(filename), status='replace', action='write')
        write (unit, '(4(I0,1X),2(A,1X))') POSE_CONT_STATS_FORMAT_VERSION, iteration, part, &
            &nparts, trim(route), trim(objective_name)
        write (unit, '(15(I0,1X))') self%attempted, self%improved, &
            &self%finite_no_improvement, self%bound_rejected, self%no_reliable_update, &
            &self%invalid_numerics, self%iteration_limit, self%invalid_preparation, &
            &self%unknown_status, self%iterations, self%proposals, self%accepted_proposals, &
            &self%bound_hits, self%objective_pairs, self%objective_increased
        write (unit, '(8(ES24.16,1X))') self%rotation_motion_sum, self%shift_motion_sum, &
            &self%rotation_motion_max, self%shift_motion_max, self%objective_before_sum, &
            &self%objective_after_sum, self%runtime_sum, self%runtime_max
        close (unit)
    end subroutine write_pose_cont_partition

    subroutine read_pose_cont_partition(filename, config, expected_iteration, expected_part, &
        &expected_nparts, stats)
        character(len=*), intent(in) :: filename
        type(pose_cont_config), intent(in) :: config
        integer, intent(in) :: expected_iteration, expected_part, expected_nparts
        type(pose_cont_run_stats), intent(out) :: stats
        character(len=32) :: route, objective_name, file_route, file_objective
        integer :: unit, ios, format_version, iteration, part, nparts

        stats = pose_cont_run_stats()
        open (newunit=unit, file=trim(filename), status='old', action='read', iostat=ios)
        if (ios /= 0) THROW_HARD('missing pose_cont partition statistics')
        read (unit, *, iostat=ios) format_version, iteration, part, nparts, &
            &file_route, file_objective
        if (ios /= 0) THROW_HARD('invalid pose_cont partition statistics header')
        call pose_cont_stats_labels(config, route, objective_name)
        if (format_version /= POSE_CONT_STATS_FORMAT_VERSION .or. &
            &iteration /= expected_iteration .or. part /= expected_part .or. &
            &nparts /= expected_nparts .or. trim(file_route) /= trim(route) .or. &
            &trim(file_objective) /= trim(objective_name)) &
            &THROW_HARD('mismatched pose_cont partition statistics metadata')
        read (unit, *, iostat=ios) stats%attempted, stats%improved, &
            &stats%finite_no_improvement, stats%bound_rejected, stats%no_reliable_update, &
            &stats%invalid_numerics, stats%iteration_limit, stats%invalid_preparation, &
            &stats%unknown_status, stats%iterations, stats%proposals, stats%accepted_proposals, &
            &stats%bound_hits, stats%objective_pairs, stats%objective_increased
        if (ios /= 0) THROW_HARD('invalid pose_cont partition statistics counts')
        read (unit, *, iostat=ios) stats%rotation_motion_sum, stats%shift_motion_sum, &
            &stats%rotation_motion_max, stats%shift_motion_max, stats%objective_before_sum, &
            &stats%objective_after_sum, stats%runtime_sum, stats%runtime_max
        if (ios /= 0) THROW_HARD('invalid pose_cont partition statistics values')
        call validate_pose_cont_stats(stats)
        close (unit)
    end subroutine read_pose_cont_partition

    subroutine delete_pose_cont_partition(filename)
        character(len=*), intent(in) :: filename
        integer :: unit, ios

        open (newunit=unit, file=trim(filename), status='old', action='read', iostat=ios)
        if (ios /= 0) THROW_HARD('failed reopening pose_cont partition statistics')
        close (unit, status='delete')
    end subroutine delete_pose_cont_partition

    pure integer function pose_cont_invalid_count(stats) result(count)
        class(pose_cont_run_stats), intent(in) :: stats

        count = stats%no_reliable_update + stats%invalid_numerics + stats%iteration_limit + &
            &stats%invalid_preparation + stats%unknown_status
    end function pose_cont_invalid_count

    subroutine validate_pose_cont_stats(stats)
        class(pose_cont_run_stats), intent(in) :: stats
        integer :: terminal_count

        terminal_count = stats%improved + stats%finite_no_improvement + stats%bound_rejected + &
            &pose_cont_invalid_count(stats)
        if (terminal_count /= stats%attempted) &
            &THROW_HARD('pose_cont terminal accounting does not balance')
        if (.not. all(ieee_is_finite([stats%rotation_motion_sum, stats%shift_motion_sum, &
            &stats%rotation_motion_max, stats%shift_motion_max, stats%objective_before_sum, &
            &stats%objective_after_sum, stats%runtime_sum, stats%runtime_max]))) &
            &THROW_HARD('pose_cont statistics contain non-finite values')
    end subroutine validate_pose_cont_stats

    pure subroutine pose_cont_stats_labels(config, route, objective_name)
        type(pose_cont_config), intent(in) :: config
        character(len=*), intent(out) :: route, objective_name

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
    end subroutine pose_cont_stats_labels

    pure logical function materially_increased(before, after) result(increased)
        real(dp), intent(in) :: before, after
        real(dp) :: tolerance

        tolerance = max(1.e-10_dp, 1.e-8_dp*abs(before))
        increased = after > before + tolerance
    end function materially_increased

end module simple_pose_cont_run_stats
