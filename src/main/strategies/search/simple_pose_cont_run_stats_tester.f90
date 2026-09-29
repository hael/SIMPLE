!@descr: unit tests for pose_cont run-statistics aggregation (simple_pose_cont_run_stats)
module simple_pose_cont_run_stats_tester
use simple_defs, only: dp
use simple_pose_cont_refine3D_adapter, only: LM_ACCEPTED_IMPROVEMENT, pose_cont_config, &
    &pose_cont_config_from_route, pose_cont_transaction_result
use simple_pose_cont_run_stats, only: pose_cont_run_stats, aggregate_pose_cont_stats_files
use simple_syslib, only: del_file, file_exists
use simple_test_utils
implicit none
private
public :: run_all_pose_cont_run_stats_tests

contains

    subroutine run_all_pose_cont_run_stats_tests()
        write (*, '(A)') '**** running all pose_cont run-statistics tests ****'
        write (*, '(A)') 'test_partition_aggregation'
        call test_partition_aggregation()
    end subroutine run_all_pose_cont_run_stats_tests

    subroutine test_partition_aggregation()
        type(pose_cont_run_stats) :: first, second
        type(pose_cont_transaction_result) :: result
        type(pose_cont_config) :: config
        real(dp) :: value

        call remove_stats_fixture()
        config = pose_cont_config_from_route('joint')
        result%status = LM_ACCEPTED_IMPROVEMENT
        result%objective_before = 0.123456789012345_dp
        result%objective_after = 0.103456789012345_dp
        result%cumulative_rotation = 4.e-4_dp
        result%cumulative_shift = 5.e-3_dp
        call first%record(result, 0.1_dp)
        result%objective_before = 0.987654321098765_dp
        result%objective_after = 0.907654321098765_dp
        result%cumulative_rotation = 6.e-4_dp
        result%cumulative_shift = 7.e-3_dp
        call second%record(result, 0.2_dp)

        call first%report(config, 997, 1, 2)
        call second%report(config, 997, 2, 2)
        call aggregate_pose_cont_stats_files(config, 997, 2)
        call read_stats_value('POSE_CONT_STATS_ITER997.txt', 'attempted', value)
        call assert_int(2, nint(value), 'aggregated pose-cont attempted count')
        call read_stats_value('POSE_CONT_STATS_ITER997.txt', 'objective_pairs', value)
        call assert_int(2, nint(value), 'aggregated pose-cont objective count')
        call read_stats_value('POSE_CONT_STATS_ITER997.txt', 'objective_increased', value)
        call assert_int(0, nint(value), 'aggregated pose-cont objective increases')
        call read_stats_value('POSE_CONT_STATS_ITER997.txt', 'mean_objective_before', value)
        call assert_double(0.555555555055555_dp, value, &
            &'aggregated pose-cont mean objective before', ulp_tol=1.e4_dp)
        call read_stats_value('POSE_CONT_STATS_ITER997.txt', 'mean_objective_after', value)
        call assert_double(0.505555555055555_dp, value, &
            &'aggregated pose-cont mean objective after', ulp_tol=1.e4_dp)
        call assert_false(file_exists('POSE_CONT_STATS_ITER997_PART001.txt'), &
            &'pose-cont aggregation retained partition 1 statistics')
        call assert_false(file_exists('POSE_CONT_STATS_ITER997_PART002.txt'), &
            &'pose-cont aggregation retained partition 2 statistics')
        call remove_stats_fixture()
    end subroutine test_partition_aggregation

    subroutine read_stats_value(filename, requested_key, value)
        character(len=*), intent(in) :: filename, requested_key
        real(dp), intent(out) :: value
        character(len=256) :: line
        integer :: unit, ios, separator

        value = 0._dp
        open (newunit=unit, file=trim(filename), status='old', action='read', iostat=ios)
        call assert_int(0, ios, 'pose-cont aggregate report opens')
        if (ios /= 0) return
        do
            read (unit, '(A)', iostat=ios) line
            if (ios /= 0) exit
            separator = index(line, '=')
            if (separator < 2) cycle
            if (trim(line(:separator - 1)) /= trim(requested_key)) cycle
            read (line(separator + 1:), *, iostat=ios) value
            call assert_int(0, ios, 'pose-cont aggregate value parses')
            close (unit)
            return
        end do
        close (unit)
        call assert_true(.false., 'pose-cont aggregate value is present')
    end subroutine read_stats_value

    subroutine remove_stats_fixture()
        call del_file('POSE_CONT_STATS_ITER997.txt')
        call del_file('POSE_CONT_STATS_ITER997_PART001.txt')
        call del_file('POSE_CONT_STATS_ITER997_PART002.txt')
    end subroutine remove_stats_fixture

end module simple_pose_cont_run_stats_tester
