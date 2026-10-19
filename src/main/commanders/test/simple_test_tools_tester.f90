!@descr: unit tests for the test tooling: profile statistics (simple_test_profile_stats), packing by CPU reservation (simple_test_dispatch), identity rules (simple_test_identity), the registry (simple_test_registry)
! The expected values are worked out by hand in the comments; files are written as tmp_test_tools_*
! and removed at the end of each test.
module simple_test_tools_tester
use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan
use simple_string,             only: string
use simple_syslib,             only: simple_mkdir, simple_rmdir, del_file
use simple_test_utils,         only: assert_true, assert_false, assert_int, assert_real
use simple_test_profile_stats, only: gate_row, metric_dist, read_gate_table, robust_spread, summarize_metric, profile_table
use simple_test_dispatch,      only: next_to_start
use simple_test_identity,      only: compare_run_dirs, same_bytes, same_manifest, manifest_line_kept, SIGMA2_DIGEST_BYTES
use simple_test_registry,      only: test_entry, get_test_registry, find_test_entry
implicit none

public :: run_all_test_tools_tests
private

contains

    subroutine run_all_test_tools_tests()
        write(*,'(A)') '**** running all test tools tests ****'
        call test_robust_spread()
        call test_summarize_metric()
        call test_gate_tables()
        call test_packing()
        call test_manifest_rules()
        call test_masked_bytes()
        call test_run_dir_identity()
        call test_registry()
    end subroutine run_all_test_tools_tests

    !> s = 1.4826 MAD, or the range over 4 when the MAD is zero
    subroutine test_robust_spread()
        write(*,'(A)') 'test_robust_spread'
        ! median 3, absolute deviations 2 1 0 1 97: MAD 1
        call assert_real(1.4826, robust_spread([1., 2., 3., 4., 100.]), 1.e-5, 'robust spread is 1.4826 times the MAD')
        ! median 1, MAD 0, range 4
        call assert_real(1.0, robust_spread([1., 1., 1., 5.]), 1.e-6, 'a zero MAD falls back to the range over 4')
        call assert_real(0.0, robust_spread([2., 2., 2.]), 0., 'a constant sample has no spread')
        call assert_real(0.0, robust_spread([7.]), 0., 'one run has no spread')
    end subroutine test_robust_spread

    !> worst in the sense of the floor; failure rate over all executions, the missing ones failing
    subroutine test_summarize_metric()
        type(metric_dist) :: d
        real :: nan
        write(*,'(A)') 'test_summarize_metric'
        d = summarize_metric('err', 'max', [0.5, 0.9, 0.7], [1., 1., 1.], [.true., .false., .true.], 4)
        call assert_real(0.7, d%med, 1.e-6, 'the median of the runs')
        call assert_real(0.9, d%worst, 1.e-6, 'the worst of a ceiling metric is its largest value')
        call assert_real(0.5, d%failure_rate, 1.e-6, 'one failed and one missing execution of four fail half')
        call assert_int(3, d%runs, 'three executions wrote it')
        d = summarize_metric('corr', 'min', [0.5, 0.9, 0.7], [0.4, 0.4, 0.4], [.true., .true., .true.], 3)
        call assert_real(0.5, d%worst, 1.e-6, 'the worst of a floor metric is its smallest value')
        call assert_real(0.0, d%failure_rate, 0., 'no failure when every execution passed')
        nan = ieee_value(0., ieee_quiet_nan)
        d = summarize_metric('corr', 'min', [0.5, nan], [0.4, 0.4], [.true., .false.], 2)
        call assert_real(0.5, d%med, 1.e-6, 'a non-finite value stays out of the statistics')
        call assert_real(0.5, d%failure_rate, 1.e-6, 'and counts as a failure')
    end subroutine test_summarize_metric

    !> gate tables as test_gate writes them, aggregated over executions
    subroutine test_gate_tables()
        character(len=*), parameter :: F1 = 'tmp_test_tools_1.tsv', F2 = 'tmp_test_tools_2.tsv'
        character(len=*), parameter :: FMISSING = 'tmp_test_tools_missing.tsv'
        type(gate_row),    allocatable :: rows(:)
        type(metric_dist), allocatable :: dists(:)
        write(*,'(A)') 'test_gate_tables'
        call write_table(F1, [character(len=80) :: 'corr'//achar(9)//' 9.000000E-01'//achar(9)//' 8.000000E-01'//achar(9)//'yes'//achar(9)//'min', &
                             &'err '//achar(9)//' 1.000000E+00'//achar(9)//' 2.000000E+00'//achar(9)//'yes'//achar(9)//'max', &
                             &'err '//achar(9)//' 3.000000E+00'//achar(9)//' 2.000000E+00'//achar(9)//'no '//achar(9)//'max'])
        call write_table(F2, [character(len=80) :: 'corr'//achar(9)//' 7.000000E-01'//achar(9)//' 8.000000E-01'//achar(9)//'no '//achar(9)//'min', &
                             &'rec '//achar(9)//' 5.000000E-01'//achar(9)//'-1.000000E+00'//achar(9)//'yes'//achar(9)//'report', &
                             &'rec '//achar(9)//' 5.000000E-01'//achar(9)//' 9.000000E-01'//achar(9)//'no '//achar(9)//'min'])
        call read_gate_table(F1, rows)
        call assert_int(3, size(rows), 'every row of the table is read')
        call assert_true(rows(1)%pass .and. .not. rows(3)%pass, 'the pass column is read')
        call assert_true(trim(rows(2)%sense) == 'max', 'the sense column is read')
        call assert_real(2.0, rows(2)%floor, 1.e-6, 'the floor column is read')
        call profile_table([character(len=32) :: F1, F2, FMISSING], [.true., .true., .true.], dists)
        call assert_int(3, size(dists), 'one distribution per metric name')
        ! corr: written by both, failed in the second, missing in the third execution
        call assert_int(2, dists(1)%runs, 'corr was written by two executions')
        call assert_real(2./3., dists(1)%failure_rate, 1.e-5, 'corr failed or was missing in two of three')
        call assert_real(0.8, dists(1)%med, 1.e-5, 'the median of 0.9 and 0.7')
        ! err: written twice by the first execution, one row failing: that execution fails it
        call assert_int(1, dists(2)%runs, 'err was written by one execution')
        call assert_real(1.0, dists(2)%failure_rate, 1.e-6, 'an execution with one failing row of a metric fails it')
        call assert_real(3.0, dists(2)%worst, 1.e-6, 'the worst of its rows')
        ! rec: reported, then held to a floor it misses
        call assert_true(trim(dists(3)%sense) == 'min', 'a reported and gated metric takes the sense of its floor')
        call assert_real(0.9, dists(3)%floor, 1.e-6, 'and its floor')
        call assert_real(1.0, dists(3)%failure_rate, 1.e-6, 'and fails where its gated row fails')
        ! an execution that failed or timed out after writing passing rows fails its metrics
        call profile_table([character(len=32) :: F1, F1], [.true., .false.], dists)
        call assert_int(2, dists(2)%runs, 'both executions wrote err')
        call assert_real(1.0, dists(2)%failure_rate, 1.e-6, 'err fails in both (a failing row; a failed execution)')
        call assert_real(0.5, dists(1)%failure_rate, 1.e-6, &
            &'corr passes only where the execution passed: a failed execution fails its passing rows')
        call del_file(F1)
        call del_file(F2)
    end subroutine test_gate_tables

    !> the largest queued reservation starts first, when its cores, memory and disk fit; smaller ones
    !! wait for it; one too large runs alone
    subroutine test_packing()
        real, parameter :: M0(4) = 0., BIG = 1.e6
        write(*,'(A)') 'test_packing'
        call assert_int(1, next_to_start([8, 8, 8, 2], M0, M0, [.true., .true., .true., .true.], 24, BIG, BIG, 0), &
            &'the first in queue order among the largest')
        call assert_int(3, next_to_start([2, 4, 8, 8], M0, M0, [.true., .true., .true., .true.], 24, BIG, BIG, 0), &
            &'the largest reservation first')
        call assert_int(0, next_to_start([8, 8, 8, 2], M0, M0, [.false., .true., .true., .true.], 6, BIG, BIG, 2), &
            &'a smaller one does not start before a larger one that waits for cores')
        call assert_int(4, next_to_start([8, 8, 8, 2], M0, M0, [.false., .false., .false., .true.], 6, BIG, BIG, 3), &
            &'it starts once every larger one has')
        call assert_int(0, next_to_start([4, 4], [2., 8.], M0(1:2), [.false., .true.], 24, 5., BIG, 1), &
            &'nothing starts when its memory does not fit')
        call assert_int(0, next_to_start([4, 4], M0(1:2), [2., 8.], [.false., .true.], 24, BIG, 5., 1), &
            &'nor when its disk does not')
        call assert_int(2, next_to_start([30, 30], M0(1:2), M0(1:2), [.false., .true.], 24, BIG, BIG, 0), &
            &'a reservation above the budget runs alone')
        call assert_int(0, next_to_start([30, 30], M0(1:2), M0(1:2), [.false., .true.], 0, BIG, BIG, 1), &
            &'but not beside another')
        call assert_int(0, next_to_start([1, 1], M0(1:2), M0(1:2), [.false., .false.], 24, BIG, BIG, 0), &
            &'nothing queued, nothing to start')
    end subroutine test_packing

    !> the run-specific manifest lines go, artifact lines lose their path
    subroutine test_manifest_rules()
        character(len=*), parameter :: M1 = 'tmp_test_tools_m1.txt', M2 = 'tmp_test_tools_m2.txt'
        character(len=*), parameter :: M3 = 'tmp_test_tools_m3.txt'
        character(len=128) :: line
        write(*,'(A)') 'test_manifest_rules'
        line = 'run_id solve3D_20261016T203315.654_271529'
        call assert_false(manifest_line_kept(line), 'the run identifier is run-specific')
        line = 'stack_digest -6807008755386900676'
        call assert_false(manifest_line_kept(line), 'the stack digest is run-specific')
        line = 'nrows 51'
        call assert_true(manifest_line_kept(line), 'other lines are compared')
        line = 'artifact vol 1 1906378816907924071 /a/run/dir/recvol_state01.mrc'
        call assert_true(manifest_line_kept(line), 'artifact lines are compared')
        call assert_true(trim(line) == 'artifact vol 1 1906378816907924071', 'without their path')
        call write_table(M1, [character(len=20) :: 'run_id a', 'nrows 51', 'artifact x 1 /a', 'checksum 1'])
        call write_table(M2, [character(len=20) :: 'run_id b', 'nrows 51', 'artifact x 1 /b', 'checksum 2'])
        call write_table(M3, [character(len=20) :: 'run_id a', 'nrows 52', 'artifact x 1 /a', 'checksum 1'])
        call assert_true(same_manifest(M1, M2), 'manifests differing only in run-specific fields are the same')
        call assert_false(same_manifest(M1, M3), 'a differing row count is a difference')
        call del_file(M1)
        call del_file(M2)
        call del_file(M3)
    end subroutine test_manifest_rules

    !> a byte range excepted from the comparison, and only that range
    subroutine test_masked_bytes()
        character(len=*), parameter :: B1 = 'tmp_test_tools_b1.bin', B2 = 'tmp_test_tools_b2.bin'
        character(len=*), parameter :: B3 = 'tmp_test_tools_b3.bin'
        character(len=200) :: bytes
        write(*,'(A)') 'test_masked_bytes'
        bytes = repeat('a', 200)
        call write_bytes(B1, bytes)
        bytes(121:128) = 'DIGEST!!'
        call write_bytes(B2, bytes)
        bytes(129:129) = 'z'
        call write_bytes(B3, bytes)
        call assert_false(same_bytes(B1, B2), 'the digest bytes differ')
        call assert_true(same_bytes(B1, B2, SIGMA2_DIGEST_BYTES), 'and are excepted')
        call assert_false(same_bytes(B1, B3, SIGMA2_DIGEST_BYTES), 'a byte beside them is not')
        call del_file(B1)
        call del_file(B2)
        call del_file(B3)
    end subroutine test_masked_bytes

    !> the declared patterns select the files; a missing file counts as a difference
    subroutine test_run_dir_identity()
        character(len=*), parameter :: D1 = 'tmp_test_tools_d1', D2 = 'tmp_test_tools_d2'
        character(len=*), parameter :: REPORT = 'tmp_test_tools_report.txt'
        integer :: nsame, ndiff
        write(*,'(A)') 'test_run_dir_identity'
        call simple_mkdir(D1)
        call simple_mkdir(D2)
        call simple_mkdir(D1//'/sub')
        call simple_mkdir(D2//'/sub')
        call write_bytes(D1//'/sub/same.mrc',  'identical')
        call write_bytes(D2//'/sub/same.mrc',  'identical')
        call write_bytes(D1//'/sub/diff.mrc',  'one')
        call write_bytes(D2//'/sub/diff.mrc',  'two')
        call write_bytes(D1//'/sub/only1.mrc', 'x')
        call write_bytes(D1//'/sub/other.txt', 'not declared')
        call compare_run_dirs(D1, D2, 'sub/*.mrc', REPORT, nsame, ndiff)
        call assert_int(1, nsame, 'one declared file is the same')
        call assert_int(2, ndiff, 'one differs and one is missing in the second directory')
        call compare_run_dirs(D1, D2, 'sub/same.mrc', REPORT, nsame, ndiff)
        call assert_true(nsame == 1 .and. ndiff == 0, 'identical when the declared outputs agree')
        call simple_rmdir(D1)
        call simple_rmdir(D2)
        call del_file(REPORT)
    end subroutine test_run_dir_identity

    !> one row per CTest entry: unique names, a deadline, a reservation, outputs for the one-thread ones
    subroutine test_registry()
        type(test_entry), allocatable :: entries(:)
        integer :: i
        write(*,'(A)') 'test_registry'
        entries = get_test_registry()
        call assert_true(size(entries) > 0, 'the registry has entries')
        do i = 1, size(entries)
            call assert_int(i, find_test_entry(entries, trim(entries(i)%name)), 'entry names are unique: '//&
                &trim(entries(i)%name))
            call assert_true(len_trim(entries(i)%prg) > 0 .and. entries(i)%deadline > 0 .and. entries(i)%ncpu >= 1, &
                &'every entry has a program, a deadline and a reservation: '//trim(entries(i)%name))
            call assert_true(entries(i)%ncpu >= entries(i)%nworkers * entries(i)%wthr, &
                &'the reservation holds the workers: '//trim(entries(i)%name))
            if( entries(i)%nthr == 1 ) call assert_true(len_trim(entries(i)%outputs) > 0, &
                &'a one-thread entry declares its outputs: '//trim(entries(i)%name))
        enddo
        call assert_int(0, find_test_entry(entries, 'no_such_entry'), 'an unknown name is not found')
    end subroutine test_registry

    subroutine write_table( fname, lines )
        character(len=*), intent(in) :: fname, lines(:)
        integer :: funit, i
        open(newunit=funit, file=fname, action='write', status='replace')
        write(funit,'(a)') 'header'
        do i = 1, size(lines)
            write(funit,'(a)') trim(lines(i))
        enddo
        close(funit)
    end subroutine write_table

    subroutine write_bytes( fname, bytes )
        character(len=*), intent(in) :: fname, bytes
        integer :: funit
        open(newunit=funit, file=fname, access='stream', form='unformatted', action='write', status='replace')
        write(funit) bytes
        close(funit)
    end subroutine write_bytes

end module simple_test_tools_tester
