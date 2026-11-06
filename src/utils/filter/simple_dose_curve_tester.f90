!@descr: unit tests for the measured dose curve: its evaluation per movie fraction, its file and its checks
! Every expectation is a closed form of the knots written here.
module simple_dose_curve_tester
use simple_test_utils, only: assert_true, assert_false
use simple_defs,       only: dp
use simple_string,     only: string
use simple_syslib,     only: del_file
use simple_dose_curve, only: dose_curve, scaled_dose
implicit none
private
public :: run_all_dose_curve_tests

character(len=*), parameter :: FNAME = 'tmp_dose_curve_tester.txt'
real(dp),         parameter :: TOL   = 1.e-12_dp
! four windows of two fractions; centres at 1, 3, 5 and 7 when a fraction carries a unit scaled dose
real(dp),         parameter :: KNOT_C(4) = [1._dp, 0.9_dp, 0.8_dp, 0.6_dp]
real(dp),         parameter :: KNOT_B(4) = [0._dp, 10._dp, 30._dp, 60._dp]
real(dp),         parameter :: BETA = 10._dp, GAMMA = 0.04_dp

contains

    subroutine run_all_dose_curve_tests()
        write(*,'(A)') '**** running all dose curve tests ****'
        call test_fractions_between_windows()
        call test_fractions_past_the_windows()
        call test_rejected_window_is_skipped()
        call test_voltage_scaling()
        call test_file_round_trip()
        call test_movie_checks()
        call test_bad_files_are_rejected()
    end subroutine run_all_dose_curve_tests

    !> the curve of the knots above at voltage kv with a physical dose dpf per fraction
    subroutine make_curve( curve, kv, dpf, beta_line, gamma_line, used3 )
        class(dose_curve), intent(inout) :: curve
        real(dp),          intent(in)    :: kv, dpf, beta_line, gamma_line
        logical,           intent(in)    :: used3
        integer :: w
        call curve%new(4, kv, dpf, 20._dp, 4._dp, -20._dp, beta_line, 0.05_dp, gamma_line)
        do w = 1, 4
            call curve%set_knot(w, 2*w-1, 2*w, real(2*w-1,dp), KNOT_C(w), KNOT_B(w), w /= 3 .or. used3)
        enddo
    end subroutine make_curve

    !> between windows ln c and B are linear in dose; before the first window the first segment continues
    subroutine test_fractions_between_windows()
        type(dose_curve)      :: curve
        real(dp), allocatable :: c(:), b(:)
        write(*,'(A)') 'test_fractions_between_windows'
        call make_curve(curve, 300._dp, 1._dp, BETA, GAMMA, .true.)
        call curve%fraction_weights(7, 1._dp, 300._dp, c, b)
        ! fraction centres 0.5, 1.5, ..., 6.5
        call assert_true(abs(c(1) - 0.9_dp**(-0.25_dp)) < TOL .and. abs(b(1) + 2.5_dp) < TOL, &
            &'before the first window the first segment continues')
        call assert_true(abs(c(2) - 0.9_dp**0.25_dp) < TOL .and. abs(b(2) - 2.5_dp) < TOL, &
            &'a quarter of the way to window 2')
        call assert_true(abs(c(4) - exp(0.75_dp*log(0.9_dp) + 0.25_dp*log(0.8_dp))) < TOL .and. &
            &abs(b(4) - 15._dp) < TOL, 'linear in ln c and B between windows 2 and 3')
        call assert_true(abs(c(7) - exp(0.25_dp*log(0.8_dp) + 0.75_dp*log(0.6_dp))) < TOL .and. &
            &abs(b(7) - 52.5_dp) < TOL, 'linear in ln c and B between windows 3 and 4')
        call curve%kill
    end subroutine test_fractions_between_windows

    !> past the last window the line's slopes continue it, clamped at zero so no fraction gains
    subroutine test_fractions_past_the_windows()
        type(dose_curve)      :: curve
        real(dp), allocatable :: c(:), b(:)
        write(*,'(A)') 'test_fractions_past_the_windows'
        call make_curve(curve, 300._dp, 1._dp, BETA, GAMMA, .true.)
        call curve%fraction_weights(10, 1._dp, 300._dp, c, b)
        call assert_true(abs(c(8) - 0.6_dp*exp(-GAMMA*0.5_dp)) < TOL .and. abs(b(8) - (60._dp + BETA*0.5_dp)) < TOL, &
            &'half a fraction past the last window')
        call assert_true(abs(c(10) - 0.6_dp*exp(-GAMMA*2.5_dp)) < TOL .and. abs(b(10) - (60._dp + BETA*2.5_dp)) < TOL, &
            &'the line slopes continue past the last window')
        call make_curve(curve, 300._dp, 1._dp, -5._dp, -0.1_dp, .true.)
        call curve%fraction_weights(10, 1._dp, 300._dp, c, b)
        call assert_true(abs(c(10) - 0.6_dp) < TOL .and. abs(b(10) - 60._dp) < TOL, &
            &'a rising line holds the last window instead')
        call curve%kill
    end subroutine test_fractions_past_the_windows

    !> a window whose fit was rejected is no knot: its neighbours bracket the fractions it held
    subroutine test_rejected_window_is_skipped()
        type(dose_curve)      :: curve
        real(dp), allocatable :: c(:), b(:)
        write(*,'(A)') 'test_rejected_window_is_skipped'
        call make_curve(curve, 300._dp, 1._dp, BETA, GAMMA, .false.)
        call curve%fraction_weights(7, 1._dp, 300._dp, c, b)
        ! fraction 5 sits at 4.5, 0.375 of the way from window 2 at 3 to window 4 at 7
        call assert_true(abs(c(5) - exp(0.625_dp*log(0.9_dp) + 0.375_dp*log(0.6_dp))) < TOL .and. &
            &abs(b(5) - 28.75_dp) < TOL, 'a rejected window is bridged by its neighbours')
        call curve%kill
    end subroutine test_rejected_window_is_skipped

    !> fractions are placed in the scaled dose: at 200 kV a physical 0.8 e/A^2 is a unit scaled dose
    subroutine test_voltage_scaling()
        type(dose_curve)      :: curve
        real(dp), allocatable :: c(:), b(:)
        write(*,'(A)') 'test_voltage_scaling'
        call assert_true(abs(scaled_dose(1._dp, 300._dp) - 1._dp) < TOL, 'no scaling at 300 kV')
        call assert_true(abs(scaled_dose(1._dp, 200._dp) - 1._dp/0.8_dp) < TOL, 'dose divided by 0.8 at 200 kV')
        call assert_true(abs(scaled_dose(1._dp, 100._dp) - 1._dp/0.64_dp) < TOL, 'dose divided by 0.64 at 100 kV')
        call make_curve(curve, 200._dp, 0.8_dp, BETA, GAMMA, .true.)
        call curve%fraction_weights(4, 0.8_dp, 200._dp, c, b)
        call assert_true(abs(c(2) - 0.9_dp**0.25_dp) < TOL .and. abs(b(2) - 2.5_dp) < TOL, &
            &'a 200 kV movie meets the curve in the scaled dose')
        call curve%kill
    end subroutine test_voltage_scaling

    !> what is written is read back to the last bit
    subroutine test_file_round_trip()
        type(dose_curve)      :: curve, readback
        type(string)          :: message
        real(dp), allocatable :: c(:), b(:), c2(:), b2(:)
        logical :: ok
        write(*,'(A)') 'test_file_round_trip'
        call make_curve(curve, 300._dp, 1._dp, BETA, GAMMA, .false.)
        call curve%write(string(FNAME))
        ok = readback%read(string(FNAME), message)
        call assert_true(ok, 'a written dose curve reads back')
        if( ok )then
            call curve%fraction_weights(12, 1._dp, 300._dp, c, b)
            call readback%fraction_weights(12, 1._dp, 300._dp, c2, b2)
            call assert_true(all(c == c2) .and. all(b == b2), 'the read curve gives the written fraction weights')
        endif
        call curve%kill
        call readback%kill
        call del_file(FNAME)
    end subroutine test_file_round_trip

    !> one voltage and one dose rate: the producer's rule
    subroutine test_movie_checks()
        use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan, ieee_positive_inf
        type(dose_curve) :: curve
        type(string)     :: message
        real(dp) :: nan, inf
        nan = ieee_value(1._dp, ieee_quiet_nan)
        inf = ieee_value(1._dp, ieee_positive_inf)
        write(*,'(A)') 'test_movie_checks'
        call make_curve(curve, 300._dp, 1._dp, BETA, GAMMA, .true.)
        call assert_true(curve%check_movie(300._dp, 1._dp, message), 'a movie at the curve voltage and dose rate')
        call assert_true(curve%check_movie(300._dp, 1._dp + 1.e-7_dp, message), 'a single-precision dose rate matches')
        call assert_false(curve%check_movie(200._dp, 1._dp, message), 'another voltage is refused')
        call assert_true(message%strlen_trim() > 0, 'a refusal says why')
        call assert_false(curve%check_movie(300._dp, 1.01_dp, message), 'another dose rate is refused')
        call assert_false(curve%check_movie(nan, 1._dp, message), 'NaN movie voltage is refused')
        call assert_false(curve%check_movie(inf, 1._dp, message), 'infinite movie voltage is refused')
        call assert_false(curve%check_movie(300._dp, nan, message), 'NaN movie dose is refused')
        call assert_false(curve%check_movie(300._dp, inf, message), 'infinite movie dose is refused')
        call assert_false(curve%check_movie(300._dp, 0._dp, message), 'zero movie dose is refused')
        call assert_false(curve%check_movie(300._dp, -1._dp, message), 'negative movie dose is refused')
        call curve%kill
    end subroutine test_movie_checks

    !> each defect of the file is reported, not read
    subroutine test_bad_files_are_rejected()
        integer, parameter :: NL = 12
        character(len=48) :: valid(NL), lines(NL)
        type(dose_curve)  :: curve
        type(string)      :: message
        write(*,'(A)') 'test_bad_files_are_rejected'
        valid = [character(len=48) :: '# SIMPLE dose curve 1', '# program dose_filter3D', '# kv 3.0E+002', &
            &'# dose_per_fraction 1.0E+000', '# band 2.0E+001 4.0E+000', '# line -2.0E+001 1.0E+001', &
            &'# contrast 5.0E-002 4.0E-002', '# nknots 3', '# columns window from to dose c b used', &
            &'1 1 2 1.0 1.0 0.0 1', '2 3 4 3.0 0.9 10.0 1', '3 5 6 5.0 0.8 30.0 1']
        call write_lines(valid, NL)
        call assert_true(curve%read(string(FNAME), message), 'the valid template reads')
        lines = valid; lines(1) = '# SIMPLE dose curve 2'
        call expect_refused(lines, NL, 'an unknown version is refused')
        lines = valid; lines(7) = '# comment'
        call expect_refused(lines, NL, 'a missing header is refused')
        call expect_refused(valid, NL-1, 'fewer windows than nknots are refused')
        lines = valid; lines(10) = '1 1 2 1.0 0.9 0.0 1'
        call expect_refused(lines, NL, 'a window 1 off the gauge is refused')
        lines = valid; lines(12) = '3 5 6 3.0 0.8 30.0 1'
        call expect_refused(lines, NL, 'doses that do not increase are refused')
        lines = valid; lines(11) = '2 3 4 3.0 0.9 10.0 0'; lines(12) = '3 5 6 5.0 0.8 30.0 0'
        call expect_refused(lines, NL, 'a single fitted window is refused')
        lines = valid; lines(11) = '2 3 4 3.0 0.9 10.0'
        call expect_refused(lines, NL, 'a short window row is refused')
        lines = valid; lines(3) = '# kv 2.5E+002'
        call expect_refused(lines, NL, 'an unknown voltage is refused')
        lines = valid; lines(3) = '# kv NaN'
        call expect_refused(lines, NL, 'NaN curve voltage is refused')
        lines = valid; lines(3) = '# kv Infinity'
        call expect_refused(lines, NL, 'infinite curve voltage is refused')
        lines = valid; lines(4) = '# dose_per_fraction NaN'
        call expect_refused(lines, NL, 'NaN curve dose is refused')
        lines = valid; lines(4) = '# dose_per_fraction 1.0,2.0'
        call expect_refused(lines, NL, 'a list-directed separator is refused')
        call curve%kill
        call del_file(FNAME)

      contains

        subroutine expect_refused( text, n, msg )
            character(len=*), intent(in) :: text(:), msg
            integer,          intent(in) :: n
            logical :: ok
            call write_lines(text, n)
            ok = curve%read(string(FNAME), message)
            call assert_false(ok, msg)
            call assert_true(message%strlen_trim() > 0, msg//', with a reason')
        end subroutine expect_refused

    end subroutine test_bad_files_are_rejected

    subroutine write_lines( text, n )
        character(len=*), intent(in) :: text(:)
        integer,          intent(in) :: n
        integer :: unit, i
        open(newunit=unit, file=FNAME, status='replace', action='write')
        do i = 1, n
            write(unit,'(A)') trim(text(i))
        enddo
        close(unit)
    end subroutine write_lines

end module simple_dose_curve_tester
