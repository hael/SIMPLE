!@descr: unit tests of the solve2D stage rules: the class FRC=0.5 particle quantile and the stage omission
! In-memory class resolutions and populations; no project file or job.
module simple_solve2D_controller_tester
use simple_solve2D_controller, only: stage_params, frc05_ptcl_quantile, solve2D_stage_behind_frc05, &
    &SKIP_MIN_STAGE, FRC05_SKIP_PTCL_FRAC
use simple_test_utils
implicit none
private
public :: run_all_solve2D_stage_rule_tests

real, parameter :: TOL = 1.e-5

contains

    subroutine run_all_solve2D_stage_rule_tests()
        write(*,'(A)') '**** running all solve2D stage rule tests ****'
        call test_quantile_median()
        call test_quantile_junk_classes()
        call test_quantile_empty()
        call test_omission_rule()
        call test_omission_never()
    end subroutine run_all_solve2D_stage_rule_tests

    !> the lower population-weighted quantile and the fraction past a limit
    subroutine test_quantile_median()
        real    :: res05_q, frac_past
        write(*,'(A)') 'test_quantile_median'
        ! the dominant class decides the median particle
        call frc05_ptcl_quantile([6.,8.,10.,12.], [10,10,10,70], 0.5, res05_q, lp=11., frac_past=frac_past)
        call assert_real(12.,  res05_q,   TOL, 'median particle in the 12 A class')
        call assert_real(0.3,  frac_past, TOL, '30% of the particles past 11 A')
        ! the best classes hold the median particle
        call frc05_ptcl_quantile([6.,8.,10.,12.], [40,20,20,20], 0.5, res05_q, lp=9., frac_past=frac_past)
        call assert_real(8.,   res05_q,   TOL, 'median particle in the 8 A class')
        call assert_real(0.6,  frac_past, TOL, '60% of the particles past 9 A')
        ! unsorted input, a lower fraction
        call frc05_ptcl_quantile([12.,6.,10.,8.], [20,40,20,20], 0.25, res05_q)
        call assert_real(6.,   res05_q,   TOL, 'quarter: the best class')
        call frc05_ptcl_quantile([12.,6.,10.,8.], [20,40,20,20], 1.0, res05_q)
        call assert_real(12.,  res05_q,   TOL, 'all: the worst populated class')
    end subroutine test_quantile_median

    !> a few overfit (Nyquist) or empty classes do not move the median particle's resolution
    subroutine test_quantile_junk_classes()
        real :: res05_q, frac_past
        write(*,'(A)') 'test_quantile_junk_classes'
        call frc05_ptcl_quantile([2.5,2.5,2.5,10.,10.,30.], [1,1,1,500,500,0], FRC05_SKIP_PTCL_FRAC, res05_q, &
            &lp=9., frac_past=frac_past)
        call assert_real(10., res05_q, TOL, 'three singleton Nyquist classes do not move the median')
        call assert_real(3./1003., frac_past, TOL, 'only the singletons are past 9 A')
    end subroutine test_quantile_junk_classes

    !> no populated class, no class, mismatched input: 0
    subroutine test_quantile_empty()
        real    :: res05_q, frac_past
        integer :: none(0)
        real    :: rnone(0)
        write(*,'(A)') 'test_quantile_empty'
        call frc05_ptcl_quantile([6.,8.], [0,0], 0.5, res05_q, lp=7., frac_past=frac_past)
        call assert_real(0., res05_q,   TOL, 'no populated class: 0')
        call assert_real(0., frac_past, TOL, 'no populated class: no fraction')
        call frc05_ptcl_quantile(rnone, none, 0.5, res05_q)
        call assert_real(0., res05_q,   TOL, 'no class: 0')
        call frc05_ptcl_quantile([6.,8.], [10], 0.5, res05_q)
        call assert_real(0., res05_q,   TOL, 'mismatched sizes: 0')
    end subroutine test_quantile_empty

    !> the six-stage ladder 20, 14, 11, 9.5, 8, gold standard
    subroutine make_ladder( stage_parms )
        type(stage_params), allocatable, intent(out) :: stage_parms(:)
        allocate(stage_parms(6))
        stage_parms(:)%lp      = [20., 14., 11., 9.5, 8., 8.]
        stage_parms(:)%l_lpset = [.true., .true., .true., .true., .true., .false.]
    end subroutine make_ladder

    !> omitted when the resolution at the ends of the last two executed stages is finer than the planned limit
    subroutine test_omission_rule()
        type(stage_params), allocatable :: stage_parms(:)
        real :: hist(6)
        write(*,'(A)') 'test_omission_rule'
        call make_ladder(stage_parms)
        hist = 0.
        ! stages 1 and 2 ended at 12 and 9 A: stage 3 (11 A) is not behind both
        hist(1:2) = [12., 9.]
        call assert_false(solve2D_stage_behind_frc05(stage_parms, 3, 1, 6, hist, 2), 'held at one end only: run')
        ! 9.4 and 9 A: stages 3 (11 A) and 4 (9.5 A) behind, stage 5 (8 A) not
        hist(1:2) = [9.4, 9.]
        call assert_true( solve2D_stage_behind_frc05(stage_parms, 3, 1, 6, hist, 2), 'stage 3 behind 9.4 and 9 A')
        call assert_true( solve2D_stage_behind_frc05(stage_parms, 4, 1, 6, hist, 2), 'stage 4 behind 9.4 and 9 A')
        call assert_false(solve2D_stage_behind_frc05(stage_parms, 5, 1, 6, hist, 2), 'stage 5 at 8 A not behind')
        ! the last two executed ends decide, not the earlier ones
        hist(1:3) = [30., 7.5, 7.]
        call assert_true( solve2D_stage_behind_frc05(stage_parms, 5, 1, 6, hist, 3), 'stage 5 behind the last two ends')
        hist(1:3) = [7., 7.5, 30.]
        call assert_false(solve2D_stage_behind_frc05(stage_parms, 5, 1, 6, hist, 3), 'a worse latest end: run')
        ! exactly at the limit is not behind
        hist(1:2) = [11., 11.]
        call assert_false(solve2D_stage_behind_frc05(stage_parms, 3, 1, 6, hist, 2), 'at the limit: run')
    end subroutine test_omission_rule

    !> never the entry stage, the first SKIP_MIN_STAGE stages, the last stage, a stage without a planned limit,
    !! or on fewer than two recorded ends or a missing estimate
    subroutine test_omission_never()
        type(stage_params), allocatable :: stage_parms(:)
        real :: hist(6)
        write(*,'(A)') 'test_omission_never'
        call make_ladder(stage_parms)
        hist = 0.
        hist(1:3) = [7., 7., 7.]
        call assert_int(2, SKIP_MIN_STAGE, 'the first two stages are never omitted')
        call assert_false(solve2D_stage_behind_frc05(stage_parms, 1, 1, 6, hist, 3), 'never stage 1')
        call assert_false(solve2D_stage_behind_frc05(stage_parms, 2, 1, 6, hist, 3), 'never stage 2')
        call assert_false(solve2D_stage_behind_frc05(stage_parms, 6, 1, 6, hist, 3), 'never the last stage')
        call assert_false(solve2D_stage_behind_frc05(stage_parms, 5, 1, 5, hist, 3), 'never the last stage of the invocation')
        call assert_false(solve2D_stage_behind_frc05(stage_parms, 3, 3, 6, hist, 3), 'never the entry stage')
        call assert_true( solve2D_stage_behind_frc05(stage_parms, 4, 3, 6, hist, 3), 'the stage after the entry stage may be')
        call assert_false(solve2D_stage_behind_frc05(stage_parms, 3, 1, 6, hist, 1), 'one recorded end: run')
        call assert_false(solve2D_stage_behind_frc05(stage_parms, 3, 1, 6, hist, 0), 'no recorded end: run')
        hist(1:2) = [0., 7.]
        call assert_false(solve2D_stage_behind_frc05(stage_parms, 3, 1, 6, hist, 2), 'a missing estimate: run')
        hist(1:2) = [7., 7.]
        stage_parms(3)%l_lpset = .false.
        call assert_false(solve2D_stage_behind_frc05(stage_parms, 3, 1, 6, hist, 2), 'no planned limit: run')
    end subroutine test_omission_never

end module simple_solve2D_controller_tester
