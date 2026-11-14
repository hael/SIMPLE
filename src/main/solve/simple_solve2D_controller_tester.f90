!@descr: unit tests of the solve2D stage rules: the class FRC=0.5 particle quantile and the stage omission
! In-memory class resolutions and populations; no project file or job.
submodule (simple_solve2D_controller) simple_solve2D_controller_tester
use simple_test_utils
implicit none

real, parameter :: TOL = 1.e-5

contains

    module subroutine run_all_solve2D_stage_rule_tests()
        write(*,'(A)') '**** running all solve2D stage rule tests ****'
        call test_quantile_median()
        call test_quantile_junk_classes()
        call test_quantile_empty()
        call test_omission_rule()
        call test_omission_never()
        call test_sigma2_annealing()
    end subroutine run_all_solve2D_stage_rule_tests

    subroutine test_sigma2_annealing()
        use simple_cmdline,       only: cmdline
        use simple_sigma2_anchor, only: sigma2_stage_weight
        class(parameters), allocatable :: params
        class(cmdline), allocatable :: cline, child
        type(stage_params), allocatable :: stages(:)
        real :: first, last, hist(2)
        integer :: nstages, stage, it
        write(*,'(A)') 'test_sigma2_annealing'
        allocate(params,cline,child)
        params%sigma2_anneal='yes'
        params%refine='prob'
        params%ext='mrc'
        params%nits_per_stage=5
        params%cc_objfun=OBJFUN_EUCLID
        call cline%set('prob_beta',0.3)
        call cline%set('prob_beta_end',0.)
        do nstages=5,6
            allocate(stages(nstages))
            stages(:)%lp=8.
            stages(:)%l_lpset=.true.
            do stage=1,nstages
                call child%set('endit',5*(stage-1))
                call child%set('prob_beta',0.3)
                call child%set('prob_beta_end',0.)
                call set_cline_refine2D_stage(child,cline,params,stages,20,stage)
                first=0.; last=0.
                if( stage==nstages-1 ) last=0.5
                if( stage==nstages )then
                    first=0.5; last=1.
                endif
                call assert_real(first,child%get_rarg('sigma2_weight_start'),TOL,'stage starts at its assigned sigma weight')
                call assert_real(last,child%get_rarg('sigma2_weight_end'),TOL,'stage ends at its assigned sigma weight')
                call assert_real(0.3,child%get_rarg('prob_beta'),TOL,'sigma scheduling preserves beta')
                call assert_real(0.,child%get_rarg('prob_beta_end'),TOL,'sigma scheduling preserves explicit beta endpoint')
                do it=1,5
                    call assert_real(first+real(it-1)*(last-first)/4., &
                        &sigma2_stage_weight(it,5,first,last),TOL,'five-step ramp spans the stage')
                enddo
                call assert_real(last,sigma2_stage_weight(1,1,first,last),TOL,'single iteration reaches endpoint')
                call assert_real(last,sigma2_stage_weight(8,5,first,last),TOL,'coverage extension holds endpoint')
            enddo
            hist=[5.,5.]
            call assert_false(solve2D_stage_behind_frc05(stages,nstages-1,1,nstages,hist,2,sigma2_anneal=.true.), &
                &'penultimate annealing stage cannot be omitted')
            call assert_true(solve2D_stage_behind_frc05(stages,nstages-1,1,nstages,hist,2,sigma2_anneal=.false.), &
                &'ordinary stage omission remains enabled')
            call set_cline_refine2D_seed_pass(child,params,stages,10,'seed.mrc')
            call assert_real(0.,child%get_rarg('sigma2_weight_end'),TOL,'seed pass keeps the bootstrap')
            deallocate(stages)
        enddo
        call assert_real(1.,sigma2_stage_weight(1,1,1.,1.),TOL,'terminal pass is measured')
        params%sigma2_anneal='no'
        allocate(stages(6))
        call child%kill()
        call set_cline_refine2D_stage(child,cline,params,stages,20,1)
        call assert_false(child%defined('sigma2_weight_start'),'off path adds no sigma schedule')
        call child%kill()
        call cline%kill()
    end subroutine test_sigma2_annealing

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

end submodule simple_solve2D_controller_tester
