!@descr: unit tests for the 3D convergence report and the refine=cont convergence rule (simple_convergence)
! N19: a Cartesian pass reports the mean corr_cart, attempts and improved percentage over its sample;
! a polar pass reports corr. N25: a refine=cont pass converges when 90% of every state moved at most
! 0.5 deg and 1 A with an accepted or no-improvement outcome, over a streak covering 90% of the active
! particles (the cont_conv_rot, cont_conv_shift and cont_conv_frac defaults; the parameters replace
! them); a discrete pass, with or without the polish, keeps its own rule. The multi-state rule
! converges on the state overlap and the scan together; with scan_exit=yes (solve3D's multi-state
! shc stage 2, fracsrch=95) the scan alone also ends it, and without the key (stage 1, fracsrch=99,
! and every other user) a low overlap holds it. Expected values are hand-counted from the fixtures.
module simple_convergence_tester
use simple_core_module_api, only: oris, del_file
use simple_cmdline,         only: cmdline
use simple_convergence,     only: convergence, cont_particle_stable
use simple_parameters,      only: parameters, PROB_BETA_TAIL
use simple_defs_fname,      only: STATS_FILE
use simple_defs_conv,       only: CONT_CONV_ROT_DEG, CONT_CONV_SHIFT_A, CONT_CONV_FRAC
use simple_cartft_pose_opt, only: CARTFT_NOT_ATTEMPTED, CARTFT_INVALID_PREPARATION, CARTFT_ACCEPTED, &
    &CARTFT_NO_IMPROVEMENT, CARTFT_NO_RELIABLE_UPDATE, CARTFT_BOUND_REJECTED, CARTFT_INVALID_NUMERICS, &
    &CARTFT_ITERATION_LIMIT
use simple_test_utils
implicit none
private
public :: run_all_convergence_tests

integer, parameter :: NPTCLS = 6
integer, parameter :: NCONT  = 10   !< particles of the refine=cont rule fixture (90% is nine)
real,    parameter :: TOL    = 1.e-5

contains

    subroutine run_all_convergence_tests()
        write(*,'(A)') '**** running all convergence tests ****'
        write(*,'(A)') 'test_cartesian_pass_scores'
        call test_cartesian_pass_scores()
        write(*,'(A)') 'test_polar_pass_scores'
        call test_polar_pass_scores()
        write(*,'(A)') 'test_cartesian_pass_statistics'
        call test_cartesian_pass_statistics()
        write(*,'(A)') 'test_cont_particle_stable'
        call test_cont_particle_stable()
        write(*,'(A)') 'test_cont_rule_full_sample'
        call test_cont_rule_full_sample()
        write(*,'(A)') 'test_cont_rule_streak'
        call test_cont_rule_streak()
        write(*,'(A)') 'test_cont_rule_tunable'
        call test_cont_rule_tunable()
        write(*,'(A)') 'test_multistate_scan_exit'
        call test_multistate_scan_exit()
        write(*,'(A)') 'test_lpexit_gain_guard'
        call test_lpexit_gain_guard()
        write(*,'(A)') 'test_posterior_assignment_bypass'
        call test_posterior_assignment_bypass()
        write(*,'(A)') 'test_posterior_anneal_tail'
        call test_posterior_anneal_tail()
        call del_file(STATS_FILE)
    end subroutine run_all_convergence_tests

    !> the annealing tail (pose_marginalization.md 19): a refine=prob run with prob_beta_end above prob_beta
    !! searches flat while the resolution gains (FSC 0.5 or 0.143 finer than the previous iteration), starts the
    !! tail at the first iteration where neither gained (tail_start = the next iteration), and converges at the
    !! tail's last iteration, PROB_BETA_TAIL after its start; a run without prob_beta_end never starts a tail
    subroutine test_posterior_anneal_tail()
        class(parameters), allocatable :: params
        type(convergence) :: conv
        type(cmdline)     :: cline
        type(oris)        :: os
        real    :: res05(8), res143(8)
        logical :: conv_at(8)
        integer :: it
        call del_file(STATS_FILE)
        res05  = [7.5, 7.0, 7.0, 6.9, 6.9, 6.9, 6.9, 6.9]  ! gain at 2, flat at 3 (stall), a late gain inside the tail
        res143 = [6.5, 6.0, 6.0, 6.0, 6.0, 6.0, 6.0, 6.0]
        allocate(params)
        params%nstates          = 1
        params%trs              = 5.
        params%refine           = 'prob'
        params%l_prob_align_mode = .true.
        params%prob_beta        = 0.3
        params%prob_beta_end    = 1.
        params%startit          = 1
        params%maxits           = 25
        params%prob_beta_tail_start = 0
        call cline%set('trs',      5.)
        call cline%set('overlap',  0.99)
        call cline%set('fracsrch', 99.)
        do it = 1, 8
            params%which_iter = it
            call frozen_field(os, res05(it), res143(it))
            conv_at(it) = conv%check_conv3D(params, cline, os, 40.)
            call os%kill
            if( it == 2 ) call assert_int(0, params%prob_beta_tail_start, 'a gaining iteration starts no tail')
            if( it == 3 ) call assert_int(4, params%prob_beta_tail_start, 'the first iteration without gain starts the tail at the next')
        enddo
        call assert_true(.not. any(conv_at(1:3)), 'the search phase does not converge')
        call assert_true(.not. any(conv_at(4:4+PROB_BETA_TAIL-2)), 'the tail does not converge before its last iteration')
        call assert_true(conv_at(4+PROB_BETA_TAIL-1), 'the tail converges at its last iteration')
        call assert_int(4, params%prob_beta_tail_start, 'a gain inside the tail does not restart it')
        ! no prob_beta_end: no tail
        call del_file(STATS_FILE)
        params%prob_beta_end        = params%prob_beta
        params%prob_beta_tail_start = 0
        do it = 1, 4
            params%which_iter = it
            call frozen_field(os, 7.0, 6.0)
            conv_at(it) = conv%check_conv3D(params, cline, os, 40.)
            call os%kill
        enddo
        call assert_int(0, params%prob_beta_tail_start, 'a flat run (prob_beta_end = prob_beta) starts no tail')
        call assert_true(.not. any(conv_at(1:4)), 'a flat run does not converge on the stall')
        call cline%kill
        deallocate(params)
        call del_file(STATS_FILE)

    contains

        !> NPTCLS frozen particles (full scan, overlap 1) at the given FSC=0.5 and FSC=0.143 resolutions
        subroutine frozen_field( os_loc, r05, r143 )
            type(oris), intent(inout) :: os_loc
            real,       intent(in)    :: r05, r143
            integer :: i
            call os_loc%new(NPTCLS, is_ptcl=.true.)
            do i = 1, NPTCLS
                call os_loc%set_state(i, 1)
                call os_loc%set(i, 'updatecnt', 1.)
                call os_loc%set(i, 'sampled',   1.)
                call os_loc%set(i, 'frac',      100.)
                call os_loc%set(i, 'mi_proj',   1.)
                call os_loc%set(i, 'corr',      0.5)
                call os_loc%set(i, 'res05',     r05)
                call os_loc%set(i, 'res',       r143)
            end do
        end subroutine frozen_field

    end subroutine test_posterior_anneal_tail

    !> the probability tables' posterior draws (pose_marginalization.md 10.3): a frozen assignment (full scan,
    !! overlap 1) converges the regular rule on a shc pass and must not converge on the dense table (refine=prob)
    !! nor on the neighbourhood table (refine=prob_neigh, every neighbourhood mode), whose committed pose is a
    !! draw; the lpexit rule still ends the stage (a stalled 7.5 A past lpexit 8 at the stage's second iteration)
    subroutine test_posterior_assignment_bypass()
        call del_file(STATS_FILE)
        call assert_true(frozen_converged('shc', 0., 1, 1), 'the regular rule did not converge on a frozen shc assignment')
        call assert_false(frozen_converged('prob', 0., 1, 1), 'the pose-overlap rule ended a dense-table stage')
        call assert_false(frozen_converged('prob', 0., 2, 1), 'the pose-overlap rule ended a dense-table stage (second iteration)')
        call assert_false(frozen_converged('prob_neigh', 0., 3, 1, neigh_mode='shc'), 'the pose-overlap rule ended a neighbourhood stage (shc)')
        call assert_false(frozen_converged('prob_neigh', 0., 3, 1, neigh_mode='state'), 'the pose-overlap rule ended a neighbourhood stage (state)')
        call assert_true(frozen_converged('shc', 0., 3, 1), 'a shc pass lost its pose-overlap rule')
        call del_file(STATS_FILE)
        call assert_false(frozen_converged('prob', 8., 1, 1), 'the stage ended at its first iteration')
        call assert_true(frozen_converged('prob', 8., 2, 1),  'the resolution exit did not end the stage')
        call del_file(STATS_FILE)
    end subroutine test_posterior_assignment_bypass

    !> check_conv3D of a single-state pass of NPTCLS frozen particles (frac 100, mi_proj 1) at FSC=0.5
    !! 7.5 A, iteration which_iter of a stage from startit, with the refine mode, lpexit and, for prob_neigh, the
    !! neighbourhood mode as given
    logical function frozen_converged( refine, lpexit, which_iter, startit, neigh_mode ) result( converged )
        character(len=*),           intent(in) :: refine
        real,                       intent(in) :: lpexit
        integer,                    intent(in) :: which_iter, startit
        character(len=*), optional, intent(in) :: neigh_mode
        class(parameters), allocatable :: params
        type(convergence) :: conv
        type(cmdline)     :: cline
        type(oris)        :: os
        integer :: i
        call os%new(NPTCLS, is_ptcl=.true.)
        do i = 1, NPTCLS
            call os%set_state(i, 1)
            call os%set(i, 'updatecnt', 1.)
            call os%set(i, 'sampled',   1.)
            call os%set(i, 'frac',      100.)
            call os%set(i, 'mi_proj',   1.)
            call os%set(i, 'corr',      0.5)
            call os%set(i, 'res05',     7.5)
        end do
        allocate(params)
        params%nstates    = 1
        params%trs        = 5.
        params%refine     = refine
        params%lpexit     = lpexit
        params%startit    = startit
        params%which_iter = which_iter
        if( present(neigh_mode) ) params%prob_neigh_mode = neigh_mode
        call cline%set('trs',      5.)
        call cline%set('overlap',  0.99)
        call cline%set('fracsrch', 99.)
        converged = conv%check_conv3D(params, cline, os, 40.)
        call cline%kill
        call os%kill
        deallocate(params)
    end function frozen_converged

    !> Six particles: corr 0.1 i, corr_cart 0.5 + 0.05 i, improved on the odd ones; the
    !! sampled mask holds particles 1, 2, 3 and 6.
    subroutine make_field( os, mask )
        type(oris), intent(inout) :: os
        logical,    intent(out)   :: mask(NPTCLS)
        integer :: i
        call os%new(NPTCLS, is_ptcl=.true.)
        do i = 1, NPTCLS
            call os%set_state(i, 1)
            call os%set(i, 'corr',               0.1*real(i))
            call os%set(i, 'corr_cart',          0.5 + 0.05*real(i))
            call os%set(i, 'pose_cont_improved', real(mod(i,2)))
        end do
        mask = [.true., .true., .true., .false., .false., .true.]
    end subroutine make_field

    subroutine test_cartesian_pass_scores()
        type(convergence) :: conv
        type(oris)        :: os
        logical :: mask(NPTCLS), cart_stats
        call make_field(os, mask)
        call conv%calc_pass_scores(os, mask, .true., cart_stats)
        call assert_true(cart_stats, 'a Cartesian pass reports pose statistics')
        ! (0.55 + 0.60 + 0.65 + 0.80) / 4
        call assert_real(0.65, conv%get('corr'), TOL, 'the Cartesian score is the mean corr_cart over the sampled')
        call assert_real(4., conv%get('pose_cont_attempts'), TOL, 'the attempts are the sampled particles')
        ! particles 1 and 3 of the four improved
        call assert_real(50., conv%get('pose_cont_improved_pct'), TOL, 'the improved percentage counts the flag')
        call os%kill
    end subroutine test_cartesian_pass_scores

    subroutine test_polar_pass_scores()
        type(convergence) :: conv
        type(oris)        :: os
        logical :: mask(NPTCLS), cart_stats
        call make_field(os, mask)
        call conv%calc_pass_scores(os, mask, .true., cart_stats)
        call conv%calc_pass_scores(os, mask, .false., cart_stats)
        call assert_false(cart_stats, 'a polar pass reports no pose statistics')
        ! (0.1 + 0.2 + 0.3 + 0.6) / 4
        call assert_real(0.3, conv%get('corr'), TOL, 'the polar score is the mean corr')
        call assert_real(0., conv%get('pose_cont_attempts'), TOL, 'a polar pass clears the attempts')
        call assert_real(0., conv%get('pose_cont_improved_pct'), TOL, 'a polar pass clears the improved percentage')
        call os%kill
    end subroutine test_polar_pass_scores

    !> with no transaction status recorded a Cartesian pass does not converge, and it still reports
    !! its motion statistics; the same field under a discrete mode, with or without the polish,
    !! converges by the discrete rule
    subroutine test_cartesian_pass_statistics()
        type(convergence) :: conv
        type(parameters)  :: params
        type(cmdline)     :: cline
        type(oris)        :: os
        integer :: i, ilim
        real    :: limits(2,3)
        logical :: converged
        ! (overlap, fracsrch): the defaults, permissive and impossible-to-miss limits
        limits = reshape([0.99, 99., 0.5, 50., 0., 0.], [2,3])
        call os%new(NPTCLS, is_ptcl=.true.)
        do i = 1, NPTCLS
            call os%set_state(i, 1)
            call os%set(i, 'updatecnt',          1.)
            call os%set(i, 'sampled',            1.)
            call os%set(i, 'frac',               100.)
            call os%set(i, 'mi_proj',            1.)
            call os%set(i, 'dist',               0.1*real(i))
            call os%set(i, 'shincarg',           0.01*real(i))
            call os%set(i, 'corr',               0.5)
            call os%set(i, 'corr_cart',          0.4)
            call os%set(i, 'pose_cont_improved', real(mod(i,2)))
        end do
        params%nstates  = 1
        params%trs      = 5.
        call cline%set('trs', 5.)
        do ilim = 1, size(limits,2)
            call cline%set('overlap',  limits(1,ilim))
            call cline%set('fracsrch', limits(2,ilim))
            params%refine        = 'cont'
            params%l_cart_refine = .true.
            converged = conv%check_conv3D(params, cline, os, 40.)
            call assert_false(converged, 'a Cartesian pass without a recorded status declared convergence')
            ! the motion statistics are computed all the same: (0.1+...+0.6)/6, (0.01+...+0.06)/6
            call assert_real(0.35, conv%get('dist'), TOL, 'the Cartesian pass mean orientation motion')
            call assert_real(50., conv%get('pose_cont_improved_pct'), TOL, 'the Cartesian pass improved fraction')
            call assert_real(real(NPTCLS), conv%get('pose_cont_attempts'), TOL, 'the Cartesian pass attempts')
            ! the same field under a discrete mode: the discrete rule is unchanged
            params%refine        = 'neigh'
            params%l_cart_refine = .false.
            converged = conv%check_conv3D(params, cline, os, 40.)
            call assert_true(converged, 'the discrete rule changed')
            ! the same discrete pass followed by the polish
            params%pose_cont = 'yes'
            converged = conv%check_conv3D(params, cline, os, 40.)
            call assert_true(converged, 'the polish changed the discrete rule')
            call assert_real(0.35, conv%get('dist'), TOL, 'the discrete motion with the polish')
            call assert_real(real(NPTCLS), conv%get('pose_cont_attempts'), TOL, 'the polish attempts')
            call assert_real(50., conv%get('pose_cont_improved_pct'), TOL, 'the polish improved fraction')
            params%pose_cont = 'no'
        end do
        call del_file(STATS_FILE)
        call cline%kill
        call os%kill
    end subroutine test_cartesian_pass_statistics

    !> three states of two particles, all updated: a mean scan of 96% with a joint overlap of 0.5 ends
    !! solve3D's multi-state stage 2 (scan_exit=yes, fracsrch=95) and not stage 1 (no key, fracsrch=99),
    !! nor the rule without the key at fracsrch=95; 94% does not end stage 2; with every state's
    !! overlap above 0.96 the present rule still converges without the key
    subroutine test_multistate_scan_exit()
        type(oris) :: os
        call make_multistate_field(os, frac=96., mi_proj=0.5)
        call assert_true(multistate_converged(os, 95., .true.),   'stage 2 did not end at a 96% scan with a low overlap')
        call assert_false(multistate_converged(os, 99., .false.), 'stage 1 ended at a 96% scan with a low overlap')
        call assert_false(multistate_converged(os, 95., .false.), 'the rule without scan_exit ended on the scan alone')
        call make_multistate_field(os, frac=94., mi_proj=0.5)
        call assert_false(multistate_converged(os, 95., .true.),  'stage 2 ended at a 94% scan')
        call make_multistate_field(os, frac=96., mi_proj=0.98)
        call assert_true(multistate_converged(os, 95., .false.),  'the overlap-and-scan rule no longer converges')
        call assert_true(multistate_converged(os, 95., .true.),   'scan_exit blocked the overlap-and-scan rule')
        call del_file(STATS_FILE)
        call os%kill
    end subroutine test_multistate_scan_exit

    !> NPTCLS particles, two per state 1-3, all sampled and updated in generation 1, each scanning
    !! frac % of the search space with orientation overlap mi_proj
    subroutine make_multistate_field( os, frac, mi_proj )
        type(oris), intent(inout) :: os
        real,       intent(in)    :: frac, mi_proj
        integer :: i
        call os%kill
        call os%new(NPTCLS, is_ptcl=.true.)
        do i = 1, NPTCLS
            call os%set_state(i, (i + 1)/2)
            call os%set(i, 'updatecnt', 1.)
            call os%set(i, 'sampled',   1.)
            call os%set(i, 'frac',      frac)
            call os%set(i, 'mi_proj',   mi_proj)
            call os%set(i, 'corr',      0.5)
        end do
    end subroutine make_multistate_field

    !> check_conv3D of a three-state shc pass on os with the stage's fracsrch, with or without scan_exit
    logical function multistate_converged( os, fracsrch, l_scan_exit ) result( converged )
        type(oris), intent(inout) :: os
        real,       intent(in)    :: fracsrch
        logical,    intent(in)    :: l_scan_exit
        class(parameters), allocatable :: params
        type(convergence) :: conv
        type(cmdline)     :: cline
        allocate(params)
        params%nstates     = 3
        params%trs         = 5.
        params%refine      = 'shc'
        params%l_scan_exit = l_scan_exit
        params%startit     = 1
        params%which_iter  = 1
        call cline%set('trs',      5.)
        call cline%set('overlap',  0.99)
        call cline%set('fracsrch', fracsrch)
        converged = conv%check_conv3D(params, cline, os, 40.)
        call cline%kill
        deallocate(params)
    end function multistate_converged

    !> lpexit=8 A over consecutive iterations of a stage from iteration 1: 9 A holds it, 7.5 A still
    !! gaining holds it, a repeated 7.5 A ends it and so does a loss to 7.6 A; the first iteration of
    !! a stage from 5 does not compare with the last of the previous stage, an earlier stage's stats
    !! are not this stage's and without lpexit a stalled 7.5 A holds; with three states the best
    !! (7.5 A, not the particle average of 8.83 A) is compared with the previous best
    subroutine test_lpexit_gain_guard()
        call del_file(STATS_FILE)
        call assert_false(lpexit_converged([9.0], 1, 1, 8.), 'the stage ended short of lpexit')
        call assert_false(lpexit_converged([7.5], 2, 1, 8.), 'the stage ended while still gaining')
        call assert_true(lpexit_converged([7.5],  3, 1, 8.), 'a stalled resolution past lpexit did not end the stage')
        call assert_true(lpexit_converged([7.6],  4, 1, 8.), 'a lost resolution past lpexit did not end the stage')
        call assert_false(lpexit_converged([7.6], 5, 5, 8.), 'the first iteration compared with the previous stage')
        call assert_false(lpexit_converged([7.6], 9, 8, 8.), 'an earlier stage''s stats were compared')
        call assert_false(lpexit_converged([7.6], 10, 8, 0.), 'the exit fired without lpexit')
        call del_file(STATS_FILE)
        call assert_false(lpexit_converged([7.5, 9.0, 10.], 1, 1, 8.), 'the first three-state iteration ended the stage')
        call assert_true(lpexit_converged([7.5, 9.0, 10.],  2, 1, 8.), 'the best state was not compared with the previous best')
        call del_file(STATS_FILE)
    end subroutine test_lpexit_gain_guard

    !> check_conv3D of a shc pass at iteration which_iter of a stage from startit, NPTCLS particles
    !! split evenly over size(res05) states at FSC=0.5 resolutions res05; a 50% scan with overlap 0.5
    !! holds the regular rule, so only the lpexit rule can converge
    logical function lpexit_converged( res05, which_iter, startit, lpexit ) result( converged )
        real,    intent(in) :: res05(:), lpexit
        integer, intent(in) :: which_iter, startit
        class(parameters), allocatable :: params
        type(convergence) :: conv
        type(cmdline)     :: cline
        type(oris)        :: os
        integer :: i, istate
        call os%new(NPTCLS, is_ptcl=.true.)
        do i = 1, NPTCLS
            istate = (i - 1)*size(res05)/NPTCLS + 1
            call os%set_state(i, istate)
            call os%set(i, 'updatecnt', 1.)
            call os%set(i, 'sampled',   1.)
            call os%set(i, 'frac',      50.)
            call os%set(i, 'mi_proj',   0.5)
            call os%set(i, 'corr',      0.5)
            call os%set(i, 'res05',     res05(istate))
        end do
        allocate(params)
        params%nstates    = size(res05)
        params%trs        = 5.
        params%refine     = 'shc'
        params%lpexit     = lpexit
        params%startit    = startit
        params%which_iter = which_iter
        call cline%set('trs',      5.)
        call cline%set('overlap',  0.99)
        call cline%set('fracsrch', 99.)
        converged = conv%check_conv3D(params, cline, os, 40.)
        call cline%kill
        call os%kill
        deallocate(params)
    end function lpexit_converged

    !> the bounds are inclusive; only an accepted or a finite no-improvement outcome is stable
    subroutine test_cont_particle_stable()
        integer :: rejected(6), i
        call assert_true(cont_particle_stable(CARTFT_ACCEPTED, 0.5, 1.0, CONT_CONV_ROT_DEG, CONT_CONV_SHIFT_A), &
            &'the bounds are inclusive')
        call assert_false(cont_particle_stable(CARTFT_ACCEPTED, 0.51, 0., CONT_CONV_ROT_DEG, CONT_CONV_SHIFT_A), &
            &'a rotation above 0.5 deg is stable')
        call assert_false(cont_particle_stable(CARTFT_ACCEPTED, 0., 1.01, CONT_CONV_ROT_DEG, CONT_CONV_SHIFT_A), &
            &'a shift above 1 A is stable')
        call assert_true(cont_particle_stable(CARTFT_NO_IMPROVEMENT, 0., 0., CONT_CONV_ROT_DEG, CONT_CONV_SHIFT_A), &
            &'a finite no-improvement is not stable')
        call assert_false(cont_particle_stable(CARTFT_ACCEPTED, 0.06, 0., 0.05, 0.1), 'a rotation above a 0.05 deg bound is stable')
        call assert_false(cont_particle_stable(CARTFT_ACCEPTED, 0., 0.11, 0.05, 0.1), 'a shift above a 0.1 A bound is stable')
        rejected = [CARTFT_NOT_ATTEMPTED, CARTFT_INVALID_PREPARATION, CARTFT_NO_RELIABLE_UPDATE, &
            &CARTFT_BOUND_REJECTED, CARTFT_INVALID_NUMERICS, CARTFT_ITERATION_LIMIT]
        do i = 1, size(rejected)
            call assert_false(cont_particle_stable(rejected(i), 0., 0., CONT_CONV_ROT_DEG, CONT_CONV_SHIFT_A), &
                &'a rejected transaction at zero motion is stable')
        enddo
    end subroutine test_cont_particle_stable

    !> every particle sampled: nine of ten stable converges, eight does not, a finite no-improvement
    !! converges, an all-rejected pass at zero motion does not, the shift bound is in A, one state
    !! below 90% vetoes, and minits holds
    subroutine test_cont_rule_full_sample()
        type(convergence) :: conv
        type(oris)        :: os
        call make_cont_field(os)
        call assert_true(cont_converged(conv, os, 1., 1, .false.), 'an all-stable pass did not converge')
        call assert_real(100., conv%get('cont_stable_pct'), TOL, 'the stable percentage of an all-stable pass')
        call os%set(1, 'dist', 0.4)
        call assert_true(cont_converged(conv, os, 1., 1, .false.), 'nine of ten stable did not converge')
        call os%set(2, 'dist', 0.4)
        call assert_false(cont_converged(conv, os, 1., 1, .false.), 'eight of ten stable converged')
        call assert_real(80., conv%get('cont_stable_pct'), TOL, 'the stable percentage of eight of ten')
        ! a finite no-improvement keeps the seed: zero motion, stable
        call make_cont_field(os, CARTFT_NO_IMPROVEMENT, 0., 0.)
        call assert_true(cont_converged(conv, os, 1., 1, .false.), 'a no-improvement pass did not converge')
        ! every transaction left its bound: zero motion, the failure signature, not stable
        call make_cont_field(os, CARTFT_BOUND_REJECTED, 0., 0.)
        call assert_false(cont_converged(conv, os, 1., 1, .false.), 'an all-rejected pass converged')
        call assert_real(0., conv%get('cont_stable_pct'), TOL, 'the stable percentage of an all-rejected pass')
        ! 3 px of shift: 3 A at smpd 1, 0.9 A at smpd 0.3
        call make_cont_field(os, CARTFT_ACCEPTED, 0., 3.)
        call assert_false(cont_converged(conv, os, 1.,  1, .false.), 'a 3 A shift converged')
        call assert_true(cont_converged(conv, os, 0.3, 1, .false.), 'a 0.9 A shift did not converge')
        ! states: 1-5 all stable, 6-10 three of five (80% overall)
        call make_cont_field(os)
        call os%set_state(6, 2)
        call os%set_state(7, 2)
        call os%set_state(8, 2)
        call os%set_state(9, 2)
        call os%set_state(10, 2)
        call os%set(9,  'dist', 0.6)
        call os%set(10, 'dist', 0.6)
        call assert_false(cont_converged(conv, os, 1., 2, .false.), 'a state with 60% stable did not veto')
        call assert_real(60., conv%get('cont_stable_pct'), TOL, 'the stable percentage is the least stable state')
        ! minits: a converging pass in iteration 1 of 3
        call make_cont_field(os)
        call assert_false(cont_converged(conv, os, 1., 1, .false., minits=3), 'minits was not honoured')
        call os%kill
    end subroutine test_cont_rule_full_sample

    !> half samples: a passing first half covers 50% (no convergence), the passing second half
    !! completes the streak (convergence); a failing iteration resets the streak, after which one
    !! passing half covers 50% again
    subroutine test_cont_rule_streak()
        type(convergence) :: conv, conv_fresh
        type(oris)        :: os
        call make_cont_field(os)
        call clear_samples(os)
        call sample_half(os, 1, 1)
        call assert_false(cont_converged(conv, os, 1., 1, .true.), 'half the particles completed the streak')
        call assert_real(50., conv%get('cont_coverage_pct'), TOL, 'the coverage of the first half')
        call sample_half(os, 2, 2)
        call assert_true(cont_converged(conv, os, 1., 1, .true.), 'two passing halves did not converge')
        call assert_real(100., conv%get('cont_coverage_pct'), TOL, 'the coverage of two halves')
        ! a fresh history: first half passes, second half fails, first half passes again
        call make_cont_field(os)
        call clear_samples(os)
        call sample_half(os, 1, 1)
        call assert_false(cont_converged(conv_fresh, os, 1., 1, .true.), 'the first half converged alone')
        call sample_half(os, 2, 2)
        call os%set(6, 'dist', 0.6)
        call os%set(7, 'dist', 0.6)
        call assert_false(cont_converged(conv_fresh, os, 1., 1, .true.), 'a failing half converged')
        call assert_real(0., conv_fresh%get('cont_coverage_pct'), TOL, 'a failing iteration did not reset the streak')
        call sample_half(os, 1, 3)
        call assert_false(cont_converged(conv_fresh, os, 1., 1, .true.), 'the streak survived a failing iteration')
        call assert_real(50., conv_fresh%get('cont_coverage_pct'), TOL, 'the coverage after the reset')
        call os%kill
    end subroutine test_cont_rule_streak

    !> the parameters cont_conv_rot, cont_conv_shift and cont_conv_frac replace the defaults: a 0.4 deg
    !! move fails a 0.3 deg bound, a 0.5 A shift fails 0.4 A, eight of ten stable pass a 0.8 fraction,
    !! and the tight tomography bounds (0.05 deg, 0.1 A) accept only small motion; fresh histories
    subroutine test_cont_rule_tunable()
        class(parameters), allocatable :: defaults
        type(oris) :: os
        allocate(defaults)
        call assert_true(defaults%cont_conv_rot == CONT_CONV_ROT_DEG .and. defaults%cont_conv_shift == CONT_CONV_SHIFT_A &
            &.and. defaults%cont_conv_frac == CONT_CONV_FRAC, 'the parameter defaults are not the simple_defs_conv constants')
        deallocate(defaults)
        call make_cont_field(os)
        call assert_true(fresh_cont_converged(os),                 'the default bounds did not converge')
        call assert_false(fresh_cont_converged(os, rot_lim=0.3),   'a 0.4 deg move passed a 0.3 deg bound')
        call assert_false(fresh_cont_converged(os, shift_lim=0.4), 'a 0.5 A shift passed a 0.4 A bound')
        call assert_true(fresh_cont_converged(os, shift_lim=0.6),  'a 0.5 A shift failed a 0.6 A bound')
        call os%set(1, 'dist', 0.4)
        call os%set(2, 'dist', 0.4)
        call assert_false(fresh_cont_converged(os),                'eight of ten stable passed the default fraction')
        call assert_true(fresh_cont_converged(os, frac_lim=0.8),   'eight of ten stable failed a 0.8 fraction')
        call make_cont_field(os, CARTFT_ACCEPTED, 0.04, 0.05)
        call assert_true(fresh_cont_converged(os, rot_lim=0.05, shift_lim=0.1), 'small motion failed the tight bounds')
        call make_cont_field(os)
        call assert_false(fresh_cont_converged(os, rot_lim=0.05, shift_lim=0.1), 'a 0.4 deg move passed the tight bounds')
        call os%kill
    end subroutine test_cont_rule_tunable

    !> check_conv3D of a refine=cont pass at smpd 1 with a fresh history and optional bounds
    logical function fresh_cont_converged( os, rot_lim, shift_lim, frac_lim ) result( converged )
        type(oris),     intent(inout) :: os
        real, optional, intent(in)    :: rot_lim, shift_lim, frac_lim
        type(convergence) :: conv
        converged = cont_converged(conv, os, 1., 1, .false., rot_lim=rot_lim, shift_lim=shift_lim, frac_lim=frac_lim)
    end function fresh_cont_converged

    !> NCONT particles of state 1, all sampled in generation 1, every transaction ending with status
    !! (default accepted), rotation rot/2 + rot/2 deg (default 0.2 + 0.2) and shift sh px (default 0.5)
    subroutine make_cont_field( os, status, rot, sh )
        type(oris),        intent(inout) :: os
        integer, optional, intent(in)    :: status
        real,    optional, intent(in)    :: rot, sh
        integer :: i, istatus
        real    :: rot_here, sh_here
        istatus  = CARTFT_ACCEPTED
        rot_here = 0.4
        sh_here  = 0.5
        if( present(status) ) istatus  = status
        if( present(rot)    ) rot_here = rot
        if( present(sh)     ) sh_here  = sh
        call os%kill
        call os%new(NCONT, is_ptcl=.true.)
        do i = 1, NCONT
            call os%set_state(i, 1)
            call os%set(i, 'updatecnt',          1.)
            call os%set(i, 'sampled',            1.)
            call os%set(i, 'frac',               100.)
            call os%set(i, 'corr_cart',          0.4)
            call os%set(i, 'pose_cont_status',   real(istatus))
            call os%set(i, 'pose_cont_improved', merge(1., 0., istatus == CARTFT_ACCEPTED))
            call os%set(i, 'dist',               0.5*rot_here)
            call os%set(i, 'dist_inpl',          0.5*rot_here)
            call os%set(i, 'shincarg',           sh_here)
        end do
    end subroutine make_cont_field

    !> no particle sampled or updated yet
    subroutine clear_samples( os )
        type(oris), intent(inout) :: os
        integer :: i
        do i = 1, NCONT
            call os%set(i, 'sampled',   0.)
            call os%set(i, 'updatecnt', 0.)
        end do
    end subroutine clear_samples

    !> half ihalf (1: particles 1-5, 2: 6-10) sampled in generation gen
    subroutine sample_half( os, ihalf, gen )
        type(oris), intent(inout) :: os
        integer,    intent(in)    :: ihalf, gen
        integer :: i
        do i = (ihalf - 1)*NCONT/2 + 1, ihalf*NCONT/2
            call os%set(i, 'sampled',   real(gen))
            call os%set(i, 'updatecnt', real(gen))
        end do
    end subroutine sample_half

    !> check_conv3D of a refine=cont pass on os, with the history conv carries
    logical function cont_converged( conv, os, smpd, nstates, l_frac, minits, rot_lim, shift_lim, frac_lim ) result( converged )
        type(convergence), intent(inout) :: conv
        type(oris),        intent(inout) :: os
        real,              intent(in)    :: smpd
        integer,           intent(in)    :: nstates
        logical,           intent(in)    :: l_frac
        integer, optional, intent(in)    :: minits
        real,    optional, intent(in)    :: rot_lim, shift_lim, frac_lim
        class(parameters), allocatable :: params
        type(cmdline) :: cline
        allocate(params)
        params%nstates       = nstates
        params%smpd          = smpd
        params%trs           = 5.
        params%refine        = 'cont'
        params%l_cart_refine = .true.
        params%l_update_frac = l_frac
        params%startit       = 1
        params%which_iter    = 1
        if( present(minits)    ) params%minits          = minits
        if( present(rot_lim)   ) params%cont_conv_rot   = rot_lim
        if( present(shift_lim) ) params%cont_conv_shift = shift_lim
        if( present(frac_lim)  ) params%cont_conv_frac  = frac_lim
        call cline%set('trs', 5.)
        converged = conv%check_conv3D(params, cline, os, 40.)
        call cline%kill
        deallocate(params)
    end function cont_converged

end module simple_convergence_tester
