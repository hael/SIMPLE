!@descr: unit tests for the 3D convergence report and the refine=cont convergence rule (simple_convergence)
! N19: a Cartesian pass reports the mean corr_cart, attempts and improved percentage over its sample;
! a polar pass reports corr. N25: a refine=cont pass converges when 90% of every state moved at most
! 0.5 deg and 1 A with an accepted or no-improvement outcome, over a streak covering 90% of the active
! particles; a discrete pass, with or without the polish, keeps its own rule. Expected values are
! hand-counted from the fixtures.
module simple_convergence_tester
use simple_core_module_api, only: oris, del_file
use simple_cmdline,         only: cmdline
use simple_convergence,     only: convergence, cont_particle_stable
use simple_parameters,      only: parameters
use simple_defs_fname,      only: STATS_FILE
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
        call del_file(STATS_FILE)
    end subroutine run_all_convergence_tests

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

    !> the bounds are inclusive; only an accepted or a finite no-improvement outcome is stable
    subroutine test_cont_particle_stable()
        integer :: rejected(6), i
        call assert_true(cont_particle_stable(CARTFT_ACCEPTED, 0.5, 1.0), 'the bounds are inclusive')
        call assert_false(cont_particle_stable(CARTFT_ACCEPTED, 0.51, 0.), 'a rotation above 0.5 deg is stable')
        call assert_false(cont_particle_stable(CARTFT_ACCEPTED, 0., 1.01), 'a shift above 1 A is stable')
        call assert_true(cont_particle_stable(CARTFT_NO_IMPROVEMENT, 0., 0.), 'a finite no-improvement is not stable')
        rejected = [CARTFT_NOT_ATTEMPTED, CARTFT_INVALID_PREPARATION, CARTFT_NO_RELIABLE_UPDATE, &
            &CARTFT_BOUND_REJECTED, CARTFT_INVALID_NUMERICS, CARTFT_ITERATION_LIMIT]
        do i = 1, size(rejected)
            call assert_false(cont_particle_stable(rejected(i), 0., 0.), 'a rejected transaction at zero motion is stable')
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
    logical function cont_converged( conv, os, smpd, nstates, l_frac, minits ) result( converged )
        type(convergence), intent(inout) :: conv
        type(oris),        intent(inout) :: os
        real,              intent(in)    :: smpd
        integer,           intent(in)    :: nstates
        logical,           intent(in)    :: l_frac
        integer, optional, intent(in)    :: minits
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
        if( present(minits) ) params%minits = minits
        call cline%set('trs', 5.)
        converged = conv%check_conv3D(params, cline, os, 40.)
        call cline%kill
        deallocate(params)
    end function cont_converged

end module simple_convergence_tester
