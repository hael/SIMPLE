!@descr: unit tests for the pass scores of the 3D convergence report (simple_convergence)
! N19 (C14, C15, O8): a Cartesian pass reports the mean of corr_cart over the particles it
! sampled, its attempts as their count and its improved percentage from the improved flag; a
! polar pass reports corr and no pose statistics. Expected values: hand-counted from the field.
! N25 (O2 interim rule, section 7.2): with every particle sampled, frac = 100 and no motion,
! check_conv3D declares no convergence under refine=cont for any overlap and search-fraction
! limit, and still computes the motion statistics (dist, shift increment, improved fraction);
! the same field under a discrete mode converges as before. Expected values: the rule of 7.2
! and hand means over the field. Phase 9: a discrete pass followed by the polish (pose_cont=yes)
! keeps the discrete rule and its motion statistics, and reports the polish's attempts and
! improved percentage over the sample. Expected values: the same hand counts.
module simple_convergence_tester
use simple_core_module_api, only: oris, del_file
use simple_cmdline,         only: cmdline
use simple_convergence,     only: convergence
use simple_parameters,      only: parameters
use simple_defs_fname,      only: STATS_FILE
use simple_test_utils
implicit none
private
public :: run_all_convergence_tests

integer, parameter :: NPTCLS = 6
real,    parameter :: TOL    = 1.e-5

contains

    subroutine run_all_convergence_tests()
        write(*,'(A)') '**** running all convergence tests ****'
        write(*,'(A)') 'test_cartesian_pass_scores'
        call test_cartesian_pass_scores()
        write(*,'(A)') 'test_polar_pass_scores'
        call test_polar_pass_scores()
        write(*,'(A)') 'test_cartesian_pass_never_converges'
        call test_cartesian_pass_never_converges()
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

    subroutine test_cartesian_pass_never_converges()
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
            call assert_false(converged, 'a Cartesian pass declared convergence')
            ! the motion statistics are computed all the same: (0.1+...+0.6)/6, (0.01+...+0.06)/6
            call assert_real(0.35, conv%get('dist'), TOL, 'the Cartesian pass mean orientation motion')
            call assert_real(50., conv%get('pose_cont_improved_pct'), TOL, 'the Cartesian pass improved fraction')
            call assert_real(real(NPTCLS), conv%get('pose_cont_attempts'), TOL, 'the Cartesian pass attempts')
            ! the same field under a discrete mode: the rule of 7.2 is unchanged
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
    end subroutine test_cartesian_pass_never_converges

end module simple_convergence_tester
