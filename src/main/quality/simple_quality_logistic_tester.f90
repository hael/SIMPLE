!@descr: unit tests for the bank-free logistic core (simple_quality_logistic) and the rank-based AUC
submodule (simple_quality_logistic) simple_quality_logistic_tester
use simple_test_utils
use simple_defs,               only: dp
use simple_cavg_quality_stats, only: auc_for_values
use simple_linalg,             only: cholesky
implicit none

integer, parameter :: NTOY = 400
character(len=6), parameter :: SOLVERS(2) = ['lbfgsb', 'lm    ']

contains

    module subroutine run_all_quality_logistic_tests()
        write(*,'(A)') '**** running all quality logistic tests ****'
        call test_toy_plane_recovered()
        call test_fisher_is_hessian()
        call test_fit_repeatable()
        call test_roles_weights_macro()
        call test_auc_trajectory_scale()
        call test_auc_matches_pairwise()
        call test_auc_eps_ties()
    end subroutine run_all_quality_logistic_tests

    !> labels from the plane x1 - 2 x2 + 0.5 = 0 with every 20th label flipped: the fitted linear part points along
    !! the plane's normal, and the fit classifies the clean rows, under both solvers; a comparison row per solver with
    !! the median seconds of NTIME identical fits
    subroutine test_toy_plane_recovered()
        use simple_stat, only: median
        integer, parameter :: NTIME = 20
        type(quality_logistic_problem) :: prob
        real(kind=8), allocatable :: sol(:)
        real    :: x(2,NTOY), obj, cosang
        logical :: good(NTOY)
        integer :: i, ncorrect, nclean, isolver, nevals, niter, irep
        integer(kind=8) :: c0, c1, crate
        real(kind=8)    :: times(NTIME)
        write(*,'(A)') 'test_toy_plane_recovered'
        do isolver = 1, size(SOLVERS)
            call toy_problem(prob, x, good, trim(SOLVERS(isolver)))
            do irep = 1, NTIME
                call system_clock(c0, crate)
                call prob%fit(sol, obj, nevals, niter)
                call system_clock(c1)
                times(irep) = real(c1 - c0, 8) / real(crate, 8)
            end do
            cosang = real((sol(2) * 1.d0 - sol(3) * 2.d0) / (sqrt(sol(2)**2 + sol(3)**2) * sqrt(5.d0)))
            call assert_true(cosang > 0.99, trim(SOLVERS(isolver))//': linear coefficients along the plane normal')
            call assert_true(abs(sol(4)) < 0.25 * sqrt(sol(2)**2 + sol(3)**2), &
                &trim(SOLVERS(isolver))//': interaction coefficient small')
            ncorrect = 0
            nclean   = 0
            do i = 1, NTOY
                if( mod(i, 20) == 0 ) cycle
                nclean = nclean + 1
                if( (stable_sigmoid_8(sol(1) + sol(2) * x(1,i) + sol(3) * x(2,i) + sol(4) * real(x(1,i) * x(2,i), 8)) &
                    &> 0.5d0) .eqv. good(i) ) ncorrect = ncorrect + 1
            end do
            call assert_true(real(ncorrect) / real(nclean) > 0.95, trim(SOLVERS(isolver))//': clean rows classified')
            call assert_true(obj > 0.0, trim(SOLVERS(isolver))//': objective positive')
            write(*,'(A,A6,A,I4,A,I4,A,ES10.3,A,F9.6,A,F9.6,A,ES16.8)') '  LM_COMPARE quality_logistic ', &
                &SOLVERS(isolver), ' evals ', nevals, ' iters ', niter, ' s median ', &
                &median(times), ' cos ', cosang, ' clean ', real(ncorrect) / real(nclean), ' obj ', obj
            call prob%kill
        end do
    end subroutine test_toy_plane_recovered

    !> the curvature is the exact Hessian: the central difference of the analytic gradient within 1e-6 of
    !! max(1, max|C|) (guessed), at the start and at a point away from it, and positive definite
    subroutine test_fisher_is_hessian()
        type(quality_logistic_problem) :: prob
        real(kind=8), allocatable :: vec(:), gp(:), gm(:), curv(:,:), hfd(:,:), chol(:,:)
        real(kind=8), parameter :: STEP = 1.d-5
        real(kind=8) :: f, err
        real    :: x(2,NTOY)
        logical :: good(NTOY), ok
        integer :: ipoint, j, d
        write(*,'(A)') 'test_fisher_is_hessian'
        call toy_problem(prob, x, good)
        d = prob%ndim
        allocate(vec(d), gp(d), gm(d), curv(d,d), hfd(d,d), chol(d,d))
        do ipoint = 1, 2
            vec = 0.d0
            if( ipoint == 2 ) vec = [0.3d0, 1.7d0, -2.9d0, 0.4d0]
            call prob%curvature(vec, curv, d)
            do j = 1, d
                vec(j) = vec(j) + STEP
                call prob%fdf(vec, f, gp, d)
                vec(j) = vec(j) - 2.d0 * STEP
                call prob%fdf(vec, f, gm, d)
                vec(j) = vec(j) + STEP
                hfd(:,j) = (gp - gm) / (2.d0 * STEP)
            end do
            err = maxval(abs(curv - hfd)) / max(1.d0, maxval(abs(curv)))
            write(*,'(A,I0,A,ES10.3)') '  Fisher oracle at point ', ipoint, ': max relative difference to the Hessian ', err
            call assert_true(err < 1.d-6, 'the Fisher matrix plus the ridge is the Hessian')
            call assert_true(all(curv == transpose(curv)), 'the curvature is symmetric')
            call cholesky(curv, chol, d, ok)
            call assert_true(ok, 'the curvature is positive definite')
        end do
        call prob%kill
    end subroutine test_fisher_is_hessian

    !> the same problem fitted twice gives the same solution, bit for bit
    subroutine test_fit_repeatable()
        type(quality_logistic_problem) :: prob
        real(kind=8), allocatable :: sol1(:), sol2(:)
        real    :: x(2,NTOY), obj1, obj2
        logical :: good(NTOY)
        write(*,'(A)') 'test_fit_repeatable'
        call toy_problem(prob, x, good)
        call prob%fit(sol1, obj1)
        call prob%fit(sol2, obj2)
        call assert_true(all(sol1 == sol2), 'solutions identical')
        call assert_true(obj1 == obj2, 'objectives identical')
        call prob%kill
        call assert_int(0, prob%nobs, 'killed')
    end subroutine test_fit_repeatable

    subroutine test_roles_weights_macro()
        real(kind=8) :: gw, bw
        write(*,'(A)') 'test_roles_weights_macro'
        call assert_int(LEARN_ROLE_BALANCED,         dataset_learn_role(3, 2, 0), 'balanced')
        call assert_int(LEARN_ROLE_RECALL_ONLY,      dataset_learn_role(3, 0, 0), 'good only')
        call assert_int(LEARN_ROLE_SPECIFICITY_ONLY, dataset_learn_role(0, 2, 0), 'bad only')
        call assert_int(LEARN_ROLE_SKIP,             dataset_learn_role(0, 2, 1), 'bad only with hard-rejected good')
        call logistic_dataset_class_weights(LEARN_ROLE_BALANCED, 4, 7, gw, bw)
        call assert_double(real(0.65, 8) / 4.d0, gw, 'balanced good weight', 4.d0)
        call assert_double(real(1.0 - 0.65, 8) / 7.d0, bw, 'balanced bad weight', 4.d0)
        ! four scores: mean 0.6, tail (one score) 0.3 -> 0.45
        call assert_real(0.45, robust_macro_score([0.9, 0.3, 0.6, 0.6]), 1.e-6, 'robust macro score')
    end subroutine test_roles_weights_macro

    !> perfect separation at the size of the first trajectory (5,336 good, 7,635 bad: 40.7 million pairs, beyond a
    !! single-precision counter) is 1, its reverse 0, a constant score 1/2
    subroutine test_auc_trajectory_scale()
        integer, parameter :: NG = 5336, NB = 7635
        real    :: v(NG+NB)
        integer :: st(NG+NB), i
        write(*,'(A)') 'test_auc_trajectory_scale'
        do i = 1, NG + NB
            if( i <= NG )then
                st(i) = 1
                v(i)  = 10. + real(i) * 1.e-3
            else
                st(i) = 0
                v(i)  = real(i - NG) * 1.e-3
            endif
        end do
        call assert_real(1.0, auc_for_values(v, st), 0.0, 'perfect separation is 1')
        call assert_real(0.0, auc_for_values(-v, st), 0.0, 'reversed separation is 0')
        v = 3.
        call assert_real(0.5, auc_for_values(v, st), 0.0, 'constant score is 1/2')
    end subroutine test_auc_trajectory_scale

    !> on a small sample with ties the rank count equals the pairwise count of the class-average definition
    subroutine test_auc_matches_pairwise()
        integer, parameter :: N = 301
        real    :: v(N), wins
        integer :: st(N), i, j, seed, ngood
        write(*,'(A)') 'test_auc_matches_pairwise'
        seed = 12345
        do i = 1, N
            seed  = mod(seed * 16807, 2147483647)
            v(i)  = real(mod(seed, 50)) * 0.1          ! coarse values: many exact ties
            seed  = mod(seed * 16807, 2147483647)
            st(i) = merge(1, 0, mod(seed, 3) == 0)
        end do
        wins  = 0.
        ngood = count(st > 0)
        do i = 1, N
            if( st(i) <= 0 ) cycle
            do j = 1, N
                if( st(j) > 0 ) cycle
                if( v(i) > v(j) )then
                    wins = wins + 1.
                else if( abs(v(i) - v(j)) <= 1.e-6 )then
                    wins = wins + 0.5
                endif
            end do
        end do
        call assert_real(wins / real(ngood * (N - ngood)), auc_for_values(v, st), 0.0, 'equals the pairwise count')
    end subroutine test_auc_matches_pairwise

    !> a bad value at most EPS above a good one counts 1/2, one further above counts 0, one below counts 1
    subroutine test_auc_eps_ties()
        write(*,'(A)') 'test_auc_eps_ties'
        call assert_real(0.5, auc_for_values([1.0, 1.0 + 5.e-7], [1, 0]), 0.0, 'within EPS above: 1/2')
        call assert_real(0.0, auc_for_values([1.0, 1.0 + 1.e-4], [1, 0]), 0.0, 'beyond EPS above: 0')
        call assert_real(1.0, auc_for_values([1.0, 1.0 - 5.e-7], [1, 0]), 0.0, 'below: 1')
    end subroutine test_auc_eps_ties

    ! the toy design: two features on a deterministic lattice in [-2, 2]^2, balanced weights, lambda 1e-4
    subroutine toy_problem( prob, x, good, solver )
        type(quality_logistic_problem), intent(inout) :: prob
        real,                           intent(out)   :: x(2,NTOY)
        logical,                        intent(out)   :: good(NTOY)
        character(len=*), optional,     intent(in)    :: solver
        real(kind=8) :: gw, bw
        integer :: i, ng
        do i = 1, NTOY
            x(1,i)  = -2. + 4. * real(mod(i * 37, 101)) / 100.
            x(2,i)  = -2. + 4. * real(mod(i * 53, 97)) / 96.
            good(i) = x(1,i) - 2. * x(2,i) + 0.5 > 0.
            if( mod(i, 20) == 0 ) good(i) = .not. good(i)
        end do
        ng = count(good)
        call logistic_dataset_class_weights(LEARN_ROLE_BALANCED, ng, NTOY - ng, gw, bw)
        call prob%new([1, 2], NTOY, .false., 1.e-4, solver)
        do i = 1, NTOY
            call prob%set_row(i, x(:,i), good(i), merge(gw, bw, good(i)))
        end do
        call prob%finalize()
    end subroutine toy_problem

end submodule simple_quality_logistic_tester
