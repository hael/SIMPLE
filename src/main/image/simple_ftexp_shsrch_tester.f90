!@descr: unit tests for the expanded-Fourier shift search (simple_ftexp_shsrch): correlator peak, gradients and optimiser
! Replaces the in-module test_ftexp_shsrch and test_ftexp_shsrch2 (called private cost functions,
! stopped at the first failure; the optimiser test never looked at the shift it found). The
! callbacks are reached through the public ospec, as the optimiser reaches them: costfun_8 returns
! minus the normalised correlation, fdfcostfun_8 minus the correlation times 1e8 (the scale L-BFGS-B
! works in) and its gradient. Convention pinned: with the reference as the unshifted image and the
! particle shifted by image%shift(s), the cost is lowest at vec = s, and minimize returns s.
module simple_ftexp_shsrch_tester
use simple_test_utils
use simple_test_image_fixture, only: draw_square
use simple_defs
use simple_rnd,          only: ran3
use simple_image,        only: image, unmemoize_mask_coords
use simple_ft_expanded,  only: ft_expanded, ftexp_transfmat_init, ftexp_transfmat_kill
use simple_ftexp_shsrch, only: ftexp_shsrch
implicit none
private
public :: run_all_ftexp_shsrch_tests

real(dp), parameter :: COST_SCALE = 1.d8 ! fdfcostfun_8 works in the correlation times this
integer,  parameter :: NTST_OPT   = 100
real    :: opt_shifts(2,NTST_OPT)            ! test_optimiser's shifts, drawn once: both solvers see the same
logical :: opt_shifts_drawn = .false.

contains

    subroutine run_all_ftexp_shsrch_tests()
        write(*,'(A)') '**** running all ftexp_shsrch tests ****'
        call test_correlator()
        call test_curvature_oracle()
        call test_optimiser('lbfgsb')
        call test_optimiser('lm')
        call ftexp_transfmat_kill
    end subroutine run_all_ftexp_shsrch_tests

    ! a square of half-width 32 in a 128 box, shifted by integer (x,y) within +/-10: an exhaustive
    ! integer search of the cost finds (x,y), where the correlation is 1, the fdf callback agrees with
    ! the cost callback and the gradient vanishes; off the peak the analytic gradient matches central
    ! differences
    subroutine test_correlator()
        integer,  parameter :: LDIM(3) = [128,128,1], SQRAD = 32, NTST = 10, ITRS = 10
        real,     parameter :: SMPD = 1.77, TRS = real(ITRS), HP = 100., LP = 8.
        real(dp), parameter :: FD_STEP = 1.d-3
        type(image)        :: img, img_shifted
        type(ft_expanded)  :: ftexp_ref, ftexp_ptcl
        type(ftexp_shsrch) :: srch
        real(dp) :: vec(2), cost, cost_best, f, fp, fm, grad(2), grad_fd(2), gp(2)
        real(dp) :: corr_min, fdf_diff_max, grad_mag_avg, fd_err_max
        integer  :: itst, x, y, xsh, ysh, xbest, ybest, nwrong, d
        write(*,'(A)') 'test_correlator'
        call img%new(LDIM, SMPD)
        call draw_square(img, SQRAD)
        call ftexp_ref%new(img, HP, LP, .true.)
        call ftexp_transfmat_init(img, LP)
        call img_shifted%new(LDIM, SMPD)
        nwrong       = 0
        corr_min     = huge(corr_min)
        fdf_diff_max = 0.d0
        grad_mag_avg = 0.d0
        fd_err_max   = 0.d0
        do itst = 1,NTST
            x = nint(ran3() * 2. * TRS - TRS)
            y = nint(ran3() * 2. * TRS - TRS)
            img_shifted = img
            call img_shifted%shift([real(x), real(y), 0.])
            call ftexp_ptcl%new(img_shifted, HP, LP, .true.)
            call srch%new(ftexp_ref, ftexp_ptcl, TRS)
            cost_best = huge(cost_best)
            xbest     = 0
            ybest     = 0
            do xsh = -ITRS,ITRS
                do ysh = -ITRS,ITRS
                    vec  = real([xsh, ysh], dp)
                    cost = srch%ospec%costfun_8(srch, vec, 2)
                    if( cost < cost_best )then
                        cost_best = cost
                        xbest     = xsh
                        ybest     = ysh
                    endif
                end do
            end do
            if( xbest /= x .or. ybest /= y ) nwrong = nwrong + 1
            corr_min = min(corr_min, -cost_best)
            ! the fdf callback at the peak
            vec = real([xbest, ybest], dp)
            call srch%ospec%fdfcostfun_8(srch, vec, f, grad, 2)
            fdf_diff_max = max(fdf_diff_max, abs(f / COST_SCALE - cost_best))
            grad_mag_avg = grad_mag_avg + sqrt(sum((grad / COST_SCALE)**2)) / real(NTST, dp)
            ! off the peak: central differences of the fdf cost
            vec = real([x, y], dp) + [0.3d0, -0.4d0]
            call srch%ospec%fdfcostfun_8(srch, vec, f, grad, 2)
            do d = 1,2
                vec(d) = vec(d) + FD_STEP
                call srch%ospec%fdfcostfun_8(srch, vec, fp, gp, 2)
                vec(d) = vec(d) - 2.d0 * FD_STEP
                call srch%ospec%fdfcostfun_8(srch, vec, fm, gp, 2)
                vec(d) = vec(d) + FD_STEP
                grad_fd(d) = (fp - fm) / (2.d0 * FD_STEP)
            end do
            fd_err_max = max(fd_err_max, sqrt(sum((grad - grad_fd)**2)) / sqrt(sum(grad**2)))
            call srch%kill
            call ftexp_ptcl%kill
        end do
        call assert_int(0, nwrong, 'the cost is lowest at the applied integer shift in every trial')
        call assert_true(corr_min > 0.999d0, 'the correlation at the peak is 1 (above 0.999)')
        call assert_true(fdf_diff_max < 1.d-5, 'fdfcostfun_8 / 1e8 equals costfun_8 at the peak')
        call assert_true(grad_mag_avg < 1.d-6, 'the gradient vanishes at the peak (mean magnitude below 1e-6)')
        call assert_true(fd_err_max < 1.d-3, 'off the peak the analytic gradient matches central differences (relative 1e-3)')
        call ftexp_ref%kill
        call img%kill
        call img_shifted%kill
    end subroutine test_correlator

    ! the Levenberg-Marquardt curvature at the peak of an integer-shifted square (zero residual): the central
    ! difference of the analytic gradient within the gradient check's relative 1e-3, and positive semi-definite
    subroutine test_curvature_oracle()
        integer,  parameter :: LDIM(3) = [128,128,1], SQRAD = 32
        real,     parameter :: SMPD = 1.77, TRS = 10., HP = 100., LP = 8.
        real(dp), parameter :: FD_STEP = 1.d-3
        type(image)        :: img, img_shifted
        type(ft_expanded)  :: ftexp_ref, ftexp_ptcl
        type(ftexp_shsrch) :: srch
        real(dp) :: vec(2), f, grad(2), gp(2), gm(2), curv(2,2), hfd(2,2), rel_err
        real     :: cxy(3)
        integer  :: d
        write(*,'(A)') 'test_curvature_oracle'
        call img%new(LDIM, SMPD)
        call draw_square(img, SQRAD)
        call ftexp_ref%new(img, HP, LP, .true.)
        call ftexp_transfmat_init(img, LP)
        call img_shifted%new(LDIM, SMPD)
        img_shifted = img
        call img_shifted%shift([3., -4., 0.])
        call ftexp_ptcl%new(img_shifted, HP, LP, .true.)
        call srch%new(ftexp_ref, ftexp_ptcl, TRS, solver='lm')
        cxy = srch%minimize()          ! prepares the solve's terms, the curvature among them
        vec = [3.d0, -4.d0]
        call srch%ospec%fdfcurvfun_8(srch, vec, f, grad, curv, 2)
        do d = 1,2
            vec(d) = vec(d) + FD_STEP
            call srch%ospec%fdfcostfun_8(srch, vec, f, gp, 2)
            vec(d) = vec(d) - 2.d0 * FD_STEP
            call srch%ospec%fdfcostfun_8(srch, vec, f, gm, 2)
            vec(d) = vec(d) + FD_STEP
            hfd(:,d) = (gp - gm) / (2.d0 * FD_STEP)
        end do
        rel_err = sqrt(sum((curv - hfd)**2)) / sqrt(sum(hfd**2))
        write(logfhandle,'(A,ES10.3,A,2ES12.4)') 'curvature oracle: relative difference to the Hessian ', rel_err, &
            &', diagonal ', curv(1,1), curv(2,2)
        call assert_true(rel_err < 1.d-3, 'at zero residual the curvature is the Hessian of the cost (relative 1e-3)')
        call assert_true(curv(1,1) >= 0.d0 .and. curv(2,2) >= 0.d0 .and. &
            &curv(1,1)*curv(2,2) - curv(1,2)*curv(2,1) >= 0.d0, 'the curvature is positive semi-definite')
        call assert_true(curv(1,2) == curv(2,1), 'the curvature is symmetric')
        call srch%kill
        call ftexp_ptcl%kill
        call ftexp_ref%kill
        call img%kill
        call img_shifted%kill
    end subroutine test_curvature_oracle

    ! a soft disc of radius 8 in a 32 box, shifted by random sub-pixel amounts within +/-5 in Fourier
    ! space: minimize finds the shift and a correlation of 1, under either solver; a comparison row per solver
    subroutine test_optimiser( solver )
        character(len=*), intent(in) :: solver
        integer, parameter :: NTST = NTST_OPT
        real,    parameter :: TRS = 5., LP = 6., HP = 100.
        type(image)        :: img_ref, img_ptcl
        type(ft_expanded)  :: ftexp_ref, ftexp_ptcl
        type(ftexp_shsrch) :: srch
        real    :: cxy(3), x, y, corr_min, sh_err_max
        integer :: i, nevals, niter
        integer(kind=8) :: c0, c1, crate, csum
        write(*,'(A)') 'test_optimiser '//solver
        call img_ref%new([32,32,1], 2.)
        call img_ptcl%new([32,32,1], 2.)
        img_ref = 1.
        call img_ref%memoize_mask_coords
        call img_ref%mask2D_soft(8., backgr=0.)
        call img_ref%fft()
        call ftexp_ref%new(img_ref, HP, LP, .true.)
        call ftexp_ptcl%new(img_ptcl, HP, LP, .false.)
        call srch%new(ftexp_ref, ftexp_ptcl, TRS, solver=solver)
        call ftexp_transfmat_init(img_ref, LP)
        corr_min   = huge(corr_min)
        sh_err_max = 0.
        nevals     = 0
        niter      = 0
        csum       = 0
        do i = 1,NTST
            if( .not. opt_shifts_drawn )then
                opt_shifts(1,i) = ran3() * 2. * TRS - TRS
                opt_shifts(2,i) = ran3() * 2. * TRS - TRS
            endif
            x = opt_shifts(1,i)
            y = opt_shifts(2,i)
            img_ptcl = img_ref
            call img_ptcl%shift([x, y, 0.])
            call ftexp_ptcl%new(img_ptcl, HP, LP, .true.)
            call system_clock(c0, crate)
            cxy        = srch%minimize()
            call system_clock(c1)
            csum       = csum + (c1 - c0)
            nevals     = nevals + srch%ospec%nevals
            niter      = niter  + srch%ospec%niter
            corr_min   = min(corr_min, cxy(1))
            sh_err_max = max(sh_err_max, maxval(abs(cxy(2:3) - [x, y])))
        end do
        opt_shifts_drawn = .true.
        write(logfhandle,'(A,F8.5,A,ES10.3)') 'optimiser: lowest correlation ', corr_min, ', largest shift error ', sh_err_max
        write(logfhandle,'(A,A6,A,F7.2,A,F6.2,A,ES10.3,A,F8.5,A,ES10.3)') '  LM_COMPARE ftexp_shsrch ', solver, &
            &' evals/solve ', real(nevals)/real(NTST), ' iters/solve ', real(niter)/real(NTST), ' s/solve ', &
            &real(csum, 8)/real(crate, 8)/real(NTST, 8), ' corr_min ', corr_min, ' sh_err_max ', sh_err_max
        call assert_true(corr_min >= 0.999, 'minimize reaches a correlation of 1 (at least 0.999) in every trial')
        call assert_true(sh_err_max < 0.05, 'minimize returns the applied sub-pixel shift (within 0.05 px)')
        call srch%kill
        call ftexp_ref%kill
        call ftexp_ptcl%kill
        call img_ref%kill
        call img_ptcl%kill
        call unmemoize_mask_coords
    end subroutine test_optimiser

end module simple_ftexp_shsrch_tester
