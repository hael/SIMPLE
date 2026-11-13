!@descr: class-average sums: the weighted accumulation of one particle plane at several in-plane angles (topk_inpl)
submodule (simple_classaverager) simple_classaverager_insert_tester
use simple_test_utils, only: assert_real, assert_true
implicit none

contains

    module subroutine run_all_cavger_insert_tests()
        write(*,'(A)') '**** running all class-average in-plane entry tests ****'
        call test_weighted_accumulation()
    end subroutine run_all_cavger_insert_tests

    !> One particle plane gridded at in-plane angles with weights: a weight of one is the unweighted insertion bit
    !! for bit; a power-of-two weight scales both sums exactly; two halves at one angle equal the whole; a set of
    !! weights summing to one is the weighted sum of the single insertions (the insertion is linear in the plane)
    !! and conserves the CTF^2 mass; the angle is periodic; a vanishing weight adds nothing. The plane is
    !! asymmetric (a real part tilted in h, an imaginary part in k) so that a wrong sign or a lost angle shows
    subroutine test_weighted_accumulation()
        integer,  parameter :: LDIM = 32
        real,     parameter :: ANGLE_A = 30., ANGLE_B = 75.
        real,     parameter :: TOL_SUM = 1.e-5  ! single-precision sums of a few hundred terms, relative to the peak
        real,     parameter :: TOL_PERIODIC = 1.e-4 ! the angle in single-precision radians moves a window position by ~1e-5 pixel
        type(stack)       :: stk
        type(fplane_type) :: fpl
        real    :: peak_c, peak_t, diff_c, diff_t, mass_a, mass_b
        integer :: lim, h, k
        write(*,'(A)') 'test_weighted_accumulation'
        call stk%new_stack([LDIM,LDIM], 8)
        lim = maxval(abs(stk%flims(1:2,:))) + 1
        allocate(fpl%cmplx_plane(-lim:lim,-lim:0), fpl%ctfsq_plane(-lim:lim,-lim:0))
        do k = -lim,0
            do h = -lim,lim
                fpl%cmplx_plane(h,k) = exp(-real(h*h+k*k)/200.) * cmplx(1.+0.5*real(h)/real(lim), 0.3*real(k)/real(lim))
                fpl%ctfsq_plane(h,k) = 0.4 + 0.3*cos(0.4*real(h))*cos(0.3*real(k))
            enddo
        enddo
        ! slice 1: unweighted at A; 2: weight 1; 3: 0.5 twice at A; 4: 0.3 at A and 0.7 at B; 5: B alone;
        ! 6: weight 0.25 at A; 7: A plus a full turn; 8: weight 0
        call stk%accumulate_fplane(ANGLE_A, fpl, 1)
        call stk%accumulate_fplane(ANGLE_A, fpl, 2, 1.0)
        call stk%accumulate_fplane(ANGLE_A, fpl, 3, 0.5)
        call stk%accumulate_fplane(ANGLE_A, fpl, 3, 0.5)
        call stk%accumulate_fplane(ANGLE_A, fpl, 4, 0.3)
        call stk%accumulate_fplane(ANGLE_B, fpl, 4, 0.7)
        call stk%accumulate_fplane(ANGLE_B, fpl, 5)
        call stk%accumulate_fplane(ANGLE_A, fpl, 6, 0.25)
        call stk%accumulate_fplane(ANGLE_A + 360., fpl, 7)
        call stk%accumulate_fplane(ANGLE_A, fpl, 8, 0.)
        peak_c = maxval(abs(stk%cmat(:,:,1)))
        peak_t = maxval(stk%ctfsq(:,:,1))
        call assert_true(peak_c > 0. .and. peak_t > 0., 'the insertion at one angle deposits a numerator and a CTF^2 sum')
        call assert_true(all(stk%cmat(:,:,1) == stk%cmat(:,:,2)) .and. all(stk%ctfsq(:,:,1) == stk%ctfsq(:,:,2)), &
            &'a weight of one is the unweighted insertion bit for bit, on the numerator and the CTF^2 sum')
        call assert_true(all(stk%cmat(:,:,6) == 0.25 * stk%cmat(:,:,1)) .and. all(stk%ctfsq(:,:,6) == 0.25 * stk%ctfsq(:,:,1)), &
            &'a weight of one quarter scales the numerator and the CTF^2 sum exactly')
        diff_c = maxval(abs(stk%cmat(:,:,3) - stk%cmat(:,:,1)))
        diff_t = maxval(abs(stk%ctfsq(:,:,3) - stk%ctfsq(:,:,1)))
        call assert_real(0., diff_c / peak_c, TOL_SUM, 'two halves at one angle equal the whole, numerator')
        call assert_real(0., diff_t / peak_t, TOL_SUM, 'two halves at one angle equal the whole, CTF^2 sum')
        diff_c = maxval(abs(stk%cmat(:,:,4)  - (0.3 * stk%cmat(:,:,1)  + 0.7 * stk%cmat(:,:,5))))
        diff_t = maxval(abs(stk%ctfsq(:,:,4) - (0.3 * stk%ctfsq(:,:,1) + 0.7 * stk%ctfsq(:,:,5))))
        call assert_real(0., diff_c / peak_c, TOL_SUM, 'weights summing to one give the weighted sum of the single insertions, numerator')
        call assert_real(0., diff_t / peak_t, TOL_SUM, 'weights summing to one give the weighted sum of the single insertions, CTF^2 sum')
        mass_a = sum(stk%ctfsq(:,:,1))
        mass_b = sum(stk%ctfsq(:,:,5))
        call assert_real(0., abs(sum(stk%ctfsq(:,:,4)) - (0.3 * mass_a + 0.7 * mass_b)) / mass_a, TOL_SUM, &
            &'the CTF^2 mass of a particle at weights summing to one is that of the single insertions')
        call assert_real(0., abs(mass_a - mass_b) / mass_a, TOL_SUM, 'the KB window conserves the CTF^2 mass at any angle')
        call assert_true(maxval(abs(stk%cmat(:,:,1) - stk%cmat(:,:,5))) > 0.01 * peak_c, &
            &'two angles 45 degrees apart deposit different numerators (the angle is used)')
        diff_c = maxval(abs(stk%cmat(:,:,7) - stk%cmat(:,:,1)))
        diff_t = maxval(abs(stk%ctfsq(:,:,7) - stk%ctfsq(:,:,1)))
        call assert_real(0., max(diff_c / peak_c, diff_t / peak_t), TOL_PERIODIC, 'the in-plane angle is periodic')
        call assert_true(all(stk%cmat(:,:,8) == 0.) .and. all(stk%ctfsq(:,:,8) == 0.), 'a vanishing weight adds nothing')
        deallocate(fpl%cmplx_plane, fpl%ctfsq_plane)
        call stk%kill_stack
    end subroutine test_weighted_accumulation

end submodule simple_classaverager_insert_tester
