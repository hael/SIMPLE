!@descr: unit test routines for the annealing and particle-sampling schedules (simple_decay_funs)
! The cosine decays that drive eps/lambda in refine3D, the sampling schedules behind the stochastic update
! fraction and the extremal-search neighbourhood fractions of the 2D/3D probabilistic searches.
module simple_decay_funs_tester
use simple_test_utils      ! assertions etc.
use simple_defs            ! NSAMPLE_MINMAX_DEFAULT, SNHC2D_INITFRAC, MAX_EXTRLIM2D
use simple_string_utils,   only: int2str
use simple_decay_funs,     only: cos_decay, inv_cos_decay, calc_nsampl_fromto, inv_nsampl_decay,&
                                &calc_update_frac, calc_update_frac_dyn, extremal_decay, extremal_decay2D
implicit none
private
public :: run_all_decay_funs_tests

integer, parameter :: MAXITS      = 100
real,    parameter :: BOUNDS(2)   = [0.05, 1.0]
real,    parameter :: EPS         = 1.0e-5

contains

    subroutine run_all_decay_funs_tests()
        write(*,'(A)') '**** running all decay schedule tests ****'
        call test_cos_decays()
        call test_calc_nsampl_fromto()
        call test_inv_nsampl_decay()
        call test_update_frac()
        call test_extremal_decays()
    end subroutine run_all_decay_funs_tests

    !---------------- cosine decays ----------------

    ! cos_decay runs from the upper bound at i = 0 to the lower bound at i = maxits, inv_cos_decay the
    ! other way; their sum is constant (mirror images) and both are monotone
    subroutine test_cos_decays()
        real    :: prev, cur
        integer :: i
        logical :: l_mono_dec, l_mono_inc
        write(*,'(A)') 'test_cos_decays'
        call assert_real(1.0,      cos_decay(0,      MAXITS, BOUNDS), EPS, 'cos_decay starts at the upper bound')
        call assert_real(0.860876, cos_decay(25,     MAXITS, BOUNDS), EPS, 'cos_decay at a quarter')
        call assert_real(0.525,    cos_decay(50,     MAXITS, BOUNDS), EPS, 'cos_decay at the midpoint is the mean of the bounds')
        call assert_real(0.189124, cos_decay(75,     MAXITS, BOUNDS), EPS, 'cos_decay at three quarters')
        call assert_real(0.05,     cos_decay(MAXITS, MAXITS, BOUNDS), EPS, 'cos_decay ends at the lower bound')
        call assert_real(0.05,     inv_cos_decay(0,      MAXITS, BOUNDS), EPS, 'inv_cos_decay starts at the lower bound')
        call assert_real(0.189124, inv_cos_decay(25,     MAXITS, BOUNDS), EPS, 'inv_cos_decay at a quarter')
        call assert_real(0.525,    inv_cos_decay(50,     MAXITS, BOUNDS), EPS, 'inv_cos_decay at the midpoint')
        call assert_real(1.0,      inv_cos_decay(MAXITS, MAXITS, BOUNDS), EPS, 'inv_cos_decay ends at the upper bound')
        l_mono_dec = .true.
        l_mono_inc = .true.
        prev = cos_decay(0, MAXITS, BOUNDS)
        do i = 1,MAXITS
            cur = cos_decay(i, MAXITS, BOUNDS)
            if( cur > prev ) l_mono_dec = .false.
            prev = cur
        end do
        prev = inv_cos_decay(0, MAXITS, BOUNDS)
        do i = 1,MAXITS
            cur = inv_cos_decay(i, MAXITS, BOUNDS)
            if( cur < prev ) l_mono_inc = .false.
            prev = cur
        end do
        call assert_true(l_mono_dec, 'cos_decay is monotone non-increasing over [0,maxits]')
        call assert_true(l_mono_inc, 'inv_cos_decay is monotone non-decreasing over [0,maxits]')
        do i = 0,MAXITS,10
            call assert_real(BOUNDS(1) + BOUNDS(2), cos_decay(i, MAXITS, BOUNDS) + inv_cos_decay(i, MAXITS, BOUNDS), EPS,&
                &'cos_decay and inv_cos_decay are mirror images (i = '//int2str(i)//')')
        end do
    end subroutine test_cos_decays

    !---------------- particle sampling ----------------

    ! small sets (half the particles below the minimum) sample from a quarter to all of the particles;
    ! large sets go from 1/20 of the upper limit (at least the minimum) to half the particles capped
    ! at the maximum
    subroutine test_calc_nsampl_fromto()
        integer :: fromto(2)
        write(*,'(A)') 'test_calc_nsampl_fromto'
        fromto = calc_nsampl_fromto(5000, NSAMPLE_MINMAX_DEFAULT)
        call assert_int(1250,  fromto(1), 'small set: lower limit is a quarter of the particles')
        call assert_int(5000,  fromto(2), 'small set: upper limit is all the particles')
        fromto = calc_nsampl_fromto(12000, NSAMPLE_MINMAX_DEFAULT)
        call assert_int(3000,  fromto(1), 'below-minimum half: lower limit is a quarter')
        call assert_int(12000, fromto(2), 'below-minimum half: upper limit is everything')
        fromto = calc_nsampl_fromto(36000, NSAMPLE_MINMAX_DEFAULT)
        call assert_int(10000, fromto(1), '36000 particles: lower limit is the minimum (900 < 10000)')
        call assert_int(18000, fromto(2), '36000 particles: upper limit is half the particles')
        fromto = calc_nsampl_fromto(80000, NSAMPLE_MINMAX_DEFAULT)
        call assert_int(10000, fromto(1), '80000 particles: lower limit is the minimum')
        call assert_int(25000, fromto(2), '80000 particles: upper limit is capped at the maximum')
        fromto = calc_nsampl_fromto(1000000, NSAMPLE_MINMAX_DEFAULT)
        call assert_int(10000, fromto(1), '1M particles: lower limit is the minimum')
        call assert_int(25000, fromto(2), '1M particles: upper limit is the maximum')
        fromto = calc_nsampl_fromto(80000, [20000, 30000])
        call assert_int(20000, fromto(1), 'custom limits: lower is max(30000/20, 20000)')
        call assert_int(30000, fromto(2), 'custom limits: upper is min(40000, 30000)')
    end subroutine test_calc_nsampl_fromto

    ! the inverse cosine schedule on the sampling window: from the lower limit at it = 0 to the upper
    ! limit at it = maxits, monotone, and held at the upper limit past maxits
    subroutine test_inv_nsampl_decay()
        integer, parameter :: NPTCLS = 36000, MITS = 20
        integer :: it, prev, cur
        logical :: l_mono
        write(*,'(A)') 'test_inv_nsampl_decay'
        call assert_int(10000, inv_nsampl_decay(0,  MITS, NPTCLS, NSAMPLE_MINMAX_DEFAULT), 'it = 0 samples the lower limit')
        call assert_int(10049, inv_nsampl_decay(1,  MITS, NPTCLS, NSAMPLE_MINMAX_DEFAULT), 'it = 1 has barely moved off the lower limit')
        call assert_int(14000, inv_nsampl_decay(10, MITS, NPTCLS, NSAMPLE_MINMAX_DEFAULT), 'it = maxits/2 samples the midpoint')
        call assert_int(18000, inv_nsampl_decay(20, MITS, NPTCLS, NSAMPLE_MINMAX_DEFAULT), 'it = maxits samples the upper limit')
        call assert_int(18000, inv_nsampl_decay(25, MITS, NPTCLS, NSAMPLE_MINMAX_DEFAULT), 'past maxits the upper limit is held')
        l_mono = .true.
        prev   = inv_nsampl_decay(0, MITS, NPTCLS, NSAMPLE_MINMAX_DEFAULT)
        do it = 1,MITS
            cur = inv_nsampl_decay(it, MITS, NPTCLS, NSAMPLE_MINMAX_DEFAULT)
            if( cur < prev ) l_mono = .false.
            prev = cur
        end do
        call assert_true(l_mono, 'the sample size never shrinks along the schedule')
    end subroutine test_inv_nsampl_decay

    ! the static fraction is the maximum sample size (scaled by the number of states, capped at the
    ! particle count) over the particle count; the dynamic one follows inv_nsampl_decay
    subroutine test_update_frac()
        write(*,'(A)') 'test_update_frac'
        call assert_real(25000./36000., calc_update_frac(36000,  1, NSAMPLE_MINMAX_DEFAULT), EPS, 'static fraction: maximum over particles')
        call assert_real(1.0,           calc_update_frac(12000,  1, NSAMPLE_MINMAX_DEFAULT), EPS, 'static fraction: fewer particles than the maximum gives 1')
        call assert_real(0.5,           calc_update_frac(100000, 2, NSAMPLE_MINMAX_DEFAULT), EPS, 'static fraction: two states double the maximum')
        call assert_real(1.0,           calc_update_frac(20000,  3, NSAMPLE_MINMAX_DEFAULT), EPS, 'static fraction: three states saturate at 1')
        call assert_real(10000./36000., calc_update_frac_dyn(36000, 1, NSAMPLE_MINMAX_DEFAULT, 0,  20), EPS, 'dynamic fraction at it = 0 is the lower limit')
        call assert_real(14000./36000., calc_update_frac_dyn(36000, 1, NSAMPLE_MINMAX_DEFAULT, 10, 20), EPS, 'dynamic fraction at the midpoint')
        call assert_real(0.5,           calc_update_frac_dyn(36000, 1, NSAMPLE_MINMAX_DEFAULT, 20, 20), EPS, 'dynamic fraction at maxits is the upper limit')
    end subroutine test_update_frac

    !---------------- extremal searches ----------------

    ! extremal_decay anneals the neighbourhood fraction from 0.95 at it = 1 to 0.3 at it = maxits and holds
    ! 0.3 afterwards; extremal_decay2D is the factorial decay of the 2D SNHC, clamped to [0, SNHC2D_INITFRAC]
    subroutine test_extremal_decays()
        real    :: prev, cur
        integer :: it
        logical :: l_mono
        write(*,'(A)') 'test_extremal_decays'
        call assert_real(0.95,     extremal_decay(1,  10), EPS, 'extremal_decay starts at 0.95')
        call assert_real(0.7875,   extremal_decay(4,  10), EPS, 'extremal_decay at it = 4 of 10')
        call assert_real(0.681436, extremal_decay(5,  10), EPS, 'extremal_decay at it = 5 of 10')
        call assert_real(0.3,      extremal_decay(10, 10), EPS, 'extremal_decay ends at 0.3')
        call assert_real(0.3,      extremal_decay(11, 10), EPS, 'extremal_decay holds 0.3 past maxits')
        call assert_real(0.95,     extremal_decay(1,  1),  EPS, 'a single iteration starts at 0.95 (no division by zero)')
        call assert_real(0.3,      extremal_decay(2,  1),  EPS, 'a single iteration holds 0.3 from it = 2')
        l_mono = .true.
        prev   = extremal_decay(1, 10)
        do it = 2,12
            cur = extremal_decay(it, 10)
            if( cur > prev ) l_mono = .false.
            prev = cur
        end do
        call assert_true(l_mono, 'extremal_decay is monotone non-increasing')
        call assert_real(SNHC2D_INITFRAC, extremal_decay2D(1,  30), EPS, 'extremal_decay2D is clamped to the initial fraction early on')
        call assert_real(SNHC2D_INITFRAC, extremal_decay2D(4,  30), EPS, 'extremal_decay2D still clamped at the greedy start')
        call assert_real(0.256,           extremal_decay2D(10, 30), EPS, 'extremal_decay2D at it = 10 of 30 (power 3)')
        call assert_real(0.027488,        extremal_decay2D(30, 30), EPS, 'extremal_decay2D at the extremal limit (power 13)')
        l_mono = .true.
        prev   = extremal_decay2D(1, 30)
        do it = 2,30
            cur = extremal_decay2D(it, 30)
            if( cur > prev ) l_mono = .false.
            prev = cur
        end do
        call assert_true(l_mono, 'extremal_decay2D is monotone non-increasing')
        call assert_true(extremal_decay2D(60, 30) >= 0., 'extremal_decay2D never goes negative')
    end subroutine test_extremal_decays

end module simple_decay_funs_tester
