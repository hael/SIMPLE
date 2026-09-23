!@descr: unit test routines for the Kaiser-Bessel interpolation kernel (simple_kbinterpol)
! The apodisation function and its polynomial fast form, the value/derivative pair, the separable 2D/3D
! stencils (exact, fast, fast with gradient) and the instrument function, pinned on closed-form values of
! the standard KBWINSZ=1.5, KBALPHA=2 kernel.
module simple_kbinterpol_tester
use simple_test_utils    ! assertions etc.
use simple_defs          ! sp, dp, PI, KBWINSZ, KBALPHA, KB_BETA_KB15_A2
use simple_string_utils, only: int2str
use simple_kbinterpol,   only: kbinterpol, kb_windim, apod_device, apod_fast_device, apod_kb15_a2
implicit none
private
public :: run_all_kbinterpol_tests

! the standard kernel: W = 3, beta = pi sqrt((W/alpha)^2 (alpha - 1/2)^2 - 0.8) = 6.48607653
real(dp), parameter :: BETA_REF    = 6.48607653051083_dp
real(dp), parameter :: I0_BETA_REF = 104.94085488702105_dp  ! I0(beta), power series to 1e-18
real(dp), parameter :: SINHC_REF   = 50.5654844798751_dp    ! sinh(beta)/beta, the instrument function at 0
integer,  parameter :: NCOEFF      = 15                     ! degree-14 polynomial of apod_fast
integer,  parameter :: IWINSZ      = 1                      ! ceiling(KBWINSZ - 0.5)
integer,  parameter :: WDIM        = 2 * IWINSZ + 1
real,     parameter :: FD_STEP     = 1.0e-3

contains

    subroutine run_all_kbinterpol_tests()
        write(*,'(A)') '**** running all Kaiser-Bessel kernel tests ****'
        call test_window_geometry()
        call test_apod_closed_form()
        call test_apod_fast_polynomial()
        call test_apod_fast_value_deriv()
        call test_device_forms()
        call test_apod_mat_2d()
        call test_apod_mat_3d()
        call test_fast_stencils()
        call test_apod_mat_3d_fast_grad()
        call test_instr()
    end subroutine run_all_kbinterpol_tests

    function standard_kernel() result( kb )
        type(kbinterpol) :: kb
        kb = kbinterpol(KBWINSZ, KBALPHA)
    end function standard_kernel

    ! I0 by its power series in double precision
    pure real(dp) function bessel_i0( x )
        real(dp), intent(in) :: x
        real(dp) :: term
        integer  :: k
        bessel_i0 = 1.0_dp
        term      = 1.0_dp
        do k = 1,60
            term      = term * (x / 2.0_dp)**2 / real(k, dp)**2
            bessel_i0 = bessel_i0 + term
            if( term < 1.0e-18_dp ) exit
        enddo
    end function bessel_i0

    ! the closed form apod(x) = I0(beta sqrt(1 - (2x/W)^2)) / W, in double precision
    pure real(dp) function apod_ref( x )
        real, intent(in) :: x
        real(dp) :: w, u
        w = 2.0_dp * real(KBWINSZ, dp)
        u = 1.0_dp - (2.0_dp * real(x, dp) / w)**2
        if( abs(x) > KBWINSZ )then
            apod_ref = 0.0_dp
        else
            apod_ref = bessel_i0(BETA_REF * sqrt(max(0.0_dp, u))) / w
        endif
    end function apod_ref

    !---------------- window geometry ----------------

    subroutine test_window_geometry()
        type(kbinterpol) :: kb
        write(*,'(A)') 'test_window_geometry'
        kb = standard_kernel()
        call assert_real(KBWINSZ, kb%get_winsz(), 1.e-6, 'get_winsz returns the half-width')
        call assert_real(KBALPHA, kb%get_alpha(), 1.e-6, 'get_alpha returns the oversampling factor')
        call assert_int(WDIM, kb%get_wdim(),             'get_wdim: 2 ceiling(W/2 - 1/2) + 1 = 3 for W/2 = 1.5')
        call assert_int(WDIM, kb_windim(KBWINSZ),        'kb_windim agrees with get_wdim')
        call assert_int(5,    kb_windim(2.5),            'kb_windim: half-width 2.5 gives a 5-point window')
        call assert_real(real(BETA_REF), KB_BETA_KB15_A2, 1.e-5, 'KB_BETA_KB15_A2 is the closed-form beta of the standard kernel')
    end subroutine test_window_geometry

    !---------------- the apodisation function ----------------

    subroutine test_apod_closed_form()
        type(kbinterpol) :: kb
        real    :: x, prev
        integer :: i
        logical :: symmetric, decreasing
        write(*,'(A)') 'test_apod_closed_form'
        kb = standard_kernel()
        call assert_real(real(I0_BETA_REF / 3.0_dp), kb%apod(0.),      2.e-4, 'apod(0) = I0(beta)/W')
        call assert_real(1.0 / 3.0,                  kb%apod(KBWINSZ), 1.e-6, 'apod at the half-width = 1/W (I0(0) = 1)')
        call assert_real(0.0,                        kb%apod(1.5001),  0.0,   'apod is exactly zero beyond the half-width')
        call assert_real(0.0,                        kb%apod(-7.),     0.0,   'apod is exactly zero far outside')
        call assert_real(real(apod_ref(0.5)),  kb%apod(0.5),  1.e-4, 'apod(0.5) against the double-precision closed form')
        call assert_real(real(apod_ref(1.0)),  kb%apod(1.0),  1.e-4, 'apod(1.0) against the double-precision closed form')
        call assert_real(real(apod_ref(1.25)), kb%apod(1.25), 1.e-4, 'apod(1.25) against the double-precision closed form')
        symmetric  = .true.
        decreasing = .true.
        prev = kb%apod(0.)
        do i = 1,300
            x = KBWINSZ * real(i) / 300.
            if( kb%apod(x) /= kb%apod(-x) ) symmetric = .false.
            if( kb%apod(x) > prev ) decreasing = .false.
            prev = kb%apod(x)
        enddo
        call assert_true(symmetric,  'apod is even')
        call assert_true(decreasing, 'apod decreases monotonically from 0 to the half-width')
    end subroutine test_apod_closed_form

    ! apod_fast is the I0 power series truncated at degree 14 in u = 1 - (2x/W)^2
    subroutine test_apod_fast_polynomial()
        type(kbinterpol) :: kb, other
        real     :: x, w_exact, w_fast, max_abs, max_rel
        real(dp) :: coeff
        integer  :: i, k
        logical  :: fallback_exact
        write(*,'(A)') 'test_apod_fast_polynomial'
        kb = standard_kernel()
        ! the coefficients are (beta^2/4)^k / (k!)^2
        coeff = 1.0_dp
        do k = 1,NCOEFF-1
            coeff = coeff * (BETA_REF**2 / 4.0_dp) / real(k, dp)**2
        enddo
        call assert_true(abs(coeff - 2.6658846e-8_dp) < 1.e-14_dp, 'apod_fast: the degree-14 coefficient is (beta^2/4)^14 / (14!)^2')
        max_abs = 0.
        max_rel = 0.
        do i = -2000,2000
            x       = KBWINSZ * real(i) / 2000.
            w_exact = kb%apod(x)
            w_fast  = kb%apod_fast(x)
            max_abs = max(max_abs, abs(w_exact - w_fast))
            max_rel = max(max_rel, abs(w_exact - w_fast) / max(abs(w_exact), TINY))
        enddo
        call assert_true(max_abs < 1.e-4, 'apod_fast: absolute deviation from apod below 1e-4 over the window')
        call assert_true(max_rel < 1.e-5, 'apod_fast: relative deviation from apod below 1e-5 over the window')
        call assert_real(0.0, kb%apod_fast(1.5001), 0.0, 'apod_fast is exactly zero beyond the half-width')
        call assert_real(kb%apod_fast(0.7), kb%apod_fast(-0.7), 0.0, 'apod_fast is even')
        ! any other kernel falls back to apod exactly
        other = kbinterpol(2.0, 2.0)
        fallback_exact = .true.
        do i = -40,40
            x = 2.0 * real(i) / 40.
            if( other%apod_fast(x) /= other%apod(x) ) fallback_exact = .false.
        enddo
        call assert_true(fallback_exact, 'apod_fast falls back to apod for a non-standard kernel')
    end subroutine test_apod_fast_polynomial

    subroutine test_apod_fast_value_deriv()
        type(kbinterpol) :: kb
        real    :: x, value, deriv, vplus, vminus, fd, max_err, max_val
        integer :: i
        logical :: value_matches, antisymmetric
        write(*,'(A)') 'test_apod_fast_value_deriv'
        kb = standard_kernel()
        value_matches = .true.
        antisymmetric = .true.
        max_err = 0.
        max_val = 0.
        do i = -140,140
            x = 1.4 * real(i) / 140.
            call kb%apod_fast_value_deriv(x, value, deriv)
            if( value /= kb%apod_fast(x) ) value_matches = .false.
            call kb%apod_fast_value_deriv(-x, vplus, vminus) ! vminus is the derivative at -x here
            if( abs(vminus + deriv) > 1.e-5 * max(1., abs(deriv)) ) antisymmetric = .false.
            vplus  = kb%apod_fast(x + FD_STEP)
            vminus = kb%apod_fast(x - FD_STEP)
            fd = (vplus - vminus) / (2. * FD_STEP)
            max_err = max(max_err, abs(fd - deriv))
            max_val = max(max_val, abs(deriv))
        enddo
        call assert_true(value_matches, 'apod_fast_value_deriv: the value is apod_fast bit for bit')
        call assert_true(antisymmetric,  'apod_fast_value_deriv: the derivative is odd')
        call assert_true(max_val > 30.,  'apod_fast_value_deriv: the kernel slope reaches tens per pixel')
        call assert_true(max_err < 2.e-3 * max_val, 'apod_fast_value_deriv: derivative matches a central difference of apod_fast')
        call kb%apod_fast_value_deriv(0., value, deriv)
        call assert_real(0.0, deriv, 0.0, 'apod_fast_value_deriv: zero slope at the centre')
        call kb%apod_fast_value_deriv(2., value, deriv)
        call assert_real(0.0, value, 0.0, 'apod_fast_value_deriv: zero value outside the window')
        call assert_real(0.0, deriv, 0.0, 'apod_fast_value_deriv: zero slope outside the window')
    end subroutine test_apod_fast_value_deriv

    ! the non-polymorphic device forms must stay identical to the type-bound ones
    subroutine test_device_forms()
        type(kbinterpol) :: kb
        real    :: x
        integer :: i
        logical :: same_apod, same_fast, same_kb15
        write(*,'(A)') 'test_device_forms'
        kb = standard_kernel()
        same_apod = .true.
        same_fast = .true.
        same_kb15 = .true.
        do i = -200,200
            x = 2.0 * real(i) / 200.
            if( apod_device(kb, x)      /= kb%apod(x)      ) same_apod = .false.
            if( apod_fast_device(kb, x) /= kb%apod_fast(x) ) same_fast = .false.
            if( apod_kb15_a2(x)         /= kb%apod_fast(x) ) same_kb15 = .false.
        enddo
        call assert_true(same_apod, 'apod_device equals apod bit for bit')
        call assert_true(same_fast, 'apod_fast_device equals apod_fast bit for bit')
        call assert_true(same_kb15, 'apod_kb15_a2 equals apod_fast of the standard kernel bit for bit')
    end subroutine test_device_forms

    !---------------- separable stencils ----------------

    ! the reference the old test computed by hand: product of the axis kernels, normalised once
    subroutine reference_2d( kb, loc, kbw )
        type(kbinterpol), intent(in)  :: kb
        real,             intent(in)  :: loc(2)
        real,             intent(out) :: kbw(WDIM,WDIM)
        integer :: win_lo(2), i
        real    :: d(2)
        win_lo = nint(loc) - IWINSZ
        kbw    = 1.
        do i = 1,WDIM
            d = real(win_lo + i - 1) - loc
            kbw(i,:) = kbw(i,:) * kb%apod(d(1))
            kbw(:,i) = kbw(:,i) * kb%apod(d(2))
        enddo
        kbw = kbw / sum(kbw)
    end subroutine reference_2d

    subroutine reference_3d( kb, loc, kbw )
        type(kbinterpol), intent(in)  :: kb
        real,             intent(in)  :: loc(3)
        real,             intent(out) :: kbw(WDIM,WDIM,WDIM)
        integer :: win_lo(3), i
        real    :: d(3)
        win_lo = nint(loc) - IWINSZ
        kbw    = 1.
        do i = 1,WDIM
            d = real(win_lo + i - 1) - loc
            kbw(i,:,:) = kbw(i,:,:) * kb%apod(d(1))
            kbw(:,i,:) = kbw(:,i,:) * kb%apod(d(2))
            kbw(:,:,i) = kbw(:,:,i) * kb%apod(d(3))
        enddo
        kbw = kbw / sum(kbw)
    end subroutine reference_3d

    subroutine test_apod_mat_2d()
        real, parameter :: LOCS(2,5) = reshape([0.0, 0.0,  0.3, -0.7,  2.5, 1.25,  -3.1, 4.9,  0.49, -0.49], [2,5])
        type(kbinterpol) :: kb
        real    :: kbw(WDIM,WDIM), ref(WDIM,WDIM)
        integer :: iloc
        write(*,'(A)') 'test_apod_mat_2d'
        kb = standard_kernel()
        do iloc = 1,size(LOCS, 2)
            call kb%apod_mat_2d(LOCS(:,iloc), IWINSZ, WDIM, kbw)
            call reference_2d(kb, LOCS(:,iloc), ref)
            call assert_true(maxval(abs(kbw - ref)) < 1.e-6, 'apod_mat_2d equals the normalised outer product at location '//int2str(iloc))
            call assert_real(1.0, sum(kbw), 1.e-6, 'apod_mat_2d sums to one at location '//int2str(iloc))
            call assert_true(all(kbw >= 0.), 'apod_mat_2d is non-negative at location '//int2str(iloc))
        enddo
        ! on a grid point the centre pixel dominates and the stencil is symmetric
        call kb%apod_mat_2d([0., 0.], IWINSZ, WDIM, kbw)
        call assert_true(kbw(2,2) > maxval(kbw) - 1.e-7, 'apod_mat_2d on a grid point peaks at the centre')
        call assert_real(kbw(1,2), kbw(3,2), 1.e-7, 'apod_mat_2d on a grid point is symmetric in x')
        call assert_real(kbw(2,1), kbw(2,3), 1.e-7, 'apod_mat_2d on a grid point is symmetric in y')
    end subroutine test_apod_mat_2d

    subroutine test_apod_mat_3d()
        real, parameter :: LOCS(3,4) = reshape([0.0, 0.0, 0.0,  0.3, -0.7, 1.2,  2.5, 1.25, -0.5,  -3.1, 4.9, 0.49], [3,4])
        type(kbinterpol) :: kb
        real    :: kbw(WDIM,WDIM,WDIM), ref(WDIM,WDIM,WDIM)
        integer :: iloc
        write(*,'(A)') 'test_apod_mat_3d'
        kb = standard_kernel()
        do iloc = 1,size(LOCS, 2)
            call kb%apod_mat_3d(LOCS(:,iloc), IWINSZ, WDIM, kbw)
            call reference_3d(kb, LOCS(:,iloc), ref)
            call assert_true(maxval(abs(kbw - ref)) < 1.e-6, 'apod_mat_3d equals the normalised outer product at location '//int2str(iloc))
            call assert_real(1.0, sum(kbw), 1.e-6, 'apod_mat_3d sums to one at location '//int2str(iloc))
        enddo
        call kb%apod_mat_3d([0., 0., 0.], IWINSZ, WDIM, kbw)
        call assert_true(kbw(2,2,2) > maxval(kbw) - 1.e-7, 'apod_mat_3d on a grid point peaks at the centre')
        call assert_real(kbw(1,2,2), kbw(3,2,2), 1.e-7, 'apod_mat_3d on a grid point is symmetric in x')
        call assert_real(kbw(2,2,1), kbw(2,2,3), 1.e-7, 'apod_mat_3d on a grid point is symmetric in z')
    end subroutine test_apod_mat_3d

    subroutine test_fast_stencils()
        real, parameter :: LOCS2(2,3) = reshape([0.3, -0.7,  2.5, 1.25,  -3.1, 4.9], [2,3])
        real, parameter :: LOCS3(3,3) = reshape([0.3, -0.7, 1.2,  2.5, 1.25, -0.5,  -3.1, 4.9, 0.49], [3,3])
        type(kbinterpol) :: kb
        real    :: k2(WDIM,WDIM), f2(WDIM,WDIM), k3(WDIM,WDIM,WDIM), f3(WDIM,WDIM,WDIM)
        integer :: iloc
        write(*,'(A)') 'test_fast_stencils'
        kb = standard_kernel()
        do iloc = 1,size(LOCS2, 2)
            call kb%apod_mat_2d(LOCS2(:,iloc), IWINSZ, WDIM, k2)
            call kb%apod_mat_2d_fast(LOCS2(:,iloc), IWINSZ, WDIM, f2)
            call assert_true(maxval(abs(k2 - f2)) < 1.e-5, 'apod_mat_2d_fast matches apod_mat_2d at location '//int2str(iloc))
            call assert_real(1.0, sum(f2), 1.e-5, 'apod_mat_2d_fast sums to one without the final renormalisation, location '//int2str(iloc))
        enddo
        do iloc = 1,size(LOCS3, 2)
            call kb%apod_mat_3d(LOCS3(:,iloc), IWINSZ, WDIM, k3)
            call kb%apod_mat_3d_fast(LOCS3(:,iloc), IWINSZ, WDIM, f3)
            call assert_true(maxval(abs(k3 - f3)) < 1.e-5, 'apod_mat_3d_fast matches apod_mat_3d at location '//int2str(iloc))
            call assert_real(1.0, sum(f3), 1.e-5, 'apod_mat_3d_fast sums to one without the final renormalisation, location '//int2str(iloc))
        enddo
    end subroutine test_fast_stencils

    ! the gradient of the fixed-cell stencil against central differences of apod_mat_3d_fast
    subroutine test_apod_mat_3d_fast_grad()
        real, parameter :: LOC(3) = [0.3, -0.7, 1.2]
        type(kbinterpol) :: kb
        real    :: kbw(WDIM,WDIM,WDIM), fast(WDIM,WDIM,WDIM), dw(WDIM,WDIM,WDIM,3)
        real    :: plus(WDIM,WDIM,WDIM), minus(WDIM,WDIM,WDIM), fd(WDIM,WDIM,WDIM), shifted(3), margin(3)
        integer :: i0(3), axis
        write(*,'(A)') 'test_apod_mat_3d_fast_grad'
        kb = standard_kernel()
        call kb%apod_mat_3d_fast_grad(LOC, IWINSZ, WDIM, i0, margin, kbw, dw)
        call assert_true(all(i0 == nint(LOC) - IWINSZ), 'apod_mat_3d_fast_grad: stencil corner is nint(loc) - iwinsz')
        call assert_true(maxval(abs(margin - (0.5 - abs(LOC - real(nint(LOC)))))) < 1.e-6, &
            &'apod_mat_3d_fast_grad: switch margin is the distance to the next nint switch')
        call kb%apod_mat_3d_fast(LOC, IWINSZ, WDIM, fast)
        call assert_true(maxval(abs(kbw - fast)) < 1.e-6, 'apod_mat_3d_fast_grad: the stencil is apod_mat_3d_fast')
        do axis = 1,3
            call assert_real(0.0, sum(dw(:,:,:,axis)), 1.e-5, 'apod_mat_3d_fast_grad: the gradient of a normalised stencil sums to zero, axis '//int2str(axis))
            shifted = LOC
            shifted(axis) = LOC(axis) + FD_STEP
            call kb%apod_mat_3d_fast(shifted, IWINSZ, WDIM, plus)
            shifted(axis) = LOC(axis) - FD_STEP
            call kb%apod_mat_3d_fast(shifted, IWINSZ, WDIM, minus)
            fd = (plus - minus) / (2. * FD_STEP)
            call assert_true(maxval(abs(fd - dw(:,:,:,axis))) < 1.e-3, 'apod_mat_3d_fast_grad: gradient matches central differences, axis '//int2str(axis))
            call assert_true(maxval(abs(dw(:,:,:,axis))) > 0.1, 'apod_mat_3d_fast_grad: the gradient is not trivially small, axis '//int2str(axis))
        enddo
        ! on a half-integer the margin is zero
        call kb%apod_mat_3d_fast_grad([0.5, 1.0, -2.5], IWINSZ, WDIM, i0, margin, kbw, dw)
        call assert_real(0.0, margin(1), 1.e-7, 'apod_mat_3d_fast_grad: zero margin at a half-integer coordinate')
        call assert_real(0.5, margin(2), 1.e-7, 'apod_mat_3d_fast_grad: half margin on a grid point')
        call assert_real(0.0, margin(3), 1.e-7, 'apod_mat_3d_fast_grad: zero margin at a negative half-integer')
    end subroutine test_apod_mat_3d_fast_grad

    !---------------- the instrument function ----------------

    subroutine test_instr()
        type(kbinterpol) :: kb
        real    :: x, prev
        integer :: i
        logical :: decreasing, symmetric
        write(*,'(A)') 'test_instr'
        kb = standard_kernel()
        call assert_real(real(SINHC_REF), kb%instr(0.), 1.e-3, 'instr(0) = sinh(beta)/beta')
        call assert_real(1.0, kb%instr(0.7),  1.e-6, 'instr is one beyond beta/(pi W) = 0.688')
        call assert_real(1.0, kb%instr(-2.0), 1.e-6, 'instr is one far outside')
        decreasing = .true.
        symmetric  = .true.
        prev = kb%instr(0.)
        do i = 1,68
            x = 0.01 * real(i)
            if( kb%instr(x) > prev + 1.e-6 ) decreasing = .false.
            if( abs(kb%instr(x) - kb%instr(-x)) > 1.e-6 ) symmetric = .false.
            prev = kb%instr(x)
        enddo
        call assert_true(decreasing, 'instr decreases from 0 to the threshold')
        call assert_true(symmetric,  'instr is even')
        call assert_true(kb%instr(0.68) > 1.0, 'instr is above one just inside the threshold')
    end subroutine test_instr

end module simple_kbinterpol_tester
