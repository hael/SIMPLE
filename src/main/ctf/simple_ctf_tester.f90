!@descr: unit test routines for the contrast transfer function (simple_ctf)
! The CTFFIND4 wavelength, defocus and phase conventions, the additive phase-shift policy (canonicalisation
! modulo 2pi, never pi), the extrema count and the zeros, the fast image-layer kernel (ctf2img, ft2img placement) and the
! H*y / |H|*y restoration contracts of gen_fplane4rec; the closed-form values are evaluated in double
! precision from the same formulas (Rohou & Grigorieff 2015).
module simple_ctf_tester
use simple_test_utils          ! assertions etc.
use simple_defs                ! sp, dp, PI, PIO2, TWOPI, OSMPL_PAD_FAC
use simple_type_defs,          only: ctfparams, ctfvars, fplane_type, canonical_phshift, CTFFLAG_YES, CTFFLAG_FLIP
use simple_string_utils,       only: int2str
use simple_image,              only: image
use simple_ctf,                only: ctf
use simple_memoize_ft_maps,    only: memoize_ft_maps
implicit none
private
public :: run_all_ctf_tests

integer, parameter :: LDIM(3) = [64,64,1]
integer, parameter :: H = 5, K = 3
real,    parameter :: SMPD = 1.0, DFX = 2.0, DFY = 2.3, ANGAST = 27.
real,    parameter :: KV = 300., CS = 2.0, AC = 0.1
real,    parameter :: TOL = 2.e-6   ! Fortran-vs-Fortran comparisons
! 3.0 rad ~ 172 deg exercises the near-pi regime typical of a laser phase plate,
! where an erroneous fold at pi would flip the sign of the transfer function
real,    parameter :: PHASES(4) = [0., PI/4., PIO2, 3.0]
! closed-form references (double precision, ctf_ref.py of the stats batch, 2026-09-23)
real,    parameter :: WL_300KV    = 0.0196808236   ! Angstrom, 12.26/sqrt(V + 0.9784 V^2/1e6)
real,    parameter :: WL_200KV    = 0.0250707883
real,    parameter :: AC_CONST    = 0.1001674212   ! atan(0.1/sqrt(1-0.01))
real,    parameter :: CTF_5_3_0   = -0.8021261     ! (h,k) = (5,3), phase shift 0, chi = 10.3 rad
real,    parameter :: CTF_5_3_Q   = -0.9894409     ! (h,k) = (5,3), phase shift pi/4
real,    parameter :: CTF_10_M7_0 = -0.2373072     ! (h,k) = (10,-7), phase shift 0, chi = 49.9 rad
real,    parameter :: CTF_10_M7_Q =  0.5191065     ! (h,k) = (10,-7), phase shift pi/4
real,    parameter :: CTF_20_12_0 = -0.0059576     ! (h,k) = (20,12), phase shift 0, chi = 160 rad
! the first three zeros along the astigmatism axis (df = 2.0 um): 20.2, 14.1 and 11.5 Angstrom
real,    parameter :: S2_ZEROS(3) = [0.00246071, 0.00500494, 0.00755167]

contains

    subroutine run_all_ctf_tests()
        write(*,'(A)') '**** running all CTF tests ****'
        call test_constructor_and_conventions()
        call test_closed_form_values()
        call test_phase_shift_policy()
        call test_extrema()
        call test_fast_kernel_and_ft2img()
        call test_restoration_contracts()
    end subroutine run_all_ctf_tests

    function make_tfun() result( tfun )
        type(ctf) :: tfun
        tfun = ctf(SMPD, KV, CS, AC)
        call tfun%init(DFX, DFY, ANGAST)
    end function make_tfun

    pure real function spafreqsq( hh, kk )
        integer, intent(in) :: hh, kk
        spafreqsq = (real(hh)/real(LDIM(1)))**2 + (real(kk)/real(LDIM(2)))**2
    end function spafreqsq

    pure real function angle( hh, kk )
        integer, intent(in) :: hh, kk
        angle = atan2(real(kk), real(hh))
    end function angle

    ! the CTF at an arbitrary (not yet canonical) phase shift: the production path canonicalises once
    ! per particle (get_ctfvars) and evaluates with eval_canonical in the pixel loops
    real function ctf_at( tfun, s2, ang, phshift )
        type(ctf), intent(in) :: tfun
        real,      intent(in) :: s2, ang, phshift
        ctf_at = tfun%eval_canonical(s2, ang, canonical_phshift(phshift))
    end function ctf_at

    ! sin(pi wl s^2 (df - 0.5 wl^2 s^2 Cs) + phshift + atan(ac/sqrt(1-ac^2))) in double precision, with
    ! df, wl and Cs in pixel units as the ctf object stores them
    pure real function closed_form( s2, df_um, phshift )
        real, intent(in) :: s2, df_um, phshift
        real(dp) :: wl, cs_pix, df_pix, chi
        wl     = real(WL_300KV, dp) / real(SMPD, dp)
        cs_pix = real(CS, dp) * 1.0d7 / real(SMPD, dp)
        df_pix = real(df_um, dp) * 1.0d4 / real(SMPD, dp)
        chi    = acos(-1.d0) * wl * real(s2,dp) * (df_pix - 0.5d0 * wl * wl * real(s2,dp) * cs_pix) + real(phshift,dp)
        closed_form = real(sin(chi + real(AC_CONST,dp)))
    end function closed_form

    !---------------- constructor, get_ctfvars, convention ----------------

    subroutine test_constructor_and_conventions()
        type(ctf)     :: tfun
        type(ctfvars) :: vars
        real          :: df_1, df_2, ang_c   ! not dfx/dfy/angast: those would hide the module constants
        write(*,'(A)') 'test_constructor_and_conventions'
        tfun = make_tfun()
        vars = tfun%get_ctfvars(0.)
        call assert_real(WL_300KV/SMPD, vars%wl, 1.e-7, 'wavelength at 300 kV in pixels (relativistic de Broglie)')
        call assert_real(AC_CONST, vars%amp_contr_const, 1.e-6, 'amplitude contrast enters as atan(ac/sqrt(1-ac^2))')
        call assert_real(DFX*1.e4/SMPD, vars%dfx, 1.e-2, 'dfx is converted from microns to pixels')
        call assert_real(DFY*1.e4/SMPD, vars%dfy, 1.e-2, 'dfy is converted from microns to pixels')
        call assert_real(ANGAST*PI/180., vars%angast, 1.e-6, 'astigmatism angle is stored in radians')
        call assert_real(KV, vars%kv, TOL, 'ctfvars carries the voltage')
        call assert_real(SMPD, vars%smpd, TOL, 'ctfvars carries the sampling distance')
        tfun = ctf(SMPD, 200., CS, AC)
        vars = tfun%get_ctfvars(0.)
        call assert_real(WL_200KV/SMPD, vars%wl, 1.e-7, 'wavelength at 200 kV')
        tfun = ctf(2.0, KV, CS, AC)
        call tfun%init(DFX, DFY, ANGAST)
        vars = tfun%get_ctfvars(0.)
        call assert_real(WL_300KV/2.0, vars%wl, 1.e-7, 'wavelength scales with 1/smpd')
        call assert_real(DFX*1.e4/2.0, vars%dfx, 1.e-2, 'defocus in pixels scales with 1/smpd')
        ! apply_convention: negative angles wrap, dfx >= dfy with the angle rotated by 90 degrees
        tfun   = make_tfun()
        df_1    = 1.0
        df_2    = 2.0
        ang_c = -30.
        call tfun%apply_convention(df_1, df_2, ang_c)
        call assert_real(2.0,  df_1,    TOL, 'apply_convention: the larger defocus becomes dfx')
        call assert_real(1.0,  df_2,    TOL, 'apply_convention: the smaller defocus becomes dfy')
        call assert_real(60.0, ang_c, 1.e-4, 'apply_convention: -30 -> 150 -> +90 for the swap -> 240 -> 60')
        df_1    = 2.5
        df_2    = 2.0
        ang_c = 100.
        call tfun%apply_convention(df_1, df_2, ang_c)
        call assert_real(2.5,   df_1,    TOL, 'apply_convention leaves an ordered pair alone')
        call assert_real(2.0,   df_2,    TOL, 'apply_convention leaves an ordered pair alone (dfy)')
        call assert_real(100.0, ang_c, 1.e-4, 'apply_convention leaves an in-range angle alone')
    end subroutine test_constructor_and_conventions

    !---------------- closed-form values ----------------

    ! the tolerances follow the size of the phase argument evaluated in single precision (10, 50 and
    ! 160 radians at the three frequencies)
    subroutine test_closed_form_values()
        type(ctf) :: tfun
        real      :: s2, ang_on, ang_off
        write(*,'(A)') 'test_closed_form_values'
        tfun = make_tfun()
        call assert_real(CTF_5_3_0,   ctf_at(tfun, spafreqsq(5,3),   angle(5,3),   0.),    1.e-5, 'CTF at (5,3), no phase shift')
        call assert_real(CTF_5_3_Q,   ctf_at(tfun, spafreqsq(5,3),   angle(5,3),   PI/4.), 1.e-5, 'CTF at (5,3), pi/4 phase shift')
        call assert_real(CTF_10_M7_0, ctf_at(tfun, spafreqsq(10,-7), angle(10,-7), 0.),    5.e-5, 'CTF at (10,-7), no phase shift')
        call assert_real(CTF_10_M7_Q, ctf_at(tfun, spafreqsq(10,-7), angle(10,-7), PI/4.), 5.e-5, 'CTF at (10,-7), pi/4 phase shift')
        call assert_real(CTF_20_12_0, ctf_at(tfun, spafreqsq(20,12), angle(20,12), 0.),    5.e-4, 'CTF at (20,12), chi of 160 radians')
        ! the effective defocus df = 0.5 (dfx + dfy + cos(2(ang - angast)) (dfx - dfy)): dfx along the
        ! astigmatism axis, dfy perpendicular to it, their mean at 45 degrees
        s2      = spafreqsq(2,0)
        ang_on  = ANGAST*PI/180.
        ang_off = ang_on + PIO2
        call assert_real(closed_form(s2, DFX, 0.), ctf_at(tfun, s2, ang_on, 0.), 2.e-5, &
            &'on the astigmatism axis the effective defocus is dfx')
        call assert_real(closed_form(s2, DFY, 0.), ctf_at(tfun, s2, ang_off, 0.), 2.e-5, &
            &'perpendicular to it the effective defocus is dfy')
        call assert_real(closed_form(s2, 0.5*(DFX+DFY), 0.), ctf_at(tfun, s2, ang_on + PI/4., 0.), 2.e-5, &
            &'at 45 degrees the effective defocus is the mean')
        call assert_real(closed_form(s2, DFX, 1.2), ctf_at(tfun, s2, ang_on, 1.2), 2.e-5, &
            &'an additive phase shift enters the argument as is')
    end subroutine test_closed_form_values

    !---------------- phase shift policy ----------------

    ! Canonicalization is modulo 2*pi, which is an identity for the transfer function.
    ! It must never be modulo pi: a pi offset negates the CTF, so folding at pi would
    ! split a physically homogeneous phase population into two sign-opposite groups.
    subroutine test_phase_shift_policy()
        type(ctf)     :: tfun
        type(ctfvars) :: ctfvals
        real          :: s2, ang
        write(*,'(A)') 'test_phase_shift_policy'
        tfun = make_tfun()
        s2   = spafreqsq(H,K)
        ang  = angle(H,K)
        call assert_real(PI, canonical_phshift(PI), TOL, &
            &'phase canonicalization leaves pi alone; it is a distinct, sign-flipped CTF')
        call assert_real(0., canonical_phshift(TWOPI), TOL, &
            &'phase canonicalization maps 2*pi to zero')
        call assert_real(TWOPI - PI/4., canonical_phshift(-PI/4.), TOL, &
            &'phase canonicalization maps negative values into [0,2pi)')
        call assert_real(3.0, canonical_phshift(3.0), TOL, &
            &'phase canonicalization leaves the near-pi phase-plate regime alone')
        ctfvals = tfun%get_ctfvars(TWOPI + PI/4.)
        call assert_real(PI/4., ctfvals%phshift, TOL, &
            &'ctfvars carries the canonical additive phase shift')
        ctfvals = tfun%get_ctfvars(-PIO2)
        call assert_real(TWOPI - PIO2, ctfvals%phshift, TOL, &
            &'ctfvars maps a negative phase shift into [0,2pi)')
        call assert_real(ctf_at(tfun, s2, ang, PI/4.), -ctf_at(tfun, s2, ang, PI + PI/4.), TOL, &
            &'a pi phase offset negates the transfer function')
        call assert_real(ctf_at(tfun, s2, ang, PI/4.), ctf_at(tfun, s2, ang, TWOPI + PI/4.), TOL, &
            &'CTF evaluation is invariant under a 2*pi phase offset')
        call assert_real(tfun%eval_canonical(s2, ang, PI/4.), tfun%eval_canonical(s2, ang, TWOPI + PI/4.), 1.e-5, &
            &'the canonical-input evaluator is itself 2*pi periodic (sine periodicity)')
    end subroutine test_phase_shift_policy

    !---------------- extrema and zeros ----------------

    ! nextrema counts the extrema of the CTF below a frequency: floor((chi + ac)/pi + 1/2) below the
    ! phase-aberration extremum; at the n-th zero exactly n extrema precede it. spafreqsqatnthzero
    ! (the fitting ranges of the CTF estimator) inverts the quadratic in s^2 for the n-th zero
    subroutine test_extrema()
        type(ctf) :: tfun
        real      :: ang, val, s2
        integer   :: n
        write(*,'(A)') 'test_extrema'
        tfun = make_tfun()
        call assert_int(3,  tfun%nextrema(spafreqsq(5,3),   angle(5,3),   0.), 'three extrema below (5,3)')
        call assert_int(16, tfun%nextrema(spafreqsq(10,-7), angle(10,-7), 0.), 'sixteen extrema below (10,-7)')
        call assert_int(51, tfun%nextrema(spafreqsq(20,12), angle(20,12), 0.), 'fifty-one extrema below (20,12)')
        call assert_int(0,  tfun%nextrema(0., 0., 0.), 'no extremum at the origin')
        ang = ANGAST*PI/180.
        do n = 1,3
            val = ctf_at(tfun, S2_ZEROS(n), ang, 0.)
            call assert_true(abs(val) < 1.e-3, 'the CTF vanishes at its zero '//int2str(n))
            call assert_int(n, tfun%nextrema(S2_ZEROS(n), ang, 0.), 'exactly '//int2str(n)//' extrema precede zero '//int2str(n))
        end do
        call assert_true(ctf_at(tfun, 0.5*(S2_ZEROS(1)+S2_ZEROS(2)), ang, 0.) * &
                        &ctf_at(tfun, 0.5*(S2_ZEROS(2)+S2_ZEROS(3)), ang, 0.) < 0., &
            &'the CTF changes sign across its second zero')
        call assert_int(3, tfun%nextrema(S2_ZEROS(2), ang, PI), 'a pi phase shift adds one extremum before the second zero')
        do n = 1,3
            s2 = tfun%spafreqsqatnthzero(n, 0., ang)
            call assert_real(S2_ZEROS(n), s2, 2.e-6, 'spafreqsqatnthzero finds zero '//int2str(n)//' (single-precision quadratic)')
            call assert_true(abs(ctf_at(tfun, s2, ang, 0.)) < 2.e-3, 'the CTF vanishes at the zero it computed ('//int2str(n)//')')
        end do
        call assert_true(tfun%spafreqsqatnthzero(1, PI/4., ang) < S2_ZEROS(1), 'a positive phase shift brings the first zero in')
        call assert_true(tfun%spafreqsqatnthzero(2, 0., ang + PIO2) < S2_ZEROS(2), 'the larger defocus across the axis brings the zeros in')
    end subroutine test_extrema

    !---------------- image-layer kernel ----------------

    ! The image-layer fast kernel must agree with the scalar CTF for conventional,
    ! intermediate, and near-quadrature phase shifts; ft2img places component (h,k) at
    ! (h+mh+1, k+mk+1) where mh, mk are the lower Fourier limits
    subroutine test_fast_kernel_and_ft2img()
        type(image) :: img, img_spec
        type(ctf)   :: tfun
        real        :: s2, ang, expected, actual, conventional, shifted
        integer     :: iph, lims(3,2), mh, mk
        write(*,'(A)') 'test_fast_kernel_and_ft2img'
        call img%new(LDIM, SMPD)
        call img_spec%new(LDIM, SMPD)
        tfun = ctf(SMPD, KV, CS, AC)
        call memoize_ft_maps(LDIM, SMPD)
        s2   = spafreqsq(H,K)
        ang  = angle(H,K)
        conventional = 0.
        shifted      = 0.
        do iph = 1,size(PHASES)
            call img%ctf2img(tfun, DFX, DFY, ANGAST, PHASES(iph))
            call tfun%init(DFX, DFY, ANGAST)
            expected = ctf_at(tfun, s2, ang, PHASES(iph))
            actual   = real(img%get_fcomp2D(H,K))
            call assert_real(expected, actual, TOL, 'fast CTF kernel agrees with scalar phase-shift CTF')
            if( iph == 1 ) conventional = actual
            if( iph == size(PHASES) ) shifted = actual
        enddo
        call img%ctf2img(tfun, DFX, DFY, ANGAST, 0.)
        call assert_real(conventional, real(img%get_fcomp2D(H,K)), TOL, &
            &'explicit zero phase shift preserves conventional CTF behavior')
        call assert_true(abs(shifted-conventional) > 0.1, &
            &'nonzero phase shift materially changes the transfer function')
        call assert_real(CTF_5_3_0,   real(img%get_fcomp2D(H,K)),   1.e-5, 'the CTF image carries the closed-form value at (5,3)')
        call assert_real(CTF_10_M7_0, real(img%get_fcomp2D(10,-7)), 5.e-5, 'the CTF image carries the closed-form value at (10,-7)')
        call assert_real(0., aimag(img%get_fcomp2D(H,K)), TOL, 'the CTF image is real')
        ! ft2img: real part, power and modulus land at the shifted (h,k) positions
        lims = img%loop_lims(3)
        mh   = abs(lims(1,1))
        mk   = abs(lims(2,1))
        call img%ft2img('real', img_spec)
        call assert_false(img_spec%is_ft(), 'ft2img produces a real-space image')
        call assert_real(real(img%get_fcomp2D(H,K)), img_spec%get_rmat_at(H+mh+1, K+mk+1, 1), TOL, &
            &'ft2img(real) places the real part of (h,k) at (h+mh+1,k+mk+1)')
        call assert_real(real(img%get_fcomp2D(10,-7)), img_spec%get_rmat_at(10+mh+1, -7+mk+1, 1), TOL, &
            &'ft2img(real) places a negative-k component below the centre')
        call img%ft2img('power', img_spec)
        call assert_real(real(img%get_fcomp2D(H,K))**2, img_spec%get_rmat_at(H+mh+1, K+mk+1, 1), TOL, &
            &'ft2img(power) stores the squared modulus')
        call img%ft2img('sqrt', img_spec)
        call assert_real(abs(real(img%get_fcomp2D(H,K))), img_spec%get_rmat_at(H+mh+1, K+mk+1, 1), TOL, &
            &'ft2img(sqrt) stores the modulus')
        call assert_true(img%is_ft(), 'ft2img leaves an input that was already in Fourier space there')
        call img%kill
        call img_spec%kill
    end subroutine test_fast_kernel_and_ft2img

    !---------------- 3D restoration contracts ----------------

    ! Reconstruction uses H*y and H^2 for raw observations, and |H|*y and H^2
    ! for phase-flipped observations.  A unit Fourier image makes both contracts
    ! directly observable at a native-grid sample.
    subroutine test_restoration_contracts()
        type(image)       :: img
        type(ctfparams)   :: ctfparms
        type(fplane_type) :: fplane
        real              :: actual
        integer           :: hp, kp
        write(*,'(A)') 'test_restoration_contracts'
        call img%new(LDIM, SMPD)
        call memoize_ft_maps(LDIM, SMPD)
        call img%set_cmat(cmplx(1.,0.))   ! puts the image in Fourier space with unit components
        ctfparms%smpd        = SMPD
        ctfparms%kv          = KV
        ctfparms%cs          = CS
        ctfparms%fraca       = AC
        ctfparms%dfx         = DFX
        ctfparms%dfy         = DFY
        ctfparms%angast      = ANGAST
        ctfparms%phshift     = PI/4.
        ctfparms%ctfflag     = CTFFLAG_YES
        call img%gen_fplane4rec([0,LDIM(1)/(2*OSMPL_PAD_FAC)], SMPD, ctfparms, [0.,0.], fplane)
        hp = OSMPL_PAD_FAC*H
        kp = -OSMPL_PAD_FAC*K
        actual = real(fplane%cmplx_plane(hp,kp))
        call assert_real(actual*actual, fplane%ctfsq_plane(hp,kp), 5.e-6, &
            &'raw 3D restoration stores H*y over H^2')
        ctfparms%ctfflag = CTFFLAG_FLIP
        call img%gen_fplane4rec([0,LDIM(1)/(2*OSMPL_PAD_FAC)], SMPD, ctfparms, [0.,0.], fplane)
        actual = real(fplane%cmplx_plane(hp,kp))
        call assert_true(actual >= 0., 'phase-flipped 3D restoration uses nonnegative |H|')
        call assert_real(actual*actual, fplane%ctfsq_plane(hp,kp), 5.e-6, &
            &'phase-flipped 3D restoration stores |H|*y over H^2')
        call img%kill
    end subroutine test_restoration_contracts

end module simple_ctf_tester
