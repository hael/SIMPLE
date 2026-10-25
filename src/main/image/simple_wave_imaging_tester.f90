!@descr: unit tests of wave_imaging: adjoint identity, finite-difference gradient, weak-phase limit, model hierarchy
! Native box 24 on the twice-finer grid, 300 kV, astigmatic underfocus, Cs 2.7 mm, amplitude contrast 0.1,
! phase shift 0.2 rad, a B-factor of 20 A^2, a random spectrum of about half a radian of phase.
! - The transpose of the linearisation at a random point: Re<J ds, q> = Re<ds, J^T q> to 1e-5 (the plan's
!   bound for single-precision chains with double-precision dot products), every model.
! - The gradient of 1/2 |y - ihat(s)|^2 by the adjoint against a fourth-order central difference at 3e-2 of
!   the point's scale, to 1e-3 (the plan's bound).
! - The weak-phase limit: at s = 0 the linearisation of a Hermitian ds (the section of a real potential) is
!   -2 sigma a CTF ds, with SIMPLE's CTF sin(eval_phase + atan(fraca/sqrt(1-fraca^2))) and the B-factor
!   envelope, to 1e-5 of its largest value, every model; at nonzero frequencies, since the unscattered beam
!   (the origin, excluded from every misfit) carries no CTF.
! - First Born is the linearisation of Rytov: at s = 0 their linearisations agree to 1e-6.
! - The lens and detection of a given exit wave (multislice's entry): the image of the exit wave projection
!   formed is projection's image, bitwise; the transpose of their linearisation to 1e-5.
! - The detector of simulate_scattering: Poisson-sampled intensities are nonnegative counts over the dose.
module simple_wave_imaging_tester
use simple_defs,         only: dp, PI
use simple_type_defs,    only: ctfparams, CTFFLAG_YES
use simple_ctf,          only: ctf
use simple_math_ctf,     only: electron_interaction_constant
use simple_wave_imaging, only: wave_imaging, WAVE_PROJECTION, WAVE_FIRST_BORN, WAVE_RYTOV, WAVE_MSLICE, &
    &detect_intensity
use simple_image,        only: image
use simple_test_utils,   only: assert_true, fix_test_seed
implicit none
private
public :: run_all_wave_imaging_tests

integer, parameter :: NAT = 24, R = NAT/2, OSF = 2
real,    parameter :: BFAC = 20.

contains

    subroutine run_all_wave_imaging_tests()
        integer :: model
        write(*,'(A)') '**** running all wave imaging tests ****'
        call fix_test_seed()
        do model = WAVE_PROJECTION, WAVE_RYTOV
            call test_adjoint_and_gradient(model)
            call test_weak_phase(model)
        end do
        call test_born_is_rytov_linearised()
        call test_exit_wave()
        call test_detector()
    end subroutine run_all_wave_imaging_tests

    function optics() result( c )
        type(ctfparams) :: c
        c%ctfflag = CTFFLAG_YES
        c%smpd    = 1.0
        c%kv      = 300.
        c%cs      = 2.7
        c%fraca   = 0.1
        c%dfx     = 1.2
        c%dfy     = 1.0
        c%angast  = 30.
        c%phshift = 0.2
    end function optics

    subroutine test_adjoint_and_gradient( model )
        integer, intent(in) :: model
        type(wave_imaging) :: wi
        complex  :: s(-R:R,-R:R), ds(-R:R,-R:R), q(-R:R,-R:R), y(-R:R,-R:R), ihat(-R:R,-R:R)
        complex  :: jds(-R:R,-R:R), jtq(-R:R,-R:R), grad(-R:R,-R:R)
        real(dp) :: lhs, rhs, fd, an
        real     :: h
        call wi%new(NAT, OSF, model)
        call wi%set_optics(optics(), BFAC, real(R))
        ! about half a radian of phase: u ~ |s|/nat^2 summed over the disk
        s  = spectrum_scaled(0.5)
        ds = random_disk()
        q  = random_disk()
        call wi%forward(s, ihat)
        call wi%linearize(ds, jds)
        call wi%adjoint(q, jtq)
        lhs = redot(jds, q)
        rhs = redot(ds, jtq)
        write(*,'(A,I2,A,ES10.3)') '    model ', model, ' adjoint relative error ', abs(lhs-rhs)/max(abs(lhs),abs(rhs))
        call assert_true(abs(lhs - rhs) <= 1.0e-5_dp * max(abs(lhs), abs(rhs)), &
            &'the adjoint is the transpose of the linearisation')
        ! the gradient of 1/2 |y - ihat|^2 along ds
        y = ihat + 0.3 * maxval(abs(ihat(1:,:))) * random_disk()
        call wi%forward(s, ihat)
        call wi%adjoint(-(y - ihat), grad)
        an = redot(grad, ds)
        ! fourth-order central difference at 3e-2 of the point's scale: truncation O(h^4) about 1e-6,
        ! the single-precision chain's rounding about 2e-5 of the derivative
        h   = 3.0e-2 * maxval(abs(s)) / maxval(abs(ds))
        fd  = 8.0_dp * (misfit2d(wi, y, s + h*ds) - misfit2d(wi, y, s - h*ds)) &
            &- (misfit2d(wi, y, s + 2.*h*ds) - misfit2d(wi, y, s - 2.*h*ds))
        fd  = fd / (12.0_dp * real(h,dp))
        write(*,'(A,I2,A,ES10.3)') '    model ', model, ' gradient relative error ', abs(fd-an)/abs(an)
        call assert_true(abs(fd - an) <= 1.0e-3_dp * abs(an), 'the adjoint gives the gradient of the misfit')
        call wi%kill
    end subroutine test_adjoint_and_gradient

    !> at s = 0: J ds = -2 sigma a CTF ds for a Hermitian ds
    subroutine test_weak_phase( model )
        integer, intent(in) :: model
        type(wave_imaging) :: wi
        type(ctfparams)    :: c
        type(ctf)          :: tfun
        complex :: s0(-R:R,-R:R), ds(-R:R,-R:R), ihat(-R:R,-R:R), jds(-R:R,-R:R), expected(-R:R,-R:R)
        real    :: sig, acc, ph, s2, ang, env, err
        integer :: h, k
        c = optics()
        call wi%new(NAT, OSF, model)
        call wi%set_optics(c, BFAC, real(R))
        s0 = cmplx(0.,0.)
        ds = hermitian(random_disk())
        call wi%forward(s0, ihat)
        call wi%linearize(ds, jds)
        jds(0,0) = cmplx(0.,0.)
        sig  = electron_interaction_constant(c%kv)
        acc  = atan(c%fraca / sqrt(1. - c%fraca**2))
        tfun = ctf(c%smpd, c%kv, c%cs, c%fraca)
        call tfun%init(c%dfx, c%dfy, c%angast)
        expected = cmplx(0.,0.)
        do k = -R, R
            do h = -R, R
                if( h*h + k*k > R*R ) cycle
                ! the unscattered beam itself: no CTF there (only absorption), and no misfit either
                if( h == 0 .and. k == 0 ) cycle
                s2  = real(h*h + k*k) / real(NAT)**2
                ang = atan2(real(k), real(h))
                ph  = tfun%eval_phase(s2, ang, c%phshift)
                ! the lens envelope multiplies both terms of the intensity: exp(-B s^2 / 4) once each
                env = exp(-0.25 * BFAC * s2 / c%smpd**2)
                expected(h,k) = -2. * sig * c%smpd * env * sin(ph + acc) * ds(h,k)
            end do
        end do
        err = maxval(abs(jds - expected)) / maxval(abs(expected))
        write(*,'(A,I2,A,ES10.3)') '    model ', model, ' weak-phase identity relative error ', err
        call assert_true(err <= 1.0e-5, 'the weak-phase limit is SIMPLE''s CTF times the section')
        call wi%kill
    end subroutine test_weak_phase

    subroutine test_born_is_rytov_linearised()
        type(wave_imaging) :: wr, wb
        complex :: s0(-R:R,-R:R), ds(-R:R,-R:R), ihat(-R:R,-R:R), jr(-R:R,-R:R), jb(-R:R,-R:R)
        call wr%new(NAT, OSF, WAVE_RYTOV)
        call wb%new(NAT, OSF, WAVE_FIRST_BORN)
        call wr%set_optics(optics(), BFAC, real(R))
        call wb%set_optics(optics(), BFAC, real(R))
        s0 = cmplx(0.,0.)
        ds = random_disk()
        call wr%forward(s0, ihat)
        call wb%forward(s0, ihat)
        call wr%linearize(ds, jr)
        call wb%linearize(ds, jb)
        call assert_true(maxval(abs(jr - jb)) <= 1.0e-6 * maxval(abs(jr)), 'first Born is the linearisation of Rytov')
        call wr%kill
        call wb%kill
    end subroutine test_born_is_rytov_linearised

    subroutine test_exit_wave()
        type(wave_imaging) :: wp, wm
        complex, allocatable :: psi(:,:), dpsi(:,:), adj(:,:)
        complex  :: s(-R:R,-R:R), ihat(-R:R,-R:R), ihat_m(-R:R,-R:R), q(-R:R,-R:R), dihat(-R:R,-R:R)
        real     :: re(NAT*OSF,NAT*OSF), im(NAT*OSF,NAT*OSF)
        real(dp) :: lhs, rhs
        call wp%new(NAT, OSF, WAVE_PROJECTION)
        call wm%new(NAT, OSF, WAVE_MSLICE)
        call wp%set_optics(optics(), BFAC, 0.)
        call wm%set_optics(optics(), BFAC, 0.)
        s = spectrum_scaled(0.5)
        call wp%forward(s, ihat)
        allocate(psi(NAT*OSF,NAT*OSF), dpsi(NAT*OSF,NAT*OSF), adj(NAT*OSF,NAT*OSF))
        call wp%get_exit_wave(psi)
        call wm%forward_wave(psi, ihat_m)
        call assert_true(all(ihat_m == ihat), 'the image of a given exit wave is the full chain''s, bitwise')
        call random_number(re)
        call random_number(im)
        dpsi = cmplx(re - 0.5, im - 0.5)
        q    = random_disk()
        call wm%linearize_wave(dpsi, dihat)
        call wm%adjoint_wave(q, adj)
        lhs = redot(dihat, q)
        rhs = sum(real(dpsi,dp)*real(adj,dp) + real(aimag(dpsi),dp)*real(aimag(adj),dp))
        call assert_true(abs(lhs - rhs) <= 1.0e-5_dp * max(abs(lhs), abs(rhs)), &
            &'the exit wave''s lens and detection: the adjoint is the transpose of the linearisation')
        call wp%kill
        call wm%kill
    end subroutine test_exit_wave

    subroutine test_detector()
        type(image)       :: img
        real, allocatable :: rmat(:,:,:), intensity(:,:,:)
        real, parameter   :: DOSE = 20.
        call img%new([NAT,NAT,1], 1.5)
        allocate(rmat(NAT,NAT,1), source=0.9)
        call img%set_rmat(rmat, .false.)
        call detect_intensity(img, DOSE, 0.)
        intensity = img%get_rmat()
        call assert_true(all(intensity >= 0.0), 'Poisson-sampled detector intensity is nonnegative')
        call assert_true(maxval(abs(intensity * (DOSE * 1.5**2) - real(nint(intensity * (DOSE * 1.5**2))))) < 2.0e-5, &
            &'dose-normalized detector values correspond to integral electron counts')
        call img%kill
    end subroutine test_detector

    ! ---- helpers ----

    real(dp) function misfit2d( wi, y, s )
        type(wave_imaging), intent(inout) :: wi
        complex,            intent(in)    :: y(-R:R,-R:R), s(-R:R,-R:R)
        complex :: ihat(-R:R,-R:R)
        call wi%forward(s, ihat)
        misfit2d = 0.5_dp * sum(real(abs(y - ihat),dp)**2)
    end function misfit2d

    !> a random spectrum whose field carries about phase radians
    function spectrum_scaled( phase ) result( s )
        real, intent(in) :: phase
        complex :: s(-R:R,-R:R)
        real    :: sig
        sig = electron_interaction_constant(300.)
        s   = random_disk()
        ! the field is sum(s)/nat^2 per pixel in magnitude about sqrt(npoints)|s|/nat^2
        s   = s * (phase / sig) * real(NAT)**2 / sqrt(real(count_disk()))
    end function spectrum_scaled

    integer function count_disk()
        integer :: h, k
        count_disk = 0
        do k = -R, R
            do h = -R, R
                if( h*h + k*k <= R*R ) count_disk = count_disk + 1
            end do
        end do
    end function count_disk

    function random_disk() result( s )
        complex :: s(-R:R,-R:R)
        real    :: re(-R:R,-R:R), im(-R:R,-R:R)
        integer :: h, k
        call random_number(re)
        call random_number(im)
        s = cmplx(re - 0.5, im - 0.5)
        do k = -R, R
            do h = -R, R
                if( h*h + k*k > R*R ) s(h,k) = cmplx(0.,0.)
            end do
        end do
    end function random_disk

    !> the Hermitian part, s(-h,-k) = conj(s(h,k)), as the section of a real potential
    function hermitian( a ) result( s )
        complex, intent(in) :: a(-R:R,-R:R)
        complex :: s(-R:R,-R:R)
        integer :: h, k
        do k = -R, R
            do h = -R, R
                s(h,k) = 0.5 * (a(h,k) + conjg(a(-h,-k)))
            end do
        end do
    end function hermitian

    real(dp) function redot( a, b )
        complex, intent(in) :: a(-R:R,-R:R), b(-R:R,-R:R)
        redot = sum(real(a,dp)*real(b,dp) + real(aimag(a),dp)*real(aimag(b),dp))
    end function redot

end module simple_wave_imaging_tester
