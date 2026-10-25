!@descr: unit tests of wave_multislice: adjoint of the sweep, finite-difference gradient, signs and midplane, analytic limits
! The chain alone, on cube spectra built here (no lattice interpolation in the way).
! - The transpose of the linearisation at a random point of about half a radian: Re<J dC, q> = Re<dC, J^T q> to
!   1e-5 (single-precision chains, double-precision dot products), both curvature signs.
! - The gradient of 1/2 |y - ihat(C)|^2 by the adjoint against a fourth-order central difference at 3e-2 of the
!   point's scale, to 1e-3, both signs.
! - A weak sheet (random in the plane, band-limited to half the native Nyquist; a Gaussian of 1.5 sections along
!   the beam), centred on the midplane, images as SIMPLE's weak-phase model of a thick object: each section j
!   at its own defocus df + z_j a, -2 sigma a sum_j sin(eval_phase + atan(fraca/sqrt(1-fraca^2))) s_j: the signs
!   of transmission, propagation and lens and the midplane reference at once. Bound 5e-3 of the largest term:
!   the second order is at most about half the projected phase's peak (5e-3 rad) relative to the first, single
!   precision through the 64 sections about 1e-3 at that contrast, the sheet's truncation to the sphere 2e-4.
!   The same sheet 16 sections downstream (ewald=neg) images at df + 16 A to the same bound and
!   at df - 16 A at least ten times worse; with the beam entering through the last section (ewald=pos) the
!   sheet is upstream and df - 16 A holds.
! - Analytic limits: zero potential and a constant real potential give the vacuum (intensity 1), a constant
!   potential with amplitude contrast attenuates as exp(-2 sigma a N fraca V).
module simple_wave_multislice_tester
use simple_defs,            only: dp, PI
use simple_type_defs,       only: ctfparams, CTFFLAG_YES
use simple_ctf,             only: ctf
use simple_math_ctf,        only: electron_interaction_constant
use simple_wave_multislice, only: wave_multislice
use simple_test_utils,      only: assert_true, fix_test_seed
implicit none
private
public :: run_all_wave_multislice_tests

real, parameter :: SMPD = 1.0, KV = 300.

contains

    subroutine run_all_wave_multislice_tests()
        write(*,'(A)') '**** running all wave multislice tests ****'
        call fix_test_seed()
        call test_adjoint_and_gradient(-1)
        call test_adjoint_and_gradient(1)
        call test_weak_sheet()
        call test_analytic_limits()
    end subroutine run_all_wave_multislice_tests

    function optics( fraca ) result( c )
        real, intent(in) :: fraca
        type(ctfparams) :: c
        c%ctfflag = CTFFLAG_YES
        c%smpd    = SMPD
        c%kv      = KV
        c%cs      = 2.7
        c%fraca   = fraca
        c%dfx     = 1.2
        c%dfy     = 1.0
        c%angast  = 30.
        c%phshift = 0.
    end function optics

    subroutine test_adjoint_and_gradient( kappa )
        integer, intent(in) :: kappa
        integer, parameter :: N = 24, R = N/2, RL = R-1
        class(wave_multislice), allocatable :: ms
        complex, allocatable :: c(:,:,:), dc(:,:,:), cadj(:,:,:), cg(:,:,:)
        complex  :: q(-R:R,-R:R), y(-R:R,-R:R), ihat(-R:R,-R:R), jdc(-R:R,-R:R)
        real(dp) :: lhs, rhs, fd, an
        real     :: h, sig
        allocate(ms)
        call ms%new(N, 2)
        call ms%set_optics(optics(0.1), 20., kappa)
        sig = electron_interaction_constant(KV) * SMPD
        ! about half a radian through the box: |u| ~ sqrt(npoints) |C| / N^3 per voxel, summed over N sections
        c  = random_cube(N) * (0.5 / sig) * real(N)**2 / sqrt(real(count_sphere(N)))
        dc = random_cube(N)
        q  = random_disk(R)
        call ms%forward(c, ihat)
        call ms%linearize(dc, jdc)
        allocate(cadj, mold=c)
        call ms%adjoint(q, cadj)
        lhs = redot2(jdc, q)
        rhs = redot3(dc, cadj)
        write(*,'(A,I3,A,ES10.3)') '    kappa', kappa, ' multislice adjoint relative error ', abs(lhs-rhs)/max(abs(lhs),abs(rhs))
        call assert_true(abs(lhs - rhs) <= 1.0e-5_dp * max(abs(lhs), abs(rhs)), &
            &'the multislice adjoint is the transpose of its linearisation')
        ! the gradient of 1/2 |y - ihat|^2 along dc
        y = ihat + 0.3 * maxval(abs(ihat(1:,:))) * random_disk(R)
        call ms%forward(c, ihat)
        allocate(cg, mold=c)
        call ms%adjoint(-(y - ihat), cg)
        an = redot3(cg, dc)
        h  = 3.0e-2 * maxval(abs(c)) / maxval(abs(dc))
        fd = 8.0_dp * (misfit(ms, y, c + h*dc) - misfit(ms, y, c - h*dc)) - (misfit(ms, y, c + 2.*h*dc) - &
            &misfit(ms, y, c - 2.*h*dc))
        fd = fd / (12.0_dp * real(h,dp))
        write(*,'(A,I3,A,ES10.3)') '    kappa', kappa, ' multislice gradient relative error ', abs(fd-an)/abs(an)
        call assert_true(abs(fd - an) <= 1.0e-3_dp * abs(an), 'the multislice adjoint gives the gradient of the misfit')
        call ms%kill
    end subroutine test_adjoint_and_gradient

    subroutine test_weak_sheet()
        integer, parameter :: N = 64, R = N/2
        real,    parameter :: PEAK_PHASE = 5.0e-3, ZOFF = 16.
        class(wave_multislice), allocatable :: ms
        complex, allocatable :: sj(:,:,:), cube(:,:,:)
        complex :: ihat(-R:R,-R:R)
        real    :: err0, errp, errm, errpos
        allocate(ms)
        call ms%new(N, 2)
        ! the sheet in the midplane, ewald=neg
        call sheet(0., sj, cube)
        call ms%set_optics(optics(0.1), 0., -1)
        call ms%forward(cube, ihat)
        err0 = rel_err(ihat, weak_phase(sj, 1.))
        write(*,'(A,ES10.3)') '    weak sheet in the midplane against the thick weak-phase model: ', err0
        call assert_true(err0 <= 5.e-3, 'a weak midplane sheet images as SIMPLE''s weak-phase model (signs, reference)')
        ! 16 sections downstream
        call sheet(ZOFF, sj, cube)
        call ms%forward(cube, ihat)
        errp = rel_err(ihat, weak_phase(sj,  1.))
        errm = rel_err(ihat, weak_phase(sj, -1.))
        write(*,'(A,2ES10.3)') '    sheet 16 A downstream, ewald=neg, against df+z and df-z: ', errp, errm
        call assert_true(errp <= 5.e-3 .and. errm >= 10. * errp, &
            &'a sheet downstream of the midplane images at the larger defocus (ewald=neg)')
        ! the beam through the last section: the same sheet is upstream
        call ms%set_optics(optics(0.1), 0., 1)
        call ms%forward(cube, ihat)
        errpos = rel_err(ihat, weak_phase(sj, -1.))
        write(*,'(A,ES10.3)') '    the same sheet, ewald=pos, against df-z: ', errpos
        call assert_true(errpos <= 5.e-3, 'with the beam reversed (ewald=pos) the sheet images at the smaller defocus')
        call ms%kill

    contains

        !> the sheet's section spectra sj(h,k,j) (native DFT units, volts) and its cube spectrum, centred z0
        !! sections downstream of the midplane
        subroutine sheet( z0, sj, cube )
            real,                 intent(in)  :: z0
            complex, allocatable, intent(out) :: sj(:,:,:), cube(:,:,:)
            complex :: s0(-R:R,-R:R), a(-R:R,-R:R)
            real    :: g(N), re, im, zj, scale
            integer :: h, k, l, j
            s0 = cmplx(0.,0.)
            do k = -R, R
                do h = -R, R
                    if( h*h + k*k > (R/2)**2 ) cycle
                    call random_number(re)
                    call random_number(im)
                    s0(h,k) = cmplx(re - 0.5, im - 0.5)
                end do
            end do
            ! Hermitian: the section of a real potential (from a copy: in place, the pairs would mix)
            a = s0
            do k = -R, R
                do h = -R, R
                    s0(h,k) = 0.5 * (a(h,k) + conjg(a(-h,-k)))
                end do
            end do
            do j = 1, N
                zj   = real(j - (N/2 + 1))
                g(j) = exp(-0.5 * ((zj - z0) / 1.5)**2)
            end do
            ! the projected phase's peak, sigma a max_x |sum_j u_j(x)|, is PEAK_PHASE: the real-space peak of
            ! s0 by a direct inverse DFT (the origin at the box centre, as SIMPLE's)
            scale = PEAK_PHASE / (electron_interaction_constant(KV) * SMPD * sum(g) * peak_of(s0))
            allocate(sj(-R:R,-R:R,N), cube(-R:R,-R:R,-(R-1):R-1))
            do j = 1, N
                sj(:,:,j) = scale * g(j) * s0
            end do
            ! C(h,k,l) = sum_j s_j(h,k) exp(-2 pi i l z_j / N), on the sphere, |l| < R
            cube = cmplx(0.,0.)
            do l = -(R-1), R-1
                do k = -R, R
                    do h = -R, R
                        if( h*h + k*k + l*l > R*R ) cycle
                        do j = 1, N
                            zj = real(j - (N/2 + 1))
                            cube(h,k,l) = cube(h,k,l) + sj(h,k,j) * cmplx(cos(2.*PI*real(l)*zj/real(N)), &
                                &-sin(2.*PI*real(l)*zj/real(N)))
                        end do
                    end do
                end do
            end do
        end subroutine sheet

        !> max_x |(1/N^2) sum_(h,k) s(h,k) exp(2 pi i (h x + k y) / N)|
        real function peak_of( s )
            complex, intent(in) :: s(-R:R,-R:R)
            complex :: v
            integer :: x, y, h, k
            peak_of = 0.
            do y = -R, R-1
                do x = -R, R-1
                    v = cmplx(0.,0.)
                    do k = -R/2, R/2
                        do h = -R/2, R/2
                            v = v + s(h,k) * cmplx(cos(2.*PI*real(h*x+k*y)/real(N)), sin(2.*PI*real(h*x+k*y)/real(N)))
                        end do
                    end do
                    peak_of = max(peak_of, abs(v) / real(N)**2)
                end do
            end do
        end function peak_of

        !> the first-order image of the sections sj, each at defocus df + side z_j a: -2 sigma a sum_j CTF_j s_j
        function weak_phase( sj, side ) result( dihat )
            complex, intent(in) :: sj(-R:R,-R:R,N)
            real,    intent(in) :: side
            complex :: dihat(-R:R,-R:R)
            type(ctf)       :: tfun
            type(ctfparams) :: c
            real    :: sig, acc, ph, ddf
            integer :: h, k, j
            c    = optics(0.1)
            sig  = electron_interaction_constant(KV) * SMPD
            acc  = atan(c%fraca / sqrt(1. - c%fraca**2))
            dihat = cmplx(0.,0.)
            do j = 1, N
                ddf  = side * real(j - (N/2 + 1)) * SMPD * 1.e-4
                tfun = ctf(c%smpd, c%kv, c%cs, c%fraca)
                call tfun%init(c%dfx + ddf, c%dfy + ddf, c%angast)
                do k = -R, R
                    do h = -R, R
                        if( h*h + k*k > R*R ) cycle
                        if( h == 0 .and. k == 0 ) cycle
                        ph = tfun%eval_phase(real(h*h + k*k) / real(N)**2, atan2(real(k), real(h)), c%phshift)
                        dihat(h,k) = dihat(h,k) - 2. * sig * sin(ph + acc) * sj(h,k,j)
                    end do
                end do
            end do
        end function weak_phase

        !> the image's departure from the vacuum (all but the origin) against the first-order model
        real function rel_err( ihat, dexp )
            complex, intent(in) :: ihat(-R:R,-R:R), dexp(-R:R,-R:R)
            complex :: d(-R:R,-R:R)
            d = ihat
            d(0,0) = cmplx(0.,0.)
            rel_err = maxval(abs(d - dexp)) / maxval(abs(dexp))
        end function rel_err

    end subroutine test_weak_sheet

    subroutine test_analytic_limits()
        integer, parameter :: N = 24, R = N/2
        real,    parameter :: V0 = 0.25
        class(wave_multislice), allocatable :: ms
        complex, allocatable :: cube(:,:,:)
        complex :: ihat(-R:R,-R:R)
        real    :: expected
        allocate(ms)
        call ms%new(N, 2)
        allocate(cube(-R:R,-R:R,-(R-1):R-1), source=cmplx(0.,0.))
        call ms%set_optics(optics(0.), 64., -1)
        call ms%forward(cube, ihat)
        call assert_true(vacuum_dev(ihat, 1.) < 2.e-5, 'zero potential gives the vacuum')
        ! a constant potential: the cube's origin holds N^3 V0
        cube(0,0,0) = cmplx(real(N)**3 * V0, 0.)
        call ms%forward(cube, ihat)
        call assert_true(vacuum_dev(ihat, 1.) < 3.e-5, 'a constant real potential is a pure phase object')
        call ms%set_optics(optics(0.1), 64., -1)
        expected = exp(-2. * electron_interaction_constant(KV) * SMPD * real(N) * 0.1 * V0)
        call ms%forward(cube, ihat)
        call assert_true(vacuum_dev(ihat, expected) < 3.e-5, &
            &'amplitude contrast attenuates a constant potential as exp(-2 sigma a N fraca V)')
        call ms%kill

    contains

        !> the largest departure from a flat intensity i0: |ihat(0,0)/N^2 - i0| and |ihat|/N^2 elsewhere
        real function vacuum_dev( ihat, i0 )
            complex, intent(in) :: ihat(-R:R,-R:R)
            real,    intent(in) :: i0
            complex :: d(-R:R,-R:R)
            d = ihat / real(N)**2
            d(0,0) = d(0,0) - i0
            vacuum_dev = maxval(abs(d))
        end function vacuum_dev

    end subroutine test_analytic_limits

    ! ---- helpers ----

    real(dp) function misfit( ms, y, c )
        class(wave_multislice), intent(inout) :: ms
        complex,                intent(in)    :: y(:,:), c(:,:,:)
        complex, allocatable :: ihat(:,:)
        allocate(ihat, mold=y)
        call ms%forward(c, ihat)
        misfit = 0.5_dp * sum(real(abs(y - ihat),dp)**2)
    end function misfit

    function random_cube( n ) result( c )
        integer, intent(in) :: n
        complex, allocatable :: c(:,:,:)
        real,    allocatable :: re(:,:,:), im(:,:,:)
        integer :: h, k, l, r
        r = n / 2
        allocate(re(-r:r,-r:r,-(r-1):r-1), im(-r:r,-r:r,-(r-1):r-1))
        call random_number(re)
        call random_number(im)
        allocate(c(-r:r,-r:r,-(r-1):r-1))
        c = cmplx(re - 0.5, im - 0.5)
        do l = -(r-1), r-1
            do k = -r, r
                do h = -r, r
                    if( h*h + k*k + l*l > r*r ) c(h,k,l) = cmplx(0.,0.)
                end do
            end do
        end do
    end function random_cube

    integer function count_sphere( n )
        integer, intent(in) :: n
        integer :: h, k, l, r
        r = n / 2
        count_sphere = 0
        do l = -(r-1), r-1
            do k = -r, r
                do h = -r, r
                    if( h*h + k*k + l*l <= r*r ) count_sphere = count_sphere + 1
                end do
            end do
        end do
    end function count_sphere

    function random_disk( r ) result( s )
        integer, intent(in) :: r
        complex :: s(-r:r,-r:r)
        real    :: re(-r:r,-r:r), im(-r:r,-r:r)
        integer :: h, k
        call random_number(re)
        call random_number(im)
        s = cmplx(re - 0.5, im - 0.5)
        do k = -r, r
            do h = -r, r
                if( h*h + k*k > r*r ) s(h,k) = cmplx(0.,0.)
            end do
        end do
    end function random_disk

    real(dp) function redot2( a, b )
        complex, intent(in) :: a(:,:), b(:,:)
        redot2 = sum(real(a,dp)*real(b,dp) + real(aimag(a),dp)*real(aimag(b),dp))
    end function redot2

    real(dp) function redot3( a, b )
        complex, intent(in) :: a(:,:,:), b(:,:,:)
        redot3 = sum(real(a,dp)*real(b,dp) + real(aimag(a),dp)*real(aimag(b),dp))
    end function redot3

end module simple_wave_multislice_tester
