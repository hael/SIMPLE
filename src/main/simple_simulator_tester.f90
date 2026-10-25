!@descr: unit tests of simulate_particles' wave mode (simple_simulator, wave_projector%plane_to_image, the scale link)
! Box 64 at 1 A, 300 kV, astigmatic underfocus of 0.4 microns (no CTF delocalisation past the box), a map of smooth Gaussian blobs at SIMPLE's scale.
! - The weak-signal limit through the program's chain (map / scale link, wave_projector%predict, plane_to_image,
!   simimg_wave without noise): the part odd in the map's sign is twice simulate_particles' flat particle of the
!   deapodised map (fproject, shift, simimg with the CTF), PCG's scale; the even part is second order.
! - Per-particle lens B-factors: a particle's prediction is the one of an operator with its B-factor for all.
! - The noise of simimg_wave in an empty region: counts of mean D a^2 and variance D a^2 + sdev^2.
module simple_simulator_tester
use simple_core_module_api
use simple_test_utils
use simple_image,          only: image
use simple_projector,      only: projector
use simple_ctf,            only: ctf
use simple_simulator,      only: simimg, simimg_wave
use simple_math_ctf,       only: electron_scale_link, electron_interaction_constant
use simple_wave_imaging,   only: WAVE_PROJECTION, WAVE_FIRST_BORN
use simple_wave_projector, only: wave_projector
use simple_memoize_ft_maps, only: memoize_ft_maps
use simple_kbinterpol,      only: kbinterpol
use simple_gridding,        only: kb_stencil_centered_crop_inv_envelope_1d
implicit none
private
public :: run_all_simulator_tests
#include "simple_local_flags.inc"

integer, parameter :: BOX = 64
real,    parameter :: SMPD = 1.0, KV = 300., CS = 2.7, FRACA = 0.1, DOSE = 40.
real,    parameter :: S_BLOB = 2.0
! the lattice's interpolation error against an exact reference (TOL_GOLDEN of the wave projector tests)
real,    parameter :: TOL_LATTICE = 1.10e-2

contains

    subroutine run_all_simulator_tests()
        write(*,'(A)') '**** running all simulator tests ****'
        call test_weak_limit(WAVE_PROJECTION)
        call test_weak_limit(WAVE_FIRST_BORN)
        call test_bfac_per_particle()
        call test_detector_noise()
    end subroutine run_all_simulator_tests

    ! three blobs of width S_BLOB voxels near the centre, amplitude amp at SIMPLE's scale
    subroutine blob_map( vol, amp )
        class(image), intent(inout) :: vol
        real,         intent(in)    :: amp
        real, parameter :: CEN(3,3) = reshape([3.,-2.,1., -3.,2.,-2., 1.,3.,3.], [3,3])
        real, allocatable :: rmat(:,:,:)
        integer :: i, j, k, b
        real    :: r2
        call vol%new([BOX,BOX,BOX], SMPD)
        allocate(rmat(BOX,BOX,BOX), source=0.)
        do k = 1, BOX
            do j = 1, BOX
                do i = 1, BOX
                    do b = 1, 3
                        r2 = (real(i - BOX/2 - 1) - CEN(1,b))**2 + (real(j - BOX/2 - 1) - CEN(2,b))**2 + &
                            &(real(k - BOX/2 - 1) - CEN(3,b))**2
                        rmat(i,j,k) = rmat(i,j,k) + amp * exp(-r2 / (2. * S_BLOB**2))
                    end do
                end do
            end do
        end do
        call vol%set_rmat(rmat, .false.)
    end subroutine blob_map

    ! a particle at a fixed pose and shift, astigmatic underfocus
    subroutine test_particle( o )
        class(ori), intent(inout) :: o
        type(ctfparams) :: ctfparms
        call o%new(is_ptcl=.true.)
        call o%set_euler([37., 68., 121.])
        call o%set('x', 0.)
        call o%set('y', 0.)
        ctfparms%ctfflag = CTFFLAG_YES
        ctfparms%smpd    = SMPD
        ctfparms%kv      = KV
        ctfparms%cs      = CS
        ctfparms%fraca   = FRACA
        ctfparms%dfx     = 0.45
        ctfparms%dfy     = 0.35
        ctfparms%angast  = 30.
        ctfparms%phshift = 0.
        call o%set_ctfvars(ctfparms)
        call o%set_state(1)
    end subroutine test_particle

    ! simulate_particles' flat chain without noise: project, shift, CTF, clip
    subroutine flat_particle( vol, o, img )
        class(image), intent(inout) :: vol, img
        class(ori),   intent(inout) :: o
        type(projector) :: vol_pad
        type(image)     :: img_pad
        type(ctf)       :: tfun
        tfun = ctf(SMPD, KV, CS, FRACA)
        call vol_pad%new([2*BOX,2*BOX,2*BOX], SMPD)
        call vol%pad(vol_pad)
        call vol_pad%fft
        call vol_pad%expand_cmat()
        call img_pad%new([2*BOX,2*BOX,1], SMPD)
        call memoize_ft_maps(img_pad%get_ldim(), img_pad%get_smpd())
        call img_pad%zero_and_flag_ft
        call vol_pad%fproject(o, img_pad)
        call img_pad%shift([o%get('x'), o%get('y'), 0.])
        call simimg(img_pad, o, tfun, 'yes', 10., apply_ctf=.true.)
        call img%new([BOX,BOX,1], SMPD)
        call img_pad%clip(img)
        call vol_pad%kill_expanded
        call vol_pad%kill
        call img_pad%kill
    end subroutine flat_particle

    ! the wave mode's chain without noise: the map in volts through the operator to SIMPLE's particle units
    subroutine wave_particle( vol, o, model, img )
        class(image), intent(inout) :: vol, img
        class(ori),   intent(inout) :: o
        integer,      intent(in)    :: model
        class(wave_projector), allocatable :: wproj
        type(oris)           :: os
        complex, allocatable :: plane(:,:,:)
        real,    allocatable :: volts(:,:,:)
        call os%new(1, is_ptcl=.true.)
        call os%set_ori(1, o)
        volts = vol%get_rmat() / electron_scale_link(DOSE, KV, SMPD, BOX)
        allocate(wproj)
        call wproj%new(BOX, SMPD, model, 0, 1.0)
        call wproj%prep_particles(os)
        call wproj%set_volume(volts)
        allocate(plane(-BOX/2:BOX/2,-BOX/2:BOX/2,1))
        call img%new([BOX,BOX,1], SMPD)
        call wproj%predict(1, 1, plane)
        call wproj%plane_to_image(plane(:,:,1), img)
        call simimg_wave(img, DOSE, 10., shot_noise=.false.)
        call wproj%kill
        call os%kill
    end subroutine wave_particle

    !> a weak map (phase below 1e-3): the part of the wave particle odd in the map's sign is the flat particle at
    !! PCG's scale to the lattice's interpolation error. The flat chain projects the 2N-padded volume without
    !! deapodisation: fed the deapodised map, twice its particle is PCG's forward model, the scale link's reference
    subroutine test_weak_limit( model )
        integer, intent(in) :: model
        real, parameter :: AMP_LIN = 0.02, AMP_SO = 0.8
        type(image) :: vol, vol_deapod, flat
        type(ori)   :: o
        real, allocatable :: f(:,:,:), odd(:,:,:), even(:,:,:), deapod(:,:,:)
        real :: err_lin, ratio_even(2), phi_lin, phi_so
        integer :: iamp
        write(*,'(A,I0)') 'test_weak_limit model ', model
        call test_particle(o)
        call blob_map(vol, AMP_LIN)
        deapod = vol%get_rmat() * deapod_envelope()
        call vol_deapod%new([BOX,BOX,BOX], SMPD)
        call vol_deapod%set_rmat(deapod, .false.)
        call flat_particle(vol_deapod, o, flat)
        f = 2. * flat%get_rmat()
        call odd_even(vol, o, model, odd, even)
        err_lin = sqrt(sum((odd - f)**2) / sum(f**2))
        ! the even part is second order: its share of the odd part doubles with the map, up to third-order terms
        ! of relative size phi; at weaker maps single-precision rounding of the unit intensity hides it
        do iamp = 1, 2
            call blob_map(vol, real(iamp) * AMP_SO)
            call odd_even(vol, o, model, odd, even)
            ratio_even(iamp) = sqrt(sum(even**2) / sum(odd**2))
        end do
        phi_lin = max_phase(AMP_LIN)
        phi_so  = max_phase(2. * AMP_SO)
        write(*,'(A,ES10.3,A,ES10.3,A,2ES10.3,A,ES10.3,A)') '  linear part vs flat at PCG''s scale ', err_lin, &
            &' (phi ', phi_lin, '); even/odd ', ratio_even, ' (phi ', phi_so, ')'
        call assert_true(phi_lin < 1.e-3, 'weak limit: the map of the linear check is weak (phase below 1e-3)')
        call assert_real(0., err_lin, TOL_LATTICE, 'weak limit: the linear part is the flat particle at PCG''s scale')
        call assert_true(phi_so < 5.e-2, 'weak limit: the maps of the second-order check are weak (phase below 0.05)')
        call assert_real(2., ratio_even(2) / ratio_even(1), 2. * phi_so, 'weak limit: the even part is second order')
        call vol%kill; call vol_deapod%kill; call flat%kill; call o%kill
    end subroutine test_weak_limit

    ! the largest phase any ray picks up from the blob map of amplitude amp: three blobs on one line,
    ! each sqrt(2 pi) s V a with V the blob's amplitude in volts
    real function max_phase( amp )
        real, intent(in) :: amp
        max_phase = 3. * sqrt(2. * PI) * S_BLOB * SMPD * electron_interaction_constant(KV) * &
            &amp / electron_scale_link(DOSE, KV, SMPD, BOX)
    end function max_phase

    ! the parts of the wave particle odd and even in the map's sign
    subroutine odd_even( vol, o, model, odd, even )
        class(image),      intent(inout) :: vol
        class(ori),        intent(inout) :: o
        integer,           intent(in)    :: model
        real, allocatable, intent(out)   :: odd(:,:,:), even(:,:,:)
        type(image) :: wplus, wminus
        real, allocatable :: p(:,:,:), m(:,:,:)
        call wave_particle(vol, o, model, wplus)
        call vol%mul(-1.)
        call wave_particle(vol, o, model, wminus)
        call vol%mul(-1.)
        p    = wplus%get_rmat()
        m    = wminus%get_rmat()
        odd  = 0.5 * (p - m)
        even = 0.5 * (p + m)
        call wplus%kill; call wminus%kill
    end subroutine odd_even

    ! the inverse of the flat projector's Kaiser-Bessel envelope on the native box (2N padding)
    function deapod_envelope() result( env )
        real, allocatable :: env(:,:,:), inv1d(:)
        type(kbinterpol)  :: kbwin
        integer :: i, j, k
        kbwin = kbinterpol(KBWINSZ, KBALPHA)
        call kb_stencil_centered_crop_inv_envelope_1d(kbwin, 2*BOX, BOX, inv1d)
        allocate(env(BOX,BOX,BOX))
        do k = 1, BOX
            do j = 1, BOX
                do i = 1, BOX
                    env(i,j,k) = inv1d(i) * inv1d(j) * inv1d(k)
                end do
            end do
        end do
    end function deapod_envelope

    !> per-particle B-factors: particle 2 of a pair with B-factors 0 and B predicts bitwise as an operator of
    !! B-factor B for every particle, and differs from particle 1 at the same pose
    subroutine test_bfac_per_particle()
        real, parameter :: BFAC = 80.
        class(wave_projector), allocatable :: wp_pp, wp_one
        type(image)          :: vol
        type(ori)            :: o
        type(oris)           :: os2, os1
        complex, allocatable :: planes2(:,:,:), plane1(:,:,:)
        real,    allocatable :: volts(:,:,:)
        write(*,'(A)') 'test_bfac_per_particle'
        call test_particle(o)
        call os2%new(2, is_ptcl=.true.)
        call os2%set_ori(1, o)
        call os2%set_ori(2, o)
        call os1%new(1, is_ptcl=.true.)
        call os1%set_ori(1, o)
        call blob_map(vol, 0.8)
        volts = vol%get_rmat() / electron_scale_link(DOSE, KV, SMPD, BOX)
        allocate(planes2(-BOX/2:BOX/2,-BOX/2:BOX/2,2), plane1(-BOX/2:BOX/2,-BOX/2:BOX/2,1))
        allocate(wp_pp, wp_one)
        call wp_pp%new(BOX, SMPD, WAVE_PROJECTION, 0, 1.0)
        call wp_pp%prep_particles(os2, bfacs=[0., BFAC])
        call wp_pp%set_volume(volts)
        call wp_pp%predict(2, 1, planes2)
        call wp_one%new(BOX, SMPD, WAVE_PROJECTION, 0, 1.0, bfac=BFAC)
        call wp_one%prep_particles(os1)
        call wp_one%set_volume(volts)
        call wp_one%predict(1, 1, plane1)
        call assert_true(all(planes2(:,:,2) == plane1(:,:,1)), 'per-particle B-factor: as the operator''s B-factor')
        call assert_true(any(planes2(:,:,1) /= planes2(:,:,2)), 'per-particle B-factor: each particle its own')
        call wp_pp%kill; call wp_one%kill; call os2%kill; call os1%kill; call vol%kill; call o%kill
    end subroutine test_bfac_per_particle

    !> simimg_wave in an empty region: before the normalisation the counts have mean D a^2 and variance
    !! D a^2 + sdev^2, sdev^2 = var(E) / snr; checked to four standard errors over NIMG images
    subroutine test_detector_noise()
        integer, parameter :: NIMG = 60, NB = 64, NC = 12
        real,    parameter :: SNR = 0.5
        type(image) :: img
        real, allocatable :: rmat(:,:,:), corner(:)
        real(dp) :: s1, s2, mean_c, var_c, dose_px, sdev2, var_i, var_exp
        integer  :: n, iimg, i, j
        write(*,'(A)') 'test_detector_noise'
        call set_fixed_seed(20261024)
        dose_px = real(DOSE, dp) * real(SMPD, dp)**2
        call img%new([NB,NB,1], SMPD)
        allocate(corner(NIMG * 4 * NC * NC))
        n = 0
        do iimg = 1, NIMG
            ! relative intensity: 1 in the corners, a disk of contrast 0.3 at the centre
            allocate(rmat(NB,NB,1), source=1.)
            do j = 1, NB
                do i = 1, NB
                    if( (i - NB/2)**2 + (j - NB/2)**2 < 100 ) rmat(i,j,1) = 1.3
                end do
            end do
            if( iimg == 1 )then
                var_i = real(sum((rmat - sum(rmat) / size(rmat))**2) / size(rmat), dp)
            endif
            call img%set_rmat(rmat, .false.)
            call simimg_wave(img, DOSE, SNR)
            rmat = img%get_rmat()
            ! back to counts: x = -(c - D a^2) / sqrt(D a^2)
            rmat = real(dose_px - sqrt(dose_px) * real(rmat, dp))
            do j = 1, NB
                do i = 1, NB
                    if( (i <= NC .or. i > NB - NC) .and. (j <= NC .or. j > NB - NC) )then
                        n = n + 1
                        corner(n) = rmat(i,j,1)
                    endif
                end do
            end do
            deallocate(rmat)
        end do
        s1 = sum(real(corner(1:n), dp))
        s2 = sum(real(corner(1:n), dp)**2)
        mean_c  = s1 / n
        var_c   = s2 / n - mean_c**2
        sdev2   = var_i / real(SNR, dp) * dose_px**2
        var_exp = dose_px + sdev2
        write(*,'(A,F10.3,A,F10.3,A,F10.3,A,F10.3,A)') '  counts mean ', mean_c, ' (', dose_px, '), variance ', var_c, &
            &' (', var_exp, ')'
        call assert_real(real(dose_px), real(mean_c), real(4._dp * sqrt(var_exp / n)), 'detector: mean counts D a^2')
        call assert_real(real(var_exp), real(var_c), real(4._dp * var_exp * sqrt(2._dp / n)), &
            &'detector: variance D a^2 + sdev^2')
        call img%kill
    end subroutine test_detector_noise

end module simple_simulator_tester
