!@descr: unit tests of wave_projector: adjoint, gradients, weak-phase identity, curvature, model hierarchy, threads
! Box 24 at 1 A, 300 kV, astigmatic underfocus, potentials of random Gaussian blobs, particles
! at random poses and shifts.
! - The transpose of the volume-level linearisation (gather, shift, 2D chain): Re<J dV, q> = <dV, J^T q> to
!   1e-5, every model and both curvature signs, with point-group replication.
! - The gradient of the misfit against a fourth-order central difference to 1e-3, every model and both signs,
!   on a potential of about a radian of phase (the nonlinear regime), data predicted at V + dV plus noise.
! - The weak-phase identity: at V = 0 the linearisation is 2 g sigma a N times PCG's forward_plane times
!   build_transfer (CTF, shift) of the same deapodised volume, to 1e-5.
! - Rytov with no curvature is projection, bitwise; first Born is the linearisation of Rytov at V = 0; at V = 0
!   every model predicts the vacuum, relative intensity 1 (p = -g at the origin, 0 elsewhere), to 1e-5.
! - line_losses at x + t d equal the misfit of a volume moved to x + t d, to 1e-5.
! - A point scatterer z0 voxels off the centre plane on the paraboloid images as the flat model at defocus
!   df - kappa z0 (ewald=neg: df + z0): see the size and sign of the curvature.
! - predict, loss_grad, curvature, line_losses, linearize_planes, adjoint_planes and gain_terms are bitwise
!   identical at 1, 11, 16 and 24 threads, twice each.
! - Oversampling 1 (a measurement setting) keeps the exact transpose.
! - The scale link of reconstruct3D_wave: at V = 0 the Gauss-Newton curvature of the misfit along d,
!   sum w |J d|^2, is kappa^2 d^T H d / boxpd^3 with H PCG's (matrix-free) data normal operator and
!   kappa = 2 g sigma a N, to 1e-4 (both exact, single-precision lattice sums): a PCG map is kappa times volts,
!   and PCG's normal operator, density and prior carry boxpd^3 beside kappa^2 (its adjoint folds without the
!   1/boxpd^3 of the transform's transpose, which its solution does not see).
! - Multislice (the rotated cube gathered from the lattice): its adjoint and gradient as above, both signs; the
!   vacuum at zero potential; on a weak potential (about 1e-3 rad) its prediction tends to first Born's on the
!   paraboloid of the same sign, both signs, to 2.2e-2 of the largest term: two lattice interpolations (first
!   Born's at the paraboloid's fractional heights, multislice's cube at whole l with the transform along z), each
!   within the lattice's error as measured against specter (TOL_GOLDEN), and the second order, 1e-3.
! - The cube's z axis and the paraboloid agree: at 0.25 A pixels and 100 kV (the paraboloid 0.44 voxel off the
!   plane at Nyquist), a smooth blob 4 voxels downstream images under multislice (ewald=neg) as under first Born of
!   the same sign to 2.2e-2 (the weak limit's bound) and at least four times worse against the opposite sign, on
!   the half disk |k| <= R/2 (towards the rim the cube's truncation to the sphere enters, see the test). A flipped
!   axis would swap the two errors; the two signs themselves differ by only about 2% in this box, so the factor
!   stays well below the separation's (about 9).
! - Thread identity as above for multislice, whose gradient passes run in chunks of one particle per thread.
! - specter's golden images (simple_wave_projector_golden_tester) to the KB lattice's accuracy: see TOL_GOLDEN.
module simple_wave_projector_tester
use simple_defs,              only: dp, PI
use simple_type_defs,         only: ctfparams, CTFFLAG_YES
use simple_ori,               only: ori
use simple_oris,              only: oris
use simple_sym,               only: sym
use simple_math_ctf,          only: electron_interaction_constant
use simple_gridding,          only: kb_stencil_centered_crop_inv_envelope_1d
use simple_reconstructor_pcg, only: reconstructor_pcg
use simple_wave_imaging,      only: WAVE_PROJECTION, WAVE_FIRST_BORN, WAVE_RYTOV, WAVE_MSLICE
use simple_wave_projector,    only: wave_projector
use simple_image,             only: image
use simple_wave_projector_golden_tester, only: CENTRES16, CENTRES24, G16_PROJECTION_PAR, G16_PROJECTION_IMG, &
    &G16_FB_NEG_PAR, G16_FB_NEG_IMG, G16_FB_POS_PAR, G16_FB_POS_IMG, G16_RY_NEG_PAR, G16_RY_NEG_IMG, G16_RY_POS_PAR, &
    &G16_RY_POS_IMG, G24_PROJECTION_PAR, G24_PROJECTION_IMG, G24_RY_NEG_PAR, G24_RY_NEG_IMG, G24_FB_POS_PAR, G24_FB_POS_IMG
use simple_test_utils,        only: assert_true, fix_test_seed
!$ use omp_lib,                only: omp_get_max_threads, omp_set_num_threads
implicit none
private
public :: run_all_wave_projector_tests

integer, parameter :: BOX = 24, R = BOX/2, NP = 4
! the golden images' bounds, the plan's tolerance: the KB lattice's error, measured at zero curvature where
! both packages are exact by another route, with a margin of 1.25 for its spread over poses (+-20%).
! Three-tap window of the PCG lattice: the largest zero-curvature (projection) interior image error of the
! parity campaign (box 64, 96, 128: 8.2e-3) and of the golden projection cases here (8.77e-3); seven-tap
! window (winsz 3): the parity campaign's 4.9e-4
real,    parameter :: TOL_GOLDEN    = 1.10e-2
real,    parameter :: TOL_GOLDEN_W3 = 6.1e-4
real,    parameter :: SMPD = 1.0, KV = 300., GAIN = 6.0

contains

    subroutine run_all_wave_projector_tests()
        integer :: model, kappa
        write(*,'(A)') '**** running all wave projector tests ****'
        call fix_test_seed()
        do model = WAVE_PROJECTION, WAVE_RYTOV
            do kappa = -1, 1, 2
                if( model == WAVE_PROJECTION .and. kappa == 1 ) cycle
                call test_adjoint_and_gradient(model, kappa)
            end do
        end do
        call test_adjoint_and_gradient(WAVE_MSLICE, -1)
        call test_adjoint_and_gradient(WAVE_MSLICE,  1)
        call test_weak_phase_identity()
        call test_model_hierarchy()
        call test_line_losses()
        call test_point_scatterer()
        call test_thread_identity(WAVE_RYTOV)
        call test_thread_identity(WAVE_MSLICE)
        call test_oversampling_one()
        call test_golden()
        call test_pcg_scale()
        call test_multislice_weak_limit(-1)
        call test_multislice_weak_limit(1)
        call test_multislice_orientation()
    end subroutine run_all_wave_projector_tests

    ! ---- tests ----

    subroutine test_adjoint_and_gradient( model, kappa )
        integer, intent(in) :: model, kappa
        class(wave_projector), allocatable :: wp
        type(oris) :: os
        type(sym)  :: c2
        real,    allocatable :: v(:,:,:), dv(:,:,:), g(:,:,:)
        complex, allocatable :: q(:,:,:), jdv(:,:,:), y(:,:,:)
        real(dp) :: lhs, rhs, fd, an
        real     :: hstep
        allocate(wp)
        call wp%new(BOX, SMPD, model, kappa, GAIN)
        call c2%new('c2')
        call wp%set_sym(c2)
        os = particles(NP)
        call wp%prep_particles(os)
        ! blobs of about half a radian of phase each (300 V peak, one voxel wide): the nonlinear regime the
        ! difference step of wave imaging's tester was chosen for
        v  = blobs(300.0)
        dv = blobs(100.0)
        call wp%set_volume(v)
        call wp%set_direction(dv)
        ! adjoint
        allocate(jdv(-R:R,-R:R,NP), q(-R:R,-R:R,NP))
        call wp%linearize_planes(NP, 1, jdv)
        q = random_planes(NP)
        call wp%begin_grad
        call wp%adjoint_planes(q, NP, 1)
        allocate(g(BOX,BOX,BOX))
        call wp%get_grad(g)
        lhs = sum(real(conjg(jdv)*q, dp))
        rhs = sum(real(dv,dp)*real(g,dp))
        write(*,'(A,2I3,A,ES10.3)') '    model, kappa', model, kappa, ': adjoint relative error ', &
            &abs(lhs-rhs)/max(abs(lhs),abs(rhs))
        call assert_true(abs(lhs - rhs) <= 1.0e-5_dp * max(abs(lhs), abs(rhs)), &
            &'the volume gradient is the transpose of the volume linearisation')
        ! gradient of the misfit against a central difference, on data predicted at v + dv plus noise: the residual
        ! has a part along J dv, so the derivative along dv cannot cancel down to the misfit's rounding
        allocate(y(-R:R,-R:R,NP))
        call wp%set_volume(v + dv)
        call wp%predict(NP, 1, y)
        call wp%set_volume(v)
        y = y + 0.3 * maxval(abs(y(1:,:,:))) * random_planes(NP)
        call wp%begin_grad
        an = wp%loss_grad(y, NP, 1)
        call wp%get_grad(g)
        an = sum(real(dv,dp)*real(g,dp))
        ! fourth-order central difference at 3e-2 of the volume's scale (see wave imaging's tester)
        hstep = 3.0e-2 * maxval(abs(v)) / maxval(abs(dv))
        fd = 8.0_dp * (misfit3d(wp, y, v + hstep*dv) - misfit3d(wp, y, v - hstep*dv)) &
            &- (misfit3d(wp, y, v + 2.*hstep*dv) - misfit3d(wp, y, v - 2.*hstep*dv))
        fd = fd / (12.0_dp * real(hstep,dp))
        write(*,'(A,2I3,A,ES10.3)') '    model, kappa', model, kappa, ': gradient relative error ', abs(fd-an)/abs(an)
        call assert_true(abs(fd - an) <= 1.0e-3_dp * abs(an), 'the volume gradient of the misfit, against a difference')
        call wp%kill
        call os%kill
        deallocate(wp)
    end subroutine test_adjoint_and_gradient

    !> J(0) dV = 2 g sigma a N * build_transfer * forward_plane of the deapodised dV
    subroutine test_weak_phase_identity()
        class(wave_projector),    allocatable :: wp
        class(reconstructor_pcg), allocatable :: rp
        type(oris) :: os
        type(ori)  :: o
        real,    allocatable :: v0(:,:,:), dv(:,:,:), inv1d(:), dvd(:,:,:)
        complex, allocatable :: jdv(:,:,:), plane(:,:), tr(:,:), expected(:,:)
        type(ctfparams)  :: cp
        real    :: scale, err, errmax
        integer :: i, j, k, model
        os = particles(NP)
        dv = blobs(1.0)
        allocate(v0(BOX,BOX,BOX), source=0.)
        ! PCG's section of the same deapodised volume, with its transfer (CTF and shift)
        allocate(rp)
        call rp%new(BOX, SMPD, 0.)
        call rp%prep_particles(os, use_ctf=.true.)
        call kb_stencil_centered_crop_inv_envelope_1d(rp%kbwin, rp%boxpd, BOX, inv1d)
        allocate(dvd(BOX,BOX,BOX))
        do k = 1, BOX
            do j = 1, BOX
                do i = 1, BOX
                    dvd(i,j,k) = dv(i,j,k) * inv1d(i) * inv1d(j) * inv1d(k)
                end do
            end do
        end do
        call rp%set_volume(dvd)
        allocate(plane(-R:R,-R:R), tr(-R:R,-R:R), expected(-R:R,-R:R), jdv(-R:R,-R:R,NP))
        scale = 2. * GAIN * electron_interaction_constant(KV) * SMPD * real(BOX)
        errmax = 0.
        do model = WAVE_PROJECTION, WAVE_RYTOV
            allocate(wp)
            call wp%new(BOX, SMPD, model, 0, GAIN)
            call wp%prep_particles(os)
            call wp%set_volume(v0)
            call wp%set_direction(dv)
            call wp%linearize_planes(NP, 1, jdv)
            do i = 1, NP
                call os%get_ori(i, o)
                cp = o%get_ctfvars()
                call rp%forward_plane(o, plane)
                tr = rp%build_transfer(cp, o%get_2Dshift())
                expected = scale * tr * plane
                err = maxval(abs(jdv(:,:,i) - expected)) / maxval(abs(expected))
                errmax = max(errmax, err)
            end do
            call wp%kill
            deallocate(wp)
        end do
        write(*,'(A,ES10.3)') '    weak-phase identity against forward_plane and build_transfer: ', errmax
        call assert_true(errmax <= 1.0e-5, 'the linearisation at zero is 2 g sigma a N CTF S of the PCG backend')
        call rp%kill
        call os%kill
        call o%kill
    end subroutine test_weak_phase_identity

    subroutine test_model_hierarchy()
        class(wave_projector), allocatable :: wr, wpj, wb
        type(oris) :: os
        real,    allocatable :: v(:,:,:), v0(:,:,:), dv(:,:,:)
        complex, allocatable :: pr(:,:,:), pp(:,:,:), jr(:,:,:), jb(:,:,:)
        os = particles(NP)
        v  = blobs(3.0)
        dv = blobs(1.0)
        allocate(v0(BOX,BOX,BOX), source=0.)
        allocate(pr(-R:R,-R:R,NP), pp(-R:R,-R:R,NP), jr(-R:R,-R:R,NP), jb(-R:R,-R:R,NP))
        allocate(wr, wpj, wb)
        call wr%new(BOX, SMPD, WAVE_RYTOV, 0, GAIN)
        call wpj%new(BOX, SMPD, WAVE_PROJECTION, 0, GAIN)
        call wr%prep_particles(os)
        call wpj%prep_particles(os)
        call wr%set_volume(v)
        call wpj%set_volume(v)
        call wr%predict(NP, 1, pr)
        call wpj%predict(NP, 1, pp)
        call assert_true(all(pr == pp), 'Rytov without curvature is projection, bitwise')
        ! the vacuum: relative intensity 1, SIMPLE's transform 1 at the origin, p = -g there
        call wpj%set_volume(v0)
        call wpj%predict(NP, 1, pp)
        pp(0,0,:) = pp(0,0,:) + GAIN
        call assert_true(maxval(abs(pp)) <= 1.0e-5 * GAIN, 'at zero potential the prediction is the vacuum, intensity 1')
        ! multislice too
        call wr%kill
        call wr%new(BOX, SMPD, WAVE_MSLICE, -1, GAIN)
        call wr%prep_particles(os)
        call wr%set_volume(v0)
        call wr%predict(NP, 1, pr)
        pr(0,0,:) = pr(0,0,:) + GAIN
        call assert_true(maxval(abs(pr)) <= 1.0e-5 * GAIN, 'at zero potential multislice predicts the vacuum')
        call wr%kill
        call wr%new(BOX, SMPD, WAVE_RYTOV, -1, GAIN)
        call wb%new(BOX, SMPD, WAVE_FIRST_BORN, -1, GAIN)
        call wr%prep_particles(os)
        call wb%prep_particles(os)
        call wr%set_volume(v0)
        call wb%set_volume(v0)
        call wr%set_direction(dv)
        call wb%set_direction(dv)
        call wr%linearize_planes(NP, 1, jr)
        call wb%linearize_planes(NP, 1, jb)
        call assert_true(maxval(abs(jr - jb)) <= 1.0e-6 * maxval(abs(jr)), &
            &'first Born is the linearisation of Rytov, on the paraboloid')
        call wr%kill
        call wb%kill
        call wpj%kill
        call os%kill
    end subroutine test_model_hierarchy

    subroutine test_line_losses()
        class(wave_projector), allocatable :: wp
        type(oris) :: os
        real,    allocatable :: v(:,:,:), dv(:,:,:)
        complex, allocatable :: y(:,:,:)
        real(dp) :: t(3), losses(3), direct
        real     :: err
        integer  :: j
        os = particles(NP)
        v  = blobs(3.0)
        dv = blobs(1.0)
        allocate(wp)
        call wp%new(BOX, SMPD, WAVE_RYTOV, -1, GAIN)
        call wp%prep_particles(os)
        call wp%set_volume(v)
        allocate(y(-R:R,-R:R,NP))
        call wp%predict(NP, 1, y)
        y = y + 0.3 * maxval(abs(y(1:,:,:))) * random_planes(NP)
        call wp%set_direction(dv)
        t = [0.0_dp, 0.5_dp, 1.0_dp]
        call wp%line_losses(y, NP, 1, t, losses)
        err = 0.
        do j = 1, 3
            call wp%set_volume(v + real(t(j)) * dv)
            call wp%begin_grad
            direct = wp%loss_grad(y, NP, 1)
            err = max(err, real(abs(direct - losses(j)) / direct))
        end do
        call assert_true(err <= 1.0e-5, 'line_losses are the misfits of the moved volume')
        call wp%kill
        call os%kill
    end subroutine test_line_losses

    !> a point z0 voxels along the beam (identity pose, no shift): on the paraboloid, first Born, it
    !! images as the flat first Born at defocus df - kappa z0; the other sign is far off
    subroutine test_point_scatterer()
        integer, parameter :: Z0 = 4
        class(wave_projector), allocatable :: wc, wf
        type(oris) :: os, osf
        real,    allocatable :: v(:,:,:)
        complex, allocatable :: pc(:,:,:), pf(:,:,:)
        real    :: err(2), dfz
        integer :: is
        allocate(v(BOX,BOX,BOX), source=0.)
        v(BOX/2+1, BOX/2+1, BOX/2+1+Z0) = 50.
        call os%new(1, .true.)
        call set_particle(os, 1, [0.,0.,0.], [0.,0.], 1.0)
        allocate(pc(-R:R,-R:R,1), pf(-R:R,-R:R,1))
        allocate(wc, wf)
        call wc%new(BOX, SMPD, WAVE_FIRST_BORN, -1, GAIN)
        call wc%prep_particles(os)
        call wc%set_volume(v)
        call wc%predict(1, 1, pc)
        do is = 1, 2
            ! is = 1: df + z0 (what ewald=neg must give), is = 2: df - z0
            dfz = 1.0 + merge(1., -1., is == 1) * real(Z0) * SMPD * 1.0e-4
            v = 0.
            v(BOX/2+1, BOX/2+1, BOX/2+1) = 50.
            call osf%new(1, .true.)
            call set_particle(osf, 1, [0.,0.,0.], [0.,0.], dfz)
            call wf%new(BOX, SMPD, WAVE_FIRST_BORN, 0, GAIN)
            call wf%prep_particles(osf)
            call wf%set_volume(v)
            call wf%predict(1, 1, pf)
            pc(0,0,1) = 0.; pf(0,0,1) = 0.
            err(is) = maxval(abs(pc - pf)) / maxval(abs(pf))
            call wf%kill
            call osf%kill
        end do
        write(*,'(A,2ES10.3)') '    point off the centre plane, ewald=neg against df+z0 and df-z0: ', err
        call assert_true(err(1) <= 2.0e-2, 'ewald=neg images a point at height z0 at defocus df + z0')
        call assert_true(err(2) >= 10. * err(1), 'and not at df - z0')
        call wc%kill
        call os%kill
    end subroutine test_point_scatterer

    subroutine test_thread_identity( model )
        integer, intent(in) :: model
        integer, parameter :: NTHR(4) = [1, 11, 16, 24]
        class(wave_projector), allocatable :: wp
        type(oris) :: os
        real,    allocatable :: v(:,:,:), g(:,:,:), gref(:,:,:), ga(:,:,:), garef(:,:,:)
        complex, allocatable :: y(:,:,:), p(:,:,:), pref(:,:,:), jd(:,:,:), jdref(:,:,:)
        real,    allocatable :: dv(:,:,:)
        real(dp) :: loss, lref, curv, cref, ll(2), llref(2), gnum(2*NP+3), gden(2*NP+3), gnref(2*NP+3), gdref(2*NP+3)
        integer  :: it, irep, nthr0
        logical  :: l_same
        os = particles(2*NP+3)
        v  = blobs(3.0)
        dv = blobs(1.0)
        allocate(wp)
        call wp%new(BOX, SMPD, model, -1, GAIN)
        call wp%prep_particles(os)
        call wp%set_volume(v)
        call wp%set_direction(dv)
        allocate(y(-R:R,-R:R,2*NP+3), p(-R:R,-R:R,2*NP+3), pref(-R:R,-R:R,2*NP+3))
        y = random_planes(2*NP+3)
        allocate(jd(-R:R,-R:R,2*NP+3), jdref(-R:R,-R:R,2*NP+3))
        allocate(g(BOX,BOX,BOX), gref(BOX,BOX,BOX), ga(BOX,BOX,BOX), garef(BOX,BOX,BOX))
        nthr0 = 1
        !$ nthr0 = omp_get_max_threads()
        l_same = .true.
        do it = 1, size(NTHR)
            do irep = 1, 2
                !$ call omp_set_num_threads(NTHR(it))
                call wp%predict(2*NP+3, 1, p)
                call wp%begin_grad
                loss = wp%loss_grad(y, 2*NP+3, 1)
                call wp%get_grad(g)
                curv = wp%curvature(2*NP+3, 1)
                call wp%line_losses(y, 2*NP+3, 1, [0.5_dp, 1.0_dp], ll)
                call wp%linearize_planes(2*NP+3, 1, jd)
                call wp%begin_grad
                call wp%adjoint_planes(y, 2*NP+3, 1)
                call wp%get_grad(ga)
                call wp%gain_terms(y, 2*NP+3, 1, gnum, gden)
                if( it == 1 .and. irep == 1 )then
                    pref = p; gref = g; lref = loss; cref = curv; llref = ll; jdref = jd; garef = ga
                    gnref = gnum; gdref = gden
                else
                    l_same = l_same .and. all(gnum == gnref) .and. all(gden == gdref)
                    l_same = l_same .and. all(p == pref) .and. all(g == gref) .and. loss == lref
                    l_same = l_same .and. curv == cref .and. all(ll == llref)
                    l_same = l_same .and. all(jd == jdref) .and. all(ga == garef)
                endif
            end do
        end do
        !$ call omp_set_num_threads(nthr0)
        call assert_true(l_same, &
            &'every output of the operator is identical at 1, 11, 16 and 24 threads')
        call wp%kill
        call os%kill
    end subroutine test_thread_identity

    !> the image formation on the native grid itself (oversampling 1, for measurement): still an exact
    !! transpose, the disk's Nyquist pairs folded onto their bins
    subroutine test_oversampling_one()
        ! a potential of a few volts: the test is of the transpose's structure on the coarse grid; at 300 V
        ! the unoversampled chain's single-precision rounding reaches 1.2e-5 (4e-7 at 30 V, 1.7e-7 at 3 V),
        ! a reason more to keep the oversampling at 2
        real, parameter :: VAMP = 3.
        class(wave_projector), allocatable :: wp
        type(oris) :: os
        real,    allocatable :: v(:,:,:), dv(:,:,:), g(:,:,:)
        complex, allocatable :: q(:,:,:), jdv(:,:,:)
        real(dp) :: lhs, rhs
        os = particles(NP)
        v  = blobs(VAMP)
        dv = blobs(VAMP/3.)
        allocate(wp)
        call wp%new(BOX, SMPD, WAVE_RYTOV, -1, GAIN, osf=1)
        call wp%prep_particles(os)
        call wp%set_volume(v)
        call wp%set_direction(dv)
        allocate(jdv(-R:R,-R:R,NP), g(BOX,BOX,BOX))
        call wp%linearize_planes(NP, 1, jdv)
        q = random_planes(NP)
        call wp%begin_grad
        call wp%adjoint_planes(q, NP, 1)
        call wp%get_grad(g)
        lhs = sum(real(conjg(jdv)*q, dp))
        rhs = sum(real(dv,dp)*real(g,dp))
        call assert_true(abs(lhs - rhs) <= 1.0e-5_dp * max(abs(lhs), abs(rhs)), &
            &'without oversampling the gradient is still the transpose of the linearisation')
        call wp%kill
        call os%kill
    end subroutine test_oversampling_one

    !> specter's golden images (simple_wave_projector_golden_tester): the interior relative intensity of the same
    !! tapered blob potential, poses and CTF, per particle, within TOL_GOLDEN
    subroutine test_golden()
        real    :: err(8)
        integer :: iw
        do iw = 1, 2
            ! iw = 1: the PCG lattice's KB window; iw = 2: seven taps, where the operator's own numerics show
            call golden_case(16, 300., CENTRES16, WAVE_PROJECTION,  0, G16_PROJECTION_PAR, G16_PROJECTION_IMG, iw, err(1))
            call golden_case(16, 300., CENTRES16, WAVE_FIRST_BORN, -1, G16_FB_NEG_PAR,     G16_FB_NEG_IMG,     iw, err(2))
            call golden_case(16, 300., CENTRES16, WAVE_FIRST_BORN,  1, G16_FB_POS_PAR,     G16_FB_POS_IMG,     iw, err(3))
            call golden_case(16, 300., CENTRES16, WAVE_RYTOV,      -1, G16_RY_NEG_PAR,     G16_RY_NEG_IMG,     iw, err(4))
            call golden_case(16, 300., CENTRES16, WAVE_RYTOV,       1, G16_RY_POS_PAR,     G16_RY_POS_IMG,     iw, err(5))
            call golden_case(24, 100., CENTRES24, WAVE_PROJECTION,  0, G24_PROJECTION_PAR, G24_PROJECTION_IMG, iw, err(6))
            call golden_case(24, 100., CENTRES24, WAVE_RYTOV,      -1, G24_RY_NEG_PAR,     G24_RY_NEG_IMG,     iw, err(7))
            call golden_case(24, 100., CENTRES24, WAVE_FIRST_BORN,  1, G24_FB_POS_PAR,     G24_FB_POS_IMG,     iw, err(8))
            write(*,'(A,I2,A,8ES10.2)') '    golden cases against specter, window ', iw, ', relative interior image error: ', err
            if( iw == 1 )then
                call assert_true(all(err <= TOL_GOLDEN), 'the predicted images are specter''s, within the lattice''s accuracy')
            else
                call assert_true(all(err <= TOL_GOLDEN_W3), 'with a seven-tap window, within its accuracy')
            endif
        end do
    end subroutine test_golden

    subroutine golden_case( n, kv_case, centres, model, kappa, par, img, iw, err )
        integer, intent(in)  :: n, model, kappa, iw
        real,    intent(in)  :: kv_case, centres(3,5), par(10,2), img(n,n,2)
        real,    intent(out) :: err
        class(wave_projector), allocatable :: wp
        type(oris)      :: os
        type(ori)       :: o
        type(ctfparams) :: cp
        type(image)     :: vol, pimg
        real,    allocatable :: v(:,:,:), ours(:,:,:)
        complex, allocatable :: p(:,:,:)
        complex :: comp
        real    :: x, y, z, r, lo, hi, t, ddf
        integer :: i, j, k, b, h, l, m, rn, phys(3), ip
        ! the blobs, then the parity taper between 0.70 and 0.85 of the native Nyquist radius
        allocate(v(n,n,n), source=0.)
        do k = 1, n
            do j = 1, n
                do i = 1, n
                    x = real(i - 1 - n/2); y = real(j - 1 - n/2); z = real(k - 1 - n/2)
                    do b = 1, 5
                        v(i,j,k) = v(i,j,k) + 3. * exp(-((x-centres(1,b))**2 + (y-centres(2,b))**2 + (z-centres(3,b))**2) / 2.)
                    end do
                end do
            end do
        end do
        call vol%new([n,n,n], 1.0)
        call vol%set_rmat(v, .false.)
        call vol%fft()
        lo = 0.70 * real(n/2); hi = 0.85 * real(n/2)
        do m = -n/2, n/2 - 1
            do l = -n/2, n/2 - 1
                do h = 0, n/2
                    r    = sqrt(real(h*h + l*l + m*m))
                    t    = min(1., max(0., (r - lo) / (hi - lo)))
                    phys = vol%comp_addr_phys(h,l,m)
                    comp = vol%get_fcomp([h,l,m], phys)
                    call vol%set_fcomp([h,l,m], phys, comp * (0.5 + 0.5 * cos(PI * t)))
                end do
            end do
        end do
        call vol%ifft()
        v = vol%get_rmat()
        call vol%kill
        ! the particles; specter's positive curvature refers its defocus one voxel off the box centre
        ddf = 0.
        if( kappa == 1 ) ddf = -1.0e-4
        call os%new(2, .true.)
        do ip = 1, 2
            call o%new(.true.)
            call o%set_euler(par(1:3,ip))
            call o%set_shift(par(4:5,ip))
            cp%ctfflag = CTFFLAG_YES; cp%smpd = 1.0; cp%kv = kv_case; cp%cs = par(9,ip); cp%fraca = 0.07
            cp%dfx = par(6,ip) + ddf; cp%dfy = par(7,ip) + ddf; cp%angast = par(8,ip); cp%phshift = par(10,ip)
            call o%set_ctfvars(cp)
            call os%set_ori(ip, o)
            call o%kill
        end do
        allocate(wp)
        if( iw == 1 )then
            call wp%new(n, 1.0, model, kappa, 1.0)
        else
            call wp%new(n, 1.0, model, kappa, 1.0, winsz=3.0)
        endif
        call wp%prep_particles(os)
        call wp%set_volume(v)
        rn = n / 2
        allocate(p(-rn:rn,-rn:rn,2), ours(n,n,2))
        call wp%predict(2, 1, p)
        ! the interior relative intensity minus one: -p (gain 1) for 0 < |k| < R - 2
        call pimg%new([n,n,1], 1.0)
        err = 0.
        do ip = 1, 2
            call pimg%zero_and_flag_ft()
            do k = -rn, rn
                do h = 0, rn
                    if( h*h + k*k == 0 .or. real(h*h + k*k) >= real(rn-2)**2 ) cycle
                    phys = pimg%comp_addr_phys(h,k,0)
                    call pimg%set_fcomp([h,k,0], phys, -p(h,k,ip))
                end do
            end do
            call pimg%ifft()
            ours = 0.
            ours(:,:,ip:ip) = pimg%get_rmat()
            err = max(err, sqrt(sum((ours(:,:,ip) - img(:,:,ip))**2) / sum(img(:,:,ip)**2)))
        end do
        call pimg%kill
        call wp%kill
        call os%kill
    end subroutine golden_case

    ! ---- helpers ----

    real(dp) function misfit3d( wp, y, v )
        class(wave_projector), intent(inout) :: wp
        complex,               intent(in)    :: y(-R:R,-R:R,NP)
        real,                  intent(in)    :: v(BOX,BOX,BOX)
        call wp%set_volume(v)
        call wp%begin_grad
        misfit3d = wp%loss_grad(y, NP, 1)
    end function misfit3d

    subroutine test_pcg_scale()
        class(wave_projector),    allocatable :: wp
        class(reconstructor_pcg), allocatable :: rp
        type(oris) :: os
        real, allocatable :: v0(:,:,:), b1(:,:,:), b2(:,:,:), d(:,:,:), hd(:,:,:), inv1d(:), invenv(:,:,:)
        real(dp) :: curv_wave, dhd, kappa, err
        integer  :: i, j, k
        os = particles(NP)
        ! d with a zero sum of d E^-1: no section of it has a DC term (the wave misfit leaves the origin out)
        allocate(rp)
        call rp%new(BOX, SMPD, 0.)
        call rp%prep_particles(os, use_ctf=.true.)
        call kb_stencil_centered_crop_inv_envelope_1d(rp%kbwin, rp%boxpd, BOX, inv1d)
        allocate(invenv(BOX,BOX,BOX))
        do k = 1, BOX
            do j = 1, BOX
                do i = 1, BOX
                    invenv(i,j,k) = inv1d(i) * inv1d(j) * inv1d(k)
                end do
            end do
        end do
        b1 = blobs(1.0)
        b2 = blobs(1.0)
        d  = b1 - (sum(b1 * invenv) / sum(b2 * invenv)) * b2
        hd = rp%apply_normal(d)
        dhd = sum(real(d,dp) * real(hd,dp))
        allocate(v0(BOX,BOX,BOX), source=0.)
        allocate(wp)
        call wp%new(BOX, SMPD, WAVE_FIRST_BORN, 0, GAIN)
        call wp%prep_particles(os)
        call wp%set_volume(v0)
        call wp%set_direction(d)
        curv_wave = wp%curvature(NP, 1)
        kappa = 2.0_dp * real(GAIN,dp) * real(electron_interaction_constant(KV),dp) * real(SMPD,dp) * real(BOX,dp)
        err   = abs(curv_wave - kappa**2 * dhd / real(rp%boxpd,dp)**3) / abs(curv_wave)
        write(*,'(A,ES10.3)') '    Gauss-Newton curvature at zero against kappa^2 times PCG''s normal operator: ', err
        call assert_true(err <= 1.0e-4_dp, 'a PCG map is kappa = 2 g sigma a N times the potential in volts')
        call wp%kill
        call rp%kill
        call os%kill
    end subroutine test_pcg_scale

    subroutine test_multislice_weak_limit( kappa )
        integer, intent(in) :: kappa
        class(wave_projector), allocatable :: wms, wfb
        type(oris) :: os
        complex, allocatable :: pms(:,:,:), pfb(:,:,:)
        real,    allocatable :: v(:,:,:)
        real :: err
        os = particles(NP)
        ! a weak, smooth potential: about 1e-3 rad through a blob
        v = blobs(0.5, 2.0)
        allocate(wms, wfb)
        call wms%new(BOX, SMPD, WAVE_MSLICE,     kappa, GAIN)
        call wfb%new(BOX, SMPD, WAVE_FIRST_BORN, kappa, GAIN)
        call wms%prep_particles(os)
        call wfb%prep_particles(os)
        call wms%set_volume(v)
        call wfb%set_volume(v)
        allocate(pms(-R:R,-R:R,NP), pfb(-R:R,-R:R,NP))
        call wms%predict(NP, 1, pms)
        call wfb%predict(NP, 1, pfb)
        ! the vacuum term at the origin differs by the absorption's second order; the contrast is the rest
        pms(0,0,:) = cmplx(0.,0.)
        pfb(0,0,:) = cmplx(0.,0.)
        err = maxval(abs(pms - pfb)) / maxval(abs(pfb))
        write(*,'(A,I3,A,ES10.3)') '    multislice of a weak potential against first Born, kappa', kappa, ': ', err
        call assert_true(err <= 2.2e-2, 'multislice tends to first Born on the paraboloid of its sign')
        call wms%kill
        call wfb%kill
        call os%kill
    end subroutine test_multislice_weak_limit

    subroutine test_multislice_orientation()
        real, parameter :: SMPD_O = 0.25, KV_O = 100.
        class(wave_projector), allocatable :: wms, wfb
        type(oris)      :: os
        type(ori)       :: o
        type(ctfparams) :: cp
        complex, allocatable :: pms(:,:,:), pfb(:,:,:)
        real,    allocatable :: v(:,:,:)
        real    :: rnd(3), err_same, err_opp, d2
        integer :: i, j, k, kappa
        call os%new(NP, .true.)
        do i = 1, NP
            call random_number(rnd)
            call o%new(.true.)
            call o%set_euler([360.*rnd(1), 30.*rnd(2), 360.*rnd(3)])
            call o%set_shift([0.,0.])
            cp%ctfflag = CTFFLAG_YES; cp%smpd = SMPD_O; cp%kv = KV_O; cp%cs = 2.7; cp%fraca = 0.1
            cp%dfx = 0.3; cp%dfy = 0.27; cp%angast = 40.; cp%phshift = 0.
            call o%set_ctfvars(cp)
            call os%set_ori(i, o)
            call o%kill
        end do
        ! a weak, smooth blob (2 voxels) four voxels downstream of the centre, off the axis, within the central
        ! half of the box: a sharp one has content near the Nyquist sphere, where the cube is truncated and the
        ! paraboloid's fractional heights are interpolated from few samples along l
        allocate(v(BOX,BOX,BOX))
        do k = 1, BOX
            do j = 1, BOX
                do i = 1, BOX
                    d2 = real(i-BOX/2-1-2)**2 + real(j-BOX/2-1+1)**2 + real(k-BOX/2-1-4)**2
                    v(i,j,k) = 0.5 * exp(-d2 / 8.)
                end do
            end do
        end do
        allocate(pms(-R:R,-R:R,NP), pfb(-R:R,-R:R,NP))
        allocate(wms, wfb)
        call wms%new(BOX, SMPD_O, WAVE_MSLICE, -1, GAIN)
        call wms%prep_particles(os)
        call wms%set_volume(v)
        call wms%predict(NP, 1, pms)
        ! the half disk, |k| <= R/2: there the cube holds at least 0.87 R along l, so its truncation to the sphere
        ! (towards the rim the paraboloid's fractional heights are interpolated from few samples along l, which first
        ! Born reads from the whole lattice) cannot enter, while the two signs differ by up to 0.4 rad there
        call interior(pms)
        do kappa = -1, 1, 2
            call wfb%new(BOX, SMPD_O, WAVE_FIRST_BORN, kappa, GAIN)
            call wfb%prep_particles(os)
            call wfb%set_volume(v)
            call wfb%predict(NP, 1, pfb)
            call interior(pfb)
            if( kappa == -1 )then
                err_same = maxval(abs(pms - pfb)) / maxval(abs(pfb))
            else
                err_opp  = maxval(abs(pms - pfb)) / maxval(abs(pfb))
            endif
            call wfb%kill
        end do
        write(*,'(A,2ES10.3)') '    multislice (ewald=neg) against first Born of the same and the opposite sign: ', err_same, err_opp
        call assert_true(err_same <= 2.2e-2 .and. err_opp >= 4. * err_same, &
            &'the multislice cube''s z axis is the paraboloid''s: the curvature sign of the data')
        call wms%kill
        call os%kill

    contains

        subroutine interior( p )
            complex, intent(inout) :: p(-R:R,-R:R,NP)
            integer :: h, k
            do k = -R, R
                do h = -R, R
                    if( 4*(h*h + k*k) > R*R .or. (h == 0 .and. k == 0) ) p(h,k,:) = cmplx(0.,0.)
                end do
            end do
        end subroutine interior
    end subroutine test_multislice_orientation

    !> n particles at random poses, shifts up to 2 pixels, astigmatic underfocus 0.8-1.6 um
    function particles( n ) result( os )
        integer, intent(in) :: n
        type(oris) :: os
        real    :: r(6)
        integer :: i
        call os%new(n, .true.)
        do i = 1, n
            call random_number(r)
            call set_particle(os, i, [360.*r(1), 180.*r(2), 360.*r(3)], 4.*(r(4:5)-0.5), 0.8 + 0.8*r(6))
        end do
    end function particles

    subroutine set_particle( os, i, euls, sh, df )
        type(oris), intent(inout) :: os
        integer,    intent(in)    :: i
        real,       intent(in)    :: euls(3), sh(2), df
        type(ori)       :: o
        type(ctfparams) :: cp
        call o%new(.true.)
        call o%set_euler(euls)
        call o%set_shift(sh)
        cp%ctfflag = CTFFLAG_YES
        cp%smpd    = SMPD
        cp%kv      = KV
        cp%cs      = 2.7
        cp%fraca   = 0.1
        cp%dfx     = df
        cp%dfy     = 0.9 * df
        cp%angast  = 40.
        cp%phshift = 0.
        call o%set_ctfvars(cp)
        call os%set_ori(i, o)
        call o%kill
    end subroutine set_particle

    !> five Gaussian blobs of the given peak (volts) and standard deviation (voxels, 1 by default) at random
    !! places within the central half
    function blobs( peak, sdev ) result( v )
        real,           intent(in) :: peak
        real, optional, intent(in) :: sdev
        real, allocatable :: v(:,:,:)
        real    :: c(3,5), d2, s2
        integer :: i, j, k, b
        s2 = 1.
        if( present(sdev) ) s2 = sdev**2
        call random_number(c)
        c = real(BOX/2+1) + (c - 0.5) * real(BOX/2)
        allocate(v(BOX,BOX,BOX), source=0.)
        do k = 1, BOX
            do j = 1, BOX
                do i = 1, BOX
                    do b = 1, 5
                        d2 = (real(i)-c(1,b))**2 + (real(j)-c(2,b))**2 + (real(k)-c(3,b))**2
                        v(i,j,k) = v(i,j,k) + peak * exp(-d2 / (2. * s2))
                    end do
                end do
            end do
        end do
    end function blobs

    function random_planes( n ) result( q )
        integer, intent(in) :: n
        complex, allocatable :: q(:,:,:)
        real,    allocatable :: re(:,:,:), im(:,:,:)
        integer :: h, k
        allocate(re(-R:R,-R:R,n), im(-R:R,-R:R,n))
        call random_number(re)
        call random_number(im)
        allocate(q(-R:R,-R:R,n))
        q(:,:,:) = cmplx(re - 0.5, im - 0.5)
        do k = -R, R
            do h = -R, R
                if( h*h + k*k > R*R ) q(h,k,:) = cmplx(0.,0.)
            end do
        end do
    end function random_planes

end module simple_wave_projector_tester
