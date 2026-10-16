!@descr: unit tests for continuous in-plane registration on the polar Fourier transform (simple_polarft_calc, simple_pftc_shsrch_grad)
! A phantom volume (written to the run directory, removed at the end) projected through the strategy3D
! toolbox gives a reference at e3 = 0 and a particle rotated and shifted by a known amount. For the raw
! Euclidean and cc continuous-angle evaluators: identity with the discrete routes at integer
! angles, analytic vs central-difference gradients, non-negative losses under near-noiseless sigma2, the cc
! penalty, the joint seed selection, recovery of the known pose, and the strategy2D routes under inpl_cont.
module simple_pftc_inplane_tester
use simple_pftc_srch_api
use simple_string,                  only: string
use simple_syslib,                  only: del_file
use simple_builder,                 only: builder
use simple_image,                   only: image
use simple_matcher_smpl_and_lplims, only: set_bp_range3D
use simple_polarft_calc,            only: polarft_calc, vol_pad2ref_pfts
use simple_pftc_shsrch_grad,        only: pftc_shsrch_grad, fit_quadratic_step, grid_stencil_exists, FIT_CURVATURE, FIT_STEP, &
    &shift_ptcl2ref, shift_ref2ptcl
use simple_strategy2D_srch,         only: strategy2D_srch, strategy2D_spec
use simple_test_utils
use, intrinsic :: ieee_arithmetic,  only: ieee_is_finite, ieee_value, ieee_quiet_nan
use, intrinsic :: iso_fortran_env,  only: int64
implicit none
private
public :: run_all_pftc_inplane_tests

! the phantom: box 64 at 1.3 A, four Gaussian blobs inside a 60 A mask
integer,          parameter :: BOX     = 64
real,             parameter :: SMPD    = 1.3
real,             parameter :: MSKDIAM = 60.
real,             parameter :: LP_LOW  = 8.
real,             parameter :: LP_FULL = 2.7
real,             parameter :: TRS     = 5.
character(len=*), parameter :: PHANTOM_FILE = 'pftc_inplane_phantom.mrc'
integer,          parameter :: NBLOBS = 4
real,             parameter :: CTRS(3,NBLOBS) = reshape([&
    &-9.0, -5.0,  3.0, &
    & 7.0,  8.0, -5.0, &
    & 0.0,-10.0, -8.0, &
    & 5.0, -3.0, 10.0], [3,NBLOBS])
real,             parameter :: SIGMAS(NBLOBS) = [3.5, 4.0, 3.0, 3.8]
real,             parameter :: AMPS(NBLOBS)   = [1.0, 0.8, 0.6, 0.5]
! the truth pose of the particle relative to the reference
real,             parameter :: TRUTH_ANGLE    = 37.
real,             parameter :: TRUTH_SHIFT(2) = [2., -1.5]
! FD noise floor: single-precision coefficient series differentiated with h=1e-3
! leave ~1e-3 absolute noise; a broken gradient errs at O(|grad|)
real(dp),         parameter :: FD_STEP        = 1.d-3
real(dp),         parameter :: GRAD_FD_RTOL   = 1.d-2
real(dp),         parameter :: GRAD_FD_AFLOOR = 3.d-3
real(dp),         parameter :: IDENTITY_ATOL  = 1.d-4
real(dp),         parameter :: PROBE_SHIFTS(2,3) = reshape([0.d0,0.d0, 1.7d0,-2.3d0, -0.6d0,0.4d0], [2,3])
real(dp),         parameter :: PROBE_ROTOFFS(4)  = [0.37d0, -1.42d0, 0.5d0, 1.93d0]

integer, parameter :: OBJ_EUCLID = 1, OBJ_CC = 2

! the pftc keeps a reference to the sigma2 array it is given: module lifetime
real(sp), allocatable, target :: sigma2_fixture(:,:)

contains

    subroutine run_all_pftc_inplane_tests()
        write(*,'(A)') '**** running all pftc in-plane tests ****'
        call test_grid_fit_step()
        call test_shift_frames()
        call write_phantom()
        call test_euclid_routes_and_recovery()
        call test_euclid_band_edges()
        call test_cc_evaluator()
        call test_strategy2D_route_flags()
        call test_grid_fit_search()
        call test_particle_bank()
        call test_grid_ptcl_search()
        call del_file(PHANTOM_FILE)
    end subroutine run_all_pftc_inplane_tests

    ! ---- fixture ---------------------------------------------------------------

    subroutine write_phantom()
        type(image) :: vol
        real, allocatable :: rmat(:,:,:)
        real    :: ctr, dx, dy, dz
        integer :: b, i, j, k
        allocate(rmat(BOX,BOX,BOX), source=0.)
        ctr = real(BOX)/2. + 1.
        do k = 1, BOX
            do j = 1, BOX
                do i = 1, BOX
                    do b = 1, NBLOBS
                        dx = real(i) - ctr - CTRS(1,b); dy = real(j) - ctr - CTRS(2,b); dz = real(k) - ctr - CTRS(3,b)
                        rmat(i,j,k) = rmat(i,j,k) + AMPS(b) * exp(-(dx*dx + dy*dy + dz*dz) / (2. * SIGMAS(b)**2))
                    enddo
                enddo
            enddo
        enddo
        call vol%new([BOX,BOX,BOX], SMPD)
        call vol%set_rmat(rmat, .false.)
        call vol%write(string(PHANTOM_FILE), del_if_exists=.true.)
        call vol%kill()
    end subroutine write_phantom

    !> reference at e3 = 0 and one particle at (truth_angle, shift) from the phantom
    subroutine build_fixture( b, p, objfun, lp, truth_angle, shift, hard_edge )
        type(builder),    target, intent(inout) :: b
        type(parameters), target, intent(inout) :: p
        integer,                  intent(in)    :: objfun
        real,                     intent(in)    :: lp, truth_angle, shift(2)
        logical,                  intent(in)    :: hard_edge
        type(cmdline) :: cline
        type(ori)     :: o_ref, o_particle
        call cline%set('vol1',    PHANTOM_FILE)
        call cline%set('mskdiam', MSKDIAM)
        call cline%set('smpd',    SMPD)
        call cline%set('lp',      lp)
        call cline%set('trs',     TRS)
        call cline%set('nptcls',  1.)
        ! the smallest even projection-direction count (build_refspiral needs an even
        ! nspace); vol_pad2ref_pfts fills nspace references, the fixture uses the first
        call cline%set('nspace',  2)
        call cline%set('ctf',     'no')
        select case(objfun)
            case(OBJ_EUCLID)
                call cline%set('objfun', 'euclid')
            case(OBJ_CC)
                call cline%set('objfun', 'cc')
        end select
        call cline%check
        call b%init_params_and_build_strategy3D_tbox(cline, p)
        call set_bp_range3D(p, b, cline)
        call b%pftc%new(p, p%nspace, [1,1], p%kfromto)
        call b%vol%read(p%vols(1))
        if( hard_edge )then
            call b%vol%mask3D_hard(p%msk)
        else
            call b%vol%mask3D_soft(p%msk)
        endif
        ! the production projector: the masked volume padded by OSMPL_PAD_FAC, since the
        ! polar coordinates of the pftc live on the padded lattice (an unpadded projector
        ! is read beyond its expanded bounds at the full band)
        call b%vol_pad%new([p%box_croppd, p%box_croppd, p%box_croppd], p%smpd_crop, wthreads=.false.)
        call b%vol%pad_fft(b%vol_pad)
        call b%vol_pad%expand_cmat()
        call b%eulspace%get_ori(1, o_ref)
        call o_ref%e3set(0.0)
        o_particle = o_ref
        call o_particle%e3set(truth_angle)
        call b%eulspace%set_ori(1, o_particle)
        call vol_pad2ref_pfts(b%pftc, b%vol_pad, b%eulspace, 1, iseven=.true.)
        call b%pftc%cp_even_ref2ptcl(1, 1)
        call b%eulspace%set_ori(1, o_ref)
        call vol_pad2ref_pfts(b%pftc, b%vol_pad, b%eulspace, 1, iseven=.true.)
        call b%pftc%set_eo(1, .true.)
        ! the particle is rotated first, then shifted by the production phase: the
        ! recovered shift is expressed in the rotated frame, R(truth_angle) * shift
        if( sum(abs(shift)) > 0. ) call b%pftc%shift_ptcl(1, shift)
        call b%pftc%memoize_refs
        call b%pftc%memoize_ptcls
        if( allocated(sigma2_fixture) ) deallocate(sigma2_fixture)
        allocate(sigma2_fixture(p%kfromto(1):p%kfromto(2),2), source=1.0_sp)
        call b%pftc%assign_sigma2_noise(sigma2_fixture)
        call b%pftc%memoize_sqsum_ptcl(1)
        call cline%kill
    end subroutine build_fixture

    subroutine kill_fixture( b )
        type(builder), intent(inout) :: b
        call b%vol_pad%kill_expanded
        call b%vol_pad%kill
        call b%kill_strategy3D_tbox
        call b%kill_general_tbox
        if( allocated(sigma2_fixture) ) deallocate(sigma2_fixture)
    end subroutine kill_fixture

    ! ---- measures --------------------------------------------------------------

    !> f and grad of one evaluator at one pose
    subroutine evaluate( b, objfun, shift, rotind, f, grad )
        type(builder), intent(inout) :: b
        integer,       intent(in)    :: objfun
        real(dp),      intent(in)    :: shift(2), rotind
        real(dp),      intent(out)   :: f, grad(3)
        select case(objfun)
            case(OBJ_EUCLID)
                call b%pftc%gen_raw_euclid_grad_at_angle(1, 1, shift, rotind, f, grad)
            case(OBJ_CC)
                call b%pftc%gen_corr_grad_at_angle(1, 1, shift, rotind, f, grad)
            case default
                f    = ieee_value(f, ieee_quiet_nan)
                grad = f
        end select
    end subroutine evaluate

    !> central-difference check of the analytic gradient at one pose
    subroutine fd_check( b, objfun, shift, rotind, max_error, tol, ok )
        type(builder), intent(inout) :: b
        integer,       intent(in)    :: objfun
        real(dp),      intent(in)    :: shift(2), rotind
        real(dp),      intent(out)   :: max_error, tol
        logical,       intent(out)   :: ok
        real(dp) :: f0, g0(3), gtmp(3), fp, fm, sh(2), ri
        integer  :: idim
        call evaluate(b, objfun, shift, rotind, f0, g0)
        ok        = ieee_is_finite(f0) .and. all(ieee_is_finite(g0))
        max_error = 0.d0
        do idim = 1, 3
            sh = shift; ri = rotind
            if( idim <= 2 )then
                sh(idim) = sh(idim) + FD_STEP
            else
                ri = ri + FD_STEP
            endif
            call evaluate(b, objfun, sh, ri, fp, gtmp)
            sh = shift; ri = rotind
            if( idim <= 2 )then
                sh(idim) = sh(idim) - FD_STEP
            else
                ri = ri - FD_STEP
            endif
            call evaluate(b, objfun, sh, ri, fm, gtmp)
            max_error = max(max_error, abs((fp - fm)/(2.d0*FD_STEP) - g0(idim)))
        enddo
        tol = max(GRAD_FD_RTOL*maxval(abs(g0)), GRAD_FD_AFLOOR)
        ok  = ok .and. ieee_is_finite(max_error) .and. max_error <= tol
    end subroutine fd_check

    !> the gradient check at the twelve search-realistic probe poses around the grid selection
    subroutine gradient_probe_set( b, objfun, igrid, label )
        type(builder),    intent(inout) :: b
        integer,          intent(in)    :: objfun, igrid
        character(len=*), intent(in)    :: label
        real(dp) :: err, tol, worst
        logical  :: ok, all_ok
        integer  :: ishift, ioff
        all_ok = .true.
        worst  = 0.d0
        do ishift = 1, size(PROBE_SHIFTS, 2)
            do ioff = 1, size(PROBE_ROTOFFS)
                call fd_check(b, objfun, PROBE_SHIFTS(:,ishift), real(igrid,dp) + PROBE_ROTOFFS(ioff), err, tol, ok)
                all_ok = all_ok .and. ok
                worst  = max(worst, err)
            enddo
        enddo
        call assert_true(all_ok, label//': analytic gradient matches central differences at 12 probe poses')
        write(*,'(A,ES10.3)') '  '//label//' worst FD gradient error: ', worst
    end subroutine gradient_probe_set

    !> at integer angles and several shifts -f reproduces the discrete scoring route
    subroutine grid_identity( b, objfun, nrots, label )
        type(builder),    intent(inout) :: b
        integer,          intent(in)    :: objfun, nrots
        character(len=*), intent(in)    :: label
        real(sp), allocatable :: scores(:)
        real(dp) :: f, grad(3), worst
        integer  :: ishift, irot
        logical  :: finite
        allocate(scores(nrots))
        worst  = 0.d0
        finite = .true.
        do ishift = 1, size(PROBE_SHIFTS, 2)
            call b%pftc%gen_objfun_vals(1, 1, real(PROBE_SHIFTS(:,ishift),sp), scores)
            do irot = 1, nrots
                call evaluate(b, objfun, PROBE_SHIFTS(:,ishift), real(irot,dp), f, grad)
                worst  = max(worst, abs(-f - real(scores(irot),dp)))
                finite = finite .and. ieee_is_finite(f) .and. all(ieee_is_finite(grad))
            enddo
        enddo
        call assert_true(finite, label//': evaluator is finite at every grid angle and probe shift')
        call assert_true(worst <= IDENTITY_ATOL, label//': -f equals the gen_objfun_vals score at grid angles (1e-4)')
    end subroutine grid_identity

    real(dp) function angle_of( pftc, rotind )
        class(polarft_calc), intent(in) :: pftc
        real(dp),            intent(in) :: rotind
        angle_of = (rotind - 1.d0) * real(pftc%get_dang(),dp)
    end function angle_of

    pure real(dp) function angular_error( angle, target )
        real(dp), intent(in) :: angle, target
        angular_error = abs(modulo(angle - target + 180.d0, 360.d0) - 180.d0)
    end function angular_error

    ! ---- tests -----------------------------------------------------------------

    !> raw Euclidean: route identities, gradient, parity, stress, seed parity, joint recovery
    subroutine test_euclid_routes_and_recovery()
        type(builder),    target :: b
        type(parameters), target :: p
        type(pftc_shsrch_grad) :: fixed_search, joint_search
        real(sp), allocatable :: scores(:), raw_losses(:)
        real(dp), parameter   :: PARITY_SHIFTS(2,3) = reshape([1.7d0,-2.3d0, -0.6d0,0.4d0, 3.1d0,1.2d0], [2,3])
        real(dp) :: f_cont, g_cont(3), f_disc, g_disc(2), scalar_loss, grad2(2)
        real(dp) :: max_legacy, max_scalar, parity_err, parity_scale, tol
        real(dp) :: series_min, stale_parity, fresh_parity
        real(dp) :: grid_loss, joint_loss, initial_loss, joint_coord, recovered(2), expected(2), theta, g3(3)
        real(dp) :: grid_angle_err, joint_angle_err, shift_rms
        real     :: cxy(3), seed(2), corr, limits(2,2), joint_limits(3,2)
        integer  :: nrots, irot, igrid, irots(3), j, m, k, selected, expected_irot, irot_fixed, joint_index
        logical  :: finite, valid, improved
        write(*,'(A)') 'test_euclid_routes_and_recovery'
        call build_fixture(b, p, OBJ_EUCLID, LP_LOW, TRUTH_ANGLE, TRUTH_SHIFT, .false.)
        nrots = b%pftc%get_nrots()
        call assert_true(nrots > 8, 'polar grid has rotations')
        allocate(scores(nrots), raw_losses(nrots))
        ! (a) the vector scoring routes share one loss: score = exp(-raw loss), and the
        !     scalar discrete route returns the same raw loss
        call b%pftc%gen_objfun_vals(1, 1, [0.0_sp, 0.0_sp], scores)
        call b%pftc%gen_raw_euclid_vals(1, 1, [0.0_sp, 0.0_sp], raw_losses)
        max_legacy = 0.d0
        max_scalar = 0.d0
        do irot = 1, nrots
            max_legacy = max(max_legacy, abs(real(scores(irot),dp) - exp(-real(raw_losses(irot),dp))))
            call b%pftc%gen_raw_euclid_grad_for_rot_8(1, 1, [0.0_dp, 0.0_dp], irot, scalar_loss, grad2)
            max_scalar = max(max_scalar, abs(scalar_loss - real(raw_losses(irot),dp)))
        enddo
        tol = 5.0d-4 * (1.d0 + maxval(abs(real(raw_losses,dp))))
        call assert_true(max_legacy <= tol, 'gen_objfun_vals = exp(-gen_raw_euclid_vals) at every angle')
        call assert_true(max_scalar <= tol, 'gen_raw_euclid_grad_for_rot_8 returns the vector route loss')
        ! (b) the continuous evaluator at integer angles equals the discrete reference
        igrid = minloc(raw_losses, dim=1)
        irots = [igrid, modulo(igrid-1+nrots/3, nrots)+1, modulo(igrid-1+(2*nrots)/3, nrots)+1]
        parity_err   = 0.d0
        parity_scale = 0.d0
        finite       = .true.
        do j = 1, 3
            do m = 1, 3
                call b%pftc%gen_raw_euclid_grad_at_angle(1, 1, PARITY_SHIFTS(:,m), real(irots(j),dp), f_cont, g_cont)
                call b%pftc%gen_raw_euclid_grad_for_rot_8(1, 1, PARITY_SHIFTS(:,m), irots(j), f_disc, g_disc)
                finite = finite .and. ieee_is_finite(f_cont) .and. all(ieee_is_finite(g_cont)) .and. &
                    &ieee_is_finite(f_disc) .and. all(ieee_is_finite(g_disc))
                parity_err   = max(parity_err, abs(f_cont-f_disc), abs(g_cont(1)-g_disc(1)), abs(g_cont(2)-g_disc(2)))
                parity_scale = max(parity_scale, abs(f_disc), maxval(abs(g_disc)))
            enddo
        enddo
        call assert_true(finite, 'continuous and discrete evaluators are finite at grid angles')
        call assert_true(parity_err <= GRAD_FD_RTOL*max(1.d0, parity_scale), &
            &'continuous evaluator equals the discrete reference (loss and x/y gradient) at grid angles')
        ! (c) analytic gradient against central differences
        call gradient_probe_set(b, OBJ_EUCLID, igrid, 'euclid')
        ! (d) seed parity: the joint route's discrete selection is one alternating-search callback at the seed
        seed = [0.5, -0.25]
        call b%pftc%gen_objfun_vals(1, 1, seed, scores)
        expected_irot = maxloc(scores, dim=1)
        joint_limits(1:2,1) = -1.; joint_limits(1:2,2) = 1.
        joint_limits(3,:)   = [1.-2., real(nrots)+2.]
        call joint_search%new_joint(b, joint_limits, 100)
        call joint_search%set_indices(1, 1)
        call joint_search%select_best_discrete_angle(seed, selected, corr)
        call assert_false(joint_search%does_opt_angle(), 'joint seed selector does not attach the alternating-search callback')
        call assert_int(expected_irot, selected, 'joint seed selector reproduces the callback index at the x/y seed')
        call assert_true(abs(corr - scores(expected_irot)) <= 5.e-6, 'joint seed selector reproduces the callback score')
        call joint_search%kill
        ! (e) recovery: the fixed (grid-angle) shift search, then the joint solve from the grid seed
        limits(:,1) = -TRS; limits(:,2) = TRS
        call fixed_search%new_fixed(b, limits, maxits=16)
        call fixed_search%set_indices(1, 1)
        irot_fixed = igrid
        cxy = fixed_search%minimize(irot_fixed, sh_rot=.false., xy_in=[0.,0.])
        call assert_true(irot_fixed > 0, 'fixed-angle shift search accepts a solution')
        call fixed_search%kill
        call b%pftc%gen_raw_euclid_grad_at_angle(1, 1, [0.d0,0.d0], real(igrid,dp), grid_loss, g3)
        joint_limits(:,1) = [-TRS, -TRS, real(igrid)-2.]
        joint_limits(:,2) = [ TRS,  TRS, real(igrid)+2.]
        call joint_search%new_joint(b, joint_limits, 100)
        call joint_search%set_indices(1, 1)
        joint_index = igrid
        cxy = joint_search%minimize_joint(joint_index, [0.,0.], sh_rot=.false., rotind_frac=joint_coord, &
            &evaluation_valid=valid, improved=improved, initial_cost_out=initial_loss, irot_in=igrid)
        call assert_true(valid, 'joint evaluation from the grid seed is valid')
        call assert_true(improved .and. joint_index > 0, 'joint solve improves on the grid seed')
        recovered = real(cxy(2:3),dp)
        call b%pftc%gen_raw_euclid_grad_at_angle(1, 1, recovered, joint_coord, joint_loss, g3)
        call joint_search%kill
        theta    = real(TRUTH_ANGLE,dp) * acos(-1.d0) / 180.d0
        expected = [cos(theta)*real(TRUTH_SHIFT(1),dp) - sin(theta)*real(TRUTH_SHIFT(2),dp), &
                   &sin(theta)*real(TRUTH_SHIFT(1),dp) + cos(theta)*real(TRUTH_SHIFT(2),dp)]
        grid_angle_err  = angular_error(angle_of(b%pftc, real(igrid,dp)), -real(TRUTH_ANGLE,dp))
        joint_angle_err = angular_error(angle_of(b%pftc, joint_coord),    -real(TRUTH_ANGLE,dp))
        shift_rms       = sqrt(sum((recovered - expected)**2)/2.d0)
        call assert_true(all(ieee_is_finite([grid_loss, joint_loss, initial_loss, joint_angle_err, shift_rms])), &
            &'recovery quantities are finite')
        call assert_true(abs(initial_loss - grid_loss) <= 16.d0*1.d-8*max(1.d0, abs(grid_loss), abs(joint_loss)), &
            &'the joint seed cost is the grid objective')
        call assert_true(joint_loss <= grid_loss + 1.d-8*max(1.d0, abs(grid_loss), abs(joint_loss)), &
            &'the joint solve does not worsen the objective')
        call assert_true(joint_angle_err < grid_angle_err - 1.d-4, 'the joint solve improves the angle over the grid')
        call assert_true(shift_rms <= 0.25d0, 'the joint solve recovers the known shift within 0.25 px rms')
        write(*,'(A,3F9.4)') '  grid/joint angle error (deg), shift rms (px): ', grid_angle_err, joint_angle_err, shift_rms
        ! (f) stress: cavgs-like near-noiseless shell-dependent sigma2; assigning the spectra refreshes the
        !     shell weights and the memoised self-terms (memoize_sigma2_weights), then re-memoising changes nothing
        do k = p%kfromto(1), p%kfromto(2)
            sigma2_fixture(k,:) = 1.e-4_sp / real(k,sp)
        enddo
        call b%pftc%assign_sigma2_noise(sigma2_fixture)
        call series_floor_scan(b, igrid, series_min, stale_parity)
        call b%pftc%gen_raw_euclid_vals(1, 1, [1.7_sp, -2.3_sp], raw_losses)
        call b%pftc%gen_objfun_vals(1, 1, [1.7_sp, -2.3_sp], scores)
        series_min = min(series_min, real(minval(raw_losses),dp))
        call assert_true(ieee_is_finite(series_min) .and. series_min > -GRAD_FD_RTOL, &
            &'loss series stays non-negative under near-noiseless sigma2 with stale square sums')
        call assert_true(all(ieee_is_finite(scores)) .and. minval(scores) >= 0._sp .and. maxval(scores) <= 1._sp, &
            &'scores stay in [0,1] under near-noiseless sigma2')
        call b%pftc%memoize_sqsum_ptcl(1)
        call series_floor_scan(b, igrid, series_min, fresh_parity)
        call assert_true(series_min > -GRAD_FD_RTOL, 'loss series stays non-negative after re-memoising')
        call assert_true(fresh_parity <= GRAD_FD_RTOL, 'continuous and discrete losses agree at grid angles under stress')
        call kill_fixture(b)
    end subroutine test_euclid_routes_and_recovery

    !> dense scan of the joint objective over the searched (shift, fractional angle) domain
    subroutine series_floor_scan( b, igrid, series_min, parity_max )
        type(builder), intent(inout) :: b
        integer,       intent(in)    :: igrid
        real(dp),      intent(out)   :: series_min, parity_max
        real(dp) :: f, g(3), f_disc, g_disc(2), shift(2)
        integer  :: ish, jsh, ith
        series_min = huge(0.d0)
        parity_max = 0.d0
        do ish = -1, 1
            do jsh = -1, 1
                shift = [real(ish,dp)*1.25d0, real(jsh,dp)*1.25d0]
                do ith = -20, 20
                    call b%pftc%gen_raw_euclid_grad_at_angle(1, 1, shift, real(igrid,dp) + real(ith,dp)*0.1d0, f, g)
                    series_min = min(series_min, f)
                enddo
                call b%pftc%gen_raw_euclid_grad_at_angle(1, 1, shift, real(igrid,dp), f, g)
                call b%pftc%gen_raw_euclid_grad_for_rot_8(1, 1, shift, igrid, f_disc, g_disc)
                parity_max = max(parity_max, abs(f - f_disc))
            enddo
        enddo
    end subroutine series_floor_scan

    !> the wrap-around angle with a hard mask edge at the full band, and the zero angle
    subroutine test_euclid_band_edges()
        type(builder),    target :: b
        type(parameters), target :: p
        real(sp), allocatable :: raw_losses(:)
        real(dp) :: f, g(3)
        integer  :: nrots, igrid
        write(*,'(A)') 'test_euclid_band_edges'
        ! full band, hard edge, truth just below 360
        call build_fixture(b, p, OBJ_EUCLID, LP_FULL, 359.375, [0.,0.], .true.)
        nrots = b%pftc%get_nrots()
        allocate(raw_losses(nrots))
        call b%pftc%gen_raw_euclid_vals(1, 1, [0._sp,0._sp], raw_losses)
        igrid = minloc(raw_losses, dim=1)
        call assert_true(all(ieee_is_finite(raw_losses)), 'full-band losses are finite')
        call assert_true(angular_error(angle_of(b%pftc, real(igrid,dp)), -359.375d0) <= 1.5d0*real(b%pftc%get_dang(),dp), &
            &'the grid selection is within one step of the wrapped truth angle')
        call gradient_probe_set(b, OBJ_EUCLID, igrid, 'euclid full band')
        call b%pftc%gen_raw_euclid_grad_at_angle(1, 1, [0.d0,0.d0], real(igrid,dp) - 0.5d0, f, g)
        call assert_true(ieee_is_finite(f) .and. all(ieee_is_finite(g)), 'fractional angle across the wrap is finite')
        call kill_fixture(b)
        deallocate(raw_losses)
        ! zero truth angle, low band
        call build_fixture(b, p, OBJ_EUCLID, LP_LOW, 0., [0.,0.], .false.)
        nrots = b%pftc%get_nrots()
        allocate(raw_losses(nrots))
        call b%pftc%gen_raw_euclid_vals(1, 1, [0._sp,0._sp], raw_losses)
        igrid = minloc(raw_losses, dim=1)
        call assert_int(1, igrid, 'an unrotated, unshifted particle selects the first grid angle')
        call assert_true(raw_losses(1) <= 1.d-3*max(1.d0, real(maxval(raw_losses),dp)), 'the loss vanishes at the truth pose')
        call gradient_probe_set(b, OBJ_EUCLID, igrid, 'euclid zero angle')
        call kill_fixture(b)
    end subroutine test_euclid_band_edges

    !> cc: grid identity, gradient, and the finite penalty for an undefined denominator
    subroutine test_cc_evaluator()
        type(builder),    target :: b
        type(parameters), target :: p
        real(sp),    allocatable :: scores(:)
        complex(sp), allocatable :: zero_pft(:,:)
        real(dp) :: f, grad(3)
        integer  :: nrots, igrid
        write(*,'(A)') 'test_cc_evaluator'
        call build_fixture(b, p, OBJ_CC, LP_LOW, TRUTH_ANGLE, TRUTH_SHIFT, .false.)
        nrots = b%pftc%get_nrots()
        allocate(scores(nrots))
        call grid_identity(b, OBJ_CC, nrots, 'cc')
        call b%pftc%gen_objfun_vals(1, 1, [0._sp,0._sp], scores)
        igrid = maxloc(scores, dim=1)
        call gradient_probe_set(b, OBJ_CC, igrid, 'cc')
        ! a zero reference gives an identically zero denominator series at every angle:
        ! the evaluator must return a finite loss above the physical range (> 1) with a
        ! zero gradient, so such a pose can never beat a seed (run last: destroys the reference)
        zero_pft = b%pftc%allocate_pft()
        zero_pft = cmplx(0._sp, 0._sp, kind=sp)
        call b%pftc%set_ref_pft(1, zero_pft, iseven=.true.)
        call b%pftc%memoize_refs
        call b%pftc%gen_corr_grad_at_angle(1, 1, [0.d0,0.d0], real(igrid,dp) + 0.3d0, f, grad)
        call assert_true(ieee_is_finite(f) .and. all(ieee_is_finite(grad)), 'degenerate denominator gives finite output')
        call assert_true(f > 1.d0, 'degenerate denominator gives the penalty loss (> 1)')
        call assert_true(all(grad == 0.d0), 'degenerate denominator gives a zero gradient')
        call kill_fixture(b)
    end subroutine test_cc_evaluator

    !> strategy2D_srch constructs the joint optimiser only under inpl_cont=yes and always
    !! keeps the alternating seed-search angle update (selection parity)
    subroutine test_strategy2D_route_flags()
        type(builder),    target :: b
        type(parameters), target :: p
        type(strategy2D_srch) :: srch
        type(strategy2D_spec) :: spec
        type(pftc_shsrch_grad) :: alternating
        real :: limits(2,2)
        write(*,'(A)') 'test_strategy2D_route_flags'
        call build_fixture(b, p, OBJ_EUCLID, LP_LOW, TRUTH_ANGLE, TRUTH_SHIFT, .false.)
        limits(:,1) = -1.; limits(:,2) = 1.
        call alternating%new_alternating(b, limits)
        call alternating%kill
        spec%iptcl       = 1
        spec%iptcl_batch = 1
        spec%iptcl_map   = 1
        p%l_prob_align_mode = .false.
        p%inpl_cont         = 'no'
        call srch%new(p, spec, b)
        call assert_false(srch%uses_continuous_refinement(), 'inpl_cont=no: no continuous polish')
        call assert_false(srch%joint_inpl_optimizer%uses_joint_inplane(), 'inpl_cont=no: no joint optimiser')
        call assert_true(srch%grad_shsrch_first_obj%does_opt_angle(), 'inpl_cont=no: alternating seed search updates the angle')
        call srch%kill
        p%inpl_cont = 'yes'
        call srch%new(p, spec, b)
        call assert_true(srch%uses_continuous_refinement(), 'inpl_cont=yes: continuous polish')
        call assert_true(srch%joint_inpl_optimizer%uses_joint_inplane(), 'inpl_cont=yes: joint optimiser constructed')
        call assert_true(srch%grad_shsrch_first_obj%does_opt_angle(), &
            &'inpl_cont=yes: alternating seed search still updates the angle (selection parity)')
        call srch%kill
        p%l_prob_align_mode = .true.
        call srch%new(p, spec, b)
        call assert_true(srch%uses_continuous_refinement(), 'probabilistic route: continuous polish of the selected candidate')
        call assert_true(srch%joint_inpl_optimizer%uses_joint_inplane(), 'probabilistic route: joint optimiser constructed')
        call srch%kill
        p%l_prob_align_mode = .false.
        p%inpl_cont         = 'no'
        call kill_fixture(b)
    end subroutine test_strategy2D_route_flags

    ! ---- quadratic shift step from the coarse grid (multi-state table) -------

    !> the pure step: exact recovery of a quadratic's stationary point, and each guard that refuses a fit
    subroutine test_grid_fit_step()
        real(dp), parameter :: HESS(2,2) = reshape([2.d0, 0.5d0, 0.5d0, 1.5d0], [2,2])
        real(dp) :: v(-1:1,-1:1), da, db
        integer  :: ia, ib, reason
        write(*,'(A)') 'test_grid_fit_step'
        ! an exact quadratic with its stationary point inside the centre cell
        call fill_quadratic([0.3d0, -0.2d0], v)
        call fit_quadratic_step(v, da, db, reason)
        call assert_int(0, reason, 'an exact quadratic is accepted')
        call assert_true(abs(da - 0.3d0) <= 1.d-10 .and. abs(db + 0.2d0) <= 1.d-10, &
            &'the stationary point of an exact quadratic is recovered (1e-10 grid units)')
        ! a stationary point outside the centre cell
        call fill_quadratic([0.8d0, 0.d0], v)
        call fit_quadratic_step(v, da, db, reason)
        call assert_int(FIT_STEP, reason, 'a stationary point outside the centre cell is refused')
        ! not positive definite: a maximum, a saddle, a plane
        do ia = -1,1
            do ib = -1,1
                v(ia,ib) = -real(ia*ia + ib*ib, dp)
            enddo
        enddo
        call fit_quadratic_step(v, da, db, reason)
        call assert_int(FIT_CURVATURE, reason, 'a maximum is refused')
        do ia = -1,1
            do ib = -1,1
                v(ia,ib) = real(ia*ia - ib*ib, dp)
            enddo
        enddo
        call fit_quadratic_step(v, da, db, reason)
        call assert_int(FIT_CURVATURE, reason, 'a saddle is refused')
        do ia = -1,1
            do ib = -1,1
                v(ia,ib) = 0.3d0 * real(ia, dp)
            enddo
        enddo
        call fit_quadratic_step(v, da, db, reason)
        call assert_int(FIT_CURVATURE, reason, 'a plane is refused')
        ! ill conditioned: eigenvalue ratio 1e-3, below the 1e-2 floor
        do ia = -1,1
            do ib = -1,1
                v(ia,ib) = real(ia*ia, dp) + 1.d-3 * real(ib*ib, dp)
            enddo
        enddo
        call fit_quadratic_step(v, da, db, reason)
        call assert_int(FIT_CURVATURE, reason, 'an ill-conditioned curvature is refused')
        ! the stencil exists away from the edge of the 5 x 5 lattice only
        call assert_true(grid_stencil_exists(3,3) .and. grid_stencil_exists(2,2) .and. grid_stencil_exists(4,4), &
            &'interior grid points have a full stencil')
        call assert_false(grid_stencil_exists(1,3) .or. grid_stencil_exists(5,3) .or. grid_stencil_exists(3,1) .or. &
            &grid_stencil_exists(3,5), 'edge grid points have none')

    contains

        subroutine fill_quadratic( x0, vv )
            real(dp), intent(in)  :: x0(2)
            real(dp), intent(out) :: vv(-1:1,-1:1)
            real(dp) :: d(2)
            integer  :: i, j
            do i = -1,1
                do j = -1,1
                    d = [real(i,dp), real(j,dp)] - x0
                    vv(i,j) = 7.d0 + 0.5d0 * dot_product(d, matmul(HESS, d))
                enddo
            enddo
        end subroutine fill_quadratic

    end subroutine test_grid_fit_step

    !> grid_fit against the L-BFGS-B path from the same grid on the phantom: a kept fit is within 0.2 px of its
    !! shift, a fallback and the option off are bit-identical to it. The printed counts say why a case did not fit.
    subroutine test_grid_fit_search()
        integer :: icase
        write(*,'(A)') 'test_grid_fit_search'
        do icase = 1, 3
            call grid_fit_case(icase)
        enddo

    contains

        !> one phantom fixture per call, so the builder and parameters are fresh every time
        subroutine grid_fit_case( icase_loc )
            integer, intent(in) :: icase_loc
            real(sp), parameter :: TRUE_SHIFTS(2,3) = reshape([0.5, 0.3, 1.2, -0.8, 4.6, 0.0], [2,3])
            type(builder),    target :: b
            type(parameters), target :: p
            type(pftc_shsrch_grad) :: base, with_fit, off
            integer(int64) :: counts(0:5), counts_off(0:5)
            real    :: limits(2,2), lims_init(2,2), cxy_a(3), cxy_b(3), cxy_c(3)
            integer :: irot_a, irot_b, irot_c
            limits(:,1)    = -TRS
            limits(:,2)    =  TRS
            lims_init(:,1) = -SHC_INPL_TRSHWDTH
            lims_init(:,2) =  SHC_INPL_TRSHWDTH
            call build_fixture(b, p, OBJ_EUCLID, LP_LOW, TRUTH_ANGLE, TRUE_SHIFTS(:,icase_loc), .false.)
            call base%new_alternating(b, limits, lims_init=lims_init, maxits=60, coarse_init=.true.)
            call with_fit%new_alternating(b, limits, lims_init=lims_init, maxits=60, coarse_init=.true., grid_fit=.true.)
            call off%new_alternating(b, limits, lims_init=lims_init, maxits=60, coarse_init=.true., grid_fit=.false.)
            call base%set_indices(1, 1)
            call with_fit%set_indices(1, 1)
            call off%set_indices(1, 1)
            irot_a = 1; irot_b = 1; irot_c = 1
            cxy_a = base%minimize(irot_a, sh_rot=.false.)
            cxy_b = with_fit%minimize(irot_b, sh_rot=.false.)
            cxy_c = off%minimize(irot_c, sh_rot=.false.)
            call with_fit%get_fit_profile(counts)
            call off%get_fit_profile(counts_off)
            call assert_int(1, int(sum(counts)), 'each shift search takes exactly one fit decision')
            call assert_true(all(counts_off == 0_int64), 'no fit decisions without the option')
            call assert_true(same_result(irot_c, cxy_c, irot_a, cxy_a), &
                &'option off: identical to the search built without it')
            if( counts(0) == 1_int64 )then
                ! irot = 0 flags a search that found nothing better than zero shift; cxy is then undefined
                if( irot_a > 0 .and. irot_b > 0 )then
                    call assert_true(all(abs(cxy_b(2:3) - cxy_a(2:3)) <= 0.2), &
                        &'a kept fit lies within 0.2 px of the L-BFGS-B shift from the same grid point')
                endif
            else
                call assert_true(same_result(irot_b, cxy_b, irot_a, cxy_a), 'a fallback is identical to the L-BFGS-B path')
            endif
            write(*,'(A,I0,A,6I4)') '  case ', icase_loc, ' fit counts (kept, edge, curvature, step, box, value): ', int(counts)
            write(*,'(A,2F8.3,A,2F8.3)') '  L-BFGS-B shift: ', cxy_a(2:3), '   with fit: ', cxy_b(2:3)
            call base%kill
            call with_fit%kill
            call off%kill
            call kill_fixture(b)
        end subroutine grid_fit_case

        logical function same_result( irot1, cxy1, irot2, cxy2 )
            integer, intent(in) :: irot1, irot2
            real,    intent(in) :: cxy1(3), cxy2(3)
            same_result = irot1 == irot2
            if( same_result .and. irot1 > 0 ) same_result = all(cxy1 == cxy2)
        end function same_result

    end subroutine test_grid_fit_search

    ! ---- particle-frame shifts (the probabilistic tables) ----------------------

    !> the two frame conversions: a quarter turn, the round trip, and a lattice corner that leaves the box
    subroutine test_shift_frames()
        real :: u(2), s(2), c(2)
        integer :: irot
        write(*,'(A)') 'test_shift_frames'
        u = shift_ref2ptcl([1., 0.], 90.)
        call assert_true(all(abs(u - [0., 1.]) <= 1.e-5), 'a reference-frame shift at 90 degrees is a quarter turn of the particle shift')
        s = shift_ptcl2ref([0., 1.], 90.)
        call assert_true(all(abs(s - [1., 0.]) <= 1.e-5), 'the particle-frame conversion is its inverse')
        s = [1.7, -2.3]
        do irot = 0, 7
            u = shift_ref2ptcl(s, 47.5 * real(irot))
            c = shift_ptcl2ref(u, 47.5 * real(irot))
            call assert_true(all(abs(c - s) <= 1.e-5), 'the round trip is the identity')
            call assert_true(abs(sum(u**2) - sum(s**2)) <= 1.e-4, 'a rotation keeps the length')
        enddo
        ! the corner of the 5 x 5 lattice, 0.8 trs in both axes, is 1.13 trs from the origin
        c = shift_ptcl2ref([0.8 * TRS, 0.8 * TRS], 45.)
        call assert_true(c(1) > TRS, 'a lattice corner can leave the search box: the start is clamped')
    end subroutine test_shift_frames

    !> a bank entry (shifted particle memo) against the direct sum at the reference-frame shift it corresponds to
    !! at every rotation, for both objectives, and the likelihood route on the same entry
    subroutine test_particle_bank()
        real(dp), parameter :: TODAY_ATOL = 5.d-2 ! the dropped Nyquist bin separates the two kernels; a frame error is O(1)
        type(builder),    target :: b
        type(parameters), target :: p
        complex(sp), allocatable :: bank(:,:)
        real(sp),    allocatable :: vals(:), vals_today(:), pvec(:)
        integer,     allocatable :: inds(:)
        real(dp) :: max_err, max_today, direct
        real(sp) :: t(2), s(2), dist, corr
        integer  :: iobj, ish, irot, irot_s, nrots, kfromto(2), irot_seed
        write(*,'(A)') 'test_particle_bank'
        do iobj = OBJ_EUCLID, OBJ_CC
            call build_fixture(b, p, iobj, LP_LOW, TRUTH_ANGLE, [2., -1.5], .false.)
            nrots   = b%pftc%get_nrots()
            kfromto = b%pftc%get_kfromto()
            allocate(bank(b%pftc%get_pftsz()+1, kfromto(2)-kfromto(1)+1), vals(nrots), vals_today(nrots), pvec(nrots), inds(nrots))
            max_err   = 0.d0
            max_today = 0.d0
            do ish = 1, size(PROBE_SHIFTS,2)
                t = real(PROBE_SHIFTS(:,ish), sp)
                call b%pftc%memoize_shifted_ptcl(1, t, bank)
                call b%pftc%gen_objfun_vals_bank(1, 1, bank, vals)
                do irot = 1, nrots
                    s = shift_ptcl2ref(t, b%pftc%get_rot(irot))
                    call b%pftc%gen_objfun_vals(1, 1, s, vals_today)
                    direct = direct_reference(iobj, s, irot, vals_today)
                    max_err = max(max_err, abs(real(vals(irot),dp) - direct))
                    max_today = max(max_today, abs(real(vals(irot) - vals_today(irot),dp)))
                enddo
            enddo
            call assert_true(max_err <= IDENTITY_ATOL, &
                &'a bank entry scores every rotation as the reference shifted by its particle shift (direct sum)')
            call assert_true(max_today <= TODAY_ATOL, 'and agrees with the shifted-reference kernel up to the Nyquist bin')
            write(*,'(A,I0,A,ES10.2,A,ES10.2)') '  objective ', iobj, ': max |bank - direct| ', max_err, &
                &'   max |bank - shifted-reference kernel| ', max_today
            ! the likelihood route with one sample takes the best rotation of the same entry
            call b%pftc%gen_likelihood_val_bank(1, 1, bank, 1, dist, corr, irot_s, pvec, inds)
            call assert_int(maxloc(vals, dim=1), irot_s, 'the likelihood route picks the best rotation of the entry')
            call assert_true(abs(corr - vals(irot_s)) <= 1.e-4, 'and reports its score')
            ! a seed from the shift search: estimated at irot_seed as a reference shift, held as one particle shift
            irot_seed = nrots / 3
            s = [1.2, -0.8]
            t = shift_ref2ptcl(s, b%pftc%get_rot(irot_seed))
            call b%pftc%memoize_shifted_ptcl(1, t, bank)
            call b%pftc%gen_objfun_vals_bank(1, 1, bank, vals)
            call b%pftc%gen_objfun_vals(1, 1, s, vals_today)
            direct = direct_reference(iobj, s, irot_seed, vals_today)
            call assert_true(abs(real(vals(irot_seed),dp) - direct) <= IDENTITY_ATOL, &
                &'a seed scores at its own rotation exactly as the reference-frame seed')
            call assert_true(abs(vals(irot_seed) - vals_today(irot_seed)) <= TODAY_ATOL, &
                &'and as today''s kernel up to the Nyquist bin')
            deallocate(bank, vals, vals_today, pvec, inds)
            call kill_fixture(b)
        enddo

    contains

        !> the reference-frame score at one rotation: the direct sum for Euclid; for cc the direct sum is the
        !! k-weighted correlation, a different objective from the FFT route, so the shifted-reference kernel is the reference
        real(dp) function direct_reference( iobj_loc, s_loc, irot_loc, today_loc )
            integer,  intent(in) :: iobj_loc, irot_loc
            real(sp), intent(in) :: s_loc(2), today_loc(:)
            if( iobj_loc == OBJ_EUCLID )then
                direct_reference = b%pftc%gen_corr_for_rot_8(1, 1, real(s_loc,dp), irot_loc)
            else
                direct_reference = real(today_loc(irot_loc), dp)
            endif
        end function direct_reference

    end subroutine test_particle_bank

    !> the particle-frame grid (multi-state table) on the phantom: option off is bit-identical, the bank is built once per particle, the result
    !! agrees with the reference-frame grid, and the quadratic step works on the lattice
    subroutine test_grid_ptcl_search()
        integer :: icase
        write(*,'(A)') 'test_grid_ptcl_search'
        ! the interior shift and the one next to the box edge
        do icase = 1, 3, 2
            call grid_ptcl_case(icase)
        enddo

    contains

        subroutine grid_ptcl_case( icase_loc )
            integer, intent(in) :: icase_loc
            real(sp), parameter :: TRUE_SHIFTS(2,3) = reshape([0.5, 0.3, 1.2, -0.8, 4.6, 0.0], [2,3])
            type(builder),    target :: b
            type(parameters), target :: p
            type(pftc_shsrch_grad) :: base, off, ptcl, ptcl_fit
            integer(int64) :: counts(0:5)
            real    :: limits(2,2), lims_init(2,2), cxy_a(3), cxy_o(3), cxy_p(3), cxy_f(3)
            integer :: irot_a, irot_o, irot_p, irot_f, nrots, drot
            limits(:,1)    = -TRS
            limits(:,2)    =  TRS
            lims_init(:,1) = -SHC_INPL_TRSHWDTH
            lims_init(:,2) =  SHC_INPL_TRSHWDTH
            call build_fixture(b, p, OBJ_EUCLID, LP_LOW, TRUTH_ANGLE, TRUE_SHIFTS(:,icase_loc), .false.)
            nrots = b%pftc%get_nrots()
            call base%new_alternating(b, limits, lims_init=lims_init, maxits=60, coarse_init=.true.)
            call off%new_alternating(b, limits, lims_init=lims_init, maxits=60, coarse_init=.true., grid_ptcl=.false.)
            call ptcl%new_alternating(b, limits, lims_init=lims_init, maxits=60, coarse_init=.true., grid_ptcl=.true.)
            call ptcl_fit%new_alternating(b, limits, lims_init=lims_init, maxits=60, coarse_init=.true., grid_ptcl=.true., &
                &grid_fit=.true.)
            call base%set_indices(1, 1)
            call off%set_indices(1, 1)
            call ptcl%set_indices(1, 1)
            call ptcl_fit%set_indices(1, 1)
            irot_a = 1; irot_o = 1; irot_p = 1; irot_f = 1
            cxy_a = base%minimize(irot_a, sh_rot=.false.)
            cxy_o = off%minimize(irot_o, sh_rot=.false.)
            cxy_p = ptcl%minimize(irot_p, sh_rot=.false.)
            cxy_f = ptcl_fit%minimize(irot_f, sh_rot=.false.)
            call assert_true(irot_o == irot_a, 'option off: same rotation as the search built without it')
            if( irot_o > 0 ) call assert_true(all(cxy_o == cxy_a), 'option off: bit-identical to the search built without it')
            call assert_true(off%get_bank_builds() == 0_int64, 'option off: no shifted particle memos')
            call assert_true(ptcl%get_bank_builds() == 25_int64, 'one memo per lattice point')
            ! a second reference of the same particle reuses the bank
            call ptcl%set_indices(2, 1)
            irot_p = 1
            cxy_p = ptcl%minimize(irot_p, sh_rot=.false.)
            call assert_true(ptcl%get_bank_builds() == 25_int64, 'the bank is built once per particle and shared by its references')
            call ptcl%set_indices(1, 1)
            irot_p = 1
            cxy_p = ptcl%minimize(irot_p, sh_rot=.false.)
            if( irot_a > 0 .and. irot_p > 0 )then
                drot = abs(irot_p - irot_a)
                call assert_true(min(drot, nrots - drot) <= 1, 'the particle-frame grid finds the same rotation')
                call assert_true(all(abs(cxy_p(2:3) - cxy_a(2:3)) <= 0.3), 'and the same shift within 0.3 px')
                call assert_true(cxy_p(1) >= cxy_a(1) - 1.e-3, 'and is not worse than the reference-frame grid')
                call assert_true(all(abs(cxy_p(2:3)) <= TRS + 1.e-3), 'and stays in the search box')
            endif
            call ptcl_fit%get_fit_profile(counts)
            call assert_int(1, int(sum(counts)), 'the quadratic step takes exactly one decision on the lattice')
            if( irot_p > 0 .and. irot_f > 0 )then
                call assert_true(all(abs(cxy_f(2:3) - cxy_p(2:3)) <= 0.3), 'with the step the shift agrees with the lattice search')
            endif
            write(*,'(A,I0,A,6I4)') '  case ', icase_loc, ' lattice fit counts (kept, edge, curvature, step, box, value): ', int(counts)
            write(*,'(A,2F8.3,A,2F8.3,A,2F8.3)') '  reference grid: ', cxy_a(2:3), '   particle grid: ', cxy_p(2:3), &
                &'   with fit: ', cxy_f(2:3)
            call base%kill
            call off%kill
            call ptcl%kill
            call ptcl_fit%kill
            call kill_fixture(b)
        end subroutine grid_ptcl_case

    end subroutine test_grid_ptcl_search

end module simple_pftc_inplane_tester
