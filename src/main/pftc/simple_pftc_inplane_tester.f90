!@descr: unit tests for continuous in-plane registration on the polar Fourier transform (simple_polarft_calc, simple_pftc_shsrch_grad)
! One synthetic phantom volume (written to the run directory for the suite, removed at
! the end) projected through the strategy3D toolbox gives a reference at e3 = 0 and a
! particle rotated and shifted by a known amount. Pinned, for the raw Euclidean, cc and
! hybrid continuous-angle evaluators: identity with the discrete scoring routes at
! integer angles, analytic against central-difference gradients in (sx, sy, theta),
! the non-negative loss series under near-noiseless sigma2 (stale and re-memoised
! square sums), the cc penalty for an undefined denominator, the joint route's seed
! parity with the legacy callback, recovery of the known angle and shift by the joint
! solve, and the strategy2D route construction under inpl_cont=no|yes.
module simple_pftc_inplane_tester
use simple_pftc_srch_api
use simple_string,                  only: string
use simple_syslib,                  only: del_file
use simple_builder,                 only: builder
use simple_image,                   only: image
use simple_matcher_smpl_and_lplims, only: set_bp_range3D
use simple_polarft_calc,            only: polarft_calc, vol_pad2ref_pfts
use simple_pftc_shsrch_grad,        only: pftc_shsrch_grad
use simple_strategy2D_srch,         only: strategy2D_srch, strategy2D_spec
use simple_test_utils
use, intrinsic :: ieee_arithmetic,  only: ieee_is_finite, ieee_value, ieee_quiet_nan
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

integer, parameter :: OBJ_EUCLID = 1, OBJ_CC = 2, OBJ_HYBRID = 3

! the pftc keeps a reference to the sigma2 array it is given: module lifetime
real(sp), allocatable, target :: sigma2_fixture(:,:)

contains

    subroutine run_all_pftc_inplane_tests()
        write(*,'(A)') '**** running all pftc in-plane tests ****'
        call write_phantom()
        call test_euclid_routes_and_recovery()
        call test_euclid_band_edges()
        call test_cc_evaluator()
        call test_hybrid_evaluator()
        call test_strategy2D_route_flags()
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
        complex(sp), allocatable :: den_pft(:,:)
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
            case(OBJ_HYBRID)
                call cline%set('objfun',       'euclid')
                call cline%set('objfun_den',   'yes')
                call cline%set('objfun_den_w', 0.37)
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
        if( objfun == OBJ_HYBRID )then
            ! the shifted raw particle doubles as the denoised fixture: isolates the
            ! hybrid calculus while exercising its storage, memoisation and normalisation
            den_pft = b%pftc%allocate_pft()
            call b%pftc%get_ptcl_pft(1, den_pft(:,p%kfromto(1):p%kfromto(2)))
            call b%pftc%set_ptcl_den_pft(1, den_pft)
        endif
        call b%pftc%memoize_refs
        call b%pftc%memoize_ptcls
        if( objfun == OBJ_HYBRID ) call b%pftc%memoize_ptcls_den
        if( allocated(sigma2_fixture) ) deallocate(sigma2_fixture)
        allocate(sigma2_fixture(p%kfromto(1):p%kfromto(2),1), source=1.0_sp)
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
            case(OBJ_HYBRID)
                call b%pftc%gen_hybrid_grad_at_angle(1, 1, shift, rotind, f, grad)
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
        ! (d) seed parity: the joint route's discrete selection is one legacy callback at the seed
        seed = [0.5, -0.25]
        call b%pftc%gen_objfun_vals(1, 1, seed, scores)
        expected_irot = maxloc(scores, dim=1)
        joint_limits(1:2,1) = -1.; joint_limits(1:2,2) = 1.
        joint_limits(3,:)   = [1.-2., real(nrots)+2.]
        call joint_search%new_joint(b, joint_limits, 100)
        call joint_search%set_indices(1, 1)
        call joint_search%select_best_discrete_angle(seed, selected, corr)
        call assert_false(joint_search%does_opt_angle(), 'joint seed selector does not attach the legacy callback')
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
        ! (f) stress: cavgs-like near-noiseless shell-dependent sigma2, first with stale
        !     memoised square sums (production updates sigma2 without re-memoising), then fresh
        do k = p%kfromto(1), p%kfromto(2)
            sigma2_fixture(k,1) = 1.e-4_sp / real(k,sp)
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

    !> hybrid (euclid with denoised term): capability, grid identity, gradient
    subroutine test_hybrid_evaluator()
        type(builder),    target :: b
        type(parameters), target :: p
        type(pftc_shsrch_grad) :: joint_search
        real(sp), allocatable :: scores(:)
        real     :: joint_limits(3,2)
        integer  :: nrots, igrid
        write(*,'(A)') 'test_hybrid_evaluator'
        call build_fixture(b, p, OBJ_HYBRID, LP_LOW, TRUTH_ANGLE, TRUTH_SHIFT, .false.)
        nrots = b%pftc%get_nrots()
        allocate(scores(nrots))
        call assert_true(b%pftc%is_hybrid_objfun(),     'objfun_den=yes is the hybrid objective')
        call assert_true(b%pftc%is_joint_grad_objfun(), 'the hybrid objective advertises joint gradients')
        call assert_false(b%pftc%is_raw_euclid_objfun(), 'the hybrid objective is not the raw Euclidean one')
        joint_limits(1:2,1) = -p%trs; joint_limits(1:2,2) = p%trs
        joint_limits(3,:)   = [1.-2., real(nrots)+2.]
        call joint_search%new_joint(b, joint_limits, p%maxits_sh)
        call joint_search%kill
        call grid_identity(b, OBJ_HYBRID, nrots, 'hybrid')
        call b%pftc%gen_objfun_vals(1, 1, [0._sp,0._sp], scores)
        igrid = maxloc(scores, dim=1)
        call gradient_probe_set(b, OBJ_HYBRID, igrid, 'hybrid')
        call kill_fixture(b)
    end subroutine test_hybrid_evaluator

    !> strategy2D_srch constructs the joint optimiser only under inpl_cont=yes and always
    !! keeps the legacy seed-search angle update (selection parity)
    subroutine test_strategy2D_route_flags()
        type(builder),    target :: b
        type(parameters), target :: p
        type(strategy2D_srch) :: srch
        type(strategy2D_spec) :: spec
        type(pftc_shsrch_grad) :: legacy
        real :: limits(2,2)
        write(*,'(A)') 'test_strategy2D_route_flags'
        call build_fixture(b, p, OBJ_EUCLID, LP_LOW, TRUTH_ANGLE, TRUTH_SHIFT, .false.)
        limits(:,1) = -1.; limits(:,2) = 1.
        call legacy%new_legacy(b, limits)
        call legacy%kill
        spec%iptcl       = 1
        spec%iptcl_batch = 1
        spec%iptcl_map   = 1
        p%l_prob_align_mode = .false.
        p%l_objfun_den      = .false.
        p%inpl_cont         = 'no'
        call srch%new(p, spec, b)
        call assert_false(srch%uses_continuous_refinement(), 'inpl_cont=no: no continuous polish')
        call assert_false(srch%joint_inpl_optimizer%uses_joint_inplane(), 'inpl_cont=no: no joint optimiser')
        call assert_true(srch%grad_shsrch_first_obj%does_opt_angle(), 'inpl_cont=no: legacy seed search updates the angle')
        call srch%kill
        p%inpl_cont = 'yes'
        call srch%new(p, spec, b)
        call assert_true(srch%uses_continuous_refinement(), 'inpl_cont=yes: continuous polish')
        call assert_true(srch%joint_inpl_optimizer%uses_joint_inplane(), 'inpl_cont=yes: joint optimiser constructed')
        call assert_true(srch%grad_shsrch_first_obj%does_opt_angle(), &
            &'inpl_cont=yes: legacy seed search still updates the angle (selection parity)')
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

end module simple_pftc_inplane_tester
