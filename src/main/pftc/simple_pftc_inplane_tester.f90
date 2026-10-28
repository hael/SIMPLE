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
use simple_eul_prob_tab_utils,      only: posterior_temperature, shift_prior_sigma2, SHIFT_PRIOR_SIGMA_MIN, SHIFT_PRIOR_LP_FRAC
use simple_oris,                    only: oris
use simple_test_utils
use, intrinsic :: ieee_arithmetic,  only: ieee_is_finite, ieee_value, ieee_quiet_nan
use, intrinsic :: iso_fortran_env,  only: int64
!$ use omp_lib, only: omp_get_max_threads, omp_get_num_procs, omp_get_thread_num, omp_set_num_threads
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
        call test_shmat_searched_shells()
        call test_search_threads()
        call test_loaded_reference()
        call test_posterior_threads()
        call test_shift_prior()
        call test_shift_prior_frame()
        call test_shift_prior_grid_fit()
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

    !> reference at e3 = 0 and one particle at (truth_angle, shift) from the phantom; with nrefs, references 2..nrefs
    !! are the phantom's projections along the other directions of the spiral
    subroutine build_fixture( b, p, objfun, lp, truth_angle, shift, hard_edge, nrefs, nthr )
        type(builder),    target, intent(inout) :: b
        type(parameters), target, intent(inout) :: p
        integer,                  intent(in)    :: objfun
        real,                     intent(in)    :: lp, truth_angle, shift(2)
        logical,                  intent(in)    :: hard_edge
        integer,        optional, intent(in)    :: nrefs, nthr
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
        if( present(nrefs) ) call cline%set('nspace', nrefs)
        if( present(nthr)  ) call cline%set('nthr',   nthr)
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
    !! at every rotation, for both objectives, and the posterior route on the same entry
    subroutine test_particle_bank()
        real(dp), parameter :: TODAY_ATOL = 5.d-2 ! the dropped Nyquist bin separates the two kernels; a frame error is O(1)
        type(builder),    target :: b
        type(parameters), target :: p
        complex(sp), allocatable :: bank(:,:)
        real(sp),    allocatable :: vals(:), vals_today(:), pvec(:)
        real(dp) :: max_err, max_today, direct
        real(sp) :: t(2), s(2), dist
        integer  :: iobj, ish, irot, irot_s, nrots, kfromto(2), irot_seed
        write(*,'(A)') 'test_particle_bank'
        do iobj = OBJ_EUCLID, OBJ_CC
            call build_fixture(b, p, iobj, LP_LOW, TRUTH_ANGLE, [2., -1.5], .false.)
            nrots   = b%pftc%get_nrots()
            kfromto = b%pftc%get_kfromto()
            allocate(bank(b%pftc%get_pftsz()+1, kfromto(2)-kfromto(1)+1), vals(nrots), vals_today(nrots), pvec(nrots))
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
            ! the posterior route (euclid): the bank row is the loss of the same entry at every rotation, and a draw
            ! at a cold temperature takes its best rotation
            if( iobj == OBJ_EUCLID )then
                call b%pftc%gen_euclid_dists_bank(1, 1, bank, pvec)
                call assert_true(maxval(abs(exp(-pvec) - vals)) <= 1.e-4, 'the bank row is -log of the bank values')
                call b%pftc%gen_posterior_val_bank(1, 1, bank, huge(1._sp), dist, irot_s)
                call assert_int(maxloc(vals, dim=1), irot_s, 'the cold posterior draw picks the best rotation of the entry')
                call assert_true(abs(exp(-dist) - vals(irot_s)) <= 1.e-4, 'and reports its loss')
            endif
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
            deallocate(bank, vals, vals_today, pvec)
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

    ! ---- exact checks: shift matrix range, loaded reference, the tables' kernels across threads ----

    !> the alignment shift matrix holds only the searched shells: the scores are the same, bit for bit, whatever the
    !! shells above kfromto(2) hold (shift_ptcl fills them with another shift; the searched shells are then reset)
    subroutine test_shmat_searched_shells()
        type(builder),    target :: b
        type(parameters), target :: p
        complex(sp), allocatable :: pft(:,:), pft_full(:,:)
        real(sp),    allocatable :: v1(:), v2(:)
        integer :: iobj, ish, nrots, kfromto(2), ndiff
        write(*,'(A)') 'test_shmat_searched_shells'
        do iobj = OBJ_EUCLID, OBJ_CC
            call build_fixture(b, p, iobj, LP_LOW, TRUTH_ANGLE, TRUTH_SHIFT, .false.)
            nrots   = b%pftc%get_nrots()
            kfromto = b%pftc%get_kfromto()
            allocate(v1(nrots), v2(nrots), pft(b%pftc%get_pftsz(),kfromto(1):kfromto(2)))
            pft_full = b%pftc%allocate_ptcl_pft()
            call b%pftc%get_ptcl_pft(1, pft)
            pft_full(:,kfromto(1):kfromto(2)) = pft
            ndiff = 0
            do ish = 1,size(PROBE_SHIFTS,2)
                call reset_particle
                call b%pftc%gen_objfun_vals(1, 1, real(PROBE_SHIFTS(:,ish),sp), v1)
                call b%pftc%shift_ptcl(1, [3.1, -2.7] * real(ish))
                call reset_particle
                call b%pftc%gen_objfun_vals(1, 1, real(PROBE_SHIFTS(:,ish),sp), v2)
                ndiff = ndiff + count(v1 /= v2)
            enddo
            call assert_int(0, ndiff, 'the scores do not depend on the shells above the searched band (exact)')
            deallocate(v1, v2, pft, pft_full)
            call kill_fixture(b)
        enddo

    contains

        !> the searched shells back, nothing above them, and the particle memoised afresh
        subroutine reset_particle
            call b%pftc%set_ptcl_pft(1, pft_full)
            call b%pftc%memoize_ptcls
        end subroutine reset_particle

    end subroutine test_shmat_searched_shells

    !> the evaluations of the shift search with the reference loaded once equal those that copy it every time and those
    !! after a fresh load, exactly, also when other evaluations and loads come in between
    subroutine test_loaded_reference()
        integer, parameter :: NREFS_LR = 6
        type(builder),    target :: b
        type(parameters), target :: p
        real(sp), allocatable :: vals(:)
        real(dp) :: f_load, f_copy, g_load(2), g_copy(2), c_load, c_copy
        integer  :: iobj, iref, ish, irot, nrots, ndiff
        write(*,'(A)') 'test_loaded_reference'
        do iobj = OBJ_EUCLID, OBJ_CC
            call build_fixture(b, p, iobj, LP_LOW, TRUTH_ANGLE, TRUTH_SHIFT, .false., nrefs=NREFS_LR)
            nrots = b%pftc%get_nrots()
            allocate(vals(nrots))
            ndiff = 0
            do iref = 1,NREFS_LR
                call b%pftc%load_ref_aln_8(iref, 1)
                do ish = 1,size(PROBE_SHIFTS,2)
                    do irot = 1,nrots,7
                        ! evaluations of another reference in between
                        call b%pftc%gen_objfun_vals(mod(iref,NREFS_LR)+1, 1, real(PROBE_SHIFTS(:,ish),sp), vals)
                        c_copy = b%pftc%gen_corr_for_rot_8(mod(iref,NREFS_LR)+1, 1, PROBE_SHIFTS(:,ish), irot)
                        c_load = b%pftc%gen_corr_for_rot_8_aln(iref, 1, PROBE_SHIFTS(:,ish), irot)
                        c_copy = b%pftc%gen_corr_for_rot_8(iref, 1, PROBE_SHIFTS(:,ish), irot)
                        call b%pftc%gen_corr_grad_for_rot_8_aln(iref, 1, PROBE_SHIFTS(:,ish), irot, f_load, g_load)
                        ! another reference loaded and evaluated, then this one loaded afresh
                        call b%pftc%load_ref_aln_8(mod(iref,NREFS_LR)+1, 1)
                        call b%pftc%gen_corr_grad_for_rot_8_aln(mod(iref,NREFS_LR)+1, 1, PROBE_SHIFTS(:,ish), irot,&
                            &f_copy, g_copy)
                        call b%pftc%load_ref_aln_8(iref, 1)
                        call b%pftc%gen_corr_grad_for_rot_8_aln(iref, 1, PROBE_SHIFTS(:,ish), irot, f_copy, g_copy)
                        if( c_load /= c_copy .or. f_load /= f_copy .or. any(g_load /= g_copy) ) ndiff = ndiff + 1
                    enddo
                enddo
            enddo
            call assert_int(0, ndiff, 'the loaded reference evaluates exactly as the copied one')
            deallocate(vals)
            call kill_fixture(b)
        enddo
    end subroutine test_loaded_reference

    !> the kernels the probability tables score with, over 24 references at one thread and 20 times at every
    !! processor: the rows of gen_euclid_dists (shifted-reference kernel) and gen_euclid_dists_bank, the posterior
    !! draws at a cold temperature (the row's argmin, deterministic) and the SHC strategies' single-sample draws, bit
    !! for bit the same; and a draw at the particle's own temperature always returns a rotation of its row at the
    !! row's value. cc has no posterior kernel: its rows are gen_objfun_vals and gen_objfun_vals_bank
    subroutine test_posterior_threads()
        integer, parameter :: NREFS_LT = 24, NRES = 10
        real(sp), parameter :: COLD = huge(1._sp)
        type(builder),    target :: b
        type(parameters), target :: p
        complex(sp), allocatable :: bank(:,:)
        real(sp),    allocatable :: res1(:,:), resn(:,:), rows1(:,:,:), rowsn(:,:,:)
        integer :: iobj, nthr_prev, nthr_full, kfromto(2), irep, ndiff, nrots, nbad
        write(*,'(A)') 'test_posterior_threads'
        nthr_prev = 1
        nthr_full = 1
        !$ nthr_prev = omp_get_max_threads()
        !$ nthr_full = omp_get_num_procs()
        allocate(res1(NRES,NREFS_LT), resn(NRES,NREFS_LT))
        do iobj = OBJ_EUCLID, OBJ_CC
            call build_fixture(b, p, iobj, LP_LOW, TRUTH_ANGLE, TRUTH_SHIFT, .false., nrefs=NREFS_LT, nthr=nthr_full)
            kfromto = b%pftc%get_kfromto()
            nrots   = b%pftc%get_nrots()
            allocate(bank(b%pftc%get_pftsz()+1, kfromto(2)-kfromto(1)+1), rows1(nrots,2,NREFS_LT), rowsn(nrots,2,NREFS_LT))
            call b%pftc%memoize_shifted_ptcl(1, [0.7, -0.4], bank)
            !$ call omp_set_num_threads(1)
            call score_refs(res1, rows1)
            !$ call omp_set_num_threads(nthr_full)
            ndiff = 0
            nbad  = 0
            do irep = 1,20
                call score_refs(resn, rowsn)
                ndiff = ndiff + count(transfer(res1,0,size(res1)) /= transfer(resn,0,size(resn)))
                ndiff = ndiff + count(transfer(rows1,0,size(rows1)) /= transfer(rowsn,0,size(rowsn)))
                if( iobj == OBJ_EUCLID ) nbad = nbad + draws_off_row(rowsn)
            enddo
            !$ call omp_set_num_threads(nthr_prev)
            call assert_int(0, ndiff, 'the rows, the cold draws and the single-sample draws at every processor equal those at one, bit for bit')
            if( iobj == OBJ_EUCLID ) call assert_int(0, nbad, 'a draw at the particle temperature returns a rotation of the row at its value')
            deallocate(bank, rows1, rowsn)
            call kill_fixture(b)
        enddo
        nthr_glob = nthr_prev

    contains

        !> res: the cold draws (dist, irot) of the shifted and the bank kernel, then (euclid) the SHC strategies'
        !! single-sample draws; rows: the two kernels' rows
        subroutine score_refs( res, rows )
            real(sp), intent(out) :: res(:,:), rows(:,:,:)
            real(sp) :: pvec(nrots), dist, corr
            integer  :: inds(nrots), iref, irot
            res = 0.
            !$omp parallel do default(shared) private(iref,pvec,inds,dist,corr,irot) schedule(dynamic)
            do iref = 1,NREFS_LT
                if( iobj == OBJ_EUCLID )then
                    call b%pftc%gen_euclid_dists(iref, 1, [1.1, -0.6], rows(:,1,iref))
                    call b%pftc%gen_euclid_dists_bank(iref, 1, bank, rows(:,2,iref))
                    call b%pftc%gen_posterior_val(iref, 1, [1.1, -0.6], COLD, dist, irot)
                    res(1:2,iref) = [dist, real(irot)]
                    call b%pftc%gen_posterior_val_bank(iref, 1, bank, COLD, dist, irot)
                    res(3:4,iref) = [dist, real(irot)]
                    call b%pftc%gen_prob_power_objfun_val(iref, 1, [1.1, -0.6], 2., 1, dist, corr, irot, pvec, inds)
                    res(5:7,iref) = [dist, corr, real(irot)]
                    call b%pftc%gen_prob_objfun_val(iref, 1, [1.1, -0.6], 0.1, 0.1, dist, irot, pvec, inds)
                    res(8:9,iref) = [dist, real(irot)]
                    res(10,iref)  = b%pftc%ptcl_posterior_temp(1)
                else
                    call b%pftc%gen_objfun_vals(iref, 1, [1.1, -0.6], rows(:,1,iref))
                    call b%pftc%gen_objfun_vals_bank(iref, 1, bank, rows(:,2,iref))
                endif
            enddo
            !$omp end parallel do
        end subroutine score_refs

        !> draws at the particle's temperature over the shifted kernel: the number whose rotation is outside the
        !! row or whose value is not the row's at that rotation
        integer function draws_off_row( rows ) result( nbad_loc )
            real(sp), intent(in) :: rows(:,:,:)
            real(sp) :: dist, temp
            integer  :: iref, irot
            temp     = b%pftc%ptcl_posterior_temp(1)
            nbad_loc = 0
            !$omp parallel do default(shared) private(iref,dist,irot) schedule(dynamic) reduction(+:nbad_loc)
            do iref = 1,NREFS_LT
                call b%pftc%gen_posterior_val(iref, 1, [1.1, -0.6], temp, dist, irot)
                if( irot < 1 .or. irot > nrots )then
                    nbad_loc = nbad_loc + 1
                else if( dist /= rows(irot,1,iref) )then
                    nbad_loc = nbad_loc + 1
                endif
            enddo
            !$omp end parallel do
        end function draws_off_row

    end subroutine test_posterior_threads

    !> The shift prior of the probability tables' searches (shift_prior=yes): off, the search is the search
    !! (bitwise, also after a prior was set and cleared); on, the correlation a search reports is the raw
    !! correlation at its shift times exp(-|u - u0|^2 / (2 sigma^2 T)) with u0 the start as a stored shift and u the
    !! candidate's stored shift at its rotation (the frame rule of test_shift_prior_frame), so the in-plane choice
    !! is the argmax of the raw objective times the rotation's factor; a prior far wider than the move leaves the
    !! free search's objective (its value, not its path: a factor within L-BFGS-B's stopping tolerance can end an
    !! alternation elsewhere) and a prior far narrower pins the stored shift at its start, on the shift route and
    !! on the joint polish alike. And the width from a field: half the mean squared increment of the sampled
    !! particles, floored by the band
    subroutine test_shift_prior()
        real, parameter :: START(2) = [1.5, -1.0]
        type(builder),    target :: b
        type(parameters), target :: p
        type(pftc_shsrch_grad)   :: srch
        type(oris) :: os
        real(dp) :: w, g, rotind_frac, cost_free, cost_prior, best
        real     :: limits(2,2), lims_init(2,2), joint_lims(3,2), temp, cxy_free(3), cxy_again(3), cxy(3)
        real     :: corr_prior, sigma2, sigma2_weak, sigma2_strong, u0(2), u0_joint(2)
        real,    allocatable :: scores(:)
        integer  :: irot_free, irot_again, irot, irot2, irot3, irot_seed, i
        logical  :: l_valid, l_improved
        write(*,'(A)') 'test_shift_prior'
        limits(:,1)    = -TRS
        limits(:,2)    =  TRS
        lims_init(:,1) = -SHC_INPL_TRSHWDTH
        lims_init(:,2) =  SHC_INPL_TRSHWDTH
        call build_fixture(b, p, OBJ_EUCLID, LP_LOW, TRUTH_ANGLE, TRUTH_SHIFT, .false.)
        temp = posterior_temperature(b%pftc%wsqsum_ptcl_now(1), b%pftc%get_pftsz())
        call assert_true(temp > 0., 'the fixture particle has a temperature')
        allocate(scores(b%pftc%get_nrots()))
        ! the widths in the fixture's own units: w = 1/(2 sigma2 T); weak: w (10 px)^2 < 1e-6 nats, strong:
        ! w (0.1 px)^2 > 5 nats
        sigma2_weak   = 1.e8 / temp
        sigma2_strong = 1.e-3 / temp
        call srch%new_alternating(b, limits, lims_init=lims_init, maxits=p%maxits_sh, coarse_init=.true.)
        call srch%set_indices(1, 1)
        call assert_false(srch%has_shift_prior(), 'a new search object has no prior')
        ! the free search from the start, twice: the second after a prior was set and cleared
        irot_free = 1
        cxy_free  = srch%minimize(irot_free, sh_rot=.false., xy_in=START)
        call assert_true(irot_free > 0, 'the free search moves from the start')
        call assert_true(norm2(cxy_free(2:3) - START) > 0.5, 'and moves it by more than half a pixel')
        call srch%set_shift_prior(1., temp)
        call assert_true(srch%has_shift_prior(), 'set_shift_prior arms the prior')
        call srch%clear_shift_prior
        irot_again = 1
        cxy_again  = srch%minimize(irot_again, sh_rot=.false., xy_in=START)
        call assert_true(irot_again == irot_free .and. all(transfer(cxy_again,0,3) == transfer(cxy_free,0,3)), &
            &'a cleared prior leaves the search bit for bit as it was')
        ! the factor identity at the prior's own result: the centre is START as a stored shift (expressed at
        ! rotation 1 on entry), so at rotation i the factor is exp(-w |xy - R(i) u0|^2); the in-plane choice is
        ! the argmax of the raw objective times that factor and the reported correlation is the product
        call srch%set_shift_prior(1., temp)
        irot = 1
        cxy  = srch%minimize(irot, sh_rot=.false., xy_in=START)
        call assert_true(irot > 0, 'the search with a unit prior moves')
        call srch%select_best_discrete_angle(cxy(2:3), irot2, corr_prior)
        call b%pftc%gen_objfun_vals(1, 1, cxy(2:3), scores)
        w     = 1.d0 / (2.d0 * 1.d0 * real(temp,dp))
        u0    = shift_ref2ptcl(START, b%pftc%get_rot(1))
        best  = -huge(best)
        irot3 = 0
        do i = 1,b%pftc%get_nrots()
            g = exp(-w * sum((real(cxy(2:3),dp) - real(shift_ptcl2ref(u0, b%pftc%get_rot(i)),dp))**2))
            if( real(scores(i),dp) * g > best )then
                best  = real(scores(i),dp) * g
                irot3 = i
            endif
        enddo
        call assert_int(irot3, irot2, 'the prior''s in-plane choice is the argmax of the raw objective times the rotation''s factor')
        call assert_true(abs(real(corr_prior,dp) - best) <= 1.d-5, &
            &'the reported correlation is the raw one times the prior factor at its rotation')
        call assert_true(abs(cxy(1) - corr_prior) <= 1.e-5, 'and minimize returns that correlation')
        ! a wide prior changes nothing that matters, a narrow one pins the stored shift at the start
        call srch%set_shift_prior(sigma2_weak, temp)
        irot = 1
        cxy  = srch%minimize(irot, sh_rot=.false., xy_in=START)
        call assert_true(irot > 0 .and. abs(cxy(1) - cxy_free(1)) <= 1.e-4 .and. norm2(cxy(2:3) - cxy_free(2:3)) < 0.25, &
            &'a prior far wider than the move leaves the free search''s objective, the shift within a quarter pixel')
        call srch%set_shift_prior(sigma2_strong, temp)
        irot = 1
        cxy  = srch%minimize(irot, sh_rot=.true., xy_in=START)
        call assert_true(irot == 0 .or. norm2(cxy(2:3) - u0) < 0.1, &
            &'a prior far narrower than the move pins the stored shift at its start, whatever the rotation')
        call srch%kill
        ! the joint route (the polish of the committed pose): the prior is the same loss term; its initial cost
        ! (the start re-scored) is untouched, a wide prior leaves the free polish's objective, a narrow one pins
        ! the start. The seed cell is the discrete selection at the start, as the polish is seeded in production
        call b%pftc%gen_objfun_vals(1, 1, START, scores)
        irot_seed = maxloc(scores, dim=1)
        u0_joint  = shift_ref2ptcl(START, b%pftc%get_rot(irot_seed))
        joint_lims(1:2,:) = limits
        joint_lims(3,:)   = [1.-2., real(b%pftc%get_nrots())+2.]
        call srch%new_joint(b, joint_lims, p%maxits_sh)
        call srch%set_indices(1, 1)
        cxy_free = srch%minimize_joint(irot_free, START, sh_rot=.false., rotind_frac=rotind_frac, &
            &evaluation_valid=l_valid, improved=l_improved, initial_cost_out=cost_free, irot_in=irot_seed)
        call assert_true(l_valid .and. l_improved, 'the free joint polish improves from the start')
        call srch%set_shift_prior(sigma2_weak, temp)
        cxy = srch%minimize_joint(irot, START, sh_rot=.false., rotind_frac=rotind_frac, &
            &evaluation_valid=l_valid, improved=l_improved, initial_cost_out=cost_prior, irot_in=irot_seed)
        call assert_true(abs(cost_prior - cost_free) <= 1.d-12 * max(1.d0, abs(cost_free)), &
            &'the joint prior leaves the cost at the start untouched')
        call assert_true(l_valid .and. l_improved .and. abs(cxy(1) - cxy_free(1)) <= 1.e-4 .and. &
            &norm2(cxy(2:3) - cxy_free(2:3)) < 0.25, &
            &'a wide prior leaves the free joint polish''s objective, the shift within a quarter pixel')
        call srch%set_shift_prior(sigma2_strong, temp)
        cxy = srch%minimize_joint(irot, START, sh_rot=.true., rotind_frac=rotind_frac, &
            &evaluation_valid=l_valid, improved=l_improved, irot_in=irot_seed)
        call assert_true((.not. l_improved) .or. norm2(cxy(2:3) - u0_joint) < 0.1, &
            &'a narrow prior pins the joint polish''s stored shift at its start')
        call srch%kill
        call kill_fixture(b)
        ! the width from a field: five sampled particles with increments 2, 0, 4 and two never searched (zero,
        ! which a particle record cannot tell from a search that kept its shift), one unsampled that must not count
        call os%new(6, is_ptcl=.true.)
        do i = 1,6
            call os%set(i, 'sampled', 1.)
        enddo
        call os%set(6, 'sampled', 0.)
        call os%set(1, 'shincarg', 2.)
        call os%set(2, 'shincarg', 0.)
        call os%set(3, 'shincarg', 4.)
        call os%set(6, 'shincarg', 100.)
        sigma2 = shift_prior_sigma2(os, 2.)
        call assert_true(abs(sigma2 - 0.5 * (4. + 0. + 16.) / 5.) <= 1.e-5, &
            &'the width is half the mean squared increment of the sampled particles')
        sigma2 = shift_prior_sigma2(os, 40.)
        call assert_true(abs(sigma2 - (SHIFT_PRIOR_LP_FRAC * 40.)**2) <= 1.e-4, &
            &'and never below a quarter of the matching wavelength')
        call os%set(1, 'shincarg', 0.5)
        call os%set(3, 'shincarg', 0.5)
        sigma2 = shift_prior_sigma2(os, 2.)
        call assert_true(abs(sigma2 - SHIFT_PRIOR_SIGMA_MIN**2) <= 1.e-6, 'nor below a pixel')
        call os%kill
        deallocate(scores)
    end subroutine test_shift_prior

    !> The frame of the shift prior (2026-10-28): the optimizer's shift x at rotation irot is the stored shift
    !! rotated into the reference's frame, x = R(irot) u, so the prior's centre is kept as the stored shift u0 and
    !! the penalty at rotation irot is w |x - R(irot) u0|^2 = w |u - u0|^2. With a nonzero seed u0 expressed at
    !! irot_seed: no penalty at the seed; none at another rotation for the x that keeps the stored shift; the
    !! stored displacement's cost, w (2 |u0| sin(dtheta/2))^2, for the seed's own x at another rotation (which
    !! paid nothing before, the centre being x itself); the particle-frame penalty of that candidate's stored
    !! shift is the same number; a zero seed is the same centre in every frame. And the search: a narrow prior
    !! keeps the stored shift at the seed whatever rotation it ends on
    subroutine test_shift_prior_frame()
        real, parameter :: U0(2) = [2.0, -1.5]
        type(builder),    target :: b
        type(parameters), target :: p
        type(pftc_shsrch_grad)   :: srch
        real     :: limits(2,2), temp, xy_seed(2), xy_other(2), rot_seed, rot_other, cxy(3), sigma2_strong
        real(dp) :: w, pen, expected, u(2), scale
        integer  :: irot_seed, irot_other, irot, nrots
        write(*,'(A)') 'test_shift_prior_frame'
        limits(:,1) = -TRS
        limits(:,2) =  TRS
        call build_fixture(b, p, OBJ_EUCLID, LP_LOW, TRUTH_ANGLE, TRUTH_SHIFT, .false.)
        temp  = posterior_temperature(b%pftc%wsqsum_ptcl_now(1), b%pftc%get_pftsz())
        nrots = b%pftc%get_nrots()
        call srch%new_alternating(b, limits, maxits=p%maxits_sh)
        call srch%set_indices(1, 1)
        call srch%set_shift_prior(1., temp)
        w     = 1.d0 / (2.d0 * real(temp,dp))
        scale = w * sum(real(U0,dp)**2) ! the penalty's natural size
        irot_seed  = 1 + nrots / 8
        irot_other = 1 + nrots / 3
        rot_seed   = b%pftc%get_rot(irot_seed)
        rot_other  = b%pftc%get_rot(irot_other)
        xy_seed    = shift_ptcl2ref(U0, rot_seed)   ! the seed as the caller passes it: expressed at its rotation
        xy_other   = shift_ptcl2ref(U0, rot_other)  ! the same stored shift at another rotation
        call srch%set_prior_centre(real(xy_seed,dp), irot_seed)
        call assert_true(srch%shift_prior_penalty(real(xy_seed,dp), irot_seed) <= 1.d-5 * scale, &
            &'no penalty at the seed')
        call assert_true(srch%shift_prior_penalty(real(xy_other,dp), irot_other) <= 1.d-5 * scale, &
            &'no penalty at another rotation for the shift that keeps the stored shift')
        pen      = srch%shift_prior_penalty(real(xy_seed,dp), irot_other)
        expected = w * (2.d0 * norm2(real(U0,dp)) * sin(0.5d0 * real(rot_other - rot_seed,dp) * DPI / 180.d0))**2
        call assert_true(pen > 0.1d0 * scale, 'the seed''s own shift at another rotation pays: the stored shift moved')
        call assert_true(abs(pen - expected) <= 1.d-5 * scale, 'and pays the stored displacement, 2 |u0| sin(dtheta/2)')
        u = real(shift_ref2ptcl(xy_seed, rot_other), dp)
        call assert_true(abs(srch%shift_prior_penalty_ptcl(u) - pen) <= 1.d-5 * scale, &
            &'the particle-frame penalty of that candidate''s stored shift is the same number')
        call assert_true(srch%shift_prior_penalty_ptcl(real(U0,dp)) <= 1.d-5 * scale, &
            &'and zero at the seed''s stored shift')
        call srch%set_prior_centre([0.d0, 0.d0], irot_seed)
        call assert_true(srch%shift_prior_penalty([0.d0, 0.d0], irot_other) <= 1.d-12, &
            &'a zero seed is the same centre at every rotation')
        pen = srch%shift_prior_penalty(real(xy_seed,dp), irot_other)
        call assert_true(abs(pen - w * sum(real(xy_seed,dp)**2)) <= 1.d-5 * scale, &
            &'and its penalty is the plain squared shift')
        ! the search from the seed with a narrow prior (w (0.1 px)^2 > 5 nats) keeps the stored shift at the seed,
        ! whatever rotation it ends on (the returned shift is the stored one with sh_rot)
        sigma2_strong = 1.e-3 / temp
        call srch%set_shift_prior(sigma2_strong, temp)
        irot = irot_seed
        cxy  = srch%minimize(irot, sh_rot=.true., xy_in=xy_seed)
        call assert_true(irot == 0 .or. norm2(cxy(2:3) - U0) < 0.1, &
            &'a narrow prior keeps the stored shift at a nonzero seed whatever the rotation the search ends on')
        call srch%kill
        call kill_fixture(b)
    end subroutine test_shift_prior_frame

    !> The shift prior through the quadratic grid step (grid_fit, the multi-state table's search), whose exact
    !! re-score at the fitted shift is its own evaluation, not the all-rotation objective's: on both lattices (the
    !! reference frame and the particle frame) the correlation a search reports is the raw correlation at its shift
    !! times the prior's factor exp(-w |s|^2) (centred at zero, a search without a start), whether the fit was kept
    !! or fell back to the alternations; and the fixture keeps a fit on each lattice, so the step's re-score is
    !! exercised. w = 0.01 per px^2 is weak on the noise-free phantom (the search moves as without it) but its
    !! factor, 0.3% to 2% at the cases' shifts, is far above the tolerance
    subroutine test_shift_prior_grid_fit()
        real(sp), parameter :: TRUE_SHIFTS(2,3) = reshape([0.5, 0.3, 1.2, -0.8, 4.6, 0.0], [2,3])
        real(dp), parameter :: W_PRIOR = 1.d-2 ! the prior's weight per px^2, in loss units
        type(builder),    target :: b
        type(parameters), target :: p
        type(pftc_shsrch_grad)   :: srch
        complex(sp), allocatable :: bank(:,:)
        real,        allocatable :: vals(:)
        integer(int64) :: counts(0:5)
        real     :: limits(2,2), lims_init(2,2), cxy(3), temp, u(2)
        real(dp) :: g, raw
        integer  :: icase, ilat, irot, nkept(2), kfromto(2)
        logical  :: l_ptcl, l_kept
        write(*,'(A)') 'test_shift_prior_grid_fit'
        limits(:,1)    = -TRS
        limits(:,2)    =  TRS
        lims_init(:,1) = -SHC_INPL_TRSHWDTH
        lims_init(:,2) =  SHC_INPL_TRSHWDTH
        nkept = 0
        do ilat = 1,2
            l_ptcl = ilat == 2
            do icase = 1,size(TRUE_SHIFTS,2)
                call build_fixture(b, p, OBJ_EUCLID, LP_LOW, TRUTH_ANGLE, TRUE_SHIFTS(:,icase), .false.)
                kfromto = b%pftc%get_kfromto()
                allocate(bank(b%pftc%get_pftsz()+1, kfromto(2)-kfromto(1)+1), vals(b%pftc%get_nrots()))
                temp = posterior_temperature(b%pftc%wsqsum_ptcl_now(1), b%pftc%get_pftsz())
                call srch%new_alternating(b, limits, lims_init=lims_init, maxits=60, coarse_init=.true., &
                    &grid_fit=.true., grid_ptcl=l_ptcl)
                call srch%set_indices(1, 1)
                ! w = 1/(2 sigma2 T) = W_PRIOR
                call srch%set_shift_prior(real(1.d0 / (2.d0 * W_PRIOR * real(temp,dp))), temp)
                irot = 1
                cxy  = srch%minimize(irot, sh_rot=.false.)
                call srch%get_fit_profile(counts)
                call assert_int(1, int(sum(counts)), 'each prior search takes exactly one fit decision')
                l_kept = counts(0) == 1_int64
                if( l_kept ) nkept(ilat) = nkept(ilat) + 1
                if( irot > 0 )then
                    ! the raw correlation at the reported shift by the kernel that scored it: a kept particle-frame
                    ! fit through the shifted particle, every other result through the shifted reference
                    if( l_kept .and. l_ptcl )then
                        u = shift_ref2ptcl(cxy(2:3), b%pftc%get_rot(irot))
                        call b%pftc%memoize_shifted_ptcl(1, u, bank)
                        call b%pftc%gen_objfun_vals_bank(1, 1, bank, vals)
                    else
                        call b%pftc%gen_objfun_vals(1, 1, cxy(2:3), vals)
                    endif
                    raw = real(vals(irot),dp)
                    g   = exp(-W_PRIOR * sum(real(cxy(2:3),dp)**2))
                    call assert_true(abs(real(cxy(1),dp) - raw * g) <= 1.d-5, &
                        &'the grid search reports the raw correlation at its shift times the prior factor')
                endif
                write(*,'(A,L1,A,I0,A,L1,A,I0,A,2F8.3)') '  particle frame ', l_ptcl, ', case ', icase, &
                    &': fit kept ', l_kept, ', irot ', irot, ', shift ', cxy(2:3)
                call srch%kill
                deallocate(bank, vals)
                call kill_fixture(b)
            enddo
        enddo
        call assert_true(nkept(1) > 0, 'a kept fit on the reference-frame lattice exercises the step''s re-score')
        call assert_true(nkept(2) > 0, 'a kept fit on the particle-frame lattice exercises the step''s re-score')
    end subroutine test_shift_prior_grid_fit

    !> the alternating shift search (as the probabilistic tables run it) over 24 references and three start shifts,
    !! euclid and cc, at one thread and 20 times at every processor: every result bit for bit the same
    subroutine test_search_threads()
        integer,  parameter :: NREFS_TH = 24, NSTARTS = 3
        real,     parameter :: STARTS(2,NSTARTS) = reshape([0.,0., 1.5,-1., -2.,0.5], [2,NSTARTS])
        type(builder),    target :: b
        type(parameters), target :: p
        type(pftc_shsrch_grad), allocatable :: srch(:)
        real,    allocatable :: cxy1(:,:,:), cxy4(:,:,:)
        integer, allocatable :: irot1(:,:), irot4(:,:)
        real    :: limits(2,2), lims_init(2,2)
        integer :: iobj, nthr_prev, ithr, nthr_full, irep, ndiff
        write(*,'(A)') 'test_search_threads'
        limits(:,1)    = -TRS
        limits(:,2)    =  TRS
        lims_init(:,1) = -SHC_INPL_TRSHWDTH
        lims_init(:,2) =  SHC_INPL_TRSHWDTH
        nthr_prev = 1
        nthr_full = 1
        !$ nthr_prev = omp_get_max_threads()
        !$ nthr_full = omp_get_num_procs()
        allocate(srch(nthr_full), cxy1(3,NSTARTS,NREFS_TH), cxy4(3,NSTARTS,NREFS_TH), irot1(NSTARTS,NREFS_TH),&
            &irot4(NSTARTS,NREFS_TH))
        do iobj = OBJ_EUCLID, OBJ_CC
            ! the pftc's per-thread work arrays are made for the full team
            call build_fixture(b, p, iobj, LP_LOW, TRUTH_ANGLE, TRUTH_SHIFT, .false., nrefs=NREFS_TH, nthr=nthr_full)
            do ithr = 1,nthr_full
                call srch(ithr)%new_alternating(b, limits, lims_init=lims_init, maxits=p%maxits_sh, coarse_init=.true.)
            enddo
            !$ call omp_set_num_threads(1)
            call run_searches(cxy1, irot1)
            !$ call omp_set_num_threads(nthr_full)
            ndiff = 0
            do irep = 1,20
                call run_searches(cxy4, irot4)
                ndiff = ndiff + count(irot1 /= irot4) + count(transfer(cxy1,0,size(cxy1)) /= transfer(cxy4,0,size(cxy4)))
            enddo
            !$ call omp_set_num_threads(nthr_prev)
            call assert_int(0, ndiff, 'the shift searches at every processor equal those at one thread, bit for bit')
            call assert_true(count(irot1 > 0) > NREFS_TH, 'the searches find solutions')
            do ithr = 1,nthr_full
                call srch(ithr)%kill
            enddo
            call kill_fixture(b)
        enddo
        nthr_glob = nthr_prev

    contains

        subroutine run_searches( cxy, irot )
            real,    intent(out) :: cxy(:,:,:)
            integer, intent(out) :: irot(:,:)
            integer :: iref, istart, ith
            cxy = 0.
            !$omp parallel do default(shared) private(iref,istart,ith) schedule(dynamic)
            do iref = 1,NREFS_TH
                ith = 1
                !$ ith = omp_get_thread_num() + 1
                call srch(ith)%set_indices(iref, 1)
                do istart = 1,NSTARTS
                    irot(istart,iref) = 1
                    if( istart == 1 )then
                        cxy(:,istart,iref) = srch(ith)%minimize(irot(istart,iref), sh_rot=.true.)
                    else
                        cxy(:,istart,iref) = srch(ith)%minimize(irot(istart,iref), sh_rot=.true., xy_in=STARTS(:,istart))
                    endif
                    if( irot(istart,iref) == 0 ) cxy(:,istart,iref) = 0.
                enddo
            enddo
            !$omp end parallel do
        end subroutine run_searches

    end subroutine test_search_threads

end module simple_pftc_inplane_tester
