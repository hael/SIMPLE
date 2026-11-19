!@descr: unit tests of the continuous view refinement on polar sections (simple_polarft_cont_calc): T1, T2 and T4 of polar_continuous_view.md
! The 4-blob phantom of simple_pftc_inplane_tester, unmasked (its transform is analytic only without the mask). T1/T2:
! the B-spline section and its normal derivative at the 20 directions of the nspace=20 spiral against the analytic
! transform, per shell, within twice the KB bank's errors frozen in Phase 0 (the plan's Progress). T4: the five-variable
! gradient against central differences, the Gauss-Newton matrix against the central-difference Hessian at zero residual,
! each with a negative control; thread consistency of the per-thread cache; ns per section, B-spline against KB.
module simple_polarft_cont_calc_tester
use simple_core_module_api,         only: dp, sp, DTWOPI, DPI, OSMPL_PAD_FAC, del_file
use simple_string,                  only: string
use simple_cmdline,                 only: cmdline
use simple_parameters,              only: parameters
use simple_builder,                 only: builder
use simple_image,                   only: image
use simple_matcher_smpl_and_lplims, only: set_bp_range3D
use simple_polarft_calc,            only: vol_pad2ref_pfts_opt
use simple_polarft_cont_calc,       only: polarft_cont_calc
use simple_ori_utils,               only: euler2m
use simple_test_polar_fixture,      only: shift_args, rotate_half_ring
use simple_test_utils
use, intrinsic :: ieee_arithmetic,  only: ieee_is_finite
use, intrinsic :: iso_fortran_env,  only: int64
!$ use omp_lib,                     only: omp_get_max_threads, omp_set_num_threads, omp_get_num_threads
implicit none
private
public :: run_all_polarft_cont_calc_tests

! the phantom of simple_pftc_inplane_tester: box 64 at 1.3 A, four Gaussian blobs
integer,          parameter :: BOX     = 64
real,             parameter :: SMPD    = 1.3
real,             parameter :: MSKDIAM = 60.
real,             parameter :: TRS     = 5.
real,             parameter :: LP_LOW  = 8.
real,             parameter :: LP_FULL = 2.7
character(len=*), parameter :: PHANTOM_FILE = 'polarft_cont_phantom.mrc'
integer,          parameter :: NBLOBS = 4
real,             parameter :: CTRS(3,NBLOBS) = reshape([&
    &-9.0, -5.0,  3.0, &
    & 7.0,  8.0, -5.0, &
    & 0.0,-10.0, -8.0, &
    & 5.0, -3.0, 10.0], [3,NBLOBS])
real,             parameter :: SIGMAS(NBLOBS) = [3.5, 4.0, 3.0, 3.8]
real,             parameter :: AMPS(NBLOBS)   = [1.0, 0.8, 0.6, 0.5]
! T1/T2 (frozen in Phase 0, polar_continuous_view.md Progress): the 20 directions of the nspace=20 spiral, the judged
! shells 2-17, and the KB bank's per-shell relative RMS errors of the section and of its normal derivative; floor 2x KB
integer,          parameter :: NSPACE_T12 = 20
integer,          parameter :: KJUDGE(2)  = [2, 17]
real(dp),         parameter :: KB_SECT_ERR(2:17) = [5.0326d-02, 6.7328d-02, 5.4585d-02, 4.9003d-02, 5.2065d-02, &
    &5.5024d-02, 5.7793d-02, 5.5900d-02, 5.4462d-02, 5.9976d-02, 6.4599d-02, 6.2582d-02, 6.3757d-02, 6.9081d-02, &
    &7.3676d-02, 7.9389d-02]
real(dp),         parameter :: KB_DNV_ERR(2:17)  = [1.2680d-01, 9.8836d-02, 9.3614d-02, 9.8526d-02, 9.7769d-02, &
    &1.0279d-01, 1.0630d-01, 9.7726d-02, 1.0088d-01, 1.1254d-01, 1.0753d-01, 1.0549d-01, 1.0019d-01, 1.0729d-01, &
    &1.1706d-01, 1.1582d-01]
real(dp),         parameter :: FLOOR_FACTOR = 2.d0
! the bank's recomputed errors reproduce Phase 0's (5 printed digits) when the setting is the same
real(dp),         parameter :: KB_REPRO_RTOL = 1.d-3
! the central-difference constants of simple_pftc_inplane_tester (T4 names them)
real(dp),         parameter :: FD_STEP        = 1.d-3
real(dp),         parameter :: GRAD_FD_RTOL   = 1.d-2
real(dp),         parameter :: GRAD_FD_AFLOOR = 3.d-3
real(dp),         parameter :: PROBE_SHIFTS(2,3) = reshape([0.d0,0.d0, 1.7d0,-2.3d0, -0.6d0,0.4d0], [2,3])
real(dp),         parameter :: PROBE_ROTOFFS(4)  = [0.37d0, -1.42d0, 0.5d0, 1.93d0]
! T4's truth pose: the drawn reference, the tilts (rad), a grid rotation index and a shift (px); the probes' tilt offset
real,             parameter :: REF_EULER(3)  = [47., 63., 0.]
real(dp),         parameter :: TRUTH_TILT(2) = [0.06d0, -0.04d0]
integer,          parameter :: TRUTH_IROT    = 11
real(dp),         parameter :: TRUTH_SH(2)   = [1.4d0, -0.6d0]
real(dp),         parameter :: PROBE_DTILT(2) = [0.025d0, -0.03d0]
! the particle is the model at the truth: the loss there is zero up to single-precision rounding
real(dp),         parameter :: ZERO_RESIDUAL_TOL = 1.d-5

! the pftc keeps a reference to the sigma2 array it is given: module lifetime
real(sp), allocatable, target :: sigma2_fixture(:,:)

contains

    subroutine run_all_polarft_cont_calc_tests()
        write(*,'(A)') '**** running all polar continuous view tests ****'
        call write_phantom()
        call test_section_accuracy()
        call test_view_gradient_and_curvature(.false.)
        call test_view_gradient_and_curvature(.true.)
        call test_threads()
        call del_file(PHANTOM_FILE)
    end subroutine run_all_polarft_cont_calc_tests

    ! ---- fixture ---------------------------------------------------------------

    subroutine write_phantom()
        type(image) :: vol
        real, allocatable :: rmat(:,:,:)
        real    :: ctr, dx, dy, dz
        integer :: ib, i, j, k
        allocate(rmat(BOX,BOX,BOX), source=0.)
        ctr = real(BOX)/2. + 1.
        do k = 1, BOX
            do j = 1, BOX
                do i = 1, BOX
                    do ib = 1, NBLOBS
                        dx = real(i) - ctr - CTRS(1,ib); dy = real(j) - ctr - CTRS(2,ib); dz = real(k) - ctr - CTRS(3,ib)
                        rmat(i,j,k) = rmat(i,j,k) + AMPS(ib) * exp(-(dx*dx + dy*dy + dz*dz) / (2. * SIGMAS(ib)**2))
                    enddo
                enddo
            enddo
        enddo
        call vol%new([BOX,BOX,BOX], SMPD)
        call vol%set_rmat(rmat, .false.)
        call vol%write(string(PHANTOM_FILE), del_if_exists=.true.)
        call vol%kill()
    end subroutine write_phantom

    !> the unmasked phantom as the KB bank of nspace spiral directions (even half) and as the B-spline references of one
    !! state (both halves), the polar grid of lp, one particle slot (odd: of the odd half, its sigma2 rising), unit sigma2
    subroutine build_fixture( b, p, pcont, lp, nspace, nthr, odd )
        class(builder),    target, intent(inout) :: b
        class(parameters), target, intent(inout) :: p
        class(polarft_cont_calc),  intent(inout) :: pcont
        real,                      intent(in)    :: lp
        integer,                   intent(in)    :: nspace
        integer, optional,         intent(in)    :: nthr
        logical, optional,         intent(in)    :: odd
        type(cmdline)     :: cline
        real, allocatable :: rmat(:,:,:)
        logical :: l_odd
        integer :: k
        l_odd = .false.
        if( present(odd) ) l_odd = odd
        call cline%set('vol1',      PHANTOM_FILE)
        call cline%set('mskdiam',   MSKDIAM)
        call cline%set('smpd',      SMPD)
        call cline%set('lp',        lp)
        call cline%set('trs',       TRS)
        call cline%set('nptcls',    1.)
        call cline%set('nspace',    nspace)
        call cline%set('ctf',       'no')
        call cline%set('objfun',    'euclid')
        call cline%set('mirr_proj', 'no')
        if( present(nthr) ) call cline%set('nthr', nthr)
        call cline%check
        call b%init_params_and_build_strategy3D_tbox(cline, p)
        call set_bp_range3D(p, b, cline)
        call b%pftc%new(p, p%nspace, [1,1], p%kfromto)
        call b%vol%read(p%vols(1))
        call b%vol_pad%new([p%box_croppd, p%box_croppd, p%box_croppd], p%smpd_crop, wthreads=.false.)
        call b%vol%pad_fft(b%vol_pad)
        call b%vol_pad%expand_cmat()
        call vol_pad2ref_pfts_opt(b%pftc, b%vol_pad, b%eulspace, 1, iseven=.true.)
        call b%pftc%set_eo(1, .not. l_odd)
        call b%pftc%memoize_refs
        rmat = b%vol%get_rmat()
        call b%cftc%new(1, p%box, 1)
        call b%cftc%set_ref(1, .true.,  rmat)
        call b%cftc%set_ref(1, .false., rmat)
        if( allocated(sigma2_fixture) ) deallocate(sigma2_fixture)
        allocate(sigma2_fixture(p%kfromto(1):p%kfromto(2),2), source=1.0_sp)
        ! the odd half's spectrum rises over the shells: the half's column and the shell weights are exercised
        if( l_odd ) sigma2_fixture(:,2) = [(1.0_sp + 0.05_sp*real(k - p%kfromto(1)), k = p%kfromto(1), p%kfromto(2))]
        call b%pftc%assign_sigma2_noise(sigma2_fixture)
        call pcont%new(b%pftc, b%cftc)
        call cline%kill
    end subroutine build_fixture

    subroutine kill_fixture( b, pcont )
        class(builder),           intent(inout) :: b
        class(polarft_cont_calc), intent(inout) :: pcont
        call pcont%kill
        call b%cftc%kill
        call b%vol_pad%kill_expanded
        call b%vol_pad%kill
        call b%kill_strategy3D_tbox
        call b%kill_general_tbox
        if( allocated(sigma2_fixture) ) deallocate(sigma2_fixture)
    end subroutine kill_fixture

    !> the analytic transform of the phantom and its derivative along n at q (native lattice; phase origin at the box
    !! centre, the scale N^-3 of a polar reference section)
    subroutine analytic( q, n, fa, fan )
        real(dp),    intent(in)  :: q(3), n(3)
        complex(dp), intent(out) :: fa, fan
        complex(dp) :: fb
        real(dp)    :: sig2, g, ph
        integer     :: ib
        fa  = (0.d0,0.d0)
        fan = (0.d0,0.d0)
        do ib = 1, NBLOBS
            sig2 = real(SIGMAS(ib),dp)**2
            g    = real(AMPS(ib),dp) * (DTWOPI*sig2)**1.5d0 * exp(-2.d0*DPI**2*sig2*sum(q*q)/real(BOX,dp)**2)
            ph   = -DTWOPI * sum(q*real(CTRS(:,ib),dp)) / real(BOX,dp)
            fb   = g * cmplx(cos(ph), sin(ph), kind=dp) / real(BOX,dp)**3
            fa   = fa  + fb
            fan  = fan + fb * cmplx(-4.d0*DPI**2*sig2*sum(q*n)/real(BOX,dp)**2, &
                &-DTWOPI*sum(real(CTRS(:,ib),dp)*n)/real(BOX,dp), kind=dp)
        enddo
    end subroutine analytic

    ! ---- T1, T2: the section and its normal derivative against the analytic transform ----------------------------

    subroutine test_section_accuracy()
        class(builder),    allocatable, target :: b
        class(parameters), allocatable, target :: p
        type(polarft_cont_calc)  :: pcont
        complex(sp), allocatable :: sect(:,:), dnv(:,:), bank(:,:)
        real(sp),    allocatable :: hc(:,:), kc(:,:)
        real(dp),    allocatable :: nb(:), nk(:), nn(:), den(:), denn(:)
        complex(dp) :: fa, fan
        real(dp)    :: frame(3,3), q(3), n(3), err_b, err_k, err_n, t_bspline, t_kb
        integer     :: pftsz, kfromto(2), iproj, k, irot, rep
        integer(int64) :: c0, c1, crate
        logical     :: ok1, ok2, okr
        integer, parameter :: NREP = 10
        write(*,'(A)') 'test_section_accuracy (T1, T2)'
        allocate(b, p)
        call build_fixture(b, p, pcont, LP_FULL, NSPACE_T12)
        pftsz   = b%pftc%get_pftsz()
        kfromto = b%pftc%get_kfromto()
        call assert_true(kfromto(1) <= KJUDGE(1) .and. kfromto(2) >= KJUDGE(2), 'the band holds the judged shells 2-17')
        allocate(sect(pftsz,kfromto(1):kfromto(2)), dnv(pftsz,kfromto(1):kfromto(2)), bank(pftsz,kfromto(1):kfromto(2)), &
            &hc(kfromto(1):kfromto(2),pftsz), kc(kfromto(1):kfromto(2),pftsz))
        allocate(nb(kfromto(1):kfromto(2)), nk(kfromto(1):kfromto(2)), nn(kfromto(1):kfromto(2)), &
            &den(kfromto(1):kfromto(2)), denn(kfromto(1):kfromto(2)), source=0.d0)
        call b%pftc%get_polar_coords(1, hc)
        call b%pftc%get_polar_coords(2, kc)
        do iproj = 1, NSPACE_T12
            frame = real(b%eulspace%get_mat(iproj), dp)
            n     = frame(3,:)
            call pcont%section(1, .true., frame, sect, dnv)
            call b%pftc%get_ref_pft(iproj, .true., bank)
            do k = kfromto(1), kfromto(2)
                do irot = 1, pftsz
                    q = real(hc(k,irot),dp)*frame(1,:) + real(kc(k,irot),dp)*frame(2,:)
                    call analytic(q, n, fa, fan)
                    nb(k)   = nb(k)   + abs(cmplx(sect(irot,k),kind=dp) - fa)**2
                    nk(k)   = nk(k)   + abs(cmplx(bank(irot,k),kind=dp) - fa)**2
                    nn(k)   = nn(k)   + abs(cmplx(dnv(irot,k), kind=dp) - fan)**2
                    den(k)  = den(k)  + abs(fa)**2
                    denn(k) = denn(k) + abs(fan)**2
                enddo
            enddo
        enddo
        ok1 = .true.; ok2 = .true.; okr = .true.
        write(*,'(A)') '  shell  KB section  B-spline section  T1 floor   KB dV/dn  B-spline dV/dn  T2 floor'
        do k = KJUDGE(1), KJUDGE(2)
            err_b = sqrt(nb(k)/den(k))
            err_k = sqrt(nk(k)/den(k))
            err_n = sqrt(nn(k)/denn(k))
            write(*,'(I7,6ES13.4)') k, err_k, err_b, FLOOR_FACTOR*KB_SECT_ERR(k), KB_DNV_ERR(k), err_n, &
                &FLOOR_FACTOR*KB_DNV_ERR(k)
            okr = okr .and. abs(err_k - KB_SECT_ERR(k)) <= KB_REPRO_RTOL*KB_SECT_ERR(k)
            ok1 = ok1 .and. err_b <= FLOOR_FACTOR*KB_SECT_ERR(k)
            ok2 = ok2 .and. err_n <= FLOOR_FACTOR*KB_DNV_ERR(k)
        enddo
        call assert_true(okr, 'T1 setting: the KB bank reproduces the errors frozen in Phase 0')
        call assert_true(ok1, 'T1: the B-spline section is within twice the KB error at every judged shell')
        call assert_true(ok2, 'T2: the B-spline dV/dn is within twice the KB error at every judged shell')
        ! cost per section: the B-spline section with its normal derivative against the KB bank (fproject over the spiral)
        call system_clock(c0, crate)
        do rep = 1, NREP
            do iproj = 1, NSPACE_T12
                frame = real(b%eulspace%get_mat(iproj), dp)
                call pcont%section(1, .true., frame, sect, dnv)
            enddo
        enddo
        call system_clock(c1)
        t_bspline = real(c1 - c0, dp)/real(crate, dp)/real(NREP*NSPACE_T12, dp)
        call system_clock(c0)
        do rep = 1, NREP
            call vol_pad2ref_pfts_opt(b%pftc, b%vol_pad, b%eulspace, 1, iseven=.true.)
        enddo
        call system_clock(c1)
        t_kb = real(c1 - c0, dp)/real(crate, dp)/real(NREP*NSPACE_T12, dp)
        write(*,'(A,I0,A,I0,A,F10.0,A,F10.0,A)') '  ns per section (', pftsz, ' x ', kfromto(2)-kfromto(1)+1, &
            &' samples): B-spline with dV/dn ', 1.d9*t_bspline, ', KB bank ', 1.d9*t_kb, ' (the team in use)'
        call kill_fixture(b, pcont)
        deallocate(b, p)
    end subroutine test_section_accuracy

    ! ---- T4: the five-variable gradient and Gauss-Newton matrix ----------------------------------------------------

    !> the particle is the model at the truth: the B-spline section at the truth tilts, shifted, rotated by a grid index
    subroutine set_truth_particle( b, pcont, rref )
        class(builder),           intent(inout) :: b
        class(polarft_cont_calc), intent(inout) :: pcont
        real(dp),                 intent(in)    :: rref(3,3)
        complex(sp), allocatable :: sect(:,:), dnv(:,:), rot(:,:), pft(:,:)
        real(dp),    allocatable :: ax(:,:), ay(:,:)
        integer :: pftsz, kfromto(2), pdim(3)
        pftsz   = b%pftc%get_pftsz()
        kfromto = b%pftc%get_kfromto()
        pdim    = b%pftc%get_pdim_interp()
        allocate(sect(pftsz,kfromto(1):kfromto(2)), dnv(pftsz,kfromto(1):kfromto(2)), rot(pftsz,kfromto(1):kfromto(2)), &
            &ax(pftsz,kfromto(1):kfromto(2)), ay(pftsz,kfromto(1):kfromto(2)), pft(pdim(1),pdim(2):pdim(3)))
        call pcont%section(1, .true., pcont%view_frame(rref, TRUTH_TILT), sect, dnv)
        call shift_args(pftsz, kfromto, BOX, ax, ay)
        sect = sect * cmplx(cos(ax*TRUTH_SH(1) + ay*TRUTH_SH(2)), sin(ax*TRUTH_SH(1) + ay*TRUTH_SH(2)), kind=sp)
        call rotate_half_ring(pftsz, kfromto(2)-kfromto(1)+1, sect, TRUTH_IROT, rot)
        pft = cmplx(0.,0.,kind=sp)
        pft(:,kfromto(1):kfromto(2)) = rot
        call b%pftc%set_ptcl_pft(1, pft)
        call b%pftc%memoize_ptcls
        call b%pftc%memoize_sqsum_ptcl(1)
    end subroutine set_truth_particle

    !> central differences of f in every variable at x: the largest error and the check's tolerance
    subroutine fd_gradient( pcont, x, g, max_error, tol )
        class(polarft_cont_calc), intent(inout) :: pcont
        real(dp),                 intent(in)    :: x(5), g(5)
        real(dp),                 intent(out)   :: max_error, tol
        real(dp) :: xp(5), fp, fm, gtmp(5)
        integer  :: d
        max_error = 0.d0
        do d = 1, 5
            xp = x; xp(d) = xp(d) + FD_STEP
            call pcont%eval(1, xp, fp, gtmp)
            xp = x; xp(d) = xp(d) - FD_STEP
            call pcont%eval(1, xp, fm, gtmp)
            max_error = max(max_error, abs((fp - fm)/(2.d0*FD_STEP) - g(d)))
        enddo
        tol = max(GRAD_FD_RTOL*maxval(abs(g)), GRAD_FD_AFLOOR)
    end subroutine fd_gradient

    !> the central-difference Hessian of the analytic gradient at x
    subroutine fd_hessian( pcont, x, h )
        class(polarft_cont_calc), intent(inout) :: pcont
        real(dp),                 intent(in)    :: x(5)
        real(dp),                 intent(out)   :: h(5,5)
        real(dp) :: xp(5), f, gp(5), gm(5)
        integer  :: d
        do d = 1, 5
            xp = x; xp(d) = xp(d) + FD_STEP
            call pcont%eval(1, xp, f, gp)
            xp = x; xp(d) = xp(d) - FD_STEP
            call pcont%eval(1, xp, f, gm)
            h(:,d) = (gp - gm)/(2.d0*FD_STEP)
        enddo
    end subroutine fd_hessian

    subroutine test_view_gradient_and_curvature( l_odd )
        use simple_linalg, only: cholesky
        logical, intent(in) :: l_odd
        class(builder),    allocatable, target :: b
        class(parameters), allocatable, target :: p
        type(polarft_cont_calc) :: pcont
        real(dp) :: rref(3,3), xt(5), x(5), f, g(5), gflip(5), c(5,5), h(5,5), cbad(5,5), l(5,5), eye(5,5)
        real(dp) :: err, tol, worst, err_flip, err_c, err_bad, err_blk, flip_ratio
        integer  :: ishift, ioff, i
        logical  :: all_ok, any_flip_fail, finite, ok
        write(*,'(A,L1)') 'test_view_gradient_and_curvature (T4), odd half with rising sigma2: ', l_odd
        allocate(b, p)
        call build_fixture(b, p, pcont, LP_LOW, 2, odd=l_odd)
        rref = real(euler2m(REF_EULER), dp)
        call set_truth_particle(b, pcont, rref)
        call pcont%set_view_ref(1, .not. l_odd, rref)
        xt = [TRUTH_SH(1), TRUTH_SH(2), real(TRUTH_IROT,dp), TRUTH_TILT(1), TRUTH_TILT(2)]
        call pcont%eval(1, xt, f, g)
        write(*,'(A,ES11.3,A,5ES11.3)') '  at the truth: loss ', f, ', gradient ', g
        call assert_true(abs(f) <= ZERO_RESIDUAL_TOL, 'T4: the particle is the model at the truth (zero residual)')
        ! the gradient at the twelve probe poses, tilts off the truth
        all_ok        = .true.
        any_flip_fail = .false.
        finite        = .true.
        worst         = 0.d0
        flip_ratio    = 0.d0
        do ishift = 1, size(PROBE_SHIFTS, 2)
            do ioff = 1, size(PROBE_ROTOFFS)
                x = [PROBE_SHIFTS(1,ishift), PROBE_SHIFTS(2,ishift), real(TRUTH_IROT,dp) + PROBE_ROTOFFS(ioff), &
                    &TRUTH_TILT(1) + PROBE_DTILT(1), TRUTH_TILT(2) + PROBE_DTILT(2)]
                call pcont%eval(1, x, f, g)
                finite = finite .and. ieee_is_finite(f) .and. all(ieee_is_finite(g))
                call fd_gradient(pcont, x, g, err, tol)
                all_ok = all_ok .and. err <= tol
                worst  = max(worst, err)
                ! negative control: dV/dn with its sign flipped negates both view components
                gflip      = g
                gflip(4:5) = -g(4:5)
                call fd_gradient(pcont, x, gflip, err_flip, tol)
                any_flip_fail = any_flip_fail .or. err_flip > tol
                flip_ratio    = max(flip_ratio, err_flip/tol)
            enddo
        enddo
        write(*,'(A,ES10.3,A,F8.1)') '  worst FD gradient error over the 12 probes: ', worst, &
            &'; sign-flipped dV/dn, largest error over tolerance ', flip_ratio
        call assert_true(finite, 'T4: the cost is finite at every probe')
        call assert_true(all_ok, 'T4: the analytic gradient matches central differences at 12 probe poses')
        call assert_true(any_flip_fail, 'T4 negative control: a sign-flipped dV/dn fails the gradient check')
        ! the Gauss-Newton matrix at zero residual is the central-difference Hessian of the analytic gradient
        call pcont%eval(1, xt, f, g, c)
        call fd_hessian(pcont, xt, h)
        err_c = maxval(abs(c - h))/maxval(abs(h))
        err_blk = maxval(abs(c(4:5,:) - h(4:5,:)))/maxval(abs(h(4:5,:)))
        write(*,'(A,ES10.3,A,ES10.3)') '  Gauss-Newton against the Hessian: max relative difference ', err_c, &
            &', view rows ', err_blk
        do i = 1, 5
            write(*,'(A,5ES11.3,A,5ES11.3)') '   C ', c(i,:), ' | H ', h(i,:)
        enddo
        call assert_true(err_c <= GRAD_FD_RTOL .and. err_blk <= GRAD_FD_RTOL, &
            &'T4: the Gauss-Newton matrix is the central-difference Hessian at zero residual')
        eye = 0.d0
        do i = 1, 5
            eye(i,i) = 1.d0
        enddo
        call cholesky(c + 1.d-12*maxval(abs(c))*eye, l, 5, ok)
        call assert_true(ok .and. all(c == transpose(c)), 'T4: the Gauss-Newton matrix is symmetric positive semi-definite')
        ! negative control: the view-shift coupling with its sign flipped fails the Hessian check
        cbad        = c
        cbad(4:5,1:3) = -c(4:5,1:3)
        cbad(1:3,4:5) = -c(1:3,4:5)
        err_bad = maxval(abs(cbad(4:5,:) - h(4:5,:)))/maxval(abs(h(4:5,:)))
        write(*,'(A,ES10.3)') '  flipped view coupling against the Hessian: view rows ', err_bad
        call assert_true(err_bad > GRAD_FD_RTOL, 'T4 negative control: a flipped view coupling fails the Hessian check')
        write(*,'(A,I0)') '  sections cut: ', pcont%get_nsections()
        call kill_fixture(b, pcont)
        deallocate(b, p)
    end subroutine test_view_gradient_and_curvature

    ! ---- threads ---------------------------------------------------------------------------------------------------

    !> a team of four evaluates the probe poses, each thread with its own drawn reference and section cache (reused: a
    !! thread's probes share one tilt), bit for bit as one thread does, which re-cuts at every probe
    subroutine test_threads()
        class(builder),    allocatable, target :: b
        class(parameters), allocatable, target :: p
        type(polarft_cont_calc) :: pcont
        integer,  parameter :: NPROBE = 12, NTEAM = 4
        real(dp) :: rref(3,3), x(5,NPROBE), fs(NPROBE), gs(5,NPROBE), cs(5,5,NPROBE), ft(NPROBE), gt(5,NPROBE), &
            &ct(5,5,NPROBE)
        integer  :: i, nthr_prev
        integer(int64) :: nsect_serial
        integer  :: team_size
        write(*,'(A)') 'test_threads'
        nthr_prev = 1
        !$ nthr_prev = omp_get_max_threads()
        allocate(b, p)
        call build_fixture(b, p, pcont, LP_LOW, 2, nthr=NTEAM)
        rref = real(euler2m(REF_EULER), dp)
        call set_truth_particle(b, pcont, rref)
        do i = 1, NPROBE
            x(:,i) = [PROBE_SHIFTS(1,mod(i,3)+1), PROBE_SHIFTS(2,mod(i,3)+1), real(TRUTH_IROT,dp) + PROBE_ROTOFFS(mod(i,4)+1), &
                &TRUTH_TILT(1) + 0.01d0*real(mod(i-1,NTEAM),dp), TRUTH_TILT(2) - 0.005d0*real(mod(i-1,NTEAM),dp)]
        enddo
        call pcont%set_view_ref(1, .true., rref)
        do i = 1, NPROBE
            call pcont%eval(1, x(:,i), fs(i), gs(:,i), cs(:,:,i))
        enddo
        nsect_serial = pcont%get_nsections()
        !$omp parallel num_threads(NTEAM) default(shared) proc_bind(close)
        call pcont%set_view_ref(1, .true., rref)
        !$omp single
        team_size = 1
        !$ team_size = omp_get_num_threads()
        !$omp end single
        !$omp do private(i) schedule(static,1)
        do i = 1, NPROBE
            call pcont%eval(1, x(:,i), ft(i), gt(:,i), ct(:,:,i))
        enddo
        !$omp end do
        !$omp end parallel
        write(*,'(A,I0,A,I0)') '  sections cut: serial ', nsect_serial, ', team ', pcont%get_nsections() - nsect_serial
        call assert_true(all(fs == ft) .and. all(gs == gt) .and. all(cs == ct), &
            &'a team of four evaluates the cost bit for bit as one thread')
        ! with a full team each thread's probes share one tilt (a thread limit below four mixes them)
        if( team_size == NTEAM )then
            call assert_true(pcont%get_nsections() - nsect_serial == int(NTEAM,int64), &
                &'each thread of the team cuts its one tilt once and reuses the section')
        endif
        call kill_fixture(b, pcont)
        deallocate(b, p)
        !$ call omp_set_num_threads(nthr_prev)
    end subroutine test_threads

end module simple_polarft_cont_calc_tester
