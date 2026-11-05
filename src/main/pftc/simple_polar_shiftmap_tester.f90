!@descr: unit tests of the polar shift map: the Fourier-Bessel mechanics against closed forms and the calculator's shifted, rotated objective, and the calculator's all-rotation kernel against the mechanics
! The mechanics must equal (1) the closed-form cosine sum of a pure phase ramp on the whole grid and off it,
! (2) gen_corr_for_rot_8 (the reference shifted in its own frame, then rotated, objfun=cc) at random grid poses
! including every branch of the rotation, and (3) locate a shifted, rotated copy at its pose, exact and at unit
! SNR. The kernel (gen_shiftmap) must equal the mechanics and the calculator's objective for cc and euclid at
! random poses, give the same grids on four threads as on one, and its cost is tabled against the coarse grid.
module simple_polar_shiftmap_tester
use simple_core_module_api
use simple_cmdline,         only: cmdline
use simple_parameters,      only: parameters
use simple_polarft_calc,    only: polarft_calc, SHIFTMAP_FULL, SHIFTMAP_SCREEN
use simple_polar_shiftmap,  only: polar_shiftmap, shiftmap_ndirections
use simple_timer,           only: tic, toc, timer_int_kind
use simple_rnd,             only: rnd_scope_key, RND_PROG_REFINE3D, RND_PASS_ALIGN
use, intrinsic :: iso_fortran_env, only: int64
use simple_test_utils
!$ use omp_lib,             only: omp_get_max_threads, omp_set_num_threads, omp_get_num_procs, omp_get_num_threads
implicit none
private
public :: run_all_polar_shiftmap_tests

contains

    subroutine run_all_polar_shiftmap_tests()
        write(*,'(A)') '**** running all polar shift map tests ****'
        call test_closed_form_ramp()
        call test_against_pftc_cc()
        call test_peak_recovery()
        call test_kernel_identity('cc')
        call test_kernel_identity('euclid')
        call test_orders_beyond_ring()
        call test_angular_nyquist()
        call test_cell_areas()
        call test_kernel_threads()
        call test_kernel_draw()
        call test_physical_width()
        call test_zero_prior_width_rounding()
        call test_map_cache_lifecycle()
        call test_screen_width_identity()
        call test_disk_within_radius()
        call test_kernel_timing()
    end subroutine run_all_polar_shiftmap_tests

    !> the shift arguments of the polar shift matrix, as the calculator forms them (argtransf):
    !! ax = k sin(theta_p) 2 pi / box, ay = -k cos(theta_p) 2 pi / box, theta_p = (p-1) 2 pi / nrots
    pure subroutine shift_args( pftsz, kfromto, box, ax, ay )
        integer,  intent(in)  :: pftsz, kfromto(2), box
        real(dp), intent(out) :: ax(pftsz,kfromto(1):kfromto(2)), ay(pftsz,kfromto(1):kfromto(2))
        real(dp) :: theta
        integer  :: p, k
        do p = 1, pftsz
            theta = real(p-1,dp) * DTWOPI / real(2*pftsz,dp)
            do k = kfromto(1), kfromto(2)
                ax(p,k) =  real(k,dp) * sin(theta) * DTWOPI / real(box,dp)
                ay(p,k) = -real(k,dp) * cos(theta) * DTWOPI / real(box,dp)
            enddo
        enddo
    end subroutine shift_args

    !> a half ring rotated by grid rotation irot, as rotate_ref_8 does it: entry p of the result is entry
    !! p - irot + 1 (modulo the turn) of the Friedel-extended ring, whose second half is the conjugate of the first
    pure subroutine rotate_half_ring( pftsz, nk, pft, irot, pft_rot )
        integer,     intent(in)  :: pftsz, nk, irot
        complex(sp), intent(in)  :: pft(pftsz,nk)
        complex(sp), intent(out) :: pft_rot(pftsz,nk)
        integer :: p, q
        do p = 1, pftsz
            q = modulo(p - irot, 2*pftsz) + 1
            if( q <= pftsz )then
                pft_rot(p,:) = pft(q,:)
            else
                pft_rot(p,:) = conjg(pft(q-pftsz,:))
            endif
        enddo
    end subroutine rotate_half_ring

    !> A pure phase ramp Q = conj(shift matrix of s0) has the closed-form cross term
    !! C(s) = sum_k w_k sum_p cos(ax (sx - s0x) + ay (sy - s0y)), largest at s = s0 where it is sum_k w_k pftsz.
    !! The map must reproduce it on every grid point and off the grid, peak at s0 when s0 is a grid point, be
    !! direction-independent at zero magnitude (order 0 only), and keep at least ceiling(z) orders at every magnitude.
    subroutine test_closed_form_ramp()
        integer,  parameter :: PFTSZ = 40, BOX = 64, NALPHA = 24, NRHO = 6
        integer,  parameter :: KFROMTO(2) = [1, 12]
        real(dp), parameter :: RHOS(NRHO) = [0.d0, 1.d0, 2.5d0, 5.d0, 8.d0, 12.d0]
        type(polar_shiftmap)     :: sm
        real(dp),    allocatable :: ax(:,:), ay(:,:), w(:), cmap(:,:)
        complex(dp), allocatable :: Q(:,:), A(:,:)
        real(dp) :: s0(2), s(2), expected, scale, maxerr, alpha, z, val
        integer  :: k, ir, j, irho, ialpha, nmax, itrial
        write(*,'(A)') 'test_closed_form_ramp'
        call set_fixed_seed(20261030)
        allocate(ax(PFTSZ,KFROMTO(1):KFROMTO(2)), ay(PFTSZ,KFROMTO(1):KFROMTO(2)), w(KFROMTO(1):KFROMTO(2)),&
                &Q(PFTSZ,KFROMTO(1):KFROMTO(2)))
        call shift_args(PFTSZ, KFROMTO, BOX, ax, ay)
        do k = KFROMTO(1), KFROMTO(2)
            w(k) = real(k,dp)
        enddo
        scale = sum(w) * real(PFTSZ,dp)
        call sm%new(PFTSZ, KFROMTO, BOX, RHOS, NALPHA)
        nmax = sm%get_nharm_max()
        allocate(A(-nmax:nmax,NRHO), cmap(NALPHA,NRHO))
        ! the truncation rule
        call assert_int(0, sm%get_nharm(1), 'zero magnitude keeps order 0 only')
        do ir = 2, NRHO
            z = DTWOPI * real(KFROMTO(2),dp) * RHOS(ir) / real(BOX,dp)
            call assert_true(sm%get_nharm(ir) >= ceiling(z), 'at least ceiling(z) orders are kept')
            call assert_true(sm%get_nharm(ir) >= sm%get_nharm(ir-1), 'the order kept grows with the magnitude')
        enddo
        write(*,'(A,6I5)') '  orders kept at the magnitudes: ', (sm%get_nharm(ir), ir=1,NRHO)
        ! s0 a grid point: the map on the grid, the peak, the zero-magnitude row
        s0 = sm%get_shift(4, 7)
        Q  = exp(cmplx(0.d0, -1.d0, kind=dp) * (ax * s0(1) + ay * s0(2)))
        call sm%coeffs(Q, w, A)
        call sm%map(A, cmap)
        maxerr = 0.d0
        do ir = 1, NRHO
            do j = 1, NALPHA
                s        = sm%get_shift(ir, j)
                expected = sum(spread(w, 1, PFTSZ) * cos(ax * (s(1) - s0(1)) + ay * (s(2) - s0(2))))
                maxerr   = max(maxerr, abs(cmap(j,ir) - expected))
            enddo
        enddo
        write(*,'(A,ES10.2,A,ES10.2)') '  grid: max |map - closed form| ', maxerr, ' at scale ', scale
        call assert_true(maxerr <= 1.d-9 * scale, 'the map equals the closed form on the grid')
        call sm%peak(cmap, irho, ialpha, val)
        call assert_int(4, irho,   'the peak is at the ramp''s magnitude')
        call assert_int(7, ialpha, 'the peak is at the ramp''s direction')
        call assert_true(abs(val - scale) <= 1.d-9 * scale, 'the peak value is the full weight')
        call assert_true(maxval(cmap(:,1)) - minval(cmap(:,1)) <= 1.d-12 * scale, 'zero magnitude is direction-independent')
        expected = sum(spread(w, 1, PFTSZ) * real(Q, dp))
        call assert_true(abs(cmap(1,1) - expected) <= 1.d-9 * scale, 'zero magnitude is the unshifted cross term')
        ! s0 off the grid, evaluated off the grid in direction
        s0 = [-3.3d0, 1.7d0]
        Q  = exp(cmplx(0.d0, -1.d0, kind=dp) * (ax * s0(1) + ay * s0(2)))
        call sm%coeffs(Q, w, A)
        maxerr = 0.d0
        do itrial = 1, 50
            ir    = 1 + int(ran3() * real(NRHO))
            ir    = min(max(ir, 1), NRHO)
            alpha = real(ran3(),dp) * DTWOPI
            s     = sm%get_rho(ir) * [cos(alpha), sin(alpha)]
            expected = sum(spread(w, 1, PFTSZ) * cos(ax * (s(1) - s0(1)) + ay * (s(2) - s0(2))))
            maxerr   = max(maxerr, abs(sm%eval(A, ir, alpha) - expected))
        enddo
        write(*,'(A,ES10.2)') '  off-grid directions: max |eval - closed form| ', maxerr
        call assert_true(maxerr <= 1.d-9 * scale, 'eval equals the closed form at any direction')
        call sm%kill
        deallocate(ax, ay, w, Q, A, cmap)
    end subroutine test_closed_form_ramp

    !> the calculator's fixture of these tests: random band-limited polar sections, objfun=cc, no CTF
    subroutine make_pftc( p, pftc, ref, ptcl, kfromto, pftsz, box )
        class(parameters), target, intent(inout) :: p
        class(polarft_calc),       intent(inout) :: pftc
        complex(sp), allocatable,  intent(out)   :: ref(:,:), ptcl(:,:)
        integer,                   intent(out)   :: kfromto(2), pftsz, box
        integer, parameter :: BOXSZ = 64, PFTSIZE = 48, KHI = 10
        type(cmdline) :: cline
        integer :: k, pp, pdim(3)
        call cline%set('box',     real(BOXSZ))
        call cline%set('smpd',    1.0)
        call cline%set('mskdiam', 48.0)
        call cline%set('pftsz',   real(PFTSIZE))
        call cline%set('nptcls',  1.0)
        call cline%set('nthr',    1.0)
        call cline%set('ctf',     'no')
        call cline%set('objfun',  'cc')
        call p%new(cline)
        kfromto = [2, KHI]
        call pftc%new(p, 1, [1,1], kfromto)
        call pftc%set_rnd_scope(rnd_scope_key(20261105_int64, RND_PROG_REFINE3D, RND_PASS_ALIGN, 1))
        pftsz = pftc%get_pftsz()
        box   = BOXSZ
        call assert_int(PFTSIZE, pftsz, 'the calculator takes the requested pftsz')
        ! the calculator's section layout (pftsz, kfromto(1):interpklim); only the band is nonzero
        pdim = pftc%get_pdim_interp()
        allocate(ref(pdim(1),pdim(2):pdim(3)), ptcl(pdim(1),pdim(2):pdim(3)), source=cmplx(0.,0.))
        do k = kfromto(1), kfromto(2)
            do pp = 1, pftsz
                ref(pp,k)  = cmplx(gasdev(), gasdev())
                ptcl(pp,k) = cmplx(gasdev(), gasdev())
            enddo
        enddo
    end subroutine make_pftc

    !> the k-weighted self term sum_k k sum_p |pft|^2 of a band (pftsz, kfromto(1):kfromto(2)); rotation and
    !! shift leave it unchanged
    real(dp) function kweighted_sqsum( pft, kfromto )
        integer,     intent(in) :: kfromto(2)
        complex(sp), intent(in) :: pft(:,kfromto(1):)
        integer :: k
        kweighted_sqsum = 0.d0
        do k = kfromto(1), kfromto(2)
            kweighted_sqsum = kweighted_sqsum + real(k,dp) * sum(real(pft(:,k)*conjg(pft(:,k)),dp))
        enddo
    end function kweighted_sqsum

    !> The map's cross term, normalised by the two k-weighted self terms, is gen_corr_for_rot_8 at the same grid
    !! shift and rotation (objfun=cc, no CTF): at the four rotations that exercise every branch of rotate_ref_8
    !! and at random rotations, magnitudes and directions. The calculator builds its shift arguments from
    !! single-precision angles, so the agreement is to 1e-6 (the double-precision exactness is the closed form's).
    subroutine test_against_pftc_cc()
        integer,  parameter :: NALPHA = 20, NRHO = 5, NTRIALS = 40
        real(dp), parameter :: RHOS(NRHO) = [0.d0, 1.5d0, 4.d0, 9.d0, 15.d0]
        class(parameters),   allocatable, target :: p
        class(polarft_calc), allocatable         :: pftc
        type(polar_shiftmap)       :: sm
        complex(sp), allocatable   :: ref(:,:), ptcl(:,:)
        real(sp),    allocatable   :: ctf(:,:)
        real(dp),    allocatable   :: w(:)
        complex(dp), allocatable   :: Q(:,:), A(:,:)
        real(dp) :: s(2), norm, cc_ref, cc_map, maxerr
        integer  :: kfromto(2), pftsz, box, nrots, nmax, k, itrial, irot, ir, j
        integer  :: branch_rots(4)
        write(*,'(A)') 'test_against_pftc_cc'
        call set_fixed_seed(20261031)
        allocate(p, pftc)
        call make_pftc(p, pftc, ref, ptcl, kfromto, pftsz, box)
        nrots = pftc%get_nrots()
        call pftc%set_ref_pft(1, ref, iseven=.true.)
        call pftc%set_ptcl_pft(1, ptcl)
        call pftc%set_eo(1, .true.)
        call sm%new(pftsz, kfromto, box, RHOS, NALPHA)
        nmax = sm%get_nharm_max()
        allocate(w(kfromto(1):kfromto(2)), ctf(pftsz,kfromto(1):kfromto(2)), Q(pftsz,kfromto(1):kfromto(2)),&
                &A(-nmax:nmax,NRHO))
        do k = kfromto(1), kfromto(2)
            w(k) = real(k,dp)
        enddo
        ctf  = 1.0
        norm = sqrt(kweighted_sqsum(ref, kfromto) * kweighted_sqsum(ptcl, kfromto))
        branch_rots = [1, pftsz, pftsz + 1, nrots]
        maxerr = 0.d0
        do itrial = 1, NTRIALS + 4
            if( itrial <= 4 )then
                irot = branch_rots(itrial)
                ir   = 1 + mod(itrial, NRHO)
                j    = 1 + mod(3 * itrial, NALPHA)
            else
                irot = min(nrots,  1 + int(ran3() * real(nrots)))
                ir   = min(NRHO,   1 + int(ran3() * real(NRHO)))
                j    = min(NALPHA, 1 + int(ran3() * real(NALPHA)))
            endif
            s      = sm%get_shift(ir, j)
            cc_ref = pftc%gen_corr_for_rot_8(1, 1, s, irot)
            call sm%cross_spectrum(ref(:,kfromto(1):kfromto(2)), ptcl(:,kfromto(1):kfromto(2)), ctf, irot, Q)
            call sm%coeffs(Q, w, A)
            cc_map = sm%eval(A, ir, sm%get_alpha(j)) / norm
            maxerr = max(maxerr, abs(cc_map - cc_ref))
        enddo
        write(*,'(A,ES10.2)') '  max |map cc - gen_corr_for_rot_8| ', maxerr
        call assert_true(maxerr <= 1.d-6, 'the map is the calculator''s shifted, rotated cc at every probe')
        call sm%kill
        call pftc%kill
        deallocate(ref, ptcl, w, ctf, Q, A, p, pftc)
    end subroutine test_against_pftc_cc

    !> A particle that is the reference shifted by a grid shift and rotated by a grid rotation: the calculator's cc
    !! at that pose is 1 (the construction follows its conventions), and the map's largest value over every
    !! rotation and grid point is at that pose, exact and with white noise of the signal's power added.
    subroutine test_peak_recovery()
        integer,  parameter :: NALPHA = 24, NRHO = 9
        real(dp), parameter :: RHOS(NRHO) = [0.d0, 1.d0, 2.d0, 3.d0, 4.d0, 6.d0, 8.d0, 11.d0, 14.d0]
        integer,  parameter :: IRHO_TRUE = 6, IALPHA_TRUE = 17, IROT_TRUE = 61
        class(parameters),   allocatable, target :: p
        class(polarft_calc), allocatable         :: pftc
        type(polar_shiftmap)       :: sm
        complex(sp), allocatable   :: ref(:,:), ptcl(:,:), shifted(:,:)
        real(sp),    allocatable   :: ctf(:,:)
        real(dp),    allocatable   :: w(:), ax(:,:), ay(:,:), cmap(:,:)
        complex(dp), allocatable   :: Q(:,:), A(:,:)
        real(dp) :: s_true(2), cc_true, best, val, sig
        integer  :: kfromto(2), pftsz, box, nrots, nmax, k, pp, irot, irho, ialpha, best_pose(3), ipass
        write(*,'(A)') 'test_peak_recovery'
        call set_fixed_seed(20261101)
        allocate(p, pftc)
        call make_pftc(p, pftc, ref, ptcl, kfromto, pftsz, box)
        nrots = pftc%get_nrots()
        call sm%new(pftsz, kfromto, box, RHOS, NALPHA)
        nmax = sm%get_nharm_max()
        allocate(w(kfromto(1):kfromto(2)), ctf(pftsz,kfromto(1):kfromto(2)), Q(pftsz,kfromto(1):kfromto(2)),&
                &A(-nmax:nmax,NRHO), cmap(NALPHA,NRHO), ax(pftsz,kfromto(1):kfromto(2)), ay(pftsz,kfromto(1):kfromto(2)),&
                &shifted(pftsz,kfromto(1):kfromto(2)))
        do k = kfromto(1), kfromto(2)
            w(k) = real(k,dp)
        enddo
        ctf = 1.0
        call shift_args(pftsz, kfromto, box, ax, ay)
        s_true = sm%get_shift(IRHO_TRUE, IALPHA_TRUE)
        ! the particle: the reference shifted in its frame by s_true, then rotated by IROT_TRUE
        shifted = cmplx(cmplx(ref(:,kfromto(1):kfromto(2)),kind=dp) * exp(cmplx(0.d0,1.d0,kind=dp) * (ax*s_true(1) + ay*s_true(2))))
        ptcl = cmplx(0.,0.)
        call rotate_half_ring(pftsz, kfromto(2)-kfromto(1)+1, shifted, IROT_TRUE, ptcl(:,kfromto(1):kfromto(2)))
        call pftc%set_ref_pft(1, ref, iseven=.true.)
        call pftc%set_eo(1, .true.)
        sig = sqrt(sum(real(ref(:,kfromto(1):kfromto(2))*conjg(ref(:,kfromto(1):kfromto(2))),dp)) &
            &/ real(pftsz*(kfromto(2)-kfromto(1)+1),dp))
        do ipass = 1, 2
            if( ipass == 2 )then
                ! white noise of the signal's power on the particle
                do k = kfromto(1), kfromto(2)
                    do pp = 1, pftsz
                        ptcl(pp,k) = ptcl(pp,k) + cmplx(real(sig)*gasdev()/sqrt(2.), real(sig)*gasdev()/sqrt(2.))
                    enddo
                enddo
            endif
            call pftc%set_ptcl_pft(1, ptcl)
            cc_true = pftc%gen_corr_for_rot_8(1, 1, s_true, IROT_TRUE)
            if( ipass == 1 )then
                call assert_true(cc_true > 0.999d0, 'the construction reproduces the calculator''s pose conventions (cc 1 at the true pose)')
            else
                write(*,'(A,F8.4)') '  cc at the true pose at unit SNR ', cc_true
            endif
            best      = -huge(best)
            best_pose = 0
            do irot = 1, nrots
                call sm%cross_spectrum(ref(:,kfromto(1):kfromto(2)), ptcl(:,kfromto(1):kfromto(2)), ctf, irot, Q)
                call sm%coeffs(Q, w, A)
                call sm%map(A, cmap)
                call sm%peak(cmap, irho, ialpha, val)
                if( val > best )then
                    best      = val
                    best_pose = [irot, irho, ialpha]
                endif
            enddo
            if( ipass == 1 )then
                call assert_int(IROT_TRUE,   best_pose(1), 'exact copy: the map peaks at the true rotation')
                call assert_int(IRHO_TRUE,   best_pose(2), 'exact copy: the map peaks at the true magnitude')
                call assert_int(IALPHA_TRUE, best_pose(3), 'exact copy: the map peaks at the true direction')
            else
                call assert_int(IROT_TRUE,   best_pose(1), 'unit SNR: the map peaks at the true rotation')
                call assert_int(IRHO_TRUE,   best_pose(2), 'unit SNR: the map peaks at the true magnitude')
                call assert_int(IALPHA_TRUE, best_pose(3), 'unit SNR: the map peaks at the true direction')
            endif
        enddo
        call sm%kill
        call pftc%kill
        deallocate(ref, ptcl, shifted, w, ctf, Q, A, cmap, ax, ay, p, pftc)
    end subroutine test_peak_recovery

    !> The kernel's fixture: nrefs random band-limited references with distinct halves, nptcls particles alternating
    !! odd and even (particle 1 odd), each with its own angularly varying |CTF| (set_ptcl_ctfmat), the memoised
    !! spectra; for euclid two different shell-dependent sigma2 spectra, kept by the caller while the calculator
    !! lives. pftsz_in 0 lets the parameters derive pftsz from the mask
    subroutine make_pftc_kernel( p, pftc, objfun, nthr, box, smpd, mskdiam, lp, pftsz_in, nrefs, nptcls, sigma2, refs,&
            &ptcls, ctfs, kfromto )
        class(parameters), target,        intent(inout) :: p
        class(polarft_calc),              intent(inout) :: pftc
        character(len=*),                 intent(in)    :: objfun
        integer,                          intent(in)    :: nthr, box, pftsz_in, nrefs, nptcls
        real,                             intent(in)    :: smpd, mskdiam, lp
        real,        allocatable, target, intent(inout) :: sigma2(:,:)
        complex(sp), allocatable,         intent(out)   :: refs(:,:,:,:), ptcls(:,:,:)
        real(sp),    allocatable,         intent(out)   :: ctfs(:,:,:)
        integer,                          intent(out)   :: kfromto(2)
        type(cmdline) :: cline
        real    :: theta
        integer :: k, pp, pdim(3), iref, iptcl, pftsz, ieo
        call cline%set('box',     real(box))
        call cline%set('smpd',    smpd)
        call cline%set('mskdiam', mskdiam)
        if( pftsz_in > 0 ) call cline%set('pftsz', real(pftsz_in))
        call cline%set('nptcls',  real(nptcls))
        call cline%set('nthr',    real(nthr))
        call cline%set('ctf',     'yes')
        call cline%set('objfun',  objfun)
        call p%new(cline)
        kfromto = [2, min(box/2 - 1, max(3, nint(real(box) * smpd / lp)))]
        call pftc%new(p, nrefs, [1,nptcls], kfromto)
        call pftc%set_rnd_scope(rnd_scope_key(20261105_int64, RND_PROG_REFINE3D, RND_PASS_ALIGN, 1))
        pftsz = pftc%get_pftsz()
        pdim  = pftc%get_pdim_interp()
        allocate(refs(pdim(1),pdim(2):pdim(3),nrefs,2), ptcls(pdim(1),pdim(2):pdim(3),nptcls), source=cmplx(0.,0.))
        allocate(ctfs(pdim(1),pdim(2):pdim(3),nptcls), source=1.0)
        do iref = 1, nrefs
            do ieo = 1, 2
                do k = kfromto(1), kfromto(2)
                    do pp = 1, pftsz
                        refs(pp,k,iref,ieo) = cmplx(gasdev(), gasdev())
                    enddo
                enddo
                call pftc%set_ref_pft(iref, refs(:,:,iref,ieo), iseven=(ieo == 1))
            enddo
        enddo
        do iptcl = 1, nptcls
            do k = kfromto(1), kfromto(2)
                do pp = 1, pftsz
                    theta = real(pp-1) * TWOPI / real(2*pftsz)
                    ptcls(pp,k,iptcl) = cmplx(gasdev(), gasdev())
                    ctfs(pp,k,iptcl)  = 0.3 + 0.7 * abs(sin(0.15 * real(k) * (1. + 0.1 * real(iptcl)) + 0.8 * cos(2. * theta + 0.3)))
                enddo
            enddo
            call pftc%set_eo(iptcl, mod(iptcl,2) == 0)
            call pftc%set_ptcl_ctfmat(iptcl, ctfs(:,:,iptcl))
            call pftc%set_ptcl_pft(iptcl, ptcls(:,:,iptcl))
        enddo
        call pftc%memoize_refs
        call pftc%memoize_ptcls
        if( trim(objfun) == 'euclid' )then
            if( allocated(sigma2) ) deallocate(sigma2)
            allocate(sigma2(kfromto(1):kfromto(2),2))
            do k = kfromto(1), kfromto(2)
                sigma2(k,1) = 0.5 + 0.1 * real(k)
                sigma2(k,2) = 0.8 + 0.05 * real(k)
            enddo
            call pftc%assign_sigma2_noise(sigma2)
        endif
    end subroutine make_pftc_kernel

    !> The kernel against the mechanics class and the calculator, particle 1 (odd, its own CTF) against reference 2's
    !! odd half: at random grid poses, including the four rotations that exercise every branch of the rotation, the
    !! cross term equals the mechanics class's with the objective's weights, and the grid's objective equals
    !! gen_corr_for_rot_8 (euclid returns exp(-L)); the zero-magnitude row is the unshifted objective at every
    !! rotation; a grid with fewer directions than orders is exact on its directions (the folding). The kernel's
    !! spectra and transforms are single precision, so the agreement is 1e-6 of the Cauchy-Schwarz scale of the cross
    !! term (the mechanics are exact in double precision); for euclid the error is also given in the sampler's
    !! units, T |dL| at the particle's temperature.
    subroutine test_kernel_identity( objfun )
        use simple_eul_prob_tab_utils, only: posterior_temperature
        character(len=*), intent(in) :: objfun
        integer,  parameter :: NALPHA = 20, NRHO = 5, NTRIALS = 60, NALPHA_FOLD = 7, KSCREEN = 6
        real(dp), parameter :: RHOS(NRHO) = [0.d0, 1.5d0, 4.d0, 9.d0, 15.d0]
        class(parameters),   allocatable, target :: p
        class(polarft_calc), allocatable         :: pftc
        type(polar_shiftmap)      :: sm
        real,        allocatable  :: sigma2(:,:)
        complex(sp), allocatable  :: refs(:,:,:,:), ptcls(:,:,:)
        real(sp),    allocatable  :: ctfs(:,:,:), vals(:,:,:), vals_before(:,:,:)
        real(dp),    allocatable  :: w(:)
        complex(dp), allocatable  :: Q(:,:), A(:,:)
        real(dp) :: scale, alpha, c_mech, c_kern, err_mech, err_grid, err_zero, err_fold, val_ref, val_map, s(2), temp
        real(dp) :: err_nats
        integer  :: kfromto(2), pftsz, nrots, nmax, k, itrial, irot, ir, j, branch_rots(4), nrho_c, nalpha_c
        write(*,'(A)') 'test_kernel_identity ('//trim(objfun)//')'
        call set_fixed_seed(20261102)
        allocate(p, pftc)
        call make_pftc_kernel(p, pftc, objfun, 1, 64, 1.0, 48.0, 6.0, 48, 2, 1, sigma2, refs, ptcls, ctfs, kfromto)
        pftsz = pftc%get_pftsz()
        nrots = pftc%get_nrots()
        call pftc%new_shiftmap(SHIFTMAP_FULL, RHOS, NALPHA, kfromto(2))
        call assert_true(pftc%has_shiftmap(SHIFTMAP_FULL), 'the calculator holds a shift map')
        call pftc%shiftmap_dims(SHIFTMAP_FULL, nrho_c, nalpha_c)
        call assert_int(NRHO,   nrho_c,   'the calculator holds the grid magnitudes')
        call assert_int(NALPHA, nalpha_c, 'the calculator holds the grid directions')
        call sm%new(pftsz, kfromto, 64, RHOS, NALPHA)
        nmax = sm%get_nharm_max()
        allocate(w(kfromto(1):kfromto(2)), Q(pftsz,kfromto(1):kfromto(2)), A(-nmax:nmax,NRHO), vals(NALPHA,nrots,NRHO))
        do k = kfromto(1), kfromto(2)
            if( trim(objfun) == 'cc' )then
                w(k) = real(k,dp)
            else
                w(k) = real(k,dp) / real(sigma2(k,2),dp)
            endif
        enddo
        ! the Cauchy-Schwarz scale of the cross term with the CTF on the particle side
        scale = 0.d0
        do k = kfromto(1), kfromto(2)
            scale = scale + w(k) * sqrt(sum(real(refs(:,k,2,2)*conjg(refs(:,k,2,2)),dp)) * &
                &sum(real(ctfs(:,k,1)**2 * ptcls(:,k,1)*conjg(ptcls(:,k,1)),dp)))
        enddo
        temp = 0.d0
        if( trim(objfun) == 'euclid' ) temp = real(posterior_temperature(pftc%wsqsum_ptcl_now(1), pftsz), dp)
        call pftc%gen_shiftmap(SHIFTMAP_FULL, 2, 1)
        call pftc%shiftmap_objfun_grid(SHIFTMAP_FULL, vals)
        if( trim(objfun) == 'cc' )then
            allocate(vals_before,source=vals)
            call pftc%memoize_sigma2_weights
            call pftc%gen_shiftmap(SHIFTMAP_FULL, 2, 1)
            call pftc%shiftmap_objfun_grid(SHIFTMAP_FULL, vals)
            call assert_true(all(vals == vals_before), 'CC map weights survive a noise refresh without sigma2')
            deallocate(vals_before)
        endif
        branch_rots = [1, pftsz, pftsz + 1, nrots]
        err_mech = 0.d0
        err_grid = 0.d0
        err_nats = 0.d0
        do itrial = 1, NTRIALS + 4
            if( itrial <= 4 )then
                irot = branch_rots(itrial)
                ir   = 1 + mod(itrial, NRHO)
                j    = 1 + mod(3 * itrial, NALPHA)
            else
                irot = min(nrots,  1 + int(ran3() * real(nrots)))
                ir   = min(NRHO,   1 + int(ran3() * real(NRHO)))
                j    = min(NALPHA, 1 + int(ran3() * real(NALPHA)))
            endif
            alpha = real(ran3(),dp) * DTWOPI
            call sm%cross_spectrum(refs(:,kfromto(1):kfromto(2),2,2), ptcls(:,kfromto(1):kfromto(2),1),&
                &ctfs(:,kfromto(1):kfromto(2),1), irot, Q)
            call sm%coeffs(Q, w, A)
            c_mech   = sm%eval(A, ir, alpha)
            c_kern   = pftc%shiftmap_cross(SHIFTMAP_FULL, irot, ir, alpha)
            err_mech = max(err_mech, abs(c_kern - c_mech))
            s       = pftc%shiftmap_shift(SHIFTMAP_FULL, ir, j)
            val_ref = pftc%gen_corr_for_rot_8(2, 1, s, irot)
            if( trim(objfun) == 'cc' )then
                val_map = real(vals(j,irot,ir), dp)
            else
                val_map  = exp(-real(vals(j,irot,ir), dp))
                err_nats = max(err_nats, temp * abs(real(vals(j,irot,ir), dp) + log(val_ref)))
            endif
            err_grid = max(err_grid, abs(val_map - val_ref))
        enddo
        write(*,'(A,ES10.2,A,ES10.2)') '  max |kernel - mechanics| ', err_mech, ' at scale ', scale
        write(*,'(A,ES10.2)')          '  max |grid objective - gen_corr_for_rot_8| ', err_grid
        if( trim(objfun) == 'euclid' ) write(*,'(A,ES10.2,A)') '  max T |dL| ', err_nats, ' nats at the particle''s temperature'
        call assert_true(err_mech <= 1.d-6 * scale, 'the kernel''s cross term is the mechanics class''s at random poses')
        call assert_true(err_grid <= 1.d-6, 'the grid''s objective is the calculator''s shifted, rotated objective')
        err_zero = 0.d0
        do irot = 1, nrots
            val_ref = pftc%gen_corr_for_rot_8(2, 1, [0.d0, 0.d0], irot)
            if( trim(objfun) == 'cc' )then
                val_map = real(vals(1,irot,1), dp)
            else
                val_map = exp(-real(vals(1,irot,1), dp))
            endif
            err_zero = max(err_zero, abs(val_map - val_ref))
        enddo
        write(*,'(A,ES10.2)') '  max |zero-magnitude row - unshifted objective| ', err_zero
        call assert_true(err_zero <= 1.d-6, 'the zero-magnitude row is the unshifted objective at every rotation')
        ! fewer directions than orders: the folded transform is exact on its own directions (the grid against the
        ! calculator's shifted, rotated objective at the grid's shifts, as above; cc)
        if( trim(objfun) == 'cc' )then
            call pftc%new_shiftmap(SHIFTMAP_FULL, RHOS, NALPHA_FOLD, kfromto(2))
            deallocate(vals)
            allocate(vals(NALPHA_FOLD,nrots,NRHO))
            call pftc%gen_shiftmap(SHIFTMAP_FULL, 2, 1)
            call pftc%shiftmap_objfun_grid(SHIFTMAP_FULL, vals)
            err_fold = 0.d0
            do itrial = 1, 24
                irot    = min(nrots, 1 + int(ran3() * real(nrots)))
                ir      = NRHO - mod(itrial, 2)
                j       = 1 + mod(itrial, NALPHA_FOLD)
                s       = pftc%shiftmap_shift(SHIFTMAP_FULL, ir, j)
                val_ref = pftc%gen_corr_for_rot_8(2, 1, s, irot)
                err_fold = max(err_fold, abs(real(vals(j,irot,ir), dp) - val_ref))
            enddo
            write(*,'(A,I2,A,ES10.2)') '  ', NALPHA_FOLD, ' directions, orders folded: max |grid - gen_corr_for_rot_8| ', err_fold
            call assert_true(err_fold <= 1.d-6, 'a grid with fewer directions than orders is exact on its directions')
            ! the screen: a second map over the lower half of the band is the mechanics over that band
            call pftc%new_shiftmap(SHIFTMAP_SCREEN, RHOS, NALPHA, KSCREEN)
            call sm%kill
            call sm%new(pftsz, [kfromto(1), KSCREEN], 64, RHOS, NALPHA)
            nmax = sm%get_nharm_max()
            deallocate(A, Q)
            allocate(A(-nmax:nmax,NRHO), Q(pftsz,kfromto(1):KSCREEN))
            call pftc%gen_shiftmap(SHIFTMAP_SCREEN, 2, 1)
            err_fold = 0.d0
            do itrial = 1, 24
                irot  = min(nrots, 1 + int(ran3() * real(nrots)))
                ir    = min(NRHO,  1 + int(ran3() * real(NRHO)))
                alpha = real(ran3(),dp) * DTWOPI
                call sm%cross_spectrum(refs(:,kfromto(1):KSCREEN,2,2), ptcls(:,kfromto(1):KSCREEN,1),&
                    &ctfs(:,kfromto(1):KSCREEN,1), irot, Q)
                call sm%coeffs(Q, w(kfromto(1):KSCREEN), A)
                err_fold = max(err_fold, abs(pftc%shiftmap_cross(SHIFTMAP_SCREEN, irot, ir, alpha) - sm%eval(A, ir, alpha)))
            enddo
            write(*,'(A,I2,A,ES10.2)') '  the screen over shells 2-', KSCREEN, ': max |kernel - mechanics| ', err_fold
            call assert_true(err_fold <= 1.d-6 * scale, 'the screen''s map is the mechanics over its band')
            call assert_true(pftc%has_shiftmap(SHIFTMAP_FULL) .and. pftc%has_shiftmap(SHIFTMAP_SCREEN), 'both maps exist')
        endif
        call sm%kill
        call pftc%kill
        deallocate(refs, ptcls, ctfs, w, Q, A, vals, p, pftc)
        if( allocated(sigma2) ) deallocate(sigma2)
    end subroutine test_kernel_identity

    !> the polar measure: the cells of a grid sum to the disk out to the outer edge, and the origin's cells to the
    !! central disk; an unsorted magnitude list gives the same areas
    !> The review's M3: a disk wide against a short ring needs Bessel orders beyond the ring's length, whose spectral
    !! index wraps modulo the ring; the grid's objective is still the calculator's shifted, rotated objective (cc, a
    !! ring of 2 x 16 samples, magnitudes out to 20 px at shell 21)
    subroutine test_orders_beyond_ring()
        integer,  parameter :: NALPHA = 24, NRHO = 4, NTRIALS = 40
        real(dp), parameter :: RHOS(NRHO) = [0.d0, 6.d0, 13.d0, 20.d0]
        class(parameters),   allocatable, target :: p
        class(polarft_calc), allocatable         :: pftc
        type(polar_shiftmap)      :: sm
        real,        allocatable  :: sigma2(:,:)
        complex(sp), allocatable  :: refs(:,:,:,:), ptcls(:,:,:)
        real(sp),    allocatable  :: ctfs(:,:,:), vals(:,:,:)
        real(dp) :: s(2), err, val_ref
        integer  :: kfromto(2), pftsz, nrots, itrial, irot, ir, j
        write(*,'(A)') 'test_orders_beyond_ring'
        call set_fixed_seed(20261104)
        allocate(p, pftc)
        call make_pftc_kernel(p, pftc, 'cc', 1, 64, 1.0, 48.0, 3.0, 16, 2, 1, sigma2, refs, ptcls, ctfs, kfromto)
        pftsz = pftc%get_pftsz()
        nrots = pftc%get_nrots()
        call sm%new(pftsz, kfromto, 64, RHOS, NALPHA)
        call assert_true(sm%get_nharm_max() > nrots, 'the fixture needs orders beyond the ring''s length')
        call pftc%new_shiftmap(SHIFTMAP_FULL, RHOS, NALPHA, kfromto(2))
        allocate(vals(NALPHA,nrots,NRHO))
        call pftc%gen_shiftmap(SHIFTMAP_FULL, 2, 1)
        call pftc%shiftmap_objfun_grid(SHIFTMAP_FULL, vals)
        err = 0.d0
        do itrial = 1, NTRIALS
            irot    = min(nrots,  1 + int(ran3() * real(nrots)))
            ir      = min(NRHO,   1 + int(ran3() * real(NRHO)))
            j       = min(NALPHA, 1 + int(ran3() * real(NALPHA)))
            s       = pftc%shiftmap_shift(SHIFTMAP_FULL, ir, j)
            val_ref = pftc%gen_corr_for_rot_8(2, 1, s, irot)
            err     = max(err, abs(real(vals(j,irot,ir), dp) - val_ref))
        enddo
        write(*,'(A,I3,A,I3,A,ES10.2)') '  orders to ', sm%get_nharm_max(), ' on a ring of ', nrots,&
            &': max |grid - gen_corr_for_rot_8| ', err
        call assert_true(err <= 1.d-6, 'orders beyond the ring wrap: the grid is the calculator''s objective')
        call sm%kill
        call pftc%kill
        deallocate(refs, ptcls, ctfs, vals, p, pftc)
        if( allocated(sigma2) ) deallocate(sigma2)
    end subroutine test_orders_beyond_ring

    !> The review's M1: every Euclidean objective carries the angular Nyquist bin. On the fixture's random spectra the
    !! score row (gen_objfun_vals, exp(-L)) is the direct residual (gen_corr_for_rot_8) at every rotation, unshifted and
    !! shifted, and the continuous evaluator is the row at the grid angles; with reference and particle the alternating
    !! half ring (-1)^(p-1), unit CTF, the loss is 0 at rotation 1 and 4 at rotation 2 in all four, the map included
    subroutine test_angular_nyquist()
        integer,  parameter :: NALPHA = 8, NRHO = 2
        real(dp), parameter :: RHOS(NRHO) = [0.d0, 2.d0], TOL = 1.d-5
        class(parameters),   allocatable, target :: p
        class(polarft_calc), allocatable         :: pftc
        real,        allocatable  :: sigma2(:,:)
        complex(sp), allocatable  :: refs(:,:,:,:), ptcls(:,:,:), alt(:,:)
        real(sp),    allocatable  :: ctfs(:,:,:), row(:), vals(:,:,:), ones(:,:)
        real(dp) :: f, grad(3), err_row, err_cont, shifts(2,2), l_dir, l_row, l_map, l_cont
        integer  :: kfromto(2), pftsz, nrots, k, pp, irot, is, pdim(3)
        write(*,'(A)') 'test_angular_nyquist'
        call set_fixed_seed(20261105)
        allocate(p, pftc)
        call make_pftc_kernel(p, pftc, 'euclid', 1, 64, 1.0, 48.0, 6.0, 36, 2, 1, sigma2, refs, ptcls, ctfs, kfromto)
        pftsz  = pftc%get_pftsz()
        nrots  = pftc%get_nrots()
        pdim   = pftc%get_pdim_interp()
        shifts = reshape([0.d0, 0.d0, 1.3d0, -0.7d0], [2,2])
        allocate(row(nrots))
        err_row  = 0.d0
        err_cont = 0.d0
        do is = 1, 2
            call pftc%gen_objfun_vals(2, 1, real(shifts(:,is)), row)
            do irot = 1, nrots
                err_row = max(err_row, abs(real(row(irot),dp) - pftc%gen_corr_for_rot_8(2, 1, shifts(:,is), irot)))
                if( is > 1 ) cycle
                call pftc%gen_raw_euclid_grad_at_angle(2, 1, [0.d0, 0.d0], real(irot,dp), f, grad)
                err_cont = max(err_cont, abs(exp(-f) - real(row(irot),dp)))
            enddo
        enddo
        write(*,'(A,ES10.2,A,ES10.2)') '  max |row - direct| ', err_row, ', max |continuous - row| ', err_cont
        call assert_true(err_row  <= TOL, 'the score row is the direct residual at every rotation, unshifted and shifted')
        call assert_true(err_cont <= TOL, 'the continuous evaluator is the score row at the grid angles')
        ! the alternating half ring: all its energy at the angular Nyquist bin
        allocate(alt(pdim(1),pdim(2):pdim(3)), source=cmplx(0.,0.))
        allocate(ones(pdim(1),pdim(2):pdim(3)), source=1.)
        do k = kfromto(1), kfromto(2)
            do pp = 1, pftsz
                alt(pp,k) = cmplx(real(1 - 2 * mod(pp-1, 2)), 0.)
            enddo
        enddo
        call pftc%set_ref_pft(2, alt, iseven=.true.)
        call pftc%set_ref_pft(2, alt, iseven=.false.)
        call pftc%set_ptcl_ctfmat(1, ones)
        call pftc%set_ptcl_pft(1, alt)
        call pftc%memoize_refs
        call pftc%memoize_ptcls
        call pftc%assign_sigma2_noise(sigma2)
        call pftc%new_shiftmap(SHIFTMAP_FULL, RHOS, NALPHA, kfromto(2))
        allocate(vals(NALPHA,nrots,NRHO))
        call pftc%gen_shiftmap(SHIFTMAP_FULL, 2, 1)
        call pftc%shiftmap_objfun_grid(SHIFTMAP_FULL, vals)
        call pftc%gen_objfun_vals(2, 1, [0., 0.], row)
        do irot = 1, 2
            l_dir  = -log(pftc%gen_corr_for_rot_8(2, 1, [0.d0, 0.d0], irot))
            l_row  = -log(real(row(irot),dp))
            l_map  = real(vals(1,irot,1),dp)
            call pftc%gen_raw_euclid_grad_at_angle(2, 1, [0.d0, 0.d0], real(irot,dp), l_cont, grad)
            write(*,'(A,I1,A,4F9.5)') '  rotation ', irot, ': direct, row, map, continuous ', l_dir, l_row, l_map, l_cont
            call assert_true(abs(l_dir - real(4*(irot-1),dp)) <= 1.d-4, 'the direct residual of the Nyquist pattern')
            call assert_true(abs(l_row - l_dir) <= 1.d-4 .and. abs(l_map - l_dir) <= 1.d-4 .and. abs(l_cont - l_dir) <= 1.d-4,&
                &'the row, the map and the continuous evaluator agree with it on the Nyquist bin')
        enddo
        call pftc%kill
        deallocate(refs, ptcls, ctfs, row, alt, ones, vals, p, pftc)
        if( allocated(sigma2) ) deallocate(sigma2)
    end subroutine test_angular_nyquist

    subroutine test_cell_areas()
        integer,  parameter :: NALPHA = 12, NRHO = 5
        real(dp), parameter :: RHOS(NRHO) = [0.d0, 2.d0, 4.d0, 6.d0, 8.d0], SHUF(NRHO) = [6.d0, 0.d0, 8.d0, 2.d0, 4.d0]
        type(polar_shiftmap) :: sm, sm2
        real(dp) :: total, origin, expected
        integer  :: ir, j
        write(*,'(A)') 'test_cell_areas'
        call sm%new(40, [1, 10], 64, RHOS, NALPHA)
        total  = 0.d0
        do ir = 1, NRHO
            total = total + real(NALPHA,dp) * sm%get_cell_area(ir)
        enddo
        expected = DPI * 9.d0**2
        call assert_true(abs(total - expected) <= 1.d-9 * expected, 'the cells sum to the disk out to the outer edge')
        origin = real(NALPHA,dp) * sm%get_cell_area(1)
        call assert_true(abs(origin - DPI) <= 1.d-9 * DPI, 'the origin''s cells sum to the central disk of half a step')
        call assert_true(sm%get_cell_area(3) > sm%get_cell_area(2), 'the cell area grows with the magnitude')
        call sm2%new(40, [1, 10], 64, SHUF, NALPHA)
        do ir = 1, NRHO
            do j = 1, NRHO
                if( abs(SHUF(ir) - RHOS(j)) > 1.d-12 ) cycle
                call assert_true(abs(sm2%get_cell_area(ir) - sm%get_cell_area(j)) <= 1.d-12, 'the areas do not depend on the order of the magnitudes')
            enddo
        enddo
        call sm%kill
        call sm2%kill
        call assert_int(48, shiftmap_ndirections(47), 'the FFT-friendly direction count rounds a prime up')
        call assert_int(35, shiftmap_ndirections(33), 'the FFT-friendly direction count keeps smooth numbers close')
        call assert_int(14, shiftmap_ndirections(13), 'the FFT-friendly direction count admits a factor 7')
    end subroutine test_cell_areas

    !> A forced team must give the serial grids bit for bit with per-thread work and one shared plan.
    subroutine test_kernel_threads()
        integer,  parameter :: NALPHA = 16, NRHO = 4
        real(dp), parameter :: RHOS(NRHO) = [0.d0, 2.d0, 5.d0, 9.d0]
        class(parameters),   allocatable, target :: p
        class(polarft_calc), allocatable         :: pftc
        real,        allocatable :: sigma2(:,:)
        complex(sp), allocatable :: refs(:,:,:,:), ptcls(:,:,:)
        real(sp),    allocatable :: ctfs(:,:,:), vals_serial(:,:,:,:), vals_par(:,:,:,:)
        real(sp), allocatable :: width_serial(:), width_par(:), loss_serial(:), loss_par(:)
        integer, allocatable :: cells_serial(:,:), cells_par(:,:)
        integer :: kfromto(2), nrots, ip, nthr_prev, nthr_test, nptcls, actual_team, expected_team
        write(*,'(A)') 'test_kernel_threads'
        call set_fixed_seed(20261103)
        allocate(p, pftc)
        nthr_test = 4
        !$ nthr_test = min(24, max(4, omp_get_num_procs()))
        nptcls = 2 * nthr_test
        call make_pftc_kernel(p, pftc, 'euclid', nthr_test, 64, 1.0, 48.0, 6.0, 48, 1, nptcls, sigma2, refs, ptcls, ctfs, kfromto)
        nrots = pftc%get_nrots()
        call pftc%new_shiftmap(SHIFTMAP_FULL, RHOS, NALPHA, kfromto(2))
        allocate(vals_serial(NALPHA,nrots,NRHO,nptcls), vals_par(NALPHA,nrots,NRHO,nptcls))
        allocate(width_serial(nptcls), width_par(nptcls), loss_serial(nptcls), loss_par(nptcls), &
            &cells_serial(3,nptcls), cells_par(3,nptcls))
        nthr_prev = 1
        !$ nthr_prev = omp_get_max_threads()
        !$ call omp_set_num_threads(1)
        do ip = 1, nptcls
            call pftc%gen_shiftmap(SHIFTMAP_FULL, 1, ip)
            call pftc%shiftmap_objfun_grid(SHIFTMAP_FULL, vals_serial(:,:,:,ip))
            call pftc%shiftmap_mode(SHIFTMAP_FULL, 3.d0, 0.05d0, &
                &cells_serial(1,ip), cells_serial(2,ip), cells_serial(3,ip), loss_serial(ip), &
                &temp_phys=10.d0, sdev=width_serial(ip), &
                &width_cutoff=merge(huge(1.),-huge(1.),mod(ip,2)==0))
        enddo
        !$ call omp_set_num_threads(nthr_test)
        actual_team = 1
        expected_team = 1
        !$ expected_team = 4
        !$omp parallel default(shared) private(ip) proc_bind(close)
        !$omp single
        !$ actual_team = omp_get_num_threads()
        !$omp end single
        !$omp do schedule(dynamic)
        do ip = 1, nptcls
            call pftc%gen_shiftmap(SHIFTMAP_FULL, 1, ip)
            call pftc%shiftmap_objfun_grid(SHIFTMAP_FULL, vals_par(:,:,:,ip))
            call pftc%shiftmap_mode(SHIFTMAP_FULL, 3.d0, 0.05d0, &
                &cells_par(1,ip), cells_par(2,ip), cells_par(3,ip), loss_par(ip), &
                &temp_phys=10.d0, sdev=width_par(ip), &
                &width_cutoff=merge(huge(1.),-huge(1.),mod(ip,2)==0))
        enddo
        !$omp end do
        !$omp end parallel
        !$ call omp_set_num_threads(nthr_prev)
        write(*,'(A,I0)') '  forced map team: ', actual_team
        call assert_true(actual_team >= expected_team, 'the map comparison actually ran a parallel team')
        write(*,'(A,ES10.2)') '  max |threaded - serial| ', real(maxval(abs(vals_par - vals_serial)),dp)
        call assert_true(all(vals_par == vals_serial), 'the forced team gives the serial grids bit for bit')
        call assert_true(all(cells_par == cells_serial), 'the forced team gives the serial posterior modes')
        call assert_true(all(loss_par == loss_serial) .and. all(width_par == width_serial), &
            &'per-thread map work preserves losses and physical widths bit for bit')
        call assert_int(nptcls/2, count(width_par < 0.), 'parallel width cutoffs skip exactly the requested particles')
        call pftc%kill
        deallocate(refs, ptcls, ctfs, vals_serial, vals_par, p, pftc)
        if( allocated(sigma2) ) deallocate(sigma2)
    end subroutine test_kernel_threads

    !> The sampler on the kernel's fixture (euclid, box 64 at 1 A/px, a 30 px disk for sigma 10 px). At temperature
    !! zero the draw is the prior under the polar measure: the mean radius of 1500 draws is the discrete prior's
    !! expectation, which sits by the continuous truncated Gaussian's (12.30 for sigma 10 and a 30 disk, the review's
    !! figure) within the grid's step, and the evidence is 0; the origin's prior mass is the central disk's share (the
    !! prior is on the cell's shift: it moves with the particle); at the particle's temperature the
    !! evidence equals the direct log-sum-exp over shiftmap_objfun_grid's losses with the same weights, and the greedy
    !! draw is the cell of the largest mass.
    !> the disk never reaches beyond its radius (review R4): below one grid step the grid is the origin and the
    !! rim at the radius, at and above a step the outermost magnitude is within the radius
    subroutine test_disk_within_radius()
        class(parameters),   allocatable :: p
        class(polarft_calc), allocatable :: pftc
        real,        allocatable :: sigma2(:,:)
        complex(sp), allocatable :: refs(:,:,:,:), ptcls(:,:,:)
        real(sp),    allocatable :: ctfs(:,:,:)
        real(dp) :: step, prior_w, radius, x(2), rmax
        integer  :: kfromto(2), nrho, nalpha, kscreen, icase, ir
        write(*,'(A)') 'test_disk_within_radius'
        allocate(p, pftc)
        call make_pftc_kernel(p, pftc, 'euclid', 1, 64, 1.0, 48.0, 6.0, 48, 1, 1, sigma2, refs, ptcls, ctfs, kfromto)
        ! the fills' disk map lives on the coarse band (kscreen = max(k1 + 1, k2 / 2)); its radial step
        step = 64.d0 / real(max(kfromto(1) + 1, kfromto(2) / 2),dp) / 4.d0
        do icase = 1, 3
            select case(icase)
                case(1); radius = 0.5d0 * step
                case(2); radius = step
                case(3); radius = 1.5d0 * step
            end select
            call pftc%setup_shiftmaps(10.d0, radius, radius, prior_w, kscreen)
            call pftc%shiftmap_dims(SHIFTMAP_SCREEN, nrho, nalpha)
            rmax = 0.d0
            do ir = 1, nrho
                x    = pftc%shiftmap_shift(SHIFTMAP_SCREEN, ir, 1)
                rmax = max(rmax, sqrt(sum(x**2)))
            enddo
            call assert_true(rmax <= radius + 1.d-9, 'the outermost magnitude is within the requested radius')
            call assert_int(2, nrho, 'two magnitudes: the origin and one ring')
            if( icase == 1 ) call assert_true(abs(rmax - radius) <= 1.d-9, 'below one step the ring is the rim at the radius')
            call pftc%kill_shiftmap
        enddo
        call pftc%kill
        deallocate(p, pftc, sigma2, refs, ptcls, ctfs)
    end subroutine test_disk_within_radius

    subroutine test_kernel_draw()
        use simple_eul_prob_tab_utils, only: posterior_temperature
        integer,  parameter :: NDRAWS = 1500
        real(dp), parameter :: SIGMA = 10.d0, RADIUS = 30.d0, EXPECTED_MEAN = 12.30d0
        class(parameters),   allocatable, target :: p
        class(polarft_calc), allocatable         :: pftc
        real,        allocatable :: sigma2(:,:)
        complex(sp), allocatable :: refs(:,:,:,:), ptcls(:,:,:)
        real(sp),    allocatable :: ctfs(:,:,:), vals(:,:,:)
        real(dp),    allocatable :: rhos(:)
        real(dp) :: prior_w, lse_prior, logp0, logz, step, mean_rho, disc_mean, wsum, wr, w
        real(dp) :: lmax, s, lse, pen, loga, rho, temp
        real(sp) :: loss
        integer  :: kfromto(2), nrots, nrho, nalpha, idraw, irot, irho, ialpha, ir, j, ro, imax(3)
        write(*,'(A)') 'test_kernel_draw'
        call set_fixed_seed(20261105)
        allocate(p, pftc)
        call make_pftc_kernel(p, pftc, 'euclid', 1, 64, 1.0, 48.0, 6.0, 48, 1, 1, sigma2, refs, ptcls, ctfs, kfromto)
        nrots   = pftc%get_nrots()
        step    = 64.d0 / real(kfromto(2),dp) / 4.d0
        call magnitudes(RADIUS, step, rhos)
        nalpha  = shiftmap_ndirections(2 * ceiling(DTWOPI * real(kfromto(2),dp) * rhos(size(rhos)) / 64.d0) + 3)
        call pftc%new_shiftmap(SHIFTMAP_FULL, rhos, nalpha, kfromto(2))
        call pftc%shiftmap_dims(SHIFTMAP_FULL, nrho, nalpha)
        prior_w = 1.d0 / (2.d0 * SIGMA**2)
        call pftc%gen_shiftmap(SHIFTMAP_FULL, 1, 1)
        ! (1) temperature 0: the prior under the polar measure
        call pftc%shiftmap_prior(SHIFTMAP_FULL, prior_w, lse_prior, logp0)
        disc_mean = 0.d0
        wsum      = 0.d0
        do ir = 1, nrho
            w  = exp(-prior_w * rhos(ir)**2) * pftc%shiftmap_cell_area(SHIFTMAP_FULL, ir) * real(nalpha,dp)
            wr = w * rhos(ir)
            wsum      = wsum + w
            disc_mean = disc_mean + wr
        enddo
        disc_mean = disc_mean / wsum
        mean_rho  = 0.d0
        do idraw = 1, NDRAWS
            call pftc%shiftmap_draw(SHIFTMAP_FULL, 0.d0, prior_w, lse_prior, .false., logz, irot, irho, ialpha, loss,&
                &ordinal=idraw)
            mean_rho = mean_rho + rhos(irho)
            if( idraw == 1 ) call assert_true(abs(logz) <= 1.d-6, 'at temperature 0 the evidence is 0: the prior is normalised')
        enddo
        mean_rho = mean_rho / real(NDRAWS,dp)
        write(*,'(A,F7.3,A,F7.3,A,F7.3)') '  mean radius of the draws ', mean_rho, ', of the discrete prior ', disc_mean,&
            &', continuous ', EXPECTED_MEAN
        call assert_true(abs(mean_rho - disc_mean) <= 0.06d0 * disc_mean, 'at temperature 0 the draws follow the prior under the polar measure')
        call assert_true(abs(disc_mean - EXPECTED_MEAN) <= step, 'the discrete prior''s mean radius is the truncated Gaussian''s within a step')
        ! the origin's mass: the central disk's share of the prior, at one rotation
        w = exp(logp0) * real(nrots,dp)
        call assert_true(abs(w - exp(-prior_w * rhos(1)**2) * pftc%shiftmap_cell_area(SHIFTMAP_FULL, 1) * real(nalpha,dp) / wsum) <= 1.d-9,&
            &'the origin''s prior mass is the central disk''s share')
        ! (2) the particle's temperature: the evidence against the direct sum over the grid, the greedy cell
        call pftc%shiftmap_prior(SHIFTMAP_FULL, prior_w, lse_prior, logp0)
        temp = real(posterior_temperature(pftc%wsqsum_ptcl_now(1), pftc%get_pftsz()), dp)
        allocate(vals(nalpha,nrots,nrho))
        call pftc%shiftmap_objfun_grid(SHIFTMAP_FULL, vals)
        lmax = -huge(lmax)
        imax = 1
        do ir = 1, nrho
            rho  = rhos(ir)
            loga = log(pftc%shiftmap_cell_area(SHIFTMAP_FULL, ir))
            pen  = prior_w * rho * rho
            do ro = 1, nrots
                do j = 1, nalpha
                    s = -temp * real(vals(j,ro,ir),dp) - pen + loga
                    if( s > lmax )then
                        lmax = s
                        imax = [j, ro, ir]
                    endif
                enddo
            enddo
        enddo
        s = 0.d0
        do ir = 1, nrho
            pen  = prior_w * rhos(ir)**2
            loga = log(pftc%shiftmap_cell_area(SHIFTMAP_FULL, ir))
            do ro = 1, nrots
                do j = 1, nalpha
                    s = s + exp(-temp * real(vals(j,ro,ir),dp) - pen + loga - lmax)
                enddo
            enddo
        enddo
        lse = lmax + log(s) - lse_prior
        call pftc%shiftmap_draw(SHIFTMAP_FULL, temp, prior_w, lse_prior, .true., logz, irot, irho, ialpha, loss)
        write(*,'(A,ES10.2,A,F10.4)') '  |evidence - direct sum| ', abs(logz - lse), ' at log Z ', lse
        call assert_true(abs(logz - lse) <= 1.d-3 * max(1.d0, abs(lse)), 'the evidence is the log-sum-exp over the grid')
        call assert_int(imax(1), ialpha, 'the greedy draw is the cell of the largest mass (direction)')
        call assert_int(imax(2), irot,   'the greedy draw is the cell of the largest mass (rotation)')
        call assert_int(imax(3), irho,   'the greedy draw is the cell of the largest mass (magnitude)')
        call assert_true(abs(real(loss,dp) - real(vals(ialpha,irot,irho),dp)) <= 1.d-6, 'the returned loss is the cell''s')
        call pftc%kill
        deallocate(refs, ptcls, ctfs, vals, rhos, p, pftc)
        if( allocated(sigma2) ) deallocate(sigma2)
    end subroutine test_kernel_draw

    !> Physical moments from a direct grid sum, including broad and concentrated posteriors in both halves.
    subroutine test_physical_width()
        real(dp), parameter :: RHOS(4) = [0.d0, 1.5d0, 3.d0, 4.5d0]
        class(parameters),   allocatable, target :: p
        class(polarft_calc), allocatable :: pftc
        real, allocatable :: sigma2(:,:)
        complex(sp), allocatable :: refs(:,:,:,:), ptcls(:,:,:)
        real(sp), allocatable :: ctfs(:,:,:), vals(:,:,:)
        real(dp) :: prior_w, lse_prior, logp0, logz, temp, tphys, expected, radial_mass, radial_second
        real(sp) :: loss, width, mode_loss, mode_width
        integer :: kb(2), nrho, na, ip, icase, irot, irho, ialpha, ir, ieo, mr, mh, ma
        write(*,'(A)') 'test_physical_width'
        call set_fixed_seed(20261103)
        allocate(p, pftc)
        call make_pftc_kernel(p, pftc, 'euclid', 1, 64, 1.0, 48.0, 6.0, 48, 2, 2, sigma2, refs, ptcls, ctfs, kb)
        call pftc%new_shiftmap(SHIFTMAP_FULL, RHOS, 16, kb(2))
        call pftc%shiftmap_dims(SHIFTMAP_FULL, nrho, na)
        allocate(vals(na,pftc%get_nrots(),nrho))
        do icase = 1, 3
            if( icase == 2 )then
                call pftc%set_ref_pft(1, 0. * refs(:,:,1,1), iseven=.true.)
                call pftc%set_ref_pft(1, 0. * refs(:,:,1,2), iseven=.false.)
                call pftc%memoize_refs
            else if( icase == 3 )then
                do ip = 1, 2
                    ieo = merge(1, 2, mod(ip,2) == 0)
                    call pftc%set_ref_pft(1, refs(:,:,1,ieo), iseven=(ieo == 1))
                    call pftc%set_ptcl_pft(ip, refs(:,:,1,ieo))
                    call pftc%set_ptcl_ctfmat(ip, 0. * ctfs(:,:,ip) + 1.)
                enddo
                call pftc%memoize_refs
                call pftc%memoize_ptcls
            endif
            prior_w = 0.d0
            if( icase == 1 ) prior_w = 0.08d0
            tphys = 3.d0
            if( icase == 3 ) tphys = 300.d0
            temp = 0.3d0 * tphys
            call pftc%shiftmap_prior(SHIFTMAP_FULL, prior_w, lse_prior, logp0)
            do ip = 1, 2
                call pftc%gen_shiftmap(SHIFTMAP_FULL, 1, ip)
                call pftc%shiftmap_objfun_grid(SHIFTMAP_FULL, vals)
                expected = direct_shift_width(pftc, RHOS, vals, tphys, prior_w)
                call pftc%shiftmap_draw(SHIFTMAP_FULL, temp, prior_w, lse_prior, .true., &
                    &logz, irot, irho, ialpha, loss, temp_phys=tphys, sdev=width)
                call assert_true(abs(real(width,dp) - expected) <= 5.d-5 * max(1.d0, expected), &
                    &'physical posterior width agrees with independent grid moments')
                call pftc%shiftmap_mode(SHIFTMAP_FULL, temp, prior_w, mr, mh, ma, mode_loss, &
                    &temp_phys=tphys, sdev=mode_width)
                call assert_true(all([mr,mh,ma] == [irot,irho,ialpha]), 'mode-only selection keeps the evidence route cell')
                call assert_true(mode_loss == loss .and. mode_width == width, 'mode-only loss and width are bit identical')
                call pftc%shiftmap_mode(SHIFTMAP_FULL, temp, prior_w, mr, mh, ma, mode_loss, &
                    &temp_phys=tphys, sdev=mode_width, width_cutoff=loss)
                call assert_true(mode_width == width, 'a candidate tied at the cutoff retains its exact width')
                call pftc%shiftmap_mode(SHIFTMAP_FULL, temp, prior_w, mr, mh, ma, mode_loss, &
                    &temp_phys=tphys, sdev=mode_width, width_cutoff=-huge(1.))
                call assert_true(mode_width == -1. .and. mode_loss == loss, 'a rejected width leaves the selected loss unchanged')
                write(*,'(A,2I3,3I5,2ES24.16)') '  width identity ', icase, ip, irot, irho, ialpha, loss, width
                if( icase == 2 )then
                    radial_mass = 0.d0
                    radial_second = 0.d0
                    do ir = 1, nrho
                        radial_mass = radial_mass + pftc%shiftmap_cell_area(SHIFTMAP_FULL, ir)
                        radial_second = radial_second + pftc%shiftmap_cell_area(SHIFTMAP_FULL, ir) * RHOS(ir)**2
                    enddo
                    call assert_true(abs(real(width,dp) - sqrt(radial_second/radial_mass)) <= 1.d-5, &
                        &'flat rotational posterior retains the full radial uncertainty')
                else if( icase == 3 )then
                    call assert_true(width < 0.1, 'a concentrated aligned posterior has subpixel width')
                endif
            enddo
        enddo
        call pftc%kill
        deallocate(p, pftc, sigma2, refs, ptcls, ctfs, vals)
    end subroutine test_physical_width

    !> A sub-resolution prior leaves these represented masses unchanged. An off-origin, concentrated
    !! posterior exposes rounding changes hidden by a broad posterior or an exactly centred peak.
    subroutine test_zero_prior_width_rounding()
        real(dp), parameter :: RHOS(4) = [0.d0, 1.5d0, 3.d0, 4.5d0]
        real(dp), parameter :: TEMPS(5) = [10.d0, 20.d0, 30.d0, 40.d0, 80.d0]
        class(parameters), allocatable, target :: p
        class(polarft_calc), allocatable :: pftc
        real, allocatable :: sigma2(:,:)
        complex(sp), allocatable :: refs(:,:,:,:), ptcls(:,:,:)
        real(sp), allocatable :: ctfs(:,:,:)
        real(dp), allocatable :: ax(:,:), ay(:,:), phase(:,:)
        real(dp) :: xy(2), temp, tphys
        real(sp) :: loss, width, tiny_loss, tiny_width
        integer :: kb(2), ip, it, ieo, ro, ir, a, tr, ti, ta
        write(*,'(A)') 'test_zero_prior_width_rounding'
        call set_fixed_seed(20261103)
        allocate(p, pftc)
        call make_pftc_kernel(p, pftc, 'euclid', 1, 64, 1.0, 48.0, 6.0, 48, 2, 2, sigma2, refs, ptcls, ctfs, kb)
        allocate(ax(48,lbound(refs,2):ubound(refs,2)), ay(48,lbound(refs,2):ubound(refs,2)), &
            &phase(48,lbound(refs,2):ubound(refs,2)))
        call shift_args(48, [lbound(refs,2),ubound(refs,2)], 64, ax, ay)
        xy = 3.d0 * [cos(DTWOPI/16.d0), sin(DTWOPI/16.d0)]
        phase = ax * xy(1) + ay * xy(2)
        do ip = 1, 2
            ieo = merge(1, 2, ip == 2)
            ptcls(:,:,ip) = refs(:,:,1,ieo) * cmplx(cos(phase), sin(phase), kind=sp)
            call pftc%set_ptcl_pft(ip, ptcls(:,:,ip))
            call pftc%set_ptcl_ctfmat(ip, 0. * ctfs(:,:,ip) + 1.)
        enddo
        call pftc%memoize_ptcls
        call pftc%new_shiftmap(SHIFTMAP_FULL, RHOS, 16, kb(2))
        do ip = 1, 2
            call pftc%gen_shiftmap(SHIFTMAP_FULL, 1, ip)
            do it = 1, size(TEMPS)
                tphys = TEMPS(it)
                temp = 0.3d0 * tphys
                call pftc%shiftmap_mode(SHIFTMAP_FULL, temp, 0.d0, ro, ir, a, loss, temp_phys=tphys, sdev=width)
                call pftc%shiftmap_mode(SHIFTMAP_FULL, temp, 1.d-300, tr, ti, ta, tiny_loss, &
                    &temp_phys=tphys, sdev=tiny_width)
                call assert_true(all([ro,ir,a] == [tr,ti,ta]) .and. loss == tiny_loss, &
                    &'a sub-resolution prior leaves the selected cell and loss unchanged')
                call assert_true(width == tiny_width, 'zero and sub-resolution priors retain identical width rounding')
            enddo
        enddo
        call pftc%kill
        deallocate(p, pftc, sigma2, refs, ptcls, ctfs, ax, ay, phase)
    end subroutine test_zero_prior_width_rounding

    !> Direct posterior moments in particle axes; independent of the sampler's log-mass work and reductions.
    real(dp) function direct_shift_width( pftc, rhos, vals, temp, prior_w ) result( width )
        class(polarft_calc), intent(in) :: pftc
        real(dp), intent(in) :: rhos(:), temp, prior_w
        real(sp), intent(in) :: vals(:,:,:)
        real(dp), allocatable :: masses(:,:,:), points(:,:,:,:)
        real(dp) :: angle, a, x(2), u(2), lmax, z, mu(2), variance
        integer :: ir, ro, j, na, nr
        na = size(vals,1)
        nr = size(vals,2)
        allocate(masses(na,nr,size(rhos)), points(2,na,nr,size(rhos)))
        do ir = 1, size(rhos)
            do ro = 1, nr
                angle = real(pftc%get_rot(ro),dp) * DPI / 180.d0
                do j = 1, na
                    a = DTWOPI * real(j-1,dp) / real(na,dp)
                    x = rhos(ir) * [cos(a), sin(a)]
                    u = [cos(angle)*x(1) - sin(angle)*x(2), sin(angle)*x(1) + cos(angle)*x(2)]
                    points(:,j,ro,ir) = u
                    masses(j,ro,ir) = -temp * real(vals(j,ro,ir),dp) - prior_w * sum(x**2) &
                        &+ log(pftc%shiftmap_cell_area(SHIFTMAP_FULL,ir))
                enddo
            enddo
        enddo
        lmax = maxval(masses)
        masses = exp(masses - lmax)
        z = sum(masses)
        mu(1) = sum(masses * points(1,:,:,:)) / z
        mu(2) = sum(masses * points(2,:,:,:)) / z
        variance = sum(masses * ((points(1,:,:,:) - mu(1))**2 + (points(2,:,:,:) - mu(2))**2)) / z
        width = sqrt(max(0.d0, variance))
    end function direct_shift_width

    !> Cache reuse across noise updates, same-sized particle remaps and short batches against direct objectives.
    subroutine test_map_cache_lifecycle()
        real(dp), parameter :: RHOS(3) = [0.d0, 1.d0, 2.d0]
        class(parameters), allocatable, target :: p
        class(polarft_calc), allocatable :: pftc
        real, allocatable :: sigma2(:,:)
        complex(sp), allocatable :: refs(:,:,:,:), ptcls(:,:,:)
        real(sp), allocatable :: ctfs(:,:,:), before(:,:,:), after(:,:,:)
        real(dp) :: times(6)
        integer(int64) :: counts(8)
        integer :: kb(2), ip, cycle
        write(*,'(A)') 'test_map_cache_lifecycle'
        call set_fixed_seed(20261103)
        allocate(p, pftc)
        call make_pftc_kernel(p, pftc, 'euclid', 4, 64, 1.0, 48.0, 6.0, 48, 2, 2, sigma2, refs, ptcls, ctfs, kb)
        allocate(before(16,pftc%get_nrots(),3), after(16,pftc%get_nrots(),3))
        call pftc%new_shiftmap(SHIFTMAP_FULL, RHOS, 16, kb(2))
        call pftc%gen_shiftmap(SHIFTMAP_FULL, 1, 1)
        call pftc%shiftmap_objfun_grid(SHIFTMAP_FULL, before)
        call pftc%gen_shiftmap(SHIFTMAP_FULL, 2, 1)
        call pftc%new_shiftmap(SHIFTMAP_FULL, RHOS, 16, kb(2))
        call pftc%gen_shiftmap(SHIFTMAP_FULL, 1, 1)
        call pftc%shiftmap_objfun_grid(SHIFTMAP_FULL, after)
        call assert_true(all(before == after), 'reused geometry preserves the grid exactly')
        call pftc%get_shiftmap_profile(times, counts)
        call assert_true(counts(7) == 1 .and. counts(8) == 1, 'unchanged geometry builds once and reuses once')
        sigma2(:,2) = 1.7 * sigma2(:,2)
        call pftc%memoize_sigma2_weights
        call check_map_objective(pftc, 1, 1)
        call check_map_objective(pftc, 2, 2)
        call pftc%set_ptcl_pft(2, 0.8 * ptcls(:,:,1))
        call pftc%set_ptcl_ctfmat(2, ctfs(:,:,1))
        call pftc%memoize_ptcls
        call check_map_objective(pftc, 1, 2)
        call pftc%reallocate_ptcls(2, [3,4])
        do ip = 3, 4
            call pftc%set_eo(ip, mod(ip,2) == 0)
            call pftc%set_ptcl_pft(ip, ptcls(:,:,ip-2))
            call pftc%set_ptcl_ctfmat(ip, ctfs(:,:,ip-2))
        enddo
        call pftc%memoize_ptcls
        call check_map_objective(pftc, 1, 3)
        call check_map_objective(pftc, 2, 4)
        call pftc%reallocate_ptcls(1, [7])
        call pftc%set_eo(7, .false.)
        call pftc%set_ptcl_pft(7, ptcls(:,:,2))
        call pftc%set_ptcl_ctfmat(7, ctfs(:,:,1))
        call pftc%memoize_ptcls
        call check_map_objective(pftc, 1, 7)
        do cycle = 1, 2
            call pftc%new_shiftmap(SHIFTMAP_FULL, 0.7d0*RHOS, 14, kb(2))
            call check_map_objective(pftc, 1, 7)
            call pftc%new_shiftmap(SHIFTMAP_FULL, RHOS, 16, kb(2)-1)
            call pftc%gen_shiftmap(SHIFTMAP_FULL, 1, 7)
            call pftc%new_shiftmap(SHIFTMAP_FULL, RHOS, 16, kb(2))
            call check_map_objective(pftc, 1, 7)
            call pftc%kill_shiftmap
            call pftc%kill_shiftmap
            call pftc%new_shiftmap(SHIFTMAP_FULL, RHOS, 16, kb(2))
            call check_map_objective(pftc, 1, 7)
        enddo
        call pftc%kill
        deallocate(p, pftc, sigma2, refs, ptcls, ctfs, before, after)
    end subroutine test_map_cache_lifecycle

    subroutine check_map_objective( pftc, iref, iptcl )
        class(polarft_calc), intent(inout) :: pftc
        integer, intent(in) :: iref, iptcl
        real(sp), allocatable :: vals(:,:,:)
        real(dp) :: xy(2), expected, error
        integer :: na, nrho, ir, j, ro
        call pftc%shiftmap_dims(SHIFTMAP_FULL, nrho, na)
        allocate(vals(na,pftc%get_nrots(),nrho))
        call pftc%gen_shiftmap(SHIFTMAP_FULL, iref, iptcl)
        call pftc%shiftmap_objfun_grid(SHIFTMAP_FULL, vals)
        error = 0.d0
        do ir = 1, nrho
            do j = 1, na, max(1,na/4)
                xy = pftc%shiftmap_shift(SHIFTMAP_FULL, ir, j)
                do ro = 1, pftc%get_nrots(), max(1,pftc%get_nrots()/4)
                    expected = pftc%gen_corr_for_rot_8(iref, iptcl, xy, ro)
                    error = max(error, abs(exp(-real(vals(j,ro,ir),dp)) - expected))
                enddo
            enddo
        enddo
        call assert_true(error <= 1.d-6, 'cached map equals the independent shifted/rotated objective')
    end subroutine check_map_objective

    subroutine test_screen_width_identity()
        use simple_eul_prob_tab_utils, only: update_screen_best
        integer, parameter :: NREF = 8, KEEP = 2
        real(dp), parameter :: RHOS(3) = [0.d0, 1.d0, 2.d0]
        class(parameters), allocatable, target :: p
        class(polarft_calc), allocatable :: pftc
        real, allocatable :: sigma2(:,:)
        complex(sp), allocatable :: refs(:,:,:,:), ptcls(:,:,:)
        real(sp), allocatable :: ctfs(:,:,:)
        real :: losses(NREF), widths(NREF), selected_widths(NREF), best(KEEP), work(NREF), loss
        real(dp) :: times(6), times0(6), full_time, selective_time
        integer(int64) :: counts0(8), counts1(8)
        integer :: kb(2), r, j, ro, ir, a, selected, ip, ieo
        write(*,'(A)') 'test_screen_width_identity'
        call set_fixed_seed(20261103)
        allocate(p, pftc)
        call make_pftc_kernel(p, pftc, 'euclid', 1, 64, 1.0, 48.0, 6.0, 48, NREF, 2, sigma2, refs, ptcls, ctfs, kb)
        do ip = 1, 2
            ieo = merge(1, 2, ip == 2)
            call pftc%set_ptcl_ctfmat(ip, 0. * ctfs(:,:,ip) + 1.)
            do r = 1, NREF
                call pftc%set_ref_pft(r, real(r+1)*ptcls(:,:,ip), iseven=(ieo == 1))
            enddo
        enddo
        call pftc%memoize_refs
        call pftc%memoize_ptcls
        call pftc%new_shiftmap(SHIFTMAP_FULL, RHOS, 16, kb(2))
        do ip = 1, 2
            call pftc%get_shiftmap_profile(times0, counts0)
            do r = 1, NREF
                call pftc%gen_shiftmap(SHIFTMAP_FULL, r, ip)
                call pftc%shiftmap_mode(SHIFTMAP_FULL, 6.d0, 0.d0, ro, ir, a, losses(r), &
                    &temp_phys=20.d0, sdev=widths(r))
            enddo
            call pftc%get_shiftmap_profile(times, counts0)
            full_time = sum(times(4:5) - times0(4:5))
            times0 = times
            best = huge(1.)
            selected_widths = -1.
            do r = 1, NREF
                call pftc%gen_shiftmap(SHIFTMAP_FULL, r, ip)
                call pftc%shiftmap_mode(SHIFTMAP_FULL, 6.d0, 0.d0, ro, ir, a, loss, &
                    &temp_phys=20.d0, sdev=selected_widths(r), width_cutoff=maxval(best))
                call assert_true(loss == losses(r), 'width cutoff leaves each candidate loss identical')
                call update_screen_best(loss, best)
            enddo
            call pftc%get_shiftmap_profile(times, counts1)
            selective_time = sum(times(4:5) - times0(4:5))
            write(*,'(A,I0,A,I0,A,2F10.4)') '  screened ', NREF, ', widths ', counts1(5)-counts0(5), &
                &', full/selective posterior ms ', 1.d3*full_time, 1.d3*selective_time
            call assert_true(counts1(5)-counts0(5) < NREF, 'selective width calculation removes measured work')
            call assert_true(counts1(3)-counts0(3) == NREF, 'selective widths keep every reference map')
            call assert_true(counts1(4)-counts0(4) == NREF, 'selective widths use one selection traversal per map')
            work = losses
            do j = 1, KEEP
                selected = minloc(work,dim=1)
                call assert_true(selected_widths(selected) == widths(selected), 'each final shortlisted width is bit identical')
                work(selected) = huge(1.)
            enddo
        enddo
        call pftc%kill
        deallocate(p, pftc, sigma2, refs, ptcls, ctfs)
    end subroutine test_screen_width_identity

    !> the magnitudes 0, step, 2 step, ... up to radius
    subroutine magnitudes( radius, step, rhos )
        real(dp),              intent(in)  :: radius, step
        real(dp), allocatable, intent(out) :: rhos(:)
        integer :: n, i
        n = floor(radius / step) + 1
        allocate(rhos(n))
        do i = 1, n
            rhos(i) = real(i-1,dp) * step
        enddo
    end subroutine magnitudes

    !> The cost of the map against the present coarse grid, per pair, at the in-plane tester's fixture sizes (box 64,
    !! 1.3 A/px, mask 60 A), objfun=euclid with CTF, at lp 8 and lp 2.7 A: the 5 x 5 lattice of
    !! coarse_search_opt_angle (25 all-rotation evaluations of a shifted reference) against the map at the window's
    !! radius (5 px, step lp/4, 2 z_max + 3 directions rounded up to an FFT-friendly count) and, at lp 8, at the
    !! centring disk (30 A). Printed, not asserted.
    subroutine test_kernel_timing()
        use simple_eul_prob_tab_utils, only: posterior_temperature
        integer, parameter :: NLP = 2
        real,    parameter :: LPS(NLP) = [8.0, 2.7], SMPD = 1.3, WINDOW = 5.0, DISK_A = 30.0
        class(parameters),   allocatable, target :: p
        class(polarft_calc), allocatable         :: pftc
        real,        allocatable :: sigma2(:,:)
        complex(sp), allocatable :: refs(:,:,:,:), ptcls(:,:,:)
        real(sp),    allocatable :: ctfs(:,:,:), vals(:,:,:), corrs(:)
        real(dp),    allocatable :: rhos(:)
        integer(timer_int_kind)  :: t0
        real(dp) :: t_coarse, t_win, t_disk, t_draw, step, lse_prior, logp0, temp, logz
        real     :: x, y, loss
        integer  :: kfromto(2), nrots, ilp, irep, ix, iy, nrho, nalpha, nrep, irot, irho, ialpha
        write(*,'(A)') 'test_kernel_timing'
        write(*,'(A)') '  per pair, ms: the coarse 5x5 lattice | the map at the 5 px window | the map at the 30 A disk'
        do ilp = 1, NLP
            call set_fixed_seed(20261104)
            allocate(p, pftc)
            call make_pftc_kernel(p, pftc, 'euclid', 1, 64, SMPD, 60.0, LPS(ilp), 0, 1, 1, sigma2, refs, ptcls, ctfs, kfromto)
            nrots = pftc%get_nrots()
            nrep  = merge(20, 5, ilp == 1)
            allocate(corrs(nrots))
            t0 = tic()
            do irep = 1, nrep
                do ix = 1, 5
                    do iy = 1, 5
                        x = -WINDOW + (real(ix) - 0.5) * 2.0 * WINDOW / 5.0
                        y = -WINDOW + (real(iy) - 0.5) * 2.0 * WINDOW / 5.0
                        call pftc%gen_objfun_vals(1, 1, [x, y], corrs)
                    enddo
                enddo
            enddo
            t_coarse = 1.d3 * real(toc(t0),dp) / real(nrep,dp)
            step = real(LPS(ilp) / SMPD, dp) / 4.d0
            call magnitudes(real(WINDOW,dp), step, rhos)
            nalpha = shiftmap_ndirections(2 * ceiling(DTWOPI * real(kfromto(2),dp) * rhos(size(rhos)) / 64.d0) + 3)
            call pftc%new_shiftmap(SHIFTMAP_FULL, rhos, nalpha, kfromto(2))
            call pftc%shiftmap_dims(SHIFTMAP_FULL, nrho, nalpha)
            allocate(vals(nalpha,nrots,nrho))
            t0 = tic()
            do irep = 1, nrep
                call pftc%gen_shiftmap(SHIFTMAP_FULL, 1, 1)
                call pftc%shiftmap_objfun_grid(SHIFTMAP_FULL, vals)
            enddo
            t_win = 1.d3 * real(toc(t0),dp) / real(nrep,dp)
            write(*,'(A,F5.1,A,I4,A,I3,A,I3,A,I4,A,F9.3,A,F9.3)') '  lp ', LPS(ilp), ' A: nrots ', nrots, ', shells ',&
                &kfromto(2)-kfromto(1)+1, ', window grid ', nrho, ' x ', nalpha, ' | coarse ', t_coarse, ' | window ', t_win
            deallocate(vals, rhos)
            if( ilp == 1 )then
                call magnitudes(real(DISK_A / SMPD,dp), step, rhos)
                nalpha = shiftmap_ndirections(2 * ceiling(DTWOPI * real(kfromto(2),dp) * rhos(size(rhos)) / 64.d0) + 3)
                call pftc%new_shiftmap(SHIFTMAP_FULL, rhos, nalpha, kfromto(2))
                call pftc%shiftmap_dims(SHIFTMAP_FULL, nrho, nalpha)
                allocate(vals(nalpha,nrots,nrho))
                t0 = tic()
                do irep = 1, nrep
                    call pftc%gen_shiftmap(SHIFTMAP_FULL, 1, 1)
                    call pftc%shiftmap_objfun_grid(SHIFTMAP_FULL, vals)
                enddo
                t_disk = 1.d3 * real(toc(t0),dp) / real(nrep,dp)
                call pftc%shiftmap_prior(SHIFTMAP_FULL, 1.d0 / (2.d0 * (10.d0 / SMPD)**2), lse_prior, logp0)
                temp = real(posterior_temperature(pftc%wsqsum_ptcl_now(1), pftc%get_pftsz()), dp)
                t0 = tic()
                do irep = 1, nrep
                    call pftc%shiftmap_draw(SHIFTMAP_FULL, temp, 1.d0 / (2.d0 * (10.d0 / SMPD)**2), lse_prior, .false.,&
                        &logz, irot, irho, ialpha, loss, ordinal=irep)
                enddo
                t_draw = 1.d3 * real(toc(t0),dp) / real(nrep,dp)
                write(*,'(A,I3,A,I4,A,F9.3,A,F9.3)') '           disk grid ', nrho, ' x ', nalpha, ' | disk ', t_disk,&
                    &' | + draw ', t_draw
                deallocate(vals, rhos)
            endif
            call pftc%kill
            deallocate(refs, ptcls, ctfs, corrs, p, pftc)
            if( allocated(sigma2) ) deallocate(sigma2)
        enddo
    end subroutine test_kernel_timing

end module simple_polar_shiftmap_tester
