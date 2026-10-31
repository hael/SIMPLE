!@descr: unit tests of the polar shift map (simple_polar_shiftmap): the Fourier-Bessel map against closed forms and against the calculator's shifted, rotated cc objective
! The map must equal (1) the closed-form cosine sum of a pure phase ramp on the whole grid and off it,
! (2) gen_corr_for_rot_8 (the reference shifted in its own frame by the polar shift matrix, then rotated,
! objfun=cc) at random grid shifts and rotations including every branch of the rotation, and (3) locate a
! particle that is the reference shifted and rotated by a known pose at that pose, exact and at unit SNR.
module simple_polar_shiftmap_tester
use simple_core_module_api
use simple_cmdline,         only: cmdline
use simple_parameters,      only: parameters
use simple_polarft_calc,    only: polarft_calc
use simple_polar_shiftmap,  only: polar_shiftmap
use simple_test_utils
implicit none
private
public :: run_all_polar_shiftmap_tests

contains

    subroutine run_all_polar_shiftmap_tests()
        write(*,'(A)') '**** running all polar shift map tests ****'
        call test_closed_form_ramp()
        call test_against_pftc_cc()
        call test_peak_recovery()
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

end module simple_polar_shiftmap_tester
