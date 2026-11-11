!@descr: unit tests of simple_polar_feats: the shift against the calculator's phase ramp, the conjugate wrap, the packed product as the objective's Hermitian product, the feature cache's contract, the network's capped half ring
module simple_polar_feats_tester
use simple_test_utils
use, intrinsic :: iso_fortran_env, only: int64
use simple_defs,               only: sp, dp, DPI, DTWOPI
use simple_string,             only: string
use simple_cmdline,            only: cmdline
use simple_parameters,         only: parameters
use simple_polarft_calc,       only: polarft_calc
use simple_polar_feats,        only: polar_feats, polar_feat, polar_feats_cache, polar_feats_contract, POLAR_FEATS_NHEAD, &
    &polar_feats_net_pftsz
use simple_rnd,                only: gasdev
use simple_syslib,             only: file_exists
use simple_test_polar_fixture, only: rotate_half_ring
implicit none
private
public :: run_all_polar_feats_tests

integer, parameter :: BOX = 48, PFTSZ = 16, NPTCLS = 3, NREFS = 2
real,    parameter :: SMPD = 1.0, MSKDIAM = 40.
integer, parameter :: KFROMTO(2) = [2, 10]

contains

    subroutine run_all_polar_feats_tests()
        write(*,'(A)') '**** running all polar_feats tests ****'
        call test_shift_ramp()
        call test_conjugate_wrap()
        call test_packed_product()
        call test_ctf_ring_mean()
        call test_cache_contract()
        call test_net_half_ring()
    end subroutine run_all_polar_feats_tests

    !> a calculator of NPTCLS random particles without a CTF (the CTF matrices are ones), the Euclidean weights of
    !! two unequal spectra memoised; the particles' rings and the weights are returned beside it
    subroutine make_fixture( p, pftc, sigma2, ptcls, wk )
        type(parameters),              intent(inout) :: p
        type(polarft_calc),            intent(inout) :: pftc
        real(sp), allocatable, target, intent(inout) :: sigma2(:,:)
        complex(sp), allocatable,      intent(inout) :: ptcls(:,:,:)
        real(sp), allocatable,         intent(inout) :: wk(:,:)
        type(cmdline) :: cline
        integer :: iptcl, iref, k, pp, ieo
        complex(sp), allocatable :: ref(:,:)
        call cline%set('box',     real(BOX))
        call cline%set('smpd',    SMPD)
        call cline%set('mskdiam', MSKDIAM)
        call cline%set('pftsz',   real(PFTSZ))
        call cline%set('nptcls',  real(NPTCLS))
        call cline%set('nthr',    1.0)
        call cline%set('ctf',     'no')
        call cline%set('objfun',  'euclid')
        call p%new(cline)
        call set_fixed_seed(20261106)
        call pftc%new(p, NREFS, [1,NPTCLS], KFROMTO)
        allocate(ref(PFTSZ,KFROMTO(1):KFROMTO(2)), ptcls(PFTSZ,KFROMTO(1):KFROMTO(2),NPTCLS))
        do iref = 1,NREFS
            do ieo = 1,2
                do k = KFROMTO(1),KFROMTO(2)
                    do pp = 1,PFTSZ
                        ref(pp,k) = cmplx(gasdev(), gasdev())
                    enddo
                enddo
                call pftc%set_ref_pft(iref, ref, iseven=(ieo == 1))
            enddo
        enddo
        do iptcl = 1,NPTCLS
            do k = KFROMTO(1),KFROMTO(2)
                do pp = 1,PFTSZ
                    ptcls(pp,k,iptcl) = cmplx(gasdev(), gasdev())
                enddo
            enddo
            call pftc%set_eo(iptcl, mod(iptcl,2) == 1)
            call pftc%set_ptcl_pft(iptcl, ptcls(:,:,iptcl))
        enddo
        if( allocated(sigma2) ) deallocate(sigma2)
        allocate(sigma2(KFROMTO(1):KFROMTO(2),2))
        do k = KFROMTO(1),KFROMTO(2)
            sigma2(k,1) = 0.5 + 0.1  * real(k)
            sigma2(k,2) = 0.8 + 0.05 * real(k)
        enddo
        call pftc%assign_sigma2_noise(sigma2)
        call pftc%memoize_refs
        call pftc%memoize_ptcls
        call pftc%memoize_sigma2_weights
        call pftc%get_sigma2_weights(wk)
        deallocate(ref)
        call cline%kill
    end subroutine make_fixture

    !> a half ring shifted through the module's ramp equals the calculator's shift of the same particle
    subroutine test_shift_ramp()
        type(parameters)    :: p
        type(polarft_calc)  :: pftc
        type(polar_feats)   :: feats
        real(sp),    allocatable, target :: sigma2(:,:)
        complex(sp), allocatable :: ptcls(:,:,:), pft_sh(:,:), pft_calc(:,:)
        real(sp),    allocatable :: wk(:,:)
        real(sp), parameter :: SHIFT(2) = [1.3, -0.7]
        real    :: maxdiff, maxval_pft
        integer :: k
        write(*,'(A)') 'test_shift_ramp'
        call make_fixture(p, pftc, sigma2, ptcls, wk)
        call feats%new(PFTSZ, KFROMTO, BOX, wk)
        allocate(pft_sh(PFTSZ,KFROMTO(1):KFROMTO(2)), pft_calc(PFTSZ,KFROMTO(1):KFROMTO(2)))
        call feats%shift_pft(ptcls(:,:,1), SHIFT, pft_sh)
        call pftc%shift_ptcl(1, SHIFT)
        call pftc%get_ptcl_pft(1, pft_calc)
        maxdiff    = maxval(abs(pft_sh - pft_calc))
        maxval_pft = maxval(abs(ptcls(:,:,1)))
        call assert_true(maxdiff <= 1.e-5 * maxval_pft, 'the module''s phase ramp shifts the ring as the calculator does')
        ! a shift of zero leaves the ring as it is, exactly
        call feats%shift_pft(ptcls(:,:,2), [0.,0.], pft_sh)
        call assert_true(all(pft_sh == ptcls(:,:,2)), 'a zero shift is the identity, bit for bit')
        ! the shift keeps every ring's power (a pure phase)
        call feats%shift_pft(ptcls(:,:,3), [5.,3.], pft_sh)
        maxdiff = 0.
        do k = KFROMTO(1),KFROMTO(2)
            maxdiff = max(maxdiff, abs(sum(abs(pft_sh(:,k))**2) - sum(abs(ptcls(:,k,3))**2)))
        enddo
        call assert_true(maxdiff <= 1.e-4 * maxval_pft**2 * real(PFTSZ), 'the shift keeps each ring''s power')
        ! the shift on the packed tensor equals the packing of the shifted ring (the augmentation on a cached tensor)
        call test_shift_feat(feats, ptcls(:,:,3), [5.,3.])
        call feats%kill
        call pftc%kill
        deallocate(ptcls, pft_sh, pft_calc, wk, sigma2)
    end subroutine test_shift_ramp

    subroutine test_shift_feat( feats, pft, shift )
        type(polar_feats), intent(in) :: feats
        complex(sp),       intent(in) :: pft(PFTSZ,KFROMTO(1):KFROMTO(2))
        real(sp),          intent(in) :: shift(2)
        type(polar_feat) :: f0, f1
        complex(sp), allocatable :: pft_sh(:,:)
        real(sp),    allocatable :: ctfmat(:,:), x_sh(:,:,:)
        real :: maxdiff, maxval_x
        allocate(pft_sh(PFTSZ,KFROMTO(1):KFROMTO(2)), ctfmat(PFTSZ,KFROMTO(1):KFROMTO(2)), x_sh(PFTSZ,KFROMTO(1):KFROMTO(2),2))
        ctfmat = 1.
        call feats%build_from_pft(pft, ctfmat, .true., f0, normalise=.false.)
        call feats%shift_pft(pft, shift, pft_sh)
        call feats%build_from_pft(pft_sh, ctfmat, .true., f1, normalise=.false.)
        call feats%shift_feat(f0%x, shift, x_sh)
        maxdiff  = maxval(abs(x_sh - f1%x))
        maxval_x = maxval(abs(f0%x))
        call assert_true(maxdiff <= 1.e-5 * maxval_x, 'the shift of the packed tensor is the packing of the shifted ring')
        call f0%kill
        call f1%kill
        deallocate(pft_sh, ctfmat, x_sh)
    end subroutine test_shift_feat

    !> the conjugate wrap: the padding is the conjugate of the opposite end, a rotation by one grid step is the
    !! Friedel-extended ring's cyclic shift, a half turn is the conjugate, a full turn the identity
    subroutine test_conjugate_wrap()
        type(parameters)    :: p
        type(polarft_calc)  :: pftc
        type(polar_feats)   :: feats
        type(polar_feat)    :: feat
        real(sp),    allocatable, target :: sigma2(:,:)
        type(polar_feat)    :: feat_rot
        complex(sp), allocatable :: ptcls(:,:,:), pft_rot(:,:)
        real(sp),    allocatable :: wk(:,:), xrot(:,:,:), xpad(:,:,:), ctfmat(:,:)
        integer, parameter :: NPAD = 3
        integer :: nk, j, irot
        logical :: ok
        write(*,'(A)') 'test_conjugate_wrap'
        call make_fixture(p, pftc, sigma2, ptcls, wk)
        call feats%new(PFTSZ, KFROMTO, BOX, wk)
        nk = KFROMTO(2) - KFROMTO(1) + 1
        allocate(ctfmat(PFTSZ,KFROMTO(1):KFROMTO(2)), source=1.)
        call feats%build_from_pft(ptcls(:,:,1), ctfmat, .true., feat, normalise=.false.)
        allocate(xrot(PFTSZ,KFROMTO(1):KFROMTO(2),2), xpad(1-NPAD:PFTSZ+NPAD,KFROMTO(1):KFROMTO(2),2))
        ! the padding
        call feats%pad(feat%x, NPAD, xpad)
        ok = all(xpad(1:PFTSZ,:,:) == feat%x)
        do j = 1,NPAD
            ok = ok .and. all(xpad(1-j,:,1) == feat%x(PFTSZ+1-j,:,1)) .and. all(xpad(1-j,:,2) == -feat%x(PFTSZ+1-j,:,2))
            ok = ok .and. all(xpad(PFTSZ+j,:,1) == feat%x(j,:,1)) .and. all(xpad(PFTSZ+j,:,2) == -feat%x(j,:,2))
        enddo
        call assert_true(ok, 'the padding at each end is the conjugate of the opposite end')
        ! the rotations
        call feats%rotate(feat%x, 1, xrot)
        call assert_true(all(xrot == feat%x), 'rotation 1 is the identity')
        call feats%rotate(feat%x, PFTSZ + 1, xrot)
        call assert_true(all(xrot(:,:,1) == feat%x(:,:,1)) .and. all(xrot(:,:,2) == -feat%x(:,:,2)), &
            &'a half turn conjugates the half ring')
        call feats%rotate(feat%x, 2 * PFTSZ + 1, xrot)
        call assert_true(all(xrot == feat%x), 'a full turn is the identity')
        ! one grid step against the Friedel-extended ring's cyclic shift of the complex ring
        allocate(pft_rot(PFTSZ,KFROMTO(1):KFROMTO(2)))
        ok = .true.
        do irot = 2,5
            call rotate_half_ring(PFTSZ, nk, ptcls(:,:,1), irot, pft_rot)
            call feats%build_from_pft(pft_rot, ctfmat, .true., feat_rot, normalise=.false.)
            call feats%rotate(feat%x, irot, xrot)
            ok = ok .and. all(xrot == feat_rot%x)
        enddo
        call assert_true(ok, 'a rotation of the packed tensor is the rotation of the complex ring with the conjugate wrap')
        ! the pooled features of the Hermitian signal: the angular mean of the real part over the half ring is
        ! unchanged by a rotation (the full ring's mean), the imaginary part's changes sign under the wrap
        call feats%rotate(feat%x, 7, xrot)
        call assert_real(sum(feat%x(:,:,1)), sum(xrot(:,:,1)), 1.e-4 * sum(abs(feat%x(:,:,1))), &
            &'the half ring''s mean real part is rotation invariant')
        call feat%kill
        call feat_rot%kill
        call feats%kill
        call pftc%kill
        deallocate(ptcls, pft_rot, wk, sigma2, xrot, xpad, ctfmat)
    end subroutine test_conjugate_wrap

    !> the Euclidean product of two packed tensors is the objective's weighted Hermitian product over the half ring,
    !! half the full ring's; the weighted power is the calculator's self-term; the normalised tensor has unit norm
    subroutine test_packed_product()
        type(parameters)    :: p
        type(polarft_calc)  :: pftc
        type(polar_feats)   :: feats
        type(polar_feat)    :: fa, fb, fn
        real(sp),    allocatable, target :: sigma2(:,:)
        complex(sp), allocatable :: ptcls(:,:,:)
        real(sp),    allocatable :: wk(:,:), ctfa(:,:), ctfb(:,:)
        real(dp) :: packed, herm, herm_full, wsq_calc, norm2
        integer  :: k, pp
        write(*,'(A)') 'test_packed_product'
        call make_fixture(p, pftc, sigma2, ptcls, wk)
        call feats%new(PFTSZ, KFROMTO, BOX, wk)
        allocate(ctfa(PFTSZ,KFROMTO(1):KFROMTO(2)), ctfb(PFTSZ,KFROMTO(1):KFROMTO(2)))
        do k = KFROMTO(1),KFROMTO(2)
            do pp = 1,PFTSZ
                ctfa(pp,k) = 0.3 + 0.7 * abs(sin(0.15 * real(k) + 0.8 * cos(2. * real(pp-1) * real(DTWOPI) / real(2*PFTSZ) + 0.3)))
                ctfb(pp,k) = 0.2 + 0.8 * abs(cos(0.11 * real(k) + 0.5 * sin(2. * real(pp-1) * real(DTWOPI) / real(2*PFTSZ))))
            enddo
        enddo
        ! both particles in the even half: the same weights
        call feats%build_from_pft(ptcls(:,:,1), ctfa, .true., fa, normalise=.false.)
        call feats%build_from_pft(ptcls(:,:,2), ctfb, .true., fb, normalise=.false.)
        packed    = sum(real(fa%x,dp) * real(fb%x,dp))
        herm      = 0.d0
        herm_full = 0.d0
        do k = KFROMTO(1),KFROMTO(2)
            herm = herm + real(wk(k,1),dp) * sum(real(conjg(ptcls(:,k,1) * ctfa(:,k)) * (ptcls(:,k,2) * ctfb(:,k)), dp))
        enddo
        ! the full ring: the half ring and its Friedel mate, whose product has the same real part
        herm_full = 2.d0 * herm
        call assert_true(abs(packed - herm) <= 1.d-5 * abs(herm) + 1.d-6, &
            &'the packed Euclidean product is the objective''s weighted Hermitian product over the half ring')
        call assert_true(abs(2.d0 * packed - herm_full) <= 1.d-5 * abs(herm_full) + 1.d-6, &
            &'twice the half ring''s product is the full ring''s')
        ! the weighted power is the calculator's self-term (the raw particle, k/sigma2 weights)
        wsq_calc = pftc%wsqsum_ptcl_now(1)
        call assert_true(abs(fa%wsq - wsq_calc) <= 1.d-5 * wsq_calc, 'the weighted power equals the calculator''s wsqsum')
        ! the normalised tensor without a CTF has unit weighted norm
        ctfa = 1.
        call feats%build_from_pft(ptcls(:,:,1), ctfa, .true., fn)
        norm2 = sum(real(fn%x,dp)**2)
        call assert_true(abs(norm2 - 1.d0) <= 1.d-5, 'the normalised tensor has unit weighted norm')
        call assert_real(real(log(wsq_calc)), fn%logpow, 1.e-5 * abs(real(log(wsq_calc))), 'the side scalar is the log of the weighted power')
        ! the odd half takes the other spectrum's weights
        call feats%build_from_pft(ptcls(:,:,2), ctfb, .false., fb, normalise=.false.)
        herm = 0.d0
        do k = KFROMTO(1),KFROMTO(2)
            herm = herm + real(wk(k,2),dp) * sum(real(ptcls(:,k,2) * ctfb(:,k) * conjg(ptcls(:,k,2) * ctfb(:,k)), dp))
        enddo
        packed = sum(real(fb%x,dp)**2)
        call assert_true(abs(packed - herm) <= 1.d-5 * abs(herm), 'an odd particle is packed by the odd half''s weights')
        call fa%kill
        call fb%kill
        call fn%kill
        call feats%kill
        call pftc%kill
        deallocate(ptcls, wk, sigma2, ctfa, ctfb)
    end subroutine test_packed_product

    !> the ring mean of the CTF squared and its azimuthal variation; a tensor built from the calculator equals one
    !! built from its arrays
    subroutine test_ctf_ring_mean()
        type(parameters)    :: p
        type(polarft_calc)  :: pftc
        type(polar_feats)   :: feats
        type(polar_feat)    :: fa, fb
        real(sp),    allocatable, target :: sigma2(:,:)
        complex(sp), allocatable :: ptcls(:,:,:)
        real(sp),    allocatable :: wk(:,:), ctfa(:,:)
        real     :: maxdiff
        real(dp) :: tm
        integer  :: k, pp
        write(*,'(A)') 'test_ctf_ring_mean'
        call make_fixture(p, pftc, sigma2, ptcls, wk)
        call feats%new(PFTSZ, KFROMTO, BOX, wk)
        allocate(ctfa(PFTSZ,KFROMTO(1):KFROMTO(2)))
        do k = KFROMTO(1),KFROMTO(2)
            do pp = 1,PFTSZ
                ctfa(pp,k) = 0.1 + 0.9 * abs(sin(0.2 * real(k) + 0.3 * real(pp)))
            enddo
        enddo
        call feats%build_from_pft(ptcls(:,:,1), ctfa, .true., fa)
        maxdiff = 0.
        do k = KFROMTO(1),KFROMTO(2)
            tm = sum(real(ctfa(:,k),dp)**2) / real(PFTSZ,dp)
            maxdiff = max(maxdiff, abs(fa%ctf2(k) - real(tm)))
        enddo
        call assert_true(maxdiff <= 1.e-6, 'ctf2 is the ring mean of the CTF squared')
        call assert_true(fa%tazim > 0., 'an azimuthally varying CTF reports a positive variation')
        ctfa = 0.7
        call feats%build_from_pft(ptcls(:,:,1), ctfa, .true., fb)
        call assert_real(0.49, fb%ctf2(KFROMTO(1)), 1.e-6, 'a flat CTF''s ring mean is its square')
        call assert_real(0., fb%tazim, 1.e-6, 'a flat CTF has no azimuthal variation')
        ! from the calculator (no CTF: ones) against the arrays
        ctfa = 1.
        call feats%build_from_pft(ptcls(:,:,3), ctfa, .true., fa)
        call feats%build(pftc, 3, fb)
        call assert_true(all(fa%x == fb%x) .and. all(fa%ctf2 == fb%ctf2) .and. fa%wsq == fb%wsq .and. fa%logpow == fb%logpow, &
            &'the tensor from the calculator equals the tensor from its arrays, bit for bit')
        call assert_true(fb%iseven, 'the calculator''s half travels with the feature')
        call fa%kill
        call fb%kill
        call feats%kill
        call pftc%kill
        deallocate(ptcls, wk, sigma2, ctfa)
    end subroutine test_ctf_ring_mean

    !> the cache: written under a token and sealed with its contract, it serves the same project, selection, geometry
    !! and weights bit for bit and refuses a change of any of them
    subroutine test_cache_contract()
        type(parameters)    :: p
        type(polarft_calc)  :: pftc
        type(polar_feats)   :: feats, feats_other
        type(polar_feat)    :: feat, back
        type(polar_feats_cache) :: cache
        type(string)        :: cwd_saved, fixture_root, cname
        real(sp),    allocatable, target :: sigma2(:,:)
        complex(sp), allocatable :: ptcls(:,:,:)
        real(sp),    allocatable :: wk(:,:), wk_other(:,:), ctfa(:,:)
        integer(int64) :: key(POLAR_FEATS_NHEAD), other(POLAR_FEATS_NHEAD)
        integer, parameter :: PINDS(3) = [2, 5, 7]
        integer :: nfail0, status, i
        logical :: ok
        write(*,'(A)') 'test_cache_contract'
        nfail0 = tests_failed
        call make_fixture(p, pftc, sigma2, ptcls, wk)
        call feats%new(PFTSZ, KFROMTO, BOX, wk)
        allocate(ctfa(PFTSZ,KFROMTO(1):KFROMTO(2)), source=1.)
        call enter_fixture('polar_feats_cache', cwd_saved, fixture_root)
        cname = 'feats_cache.bin'
        key = polar_feats_contract('/data/a.simple', 101_int64, 12, BOX, BOX, SMPD, PINDS, PFTSZ, KFROMTO, wk, .true.)
        call cache%new(cname)
        call assert_false(cache%probe(feats, key), 'no cache before one is written')
        call cache%create(feats, maxval(PINDS), status)
        call assert_int(0, status, 'the cache is created when the space is there')
        call assert_false(file_exists(cname), 'an unsealed cache is not in place')
        do i = 1,size(PINDS)
            call feats%build_from_pft(ptcls(:,:,i), ctfa, mod(i,2) == 1, feat)
            feat%tazim = real(i)
            call cache%put(PINDS(i), feat)
        enddo
        call cache%seal(key)
        call assert_true(file_exists(cname), 'the sealed cache is in place')
        call assert_true(cache%probe(feats, key), 'the sealed cache serves its contract')
        ! the records back, bit for bit
        call cache%open_read(feats, key)
        ok = .true.
        call feats%alloc_feat(back)
        do i = 1,size(PINDS)
            call feats%build_from_pft(ptcls(:,:,i), ctfa, mod(i,2) == 1, feat)
            feat%tazim = real(i)
            call cache%get(PINDS(i), back)
            ok = ok .and. all(back%x == feat%x) .and. all(back%ctf2 == feat%ctf2) .and. back%tazim == feat%tazim &
                &.and. back%logpow == feat%logpow .and. back%wsq == feat%wsq .and. (back%iseven .eqv. feat%iseven)
        enddo
        call cache%close
        call assert_true(ok, 'the cached features come back bit for bit')
        ! the refusals
        other = polar_feats_contract('/data/b.simple', 101_int64, 12, BOX, BOX, SMPD, PINDS, PFTSZ, KFROMTO, wk, .true.)
        call assert_false(cache%probe(feats, other), 'a changed project path is refused')
        other = polar_feats_contract('/data/a.simple', 102_int64, 12, BOX, BOX, SMPD, PINDS, PFTSZ, KFROMTO, wk, .true.)
        call assert_false(cache%probe(feats, other), 'a replaced project (its time stamp) is refused')
        other = polar_feats_contract('/data/a.simple', 101_int64, 13, BOX, BOX, SMPD, PINDS, PFTSZ, KFROMTO, wk, .true.)
        call assert_false(cache%probe(feats, other), 'a project with other rows is refused')
        other = polar_feats_contract('/data/a.simple', 101_int64, 12, BOX, BOX, SMPD, [2, 5, 8], PFTSZ, KFROMTO, wk, .true.)
        call assert_false(cache%probe(feats, other), 'a changed selection is refused')
        other = polar_feats_contract('/data/a.simple', 101_int64, 12, BOX, BOX, SMPD, PINDS, PFTSZ, [2, 9], wk, .true.)
        call assert_false(cache%probe(feats, other), 'a changed band is refused')
        other = polar_feats_contract('/data/a.simple', 101_int64, 12, BOX, BOX/2, SMPD, PINDS, PFTSZ, KFROMTO, wk, .true.)
        call assert_false(cache%probe(feats, other), 'a changed crop is refused')
        other = polar_feats_contract('/data/a.simple', 101_int64, 12, BOX, BOX, SMPD, PINDS, PFTSZ, KFROMTO, wk, .false.)
        call assert_false(cache%probe(feats, other), 'a changed CTF treatment is refused')
        allocate(wk_other(KFROMTO(1):KFROMTO(2),2), source=wk)
        wk_other(KFROMTO(2),2) = 1.001 * wk_other(KFROMTO(2),2)
        other = polar_feats_contract('/data/a.simple', 101_int64, 12, BOX, BOX, SMPD, PINDS, PFTSZ, KFROMTO, wk_other, .true.)
        call assert_false(cache%probe(feats, other), 'changed shell weights are refused')
        other = polar_feats_contract('/data/a.simple', 101_int64, 12, BOX, BOX, SMPD, [2, 5, 7, 9], PFTSZ, KFROMTO, wk, .true.)
        call assert_false(cache%probe(feats, other), 'a selection beyond the cached rows is refused')
        ! a cache of another geometry does not even read
        call feats_other%new(PFTSZ, [2, 9], BOX, wk(2:9,:))
        other = polar_feats_contract('/data/a.simple', 101_int64, 12, BOX, BOX, SMPD, PINDS, PFTSZ, [2, 9], wk(2:9,:), .true.)
        call assert_false(cache%probe(feats_other, other), 'another geometry''s contract is refused')
        call cache%delete
        call assert_false(file_exists(cname), 'delete removes the cache')
        call cache%kill
        call feat%kill
        call back%kill
        call feats%kill
        call feats_other%kill
        call cname%kill
        call leave_fixture(cwd_saved, fixture_root, nfail0)
        call pftc%kill
        deallocate(ptcls, wk, wk_other, sigma2, ctfa)
    end subroutine test_cache_contract

    !> The network's half ring: at most twice pi k_max angles, even, never above the calculator's (bgal at 15 A: 220
    !! to 140; the medium assembly: 160 kept). A Friedel-symmetric ring whose angular harmonics lie inside the kept
    !! band, given on the calculator's angles, is resampled exactly onto the network's angles, in both channels and
    !! on every ring; the band's power scales with the angles, so the normalised tensor keeps unit weighted norm
    subroutine test_net_half_ring()
        integer, parameter :: N_IN = 72, KB(2) = [1, 4], JMAX = 10
        type(polar_feats)  :: fin, fnet
        type(polar_feat)   :: fa, fb
        complex(sp), allocatable :: pft_in(:,:), pft_net(:,:)
        real(sp),    allocatable :: ctf_in(:,:), ctf_net(:,:), wk(:,:)
        complex(dp) :: a(-JMAX:JMAX,KB(1):KB(2))
        real(dp)    :: re, im, err, big, wsum
        integer     :: m, k, j
        write(*,'(A)') 'test_net_half_ring'
        call assert_int(140, polar_feats_net_pftsz(220, 22), 'bgal at 15 A: 220 angles capped to 140')
        call assert_int(160, polar_feats_net_pftsz(160, 26), 'the medium assembly at 15 A: 160 kept (the cap is 164)')
        call assert_int(170, polar_feats_net_pftsz(480, 27), '400 px at 1 A and 15 A: 480 capped to 170')
        allocate(wk(KB(1):KB(2),2), source=1.)
        do k = KB(1),KB(2)
            wk(k,1) = real(k)
            wk(k,2) = 0.5 * real(k)
        enddo
        call fin%new(N_IN, KB, BOX, wk)
        m = fin%get_pftsz()
        call assert_int(polar_feats_net_pftsz(N_IN, KB(2)), m, 'the features take the cap')
        call assert_true(m < N_IN .and. m > JMAX, 'the fixture resamples and keeps its harmonics')
        call assert_int(N_IN, fin%get_pftsz_in(), 'the calculator''s half ring is kept beside it')
        call fnet%new(m, KB, BOX, wk)
        call assert_int(m, fnet%get_pftsz(), 'a calculator at the cap is not resampled')
        ! Friedel symmetry of the full ring: a(-j) = (-1)^j conj(a(j)), a(0) real
        do k = KB(1),KB(2)
            a(0,k) = cmplx(gasdev(0.,1.), 0., kind=dp)
            do j = 1,JMAX
                re = gasdev(0.,1.)
                im = gasdev(0.,1.)
                a(j,k)  = cmplx(re, im, kind=dp)
                a(-j,k) = real((-1)**j,dp) * conjg(a(j,k))
            enddo
        enddo
        allocate(pft_in(N_IN,KB(1):KB(2)), pft_net(m,KB(1):KB(2)))
        allocate(ctf_in(N_IN,KB(1):KB(2)), ctf_net(m,KB(1):KB(2)), source=1.)
        call sample(N_IN, pft_in)
        call sample(m, pft_net)
        call fin%build_from_pft(pft_in, ctf_in, .true., fa, normalise=.false.)
        call fnet%build_from_pft(pft_net, ctf_net, .true., fb, normalise=.false.)
        big = maxval(abs(real(fb%x,dp)))
        err = maxval(abs(real(fa%x,dp) - real(fb%x,dp)))
        call assert_true(err < 1.e-5_dp * big, 'a ring inside the band resamples to its direct sampling')
        ! the band's power on the network's ring is the calculator's times m / N_IN (normalise scales by it)
        call fin%build_from_pft(pft_in, ctf_in, .false., fa)
        wsum = 0._dp
        do k = KB(1),KB(2)
            wsum = wsum + sum(real(fa%x(:,k,:),dp)**2)
        enddo
        call assert_true(abs(wsum - 1._dp) < 1.e-4_dp, 'the normalised resampled tensor has unit weighted norm')
        call fa%kill
        call fb%kill
        call fin%kill
        call fnet%kill

      contains

        subroutine sample( n, pft )
            integer,     intent(in)  :: n
            complex(sp), intent(out) :: pft(n,KB(1):KB(2))
            complex(dp) :: v
            real(dp)    :: theta
            integer     :: p, kk, jj
            do kk = KB(1),KB(2)
                do p = 1,n
                    theta = real(p-1,dp) * DPI / real(n,dp)
                    v = (0._dp, 0._dp)
                    do jj = -JMAX,JMAX
                        v = v + a(jj,kk) * exp(cmplx(0._dp, real(jj,dp) * theta, kind=dp))
                    enddo
                    pft(p,kk) = cmplx(v, kind=sp)
                enddo
            enddo
        end subroutine sample

    end subroutine test_net_half_ring

end module simple_polar_feats_tester
