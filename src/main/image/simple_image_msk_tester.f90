!@descr: unit test routines for masks: mask bounds, graphene shells, real-space masks and binary images
module simple_image_msk_tester
use simple_test_utils ! assertions etc.
use simple_defs       ! COSMSKHALFWIDTH, GRAPHENE_BAND1/2, TINY
use simple_image,     only: image, unmemoize_mask_coords
use simple_image_bin, only: image_bin
use simple_math,      only: bounds_from_mask3D
use simple_math_ft,   only: get_resarr, calc_graphene_mask
implicit none
private
public :: run_all_mask_tests, run_all_image_bin_tests

real, parameter :: EPS = 1.0e-5

contains

    !==================================================================
    ! sub-suite 'masks'
    !==================================================================

    subroutine run_all_mask_tests()
        write(*,'(A)') '**** running all mask tests ****'
        call test_bounds_from_mask3D()
        call test_graphene_mask()
        call test_disc_and_cos_edge()
        call test_mask2D_semantics()
        call test_mask3D_semantics()
        call unmemoize_mask_coords ! leave no module state behind
    end subroutine run_all_mask_tests

    !---------------- bounds_from_mask3D ----------------

    ! the bounds must contain every set voxel and be tight for a centred cube;
    ! the search runs inwards from the faces to the centre, so a mask that does
    ! not straddle the centre is still contained, just not tightly
    subroutine test_bounds_from_mask3D()
        integer, parameter :: BOX = 64, RAD = 10
        type(image) :: cube
        logical, allocatable :: mask(:,:,:)
        integer :: lb(3), ub(3), lo(3), hi(3), i
        write(*,'(A)') 'test_bounds_from_mask3D'
        call cube%new([BOX,BOX,BOX], 1.0)
        call cube%square(RAD)
        mask = cube%bin2logical()
        call assert_int((2*RAD)**3, count(mask), 'square() sets (2*rad)^3 voxels')
        call bounds_from_mask3D(mask, lb, ub)
        call brute_force_bounds(mask, lo, hi)
        do i = 1,3
            call assert_int(lo(i), lb(i), 'bounds_from_mask3D: lower bound of a centred cube is tight')
            call assert_int(hi(i), ub(i), 'bounds_from_mask3D: upper bound of a centred cube is tight')
        end do
        call assert_int(BOX/2 - RAD + 1, lb(1), 'bounds_from_mask3D: lower bound is box/2-rad+1')
        call assert_int(BOX/2 + RAD,     ub(1), 'bounds_from_mask3D: upper bound is box/2+rad')
        ! empty mask: the bounds cross
        mask = .false.
        call bounds_from_mask3D(mask, lb, ub)
        call assert_true(all(lb > ub), 'bounds_from_mask3D: empty mask gives crossed bounds')
        ! one voxel off centre: contained
        mask(10, 20, 40) = .true.
        call bounds_from_mask3D(mask, lb, ub)
        call assert_true(lb(1) <= 10 .and. ub(1) >= 10, 'bounds_from_mask3D: off-centre voxel contained (x)')
        call assert_true(lb(2) <= 20 .and. ub(2) >= 20, 'bounds_from_mask3D: off-centre voxel contained (y)')
        call assert_true(lb(3) <= 40 .and. ub(3) >= 40, 'bounds_from_mask3D: off-centre voxel contained (z)')
        call cube%kill
        deallocate(mask)

        contains

            subroutine brute_force_bounds( m, lo, hi )
                logical, intent(in)  :: m(:,:,:)
                integer, intent(out) :: lo(3), hi(3)
                integer :: i
                lo = huge(1); hi = 0
                do i = 1,size(m,1)
                    if( any(m(i,:,:)) )then
                        lo(1) = min(lo(1), i); hi(1) = max(hi(1), i)
                    endif
                end do
                do i = 1,size(m,2)
                    if( any(m(:,i,:)) )then
                        lo(2) = min(lo(2), i); hi(2) = max(hi(2), i)
                    endif
                end do
                do i = 1,size(m,3)
                    if( any(m(:,:,i)) )then
                        lo(3) = min(lo(3), i); hi(3) = max(hi(3), i)
                    endif
                end do
            end subroutine brute_force_bounds

    end subroutine test_bounds_from_mask3D

    !---------------- graphene mask ----------------

    ! the three shells nearest each graphene band are excluded, every other shell kept
    subroutine test_graphene_mask()
        integer, parameter :: BOX = 160
        real,    parameter :: SMPD = 0.358
        real,    allocatable :: res(:), d(:)
        logical, allocatable :: gmask(:), expected(:)
        integer :: n, i, k, loc
        write(*,'(A)') 'test_graphene_mask'
        res   = get_resarr(BOX, SMPD)
        gmask = calc_graphene_mask(BOX, SMPD)
        n     = size(res)
        call assert_int(n, size(gmask), 'graphene mask has one entry per resolution shell')
        call assert_true(res(n) < GRAPHENE_BAND2, 'fixture: both graphene bands lie inside Nyquist')
        allocate(expected(n), source=.true.)
        allocate(d(n))
        do k = 1,2
            d = abs(res - merge(GRAPHENE_BAND1, GRAPHENE_BAND2, k == 1))
            do i = 1,3
                loc = minloc(d, dim=1)
                expected(loc) = .false.
                d(loc) = huge(1.0)
            end do
        end do
        call assert_int(6, count(.not. gmask),          'graphene mask excludes six shells')
        call assert_true(all(gmask .eqv. expected),      'graphene mask excludes the three shells nearest each band')
        deallocate(res, gmask, expected, d)
    end subroutine test_graphene_mask

    !---------------- disc, transfer2bimg, cos_edge ----------------

    subroutine test_disc_and_cos_edge()
        integer, parameter :: BOX = 48, EDGE = 4
        real,    parameter :: RAD = 16.0
        type(image)     :: sph
        type(image_bin) :: bsph
        integer, allocatable :: imat(:,:,:)
        integer :: npix, c, r
        real    :: v, vprev, sphere_vol
        logical :: monotone
        write(*,'(A)') 'test_disc_and_cos_edge'
        call sph%disc([BOX,BOX,BOX], 1.0, RAD, npix)
        sphere_vol = 4.0/3.0 * PI * RAD**3
        call assert_true(abs(real(npix) - sphere_vol) < 0.02 * sphere_vol, 'disc: voxel count within 2% of the sphere volume')
        c = BOX/2 + 1
        call assert_real(1.0, sph%get([c,c,c]),               EPS, 'disc: centre voxel is 1')
        call assert_real(1.0, sph%get([c+nint(RAD),c,c]),     EPS, 'disc: voxel on the radius is 1')
        call assert_real(0.0, sph%get([c+nint(RAD)+1,c,c]),   EPS, 'disc: voxel beyond the radius is 0')
        ! transfer to a binary image preserves the count
        call bsph%transfer2bimg(sph)
        call bsph%get_imat(imat)
        call assert_int(npix, sum(imat), 'transfer2bimg preserves the voxel count')
        ! cosine edge: 1 inside, falling to 0 over EDGE voxels outside the surface
        call bsph%cos_edge(EDGE)
        call assert_real(1.0, bsph%get([c,c,c]),             EPS, 'cos_edge: centre stays 1')
        call assert_real(1.0, bsph%get([c+nint(RAD),c,c]),   EPS, 'cos_edge: surface voxel stays 1')
        call assert_real(0.0, bsph%get([c+nint(RAD)+EDGE,c,c]),   EPS, 'cos_edge: value at the fall-off distance is 0')
        call assert_real(0.0, bsph%get([c+nint(RAD)+EDGE+1,c,c]), EPS, 'cos_edge: value beyond the fall-off is 0')
        monotone = .true.
        vprev    = 1.0
        do r = nint(RAD)+1, nint(RAD)+EDGE
            v = bsph%get([c+r,c,c])
            if( v > vprev + EPS ) monotone = .false.
            vprev = v
        end do
        call assert_true(monotone,                                    'cos_edge: values fall monotonically across the edge')
        call assert_real(0.5, bsph%get([c+nint(RAD)+EDGE/2,c,c]), 1.0e-4, 'cos_edge: half way across the edge is 0.5')
        call sph%kill
        call bsph%kill_bimg
        deallocate(imat)
    end subroutine test_disc_and_cos_edge

    !---------------- 2D masks: hard, soft, softavg ----------------

    ! with mskrad and box such that 2*(mskrad+COSMSKHALFWIDTH) < box the soft mask is
    ! 1 for r <= mskrad-COSMSKHALFWIDTH, 0 for r >= mskrad+COSMSKHALFWIDTH, a cosine between
    subroutine test_mask2D_semantics()
        integer, parameter :: BOX = 64
        real,    parameter :: MSKRAD = 20.0
        type(image) :: img
        integer :: c, r, rin, rout, i, j
        real    :: v, vprev
        logical :: ok
        write(*,'(A)') 'test_mask2D_semantics'
        c    = BOX/2 + 1
        rin  = nint(MSKRAD - COSMSKHALFWIDTH)
        rout = nint(MSKRAD + COSMSKHALFWIDTH)
        ! hard mask
        call img%new([BOX,BOX,1], 1.0)
        img = 1.0
        call img%mask2D_hard(MSKRAD)
        call assert_real(1.0, img%get([c,c,1]),                 EPS, 'mask2D_hard: centre is 1')
        call assert_real(1.0, img%get([c+nint(MSKRAD),c,1]),    EPS, 'mask2D_hard: pixel on the radius is 1')
        call assert_real(0.0, img%get([c+nint(MSKRAD)+1,c,1]),  EPS, 'mask2D_hard: pixel beyond the radius is 0')
        call assert_real(0.0, img%get([1,1,1]),                 EPS, 'mask2D_hard: corner is 0')
        call assert_true(symmetric_2d(img, c),                       'mask2D_hard: symmetric about the origin pixel on both axes')
        ! soft mask
        img = 1.0
        call img%mask2D_soft(MSKRAD, backgr=0.0)
        call assert_real(1.0, img%get([c,c,1]),        EPS, 'mask2D_soft: centre is 1')
        call assert_real(1.0, img%get([c+rin,c,1]),    EPS, 'mask2D_soft: pixel at mskrad-halfwidth is 1')
        call assert_real(0.0, img%get([c+rout,c,1]),   EPS, 'mask2D_soft: pixel at mskrad+halfwidth is 0')
        call assert_real(0.0, img%get([1,1,1]),        EPS, 'mask2D_soft: corner is 0')
        ok    = .true.
        vprev = 1.0
        do r = rin, rout
            v = img%get([c+r,c,1])
            if( v > vprev + EPS .or. v < -EPS .or. v > 1.0 + EPS ) ok = .false.
            vprev = v
        end do
        call assert_true(ok,                                              'mask2D_soft: cosine edge is monotone within [0,1]')
        call assert_real(0.5, img%get([c+nint(MSKRAD),c,1]),   1.0e-3, 'mask2D_soft: pixel on the radius is 0.5')
        call assert_true(symmetric_2d(img, c),                'mask2D_soft: symmetric about the origin pixel on both axes')
        ! soft mask with the background estimated and subtracted: a constant image goes to 0 everywhere
        img = 2.0
        call img%mask2D_soft(MSKRAD)
        call assert_real(0.0, img%get([c,c,1]),        1.0e-4, 'mask2D_soft: constant image minus its background is 0')
        ! softavg: outside is filled with the average of the pixels beyond mskrad
        img = 1.0
        do j = 1,BOX
            do i = 1,BOX
                if( real(i-c)**2 + real(j-c)**2 > MSKRAD**2 ) call img%set([i,j,1], 3.0)
            end do
        end do
        call img%mask2D_softavg(MSKRAD)
        call assert_real(1.0, img%get([c,c,1]),        EPS, 'mask2D_softavg: centre is untouched')
        call assert_real(3.0, img%get([1,1,1]),        EPS, 'mask2D_softavg: corner is the outside average')
        call assert_real(3.0, img%get([c+rout,c,1]),   EPS, 'mask2D_softavg: pixel at mskrad+halfwidth is the average')
        v = img%get([c+nint(MSKRAD),c,1])
        call assert_true(v > 1.0 + 1.0e-3 .and. v < 3.0 - 1.0e-3, 'mask2D_softavg: pixel on the radius is blended')
        call img%kill
    end subroutine test_mask2D_semantics

    !---------------- 3D masks: hard, soft, softavg ----------------

    subroutine test_mask3D_semantics()
        integer, parameter :: BOX = 48
        real,    parameter :: MSKRAD = 12.0
        type(image) :: vol
        integer :: c, r, rin, rout
        real    :: v, vprev
        logical :: ok
        write(*,'(A)') 'test_mask3D_semantics'
        c    = BOX/2 + 1
        rin  = nint(MSKRAD - COSMSKHALFWIDTH)
        rout = nint(MSKRAD + COSMSKHALFWIDTH)
        call vol%new([BOX,BOX,BOX], 1.0)
        vol = 1.0
        call vol%mask3D_hard(MSKRAD)
        call assert_real(1.0, vol%get([c,c,c]),                EPS, 'mask3D_hard: centre is 1')
        call assert_real(1.0, vol%get([c,c+nint(MSKRAD),c]),   EPS, 'mask3D_hard: voxel on the radius is 1')
        call assert_real(0.0, vol%get([c,c+nint(MSKRAD)+1,c]), EPS, 'mask3D_hard: voxel beyond the radius is 0')
        call assert_true(symmetric_3d(vol, c),                       'mask3D_hard: symmetric about the origin voxel on all axes')
        vol = 1.0
        call vol%mask3D_soft(MSKRAD, backgr=0.0)
        call assert_real(1.0, vol%get([c,c,c]),       EPS, 'mask3D_soft: centre is 1')
        call assert_real(1.0, vol%get([c,c,c+rin]),   EPS, 'mask3D_soft: voxel at mskrad-halfwidth is 1')
        call assert_real(0.0, vol%get([c,c,c+rout]),  EPS, 'mask3D_soft: voxel at mskrad+halfwidth is 0')
        call assert_real(0.0, vol%get([1,1,1]),       EPS, 'mask3D_soft: corner is 0')
        ok    = .true.
        vprev = 1.0
        do r = rin, rout
            v = vol%get([c,c,c+r])
            if( v > vprev + EPS .or. v < -EPS .or. v > 1.0 + EPS ) ok = .false.
            vprev = v
        end do
        call assert_true(ok, 'mask3D_soft: cosine edge is monotone within [0,1]')
        call assert_true(symmetric_3d(vol, c), 'mask3D_soft: symmetric about the origin voxel on all axes')
        vol = 1.0
        call vol%mask3D_softavg(MSKRAD)
        call assert_real(1.0, vol%get([c,c,c]),   EPS, 'mask3D_softavg: constant volume stays constant inside')
        call assert_real(1.0, vol%get([1,1,1]),   EPS, 'mask3D_softavg: constant volume stays constant outside')
        call vol%kill
    end subroutine test_mask3D_semantics


    !---------------- helpers ----------------

    ! a mask must read the same at +r and -r from the origin pixel c along every axis
    ! (and the same along x and y): this is what the mirrored loops promise
    logical function symmetric_2d( img, c )
        type(image), intent(in) :: img
        integer,     intent(in) :: c
        integer :: r
        symmetric_2d = .true.
        do r = 1, c-2
            if( abs(img%get([c+r,c,1]) - img%get([c-r,c,1])) > EPS ) symmetric_2d = .false.
            if( abs(img%get([c,c+r,1]) - img%get([c,c-r,1])) > EPS ) symmetric_2d = .false.
            if( abs(img%get([c+r,c,1]) - img%get([c,c+r,1])) > EPS ) symmetric_2d = .false.
            if( abs(img%get([c+r,c+r,1]) - img%get([c-r,c-r,1])) > EPS ) symmetric_2d = .false.
        end do
    end function symmetric_2d

    logical function symmetric_3d( vol, c )
        type(image), intent(in) :: vol
        integer,     intent(in) :: c
        integer :: r
        symmetric_3d = .true.
        do r = 1, c-2
            if( abs(vol%get([c+r,c,c]) - vol%get([c-r,c,c])) > EPS ) symmetric_3d = .false.
            if( abs(vol%get([c,c+r,c]) - vol%get([c,c-r,c])) > EPS ) symmetric_3d = .false.
            if( abs(vol%get([c,c,c+r]) - vol%get([c,c,c-r])) > EPS ) symmetric_3d = .false.
            if( abs(vol%get([c+r,c,c]) - vol%get([c,c,c+r])) > EPS ) symmetric_3d = .false.
            if( abs(vol%get([c+r,c+r,c+r]) - vol%get([c-r,c-r,c-r])) > EPS ) symmetric_3d = .false.
        end do
    end function symmetric_3d

    !==================================================================
    ! sub-suite 'binary image'
    !==================================================================

    subroutine run_all_image_bin_tests()
        write(*,'(A)') '**** running all binary image tests ****'
        call test_image_bin_examples()
        call test_ccs_connectivity_3D()
        call test_ccs_connectivity_2D()
        call test_ccs_degenerate()
        call test_morphology()
        call test_cc_bookkeeping()
        call test_holes_and_feret()
    end subroutine run_all_image_bin_tests

    !---------------- the original image_bin examples, with their answers ----------------

    subroutine test_image_bin_examples()
        type(image_bin) :: bimg, bimg3D, ccimg
        integer, allocatable :: labels(:,:,:)
        integer :: imat(4,4,1), imat3D(3,3,2), nccs
        real    :: dist
        write(*,'(A)') 'test_image_bin_examples'
        call bimg%new_bimg([4,4,1], 1.0)
        imat = 0
        imat(1,1,1) = 1
        imat(2,1,1) = 1
        imat(2,3,1) = 1
        imat(4,4,1) = 1
        call bimg%set_imat(imat)
        call bimg%max_dist(dist)
        call assert_real(sqrt(8.0), dist, 1.0e-5,        'max_dist: farthest set pixel from the centre (2,2)')
        call bimg%find_ccs(ccimg)
        call ccimg%get_nccs(nccs)
        call assert_int(3, nccs,                          'find_ccs: three 8-connected components')
        call ccimg%get_imat(labels)
        call assert_int(4, count(labels > 0),             'find_ccs: every set pixel is labelled')
        call assert_int(3, maxval(labels),                'find_ccs: labels run 1..nccs')
        call assert_int(labels(1,1,1), labels(2,1,1),     'find_ccs: adjacent pixels share a label')
        call assert_true(labels(2,3,1) /= labels(2,1,1),  'find_ccs: pixels two apart are separate components')
        call assert_true(labels(4,4,1) /= labels(2,3,1),  'find_ccs: the isolated corner is its own component')
        ! all ones
        imat = 1
        call bimg%set_imat(imat)
        call bimg%max_dist(dist)
        call assert_real(sqrt(8.0), dist, 1.0e-5,        'max_dist: full 4x4 image')
        call bimg%find_ccs(ccimg)
        call ccimg%get_nccs(nccs)
        call assert_int(1, nccs,                          'find_ccs: a full image is one component')
        ! all zeros
        imat = 0
        call bimg%set_imat(imat)
        call bimg%max_dist(dist)
        call assert_real(0.0, dist, 1.0e-5,              'max_dist: empty image is 0')
        call bimg%find_ccs(ccimg)
        call ccimg%get_nccs(nccs)
        call assert_int(0, nccs,                          'find_ccs: an empty image has no components')
        ! 3D: one voxel in the far corner of a 3x3x2 box, centre (2,2,1)
        call bimg3D%new_bimg([3,3,2], 1.0)
        imat3D = 0
        imat3D(3,3,2) = 1
        call bimg3D%set_imat(imat3D)
        call bimg3D%max_dist(dist)
        call assert_real(sqrt(3.0), dist, 1.0e-5,        'max_dist: 3D far corner')
        call bimg%kill_bimg
        call bimg3D%kill_bimg
        call ccimg%kill_bimg
        deallocate(labels)
    end subroutine test_image_bin_examples

    !---------------- 26-connectivity in 3D (the regression of the missing axial neighbours) ----------------

    subroutine test_ccs_connectivity_3D()
        integer, parameter :: N = 7, C = 4
        type(image_bin) :: bimg, ccimg
        integer, allocatable :: sz(:)
        integer :: imat(N,N,N), nccs, di, dj, dk, nfail
        write(*,'(A)') 'test_ccs_connectivity_3D'
        call bimg%new_bimg([N,N,N], 1.0)
        ! every one of the 26 offsets connects
        nfail = 0
        do dk = -1,1
            do dj = -1,1
                do di = -1,1
                    if( di == 0 .and. dj == 0 .and. dk == 0 ) cycle
                    imat = 0
                    imat(C,C,C)          = 1
                    imat(C+di,C+dj,C+dk) = 1
                    call bimg%set_imat(imat)
                    call bimg%find_ccs(ccimg)
                    call ccimg%get_nccs(nccs)
                    if( nccs /= 1 ) nfail = nfail + 1
                end do
            end do
        end do
        call assert_int(0, nfail, 'find_ccs is 26-connected in 3D (offsets that failed to connect)')
        ! Chebyshev distance 2 stays separate
        nfail = 0
        do dk = -2,2,2
            do dj = -2,2,2
                do di = -2,2,2
                    if( di == 0 .and. dj == 0 .and. dk == 0 ) cycle
                    imat = 0
                    imat(C,C,C)          = 1
                    imat(C+di,C+dj,C+dk) = 1
                    call bimg%set_imat(imat)
                    call bimg%find_ccs(ccimg)
                    call ccimg%get_nccs(nccs)
                    if( nccs /= 2 ) nfail = nfail + 1
                end do
            end do
        end do
        call assert_int(0, nfail, 'find_ccs keeps voxels two steps apart separate (offsets that merged)')
        ! the regression case: a 1x1xN column along z, and the controls along x and y
        imat = 0; imat(C,C,:) = 1
        call bimg%set_imat(imat)
        call bimg%find_ccs(ccimg)
        call ccimg%get_nccs(nccs)
        call assert_int(1, nccs, 'a 1x1xN column along z is one component')
        sz = ccimg%size_ccs()
        call assert_int(1, size(sz), 'column along z: one component size')
        call assert_int(N, sz(1),    'column along z: N voxels')
        imat = 0; imat(:,C,C) = 1
        call bimg%set_imat(imat)
        call bimg%find_ccs(ccimg)
        call ccimg%get_nccs(nccs)
        call assert_int(1, nccs, 'a 1x1xN column along x is one component')
        imat = 0; imat(C,:,C) = 1
        call bimg%set_imat(imat)
        call bimg%find_ccs(ccimg)
        call ccimg%get_nccs(nccs)
        call assert_int(1, nccs, 'a 1x1xN column along y is one component')
        call bimg%kill_bimg
        call ccimg%kill_bimg
        if( allocated(sz) ) deallocate(sz)
    end subroutine test_ccs_connectivity_3D

    !---------------- 8-connectivity in 2D ----------------

    subroutine test_ccs_connectivity_2D()
        integer, parameter :: N = 7, C = 4
        type(image_bin) :: bimg, ccimg
        integer :: imat(N,N,1), nccs, di, dj, nfail
        write(*,'(A)') 'test_ccs_connectivity_2D'
        call bimg%new_bimg([N,N,1], 1.0)
        nfail = 0
        do dj = -1,1
            do di = -1,1
                if( di == 0 .and. dj == 0 ) cycle
                imat = 0
                imat(C,C,1)       = 1
                imat(C+di,C+dj,1) = 1
                call bimg%set_imat(imat)
                call bimg%find_ccs(ccimg)
                call ccimg%get_nccs(nccs)
                if( nccs /= 1 ) nfail = nfail + 1
            end do
        end do
        call assert_int(0, nfail, 'find_ccs is 8-connected in 2D (offsets that failed to connect)')
        imat = 0
        imat(C,C,1)   = 1
        imat(C+2,C,1) = 1
        call bimg%set_imat(imat)
        call bimg%find_ccs(ccimg)
        call ccimg%get_nccs(nccs)
        call assert_int(2, nccs, 'find_ccs keeps pixels two apart separate in 2D')
        call bimg%kill_bimg
        call ccimg%kill_bimg
    end subroutine test_ccs_connectivity_2D

    !---------------- degenerate volumes ----------------

    subroutine test_ccs_degenerate()
        integer, parameter :: N = 5
        type(image_bin) :: bimg, ccimg
        integer, allocatable :: sz(:)
        integer :: imat(N,N,N), nccs
        write(*,'(A)') 'test_ccs_degenerate'
        call bimg%new_bimg([N,N,N], 1.0)
        imat = 0
        call bimg%set_imat(imat)
        call bimg%find_ccs(ccimg)
        call ccimg%get_nccs(nccs)
        call assert_int(0, nccs, 'an empty volume has no components')
        imat = 1
        call bimg%set_imat(imat)
        call bimg%find_ccs(ccimg)
        call ccimg%get_nccs(nccs)
        call assert_int(1, nccs, 'a full volume is one component')
        sz = ccimg%size_ccs()
        call assert_int(N**3, sz(1), 'full volume: N^3 voxels in the component')
        call bimg%kill_bimg
        call ccimg%kill_bimg
        if( allocated(sz) ) deallocate(sz)
    end subroutine test_ccs_degenerate

    !---------------- morphology: dilate, erode, grow_bins ----------------

    ! a 10x10 square: erode strips the outer layer (8x8), dilate puts it back exactly,
    ! dilate once more gives 12x12; grow_bins uses a disc template, so corners are not added
    subroutine test_morphology()
        integer, parameter :: BOX = 32, SIDE = 10
        type(image_bin) :: bimg
        integer, allocatable :: imat(:,:,:), imat0(:,:,:)
        integer :: i0
        write(*,'(A)') 'test_morphology'
        i0 = BOX/2 - SIDE/2 + 1
        call bimg%new_bimg([BOX,BOX,1], 1.0)
        allocate(imat0(BOX,BOX,1), source=0)
        imat0(i0:i0+SIDE-1, i0:i0+SIDE-1, 1) = 1
        call bimg%set_imat(imat0)
        call bimg%erode
        call bimg%get_imat(imat)
        call assert_int((SIDE-2)**2, sum(imat),                      'erode: the outer layer is removed')
        call assert_true(all(imat(i0+1:i0+SIDE-2, i0+1:i0+SIDE-2, 1) == 1), 'erode: the interior survives')
        call bimg%dilate
        call bimg%get_imat(imat)
        call assert_int(SIDE**2, sum(imat),                          'dilate after erode restores the square')
        call assert_true(all(imat == imat0),                         'dilate after erode restores it exactly')
        call bimg%dilate
        call bimg%get_imat(imat)
        call assert_int((SIDE+2)**2, sum(imat),                      'dilate adds one layer all round, corners included')
        ! grow_bins(1): a cross template, so the corners are not added
        call bimg%set_imat(imat0)
        call bimg%grow_bins(1)
        call bimg%get_imat(imat)
        call assert_int((SIDE+2)**2 - 4, sum(imat),                  'grow_bins(1) adds one layer without the corners')
        ! grow_bins(2) on a single pixel: the 13-pixel digital disc of radius 2
        imat0 = 0
        imat0(BOX/2, BOX/2, 1) = 1
        call bimg%set_imat(imat0)
        call bimg%grow_bins(2)
        call bimg%get_imat(imat)
        call assert_int(13, sum(imat),                               'grow_bins(2) on a pixel is the digital disc of radius 2')
        call assert_int(1, imat(BOX/2+2, BOX/2, 1),                  'grow_bins(2): two pixels along an axis are in')
        call assert_int(0, imat(BOX/2+2, BOX/2+1, 1),                'grow_bins(2): (2,1) is outside the disc')
        call bimg%kill_bimg
        deallocate(imat, imat0)
    end subroutine test_morphology

    !---------------- connected-component bookkeeping ----------------

    ! three blobs of known size and position: labels, sizes, centres, diameters, and the
    ! elimination/relabelling/extraction routines that the pickers and cavg tools use
    subroutine test_cc_bookkeeping()
        integer, parameter :: BOX = 48
        real,    parameter :: SMPD = 2.0
        type(image_bin) :: bimg, ccimg
        integer, allocatable :: imat(:,:,:), sz(:)
        integer :: nccs, lab_a, lab_b, lab_c
        real    :: xy(2), diam
        write(*,'(A)') 'test_cc_bookkeeping'
        call bimg%new_bimg([BOX,BOX,1], SMPD)
        allocate(imat(BOX,BOX,1), source=0)
        imat( 5: 6,  5: 6, 1) = 1     ! A: 2x2 = 4 pixels
        imat(10:14, 20:24, 1) = 1     ! B: 5x5 = 25 pixels, centre (12,22)
        imat(30:39, 30:39, 1) = 1     ! C: 10x10 = 100 pixels
        call bimg%set_imat(imat)
        call bimg%find_ccs(ccimg)
        call ccimg%get_nccs(nccs)
        call assert_int(3, nccs, 'find_ccs: three blobs')
        call ccimg%get_imat(imat)
        lab_a = imat(5,5,1); lab_b = imat(12,22,1); lab_c = imat(35,35,1)
        call assert_true(lab_a /= lab_b .and. lab_b /= lab_c .and. lab_a /= lab_c, 'find_ccs: distinct labels')
        sz = ccimg%size_ccs()
        call assert_int(3, size(sz),      'size_ccs: one size per component')
        call assert_int(4,   sz(lab_a),   'size_ccs: blob A')
        call assert_int(25,  sz(lab_b),   'size_ccs: blob B')
        call assert_int(100, sz(lab_c),   'size_ccs: blob C')
        ! centre of mass of B relative to the image centre (box/2+1 = 25)
        call ccimg%masscen_cc(lab_b, xy)
        call assert_real(12.0 - 25.0, xy(1), 1.0e-5, 'masscen_cc: x offset from the image centre')
        call assert_real(22.0 - 25.0, xy(2), 1.0e-5, 'masscen_cc: y offset from the image centre')
        ! diameter of B: twice the farthest pixel from its centre of mass, in Angstroms
        call ccimg%diameter_cc(lab_b, diam)
        call assert_real(2.0 * sqrt(8.0) * SMPD, diam, 1.0e-3, 'diameter_cc: 5x5 block, corner at sqrt(8) pixels')
        ! cc2bin keeps one component as a binary image
        call ccimg%cc2bin(lab_b)
        call ccimg%get_imat(imat)
        call assert_int(25, sum(imat),              'cc2bin: only the chosen component remains')
        call assert_int(1,  maxval(imat),           'cc2bin: as a 0/1 image')
        call assert_int(1,  imat(12,22,1),          'cc2bin: the chosen component is where it was')
        ! elim_ccs by size: keep the sizes in [10,50], relabel contiguously
        call bimg%find_ccs(ccimg)
        call ccimg%elim_ccs([10, 50])
        call ccimg%get_nccs(nccs)
        call assert_int(1, nccs,                    'elim_ccs: one component survives the size window')
        sz = ccimg%size_ccs()
        call assert_int(1,  size(sz),               'elim_ccs: sizes of the survivors only')
        call assert_int(25, sz(1),                  'elim_ccs: the survivor is blob B')
        call ccimg%get_imat(imat)
        call assert_int(1, imat(12,22,1),           'elim_ccs: the survivor is relabelled 1')
        ! order_ccs closes the gaps in the labelling
        imat = 0
        imat( 5: 6,  5: 6, 1) = 3
        imat(10:14, 20:24, 1) = 7
        call ccimg%set_imat(imat)
        call ccimg%order_ccs
        call ccimg%get_imat(imat)
        call assert_int(1, imat(5,5,1),             'order_ccs: the lowest label becomes 1')
        call assert_int(2, imat(12,22,1),           'order_ccs: the next label becomes 2')
        call ccimg%get_nccs(nccs)
        call assert_int(2, nccs,                    'order_ccs: nccs is the number of labels')
        call bimg%kill_bimg
        call ccimg%kill_bimg
        deallocate(imat)
        if( allocated(sz) ) deallocate(sz)
    end subroutine test_cc_bookkeeping

    !---------------- hole filling and Feret diameters ----------------

    subroutine test_holes_and_feret()
        integer, parameter :: BOX = 40
        type(image_bin) :: bimg
        integer, allocatable :: imat(:,:,:)
        real :: fmin, fmax
        write(*,'(A)') 'test_holes_and_feret'
        call bimg%new_bimg([BOX,BOX,1], 1.0)
        allocate(imat(BOX,BOX,1), source=0)
        ! a square ring: 20x20 with a 10x10 hole
        imat(11:30, 11:30, 1) = 1
        imat(16:25, 16:25, 1) = 0
        call bimg%set_imat(imat)
        call bimg%set_edgecc2background
        call bimg%get_imat(imat)
        call assert_int(400, sum(imat),                       'set_edgecc2background fills the hole')
        call assert_int(1,   imat(20,20,1),                   'set_edgecc2background: the hole centre is foreground')
        call assert_int(0,   imat(1,1,1),                     'set_edgecc2background: the outside stays background')
        ! Feret diameters of a 5 x 21 axis-aligned bar (pixel centres plus one pixel)
        imat = 0
        imat(18:22, 10:30, 1) = 1
        call bimg%set_imat(imat)
        call bimg%feret_minmax(fmin, fmax)
        call assert_real(5.0, fmin, 1.0e-3,                   'feret_minmax: minimum Feret of a 5-wide bar is 5')
        call assert_true(fmax > 21.0 .and. fmax < 21.5,       'feret_minmax: maximum Feret is the bar diagonal (21.4)')
        call bimg%kill_bimg
        deallocate(imat)
    end subroutine test_holes_and_feret

end module simple_image_msk_tester
