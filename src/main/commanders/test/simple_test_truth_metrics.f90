!@descr: workflow-gate metrics against a simulation truth: maps docked in both hands, correlation, masked FSC and frame-free pose error
! De novo maps (solve3D) have an arbitrary orientation and hand, so a map is docked in both hands
! (dock_both_hands) before compare_to_truth scores it: whole-volume or, with corr_lp, common-mask
! correlation, plus masked FSC. pair_pose_error needs no docking: relative rotation angles over
! particle pairs are invariant to a global rotation or reflection of the frame.
! The masked FSC is the production compare_volpair (simple_volpair_metrics).
module simple_test_truth_metrics
use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
use simple_core_module_api
use simple_image,           only: image
use simple_ori,             only: ori
use simple_oris,            only: oris
use simple_dock_vols,       only: dock_vols
use simple_imghead,         only: find_ldim_nptcls, find_img_smpd
use simple_volpair_metrics, only: compare_volpair
use simple_test_utils,      only: set_fixed_seed
implicit none

public :: dock_both_hands, compare_to_truth, validate_reconstructed_volume, pair_pose_error, add_gaussian_blob
private
#include "simple_local_flags.inc"

contains

    !> Dock map_fname and its x-mirror onto ref_fname (band hp..lp, mask
    !! diameter mskdiam in A) and write the better of the two, in the
    !! reference frame, to docked_fname; the intermediate maps are named after
    !! tag
    subroutine dock_both_hands( ref_fname, map_fname, mskdiam, hp, lp, tag, docked_fname, cc_direct, cc_mirror )
        class(string),    intent(in)  :: ref_fname, map_fname, docked_fname
        real,             intent(in)  :: mskdiam, hp, lp
        character(len=*), intent(in)  :: tag
        real,             intent(out) :: cc_direct, cc_mirror
        type(dock_vols) :: docker
        type(image)     :: vol
        type(string)    :: mirror, docked_direct, docked_mirror
        real    :: eulers(3), shifts(3), smpd
        integer :: ldim(3), nsections
        mirror        = string(tag//'_mirror.mrc')
        docked_direct = string(tag//'_docked_direct.mrc')
        docked_mirror = string(tag//'_docked_mirror.mrc')
        call find_ldim_nptcls(map_fname, ldim, nsections)
        smpd = find_img_smpd(map_fname)
        call vol%new(ldim, smpd, wthreads=.false.)
        call vol%read(map_fname)
        call vol%mirror('x')
        call vol%write(mirror)
        call vol%kill
        call docker%new(ref_fname, map_fname, smpd, hp, lp, mskdiam)
        call docker%srch()
        call docker%get_dock_info(eulers, shifts, cc_direct)
        call docker%rotate_target(map_fname, docked_direct)
        call docker%kill()
        call docker%new(ref_fname, mirror, smpd, hp, lp, mskdiam)
        call docker%srch()
        call docker%get_dock_info(eulers, shifts, cc_mirror)
        call docker%rotate_target(mirror, docked_mirror)
        call docker%kill()
        write(logfhandle,'(a,f7.4,a,f7.4)') '>>> '//tag//' docking correlation: direct=', cc_direct, &
            &', mirrored=', cc_mirror
        if( cc_direct >= cc_mirror )then
            call simple_copy_file(docked_direct, docked_fname)
        else
            call simple_copy_file(docked_mirror, docked_fname)
        endif
        call mirror%kill
        call docked_direct%kill
        call docked_mirror%kill
    end subroutine dock_both_hands

    !> Two maps in one frame: the whole-volume Pearson correlation by default,
    !! or the common-mask correlation through corr_lp when supplied, and the
    !! FSC=0.5 and FSC=0.143 resolutions (A) of the maps soft-masked to
    !! mskdiam, no finer than Nyquist. A missing map or a non-finite FSC gives
    !! corr 0 and resolutions -1.
    subroutine compare_to_truth( truth_fname, map_fname, mskdiam, corr, fsc05, fsc0143, corr_lp )
        class(string), intent(in)  :: truth_fname, map_fname
        real,          intent(in)  :: mskdiam
        real,          intent(out) :: corr, fsc05, fsc0143
        real, optional, intent(in) :: corr_lp
        type(image) :: truth, map
        real, allocatable :: fsc(:), res(:)
        real    :: smpd, corr_band
        integer :: ldim(3), nsections
        logical :: ok
        corr    = 0.
        fsc05   = -1.
        fsc0143 = -1.
        if( .not. file_exists(truth_fname) .or. .not. file_exists(map_fname) ) return
        call find_ldim_nptcls(map_fname, ldim, nsections)
        smpd = find_img_smpd(map_fname)
        if( present(corr_lp) )then
            call compare_volpair(truth_fname, map_fname, mskdiam, corr_lp, corr_band, fsc, res, ok)
        else
            call truth%new(ldim, smpd, wthreads=.false.)
            call map%new(ldim, smpd, wthreads=.false.)
            call truth%read(truth_fname)
            call map%read(map_fname)
            corr = truth%real_corr(map)
            call truth%kill
            call map%kill
            call compare_volpair(truth_fname, map_fname, mskdiam, 0., corr_band, fsc, res, ok)
        endif
        if( ok )then
            if( present(corr_lp) ) corr = corr_band
            call get_resolution(fsc, res, fsc05, fsc0143)
            if( fsc05   > 0. ) fsc05   = max(fsc05,   2. * smpd)
            if( fsc0143 > 0. ) fsc0143 = max(fsc0143, 2. * smpd)
        endif
    end subroutine compare_to_truth

    !> The final map of a workflow against its simulation truth: the expected
    !! cubic box and sampling (within smpd_tol), docking in both hands (band
    !! dock_hp..dock_lp), a whole-volume correlation of at least min_corr and a
    !! masked FSC=0.143 resolution no worse than max_fsc0143. When corr_lp is
    !! present, the normalized correlation uses maps identically
    !! soft-masked to mask_diameter and low-pass filtered to corr_lp.
    subroutine validate_reconstructed_volume( truth_fname, reconstruction_fname, expected_smpd, expected_box, &
        &smpd_tol, mask_diameter, dock_hp, dock_lp, min_corr, max_fsc0143, corr, fsc0143, &
        &dock_corr_direct, dock_corr_mirrored, dock_corr_selected, passed, corr_lp )
        class(string), intent(in)  :: truth_fname, reconstruction_fname
        real,          intent(in)  :: expected_smpd, smpd_tol, mask_diameter, dock_hp, dock_lp, min_corr, max_fsc0143
        integer,       intent(in)  :: expected_box
        real,          intent(out) :: corr, fsc0143
        real,          intent(out) :: dock_corr_direct, dock_corr_mirrored, dock_corr_selected
        logical,       intent(out) :: passed
        real, optional, intent(in)  :: corr_lp
        character(len=*), parameter :: DOCKED = 'workflow_reconstruction_docked.mrc'
        integer :: truth_ldim(3), reconstruction_ldim(3), nsections
        real    :: truth_smpd, reconstruction_smpd, fsc05
        logical :: corr_ok, fsc_ok
        passed  = .false.
        corr    = 0.
        fsc0143 = 0.
        dock_corr_direct   = 0.
        dock_corr_mirrored = 0.
        dock_corr_selected = 0.
        if( .not. file_exists(truth_fname) )then
            write(logfhandle,'(a)') '    FAIL: simulated truth volume was not generated'
            return
        endif
        if( .not. file_exists(reconstruction_fname) )then
            write(logfhandle,'(a,a)') '    FAIL: final reconstruction was not generated: ', reconstruction_fname%to_char()
            return
        endif
        call find_ldim_nptcls(truth_fname, truth_ldim, nsections)
        truth_smpd = find_img_smpd(truth_fname)
        call find_ldim_nptcls(reconstruction_fname, reconstruction_ldim, nsections)
        reconstruction_smpd = find_img_smpd(reconstruction_fname)
        write(logfhandle,'(a,3(i0,1x),a,f7.3)') '>>> Simulated truth dimensions/sampling: ', truth_ldim, &
            &' / ', truth_smpd
        write(logfhandle,'(a,3(i0,1x),a,f7.3)') '>>> Final volume dimensions/sampling:    ', reconstruction_ldim, &
            &' / ', reconstruction_smpd
        if( any(reconstruction_ldim /= [expected_box, expected_box, expected_box]) )then
            write(logfhandle,'(a,i0)') '    FAIL: final volume does not have the expected cubic box ', expected_box
            return
        endif
        if( abs(reconstruction_smpd - expected_smpd) > smpd_tol )then
            write(logfhandle,'(a,f7.3)') '    FAIL: final volume has incorrect sampling; expected ', expected_smpd
            return
        endif
        if( any(truth_ldim /= reconstruction_ldim) )then
            write(logfhandle,'(a)') '    FAIL: simulated truth and final volume dimensions do not match'
            return
        endif
        if( abs(truth_smpd - reconstruction_smpd) > smpd_tol )then
            write(logfhandle,'(a)') '    FAIL: simulated truth and final volume sampling do not match'
            return
        endif
        call dock_both_hands(truth_fname, reconstruction_fname, mask_diameter, dock_hp, dock_lp, 'workflow_reconstruction', &
            &string(DOCKED), dock_corr_direct, dock_corr_mirrored)
        dock_corr_selected = max(dock_corr_direct, dock_corr_mirrored)
        write(logfhandle,'(a,f7.4,a,f7.2,a,f7.2,a)') '>>> Selected docking correlation: ', dock_corr_selected, &
            &'; band ', dock_hp, '-', dock_lp, ' A'
        if( present(corr_lp) )then
            call compare_to_truth(truth_fname, string(DOCKED), mask_diameter, corr, fsc05, fsc0143, corr_lp)
            write(logfhandle,'(a,f7.2,a,f7.4,a,f7.4)') '>>> Registered soft-masked band correlation to ', &
                &corr_lp, ' A: ', corr, '; minimum ', min_corr
        else
            call compare_to_truth(truth_fname, string(DOCKED), mask_diameter, corr, fsc05, fsc0143)
            write(logfhandle,'(a,f7.4,a,f7.4)') '>>> Registered whole-volume Pearson correlation: ', corr, &
                &'; minimum ', min_corr
        endif
        corr_ok = ieee_is_finite(corr) .and. corr >= min_corr
        if( .not. corr_ok )then
            write(logfhandle,'(a)') '    FAIL: final-volume correlation is below the required minimum'
        endif
        write(logfhandle,'(a,f7.2,a,f7.2,a,f7.2,a)') '>>> Masked truth FSC: 0.500 at ', fsc05, &
            &' A; 0.143 at ', fsc0143, ' A; maximum ', max_fsc0143, ' A'
        fsc_ok = ieee_is_finite(fsc0143) .and. fsc0143 > 0. .and. fsc0143 <= max_fsc0143
        if( .not. fsc_ok ) write(logfhandle,'(a)') '    FAIL: final-volume FSC resolution is outside the accepted range'
        passed = corr_ok .and. fsc_ok
    end subroutine validate_reconstructed_volume

    !> Median over npairs seeded pairs (i from a, j from b, i /= j) of the
    !! difference, in degrees, between the relative rotation angle of the
    !! estimated poses and that of the true poses; frac5 is the fraction of
    !! pairs within 5 degrees
    real function pair_pose_error( estimated, truth, a, b, npairs, seed, frac5 ) result( err )
        class(oris), intent(in)  :: estimated, truth
        integer,     intent(in)  :: a(:), b(:), npairs, seed
        real,        intent(out) :: frac5
        real, allocatable :: errs(:)
        type(ori) :: ea, eb, ta, tb
        real    :: r1, r2
        integer :: k, i1, i2, n
        allocate(errs(npairs))
        call ea%new_ori(.false.)
        call eb%new_ori(.false.)
        call ta%new_ori(.false.)
        call tb%new_ori(.false.)
        call set_fixed_seed(seed)
        n = 0
        do k = 1, npairs
            call random_number(r1)
            call random_number(r2)
            i1 = a(1 + int(r1*real(size(a))))
            i2 = b(1 + int(r2*real(size(b))))
            if( i1 == i2 ) cycle
            call ea%set_euler(estimated%get_euler(i1))
            call eb%set_euler(estimated%get_euler(i2))
            call ta%set_euler(truth%get_euler(i1))
            call tb%set_euler(truth%get_euler(i2))
            n = n + 1
            errs(n) = abs(rad2deg(ea%geodesic_dist_trace(eb)) - rad2deg(ta%geodesic_dist_trace(tb)))
        enddo
        err   = median(errs(1:n))
        frac5 = real(count(errs(1:n) < 5.)) / real(max(1,n))
        call ea%kill
        call eb%kill
        call ta%kill
        call tb%kill
    end function pair_pose_error

    !> Add a Gaussian blob to the map in fname (in place): centre pos (A,
    !! relative to the box centre), width sigma (A) and peak amp times the
    !! map maximum
    subroutine add_gaussian_blob( fname, pos, sigma, amp )
        class(string), intent(in) :: fname
        real,          intent(in) :: pos(3), sigma, amp
        type(image) :: vol
        real, allocatable :: rmat(:,:,:)
        real    :: ctr(3), d2, vmax, smpd
        integer :: ldim(3), nsections, ix, iy, iz
        call find_ldim_nptcls(fname, ldim, nsections)
        smpd = find_img_smpd(fname)
        call vol%new(ldim, smpd)
        call vol%read(fname)
        rmat = vol%get_rmat()
        vmax = maxval(rmat)
        ctr  = real(ldim)/2. + 1.
        do iz = 1, ldim(3)
            do iy = 1, ldim(2)
                do ix = 1, ldim(1)
                    d2 = ((real(ix)-ctr(1))*smpd - pos(1))**2 + ((real(iy)-ctr(2))*smpd - pos(2))**2 + &
                        &((real(iz)-ctr(3))*smpd - pos(3))**2
                    rmat(ix,iy,iz) = rmat(ix,iy,iz) + amp * vmax * exp(-0.5 * d2 / sigma**2)
                enddo
            enddo
        enddo
        call vol%set_rmat(rmat, .false.)
        call vol%write(fname, del_if_exists=.true.)
        call vol%kill
    end subroutine add_gaussian_blob

end module simple_test_truth_metrics
