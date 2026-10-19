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
use simple_test_gate,       only: test_gate
implicit none

public :: dock_both_hands, compare_to_truth, validate_reconstructed_volume, pair_pose_error, add_gaussian_blob
public :: match_states, adjusted_rand_index, mapped_accuracy, accuracy_null, collapsed_accuracy, population_errors
public :: match_maps, rot_angle, carry_pose, group_errors, group_operators, axis_operator, discrimination
public :: finite_min, offdiag, maps_distinct, state_recalls
private
#include "simple_local_flags.inc"

contains

    !> Dock map_fname and its x-mirror onto ref_fname (band hp..lp, mask
    !! diameter mskdiam in A) and write the better of the two, in the
    !! reference frame, to docked_fname; the intermediate maps are named after
    !! tag
    subroutine dock_both_hands( ref_fname, map_fname, mskdiam, hp, lp, tag, docked_fname, cc_direct, cc_mirror, &
        &eul, shift, mirrored )
        class(string),    intent(in)  :: ref_fname, map_fname, docked_fname
        real,             intent(in)  :: mskdiam, hp, lp
        character(len=*), intent(in)  :: tag
        real,             intent(out) :: cc_direct, cc_mirror
        real,    optional, intent(out) :: eul(3), shift(3) !< the transformation of the better hand
        logical, optional, intent(out) :: mirrored         !< the better hand is the x-mirror
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
        if( present(eul)   ) eul   = eulers
        if( present(shift) ) shift = shifts
        call docker%rotate_target(map_fname, docked_direct)
        call docker%kill()
        call docker%new(ref_fname, mirror, smpd, hp, lp, mskdiam)
        call docker%srch()
        call docker%get_dock_info(eulers, shifts, cc_mirror)
        if( cc_mirror > cc_direct )then
            if( present(eul)   ) eul   = eulers
            if( present(shift) ) shift = shifts
        endif
        if( present(mirrored) ) mirrored = cc_mirror > cc_direct
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
    !! soft-masked to mask_diameter and low-pass filtered to corr_lp. Every
    !! decision is also a row of gate, named after tag.
    subroutine validate_reconstructed_volume( truth_fname, reconstruction_fname, expected_smpd, expected_box, &
        &smpd_tol, mask_diameter, dock_hp, dock_lp, min_corr, max_fsc0143, gate, tag, corr, fsc0143, &
        &dock_corr_direct, dock_corr_mirrored, dock_corr_selected, passed, corr_lp, corr_gated, fsc_gated )
        class(string),    intent(in)    :: truth_fname, reconstruction_fname
        real,             intent(in)    :: expected_smpd, smpd_tol, mask_diameter, dock_hp, dock_lp, min_corr, max_fsc0143
        integer,          intent(in)    :: expected_box
        class(test_gate), intent(inout) :: gate
        character(len=*), intent(in)    :: tag
        real,             intent(out)   :: corr, fsc0143
        real,             intent(out)   :: dock_corr_direct, dock_corr_mirrored, dock_corr_selected
        logical,          intent(out)   :: passed
        real, optional,   intent(in)    :: corr_lp
        !> the correlation and the FSC held to their bounds (the default), or reported only (a robustness
        !! profile, where the profile admits no floor) and left out of passed
        logical, optional, intent(in)   :: corr_gated, fsc_gated
        character(len=*), parameter :: DOCKED = 'workflow_reconstruction_docked.mrc'
        real    :: fsc05
        logical :: geom_ok, corr_ok, fsc_ok, l_corr_gated, l_fsc_gated
        l_corr_gated = .true.
        if( present(corr_gated) ) l_corr_gated = corr_gated
        l_fsc_gated = .true.
        if( present(fsc_gated) ) l_fsc_gated = fsc_gated
        passed  = .false.
        corr    = 0.
        fsc0143 = 0.
        dock_corr_direct   = 0.
        dock_corr_mirrored = 0.
        dock_corr_selected = 0.
        corr_ok = .false.
        fsc_ok  = .false.
        geom_ok = geometry_matches()
        if( geom_ok )then
            call dock_both_hands(truth_fname, reconstruction_fname, mask_diameter, dock_hp, dock_lp, &
                &'workflow_reconstruction', string(DOCKED), dock_corr_direct, dock_corr_mirrored)
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
            passed = (corr_ok .or. .not. l_corr_gated) .and. (fsc_ok .or. .not. l_fsc_gated)
        endif
        call gate%check(tag//'_map_written_with_expected_box_and_sampling', geom_ok)
        if( geom_ok )then
            if( l_corr_gated )then
                call gate%metric(tag//'_map_truth_corr', corr, min_corr, corr_ok)
            else
                call gate%report(tag//'_map_truth_corr', corr)
            endif
            if( l_fsc_gated )then
                call gate%metric(tag//'_map_truth_fsc0143_A', fsc0143, max_fsc0143, fsc_ok, upper=.true.)
            else
                call gate%report(tag//'_map_truth_fsc0143_A', fsc0143)
            endif
            call gate%report(tag//'_dock_corr_direct',   dock_corr_direct)
            call gate%report(tag//'_dock_corr_mirrored', dock_corr_mirrored)
        endif

    contains

        logical function geometry_matches()
            integer :: truth_ldim(3), reconstruction_ldim(3), nsections
            real    :: truth_smpd, reconstruction_smpd
            geometry_matches = .false.
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
            if( .not. (abs(reconstruction_smpd - expected_smpd) <= smpd_tol) )then
                write(logfhandle,'(a,f7.3)') '    FAIL: final volume has incorrect sampling; expected ', expected_smpd
                return
            endif
            if( any(truth_ldim /= reconstruction_ldim) )then
                write(logfhandle,'(a)') '    FAIL: simulated truth and final volume dimensions do not match'
                return
            endif
            if( .not. (abs(truth_smpd - reconstruction_smpd) <= smpd_tol) )then
                write(logfhandle,'(a)') '    FAIL: simulated truth and final volume sampling do not match'
                return
            endif
            geometry_matches = .true.
        end function geometry_matches

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

    ! ---- state labels ----------------------------------------------------------------------------

    !> The one-to-one mapping of estimated states (1 to nest) onto true states (1 to ntrue) that
    !! maximises the number of particles whose mapped label is their true one, over every
    !! permutation (at most eight states); map(e) = 0 for an estimated state left unmatched. A label
    !! outside its range (0: no state) matches nothing.
    function match_states( truth, est, ntrue, nest ) result( map )
        integer, intent(in)  :: truth(:), est(:), ntrue, nest
        integer, allocatable :: map(:)
        integer :: conf(ntrue, nest), perm(max(ntrue, nest)), best(max(ntrue, nest)), k, i, bestscore
        logical :: used(max(ntrue, nest))
        k = max(ntrue, nest)
        if( k > 8 ) THROW_HARD('match_states compares at most eight states')
        if( size(truth) /= size(est) ) THROW_HARD('match_states: label arrays differ in length')
        conf = confusion(truth, est, ntrue, nest)
        bestscore = -1
        used = .false.
        call try(1, 0)
        allocate(map(nest))
        do i = 1, nest
            map(i) = best(i)
            if( map(i) > ntrue ) map(i) = 0
        enddo

      contains

        recursive subroutine try( pos, score )
            integer, intent(in) :: pos, score
            integer :: t, add
            if( pos > k )then
                if( score > bestscore )then
                    bestscore = score
                    best      = perm
                endif
                return
            endif
            do t = 1, k
                if( used(t) ) cycle
                used(t)   = .true.
                perm(pos) = t
                add = 0
                if( pos <= nest .and. t <= ntrue ) add = conf(t, pos)
                call try(pos + 1, score + add)
                used(t) = .false.
            enddo
        end subroutine try

    end function match_states

    !> counts of (true state, estimated state) pairs; labels outside their range are left out
    function confusion( truth, est, ntrue, nest ) result( conf )
        integer, intent(in) :: truth(:), est(:), ntrue, nest
        integer :: conf(ntrue, nest), i
        conf = 0
        do i = 1, size(truth)
            if( truth(i) < 1 .or. truth(i) > ntrue .or. est(i) < 1 .or. est(i) > nest ) cycle
            conf(truth(i), est(i)) = conf(truth(i), est(i)) + 1
        enddo
    end function confusion

    !> The adjusted Rand index of two partitions (Hubert and Arabie): 1 for identical partitions up
    !! to a renaming, 0 in expectation for independent ones whatever the populations. Labels below
    !! 1 count as one more class.
    real function adjusted_rand_index( a, b ) result( ari )
        integer, intent(in) :: a(:), b(:)
        integer, allocatable :: ia(:), ib(:)
        real(dp), allocatable :: n(:,:), ra(:), rb(:)
        real(dp) :: index, expected, maxindex, pairs
        integer :: na, nb, i
        ia = max(a, 0) + 1
        ib = max(b, 0) + 1
        na = maxval(ia)
        nb = maxval(ib)
        allocate(n(na, nb), source=0._dp)
        do i = 1, size(ia)
            n(ia(i), ib(i)) = n(ia(i), ib(i)) + 1._dp
        enddo
        ra       = sum(n, dim=2)
        rb       = sum(n, dim=1)
        index    = sum(n * (n - 1._dp)) / 2._dp
        pairs    = real(size(ia), dp) * real(size(ia) - 1, dp) / 2._dp
        expected = (sum(ra * (ra - 1._dp)) / 2._dp) * (sum(rb * (rb - 1._dp)) / 2._dp) / max(pairs, 1._dp)
        maxindex = (sum(ra * (ra - 1._dp)) / 2._dp + sum(rb * (rb - 1._dp)) / 2._dp) / 2._dp
        if( abs(maxindex - expected) < 1.e-12_dp )then
            ! both partitions trivial (one class each, or all singletons): agreement or nothing
            ari = merge(1., 0., abs(index - maxindex) < 1.e-12_dp)
        else
            ari = real((index - expected) / (maxindex - expected))
        endif
    end function adjusted_rand_index

    !> the fraction of particles whose estimated state, under map, is their true one
    real function mapped_accuracy( truth, est, map ) result( acc )
        integer, intent(in) :: truth(:), est(:), map(:)
        integer :: i, ncorrect
        ncorrect = 0
        do i = 1, size(truth)
            if( est(i) < 1 .or. est(i) > size(map) ) cycle
            if( map(est(i)) == truth(i) .and. truth(i) > 0 ) ncorrect = ncorrect + 1
        enddo
        acc = real(ncorrect) / real(max(1, size(truth)))
    end function mapped_accuracy

    !> The null distribution of the matched accuracy: the true labels shuffled across particles
    !! (population counts kept) nshuffle times, each scored after its own best matching; its median
    !! and its 99.9th percentile
    subroutine accuracy_null( truth, est, ntrue, nest, nshuffle, seed, med, p999 )
        use simple_rnd, only: shuffle
        integer, intent(in)  :: truth(:), est(:), ntrue, nest, nshuffle, seed
        real,    intent(out) :: med, p999
        integer, allocatable :: shuffled(:), map(:)
        real,    allocatable :: accs(:)
        integer :: i
        allocate(accs(nshuffle))
        call set_fixed_seed(seed)
        shuffled = truth
        do i = 1, nshuffle
            call shuffle(shuffled)
            map     = match_states(shuffled, est, ntrue, nest)
            accs(i) = mapped_accuracy(shuffled, est, map)
        enddo
        med  = median(accs)
        call hpsort(accs)
        p999 = accs(min(nshuffle, ceiling(0.999 * real(nshuffle))))
    end subroutine accuracy_null

    !> the accuracy of a predictor that puts every particle in one state: the largest population
    real function collapsed_accuracy( truth, ntrue ) result( acc )
        integer, intent(in) :: truth(:), ntrue
        integer :: k, nmax
        nmax = 0
        do k = 1, ntrue
            nmax = max(nmax, count(truth == k))
        enddo
        acc = real(nmax) / real(max(1, size(truth)))
    end function collapsed_accuracy

    !> per true state, the fraction of its particles found in the estimated state mapped onto it (0
    !! when no estimated state maps onto it): a merge of two true states leaves one of them with a
    !! recall near 0 however well the labels otherwise agree
    function state_recalls( truth, est, map, ntrue ) result( rec )
        integer, intent(in) :: truth(:), est(:), map(:), ntrue
        real :: rec(ntrue)
        integer :: k, e, ntrue_k
        rec = 0.
        do k = 1, ntrue
            ntrue_k = count(truth == k)
            if( ntrue_k == 0 ) cycle
            e = findloc(map, k, dim=1)
            if( e < 1 ) cycle
            rec(k) = real(count(truth == k .and. est == e)) / real(ntrue_k)
        enddo
    end function state_recalls

    !> per true state, the absolute difference between its population fraction and that of the
    !! estimated states mapped onto it
    function population_errors( truth, est, map, ntrue ) result( errs )
        integer, intent(in) :: truth(:), est(:), map(:), ntrue
        real :: errs(ntrue)
        integer :: k, e, nest_k
        do k = 1, ntrue
            nest_k = 0
            do e = 1, size(map)
                if( map(e) == k ) nest_k = nest_k + count(est == e)
            enddo
            errs(k) = abs(real(nest_k) - real(count(truth == k))) / real(max(1, size(truth)))
        enddo
    end function population_errors

    ! ---- maps ------------------------------------------------------------------------------------

    !> The one-to-one mapping of estimated maps onto true maps that maximises the summed
    !! correlation cc(estimated, true), by the same search as match_states (at most eight)
    function match_maps( cc ) result( map )
        real, intent(in)     :: cc(:,:)  !< (estimated, true)
        integer, allocatable :: map(:)
        integer :: perm(max(size(cc,1), size(cc,2))), best(max(size(cc,1), size(cc,2))), k, ne, nt, i
        logical :: used(max(size(cc,1), size(cc,2)))
        real    :: bestscore
        ne = size(cc,1)
        nt = size(cc,2)
        k  = max(ne, nt)
        if( k > 8 ) THROW_HARD('match_maps compares at most eight maps')
        best      = [(i, i = 1, k)]
        bestscore = -huge(1.)
        used = .false.
        call try(1, 0.)
        allocate(map(ne))
        do i = 1, ne
            map(i) = best(i)
            if( map(i) > nt ) map(i) = 0
        enddo

      contains

        recursive subroutine try( pos, score )
            integer, intent(in) :: pos
            real,    intent(in) :: score
            integer :: t
            real    :: add
            if( pos > k )then
                if( score > bestscore )then
                    bestscore = score
                    best      = perm
                endif
                return
            endif
            do t = 1, k
                if( used(t) ) cycle
                used(t)   = .true.
                perm(pos) = t
                add = 0.
                ! a non-finite correlation counts as no agreement
                if( pos <= ne .and. t <= nt )then
                    if( ieee_is_finite(cc(pos, t)) ) add = cc(pos, t)
                endif
                call try(pos + 1, score + add)
                used(t) = .false.
            enddo
        end subroutine try

    end function match_maps

    !> The discrimination of estimated maps est (docked into the frame of the truths) between the
    !! true maps truth, inside the difference mask: where some pair of truths differs by more than
    !! frac of the larger maximum, dilated by dilate pixels. margins(k) is the correlation of est(k)
    !! with truth(k) less its largest correlation with another truth; diffcorr(k,l) correlates
    !! est(k)-est(l) with truth(k)-truth(l), and is not finite when the estimated maps are equal.
    !! corr(k,l), when present, is the correlation of est(k) with truth(l) inside the mask.
    subroutine discrimination( est, truth, frac, dilate, margins, diffcorr, nmask, corr )
        class(image),   intent(inout) :: est(:), truth(:)
        real,           intent(in)    :: frac
        integer,        intent(in)    :: dilate
        real,           intent(out)   :: margins(size(truth)), diffcorr(size(truth), size(truth))
        integer,        intent(out)   :: nmask
        real, optional, intent(out)   :: corr(size(truth), size(truth))
        real,    allocatable :: t(:,:,:,:), e(:,:,:,:)
        logical, allocatable :: diff(:,:,:), msk(:,:,:)
        real    :: cc(size(truth), size(truth))
        integer :: ldim(3), k, l, n, i, j, z, di, dj, dz
        n    = size(truth)
        if( n < 2 .or. size(est) /= n ) THROW_HARD('discrimination needs two or more truths and one map per truth')
        ldim = truth(1)%get_ldim()
        allocate(t(ldim(1),ldim(2),ldim(3),n), e(ldim(1),ldim(2),ldim(3),n))
        do k = 1, n
            t(:,:,:,k) = truth(k)%get_rmat()
            e(:,:,:,k) = est(k)%get_rmat()
        enddo
        allocate(diff(ldim(1),ldim(2),ldim(3)), source=.false.)
        do k = 1, n
            do l = k + 1, n
                diff = diff .or. abs(t(:,:,:,k) - t(:,:,:,l)) > frac * max(maxval(abs(t(:,:,:,k))), maxval(abs(t(:,:,:,l))))
            enddo
        enddo
        ! dilation by a cube of half-width dilate
        allocate(msk(ldim(1),ldim(2),ldim(3)), source=.false.)
        do z = 1, ldim(3)
            do j = 1, ldim(2)
                do i = 1, ldim(1)
                    if( .not. diff(i,j,z) ) cycle
                    do dz = max(1, z-dilate), min(ldim(3), z+dilate)
                        do dj = max(1, j-dilate), min(ldim(2), j+dilate)
                            do di = max(1, i-dilate), min(ldim(1), i+dilate)
                                msk(di,dj,dz) = .true.
                            enddo
                        enddo
                    enddo
                enddo
            enddo
        enddo
        nmask = count(msk)
        do k = 1, n
            do l = 1, n
                cc(k,l) = masked_corr(e(:,:,:,k), t(:,:,:,l), msk)
            enddo
        enddo
        do k = 1, n
            margins(k) = cc(k,k) - maxval(cc(k,:), mask=[(l /= k, l = 1, n)])
        enddo
        if( present(corr) ) corr = cc
        diffcorr = 0.
        do k = 1, n
            do l = 1, n
                if( l == k ) cycle
                diffcorr(k,l) = masked_corr(e(:,:,:,k) - e(:,:,:,l), t(:,:,:,k) - t(:,:,:,l), msk)
            enddo
        enddo
    end subroutine discrimination

    !> the smallest value, not finite when any value is not (a non-finite value fails its metric)
    real function finite_min( x )
        use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan
        real, intent(in) :: x(:)
        if( size(x) == 0 )then
            finite_min = ieee_value(finite_min, ieee_quiet_nan)
        else if( .not. all(ieee_is_finite(x)) )then
            finite_min = ieee_value(finite_min, ieee_quiet_nan)
        else
            finite_min = minval(x)
        endif
    end function finite_min

    !> the off-diagonal elements of a square matrix, column by column
    function offdiag( a ) result( x )
        real, intent(in)  :: a(:,:)
        real, allocatable :: x(:)
        integer :: k, l, n
        allocate(x(size(a,1) * size(a,2) - min(size(a,1), size(a,2))))
        n = 0
        do l = 1, size(a,2)
            do k = 1, size(a,1)
                if( k == l ) cycle
                n = n + 1
                x(n) = a(k,l)
            enddo
        enddo
    end function offdiag

    !> whether maps count as distinct states by the output of discrimination: every margin and
    !! every difference correlation finite and positive (copies of one map, or of their mean, are not)
    logical function maps_distinct( margins, diffcorr )
        real, intent(in) :: margins(:), diffcorr(:,:)
        real :: m, d
        m = finite_min(margins)
        d = finite_min(offdiag(diffcorr))
        maps_distinct = ieee_is_finite(m) .and. ieee_is_finite(d)
        if( maps_distinct ) maps_distinct = m > 0. .and. d > 0.
    end function maps_distinct

    !> Pearson correlation inside mask; not finite when either side is constant there
    real function masked_corr( a, b, msk ) result( cc )
        use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan
        real,    intent(in) :: a(:,:,:), b(:,:,:)
        logical, intent(in) :: msk(:,:,:)
        real(dp) :: ma, mb, sab, saa, sbb
        integer  :: n
        n = count(msk)
        if( n < 2 )then
            cc = ieee_value(cc, ieee_quiet_nan)
            return
        endif
        ma  = sum(real(a, dp), mask=msk) / real(n, dp)
        mb  = sum(real(b, dp), mask=msk) / real(n, dp)
        sab = sum((real(a, dp) - ma) * (real(b, dp) - mb), mask=msk)
        saa = sum((real(a, dp) - ma)**2, mask=msk)
        sbb = sum((real(b, dp) - mb)**2, mask=msk)
        if( saa <= 1.e-30_dp .or. sbb <= 1.e-30_dp )then
            cc = ieee_value(cc, ieee_quiet_nan)
        else
            cc = real(sab / sqrt(saa * sbb))
        endif
    end function masked_corr

    ! ---- poses -----------------------------------------------------------------------------------

    !> the rotation angle (deg) between two rotation matrices
    real function rot_angle( a, b )
        real, intent(in) :: a(3,3), b(3,3)
        rot_angle = rad2deg(acos(max(-1., min(1., (sum(a*b) - 1.)/2.))))
    end function rot_angle

    !> A particle's pose r and shift sh (pixels) relative to an estimated map, carried into the frame
    !! of the truth the map was docked onto: the docked map is rotvol of the map, or of its x-mirror
    !! when mirrored, by the rotation a and the shift s3 (pixels). The x-mirror of a box of even size
    !! box also moves the map by one pixel. The rules are pinned by the frame test of the tester.
    subroutine carry_pose( r, sh, mirrored, box, a, s3, r_out, sh_out )
        real,    intent(in)  :: r(3,3), sh(2), a(3,3), s3(3)
        logical, intent(in)  :: mirrored
        integer, intent(in)  :: box
        real,    intent(out) :: r_out(3,3), sh_out(2)
        real, parameter :: MZ(3,3) = reshape([1.,0.,0., 0.,1.,0., 0.,0.,-1.], [3,3])
        real, parameter :: MX(3,3) = reshape([-1.,0.,0., 0.,1.,0., 0.,0.,1.], [3,3])
        real :: r1(3,3), sh1(2), t(3)
        r1  = r
        sh1 = sh
        if( mirrored )then
            r1 = matmul(MZ, matmul(r, MX))
            if( mod(box, 2) == 0 )then
                t   = matmul(r, [-1., 0., 0.])
                sh1 = sh - t(1:2)
            endif
        endif
        r_out  = matmul(r1, transpose(a))
        t      = matmul(r1, s3)
        sh_out = sh1 - t(1:2)
    end subroutine carry_pose

    !> Per particle, the rotation error (deg) of the estimated pose against the true one up to the
    !! operators ops of the truth's declared point group: min over g of angle(r_est, r_true g).
    !! With projection_inversion, also up to the in-plane rotation by 180 degrees (Rz180 r_true g):
    !! the image of a centre-inverted object, so a nearly centrosymmetric object (a crystal lattice)
    !! is posed only up to it
    function group_errors( r_est, r_true, ops, projection_inversion ) result( errs )
        real,              intent(in) :: r_est(:,:,:), r_true(:,:,:), ops(:,:,:)
        logical, optional, intent(in) :: projection_inversion
        real, parameter :: RZ180(3,3) = reshape([-1.,0.,0., 0.,-1.,0., 0.,0.,1.], [3,3])
        real, allocatable :: errs(:)
        real    :: rtg(3,3)
        integer :: j, g
        logical :: l_inv
        l_inv = .false.
        if( present(projection_inversion) ) l_inv = projection_inversion
        allocate(errs(size(r_est, 3)), source=huge(1.))
        do j = 1, size(r_est, 3)
            do g = 1, size(ops, 3)
                rtg     = matmul(r_true(:,:,j), ops(:,:,g))
                errs(j) = min(errs(j), rot_angle(r_est(:,:,j), rtg))
                if( l_inv ) errs(j) = min(errs(j), rot_angle(r_est(:,:,j), matmul(RZ180, rtg)))
            enddo
        enddo
    end function group_errors

    !> the rotation matrices of SIMPLE's point group pgrp, in its standard frame
    function group_operators( pgrp ) result( ops )
        use simple_sym, only: sym
        character(len=*), intent(in) :: pgrp
        real, allocatable :: ops(:,:,:)
        type(sym) :: se
        integer   :: i
        call se%new(pgrp)
        allocate(ops(3,3,se%get_nsym()))
        do i = 1, se%get_nsym()
            call se%get_sym_rmat(i, ops(:,:,i))
        enddo
        call se%kill
    end function group_operators

    !> the rotation by ang degrees about the unit axis u
    function axis_operator( u, ang ) result( r )
        real, intent(in) :: u(3), ang
        real :: r(3,3), c, s, k(3,3)
        integer :: i
        c = cos(deg2rad(ang))
        s = sin(deg2rad(ang))
        k = reshape([0., u(3), -u(2), -u(3), 0., u(1), u(2), -u(1), 0.], [3,3])
        r = s * k + (1. - c) * matmul(k, k)
        do i = 1, 3
            r(i,i) = r(i,i) + 1.
        enddo
    end function axis_operator

end module simple_test_truth_metrics
