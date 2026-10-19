!@descr: unit tests for the truth metrics of the multi-state tests (simple_test_truth_metrics): state matching, adjusted Rand index, null accuracy, frame and pose conventions
module simple_test_truth_metrics_tester
use simple_core_module_api
use simple_image,              only: image
use simple_ori,                only: ori
use simple_oris,               only: oris
use simple_test_utils,         only: assert_true, assert_int, assert_real, set_fixed_seed
use simple_test_truth_metrics, only: match_states, adjusted_rand_index, mapped_accuracy, accuracy_null, &
    &collapsed_accuracy, population_errors, match_maps, discrimination, carry_pose, group_errors, group_operators, &
    &axis_operator, dock_both_hands, finite_min, offdiag, maps_distinct, state_recalls
implicit none

public :: run_all_truth_metrics_tests, run_all_truth_docking_tests
private

integer, parameter :: B = 32
real,    parameter :: SMPD_T = 2.

contains

    subroutine run_all_truth_metrics_tests()
        write(*,'(A)') '**** running all truth metrics tests ****'
        call test_matching()
        call test_ari()
        call test_null()
        call test_map_mapping()
        call test_discrimination()
        call test_frame_rule()
        call test_group_error()
        call test_nonfinite()
    end subroutine run_all_truth_metrics_tests

    subroutine run_all_truth_docking_tests()
        write(*,'(A)') '**** running all truth docking tests ****'
        call test_docking_frame()
    end subroutine run_all_truth_docking_tests

    !> End to end: an estimated map docked in both hands onto a truth made from its x-mirror by a
    !! known rotation and shift; dock_both_hands finds the mirrored hand and a transformation under
    !! which carry_pose takes every particle's pose and shift into the truth's frame
    subroutine test_docking_frame()
        use simple_simple_volinterp, only: rotvol
        type(image) :: t, tm, w
        type(ori)   :: od
        real    :: pos(3,4), a(3,3), s3(3), r(3,3), sh(2), rc(3,3), shc(2), eul(3), shift(3), ccd, ccm
        logical :: mirrored
        integer :: i
        write(*,'(A)') 'test_docking_frame'
        pos = reshape([10.,4.,-6., -5.,9.,3., 2.,-8.,7., -9.,-3.,-4.], [3,4])
        call blob_volume(t, pos, [3.,2.5,2.,3.5], [1.,.8,.7,.5])
        call t%write(string('tmp_truth_docking_map.mrc'))
        tm = t
        call tm%mirror('x')
        a  = euler2m([20., 50., 75.])
        s3 = [1.5, -2., 1.]
        call od%new(.false.)
        call od%set_euler(m2euler(a))
        w = rotvol(tm, od, s3)
        call w%write(string('tmp_truth_docking_truth.mrc'))
        call dock_both_hands(string('tmp_truth_docking_truth.mrc'), string('tmp_truth_docking_map.mrc'), 56., 40., 6., &
            &'tmp_truth_docking', string('tmp_truth_docking_docked.mrc'), ccd, ccm, eul, shift, mirrored)
        call assert_true(mirrored, 'docking finds the mirrored hand')
        call assert_true(ccm > 0.95 .and. ccm > ccd, 'the mirrored hand docks better')
        do i = 1, 3
            r  = euler2m([37.*real(i), 63.+10.*real(i), 112.-20.*real(i)])
            sh = [1.-0.5*real(i), -1.5+real(i)]
            call carry_pose(r, sh, mirrored, B, euler2m(eul), shift, rc, shc)
            call assert_true(proj_corr(t, r, sh, w, rc, shc) > 0.95, &
                &'an image keeps its pose and shift carried by the docking found (particle '//int2str(i)//')')
        enddo
        call t%kill
        call tm%kill
        call w%kill
        call od%kill
        call del_file('tmp_truth_docking_map.mrc')
        call del_file('tmp_truth_docking_truth.mrc')
        call del_file('tmp_truth_docking_docked.mrc')
        call del_file('tmp_truth_docking_mirror.mrc')
        call del_file('tmp_truth_docking_docked_direct.mrc')
        call del_file('tmp_truth_docking_docked_mirror.mrc')
    end subroutine test_docking_frame

    !> a volume of Gaussian blobs at pos (A from the centre)
    subroutine blob_volume( vol, pos, sig, amp )
        type(image), intent(inout) :: vol
        real,        intent(in)    :: pos(:,:), sig(:), amp(:)
        real, allocatable :: rmat(:,:,:)
        real    :: c, d2
        integer :: i, j, k, n
        allocate(rmat(B,B,B), source=0.)
        c = real(B)/2. + 1.
        do k = 1, B
            do j = 1, B
                do i = 1, B
                    do n = 1, size(sig)
                        d2 = ((real(i)-c)*SMPD_T - pos(1,n))**2 + ((real(j)-c)*SMPD_T - pos(2,n))**2 + ((real(k)-c)*SMPD_T - pos(3,n))**2
                        rmat(i,j,k) = rmat(i,j,k) + amp(n) * exp(-d2/(2.*sig(n)**2))
                    enddo
                enddo
            enddo
        enddo
        call vol%new([B,B,B], SMPD_T, wthreads=.false.)
        call vol%set_rmat(rmat, .false.)
    end subroutine blob_volume

    real function proj_corr( vola, ea, sha, volb, eb, shb )
        type(image), intent(inout) :: vola, volb
        real,        intent(in)    :: ea(3,3), eb(3,3), sha(2), shb(2)
        type(image) :: ia, ib
        call project(vola, ea, sha, ia)
        call project(volb, eb, shb, ib)
        proj_corr = ia%real_corr(ib)
        call ia%kill
        call ib%kill
    end function proj_corr

    !> the projection of vol at rotation r, shifted by sh pixels
    subroutine project( vol, r, sh, img )
        use simple_projector, only: projector
        type(image), intent(inout) :: vol
        real,        intent(in)    :: r(3,3), sh(2)
        type(image), intent(inout) :: img
        type(projector) :: vproj
        type(image)     :: plane
        type(ori)       :: e
        call vproj%new(OSMPL_PAD_FAC*[B,B,B], SMPD_T, wthreads=.false.)
        call vol%pad(vproj)
        call vproj%fft()
        call vproj%expand_cmat()
        call plane%new([OSMPL_PAD_FAC*B,OSMPL_PAD_FAC*B,1], SMPD_T, wthreads=.false.)
        call e%new(.false.)
        call e%set_euler(m2euler(r))
        call vproj%fproject_serial(e, plane)
        call plane%shift([sh(1), sh(2), 0.])
        call plane%ifft()
        call img%new([B,B,1], SMPD_T, wthreads=.false.)
        call plane%clip(img)
        call e%kill
        call plane%kill
        call vproj%kill_expanded()
        call vproj%kill()
    end subroutine project

    !> the matching of estimated onto true states, and its invariance under a renaming
    subroutine test_matching()
        integer, parameter :: RENAME(3) = [2,3,1]
        integer :: truth(10), est(10), renamed(10), i
        integer, allocatable :: map(:), map2(:)
        write(*,'(A)') 'test_matching'
        truth = [1,1,1,1,1,2,2,2,3,3]
        ! estimated states named 3, 1, 2 for true 1, 2, 3, with one particle of true state 1 wrong
        est   = [3,3,3,3,1,1,1,1,2,2]
        map   = match_states(truth, est, 3, 3)
        call assert_true(all(map == [2,3,1]), 'estimated states 3, 1, 2 map onto true 1, 2, 3')
        call assert_real(0.9, mapped_accuracy(truth, est, map), 1.e-6, 'nine of ten particles are right under the mapping')
        ! a global renaming of the estimated states changes neither the accuracy nor the ARI
        do i = 1, 10
            renamed(i) = RENAME(est(i))
        enddo
        map2 = match_states(truth, renamed, 3, 3)
        call assert_real(mapped_accuracy(truth, est, map), mapped_accuracy(truth, renamed, map2), 0., &
            &'the accuracy does not change when the estimated states are renamed')
        call assert_real(adjusted_rand_index(truth, est), adjusted_rand_index(truth, renamed), 0., &
            &'nor does the adjusted Rand index')
        ! more estimated than true states: the extra one maps onto none
        map = match_states([1,1,2,2], [1,1,2,3], 2, 3)
        call assert_true(map(1) == 1 .and. map(2) == 2 .and. map(3) == 0, 'an extra estimated state maps onto none')
        call assert_int(3, size(population_errors([1,1,2,2,3,3], [1,1,2,2,3,3], [1,2,3], 3)), &
            &'one population error per true state')
        call assert_real(0., maxval(population_errors([1,1,2,2,3,3], [2,2,3,3,1,1], [3,1,2], 3)), 1.e-6, &
            &'the populations of a renamed perfect partition agree')
        call assert_real(1., minval(state_recalls([1,1,2,2,3,3], [2,2,3,3,1,1], [3,1,2], 3)), 1.e-6, &
            &'every state of a renamed perfect partition is recalled in full')
        ! true states 2 and 3 merged into estimated state 2, estimated state 3 empty
        map = match_states([1,1,2,2,3,3], [1,1,2,2,2,2], 3, 3)
        call assert_real(0., minval(state_recalls([1,1,2,2,3,3], [1,1,2,2,2,2], map, 3)), 1.e-6, &
            &'a merge leaves one true state unrecalled')
    end subroutine test_matching

    !> adjusted Rand indices against pair counting, done independently in Python (itertools over all
    !! pairs: (both - expected)/(max - expected)); 1 for identical partitions up to renaming
    subroutine test_ari()
        write(*,'(A)') 'test_ari'
        call assert_real(0.24242424, adjusted_rand_index([1,1,1,2,2,2], [1,1,2,2,3,3]), 1.e-6, 'ARI of a published 6-item example')
        call assert_real(-0.5, adjusted_rand_index([1,1,2,2], [1,2,1,2]), 1.e-6, 'ARI of two crossed 2x2 partitions')
        call assert_real(0.39114391, adjusted_rand_index([1,1,1,1,2,2,2,3,3,3], [1,1,1,2,2,2,2,3,3,1]), 1.e-6, &
            &'ARI of a 10-item example with two errors')
        call assert_real(1., adjusted_rand_index([1,1,2,2,3], [3,3,1,1,2]), 1.e-6, 'a renamed identical partition scores 1')
        call assert_real(0., adjusted_rand_index([1,1,1,2,2,3], [1,1,1,1,1,1]), 1.e-6, 'a collapsed prediction scores 0')
    end subroutine test_ari

    !> the null distribution of the matched accuracy for populations of 50, 30 and 20 %, and the
    !! collapsed predictor; the null's median and 99.9th percentile (0.383 and 0.426) come from a
    !! Python emulation of 2000 shuffles
    subroutine test_null()
        integer :: truth(1000)
        integer, allocatable :: est(:)
        real :: med, p999
        write(*,'(A)') 'test_null'
        truth(1:500) = 1
        truth(501:800) = 2
        truth(801:1000) = 3
        est = truth
        call accuracy_null(truth, est, 3, 3, 1000, 20261017, med, p999)
        call assert_real(0.383, med, 0.01, 'the null median of the matched accuracy for 50/30/20 % populations')
        call assert_real(0.426, p999, 0.02, 'its 99.9th percentile')
        call assert_real(0.5, collapsed_accuracy(truth, 3), 1.e-6, 'the collapsed predictor scores the largest population')
        call assert_true(med < 0.5, 'the null sits below the collapsed predictor')
    end subroutine test_null

    !> the map mapping agrees with the label mapping, and fails on correct labels with swapped maps
    subroutine test_map_mapping()
        real :: cc(3,3)
        integer, allocatable :: maps(:), labels(:)
        write(*,'(A)') 'test_map_mapping'
        labels = match_states([1,1,2,2,3,3], [1,1,2,2,3,3], 3, 3)
        cc = reshape([0.9,0.5,0.4, 0.5,0.9,0.4, 0.4,0.5,0.9], [3,3])
        maps = match_maps(cc)
        call assert_true(all(maps == labels), 'maps that match their own truths agree with the labels')
        ! the maps of states 1 and 2 swapped: labels right, maps wrong
        cc = reshape([0.5,0.9,0.4, 0.9,0.5,0.4, 0.4,0.5,0.9], [3,3])
        maps = match_maps(cc)
        call assert_true(any(maps /= labels), 'swapped maps disagree with correct labels')
    end subroutine test_map_mapping

    !> the discrimination of three blob states: the truths themselves pass; copies of one map and
    !! their mean do not
    subroutine test_discrimination()
        use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
        type(image) :: truth(3), est(3)
        real    :: margins(3), diffcorr(3,3), body(3,3), lobes(3,3), corr(3,3)
        integer :: k, nmask
        write(*,'(A)') 'test_discrimination'
        body  = reshape([0.,0.,0., 8.,0.,0., -4.,6.,0.], [3,3])
        lobes = reshape([0.,10.,6., 10.,-8.,-4., -9.,-6.,5.], [3,3])
        do k = 1, 3
            call blob_volume(truth(k), reshape([body, lobes(:,k)], [3,4]), [5.,4.,4.,3.], [1.,.8,.8,.7])
        enddo
        do k = 1, 3
            est(k) = truth(k)
        enddo
        call discrimination(est, truth, 0.1, 1, margins, diffcorr, nmask, corr)
        call assert_true(nmask > 0, 'the truths differ somewhere')
        call assert_true(all([(abs(corr(k,k) - 1.) < 1.e-4, k = 1, 3)]), 'a true map correlates fully with its truth')
        call assert_true(abs(margins(2) - (corr(2,2) - max(corr(2,1), corr(2,3)))) < 1.e-6, &
            &'the margin is read off the masked correlations')
        ! a perfect map's margin is 1 less the masked correlation of two truths (about 0.18 here)
        call assert_true(all(margins > 0.1), 'each true map is closest to its own truth by a margin')
        call assert_true(diffcorr(1,2) > 0.99 .and. diffcorr(2,3) > 0.99, 'the differences of the true maps are the true differences')
        call assert_true(maps_distinct(margins, diffcorr), 'the true maps count as distinct states')
        call assert_true(size(offdiag(diffcorr)) == 6 .and. abs(finite_min(offdiag(diffcorr)) - minval(diffcorr, &
            &mask=reshape([(k /= 1, k = 1, 3), (k /= 2, k = 1, 3), (k /= 3, k = 1, 3)], [3,3]))) < 1.e-6, &
            &'the smallest off-diagonal difference correlation')
        ! copies of one map
        do k = 1, 3
            est(k) = truth(1)
        enddo
        call discrimination(est, truth, 0.1, 1, margins, diffcorr, nmask)
        call assert_true(.not. all(margins > 0.), 'copies of one map do not all discriminate')
        call assert_true(.not. ieee_is_finite(diffcorr(1,2)), 'the difference of copies has no correlation')
        call assert_true(.not. maps_distinct(margins, diffcorr), 'copies of one map do not count as distinct states')
        call assert_true(.not. ieee_is_finite(finite_min([1., diffcorr(1,2)])), 'a non-finite value makes the minimum non-finite')
        ! their mean
        call est(1)%zero_and_unflag_ft
        do k = 1, 3
            call est(1)%add(truth(k))
        enddo
        do k = 2, 3
            est(k) = est(1)
        enddo
        call discrimination(est, truth, 0.1, 1, margins, diffcorr, nmask)
        call assert_true(.not. all(margins > 0.), 'copies of the mean map do not all discriminate')
        call assert_true(.not. ieee_is_finite(diffcorr(2,3)), 'nor do their differences')
        do k = 1, 3
            call truth(k)%kill
            call est(k)%kill
        enddo
    end subroutine test_discrimination

    !> The frame and hand rule of carry_pose: a map t, a truth w made from it by a known rotation
    !! and shift (rotvol, as docking rotates), directly and from its x-mirror; every particle image
    !! of t at (r, sh) must equal the image of w at the carried pose and shift
    subroutine test_frame_rule()
        use simple_simple_volinterp, only: rotvol
        type(image) :: t, tm, w
        type(ori)   :: od
        real :: pos(3,4), a(3,3), s3(3), r(3,3), sh(2), rc(3,3), shc(2), wrong(3,3)
        integer :: i, ihand
        write(*,'(A)') 'test_frame_rule'
        pos = reshape([10.,4.,-6., -5.,9.,3., 2.,-8.,7., -9.,-3.,-4.], [3,4])
        call blob_volume(t, pos, [3.,2.5,2.,3.5], [1.,.8,.7,.5])
        a  = euler2m([20., 50., 75.])
        s3 = [1.5, -2., 1.]
        call od%new(.false.)
        call od%set_euler(m2euler(a))
        do ihand = 1, 2
            if( ihand == 1 )then
                w = rotvol(t, od, s3)
            else
                tm = t
                call tm%mirror('x')
                w = rotvol(tm, od, s3)
            endif
            do i = 1, 3
                r  = euler2m([37.*real(i), 63.+10.*real(i), 112.-20.*real(i)])
                sh = [1.-0.5*real(i), -1.5+real(i)]
                call carry_pose(r, sh, ihand == 2, B, a, s3, rc, shc)
                call assert_true(proj_corr(t, r, sh, w, rc, shc) > 0.995, &
                    &'an image keeps its pose and shift carried into the docked frame (hand '//int2str(ihand)//')')
            enddo
            ! the transposed docking rotation would not
            wrong = matmul(r, a)
            call assert_true(proj_corr(t, r, sh, w, wrong, shc) < 0.9, 'a wrong composition is detected')
        enddo
        call t%kill
        call tm%kill
        call w%kill
        call od%kill
    end subroutine test_frame_rule

    !> The pose error under a declared group: poses equal to the truth up to the group's operators
    !! (r g, the composition under which a symmetric map projects the same, pinned here by
    !! projection) score zero; a global frame error that a learnt cluster would accept does not
    subroutine test_group_error()
        type(image) :: vs
        real, allocatable :: ops(:,:,:), errs(:)
        real :: rt(3,3,12), re(3,3,12), w(3,3)
        integer :: j, g
        write(*,'(A)') 'test_group_error'
        ops = group_operators('d2')
        call assert_int(4, size(ops,3), 'D2 has four operators')
        call blob_volume(vs, reshape([10.,4.,-6., 10.,-4.,6., -10.,4.,6., -10.,-4.,-6.], [3,4]), [3.,3.,3.,3.], [1.,1.,1.,1.])
        do g = 2, 4
            call assert_true(proj_corr(vs, euler2m([37.,63.,112.]), [0.,0.], vs, matmul(euler2m([37.,63.,112.]), ops(:,:,g)), &
                &[0.,0.]) > 0.999, 'a D2 map projects the same at r and at r g (operator '//int2str(g)//')')
        enddo
        do j = 1, 12
            rt(:,:,j) = euler2m([30.*real(j), 15.*real(j), 360. - 25.*real(j)])
            re(:,:,j) = matmul(rt(:,:,j), ops(:,:,1 + mod(j, 4)))
        enddo
        errs = group_errors(re, rt, ops)
        ! single-precision acos near 1 resolves angles to about 0.03 degrees
        call assert_real(0., maxval(errs), 0.1, 'poses equal up to the group score zero')
        ! one global rotation of 30 degrees: every particle's frame difference is the same, which a
        ! cluster learnt from the data would take for a symmetry mate
        w = axis_operator([0.6, 0.8, 0.], 30.)
        do j = 1, 12
            re(:,:,j) = matmul(rt(:,:,j), w)
        enddo
        errs = group_errors(re, rt, ops)
        call assert_true(median(errs) > 25., 'a wrong frame is an error under the declared group')
        ! the in-plane partner Rz180 r g: an error (twice the angle between the beam axis in the
        ! particle's frame and the nearest two-fold, zero where they coincide), zero under projection inversion
        do j = 1, 12
            re(:,:,j) = matmul(reshape([-1.,0.,0., 0.,-1.,0., 0.,0.,1.], [3,3]), matmul(rt(:,:,j), ops(:,:,1 + mod(j, 4))))
        enddo
        errs = group_errors(re, rt, ops)
        call assert_true(median(errs) > 1., 'the in-plane partner is an error without projection inversion')
        errs = group_errors(re, rt, ops, projection_inversion=.true.)
        call assert_real(0., maxval(errs), 0.1, 'the in-plane partner scores zero under projection inversion')
        call vs%kill
    end subroutine test_group_error

    !> a non-finite metric fails its gate row
    subroutine test_nonfinite()
        use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan
        use simple_test_gate, only: test_gate
        type(test_gate) :: gate
        write(*,'(A)') 'test_nonfinite'
        call gate%new(string('tmp_truth_metrics_gate.tsv'))
        call gate%metric('nan', ieee_value(1., ieee_quiet_nan), 0., .true.)
        call assert_true(.not. gate%passed(), 'a non-finite metric fails, even with a passing verdict')
        call gate%kill
        call del_file('tmp_truth_metrics_gate.tsv')
    end subroutine test_nonfinite

end module simple_test_truth_metrics_tester
