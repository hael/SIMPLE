!@descr: unit tests of simple_pose_student: the soft targets keep a row's tails, the network classifies the mirror pairs, the cone on the folded sphere (a direction and its mirror partner one class) caps at the dense search's marker, the model file round trip is exact
module simple_pose_student_tester
use simple_test_utils
use, intrinsic :: iso_fortran_env, only: int64
use simple_defs,         only: sp
use simple_string,       only: string
use simple_polar_feats,  only: polar_feats, polar_feat
use simple_pose_student, only: pose_student, student_pred, pose_student_targets, KMAX_STUDENT
use simple_rnd,          only: gasdev, ran3
implicit none
private
public :: run_all_pose_student_tests

integer, parameter :: PFTSZ = 12, KFROMTO(2) = [2, 9], NK = 8, BOX = 32, NSPACE = 100, WIDTH = 2
real,    parameter :: SMPD = 1.0, MSKDIAM = 24., LP = 4.

contains

    subroutine run_all_pose_student_tests()
        write(*,'(A)') '**** running all pose_student tests ****'
        call test_targets()
        call test_pairs()
        call test_cone()
        call test_file_roundtrip()
    end subroutine run_all_pose_student_tests

    !> a student of the tiny geometry with random positive shell weights
    subroutine make_student( student, seed )
        type(pose_student), intent(inout) :: student
        integer(int64),     intent(in)    :: seed
        real(sp) :: wk(KFROMTO(1):KFROMTO(2),2)
        integer  :: k
        do k = KFROMTO(1),KFROMTO(2)
            wk(k,1) = 0.5 + ran3()
            wk(k,2) = 0.5 + ran3()
        enddo
        call student%new(PFTSZ, KFROMTO, BOX, BOX, SMPD, MSKDIAM, 'c1', NSPACE, LP, wk, .false., WIDTH, seed)
    end subroutine make_student

    !> The soft targets of a row: every finite entry keeps its mass (the tails are not cut), the mode is the
    !! smallest loss, invalid entries get none, and a zero temperature gives the mode alone
    subroutine test_targets()
        integer, parameter :: N = 10
        real(sp) :: dists(N), t(N)
        integer  :: i
        write(*,'(A)') 'test_targets'
        do i = 1,N
            dists(i) = real(i - 1)
        enddo
        dists(4) = huge(1.)
        call pose_student_targets(N, dists, 1., t)
        call assert_real(1., sum(t), 1.e-5, 'the targets sum to one')
        call assert_int(1, maxloc(t, 1), 'the mode is the smallest loss')
        call assert_true(all(t([1,2,3,5,6,7,8,9,10]) > 0.), 'every finite entry keeps its mass: the tails stay')
        call assert_real(0., t(4), 0., 'an invalid entry gets no mass')
        call assert_real(exp(-9.) / exp(-1.), t(10) / t(2), 1.e-5, 'the masses follow exp(-temperature times the loss above the minimum)')
        call pose_student_targets(N, dists, 0., t)
        call assert_real(1., t(1), 0., 'a zero temperature gives the mode alone')
        call assert_real(0., sum(t(2:N)), 0., 'and nothing else')
    end subroutine test_targets

    !> The network classifies the grid's mirror pairs: a prediction puts all of its mass on the pairs' unmirrored
    !! directions (the spiral's first half), none on their mirrors, which the search tests
    subroutine test_pairs()
        type(pose_student) :: student
        real(sp) :: x(PFTSZ, NK, 2), ctf2(NK), prob(NSPACE)
        integer  :: p, k, part, i
        logical  :: l_halves
        write(*,'(A)') 'test_pairs'
        call make_student(student, 20261109_int64)
        do part = 1,2
            do k = 1,NK
                do p = 1,PFTSZ
                    x(p,k,part) = gasdev()
                enddo
            enddo
        enddo
        ctf2 = 1.
        call student%prepare_predict(1)
        call student%predict_tensor(1, x, ctf2, 0., prob)
        call assert_real(1., sum(prob), 1.e-5, 'the posterior over the grid sums to one')
        call assert_true(all(prob(NSPACE/2+1:) == 0.), 'the mirrored half of the grid carries no mass')
        call assert_true(all(prob(1:NSPACE/2) > 0.), 'every pair''s unmirrored direction carries its class''s mass')
        l_halves = .true.
        do i = 1,NSPACE/2
            if( student%mirror_of(i) <= NSPACE/2 ) l_halves = .false.
        enddo
        call assert_true(l_halves, 'each class is one direction of the first half and its mirror in the second')
        call student%kill
    end subroutine test_pairs

    !> The cone on the folded sphere (ruling R1): a flat posterior folds to every mirror pair, its radius spanning the
    !! folded sphere; a peaked one gives the mode alone at full confidence and zero radius; two neighbours give a cone
    !! of two whose radius is their angle on the folded sphere; a direction and its mirror partner are one class, the
    !! heavier partner named, with the full mass and zero radius
    subroutine test_cone()
        type(pose_student) :: student
        type(student_pred) :: pred
        real(sp) :: prob(NSPACE), euls(3), ang
        integer  :: j, m7, m30
        write(*,'(A)') 'test_cone'
        call make_student(student, 20261113_int64)
        m7  = student%mirror_of(7)
        m30 = student%mirror_of(30)
        call assert_true(m7 /= 7 .and. student%mirror_of(m7) == 7, 'the mirror partner is an involution')
        call assert_int(min(7, m7), student%fold_index(7), 'the folded class is the smaller index of the pair')
        call assert_int(student%fold_index(7), student%fold_index(m7), 'shared by both partners')
        prob = 1. / real(NSPACE)
        call student%cone_of(prob, pred)
        call assert_int(NSPACE / 2, pred%k99, 'a flat posterior needs every folded pair')
        call assert_int(min(NSPACE / 2, KMAX_STUDENT), pred%ncone, 'and the cone is the cap or every pair')
        call assert_true(pred%radius > 45., 'with a radius spanning the folded sphere')
        prob    = 0.0001
        prob(7) = 1. - 0.0001 * real(NSPACE - 1)
        call student%cone_of(prob, pred)
        call assert_int(7, pred%mode, 'the mode of a peaked posterior')
        call assert_int(1, pred%k99, 'one folded direction holds 99%')
        call assert_int(1, pred%ncone, 'the cone is the mode alone')
        call assert_true(pred%conf >= 0.99, 'at full confidence')
        call assert_real(0., pred%radius, 0., 'with zero radius')
        prob     = 0.
        prob(7)  = 0.6
        prob(30) = 0.4
        call student%cone_of(prob, pred)
        call assert_int(2, pred%ncone, 'two directions sharing the mass give a cone of two')
        call assert_true(pred%cone(1) == 7 .and. pred%cone(2) == 30, 'ordered by mass')
        euls = student%grid_euler(7)
        ang  = min(student%angle_between(30, euls), student%angle_between(m30, euls))
        call assert_real(ang, pred%radius, 1.e-4, 'whose radius is their angle on the folded sphere')
        do j = 3,KMAX_STUDENT
            if( pred%cone(j) /= 0 ) call assert_true(.false., 'the cone beyond its size is empty')
        enddo
        prob     = 0.
        prob(7)  = 0.4
        prob(m7) = 0.6
        call student%cone_of(prob, pred)
        call assert_int(1, pred%k99, 'a direction and its mirror partner are one folded class')
        call assert_int(m7, pred%mode, 'named by the heavier partner')
        call assert_real(1., pred%conf, 1.e-6, 'holding their summed mass')
        call assert_real(0., pred%radius, 0., 'at zero radius')
        call student%kill
    end subroutine test_cone

    !> The model file: the header and the weights written and read back give the same geometry, grid, weights
    !! and, bit for bit, the same prediction
    subroutine test_file_roundtrip()
        type(pose_student) :: student, back
        type(polar_feats)  :: feats
        type(polar_feat)   :: feat
        type(string)       :: cwd_saved, fixture_root, fname
        complex(sp) :: pft(PFTSZ, KFROMTO(1):KFROMTO(2))
        real(sp)    :: ctfmat(PFTSZ, KFROMTO(1):KFROMTO(2)), prob0(NSPACE), prob1(NSPACE)
        real(sp), allocatable :: wk0(:,:), wk1(:,:)
        integer  :: nfail0, p, k, kf(2)
        write(*,'(A)') 'test_file_roundtrip'
        nfail0 = tests_failed
        call enter_fixture('pose_student_file', cwd_saved, fixture_root)
        call make_student(student, 20261114_int64)
        fname = 'student.txt'
        call student%write(fname)
        call back%read(fname)
        call assert_int(student%get_nspace(), back%get_nspace(), 'the grid size travels')
        call assert_int(student%get_pftsz(), back%get_pftsz(), 'the half ring travels')
        kf = back%get_kfromto()
        call assert_true(all(kf == KFROMTO), 'the band travels')
        call assert_real(student%get_lp(), back%get_lp(), 0., 'the low-pass travels')
        call assert_int(student%get_box(), back%get_box(), 'the box travels')
        call assert_true(student%get_with_ctf() .eqv. back%get_with_ctf(), 'the CTF treatment travels')
        call student%get_weights(wk0)
        call back%get_weights(wk1)
        call assert_true(all(wk0 == wk1), 'the shell weights travel bit for bit')
        ! the same prediction on a random ring
        do k = KFROMTO(1),KFROMTO(2)
            do p = 1,PFTSZ
                pft(p,k) = cmplx(gasdev(), gasdev())
            enddo
        enddo
        ctfmat = 1.
        call feats%new(PFTSZ, KFROMTO, BOX, wk0)
        call feats%build_from_pft(pft, ctfmat, .true., feat)
        call student%prepare_predict(1)
        call back%prepare_predict(1)
        call student%predict_feat(1, feat, prob0)
        call back%predict_feat(1, feat, prob1)
        call assert_true(all(prob0 == prob1), 'the read model predicts as the written one, bit for bit')
        call assert_real(1., sum(prob0), 1.e-4, 'the posterior sums to one')
        call fname%kill
        call leave_fixture(cwd_saved, fixture_root, nfail0)
        call feats%kill
        call feat%kill
        call student%kill
        call back%kill
        deallocate(wk0, wk1)
    end subroutine test_file_roundtrip

end module simple_pose_student_tester
