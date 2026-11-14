!@descr: unit tests of simple_class_student: the candidates hold the posterior's mass in descending order and cap at KMAX_CLASS, the model file round trip keeps the classes, their states and the prediction bit for bit, the read-out angle lies on the circle
submodule (simple_class_student) simple_class_student_tester
use simple_test_utils
use simple_rnd, only: gasdev, ran3
implicit none

integer, parameter :: PFTSZ = 12, KFROMTO(2) = [2, 9], NK = 8, BOX = 32, NCLS = 20, WIDTH = 2
real,    parameter :: SMPD = 1.0, MSKDIAM = 24., LP = 4.

contains

    module subroutine run_all_class_student_tests()
        write(*,'(A)') '**** running all class_student tests ****'
        call test_cands()
        call test_file_roundtrip()
    end subroutine run_all_class_student_tests

    subroutine make_student( student, seed )
        type(class_student), intent(inout) :: student
        integer(int64),      intent(in)    :: seed
        real(sp), allocatable :: wk(:,:), xc(:,:,:,:)
        integer :: k, states(NCLS), i, p, part
        allocate(wk(KFROMTO(1):KFROMTO(2),2))
        do k = KFROMTO(1),KFROMTO(2)
            wk(k,1) = 0.5 + ran3()
            wk(k,2) = 0.5 + ran3()
        enddo
        states = 1
        states(3) = 0
        call student%new(PFTSZ, KFROMTO, BOX, BOX, SMPD, MSKDIAM, LP, wk, .false., NCLS, states, WIDTH, seed)
        allocate(xc(PFTSZ, NK, 2, NCLS))
        do i = 1,NCLS
            do part = 1,2
                do k = 1,NK
                    do p = 1,PFTSZ
                        xc(p,k,part,i) = gasdev()
                    enddo
                enddo
            enddo
        enddo
        call student%init_readout(xc)
        deallocate(wk, xc)
    end subroutine make_student

    !> The candidates: the mode first, then the classes in descending mass until CAND_MASS, capped at KMAX_CLASS; a
    !! peaked posterior gives few, a flat one the cap
    subroutine test_cands()
        type(class_student) :: student
        type(class_pred)    :: pred
        real(sp) :: prob(NCLS)
        write(*,'(A)') 'test_cands'
        call make_student(student, 20261113_int64)
        prob = 0.0005
        prob(7) = 0.9
        prob(2) = 0.09
        prob = prob / sum(prob)
        call student%cands_of(prob, pred)
        call assert_int(7, pred%mode, 'the mode is the class of largest mass')
        call assert_true(pred%ncand >= 2 .and. pred%cands(1) == 7 .and. pred%cands(2) == 2, 'the candidates in descending mass')
        call assert_true(pred%conf >= CAND_MASS - 1.e-5, 'the candidates hold the declared mass')
        prob = 1. / real(NCLS)
        call student%cands_of(prob, pred)
        call assert_int(KMAX_CLASS, pred%ncand, 'a flat posterior gives the cap')
        call assert_true(pred%k99 > KMAX_CLASS, 'and says that its mass needs more classes')
        call student%kill
    end subroutine test_cands

    !> The model file: the classes, their states and the class averages' hash travel, and the read model predicts as
    !! the written one, bit for bit; the read-out angle of every class lies in [0, 360)
    subroutine test_file_roundtrip()
        type(class_student) :: student, back
        type(string)        :: cwd_saved, fixture_root
        real(sp) :: x(PFTSZ, NK, 2), ctf2(NK), p0(NCLS), p1(NCLS), e3
        integer  :: nfail0, k, p, part, ic
        logical  :: l_angles
        write(*,'(A)') 'test_file_roundtrip'
        nfail0 = tests_failed
        call enter_fixture('class_student_file', cwd_saved, fixture_root)
        call make_student(student, 20261114_int64)
        call student%set_provenance('teacher.simple', 7_int64, 4242_int64, 'every_10th_selected_stack')
        call student%write(string('student.txt'))
        call back%read(string('student.txt'))
        call assert_int(NCLS, back%get_ncls(), 'the classes travel')
        call assert_int(0, back%get_state(3), 'a rejected class stays rejected')
        call assert_int(1, back%get_state(4), 'an accepted class stays accepted')
        call assert_true(back%cavgs_hash == 4242_int64, 'the class averages'' hash travels')
        do part = 1,2
            do k = 1,NK
                do p = 1,PFTSZ
                    x(p,k,part) = gasdev()
                enddo
            enddo
        enddo
        ctf2 = 1.
        call student%prepare_predict(1)
        call back%prepare_predict(1)
        call student%predict_tensor(1, x, ctf2, 0., p0)
        call back%predict_tensor(1, x, ctf2, 0., p1)
        call assert_true(all(p0 == p1), 'the read model predicts as the written one, bit for bit')
        call assert_real(1., sum(p0), 1.e-4, 'the posterior over the classes sums to one')
        l_angles = .true.
        do ic = 1,NCLS
            e3 = back%rotation_of(1, ic)
            if( e3 < 0. .or. e3 >= 360. ) l_angles = .false.
        enddo
        call assert_true(l_angles, 'every class''s read-out angle lies in [0, 360)')
        call leave_fixture(cwd_saved, fixture_root, nfail0)
        call student%kill
        call back%kill
    end subroutine test_file_roundtrip

end submodule simple_class_student_tester
