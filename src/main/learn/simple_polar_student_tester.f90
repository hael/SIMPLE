!@descr: unit tests of simple_polar_student: the soft targets keep a row's tails, the weights resolve beside the header and must be the header's network, the core's header round trip keeps its owner's keys and reloads to the same prediction
submodule (simple_polar_student) simple_polar_student_tester
use simple_test_utils
use simple_syslib, only: simple_mkdir, simple_chdir
use simple_rnd,    only: gasdev, ran3
implicit none

integer, parameter :: PFTSZ = 12, KFROMTO(2) = [2, 9], NK = 8, BOX = 32, NCLASS = 6, WIDTH = 2
real,    parameter :: SMPD = 1.0, MSKDIAM = 24., LP = 4.

contains

    module subroutine run_all_polar_student_tests()
        write(*,'(A)') '**** running all polar_student tests ****'
        call test_targets()
        call test_weights_helpers()
        call test_header_roundtrip()
    end subroutine run_all_polar_student_tests

    !> a core of the tiny geometry with random positive shell weights and random templates
    subroutine make_core( core, seed )
        type(polar_student), intent(inout) :: core
        integer(int64),      intent(in)    :: seed
        real(sp), allocatable :: wk(:,:), xtmpl(:,:,:,:)
        integer :: k, i
        allocate(wk(KFROMTO(1):KFROMTO(2),2))
        do k = KFROMTO(1),KFROMTO(2)
            wk(k,1) = 0.5 + ran3()
            wk(k,2) = 0.5 + ran3()
        enddo
        call core%new(PFTSZ, KFROMTO, BOX, BOX, SMPD, MSKDIAM, LP, wk, .false., NCLASS, .true., WIDTH, seed)
        allocate(xtmpl(core%get_pftsz_net(), NK, 2, NCLASS))
        do i = 1,NCLASS
            call random_tensor(xtmpl(:,:,:,i))
        enddo
        call core%init_templates(xtmpl)
        deallocate(wk, xtmpl)
    end subroutine make_core

    subroutine random_tensor( x )
        real(sp), intent(out) :: x(:,:,:)
        integer :: p, k, part
        do part = 1,size(x,3)
            do k = 1,size(x,2)
                do p = 1,size(x,1)
                    x(p,k,part) = gasdev()
                enddo
            enddo
        enddo
    end subroutine random_tensor

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
        call polar_student_targets(N, dists, 1., t)
        call assert_real(1., sum(t), 1.e-5, 'the targets sum to one')
        call assert_int(1, maxloc(t, 1), 'the mode is the smallest loss')
        call assert_true(all(t([1,2,3,5,6,7,8,9,10]) > 0.), 'every finite entry keeps its mass: the tails stay')
        call assert_real(0., t(4), 0., 'an invalid entry gets no mass')
        call assert_real(exp(-9.) / exp(-1.), t(10) / t(2), 1.e-5, 'the masses follow exp(-temperature times the loss above the minimum)')
        call polar_student_targets(N, dists, 0., t)
        call assert_real(1., t(1), 0., 'a zero temperature gives the mode alone')
        call assert_real(0., sum(t(2:N)), 0., 'and nothing else')
    end subroutine test_targets

    !> A weights name resolves in the header's directory; a name with a directory, or none, is refused; a network
    !! of another seed, width or weight count than the header's is refused
    subroutine test_weights_helpers()
        type(polar_student) :: core
        type(hhcnn)         :: net
        type(string)        :: cwd_saved, fixture_root, wfile
        character(len=:), allocatable :: msg
        integer :: nfail0
        write(*,'(A)') 'test_weights_helpers'
        nfail0 = tests_failed
        call enter_fixture('polar_student_weights', cwd_saved, fixture_root)
        call make_core(core, 20261122_int64)
        call core%write(string('model.txt'), 'test model', 1)
        call header_weights_path(string('dirC/model.txt'), 'model_weights.bin', wfile, msg)
        call assert_true(len(msg) == 0 .and. wfile%to_char() == 'dirC/model_weights.bin', 'a weights name resolves in the header''s directory')
        call header_weights_path(string('dirC/model.txt'), 'dirB/model_weights.bin', wfile, msg)
        call assert_true(len(msg) > 0, 'a weights name with a directory is refused')
        call header_weights_path(string('dirC/model.txt'), '', wfile, msg)
        call assert_true(len(msg) > 0, 'a header without a weights name is refused')
        call net%read_weights(string('model_weights.bin'))
        msg = weights_header_mismatch(net, core%pftsz_net, NK, NCLASS, WIDTH, 20261122_int64, core%nweights)
        call assert_true(len(msg) == 0, 'the header''s own network passes the check')
        msg = weights_header_mismatch(net, core%pftsz_net, NK, NCLASS, WIDTH, 20261121_int64, core%nweights)
        call assert_true(len(msg) > 0, 'a network of another seed than the header''s is refused')
        msg = weights_header_mismatch(net, core%pftsz_net, NK, NCLASS, 2*WIDTH, 20261122_int64, core%nweights)
        call assert_true(len(msg) > 0, 'a network of another width than the header''s is refused')
        msg = weights_header_mismatch(net, core%pftsz_net, NK, NCLASS, WIDTH, 20261122_int64, core%nweights + 1)
        call assert_true(len(msg) > 0, 'another weight count than the header''s is refused')
        call leave_fixture(cwd_saved, fixture_root, nfail0)
        call net%kill
        call core%kill
        call wfile%kill
    end subroutine test_weights_helpers

    !> The core's model file: the owner's keys come back in order, the core's geometry and records travel, and a core
    !! rebuilt from the header and loaded predicts bit for bit as the written one
    subroutine test_header_roundtrip()
        type(polar_student)        :: core, back
        type(polar_student_header) :: hdr
        type(string)               :: cwd_saved, fixture_root
        real(sp), allocatable :: x(:,:,:)
        real(sp) :: ctf2(NK), p0(NCLASS), p1(NCLASS), t, cf, hs, hr
        integer  :: nfail0
        write(*,'(A)') 'test_header_roundtrip'
        nfail0 = tests_failed
        call enter_fixture('polar_student_header', cwd_saved, fixture_root)
        call make_core(core, 20261126_int64)
        call core%set_threshold(0.9900049, 0.95, 0.5, 1.)
        call core%set_provenance('teacher.simple', 12345_int64, 'every_10th_selected_stack')
        call core%write(string('model.txt'), 'test model', 7, ['first ', 'second'], ['one  ', 'two 2'])
        call hdr%read(string('model.txt'))
        call assert_int(7, hdr%version, 'the version travels')
        call assert_int(2, hdr%nextra, 'the owner''s keys are kept')
        call assert_true(hdr%keys(1)%to_char() == 'first' .and. hdr%keys(2)%to_char() == 'second', 'in their order')
        call assert_true(hdr%vals(2)%to_char() == 'two 2', 'with their values')
        call assert_int(PFTSZ, hdr%pftsz, 'the calculator''s half ring travels')
        call assert_true(all(hdr%kfromto == KFROMTO), 'the band travels')
        call back%new(hdr%pftsz, hdr%kfromto, hdr%box, hdr%box_crop, hdr%smpd, hdr%mskdiam, hdr%lp, hdr%wk, &
            &hdr%with_ctf, NCLASS, .true., hdr%width, hdr%seed)
        call back%load(hdr, string('model.txt'))
        call back%get_threshold(t, cf, hs, hr)
        call assert_true(transfer(t, 0) == transfer(0.9900049, 0), 'the threshold reads back bit for bit')
        allocate(x(core%get_pftsz_net(), NK, 2))
        call random_tensor(x)
        ctf2 = 1.
        call core%prepare_predict(1)
        call back%prepare_predict(1)
        call core%predict_tensor(1, x, ctf2, 0., p0)
        call back%predict_tensor(1, x, ctf2, 0., p1)
        call assert_true(all(p0 == p1), 'the reloaded core predicts as the written one, bit for bit')
        call assert_real(1., sum(p0), 1.e-4, 'the class posterior sums to one')
        call leave_fixture(cwd_saved, fixture_root, nfail0)
        call hdr%kill
        call core%kill
        call back%kill
        deallocate(x)
    end subroutine test_header_roundtrip

end submodule simple_polar_student_tester
