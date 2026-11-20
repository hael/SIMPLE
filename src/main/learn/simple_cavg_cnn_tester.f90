!@descr: unit tests of simple_cavg_cnn: its convolutions against a brute force of PyTorch's formula, GroupNorm, the symmetries of the square, the forward pass against PyTorch's logits on the committed fixture with either product, the symmetry average's invariance, predict's probabilities and thread identity, and the weights reader's refusals
submodule (simple_cavg_cnn) simple_cavg_cnn_tester
use simple_test_utils
use simple_defs,         only: SIMPLE_SOURCE_PATH
use simple_rnd,          only: gasdev
use simple_syslib,       only: del_file
use simple_string_utils, only: int2str
!$ use omp_lib,          only: omp_set_num_threads, omp_get_max_threads, omp_get_num_threads, omp_get_thread_limit
implicit none
#include "simple_local_flags.inc"

character(len=*), parameter :: FIXTURE_DIR = SIMPLE_SOURCE_PATH//'/testdata/cavg_cnn/'
! ten times the largest float32-float64 difference of PyTorch's logits on the fixture, as logits.txt records
real(sp),         parameter :: T1_LOGITS   = 3.361e-6
! a convolution output sums at most 9 cin products of order-one values in single precision
real(sp),         parameter :: TOL_CONV    = 1.e-5
integer,          parameter :: NFIX        = 3       !< the fixture's inputs

contains

    module subroutine run_all_cavg_cnn_tests()
        write(*,'(A)') '**** running all cavg_cnn tests ****'
        call test_conv_brute_force(.true.)
        call test_conv_brute_force(.false.)
        call test_conv_delta()
        call test_group_norm()
        call test_symmetries()
        call test_fixture_parity(.true.)
        call test_fixture_parity(.false.)
        call test_symmetry_average_invariance()
        call test_predict()
        call test_reader_refusals()
    end subroutine run_all_cavg_cnn_tests

    ! ===== FIXTURE =====

    !> the fixture's inputs (64, 64 and 96 pixels, 3 channels) and PyTorch's logits: plain, averaged over the symmetries
    subroutine read_fixture( x1, x2, x3, lplain, lsym )
        real(sp), allocatable, intent(out) :: x1(:,:,:), x2(:,:,:), x3(:,:,:)
        real(sp),              intent(out) :: lplain(2,NFIX), lsym(2,NFIX)
        character(len=512) :: line
        integer  :: funit, ios, nin, nx, ny, nc, i, idx
        real(sp) :: v(4)
        open(newunit=funit, file=FIXTURE_DIR//'inputs.raw', access='stream', form='unformatted', status='old', iostat=ios)
        if( ios /= 0 ) THROW_HARD('missing fixture testdata/cavg_cnn/inputs.raw')
        read(funit) nin
        if( nin /= NFIX ) THROW_HARD('the fixture holds another input count')
        do i = 1,NFIX
            read(funit) nx, ny, nc
            select case(i)
                case(1)
                    allocate(x1(nx,ny,nc))
                    read(funit) x1
                case(2)
                    allocate(x2(nx,ny,nc))
                    read(funit) x2
                case(3)
                    allocate(x3(nx,ny,nc))
                    read(funit) x3
            end select
        enddo
        close(funit)
        open(newunit=funit, file=FIXTURE_DIR//'logits.txt', status='old', action='read', iostat=ios)
        if( ios /= 0 ) THROW_HARD('missing fixture testdata/cavg_cnn/logits.txt')
        i = 0
        do
            read(funit, '(A)', iostat=ios) line
            if( ios /= 0 ) exit
            if( line(1:1) == '#' .or. len_trim(line) == 0 ) cycle
            read(line, *) idx, v
            i = i + 1
            lplain(:,idx) = v(1:2)
            lsym(:,idx)   = v(3:4)
        enddo
        close(funit)
        if( i /= NFIX ) THROW_HARD('logits.txt holds another input count')
    end subroutine read_fixture

    ! ===== THE LAYERS =====

    !> the 3x3 convolution at strides 1 and 2 and the 1x1 at stride 2 against PyTorch's formula evaluated directly
    subroutine test_conv_brute_force( l_blas )
        logical, intent(in) :: l_blas
        integer, parameter :: N = 9, CIN = 2, COUT = 3
        real(sp) :: x(N*N, CIN), w(9*CIN, COUT), w1(CIN, COUT)
        real(sp), allocatable :: y(:,:)
        real(sp) :: ref, err
        integer  :: stride, ny, co, oh0, ow0, ci, kh, kw, ih0, iw0, i
        write(*,'(A,L1)') 'test_conv_brute_force, sgemm ', l_blas
        call set_fixed_seed(20261120)
        do ci = 1,CIN
            do i = 1,N*N
                x(i,ci) = gasdev()
            enddo
        enddo
        do co = 1,COUT
            do i = 1,9*CIN
                w(i,co) = gasdev()
            enddo
            do ci = 1,CIN
                w1(ci,co) = gasdev()
            enddo
        enddo
        do stride = 1,2
            call conv3x3(x, N, CIN, w, COUT, stride, l_blas, y, ny)
            call assert_int((N - 1) / stride + 1, ny, 'the output size is PyTorch''s floor((n + 2 - 3) / stride) + 1')
            err = 0.
            do co = 1,COUT
                do oh0 = 0,ny-1
                    do ow0 = 0,ny-1
                        ref = 0.
                        do ci = 0,CIN-1
                            do kh = 0,2
                                do kw = 0,2
                                    ih0 = oh0 * stride - 1 + kh
                                    iw0 = ow0 * stride - 1 + kw
                                    if( ih0 < 0 .or. ih0 >= N .or. iw0 < 0 .or. iw0 >= N ) cycle
                                    ref = ref + w(kw + 3*kh + 9*ci + 1, co) * x(ih0*N + iw0 + 1, ci + 1)
                                enddo
                            enddo
                        enddo
                        err = max(err, abs(ref - y(oh0*ny + ow0 + 1, co)))
                    enddo
                enddo
            enddo
            call assert_true(err <= TOL_CONV, 'the 3x3 convolution at stride '//int2str(stride)//&
                &' is PyTorch''s cross-correlation with zero padding 1, taps kw fastest, then kh, then the input channel')
        enddo
        call conv1x1_s2(x, N, CIN, w1, COUT, l_blas, y, ny)
        err = 0.
        do co = 1,COUT
            do oh0 = 0,ny-1
                do ow0 = 0,ny-1
                    ref = 0.
                    do ci = 0,CIN-1
                        ref = ref + w1(ci + 1, co) * x(2*oh0*N + 2*ow0 + 1, ci + 1)
                    enddo
                    err = max(err, abs(ref - y(oh0*ny + ow0 + 1, co)))
                enddo
            enddo
        enddo
        call assert_true(err <= TOL_CONV, 'the 1x1 convolution at stride 2 samples the even pixels (0-based) of the input')
    end subroutine test_conv_brute_force

    !> a delta at the centre returns the flipped kernel around it, exactly
    subroutine test_conv_delta()
        integer, parameter :: N = 7, D = 4
        real(sp) :: x(N*N,1), w(9,1), expected(N*N)
        real(sp), allocatable :: y(:,:)
        integer  :: ny, oh, ow
        write(*,'(A)') 'test_conv_delta'
        x = 0.
        x((D-1)*N + D, 1) = 1.
        w(:,1) = [(real(oh), oh = 1,9)]
        expected = 0.
        do oh = D-1,D+1
            do ow = D-1,D+1
                expected((oh-1)*N + ow) = w((D + 1 - ow) + 3*(D + 1 - oh) + 1, 1)
            enddo
        enddo
        call conv3x3(x, N, 1, w, 1, 1, .false., y, ny)
        call assert_true(all(y(:,1) == expected), 'a delta convolved at stride 1 is the kernel flipped in both axes')
    end subroutine test_conv_delta

    !> constant channels normalise to the bias; with unit weight and zero bias each group has mean 0 and variance 1
    subroutine test_group_norm()
        integer, parameter :: NP = 50, NC = 16
        real(sp) :: a(NP,NC), g(NC), b(NC)
        real(dp) :: mean, var
        integer  :: c, i, grp, cg
        write(*,'(A)') 'test_group_norm'
        call set_fixed_seed(20261121)
        cg = NC / CAVG_CNN_GROUPS
        do c = 1,NC
            a(:,c) = real((c - 1) / cg + 1)
            g(c)   = 1. + gasdev()
            b(c)   = gasdev()
        enddo
        call group_norm(a, g, b)
        call assert_true(maxval(abs(a - spread(b, 1, NP))) <= 1.e-6, 'a constant group normalises to the bias')
        do c = 1,NC
            do i = 1,NP
                a(i,c) = 3. + 2. * gasdev()
            enddo
        enddo
        g = 1.
        b = 0.
        call group_norm(a, g, b)
        do grp = 1,CAVG_CNN_GROUPS
            mean = sum(real(a(:,(grp-1)*cg+1:grp*cg), dp)) / real(NP*cg, dp)
            var  = sum((real(a(:,(grp-1)*cg+1:grp*cg), dp) - mean)**2) / real(NP*cg, dp)
            call assert_true(abs(mean) <= 1.d-5 .and. abs(var - 1.d0) <= 1.d-4, &
                &'GroupNorm leaves each group of '//int2str(cg)//' channels at mean 0 and variance 1 (biased)')
        enddo
    end subroutine test_group_norm

    !> the 8 transforms are distinct, include the identity, form a set closed under transposition, and the quarter
    !! turn applied four times is the identity
    subroutine test_symmetries()
        integer, parameter :: N = 5
        real(sp) :: x(N,N,1), y(N,N,1,0:7), t(N,N,1), r(N,N,1)
        integer  :: i, j, s, s2, nfound
        logical  :: distinct
        write(*,'(A)') 'test_symmetries'
        do j = 1,N
            do i = 1,N
                x(i,j,1) = real(i + 10*j)
            enddo
        enddo
        do s = 0,7
            call symmetry(x, s, y(:,:,:,s))
        enddo
        call assert_true(all(y(:,:,:,0) == x), 'symmetry 0 is the identity')
        distinct = .true.
        do s = 0,7
            do s2 = s+1,7
                if( all(y(:,:,:,s) == y(:,:,:,s2)) ) distinct = .false.
            enddo
        enddo
        call assert_true(distinct, 'the 8 symmetries give 8 different images of an asymmetric input')
        nfound = 0
        do s = 0,7
            t(:,:,1) = transpose(y(:,:,1,s))
            do s2 = 0,7
                if( all(t == y(:,:,:,s2)) )then
                    nfound = nfound + 1
                    exit
                endif
            enddo
        enddo
        call assert_int(8, nfound, 'the transpose of every symmetry is one of the 8: the set is the whole group')
        r = x
        do s = 1,4
            call symmetry(r, 1, t)
            r = t
        enddo
        call assert_true(all(r == x), 'four quarter turns are the identity')
    end subroutine test_symmetries

    ! ===== THE NETWORK =====

    !> PyTorch's logits of the committed fixture, plain and averaged over the symmetries, with either product
    subroutine test_fixture_parity( l_blas )
        logical, intent(in) :: l_blas
        type(cavg_cnn) :: net
        real(sp), allocatable :: x1(:,:,:), x2(:,:,:), x3(:,:,:)
        real(sp) :: lplain(2,NFIX), lsym(2,NFIX), l(2), dplain, dsym
        integer  :: i
        write(*,'(A,L1)') 'test_fixture_parity, sgemm ', l_blas
        call read_fixture(x1, x2, x3, lplain, lsym)
        call net%new(string(FIXTURE_DIR//'weights.raw'))
        net%l_blas = l_blas
        dplain = 0.
        dsym   = 0.
        do i = 1,NFIX
            select case(i)
                case(1)
                    call net%forward(x1, l)
                case(2)
                    call net%forward(x2, l)
                case(3)
                    call net%forward(x3, l)
            end select
            dplain = max(dplain, maxval(abs(l - lplain(:,i))))
            select case(i)
                case(1)
                    call net%logits_sym(x1, l)
                case(2)
                    call net%logits_sym(x2, l)
                case(3)
                    call net%logits_sym(x3, l)
            end select
            dsym = max(dsym, maxval(abs(l - lsym(:,i))))
        enddo
        write(*,'(A,ES10.3,A,ES10.3)') '    largest |logit - PyTorch|: plain ', dplain, ', symmetry-averaged ', dsym
        call assert_true(dplain <= T1_LOGITS, 'the forward pass gives PyTorch''s logits of the fixture within T1')
        call assert_true(dsym   <= T1_LOGITS, 'the symmetry-averaged logits are PyTorch''s within T1')
        call net%kill
    end subroutine test_fixture_parity

    !> the averaged logits do not depend on the input's orientation; the plain logits do (a transposed input misses
    !! PyTorch's by far more than T1, so the parity above checks the orientation)
    subroutine test_symmetry_average_invariance()
        type(cavg_cnn) :: net
        real(sp), allocatable :: x1(:,:,:), x2(:,:,:), x3(:,:,:), xt(:,:,:)
        real(sp) :: lplain(2,NFIX), lsym(2,NFIX), l(2)
        integer  :: c
        write(*,'(A)') 'test_symmetry_average_invariance'
        call read_fixture(x1, x2, x3, lplain, lsym)
        call net%new(string(FIXTURE_DIR//'weights.raw'))
        allocate(xt, mold=x3)
        do c = 1,size(x3,3)
            xt(:,:,c) = transpose(x3(:,:,c))
        enddo
        call net%logits_sym(xt, l)
        call assert_true(maxval(abs(l - lsym(:,3))) <= T1_LOGITS, 'the symmetry-averaged logits of a transposed input are unchanged')
        call net%forward(xt, l)
        call assert_true(maxval(abs(l - lplain(:,3))) > 100. * T1_LOGITS, 'the plain logits of a transposed input differ')
        call net%kill
    end subroutine test_symmetry_average_invariance

    !> predict's probabilities follow from the averaged logits (P(protein) = sigmoid(l1), P(ss | protein) =
    !! sigmoid(l2)), and are the same bit for bit at one thread and in a team of three
    subroutine test_predict()
        integer, parameter :: NTEAM = 3
        type(cavg_cnn) :: net
        real(sp), allocatable :: x1(:,:,:), x2(:,:,:), x3(:,:,:), xb(:,:,:,:)
        real(sp) :: lplain(2,NFIX), lsym(2,NFIX)
        real(dp) :: p1(3,2), pn(3,2), pp, ps, expected(3)
        integer  :: i, nthr_prev, limit, nteam_run
        write(*,'(A)') 'test_predict'
        call read_fixture(x1, x2, x3, lplain, lsym)
        call net%new(string(FIXTURE_DIR//'weights.raw'))
        allocate(xb(size(x1,1), size(x1,2), size(x1,3), 2))
        xb(:,:,:,1) = x1
        xb(:,:,:,2) = x2
        nthr_prev = 1
        !$ nthr_prev = omp_get_max_threads()
        limit = NTEAM
        !$ limit = omp_get_thread_limit()
        !$ call omp_set_num_threads(1)
        call net%predict(xb, p1)
        do i = 1,2
            pp = 1.d0 / (1.d0 + exp(-real(lsym(1,i), dp)))
            ps = 1.d0 / (1.d0 + exp(-real(lsym(2,i), dp)))
            expected = [1.d0 - pp, pp * ps, pp * (1.d0 - ps)]
            call assert_true(maxval(abs(p1(:,i) - expected)) <= real(T1_LOGITS, dp), &
                &'predict gives junk 1 - P(protein), ss P(protein) P(ss | protein), likely the rest')
            call assert_true(abs(sum(p1(:,i)) - 1.d0) <= 1.d-12, 'the three probabilities sum to 1')
        enddo
        if( limit >= NTEAM )then
            !$ call omp_set_num_threads(NTEAM)
            nteam_run = NTEAM
            !$omp parallel default(shared)
            !$omp masked
            !$ nteam_run = omp_get_num_threads()
            !$omp end masked
            !$omp end parallel
            call net%predict(xb, pn)
            !$ call omp_set_num_threads(1)
            call assert_int(NTEAM, nteam_run, 'the identity check''s team of three formed')
            call assert_true(all(transfer(p1, 0_int64, size(p1)) == transfer(pn, 0_int64, size(pn))), &
                &'predict is the same bit for bit at one thread and in a team of three')
        else
            write(*,'(A)') 'test_predict: the team of three skipped under the thread limit'
        endif
        !$ call omp_set_num_threads(nthr_prev)
        call net%kill
    end subroutine test_predict

    ! ===== THE READER =====

    !> the reader refuses a file cut short, a wrong version, a wrong magic and data after the last tensor
    subroutine test_reader_refusals()
        character(len=*), parameter :: TMP = 'tmp_cavg_cnn_tester_weights.bin'
        type(cavg_cnn) :: net
        character(len=:), allocatable :: msg
        integer(int32), allocatable :: bytes(:)
        integer(int64) :: word
        integer :: funit, ios, nbytes, nwords
        logical :: ok
        write(*,'(A)') 'test_reader_refusals'
        inquire(file=FIXTURE_DIR//'weights.raw', size=nbytes)
        nwords = nbytes / 4
        allocate(bytes(nwords))
        open(newunit=funit, file=FIXTURE_DIR//'weights.raw', access='stream', form='unformatted', status='old', iostat=ios)
        if( ios /= 0 ) THROW_HARD('missing fixture testdata/cavg_cnn/weights.raw')
        read(funit) bytes
        close(funit)
        call net%read_weights(string(FIXTURE_DIR//'weights.raw'), ok, msg)
        call assert_true(ok, 'the fixture''s weights file is read')
        call write_words(bytes(1:nwords/2))
        call net%read_weights(string(TMP), ok, msg)
        call assert_false(ok, 'a weights file cut in half is refused')
        call write_words([bytes(1:2), 2, 0, bytes(5:nwords)])
        call net%read_weights(string(TMP), ok, msg)
        call assert_false(ok, 'a weights file of another version is refused')
        word = 12345_int64
        call write_words([transfer(word, 0_int32, 2), bytes(3:nwords)])
        call net%read_weights(string(TMP), ok, msg)
        call assert_false(ok, 'a file without the magic is refused')
        call write_words([bytes, 0])
        call net%read_weights(string(TMP), ok, msg)
        call assert_false(ok, 'a weights file with data after the last tensor is refused')
        call del_file(TMP)
        call net%kill

      contains

        subroutine write_words( words )
            integer(int32), intent(in) :: words(:)
            integer :: u
            open(newunit=u, file=TMP, access='stream', form='unformatted', status='replace')
            write(u) words
            close(u)
        end subroutine write_words

    end subroutine test_reader_refusals

end submodule simple_cavg_cnn_tester
