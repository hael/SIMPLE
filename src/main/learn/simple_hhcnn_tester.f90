!@descr: unit tests of simple_hhcnn: the Hermitian contract of every layer, the rotation test of the template readout, the mirror invariance against the calculator's mirrored references, the gradient check, the tiled blocks against the whole im2col matrix, the toy problem, the thread identity of the chunked reduction, the weights' round trip
submodule (simple_hhcnn) simple_hhcnn_tester
use simple_test_utils
use simple_cmdline,     only: cmdline
use simple_parameters,  only: parameters
use simple_image,       only: image
use simple_projector,   only: projector
use simple_oris,        only: oris
use simple_ori,         only: ori
use simple_sym,         only: sym
use simple_polar_feats, only: polar_feats, polar_feat
use simple_polarft_calc, only: polarft_calc
use simple_string_utils, only: int2str
use simple_rnd,         only: ran3, gasdev
use simple_timer,       only: tic, toc, timer_int_kind
!$ use omp_lib,         only: omp_set_num_threads, omp_get_num_threads, omp_get_thread_limit
implicit none

! the tiny network of the contract and gradient checks
integer, parameter :: T_PFTSZ = 12, T_NK = 8, T_NCLS = 5, T_WIDTH = 2

contains

    module subroutine run_all_hhcnn_tests()
        write(*,'(A)') '**** running all hhcnn tests ****'
        call test_hermitian_contract()
        call test_mirror_invariance()
        call test_ctf_normalisation()
        call test_gradient_check()
        call test_tiled_blocks()
        call test_thread_identity()
        call test_weights_roundtrip()
        call test_toy_problem()
    end subroutine run_all_hhcnn_tests

    !> a random packed input and a positive weight profile from the tester's own generator
    subroutine random_input( x, g, side )
        real(sp), intent(out) :: x(T_PFTSZ, T_NK, 2), g(T_NK), side
        integer :: p, k, part
        do part = 1,2
            do k = 1,T_NK
                do p = 1,T_PFTSZ
                    x(p,k,part) = gasdev()
                enddo
            enddo
        enddo
        do k = 1,T_NK
            g(k) = 0.5 + ran3()
        enddo
        side = 0.3
    end subroutine random_input

    !> the templates of a tiny network from random class tensors (the geometric initialisation on random "reprojections")
    subroutine random_templates( net )
        type(hhcnn), intent(inout) :: net
        real(sp) :: xr(T_PFTSZ, T_NK, 2, T_NCLS), g(T_NK), side
        integer  :: ic
        do ic = 1,T_NCLS
            call random_input(xr(:,:,:,ic), g, side)
        enddo
        g = 1.
        call net%prepare_work(1)
        call net%init_templates(xr, g)
    end subroutine random_templates

    !> the half ring rotated by one grid step with the conjugate wrap: entry p from entry p-1, entry 1 from the
    !! conjugate of the last
    subroutine rotate_one( n1, nang, nring, a, arot )
        integer,  intent(in)  :: n1, nang, nring
        real(sp), intent(in)  :: a(n1, nang, nring, 2)
        real(sp), intent(out) :: arot(n1, nang, nring, 2)
        integer :: p
        do p = 2,nang
            arot(:,p,:,:) = a(:,p-1,:,:)
        enddo
        arot(:,1,:,1) =  a(:,nang,:,1)
        arot(:,1,:,2) = -a(:,nang,:,2)
    end subroutine rotate_one

    logical function bits_equal4( a, b )
        real(sp), intent(in) :: a(:,:,:,:), b(:,:,:,:)
        bits_equal4 = all(transfer(a, [0_int32], size(a)) == transfer(b, [0_int32], size(b)))
    end function bits_equal4

    logical function bits_equal2( a, b )
        real(sp), intent(in) :: a(:,:), b(:,:)
        bits_equal2 = all(transfer(a, [0_int32], size(a)) == transfer(b, [0_int32], size(b)))
    end function bits_equal2

    logical function bits_equal1( a, b )
        real(sp), intent(in) :: a(:), b(:)
        bits_equal1 = all(transfer(a, [0_int32], size(a)) == transfer(b, [0_int32], size(b)))
    end function bits_equal1

    !> the reference convolution of a column matrix with plain loops in a fixed order
    subroutine conv_ref( cout, nrow, ncol, w, col, out )
        integer,  intent(in)  :: cout, nrow, ncol
        real(sp), intent(in)  :: w(cout, nrow), col(nrow, ncol)
        real(sp), intent(out) :: out(cout, ncol)
        integer :: c, q, r
        real(sp) :: s
        do q = 1,ncol
            do c = 1,cout
                s = 0.
                do r = 1,nrow
                    s = s + w(c,r) * col(r,q)
                enddo
                out(c,q) = s
            enddo
        enddo
    end subroutine conv_ref

    !> every layer commutes with conjugation (and the readout's maxima with it, the maximising rotation moving by
    !! half a turn); a rotation by one grid step rotates every feature map by one step (with the wrap) bit for bit,
    !! leaves every class's maximised score unchanged and moves its maximising rotation by one step (ruling R3's
    !! rotation test); the half ring with the conjugate wrap equals the explicit full ring bit for bit
    subroutine test_hermitian_contract()
        type(hhcnn) :: net
        real(sp) :: x(T_PFTSZ, T_NK, 2), xc(T_PFTSZ, T_NK, 2), xr(T_PFTSZ, T_NK, 2), g(T_NK), side
        real(sp) :: xf(1, 2*T_PFTSZ, T_NK, 2)
        real(sp), allocatable :: pre0(:,:,:,:,:), gated0(:,:,:,:,:), pooled0(:,:,:,:,:), rot(:,:,:,:)
        real(sp), allocatable :: col_h(:,:), col_f(:,:), out_h(:,:), out_f(:,:), raw0(:,:,:,:)
        real(sp) :: prob0(T_NCLS), smax0(T_NCLS), smax1(T_NCLS)
        integer  :: rbest0(T_NCLS), rbest1(T_NCLS), ic
        integer  :: b, cout, nkin, nkout, p, k, part, qh, qf, nposf
        logical  :: ok_re, ok_im, ok_rot, ok_col, ok_out, ok_shift
        real(sp) :: maxdiff
        write(*,'(A)') 'test_hermitian_contract'
        call net%new(T_PFTSZ, T_NK, T_NCLS, T_WIDTH, 20261106_int64)
        call random_templates(net)
        call random_input(x, g, side)
        ! the reference pass and its layers
        call net%forward_sample(net%work(1), x, g, side)
        prob0 = net%work(1)%prob
        call net%get_scores(1, smax0, rbest0)
        cout = net%conv(CNN_NBLK)%cout
        allocate(pre0(cout, T_PFTSZ, T_NK, 2, CNN_NBLK), gated0(cout, T_PFTSZ, T_NK, 2, CNN_NBLK), &
            &pooled0(cout, T_PFTSZ, T_NK, 2, CNN_NBLK), source=0.)
        do b = 1,CNN_NBLK
            cout  = net%conv(b)%cout
            nkin  = net%conv(b)%nkin
            nkout = net%conv(b)%nkout
            pre0(1:cout, :, 1:nkin, :, b)     = net%work(1)%blk(b)%pre
            gated0(1:cout, :, 1:nkin, :, b)   = net%work(1)%blk(b)%gated
            pooled0(1:cout, :, 1:nkout, :, b) = net%work(1)%blk(b)%pooled
        enddo
        raw0 = net%work(1)%blk(1)%raw
        ! conjugation: the conjugated input gives every layer's output conjugated, and the same probabilities
        xc(:,:,1) =  x(:,:,1)
        xc(:,:,2) = -x(:,:,2)
        call net%forward_sample(net%work(1), xc, g, side)
        ok_re = .true.
        ok_im = .true.
        do b = 1,CNN_NBLK
            cout  = net%conv(b)%cout
            nkin  = net%conv(b)%nkin
            nkout = net%conv(b)%nkout
            ok_re = ok_re .and. bits_equal2(reshape(net%work(1)%blk(b)%pre(:,:,:,1), [cout, T_PFTSZ*nkin]), &
                &reshape(pre0(1:cout, :, 1:nkin, 1, b), [cout, T_PFTSZ*nkin]))
            ok_im = ok_im .and. bits_equal2(reshape(net%work(1)%blk(b)%pre(:,:,:,2), [cout, T_PFTSZ*nkin]), &
                &reshape(-pre0(1:cout, :, 1:nkin, 2, b), [cout, T_PFTSZ*nkin]))
            ok_re = ok_re .and. bits_equal2(reshape(net%work(1)%blk(b)%gated(:,:,:,1), [cout, T_PFTSZ*nkin]), &
                &reshape(gated0(1:cout, :, 1:nkin, 1, b), [cout, T_PFTSZ*nkin]))
            ok_im = ok_im .and. bits_equal2(reshape(net%work(1)%blk(b)%gated(:,:,:,2), [cout, T_PFTSZ*nkin]), &
                &reshape(-gated0(1:cout, :, 1:nkin, 2, b), [cout, T_PFTSZ*nkin]))
            ok_re = ok_re .and. bits_equal2(reshape(net%work(1)%blk(b)%pooled(:,:,:,1), [cout, T_PFTSZ*nkout]), &
                &reshape(pooled0(1:cout, :, 1:nkout, 1, b), [cout, T_PFTSZ*nkout]))
            ok_im = ok_im .and. bits_equal2(reshape(net%work(1)%blk(b)%pooled(:,:,:,2), [cout, T_PFTSZ*nkout]), &
                &reshape(-pooled0(1:cout, :, 1:nkout, 2, b), [cout, T_PFTSZ*nkout]))
        enddo
        call assert_true(ok_re, 'conjugating the input leaves every layer''s real part unchanged, bit for bit')
        call assert_true(ok_im, 'conjugating the input negates every layer''s imaginary part, bit for bit')
        call net%get_scores(1, smax1, rbest1)
        maxdiff = maxval(abs(smax1 - smax0))
        call assert_true(maxdiff <= 1.e-5 * max(1., maxval(abs(smax0))), 'conjugating the input leaves every class''s maximised score unchanged')
        ok_shift = .true.
        do ic = 1,T_NCLS
            ok_shift = ok_shift .and. modulo(rbest1(ic) - rbest0(ic), 2*T_PFTSZ) == T_PFTSZ
        enddo
        call assert_true(ok_shift, 'and moves every maximising rotation by half a turn')
        maxdiff = maxval(abs(net%work(1)%prob - prob0))
        call assert_true(maxdiff <= 1.e-5, 'conjugating the input leaves the probabilities unchanged')
        ! rotation by one grid step: every feature map rotates by one step with the wrap, bit for bit
        call rotate_one(1, T_PFTSZ, T_NK, x, xr)
        call net%forward_sample(net%work(1), xr, g, side)
        ok_rot = .true.
        do b = 1,CNN_NBLK
            cout  = net%conv(b)%cout
            nkin  = net%conv(b)%nkin
            nkout = net%conv(b)%nkout
            allocate(rot(cout, T_PFTSZ, nkin, 2))
            call rotate_one(cout, T_PFTSZ, nkin, pre0(1:cout, :, 1:nkin, :, b), rot)
            ok_rot = ok_rot .and. bits_equal4(net%work(1)%blk(b)%pre, rot)
            call rotate_one(cout, T_PFTSZ, nkin, gated0(1:cout, :, 1:nkin, :, b), rot)
            ok_rot = ok_rot .and. bits_equal4(net%work(1)%blk(b)%gated, rot)
            deallocate(rot)
            allocate(rot(cout, T_PFTSZ, nkout, 2))
            call rotate_one(cout, T_PFTSZ, nkout, pooled0(1:cout, :, 1:nkout, :, b), rot)
            ok_rot = ok_rot .and. bits_equal4(net%work(1)%blk(b)%pooled, rot)
            deallocate(rot)
        enddo
        call assert_true(ok_rot, 'rotating the input by one grid step rotates every feature map by one step with the wrap, bit for bit')
        call net%get_scores(1, smax1, rbest1)
        maxdiff = maxval(abs(smax1 - smax0))
        call assert_true(maxdiff <= 1.e-5 * max(1., maxval(abs(smax0))), 'every class''s maximised score is unchanged by the rotation')
        ok_shift = .true.
        do ic = 1,T_NCLS
            ! a class scored by the reflected map turns the other way
            if( net%work(1)%refl(ic) )then
                ok_shift = ok_shift .and. modulo(rbest1(ic) - rbest0(ic), 2*T_PFTSZ) == 2*T_PFTSZ - 1
            else
                ok_shift = ok_shift .and. modulo(rbest1(ic) - rbest0(ic), 2*T_PFTSZ) == 1
            endif
        enddo
        call assert_true(ok_shift, 'and every maximising rotation moves by one step, backwards where the reflection scored')
        maxdiff = maxval(abs(net%work(1)%prob - prob0))
        call assert_true(maxdiff <= 1.e-5, 'the probabilities are unchanged by the rotation')
        ! the half ring with the conjugate wrap against the explicit full ring, bit for bit: the column matrices
        ! and the reference convolution of the first block
        nposf = 2 * T_PFTSZ * T_NK
        allocate(col_h(9, 2*T_PFTSZ*T_NK), col_f(9, 2*nposf), out_h(T_WIDTH, 2*T_PFTSZ*T_NK), out_f(T_WIDTH, 2*nposf))
        xf(1, 1:T_PFTSZ, :, :)            =  x
        xf(1, T_PFTSZ+1:2*T_PFTSZ, :, 1)  =  x(:,:,1)
        xf(1, T_PFTSZ+1:2*T_PFTSZ, :, 2)  = -x(:,:,2)
        call net%im2col(1, T_PFTSZ, T_NK, x, .true., col_h)
        call net%im2col(1, 2*T_PFTSZ, T_NK, xf, .false., col_f)
        call conv_ref(T_WIDTH, 9, 2*T_PFTSZ*T_NK, net%conv(1)%w, col_h, out_h)
        call conv_ref(T_WIDTH, 9, 2*nposf, net%conv(1)%w, col_f, out_f)
        ok_col = .true.
        ok_out = .true.
        do part = 1,2
            do k = 1,T_NK
                do p = 1,T_PFTSZ
                    qh = p + T_PFTSZ*(k-1) + T_PFTSZ*T_NK*(part-1)
                    qf = p + 2*T_PFTSZ*(k-1) + nposf*(part-1)
                    ok_col = ok_col .and. bits_equal1(col_h(:,qh), col_f(:,qf))
                    ok_out = ok_out .and. bits_equal1(out_h(:,qh), out_f(:,qf))
                enddo
            enddo
        enddo
        call assert_true(ok_col, 'the half ring''s taps with the conjugate wrap are the full ring''s taps, bit for bit')
        call assert_true(ok_out, 'the half ring''s convolution with the conjugate wrap is the full ring''s, bit for bit')
        maxdiff = maxval(abs(reshape(raw0, [T_WIDTH, 2*T_PFTSZ*T_NK]) - out_h))
        call assert_true(maxdiff <= 1.e-5 * max(1., maxval(abs(out_h))), 'the BLAS convolution agrees with the reference loops')
        deallocate(pre0, gated0, pooled0, raw0, col_h, col_f, out_h, out_f)
        call net%kill
    end subroutine test_hermitian_contract

    !> The mirror invariance: a sample mirrored by the calculator's mirror_ref_bank, at a frame of 0, 180 and a general
    !! angle on the rotation grid (the search absorbs the rotation), scores against every class as the sample does; the
    !! kernels are symmetric along the angle from the start and through training
    subroutine test_mirror_invariance()
        integer,  parameter :: NB = 16
        real,     parameter :: ANG(3) = [0., 180., 75.]  ! 75: a general frame on the rotation grid (5 steps of 15)
        type(hhcnn) :: net
        class(parameters),   allocatable, target :: p
        class(cmdline),      allocatable :: cline
        class(polarft_calc), allocatable :: calc
        complex(sp), allocatable :: bank(:,:,:)
        real(sp) :: x(T_PFTSZ, T_NK, 2), xm(T_PFTSZ, T_NK, 2), g(T_NK), side, s0(T_NCLS), s1(T_NCLS), smx, loss, fm
        real(sp), allocatable :: xb(:,:,:,:), gb(:,:), sb(:), tb(:,:)
        integer  :: r0(T_NCLS), r1(T_NCLS), ia, b, ic, istep, nrots, nsteps_psi, nchk
        logical  :: refl0(T_NCLS), ok_rot
        write(*,'(A)') 'test_mirror_invariance'
        call net%new(T_PFTSZ, T_NK, T_NCLS, T_WIDTH, 20261109_int64)
        call random_templates(net)
        call assert_true(kernels_symmetric(net), 'the kernels are symmetric along the angle')
        allocate(p, cline, calc)
        call cline%set('box',     32.)
        call cline%set('smpd',    1.)
        call cline%set('mskdiam', 24.)
        call cline%set('pftsz',   real(T_PFTSZ))
        call cline%set('nptcls',  2.)
        call cline%set('nthr',    1.)
        call cline%set('ctf',     'no')
        call cline%set('objfun',  'euclid')
        call p%new(cline, silent=.true.)
        call calc%new(p, 2, [1,2], [2,T_NK+1])
        allocate(bank(T_PFTSZ, T_NK, 2))
        call random_input(x, g, side)
        g = 1.
        call net%forward_sample(net%work(1), x, g, side)
        call net%get_scores(1, s0, r0)
        refl0 = net%work(1)%refl
        smx   = maxval(abs(s0))
        nrots = 2 * T_PFTSZ
        do ia = 1,size(ANG)
            bank(:,:,1) = cmplx(x(:,:,1), x(:,:,2), kind=sp)
            bank(:,:,2) = (0., 0.)
            call calc%mirror_ref_bank(bank, [2], [ANG(ia)])
            xm(:,:,1) = real(bank(:,:,2))
            xm(:,:,2) = aimag(bank(:,:,2))
            call net%forward_sample(net%work(1), xm, g, side)
            call net%get_scores(1, s1, r1)
            call assert_true(maxval(abs(s1 - s0)) <= 1.e-4 * smx, 'a mirrored sample scores against every class as the sample')
            ! the pose read-out's reflected case (pose_student%pose_of): the mirror's maximising rotation is the
            ! sample's less the frame angle, the reflection flag swapped
            nsteps_psi = nint(ANG(ia) * real(nrots) / 360.)
            ok_rot = .true.
            nchk   = 0
            do ic = 1,T_NCLS
                if( refl0(ic) .eqv. net%work(1)%refl(ic) )then
                    ok_rot = .false.
                else if( refl0(ic) )then
                    ok_rot = ok_rot .and. modulo(r1(ic) - r0(ic) - nsteps_psi, nrots) == 0
                else
                    ok_rot = ok_rot .and. modulo(r0(ic) - r1(ic) - nsteps_psi, nrots) == 0
                endif
                nchk = nchk + 1
            enddo
            call assert_true(ok_rot .and. nchk == T_NCLS, 'the mirror swaps every class''s reflection and turns its rotation by the frame angle')
        enddo
        ! a few training steps keep the kernels symmetric, bit for bit
        allocate(xb(T_PFTSZ, T_NK, 2, NB), gb(T_NK, NB), sb(NB), tb(T_NCLS, NB))
        do b = 1,NB
            call random_input(xb(:,:,:,b), gb(:,b), sb(b))
            do ic = 1,T_NCLS
                tb(ic,b) = ran3()
            enddo
            tb(:,b) = tb(:,b) / sum(tb(:,b))
        enddo
        do istep = 1,5
            call net%train_batch(NB, xb, gb, sb, tb, 1.e-2, 1.e-5, loss, fm)
        enddo
        call assert_true(kernels_symmetric(net), 'training keeps the kernels symmetric, bit for bit')
        call calc%kill
        call cline%kill
        call net%kill
    end subroutine test_mirror_invariance

    logical function kernels_symmetric( net )
        type(hhcnn), intent(in) :: net
        integer :: b, cin, ci, ik, rlo
        kernels_symmetric = .true.
        do b = 1,CNN_NBLK
            cin = net%conv(b)%cin
            do ik = -1,1
                do ci = 1,cin
                    rlo = ci + 3*cin*(ik+1)
                    if( any(net%conv(b)%w(:,rlo) /= net%conv(b)%w(:,rlo+2*cin)) ) kernels_symmetric = .false.
                enddo
            enddo
        enddo
    end function kernels_symmetric

    !> the loss of one sample at the current weights
    subroutine sample_loss_now( net, x, g, side, t, loss )
        type(hhcnn), intent(inout) :: net
        real(sp),    intent(in)    :: x(T_PFTSZ, T_NK, 2), g(T_NK), side, t(T_NCLS)
        real(dp),    intent(out)   :: loss
        call net%forward_sample(net%work(1), x, g, side)
        call net%sample_loss(T_NCLS, net%work(1)%prob, t, loss)
    end subroutine sample_loss_now

    !> One weight's central difference against its analytic gradient, in single precision. The step is 1e-2: the
    !! rounding of a loss of order one divided by 2h is then about 1e-5, well below the absolute floor, and the
    !! truncation error of order h^2 as well; a kink of the gate or the rectifier inside the step is caught by the
    !! half step, which must agree instead
    subroutine check_entry( net, x, g, side, t, w, grad, nbad, nchk, label )
        type(hhcnn),      intent(inout) :: net
        real(sp),         intent(in)    :: x(T_PFTSZ, T_NK, 2), g(T_NK), side, t(T_NCLS)
        real(sp),         intent(inout) :: w
        real(sp),         intent(in)    :: grad
        integer,          intent(inout) :: nbad, nchk
        character(len=*), intent(in)    :: label
        real(sp), parameter :: TOL_REL = 1.e-2, TOL_ABS = 2.e-4, STEP = 1.e-2
        real(dp) :: num
        real(sp) :: step_here
        logical  :: ok
        integer  :: ihalf
        ! the step halved up to three times: a harmonic magnitude or modulus near zero within the step has a
        ! curvature the central difference at STEP cannot follow, below STEP/8 the single-precision loss is noise
        step_here = STEP
        do ihalf = 0,3
            num = central_difference(net, x, g, side, t, w, step_here)
            ok  = abs(num - real(grad,dp)) <= TOL_REL * max(abs(num), abs(real(grad,dp))) + TOL_ABS
            if( ok ) exit
            step_here = 0.5 * step_here
        enddo
        nchk = nchk + 1
        if( .not. ok )then
            nbad = nbad + 1
            write(*,'(A,A,2(1X,ES12.4))') '  gradient mismatch ', label, num, grad
        endif
    end subroutine check_entry

    real(dp) function central_difference( net, x, g, side, t, w, h ) result( num )
        type(hhcnn), intent(inout) :: net
        real(sp),    intent(in)    :: x(T_PFTSZ, T_NK, 2), g(T_NK), side, t(T_NCLS), h
        real(sp),    intent(inout) :: w
        real(dp) :: lp, lm
        real(sp) :: w0
        w0 = w
        w  = w0 + h
        call sample_loss_now(net, x, g, side, t, lp)
        w = w0 - h
        call sample_loss_now(net, x, g, side, t, lm)
        w  = w0
        num = (lp - lm) / (2.d0 * real(h,dp))
    end function central_difference

    !> The first block's CTF normalisation: the normalised response to a filter's own CTF-modulated template (the
    !! observation CTF u for the unit-norm filter u, so the packed input CTF^2 u) is one at every defocus where the
    !! taps' rings carry signal, and the raw response is the packed Euclidean product with the filter's footprint
    !! (the objective's cross term, by the polar features' contract)
    subroutine test_ctf_normalisation()
        integer,  parameter :: NDEFOCUS = 4, P0 = T_PFTSZ / 2
        real(sp), parameter :: DEFOCUS(NDEFOCUS) = [0.6, 1.0, 1.5, 2.2]
        type(hhcnn) :: net
        real(sp) :: x(T_PFTSZ, T_NK, 2), g(T_NK), ctf(T_NK), side, wnorm, resp, gbar, raw_re, raw_im
        real(sp) :: xr(T_PFTSZ, T_NK, 2), gr(T_NK)
        integer  :: idef, c, k0, ia, ik, ntested, nbad_one, nbad_raw
        logical  :: l_signal
        write(*,'(A)') 'test_ctf_normalisation'
        call net%new(T_PFTSZ, T_NK, T_NCLS, T_WIDTH, 20261112_int64)
        call random_templates(net)
        side     = 0.
        ntested  = 0
        nbad_one = 0
        do idef = 1,NDEFOCUS
            ! a CTF over the rings with zeros inside the band, its ring-averaged square the normaliser's input
            do k0 = 1,T_NK
                ctf(k0) = sin(0.11 * DEFOCUS(idef) * real(k0 + 1)**2)
            enddo
            g    = ctf**2
            gbar = sum(g) / real(T_NK)
            do c = 1,net%conv(1)%cout
                wnorm = sqrt(sum(net%conv(1)%w(c,:)**2))
                do k0 = 2,T_NK-1
                    ! the template of filter c at (P0, k0): CTF^2 times the unit-norm kernel, real channel
                    x = 0.
                    l_signal = .true.
                    do ik = -1,1
                        l_signal = l_signal .and. g(k0+ik) >= 0.2 * gbar
                        do ia = -1,1
                            x(P0+ia, k0+ik, 1) = g(k0+ik) * net%conv(1)%w(c, 1 + (ia+1) + 3*(ik+1)) / wnorm
                        enddo
                    enddo
                    if( .not. l_signal ) cycle
                    call net%forward_sample(net%work(1), x, g, side)
                    resp    = net%work(1)%blk(1)%pre(c, P0, k0, 1)
                    ntested = ntested + 1
                    if( abs(resp - 1.) > 2.e-2 ) nbad_one = nbad_one + 1
                enddo
            enddo
        enddo
        call assert_true(ntested >= NDEFOCUS * net%conv(1)%cout, 'templates tested at every defocus')
        call assert_int(0, nbad_one, 'the normalised response to the filter''s own CTF-modulated template is one at every defocus')
        ! the raw response to a random input is the packed product with the footprint, both channels
        call random_input(xr, gr, side)
        call net%forward_sample(net%work(1), xr, gr, side)
        nbad_raw = 0
        do c = 1,net%conv(1)%cout
            do k0 = 2,T_NK-1
                raw_re = 0.
                raw_im = 0.
                do ik = -1,1
                    do ia = -1,1
                        raw_re = raw_re + net%conv(1)%w(c, 1 + (ia+1) + 3*(ik+1)) * xr(P0+ia, k0+ik, 1)
                        raw_im = raw_im + net%conv(1)%w(c, 1 + (ia+1) + 3*(ik+1)) * xr(P0+ia, k0+ik, 2)
                    enddo
                enddo
                if( abs(net%work(1)%blk(1)%raw(c, P0, k0, 1) - raw_re) > 1.e-5 * max(1., abs(raw_re)) ) nbad_raw = nbad_raw + 1
                if( abs(net%work(1)%blk(1)%raw(c, P0, k0, 2) - raw_im) > 1.e-5 * max(1., abs(raw_im)) ) nbad_raw = nbad_raw + 1
            enddo
        enddo
        call assert_int(0, nbad_raw, 'the raw response is the packed Euclidean product with the filter''s footprint')
        call net%kill
    end subroutine test_ctf_normalisation

    !> the gradient check of every layer and of the readout (the pointwise map, the templates, the scale) in single
    !! precision
    subroutine test_gradient_check()
        type(hhcnn) :: net
        type(cnn_grad) :: gr
        real(sp) :: x(T_PFTSZ, T_NK, 2), g(T_NK), side, t(T_NCLS)
        real(dp) :: loss
        integer  :: b, i, nbad, nchk, c, r, n, kb, m, ic, part, it, nch
        character(len=40) :: label
        write(*,'(A)') 'test_gradient_check'
        call net%new(T_PFTSZ, T_NK, T_NCLS, T_WIDTH, 20261107_int64)
        call random_templates(net)
        call random_input(x, g, side)
        do i = 1,T_NCLS
            t(i) = 0.1 + ran3()
        enddo
        t = t / sum(t)
        ! a positive gate bias so that its branch is exercised without a kink inside the step
        do b = 1,CNN_NBLK
            net%conv(b)%b = 0.05
        enddo
        call net%alloc_grad(gr)
        call net%zero_grad(gr)
        call net%ensure_chunks(1, nch)
        call net%forward_sample(net%work(1), x, g, side)
        call net%backward_sample(net%work(1), t, gr, loss, 1)
        call net%template_grads(1)
        nbad = 0
        nchk = 0
        do b = 1,CNN_NBLK
            do r = 1,size(net%conv(b)%w,2)
                do c = 1,net%conv(b)%cout
                    if( b > 1 .and. mod(c + r, 5) /= 0 ) cycle
                    write(label,'(A,I0,A,I0,A,I0)') 'conv', b, ' w(', c, ',', r
                    call check_entry(net, x, g, side, t, net%conv(b)%w(c,r), gr%blk(b)%dw(c,r), nbad, nchk, trim(label))
                enddo
            enddo
            do c = 1,net%conv(b)%cout
                write(label,'(A,I0,A,I0)') 'conv', b, ' b(', c
                call check_entry(net, x, g, side, t, net%conv(b)%b(c), gr%blk(b)%db(c), nbad, nchk, trim(label))
            enddo
        enddo
        do c = 1,net%rd%cout
            do n = 1,CNN_NREAD
                write(label,'(A,I0,A,I0)') 'readout wread(', n, ',', c
                call check_entry(net, x, g, side, t, net%rd%wread(n,c), gr%dwread(n,c), nbad, nchk, trim(label))
            enddo
        enddo
        do ic = 1,T_NCLS
            do m = 0,CNN_MREAD
                do kb = 1,net%rd%kread
                    do n = 1,CNN_NREAD
                        do part = 1,2
                            if( m == 0 .and. part == 2 ) cycle
                            if( mod(n + kb + m + ic + part, 7) /= 0 ) cycle
                            it = net%tix(n, kb, m, ic, part)
                            write(label,'(A,I0,A,I0,A,I0,A,I0,A,I0)') 'template(', n, ',', kb, ',', m, ',', ic, ',', part
                            call check_entry(net, x, g, side, t, net%rd%tmpl(it), net%gtmpl(it), nbad, nchk, trim(label))
                        enddo
                    enddo
                enddo
            enddo
        enddo
        call check_entry(net, x, g, side, t, net%rd%scale, gr%dscale, nbad, nchk, 'readout scale')
        write(*,'(A,I0,A,I0)') '  single-precision gradient entries checked: ', nchk, ', mismatches: ', nbad
        call assert_int(0, nbad, 'every layer''s and the readout''s analytic gradient matches the central difference at 1e-2 in single precision')
        call net%grad_kill(gr)
        call net%kill
    end subroutine test_gradient_check

    !> blocks 2 and 3 run in tiles of CNN_TILE columns (several, one across the real/imaginary boundary, the last a
    !! part tile); the responses, the weights' gradient and the input's equal the whole im2col matrix's products
    subroutine test_tiled_blocks()
        ! block 2: 384 columns, three tiles with the boundary at 192 inside the second; block 3: 192, a tile and a half
        integer,  parameter :: PFTSZ = 48
        ! single-precision sums of at most 384 terms in another order, relative to the largest entry
        real(sp), parameter :: TOL = 1.e-5
        type(hhcnn)    :: net
        type(cnn_grad) :: gr
        real(sp) :: x(PFTSZ, T_NK, 2), xr(PFTSZ, T_NK, 2, T_NCLS), g(T_NK), side, t(T_NCLS)
        real(sp), allocatable :: col(:,:), out(:,:), dpre(:,:), dcol(:,:), dpooled(:,:,:,:), dw(:,:)
        real(dp) :: loss
        integer  :: b, cin, cout, nkin, npos, nch, i, ic, p, k, part
        write(*,'(A)') 'test_tiled_blocks'
        call assert_true(2 * PFTSZ * (T_NK / 2) > 2 * CNN_TILE, 'the second block runs in three tiles or more')
        do ic = 0,T_NCLS
            do part = 1,2
                do k = 1,T_NK
                    do p = 1,PFTSZ
                        if( ic == 0 )then
                            x(p,k,part) = gasdev()
                        else
                            xr(p,k,part,ic) = gasdev()
                        endif
                    enddo
                enddo
            enddo
        enddo
        do i = 1,T_NK
            g(i) = 0.5 + ran3()
        enddo
        side = 0.3
        do i = 1,T_NCLS
            t(i) = 0.1 + ran3()
        enddo
        t = t / sum(t)
        call net%new(PFTSZ, T_NK, T_NCLS, T_WIDTH, 20261110_int64)
        call net%prepare_work(1)
        call net%init_templates(xr, [(1., i = 1,T_NK)])
        call net%alloc_grad(gr)
        call net%zero_grad(gr)
        call net%ensure_chunks(1, nch)
        call net%forward_sample(net%work(1), x, g, side)
        call net%backward_sample(net%work(1), t, gr, loss, 1)
        do b = 2,CNN_NBLK
            cin  = net%conv(b)%cin
            cout = net%conv(b)%cout
            nkin = net%conv(b)%nkin
            npos = net%conv(b)%npos
            allocate(col(9*cin, 2*npos), out(cout, 2*npos), dcol(9*cin, 2*npos), dpooled(cin, PFTSZ, nkin, 2))
            call im2col(cin, PFTSZ, nkin, net%work(1)%blk(b-1)%pooled, .true., col)
            call conv_ref(cout, 9*cin, 2*npos, net%conv(b)%w, col, out)
            call assert_true(maxval(abs(out - reshape(net%work(1)%blk(b)%pre, [cout, 2*npos]))) <= TOL * maxval(abs(out)), &
                &'the tiled responses of block '//int2str(b)//' are the whole matrix products')
            dpre = reshape(net%work(1)%blk(b)%dpre, [cout, 2*npos])
            dw   = matmul(dpre, transpose(col))
            call assert_true(maxval(abs(dw - gr%blk(b)%dw)) <= TOL * maxval(abs(dw)), &
                &'the tiled weight gradient of block '//int2str(b)//' is the whole matrix product')
            dcol = matmul(transpose(net%conv(b)%w), dpre)
            call col2im(cin, PFTSZ, nkin, dcol, .true., dpooled)
            call assert_true(maxval(abs(dpooled - net%work(1)%blk(b-1)%dpooled)) <= TOL * maxval(abs(dpooled)), &
                &'the tiled input gradient of block '//int2str(b)//' is the whole matrix product')
            deallocate(col, out, dpre, dcol, dpooled, dw)
        enddo
        call net%grad_kill(gr)
        call net%kill
    end subroutine test_tiled_blocks

    !> the same training steps at one thread and at an actual team of four (the fast gate's cap) give the same weights,
    !! bit for bit, on the build's arithmetic path; a team of sixteen too where the thread limit allows it (by hand: a
    !! serial run fails the check that the team formed)
    subroutine test_thread_identity()
        integer, parameter :: N = 40, NSTEPS = 3, NTEAMS(2) = [4, 16]
        type(hhcnn) :: net1, netn
        real(sp), allocatable :: x(:,:,:,:), g(:,:), side(:), t(:,:)
        real(sp) :: loss1, lossn, fm1, fmn
        integer  :: i, istep, nthr_prev, nteam_run, it, limit
        write(*,'(A)') 'test_thread_identity'
        allocate(x(T_PFTSZ, T_NK, 2, N), g(T_NK, N), side(N), t(T_NCLS, N))
        do i = 1,N
            call random_input(x(:,:,:,i), g(:,i), side(i))
            t(:,i) = 0.
            t(1 + mod(i, T_NCLS), i) = 1.
        enddo
        nthr_prev = 1
        !$ nthr_prev = omp_get_max_threads()
        limit = NTEAMS(1)
        !$ limit = omp_get_thread_limit()
        call net1%new(T_PFTSZ, T_NK, T_NCLS, T_WIDTH, 20261108_int64)
        call net1%prepare_work(1)
        call net1%init_templates(x(:,:,:,1:T_NCLS), g(:,1))
        !$ call omp_set_num_threads(1)
        do istep = 1,NSTEPS
            call net1%train_batch(N, x, g, side, t, 1.e-3, 1.e-5, loss1, fm1)
        enddo
        do it = 1,size(NTEAMS)
            if( it > 1 .and. limit < NTEAMS(it) )then
                write(*,'(A,I0,A)') 'test_thread_identity: the team of ', NTEAMS(it), ' skipped under the thread limit'
                cycle
            endif
            call netn%new(T_PFTSZ, T_NK, T_NCLS, T_WIDTH, 20261108_int64)
            call netn%prepare_work(1)
            call netn%init_templates(x(:,:,:,1:T_NCLS), g(:,1))
            !$ call omp_set_num_threads(NTEAMS(it))
            do istep = 1,NSTEPS
                call netn%train_batch(N, x, g, side, t, 1.e-3, 1.e-5, lossn, fmn)
            enddo
            nteam_run = NTEAMS(it)
            !$omp parallel default(shared)
            !$omp master
            !$ nteam_run = omp_get_num_threads()
            !$omp end master
            !$omp end parallel
            !$ call omp_set_num_threads(1)
            call assert_int(NTEAMS(it), nteam_run, 'the identity check''s team of '//int2str(NTEAMS(it))//' formed: a serial run proves nothing')
            call assert_true(net1%weights_identical(netn), 'the chunked reduction gives the same weights at one thread and at a team of '//int2str(NTEAMS(it))//', bit for bit')
            call assert_true(transfer(loss1, 0_int32) == transfer(lossn, 0_int32), 'the batch loss is the same at one thread and at a team, bit for bit')
            call netn%kill
        enddo
        !$ call omp_set_num_threads(nthr_prev)
        call net1%kill
        deallocate(x, g, side, t)
    end subroutine test_thread_identity

    !> the weights written and read back are the same bit for bit, and give the same probabilities
    subroutine test_weights_roundtrip()
        type(hhcnn) :: net, back
        type(string)   :: cwd_saved, fixture_root, fname
        real(sp) :: x(T_PFTSZ, T_NK, 2), g(T_NK), side, prob0(T_NCLS), prob1(T_NCLS)
        integer  :: nfail0
        write(*,'(A)') 'test_weights_roundtrip'
        nfail0 = tests_failed
        call enter_fixture('hhcnn_weights', cwd_saved, fixture_root)
        call net%new(T_PFTSZ, T_NK, T_NCLS, T_WIDTH, 20261109_int64)
        call random_templates(net)
        fname = 'weights.bin'
        call net%write_weights(fname)
        call back%read_weights(fname)
        call assert_true(net%weights_identical(back), 'the weights file round trip is exact')
        call assert_true(back%get_seed() == 20261109_int64, 'the seed travels with the weights')
        call random_input(x, g, side)
        call net%prepare_work(1)
        call back%prepare_work(1)
        call net%predict(1, x, g, side, prob0)
        call back%predict(1, x, g, side, prob1)
        call assert_true(bits_equal1(prob0, prob1), 'the read network predicts as the written one, bit for bit')
        call fname%kill
        call leave_fixture(cwd_saved, fixture_root, nfail0)
        call net%kill
        call back%kill
    end subroutine test_weights_roundtrip

    !> The toy problem: projections of a three-blob volume over 20 directions, each image a rotated, shifted
    !! and noisy copy at SNR 1 (add_gauran's rule: noise variance = the projection's); the network must classify
    !! held-out images at 0.95 or better
    subroutine test_toy_problem()
        real :: acc
        write(*,'(A)') 'test_toy_problem'
        call run_toy(8, 20, 1.e-2, 2.0, 1.0, 150, 32, acc)
        call assert_true(acc >= 0.95, 'the toy problem''s held-out accuracy is at least 0.95')
    end subroutine test_toy_problem

    !> one toy run: the network of width channels, nepoch epochs at rate lr in batches of nbatch, ntrain training
    !! images per direction, each shifted uniformly within shift_max pixels and noised at the image SNR snr; acc
    !! the held-out accuracy
    subroutine run_toy( width, nepoch, lr, shift_max, snr, ntrain, nbatch, acc )
        integer, intent(in)  :: width, nepoch, ntrain, nbatch
        real,    intent(in)  :: lr, shift_max, snr
        real,    intent(out) :: acc
        integer, parameter :: BOX = 32, PFTSZ = 16, NDIR = 20, NTEST = 20, KFROMTO(2) = [2, 9]
        integer, parameter :: NK = KFROMTO(2) - KFROMTO(1) + 1, PDIM(3) = [PFTSZ, KFROMTO(1), KFROMTO(2)]
        real,    parameter :: SMPD = 1.0, MSKDIAM = 24.
        type(parameters)   :: p
        type(cmdline)      :: cline
        type(image)        :: vol, img_pad, img, noise
        type(projector)    :: vol_pad
        type(oris)         :: dirs
        type(ori)          :: o
        type(sym)          :: symop
        type(polar_feats)  :: feats
        type(polar_feat)   :: feat
        type(hhcnn)     :: net
        complex(sp), allocatable :: refs(:,:,:), pft(:,:), pft_rot(:,:), pft_sh(:,:)
        real(sp),    allocatable :: rmat(:,:,:), xtr(:,:,:,:), gtr(:,:), str(:), ttr(:,:), xte(:,:,:,:), gte(:,:), ste(:), tte(:,:)
        real(sp),    allocatable :: ctfmat(:,:), wk(:,:), noise_sd(:), xref(:,:,:,:)
        real(sp) :: loss, fm, shift(2), side_mean, side_sd
        real(dp) :: lsum, fsum
        integer  :: i, j, k, idir, isample, ntr, nte, irot, pp, q, nrots, ib, i0, i1, nb
        integer, allocatable :: order(:)
        integer(timer_int_kind) :: t0
        real(timer_int_kind)    :: t_data, t_train
        t0 = tic()
        call cline%set('box',     real(BOX))
        call cline%set('smpd',    SMPD)
        call cline%set('mskdiam', MSKDIAM)
        call cline%set('pftsz',   real(PFTSZ))
        call cline%set('nptcls',  1.)
        call cline%set('nthr',    1.)
        call cline%set('ctf',     'no')
        call cline%set('objfun',  'euclid')
        call p%new(cline)
        call set_fixed_seed(20261110)
        ! the volume: six Gaussian blobs off centre, sharp enough to fill the band
        allocate(rmat(BOX,BOX,BOX))
        do k = 1,BOX
            do j = 1,BOX
                do i = 1,BOX
                    rmat(i,j,k) = exp(-real((i-12)**2 + (j-17)**2 + (k-14)**2) / 8.) + &
                        &0.8 * exp(-real((i-21)**2 + (j-13)**2 + (k-19)**2) / 5.) + &
                        &0.7 * exp(-real((i-16)**2 + (j-22)**2 + (k-10)**2) / 6.) + &
                        &0.6 * exp(-real((i-10)**2 + (j-12)**2 + (k-20)**2) / 5.) + &
                        &0.9 * exp(-real((i-22)**2 + (j-21)**2 + (k-15)**2) / 4.) + &
                        &0.5 * exp(-real((i-15)**2 + (j-10)**2 + (k-13)**2) / 6.)
                enddo
            enddo
        enddo
        call vol%new([BOX,BOX,BOX], SMPD)
        call vol%set_rmat(rmat, .false.)
        call vol%mask3D_soft(p%msk)
        call vol_pad%new([p%box_croppd, p%box_croppd, p%box_croppd], p%smpd_crop, wthreads=.false.)
        call vol%pad_fft(vol_pad)
        call vol_pad%expand_cmat()
        ! the directions and the clean sections, polarized from the real-space projections as particles are
        call dirs%new(NDIR, is_ptcl=.false.)
        call symop%new('c1')
        call symop%build_refspiral(dirs)
        call img_pad%new([p%box_croppd, p%box_croppd, 1], p%smpd_crop, wthreads=.false.)
        call img%new([BOX,BOX,1], SMPD, wthreads=.false.)
        call noise%new([BOX,BOX,1], SMPD, wthreads=.false.)
        call img%memoize4polarize(PDIM)
        nrots = 2 * PFTSZ
        allocate(refs(PFTSZ, KFROMTO(1):KFROMTO(2), NDIR), pft(PFTSZ, KFROMTO(1):KFROMTO(2)), noise_sd(NDIR))
        allocate(pft_rot(PFTSZ, KFROMTO(1):KFROMTO(2)), pft_sh(PFTSZ, KFROMTO(1):KFROMTO(2)))
        do idir = 1,NDIR
            call dirs%get_ori(idir, o)
            call img_pad%zero_and_flag_ft
            call vol_pad%fproject(o, img_pad)
            call img_pad%ifft
            call img_pad%clip(img)
            ! the noise of SNR snr: its variance the clean projection's over snr (add_gauran's rule)
            noise_sd(idir) = sqrt(img%variance() / snr)
            call img%fft
            call img%polarize(PDIM, refs(:,:,idir))
        enddo
        ! the features' geometry: unit shell weights, no CTF
        allocate(wk(KFROMTO(1):KFROMTO(2),2), source=1.)
        allocate(ctfmat(PFTSZ, KFROMTO(1):KFROMTO(2)), source=1.)
        call feats%new(PFTSZ, KFROMTO, BOX, wk)
        ntr = NDIR * ntrain
        nte = NDIR * NTEST
        ! the classes are the spiral's mirror pairs: the first half's direction and its mirror (the network is
        ! invariant to the mirror)
        allocate(xtr(PFTSZ, NK, 2, ntr), gtr(NK, ntr), str(ntr), ttr(NDIR/2, ntr), source=0.)
        allocate(xte(PFTSZ, NK, 2, nte), gte(NK, nte), ste(nte), tte(NDIR/2, nte), source=0.)
        do idir = 1,NDIR
            do isample = 1,ntrain + NTEST
                ! a rotated, shifted and noisy copy of the section
                irot = 1 + int(ran3() * real(nrots))
                if( irot > nrots ) irot = nrots
                do pp = 1,PFTSZ
                    q = modulo(pp - irot, nrots) + 1
                    if( q <= PFTSZ )then
                        pft_rot(pp,:) = refs(q,:,idir)
                    else
                        pft_rot(pp,:) = conjg(refs(q-PFTSZ,:,idir))
                    endif
                enddo
                shift = [shift_max * (2. * ran3() - 1.), shift_max * (2. * ran3() - 1.)]
                call feats%shift_pft(pft_rot, shift, pft_sh)
                call noise%gauran(0., noise_sd(idir))
                call noise%fft
                call noise%polarize(PDIM, pft)
                pft = pft_sh + pft
                call feats%build_from_pft(pft, ctfmat, .true., feat)
                if( isample <= ntrain )then
                    i = (idir-1) * ntrain + isample
                    xtr(:,:,:,i) = feat%x
                    gtr(:,i)     = 1.
                    str(i)       = feat%logpow
                    ttr(min(idir, dirs%get_int(idir,'mirr')),i) = 1.
                else
                    i = (idir-1) * NTEST + isample - ntrain
                    xte(:,:,:,i) = feat%x
                    gte(:,i)     = 1.
                    ste(i)       = feat%logpow
                    tte(min(idir, dirs%get_int(idir,'mirr')),i) = 1.
                endif
            enddo
        enddo
        ! the side scalar standardised by the training set's mean and spread
        side_mean = real(sum(real(str,dp)) / real(ntr,dp))
        side_sd   = real(sqrt(sum((real(str,dp) - real(side_mean,dp))**2) / real(ntr,dp)))
        if( side_sd <= 0. ) side_sd = 1.
        str = (str - side_mean) / side_sd
        ste = (ste - side_mean) / side_sd
        t_data = toc(t0)
        t0     = tic()
        ! the templates from the clean sections (the classes' reprojections), then the training: shuffled
        ! mini-batches, a fixed order per epoch from the tester's generator
        allocate(xref(PFTSZ, NK, 2, NDIR/2), source=0.)
        do idir = 1,NDIR/2
            call feats%build_from_pft(refs(:,:,idir), ctfmat, .true., feat)
            xref(:,:,:,idir) = feat%x
        enddo
        call net%new(PFTSZ, NK, NDIR/2, width, 20261111_int64)
        call net%prepare_work(1)
        call net%init_templates(xref, gtr(:,1))
        allocate(order(ntr))
        nb = (ntr + nbatch - 1) / nbatch
        do i = 1,nepoch
            call shuffle(ntr, order)
            lsum = 0.d0
            fsum = 0.d0
            do ib = 1,nb
                i0 = (ib-1) * nbatch + 1
                i1 = min(ntr, ib * nbatch)
                call net%train_batch(i1 - i0 + 1, xtr(:,:,:,order(i0:i1)), gtr(:,order(i0:i1)), str(order(i0:i1)), &
                    &ttr(:,order(i0:i1)), lr, 1.e-5, loss, fm)
                lsum = lsum + real(loss,dp) * real(i1 - i0 + 1,dp)
                fsum = fsum + real(fm,dp)   * real(i1 - i0 + 1,dp)
            enddo
            if( mod(i, 5) == 0 .or. i == nepoch ) write(*,'(A,I3,A,F7.4,A,F6.3)') '  epoch ', i, ': training loss ', &
                &real(lsum / real(ntr,dp)), ', mode ', real(fsum / real(ntr,dp))
        enddo
        t_train = toc(t0)
        call net%eval_batch(nte, xte, gte, ste, tte, loss, acc)
        write(*,'(A,I0,A,I0,A,ES9.2,A,I0,A,F4.1,A,F4.1,A,I0,A,F6.3,A,F7.4,A,F6.2,A,F6.2,A)') '  toy width ', width, &
            &', epochs ', nepoch, ', lr ', lr, ', batch ', nbatch, ', shifts within ', shift_max, ' px, SNR ', snr, ', ', &
            &ntrain, ' per direction: held-out accuracy ', acc, ', loss ', loss, ' (data ', t_data, ' s, training ', t_train, ' s)'
        call net%kill
        call feats%kill
        call feat%kill
        call img_pad%kill
        call img%kill
        call noise%kill
        call o%kill
        call vol_pad%kill_expanded
        call vol_pad%kill
        call vol%kill
        call dirs%kill
        call symop%kill
        call cline%kill
        deallocate(refs, pft, pft_rot, pft_sh, rmat, xtr, gtr, str, ttr, xte, gte, ste, tte, ctfmat, wk, order, noise_sd, xref)

    contains

        !> a random permutation of 1..n from the tester's generator
        subroutine shuffle( n, perm )
            integer, intent(in)  :: n
            integer, intent(out) :: perm(n)
            integer :: ii, jj, tmp
            do ii = 1,n
                perm(ii) = ii
            enddo
            do ii = n,2,-1
                jj = 1 + int(ran3() * real(ii))
                if( jj > ii ) jj = ii
                tmp      = perm(ii)
                perm(ii) = perm(jj)
                perm(jj) = tmp
            enddo
        end subroutine shuffle

    end subroutine run_toy

end submodule simple_hhcnn_tester
