!@descr: the HHCNN kernel benchmark: per band (ring count) and width on synthetic tensors, the single-thread forward pass split into the blocks and the readout, the backward pass, a training step at one thread and at a team, the team's prediction throughput, the thread identity of a training step, and the arithmetic rates of the kernels with the step's occupancy of the team
submodule (simple_hhcnn) simple_hhcnn_bench
use simple_rnd, only: gasdev, ran3
!$ use omp_lib, only: omp_set_num_threads
implicit none

integer, parameter :: BENCH_BATCH = 256    !< the training step's mini-batch, the production size
integer, parameter :: BENCH_NONE  = 32     !< samples of each single-thread pass measurement
integer, parameter :: BENCH_NPRED = 2048   !< samples of the team's prediction throughput

contains

    !> One table row per ring count and width: pftsz angles, nks rings, nclasses templates; nthr_team the team of
    !! the multi-threaded measurements, nrep the repetitions averaged after one warm-up (the single-thread step is the
    !! identity check's, timed once); rows to funit and the log. ok_identity: a training step from the same weights
    !! gives the same weights at one thread and at the team. A second table gives per row the operation counts from
    !! the shapes (the blocks' products 2 cout 9 cin 2 npos per block and pass, the readout's class products and
    !! rotation scores; the gates, the pooling, im2col and the softmax are not counted), the arithmetic rates of the
    !! kernels timed alone at one thread (the blocks' products, the class products, the rotation scores) and of the
    !! passes and steps, and the training step's chunks against the team
    module subroutine run_hhcnn_bench( pftsz, nks, widths, nclasses, nthr_team, nrep, funit, ok_identity )
        integer, intent(in)  :: pftsz, nks(:), widths(:), nclasses, nthr_team, nrep, funit
        logical, intent(out) :: ok_identity
        type(hhcnn) :: net, net1
        real(sp), allocatable :: x(:,:,:,:), g(:,:), side(:), t(:,:), xr(:,:,:,:)
        character(len=256), allocatable :: rate_rows(:), other_rows(:)
        real(sp) :: loss, fm
        real(dp) :: rt_fwd, rt_rd, rt_bwd, rt_step1, rt_stepn, rt_predn, us_fwd, us_rd, us_bwd
        real(dp) :: us_gemm_f, us_gemm_b, us_spec, us_prod, us_rot
        real(dp) :: f_blk_f, f_blk_b, f_spec, f_prod, f_rot, f_rd_b, f_fwd, f_bwd, f_step
        integer(timer_int_kind) :: t0
        integer  :: ik, iw, nk, width, ns, i, irep, ic, nch, b, nchunk_step, nbusy, nrow
        logical  :: l_id
        character(len=*), parameter :: HEAD = '# nk width weights | fwd_us blocks_us readout_us bwd_us | '//&
            &'step1_ms stepN_ms steps_per_s_N samples_per_s_N | predN_us_per_ptcl | identity'
        character(len=*), parameter :: HEAD_OTHER = '# nk width | forward: im2col gate pool rest | backward: '//&
            &'pool gate col2im rest (rest: the first block''s normalisation, the copies, the readout''s backward)'
        real(dp) :: us_k(6)
        integer  :: ik_k
        character(len=*), parameter :: HEAD_RATES = '# nk width | MFLOP per sample: blocks_fwd blocks_bwd products '//&
            &'rotations; per step: step | GFLOP/s at one thread: blocks_gemm_fwd blocks_gemm_bwd blocks_fwd products '//&
            &'rotations fwd bwd step1; at the team: stepN | step: chunks threads_busy occupancy'
        ok_identity = .true.
        write(funit,'(A)') HEAD
        write(logfhandle,'(A)') HEAD
        ns = max(BENCH_BATCH, BENCH_NPRED)
        allocate(rate_rows(size(nks) * size(widths)), other_rows(size(nks) * size(widths)))
        nrow = 0
        do ik = 1,size(nks)
            nk = nks(ik)
            allocate(x(pftsz, nk, 2, ns), g(nk, ns), side(ns), t(nclasses, ns), xr(pftsz, nk, 2, nclasses))
            call fill_inputs
            do iw = 1,size(widths)
                width = widths(iw)
                call net%new(pftsz, nk, nclasses, width, 20261110_int64)
                call net%prepare_work(max(1, nthr_team))
                call net%init_templates(xr, g(:,1))
                ! the single-thread forward pass, the readout alone (idempotent on the last block's maps) and the backward;
                ! one thread outside a region also keeps an OpenMP BLAS on one core
                !$ call omp_set_num_threads(1)
                call net%forward_sample(net%work(1), x(:,:,:,1), g(:,1), side(1))
                t0 = tic()
                do irep = 1,nrep
                    do i = 1,BENCH_NONE
                        call net%forward_sample(net%work(1), x(:,:,:,i), g(:,i), side(i))
                    enddo
                enddo
                rt_fwd = real(toc(t0),dp)
                t0 = tic()
                do irep = 1,nrep
                    do i = 1,BENCH_NONE
                        call net%readout_forward(net%work(1))
                    enddo
                enddo
                rt_rd = real(toc(t0),dp)
                call net%ensure_chunks(1, nch)   ! one chunk's gradient buffer for the backward timing
                t0 = tic()
                do irep = 1,nrep
                    do i = 1,BENCH_NONE
                        call net%forward_sample(net%work(1), x(:,:,:,i), g(:,i), side(i))
                        call bwd_one(i)
                    enddo
                enddo
                rt_bwd = real(toc(t0),dp) - rt_fwd
                us_fwd = 1.d6 * rt_fwd / real(nrep * BENCH_NONE,dp)
                us_rd  = 1.d6 * rt_rd  / real(nrep * BENCH_NONE,dp)
                us_bwd = 1.d6 * rt_bwd / real(nrep * BENCH_NONE,dp)
                ! the kernels alone on the last sample's buffers: the blocks' products forward and backward, the
                ! readout's spectrum, its class products; the rotation scores (and the maximum, the softmax) are the remainder
                t0 = tic()
                do irep = 1,nrep
                    do i = 1,BENCH_NONE
                        call blocks_gemm_fwd
                    enddo
                enddo
                us_gemm_f = 1.d6 * real(toc(t0),dp) / real(nrep * BENCH_NONE,dp)
                t0 = tic()
                do irep = 1,nrep
                    do i = 1,BENCH_NONE
                        call blocks_gemm_bwd
                    enddo
                enddo
                us_gemm_b = 1.d6 * real(toc(t0),dp) / real(nrep * BENCH_NONE,dp)
                t0 = tic()
                do irep = 1,nrep
                    do i = 1,BENCH_NONE
                        call net%readout_spectrum(net%work(1))
                    enddo
                enddo
                us_spec = 1.d6 * real(toc(t0),dp) / real(nrep * BENCH_NONE,dp)
                t0 = tic()
                do irep = 1,nrep
                    do i = 1,BENCH_NONE
                        call net%readout_products(net%work(1))
                    enddo
                enddo
                us_prod = 1.d6 * real(toc(t0),dp) / real(nrep * BENCH_NONE,dp)
                us_rot  = us_rd - us_spec - us_prod
                ! the blocks' other kernels alone on the last sample's buffers
                us_k = 0.d0
                do ik_k = 1,size(us_k)
                    t0 = tic()
                    do irep = 1,nrep
                        do i = 1,BENCH_NONE
                            call blocks_other(ik_k)
                        enddo
                    enddo
                    us_k(ik_k) = 1.d6 * real(toc(t0),dp) / real(nrep * BENCH_NONE,dp)
                enddo
                !$ call omp_set_num_threads(nthr_team)
                ! a training step at one thread, then at the team, from the same weights: the identity of the chunked
                ! reduction, and the step times
                call net1%new(pftsz, nk, nclasses, width, 20261110_int64)
                call net1%prepare_work(max(1, nthr_team))
                call net1%init_templates(xr, g(:,1))
                !$ call omp_set_num_threads(1)
                t0 = tic()
                call net1%train_batch(BENCH_BATCH, x, g, side, t, 1.e-3, 1.e-5, loss, fm)
                rt_step1 = real(toc(t0),dp)
                !$ call omp_set_num_threads(nthr_team)
                call net%copy_weights_from(net1)
                call net1%new(pftsz, nk, nclasses, width, 20261110_int64)
                call net1%prepare_work(max(1, nthr_team))
                call net1%init_templates(xr, g(:,1))
                call net1%train_batch(BENCH_BATCH, x, g, side, t, 1.e-3, 1.e-5, loss, fm)
                l_id = net1%weights_identical(net)
                if( .not. l_id ) ok_identity = .false.
                t0 = tic()
                do irep = 1,nrep
                    call net1%train_batch(BENCH_BATCH, x, g, side, t, 1.e-3, 1.e-5, loss, fm)
                enddo
                rt_stepn = real(toc(t0),dp) / real(nrep,dp)
                ! the team's prediction throughput
                t0 = tic()
                do irep = 1,nrep
                    !$omp parallel do default(shared) private(i,ic) schedule(static) proc_bind(close)
                    do i = 1,BENCH_NPRED
                        ic = 1
                        !$ ic = omp_get_thread_num() + 1
                        call net1%forward_sample(net1%work(ic), x(:,:,:,i), g(:,i), side(i))
                    enddo
                    !$omp end parallel do
                enddo
                rt_predn = real(toc(t0),dp) / real(nrep * BENCH_NPRED,dp)
                call write_row
                call count_flops
                call rate_row
                call other_row
                call net%kill
                call net1%kill
            enddo
            deallocate(x, g, side, t, xr)
        enddo
        write(funit,'(A)') ''
        write(funit,'(A)') '*** ARITHMETIC RATES (operations counted from the shapes) ***'
        write(funit,'(A)') HEAD_RATES
        write(logfhandle,'(A)') HEAD_RATES
        do i = 1,nrow
            write(funit,'(A)') trim(rate_rows(i))
            write(logfhandle,'(A)') trim(rate_rows(i))
        enddo
        write(funit,'(A)') ''
        write(funit,'(A)') '*** THE BLOCKS'' OTHER KERNELS (us per sample at one thread, all blocks) ***'
        write(funit,'(A)') HEAD_OTHER
        write(logfhandle,'(A)') HEAD_OTHER
        do i = 1,nrow
            write(funit,'(A)') trim(other_rows(i))
            write(logfhandle,'(A)') trim(other_rows(i))
        enddo
        flush(funit)
        deallocate(rate_rows, other_rows)

    contains

        subroutine fill_inputs
            integer :: ii, p, k, part, c
            do ii = 1,ns
                do part = 1,2
                    do k = 1,nk
                        do p = 1,pftsz
                            x(p,k,part,ii) = gasdev()
                        enddo
                    enddo
                enddo
                g(:,ii)  = 0.5 + 0.5 * ran3()
                side(ii) = 0.
                do c = 1,nclasses
                    t(c,ii) = ran3()
                enddo
                t(:,ii) = t(:,ii) / sum(t(:,ii))
            enddo
            do c = 1,nclasses
                xr(:,:,:,c) = x(:,:,:,1 + mod(c-1, ns))
            enddo
        end subroutine fill_inputs

        subroutine bwd_one( ii )
            integer, intent(in) :: ii
            real(dp) :: l
            call net%zero_grad(net%gchunk(1))
            call net%backward_sample(net%work(1), t(:,ii), net%gchunk(1), l, 1)
        end subroutine bwd_one

        !> the blocks' forward products of the sample in net's first workspace, as forward_sample computes them (blocks 2
        !! and 3 tile by tile)
        subroutine blocks_gemm_fwd
            integer :: bb, cin, cout, npos, nq, q0
            do bb = 1,CNN_NBLK
                cin  = net%conv(bb)%cin
                cout = net%conv(bb)%cout
                npos = net%conv(bb)%npos
                if( bb == 1 )then
                    call gemm_nn(cout, 2*npos, 9*cin, net%conv(bb)%w, cout, net%work(1)%blk(bb)%col, 9*cin, &
                        &net%work(1)%blk(bb)%pre, cout, 0.0)
                else
                    do q0 = 1,2*npos,size(net%work(1)%blk(bb)%col,2)
                        nq = min(size(net%work(1)%blk(bb)%col,2), 2*npos - q0 + 1)
                        call gemm_nn(cout, nq, 9*cin, net%conv(bb)%w, cout, net%work(1)%blk(bb)%col, 9*cin, &
                            &net%work(1)%blk(bb)%pre, cout, 0.0)
                    enddo
                endif
            enddo
        end subroutine blocks_gemm_fwd

        !> one of the blocks' other kernels over all blocks, as forward_sample and backward_sample call them: 1 im2col,
        !! 2 the gate's forward, 3 the pooling, 4 the pooling's backward, 5 the gate's backward, 6 col2im
        subroutine blocks_other( which )
            integer, intent(in) :: which
            integer :: bb, cin, cout, nkin, nkout
            do bb = 1,CNN_NBLK
                cin   = net%conv(bb)%cin
                cout  = net%conv(bb)%cout
                nkin  = net%conv(bb)%nkin
                nkout = net%conv(bb)%nkout
                select case(which)
                    case(1)
                        if( bb == 1 )then
                            call im2col(1, pftsz, nkin, x(:,:,:,BENCH_NONE), .true., net%work(1)%blk(1)%col)
                        else
                            call im2col_cols(cin, pftsz, nkin, net%work(1)%blk(bb-1)%pooled, .true., 1, &
                                &size(net%work(1)%blk(bb)%col,2), net%work(1)%blk(bb)%col)
                        endif
                    case(2)
                        call gate_forward(cout, pftsz, nkin, net%conv(bb)%b, net%work(1)%blk(bb)%pre, net%work(1)%blk(bb)%gated)
                    case(3)
                        call pool_forward(cout, pftsz, nkin, nkout, net%work(1)%blk(bb)%gated, net%work(1)%blk(bb)%pooled)
                    case(4)
                        call pool_backward(cout, pftsz, nkin, nkout, net%work(1)%blk(bb)%dpooled, net%work(1)%blk(bb)%dgated)
                    case(5)
                        call gate_backward(cout, pftsz, nkin, net%conv(bb)%b, net%work(1)%blk(bb)%pre, &
                            &net%work(1)%blk(bb)%dgated, net%work(1)%blk(bb)%dpre, net%gchunk(1)%blk(bb)%db)
                    case(6)
                        if( bb > 1 ) call col2im_cols(cin, pftsz, nkin, net%work(1)%blk(bb)%dcol, .true., 1, &
                            &size(net%work(1)%blk(bb)%dcol,2), net%work(1)%blk(bb-1)%dpooled)
                end select
            enddo
        end subroutine blocks_other

        !> the blocks' other kernels' row: each kernel, and the remainders of the blocks' forward and backward
        !! (normalisation, copies, the readout's backward)
        subroutine other_row
            character(len=256) :: line
            write(line,'(I4,1X,I5,A,4(1X,F8.1),A,4(1X,F8.1))') nk, width, ' |', us_k(1), us_k(2), us_k(3), &
                &(us_fwd - us_rd - us_gemm_f) - us_k(1) - us_k(2) - us_k(3), ' |', us_k(4), us_k(5), us_k(6), &
                &(us_bwd - us_gemm_b) - us_k(4) - us_k(5) - us_k(6)
            other_rows(nrow) = line
        end subroutine other_row

        !> the blocks' backward products (the weights' gradient, and the input's except in the first block), as
        !! backward_sample computes them, into the first chunk's buffer
        subroutine blocks_gemm_bwd
            integer :: bb, cin, cout, npos, nq, q0
            do bb = 1,CNN_NBLK
                cin  = net%conv(bb)%cin
                cout = net%conv(bb)%cout
                npos = net%conv(bb)%npos
                if( bb == 1 )then
                    call gemm_nt_acc(cout, 9*cin, 2*npos, net%work(1)%blk(bb)%dpre, cout, net%work(1)%blk(bb)%col, 9*cin, &
                        &net%gchunk(1)%blk(bb)%dw, cout)
                else
                    do q0 = 1,2*npos,size(net%work(1)%blk(bb)%col,2)
                        nq = min(size(net%work(1)%blk(bb)%col,2), 2*npos - q0 + 1)
                        call gemm_nt_acc(cout, 9*cin, nq, net%work(1)%blk(bb)%dpre, cout, net%work(1)%blk(bb)%col, 9*cin, &
                            &net%gchunk(1)%blk(bb)%dw, cout)
                        call gemm_tn(9*cin, nq, cout, net%conv(bb)%w, cout, net%work(1)%blk(bb)%dpre, cout, &
                            &net%work(1)%blk(bb)%dcol, 9*cin)
                    enddo
                endif
            enddo
        end subroutine blocks_gemm_bwd

        !> the operation counts per sample (and per step) from the shapes
        subroutine count_flops
            real(dp) :: f_one, nterm
            integer  :: nkout_last
            f_blk_f = 0.d0
            f_blk_b = 0.d0
            do b = 1,CNN_NBLK
                f_one   = 2.d0 * real(net%conv(b)%cout,dp) * 9.d0 * real(net%conv(b)%cin,dp) * 2.d0 * real(net%conv(b)%npos,dp)
                f_blk_f = f_blk_f + f_one
                f_blk_b = f_blk_b + merge(2.d0, 1.d0, b > 1) * f_one
            enddo
            nkout_last = net%conv(CNN_NBLK)%nkout
            ! the readout map (2 parts, CNN_NREAD channels from cout at every angle and ring) and its spectrum
            f_spec = 2.d0 * 2.d0 * real(CNN_NREAD,dp) * real(net%rd%cout,dp) * real(pftsz,dp) * real(nkout_last,dp) + &
                &4.d0 * real((CNN_MREAD + 1) * net%rd%kread * CNN_NREAD,dp) * real(pftsz,dp)
            ! the readout's matrix products: per order the spectrum's four columns against the templates (forward),
            ! as much again in the backward per sample and once more per sample in the templates' gradient over the
            ! batch; the rotation basis against the map's and the reflection's coefficients
            nterm  = real((CNN_MREAD + 1) * 2 * CNN_NREAD * net%rd%kread,dp) * real(nclasses,dp)
            f_prod = 2.d0 * 4.d0 * nterm
            f_rot  = 2.d0 * real(net%nrots,dp) * real(CNN_NQ,dp) * 2.d0 * real(nclasses,dp)
            f_rd_b = f_prod + f_spec
            f_fwd  = f_blk_f + f_spec + f_prod + f_rot
            f_bwd  = f_blk_b + f_rd_b
            f_step = real(BENCH_BATCH,dp) * (f_fwd + f_bwd + f_prod)   ! the templates' gradient over the batch
            nchunk_step = (BENCH_BATCH + chunk_size(BENCH_BATCH) - 1) / chunk_size(BENCH_BATCH)
            nbusy       = min(nchunk_step, max(1, nthr_team))
        end subroutine count_flops

        subroutine rate_row
            character(len=256) :: line
            real(dp) :: occ
            occ = real(nchunk_step,dp) / real(max(1, nthr_team) * ((nchunk_step + max(1, nthr_team) - 1) / max(1, nthr_team)),dp)
            write(line,'(I4,1X,I5,A,4(1X,F8.2),1X,F9.1,A,8(1X,F7.2),A,1X,F7.2,A,2(1X,I5),1X,F5.2)') nk, width, ' |', &
                &1.d-6 * f_blk_f, 1.d-6 * f_blk_b, 1.d-6 * f_prod, 1.d-6 * f_rot, 1.d-6 * f_step, ' |', &
                &gflops(f_blk_f, us_gemm_f), gflops(f_blk_b, us_gemm_b), gflops(f_blk_f, us_fwd - us_rd), &
                &gflops(f_prod, us_prod), gflops(f_rot, us_rot), gflops(f_fwd, us_fwd), gflops(f_bwd, us_bwd), &
                &gflops(f_step, 1.d6 * rt_step1), ' |', gflops(f_step, 1.d6 * rt_stepn), ' |', nchunk_step, nbusy, occ
            nrow = nrow + 1
            rate_rows(nrow) = line
        end subroutine rate_row

        !> GFLOP/s of f operations in us microseconds
        real(dp) function gflops( f, us )
            real(dp), intent(in) :: f, us
            gflops = 0.d0
            if( us > 0.d0 ) gflops = 1.d-3 * f / us
        end function gflops

        subroutine write_row
            character(len=256) :: line
            write(line,'(I4,1X,I5,1X,I9,A,4(1X,F9.1),A,2(1X,F9.2),2(1X,F10.1),A,1X,F9.1,A,1X,L1)') nk, width, &
                &net1%get_nweights(), ' |', us_fwd, us_fwd - us_rd, us_rd, us_bwd, ' |', 1.d3 * rt_step1, 1.d3 * rt_stepn, &
                &1.d0 / rt_stepn, real(BENCH_BATCH,dp) / rt_stepn, ' |', 1.d6 * rt_predn, ' |', l_id
            write(funit,'(A)') trim(line)
            write(logfhandle,'(A)') trim(line)
            flush(funit)
            flush(logfhandle)
        end subroutine write_row

    end subroutine run_hhcnn_bench

end submodule simple_hhcnn_bench
