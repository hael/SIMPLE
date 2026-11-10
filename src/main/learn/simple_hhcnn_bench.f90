!@descr: the HHCNN kernel benchmark: per band (ring count) and width on synthetic tensors, the single-thread forward pass split into the blocks and the readout, the backward pass, a training step at one thread and at a team, the team's prediction throughput, and the thread identity of a training step
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
    !! gives the same weights at one thread and at the team
    module subroutine run_hhcnn_bench( pftsz, nks, widths, nclasses, nthr_team, nrep, funit, ok_identity )
        integer, intent(in)  :: pftsz, nks(:), widths(:), nclasses, nthr_team, nrep, funit
        logical, intent(out) :: ok_identity
        type(hhcnn) :: net, net1
        real(sp), allocatable :: x(:,:,:,:), g(:,:), side(:), t(:,:), xr(:,:,:,:)
        real(sp) :: loss, fm
        real(dp) :: rt_fwd, rt_rd, rt_bwd, rt_step1, rt_stepn, rt_predn, us_fwd, us_rd, us_bwd
        integer(timer_int_kind) :: t0
        integer  :: ik, iw, nk, width, ns, i, irep, ic, nch
        logical  :: l_id
        character(len=*), parameter :: HEAD = '# nk width weights | fwd_us blocks_us readout_us bwd_us | '//&
            &'step1_ms stepN_ms steps_per_s_N samples_per_s_N | predN_us_per_ptcl | identity'
        ok_identity = .true.
        write(funit,'(A)') HEAD
        write(logfhandle,'(A)') HEAD
        ns = max(BENCH_BATCH, BENCH_NPRED)
        do ik = 1,size(nks)
            nk = nks(ik)
            allocate(x(pftsz, nk, 2, ns), g(nk, ns), side(ns), t(nclasses, ns), xr(pftsz, nk, 2, nclasses))
            call fill_inputs
            do iw = 1,size(widths)
                width = widths(iw)
                call net%new(pftsz, nk, nclasses, width, 20261110_int64)
                call net%prepare_work(max(1, nthr_team))
                call net%init_templates(xr, g(:,1))
                ! the single-thread forward pass, the readout alone (idempotent on the last block's maps) and the backward
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
                call net%ensure_chunks(CNN_CHUNK, nch)   ! one chunk's gradient buffer for the backward timing
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
                call net%kill
                call net1%kill
            enddo
            deallocate(x, g, side, t, xr)
        enddo

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
            call net%backward_sample(net%work(1), t(:,ii), net%gchunk(1), l)
        end subroutine bwd_one

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
