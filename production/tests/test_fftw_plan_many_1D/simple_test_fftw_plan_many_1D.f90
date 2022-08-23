program simple_test_fftw_plan_many_1D
    include 'simple_lib.f08'
    use simple_fftw3
    use simple_parameters, only: params_glob
    
    implicit none

    integer, parameter        :: npoints = 1000, narrays = 100, niter = 100
    real, parameter           :: bound = 4
    ! Print all stages of fft to test for correctness
    character(*), parameter   :: fn_in = '/Users/wietfeldthc/Documents/simpleTestImages/testBatchFFT1DAug17/output/in.txt'
    character(*), parameter   :: fn_fft_curr = '/Users/wietfeldthc/Documents/simpleTestImages/testBatchFFT1DAug17/output/fft_curr.txt'
    character(*), parameter   :: fn_out_curr = '/Users/wietfeldthc/Documents/simpleTestImages/testBatchFFT1DAug17/output/out_curr.txt'
    character(*), parameter   :: fn_fft_batch = '/Users/wietfeldthc/Documents/simpleTestImages/testBatchFFT1DAug17/output/fft_batch.txt'
    character(*), parameter   :: fn_out_batch = '/Users/wietfeldthc/Documents/simpleTestImages/testBatchFFT1DAug17/output/out_batch.txt'
    character(*), parameter   :: fn_bench = '/Users/wietfeldthc/Documents/simpleTestImages/testBatchFFT1DAug17/output/BENCH.txt'
    type(c_ptr)               :: plan_fwd, plan_bwd
    complex(sp), allocatable  :: indat(:, :), indat_temp(:, :), cdat(:, :), cdat_temp(:, :)
    real(sp), allocatable     :: rdat(:, :) ! Use different indat and rdat to replicate calc_corrs_over_k
    integer                   :: i, j, k, fnr
    integer(timer_int_kind)   :: t_curr_fft, t_curr_ifft, t_batch_fft, t_batch_ifft
    real(timer_int_kind)      :: rt_curr_fft, rt_curr_ifft, rt_batch_fft, rt_batch_ifft

    ! We could declare these sizes but allocating data replicates how the real data is handled
    allocate(indat(npoints, narrays))
    allocate(indat_temp(npoints, narrays))
    allocate(cdat(npoints, narrays))
    allocate(cdat_temp(npoints, narrays))
    allocate(rdat(npoints, narrays))

    if (L_BENCH_GLOB) then
        rt_curr_fft = 0
        rt_curr_ifft = 0
        rt_batch_fft = 0
        rt_batch_ifft = 0
    end if
    ! Generate arrays
    do j = 1, narrays
        do i = 1, npoints
            indat(i, j) = 0.1 * j * sinc(-1 * bound + i * 2 * bound / npoints)
        end do
    end do
    call fopen(fnr, FILE=trim(fn_in), STATUS='REPLACE', action='WRITE')
    do j = 1, narrays
        do i = 1, npoints
            write (fnr, '(f10.3)') real(indat(i, j))
        end do
    end do
    call fclose(fnr)

    ! Current FFT
    print *, 'Current FFT'
    indat_temp = indat
    call fftwf_plan_with_nthreads(nthr_glob)
    plan_fwd = fftwf_plan_dft_1d(npoints, indat_temp(:,1),  cdat(:,1),  FFTW_FORWARD, FFTW_PATIENT)
    call fftwf_plan_with_nthreads(1)
    do k = 1, niter
        do j = 1, narrays
            if (L_BENCH_GLOB) t_curr_fft = tic()
            call fftwf_execute_dft(plan_fwd,  indat(:,j),  cdat(:,j))
            if (L_BENCH_GLOB) rt_curr_fft = rt_curr_fft + toc(t_curr_fft)
        end do
    end do
    call fopen(fnr, FILE=trim(fn_fft_curr), STATUS='REPLACE', action='WRITE')
    do j = 1, narrays
        do i = 1, npoints
            write(fnr, '(f10.3,SP,f10.3,"i")') cdat(i, j)
        end do
    end do
    call fclose(fnr)

    ! Current iFFT
    print *, 'Current iFFT'
    cdat_temp(:, 1) = cdat(:, 1)
    call fftwf_plan_with_nthreads(nthr_glob)
    plan_bwd = fftwf_plan_dft_c2r_1d(npoints, cdat_temp(:,1),  rdat(:,1), FFTW_PATIENT)
    call fftwf_plan_with_nthreads(1)
    do k = 1, niter
        do j = 1, narrays
            cdat_temp(:, j) = cdat(:, j) ! dft_c2r destroys input data so use temp cdat if niter > 1
            if (L_BENCH_GLOB) t_curr_ifft = tic()
            call fftwf_execute_dft_c2r(plan_bwd,  cdat_temp(:,j),  rdat(:,j))
            if (L_BENCH_GLOB) rt_curr_ifft = rt_curr_ifft + toc(t_curr_ifft)
            rdat(:, j) = 1. / npoints * rdat(:, j) ! So that iFFT(FFT(f(t))) = f(t) for testing correctness
        end do
    end do
    call fopen(fnr, FILE=trim(fn_out_curr), STATUS='REPLACE', action='WRITE')
    do j = 1, narrays
        do i = 1, npoints
            write(fnr, '(f10.3)') rdat(i, j)
        end do
    end do
    call fclose(fnr)

    ! Batch FFT
    print *, 'BATCH FFT'
    indat_temp = indat
    call fftwf_plan_with_nthreads(nthr_glob)
    plan_fwd = fftwf_plan_many_dft(1, [npoints], narrays, indat_temp, [npoints], 1, npoints, cdat, [npoints], 1, npoints, FFTW_FORWARD, FFTW_PATIENT)
    call fftwf_plan_with_nthreads(1)
    do k = 1, niter
        if (L_BENCH_GLOB) t_batch_fft = tic()
        call fftwf_execute_dft(plan_fwd, indat, cdat)
        if (L_BENCH_GLOB) rt_batch_fft = rt_batch_fft + toc(t_batch_fft)
    end do
    call fopen(fnr, FILE=trim(fn_fft_batch), STATUS='REPLACE', action='WRITE')
    do j = 1, narrays
        do i = 1, npoints
            write(fnr, '(f10.3,SP,f10.3,"i")') cdat(i, j)
        end do
    end do
    call fclose(fnr)

    ! Batch iFFT
    print *, 'BATCH IFFT'
    cdat_temp = cdat
    call fftwf_plan_with_nthreads(nthr_glob)
    plan_bwd = fftwf_plan_many_dft_c2r(1, [npoints], narrays, cdat_temp, [npoints], 1, npoints, rdat, [npoints], 1, npoints, FFTW_PATIENT)
    call fftwf_plan_with_nthreads(1)
    do k = 1, niter
        cdat_temp = cdat  ! dft_c2r destroys input data so use temp cdat if niter > 1
        if (L_BENCH_GLOB) t_batch_ifft = tic()
        call fftwf_execute_dft_c2r(plan_bwd, cdat_temp, rdat)
        if (L_BENCH_GLOB) rt_batch_ifft = rt_batch_ifft + toc(t_batch_ifft)
    end do
    rdat = 1. / npoints * rdat ! So that iFFT(FFT(f(t))) = f(t) for testing correctness
    call fopen(fnr, FILE=trim(fn_out_batch), STATUS='REPLACE', action='WRITE')
    do j = 1, narrays
        do i = 1, npoints
            write(fnr, '(f10.3)') rdat(i, j)
        end do
    end do
    call fclose(fnr)

    if( L_BENCH_GLOB )then
        call fopen(fnr, FILE=trim(fn_bench), STATUS='REPLACE', action='WRITE')
        write(fnr,'(a)') '*** TIMINGS (s) ***'
        write(fnr,'(a,1x,f9.4)') 'current_fft    : ', rt_curr_fft
        write(fnr,'(a,1x,f9.4)') 'batch_fft      : ', rt_batch_fft
        write(fnr,'(a,1x,f9.4)') 'current_ifft   : ', rt_curr_ifft
        write(fnr,'(a,1x,f9.4)') 'batch_ifft     : ', rt_batch_ifft
        write(fnr,'(a,1x,f9.4)') 'current_tot    : ', rt_curr_fft + rt_curr_ifft
        write(fnr,'(a,1x,f9.4)') 'batch_tot      : ', rt_batch_fft + rt_batch_ifft
        write(fnr,'(a)') '*** RELATIVE TIMING (%) ***'
        write(fnr,'(a,1x,f9.4)') 'batch_tot / current_tot    : ', (rt_batch_fft + rt_batch_ifft) / (rt_curr_fft + rt_curr_ifft) 
        call fclose(fnr)
    endif

end program simple_test_fftw_plan_many_1D