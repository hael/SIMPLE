!@descr: flex_pca EM: environment overrides, memory/dimension budgets and run-stage subsampling
submodule (simple_flex_pca_em) simple_flex_pca_em_env
use simple_flex_gpu,        only: flex_gpu_available, flex_gpu_prep_begin_f, flex_gpu_prep_free_f,&
    &flex_gpu_prep_ready
implicit none
#include "simple_local_flags.inc"

! Width of the RIGHT kernel -- the one that reads each image's value at the column frequency.
! Zero uses the shared three-tap KB backprojection stencil for both sides.
real    :: COV_RIGHT_KERNEL_W = 0.0
logical :: cov_rkw_read = .false.

contains


    !> Pearson correlation of two double vectors.
    real(dp) module function corr_dp( a, b, n ) result( r )
        integer,  intent(in) :: n
        real(dp), intent(in) :: a(n), b(n)
        real(dp) :: ma, mb, sa, sb, sab
        integer  :: i
        r  = 0.d0
        if( n < 3 ) return
        ma = sum(a)/real(n,dp); mb = sum(b)/real(n,dp)
        sa = 0.d0; sb = 0.d0; sab = 0.d0
        do i = 1, n
            sa  = sa  + (a(i)-ma)**2
            sb  = sb  + (b(i)-mb)**2
            sab = sab + (a(i)-ma)*(b(i)-mb)
        end do
        if( sa <= DTINY .or. sb <= DTINY ) return
        r = sab / sqrt(sa*sb)
    end function corr_dp

    !>  True only when an environment variable is explicitly set to zero (an opt-OUT switch).
    logical module function cov_env_int_off( name ) result(off)
        character(len=*), intent(in) :: name
        character(len=32) :: envval
        integer :: stat, ln, ival
        off = .false.
        call get_environment_variable(name, envval, ln, stat)
        if( stat /= 0 .or. ln < 1 ) return
        read(envval(:ln), *, iostat=stat) ival
        if( stat == 0 ) off = ival == 0
    end function cov_env_int_off

    ! True only when set AND non-zero. Not the complement of cov_env_int_off: an opt-in reads unset as OFF.
    logical module function cov_env_int_on( name ) result(on)
        character(len=*), intent(in) :: name
        character(len=32) :: envval
        integer :: stat, ln, ival
        on = .false.
        call get_environment_variable(name, envval, ln, stat)
        if( stat /= 0 .or. ln < 1 ) return
        read(envval(:ln), *, iostat=stat) ival
        if( stat == 0 ) on = ival /= 0
    end function cov_env_int_on

    ! Rank at which the Gram spectrum enters its noise bulk. Noise level = median of the lower half,
    ! so the leading signal directions cannot inflate it. Scale-free.
    pure integer module function cov_signal_rank( eval, n ) result( d )
        integer,  intent(in) :: n
        real(dp), intent(in) :: eval(n)          !< DESCENDING eigenvalues
        real(dp) :: noise
        integer  :: lo, m
        d = 1
        if( n < 4 ) return
        lo    = n/2 + 1
        m     = n - lo + 1
        noise = eval(lo + m/2)
        if( noise <= DTINY )then
            d = n
            return
        endif
        d = 0
        do while( d < n )
            if( eval(d+1) <= COV_SIGNAL_FACTOR*noise ) exit
            d = d + 1
        end do
        d = max(1, min(n, d))
    end function cov_signal_rank

    ! Halfset-safe capped subsample, shared by the column-subspace initialiser and the probe EM.
    ! `eo` alternates strictly by particle index, so a plain stride of 2 selects one halfset entirely
    ! and the even/odd FSC that regularises every M-step is then computed against nothing; stride
    ! WITHIN each halfset instead. `maxtot` is a total across processes, so only a WORKER passes
    ! nparts -- the master holds every particle and dividing there inflates the stride by nparts.
    module subroutine cov_stage_subsample( build, pinds, nptcls, nparts, maxtot, env_stride, env_max, &
        &label, spinds, nsel )
        type(builder),        intent(inout) :: build
        integer,              intent(in)    :: pinds(:), nptcls, nparts, maxtot
        character(len=*),     intent(in)    :: env_stride, env_max, label
        integer, allocatable, intent(out)   :: spinds(:)
        integer,              intent(out)   :: nsel
        integer :: stride, nmax_tot, nmax_part, ihalf, i, ii, nkept, n_half, ntgt, vstrat
        logical :: l_stride, l_ostrat
        integer, allocatable :: worder(:)
        real,    allocatable :: okey(:)
        real :: th_os, ph_os
        nmax_tot = maxtot
        call cov_env_int(env_max, nmax_tot)
        ! cap off (the default): hand back every particle, in project order
        if( nmax_tot < 1 )then
            allocate(spinds(nptcls), source=pinds(:nptcls))
            nsel = nptcls
            call hpsort(spinds)
            return
        endif
        nmax_part = max(1, nmax_tot / max(1, nparts))
        ! integer stride quantises the budget badly (100000 against a cap of 80000 wants 1.25, gets 2,
        ! keeps 50000), so default to exact-count Bresenham. env_stride restores it; =1 keeps all.
        stride   = 0
        call cov_env_int(env_stride, stride)
        l_stride = stride > 0
        stride   = max(1, stride)
        ! ---- SIMPLE_COV_OSTRAT=1: walk each halfset in ORIENTATION order instead of project
        ! order before the Bresenham pick. Project order is micrograph order, and orientation
        ! correlates with micrograph (ice, support, charging), so an acquisition-order stride
        ! can over/under-hit view clusters; the sorted walk gives low-discrepancy proportional
        ! coverage of the view sphere -- no colatitude band missed by ordering luck. Key: 48
        ! equal-area colatitude bands (1-cos theta), azimuth spreading within each band.
        ! Measured on the polar branch: +5.4 pts GT basis capture on preferred-orientation
        ! data, nothing on isotropic. Selection happens ONCE per stage; per-iteration
        ! resampling remains measured-negative (this is not an incremental EM).
        vstrat = 0
        call cov_env_int('SIMPLE_COV_OSTRAT', vstrat)
        l_ostrat = vstrat > 0
        allocate(worder(nptcls))
        do i = 1, nptcls
            worder(i) = i
        end do
        if( l_ostrat )then
            allocate(okey(nptcls))
            do i = 1, nptcls
                th_os = build%spproj_field%e2get(pinds(i))
                ph_os = build%spproj_field%e1get(pinds(i))
                okey(i) = real(int((1.0 - cos(th_os*PI/180.0))*24.0))*1000.0 + ph_os + 180.0
            end do
            call hpsort(okey, worder)
            deallocate(okey)
            write(logfhandle,'(A,A,A)') '>>> FLEX_PCA ',trim(label), &
                &' subsample: ORIENTATION-STRATIFIED walk (48 equal-area colatitude bands)'
            call flush(logfhandle)
        endif
        allocate(spinds(nptcls))
        nsel = 0
        do ihalf = 0, 1
            n_half = 0
            do i = 1, nptcls
                if( build%spproj_field%get_eo(pinds(i)) == ihalf ) n_half = n_half + 1
            end do
            if( n_half < 1 ) cycle
            ! split the per-part budget evenly between halfsets, never starving one
            ntgt  = min(n_half, max(1, (nmax_part + 1 - ihalf)/2))
            nkept = 0
            do ii = 1, nptcls
                i = worder(ii)
                if( build%spproj_field%get_eo(pinds(i)) /= ihalf ) cycle
                if( l_stride )then
                    if( mod(nkept, stride) == 0 )then
                        nsel = nsel + 1
                        spinds(nsel) = pinds(i)
                    endif
                else
                    ! real(dp) rather than integer products: nkept*ntgt overflows int32 at these sizes
                    if( int(real(nkept+1,dp)*real(ntgt,dp)/real(n_half,dp)) > &
                       &int(real(nkept,  dp)*real(ntgt,dp)/real(n_half,dp)) )then
                        nsel = nsel + 1
                        spinds(nsel) = pinds(i)
                    endif
                endif
                nkept = nkept + 1
            end do
        end do
        deallocate(worder)
        if( nsel < 2 ) THROW_HARD('stage subsample left too few particles; raise '//trim(env_max))
        if( l_ostrat .and. vstrat >= 2 )then
            ! SIMPLE_COV_OSTRAT=2: KEEP a stratified order instead of restoring project order,
            ! so that under ONLINE windows every contiguous window is view-balanced AND
            ! eo-balanced: deal the selection round-robin over (halfset x colatitude band).
            ! Costs scattered reads; on preferred-orientation data (cFAR ~0.06) the window
            ! composition is what the basis sees each iteration, and acquisition-ordered
            ! windows inherit micrograph-correlated view skew.
            block
                integer, allocatable :: sp2(:), bkey(:), bcnt(:), boff(:), bcur(:), bmem(:)
                integer :: kdeal, nout, ib2, nb2, maxcnt
                real    :: th2
                allocate(sp2(nsel), bkey(nsel))
                do i = 1, nsel
                    th2 = build%spproj_field%e2get(spinds(i))
                    bkey(i) = int((1.0 - cos(th2*PI/180.0))*24.0)               ! 48 equal-area bands
                    bkey(i) = bkey(i)*2 + build%spproj_field%get_eo(spinds(i))  ! x eo
                end do
                nb2 = maxval(bkey) + 1
                allocate(bcnt(0:nb2-1), boff(0:nb2-1), bcur(0:nb2-1), bmem(nsel))
                bcnt = 0
                do i = 1, nsel
                    bcnt(bkey(i)) = bcnt(bkey(i)) + 1
                end do
                boff(0) = 0
                do ib2 = 1, nb2 - 1
                    boff(ib2) = boff(ib2-1) + bcnt(ib2-1)
                end do
                bcur = 0
                do i = 1, nsel
                    ib2 = bkey(i)
                    bcur(ib2) = bcur(ib2) + 1
                    bmem(boff(ib2) + bcur(ib2)) = spinds(i)
                end do
                maxcnt = maxval(bcnt)
                nout = 0
                do kdeal = 1, maxcnt
                    do ib2 = 0, nb2 - 1
                        if( kdeal <= bcnt(ib2) )then
                            nout = nout + 1
                            sp2(nout) = bmem(boff(ib2) + kdeal)
                        endif
                    end do
                end do
                spinds(:nsel) = sp2(:nsel)
                deallocate(sp2, bkey, bcnt, boff, bcur, bmem)
                write(logfhandle,'(A,A,A,I0,A)') '>>> FLEX_PCA ',trim(label), &
                    &' subsample: STRATIFIED ORDER retained (', nb2, ' eo x band buckets, round-robin)'
                call flush(logfhandle)
            end block
        else
            call hpsort(spinds(:nsel))   ! restore project order so batched image reads stay sequential
        endif
        if( nsel < nptcls )then
            write(logfhandle,'(A,A,A,I0,A,I0,A,I0,A)') '>>> FLEX_PCA ',trim(label),' subsample (', &
                &merge(1,0,l_stride),'=stride): using ',nsel,' of ',nptcls,' particles'
            call flush(logfhandle)
        endif
    end subroutine cov_stage_subsample

    !>  Write the two half-data latent solves so the per-particle error model can be calibrated
    !!  against the error actually observed. The halves are disjoint checkerboard subsets of one
    !!  particle's own Fourier samples (see cov_herm_inner's `half` argument), so var(z1 - z2)
    !!  measures that particle's estimation error directly, including the model misspecification the
    !!  analytic posterior covariance cannot express. `prior` is written alongside because the half
    !!  solves are shrunk by it. Gated on SIMPLE_COV_ZHALF.
    module subroutine write_zhalf_replicates( zhalf, prior, nptcls, ncomp )
        real(dp), intent(in) :: zhalf(nptcls,ncomp,2), prior(ncomp)
        integer,  intent(in) :: nptcls, ncomp
        integer :: enable, funit, io_stat
        enable = 0
        call cov_env_int('SIMPLE_COV_ZHALF', enable)
        if( enable <= 0 ) return
        call fopen(funit, file=string('flex_pca_zhalf.bin'), access='stream', action='WRITE', &
            &status='REPLACE', iostat=io_stat)
        if( io_stat /= 0 )then
            THROW_WARN('could not open flex_pca_zhalf.bin; skipping half-solve export')
            return
        endif
        write(funit) nptcls, ncomp
        write(funit) zhalf
        write(funit) prior
        call fclose(funit)
        write(logfhandle,'(A,I0,A,I0,A)') '>>> FLEX_PCA wrote flex_pca_zhalf.bin (nptcls=', &
            &nptcls,' ncomp=',ncomp,')'
        call flush(logfhandle)
    end subroutine write_zhalf_replicates

    !>  Override an integer from the environment, if the variable is set and parses.
    module subroutine cov_env_int_pub( name, val )
        character(len=*), intent(in)    :: name
        integer,          intent(inout) :: val
        call cov_env_int(name, val)
    end subroutine cov_env_int_pub

    !> Is an environment variable present at all? Used where the DEFAULT must be "behave exactly as
    !! before", not "behave as if the variable were at its documented default".
    logical module function cov_env_is_set( name )
        character(len=*), intent(in) :: name
        character(len=32) :: envval
        integer :: stat, ln
        call get_environment_variable(name, envval, ln, stat)
        cov_env_is_set = (stat == 0 .and. ln >= 1)
    end function cov_env_is_set

    module subroutine cov_env_int( name, val )
        character(len=*), intent(in)    :: name
        integer,          intent(inout) :: val
        character(len=32) :: envval
        integer :: stat, ln, ival
        call get_environment_variable(name, envval, ln, stat)
        if( stat /= 0 .or. ln < 1 ) return
        read(envval(:ln), *, iostat=stat) ival
        if( stat == 0 .and. ival > 0 )then
            val = ival
            write(logfhandle,'(A,A,A,I0)') '>>> FLEX_PCA ',trim(name),' override: ',ival
            call flush(logfhandle)
        endif
    end subroutine cov_env_int


    !>  Packed accumulation + matrix-free CG is the DEFAULT reduced solve, so an UNSET
    !!  SIMPLE_COV_CGSOLVE selects packed. SIMPLE_COV_CGSOLVE=0 is the documented escape hatch back to
    !!  the dense d^2 x d^2 accumulator + Cholesky; cov_env_int cannot express that (it ignores every
    !!  value <= 0), which is why this goes through the presence-and-zero test in cov_env_int_off.
    !!  EVERY site whose memory model depends on the choice must call this -- if the dimension budget
    !!  and the solve disagree, d_tilde is sized against an accumulator that is never allocated.
    logical module function cov_packed_cgsolve() result( packed )
        packed = .not. cov_env_int_off('SIMPLE_COV_CGSOLVE')
    end function cov_packed_cgsolve

    !>  Bytes in the reduced solve's ONE shared accumulator at column dimension d.
    pure real(dp) module function cov_accum_bytes( d, packed ) result( nbytes )
        integer, intent(in) :: d
        logical, intent(in) :: packed
        real(dp) :: n
        if( packed )then
            n = real(d,dp)*real(d+1,dp)/2.d0   ! Mspk(npk,npk), npk = d(d+1)/2
        else
            n = real(d,dp)**2                  ! A(d^2,d^2)
        endif
        nbytes = 8.d0*n*n
    end function cov_accum_bytes

    !>  Largest d whose accumulator fits COV_ATHR_BUDGET under the model the solve will ACTUALLY use,
    !!  i.e. cov_accum_bytes(d, packed) <= COV_ATHR_BUDGET.
    pure integer module function cov_dim_budget( packed ) result( d )
        logical, intent(in) :: packed
        if( packed )then
            ! d(d+1)/2 = sqrt(BUDGET/8)  =>  d = (-1 + sqrt(1 + 8*sqrt(BUDGET/8)))/2
            d = max(1, int((-1.d0 + sqrt(1.d0 + 8.d0*sqrt(COV_ATHR_BUDGET/8.d0)))/2.d0))
        else
            d = max(1, int((COV_ATHR_BUDGET/8.d0)**0.25d0))
        endif
    end function cov_dim_budget

    !>  One-shot read of the SIMPLE_COV_RKW override for COV_RIGHT_KERNEL_W.

    !> begin the device prep lifecycle for a stage batch loop when the GPU is present and
    !! enabled; the shared prep funnel (prep_imgs4projected_model) then takes its device
    !! branch for every batch. No-op (l_on=.false.) when an outer lifecycle already owns it.
    module subroutine cov_dev_prep_start( params, build, l_on )
        class(parameters), intent(in)    :: params
        class(builder),    intent(inout) :: build
        logical,           intent(out)   :: l_on
        integer :: vprep
        l_on = .false.
        if( .not. flex_gpu_available() ) return
        if( flex_gpu_prep_ready() )      return   ! outer owner
        ! OPT-IN (SIMPLE_COV_GPU_PREP_STAGES=1): measured on the 4-worker Ribosembly arm, the
        ! fetch+unpack funnel is slower than the 6.4 s threaded CPU prep at these stages (8.4 s
        ! under device contention), and its 1e-6-level numerics nudged the probe convergence
        ! rule from 3 to 5 rounds (+31 s wall). The probe's own resident prep (no fetch) and
        ! the STATEREC resident hand-off are the configurations that pay.
        vprep = 0
        call cov_env_int('SIMPLE_COV_GPU_PREP_STAGES', vprep)
        if( vprep <= 0 ) return
        if( params%l_ml_reg )then
            ! whitening path needs the loaded sigma2 spectra (CPU prep THROWs without them)
            if( .not. allocated(build%esig%sigma2_noise) ) return
        endif
        if( cov_image_mask_radius(params) > 0. ) return   ! mask variant stays on the CPU
        call flex_gpu_prep_begin_f(build%lmsk, params%box, params%boxpd, MAXIMGBATCHSZ, &
            &0.0, .true.)
        l_on = .true.
    end subroutine cov_dev_prep_start

    module subroutine cov_dev_prep_stop( l_on )
        logical, intent(in) :: l_on
        if( l_on ) call flex_gpu_prep_free_f
    end subroutine cov_dev_prep_stop

    module subroutine cov_init_right_kernel_width
        character(len=32) :: envval
        integer :: stat, ln
        real    :: rval
        if( cov_rkw_read ) return
        cov_rkw_read = .true.
        call get_environment_variable('SIMPLE_COV_RKW', envval, ln, stat)
        if( stat == 0 .and. ln > 0 )then
            read(envval(:ln), *, iostat=stat) rval
            if( stat == 0 ) COV_RIGHT_KERNEL_W = rval
        endif
        if( COV_RIGHT_KERNEL_W > 0. )then
            write(logfhandle,'(A,F6.3)') '>>> FLEX_PCA right (column-gather) kernel: triangular, width ', &
                &COV_RIGHT_KERNEL_W
        else
            write(logfhandle,'(A)') '>>> FLEX_PCA right (column-gather) kernel: shared KB backprojection stencil'
        endif
        call flush(logfhandle)
    end subroutine cov_init_right_kernel_width

    !> Sampling precision of the MAP latent estimate, Q = A*Gtil^+*A with A = Gtil + diag(prior). This is
    !! the precision of the ESTIMATOR z_hat, not the posterior precision A, so distances measured with it
    !! reflect how well each component was actually determined for the particle.
    module subroutine map_sampling_precision( Gtil, prior, n, Qout )
        integer,  intent(in)  :: n
        real(dp), intent(in)  :: Gtil(n,n), prior(n)
        real(dp), intent(out) :: Qout(n,n)
        real(dp) :: Amat(n,n), Gpinv(n,n), Vmat(n,n), Awork(n,n), ev(n), thresh
        integer  :: ii, jj, kk, nrot
        Amat = Gtil
        do ii = 1, n
            Amat(ii,ii) = Amat(ii,ii) + prior(ii)
        end do
        Awork = Gtil
        call jacobi(Awork, n, n, ev, Vmat, nrot)   ! symmetric eigendecomposition (LAPACK dsyev)
        thresh = COV_PINV_RCOND * maxval(abs(ev))
        Gpinv  = 0.d0
        do kk = 1, n
            if( abs(ev(kk)) <= thresh ) cycle      ! drop the null space, as pinv does
            do jj = 1, n
                do ii = 1, n
                    Gpinv(ii,jj) = Gpinv(ii,jj) + Vmat(ii,kk)*Vmat(jj,kk)/ev(kk)
                end do
            end do
        end do
        Qout = matmul(Amat, matmul(Gpinv, Amat))
        Qout = 0.5d0*(Qout + transpose(Qout))      ! symmetrise away round-off
    end subroutine map_sampling_precision

end submodule simple_flex_pca_em_env
