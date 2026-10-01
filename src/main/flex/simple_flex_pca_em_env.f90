!@descr: flex_pca EM: environment overrides, memory/dimension budgets and run-stage subsampling
submodule (simple_flex_pca_em) simple_flex_pca_em_env
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
    module subroutine cov_stage_subsample( build, pinds, nptcls, nparts, maxtot, env_max, &
        &label, spinds, nsel )
        type(builder),        intent(inout) :: build
        integer,              intent(in)    :: pinds(:), nptcls, nparts, maxtot
        character(len=*),     intent(in)    :: env_max, label
        integer, allocatable, intent(out)   :: spinds(:)
        integer,              intent(out)   :: nsel
        integer :: nmax_tot, nmax_part, ihalf, i, nkept, n_half, ntgt
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
            do i = 1, nptcls
                if( build%spproj_field%get_eo(pinds(i)) /= ihalf ) cycle
                ! real(dp) rather than integer products: nkept*ntgt overflows int32 at these sizes
                if( int(real(nkept+1,dp)*real(ntgt,dp)/real(n_half,dp)) > &
                   &int(real(nkept,  dp)*real(ntgt,dp)/real(n_half,dp)) )then
                    nsel = nsel + 1
                    spinds(nsel) = pinds(i)
                endif
                nkept = nkept + 1
            end do
        end do
        if( nsel < 2 ) THROW_HARD('stage subsample left too few particles; raise '//trim(env_max))
        call hpsort(spinds(:nsel))   ! restore project order so batched image reads stay sequential
        if( nsel < nptcls )then
            write(logfhandle,'(A,A,A,I0,A,I0,A)') '>>> FLEX_PCA ',trim(label),' subsample: using ', &
                &nsel,' of ',nptcls,' particles'
            call flush(logfhandle)
        endif
    end subroutine cov_stage_subsample


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



    !>  Bytes of the packed [d(d+1)/2]^2 array model that cov_dim_budget sizes d against.
    pure real(dp) module function cov_accum_bytes( d ) result( nbytes )
        integer, intent(in) :: d
        real(dp) :: n
        n = real(d,dp)*real(d+1,dp)/2.d0   ! Mspk(npk,npk), npk = d(d+1)/2
        nbytes = 8.d0*n*n
    end function cov_accum_bytes

    !>  Largest d with cov_accum_bytes(d) <= COV_ATHR_BUDGET.
    pure integer module function cov_dim_budget() result( d )
        ! d(d+1)/2 = sqrt(BUDGET/8)  =>  d = (-1 + sqrt(1 + 8*sqrt(BUDGET/8)))/2
        d = max(1, int((-1.d0 + sqrt(1.d0 + 8.d0*sqrt(COV_ATHR_BUDGET/8.d0)))/2.d0))
    end function cov_dim_budget

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
