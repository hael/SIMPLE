!@descr: flex_pca shared helpers: chi-squared median, unimodality test, kernel weights, delivery naming
module simple_flex_pca_util
use simple_core_module_api, only: dp, dtiny, fname2ext, get_fbody, hpsort, logfhandle, mrc_ext, oris, &
    &simple_exception, string, tiny
use simple_defs_flex, only: FLEX_MAX_BW_GROW, FLEX_ACCUM_BYTE_BUDGET
use simple_image,     only: image
use simple_oris,      only: oris
implicit none
private
#include "simple_local_flags.inc"

public :: chi2_median, two_gauss_unimodal
public :: kernel_weights_at_bandwidth, project_onto_target_polyline, dilation_template
public :: flex_pca_write_state
public :: corr_dp, cov_signal_rank, cov_stage_subsample, cov_accum_bytes, cov_dim_budget

real(dp), parameter :: COV_SIGNAL_FACTOR = 4.0d0

contains

    !> Median of chi-squared with k dof, Wilson-Hilferty: k*(1 - 2/(9k))^3. Good to 3 % at k=1, which is
    !! far inside the tolerance of a bandwidth FLOOR and needs no gamma inverse.
    pure real(dp) function chi2_median( k )
        integer, intent(in) :: k
        real(dp) :: kk
        kk = real(max(k,1),dp)
        chi2_median = kk * (1.d0 - 2.d0/(9.d0*kk))**3
    end function chi2_median

    !> Unimodality of the 1-D two-component equal-variance mixture p1*N(0,1) + p2*N(d,1).
    !! With a tied covariance the density along the line between two component means IS this
    !! mixture, so the test is exact, and it is parameter-free: adjacent tiles of one continuum
    !! overlap into a single mode, a discrete island keeps a dip between the modes.
    logical function two_gauss_unimodal( d, p1, p2 )
        real(dp), intent(in) :: d, p1, p2
        real(dp) :: g(0:200), t
        integer  :: j, j1, j2
        if( d <= 2.d0 )then
            ! equal-variance pair is unimodal for ANY mixing weights at separation <= 2 sigma
            two_gauss_unimodal = .true.
            return
        endif
        do j = 0, 200
            t    = d*real(j,dp)/200.d0
            g(j) = p1*exp(-0.5d0*t*t) + p2*exp(-0.5d0*(t - d)*(t - d))
        end do
        j1 = maxloc(g(0:100),   dim=1) - 1
        j2 = maxloc(g(100:200), dim=1) + 99
        if( j2 <= j1 + 1 )then
            two_gauss_unimodal = .true.
            return
        endif
        two_gauss_unimodal = minval(g(j1:j2)) >= (1.d0 - 1.d-9)*min(g(j1), g(j2))
    end function two_gauss_unimodal


    !> Epanechnikov weights at bandwidth h_in on a squared-distance vector; the bandwidth grows
    !! (safety only) until at least min_neff particles have non-zero weight.
    subroutine kernel_weights_at_bandwidth( dist, nptcls, h_in, min_neff, w, h_out, neff_out )
        integer,  intent(in)  :: nptcls, min_neff
        real(dp), intent(in)  :: dist(nptcls), h_in
        real,     intent(out) :: w(nptcls)
        real(dp), intent(out) :: h_out
        real,     intent(out) :: neff_out
        real(dp) :: h, u2, sumw, sumw2
        integer  :: i, grow, nsupp
        h = h_in
        nsupp = 0
        do grow = 0, FLEX_MAX_BW_GROW
            sumw = 0.d0; sumw2 = 0.d0; nsupp = 0
            !$omp parallel do default(shared) private(i,u2) schedule(static) &
            !$omp& reduction(+:sumw,sumw2,nsupp)
            do i = 1, nptcls
                u2 = dist(i) / (h*h)
                w(i) = real(max(0.d0, 1.d0 - u2))
                sumw  = sumw  + real(w(i),dp)
                sumw2 = sumw2 + real(w(i),dp)**2
                if( w(i) > 0. ) nsupp = nsupp + 1
            end do
            !$omp end parallel do
            if( nsupp >= min(min_neff, nptcls) ) exit
            if( grow >= FLEX_MAX_BW_GROW      ) exit
            h = 1.3d0*h
        end do
        if( maxval(w) > TINY ) w = w / maxval(w)
        sumw  = sum(real(w,dp))
        sumw2 = sum(real(w,dp)**2)
        h_out    = h
        neff_out = real(sumw*sumw/max(sumw2,DTINY))
    end subroutine kernel_weights_at_bandwidth

    !> Arc-length coordinate of every particle along the polyline through the ordered targets,
    !! and the targets' own coordinates on it.
    subroutine project_onto_target_polyline( z, nptcls, nk, tcen, nstates, ppath, tpath )
        integer,  intent(in)  :: nptcls, nk, nstates
        real(dp), intent(in)  :: z(nptcls,nk), tcen(nk,nstates)
        real(dp), intent(out) :: ppath(nptcls), tpath(nstates)
        real(dp) :: seg(nk,max(1,nstates-1)), sl2(max(1,nstates-1)), seglen(max(1,nstates-1))
        real(dp) :: dz(nk), t, d2, best_d2, best_c
        integer  :: i, s, q
        tpath(1) = 0.d0
        do s = 1, nstates-1
            seg(:,s)   = tcen(:,s+1) - tcen(:,s)
            sl2(s)     = sum(seg(:,s)**2)
            seglen(s)  = sqrt(sl2(s))
            tpath(s+1) = tpath(s) + seglen(s)
        end do
        !$omp parallel do default(shared) private(i,s,q,dz,t,d2,best_d2,best_c) &
        !$omp& schedule(static) proc_bind(close)
        do i = 1, nptcls
            best_d2 = huge(0.d0)
            best_c  = 0.d0
            do s = 1, nstates-1
                if( sl2(s) <= DTINY ) cycle
                do q = 1, nk
                    dz(q) = z(i,q) - tcen(q,s)
                end do
                t  = max(0.d0, min(1.d0, sum(dz*seg(:,s))/sl2(s)))
                d2 = 0.d0
                do q = 1, nk
                    d2 = d2 + (dz(q) - t*seg(q,s))**2
                end do
                if( d2 < best_d2 )then
                    best_d2 = d2
                    best_c  = tpath(s) + t*seglen(s)
                endif
            end do
            ppath(i) = best_c
        end do
        !$omp end parallel do
    end subroutine project_onto_target_polyline


    !> Replace a real-space volume by its dilation mode (x - c) . grad rho on the box: the direction
    !! a uniform magnification/defocus scatter ('breathing') moves the density along. Central
    !! differences; the outermost layer is zeroed.
    subroutine dilation_template( img, box )
        class(image), intent(inout) :: img
        integer,      intent(in)    :: box
        real, pointer :: r(:,:,:)
        real, allocatable :: d(:,:,:)
        real    :: c, gx, gy, gz
        integer :: i, j, k
        call img%get_rmat_ptr(r)
        allocate(d(box,box,box), source=0.)
        c = real(box/2 + 1)
        do k = 2, box - 1
            do j = 2, box - 1
                do i = 2, box - 1
                    gx = 0.5*(r(i+1,j,k) - r(i-1,j,k))
                    gy = 0.5*(r(i,j+1,k) - r(i,j-1,k))
                    gz = 0.5*(r(i,j,k+1) - r(i,j,k-1))
                    d(i,j,k) = (real(i)-c)*gx + (real(j)-c)*gy + (real(k)-c)*gz
                end do
            end do
        end do
        r = 0.
        r(1:box,1:box,1:box) = d
        deallocate(d)
    end subroutine dilation_template

    !> delivered state map name from outvol: state 1 keeps the name, others get _NNN
    subroutine flex_pca_write_state( outvol, img, state, vol_fname )
        type(string),  intent(in)    :: outvol
        class(image),  intent(inout) :: img
        integer,       intent(in)    :: state
        class(string), intent(inout) :: vol_fname
        type(string) :: prefix, ext
        character(len=:), allocatable :: stem
        character(len=3) :: tag
        if( state==1 )then
            vol_fname = outvol
        else
            ext=fname2ext(outvol)
            prefix=get_fbody(outvol,ext)
            stem=prefix%to_char()
            if( len_trim(stem)>4 )then
                if( stem(len_trim(stem)-3:len_trim(stem))=='_001' ) stem=stem(:len_trim(stem)-4)
            endif
            prefix=string(stem)
            write(tag,'(I3.3)') state
            vol_fname = prefix//'_'//tag//MRC_EXT
        endif
        call img%write(vol_fname,del_if_exists=.true.)
        write(logfhandle,'(A,I0,A,A)') '>>> FLEX DIFFMAP NYSTROM PRE-IMAGE ',state,': ',vol_fname%to_char()
        call prefix%kill
        call ext%kill
    end subroutine flex_pca_write_state

    !> Pearson correlation of two double vectors.
    real(dp) function corr_dp( a, b, n ) result( r )
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

    ! Rank at which the Gram spectrum enters its noise bulk. Noise level = median of the lower half,
    ! so the leading signal directions cannot inflate it. Scale-free.
    pure integer function cov_signal_rank( eval, n ) result( d )
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
    subroutine cov_stage_subsample( proj_field, pinds, nptcls, nparts, maxtot, label, spinds, nsel )
        class(oris),          intent(in)  :: proj_field
        integer,              intent(in)  :: pinds(:), nptcls, nparts, maxtot
        character(len=*),     intent(in)  :: label
        integer, allocatable, intent(out) :: spinds(:)
        integer,              intent(out) :: nsel
        integer :: nmax_tot, nmax_part, ihalf, i, nkept, n_half, ntgt
        nmax_tot = maxtot
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
                if( proj_field%get_eo(pinds(i)) == ihalf ) n_half = n_half + 1
            end do
            if( n_half < 1 ) cycle
            ! split the per-part budget evenly between halfsets, never starving one
            ntgt  = min(n_half, max(1, (nmax_part + 1 - ihalf)/2))
            nkept = 0
            do i = 1, nptcls
                if( proj_field%get_eo(pinds(i)) /= ihalf ) cycle
                ! real(dp) rather than integer products: nkept*ntgt overflows int32 at these sizes
                if( int(real(nkept+1,dp)*real(ntgt,dp)/real(n_half,dp)) > &
                   &int(real(nkept,  dp)*real(ntgt,dp)/real(n_half,dp)) )then
                    nsel = nsel + 1
                    spinds(nsel) = pinds(i)
                endif
                nkept = nkept + 1
            end do
        end do
        if( nsel < 2 ) THROW_HARD('stage subsample left too few particles; raise the '//trim(label)//' particle cap')
        call hpsort(spinds(:nsel))   ! restore project order so batched image reads stay sequential
        if( nsel < nptcls )then
            write(logfhandle,'(A,A,A,I0,A,I0,A)') '>>> FLEX_PCA ',trim(label),' subsample: using ', &
                &nsel,' of ',nptcls,' particles'
            call flush(logfhandle)
        endif
    end subroutine cov_stage_subsample

    !>  Bytes of the packed [d(d+1)/2]^2 array model that cov_dim_budget sizes d against.
    pure real(dp) function cov_accum_bytes( d ) result( nbytes )
        integer, intent(in) :: d
        real(dp) :: n
        n = real(d,dp)*real(d+1,dp)/2.d0   ! Mspk(npk,npk), npk = d(d+1)/2
        nbytes = 8.d0*n*n
    end function cov_accum_bytes

    !> Largest d with cov_accum_bytes(d) <= FLEX_ACCUM_BYTE_BUDGET.
    pure integer function cov_dim_budget() result( d )
        ! d(d+1)/2 = sqrt(BUDGET/8)  =>  d = (-1 + sqrt(1 + 8*sqrt(BUDGET/8)))/2
        d = max(1, int((-1.d0 + sqrt(1.d0 + 8.d0*sqrt(FLEX_ACCUM_BYTE_BUDGET/8.d0)))/2.d0))
    end function cov_dim_budget

end module simple_flex_pca_util
