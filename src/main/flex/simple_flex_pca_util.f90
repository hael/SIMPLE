!@descr: flex_pca shared helpers: environment switches, chi-squared median, unimodality test
module simple_flex_pca_util
use simple_core_module_api
use simple_image, only: image
use simple_parameters, only: parameters
implicit none
private
#include "simple_local_flags.inc"

public :: cov_env_int, cov_env_int_pub, cov_env_flag_on, cov_env_flag_off, cov_env_dp
public :: chi2_median, punit, two_gauss_unimodal
public :: kernel_weights_at_bandwidth, project_onto_target_polyline, dilation_template
public :: flex_pca_write_state

!> safety cap on kernel bandwidth growth when a state's support falls below min_neff
integer, parameter :: COV_MAX_BW_GROW = 4


contains

    !>  Override an integer from the environment, if the variable is set and parses (values > 0 only).
    subroutine cov_env_int( name, val )
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

    !>  Override an integer from the environment, if the variable is set and parses.
    subroutine cov_env_int_pub( name, val )
        character(len=*), intent(in)    :: name
        integer,          intent(inout) :: val
        call cov_env_int(name, val)
    end subroutine cov_env_int_pub


    !>  True only when an environment flag is set to a nonzero integer (an opt-IN switch).
    logical function cov_env_flag_on( name ) result(on)
        character(len=*), intent(in) :: name
        character(len=32) :: envval
        integer :: stat, ln, ival
        on = .false.
        call get_environment_variable(name, envval, ln, stat)
        if( stat /= 0 .or. ln < 1 ) return
        read(envval(:ln), *, iostat=stat) ival
        if( stat == 0 ) on = ival /= 0
    end function cov_env_flag_on

    !>  True only when an environment flag is explicitly set to zero (an opt-OUT switch).
    logical function cov_env_flag_off( name ) result(off)
        character(len=*), intent(in) :: name
        character(len=32) :: envval
        integer :: stat, ln, ival
        off = .false.
        call get_environment_variable(name, envval, ln, stat)
        if( stat /= 0 .or. ln < 1 ) return
        read(envval(:ln), *, iostat=stat) ival
        if( stat == 0 ) off = ival == 0
    end function cov_env_flag_off

    !>  Real-valued environment override, leaving `val` untouched when unset. Companion to cov_env_int.
    subroutine cov_env_dp( name, val )
        character(len=*), intent(in)    :: name
        real(dp),         intent(inout) :: val
        character(len=32) :: envval
        integer  :: stat, ln_env
        real(dp) :: rval
        call get_environment_variable(name, envval, ln_env, stat)
        if( stat /= 0 .or. ln_env < 1 ) return
        read(envval(:ln_env), *, iostat=stat) rval
        if( stat == 0 )then
            val = rval
            write(logfhandle,'(A,A,A,ES12.4)') '>>> FLEX_PCA ',trim(name),' override: ',rval
            call flush(logfhandle)
        endif
    end subroutine cov_env_dp

    !> Median of chi-squared with k dof, Wilson-Hilferty: k*(1 - 2/(9k))^3. Good to 3 % at k=1, which is
    !! far inside the tolerance of a bandwidth FLOOR and needs no gamma inverse.
    pure real(dp) function chi2_median( k )
        integer, intent(in) :: k
        real(dp) :: kk
        kk = real(max(k,1),dp)
        chi2_median = kk * (1.d0 - 2.d0/(9.d0*kk))**3
    end function chi2_median

    !> deterministic centred unit-variance pseudo-random draw, so the tests do not depend on an RNG
    real(dp) function punit( n ) result( u )
        integer, intent(in) :: n
        real(dp) :: t
        t = sin(real(n,dp)*12.9898d0)*43758.5453d0
        u = t - floor(t)                               ! uniform(0,1)
        u = (2.d0*u - 1.d0)*sqrt(3.d0)                 ! centred, unit variance
    end function punit

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
        do grow = 0, COV_MAX_BW_GROW
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
            if( grow >= COV_MAX_BW_GROW      ) exit
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

    !> delivered state map name from params%outvol: state 1 keeps the name, others get _NNN
    subroutine flex_pca_write_state( params, img, state, vol_fname )
        class(parameters), intent(in)    :: params
        class(image),      intent(inout) :: img
        integer,           intent(in)    :: state
        class(string),     intent(inout) :: vol_fname
        type(string) :: prefix, ext
        character(len=:), allocatable :: stem
        character(len=3) :: tag
        if( state==1 )then
            vol_fname = params%outvol
        else
            ext=fname2ext(params%outvol)
            prefix=get_fbody(params%outvol,ext)
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

end module simple_flex_pca_util
