!@descr: the trajectory's noise floor (two-component mixture on the 5-frame median of the peak), its guard, and the floor normalization of the frame-quality bank
module single_frame_quality_floor
use simple_defs,                only: dp
use simple_error,               only: simple_exception
use simple_xd_gmm,              only: xd_gmm
use single_frame_quality_types, only: frame_quality_table, FQ_NFEATS, FQ_NLOC, FQ_FLOOR_K, FQ_GUARD_WEIGHT, &
    FQ_GUARD_SEP, I_LATTICE_PEAK, I_BAND_MAX, I_PEAK_MEDIAN, I_SPOT_COUNT
implicit none
private
#include "simple_local_flags.inc"

public :: fq_running_median, fq_fit_floor, fq_floor_normalize, fq_median, fq_percentile

integer,  parameter :: MAXITS_EM = 200
real(dp), parameter :: TOL_EM    = 1.d-8
real(dp), parameter :: PI_MIN_EM = 1.d-6
real,     parameter :: MAD_GAU   = 1.4826

contains

    !> running median over [i - hw, i + hw] clipped to the trajectory (the mean of the two middle values for an even
    !! count); hw = 0 returns x
    function fq_running_median( x, hw ) result( y )
        real,    intent(in) :: x(:)
        integer, intent(in) :: hw
        real    :: y(size(x)), win(2*max(hw,0)+1)
        integer :: n, i, j0, j1, m
        n = size(x)
        if( hw <= 0 )then
            y = x
            return
        endif
        do i = 1, n
            j0 = max(1, i - hw)
            j1 = min(n, i + hw)
            m  = j1 - j0 + 1
            win(1:m) = x(j0:j1)
            y(i) = fq_median(win(1:m))
        end do
    end function fq_running_median

    !> the two-component mixture of a one-dimensional sample, as fit_k of simple_nano_species: the better of the
    !! equal-mass quantile start and the largest-gap start; components in increasing mean, resp the floor's posterior
    subroutine fq_fit_floor( x, w, mu, sd, resp_floor )
        real, intent(in)  :: x(:)
        real, intent(out) :: w(2), mu(2), sd(2), resp_floor(size(x))
        real(dp), allocatable :: xx(:,:), rr(:,:,:), nz(:,:,:), post(:,:), post2(:,:)
        real(dp) :: means0(1,2), m(2), v(2), pk(2), lnl, m2(2), v2(2), pk2(2), lnl2
        real     :: xs(size(x))
        integer  :: n, g, i
        n = size(x)
        if( n < 4 ) THROW_HARD('too few frames for the floor; fq_fit_floor')
        allocate(xx(n,1), rr(1,1,n), nz(1,1,n), post(n,2), post2(n,2))
        xx(:,1) = real(x, dp)
        rr = 1._dp
        nz = 0._dp
        call fit_from(m, v, pk, post, lnl)
        xs = x
        call sort_asc(xs)
        g = 1
        do i = 2, n - 1
            if( xs(i+1) - xs(i) > xs(g+1) - xs(g) ) g = i
        end do
        means0(1,1) = sum(real(xs(1:g), dp)) / real(g, dp)
        means0(1,2) = sum(real(xs(g+1:n), dp)) / real(n - g, dp)
        call fit_from(m2, v2, pk2, post2, lnl2, means0)
        if( lnl2 > lnl )then
            m    = m2
            v    = v2
            pk   = pk2
            post = post2
        endif
        if( m(1) <= m(2) )then
            mu = real(m); sd = real(sqrt(v)); w = real(pk); resp_floor = real(post(:,1))
        else
            mu = real(m(2:1:-1)); sd = real(sqrt(v(2:1:-1))); w = real(pk(2:1:-1)); resp_floor = real(post(:,2))
        endif

        contains

            subroutine fit_from( m_o, v_o, pk_o, post_o, lnl_o, means0_o )
                real(dp),           intent(out) :: m_o(2), v_o(2), pk_o(2), post_o(:,:), lnl_o
                real(dp), optional, intent(in)  :: means0_o(1,2)
                type(xd_gmm) :: gm
                real(dp) :: mm(1,2), ss(1,1,2), xhat(n,1), xcov(1,1,n)
                call gm%new(1, 2, MAXITS_EM, TOL_EM, PI_MIN_EM)
                if( present(means0_o) )then
                    call gm%init(xx, nz, means0_o)
                else
                    call gm%init(xx, nz)
                endif
                call gm%fit(xx, rr, nz)
                lnl_o = gm%get_loglik()
                call gm%get_means(mm)
                call gm%get_covs(ss)
                call gm%get_pi(pk_o)
                call gm%posterior(xx, rr, nz, xhat, xcov, post_o)
                call gm%kill
                m_o = mm(1,:)
                v_o = ss(1,1,:)
            end subroutine fit_from

    end subroutine fq_fit_floor

    !> peak_median, the floors of lattice_peak and band_max, the guard, spot_count and the normalized features
    !! (trajectory_model_frames_rejection.md section 4); the raw lattice_peak, band_max, local maxima and the other raw features
    !! must be in place. Without a floor, every feature is normalized by its trajectory median and MAD.
    !! The floor that defines the trajectory (the guard, the spot threshold, peak_median's scale) is fitted on the
    !! 5-frame median of lattice_peak; the single-frame lattice_peak is normalized by its own, wider floor, because a
    !! frame's own peak decides too: lattice visibility flickers frame to frame as the particle tumbles, and a 21-frame
    !! median alone rejects the flashes (seen on 71-17/traj02, 2026-10-27).
    subroutine fq_floor_normalize( tab )
        type(frame_quality_table), intent(inout) :: tab
        real    :: w(2), mu(2), sd(2), wb(2), mub(2), sdb(2), ws(2), mus(2), sds(2), med, mad, mins, scale
        real, allocatable :: resp(:), respb(:), resps(:), vals(:)
        logical, allocatable :: ok(:)
        integer :: i, ifeat
        allocate(resp(tab%n), respb(tab%n), resps(tab%n))
        ok = .not. tab%invalid
        tab%raw(:,I_PEAK_MEDIAN) = fq_running_median(tab%raw(:,I_LATTICE_PEAK), tab%winsz)
        call fq_fit_floor(fq_running_median(tab%raw(:,I_LATTICE_PEAK), 2), w, mu, sd, resp)
        call fq_fit_floor(tab%raw(:,I_LATTICE_PEAK), ws, mus, sds, resps)
        call fq_fit_floor(fq_running_median(tab%raw(:,I_BAND_MAX), 2), wb, mub, sdb, respb)
        tab%floor_mean    = mu(1)
        tab%floor_sd      = sd(1)
        tab%floor_weight  = w(1)
        tab%upper_mean    = mu(2)
        tab%upper_sd      = sd(2)
        tab%sf_floor_mean = mus(1)
        tab%sf_floor_sd   = sds(1)
        tab%bm_floor_mean = mub(1)
        tab%bm_floor_sd   = sdb(1)
        tab%floor_frames  = resp > 0.5
        tab%has_floor     = w(1) >= FQ_GUARD_WEIGHT .and. (mu(2) - mu(1)) >= FQ_GUARD_SEP * sd(1)
        tab%spot_threshold = mu(1) + FQ_FLOOR_K * sd(1)
        do i = 1, tab%n
            tab%raw(i,I_SPOT_COUNT) = real(count(tab%locmax(i,:) > 10.**tab%spot_threshold))
        end do
        if( tab%has_floor )then
            tab%z(:,I_LATTICE_PEAK) = (tab%raw(:,I_LATTICE_PEAK) - mus(1)) / sds(1)
            tab%z(:,I_PEAK_MEDIAN)  = (tab%raw(:,I_PEAK_MEDIAN)  - mu(1)) / sd(1)
            tab%z(:,I_BAND_MAX)     = (tab%raw(:,I_BAND_MAX)     - mub(1)) / sdb(1)
            do ifeat = 1, FQ_NFEATS
                if( ifeat == I_LATTICE_PEAK .or. ifeat == I_PEAK_MEDIAN .or. ifeat == I_BAND_MAX ) cycle
                vals = pack(tab%raw(:,ifeat), tab%floor_frames .and. ok)
                med  = fq_median(vals)
                mad  = MAD_GAU * fq_median(abs(vals - med))
                if( ifeat == I_SPOT_COUNT )then
                    mins = 1.
                else
                    vals = pack(tab%raw(:,ifeat), ok)
                    mins = 0.1 * (fq_percentile(vals, 75.) - fq_percentile(vals, 25.))
                endif
                scale = max(mad, mins)
                call set_z(ifeat, med, scale)
            end do
        else
            do ifeat = 1, FQ_NFEATS
                vals = pack(tab%raw(:,ifeat), ok)
                med  = fq_median(vals)
                mad  = MAD_GAU * fq_median(abs(vals - med))
                call set_z(ifeat, med, mad)
            end do
        endif

        contains

            ! a feature without spread carries nothing: zero
            subroutine set_z( jfeat, center, sc )
                integer, intent(in) :: jfeat
                real,    intent(in) :: center, sc
                if( sc > tiny(sc) )then
                    tab%z(:,jfeat) = (tab%raw(:,jfeat) - center) / sc
                else
                    tab%z(:,jfeat) = 0.
                endif
            end subroutine set_z

    end subroutine fq_floor_normalize

    !> the median (the mean of the two middle values for an even count)
    real function fq_median( x )
        real, intent(in) :: x(:)
        real    :: xs(size(x))
        integer :: m
        m = size(x)
        if( m == 0 )then
            fq_median = 0.
            return
        endif
        xs = x
        call sort_asc(xs)
        if( mod(m, 2) == 1 )then
            fq_median = xs(m/2 + 1)
        else
            fq_median = 0.5 * (xs(m/2) + xs(m/2 + 1))
        endif
    end function fq_median

    !> the q-th percentile by linear interpolation between order statistics (numpy's default)
    real function fq_percentile( x, q )
        real, intent(in) :: x(:), q
        real    :: xs(size(x)), pos
        integer :: m, lo
        m = size(x)
        if( m == 0 )then
            fq_percentile = 0.
            return
        endif
        xs = x
        call sort_asc(xs)
        pos = real(m - 1) * q / 100.
        lo  = min(m - 1, int(pos))
        if( m == 1 )then
            fq_percentile = xs(1)
        else
            fq_percentile = xs(lo + 1) + (pos - real(lo)) * (xs(lo + 2) - xs(lo + 1))
        endif
    end function fq_percentile

    !> ascending heapsort
    subroutine sort_asc( v )
        real, intent(inout) :: v(:)
        integer :: n, i
        real    :: t
        n = size(v)
        do i = n / 2, 1, -1
            call sift(i, n)
        end do
        do i = n, 2, -1
            t = v(1); v(1) = v(i); v(i) = t
            call sift(1, i - 1)
        end do

        contains

            subroutine sift( start, last )
                integer, intent(in) :: start, last
                integer :: root, child
                real    :: tmp
                root = start
                do while( 2 * root <= last )
                    child = 2 * root
                    if( child < last )then
                        if( v(child) < v(child + 1) ) child = child + 1
                    endif
                    if( v(root) >= v(child) ) exit
                    tmp = v(root); v(root) = v(child); v(child) = tmp
                    root = child
                end do
            end subroutine sift

    end subroutine sort_asc

end module single_frame_quality_floor
