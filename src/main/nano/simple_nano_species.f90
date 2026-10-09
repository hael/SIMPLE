!@descr: species discovery numerics for SINGLE: intensity mixture and class choice, enclosed Gaussian fraction, detection-threshold calibration, Gaussian width fit, filters
! Array numerics with no image dependence; the design is doc/implementation_notes/planned/species_discovery.md
! (sections 3.2 to 3.7).
module simple_nano_species
use simple_defs,          only: dp, DPI
use simple_error,         only: simple_exception
use simple_srch_sort_loc, only: hpsort
use simple_stat,          only: median
use simple_xd_gmm,        only: xd_gmm
implicit none

public :: fit_species_mixture, class_separation, enclosed_fraction, maxima_shape, calibrate_threshold, expected_false_count
public :: fit_gauss_width, gauss_filter3D, local_maxima, robust_spread
public :: MAX_NSPECIES, MIN_SEPARATION, MIN_CLASS_ATOMS, MIN_CLASS_FRAC
private
#include "simple_local_flags.inc"

integer,  parameter :: MAX_NSPECIES    = 3      ! largest number of classes tried
real,     parameter :: MIN_SEPARATION  = 3.     ! adjacent classes need a separation D of at least this
integer,  parameter :: MIN_CLASS_ATOMS = 8      ! a class holds at least max(MIN_CLASS_ATOMS, MIN_CLASS_FRAC * N) atoms
real,     parameter :: MIN_CLASS_FRAC  = 0.02
integer,  parameter :: MAXITS_EM       = 1000
real(dp), parameter :: TOL_EM          = 1.d-10 ! relative change of the log-likelihood that ends EM
real(dp), parameter :: PI_MIN_EM       = 1.d-6  ! mixing floor of the deconvolution fit: a component below it freezes
real,     parameter :: U_MAX           = 50.    ! upper end of the threshold search
real,     parameter :: BFAC_MIN        = 1.     ! range of the width search, B in A**2
real,     parameter :: BFAC_MAX        = 200.
integer,  parameter :: NGRID_B         = 100    ! log-spaced grid of the width search before golden-section refinement
integer,  parameter :: NGOLD_B         = 40

contains

    ! One-dimensional Gaussian mixture of the intensities x for K = 1..kmax, fitted as an extreme deconvolution
    ! (simple_xd_gmm): every intensity is a class value plus measurement noise of spread s_meas, so the class
    ! variances are the intrinsic spreads and the observed variance of a class is var_k + s_meas**2 (var returns
    ! the observed one). Best of two deterministic starts, sorted quantiles and largest gaps. With nspecies_in
    ! positive, kmax = K = nspecies_in; otherwise kmax = MAX_NSPECIES and K is the admissible K of lowest BIC.
    ! Classes are numbered by decreasing mean; post(i,k) is the posterior. A K fixed by nspecies_in is taken even
    ! when inadmissible, so a class may then hold no atom. bic and admissible are of size kmax.
    subroutine fit_species_mixture( x, s_meas, nspecies_in, K, labels, post, mu, var, bic, admissible )
        real,                           intent(in)  :: x(:)
        real,                           intent(in)  :: s_meas
        integer,                        intent(in)  :: nspecies_in
        integer,                        intent(out) :: K
        integer,                        intent(out) :: labels(size(x))
        real,              allocatable, intent(out) :: post(:,:), mu(:), var(:)
        real,    optional, allocatable, intent(out) :: bic(:)
        logical, optional, allocatable, intent(out) :: admissible(:)
        real(dp), allocatable :: mus(:,:), vars(:,:), ws(:,:), posts(:,:,:), bics(:)
        logical,  allocatable :: adm(:)
        real(dp) :: X1(size(x),1), R(1,1,size(x)), Nz(1,1,size(x)), lnl
        integer  :: n, kk, kmax
        n = size(x)
        if( n < 1 )               THROW_HARD('no intensities; fit_species_mixture')
        if( nspecies_in < 0 )     THROW_HARD('nspecies must be 0 (automatic) or positive; fit_species_mixture')
        if( nspecies_in > n )     THROW_HARD('fewer intensities than nspecies; fit_species_mixture')
        if( s_meas <= 0. )        THROW_HARD('measurement spread must be positive; fit_species_mixture')
        X1(:,1) = real(x, dp)
        R       = 1._dp
        Nz      = real(s_meas, dp)**2
        if( nspecies_in > 0 )then
            kmax = nspecies_in
        else
            kmax = min(MAX_NSPECIES, n)
        endif
        allocate(mus(kmax,kmax), vars(kmax,kmax), ws(kmax,kmax), posts(n,kmax,kmax), source=0._dp)
        allocate(bics(kmax), source=huge(1._dp))
        allocate(adm(kmax),  source=.false.)
        do kk = 1,kmax
            call fit_k(X1, R, Nz, kk, mus(:kk,kk), vars(:kk,kk), ws(:kk,kk), posts(:,:kk,kk), lnl)
            bics(kk) = -2._dp * lnl + real(3*kk-1, dp) * log(real(n, dp))
            adm(kk)  = is_admissible(posts(:,:kk,kk), mus(:kk,kk), vars(:kk,kk))
        enddo
        if( nspecies_in > 0 )then
            K = nspecies_in
        else
            K = 1
            do kk = 2,kmax
                if( adm(kk) .and. bics(kk) < bics(K) ) K = kk
            enddo
        endif
        labels = maxloc(posts(:,:K,K), dim=2)
        post   = real(posts(:,:K,K))
        mu     = real(mus(:K,K))
        var    = real(vars(:K,K))
        if( present(bic) )        bic        = real(min(bics, real(huge(1.), dp)))
        if( present(admissible) ) admissible = adm
    end subroutine fit_species_mixture

    ! D = |mu_a - mu_b| / sqrt((var_a + var_b) / 2)
    elemental real function class_separation( mu_a, var_a, mu_b, var_b )
        real, intent(in) :: mu_a, var_a, mu_b, var_b
        class_separation = abs(mu_a - mu_b) / sqrt(0.5 * (var_a + var_b))
    end function class_separation

    ! fraction of an isotropic 3D Gaussian inside t standard deviations
    elemental real function enclosed_fraction( t )
        real, intent(in) :: t
        if( t <= 0. )then
            enclosed_fraction = 0.
        else
            enclosed_fraction = erf(t / sqrt(2.)) - sqrt(2. / real(DPI)) * t * exp(-0.5 * t * t)
        endif
    end function enclosed_fraction

    ! u dependence of the expected number of maxima above u of a smooth Gaussian field in 3D
    ! (Euler characteristic density): C (u**2 - 1) exp(-u**2 / 2)
    elemental real function maxima_shape( u )
        real, intent(in) :: u
        maxima_shape = (u * u - 1.) * exp(-0.5 * u * u)
    end function maxima_shape

    ! counts(i) maxima above u(i) in the noise region give C = sum(counts) / sum(maxima_shape(u)); scaled by
    ! the search-to-region volume ratio, k is the threshold at which target noise maxima are expected
    subroutine calibrate_threshold( u, counts, vol_ratio, target, k, c_search )
        real, intent(in)  :: u(:), counts(:)
        real, intent(in)  :: vol_ratio, target
        real, intent(out) :: k, c_search
        real    :: lo, hi, mid
        integer :: it
        if( size(u) < 1 .or. size(u) /= size(counts) ) THROW_HARD('levels and counts differ in size; calibrate_threshold')
        if( any(u <= 1.) )                             THROW_HARD('calibration levels must exceed 1; calibrate_threshold')
        if( any(counts < 0.) .or. sum(counts) <= 0. )  THROW_HARD('no region maxima above the calibration levels; calibrate_threshold')
        if( vol_ratio <= 0. .or. target <= 0. )        THROW_HARD('volume ratio and target must be positive; calibrate_threshold')
        c_search = vol_ratio * sum(counts) / sum(maxima_shape(u))
        ! maxima_shape falls monotonically above sqrt(3)
        lo = sqrt(3.)
        hi = U_MAX
        if( c_search * maxima_shape(lo) <= target )then
            k = lo
            return
        endif
        do it = 1,60
            mid = 0.5 * (lo + hi)
            if( c_search * maxima_shape(mid) > target )then
                lo = mid
            else
                hi = mid
            endif
        enddo
        k = 0.5 * (lo + hi)
    end subroutine calibrate_threshold

    elemental real function expected_false_count( c_search, k )
        real, intent(in) :: c_search, k
        expected_false_count = c_search * maxima_shape(k)
    end function expected_false_count

    ! Amplitude A and B factor of A exp(-4 pi**2 r2 / B) fitted to the samples y at squared distances r2 (A**2),
    ! noise variance sig2n. Free amplitude with the prior ((ln B - lnb_mean) / lnb_sd)**2 (stage 1), or, with
    ! i_class and tau_class, the amplitude tied to the class: ((A (B / 4 pi)**1.5 - i_class) / tau_class)**2 (stage 2).
    subroutine fit_gauss_width( y, r2, sig2n, lnb_mean, lnb_sd, amp, bfac, i_class, tau_class )
        real,           intent(in)  :: y(:), r2(:)
        real,           intent(in)  :: sig2n, lnb_mean, lnb_sd
        real,           intent(out) :: amp, bfac
        real, optional, intent(in)  :: i_class, tau_class
        real(dp), parameter :: GOLD = (sqrt(5._dp) - 1._dp) / 2._dp
        real(dp) :: t(NGRID_B), c(NGRID_B), lo, hi, x1, x2, f1, f2, syy
        logical  :: l_tied
        integer  :: i, ibest
        if( size(y) < 1 .or. size(y) /= size(r2) ) THROW_HARD('samples and distances differ in size; fit_gauss_width')
        l_tied = present(i_class) .and. present(tau_class)
        syy    = sum(real(y, dp)**2)
        do i = 1,NGRID_B
            t(i) = log(real(BFAC_MIN, dp)) + real(i-1, dp) * log(real(BFAC_MAX/BFAC_MIN, dp)) / real(NGRID_B-1, dp)
            c(i) = cost(t(i))
        enddo
        ibest = minloc(c, dim=1)
        lo = t(max(ibest-1, 1))
        hi = t(min(ibest+1, NGRID_B))
        x1 = hi - GOLD * (hi - lo)
        x2 = lo + GOLD * (hi - lo)
        f1 = cost(x1)
        f2 = cost(x2)
        do i = 1,NGOLD_B
            if( f1 < f2 )then
                hi = x2
                x2 = x1
                f2 = f1
                x1 = hi - GOLD * (hi - lo)
                f1 = cost(x1)
            else
                lo = x1
                x1 = x2
                f1 = f2
                x2 = lo + GOLD * (hi - lo)
                f2 = cost(x2)
            endif
        enddo
        if( min(f1, f2) > c(ibest) )then
            bfac = real(exp(t(ibest)))
        else
            bfac = real(exp(0.5_dp * (lo + hi)))
        endif
        amp = real(amp_at(log(real(bfac, dp))))

    contains

        ! amplitude that minimises the objective at ln B = tb
        real(dp) function amp_at( tb )
            real(dp), intent(in) :: tb
            real(dp) :: g(size(y)), sgg, syg, cb
            g   = exp(-4._dp * DPI**2 * real(r2, dp) / exp(tb))
            sgg = sum(g * g)
            syg = sum(real(y, dp) * g)
            if( l_tied )then
                cb     = (exp(tb) / (4._dp * DPI))**1.5_dp
                amp_at = (syg / sig2n + cb * i_class / tau_class**2) / (sgg / sig2n + cb * cb / tau_class**2)
            else
                amp_at = syg / max(sgg, tiny(1._dp))
            endif
            amp_at = max(amp_at, 0._dp)
        end function amp_at

        real(dp) function cost( tb )
            real(dp), intent(in) :: tb
            real(dp) :: g(size(y)), sgg, syg, a, cb
            g    = exp(-4._dp * DPI**2 * real(r2, dp) / exp(tb))
            sgg  = sum(g * g)
            syg  = sum(real(y, dp) * g)
            a    = amp_at(tb)
            cost = (syy - 2._dp * a * syg + a * a * sgg) / sig2n
            if( l_tied )then
                cb   = (exp(tb) / (4._dp * DPI))**1.5_dp
                cost = cost + ((a * cb - i_class) / tau_class)**2
            else
                cost = cost + ((tb - lnb_mean) / lnb_sd)**2
            endif
        end function cost

    end subroutine fit_gauss_width

    ! b = a convolved with the normalised 3D Gaussian of standard deviation sigma_vox voxels, zero outside the box
    subroutine gauss_filter3D( a, sigma_vox, b )
        real, intent(in)  :: a(:,:,:)
        real, intent(in)  :: sigma_vox
        real, intent(out) :: b(:,:,:)
        real, allocatable :: tmp(:,:,:), w(:)
        integer :: n(3), hw, i, j, k, l
        n  = shape(a)
        hw = max(1, ceiling(4. * sigma_vox))
        allocate(w(-hw:hw))
        do l = -hw,hw
            w(l) = exp(-0.5 * (real(l) / sigma_vox)**2)
        enddo
        w = w / sum(w)
        allocate(tmp(n(1),n(2),n(3)))
        !$omp parallel do default(shared) private(i,j,k,l) schedule(static) proc_bind(close)
        do k = 1,n(3)
            do j = 1,n(2)
                do i = 1,n(1)
                    b(i,j,k) = 0.
                    do l = max(-hw, 1-i), min(hw, n(1)-i)
                        b(i,j,k) = b(i,j,k) + w(l) * a(i+l,j,k)
                    enddo
                enddo
            enddo
        enddo
        !$omp end parallel do
        !$omp parallel do default(shared) private(i,j,k,l) schedule(static) proc_bind(close)
        do k = 1,n(3)
            do j = 1,n(2)
                do i = 1,n(1)
                    tmp(i,j,k) = 0.
                    do l = max(-hw, 1-j), min(hw, n(2)-j)
                        tmp(i,j,k) = tmp(i,j,k) + w(l) * b(i,j+l,k)
                    enddo
                enddo
            enddo
        enddo
        !$omp end parallel do
        !$omp parallel do default(shared) private(i,j,k,l) schedule(static) proc_bind(close)
        do k = 1,n(3)
            do j = 1,n(2)
                do i = 1,n(1)
                    b(i,j,k) = 0.
                    do l = max(-hw, 1-k), min(hw, n(3)-k)
                        b(i,j,k) = b(i,j,k) + w(l) * tmp(i,j,k+l)
                    enddo
                enddo
            enddo
        enddo
        !$omp end parallel do
    end subroutine gauss_filter3D

    ! voxels of z inside mask that exceed all 26 neighbours and thres, with at least nmin_above neighbours
    ! above thres; ijk(3,n) and vals(n) by decreasing value
    subroutine local_maxima( z, mask, thres, nmin_above, ijk, vals )
        real,                 intent(in)  :: z(:,:,:)
        logical,              intent(in)  :: mask(:,:,:)
        real,                 intent(in)  :: thres
        integer,              intent(in)  :: nmin_above
        integer, allocatable, intent(out) :: ijk(:,:)
        real,    allocatable, intent(out) :: vals(:)
        integer, allocatable :: tmp_ijk(:,:), order(:)
        real,    allocatable :: tmp_vals(:)
        integer :: n(3), i, j, k, cnt, nabove, nmax
        logical :: l_max
        n    = shape(z)
        nmax = count(mask .and. z > thres)
        allocate(tmp_ijk(3,max(nmax,1)), tmp_vals(max(nmax,1)))
        cnt = 0
        do k = 2,n(3)-1
            do j = 2,n(2)-1
                do i = 2,n(1)-1
                    if( .not. mask(i,j,k) ) cycle
                    if( z(i,j,k) <= thres ) cycle
                    l_max  = z(i,j,k) > maxval(z(i-1:i+1,j-1:j+1,k-1:k+1), mask=neighbour_mask())
                    if( .not. l_max ) cycle
                    nabove = count(z(i-1:i+1,j-1:j+1,k-1:k+1) > thres) - 1
                    if( nabove < nmin_above ) cycle
                    cnt = cnt + 1
                    tmp_ijk(:,cnt) = [i,j,k]
                    tmp_vals(cnt)  = z(i,j,k)
                enddo
            enddo
        enddo
        allocate(ijk(3,cnt), vals(cnt), order(cnt))
        if( cnt == 0 ) return
        vals  = -tmp_vals(:cnt)
        order = [(i, i=1,cnt)]
        call hpsort(vals, order)
        vals = -vals
        ijk  = tmp_ijk(:,order)

    contains

        pure function neighbour_mask() result( m )
            logical :: m(3,3,3)
            m        = .true.
            m(2,2,2) = .false.
        end function neighbour_mask

    end subroutine local_maxima

    ! 1.4826 times the median absolute deviation from the median
    real function robust_spread( x )
        real, intent(in) :: x(:)
        real :: med
        med           = median(x)
        robust_spread = 1.4826 * median(abs(x - med))
    end function robust_spread

    ! PRIVATE

    ! Extreme-deconvolution fits of kk classes from the sorted-quantile and the largest-gap starts; the fit of
    ! higher likelihood, by decreasing mean. var is the observed class variance, intrinsic plus the noise.
    subroutine fit_k( X1, R, Nz, kk, mu, var, w, post, lnl )
        real(dp), intent(in)  :: X1(:,:), R(:,:,:), Nz(:,:,:)
        integer,  intent(in)  :: kk
        real(dp), intent(out) :: mu(kk), var(kk), w(kk), post(size(X1,1),kk), lnl
        real(dp) :: mu2(kk), var2(kk), w2(kk), post2(size(X1,1),kk), lnl2, means0(1,kk)
        real     :: xs(size(X1,1)), gaps(size(X1,1)-1)
        integer  :: idx(size(X1,1)), gidx(size(X1,1)-1), bnd(kk+1), n, i, g, order(kk)
        n = size(X1,1)
        call fit_from(mu, var, w, post, lnl)
        if( kk > 1 )then
            ! the kk-1 largest gaps of the sorted sample bound the start groups
            xs  = real(X1(:,1))
            idx = [(i, i=1,n)]
            call hpsort(xs, idx)
            gaps = xs(2:n) - xs(1:n-1)
            gidx = [(i, i=1,n-1)]
            call hpsort(gaps, gidx)
            bnd(1)    = 0
            bnd(kk+1) = n
            bnd(2:kk) = gidx(n-kk+1:n-1)
            call hpsort(bnd(2:kk))
            do g = 1,kk
                means0(1,g) = sum(X1(idx(bnd(g)+1:bnd(g+1)),1)) / real(bnd(g+1) - bnd(g), dp)
            enddo
            call fit_from(mu2, var2, w2, post2, lnl2, means0)
            if( lnl2 > lnl )then
                mu   = mu2
                var  = var2
                w    = w2
                post = post2
                lnl  = lnl2
            endif
        endif
        ! decreasing mean
        order = sort_desc(mu)
        mu   = mu(order)
        var  = var(order)
        w    = w(order)
        post = post(:,order)

        contains

            subroutine fit_from( mu_o, var_o, w_o, post_o, lnl_o, means0_o )
                real(dp),           intent(out) :: mu_o(kk), var_o(kk), w_o(kk), post_o(size(X1,1),kk), lnl_o
                real(dp), optional, intent(in)  :: means0_o(1,kk)
                type(xd_gmm) :: g
                real(dp) :: m(1,kk), S(1,1,kk), xhat(size(X1,1),1), xcov(1,1,size(X1,1))
                call g%new(1, kk, MAXITS_EM, TOL_EM, PI_MIN_EM)
                if( present(means0_o) )then
                    call g%init(X1, Nz, means0_o)
                else
                    call g%init(X1, Nz)
                endif
                call g%fit(X1, R, Nz)
                lnl_o = g%get_loglik()
                call g%get_means(m)
                call g%get_covs(S)
                call g%get_pi(w_o)
                call g%posterior(X1, R, Nz, xhat, xcov, post_o)
                call g%kill
                mu_o  = m(1,:)
                var_o = S(1,1,:) + Nz(1,1,1)
            end subroutine fit_from

    end subroutine fit_k

    ! permutation that sorts a short array into decreasing order
    function sort_desc( v ) result( p )
        real(dp), intent(in) :: v(:)
        integer :: p(size(v)), i, j, t
        p = [(i, i=1,size(v))]
        do i = 2,size(v)
            j = i
            do while( j > 1 )
                if( v(p(j)) <= v(p(j-1)) ) exit
                t = p(j); p(j) = p(j-1); p(j-1) = t
                j = j - 1
            enddo
        enddo
    end function sort_desc

    ! every class holds enough atoms (by largest posterior) and adjacent classes are separated by D >= MIN_SEPARATION
    logical function is_admissible( post, mu, var )
        real(dp), intent(in) :: post(:,:), mu(:), var(:)
        integer :: lab(size(post,1)), kk, n
        n = size(post,1)
        is_admissible = .true.
        if( size(mu) == 1 ) return
        lab = maxloc(post, dim=2)
        do kk = 1,size(mu)
            if( real(count(lab == kk)) < max(real(MIN_CLASS_ATOMS), MIN_CLASS_FRAC * real(n)) ) is_admissible = .false.
        enddo
        do kk = 1,size(mu)-1
            if( class_separation(real(mu(kk)), real(var(kk)), real(mu(kk+1)), real(var(kk+1))) < MIN_SEPARATION )then
                is_admissible = .false.
            endif
        enddo
    end function is_admissible

end module simple_nano_species
