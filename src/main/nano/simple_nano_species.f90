!@descr: species numerics for SINGLE: the class mixture with a shared coordination term and its posteriors, enclosed Gaussian fraction, detection-threshold calibration, Gaussian width fit, local maxima
! Array numerics with no image dependence; the design is doc/implementation_notes/completed/species_discovery_completed.md
! (sections 3.2 to 3.7).
module simple_nano_species
use simple_defs,          only: dp, DPI, DTINY
use simple_error,         only: simple_exception
use simple_srch_sort_loc, only: hpsort
use simple_xd_gmm,        only: xd_gmm
implicit none

public :: fit_species_mixture, class_posterior, class_separation, enclosed_fraction, maxima_shape, calibrate_threshold
public :: expected_false_count
public :: fit_gauss_width, local_maxima
public :: MAX_NSPECIES, MIN_SEPARATION, MIN_CLASS_ATOMS, MIN_CLASS_FRAC
private
#include "simple_local_flags.inc"

integer,  parameter :: MAX_NSPECIES    = 3      ! largest number of classes tried
real,     parameter :: MIN_SEPARATION  = 3.     ! adjacent classes need a separation D of at least this
integer,  parameter :: MIN_CLASS_ATOMS = 8      ! a class holds at least max(MIN_CLASS_ATOMS, MIN_CLASS_FRAC * N) atoms
real,     parameter :: MIN_CLASS_FRAC  = 0.02
integer,  parameter :: MAXITS_EM       = 200    ! a fit still moving after this many iterations is crawling on a flat likelihood
real(dp), parameter :: TOL_EM          = 1.d-8  ! relative change of the log-likelihood that ends EM
real(dp), parameter :: PI_MIN_EM       = 1.d-6  ! mixing floor of the deconvolution fit: a component below it freezes
real,     parameter :: U_MAX           = 50.    ! upper end of the threshold search
real,     parameter :: BFAC_MIN        = 1.     ! range of the width search, B in A**2
real,     parameter :: BFAC_MAX        = 200.
integer,  parameter :: NGRID_B         = 100    ! log-spaced grid of the width search before golden-section refinement
integer,  parameter :: NGOLD_B         = 40
integer,  parameter :: NITER_SLOPE     = 20     ! alternations of the class fit and the shared covariate slope
real(dp), parameter :: SLOPE_TOL_REL   = 1.d-3  ! the slope has converged when it moves by less than this mean intensity
real(dp), parameter :: TOL_SLOPE_EM    = 1.d-4  ! EM tolerance of the class fits inside the alternation
real(dp), parameter :: TOL_SLOPE_FIN   = 1.d-6  ! EM tolerance of the class fit at the slope chosen

contains

    ! One-dimensional Gaussian mixture of the values x for K = 1..kmax, fitted as an extreme deconvolution
    ! (simple_xd_gmm): every value is a class value plus noise of its own variance s2(i), so the class variances var
    ! are the intrinsic spreads beyond it (a value's class likelihood has the variance var_k + s2(i)). Best of two
    ! deterministic starts, sorted quantiles and largest gaps. With nspecies_in positive, kmax = K = nspecies_in;
    ! otherwise kmax = MAX_NSPECIES and K is the admissible K of lowest BIC. Classes are numbered by decreasing mean;
    ! post(i,k) is the posterior, w the class fractions. A K fixed by nspecies_in is taken even when inadmissible, so a
    ! class may then hold no atom. bic and admissible are of size kmax. With cov, every value carries one term
    ! slope * cov(i) shared by the classes (3.6): mu is then the class value at cov = 0, the slope is fitted for each K
    ! and counts as one more parameter, and slope returns that of K. detect_atomic_species fits x = ln I with
    ! cov = c - 12, the multiplicative form ln I = ln S_k + beta (c - 12) (atomic_density_model.md, section 2)
    subroutine fit_species_mixture( x, s2, nspecies_in, K, labels, post, mu, var, bic, admissible, cov, slope, w )
        real,                           intent(in)  :: x(:)
        real,                           intent(in)  :: s2(:)
        integer,                        intent(in)  :: nspecies_in
        integer,                        intent(out) :: K
        integer,                        intent(out) :: labels(size(x))
        real,              allocatable, intent(out) :: post(:,:), mu(:), var(:)
        real,    optional, allocatable, intent(out) :: bic(:)
        logical, optional, allocatable, intent(out) :: admissible(:)
        real,    optional,              intent(in)  :: cov(:)
        real,    optional,              intent(out) :: slope
        real,    optional, allocatable, intent(out) :: w(:)
        real(dp), allocatable :: mus(:,:), vars(:,:), ws(:,:), posts(:,:,:), bics(:), slopes(:)
        logical,  allocatable :: adm(:)
        real(dp) :: X1(size(x),1), R(1,1,size(x)), Nz(1,1,size(x)), lnl
        integer  :: n, kk, kmax, npar
        n = size(x)
        if( n < 1 )               THROW_HARD('no intensities; fit_species_mixture')
        if( nspecies_in < 0 )     THROW_HARD('nspecies must be 0 (automatic) or positive; fit_species_mixture')
        if( nspecies_in > n )     THROW_HARD('fewer intensities than nspecies; fit_species_mixture')
        if( size(s2) /= n )       THROW_HARD('one noise variance per value; fit_species_mixture')
        if( any(.not. s2 > 0.) )  THROW_HARD('noise variances must be positive; fit_species_mixture')
        if( present(cov) )then
            if( size(cov) /= n )  THROW_HARD('one covariate per intensity; fit_species_mixture')
        endif
        X1(:,1) = real(x, dp)
        R       = 1._dp
        Nz(1,1,:) = real(s2, dp)
        if( nspecies_in > 0 )then
            kmax = nspecies_in
        else
            kmax = min(MAX_NSPECIES, n)
        endif
        allocate(mus(kmax,kmax), vars(kmax,kmax), ws(kmax,kmax), posts(n,kmax,kmax), source=0._dp)
        allocate(bics(kmax), source=huge(1._dp))
        allocate(slopes(kmax), source=0._dp)
        allocate(adm(kmax),  source=.false.)
        do kk = 1,kmax
            if( present(cov) )then
                call fit_k_cov(X1, real(cov, dp), R, Nz, kk, mus(:kk,kk), vars(:kk,kk), ws(:kk,kk), posts(:,:kk,kk),&
                    &lnl, slopes(kk))
                npar = 3 * kk
            else
                call fit_k(X1, R, Nz, kk, mus(:kk,kk), vars(:kk,kk), ws(:kk,kk), posts(:,:kk,kk), lnl)
                npar = 3 * kk - 1
            endif
            bics(kk) = -2._dp * lnl + real(npar, dp) * log(real(n, dp))
            ! admissibility on the observed spreads: the intrinsic one and the typical noise
            adm(kk)  = is_admissible(posts(:,:kk,kk), mus(:kk,kk), vars(:kk,kk) + median_dp(Nz(1,1,:)))
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
        if( present(slope) )      slope      = real(slopes(K))
        if( present(w) )          w          = real(ws(:K,K))
    end subroutine fit_species_mixture

    ! the class posterior of one value x of noise variance nz under classes of means m (the shared covariate term
    ! included), intrinsic variances var and fractions w: post(k) ~ w_k N(x; m_k, var_k + nz); and the standard error
    ! of the posterior of the most probable class from the measurement noise v of x (part of nz), by propagation:
    ! sd = sqrt(v) |dP_k/dx|, dP_k/dx = P_k (g_k - sum_j P_j g_j), g_j = (m_j - x) / (var_j + nz). For two classes of
    ! one variance s2 this is P (1 - P) |m_1 - m_2| / s2 sqrt(v)
    subroutine class_posterior( x, nz, v, m, var, w, post, label, sd_post )
        real,    intent(in)  :: x, nz, v, m(:), var(:), w(:)
        real,    intent(out) :: post(size(m)), sd_post
        integer, intent(out) :: label
        real(dp) :: lp(size(m)), g(size(m)), p(size(m)), s2k, gbar
        integer  :: k
        do k = 1,size(m)
            s2k = real(var(k), dp) + real(nz, dp)
            if( w(k) > 0. )then
                lp(k) = log(real(w(k), dp)) - 0.5_dp * log(s2k) - 0.5_dp * (real(x, dp) - real(m(k), dp))**2 / s2k
            else
                lp(k) = -huge(1._dp)
            endif
            g(k) = (real(m(k), dp) - real(x, dp)) / s2k
        enddo
        p    = exp(lp - maxval(lp))
        p    = p / sum(p)
        post = real(p)
        label = maxloc(p, dim=1)
        gbar    = sum(p * g)
        sd_post = real(sqrt(max(real(v, dp), 0._dp)) * abs(p(label) * (g(label) - gbar)))
    end subroutine class_posterior

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

    ! PRIVATE

    ! Extreme-deconvolution fits of kk classes from the sorted-quantile and the largest-gap starts; the fit of
    ! higher likelihood, by decreasing mean. var is the intrinsic class variance (the noise of each value apart).
    subroutine fit_k( X1, R, Nz, kk, mu, var, w, post, lnl, tol, mstart )
        real(dp),           intent(in)  :: X1(:,:), R(:,:,:), Nz(:,:,:)
        integer,            intent(in)  :: kk
        real(dp),           intent(out) :: mu(kk), var(kk), w(kk), post(size(X1,1),kk), lnl
        real(dp), optional, intent(in)  :: tol, mstart(kk)
        real(dp) :: mu2(kk), var2(kk), w2(kk), post2(size(X1,1),kk), lnl2, means0(1,kk)
        real     :: xs(size(X1,1)), gaps(size(X1,1)-1)
        integer  :: idx(size(X1,1)), gidx(size(X1,1)-1), bnd(kk+1), n, i, g, order(kk)
        n = size(X1,1)
        if( present(mstart) )then
            ! one start, from the given means
            means0(1,:) = mstart
            call fit_from(mu, var, w, post, lnl, means0)
        else
            call fit_from(mu, var, w, post, lnl)
        endif
        if( kk > 1 .and. .not. present(mstart) )then
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
                if( present(tol) )then
                    call g%new(1, kk, MAXITS_EM, tol, PI_MIN_EM)
                else
                    call g%new(1, kk, MAXITS_EM, TOL_EM, PI_MIN_EM)
                endif
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
                var_o = S(1,1,:)
            end subroutine fit_from

    end subroutine fit_k

    ! fit_k with one slope on the covariate cv shared by the classes; mu is the class value at cv = 0. One class: the
    ! least-squares slope, its maximum likelihood. Several: alternation of the classes fitted to x - slope cv (warm
    ! started, at the loose TOL_SLOPE_EM) and the slope refitted as the within-class weighted regression of x on cv,
    ! until it moves by less than SLOPE_TOL_REL of the mean intensity, from slope 0 and from the least-squares slope;
    ! the likelier slope is kept (a species that follows cv, such as a core, pulls the least-squares start away) and
    ! the classes refitted at it from their means at TOL_SLOPE_FIN
    subroutine fit_k_cov( X1, cv, R, Nz, kk, mu, var, w, post, lnl, slope )
        real(dp), intent(in)  :: X1(:,:), cv(:), R(:,:,:), Nz(:,:,:)
        integer,  intent(in)  :: kk
        real(dp), intent(out) :: mu(kk), var(kk), w(kk), post(size(X1,1),kk), lnl, slope
        real(dp) :: mu2(kk), var2(kk), w2(kk), post2(size(X1,1),kk), mstart(kk), lnl2, slope2, s0, xm, cm, sxx, smean, nzm
        integer  :: n
        n     = size(X1,1)
        nzm   = median_dp(Nz(1,1,:))
        xm    = sum(X1(:,1)) / real(n, dp)
        cm    = sum(cv) / real(n, dp)
        sxx   = sum((cv - cm)**2)
        smean = sum(abs(X1(:,1))) / real(n, dp)
        s0    = 0._dp
        if( sxx > DTINY ) s0 = sum((cv - cm) * (X1(:,1) - xm)) / sxx
        if( kk == 1 )then
            slope = s0
            call fit_k(shifted(slope), R, Nz, kk, mu, var, w, post, lnl)
            return
        endif
        call alternate(0._dp, mu,  var,  w,  post,  lnl,  slope)
        call alternate(s0,    mu2, var2, w2, post2, lnl2, slope2)
        if( lnl2 > lnl )then
            mu    = mu2
            slope = slope2
        endif
        mstart = mu
        call fit_k(shifted(slope), R, Nz, kk, mu, var, w, post, lnl, tol=TOL_SLOPE_FIN, mstart=mstart)

        contains

            subroutine alternate( slope0, mu_o, var_o, w_o, post_o, lnl_o, slope_o )
                real(dp), intent(in)  :: slope0
                real(dp), intent(out) :: mu_o(kk), var_o(kk), w_o(kk), post_o(size(X1,1),kk), lnl_o, slope_o
                real(dp) :: mprev(kk), num, den, wk, xb, cb, snew
                integer  :: it, k
                logical  :: l_conv
                slope_o = slope0
                call fit_k(shifted(slope_o), R, Nz, kk, mu_o, var_o, w_o, post_o, lnl_o, tol=TOL_SLOPE_EM)
                do it = 1,NITER_SLOPE
                    num = 0._dp
                    den = 0._dp
                    do k = 1,kk
                        wk = sum(post_o(:,k))
                        if( wk < DTINY ) cycle
                        xb  = sum(post_o(:,k) * X1(:,1)) / wk
                        cb  = sum(post_o(:,k) * cv) / wk
                        num = num + sum(post_o(:,k) * (X1(:,1) - xb) * (cv - cb)) / (var_o(k) + nzm)
                        den = den + sum(post_o(:,k) * (cv - cb)**2) / (var_o(k) + nzm)
                    enddo
                    snew = slope_o
                    if( den > DTINY ) snew = num / den
                    l_conv  = abs(snew - slope_o) <= SLOPE_TOL_REL * smean
                    slope_o = snew
                    mprev   = mu_o
                    call fit_k(shifted(slope_o), R, Nz, kk, mu_o, var_o, w_o, post_o, lnl_o, tol=TOL_SLOPE_EM, mstart=mprev)
                    if( l_conv ) exit
                enddo
            end subroutine alternate

            function shifted( b ) result( xs )
                real(dp), intent(in) :: b
                real(dp) :: xs(size(X1,1),1)
                xs(:,1) = X1(:,1) - b * cv
            end function shifted

    end subroutine fit_k_cov

    ! median of a vector
    real(dp) function median_dp( v )
        real(dp), intent(in) :: v(:)
        real :: vs(size(v))
        vs = real(v)
        call hpsort(vs)
        if( mod(size(vs), 2) == 0 )then
            median_dp = 0.5_dp * real(vs(size(vs)/2) + vs(size(vs)/2 + 1), dp)
        else
            median_dp = real(vs(size(vs)/2 + 1), dp)
        endif
    end function median_dp

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
