!@descr: dose-series attenuation model: one B-factor per window on a line in dose, and the matched-filter accumulation
! Model, equations and justification: doc/implementation_notes/planned/dose_filter3D.md
module simple_dose_series
use simple_defs
use simple_error,        only: simple_exception
use simple_fileio,       only: fopen, fclose, fileiochk
use simple_stat,         only: kish_ess
use simple_string_utils, only: int2str
use simple_string,       only: string
use simple_dose_curve,   only: dose_curve, scaled_dose
implicit none
#include "simple_local_flags.inc"

public :: dose_series, shell_independence, shell_correlation, line_tau2, DOSE_WEIGHTS_MODEL, DOSE_WEIGHTS_CONTROL, DOSE_WEIGHTS_FLAT
private

!> the accumulation schemes: the measured model, the Grant & Grigorieff control, the flat mean (a = 1)
integer, parameter :: DOSE_WEIGHTS_MODEL = 1, DOSE_WEIGHTS_CONTROL = 2, DOSE_WEIGHTS_FLAT = 3
integer, parameter :: N_SCHEMES = 3

integer,  parameter :: MIN_BAND_SHELLS = 6       !< a line in k^2 needs a band to be a line over
integer,  parameter :: N_ALT_PASSES    = 3       !< alternating passes of the rank-one fit
integer,  parameter :: N_LM_ITER       = 50      !< Levenberg-Marquardt iterations of the Guinier fit
integer,  parameter :: N_LM_TRIES      = 12      !< damping increases tried before a step is given up
integer,  parameter :: N_BANDS         = 8       !< resolution bands of the readout
real(dp), parameter :: OWN_FSC_FLOOR   = 0.1_dp
real(dp), parameter :: SE_RATIO        = 2.0_dp  !< shells entering the log-domain start of the Guinier fit
real(dp), parameter :: C_DETECT        = 5.0_dp  !< errors a window's scale must clear to count as measured

type :: dose_series
    private
    integer  :: nwin = 0, nsh = 0, klo = 0, khi = 0
    real(dp) :: kv = 300._dp
    real(dp) :: rscale = 1._dp, mask_frac = 1._dp
    real(dp) :: shell_corr = 1._dp  !< variance inflation of a fit across shells, from their correlation
    logical  :: l_line   = .true.
    logical  :: l_paired = .false.
    logical  :: l_exists = .false.
    integer,  allocatable :: nshell(:)       !< (nsh) coefficients per shell in the statistics
    real(dp), allocatable :: kfreq(:)        !< (nsh) spatial frequency of each shell, 1/A
    real(dp), allocatable :: dose(:)         !< (nwin) scaled accumulated dose at the window centre
    real(dp), allocatable :: ddose(:)        !< (nwin) scaled dose span of the window
    real(dp), allocatable :: sigma2(:,:)     !< (nwin,nsh) one half's noise variance
    real(dp), allocatable :: prec(:,:)       !< (nwin,nsh) 1/sigma2
    real(dp), allocatable :: gram(:,:,:)     !< (2*nwin,2*nwin,nsh) Gram matrix of the halves
    real(dp), allocatable :: cvar(:,:,:)     !< (nwin,nwin,nsh) variance of the cross-half products
    real(dp), allocatable :: atten(:,:)      !< (nwin,nsh) model-free relative attenuation
    real(dp), allocatable :: atten_se(:,:)   !< (nwin,nsh) its standard error
    real(dp), allocatable :: atten_cov(:,:,:) !< (nwin,nwin,nsh) joint attenuation covariance
    real(dp), allocatable :: bw_cov(:,:)      !< (nwin,nwin) covariance of the measured B-factors
    real(dp), allocatable :: sigpow(:)       !< (nsh) signal power of window 1
    real(dp), allocatable :: cw(:), cw_se(:), bw(:), bw_se(:), bw_star(:), zw(:), grms(:)
    logical,  allocatable :: l_rejected(:)   !< (nwin) no fit of the window's attenuation
    real(dp), allocatable :: weight(:,:,:)   !< (nwin,nsh,N_SCHEMES) accumulation weights
    real(dp), allocatable :: ess(:), abar(:) !< (nsh) of the model's weights
    real(dp), allocatable :: cv(:,:)         !< (N_SCHEMES,nsh) mean whitened cross-half prediction error
    real(dp), allocatable :: dsum(:)         !< (nsh) sum of paired model-minus-control errors
    real(dp), allocatable :: paired_delta(:) !< (nsh) paired mean error difference; descriptive only
    real(dp), allocatable :: ne_shell(:), ne_se(:) !< (nsh) model-free critical exposure and its error
    real(dp), allocatable :: ne_model_est(:)       !< (nsh) the per-window model's, through the same estimator
    real(dp) :: b0 = 0._dp, beta = 0._dp, beta_se = 0._dp, tau = 0._dp
    real(dp) :: lnc0 = 0._dp, gamma = 0._dp, gamma_se = 0._dp, tau_c = 0._dp  !< ln c = lnc0 - gamma D
    logical  :: l_line_fitted = .false.
  contains
    procedure :: new
    procedure :: set_window
    procedure :: accumulate_gram
    procedure :: set_band
    procedure :: fit
    procedure :: paired_cv
    procedure :: accumulate
    procedure :: zerodose
    procedure :: write_model_table
    procedure :: write_attenuation_table
    procedure :: write_weights_table
    procedure :: write_cv_table
    procedure :: write_decay_table
    procedure :: write_step_exposure
    procedure :: write_resolution_weights
    procedure :: export_curve
    procedure :: log_readout
    procedure :: get_beta
    procedure :: get_band
    procedure :: get_weights
    procedure :: get_cv
    procedure :: get_attenuation
    procedure :: kill
    procedure, private :: shell_buckets
    procedure, private :: noise_spectra
    procedure, private :: derive_band
    procedure, private :: fit_attenuation
    procedure, private :: fit_bfactors
    procedure, private :: fit_dose_curve
    procedure, private :: fit_shell_exposure
    procedure, private :: shell_exposure
    procedure, private :: step_exposure
    procedure, private :: set_weights
    procedure, private :: cross_validate
    procedure, private :: own_fsc
    procedure, private :: atten_star
    procedure, private :: scheme_atten
    procedure, private :: band_difference
    procedure, private :: sharpest
    procedure, private :: rel_se
end type dose_series

contains

    !> allocate for nwin windows over nsh shells; kfreq is the shell frequency grid in 1/A
    subroutine new( self, nwin, kfreq, kv, l_line )
        class(dose_series), intent(inout) :: self
        integer,            intent(in)    :: nwin
        real(dp),           intent(in)    :: kfreq(:), kv
        logical,            intent(in)    :: l_line
        call self%kill
        if( nwin < 2 )        THROW_HARD('A dose series needs at least two windows')
        if( size(kfreq) < 1 ) THROW_HARD('Empty shell grid')
        self%nwin   = nwin
        self%nsh    = size(kfreq)
        self%kv     = kv
        self%l_line = l_line
        allocate(self%kfreq(self%nsh), source=kfreq)
        allocate(self%nshell(self%nsh), source=0)
        allocate(self%dose(nwin), self%ddose(nwin), source=0._dp)
        allocate(self%sigma2(nwin,self%nsh), self%prec(nwin,self%nsh), source=0._dp)
        allocate(self%gram(2*nwin,2*nwin,self%nsh), source=0._dp)
        allocate(self%cvar(nwin,nwin,self%nsh), source=0._dp)
        allocate(self%atten(nwin,self%nsh), self%atten_se(nwin,self%nsh), source=0._dp)
        allocate(self%atten_cov(nwin,nwin,self%nsh), self%bw_cov(nwin,nwin), source=0._dp)
        allocate(self%sigpow(self%nsh), self%ess(self%nsh), self%abar(self%nsh), source=0._dp)
        allocate(self%cw(nwin), self%cw_se(nwin), self%bw(nwin), self%bw_se(nwin), source=0._dp)
        allocate(self%bw_star(nwin), self%zw(nwin), self%grms(nwin), source=0._dp)
        allocate(self%l_rejected(nwin), source=.false.)
        allocate(self%weight(nwin,self%nsh,N_SCHEMES), source=0._dp)
        allocate(self%cv(N_SCHEMES,self%nsh), source=0._dp)
        allocate(self%dsum(self%nsh), self%paired_delta(self%nsh), source=0._dp)
        allocate(self%ne_shell(self%nsh), self%ne_se(self%nsh), self%ne_model_est(self%nsh), source=0._dp)
        self%l_exists = .true.
    end subroutine new

    !> Fraction of a shell's Fourier coefficients that are independent under a real-space mask m, the
    !! mask_frac of accumulate_gram. Masked noise is correlated in Fourier space by a kernel of m^2, and a
    !! shell cuts only a plane of it: by the projection-slice theorem the plane's sum of squared
    !! correlations is N^2 sum_xy (sum_z m^2)^2 / (sum m^2)^2, averaged here over the three axes.
    pure real(dp) function shell_independence( m ) result( frac )
        real, intent(in) :: m(:,:,:)
        real(dp), allocatable :: m2(:,:,:)
        real(dp) :: total, g
        integer  :: n(3)
        n     = shape(m)
        m2    = real(m, dp)**2
        total = sum(m2)
        frac  = 0._dp
        if( total <= 0._dp ) return
        g = real(n(1),dp)*real(n(2),dp) * sum(sum(m2, dim=3)**2) + &
           &real(n(1),dp)*real(n(3),dp) * sum(sum(m2, dim=2)**2) + &
           &real(n(2),dp)*real(n(3),dp) * sum(sum(m2, dim=1)**2)
        frac = 3._dp * total**2 / g
    end function shell_independence

    !> Variance inflation of a smooth fit across shells under a real-space mask m, the shell_corr of
    !! accumulate_gram: the |FT(m^2)|^2 kernel that correlates coefficients reaches into neighbouring shells.
    !! A shell holds its plane through the kernel, n_a sum_xy (sum_a m^2)^2 against all of it, N^3 sum m^4;
    !! the ratio, averaged over the three axes, is the sum of the correlations between shells. 1 without a mask.
    pure real(dp) function shell_correlation( m ) result( f )
        real, intent(in) :: m(:,:,:)
        real(dp), allocatable :: m2(:,:,:)
        real(dp) :: quart, p(3)
        integer  :: n(3)
        n     = shape(m)
        m2    = real(m, dp)**2
        quart = sum(m2 * m2)
        f     = 1._dp
        p(1)  = sum(sum(m2, dim=1)**2)
        p(2)  = sum(sum(m2, dim=2)**2)
        p(3)  = sum(sum(m2, dim=3)**2)
        if( quart <= 0._dp .or. any(p <= 0._dp) ) return
        f = max(1._dp, (real(n(1),dp)/p(1) + real(n(2),dp)/p(2) + real(n(3),dp)/p(3)) * quart / 3._dp)
    end function shell_correlation

    !> Method-of-moments (DerSimonian-Laird) excess variance tau^2 of points b(d) +- se about their weighted
    !! line, untruncated. E[Q] = (k-2) + tau^2 tr(P), P = W - WX(X'WX)^-1 X'W, X = [1, d], W = diag(1/se^2):
    !! the fitted slope takes a second degree of freedom out of the residuals and out of tr(P) alike.
    pure real(dp) function line_tau2( b, d, se ) result( tau2 )
        real(dp), intent(in) :: b(:), d(:), se(:)
        real(dp) :: w(size(b)), s0, s1, s2, t0, t1, t2, sy, sxy, det, b0, slope, q, trp
        integer  :: k
        tau2 = 0._dp
        k    = size(b)
        if( k < 3 ) return
        w   = 1._dp / se**2
        s0  = sum(w);      s1  = sum(w*d);      s2  = sum(w*d*d)
        t0  = sum(w*w);    t1  = sum(w*w*d);    t2  = sum(w*w*d*d)
        sy  = sum(w*b);    sxy = sum(w*d*b)
        det = s0*s2 - s1*s1
        if( det <= DTINY ) return
        b0    = ( s2*sy - s1*sxy) / det
        slope = (-s1*sy + s0*sxy) / det
        q     = sum(w * (b - b0 - slope*d)**2)
        trp   = s0 - (s2*t0 - 2._dp*s1*t1 + s0*t2) / det
        if( trp <= DTINY ) return
        tau2  = (q - real(k-2, dp)) / trp
    end function line_tau2

    !> The weighted line y = a + b d through points with errors se, in two passes: weights 1/se^2, then
    !! 1/(se^2 + tau^2) with the excess scatter tau of line_tau2. ok is false when d does not vary.
    pure subroutine scatter_line( y, d, se, a, b, b_se, tau, ok )
        real(dp), intent(in)  :: y(:), d(:), se(:)
        real(dp), intent(out) :: a, b, b_se, tau
        logical,  intent(out) :: ok
        real(dp) :: wgt(size(y)), s0, s1, s2, sy, sxy, det
        integer  :: it
        a = 0._dp; b = 0._dp; b_se = 0._dp; tau = 0._dp
        ok = .false.
        do it = 1, 2
            wgt = 1._dp / (se**2 + tau**2)
            s0  = sum(wgt);   s1  = sum(wgt*d);   s2 = sum(wgt*d*d)
            sy  = sum(wgt*y); sxy = sum(wgt*d*y)
            det = s0*s2 - s1*s1
            if( abs(det) <= DTINY ) exit
            a    = ( s2*sy - s1*sxy) / det
            b    = (-s1*sy + s0*sxy) / det
            b_se = sqrt(max(DTINY, s0 / det))
            ok   = .true.
            if( it == 1 ) tau = sqrt(max(0._dp, line_tau2(y, d, se)))
        enddo
    end subroutine scatter_line

    !> physical dose bounds of one window, kept in the scaled dose of calc_dose_weights
    subroutine set_window( self, iwin, dose_lo, dose_hi )
        class(dose_series), intent(inout) :: self
        integer,            intent(in)    :: iwin
        real(dp),           intent(in)    :: dose_lo, dose_hi
        if( .not. self%l_exists )            THROW_HARD('Dose series not constructed')
        if( iwin < 1 .or. iwin > self%nwin ) THROW_HARD('Window index out of range')
        if( dose_hi <= dose_lo )             THROW_HARD('Window dose span must be positive')
        self%dose(iwin)  = scaled_dose(0.5_dp * (dose_lo + dose_hi), self%kv)
        self%ddose(iwin) = scaled_dose(dose_hi - dose_lo, self%kv)
    end subroutine set_window

    !> order the in-range coefficients by shell, so a loop over shells sums each one in a fixed order
    subroutine shell_buckets( self, shell, offset, order )
        class(dose_series),   intent(in)    :: self
        integer,              intent(in)    :: shell(:)
        integer, allocatable, intent(inout) :: offset(:), order(:)
        integer, allocatable :: cnt(:)
        integer :: i, s
        if( allocated(offset) ) deallocate(offset)
        if( allocated(order)  ) deallocate(order)
        allocate(cnt(self%nsh), source=0)
        do i = 1, size(shell)
            s = shell(i)
            if( s >= 1 .and. s <= self%nsh ) cnt(s) = cnt(s) + 1
        enddo
        allocate(offset(self%nsh+1))
        offset(1) = 1
        do s = 1, self%nsh
            offset(s+1) = offset(s) + cnt(s)
        enddo
        allocate(order(max(1, offset(self%nsh+1) - 1)), source=0)
        cnt = 0
        do i = 1, size(shell)
            s = shell(i)
            if( s < 1 .or. s > self%nsh ) cycle
            order(offset(s) + cnt(s)) = i
            cnt(s) = cnt(s) + 1
        enddo
    end subroutine shell_buckets

    !> one pass over the masked coefficients: noise spectra, Gram matrices and cross-product variances
    subroutine accumulate_gram( self, em, om, shell, mask_frac, shell_corr )
        class(dose_series), intent(inout) :: self
        complex(sp),        intent(in)    :: em(:,:), om(:,:)   !< (ncoeff,nwin) masked halves
        integer,            intent(in)    :: shell(:)           !< (ncoeff) shell of each coefficient
        real(dp),           intent(in)    :: mask_frac          !< fraction of a shell's coefficients the mask leaves independent
        real(dp),           intent(in)    :: shell_corr         !< variance inflation of a fit across shells (shell_correlation)
        integer,  allocatable :: offset(:), order(:)
        complex(dp) :: y(2*self%nwin)
        real(dp)    :: x, invn, invneff, pwr, nvals, vx(self%nwin,self%nwin)
        integer     :: ncoeff, i, s, w, v, a, nn, lo, hi
        if( .not. self%l_exists ) THROW_HARD('Dose series not constructed')
        ncoeff = size(shell)
        if( size(em,1) /= ncoeff .or. size(om,1) /= ncoeff ) THROW_HARD('Coefficient count mismatch')
        if( size(em,2) /= self%nwin .or. size(om,2) /= self%nwin ) THROW_HARD('Window count mismatch')
        if( mask_frac <= 0._dp .or. mask_frac > 1._dp ) THROW_HARD('Independent shell fraction out of range')
        if( shell_corr < 1._dp ) THROW_HARD('Shell correlation below one')
        call self%shell_buckets(shell, offset, order)
        ! every statistic downstream is a ratio, so one global scale is free; fixing it at unit mean
        ! power keeps the guards meaningful whatever the FFT normalisation of the caller
        pwr = 0._dp
        do a = 1, offset(self%nsh+1) - 1
            i = order(a)
            do w = 1, self%nwin
                pwr = pwr + real(abs(em(i,w)),dp)**2 + real(abs(om(i,w)),dp)**2
            enddo
        enddo
        ! the count of values summed passes 2^31 at box ~590 with 20 windows: form it in double precision
        nvals = 2._dp * real(self%nwin, dp) * real(offset(self%nsh+1) - 1, dp)
        if( nvals < 1._dp .or. pwr <= 0._dp ) THROW_HARD('No signal in the dose series coefficients')
        self%rscale    = 1._dp / sqrt(pwr / nvals)
        self%mask_frac = mask_frac
        self%shell_corr = shell_corr
        do s = 1, self%nsh
            self%nshell(s) = offset(s+1) - offset(s)
        enddo
        self%gram   = 0._dp
        self%sigma2 = 0._dp
        !$omp parallel do default(shared) private(s,lo,hi,nn,i,w,v,a,y,x,invn,invneff,vx) schedule(dynamic)
        do s = 1, self%nsh
            lo = offset(s); hi = offset(s+1) - 1
            nn = hi - lo + 1
            if( nn < 1 ) cycle
            do a = lo, hi
                i = order(a)
                do w = 1, self%nwin
                    y(w)           = self%rscale * cmplx(em(i,w), kind=dp)
                    y(self%nwin+w) = self%rscale * cmplx(om(i,w), kind=dp)
                enddo
                do w = 1, 2*self%nwin
                    do v = w, 2*self%nwin
                        self%gram(w,v,s) = self%gram(w,v,s) + real(y(w)*conjg(y(v)), dp)
                    enddo
                enddo
                do w = 1, self%nwin
                    x = real(abs(y(w) - y(self%nwin+w)), dp)
                    self%sigma2(w,s) = self%sigma2(w,s) + x*x
                enddo
            enddo
            invn = 1._dp / real(nn, dp)
            do w = 1, 2*self%nwin
                do v = w, 2*self%nwin
                    self%gram(w,v,s) = self%gram(w,v,s) * invn
                    self%gram(v,w,s) = self%gram(w,v,s)
                enddo
            enddo
            do w = 1, self%nwin
                self%sigma2(w,s) = 0.5_dp * self%sigma2(w,s) * invn
            enddo
            ! The error of X_wv as an estimate of a_w a_v P is its noise alone, given the structure:
            ! (|a_w S|^2 sigma_v^2 + |a_v S|^2 sigma_w^2 + sigma_w^2 sigma_v^2) / 2 per coefficient. The
            ! spread of the products themselves would add the structure's own power fluctuation,
            ! which is common to every pair and cancels in the attenuations.
            ! The mask correlates neighbouring coefficients: the shell holds fewer independent samples.
            invneff = 1._dp / max(1._dp, mask_frac * real(nn, dp))
            do w = 1, self%nwin
                do v = 1, self%nwin
                    vx(w,v) = 0.5_dp * invneff * (self%gram(w,w,s) * self%sigma2(v,s) + &
                        &self%gram(self%nwin+v,self%nwin+v,s) * self%sigma2(w,s) - self%sigma2(w,s) * self%sigma2(v,s))
                enddo
            enddo
            ! C = (X + X^T)/2 averages two estimates with independent noise
            do w = 1, self%nwin
                self%cvar(w,w,s) = max(DTINY, vx(w,w))
                do v = 1, self%nwin
                    if( v /= w ) self%cvar(w,v,s) = max(DTINY, 0.25_dp * (vx(w,v) + vx(v,w)))
                enddo
            enddo
        enddo
        !$omp end parallel do
        call self%noise_spectra
    end subroutine accumulate_gram

    !> precision from the frozen noise spectra; a dead shell gets zero weight rather than a huge one
    subroutine noise_spectra( self )
        class(dose_series), intent(inout) :: self
        integer :: w, s
        do s = 1, self%nsh
            do w = 1, self%nwin
                if( self%sigma2(w,s) > DTINY )then
                    self%prec(w,s) = 1._dp / self%sigma2(w,s)
                else
                    self%prec(w,s) = 0._dp
                endif
            enddo
        enddo
    end subroutine noise_spectra

    !> the band the model is fitted over: from the high-pass shell, which keeps the mask and solvent
    !! artefacts of the lowest shells out of the Guinier fit, to lp_shell; lp_shell <= 0 derives it
    !! from the windows' own FSCs
    subroutine set_band( self, hp_shell, lp_shell )
        class(dose_series), intent(inout) :: self
        integer,            intent(in)    :: hp_shell, lp_shell
        if( .not. self%l_exists ) THROW_HARD('Dose series not constructed')
        self%klo = max(1, hp_shell)
        if( lp_shell > 0 )then
            self%khi = min(lp_shell, self%nsh)
        else
            self%khi = self%derive_band()
        endif
        if( self%khi - self%klo + 1 < MIN_BAND_SHELLS ) THROW_HARD('Fitting band is shorter than six shells')
    end subroutine set_band

    !> last shell where any window's own FSC still holds; the early windows reach past the flat sum
    integer function derive_band( self ) result( khi )
        class(dose_series), intent(in) :: self
        integer :: s, w
        khi = self%klo
        do s = 1, self%nsh
            do w = 1, self%nwin
                if( self%own_fsc(w,s) > OWN_FSC_FLOOR )then
                    khi = s
                    exit
                endif
            enddo
        enddo
    end function derive_band

    !> FSC of one window against itself: its cross-half product over its two half powers
    real(dp) function own_fsc( self, w, s ) result( fsc )
        class(dose_series), intent(in) :: self
        integer,            intent(in) :: w, s
        real(dp) :: denom
        denom = self%gram(w,w,s) * self%gram(self%nwin+w,self%nwin+w,s)
        if( denom <= DTINY )then
            fsc = 0._dp
        else
            fsc = self%gram(w,self%nwin+w,s) / sqrt(denom)
        endif
    end function own_fsc

    !> steps 3 to 6 of the note, in order; stops when nothing could be measured rather than accumulating blind
    subroutine fit( self )
        class(dose_series), intent(inout) :: self
        real(dp) :: noise, power
        integer  :: s, w
        if( .not. self%l_exists ) THROW_HARD('Dose series not constructed')
        if( self%khi < self%klo ) THROW_HARD('Band not set')
        ! the errors are built on the half difference: identical halves leave every one of them undefined
        noise = 0._dp; power = 0._dp
        do s = self%klo, self%khi
            do w = 1, self%nwin
                noise = noise + self%sigma2(w,s)
                power = power + 0.5_dp * (self%gram(w,w,s) + self%gram(self%nwin+w,self%nwin+w,s))
            enddo
        enddo
        if( noise <= 1.e-6_dp * power ) THROW_HARD('The even and odd halves carry no measurable noise: are they the same map?')
        self%l_paired = .false.
        call self%fit_attenuation
        call self%fit_bfactors
        if( all(self%l_rejected(2:)) ) THROW_HARD('No window after the first could be fitted: no dose curve to measure')
        call self%fit_dose_curve
        call self%fit_shell_exposure
        call self%set_weights
        call self%cross_validate
    end subroutine fit

    !> rank-one weighted least squares of the cross-half block, per shell, gauged at a_1 = 1: objective
    !! (9) over unordered pairs, with the errors of attenuation_se
    subroutine fit_attenuation( self )
        class(dose_series), intent(inout) :: self
        real(dp) :: a(self%nwin), c(self%nwin,self%nwin), se(self%nwin), p, num, den, aold, neff
        integer  :: s, w, v, it
        !$omp parallel do default(shared) private(s,w,v,it,a,c,se,p,num,den,aold,neff) schedule(dynamic)
        do s = self%klo, self%khi
            do w = 1, self%nwin
                do v = 1, self%nwin
                    c(w,v) = 0.5_dp * (self%gram(w,self%nwin+v,s) + self%gram(v,self%nwin+w,s))
                enddo
            enddo
            ! start from the leading eigenvector by power iteration, then alternate
            a = 1._dp
            do it = 1, 8
                a = matmul(c, a)
                den = sqrt(max(DTINY, sum(a*a)))
                a = a / den
            enddo
            if( abs(a(1)) < 1.e-8_dp )then
                a = 1._dp
            else
                a = a / a(1)
            endif
            p = max(DTINY, self%gram(1,self%nwin+1,s))
            do it = 1, N_ALT_PASSES
                ! each a_w given the others; its own pair enters (9) once but holds a_w twice
                do w = 2, self%nwin
                    aold = a(w)
                    num  = 2._dp * c(w,w) * aold / self%cvar(w,w,s)
                    den  = 2._dp * aold * aold / self%cvar(w,w,s)
                    do v = 1, self%nwin
                        if( v == w ) cycle
                        num = num + c(w,v) * a(v) / self%cvar(w,v,s)
                        den = den + a(v) * a(v) / self%cvar(w,v,s)
                    enddo
                    den = p * den
                    if( abs(den) > DTINY ) a(w) = num / den
                enddo
                num = 0._dp; den = 0._dp
                do w = 1, self%nwin
                    do v = w, self%nwin
                        num = num + c(w,v) * a(w) * a(v) / self%cvar(w,v,s)
                        den = den + (a(w)*a(v))**2 / self%cvar(w,v,s)
                    enddo
                enddo
                if( den > DTINY ) p = num / den
                if( abs(a(1)) > 1.e-8_dp )then
                    p = p * a(1)*a(1)
                    a = a / a(1)
                endif
            enddo
            self%sigpow(s) = p
            neff = max(1._dp, self%mask_frac * real(self%nshell(s), dp))
            call attenuation_se(a, p, self%cvar(:,:,s), self%sigma2(:,s), neff, se, self%atten_cov(:,:,s))
            self%atten(:,s)    = a
            self%atten_se(:,s) = se
            self%atten(1,s)    = 1._dp
            self%atten_se(1,s) = 0._dp
        enddo
        !$omp end parallel do
    end subroutine fit_attenuation

    !> Standard errors of a_2..a_W from the sandwich covariance of (9). To first order the error of C_wv is
    !! a_w g_v + a_v g_w + r_wv: g_k, the projection of window k's half noise onto the structure, is shared
    !! by every pair holding k and does not average down over them, which weights alone would assume.
    pure subroutine attenuation_se( a, p, cvar, sigma2, neff, se, acov )
        use simple_linalg, only: spd_inverse
        real(dp), intent(in)  :: a(:), p, cvar(:,:), sigma2(:), neff
        real(dp), intent(out) :: se(:), acov(:,:)
        real(dp), allocatable :: jac(:,:), lmat(:,:), dw(:), rv(:), lam(:), h(:,:), hinv(:,:), b(:,:), jdl(:,:), cov(:,:)
        integer :: nw, m, ip, w, v, k
        logical :: ok
        nw = size(a)
        m  = nw*(nw+1)/2
        ! parameters: a_2..a_nw in columns 1..nw-1, P in column nw
        allocate(jac(m,nw), lmat(m,nw), dw(m), rv(m), lam(nw), h(nw,nw), hinv(nw,nw), b(nw,nw), jdl(nw,nw), cov(nw,nw))
        jac  = 0._dp
        lmat = 0._dp
        ip   = 0
        do w = 1, nw
            do v = w, nw
                ip     = ip + 1
                dw(ip) = 1._dp / cvar(w,v)
                if( w == v )then
                    if( w > 1 ) jac(ip,w-1) = 2._dp * p * a(w)
                    lmat(ip,w) = 2._dp * a(w)
                    rv(ip)     = sigma2(w)**2 / (2._dp * neff)
                else
                    if( w > 1 ) jac(ip,w-1) = p * a(v)
                    jac(ip,v-1) = p * a(w)
                    lmat(ip,v)  = a(w)
                    lmat(ip,w)  = a(v)
                    rv(ip)      = sigma2(w) * sigma2(v) / (4._dp * neff)
                endif
                jac(ip,nw) = a(w) * a(v)
            enddo
        enddo
        do k = 1, nw
            lam(k) = max(0._dp, p) * sigma2(k) / (4._dp * neff)
        enddo
        h   = matmul(transpose(jac), spread(dw, 2, nw) * jac)
        jdl = matmul(transpose(jac), spread(dw, 2, nw) * lmat)
        b   = matmul(jdl, spread(lam, 2, nw) * transpose(jdl)) + &
            &matmul(transpose(jac), spread(dw*dw*rv, 2, nw) * jac)
        call spd_inverse(h, hinv, nw, ok)
        se = huge(1._dp)
        se(1) = 0._dp
        acov = 0._dp
        if( .not. ok ) return
        cov = matmul(hinv, matmul(b, hinv))
        acov(2:,2:) = 0.5_dp * (cov(:nw-1,:nw-1) + transpose(cov(:nw-1,:nw-1)))
        do k = 2, nw
            if( cov(k-1,k-1) > 0._dp ) se(k) = sqrt(cov(k-1,k-1))
        enddo
    end subroutine attenuation_se

    !> One scale and one B-factor per window, refined in the linear domain so faded shells count right.
    !! Levenberg-Marquardt keeps a step only when it lowers (10). A fit that fails, or whose scale does not clear
    !! C_DETECT errors, is rejected: its window is left out of the line and of the model's weights.
    subroutine fit_bfactors( self )
        use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
        class(dose_series), intent(inout) :: self
        character(len=:), allocatable :: msg
        real(dp) :: k2(self%nsh), wgt, lsum(4), sxy, rhs(2), jtj(3), det, c, b, cn, bn, obj, objn, lam, varc
        real(dp) :: db_da(self%nwin,self%nsh)
        integer  :: w, v, s, it, ntry, nused
        logical  :: l_step
        k2 = self%kfreq * self%kfreq
        self%l_rejected = .false.
        db_da = 0._dp
        self%bw_cov = 0._dp
        self%cw(1)    = 1._dp
        self%cw_se(1) = 0._dp
        self%bw(1)    = 0._dp
        self%bw_se(1) = 0._dp
        self%grms(1)  = 0._dp
        do w = 2, self%nwin
            ! weighted regression of log a on k^2 over the shells the window still resolves
            lsum = 0._dp; sxy = 0._dp; nused = 0
            do s = self%klo, self%khi
                if( self%atten_se(w,s) >= huge(1._dp) ) cycle
                if( self%atten(w,s) <= 0._dp ) cycle
                if( self%atten(w,s) <= SE_RATIO * self%atten_se(w,s) ) cycle
                wgt     = (self%atten(w,s) / max(DTINY, self%atten_se(w,s)))**2
                lsum(1) = lsum(1) + wgt
                lsum(2) = lsum(2) + wgt * k2(s)
                lsum(3) = lsum(3) + wgt * k2(s) * k2(s)
                lsum(4) = lsum(4) + wgt * log(self%atten(w,s))
                sxy     = sxy     + wgt * k2(s) * log(self%atten(w,s))
                nused   = nused + 1
            enddo
            det = 0._dp
            if( nused >= 2 ) det = lsum(1)*lsum(3) - lsum(2)*lsum(2)
            if( abs(det) > DTINY )then
                ! the slope of log a against k^2 is -B/4
                b = -4._dp * (lsum(1)*sxy - lsum(2)*lsum(4)) / det
                c = exp((lsum(3)*lsum(4) - lsum(2)*sxy) / det)
            else
                b = 0._dp
                c = 1._dp
            endif
            ! Levenberg-Marquardt in the linear domain
            obj = guinier_obj(c, b)
            lam = 1.e-3_dp
            do it = 1, N_LM_ITER
                call guinier_normal(c, b, jtj, rhs)
                l_step = .false.
                do ntry = 1, N_LM_TRIES
                    det = jtj(1)*(1._dp+lam) * jtj(3)*(1._dp+lam) - jtj(2)*jtj(2)
                    if( det > DTINY )then
                        cn   = c + ( jtj(3)*(1._dp+lam)*rhs(1) - jtj(2)*rhs(2)) / det
                        bn   = b + (-jtj(2)*rhs(1) + jtj(1)*(1._dp+lam)*rhs(2)) / det
                        objn = guinier_obj(cn, bn)
                        if( ieee_is_finite(objn) )then
                            if( objn < obj )then
                                l_step = .true.
                                exit
                            endif
                        endif
                    endif
                    lam = 10._dp * lam
                enddo
                if( .not. l_step ) exit
                c   = cn
                b   = bn
                lam = max(1.e-12_dp, 0.1_dp * lam)
                if( obj - objn <= 1.e-10_dp * obj )then
                    obj = objn
                    exit
                endif
                obj = objn
            enddo
            self%cw(w) = c
            self%bw(w) = b
            call guinier_normal(c, b, jtj, rhs)
            det = jtj(1)*jtj(3) - jtj(2)*jtj(2)
            varc = huge(1._dp)
            self%bw_se(w) = huge(1._dp)
            if( det > DTINY .and. jtj(1) > DTINY )then
                ! the shells are not independent under the mask: their correlation inflates both variances
                self%bw_se(w) = sqrt(max(DTINY, self%shell_corr * jtj(1) / det))
                varc          = max(0._dp, self%shell_corr * jtj(3) / det)
            endif
            self%cw_se(w) = sqrt(max(DTINY, varc))
            self%grms(w) = huge(1._dp)
            if( ieee_is_finite(obj) .and. nused_shells() > 0 ) self%grms(w) = sqrt(obj / real(nused_shells(), dp))
            ! no curvature, a non-finite fit, or a scale lost in its own error: no measurable signal to weight
            if( .not. ieee_is_finite(obj) .or. .not. ieee_is_finite(c) .or. .not. ieee_is_finite(b) .or. &
                &self%bw_se(w) >= huge(1._dp) )then
                self%l_rejected(w) = .true.
            else if( c <= C_DETECT * sqrt(varc) )then
                self%l_rejected(w) = .true.
            endif
            if( self%l_rejected(w) )then
                self%bw_se(w) = huge(1._dp)
                msg = 'Dose series window '//int2str(w)//' has no fit of its attenuation: left out of the model weights'
                THROW_WARN(msg)
            else
                do s = self%klo, self%khi
                    if( self%atten_se(w,s) >= huge(1._dp) ) cycle
                    wgt = 1._dp / max(DTINY, self%atten_se(w,s)**2)
                    db_da(w,s) = wgt * bfac_envelope(b, k2(s)) * (-jtj(2) - jtj(1)*c*k2(s)/4._dp) / det
                enddo
            endif
        enddo
        ! The fits share their noisy reference; retain that covariance when changing the reporting gauge.
        do w = 2, self%nwin
            do v = 2, self%nwin
                self%bw_cov(w,v) = self%shell_corr * sum(db_da(w,:) * db_da(v,:) * self%atten_cov(w,v,:))
            enddo
        enddo

      contains

        !> the objective (10) of window w at (cc, bb)
        real(dp) function guinier_obj( cc, bb ) result( f )
            real(dp), intent(in) :: cc, bb
            integer :: ss
            f = 0._dp
            do ss = self%klo, self%khi
                if( self%atten_se(w,ss) >= huge(1._dp) ) cycle
                f = f + ((self%atten(w,ss) - cc * bfac_envelope(bb, k2(ss))) / max(DTINY, self%atten_se(w,ss)))**2
            enddo
        end function guinier_obj

        !> the Gauss-Newton normal equations of (10) at (cc, bb): jtj holds (cc,cc), (cc,bb), (bb,bb)
        subroutine guinier_normal( cc, bb, jtj_out, rhs_out )
            real(dp), intent(in)  :: cc, bb
            real(dp), intent(out) :: jtj_out(3), rhs_out(2)
            real(dp) :: ww, ee, rr, dcc, dbb
            integer  :: ss
            jtj_out = 0._dp; rhs_out = 0._dp
            do ss = self%klo, self%khi
                if( self%atten_se(w,ss) >= huge(1._dp) ) cycle
                ww  = 1._dp / max(DTINY, self%atten_se(w,ss)**2)
                ee  = bfac_envelope(bb, k2(ss))
                rr  = self%atten(w,ss) - cc * ee
                dcc = ee
                dbb = -cc * ee * k2(ss) / 4._dp
                jtj_out(1) = jtj_out(1) + ww * dcc * dcc
                jtj_out(2) = jtj_out(2) + ww * dcc * dbb
                jtj_out(3) = jtj_out(3) + ww * dbb * dbb
                rhs_out(1) = rhs_out(1) + ww * dcc * rr
                rhs_out(2) = rhs_out(2) + ww * dbb * rr
            enddo
        end subroutine guinier_normal

        integer function nused_shells() result( m )
            integer :: ss
            m = 0
            do ss = self%klo, self%khi
                if( self%atten_se(w,ss) < huge(1._dp) ) m = m + 1
            enddo
        end function nused_shells

    end subroutine fit_bfactors

    !> the weighted lines of B and ln c in dose through the measured windows, their excess scatter, and the
    !! shrunk B-factors
    subroutine fit_dose_curve( self )
        class(dose_series), intent(inout) :: self
        real(dp) :: wgt, det, pred, vb, slope
        logical  :: free(self%nwin-1), l_fitted
        integer  :: w
        self%b0   = 0._dp; self%beta  = 0._dp; self%beta_se  = 0._dp; self%tau   = 0._dp
        self%lnc0 = 0._dp; self%gamma = 0._dp; self%gamma_se = 0._dp; self%tau_c = 0._dp
        self%l_line_fitted = .false.
        free = self%bw_se(2:) < huge(1._dp)
        if( count(free) >= 2 )then
            call scatter_line(pack(self%bw(2:), free), pack(self%dose(2:), free), pack(self%bw_se(2:), free), &
                &self%b0, self%beta, self%beta_se, self%tau, self%l_line_fitted)
            ! a measured window has c > 2 se_c > 0
            call scatter_line(log(pack(self%cw(2:), free)), pack(self%dose(2:), free), &
                &pack(self%cw_se(2:), free) / pack(self%cw(2:), free), self%lnc0, slope, self%gamma_se, self%tau_c, l_fitted)
            self%gamma = -slope
        else
            ! only one measured window: the lines through the gauge point and it
            do w = 2, self%nwin
                if( .not. free(w-1) ) cycle
                det = self%dose(w) - self%dose(1)
                if( abs(det) > DTINY )then
                    self%beta     = self%bw(w) / det
                    self%beta_se  = self%bw_se(w) / abs(det)
                    self%b0       = -self%beta * self%dose(1)
                    self%gamma    = -log(self%cw(w)) / det
                    self%gamma_se = self%cw_se(w) / (self%cw(w) * abs(det))
                    self%lnc0     = self%gamma * self%dose(1)
                    self%l_line_fitted = .true.
                endif
                exit
            enddo
        endif
        ! window 1 is the gauge: its diagnostic is the predicted B(D_1), not a residual
        self%zw(1) = 0._dp
        do w = 2, self%nwin
            pred = self%b0 + self%beta * self%dose(w)
            self%zw(w) = 0._dp
            if( self%bw_se(w) >= huge(1._dp) ) cycle
            vb = sqrt(self%bw_se(w)**2 + self%tau**2)
            if( vb > DTINY ) self%zw(w) = (self%bw(w) - pred) / vb
        enddo
        ! window 1 is the gauge, not a measurement: it keeps B = 0 whatever the line says
        self%bw_star(1) = 0._dp
        do w = 2, self%nwin
            if( self%l_rejected(w) .or. .not. self%l_line .or. .not. self%l_line_fitted )then
                self%bw_star(w) = self%bw(w)
            else if( self%bw_se(w) >= huge(1._dp) )then
                self%bw_star(w) = self%b0 + self%beta * self%dose(w)
            else if( self%tau <= DTINY )then
                self%bw_star(w) = self%b0 + self%beta * self%dose(w)
            else
                wgt = 1._dp/self%bw_se(w)**2 + 1._dp/self%tau**2
                self%bw_star(w) = (self%bw(w)/self%bw_se(w)**2 + &
                    &(self%b0 + self%beta*self%dose(w))/self%tau**2) / wgt
            endif
        enddo
    end subroutine fit_dose_curve

    !> Critical exposure per shell, measured and from the per-window model, through one estimator
    subroutine fit_shell_exposure( self )
        class(dose_series), intent(inout) :: self
        real(dp) :: model_se
        integer  :: s
        do s = 1, self%nsh
            call self%shell_exposure(s, .false., self%ne_shell(s), self%ne_se(s))
            call self%shell_exposure(s, .true.,  self%ne_model_est(s), model_se)
        enddo
    end subroutine fit_shell_exposure

    !> Critical exposure at shell s with no model of its dose dependence: ln a_w(s) regressed through the origin on
    !! D_w - D_1 over the windows the shell still resolves, weighted by (a_w/se_w)^2, Ne = -1/(2 slope) as
    !! exp(-D/2Ne) has it. l_model regresses the per-window model's attenuation instead, with the same windows and
    !! weights. NaN outside the band and where the shell shows no decay beyond twice its error.
    subroutine shell_exposure( self, s, l_model, ne, ne_se )
        use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan
        class(dose_series), intent(in)  :: self
        integer,            intent(in)  :: s
        logical,            intent(in)  :: l_model
        real(dp),           intent(out) :: ne, ne_se
        real(dp) :: sxx, sxy, x, y, wgt, m, sem, gradient(self%nwin)
        integer  :: w, nused
        ne    = ieee_value(1._dp, ieee_quiet_nan)
        ne_se = ne
        if( s < self%klo .or. s > self%khi ) return
        sxx = 0._dp; sxy = 0._dp; nused = 0
        gradient = 0._dp
        do w = 2, self%nwin
            if( self%atten_se(w,s) >= huge(1._dp) ) cycle
            if( self%atten(w,s) <= SE_RATIO * self%atten_se(w,s) ) cycle
            if( l_model )then
                ! a rejected window has no model attenuation
                if( self%l_rejected(w) ) cycle
                y = self%atten_star(w,s)
            else
                y = self%atten(w,s)
            endif
            x   = self%dose(w) - self%dose(1)
            wgt = (self%atten(w,s) / max(DTINY, self%atten_se(w,s)))**2
            sxx = sxx + wgt * x * x
            sxy = sxy + wgt * x * log(y)
            gradient(w) = wgt * x / self%atten(w,s)
            nused = nused + 1
        enddo
        if( nused < 2 .or. sxx <= DTINY ) return
        m   = sxy / sxx
        gradient = gradient / sxx
        sem = sqrt(max(0._dp, dot_product(gradient, matmul(self%atten_cov(:,:,s), gradient))))
        if( m >= -SE_RATIO * sem ) return
        ne    = -0.5_dp / m
        ne_se = 0.5_dp * sem / (m * m)
    end subroutine shell_exposure

    !> Critical exposure between neighbouring windows w and w+1 at shell s, with no model:
    !! Ne = (D_{w+1} - D_w) / (2 ln(a_w/a_{w+1})), with the joint attenuation covariance. Fading
    !! exponential in dose gives the same value for every pair. NaN unless both windows resolve the shell
    !! and it fades by more than twice its error.
    subroutine step_exposure( self, w, s, ne, ne_se )
        use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan
        class(dose_series), intent(in)  :: self
        integer,            intent(in)  :: w, s
        real(dp),           intent(out) :: ne, ne_se
        real(dp) :: r1, r2, dl, sdl
        integer  :: v
        ne    = ieee_value(1._dp, ieee_quiet_nan)
        ne_se = ne
        if( s < self%klo .or. s > self%khi ) return
        if( w < 1 .or. w >= self%nwin ) return
        ! window 1 is the gauge, a_1 = 1 exactly
        do v = max(2, w), w + 1
            if( self%atten_se(v,s) >= huge(1._dp) ) return
            if( self%atten(v,s) <= SE_RATIO * self%atten_se(v,s) ) return
        enddo
        r1 = self%atten_se(w,s)   / self%atten(w,s)
        r2 = self%atten_se(w+1,s) / self%atten(w+1,s)
        dl  = log(self%atten(w,s) / self%atten(w+1,s))
        sdl = sqrt(max(0._dp, r1*r1 + r2*r2 - &
            &2._dp*self%atten_cov(w,w+1,s)/(self%atten(w,s)*self%atten(w+1,s))))
        if( dl <= SE_RATIO * sdl ) return
        ne    = 0.5_dp * (self%dose(w+1) - self%dose(w)) / dl
        ne_se = ne * sdl / dl
    end subroutine step_exposure

    !> attenuation of a window at a shell under the shrunk parameters; none for a rejected window
    real(dp) function atten_star( self, w, s ) result( a )
        class(dose_series), intent(in) :: self
        integer,            intent(in) :: w, s
        if( self%l_rejected(w) )then
            a = 0._dp
        else
            a = max(DTINY, self%cw(w)) * bfac_envelope(self%bw_star(w), self%kfreq(s)**2)
        endif
    end function atten_star

    !> attenuation of window w at shell s under one accumulation scheme; the control is gauged at
    !! window 1 like the model, and the flat mean is the model at beta = 0
    real(dp) function scheme_atten( self, scheme, w, s ) result( a )
        class(dose_series), intent(in) :: self
        integer,            intent(in) :: scheme, w, s
        real(dp) :: ne
        select case(scheme)
            case(DOSE_WEIGHTS_MODEL)
                a = self%atten_star(w,s)
            case(DOSE_WEIGHTS_CONTROL)
                ne = 0.245_dp * self%kfreq(s)**(-1.665_dp) + 2.81_dp
                a  = exp(-(self%dose(w) - self%dose(1)) / (2._dp*ne))
            case default
                a = 1._dp
        end select
    end function scheme_atten

    !> exp(-B k^2/4), with the argument clamped so a diverging fit cannot overflow
    pure real(dp) function bfac_envelope( b, k2 ) result( e )
        real(dp), intent(in) :: b, k2
        e = exp(max(-700._dp, min(700._dp, -b * k2 / 4._dp)))
    end function bfac_envelope

    !> the matched filter of every scheme under its attenuation and the windows' noise, normalised to sum one
    subroutine set_weights( self )
        class(dose_series), intent(inout) :: self
        real(dp) :: a(self%nwin), den
        integer  :: s, w, k
        do k = 1, N_SCHEMES
            do s = self%klo, self%khi
                den = 0._dp
                do w = 1, self%nwin
                    a(w) = self%scheme_atten(k, w, s) * self%prec(w,s)
                    den  = den + a(w)
                enddo
                if( den > DTINY )then
                    self%weight(:,s,k) = a / den
                else
                    self%weight(:,s,k) = 1._dp / real(self%nwin, dp)
                endif
            enddo
            ! beyond the band the last in-band shell's weights continue
            do s = 1, self%klo - 1
                self%weight(:,s,k) = self%weight(:,self%klo,k)
            enddo
            do s = self%khi + 1, self%nsh
                self%weight(:,s,k) = self%weight(:,self%khi,k)
            enddo
        enddo
        do s = self%klo, self%khi
            self%ess(s)  = kish_ess(self%weight(:,s,DOSE_WEIGHTS_MODEL))
            self%abar(s) = 0._dp
            do w = 1, self%nwin
                self%abar(s) = self%abar(s) + self%weight(w,s,DOSE_WEIGHTS_MODEL) * self%atten_star(w,s)
            enddo
        enddo
        do s = 1, self%klo - 1
            self%ess(s)  = self%ess(self%klo)
            self%abar(s) = self%abar(self%klo)
        enddo
        do s = self%khi + 1, self%nsh
            self%ess(s)  = self%ess(self%khi)
            self%abar(s) = self%abar(self%khi)
        enddo
    end subroutine set_weights

    !> Cross-half prediction error under each attenuation scheme, whitened by the target noise.
    subroutine cross_validate( self )
        class(dose_series), intent(inout) :: self
        real(dp) :: a(self%nwin), g(self%nwin), errs, den, eo, ee, at
        integer  :: s, w, v, t, icv
        do s = self%klo, self%khi
            do icv = 1, N_SCHEMES
                do w = 1, self%nwin
                    a(w) = self%scheme_atten(icv, w, s)
                enddo
                den = 0._dp
                do v = 1, self%nwin
                    den = den + a(v)*a(v) * self%prec(v,s)
                enddo
                if( den <= DTINY )then
                    self%cv(icv,s) = 0._dp
                    cycle
                endif
                do v = 1, self%nwin
                    g(v) = a(v) * self%prec(v,s) / den
                enddo
                errs = 0._dp
                do t = 1, self%nwin
                    at = a(t)
                    ! predicting the odd half from the even halves, and the reverse
                    eo = self%gram(self%nwin+t,self%nwin+t,s)
                    ee = self%gram(t,t,s)
                    do v = 1, self%nwin
                        eo = eo - 2._dp * at * g(v) * self%gram(v,self%nwin+t,s)
                        ee = ee - 2._dp * at * g(v) * self%gram(t,self%nwin+v,s)
                    enddo
                    do v = 1, self%nwin
                        do w = 1, self%nwin
                            eo = eo + at*at * g(v) * g(w) * self%gram(v,w,s)
                            ee = ee + at*at * g(v) * g(w) * self%gram(self%nwin+v,self%nwin+w,s)
                        enddo
                    enddo
                    errs = errs + 0.5_dp * self%prec(t,s) * (eo + ee)
                enddo
                self%cv(icv,s) = errs / real(self%nwin, dp)
            enddo
        enddo
    end subroutine cross_validate

    !> the model against the control, paired coefficient by coefficient over the masked halves; the
    !! per-shell mean of the difference is the difference of their Gram-based held-out errors
    subroutine paired_cv( self, em, om, shell )
        class(dose_series), intent(inout) :: self
        complex(sp),        intent(in)    :: em(:,:), om(:,:)   !< (ncoeff,nwin) masked halves
        integer,            intent(in)    :: shell(:)
        integer,     allocatable :: offset(:), order(:)
        real(dp)    :: am(self%nwin), ac(self%nwin), gm(self%nwin), gc(self%nwin), d, rm, rc, dm, dc
        complex(dp) :: e(self%nwin), o(self%nwin), sem, som, sec, soc
        integer     :: s, a, i, w, t
        if( .not. self%l_exists ) THROW_HARD('Dose series not constructed')
        if( self%khi < self%klo ) THROW_HARD('Fit the series before pairing its cross-validation')
        if( size(em,1) /= size(shell) .or. size(om,1) /= size(shell) ) THROW_HARD('Coefficient count mismatch')
        if( size(em,2) /= self%nwin .or. size(om,2) /= self%nwin ) THROW_HARD('Window count mismatch')
        call self%shell_buckets(shell, offset, order)
        self%dsum = 0._dp
        self%paired_delta = 0._dp
        !$omp parallel do default(shared) private(s,a,i,w,t,am,ac,gm,gc,d,rm,rc,dm,dc,e,o,sem,som,sec,soc)&
        !$omp schedule(dynamic)
        do s = self%klo, self%khi
            ! the signal estimates of (18) for both schemes
            dm = 0._dp; dc = 0._dp
            do w = 1, self%nwin
                am(w) = self%scheme_atten(DOSE_WEIGHTS_MODEL,   w, s)
                ac(w) = self%scheme_atten(DOSE_WEIGHTS_CONTROL, w, s)
                dm    = dm + am(w)*am(w)*self%prec(w,s)
                dc    = dc + ac(w)*ac(w)*self%prec(w,s)
            enddo
            if( dm <= DTINY .or. dc <= DTINY ) cycle
            do w = 1, self%nwin
                gm(w) = am(w) * self%prec(w,s) / dm
                gc(w) = ac(w) * self%prec(w,s) / dc
            enddo
            do a = offset(s), offset(s+1) - 1
                i = order(a)
                sem = (0._dp,0._dp); som = (0._dp,0._dp); sec = (0._dp,0._dp); soc = (0._dp,0._dp)
                do w = 1, self%nwin
                    e(w) = self%rscale * cmplx(em(i,w), kind=dp)
                    o(w) = self%rscale * cmplx(om(i,w), kind=dp)
                    sem  = sem + gm(w) * e(w)
                    som  = som + gm(w) * o(w)
                    sec  = sec + gc(w) * e(w)
                    soc  = soc + gc(w) * o(w)
                enddo
                rm = 0._dp; rc = 0._dp
                do t = 1, self%nwin
                    rm = rm + 0.5_dp*self%prec(t,s) * (abs(o(t) - am(t)*sem)**2 + abs(e(t) - am(t)*som)**2)
                    rc = rc + 0.5_dp*self%prec(t,s) * (abs(o(t) - ac(t)*sec)**2 + abs(e(t) - ac(t)*soc)**2)
                enddo
                d = (rm - rc) / real(self%nwin, dp)
                self%dsum(s) = self%dsum(s) + d
            enddo
        enddo
        !$omp end parallel do
        do s = self%klo, self%khi
            self%paired_delta(s) = self%band_difference(s, s)
        enddo
        self%l_paired = .true.
    end subroutine paired_cv

    !> Mean paired prediction-error difference over a band; negative favours the model descriptively.
    real(dp) function band_difference( self, k1, k2 ) result( delta )
        class(dose_series), intent(in) :: self
        integer,            intent(in) :: k1, k2
        real(dp) :: n
        delta = 0._dp
        n = sum(real(self%nshell(k1:k2), dp))
        if( n > 0._dp ) delta = sum(self%dsum(k1:k2)) / n
    end function band_difference

    !> apply one scheme's weights to the halves as given; a coefficient of the output is a weighted mean
    subroutine accumulate( self, eu, ou, shell, scheme, ae, ao )
        class(dose_series), intent(in)  :: self
        complex(sp),        intent(in)  :: eu(:,:), ou(:,:)   !< (ncoeff,nwin) halves as given
        integer,            intent(in)  :: shell(:)
        integer,            intent(in)  :: scheme
        complex(sp),        intent(out) :: ae(:), ao(:)
        complex(dp) :: se, so
        integer     :: i, s, w
        if( .not. self%l_exists ) THROW_HARD('Dose series not constructed')
        if( scheme < 1 .or. scheme > N_SCHEMES ) THROW_HARD('Unknown accumulation scheme')
        if( size(eu,2) /= self%nwin .or. size(ou,2) /= self%nwin ) THROW_HARD('Window count mismatch in accumulation')
        !$omp parallel do default(shared) private(i,s,w,se,so) schedule(static)
        do i = 1, size(shell)
            s  = min(max(shell(i), 1), self%nsh)
            se = (0._dp, 0._dp)
            so = (0._dp, 0._dp)
            do w = 1, self%nwin
                se = se + self%weight(w,s,scheme) * cmplx(eu(i,w), kind=dp)
                so = so + self%weight(w,s,scheme) * cmplx(ou(i,w), kind=dp)
            enddo
            ae(i) = cmplx(se, kind=sp)
            ao(i) = cmplx(so, kind=sp)
        enddo
        !$omp end parallel do
    end subroutine accumulate

    !> the zero-dose extrapolation (15), a per-shell rescaling of the accumulated map; past the band the
    !! factor of its last shell is held, so the noise there is not amplified without bound
    subroutine zerodose( self, a, shell, s0 )
        class(dose_series), intent(in)  :: self
        complex(sp),        intent(in)  :: a(:)
        integer,            intent(in)  :: shell(:)
        complex(sp),        intent(out) :: s0(:)
        real(dp) :: fac(self%nsh)
        integer  :: i, s, ss
        if( .not. self%l_exists ) THROW_HARD('Dose series not constructed')
        if( self%khi < self%klo ) THROW_HARD('Fit the series before extrapolating it')
        do s = 1, self%nsh
            ss     = min(s, self%khi)
            fac(s) = 1._dp / (max(DTINY, self%abar(ss)) * bfac_envelope(-self%b0, self%kfreq(ss)**2))
        enddo
        !$omp parallel do default(shared) private(i,s) schedule(static)
        do i = 1, size(shell)
            s     = min(max(shell(i), 1), self%nsh)
            s0(i) = cmplx(fac(s) * cmplx(a(i), kind=dp), kind=sp)
        enddo
        !$omp end parallel do
    end subroutine zerodose

    real(dp) function get_beta( self ) result( b )
        class(dose_series), intent(in) :: self
        b = self%beta
    end function get_beta

    subroutine get_band( self, klo, khi )
        class(dose_series), intent(in)  :: self
        integer,            intent(out) :: klo, khi
        klo = self%klo
        khi = self%khi
    end subroutine get_band

    !> per-shell weights of one scheme
    subroutine get_weights( self, w, scheme )
        class(dose_series),    intent(in)    :: self
        real(dp), allocatable, intent(inout) :: w(:,:)
        integer,               intent(in)    :: scheme
        if( scheme < 1 .or. scheme > N_SCHEMES ) THROW_HARD('Unknown accumulation scheme')
        if( allocated(w) ) deallocate(w)
        allocate(w(self%nwin,self%nsh), source=self%weight(:,:,scheme))
    end subroutine get_weights

    !> Per-shell cross-half errors and the paired model-minus-control mean difference; not a significance test.
    subroutine get_cv( self, cv, delta )
        class(dose_series),    intent(in)    :: self
        real(dp), allocatable, intent(inout) :: cv(:,:), delta(:)
        if( allocated(cv) ) deallocate(cv)
        if( allocated(delta) ) deallocate(delta)
        allocate(cv(N_SCHEMES,self%nsh), source=self%cv)
        allocate(delta(self%nsh), source=self%paired_delta)
    end subroutine get_cv

    !> the model-free attenuation of every window and shell with its standard error, and each window's
    !! Guinier rms: the residual of its B-factor curve in units of those errors
    subroutine get_attenuation( self, atten, se, grms )
        class(dose_series),    intent(in)    :: self
        real(dp), allocatable, intent(inout) :: atten(:,:), se(:,:), grms(:)
        if( allocated(atten) ) deallocate(atten)
        if( allocated(se)    ) deallocate(se)
        if( allocated(grms)  ) deallocate(grms)
        allocate(atten(self%nwin,self%nsh), source=self%atten)
        allocate(se(self%nwin,self%nsh),    source=self%atten_se)
        allocate(grms(self%nwin),           source=self%grms)
    end subroutine get_attenuation

    !> The windows' fit, reported against the sharpest window: c and B divide and shift by its values, which
    !! leaves every ratio between windows, the line's slopes and the weights as they are
    subroutine write_model_table( self, fname )
        class(dose_series), intent(in) :: self
        class(string),      intent(in) :: fname
        character, parameter :: TAB = achar(9)
        real(dp) :: bref, cref
        integer  :: funit, w, io_stat, iref
        iref = self%sharpest()
        bref = self%bw_star(iref)
        cref = self%cw(iref)
        call fopen(funit, file=fname, status='replace', action='write', form='formatted')
        write(funit,'(A,I0,A)',iostat=io_stat) '# dose_series model; B in A^2 and c relative to window ', iref, &
            &', the sharpest; dose in the scaled coordinate'
        call fileiochk('Writing dose series model header', io_stat)
        write(funit,'(A)',iostat=io_stat) 'window'//TAB//'dose'//TAB//'c_w'//TAB//'B_w'//TAB//'se_B'// &
            &TAB//'B_line'//TAB//'B_star'//TAB//'z'//TAB//'guinier_rms'//TAB//'used'
        call fileiochk('Writing dose series model columns', io_stat)
        do w = 1, self%nwin
            write(funit,'(I0,8(A,ES16.8E3),A,I0)',iostat=io_stat) w, TAB, self%dose(w), TAB, self%cw(w)/cref, &
                &TAB, self%bw(w) - bref, TAB, self%rel_se(w, iref), TAB, self%b0 + self%beta*self%dose(w) - bref, &
                &TAB, self%bw_star(w) - bref, TAB, self%zw(w), TAB, self%grms(w), TAB, merge(0, 1, self%l_rejected(w))
            call fileiochk('Writing dose series model row', io_stat)
        enddo
        write(funit,'(A,4(A,ES16.8E3))',iostat=io_stat) '# line B0 beta se_beta tau', TAB, self%b0 - bref, &
            &TAB, self%beta, TAB, self%beta_se, TAB, self%tau
        call fileiochk('Writing dose series line', io_stat)
        write(funit,'(A,4(A,ES16.8E3))',iostat=io_stat) '# contrast lnc0 gamma se_gamma tau', TAB, self%lnc0 - log(cref), &
            &TAB, self%gamma, TAB, self%gamma_se, TAB, self%tau_c
        call fileiochk('Writing dose series contrast line', io_stat)
        call fclose(funit)
    end subroutine write_model_table

    subroutine write_attenuation_table( self, fname )
        class(dose_series), intent(in) :: self
        class(string),      intent(in) :: fname
        character, parameter :: TAB = achar(9)
        integer :: funit, s, w, io_stat
        call fopen(funit, file=fname, status='replace', action='write', form='formatted')
        write(funit,'(A)',advance='no',iostat=io_stat) 'shell'//TAB//'k_inv_A'
        call fileiochk('Writing attenuation header', io_stat)
        do w = 1, self%nwin
            write(funit,'(A)',advance='no',iostat=io_stat) TAB//'a_'//int2str(w)//TAB//'se_'//int2str(w)
            call fileiochk('Writing attenuation columns', io_stat)
        enddo
        write(funit,'(A)',iostat=io_stat) ''
        call fileiochk('Writing attenuation newline', io_stat)
        do s = self%klo, self%khi
            write(funit,'(I0,A,ES16.8E3)',advance='no',iostat=io_stat) s, TAB, self%kfreq(s)
            call fileiochk('Writing attenuation row', io_stat)
            do w = 1, self%nwin
                write(funit,'(2(A,ES16.8E3))',advance='no',iostat=io_stat) TAB, self%atten(w,s), &
                    &TAB, self%atten_se(w,s)
                call fileiochk('Writing attenuation cell', io_stat)
            enddo
            write(funit,'(A)',iostat=io_stat) ''
            call fileiochk('Writing attenuation row end', io_stat)
        enddo
        call fclose(funit)
    end subroutine write_attenuation_table

    subroutine write_weights_table( self, fname )
        class(dose_series), intent(in) :: self
        class(string),      intent(in) :: fname
        character, parameter :: TAB = achar(9)
        integer :: funit, s, w, io_stat
        call fopen(funit, file=fname, status='replace', action='write', form='formatted')
        write(funit,'(A)',advance='no',iostat=io_stat) 'shell'//TAB//'k_inv_A'//TAB//'neff'//TAB//'abar'
        call fileiochk('Writing weights header', io_stat)
        do w = 1, self%nwin
            write(funit,'(A)',advance='no',iostat=io_stat) TAB//'w_'//int2str(w)
            call fileiochk('Writing weights columns', io_stat)
        enddo
        write(funit,'(A)',iostat=io_stat) ''
        call fileiochk('Writing weights newline', io_stat)
        do s = 1, self%nsh
            write(funit,'(I0,3(A,ES16.8E3))',advance='no',iostat=io_stat) s, TAB, self%kfreq(s), &
                &TAB, self%ess(s), TAB, self%abar(s)
            call fileiochk('Writing weights row', io_stat)
            do w = 1, self%nwin
                write(funit,'(A,ES16.8E3)',advance='no',iostat=io_stat) TAB, self%weight(w,s,DOSE_WEIGHTS_MODEL)
                call fileiochk('Writing weights cell', io_stat)
            enddo
            write(funit,'(A)',iostat=io_stat) ''
            call fileiochk('Writing weights row end', io_stat)
        enddo
        call fclose(funit)
    end subroutine write_weights_table

    subroutine write_cv_table( self, fname )
        class(dose_series), intent(in) :: self
        class(string),      intent(in) :: fname
        character, parameter :: TAB = achar(9)
        integer :: funit, s, io_stat
        call fopen(funit, file=fname, status='replace', action='write', form='formatted')
        write(funit,'(A)',iostat=io_stat) '# cross-half prediction error: 1 is noise-free, 2 the window alone; parameters fitted on both halves;'// &
            &' paired differences are descriptive, with no calibrated significance or winner'
        call fileiochk('Writing cv header', io_stat)
        write(funit,'(A)',iostat=io_stat) 'shell'//TAB//'k_inv_A'//TAB//'model'//TAB//'grant_grigorieff'//TAB//'flat'// &
            &TAB//'error_model_minus_control'
        call fileiochk('Writing cv columns', io_stat)
        do s = self%klo, self%khi
            write(funit,'(I0,5(A,ES16.8E3))',iostat=io_stat) s, TAB, self%kfreq(s), &
                &TAB, self%cv(DOSE_WEIGHTS_MODEL,s), TAB, self%cv(DOSE_WEIGHTS_CONTROL,s), &
                &TAB, self%cv(DOSE_WEIGHTS_FLAT,s), TAB, self%paired_delta(s)
            call fileiochk('Writing cv row', io_stat)
        enddo
        call fclose(funit)
    end subroutine write_cv_table

    !> Per band shell: the model-free critical exposure and its error, the lines' and Grant & Grigorieff's, the
    !! per-window model's through the same estimator, then per window the measured attenuation, its error, the
    !! model's and theirs. The lines fade as exp(-(gamma + beta k^2/4) D): critical exposure 1/(2 gamma + beta k^2/2).
    subroutine write_decay_table( self, fname )
        use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan
        class(dose_series), intent(in) :: self
        class(string),      intent(in) :: fname
        character, parameter :: TAB = achar(9)
        real(dp) :: nan, ne_model, ne_gg, se, am, rate
        integer  :: funit, s, w, io_stat
        nan = ieee_value(1._dp, ieee_quiet_nan)
        call fopen(funit, file=fname, status='replace', action='write', form='formatted')
        write(funit,'(A)',iostat=io_stat) '# critical exposure in the scaled dose coordinate (e/A^2); NaN where not measured'
        call fileiochk('Writing decay header', io_stat)
        write(funit,'(A)',advance='no',iostat=io_stat) 'shell'//TAB//'k_inv_A'//TAB//'ne_measured'//TAB//'ne_se'// &
            &TAB//'ne_model'//TAB//'ne_gg'//TAB//'ne_model_est'
        call fileiochk('Writing decay columns', io_stat)
        do w = 1, self%nwin
            write(funit,'(A)',advance='no',iostat=io_stat) TAB//'a_'//int2str(w)//TAB//'se_'//int2str(w)// &
                &TAB//'model_'//int2str(w)//TAB//'gg_'//int2str(w)
            call fileiochk('Writing decay window columns', io_stat)
        enddo
        write(funit,'(A)',iostat=io_stat) ''
        call fileiochk('Writing decay newline', io_stat)
        do s = self%klo, self%khi
            ne_model = nan
            rate     = self%gamma + self%beta * self%kfreq(s)**2 / 4._dp
            if( rate > DTINY ) ne_model = 0.5_dp / rate
            ne_gg = 0.245_dp * self%kfreq(s)**(-1.665_dp) + 2.81_dp
            write(funit,'(I0,6(A,ES16.8E3))',advance='no',iostat=io_stat) s, TAB, self%kfreq(s), &
                &TAB, self%ne_shell(s), TAB, self%ne_se(s), TAB, ne_model, TAB, ne_gg, TAB, self%ne_model_est(s)
            call fileiochk('Writing decay row', io_stat)
            do w = 1, self%nwin
                se = self%atten_se(w,s)
                if( se >= huge(1._dp) ) se = nan
                am = self%atten_star(w,s)
                if( self%l_rejected(w) ) am = nan
                write(funit,'(4(A,ES16.8E3))',advance='no',iostat=io_stat) TAB, self%atten(w,s), TAB, se, &
                    &TAB, am, TAB, self%scheme_atten(DOSE_WEIGHTS_CONTROL, w, s)
                call fileiochk('Writing decay cell', io_stat)
            enddo
            write(funit,'(A)',iostat=io_stat) ''
            call fileiochk('Writing decay row end', io_stat)
        enddo
        call fclose(funit)
    end subroutine write_decay_table

    !> Per window: its dose and, at the shell nearest each requested resolution, the model's and the Grant &
    !! Grigorieff control's accumulation weights. A resolution past Nyquist gets NaN columns.
    !> Per pair of neighbouring windows: its mean dose and, at the band shells nearest res (A), the critical
    !! exposure between the two with its error and Grant & Grigorieff's; NaN where not measured
    subroutine write_step_exposure( self, fname, res )
        use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan
        class(dose_series), intent(in) :: self
        class(string),      intent(in) :: fname
        real(dp),           intent(in) :: res(:)
        character, parameter :: TAB = achar(9)
        real(dp) :: nan, ne, se, ne_gg
        integer  :: funit, w, i, io_stat, sres(size(res))
        nan = ieee_value(1._dp, ieee_quiet_nan)
        do i = 1, size(res)
            sres(i) = nint(1._dp / (res(i) * self%kfreq(1)))
            if( sres(i) < self%klo .or. sres(i) > self%khi ) sres(i) = 0
        enddo
        call fopen(funit, file=fname, status='replace', action='write', form='formatted')
        write(funit,'(A)',advance='no',iostat=io_stat) 'pair'//TAB//'dose'
        call fileiochk('Writing step exposure header', io_stat)
        do i = 1, size(res)
            write(funit,'(A)',advance='no',iostat=io_stat) TAB//'ne_'//int2str(nint(res(i)))//'A'// &
                &TAB//'se_'//int2str(nint(res(i)))//'A'//TAB//'gg_'//int2str(nint(res(i)))//'A'
            call fileiochk('Writing step exposure columns', io_stat)
        enddo
        write(funit,'(A)',iostat=io_stat) ''
        call fileiochk('Writing step exposure newline', io_stat)
        do w = 1, self%nwin - 1
            write(funit,'(I0,A,ES16.8E3)',advance='no',iostat=io_stat) w, TAB, 0.5_dp * (self%dose(w) + self%dose(w+1))
            call fileiochk('Writing step exposure row', io_stat)
            do i = 1, size(res)
                ne = nan; se = nan; ne_gg = nan
                if( sres(i) > 0 )then
                    call self%step_exposure(w, sres(i), ne, se)
                    ne_gg = 0.245_dp * self%kfreq(sres(i))**(-1.665_dp) + 2.81_dp
                endif
                write(funit,'(3(A,ES16.8E3))',advance='no',iostat=io_stat) TAB, ne, TAB, se, TAB, ne_gg
                call fileiochk('Writing step exposure cell', io_stat)
            enddo
            write(funit,'(A)',iostat=io_stat) ''
            call fileiochk('Writing step exposure row end', io_stat)
        enddo
        call fclose(funit)
    end subroutine write_step_exposure

    subroutine write_resolution_weights( self, fname, res )
        use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan
        class(dose_series), intent(in) :: self
        class(string),      intent(in) :: fname
        real(dp),           intent(in) :: res(:)
        character, parameter :: TAB = achar(9)
        real(dp) :: nan
        integer  :: funit, w, i, io_stat, sres(size(res))
        nan = ieee_value(1._dp, ieee_quiet_nan)
        do i = 1, size(res)
            sres(i) = nint(1._dp / (res(i) * self%kfreq(1)))
            if( sres(i) < 1 .or. sres(i) > self%nsh ) sres(i) = 0
        enddo
        call fopen(funit, file=fname, status='replace', action='write', form='formatted')
        write(funit,'(A)',advance='no',iostat=io_stat) 'window'//TAB//'dose'
        call fileiochk('Writing resolution weights header', io_stat)
        do i = 1, size(res)
            write(funit,'(A)',advance='no',iostat=io_stat) TAB//'model_'//int2str(nint(res(i)))//'A'// &
                &TAB//'gg_'//int2str(nint(res(i)))//'A'
            call fileiochk('Writing resolution weights columns', io_stat)
        enddo
        write(funit,'(A)',iostat=io_stat) ''
        call fileiochk('Writing resolution weights newline', io_stat)
        do w = 1, self%nwin
            write(funit,'(I0,A,ES16.8E3)',advance='no',iostat=io_stat) w, TAB, self%dose(w)
            call fileiochk('Writing resolution weights row', io_stat)
            do i = 1, size(res)
                if( sres(i) == 0 )then
                    write(funit,'(2(A,ES16.8E3))',advance='no',iostat=io_stat) TAB, nan, TAB, nan
                else
                    write(funit,'(2(A,ES16.8E3))',advance='no',iostat=io_stat) TAB, &
                        &self%weight(w,sres(i),DOSE_WEIGHTS_MODEL), TAB, self%weight(w,sres(i),DOSE_WEIGHTS_CONTROL)
                endif
                call fileiochk('Writing resolution weights cell', io_stat)
            enddo
            write(funit,'(A)',iostat=io_stat) ''
            call fileiochk('Writing resolution weights row end', io_stat)
        enddo
        call fclose(funit)
    end subroutine write_resolution_weights

    !> the fitted model as the dose curve the movie re-extraction reads: the windows' contrast and shrunk
    !! B-factor, a rejected window flagged and zeroed; fractions and physical dose rate from the manifest
    subroutine export_curve( self, curve, ffrom, fto, dose_per_fraction )
        class(dose_series), intent(in)    :: self
        class(dose_curve),  intent(inout) :: curve
        integer,            intent(in)    :: ffrom(:), fto(:)
        real(dp),           intent(in)    :: dose_per_fraction
        integer :: w
        if( .not. self%l_exists ) THROW_HARD('Dose series not constructed')
        if( size(ffrom) /= self%nwin .or. size(fto) /= self%nwin ) THROW_HARD('One fraction range per window')
        call curve%new(self%nwin, self%kv, dose_per_fraction, 1._dp/self%kfreq(self%klo), 1._dp/self%kfreq(self%khi), &
            &self%b0, self%beta, self%lnc0, self%gamma)
        do w = 1, self%nwin
            if( self%l_rejected(w) )then
                call curve%set_knot(w, ffrom(w), fto(w), self%dose(w), 0._dp, 0._dp, .false.)
            else
                call curve%set_knot(w, ffrom(w), fto(w), self%dose(w), self%cw(w), self%bw_star(w), .true.)
            endif
        enddo
    end subroutine export_curve

    !> the Dose series: lines of the log, per window, for the line, and per resolution band
    subroutine log_readout( self )
        class(dose_series), intent(in) :: self
        real(dp), parameter :: RES_PROBE(3) = [2._dp, 3._dp, 10._dp]
        real(dp) :: ne, equiv(3), n, nsw, ess_b, abar_b, cvb(N_SCHEMES), delta, bref, cref
        integer  :: w, i, j, nb, nsb, k1, k2, s, iref
        iref = self%sharpest()
        bref = self%bw_star(iref)
        cref = self%cw(iref)
        write(logfhandle,'(A,I0,A)') '>>> DOSE SERIES: PER WINDOW, c AND B RELATIVE TO WINDOW ', iref, ', THE SHARPEST'
        do w = 1, self%nwin
            if( self%l_rejected(w) )then
                write(logfhandle,'(A,I3,A,F11.4,A)') '>>> Dose series: window ', w, '  D=', self%dose(w), &
                    &'  REJECTED: no fit of its attenuation, left out of the model weights'
                cycle
            endif
            write(logfhandle,'(A,I3,8(A,F11.4))') '>>> Dose series: window ', w, '  D=', self%dose(w), &
                &'  c=', self%cw(w)/cref, '  B=', self%bw(w) - bref, ' +-', self%rel_se(w, iref), &
                &'  B(D)=', self%b0 + self%beta*self%dose(w) - bref, '  B*=', self%bw_star(w) - bref, '  z=', self%zw(w), &
                &'  guinier_rms=', self%grms(w)
        enddo
        do i = 1, 3
            ne       = 0.245_dp * (1._dp/RES_PROBE(i))**(-1.665_dp) + 2.81_dp
            equiv(i) = 2._dp / (ne * (1._dp/RES_PROBE(i))**2)
        enddo
        write(logfhandle,'(A,4(A,F11.4))') '>>> Dose series: line', '  B0=', self%b0 - bref, &
            &'  beta=', self%beta, ' +-', self%beta_se, '  tau=', self%tau
        write(logfhandle,'(A,4(A,F11.6))') '>>> Dose series: contrast, ln c = lnc0 - gamma D', '  lnc0=', self%lnc0 - log(cref), &
            &'  gamma=', self%gamma, ' +-', self%gamma_se, '  tau=', self%tau_c
        write(logfhandle,'(A,3(A,F10.3))') '>>> Dose series: Grant & Grigorieff beta at 2/3/10 A', &
            &'  ', equiv(1), '  ', equiv(2), '  ', equiv(3)
        write(logfhandle,'(A,2(A,I0))') '>>> Dose series: band shells', '  ', self%klo, ' to ', self%khi
        nsb = self%khi - self%klo + 1
        nb  = min(N_BANDS, nsb)
        write(logfhandle,'(A)') '>>> DOSE SERIES: PER RESOLUTION BAND (cv 1 is noise-free, 2 the window alone;'// &
            &' paired mean differences are descriptive, not significance tests)'
        do j = 1, nb
            k1 = self%klo + ((j-1)*nsb)/nb
            k2 = self%klo + (j*nsb)/nb - 1
            n = 0._dp; ess_b = 0._dp; abar_b = 0._dp; cvb = 0._dp
            do s = k1, k2
                nsw    = real(self%nshell(s), dp)
                n      = n      + nsw
                ess_b  = ess_b  + nsw * self%ess(s)
                abar_b = abar_b + nsw * self%abar(s)
                cvb    = cvb    + nsw * self%cv(:,s)
            enddo
            if( n > 0._dp )then
                ess_b  = ess_b  / n
                abar_b = abar_b / n
                cvb    = cvb    / n
            endif
            delta = 0._dp
            if( self%l_paired ) delta = self%band_difference(k1, k2)
            write(logfhandle,'(A,F7.2,A,F7.2,A,F7.3,A,F7.3,3(A,F7.4),A,ES12.4)') '>>> Dose series: band ', &
                &1._dp/self%kfreq(k1), '-', 1._dp/self%kfreq(k2), ' A  neff=', ess_b, '  abar=', abar_b, &
                &'  cv model=', cvb(DOSE_WEIGHTS_MODEL), ' control=', cvb(DOSE_WEIGHTS_CONTROL), &
                &' flat=', cvb(DOSE_WEIGHTS_FLAT), '  paired model-control=', delta
        enddo
    end subroutine log_readout

    !> The fitted window with the least B*: the reference of the reports. Window 1 stays the gauge of the fit;
    !! early motion can leave it blurrier than a later window, whose B would then read negative against it.
    integer function sharpest( self ) result( iref )
        class(dose_series), intent(in) :: self
        integer :: w
        iref = 1
        do w = 2, self%nwin
            if( self%l_rejected(w) ) cycle
            if( self%bw_star(w) < self%bw_star(iref) ) iref = w
        enddo
    end function sharpest

    !> Error of B_w minus B*_iref, propagating the joint fit covariance with line weights and tau fixed.
    real(dp) function rel_se( self, w, iref ) result( se )
        class(dose_series), intent(in) :: self
        integer,            intent(in) :: w, iref
        real(dp) :: gradient(self%nwin), weight(self%nwin), s0, s1, s2, det, alpha
        integer :: v
        se = huge(1._dp)
        if( self%l_rejected(w) ) return
        gradient = 0._dp
        if( iref /= 1 )then
            gradient(iref) = 1._dp
            if( self%l_line .and. self%l_line_fitted .and. count(.not. self%l_rejected(2:)) >= 2 )then
                weight = 0._dp
                do v = 2, self%nwin
                    if( self%l_rejected(v) ) cycle
                    weight(v) = 1._dp / (self%bw_se(v)**2 + self%tau**2)
                enddo
                s0 = sum(weight); s1 = sum(weight*self%dose); s2 = sum(weight*self%dose**2)
                det = s0*s2 - s1*s1
                if( det > DTINY )then
                    alpha = self%tau**2 / (self%tau**2 + self%bw_se(iref)**2)
                    gradient = (1._dp-alpha) * weight * &
                        &(s2 - s1*self%dose(iref) + (s0*self%dose(iref)-s1)*self%dose) / det
                    gradient(iref) = gradient(iref) + alpha
                endif
            endif
        endif
        gradient = -gradient
        gradient(w) = gradient(w) + 1._dp
        se = sqrt(max(0._dp, dot_product(gradient, matmul(self%bw_cov, gradient))))
    end function rel_se

    subroutine kill( self )
        class(dose_series), intent(inout) :: self
        if( .not. self%l_exists ) return
        if( allocated(self%nshell)   ) deallocate(self%nshell)
        if( allocated(self%kfreq)    ) deallocate(self%kfreq)
        if( allocated(self%dose)     ) deallocate(self%dose)
        if( allocated(self%ddose)    ) deallocate(self%ddose)
        if( allocated(self%sigma2)   ) deallocate(self%sigma2)
        if( allocated(self%prec)     ) deallocate(self%prec)
        if( allocated(self%gram)     ) deallocate(self%gram)
        if( allocated(self%cvar)     ) deallocate(self%cvar)
        if( allocated(self%atten)    ) deallocate(self%atten)
        if( allocated(self%atten_se) ) deallocate(self%atten_se)
        if( allocated(self%atten_cov) ) deallocate(self%atten_cov)
        if( allocated(self%bw_cov)   ) deallocate(self%bw_cov)
        if( allocated(self%sigpow)   ) deallocate(self%sigpow)
        if( allocated(self%cw)       ) deallocate(self%cw)
        if( allocated(self%cw_se)    ) deallocate(self%cw_se)
        if( allocated(self%bw)       ) deallocate(self%bw)
        if( allocated(self%bw_se)    ) deallocate(self%bw_se)
        if( allocated(self%bw_star)  ) deallocate(self%bw_star)
        if( allocated(self%zw)       ) deallocate(self%zw)
        if( allocated(self%grms)     ) deallocate(self%grms)
        if( allocated(self%l_rejected) ) deallocate(self%l_rejected)
        if( allocated(self%weight)   ) deallocate(self%weight)
        if( allocated(self%ess)      ) deallocate(self%ess)
        if( allocated(self%abar)     ) deallocate(self%abar)
        if( allocated(self%cv)       ) deallocate(self%cv)
        if( allocated(self%dsum)     ) deallocate(self%dsum)
        if( allocated(self%paired_delta) ) deallocate(self%paired_delta)
        if( allocated(self%ne_shell) ) deallocate(self%ne_shell)
        if( allocated(self%ne_se)    ) deallocate(self%ne_se)
        if( allocated(self%ne_model_est) ) deallocate(self%ne_model_est)
        self%nwin = 0; self%nsh = 0; self%klo = 0; self%khi = 0
        self%rscale = 1._dp; self%mask_frac = 1._dp; self%shell_corr = 1._dp
        self%b0 = 0._dp; self%beta = 0._dp; self%beta_se = 0._dp; self%tau = 0._dp
        self%lnc0 = 0._dp; self%gamma = 0._dp; self%gamma_se = 0._dp; self%tau_c = 0._dp
        self%l_line_fitted = .false.
        self%l_paired      = .false.
        self%l_exists      = .false.
    end subroutine kill

end module simple_dose_series
