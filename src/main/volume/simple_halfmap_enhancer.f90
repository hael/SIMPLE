!@descr: model-free local Wiener sharpening of independent cryo-EM half-maps with tapered overlap-add
module simple_halfmap_enhancer
use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
use simple_defs,          only: DP, DPI, HPLIM_GUINIER, TINY
use simple_image,         only: image
use simple_estimate_ssnr, only: fsc2optlp
use simple_butterworth,   only: butterworth_filter
use simple_math_ft,       only: calc_fourier_index
use simple_stat,          only: median
use simple_error,         only: simple_exception
implicit none
private
#include "simple_local_flags.inc"

public :: halfmap_enhancer

real, parameter :: FSC_CUTOFF            = 0.143
real, parameter :: ENVELOPE_FULL_SUPPORT = 0.05
real, parameter :: MAX_SHARPEN_B          = 200.
real, parameter :: MAX_SPECTRAL_GAIN      = 8.
real, parameter :: TARGET_WINDOW_A        = 48.
real, parameter :: CORRELATED_NOISE_RATIO = 1.10
integer, parameter :: MIN_WINDOW          = 16
integer, parameter :: MAX_WINDOW          = 64

type :: halfmap_enhancer
    private
    logical :: l_exists = .false.
    integer :: ldim(3)  = 0
    integer :: window   = 0
    integer :: stride   = 0
    real    :: smpd     = 0.
    real    :: noise_mean = 0.
    real    :: noise_sdev = 0.
    real    :: correlated_sdev = 0.
    real, allocatable :: taper(:,:,:)
  contains
    procedure :: new
    procedure :: enhance
    procedure :: get_noise_stats
    procedure :: get_window
    procedure :: kill
end type halfmap_enhancer

contains

    subroutine new( self, ldim, smpd, window, stride )
        class(halfmap_enhancer), intent(inout) :: self
        integer,                 intent(in)    :: ldim(3)
        real,                    intent(in)    :: smpd
        integer, optional,       intent(in)    :: window, stride
        real, allocatable :: w1(:)
        integer :: i, j, k, win
        call self%kill()
        if( any(ldim < MIN_WINDOW) ) THROW_HARD('half-map dimensions must be at least 16 voxels; halfmap_enhancer%new')
        if( smpd <= 0. ) THROW_HARD('sampling distance must be positive; halfmap_enhancer%new')
        win = 2 * nint(0.5 * TARGET_WINDOW_A / smpd)
        win = max(MIN_WINDOW, min(MAX_WINDOW, win, minval(ldim)))
        if( mod(win, 2) /= 0 ) win = win - 1
        if( present(window) ) win = window
        if( win < MIN_WINDOW .or. win > minval(ldim) .or. mod(win, 2) /= 0 ) &
            &THROW_HARD('window must be even and fit every volume dimension; halfmap_enhancer%new')
        self%window = win
        self%stride = max(1, win / 2)
        if( present(stride) ) self%stride = stride
        if( self%stride < 1 .or. self%stride > win ) &
            &THROW_HARD('stride must lie between one and the window size; halfmap_enhancer%new')
        self%ldim = ldim
        self%smpd = smpd
        allocate(w1(win), self%taper(win,win,win))
        do i = 1, win
            w1(i) = sin(real(DPI) * (real(i) - 0.5) / real(win))
        enddo
        do k = 1, win
            do j = 1, win
                do i = 1, win
                    self%taper(i,j,k) = w1(i) * w1(j) * w1(k)
                enddo
            enddo
        enddo
        deallocate(w1)
        self%l_exists = .true.
    end subroutine new

    subroutine enhance( self, even, odd, enhanced, confidence, local_resolution, envelope )
        class(halfmap_enhancer), intent(inout) :: self
        class(image), target,     intent(in)    :: even, odd
        class(image),             intent(inout) :: enhanced, confidence, local_resolution, envelope
        type(image) :: even_win, odd_win, merged_win
        real, pointer :: even_r(:,:,:), odd_r(:,:,:)
        real, allocatable :: even_patch(:,:,:), odd_patch(:,:,:), merged_patch(:,:,:), filtered_patch(:,:,:)
        real, allocatable :: accum(:,:,:), norm(:,:,:), conf_accum(:,:,:), res_accum(:,:,:), env_accum(:,:,:)
        real, allocatable :: out(:,:,:), conf(:,:,:), locres(:,:,:), env(:,:,:)
        real, allocatable :: fsc(:), fsc_reg(:), wiener(:), lowpass(:), cap(:), res(:)
        integer, allocatable :: xs(:), ys(:), zs(:)
        real :: local_res, local_conf, env_conf, bfac
        integer :: ix, iy, iz, i0, j0, k0, i1, j1, k1, k, kcut, fromk, nyq
        if( .not. self%l_exists ) THROW_HARD('halfmap_enhancer must be constructed before enhance')
        if( even%is_ft() .or. odd%is_ft() ) THROW_HARD('halfmap_enhancer expects real-space half-maps')
        if( any(even%get_ldim() /= self%ldim) .or. any(odd%get_ldim() /= self%ldim) ) &
            &THROW_HARD('half-map dimensions do not match the enhancer')
        call even%get_rmat_ptr(even_r)
        call odd%get_rmat_ptr(odd_r)
        call estimate_noise(self, even_r, odd_r)
        xs = window_starts(self%ldim(1), self%window, self%stride)
        ys = window_starts(self%ldim(2), self%window, self%stride)
        zs = window_starts(self%ldim(3), self%window, self%stride)
        allocate(accum(self%ldim(1),self%ldim(2),self%ldim(3)), source=0.)
        allocate(norm(self%ldim(1),self%ldim(2),self%ldim(3)), source=0.)
        allocate(conf_accum(self%ldim(1),self%ldim(2),self%ldim(3)), source=0.)
        allocate(res_accum(self%ldim(1),self%ldim(2),self%ldim(3)), source=0.)
        allocate(env_accum(self%ldim(1),self%ldim(2),self%ldim(3)), source=0.)
        allocate(even_patch(self%window,self%window,self%window))
        allocate(odd_patch(self%window,self%window,self%window))
        allocate(merged_patch(self%window,self%window,self%window))
        allocate(filtered_patch(self%window,self%window,self%window))
        call even_win%new([self%window,self%window,self%window], self%smpd, wthreads=.false.)
        call odd_win%new([self%window,self%window,self%window], self%smpd, wthreads=.false.)
        call merged_win%new([self%window,self%window,self%window], self%smpd, wthreads=.false.)
        nyq = merged_win%get_filtsz()
        allocate(fsc(nyq), fsc_reg(nyq), wiener(nyq), lowpass(nyq), cap(nyq), res(nyq))
        res = merged_win%get_res()
        do iz = 1, size(zs)
            k0 = zs(iz)
            k1 = k0 + self%window - 1
            do iy = 1, size(ys)
                j0 = ys(iy)
                j1 = j0 + self%window - 1
                do ix = 1, size(xs)
                    i0 = xs(ix)
                    i1 = i0 + self%window - 1
                    even_patch  = (even_r(i0:i1,j0:j1,k0:k1) - self%noise_mean) * self%taper
                    odd_patch   = (odd_r(i0:i1,j0:j1,k0:k1)  - self%noise_mean) * self%taper
                    merged_patch = 0.5 * (even_patch + odd_patch)
                    env_conf = local_envelope_confidence(even_patch, odd_patch, self%noise_sdev, &
                        &self%correlated_sdev, self%taper)
                    call even_win%set_rmat(even_patch, .false.)
                    call odd_win%set_rmat(odd_patch, .false.)
                    call even_win%fft()
                    call odd_win%fft()
                    call even_win%fsc(odd_win, fsc)
                    call regularize_fsc(fsc, fsc_reg, kcut)
                    wiener = fsc2optlp(fsc_reg)
                    if( kcut > 0 )then
                        call butterworth_filter(kcut, lowpass)
                        wiener = wiener * lowpass
                        local_res = res(kcut)
                    else
                        wiener = 0.
                        local_res = 0.
                    endif
                    local_conf = env_conf * sum(wiener) / real(nyq)
                    call merged_win%set_rmat(merged_patch, .false.)
                    bfac = 0.
                    if( kcut > 0 )then
                        fromk = calc_fourier_index(HPLIM_GUINIER, self%window, self%smpd)
                        if( local_res < 5. .and. kcut - fromk >= 2 )then
                            bfac = merged_win%guinier_bfac(HPLIM_GUINIER, local_res)
                            if( .not. ieee_is_finite(bfac) ) bfac = 0.
                            bfac = max(-MAX_SHARPEN_B, min(0., bfac))
                        endif
                    endif
                    if( bfac < 0. )then
                        do k = 1, nyq
                            cap(k) = exp(-(bfac / 4.) * (real(min(k,kcut)) / &
                                &(real(self%window) * self%smpd))**2)
                        enddo
                        cap = min(MAX_SPECTRAL_GAIN, cap)
                        wiener = wiener * cap
                    endif
                    call merged_win%apply_filter(wiener)
                    if( merged_win%is_ft() ) call merged_win%ifft()
                    filtered_patch = merged_win%get_rmat()
                    accum(i0:i1,j0:j1,k0:k1) = accum(i0:i1,j0:j1,k0:k1) + filtered_patch * self%taper
                    norm(i0:i1,j0:j1,k0:k1) = norm(i0:i1,j0:j1,k0:k1) + self%taper**2
                    conf_accum(i0:i1,j0:j1,k0:k1) = conf_accum(i0:i1,j0:j1,k0:k1) + local_conf * self%taper**2
                    res_accum(i0:i1,j0:j1,k0:k1) = res_accum(i0:i1,j0:j1,k0:k1) + local_res * self%taper**2
                    env_accum(i0:i1,j0:j1,k0:k1) = env_accum(i0:i1,j0:j1,k0:k1) + env_conf * self%taper**2
                enddo
            enddo
        enddo
        allocate(out, conf, locres, env, mold=accum)
        where( norm > TINY )
            out    = accum / norm
            conf   = conf_accum / norm
            locres = res_accum / norm
            env    = min(1., env_accum / (norm * ENVELOPE_FULL_SUPPORT))
        elsewhere
            out    = 0.
            conf   = 0.
            locres = 0.
            env    = 0.
        end where
        out  = self%noise_mean + env * out
        conf = min(1., max(0., conf * env))
        where( env <= TINY ) locres = 0.
        call enhanced%new(self%ldim, self%smpd)
        call confidence%new(self%ldim, self%smpd)
        call local_resolution%new(self%ldim, self%smpd)
        call envelope%new(self%ldim, self%smpd)
        call enhanced%set_rmat(out, .false.)
        call confidence%set_rmat(conf, .false.)
        call local_resolution%set_rmat(locres, .false.)
        call envelope%set_rmat(env, .false.)
        call even_win%kill()
        call odd_win%kill()
        call merged_win%kill()
        deallocate(even_patch, odd_patch, merged_patch, filtered_patch)
        deallocate(accum, norm, conf_accum, res_accum, env_accum, out, conf, locres, env)
        deallocate(fsc, fsc_reg, wiener, lowpass, cap, res, xs, ys, zs)
    end subroutine enhance

    subroutine estimate_noise( self, even, odd )
        class(halfmap_enhancer), intent(inout) :: self
        real,                    intent(in)    :: even(:,:,:), odd(:,:,:)
        real, allocatable :: merged_outer(:), diff_outer(:)
        real :: med_diff, merged_sdev, diff_sdev
        integer :: i, j, k, n, edge
        edge = max(2, minval(self%ldim) / 16)
        n = count_outer(self%ldim, edge)
        allocate(merged_outer(n), diff_outer(n))
        n = 0
        do k = 1, self%ldim(3)
            do j = 1, self%ldim(2)
                do i = 1, self%ldim(1)
                    if( min(i-1,self%ldim(1)-i,j-1,self%ldim(2)-j,k-1,self%ldim(3)-k) < edge )then
                        n = n + 1
                        merged_outer(n) = 0.5 * (even(i,j,k) + odd(i,j,k))
                        diff_outer(n)   = 0.5 * (even(i,j,k) - odd(i,j,k))
                    endif
                enddo
            enddo
        enddo
        self%noise_mean = median(merged_outer)
        med_diff = median(diff_outer)
        merged_sdev = 1.4826 * median(abs(merged_outer - self%noise_mean))
        diff_sdev   = 1.4826 * median(abs(diff_outer - med_diff))
        ! The half-map difference measures independent noise.  Excess variance
        ! in the outer merged map estimates correlated processing background;
        ! keep the terms separate to avoid counting independent noise twice.
        self%noise_sdev = diff_sdev
        self%correlated_sdev = 0.
        if( merged_sdev > CORRELATED_NOISE_RATIO * diff_sdev ) &
            &self%correlated_sdev = sqrt(max(0., merged_sdev**2 - diff_sdev**2))
        if( self%noise_sdev <= sqrt(epsilon(1.)) ) &
            &self%noise_sdev = sqrt(sum((diff_outer - med_diff)**2) / real(max(1,n)))
        self%noise_sdev = max(self%noise_sdev, sqrt(epsilon(1.)))
        deallocate(merged_outer, diff_outer)
    end subroutine estimate_noise

    integer function count_outer( ldim, edge ) result( n )
        integer, intent(in) :: ldim(3), edge
        n = product(ldim) - product(max(0, ldim - 2 * edge))
    end function count_outer

    real function local_envelope_confidence( even, odd, noise_sdev, correlated_sdev, taper ) result( confidence )
        real, intent(in) :: even(:,:,:), odd(:,:,:), noise_sdev, correlated_sdev, taper(:,:,:)
        real(DP) :: signal_power, noise_power, norm
        norm = sum(real(taper,DP)**2)
        signal_power = sum(real(0.5 * (even + odd),DP)**2) / norm
        noise_power  = sum(real(0.5 * (even - odd),DP)**2) / norm
        noise_power  = max(noise_power, real(noise_sdev,DP)**2) + real(correlated_sdev,DP)**2
        signal_power = max(0._DP, signal_power - noise_power)
        confidence = real(signal_power / max(signal_power + noise_power, real(TINY,DP)))
        confidence = min(1., max(0., confidence))
    end function local_envelope_confidence

    subroutine regularize_fsc( raw, regularized, kcut )
        real,    intent(in)  :: raw(:)
        real,    intent(out) :: regularized(size(raw))
        integer, intent(out) :: kcut
        real :: work(size(raw))
        integer :: k
        work = min(0.99999, max(-1., raw))
        regularized = work
        do k = 2, size(raw) - 1
            regularized(k) = 0.25 * work(k-1) + 0.5 * work(k) + 0.25 * work(k+1)
        enddo
        do k = 2, size(raw)
            regularized(k) = min(regularized(k), regularized(k-1))
        enddo
        kcut = 0
        do k = 1, size(raw)
            if( regularized(k) < FSC_CUTOFF ) exit
            kcut = k
        enddo
        if( kcut < size(raw) ) regularized(kcut+1:) = min(regularized(kcut+1:), 0.)
    end subroutine regularize_fsc

    function window_starts( n, window, stride ) result( starts )
        integer, intent(in) :: n, window, stride
        integer, allocatable :: starts(:), work(:)
        integer :: current, next, count
        allocate(work(n), source=0)
        count = 1
        current = 1
        work(count) = current
        do while( current < n - window + 1 )
            next = min(current + stride, n - window + 1)
            if( next == current ) exit
            count = count + 1
            work(count) = next
            current = next
        enddo
        allocate(starts(count), source=work(1:count))
        deallocate(work)
    end function window_starts

    subroutine get_noise_stats( self, mean, sdev )
        class(halfmap_enhancer), intent(in) :: self
        real, intent(out) :: mean, sdev
        mean = self%noise_mean
        sdev = sqrt(self%noise_sdev**2 + self%correlated_sdev**2)
    end subroutine get_noise_stats

    integer function get_window( self ) result( window )
        class(halfmap_enhancer), intent(in) :: self
        window = self%window
    end function get_window

    subroutine kill( self )
        class(halfmap_enhancer), intent(inout) :: self
        if( allocated(self%taper) ) deallocate(self%taper)
        self%l_exists = .false.
        self%ldim = 0
        self%window = 0
        self%stride = 0
        self%smpd = 0.
        self%noise_mean = 0.
        self%noise_sdev = 0.
        self%correlated_sdev = 0.
    end subroutine kill

end module simple_halfmap_enhancer
