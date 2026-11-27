!@descr: display state of a stack of 2D images: intensity window, brightness and contrast, histogram and 8-bit tiles
! The numbers of select2d.py (doc/implementation_notes/planned/select2D.md, decision 9; limits in simple_defs_gtk):
! stack range over finite pixels, an automatic window from the medians of the visible images' percentiles, a
! histogram, adjustment in double precision, resampling as Pillow's (Lanczos-3 when shrinking, bilinear when
! enlarging) before the one conversion to 8 bits. Used by GTK views; GTK-free itself.
module simple_stack_display
use, intrinsic :: ieee_arithmetic, only: ieee_is_finite, ieee_rint
use simple_defs,     only: dp
use simple_defs_gtk, only: GTKUI_BRIGHTNESS_LIM, GTKUI_CONTRAST_LIM, GTKUI_AUTO_PCT, GTKUI_HIST_NBINS
use simple_error,    only: simple_exception
use simple_image,    only: image
use simple_stat,     only: median, percentile
implicit none

public :: stack_display
private
#include "simple_local_flags.inc"

real(dp), parameter :: LANCZOS_SUPPORT    = 3._dp
real(dp), parameter :: BILINEAR_SUPPORT   = 1._dp
real(dp), parameter :: PI_DP              = acos(-1._dp)

type :: stack_display
    private
    real,            allocatable :: pix(:,:,:)        ! nx x ny x ncls raw intensities
    integer(kind=8), allocatable :: counts(:)         ! histogram of the finite pixels of every class
    integer  :: nx = 0, ny = 0, ncls = 0
    integer  :: revision = 0                          ! changes with every display setting
    real(dp) :: range_lo = 0._dp, range_hi = 1._dp    ! manual bounds and histogram domain
    real(dp) :: auto_lo  = 0._dp, auto_hi  = 1._dp    ! automatic window
    real(dp) :: brightness = 0._dp, contrast = 1._dp
    real(dp) :: wmin = 0._dp, wmax = 1._dp            ! current window
    logical  :: usable   = .false.                    ! the stack has finite pixels
    logical  :: l_exists = .false.
  contains
    procedure :: new
    procedure :: kill
    procedure :: get_ncls
    procedure :: get_box
    procedure :: get_revision
    procedure :: is_usable
    procedure :: get_brightness
    procedure :: get_contrast
    procedure :: get_window
    procedure :: get_range
    procedure :: get_counts
    procedure :: set_brightness
    procedure :: set_contrast
    procedure :: set_window_limit
    procedure :: reset
    procedure :: mapping_box
    procedure :: render
end type stack_display

contains

    ! the raw images; visible marks the images a view shows (they set the automatic window)
    subroutine new( self, imgs, visible, auto_percentiles )
        class(stack_display), intent(inout) :: self
        class(image),        intent(in)    :: imgs(:)
        logical,             intent(in)    :: visible(:)
        real(dp), optional,  intent(in)    :: auto_percentiles(2)
        real, allocatable :: rmat(:,:,:)
        integer :: icls, ldim(3)
        call self%kill
        self%ncls = size(imgs)
        if( self%ncls == 0 ) THROW_HARD('no images; stack_display%new')
        if( size(visible) /= self%ncls ) THROW_HARD('visible mask does not match the images; stack_display%new')
        ldim    = imgs(1)%get_ldim()
        self%nx = ldim(1)
        self%ny = ldim(2)
        allocate(self%pix(self%nx,self%ny,self%ncls))
        do icls = 1,self%ncls
            rmat = imgs(icls)%get_rmat()
            if( size(rmat,1) /= self%nx .or. size(rmat,2) /= self%ny ) THROW_HARD('images differ in size; stack_display%new')
            self%pix(:,:,icls) = rmat(:,:,1)
        end do
        call stack_range(self)
        call auto_window(self, visible, auto_percentiles)
        call histogram(self)
        self%brightness = 0._dp
        self%contrast   = 1._dp
        self%wmin       = self%auto_lo
        self%wmax       = self%auto_hi
        self%revision   = 1
        self%l_exists   = .true.
    end subroutine new

    subroutine kill( self )
        class(stack_display), intent(inout) :: self
        if( allocated(self%pix)    ) deallocate(self%pix)
        if( allocated(self%counts) ) deallocate(self%counts)
        self%nx = 0
        self%ny = 0
        self%ncls = 0
        self%l_exists = .false.
    end subroutine kill

    integer function get_ncls( self )
        class(stack_display), intent(in) :: self
        get_ncls = self%ncls
    end function get_ncls

    function get_box( self ) result( box )
        class(stack_display), intent(in) :: self
        integer :: box(2)
        box = [self%nx, self%ny]
    end function get_box

    integer function get_revision( self )
        class(stack_display), intent(in) :: self
        get_revision = self%revision
    end function get_revision

    logical function is_usable( self )
        class(stack_display), intent(in) :: self
        is_usable = self%usable
    end function is_usable

    real(dp) function get_brightness( self )
        class(stack_display), intent(in) :: self
        get_brightness = self%brightness
    end function get_brightness

    real(dp) function get_contrast( self )
        class(stack_display), intent(in) :: self
        get_contrast = self%contrast
    end function get_contrast

    function get_window( self ) result( lims )
        class(stack_display), intent(in) :: self
        real(dp) :: lims(2)
        lims = [self%wmin, self%wmax]
    end function get_window

    ! the histogram's domain, which is also the range Min and Max can take
    function get_range( self ) result( lims )
        class(stack_display), intent(in) :: self
        real(dp) :: lims(2)
        lims = [self%range_lo, self%range_hi]
    end function get_range

    function get_counts( self ) result( counts )
        class(stack_display), intent(in) :: self
        integer(kind=8) :: counts(GTKUI_HIST_NBINS)
        counts = self%counts
    end function get_counts

    ! brightness in percent, clamped to [-100, 100] and rounded to an integer; true when it changed
    logical function set_brightness( self, val ) result( changed )
        class(stack_display), intent(inout) :: self
        real(dp),            intent(in)    :: val
        real(dp) :: v
        changed = .false.
        if( .not. ieee_is_finite(val) ) return
        v = anint(min(max(val, GTKUI_BRIGHTNESS_LIM(1)), GTKUI_BRIGHTNESS_LIM(2)))
        changed = v /= self%brightness
        if( changed )then
            self%brightness = v
            self%revision   = self%revision + 1
        endif
    end function set_brightness

    ! contrast factor, clamped to [0.1, 5] and rounded to two decimals; true when it changed
    logical function set_contrast( self, val ) result( changed )
        class(stack_display), intent(inout) :: self
        real(dp),            intent(in)    :: val
        real(dp) :: v
        changed = .false.
        if( .not. ieee_is_finite(val) ) return
        v = anint(min(max(val, GTKUI_CONTRAST_LIM(1)), GTKUI_CONTRAST_LIM(2)) * 100._dp) / 100._dp
        changed = v /= self%contrast
        if( changed )then
            self%contrast = v
            self%revision = self%revision + 1
        endif
    end function set_contrast

    ! Min (which = 1) or Max (which = 2) in raw units, clamped to the range and short of the other limit
    logical function set_window_limit( self, which, val ) result( changed )
        class(stack_display), intent(inout) :: self
        integer,             intent(in)    :: which
        real(dp),            intent(in)    :: val
        real(dp) :: v
        changed = .false.
        if( .not. self%usable ) return
        if( .not. ieee_is_finite(val) ) return
        if( which == 1 )then
            v = min(max(val, self%range_lo), nearest(self%wmax, -1._dp))
            changed = v /= self%wmin
            if( changed ) self%wmin = v
        else
            v = min(max(val, nearest(self%wmin, 1._dp)), self%range_hi)
            changed = v /= self%wmax
            if( changed ) self%wmax = v
        endif
        if( changed ) self%revision = self%revision + 1
    end function set_window_limit

    ! brightness 0, contrast 1 and the automatic window
    subroutine reset( self )
        class(stack_display), intent(inout) :: self
        self%brightness = 0._dp
        self%contrast   = 1._dp
        self%wmin       = self%auto_lo
        self%wmax       = self%auto_hi
        self%revision   = self%revision + 1
    end subroutine reset

    ! the raw values that map to black and to white under the current settings
    subroutine mapping_box( self, black, white )
        class(stack_display), intent(in)  :: self
        real(dp),            intent(out) :: black, white
        real(dp) :: mid, span
        span  = self%wmax - self%wmin
        mid   = 0.5_dp * (self%wmin + self%wmax)
        black = mid + (0._dp - 0.5_dp - self%brightness / 100._dp) * span / self%contrast
        white = mid + (1._dp - 0.5_dp - self%brightness / 100._dp) * span / self%contrast
    end subroutine mapping_box

    ! image icls as an 8-bit tile of size(gray,1) x size(gray,2): adjusted and windowed at full
    ! resolution, resampled, clipped, then rounded; non-finite pixels are black
    subroutine render( self, icls, gray )
        class(stack_display), intent(in)  :: self
        integer,             intent(in)  :: icls
        integer,             intent(out) :: gray(:,:)
        real(dp), allocatable :: win(:,:), out(:,:)
        real(dp) :: mid, span, shift, v
        integer  :: i, j, ow, oh
        logical  :: l_lanczos
        if( icls < 1 .or. icls > self%ncls ) THROW_HARD('image index out of range; stack_display%render')
        ow = size(gray,1)
        oh = size(gray,2)
        span  = self%wmax - self%wmin
        mid   = 0.5_dp * (self%wmin + self%wmax)
        shift = self%brightness / 100._dp * span
        allocate(win(self%nx,self%ny))
        do j = 1,self%ny
            do i = 1,self%nx
                if( ieee_is_finite(self%pix(i,j,icls)) )then
                    v = mid + self%contrast * (real(self%pix(i,j,icls), dp) - mid) + shift
                    win(i,j) = min(max((v - self%wmin) / span, 0._dp), 1._dp)
                else
                    win(i,j) = 0._dp
                endif
            end do
        end do
        l_lanczos = min(ow, oh) < min(self%nx, self%ny)
        call resample(win, ow, oh, l_lanczos, out)
        do j = 1,oh
            do i = 1,ow
                gray(i,j) = int(ieee_rint(min(max(out(i,j), 0._dp), 1._dp) * 255._dp))
            end do
        end do
    end subroutine render

    ! ---- private -----------------------------------------------------------

    ! finite range of the stack; a constant stack gets a small nonzero window
    subroutine stack_range( self )
        type(stack_display), intent(inout) :: self
        real(dp) :: lo, hi, v
        integer  :: i, j, icls
        logical  :: l_any
        l_any = .false.
        lo    = 0._dp
        hi    = 0._dp
        do icls = 1,self%ncls
            do j = 1,self%ny
                do i = 1,self%nx
                    if( .not. ieee_is_finite(self%pix(i,j,icls)) ) cycle
                    v = real(self%pix(i,j,icls), dp)
                    if( l_any )then
                        lo = min(lo, v)
                        hi = max(hi, v)
                    else
                        lo    = v
                        hi    = v
                        l_any = .true.
                    endif
                end do
            end do
        end do
        self%usable = l_any
        if( .not. l_any )then
            self%range_lo = 0._dp
            self%range_hi = 1._dp
        else
            if( hi == lo ) hi = lo + max(1._dp, abs(lo) * 1.e-6_dp)
            self%range_lo = lo
            self%range_hi = hi
        endif
    end subroutine stack_range

    ! medians, taken separately, of the visible images' low and high percentiles of their finite pixels
    subroutine auto_window( self, visible, auto_percentiles )
        type(stack_display), intent(inout) :: self
        logical,            intent(in)    :: visible(:)
        real(dp), optional, intent(in)    :: auto_percentiles(2)
        real(dp), allocatable :: vals(:), lows(:), highs(:)
        real(dp) :: lo, hi, pct(2)
        integer  :: icls, n, nlim
        allocate(lows(self%ncls), highs(self%ncls), vals(self%nx * self%ny))
        pct = GTKUI_AUTO_PCT
        if( present(auto_percentiles) ) pct = auto_percentiles
        if( .not. all(ieee_is_finite(pct)) ) THROW_HARD('nonfinite display percentiles')
        if( pct(1) < 0._dp .or. pct(2) > 100._dp .or. pct(1) >= pct(2) ) THROW_HARD('invalid display percentiles')
        nlim = 0
        do icls = 1,self%ncls
            if( .not. visible(icls) ) cycle
            n = count(ieee_is_finite(self%pix(:,:,icls)))
            if( n == 0 ) cycle
            vals(1:n) = real(pack(self%pix(:,:,icls), ieee_is_finite(self%pix(:,:,icls))), dp)
            lo = percentile(vals(1:n), pct(1))
            hi = percentile(vals(1:n), pct(2))
            if( lo < hi )then
                nlim        = nlim + 1
                lows(nlim)  = lo
                highs(nlim) = hi
            endif
        end do
        self%auto_lo = self%range_lo
        self%auto_hi = self%range_hi
        if( nlim > 0 )then
            lo = median(lows(1:nlim))
            hi = median(highs(1:nlim))
            if( lo < hi )then
                self%auto_lo = lo
                self%auto_hi = hi
            endif
        endif
    end subroutine auto_window

    ! equal bins over the stack range, as numpy.histogram (the last bin closed)
    subroutine histogram( self )
        type(stack_display), intent(inout) :: self
        real(dp) :: edges(0:GTKUI_HIST_NBINS), v, norm
        integer  :: i, j, icls, k
        allocate(self%counts(GTKUI_HIST_NBINS), source=0_8)
        if( .not. self%usable ) return
        do k = 0,GTKUI_HIST_NBINS
            edges(k) = self%range_lo + (self%range_hi - self%range_lo) * real(k, dp) / real(GTKUI_HIST_NBINS, dp)
        end do
        edges(GTKUI_HIST_NBINS) = self%range_hi
        norm = real(GTKUI_HIST_NBINS, dp) / (self%range_hi - self%range_lo)
        do icls = 1,self%ncls
            do j = 1,self%ny
                do i = 1,self%nx
                    if( .not. ieee_is_finite(self%pix(i,j,icls)) ) cycle
                    v = real(self%pix(i,j,icls), dp)
                    if( v < edges(0) .or. v > edges(GTKUI_HIST_NBINS) ) cycle
                    k = min(int((v - edges(0)) * norm), GTKUI_HIST_NBINS - 1)
                    if( v < edges(k) ) k = max(k - 1, 0)
                    if( k < GTKUI_HIST_NBINS - 1 )then
                        if( v >= edges(k + 1) ) k = k + 1
                    endif
                    self%counts(k + 1) = self%counts(k + 1) + 1_8
                end do
            end do
        end do
    end subroutine histogram

    ! separable resampling of img to ow x oh as Pillow's ImagingResample: horizontal, then vertical
    subroutine resample( img, ow, oh, l_lanczos, out )
        real(dp),              intent(in)  :: img(:,:)
        integer,               intent(in)  :: ow, oh
        logical,               intent(in)  :: l_lanczos
        real(dp), allocatable, intent(out) :: out(:,:)
        real(dp), allocatable :: tmp(:,:), wx(:,:), wy(:,:)
        integer,  allocatable :: x0(:), nxk(:), y0(:), nyk(:)
        integer :: nx, ny, i, j, k
        nx = size(img,1)
        ny = size(img,2)
        call resample_coeffs(nx, ow, l_lanczos, x0, nxk, wx)
        call resample_coeffs(ny, oh, l_lanczos, y0, nyk, wy)
        allocate(tmp(ow,ny), source=0._dp)
        do j = 1,ny
            do i = 1,ow
                do k = 1,nxk(i)
                    tmp(i,j) = tmp(i,j) + img(x0(i) + k, j) * wx(k,i)
                end do
            end do
        end do
        allocate(out(ow,oh), source=0._dp)
        do j = 1,oh
            do k = 1,nyk(j)
                out(:,j) = out(:,j) + tmp(:, y0(j) + k) * wy(k,j)
            end do
        end do
    end subroutine resample

    ! Pillow's precompute_coeffs: output pixel i reads inputs x0(i)+1 .. x0(i)+n(i) with weights w(:,i)
    subroutine resample_coeffs( insize, outsize, l_lanczos, x0, n, w )
        integer,               intent(in)  :: insize, outsize
        logical,               intent(in)  :: l_lanczos
        integer,  allocatable, intent(out) :: x0(:), n(:)
        real(dp), allocatable, intent(out) :: w(:,:)
        real(dp) :: scale, filterscale, support, center, ss, ww
        integer  :: ksize, i, k, xmin, xmax
        scale       = real(insize, dp) / real(outsize, dp)
        filterscale = max(scale, 1._dp)
        if( l_lanczos )then
            support = LANCZOS_SUPPORT * filterscale
        else
            support = BILINEAR_SUPPORT * filterscale
        endif
        ksize = 2 * ceiling(support) + 1
        allocate(x0(outsize), n(outsize), w(ksize,outsize))
        w  = 0._dp
        ss = 1._dp / filterscale
        do i = 1,outsize
            center = (real(i - 1, dp) + 0.5_dp) * scale
            xmin   = max(int(center - support + 0.5_dp), 0)
            xmax   = min(int(center + support + 0.5_dp), insize) - xmin
            xmax   = min(xmax, ksize)
            ww     = 0._dp
            do k = 1,xmax
                w(k,i) = filter_weight((real(k - 1 + xmin, dp) - center + 0.5_dp) * ss, l_lanczos)
                ww     = ww + w(k,i)
            end do
            if( ww /= 0._dp ) w(1:xmax,i) = w(1:xmax,i) / ww
            x0(i) = xmin
            n(i)  = max(xmax, 0)
        end do
    end subroutine resample_coeffs

    real(dp) function filter_weight( x, l_lanczos )
        real(dp), intent(in) :: x
        logical,  intent(in) :: l_lanczos
        if( l_lanczos )then
            if( x >= -LANCZOS_SUPPORT .and. x < LANCZOS_SUPPORT )then
                filter_weight = sinc(x) * sinc(x / LANCZOS_SUPPORT)
            else
                filter_weight = 0._dp
            endif
        else
            filter_weight = max(1._dp - abs(x), 0._dp)
        endif
    end function filter_weight

    real(dp) function sinc( x )
        real(dp), intent(in) :: x
        if( x == 0._dp )then
            sinc = 1._dp
        else
            sinc = sin(PI_DP * x) / (PI_DP * x)
        endif
    end function sinc

end module simple_stack_display
