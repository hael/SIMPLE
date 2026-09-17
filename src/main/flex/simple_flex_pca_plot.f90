!@descr: flex_pca latent figures rendered in-engine (no external plotting): a three-panel JPEG with
!! the log-density UMAP core, the UMAP coloured by label, and the latent z1 vs z2 plane coloured by
!! label. Axis limits are the 0.5-99.5 percentiles of the points, so outliers never set the scale.
module simple_flex_pca_plot
use simple_core_module_api
use simple_jpg, only: write_rgb_jpeg
implicit none
private
#include "simple_local_flags.inc"
public :: flex_plot_latent_jpg

integer, parameter :: PANEL   = 640     !< panel edge in pixels
integer, parameter :: GUTTER  = 24      !< white margin around panels
integer, parameter :: NBINS   = 160     !< density histogram bins per axis
integer, parameter :: DOT     = 2       !< scatter dot edge in pixels
integer, parameter :: LEGEND  = 14      !< legend square edge
real,    parameter :: PCT     = 0.5     !< percentile clipped at each end
!> matplotlib tab10
real, parameter :: TAB10(3,10) = reshape([ &
    &0.122,0.467,0.706,  1.000,0.498,0.055,  0.173,0.627,0.173,  0.839,0.153,0.157,  0.580,0.404,0.741, &
    &0.549,0.337,0.294,  0.890,0.467,0.761,  0.498,0.498,0.498,  0.737,0.741,0.133,  0.090,0.745,0.812 ], [3,10])
real, parameter :: GREY_NONE(3) = [0.83, 0.83, 0.83]

contains

    !> fname: output .jpg; (x,y): UMAP coordinates; (z1,z2): latent plane; labels: 0 = unlabelled,
    !! 1..nlab coloured. All arrays are over the same n points. Points with a non-finite coordinate
    !! are skipped in that panel.
    subroutine flex_plot_latent_jpg( fname, x, y, z1, z2, n, labels, nlab )
        character(len=*),  intent(in) :: fname
        real,              intent(in) :: x(:), y(:), z1(:), z2(:)
        integer,           intent(in) :: n
        integer, optional, intent(in) :: labels(:)
        integer, optional, intent(in) :: nlab
        real,    allocatable :: rgb(:,:,:)
        integer, allocatable :: lab(:)
        integer :: W, H, nl, x0
        if( n < 10 ) return
        allocate(lab(n), source=0)
        nl = 0
        if( present(labels) )then
            lab(1:n) = labels(1:n)
            nl = maxval(lab); if( present(nlab) ) nl = nlab
        endif
        W = 3*PANEL + 4*GUTTER; H = PANEL + 2*GUTTER
        allocate(rgb(3,W,H), source=1.0)
        x0 = GUTTER
        call density_panel(rgb, x0, GUTTER, x, y, n)
        x0 = x0 + PANEL + GUTTER
        call scatter_panel(rgb, x0, GUTTER, x, y, n, lab, nl)
        x0 = x0 + PANEL + GUTTER
        call scatter_panel(rgb, x0, GUTTER, z1, z2, n, lab, nl)
        call write_rgb_jpeg(fname, rgb, quality=92)
        deallocate(rgb, lab)
    end subroutine flex_plot_latent_jpg

    ! ---- panels ----------------------------------------------------------------------------------

    subroutine density_panel( rgb, px, py, a, b, n )
        real,    intent(inout) :: rgb(:,:,:)
        integer, intent(in)    :: px, py, n
        real,    intent(in)    :: a(:), b(:)
        real    :: alo, ahi, blo, bhi, cnt(NBINS,NBINS), t, g, cmax
        integer :: i, ia, ib, bi, bj, i0, j0, i1, j1, bpx
        call limits(a, n, alo, ahi); call limits(b, n, blo, bhi)
        cnt = 0.0
        do i = 1, n
            if( .not. inside(a(i), b(i), alo, ahi, blo, bhi) ) cycle
            ia = 1 + int((a(i)-alo)/(ahi-alo)*real(NBINS)); ia = max(1, min(NBINS, ia))
            ib = 1 + int((b(i)-blo)/(bhi-blo)*real(NBINS)); ib = max(1, min(NBINS, ib))
            cnt(ia,ib) = cnt(ia,ib) + 1.0
        end do
        cmax = log(1.0 + maxval(cnt)); if( cmax <= 0.0 ) cmax = 1.0
        bpx  = PANEL / NBINS
        do bj = 1, NBINS
            do bi = 1, NBINS
                if( cnt(bi,bj) <= 0.0 ) cycle
                t = log(1.0 + cnt(bi,bj)) / cmax
                g = 1.0 - 0.92*t
                i0 = px + (bi-1)*bpx; i1 = i0 + bpx - 1
                j0 = py + PANEL - bj*bpx; j1 = j0 + bpx - 1     ! b increases upwards
                call fill(rgb, i0, i1, j0, j1, [g,g,g])
            end do
        end do
        call frame(rgb, px, py)
    end subroutine density_panel

    subroutine scatter_panel( rgb, px, py, a, b, n, lab, nl )
        real,    intent(inout) :: rgb(:,:,:)
        integer, intent(in)    :: px, py, n, nl
        real,    intent(in)    :: a(:), b(:)
        integer, intent(in)    :: lab(:)
        real    :: alo, ahi, blo, bhi, col(3)
        integer :: pass, k, i, ii, stride, ix, iy, l
        call limits(a, n, alo, ahi); call limits(b, n, blo, bhi)
        ! unlabelled points first (grey), then the labels in an interleaved order so no class hides
        ! another by draw order; stride through the points with a step coprime to n
        stride = max(1, int(sqrt(real(n))) ); do while( gcd(stride, n) /= 1 ); stride = stride + 1; end do
        do pass = 0, 1
            ii = 0
            do k = 1, n
                ii = mod(ii + stride, n); i = ii + 1
                l = lab(i)
                if( (pass == 0 .and. l /= 0) .or. (pass == 1 .and. l == 0) ) cycle
                if( .not. inside(a(i), b(i), alo, ahi, blo, bhi) ) cycle
                ix = px + int((a(i)-alo)/(ahi-alo)*real(PANEL-DOT))
                iy = py + PANEL - DOT - int((b(i)-blo)/(bhi-blo)*real(PANEL-DOT))
                if( l == 0 )then
                    col = GREY_NONE
                else
                    col = TAB10(:, 1 + mod(l-1, 10))
                endif
                call fill(rgb, ix, ix+DOT-1, iy, iy+DOT-1, col)
            end do
        end do
        ! legend: one square per label, top-left, in label order
        do l = 1, min(nl, 20)
            ix = px + 8; iy = py + 8 + (l-1)*(LEGEND+4)
            call fill(rgb, ix, ix+LEGEND-1, iy, iy+LEGEND-1, TAB10(:, 1 + mod(l-1, 10)))
            call fill(rgb, ix-1, ix+LEGEND, iy-1, iy-1, [0.2,0.2,0.2]); call fill(rgb, ix-1, ix+LEGEND, iy+LEGEND, iy+LEGEND, [0.2,0.2,0.2])
        end do
        call frame(rgb, px, py)
    end subroutine scatter_panel

    ! ---- helpers ---------------------------------------------------------------------------------

    pure logical function inside( u, v, ulo, uhi, vlo, vhi )
        real, intent(in) :: u, v, ulo, uhi, vlo, vhi
        inside = .false.
        if( u /= u .or. v /= v ) return
        inside = u >= ulo .and. u <= uhi .and. v >= vlo .and. v <= vhi
    end function inside

    subroutine fill( rgb, i0, i1, j0, j1, col )
        real,    intent(inout) :: rgb(:,:,:)
        integer, intent(in)    :: i0, i1, j0, j1
        real,    intent(in)    :: col(3)
        integer :: i, j
        do j = max(1,j0), min(size(rgb,3),j1)
            do i = max(1,i0), min(size(rgb,2),i1)
                rgb(:,i,j) = col
            end do
        end do
    end subroutine fill

    subroutine frame( rgb, px, py )
        real,    intent(inout) :: rgb(:,:,:)
        integer, intent(in)    :: px, py
        real, parameter :: dark(3) = [0.15,0.15,0.15]
        call fill(rgb, px-1, px+PANEL, py-1, py-1, dark); call fill(rgb, px-1, px+PANEL, py+PANEL, py+PANEL, dark)
        call fill(rgb, px-1, px-1, py-1, py+PANEL, dark); call fill(rgb, px+PANEL, px+PANEL, py-1, py+PANEL, dark)
    end subroutine frame

    !> [PCT, 100-PCT] percentile limits of the finite values, widened slightly; degenerate -> +-1
    subroutine limits( v, n, lo, hi )
        real,    intent(in)  :: v(:)
        integer, intent(in)  :: n
        real,    intent(out) :: lo, hi
        real, allocatable :: w(:)
        integer :: m, i, klo, khi
        real    :: pad
        allocate(w(n)); m = 0
        do i = 1, n
            if( v(i) == v(i) )then; m = m + 1; w(m) = v(i); endif
        end do
        if( m < 2 )then; lo = -1.0; hi = 1.0; return; endif
        klo = max(1, nint(real(m)*PCT/100.0)); khi = min(m, nint(real(m)*(100.0-PCT)/100.0))
        lo = kth_smallest(w(1:m), klo); hi = kth_smallest(w(1:m), khi)
        if( hi <= lo )then; lo = lo - 1.0; hi = hi + 1.0; endif
        pad = 0.02*(hi-lo); lo = lo - pad; hi = hi + pad
        deallocate(w)
    end subroutine limits

    !> quickselect on a copy (Hoare partition), O(n) expected
    function kth_smallest( arr, k ) result( val )
        real,    intent(in) :: arr(:)
        integer, intent(in) :: k
        real :: val
        real, allocatable :: a(:)
        real    :: piv, tmp
        integer :: lo, hi, i, j, m
        allocate(a, source=arr); lo = 1; hi = size(a)
        do while( lo < hi )
            m = (lo + hi)/2; piv = a(m); i = lo; j = hi
            do
                do while( a(i) < piv ); i = i + 1; end do
                do while( a(j) > piv ); j = j - 1; end do
                if( i <= j )then
                    tmp = a(i); a(i) = a(j); a(j) = tmp; i = i + 1; j = j - 1
                endif
                if( i > j ) exit
            end do
            if( k <= j )then
                hi = j
            else if( k >= i )then
                lo = i
            else
                exit
            endif
        end do
        val = a(k); deallocate(a)
    end function kth_smallest

    pure integer function gcd( a, b )
        integer, intent(in) :: a, b
        integer :: x, y, t
        x = abs(a); y = abs(b)
        do while( y /= 0 ); t = mod(x, y); x = y; y = t; end do
        gcd = max(x, 1)
    end function gcd

end module simple_flex_pca_plot
