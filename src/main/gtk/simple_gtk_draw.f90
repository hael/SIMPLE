!@descr: stateless cairo and Pango helpers of SIMPLE's GTK code: colours, fonts, anchored text, text fitted to a width, marks, 8-bit images
module simple_gtk_draw
use, intrinsic :: iso_c_binding, only: c_ptr, c_int, c_int32_t, c_double, c_null_char, c_f_pointer, c_loc
use simple_defs,         only: dp
use simple_defs_gtk,     only: GTKUI_FONT, GTKUI_PANGO_UNITS, GTKUI_U_ELLIPSIS, GTKUI_COLOR_ACCENT, GTKUI_COLOR_WHITE
use simple_string_utils, only: real2str_sig, real2str_exp
use gtk,   only: CAIRO_FORMAT_RGB24, CAIRO_LINE_CAP_ROUND, CAIRO_LINE_JOIN_ROUND
use cairo, only: cairo_set_source_rgb, cairo_move_to, cairo_line_to, cairo_new_path, cairo_arc, cairo_fill, &
    cairo_stroke, cairo_set_line_width, cairo_set_line_cap, cairo_set_line_join, cairo_image_surface_create, &
    cairo_image_surface_get_data, cairo_image_surface_get_stride, cairo_surface_flush, cairo_surface_mark_dirty
use pango, only: pango_layout_set_text, pango_layout_get_pixel_size, pango_layout_set_font_description, &
    pango_font_description_from_string, pango_font_description_free, pango_font_description_set_absolute_size, &
    pango_cairo_show_layout
implicit none

public :: gtkui_font_new, gtkui_font_free, gtkui_set_rgb, gtkui_text_size, gtkui_draw_text, gtkui_fit_text
public :: gtkui_draw_check, gtkui_gray_surface, gtkui_surface_pixels
private

contains

    ! a Pango font description of GTKUI_FONT; free it with gtkui_font_free
    function gtkui_font_new() result( desc )
        type(c_ptr) :: desc
        desc = pango_font_description_from_string(GTKUI_FONT//c_null_char)
    end function gtkui_font_new

    subroutine gtkui_font_free( desc )
        type(c_ptr), intent(in) :: desc
        call pango_font_description_free(desc)
    end subroutine gtkui_font_free

    subroutine gtkui_set_rgb( cr, rgb )
        type(c_ptr), intent(in) :: cr
        real(dp),    intent(in) :: rgb(3)
        call cairo_set_source_rgb(cr, real(rgb(1), c_double), real(rgb(2), c_double), real(rgb(3), c_double))
    end subroutine gtkui_set_rgb

    ! sets the layout's font (desc at px pixels) and text, and returns its pixel size
    subroutine gtkui_text_size( layout, desc, px, text, tw, th )
        type(c_ptr),      intent(in)  :: layout, desc
        real(dp),         intent(in)  :: px
        character(len=*), intent(in)  :: text
        integer,          intent(out) :: tw, th
        integer(c_int), target :: cw, ch
        call pango_font_description_set_absolute_size(desc, real(px * GTKUI_PANGO_UNITS, c_double))
        call pango_layout_set_font_description(layout, desc)
        call pango_layout_set_text(layout, text//c_null_char, -1_c_int)
        call pango_layout_get_pixel_size(layout, c_loc(cw), c_loc(ch))
        tw = int(cw)
        th = int(ch)
    end subroutine gtkui_text_size

    ! text anchored at (x, y) as Tk anchors it: n, ne, e, se, s, sw, w, nw or c(entre)
    subroutine gtkui_draw_text( cr, layout, desc, px, text, x, y, anchor, rgb )
        type(c_ptr),      intent(in) :: cr, layout, desc
        real(dp),         intent(in) :: px, x, y, rgb(3)
        character(len=*), intent(in) :: text, anchor
        real(dp) :: x0, y0
        integer  :: tw, th
        call gtkui_text_size(layout, desc, px, text, tw, th)
        x0 = x
        y0 = y
        if( index(anchor, 'e') > 0 ) x0 = x - real(tw, dp)
        if( anchor == 'n' .or. anchor == 's' .or. anchor == 'c' ) x0 = x - 0.5_dp * real(tw, dp)
        if( index(anchor, 's') > 0 ) y0 = y - real(th, dp)
        if( anchor == 'e' .or. anchor == 'w' .or. anchor == 'c' ) y0 = y - 0.5_dp * real(th, dp)
        call gtkui_set_rgb(cr, rgb)
        call cairo_move_to(cr, anint(x0), anint(y0))
        call pango_cairo_show_layout(cr, layout)
    end subroutine gtkui_draw_text

    ! text if it fits avail pixels; else, for a number, 3, 2 or 1 significant digits in fixed or exponent form;
    ! else the longest prefix with an ellipsis, or nothing
    function gtkui_fit_text( layout, desc, px, text, avail ) result( txt )
        type(c_ptr),      intent(in) :: layout, desc
        real(dp),         intent(in) :: px
        character(len=*), intent(in) :: text
        integer,          intent(in) :: avail
        character(len=:), allocatable :: txt, abbr, digits
        real(dp) :: num
        integer  :: tw, th, p, k, ios, lo, hi, mid
        call gtkui_text_size(layout, desc, px, text, tw, th)
        if( tw <= avail )then
            txt = text
            return
        endif
        abbr   = text
        digits = ''
        do k = 1,len(text)
            if( text(k:k) /= ',' ) digits = digits//text(k:k)
        end do
        read(digits, *, iostat=ios) num
        if( ios == 0 )then
            do p = 3,1,-1
                do k = 1,2
                    if( k == 1 )then
                        abbr = compact_exponent(real2str_sig(num, p))
                    else
                        abbr = compact_exponent(real2str_exp(num, p - 1))
                    endif
                    call gtkui_text_size(layout, desc, px, abbr, tw, th)
                    if( tw <= avail )then
                        txt = abbr
                        return
                    endif
                end do
            end do
        endif
        call gtkui_text_size(layout, desc, px, GTKUI_U_ELLIPSIS, tw, th)
        if( tw > avail )then
            txt = ''
            return
        endif
        lo = 0
        hi = len(abbr)
        do while( lo < hi )
            mid = (lo + hi + 1) / 2
            call gtkui_text_size(layout, desc, px, abbr(1:mid)//GTKUI_U_ELLIPSIS, tw, th)
            if( tw <= avail )then
                lo = mid
            else
                hi = mid - 1
            endif
        end do
        txt = trim(abbr(1:lo))//GTKUI_U_ELLIPSIS
    end function gtkui_fit_text

    ! the selection mark: a check in an accent disc in the top-right corner of a tile at (x, y)
    subroutine gtkui_draw_check( cr, x, y, tile )
        type(c_ptr), intent(in) :: cr
        real(dp),    intent(in) :: x, y
        integer,     intent(in) :: tile
        real(dp) :: d, right, top
        d     = merge(18._dp, 14._dp, tile >= 90)
        right = x + real(tile - 8, dp)
        top   = y + 8._dp
        call cairo_new_path(cr)
        call cairo_arc(cr, right - 0.5_dp * d, top + 0.5_dp * d, 0.5_dp * d, 0._dp, 2._dp * acos(-1._dp))
        call gtkui_set_rgb(cr, GTKUI_COLOR_ACCENT)
        call cairo_fill(cr)
        call cairo_move_to(cr, right - d * 0.76_dp, top + d * 0.50_dp)
        call cairo_line_to(cr, right - d * 0.55_dp, top + d * 0.71_dp)
        call cairo_line_to(cr, right - d * 0.23_dp, top + d * 0.32_dp)
        call cairo_set_line_width(cr, 2._dp)
        call cairo_set_line_cap(cr, CAIRO_LINE_CAP_ROUND)
        call cairo_set_line_join(cr, CAIRO_LINE_JOIN_ROUND)
        call gtkui_set_rgb(cr, GTKUI_COLOR_WHITE)
        call cairo_stroke(cr)
    end subroutine gtkui_draw_check

    ! an 8-bit image (0..255, column = x) as a new RGB24 cairo surface with R = G = B; the caller destroys it
    function gtkui_gray_surface( gray ) result( surf )
        integer, intent(in) :: gray(:,:)
        type(c_ptr) :: surf
        integer(c_int32_t), pointer :: pix(:,:)
        integer :: j, nx, ny
        nx   = size(gray,1)
        ny   = size(gray,2)
        surf = cairo_image_surface_create(CAIRO_FORMAT_RGB24, int(nx, c_int), int(ny, c_int))
        call gtkui_surface_pixels(surf, ny, pix)
        do j = 1,ny
            pix(1:nx,j) = int(min(max(gray(:,j), 0), 255) * 65793, c_int32_t)
        end do
        call cairo_surface_mark_dirty(surf)
    end function gtkui_gray_surface

    ! the pixels of an RGB24 image surface of nrows rows: rows stride bytes apart, a pixel a native-endian
    ! 0x00RRGGBB word; call cairo_surface_mark_dirty after writing them
    subroutine gtkui_surface_pixels( surf, nrows, pix )
        type(c_ptr),                 intent(in)  :: surf
        integer,                     intent(in)  :: nrows
        integer(c_int32_t), pointer, intent(out) :: pix(:,:)
        integer :: stride
        call cairo_surface_flush(surf)
        stride = int(cairo_image_surface_get_stride(surf))
        call c_f_pointer(cairo_image_surface_get_data(surf), pix, [stride / 4, nrows])
    end subroutine gtkui_surface_pixels

    ! ---- private -----------------------------------------------------------

    ! e+05 -> e5, e-05 -> e-5, as the script abbreviated numbers
    function compact_exponent( str_in ) result( str )
        character(len=*), intent(in) :: str_in
        character(len=:), allocatable :: str, digits, sgn
        integer :: k, i
        str = str_in
        k = index(str, 'e')
        if( k == 0 ) return
        digits = str(k+1:)
        sgn    = ''
        if( len(digits) > 0 )then
            if( digits(1:1) == '-' ) sgn = '-'
            if( digits(1:1) == '-' .or. digits(1:1) == '+' ) digits = digits(2:)
        endif
        i = 1
        do while( i < len(digits) )
            if( digits(i:i) /= '0' ) exit
            i = i + 1
        end do
        str = str(1:k)//sgn//digits(i:)
    end function compact_exponent

end module simple_gtk_draw
