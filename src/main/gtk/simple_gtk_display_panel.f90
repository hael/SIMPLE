!@descr: a GTK 3 display panel for an image stack: brightness, contrast, reset and the histogram with draggable Min and Max
! Bound to a stack_display (limits in simple_defs_gtk): the sliders and number fields share one adjustment each,
! the histogram shows the raw distribution, the Min/Max lines and the raw range mapped from black to white.
! Every change asks the notify widget (the view of the images) to redraw.
module simple_gtk_display_panel
use, intrinsic :: iso_c_binding, only: c_ptr, c_int, c_double, c_float, c_null_char, c_null_ptr, c_funloc, c_loc, &
    c_f_pointer, c_associated
use simple_defs,          only: dp
use simple_defs_gtk,      only: GTKUI_BRIGHTNESS_LIM, GTKUI_CONTRAST_LIM, GTKUI_HIST_NBINS, GTKUI_HIST_W, GTKUI_HIST_H, &
    GTKUI_HIT_RADIUS, GTKUI_SLIDER_W, GTKUI_SMALL_PX, GTKUI_MESSAGE_PX, GTKUI_U_MIDDOT, GTKUI_U_TIMES, GTKUI_U_UP, &
    GTKUI_U_DOWN, GTKUI_U_LEFT, GTKUI_U_RIGHT, GTKUI_KEY_LEFT, GTKUI_KEY_UP, GTKUI_KEY_RIGHT, GTKUI_KEY_DOWN, &
    GTKUI_COLOR_FG, GTKUI_COLOR_BARS, GTKUI_COLOR_MIN, GTKUI_COLOR_MAX, GTKUI_COLOR_RANGE
use simple_stack_display, only: stack_display
use simple_string_utils,  only: int2str_sep, real2str_sig
use simple_gtk_utils,     only: gtkui_event_xy, gtkui_event_keyval, gtkui_set_cursor
use simple_gtk_draw,      only: gtkui_font_new, gtkui_font_free, gtkui_set_rgb, gtkui_draw_text
use gtk,   only: gtk_grid_new, gtk_grid_attach, gtk_grid_set_row_spacing, gtk_grid_set_column_spacing, gtk_label_new, &
    gtk_label_set_xalign, gtk_scale_new, gtk_scale_set_draw_value, gtk_spin_button_new, gtk_widget_set_size_request, &
    gtk_widget_set_hexpand, gtk_button_new_with_label, gtk_drawing_area_new, gtk_widget_set_can_focus, &
    gtk_widget_add_events, gtk_widget_queue_draw, gtk_widget_grab_focus, gtk_widget_has_focus, &
    gtk_widget_get_allocated_width, gtk_widget_get_allocated_height, gtk_widget_create_pango_layout, &
    gtk_widget_get_display, gtk_adjustment_new, gtk_adjustment_get_value, gtk_adjustment_set_value, g_signal_connect, &
    TRUE, FALSE, GTK_ORIENTATION_HORIZONTAL, GDK_BUTTON_PRESS_MASK, GDK_BUTTON_RELEASE_MASK, GDK_POINTER_MOTION_MASK, &
    GDK_KEY_PRESS_MASK, GDK_FOCUS_CHANGE_MASK
use g,     only: g_object_unref
use gdk,   only: gdk_cursor_new_from_name
use cairo, only: cairo_rectangle, cairo_fill, cairo_stroke, cairo_set_line_width, cairo_move_to, cairo_line_to, &
    cairo_set_dash, cairo_close_path, cairo_fill_preserve
implicit none

public :: gtkui_display_panel
private

type :: gtkui_display_panel
    private
    type(stack_display), pointer :: disp => null()
    type(c_ptr) :: root = c_null_ptr, hist = c_null_ptr, adj_bright = c_null_ptr, adj_contrast = c_null_ptr
    type(c_ptr) :: notify = c_null_ptr, font = c_null_ptr, cursor_ew = c_null_ptr
    real(dp) :: hb(4)    = 0._dp         ! histogram plot: left, top, right, bottom
    logical  :: hist_ok  = .false.
    logical  :: syncing  = .false.       ! adjustments set from the model, not by the user
    integer  :: drag     = 0             ! dragged limit: 1 Min, 2 Max, 0 none
    integer  :: active   = 1             ! the limit the arrow keys move
    real(dp) :: drag_off = 0._dp
  contains
    procedure :: new
    procedure :: get_widget
    procedure :: sync
    procedure :: kill
end type gtkui_display_panel

contains

    ! the widgets for disp; notify is redrawn on every change. The panel must stay at its address while
    ! the widgets live (their handlers keep it)
    subroutine new( self, disp, notify )
        class(gtkui_display_panel), target, intent(inout) :: self
        type(stack_display),        target, intent(inout) :: disp
        type(c_ptr),                        intent(in)    :: notify
        type(c_ptr) :: data, label, scale, spin, unit, button, adj
        integer     :: k
        call self%kill
        self%disp   => disp
        self%notify = notify
        self%font   = gtkui_font_new()
        data        = panel_cptr(self)
        self%root   = gtk_grid_new()
        call gtk_grid_set_row_spacing(self%root, 4_c_int)
        call gtk_grid_set_column_spacing(self%root, 8_c_int)
        ! brightness and contrast: a slider and a number field on one adjustment
        self%adj_bright = gtk_adjustment_new(real(disp%get_brightness(), c_double), &
            &real(GTKUI_BRIGHTNESS_LIM(1), c_double), real(GTKUI_BRIGHTNESS_LIM(2), c_double), &
            &1._c_double, 10._c_double, 0._c_double)
        self%adj_contrast = gtk_adjustment_new(real(disp%get_contrast(), c_double), &
            &real(GTKUI_CONTRAST_LIM(1), c_double), real(GTKUI_CONTRAST_LIM(2), c_double), &
            &0.01_c_double, 0.1_c_double, 0._c_double)
        call g_signal_connect(self%adj_bright,   'value-changed'//c_null_char, c_funloc(on_brightness), data)
        call g_signal_connect(self%adj_contrast, 'value-changed'//c_null_char, c_funloc(on_contrast),   data)
        do k = 1,2
            if( k == 1 )then
                adj   = self%adj_bright
                label = gtk_label_new('Brightness'//c_null_char)
                spin  = gtk_spin_button_new(adj, 0._c_double, 0_c_int)
                unit  = gtk_label_new('%'//c_null_char)
            else
                adj   = self%adj_contrast
                label = gtk_label_new('Contrast'//c_null_char)
                spin  = gtk_spin_button_new(adj, 0._c_double, 2_c_int)
                unit  = gtk_label_new(GTKUI_U_TIMES//c_null_char)
            endif
            scale = gtk_scale_new(GTK_ORIENTATION_HORIZONTAL, adj)
            call gtk_label_set_xalign(label, 0._c_float)
            call gtk_scale_set_draw_value(scale, FALSE)
            call gtk_widget_set_size_request(scale, int(GTKUI_SLIDER_W, c_int), -1_c_int)
            call gtk_widget_set_hexpand(scale, TRUE)
            call gtk_grid_attach(self%root, label, 0_c_int, int(k - 1, c_int), 1_c_int, 1_c_int)
            call gtk_grid_attach(self%root, scale, 1_c_int, int(k - 1, c_int), 1_c_int, 1_c_int)
            call gtk_grid_attach(self%root, spin,  2_c_int, int(k - 1, c_int), 1_c_int, 1_c_int)
            call gtk_grid_attach(self%root, unit,  3_c_int, int(k - 1, c_int), 1_c_int, 1_c_int)
        end do
        button = gtk_button_new_with_label('Reset'//c_null_char)
        call g_signal_connect(button, 'clicked'//c_null_char, c_funloc(on_reset), data)
        call gtk_grid_attach(self%root, button, 0_c_int, 2_c_int, 4_c_int, 1_c_int)
        ! the histogram with the draggable Min and Max
        self%hist = gtk_drawing_area_new()
        call gtk_widget_set_size_request(self%hist, int(GTKUI_HIST_W, c_int), int(GTKUI_HIST_H, c_int))
        call gtk_widget_set_can_focus(self%hist, TRUE)
        call gtk_widget_add_events(self%hist, ior(ior(ior(GDK_BUTTON_PRESS_MASK, GDK_BUTTON_RELEASE_MASK), &
            &ior(GDK_POINTER_MOTION_MASK, GDK_KEY_PRESS_MASK)), GDK_FOCUS_CHANGE_MASK))
        call g_signal_connect(self%hist, 'draw'//c_null_char,                 c_funloc(on_hist_draw),    data)
        call g_signal_connect(self%hist, 'button-press-event'//c_null_char,   c_funloc(on_hist_press),   data)
        call g_signal_connect(self%hist, 'motion-notify-event'//c_null_char,  c_funloc(on_hist_motion),  data)
        call g_signal_connect(self%hist, 'button-release-event'//c_null_char, c_funloc(on_hist_release), data)
        call g_signal_connect(self%hist, 'key-press-event'//c_null_char,      c_funloc(on_hist_key),     data)
        call g_signal_connect(self%hist, 'focus-in-event'//c_null_char,       c_funloc(on_hist_focus),   data)
        call g_signal_connect(self%hist, 'focus-out-event'//c_null_char,      c_funloc(on_hist_focus),   data)
        call gtk_grid_attach(self%root, self%hist, 0_c_int, 3_c_int, 4_c_int, 1_c_int)
    end subroutine new

    ! the panel, to attach into a popover or a window
    function get_widget( self ) result( widget )
        class(gtkui_display_panel), intent(in) :: self
        type(c_ptr) :: widget
        widget = self%root
    end function get_widget

    ! controls set from the model (before the panel is shown)
    subroutine sync( self )
        class(gtkui_display_panel), intent(inout) :: self
        call set_adjustment(self, self%adj_bright,   self%disp%get_brightness())
        call set_adjustment(self, self%adj_contrast, self%disp%get_contrast())
        call gtk_widget_queue_draw(self%hist)
    end subroutine sync

    ! the font and the cursor; the widgets belong to their window
    subroutine kill( self )
        class(gtkui_display_panel), intent(inout) :: self
        if( c_associated(self%font)      ) call gtkui_font_free(self%font)
        if( c_associated(self%cursor_ew) ) call g_object_unref(self%cursor_ew)
        self%font      = c_null_ptr
        self%cursor_ew = c_null_ptr
        self%root      = c_null_ptr
        self%hist      = c_null_ptr
        self%disp      => null()
    end subroutine kill

    ! ---- private -----------------------------------------------------------

    function panel_cptr( p ) result( ptr )
        type(gtkui_display_panel), target, intent(inout) :: p
        type(c_ptr) :: ptr
        ptr = c_loc(p)
    end function panel_cptr

    subroutine set_adjustment( p, adj, val )
        type(gtkui_display_panel), intent(inout) :: p
        type(c_ptr),               intent(in)    :: adj
        real(dp),                  intent(in)    :: val
        p%syncing = .true.
        call gtk_adjustment_set_value(adj, real(val, c_double))
        p%syncing = .false.
    end subroutine set_adjustment

    subroutine changed( p )
        type(gtkui_display_panel), intent(inout) :: p
        if( c_associated(p%notify) ) call gtk_widget_queue_draw(p%notify)
        call gtk_widget_queue_draw(p%hist)
    end subroutine changed

    real(dp) function hist_x( p, value )
        type(gtkui_display_panel), intent(in) :: p
        real(dp),                  intent(in) :: value
        real(dp) :: lims(2)
        lims   = p%disp%get_range()
        hist_x = p%hb(1) + (value - lims(1)) / (lims(2) - lims(1)) * (p%hb(3) - p%hb(1))
    end function hist_x

    ! the limit under (x, y): 1 Min, 2 Max, 0 none; the handles at the two ends keep coincident limits apart
    integer function hist_hit( p, x, y )
        type(gtkui_display_panel), intent(in) :: p
        real(dp),                  intent(in) :: x, y
        real(dp) :: win(2), d(2)
        logical  :: hit(2)
        hist_hit = 0
        if( .not. p%hist_ok ) return
        if( .not. p%disp%is_usable() ) return
        if( y < p%hb(2) - GTKUI_HIT_RADIUS .or. y > p%hb(4) + GTKUI_HIT_RADIUS ) return
        win = p%disp%get_window()
        d   = [abs(x - hist_x(p, win(1))), abs(x - hist_x(p, win(2)))]
        hit = d <= real(GTKUI_HIT_RADIUS, dp)
        if( hit(1) .and. y <= p%hb(2) + 12._dp )then
            hist_hit = 1
        else if( hit(2) .and. y >= p%hb(4) - 12._dp )then
            hist_hit = 2
        else if( hit(1) .and. hit(2) )then
            hist_hit = p%active
            if( d(1) < d(2) ) hist_hit = 1
            if( d(2) < d(1) ) hist_hit = 2
        else if( hit(1) )then
            hist_hit = 1
        else if( hit(2) )then
            hist_hit = 2
        endif
    end function hist_hit

    subroutine drag_to( p, x )
        type(gtkui_display_panel), intent(inout) :: p
        real(dp),                  intent(in)    :: x
        real(dp) :: lims(2), v
        if( .not. p%hist_ok ) return
        lims = p%disp%get_range()
        v    = lims(1) + (lims(2) - lims(1)) * ((x - p%drag_off - p%hb(1)) / (p%hb(3) - p%hb(1)))
        if( p%disp%set_window_limit(p%drag, v) ) call changed(p)
    end subroutine drag_to

    subroutine draw_markers( p, cr, layout )
        type(gtkui_display_panel), intent(inout) :: p
        type(c_ptr),               intent(in)    :: cr, layout
        real(c_double), target :: dashes(2)
        real(dp) :: win(2), x, y, dir, rgb(3), top, bottom
        logical  :: l_focus, l_active
        integer  :: k
        top     = p%hb(2)
        bottom  = p%hb(4)
        win     = p%disp%get_window()
        l_focus = gtk_widget_has_focus(p%hist) /= FALSE
        dashes  = [3._c_double, 2._c_double]
        do k = 1,2
            x        = hist_x(p, win(k))
            l_active = l_focus .and. k == p%active
            if( k == 1 )then
                rgb = GTKUI_COLOR_MIN
            else
                rgb = GTKUI_COLOR_MAX
            endif
            call gtkui_set_rgb(cr, rgb)
            call cairo_set_dash(cr, c_loc(dashes), 2_c_int, 0._c_double)
            call cairo_set_line_width(cr, merge(3._dp, 1._dp, l_active))
            call cairo_move_to(cr, x, top)
            call cairo_line_to(cr, x, bottom)
            call cairo_stroke(cr)
            call cairo_set_dash(cr, c_null_ptr, 0_c_int, 0._c_double)
            if( k == 1 )then
                call gtkui_draw_text(cr, layout, p%font, GTKUI_SMALL_PX, 'Min '//real2str_sig(win(k), 6), x, top - 3._dp, &
                    &'sw', rgb)
                y   = top
                dir = 1._dp
            else
                call gtkui_draw_text(cr, layout, p%font, GTKUI_SMALL_PX, 'Max '//real2str_sig(win(k), 6), x, top - 3._dp, &
                    &'se', rgb)
                y   = bottom
                dir = -1._dp
            endif
            call gtkui_set_rgb(cr, rgb)
            call cairo_move_to(cr, x - 5._dp, y)
            call cairo_line_to(cr, x + 5._dp, y)
            call cairo_line_to(cr, x, y + dir * 9._dp)
            call cairo_close_path(cr)
            call cairo_fill_preserve(cr)
            call cairo_set_line_width(cr, merge(2._dp, 1._dp, l_active))
            call cairo_stroke(cr)
        end do
    end subroutine draw_markers

    ! ---- GTK handlers (user data: the panel) -------------------------------

    subroutine on_brightness( adj, data ) bind(c)
        type(c_ptr), value :: adj, data
        type(gtkui_display_panel), pointer :: p
        real(dp) :: v
        call c_f_pointer(data, p)
        if( p%syncing ) return
        v = real(gtk_adjustment_get_value(adj), dp)
        if( p%disp%set_brightness(v) ) call changed(p)
        if( p%disp%get_brightness() /= v ) call set_adjustment(p, adj, p%disp%get_brightness())
    end subroutine on_brightness

    subroutine on_contrast( adj, data ) bind(c)
        type(c_ptr), value :: adj, data
        type(gtkui_display_panel), pointer :: p
        real(dp) :: v
        call c_f_pointer(data, p)
        if( p%syncing ) return
        v = real(gtk_adjustment_get_value(adj), dp)
        if( p%disp%set_contrast(v) ) call changed(p)
        if( p%disp%get_contrast() /= v ) call set_adjustment(p, adj, p%disp%get_contrast())
    end subroutine on_contrast

    subroutine on_reset( button, data ) bind(c)
        type(c_ptr), value :: button, data
        type(gtkui_display_panel), pointer :: p
        call c_f_pointer(data, p)
        call p%disp%reset()
        call set_adjustment(p, p%adj_bright,   p%disp%get_brightness())
        call set_adjustment(p, p%adj_contrast, p%disp%get_contrast())
        call changed(p)
    end subroutine on_reset

    function on_hist_draw( widget, cr, data ) bind(c) result( ret )
        type(c_ptr), value :: widget, cr, data
        integer(c_int)     :: ret
        type(gtkui_display_panel), pointer :: p
        type(c_ptr)     :: layout
        integer(kind=8) :: counts(GTKUI_HIST_NBINS), cmax
        real(dp)        :: width, height, left, top, right, bottom, x0, x1, lims(2), black, white
        integer         :: k
        call c_f_pointer(data, p)
        ret       = FALSE
        p%hist_ok = .false.
        width  = real(gtk_widget_get_allocated_width(widget),  dp)
        height = real(gtk_widget_get_allocated_height(widget), dp)
        layout = gtk_widget_create_pango_layout(widget, c_null_char)
        call gtkui_draw_text(cr, layout, p%font, GTKUI_SMALL_PX, 'Pixels', 8._dp, 4._dp, 'nw', GTKUI_COLOR_FG)
        call gtkui_draw_text(cr, layout, p%font, GTKUI_SMALL_PX, 'Raw intensity '//GTKUI_U_MIDDOT//' box: display range', &
            &0.5_dp * width, height - 3._dp, 's', GTKUI_COLOR_FG)
        counts = p%disp%get_counts()
        cmax   = maxval(counts)
        left   = 53._dp
        top    = 33._dp
        right  = width - 15._dp
        bottom = height - 43._dp
        if( cmax == 0_8 )then
            call gtkui_draw_text(cr, layout, p%font, GTKUI_MESSAGE_PX, 'No finite source pixels', 0.5_dp * width, &
                &0.5_dp * height, 'c', GTKUI_COLOR_FG)
        else if( right > left .and. bottom > top )then
            p%hb      = [left, top, right, bottom]
            p%hist_ok = .true.
            call gtkui_draw_text(cr, layout, p%font, GTKUI_SMALL_PX, 'Drag limits '//GTKUI_U_MIDDOT//' '//GTKUI_U_UP//'/'//&
                &GTKUI_U_DOWN//' select '//GTKUI_U_MIDDOT//' '//GTKUI_U_LEFT//'/'//GTKUI_U_RIGHT//' adjust', width - 8._dp, &
                &4._dp, 'ne', GTKUI_COLOR_FG)
            call gtkui_set_rgb(cr, GTKUI_COLOR_BARS)
            do k = 1,GTKUI_HIST_NBINS
                if( counts(k) == 0_8 ) cycle
                x0 = left + real(k - 1, dp) / real(GTKUI_HIST_NBINS, dp) * (right - left)
                x1 = left + real(k, dp)     / real(GTKUI_HIST_NBINS, dp) * (right - left)
                call cairo_rectangle(cr, x0, bottom - real(counts(k), dp) / real(cmax, dp) * (bottom - top), &
                    &max(x0 + 1._dp, x1) - x0, real(counts(k), dp) / real(cmax, dp) * (bottom - top))
            end do
            call cairo_fill(cr)
            call gtkui_set_rgb(cr, GTKUI_COLOR_FG)
            call cairo_set_line_width(cr, 1._dp)
            call cairo_move_to(cr, left + 0.5_dp, top)
            call cairo_line_to(cr, left + 0.5_dp, bottom + 0.5_dp)
            call cairo_line_to(cr, right, bottom + 0.5_dp)
            call cairo_stroke(cr)
            call gtkui_draw_text(cr, layout, p%font, GTKUI_SMALL_PX, int2str_sep(cmax), left - 5._dp, top, 'e', GTKUI_COLOR_FG)
            call gtkui_draw_text(cr, layout, p%font, GTKUI_SMALL_PX, '0', left - 5._dp, bottom, 'e', GTKUI_COLOR_FG)
            lims = p%disp%get_range()
            call gtkui_draw_text(cr, layout, p%font, GTKUI_SMALL_PX, real2str_sig(lims(1), 6), left,  bottom + 5._dp, 'nw', &
                &GTKUI_COLOR_FG)
            call gtkui_draw_text(cr, layout, p%font, GTKUI_SMALL_PX, real2str_sig(lims(2), 6), right, bottom + 5._dp, 'ne', &
                &GTKUI_COLOR_FG)
            ! the raw range mapped from black to white
            call p%disp%mapping_box(black, white)
            if( white > lims(1) .and. black < lims(2) )then
                x0 = hist_x(p, max(lims(1), black))
                x1 = hist_x(p, min(lims(2), white))
                call gtkui_set_rgb(cr, GTKUI_COLOR_RANGE)
                call cairo_set_line_width(cr, 2._dp)
                call cairo_rectangle(cr, x0, top, x1 - x0, bottom - top)
                call cairo_stroke(cr)
            endif
            call draw_markers(p, cr, layout)
        endif
        call g_object_unref(layout)
    end function on_hist_draw

    function on_hist_press( widget, event, data ) bind(c) result( handled )
        type(c_ptr), value :: widget, event, data
        integer(c_int)     :: handled
        type(gtkui_display_panel), pointer :: p
        real(dp) :: x, y, win(2)
        integer  :: k
        call c_f_pointer(data, p)
        handled = FALSE
        call gtk_widget_grab_focus(widget)
        if( .not. gtkui_event_xy(event, x, y) ) return
        k = hist_hit(p, x, y)
        if( k == 0 ) return
        win        = p%disp%get_window()
        p%active   = k
        p%drag     = k
        p%drag_off = x - hist_x(p, win(k))
        handled    = TRUE
        call gtk_widget_queue_draw(p%hist)
    end function on_hist_press

    function on_hist_motion( widget, event, data ) bind(c) result( handled )
        type(c_ptr), value :: widget, event, data
        integer(c_int)     :: handled
        type(gtkui_display_panel), pointer :: p
        real(dp) :: x, y
        call c_f_pointer(data, p)
        handled = FALSE
        if( .not. gtkui_event_xy(event, x, y) ) return
        if( p%drag == 0 )then
            if( .not. c_associated(p%cursor_ew) ) &
                &p%cursor_ew = gdk_cursor_new_from_name(gtk_widget_get_display(widget), 'ew-resize'//c_null_char)
            if( hist_hit(p, x, y) > 0 )then
                call gtkui_set_cursor(widget, p%cursor_ew)
            else
                call gtkui_set_cursor(widget, c_null_ptr)
            endif
            return
        endif
        call drag_to(p, x)
        handled = TRUE
    end function on_hist_motion

    function on_hist_release( widget, event, data ) bind(c) result( handled )
        type(c_ptr), value :: widget, event, data
        integer(c_int)     :: handled
        type(gtkui_display_panel), pointer :: p
        real(dp) :: x, y
        call c_f_pointer(data, p)
        handled = FALSE
        if( p%drag == 0 ) return
        if( gtkui_event_xy(event, x, y) ) call drag_to(p, x)
        p%drag  = 0
        handled = TRUE
    end function on_hist_release

    ! Up/Down choose Min/Max, Left/Right move it by a thousandth of the range
    function on_hist_key( widget, event, data ) bind(c) result( handled )
        type(c_ptr), value :: widget, event, data
        integer(c_int)     :: handled
        type(gtkui_display_panel), pointer :: p
        integer(c_int) :: keyval
        real(dp)       :: lims(2), win(2), step
        call c_f_pointer(data, p)
        handled = FALSE
        keyval  = gtkui_event_keyval(event)
        lims    = p%disp%get_range()
        win     = p%disp%get_window()
        step    = (lims(2) - lims(1)) / 1000._dp
        if( keyval == GTKUI_KEY_UP )then
            p%active = 1
        else if( keyval == GTKUI_KEY_DOWN )then
            p%active = 2
        else if( keyval == GTKUI_KEY_LEFT )then
            if( p%disp%set_window_limit(p%active, win(p%active) - step) ) call changed(p)
        else if( keyval == GTKUI_KEY_RIGHT )then
            if( p%disp%set_window_limit(p%active, win(p%active) + step) ) call changed(p)
        else
            return
        endif
        call gtk_widget_queue_draw(p%hist)
        handled = TRUE
    end function on_hist_key

    function on_hist_focus( widget, event, data ) bind(c) result( handled )
        type(c_ptr), value :: widget, event, data
        integer(c_int)     :: handled
        call gtk_widget_queue_draw(widget)
        handled = FALSE
    end function on_hist_focus

end module simple_gtk_display_panel
