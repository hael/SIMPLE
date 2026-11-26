!@descr: a GTK 3 scrolling grid of image tiles with captions, selection marks, zoom and hover, fed by a delegate
! One drawing area draws only the visible rows, beside its own scrollbar (select2d.py's layout; geometry and
! colours in simple_defs_gtk). An app extends gtkui_grid_delegate to say what is shown and to receive clicks,
! and attaches its context popover so the popover's grab does not swallow the grid's clicks while it is up.
module simple_gtk_image_grid
use, intrinsic :: iso_c_binding, only: c_ptr, c_int, c_double, c_null_char, c_null_ptr, c_funloc, c_loc, &
    c_f_pointer, c_associated
use simple_defs,      only: dp
use simple_defs_gtk,  only: GTKUI_GRID_MARGIN, GTKUI_GRID_SPACING, GTKUI_GRID_BORDER, GTKUI_TILE_START, &
    GTKUI_TILE_MIN, GTKUI_TILE_MAX, GTKUI_ZOOM_FACTOR, GTKUI_PT2PX, GTKUI_MESSAGE_PX, GTKUI_COLOR_BG, GTKUI_COLOR_FG, &
    GTKUI_COLOR_ACCENT, GTKUI_COLOR_BACKDROP, GTKUI_COLOR_SURFACE, GTKUI_COLOR_BORDER, GTKUI_COLOR_HOVER, &
    GTKUI_COLOR_MUTED, GTKUI_COLOR_SELECTED
use simple_gtk_utils, only: gtkui_bool, gtkui_event_xy, gtkui_event_state, gtkui_event_button, gtkui_event_scroll, &
    gtkui_event_is_press, gtkui_event_context_menu
use simple_gtk_draw,  only: gtkui_font_new, gtkui_font_free, gtkui_set_rgb, gtkui_text_size, gtkui_draw_text, &
    gtkui_fit_text, gtkui_draw_check, gtkui_gray_surface
use gtk,   only: gtk_box_new, gtk_box_pack_start, gtk_drawing_area_new, gtk_scrollbar_new, gtk_adjustment_new, &
    gtk_adjustment_get_value, gtk_adjustment_set_value, gtk_adjustment_configure, gtk_adjustment_get_page_size, &
    gtk_adjustment_get_upper, gtk_widget_set_can_focus, gtk_widget_add_events, gtk_widget_grab_focus, &
    gtk_widget_queue_draw, gtk_widget_get_allocated_width, gtk_widget_get_allocated_height, &
    gtk_widget_create_pango_layout, gtk_get_event_widget, gtk_popover_popdown, g_signal_connect, TRUE, FALSE, &
    GTK_ORIENTATION_HORIZONTAL, GTK_ORIENTATION_VERTICAL, GDK_SHIFT_MASK, GDK_CONTROL_MASK, GDK_META_MASK, GDK_MOD2_MASK, &
    GDK_BUTTON_PRESS_MASK, GDK_POINTER_MOTION_MASK, GDK_LEAVE_NOTIFY_MASK, GDK_SCROLL_MASK, GDK_SMOOTH_SCROLL_MASK
use g,     only: g_object_unref
use cairo, only: cairo_paint, cairo_rectangle, cairo_fill, cairo_stroke, cairo_set_line_width, &
    cairo_set_source_surface, cairo_surface_destroy
implicit none

public :: gtkui_image_grid, gtkui_grid_delegate
private

! what the grid shows and whom it tells: items 1..nids, nitems of them in display order
type, abstract :: gtkui_grid_delegate
  contains
    procedure(grid_count),   deferred :: nitems       ! items shown
    procedure(grid_count),   deferred :: nids         ! item identifiers are 1..nids
    procedure(grid_item_at), deferred :: item_at      ! the item at display position pos (1..nitems)
    procedure(grid_flag),    deferred :: is_selected
    procedure(grid_box),     deferred :: image_box    ! the items' image size (nx, ny)
    procedure(grid_count),   deferred :: revision     ! changes whenever the rendered images change
    procedure(grid_render),  deferred :: render       ! an item's image as 8 bits at size(gray)
    procedure(grid_count),   deferred :: ncaptions    ! caption fields below each image
    procedure(grid_caption), deferred :: caption      ! field k of an item (item 0: the widest value expected)
    procedure(grid_click),   deferred :: click        ! a click on an item, with or without Shift
    procedure(grid_context), deferred :: context      ! a context-menu request at (x, y) of the grid, also while
                                                      ! the attached popover is up (re-point it there)
end type gtkui_grid_delegate

abstract interface
    integer function grid_count( self )
        import :: gtkui_grid_delegate
        class(gtkui_grid_delegate), intent(in) :: self
    end function grid_count
    integer function grid_item_at( self, pos )
        import :: gtkui_grid_delegate
        class(gtkui_grid_delegate), intent(in) :: self
        integer,                    intent(in) :: pos
    end function grid_item_at
    logical function grid_flag( self, item )
        import :: gtkui_grid_delegate
        class(gtkui_grid_delegate), intent(in) :: self
        integer,                    intent(in) :: item
    end function grid_flag
    function grid_box( self ) result( box )
        import :: gtkui_grid_delegate
        class(gtkui_grid_delegate), intent(in) :: self
        integer :: box(2)
    end function grid_box
    subroutine grid_render( self, item, gray )
        import :: gtkui_grid_delegate
        class(gtkui_grid_delegate), intent(in)  :: self
        integer,                    intent(in)  :: item
        integer,                    intent(out) :: gray(:,:)
    end subroutine grid_render
    subroutine grid_caption( self, k, item, label, short, value )
        import :: gtkui_grid_delegate
        class(gtkui_grid_delegate),    intent(in)  :: self
        integer,                       intent(in)  :: k, item
        character(len=:), allocatable, intent(out) :: label, short, value
    end subroutine grid_caption
    subroutine grid_click( self, item, shift )
        import :: gtkui_grid_delegate
        class(gtkui_grid_delegate), intent(inout) :: self
        integer,                    intent(in)    :: item
        logical,                    intent(in)    :: shift
    end subroutine grid_click
    subroutine grid_context( self, x, y )
        import :: gtkui_grid_delegate, dp
        class(gtkui_grid_delegate), intent(inout) :: self
        real(dp),                   intent(in)    :: x, y
    end subroutine grid_context
end interface

! a rendered tile, valid for one tile size and one revision of the delegate's images
type :: tile_cache
    type(c_ptr) :: surf = c_null_ptr
    integer     :: tile = 0
    integer     :: rev  = -1
end type tile_cache

! caption geometry for the current tile size
type :: caption_layout
    logical  :: rows     = .false.       ! label and value side by side, one field per row
    real(dp) :: label_px = 10._dp
    real(dp) :: value_px = 12.5_dp
    integer  :: pad = 10, vpad = 6, gap = 6, row_h = 0, height = 0
end type caption_layout

type :: gtkui_image_grid
    private
    class(gtkui_grid_delegate), pointer :: dlg => null()
    type(c_ptr) :: root = c_null_ptr, area = c_null_ptr, vadj = c_null_ptr, font = c_null_ptr
    type(tile_cache), allocatable :: cache(:)
    type(caption_layout) :: cap
    integer  :: tile = GTKUI_TILE_START, cols = 1, card_h = 0, width = 0, height = 0
    integer  :: hover = 0                ! item under the pointer
    real(dp) :: zoom_acc = 0._dp         ! smooth-scroll steps not yet applied as zoom
  contains
    procedure :: new
    procedure :: attach_popover
    procedure :: get_widget
    procedure :: get_area
    procedure :: relayout
    procedure :: redraw
    procedure :: kill
end type gtkui_image_grid

contains

    ! the widgets; the grid object must stay at its address while they live (its GTK handlers keep it)
    subroutine new( self, dlg )
        class(gtkui_image_grid),            target, intent(inout) :: self
        class(gtkui_grid_delegate),         target, intent(inout) :: dlg
        type(c_ptr) :: data, sbar
        call self%kill
        self%dlg  => dlg
        allocate(self%cache(dlg%nids()))
        self%font = gtkui_font_new()
        data      = grid_cptr(self)
        self%root = gtk_box_new(GTK_ORIENTATION_HORIZONTAL, 0_c_int)
        self%area = gtk_drawing_area_new()
        call gtk_widget_set_can_focus(self%area, TRUE)
        call gtk_widget_add_events(self%area, ior(ior(ior(GDK_BUTTON_PRESS_MASK, GDK_POINTER_MOTION_MASK), &
            &ior(GDK_LEAVE_NOTIFY_MASK, GDK_SCROLL_MASK)), GDK_SMOOTH_SCROLL_MASK))
        call g_signal_connect(self%area, 'draw'//c_null_char,                c_funloc(on_draw),     data)
        call g_signal_connect(self%area, 'size-allocate'//c_null_char,       c_funloc(on_allocate), data)
        call g_signal_connect(self%area, 'button-press-event'//c_null_char,  c_funloc(on_press),    data)
        call g_signal_connect(self%area, 'motion-notify-event'//c_null_char, c_funloc(on_motion),   data)
        call g_signal_connect(self%area, 'leave-notify-event'//c_null_char,  c_funloc(on_leave),    data)
        call g_signal_connect(self%area, 'scroll-event'//c_null_char,        c_funloc(on_scroll),   data)
        call gtk_box_pack_start(self%root, self%area, TRUE, TRUE, 0_c_int)
        self%vadj = gtk_adjustment_new(0._c_double, 0._c_double, 1._c_double, 30._c_double, 100._c_double, 1._c_double)
        call g_signal_connect(self%vadj, 'value-changed'//c_null_char, c_funloc(on_vadj), data)
        sbar = gtk_scrollbar_new(GTK_ORIENTATION_VERTICAL, self%vadj)
        call gtk_box_pack_start(self%root, sbar, FALSE, FALSE, 0_c_int)
    end subroutine new

    ! the app's context popover over the grid. While a modal popover (GTK's default) is up it holds GTK's grab,
    ! so the grid's clicks go to the popover, not to on_press: they are taken there by on_popover_press
    subroutine attach_popover( self, popover )
        class(gtkui_image_grid), target, intent(inout) :: self
        type(c_ptr),                     intent(in)    :: popover
        type(c_ptr) :: data
        data = grid_cptr(self)
        call g_signal_connect(popover, 'button-press-event'//c_null_char,   c_funloc(on_popover_press),   data)
        call g_signal_connect(popover, 'button-release-event'//c_null_char, c_funloc(on_popover_release), data)
    end subroutine attach_popover

    ! the grid with its scrollbar, to pack into a window
    function get_widget( self ) result( widget )
        class(gtkui_image_grid), intent(in) :: self
        type(c_ptr) :: widget
        widget = self%root
    end function get_widget

    ! the drawing area: popovers point into it in its coordinates
    function get_area( self ) result( widget )
        class(gtkui_image_grid), intent(in) :: self
        type(c_ptr) :: widget
        widget = self%area
    end function get_area

    ! after the number of items changed: layout, scroll range and redraw
    subroutine relayout( self )
        class(gtkui_image_grid), intent(inout) :: self
        call update_layout(self)
    end subroutine relayout

    ! after a change of order, selection or images
    subroutine redraw( self )
        class(gtkui_image_grid), intent(inout) :: self
        if( c_associated(self%area) ) call gtk_widget_queue_draw(self%area)
    end subroutine redraw

    ! the rendered tiles and the font; the widgets belong to their window
    subroutine kill( self )
        class(gtkui_image_grid), intent(inout) :: self
        integer :: i
        if( allocated(self%cache) )then
            do i = 1,size(self%cache)
                if( c_associated(self%cache(i)%surf) ) call cairo_surface_destroy(self%cache(i)%surf)
            end do
            deallocate(self%cache)
        endif
        if( c_associated(self%font) ) call gtkui_font_free(self%font)
        self%font = c_null_ptr
        self%root = c_null_ptr
        self%area = c_null_ptr
        self%vadj = c_null_ptr
        self%dlg  => null()
    end subroutine kill

    ! ---- private -----------------------------------------------------------

    ! the C address of a grid, for the user data of its handlers
    function grid_cptr( gr ) result( ptr )
        type(gtkui_image_grid), target, intent(inout) :: gr
        type(c_ptr) :: ptr
        ptr = c_loc(gr)
    end function grid_cptr

    ! tile size within [GTKUI_TILE_MIN, min(GTKUI_TILE_MAX, width - 2 margins)], then layout when it changed
    subroutine set_tile( gr, tile, l_force )
        type(gtkui_image_grid), intent(inout) :: gr
        integer,                intent(in)    :: tile
        logical,                intent(in)    :: l_force
        integer :: t, tmax
        logical :: l_changed
        tmax = max(GTKUI_TILE_MIN, min(GTKUI_TILE_MAX, gr%width - 2 * GTKUI_GRID_MARGIN))
        t    = max(GTKUI_TILE_MIN, min(tmax, tile))
        l_changed = t /= gr%tile
        gr%tile   = t
        if( l_changed .or. l_force ) call update_layout(gr)
    end subroutine set_tile

    subroutine update_layout( gr )
        type(gtkui_image_grid), intent(inout) :: gr
        type(c_ptr) :: layout
        real(dp)    :: page, value, content
        integer     :: nitems, nrows
        layout = gtk_widget_create_pango_layout(gr%area, c_null_char)
        call caption_metrics(gr, layout)
        call g_object_unref(layout)
        gr%card_h = gr%tile + gr%cap%height
        nitems    = gr%dlg%nitems()
        gr%cols   = max(1, min((max(1, gr%width - 2 * GTKUI_GRID_MARGIN) + GTKUI_GRID_SPACING) / &
            &(gr%tile + GTKUI_GRID_SPACING), max(nitems, 1)))
        nrows     = (nitems + gr%cols - 1) / gr%cols
        content   = 0._dp
        if( nrows > 0 ) content = real(2 * GTKUI_GRID_MARGIN + nrows * gr%card_h + (nrows - 1) * GTKUI_GRID_SPACING, dp)
        page      = real(max(1, gr%height), dp)
        value     = min(real(gtk_adjustment_get_value(gr%vadj), dp), max(0._dp, content - page))
        call gtk_adjustment_configure(gr%vadj, real(value, c_double), 0._c_double, real(max(content, page), c_double), &
            &30._c_double, real(0.9_dp * page, c_double), real(page, c_double))
        call gtk_widget_queue_draw(gr%area)
    end subroutine update_layout

    ! caption fonts and height for the current tile size: select2d.py's sizes, shrunk on small tiles until the
    ! widest expected values fit
    subroutine caption_metrics( gr, layout )
        type(gtkui_image_grid), intent(inout) :: gr
        type(c_ptr),            intent(in)    :: layout
        character(len=:), allocatable :: label, short, sample
        logical :: compact, tiny, fits
        integer :: tinyfs, label_fs, value_fs, fs, tw, lw, vw, lh, vh, k, ncap
        ncap     = gr%dlg%ncaptions()
        compact  = gr%tile < 180
        tiny     = gr%tile < 90
        tinyfs   = merge(6, 7, gr%tile < 65)
        label_fs = merge(tinyfs, 8, tiny)
        value_fs = merge(tinyfs, merge(8, 10, compact), tiny)
        if( tiny )then
            gr%cap%pad = merge(3, 4, gr%tile < 65)
        else
            gr%cap%pad = merge(7, 10, compact)
        endif
        gr%cap%vpad = merge(4, 6, tiny)
        gr%cap%gap  = merge(3, 6, tiny)
        gr%cap%rows = compact
        if( compact )then
            fs = value_fs
            do while( fs > 1 )
                fits = .true.
                do k = 1,ncap
                    call gr%dlg%caption(k, 0, label, short, sample)
                    call gtkui_text_size(layout, gr%font, real(fs, dp) * GTKUI_PT2PX, short,  lw, lh)
                    call gtkui_text_size(layout, gr%font, real(fs, dp) * GTKUI_PT2PX, sample, vw, vh)
                    ! Tk's text bounds reach one pixel beyond each measured edge
                    if( lw + gr%cap%gap + vw + 2 > gr%tile - 2 * gr%cap%pad ) fits = .false.
                end do
                if( fits ) exit
                fs       = fs - 1
                label_fs = fs
                value_fs = fs
            end do
        endif
        gr%cap%label_px = real(label_fs, dp) * GTKUI_PT2PX
        gr%cap%value_px = real(value_fs, dp) * GTKUI_PT2PX
        call gtkui_text_size(layout, gr%font, gr%cap%label_px, 'Ag', tw, lh)
        call gtkui_text_size(layout, gr%font, gr%cap%value_px, 'Ag', tw, vh)
        gr%cap%row_h  = max(lh, vh)
        gr%cap%height = 0
        if( ncap > 0 ) gr%cap%height = merge(ncap, 2, gr%cap%rows) * gr%cap%row_h + 2 * gr%cap%vpad
    end subroutine caption_metrics

    ! the display position under widget point (x, y), 0 for none
    integer function hit_position( gr, x, y )
        type(gtkui_image_grid), intent(in) :: gr
        real(dp),               intent(in) :: x, y
        real(dp) :: yc, rx, ry
        integer  :: col, row, pos
        hit_position = 0
        if( gr%card_h == 0 ) return
        yc  = y + anint(real(gtk_adjustment_get_value(gr%vadj), dp))
        col = floor((x - real(GTKUI_GRID_MARGIN, dp)) / real(gr%tile + GTKUI_GRID_SPACING, dp))
        row = floor((yc - real(GTKUI_GRID_MARGIN, dp)) / real(gr%card_h + GTKUI_GRID_SPACING, dp))
        if( col < 0 .or. col >= gr%cols .or. row < 0 ) return
        rx = x  - real(GTKUI_GRID_MARGIN + col * (gr%tile + GTKUI_GRID_SPACING), dp)
        ry = yc - real(GTKUI_GRID_MARGIN + row * (gr%card_h + GTKUI_GRID_SPACING), dp)
        if( rx >= real(gr%tile, dp) .or. ry >= real(gr%card_h, dp) ) return
        pos = row * gr%cols + col + 1
        if( pos <= gr%dlg%nitems() ) hit_position = pos
    end function hit_position

    ! one item: card, image, caption, outline and selection mark
    subroutine draw_card( gr, cr, layout, item, x, y )
        type(gtkui_image_grid), intent(inout) :: gr
        type(c_ptr),            intent(in)    :: cr, layout
        integer,                intent(in)    :: item
        real(dp),               intent(in)    :: x, y
        type(c_ptr) :: surf
        integer     :: ow, oh
        real(dp)    :: t, h, px, py
        logical     :: l_sel
        t     = real(gr%tile, dp)
        h     = real(gr%card_h, dp)
        l_sel = gr%dlg%is_selected(item)
        if( l_sel )then
            call gtkui_set_rgb(cr, GTKUI_COLOR_SELECTED)
        else
            call gtkui_set_rgb(cr, GTKUI_COLOR_SURFACE)
        endif
        call cairo_rectangle(cr, x, y, t, h)
        call cairo_fill(cr)
        call gtkui_set_rgb(cr, GTKUI_COLOR_BACKDROP)
        call cairo_rectangle(cr, x + GTKUI_GRID_BORDER, y + GTKUI_GRID_BORDER, t - 2 * GTKUI_GRID_BORDER, &
            &t - 2 * GTKUI_GRID_BORDER)
        call cairo_fill(cr)
        surf = tile_surface(gr, item, ow, oh)
        px   = x + real((gr%tile - ow) / 2, dp)
        py   = y + real((gr%tile - oh) / 2, dp)
        call cairo_set_source_surface(cr, surf, px, py)
        call cairo_rectangle(cr, px, py, real(ow, dp), real(oh, dp))
        call cairo_fill(cr)
        if( gr%cap%height > 0 ) call draw_caption(gr, cr, layout, item, x, y)
        if( l_sel )then
            call gtkui_set_rgb(cr, GTKUI_COLOR_ACCENT)
            call cairo_set_line_width(cr, 2._dp)
            call cairo_rectangle(cr, x + 1._dp, y + 1._dp, t - 2._dp, h - 2._dp)
        else
            if( item == gr%hover )then
                call gtkui_set_rgb(cr, GTKUI_COLOR_HOVER)
            else
                call gtkui_set_rgb(cr, GTKUI_COLOR_BORDER)
            endif
            call cairo_set_line_width(cr, 1._dp)
            call cairo_rectangle(cr, x + 0.5_dp, y + 0.5_dp, t - 1._dp, h - 1._dp)
        endif
        call cairo_stroke(cr)
        if( l_sel ) call gtkui_draw_check(cr, x, y, gr%tile)
    end subroutine draw_card

    ! the item's tile at the current size and image revision, rendered when stale
    function tile_surface( gr, item, ow, oh ) result( surf )
        type(gtkui_image_grid), intent(inout) :: gr
        integer,                intent(in)    :: item
        integer,                intent(out)   :: ow, oh
        type(c_ptr) :: surf
        integer, allocatable :: gray(:,:)
        integer  :: box(2), inner
        real(dp) :: scale
        box   = gr%dlg%image_box()
        inner = max(1, gr%tile - 2 * GTKUI_GRID_BORDER)
        scale = min(real(inner, dp) / real(box(1), dp), real(inner, dp) / real(box(2), dp))
        ow    = max(1, nint(real(box(1), dp) * scale))
        oh    = max(1, nint(real(box(2), dp) * scale))
        surf  = gr%cache(item)%surf
        if( c_associated(surf) )then
            if( gr%cache(item)%tile == gr%tile .and. gr%cache(item)%rev == gr%dlg%revision() ) return
            call cairo_surface_destroy(surf)
        endif
        allocate(gray(ow,oh))
        call gr%dlg%render(item, gray)
        surf = gtkui_gray_surface(gray)
        gr%cache(item)%surf = surf
        gr%cache(item)%tile = gr%tile
        gr%cache(item)%rev  = gr%dlg%revision()
    end function tile_surface

    ! the caption fields below the image: one per row on small tiles, side by side on large ones
    subroutine draw_caption( gr, cr, layout, item, x, y )
        type(gtkui_image_grid), intent(inout) :: gr
        type(c_ptr),            intent(in)    :: cr, layout
        integer,                intent(in)    :: item
        real(dp),               intent(in)    :: x, y
        character(len=:), allocatable :: label, short, value, txt
        integer  :: k, ncap, lw, th
        real(dp) :: top, yk, fw, cx
        ncap = gr%dlg%ncaptions()
        top  = y + real(gr%tile + gr%cap%vpad, dp)
        fw   = real(gr%tile - 2 * gr%cap%pad, dp) / real(max(ncap, 1), dp)
        do k = 1,ncap
            call gr%dlg%caption(k, item, label, short, value)
            if( gr%cap%rows )then
                yk = top + real((k - 1) * gr%cap%row_h, dp)
                call gtkui_draw_text(cr, layout, gr%font, gr%cap%label_px, short, x + gr%cap%pad, yk, 'nw', GTKUI_COLOR_MUTED)
                call gtkui_text_size(layout, gr%font, gr%cap%label_px, short, lw, th)
                txt = gtkui_fit_text(layout, gr%font, gr%cap%value_px, value, gr%tile - 2 * gr%cap%pad - lw - gr%cap%gap)
                call gtkui_draw_text(cr, layout, gr%font, gr%cap%value_px, txt, x + real(gr%tile - gr%cap%pad, dp), yk, &
                    &'ne', GTKUI_COLOR_FG)
            else
                cx = x + real(gr%cap%pad, dp) + (real(k, dp) - 0.5_dp) * fw
                call gtkui_draw_text(cr, layout, gr%font, gr%cap%label_px, label, cx, top, 'n', GTKUI_COLOR_MUTED)
                txt = gtkui_fit_text(layout, gr%font, gr%cap%value_px, value, int(fw) - 8)
                call gtkui_draw_text(cr, layout, gr%font, gr%cap%value_px, txt, cx, top + real(gr%cap%row_h, dp), 'n', &
                    &GTKUI_COLOR_FG)
            endif
        end do
    end subroutine draw_caption

    ! ---- GTK handlers (user data: the grid) --------------------------------

    function on_draw( widget, cr, data ) bind(c) result( ret )
        type(c_ptr), value :: widget, cr, data
        integer(c_int)     :: ret
        type(gtkui_image_grid), pointer :: gr
        type(c_ptr) :: layout
        real(dp)    :: offset, x, y
        integer     :: nitems, pitch_y, r0, r1, pos, row, col
        call c_f_pointer(data, gr)
        ret = FALSE
        if( gr%card_h == 0 ) call update_layout(gr)
        call gtkui_set_rgb(cr, GTKUI_COLOR_BG)
        call cairo_paint(cr)
        layout = gtk_widget_create_pango_layout(widget, c_null_char)
        nitems = gr%dlg%nitems()
        if( nitems == 0 )then
            call gtkui_draw_text(cr, layout, gr%font, GTKUI_MESSAGE_PX, 'Nothing to show', &
                &0.5_dp * real(gr%width, dp), 0.5_dp * real(gr%height, dp), 'c', GTKUI_COLOR_MUTED)
            call g_object_unref(layout)
            return
        endif
        offset  = anint(real(gtk_adjustment_get_value(gr%vadj), dp))
        pitch_y = gr%card_h + GTKUI_GRID_SPACING
        r0 = max(0, floor((offset - real(GTKUI_GRID_MARGIN, dp)) / real(pitch_y, dp)))
        r1 = floor((offset + real(gr%height - GTKUI_GRID_MARGIN, dp)) / real(pitch_y, dp))
        do pos = r0 * gr%cols, min(nitems, (r1 + 1) * gr%cols) - 1
            row = pos / gr%cols
            col = mod(pos, gr%cols)
            x   = real(GTKUI_GRID_MARGIN + col * (gr%tile + GTKUI_GRID_SPACING), dp)
            y   = real(GTKUI_GRID_MARGIN + row * pitch_y, dp) - offset
            call draw_card(gr, cr, layout, gr%dlg%item_at(pos + 1), x, y)
        end do
        call g_object_unref(layout)
    end function on_draw

    subroutine on_allocate( widget, allocation, data ) bind(c)
        type(c_ptr), value :: widget, allocation, data
        type(gtkui_image_grid), pointer :: gr
        call c_f_pointer(data, gr)
        gr%width  = int(gtk_widget_get_allocated_width(widget))
        gr%height = int(gtk_widget_get_allocated_height(widget))
        call set_tile(gr, gr%tile, .true.)
    end subroutine on_allocate

    subroutine on_vadj( adj, data ) bind(c)
        type(c_ptr), value :: adj, data
        type(gtkui_image_grid), pointer :: gr
        call c_f_pointer(data, gr)
        call gtk_widget_queue_draw(gr%area)
    end subroutine on_vadj

    ! a left click goes to the delegate's click, a context-menu press to its context (while an attached popover
    ! is up, on_popover_press gets them instead)
    function on_press( widget, event, data ) bind(c) result( handled )
        type(c_ptr), value :: widget, event, data
        integer(c_int)     :: handled
        type(gtkui_image_grid), pointer :: gr
        real(dp) :: x, y
        integer  :: pos
        call c_f_pointer(data, gr)
        handled = TRUE
        call gtk_widget_grab_focus(widget)
        if( .not. gtkui_event_xy(event, x, y) ) return
        if( gtkui_event_context_menu(event) )then
            call gr%dlg%context(x, y)
            return
        endif
        if( .not. gtkui_event_is_press(event) ) return
        if( gtkui_event_button(event) /= 1 ) return
        pos = hit_position(gr, x, y)
        if( pos == 0 ) return
        call gr%dlg%click(gr%dlg%item_at(pos), iand(gtkui_event_state(event), GDK_SHIFT_MASK) /= 0)
        call gtk_widget_queue_draw(gr%area)
    end function on_press

    ! a press on the grid while the attached popover is up: a context-menu press moves the popover there (the
    ! delegate's context), any other press only closes it, the extra events of a double or triple click are
    ! dropped; presses inside the popover go on to its own handlers
    function on_popover_press( popover, event, data ) bind(c) result( handled )
        type(c_ptr), value :: popover, event, data
        integer(c_int)     :: handled
        type(gtkui_image_grid), pointer :: gr
        real(dp) :: x, y
        call c_f_pointer(data, gr)
        handled = FALSE
        if( .not. c_associated(gtk_get_event_widget(event), gr%area) ) return
        handled = TRUE
        if( .not. gtkui_event_is_press(event) ) return
        if( gtkui_event_context_menu(event) )then
            if( gtkui_event_xy(event, x, y) ) call gr%dlg%context(x, y)
        else
            call gtk_popover_popdown(popover)
        endif
    end function on_popover_press

    ! the release of a grid click taken by on_popover_press, kept from the popover's own release handler, which
    ! would close a popover that has just been moved
    function on_popover_release( popover, event, data ) bind(c) result( handled )
        type(c_ptr), value :: popover, event, data
        integer(c_int)     :: handled
        type(gtkui_image_grid), pointer :: gr
        call c_f_pointer(data, gr)
        handled = gtkui_bool(c_associated(gtk_get_event_widget(event), gr%area))
    end function on_popover_release

    function on_motion( widget, event, data ) bind(c) result( handled )
        type(c_ptr), value :: widget, event, data
        integer(c_int)     :: handled
        type(gtkui_image_grid), pointer :: gr
        real(dp) :: x, y
        integer  :: pos, item
        call c_f_pointer(data, gr)
        handled = FALSE
        if( .not. gtkui_event_xy(event, x, y) ) return
        pos  = hit_position(gr, x, y)
        item = 0
        if( pos > 0 ) item = gr%dlg%item_at(pos)
        if( item /= gr%hover )then
            gr%hover = item
            call gtk_widget_queue_draw(gr%area)
        endif
    end function on_motion

    function on_leave( widget, event, data ) bind(c) result( handled )
        type(c_ptr), value :: widget, event, data
        integer(c_int)     :: handled
        type(gtkui_image_grid), pointer :: gr
        call c_f_pointer(data, gr)
        handled = FALSE
        if( gr%hover /= 0 )then
            gr%hover = 0
            call gtk_widget_queue_draw(gr%area)
        endif
    end function on_leave

    ! Ctrl (Cmd on macOS) with the wheel zooms by GTKUI_ZOOM_FACTOR per step; the wheel alone scrolls
    function on_scroll( widget, event, data ) bind(c) result( handled )
        type(c_ptr), value :: widget, event, data
        integer(c_int)     :: handled
        type(gtkui_image_grid), pointer :: gr
        integer(c_int) :: zoom_mask
        real(dp) :: dy, page, value
        integer  :: steps
        call c_f_pointer(data, gr)
        handled = TRUE
        if( .not. gtkui_event_scroll(event, dy) ) return
        zoom_mask = GDK_CONTROL_MASK
#ifdef MACOSX
        zoom_mask = ior(zoom_mask, ior(GDK_META_MASK, GDK_MOD2_MASK))
#endif
        if( iand(gtkui_event_state(event), zoom_mask) /= 0 )then
            gr%zoom_acc = gr%zoom_acc + dy
            steps       = int(gr%zoom_acc)
            if( steps /= 0 )then
                gr%zoom_acc = gr%zoom_acc - real(steps, dp)
                call set_tile(gr, int(real(gr%tile, dp) * GTKUI_ZOOM_FACTOR**(-steps)), .false.)
            endif
        else
            page  = real(gtk_adjustment_get_page_size(gr%vadj), dp)
            value = real(gtk_adjustment_get_value(gr%vadj), dp) + dy * page**(2._dp / 3._dp)
            value = min(max(value, 0._dp), max(0._dp, real(gtk_adjustment_get_upper(gr%vadj), dp) - page))
            call gtk_adjustment_set_value(gr%vadj, real(value, c_double))
        endif
    end function on_scroll

end module simple_gtk_image_grid
