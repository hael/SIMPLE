!@descr: GTK manual particle picker: movable boxes, box-size slider, pan/zoom, display controls and tracking-box output
module simple_pick_positions_window
use, intrinsic :: iso_c_binding, only: c_ptr, c_int, c_double, c_null_ptr, c_null_char, c_funloc, c_loc, &
    c_f_pointer, c_associated
use simple_defs,              only: dp, logfhandle
use simple_defs_gtk,          only: GTKUI_COLOR_BACKDROP, GTKUI_ZOOM_FACTOR
use simple_manual_picker,     only: manual_picker
use simple_stack_display,     only: stack_display
use simple_string_utils,      only: int2str, real2str_sig
use simple_gtk_utils,         only: gtkui_run_app, gtkui_margins, gtkui_event_xy, gtkui_event_button, &
    gtkui_event_scroll, gtkui_event_keyval, gtkui_choice_dialog, gtkui_event_is_press
use simple_gtk_draw,          only: gtkui_gray_surface, gtkui_set_rgb
use simple_gtk_display_panel, only: gtkui_display_panel
use gtk, only: gtk_application_window_new, gtk_window_set_title, gtk_window_set_default_size, gtk_container_add, &
    gtk_widget_show_all, gtk_box_new, gtk_box_pack_start, gtk_box_pack_end, gtk_label_new, gtk_label_set_text, &
    gtk_button_new_with_label, gtk_drawing_area_new, gtk_widget_set_can_focus, gtk_widget_add_events, &
    gtk_label_set_max_width_chars, gtk_label_set_ellipsize, gtk_widget_set_tooltip_text, PANGO_ELLIPSIZE_MIDDLE, &
    gtk_widget_queue_draw, gtk_widget_grab_focus, gtk_widget_get_allocated_width, gtk_widget_get_allocated_height, &
    gtk_scale_new, gtk_scale_set_draw_value, gtk_spin_button_new, gtk_adjustment_new, gtk_adjustment_get_value, &
    gtk_widget_set_hexpand, g_signal_connect, TRUE, FALSE, GTK_ORIENTATION_HORIZONTAL, GTK_ORIENTATION_VERTICAL, &
    GDK_BUTTON_PRESS_MASK, GDK_BUTTON_RELEASE_MASK, GDK_POINTER_MOTION_MASK, GDK_SCROLL_MASK, GDK_SMOOTH_SCROLL_MASK, &
    GDK_KEY_PRESS_MASK, CAIRO_FILTER_NEAREST
use cairo, only: cairo_surface_destroy, cairo_save, cairo_restore, cairo_translate, cairo_scale, &
    cairo_set_source_surface, cairo_paint, cairo_pattern_set_filter, cairo_get_source, cairo_rectangle, &
    cairo_set_source_rgb, cairo_stroke, cairo_set_line_width, cairo_fill
implicit none
private
public :: run_pick_positions_window

type :: picking_window
    type(manual_picker), pointer :: picks => null()
    type(stack_display), pointer :: disp => null()
    type(gtkui_display_panel) :: panel
    type(c_ptr) :: window = c_null_ptr, area = c_null_ptr, status = c_null_ptr, surface = c_null_ptr
    character(len=:), allocatable :: outfile, title
    integer :: revision = -1
    logical :: fitted = .false., dragging = .false., panning = .false.
    real(dp) :: last_pointer(2) = 0._dp, drag_offset(2) = 0._dp
end type picking_window

contains

    subroutine run_pick_positions_window( picks, disp, filename, frame, outfile )
        type(manual_picker), target, intent(inout) :: picks
        type(stack_display), target, intent(inout) :: disp
        character(len=*), intent(in) :: filename, outfile
        integer, intent(in) :: frame
        type(picking_window), target :: w
        integer :: status
        w%picks => picks
        w%disp => disp
        w%outfile = outfile
        w%title = 'Pick tracking positions - '//trim(filename)//' (frame '//int2str(frame)//')'
        write(logfhandle,'(a)') '>>> PICKING '//trim(filename)//' frame '//int2str(frame)
        status = gtkui_run_app('org.simple.pick_positions', c_funloc(on_activate), c_loc(w))
        call w%panel%kill
        if( c_associated(w%surface) ) call cairo_surface_destroy(w%surface)
        if( status /= 0 ) write(logfhandle,'(a,i0)') 'pick_positions: GTK application status ', status
    end subroutine run_pick_positions_window

    subroutine on_activate( app, data ) bind(c, name='pick_positions_on_activate')
        type(c_ptr), value :: app, data
        type(picking_window), pointer :: w
        type(c_ptr) :: root, toolbar, controls, label, button, adj, scale, spin
        integer :: dims(2)
        call c_f_pointer(data, w)
        dims = w%disp%get_box()
        w%window = gtk_application_window_new(app)
        call gtk_window_set_title(w%window, w%title//c_null_char)
        call gtk_window_set_default_size(w%window, 1000_c_int, 850_c_int)
        call g_signal_connect(w%window, 'delete-event'//c_null_char, c_funloc(on_close), data)
        root = gtk_box_new(GTK_ORIENTATION_VERTICAL, 8_c_int)
        toolbar = gtk_box_new(GTK_ORIENTATION_HORIZONTAL, 8_c_int)
        call gtkui_margins(toolbar, 12, 12, 8, 0)
        label = gtk_label_new('Box size (pixels)'//c_null_char)
        call gtk_box_pack_start(toolbar, label, FALSE, FALSE, 0_c_int)
        adj = gtk_adjustment_new(real(w%picks%get_box(),c_double), 1._c_double, &
            &real(minval(dims),c_double), 1._c_double, 10._c_double, 0._c_double)
        scale = gtk_scale_new(GTK_ORIENTATION_HORIZONTAL, adj)
        call gtk_scale_set_draw_value(scale, FALSE)
        call gtk_widget_set_hexpand(scale, TRUE)
        spin = gtk_spin_button_new(adj, 1._c_double, 0_c_int)
        call g_signal_connect(adj, 'value-changed'//c_null_char, c_funloc(on_box_size), data)
        call gtk_box_pack_start(toolbar, scale, TRUE, TRUE, 0_c_int)
        call gtk_box_pack_start(toolbar, spin, FALSE, FALSE, 0_c_int)
        button = gtk_button_new_with_label('Delete selected'//c_null_char)
        call g_signal_connect(button, 'clicked'//c_null_char, c_funloc(on_remove), data)
        call gtk_box_pack_start(toolbar, button, FALSE, FALSE, 0_c_int)
        button = gtk_button_new_with_label('Fit image'//c_null_char)
        call g_signal_connect(button, 'clicked'//c_null_char, c_funloc(on_fit), data)
        call gtk_box_pack_start(toolbar, button, FALSE, FALSE, 0_c_int)
        button = gtk_button_new_with_label('Save positions'//c_null_char)
        call g_signal_connect(button, 'clicked'//c_null_char, c_funloc(on_save), data)
        call gtk_box_pack_end(toolbar, button, FALSE, FALSE, 0_c_int)
        call gtk_box_pack_start(root, toolbar, FALSE, FALSE, 0_c_int)
        label = gtk_label_new('Left click/drag: add or move | Right drag: pan | Wheel: zoom | Delete: remove'//c_null_char)
        call gtk_box_pack_start(root, label, FALSE, FALSE, 0_c_int)
        w%area = gtk_drawing_area_new()
        call gtk_widget_set_can_focus(w%area, TRUE)
        call gtk_widget_add_events(w%area, ior(ior(ior(GDK_BUTTON_PRESS_MASK, GDK_BUTTON_RELEASE_MASK), &
            &ior(GDK_POINTER_MOTION_MASK, GDK_KEY_PRESS_MASK)), ior(GDK_SCROLL_MASK, GDK_SMOOTH_SCROLL_MASK)))
        call g_signal_connect(w%area, 'draw'//c_null_char, c_funloc(on_draw), data)
        call g_signal_connect(w%area, 'button-press-event'//c_null_char, c_funloc(on_press), data)
        call g_signal_connect(w%area, 'button-release-event'//c_null_char, c_funloc(on_release), data)
        call g_signal_connect(w%area, 'motion-notify-event'//c_null_char, c_funloc(on_motion), data)
        call g_signal_connect(w%area, 'scroll-event'//c_null_char, c_funloc(on_scroll), data)
        call g_signal_connect(w%area, 'key-press-event'//c_null_char, c_funloc(on_key), data)
        call gtk_box_pack_start(root, w%area, TRUE, TRUE, 0_c_int)
        controls = gtk_box_new(GTK_ORIENTATION_HORIZONTAL, 12_c_int)
        call gtkui_margins(controls, 12, 12, 0, 8)
        call w%panel%new(w%disp, w%area)
        call gtk_box_pack_start(controls, w%panel%get_widget(), TRUE, TRUE, 0_c_int)
        w%status = gtk_label_new(c_null_char)
        call gtk_label_set_max_width_chars(w%status, 40_c_int)
        call gtk_label_set_ellipsize(w%status, PANGO_ELLIPSIZE_MIDDLE)
        call gtk_widget_set_tooltip_text(w%status, w%outfile//c_null_char)
        call gtk_box_pack_end(controls, w%status, FALSE, FALSE, 0_c_int)
        call gtk_box_pack_start(root, controls, FALSE, FALSE, 0_c_int)
        call gtk_container_add(w%window, root)
        call refresh(w)
        call gtk_widget_show_all(w%window)
        call gtk_widget_grab_focus(w%area)
    end subroutine on_activate

    function on_draw( widget, cr, data ) bind(c, name='pick_positions_on_draw') result( handled )
        type(c_ptr), value :: widget, cr, data
        integer(c_int) :: handled
        type(picking_window), pointer :: w
        integer, allocatable :: gray(:,:)
        real(dp), allocatable :: xy(:,:)
        real(dp) :: top(2), bottom(2), corner(2), other(2), factor, half
        integer :: dims(2), i
        call c_f_pointer(data, w)
        if( .not. w%fitted )then
            call w%picks%fit(int(gtk_widget_get_allocated_width(widget)), int(gtk_widget_get_allocated_height(widget)))
            w%fitted = .true.
        endif
        dims = w%disp%get_box()
        if( w%revision /= w%disp%get_revision() )then
            if( c_associated(w%surface) ) call cairo_surface_destroy(w%surface)
            allocate(gray(dims(1),dims(2)))
            call w%disp%render(1, gray)
            ! Cairo rows run down; MRC rows and box-file y run up.
            w%surface = gtkui_gray_surface(gray(:,dims(2):1:-1))
            w%revision = w%disp%get_revision()
        endif
        call gtkui_set_rgb(cr, GTKUI_COLOR_BACKDROP)
        call cairo_paint(cr)
        top = w%picks%image_to_screen([0._dp, real(dims(2),dp)])
        bottom = w%picks%image_to_screen([real(dims(1),dp), 0._dp])
        factor = (bottom(1)-top(1))/real(dims(1),dp)
        call cairo_save(cr)
        call cairo_translate(cr, top(1), top(2))
        call cairo_scale(cr, factor, factor)
        call cairo_set_source_surface(cr, w%surface, 0._c_double, 0._c_double)
        call cairo_pattern_set_filter(cairo_get_source(cr), CAIRO_FILTER_NEAREST)
        call cairo_rectangle(cr, 0._c_double, 0._c_double, real(dims(1),c_double), real(dims(2),c_double))
        call cairo_fill(cr)
        call cairo_restore(cr)
        xy = w%picks%get_centres()
        half = real(w%picks%get_box(),dp)/2._dp
        do i = 1,size(xy,2)
            corner = w%picks%image_to_screen(xy(:,i)+[-half,half])
            other = w%picks%image_to_screen(xy(:,i)+[half,-half])
            if( i == w%picks%get_selected() )then
                call cairo_set_source_rgb(cr, 1._c_double, 0._c_double, 0._c_double)
                call cairo_set_line_width(cr, 2.2_c_double)
            else
                call cairo_set_source_rgb(cr, 1._c_double, 1._c_double, 0._c_double)
                call cairo_set_line_width(cr, 1.5_c_double)
            endif
            call cairo_rectangle(cr, corner(1), corner(2), other(1)-corner(1), other(2)-corner(2))
            call cairo_stroke(cr)
        end do
        handled = FALSE
    end function on_draw

    subroutine refresh( w )
        type(picking_window), intent(inout) :: w
        real(dp), allocatable :: xy(:,:)
        character(len=:), allocatable :: text
        integer :: i
        text = int2str(w%picks%count_boxes())//' picks'
        i = w%picks%get_selected()
        if( i > 0 )then
            xy = w%picks%get_centres()
            text = text//' | Selected '//int2str(i)//': '//real2str_sig(xy(1,i),5)//', '//real2str_sig(xy(2,i),5)
        endif
        call gtk_label_set_text(w%status, text//c_null_char)
        call gtk_widget_queue_draw(w%area)
    end subroutine refresh

    function on_press( widget, event, data ) bind(c, name='pick_positions_on_press') result( handled )
        type(c_ptr), value :: widget, event, data
        integer(c_int) :: handled
        type(picking_window), pointer :: w
        real(dp) :: screen(2), xy(2)
        real(dp), allocatable :: centres(:,:)
        integer :: dims(2)
        call c_f_pointer(data, w)
        handled = FALSE
        if( .not. gtkui_event_is_press(event) ) return
        if( .not. gtkui_event_xy(event, screen(1), screen(2)) ) return
        call gtk_widget_grab_focus(widget)
        select case( gtkui_event_button(event) )
            case(1)
                xy = w%picks%screen_to_image(screen)
                dims = w%disp%get_box()
                if( any(xy < 0._dp) .or. any(xy > real(dims,dp)) ) return
                call w%picks%click(xy)
                if( w%picks%get_selected() > 0 )then
                    centres = w%picks%get_centres()
                    w%drag_offset = centres(:,w%picks%get_selected())-xy
                    w%dragging = .true.
                endif
            case(3)
                w%last_pointer = screen
                w%panning = .true.
        end select
        call refresh(w)
        handled = TRUE
    end function on_press

    function on_release( widget, event, data ) bind(c, name='pick_positions_on_release') result( handled )
        type(c_ptr), value :: widget, event, data
        integer(c_int) :: handled
        type(picking_window), pointer :: w
        call c_f_pointer(data, w)
        if( gtkui_event_button(event) == 1 ) w%dragging = .false.
        if( gtkui_event_button(event) == 3 ) w%panning = .false.
        handled = TRUE
    end function on_release

    function on_motion( widget, event, data ) bind(c, name='pick_positions_on_motion') result( handled )
        type(c_ptr), value :: widget, event, data
        integer(c_int) :: handled
        type(picking_window), pointer :: w
        real(dp) :: screen(2)
        call c_f_pointer(data, w)
        handled = FALSE
        if( .not. gtkui_event_xy(event, screen(1), screen(2)) ) return
        if( w%panning )then
            call w%picks%pan(screen-w%last_pointer)
            w%last_pointer = screen
        else if( w%dragging )then
            call w%picks%move_selected(w%picks%screen_to_image(screen)+w%drag_offset)
        else
            return
        endif
        call refresh(w)
        handled = TRUE
    end function on_motion

    function on_scroll( widget, event, data ) bind(c, name='pick_positions_on_scroll') result( handled )
        type(c_ptr), value :: widget, event, data
        integer(c_int) :: handled
        type(picking_window), pointer :: w
        real(dp) :: screen(2), dy
        call c_f_pointer(data, w)
        handled = FALSE
        if( .not. gtkui_event_xy(event, screen(1), screen(2)) ) return
        if( .not. gtkui_event_scroll(event, dy) ) return
        call w%picks%zoom(screen, GTKUI_ZOOM_FACTOR**(-max(-10._dp,min(10._dp,dy))))
        call refresh(w)
        handled = TRUE
    end function on_scroll

    function on_key( widget, event, data ) bind(c, name='pick_positions_on_key') result( handled )
        type(c_ptr), value :: widget, event, data
        integer(c_int) :: handled
        type(picking_window), pointer :: w
        call c_f_pointer(data, w)
        handled = FALSE
        select case(gtkui_event_keyval(event))
            case(int(z'FFFF',c_int), int(z'FF08',c_int))
                call w%picks%delete_selected()
                call refresh(w)
                handled = TRUE
        end select
    end function on_key

    subroutine on_box_size( adjustment, data ) bind(c, name='pick_positions_on_box_size')
        type(c_ptr), value :: adjustment, data
        type(picking_window), pointer :: w
        call c_f_pointer(data, w)
        w%dragging = .false.
        call w%picks%set_box(nint(gtk_adjustment_get_value(adjustment)))
        call refresh(w)
    end subroutine on_box_size

    subroutine on_remove( button, data ) bind(c, name='pick_positions_on_remove')
        type(c_ptr), value :: button, data
        type(picking_window), pointer :: w
        call c_f_pointer(data, w)
        call w%picks%delete_selected()
        call refresh(w)
    end subroutine on_remove

    subroutine on_fit( button, data ) bind(c, name='pick_positions_on_fit')
        type(c_ptr), value :: button, data
        type(picking_window), pointer :: w
        call c_f_pointer(data, w)
        w%fitted = .false.
        call gtk_widget_queue_draw(w%area)
    end subroutine on_fit

    logical function save_positions( w ) result( saved )
        type(picking_window), intent(inout) :: w
        integer :: status, choice
        character(len=1024) :: message
        call w%picks%save(w%outfile, status, message)
        saved = status == 0
        if( saved )then
            write(logfhandle,'(a)') '>>> SAVED '//int2str(w%picks%count_boxes())//' positions to '//w%outfile
            call gtk_label_set_text(w%status, 'Saved '//int2str(w%picks%count_boxes())//' picks to '//w%outfile//c_null_char)
        else
            choice = gtkui_choice_dialog(w%window, 'Could not save positions', trim(message), &
                &[character(len=2) :: 'OK'])
        endif
    end function save_positions

    subroutine on_save( button, data ) bind(c, name='pick_positions_on_save')
        type(c_ptr), value :: button, data
        type(picking_window), pointer :: w
        logical :: saved
        call c_f_pointer(data, w)
        saved = save_positions(w)
    end subroutine on_save

    function on_close( widget, event, data ) bind(c, name='pick_positions_on_close') result( keep_open )
        type(c_ptr), value :: widget, event, data
        integer(c_int) :: keep_open
        type(picking_window), pointer :: w
        call c_f_pointer(data, w)
        keep_open = TRUE
        if( save_positions(w) ) keep_open = FALSE
    end function on_close

end module simple_pick_positions_window
