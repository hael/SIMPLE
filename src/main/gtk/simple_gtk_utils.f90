!@descr: stateless GTK 3 helpers of SIMPLE's GTK code: display check, application run, margins, event readers, choice dialog
! The only place besides the widgets that calls gtk-fortran for glue. Event readers return Fortran values
! (the GDK accessors take C pointers); gtkui_choice_dialog stands in for the variadic GTK dialog constructors.
module simple_gtk_utils
use, intrinsic :: iso_c_binding, only: c_ptr, c_funptr, c_int, c_double, c_null_char, c_null_ptr, c_loc
use simple_defs, only: dp
use gtk,         only: gtk_init_check, gtk_application_new, gtk_widget_set_margin_start, gtk_widget_set_margin_end, &
    gtk_widget_set_margin_top, gtk_widget_set_margin_bottom, gtk_widget_get_window, gtk_widget_show_all, &
    gtk_widget_destroy, gtk_window_set_title, gtk_window_set_transient_for, gtk_window_set_modal, gtk_dialog_new, &
    gtk_dialog_add_button, gtk_dialog_get_content_area, gtk_dialog_run, gtk_container_add, gtk_label_new, &
    g_signal_connect, TRUE, FALSE, G_APPLICATION_NON_UNIQUE, GDK_BUTTON_PRESS, GDK_SCROLL_UP, GDK_SCROLL_DOWN
use g,           only: g_application_run, g_object_unref
use gdk,         only: gdk_event_get_coords, gdk_event_get_state, gdk_event_get_button, gdk_event_get_scroll_deltas, &
    gdk_event_get_scroll_direction, gdk_event_get_keyval, gdk_event_get_event_type, gdk_event_triggers_context_menu, &
    gdk_window_set_cursor
implicit none

public :: gtkui_display_available, gtkui_run_app, gtkui_bool, gtkui_margins, gtkui_set_cursor, gtkui_choice_dialog
public :: gtkui_event_xy, gtkui_event_state, gtkui_event_button, gtkui_event_keyval, gtkui_event_scroll
public :: gtkui_event_is_press, gtkui_event_context_menu
private

contains

    ! true when GTK can open a window here (a desktop session or ssh -X)
    logical function gtkui_display_available()
        gtkui_display_available = gtk_init_check(c_null_ptr, [c_null_ptr]) /= FALSE
    end function gtkui_display_available

    ! runs a GtkApplication whose 'activate' handler (activate, called with data) builds the windows;
    ! returns GApplication's exit status once the last window has closed
    integer function gtkui_run_app( app_id, activate, data ) result( status )
        character(len=*), intent(in) :: app_id
        type(c_funptr),   intent(in) :: activate
        type(c_ptr),      intent(in) :: data
        type(c_ptr) :: app
        app = gtk_application_new(app_id//c_null_char, G_APPLICATION_NON_UNIQUE)
        call g_signal_connect(app, 'activate'//c_null_char, activate, data)
        status = int(g_application_run(app, 0_c_int, [c_null_ptr]))
        call g_object_unref(app)
    end function gtkui_run_app

    integer(c_int) function gtkui_bool( l )
        logical, intent(in) :: l
        gtkui_bool = FALSE
        if( l ) gtkui_bool = TRUE
    end function gtkui_bool

    subroutine gtkui_margins( widget, mstart, mend, mtop, mbottom )
        type(c_ptr), intent(in) :: widget
        integer,     intent(in) :: mstart, mend, mtop, mbottom
        call gtk_widget_set_margin_start(widget,  int(mstart,  c_int))
        call gtk_widget_set_margin_end(widget,    int(mend,    c_int))
        call gtk_widget_set_margin_top(widget,    int(mtop,    c_int))
        call gtk_widget_set_margin_bottom(widget, int(mbottom, c_int))
    end subroutine gtkui_margins

    ! the pointer over a realized widget takes cursor (c_null_ptr: the default)
    subroutine gtkui_set_cursor( widget, cursor )
        type(c_ptr), intent(in) :: widget, cursor
        call gdk_window_set_cursor(gtk_widget_get_window(widget), cursor)
    end subroutine gtkui_set_cursor

    ! a modal message with one button per entry of buttons; the chosen index, 0 when closed otherwise
    integer function gtkui_choice_dialog( parent, title, text, buttons ) result( ichoice )
        type(c_ptr),      intent(in) :: parent
        character(len=*), intent(in) :: title, text, buttons(:)
        type(c_ptr)    :: dialog, label, button
        integer(c_int) :: resp
        integer        :: i
        dialog = gtk_dialog_new()
        call gtk_window_set_title(dialog, trim(title)//c_null_char)
        call gtk_window_set_transient_for(dialog, parent)
        call gtk_window_set_modal(dialog, TRUE)
        do i = 1,size(buttons)
            button = gtk_dialog_add_button(dialog, trim(buttons(i))//c_null_char, int(i, c_int))
        end do
        label = gtk_label_new(trim(text)//c_null_char)
        call gtkui_margins(label, 16, 16, 16, 16)
        call gtk_container_add(gtk_dialog_get_content_area(dialog), label)
        call gtk_widget_show_all(dialog)
        resp = gtk_dialog_run(dialog)
        call gtk_widget_destroy(dialog)
        ichoice = 0
        if( resp >= 1_c_int .and. resp <= int(size(buttons), c_int) ) ichoice = int(resp)
    end function gtkui_choice_dialog

    ! the event's position in the widget's coordinates; false when it has none
    logical function gtkui_event_xy( event, x, y )
        type(c_ptr), intent(in)  :: event
        real(dp),    intent(out) :: x, y
        real(c_double), target :: ex, ey
        gtkui_event_xy = gdk_event_get_coords(event, c_loc(ex), c_loc(ey)) /= FALSE
        x = 0._dp
        y = 0._dp
        if( gtkui_event_xy )then
            x = real(ex, dp)
            y = real(ey, dp)
        endif
    end function gtkui_event_xy

    ! the modifier state (GdkModifierType bits), 0 when the event has none
    integer(c_int) function gtkui_event_state( event )
        type(c_ptr), intent(in) :: event
        integer(c_int), target :: state
        gtkui_event_state = 0_c_int
        if( gdk_event_get_state(event, c_loc(state)) /= FALSE ) gtkui_event_state = state
    end function gtkui_event_state

    ! the mouse button (1 left), 0 when the event has none
    integer function gtkui_event_button( event )
        type(c_ptr), intent(in) :: event
        integer(c_int), target :: button
        gtkui_event_button = 0
        if( gdk_event_get_button(event, c_loc(button)) /= FALSE ) gtkui_event_button = int(button)
    end function gtkui_event_button

    ! the key value (GDK keysym), 0 when the event has none
    integer(c_int) function gtkui_event_keyval( event )
        type(c_ptr), intent(in) :: event
        integer(c_int), target :: keyval
        gtkui_event_keyval = 0_c_int
        if( gdk_event_get_keyval(event, c_loc(keyval)) /= FALSE ) gtkui_event_keyval = keyval
    end function gtkui_event_keyval

    ! the vertical scroll amount (smooth delta, or -1 up and +1 down for a wheel click); false for none
    logical function gtkui_event_scroll( event, dy )
        type(c_ptr), intent(in)  :: event
        real(dp),    intent(out) :: dy
        real(c_double), target :: sdx, sdy
        integer(c_int), target :: direction
        gtkui_event_scroll = .true.
        dy = 0._dp
        if( gdk_event_get_scroll_deltas(event, c_loc(sdx), c_loc(sdy)) /= FALSE )then
            dy = real(sdy, dp)
            return
        endif
        gtkui_event_scroll = .false.
        if( gdk_event_get_scroll_direction(event, c_loc(direction)) == FALSE ) return
        if( direction == GDK_SCROLL_UP )then
            dy = -1._dp
            gtkui_event_scroll = .true.
        else if( direction == GDK_SCROLL_DOWN )then
            dy = 1._dp
            gtkui_event_scroll = .true.
        endif
    end function gtkui_event_scroll

    ! a single press (not the extra events of a double or triple click)
    logical function gtkui_event_is_press( event )
        type(c_ptr), intent(in) :: event
        gtkui_event_is_press = gdk_event_get_event_type(event) == GDK_BUTTON_PRESS
    end function gtkui_event_is_press

    ! a press that opens a context menu here (right click; Ctrl-click on macOS)
    logical function gtkui_event_context_menu( event )
        type(c_ptr), intent(in) :: event
        gtkui_event_context_menu = gdk_event_triggers_context_menu(event) /= FALSE
    end function gtkui_event_context_menu

end module simple_gtk_utils
