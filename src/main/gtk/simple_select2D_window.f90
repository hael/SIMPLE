!@descr: select2D's GTK 3 window: the class-average grid, sorting and selection, the display panel and saving the selection
! An app on the shared GTK modules (doc/implementation_notes/planned/select2D.md, decisions 10-12): an image
! grid fed by a delegate over the cls2D_selection and the stack_display, and a popover with sorting,
! select all/none and the display panel.
module simple_select2D_window
use, intrinsic :: iso_c_binding, only: c_ptr, c_int, c_float, c_null_char, c_null_ptr, c_funloc, c_loc, &
    c_f_pointer, c_associated
use simple_defs,              only: dp, logfhandle
use simple_defs_gtk,          only: GTKUI_GRID_MARGIN, GTKUI_U_MIDDOT, GTKUI_U_UP, GTKUI_U_DOWN, GTKUI_U_RIGHT, GTKUI_U_CMD, &
    GTKUI_KEY_F10, GTKUI_KEY_MENU
use simple_string,            only: string
use simple_string_utils,      only: int2str, int2str_sep, real2str_sig
use simple_stack_display,     only: stack_display
use simple_cls2D_selection,   only: cls2D_selection
use simple_gtk_utils,         only: gtkui_run_app, gtkui_bool, gtkui_margins, gtkui_choice_dialog, gtkui_event_keyval, &
    gtkui_event_state
use simple_gtk_image_grid,    only: gtkui_image_grid, gtkui_grid_delegate
use simple_gtk_display_panel, only: gtkui_display_panel
use gtk, only: gtk_application_window_new, gtk_window_set_title, gtk_window_set_default_size, gtk_container_add, &
    gtk_widget_show_all, gtk_widget_set_size_request, gtk_widget_set_sensitive, gtk_widget_grab_focus, gtk_box_new, &
    gtk_box_pack_start, gtk_box_pack_end, gtk_separator_new, gtk_label_new, gtk_label_set_text, gtk_label_set_xalign, &
    gtk_label_set_ellipsize, gtk_button_new_with_label, gtk_grid_new, gtk_grid_attach, gtk_grid_set_row_spacing, &
    gtk_radio_button_new_with_label_from_widget, gtk_toggle_button_get_active, gtk_toggle_button_set_active, &
    gtk_popover_new, gtk_popover_set_pointing_to, gtk_popover_popup, gtk_popover_popdown, g_signal_connect, TRUE, FALSE, &
    GTK_ORIENTATION_HORIZONTAL, GTK_ORIENTATION_VERTICAL, PANGO_ELLIPSIZE_MIDDLE, GDK_SHIFT_MASK
use gdk_events, only: GdkRectangle
implicit none

public :: run_select2D_window
private

integer, parameter :: FIELD_POP = 1, FIELD_RES = 2    ! caption fields

! the grid's view of the selection session and the display
type, extends(gtkui_grid_delegate) :: s2d_delegate
    type(s2d_window), pointer :: w => null()
  contains
    procedure :: nitems      => s2d_nitems
    procedure :: nids        => s2d_nids
    procedure :: item_at     => s2d_item_at
    procedure :: is_selected => s2d_is_selected
    procedure :: image_box   => s2d_image_box
    procedure :: revision    => s2d_revision
    procedure :: render      => s2d_render
    procedure :: ncaptions   => s2d_ncaptions
    procedure :: caption     => s2d_caption
    procedure :: click       => s2d_click
    procedure :: context     => s2d_context
end type s2d_delegate

! one window run, shared with the GTK handlers through their user data
type :: s2d_window
    type(cls2D_selection), pointer :: sel  => null()
    type(stack_display),   pointer :: disp => null()
    character(len=:), allocatable  :: projname, projfile
    type(s2d_delegate)        :: dlg
    type(gtkui_image_grid)    :: grid
    type(gtkui_display_panel) :: panel
    type(c_ptr) :: window = c_null_ptr, path_label = c_null_ptr, footer = c_null_ptr, popover = c_null_ptr
    type(c_ptr) :: radio(3) = c_null_ptr, btn_select = c_null_ptr, btn_deselect = c_null_ptr
    integer     :: nfields = 0, fields(2) = 0
    logical     :: syncing = .false.     ! radio buttons set from the model, not by the user
end type s2d_window

contains

    ! runs the window until it closes; the selection is written only by Save (or on closing, when asked)
    subroutine run_select2D_window( sel, disp, projname, projfile )
        type(cls2D_selection), target, intent(inout) :: sel
        type(stack_display),   target, intent(inout) :: disp
        character(len=*),              intent(in)    :: projname, projfile
        type(s2d_window), target :: w
        integer :: status
        w%sel      => sel
        w%disp     => disp
        w%dlg%w    => w
        w%projname = projname
        w%projfile = projfile
        if( sel%has_field('pop') )then
            w%nfields           = w%nfields + 1
            w%fields(w%nfields) = FIELD_POP
        endif
        if( sel%has_field('res') )then
            w%nfields           = w%nfields + 1
            w%fields(w%nfields) = FIELD_RES
        endif
        status = gtkui_run_app('org.simple.select2D', c_funloc(on_activate), c_loc(w))
        call w%grid%kill
        call w%panel%kill
        if( status /= 0 ) write(logfhandle,'(a,i0)') 'select2D: the GTK application ended with status ', status
    end subroutine run_select2D_window

    ! ---- the window --------------------------------------------------------

    subroutine on_activate( app, data ) bind(c)
        type(c_ptr), value :: app, data
        type(s2d_window), pointer :: w
        type(c_ptr) :: vbox, bar, button, foot, hint
        character(len=:), allocatable :: modifier
        call c_f_pointer(data, w)
        w%window = gtk_application_window_new(app)
        call gtk_window_set_default_size(w%window, 800_c_int, 600_c_int)
        call gtk_widget_set_size_request(w%window, 650_c_int, 350_c_int)
        call g_signal_connect(w%window, 'delete-event'//c_null_char,    c_funloc(on_delete),     data)
        call g_signal_connect(w%window, 'destroy'//c_null_char,         c_funloc(on_destroy),    data)
        call g_signal_connect(w%window, 'key-press-event'//c_null_char, c_funloc(on_window_key), data)
        vbox = gtk_box_new(GTK_ORIENTATION_VERTICAL, 0_c_int)
        ! toolbar: the project (and, once saved, where the selection went) and Save
        bar = gtk_box_new(GTK_ORIENTATION_HORIZONTAL, 8_c_int)
        call gtkui_margins(bar, 12, 12, 12, 11)
        w%path_label = gtk_label_new(w%projfile//c_null_char)
        call gtk_label_set_xalign(w%path_label, 0._c_float)
        call gtk_label_set_ellipsize(w%path_label, PANGO_ELLIPSIZE_MIDDLE)
        call gtk_box_pack_start(bar, w%path_label, TRUE, TRUE, 0_c_int)
        button = gtk_button_new_with_label('Save selection'//c_null_char)
        call g_signal_connect(button, 'clicked'//c_null_char, c_funloc(on_save), data)
        call gtk_box_pack_end(bar, button, FALSE, FALSE, 0_c_int)
        call gtk_box_pack_start(vbox, bar, FALSE, FALSE, 0_c_int)
        call gtk_box_pack_start(vbox, gtk_separator_new(GTK_ORIENTATION_HORIZONTAL), FALSE, FALSE, 0_c_int)
        ! the class averages
        call w%grid%new(w%dlg)
        call gtk_box_pack_start(vbox, w%grid%get_widget(), TRUE, TRUE, 0_c_int)
        call gtk_box_pack_start(vbox, gtk_separator_new(GTK_ORIENTATION_HORIZONTAL), FALSE, FALSE, 0_c_int)
        ! footer: gestures and counts
#ifdef MACOSX
        modifier = GTKUI_U_CMD
#else
        modifier = 'Ctrl'
#endif
        foot = gtk_box_new(GTK_ORIENTATION_HORIZONTAL, 8_c_int)
        call gtkui_margins(foot, 12, 12, 10, 10)
        hint = gtk_label_new('Click to select  '//GTKUI_U_MIDDOT//'  Shift-click for range  '//GTKUI_U_MIDDOT//'  '//&
            &modifier//' + scroll to zoom  '//GTKUI_U_MIDDOT//'  right-click for display and sorting'//c_null_char)
        call gtk_label_set_xalign(hint, 0._c_float)
        call gtk_label_set_ellipsize(hint, PANGO_ELLIPSIZE_MIDDLE)
        call gtk_box_pack_start(foot, hint, TRUE, TRUE, 0_c_int)
        w%footer = gtk_label_new(c_null_char)
        call gtk_box_pack_end(foot, w%footer, FALSE, FALSE, 0_c_int)
        call gtk_box_pack_start(vbox, foot, FALSE, FALSE, 0_c_int)
        call gtk_container_add(w%window, vbox)
        call build_popover(w, data)
        call update_status(w)
        call gtk_widget_show_all(w%window)
        call gtk_widget_grab_focus(w%grid%get_area())
    end subroutine on_activate

    ! the popover: sort, select all/none and the display panel
    subroutine build_popover( w, data )
        type(s2d_window), target, intent(inout) :: w
        type(c_ptr),              intent(in)    :: data
        type(c_ptr) :: pg, label
        integer     :: k
        w%popover = gtk_popover_new(w%grid%get_area())
        pg = gtk_grid_new()
        call gtkui_margins(pg, 8, 8, 8, 8)
        call gtk_grid_set_row_spacing(pg, 4_c_int)
        label = gtk_label_new('Sort by'//c_null_char)
        call gtk_label_set_xalign(label, 0._c_float)
        call gtk_grid_attach(pg, label, 0_c_int, 0_c_int, 1_c_int, 1_c_int)
        w%radio(1) = gtk_radio_button_new_with_label_from_widget(c_null_ptr, 'Resolution '//GTKUI_U_UP//c_null_char)
        w%radio(2) = gtk_radio_button_new_with_label_from_widget(w%radio(1), 'Population '//GTKUI_U_DOWN//c_null_char)
        w%radio(3) = gtk_radio_button_new_with_label_from_widget(w%radio(1), 'Original class ID '//GTKUI_U_UP//c_null_char)
        do k = 1,3
            call g_signal_connect(w%radio(k), 'toggled'//c_null_char, c_funloc(on_sort_toggled), data)
            call gtk_grid_attach(pg, w%radio(k), 0_c_int, int(k, c_int), 1_c_int, 1_c_int)
        end do
        w%btn_select = gtk_button_new_with_label('Select All'//c_null_char)
        call g_signal_connect(w%btn_select, 'clicked'//c_null_char, c_funloc(on_select_all), data)
        call gtk_grid_attach(pg, w%btn_select, 0_c_int, 4_c_int, 1_c_int, 1_c_int)
        w%btn_deselect = gtk_button_new_with_label('Deselect All'//c_null_char)
        call g_signal_connect(w%btn_deselect, 'clicked'//c_null_char, c_funloc(on_deselect_all), data)
        call gtk_grid_attach(pg, w%btn_deselect, 0_c_int, 5_c_int, 1_c_int, 1_c_int)
        call gtk_grid_attach(pg, gtk_separator_new(GTK_ORIENTATION_HORIZONTAL), 0_c_int, 6_c_int, 1_c_int, 1_c_int)
        call w%panel%new(w%disp, w%grid%get_area())
        call gtk_grid_attach(pg, w%panel%get_widget(), 0_c_int, 7_c_int, 1_c_int, 1_c_int)
        call gtk_container_add(w%popover, pg)
        call gtk_widget_show_all(pg)
    end subroutine build_popover

    ! the popover at (x, y) of the grid, its controls set from the models
    subroutine show_popover( w, x, y )
        type(s2d_window), intent(inout) :: w
        real(dp),         intent(in)    :: x, y
        type(GdkRectangle), target :: rect
        integer(kind=8) :: psel, ptot
        integer :: nsel, nvis, k
        call w%sel%counts(nsel, nvis, psel, ptot)
        do k = 1,3
            call gtk_widget_set_sensitive(w%radio(k), gtkui_bool(nvis > 0))
        end do
        call gtk_widget_set_sensitive(w%btn_select,   gtkui_bool(nvis > 0))
        call gtk_widget_set_sensitive(w%btn_deselect, gtkui_bool(nsel > 0))
        w%syncing = .true.
        call gtk_toggle_button_set_active(w%radio(w%sel%get_sort()), TRUE)
        w%syncing = .false.
        call w%panel%sync()
        rect%x      = int(x, c_int)
        rect%y      = int(y, c_int)
        rect%width  = 1_c_int
        rect%height = 1_c_int
        call gtk_popover_set_pointing_to(w%popover, c_loc(rect))
        call gtk_popover_popup(w%popover)
    end subroutine show_popover

    ! footer counts and the window title (a leading * while the selection is not saved)
    subroutine update_status( w )
        type(s2d_window), intent(inout) :: w
        integer(kind=8) :: psel, ptot
        integer :: nsel, nvis
        character(len=:), allocatable :: title
        call w%sel%counts(nsel, nvis, psel, ptot)
        call gtk_label_set_text(w%footer, 'Selected: '//int2str(nsel)//'/'//int2str(nvis)//' classes  '//&
            &int2str_sep(psel)//'/'//int2str_sep(ptot)//' particles'//c_null_char)
        title = 'Select 2D - '//w%projname
        if( w%sel%is_dirty() ) title = '*'//title
        call gtk_window_set_title(w%window, title//c_null_char)
    end subroutine update_status

    ! ---- window handlers (user data: the s2d_window) -----------------------

    function on_delete( widget, event, data ) bind(c) result( keep_open )
        type(c_ptr), value :: widget, event, data
        integer(c_int)     :: keep_open
        type(s2d_window), pointer :: w
        type(string)              :: saved
        call c_f_pointer(data, w)
        keep_open = FALSE
        if( .not. w%sel%is_dirty() ) return
        select case( gtkui_choice_dialog(w%window, 'Unsaved selection', &
            &'The selection has changes that are not saved.', &
            &[character(len=8) :: 'Save', 'Discard', 'Cancel']) )
            case( 1 )
                call w%sel%save()
            case( 2 )
                saved = w%sel%get_saved_projfile()
                if( saved%strlen_trim() == 0 )then
                    write(logfhandle,'(a)') '>>> SELECTION DISCARDED; nothing was written'
                else
                    write(logfhandle,'(a)') '>>> SELECTION DISCARDED; '//saved%to_char()//' keeps the last saved selection'
                endif
            case DEFAULT
                keep_open = TRUE
        end select
    end function on_delete

    subroutine on_destroy( widget, data ) bind(c)
        type(c_ptr), value :: widget, data
        type(s2d_window), pointer :: w
        call c_f_pointer(data, w)
        w%window = c_null_ptr
    end subroutine on_destroy

    ! Shift+F10 or the Menu key opens the popover at the grid's corner
    function on_window_key( widget, event, data ) bind(c) result( handled )
        type(c_ptr), value :: widget, event, data
        integer(c_int)     :: handled
        type(s2d_window), pointer :: w
        integer(c_int) :: keyval
        call c_f_pointer(data, w)
        handled = FALSE
        keyval  = gtkui_event_keyval(event)
        if( keyval == GTKUI_KEY_MENU .or. (keyval == GTKUI_KEY_F10 .and. iand(gtkui_event_state(event), GDK_SHIFT_MASK) /= 0) )then
            call show_popover(w, real(GTKUI_GRID_MARGIN, dp), real(GTKUI_GRID_MARGIN, dp))
            handled = TRUE
        endif
    end function on_window_key

    subroutine on_save( button, data ) bind(c)
        type(c_ptr), value :: button, data
        type(s2d_window), pointer :: w
        type(string)              :: saved
        call c_f_pointer(data, w)
        call w%sel%save()
        saved = w%sel%get_saved_projfile()
        call gtk_label_set_text(w%path_label, w%projfile//'  '//GTKUI_U_RIGHT//'  '//saved%to_char()//c_null_char)
        call update_status(w)
    end subroutine on_save

    subroutine on_sort_toggled( button, data ) bind(c)
        type(c_ptr), value :: button, data
        type(s2d_window), pointer :: w
        integer :: k
        call c_f_pointer(data, w)
        if( w%syncing ) return
        if( gtk_toggle_button_get_active(button) == FALSE ) return
        do k = 1,3
            if( c_associated(button, w%radio(k)) )then
                call w%sel%sort(k)
                exit
            endif
        end do
        call gtk_popover_popdown(w%popover)
        call w%grid%redraw()
    end subroutine on_sort_toggled

    subroutine on_select_all( button, data ) bind(c)
        type(c_ptr), value :: button, data
        type(s2d_window), pointer :: w
        call c_f_pointer(data, w)
        call gtk_popover_popdown(w%popover)
        call w%sel%select_all()
        call update_status(w)
        call w%grid%redraw()
    end subroutine on_select_all

    subroutine on_deselect_all( button, data ) bind(c)
        type(c_ptr), value :: button, data
        type(s2d_window), pointer :: w
        call c_f_pointer(data, w)
        call gtk_popover_popdown(w%popover)
        call w%sel%deselect_all()
        call update_status(w)
        call w%grid%redraw()
    end subroutine on_deselect_all

    ! ---- the grid delegate -------------------------------------------------

    integer function s2d_nitems( self )
        class(s2d_delegate), intent(in) :: self
        s2d_nitems = self%w%sel%get_nvisible()
    end function s2d_nitems

    integer function s2d_nids( self )
        class(s2d_delegate), intent(in) :: self
        s2d_nids = self%w%sel%get_ncls()
    end function s2d_nids

    integer function s2d_item_at( self, pos )
        class(s2d_delegate), intent(in) :: self
        integer,             intent(in) :: pos
        s2d_item_at = self%w%sel%get_class(pos)
    end function s2d_item_at

    logical function s2d_is_selected( self, item )
        class(s2d_delegate), intent(in) :: self
        integer,             intent(in) :: item
        s2d_is_selected = self%w%sel%is_selected(item)
    end function s2d_is_selected

    function s2d_image_box( self ) result( box )
        class(s2d_delegate), intent(in) :: self
        integer :: box(2)
        box = self%w%disp%get_box()
    end function s2d_image_box

    integer function s2d_revision( self )
        class(s2d_delegate), intent(in) :: self
        s2d_revision = self%w%disp%get_revision()
    end function s2d_revision

    subroutine s2d_render( self, item, gray )
        class(s2d_delegate), intent(in)  :: self
        integer,             intent(in)  :: item
        integer,             intent(out) :: gray(:,:)
        call self%w%disp%render(item, gray)
    end subroutine s2d_render

    integer function s2d_ncaptions( self )
        class(s2d_delegate), intent(in) :: self
        s2d_ncaptions = self%w%nfields
    end function s2d_ncaptions

    ! population and resolution; for item 0 the widest values the captions should fit
    subroutine s2d_caption( self, k, item, label, short, value )
        class(s2d_delegate),           intent(in)  :: self
        integer,                       intent(in)  :: k, item
        character(len=:), allocatable, intent(out) :: label, short, value
        real(dp) :: pop
        value = '-'
        select case( self%w%fields(k) )
            case( FIELD_POP )
                label = 'Population'
                short = 'Pop'
                if( item == 0 )then
                    value = '999,999'
                else if( self%w%sel%has_pop(item) )then
                    pop = self%w%sel%get_pop(item)
                    if( pop == anint(pop) .and. abs(pop) < 1.e15_dp )then
                        value = int2str_sep(nint(pop, kind=8))
                    else
                        value = real2str_sig(pop, 3)
                    endif
                endif
            case DEFAULT
                label = 'Resolution'
                short = 'Res'
                if( item == 0 )then
                    value = '99.9'
                else if( self%w%sel%has_res(item) )then
                    value = real2str_sig(self%w%sel%get_res(item), 3)
                endif
        end select
    end subroutine s2d_caption

    subroutine s2d_click( self, item, shift )
        class(s2d_delegate), intent(inout) :: self
        integer,             intent(in)    :: item
        logical,             intent(in)    :: shift
        call self%w%sel%click(item, shift)
        call update_status(self%w)
    end subroutine s2d_click

    subroutine s2d_context( self, x, y )
        class(s2d_delegate), intent(inout) :: self
        real(dp),            intent(in)    :: x, y
        call show_popover(self%w, x, y)
    end subroutine s2d_context

end module simple_select2D_window
