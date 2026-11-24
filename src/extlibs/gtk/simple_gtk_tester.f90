!@descr: platform test of GTK 4 through gtk-fortran: an 8-bit tile painted with cairo, and a window that opens, draws and closes
! The tile takes the path select2D will (real image -> 8-bit -> cairo surface) and is checked
! offscreen, without a display; the window needs one (otherwise it is skipped) and closes itself
! WINDOW_MS after its first draw.
module simple_gtk_tester
use, intrinsic :: iso_c_binding, only: c_ptr, c_int, c_int32_t, c_double, c_null_char, c_null_ptr, &
    c_funloc, c_loc, c_f_pointer, c_associated
use gtk,               only: gtk_init_check, gtk_application_new, gtk_application_window_new, &
    gtk_window_set_title, gtk_window_set_default_size, gtk_window_set_child, gtk_window_present, &
    gtk_window_destroy, gtk_box_new, gtk_box_append, gtk_label_new, gtk_button_new_with_label, &
    gtk_drawing_area_new, gtk_drawing_area_set_content_width, gtk_drawing_area_set_content_height, &
    gtk_drawing_area_set_draw_func, gtk_widget_set_hexpand, gtk_widget_set_vexpand, &
    gtk_widget_set_margin_start, gtk_widget_set_margin_end, gtk_widget_set_margin_top, &
    gtk_widget_set_margin_bottom, gtk_get_major_version, gtk_get_minor_version, &
    gtk_get_micro_version, g_signal_connect, TRUE, FALSE, G_APPLICATION_NON_UNIQUE, &
    GTK_ORIENTATION_VERTICAL, CAIRO_FORMAT_RGB24, CAIRO_FILTER_GOOD, CAIRO_STATUS_SUCCESS
use g,                 only: g_application_run, g_object_unref, g_timeout_add, g_source_remove
use cairo,             only: cairo_image_surface_create, cairo_image_surface_get_data, &
    cairo_image_surface_get_stride, cairo_surface_flush, cairo_surface_mark_dirty, &
    cairo_surface_destroy, cairo_create, cairo_destroy, cairo_status, cairo_set_source_surface, &
    cairo_set_source_rgb, cairo_paint, cairo_scale, cairo_translate, cairo_save, cairo_restore, &
    cairo_pattern_set_filter, cairo_get_source
use simple_test_utils, only: assert_true, assert_int
use simple_defs,       only: logfhandle
implicit none

public :: run_all_gtk_tile_tests, run_all_gtk_window_tests
private

integer(c_int), parameter :: TILE_BOX  = 128_c_int  ! side of the synthetic tile, pixels
integer(c_int), parameter :: MARGIN    = 12_c_int   ! window padding, pixels
integer(c_int), parameter :: WINDOW_MS = 1500_c_int ! the window closes itself this long after its first draw
integer(c_int), parameter :: PAINT_W   = 300_c_int  ! offscreen canvas of the paint check, pixels
integer(c_int), parameter :: PAINT_H   = 200_c_int

! one window run, shared with the GTK callbacks through their user data
type :: gtk_run
    type(c_ptr)    :: window     = c_null_ptr
    type(c_ptr)    :: tile       = c_null_ptr
    integer(c_int) :: timeout_id = 0_c_int
    integer        :: ndraws     = 0
    integer        :: width      = 0
    integer        :: height     = 0
    logical        :: activated  = .false.
end type gtk_run

contains

    subroutine run_all_gtk_tile_tests()
        write(logfhandle,'(a)') '**** running all gtk tile tests ****'
        call test_tile_pixels
        call test_tile_paint
    end subroutine run_all_gtk_tile_tests

    subroutine run_all_gtk_window_tests()
        write(logfhandle,'(a)') '**** running all gtk window tests ****'
        call test_window
    end subroutine run_all_gtk_window_tests

    ! ---- tests -------------------------------------------------------------

    ! the tile is gray, spans 0-255 and keeps the model: the blob brightest, the ring at 1/1.5 of it
    subroutine test_tile_pixels()
        type(c_ptr)                 :: tile
        integer(c_int32_t), pointer :: pix(:,:)
        integer :: gray(TILE_BOX,TILE_BOX), loc(2), i, j
        logical :: l_gray
        real    :: centre
        write(logfhandle,'(a)') 'test_tile_pixels'
        tile = make_tile()
        call surface_pixels(tile, TILE_BOX, pix)
        l_gray = .true.
        do j = 1,TILE_BOX
            do i = 1,TILE_BOX
                gray(i,j) = int(ibits(pix(i,j), 0, 8))
                if( int(ibits(pix(i,j), 8, 8)) /= gray(i,j) .or. int(ibits(pix(i,j), 16, 8)) /= gray(i,j) ) l_gray = .false.
            end do
        end do
        call assert_true(l_gray, 'every tile pixel is gray (R = G = B)')
        call assert_int(0,   minval(gray), 'the darkest tile pixel is 0')
        call assert_int(255, maxval(gray), 'the brightest tile pixel is 255')
        centre = 0.5 * real(TILE_BOX + 1)
        loc    = maxloc(gray)
        call assert_true(abs(real(loc(1)) - (centre + 0.1 * real(TILE_BOX))) <= 1. .and. abs(real(loc(2)) - centre) <= 1., &
            &'the brightest tile pixel is the blob centre, 0.1 box right of the tile centre')
        ! opposite the blob the ring is 1 against the blob's 1.5; sampling moves either by under 0.5 %
        i = nint(centre - 0.3 * real(TILE_BOX))
        j = nint(centre)
        call assert_true(abs(gray(i,j) - 170) <= 3, 'a ring pixel maps linearly to 255 / 1.5')
        call cairo_surface_destroy(tile)
    end subroutine test_tile_pixels

    ! paint_tile on an offscreen canvas: the tile scaled to the canvas height and centred, background beside it
    subroutine test_tile_paint()
        type(c_ptr)                 :: tile, canvas, cr
        integer(c_int32_t), pointer :: pix(:,:)
        integer :: blue(PAINT_W,PAINT_H), loc(2), i, j
        logical :: l_gray
        real    :: scale, x0
        write(logfhandle,'(a)') 'test_tile_paint'
        tile   = make_tile()
        canvas = cairo_image_surface_create(CAIRO_FORMAT_RGB24, PAINT_W, PAINT_H)
        cr     = cairo_create(canvas)
        call paint_tile(cr, PAINT_W, PAINT_H, tile)
        call assert_int(int(CAIRO_STATUS_SUCCESS), int(cairo_status(cr)), 'cairo paints the tile without error')
        call cairo_destroy(cr)
        call surface_pixels(canvas, PAINT_H, pix)
        scale = real(PAINT_H) / real(TILE_BOX)
        x0    = 0.5 * (real(PAINT_W) - scale * real(TILE_BOX))
        ! background (0.13, 0.13, 0.15) is (33, 33, 38) in 8 bits, left and right of the tile
        call assert_true(is_rgb(pix(11,101), 33, 33, 38) .and. is_rgb(pix(PAINT_W-10,101), 33, 33, 38), &
            &'the background shows beside the scaled tile')
        ! inside the tile, two pixels clear of its edges where the filter meets the background
        l_gray = .true.
        blue   = 0
        do j = 3,PAINT_H-2
            do i = nint(x0)+3,nint(x0+scale*real(TILE_BOX))-2
                blue(i,j) = int(ibits(pix(i,j), 0, 8))
                if( int(ibits(pix(i,j), 8, 8)) /= blue(i,j) .or. int(ibits(pix(i,j), 16, 8)) /= blue(i,j) ) l_gray = .false.
            end do
        end do
        call assert_true(l_gray, 'the scaled tile stays gray')
        ! the blob centre, 76.8 px into the tile, lands at x0 + 76.8 * scale; pixel k covers [k-1, k)
        loc = maxloc(blue)
        call assert_true(abs(real(loc(1)) - 0.5 - (x0 + 76.8 * scale)) <= 2. .and. abs(real(loc(2)) - 0.5 - 64. * scale) <= 2., &
            &'the brightest painted pixel is the scaled blob centre')
        call cairo_surface_destroy(canvas)
        call cairo_surface_destroy(tile)
    end subroutine test_tile_paint

    ! a GtkApplication window opens, draws the tile area at least at its content size and closes itself
    subroutine test_window()
        type(gtk_run), target :: run
        type(c_ptr)           :: app
        integer(c_int)        :: status
        write(logfhandle,'(a)') 'test_window'
        if( gtk_init_check() == FALSE )then
            write(logfhandle,'(a)') 'skipped: no display to open a GTK window on'
            return
        endif
        app = gtk_application_new('org.simple.gtk_tester'//c_null_char, G_APPLICATION_NON_UNIQUE)
        call g_signal_connect(app, 'activate'//c_null_char, c_funloc(activate), c_loc(run))
        status = g_application_run(app, 0_c_int, [c_null_ptr])
        call g_object_unref(app)
        if( c_associated(run%tile) ) call cairo_surface_destroy(run%tile)
        call assert_int(0, int(status), 'the GTK application exits with status 0')
        call assert_true(run%activated, 'the application activates')
        call assert_true(run%ndraws > 0, 'GTK draws the tile area')
        call assert_true(run%width >= TILE_BOX .and. run%height >= TILE_BOX, 'the tile area is drawn at least at its content size')
        call assert_true(.not. c_associated(run%window), 'the window is gone once it closes')
    end subroutine test_window

    ! ---- the tile ----------------------------------------------------------

    ! a synthetic class average (a ring and an off-centre blob) as an RGB24 surface, R = G = B
    function make_tile() result( surface )
        type(c_ptr)                 :: surface
        integer(c_int32_t), pointer :: pix(:,:)
        real    :: img(TILE_BOX,TILE_BOX), x, y, r, box, lo, hi
        integer :: i, j, gray
        box = real(TILE_BOX)
        do j = 1,TILE_BOX
            y = real(j) - 0.5 * (box + 1.)
            do i = 1,TILE_BOX
                x = real(i) - 0.5 * (box + 1.)
                r = sqrt(x * x + y * y)
                img(i,j) = exp(-((r - 0.3 * box) / (0.05 * box))**2) + &
                    &1.5 * exp(-((x - 0.1 * box)**2 + y**2) / (0.08 * box)**2)
            end do
        end do
        lo = minval(img)
        hi = maxval(img)
        surface = cairo_image_surface_create(CAIRO_FORMAT_RGB24, TILE_BOX, TILE_BOX)
        call surface_pixels(surface, TILE_BOX, pix)
        do j = 1,TILE_BOX
            do i = 1,TILE_BOX
                gray     = nint(255. * (img(i,j) - lo) / (hi - lo))
                pix(i,j) = int(gray * 65793, c_int32_t)
            end do
        end do
        call cairo_surface_mark_dirty(surface)
    end function make_tile

    ! the tile scaled to fit width x height and centred on a dark background
    subroutine paint_tile( cr, width, height, tile )
        type(c_ptr),    intent(in) :: cr, tile
        integer(c_int), intent(in) :: width, height
        real(c_double) :: w, h, side, scale
        w     = real(width,    c_double)
        h     = real(height,   c_double)
        side  = real(TILE_BOX, c_double)
        scale = min(w, h) / side
        call cairo_set_source_rgb(cr, 0.13_c_double, 0.13_c_double, 0.15_c_double)
        call cairo_paint(cr)
        call cairo_save(cr)
        call cairo_translate(cr, 0.5_c_double * (w - scale * side), 0.5_c_double * (h - scale * side))
        call cairo_scale(cr, scale, scale)
        call cairo_set_source_surface(cr, tile, 0._c_double, 0._c_double)
        call cairo_pattern_set_filter(cairo_get_source(cr), CAIRO_FILTER_GOOD)
        call cairo_paint(cr)
        call cairo_restore(cr)
    end subroutine paint_tile

    ! the pixels of an RGB24 image surface: rows stride bytes apart, a pixel a native-endian 0x00RRGGBB word
    subroutine surface_pixels( surface, nrows, pix )
        type(c_ptr),                 intent(in)  :: surface
        integer(c_int),              intent(in)  :: nrows
        integer(c_int32_t), pointer, intent(out) :: pix(:,:)
        integer :: stride
        call cairo_surface_flush(surface)
        stride = int(cairo_image_surface_get_stride(surface))
        call c_f_pointer(cairo_image_surface_get_data(surface), pix, [stride / 4, int(nrows)])
    end subroutine surface_pixels

    ! a pixel equal to (r, g, b) up to cairo's rounding of colours to 8 bits
    logical function is_rgb( pixel, r, g, b )
        integer(c_int32_t), intent(in) :: pixel
        integer,            intent(in) :: r, g, b
        is_rgb = abs(int(ibits(pixel, 16, 8)) - r) <= 1 .and. abs(int(ibits(pixel, 8, 8)) - g) <= 1 .and. &
            &abs(int(ibits(pixel, 0, 8)) - b) <= 1
    end function is_rgb

    ! ---- GTK callbacks (user data: the gtk_run of test_window) -------------

    subroutine activate( app, user_data ) bind(c)
        type(c_ptr), value     :: app, user_data
        type(gtk_run), pointer :: run
        type(c_ptr)            :: window, box, label, area, button
        character(len=64)      :: text
        call c_f_pointer(user_data, run)
        run%activated = .true.
        write(text,'(a,i0,a,i0,a,i0)') 'SIMPLE + gtk-fortran on GTK ', gtk_get_major_version(), '.', &
            gtk_get_minor_version(), '.', gtk_get_micro_version()
        write(logfhandle,'(a)') '>>> '//trim(text)
        window     = gtk_application_window_new(app)
        run%window = window
        call g_signal_connect(window, 'destroy'//c_null_char, c_funloc(on_destroy), user_data)
        call gtk_window_set_title(window, 'SIMPLE GTK test'//c_null_char)
        call gtk_window_set_default_size(window, 360_c_int, 420_c_int)
        box = gtk_box_new(GTK_ORIENTATION_VERTICAL, 8_c_int)
        call gtk_widget_set_margin_start(box,  MARGIN)
        call gtk_widget_set_margin_end(box,    MARGIN)
        call gtk_widget_set_margin_top(box,    MARGIN)
        call gtk_widget_set_margin_bottom(box, MARGIN)
        label = gtk_label_new(trim(text)//c_null_char)
        call gtk_box_append(box, label)
        run%tile = make_tile()
        area     = gtk_drawing_area_new()
        call gtk_drawing_area_set_content_width(area,  TILE_BOX)
        call gtk_drawing_area_set_content_height(area, TILE_BOX)
        call gtk_widget_set_hexpand(area, TRUE)
        call gtk_widget_set_vexpand(area, TRUE)
        call gtk_drawing_area_set_draw_func(area, c_funloc(draw_area), user_data, c_funloc(release_tile))
        call gtk_box_append(box, area)
        button = gtk_button_new_with_label('Close'//c_null_char)
        call g_signal_connect(button, 'clicked'//c_null_char, c_funloc(close_window), window)
        call gtk_box_append(box, button)
        call gtk_window_set_child(window, box)
        call gtk_window_present(window)
    end subroutine activate

    ! records the draw; the first one starts the timer that closes the window
    subroutine draw_area( area, cr, width, height, user_data ) bind(c)
        type(c_ptr),    value  :: area, cr, user_data
        integer(c_int), value  :: width, height
        type(gtk_run), pointer :: run
        call c_f_pointer(user_data, run)
        call paint_tile(cr, width, height, run%tile)
        run%ndraws = run%ndraws + 1
        run%width  = int(width)
        run%height = int(height)
        if( run%ndraws == 1 ) run%timeout_id = g_timeout_add(WINDOW_MS, c_funloc(close_on_timeout), user_data)
    end subroutine draw_area

    function close_on_timeout( user_data ) bind(c) result( again )
        type(c_ptr), value     :: user_data
        integer(c_int)         :: again
        type(gtk_run), pointer :: run
        call c_f_pointer(user_data, run)
        run%timeout_id = 0_c_int
        if( c_associated(run%window) ) call gtk_window_destroy(run%window)
        again = FALSE
    end function close_on_timeout

    ! the window is gone (timer or Close): drop a pending timer, it would outlive the run
    subroutine on_destroy( widget, user_data ) bind(c)
        type(c_ptr), value     :: widget, user_data
        type(gtk_run), pointer :: run
        call c_f_pointer(user_data, run)
        run%window = c_null_ptr
        if( run%timeout_id /= 0_c_int )then
            if( g_source_remove(run%timeout_id) == FALSE ) write(logfhandle,'(a)') 'gtk tester: the close timer had already gone'
            run%timeout_id = 0_c_int
        endif
    end subroutine on_destroy

    ! GTK releases the draw function's data with the drawing area
    subroutine release_tile( user_data ) bind(c)
        type(c_ptr), value     :: user_data
        type(gtk_run), pointer :: run
        call c_f_pointer(user_data, run)
        if( c_associated(run%tile) ) call cairo_surface_destroy(run%tile)
        run%tile = c_null_ptr
    end subroutine release_tile

    subroutine close_window( button, window ) bind(c)
        type(c_ptr), value :: button, window
        call gtk_window_destroy(window)
    end subroutine close_window

end module simple_gtk_tester
