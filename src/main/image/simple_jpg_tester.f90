!@descr: self-contained MRC-to-JPEG conversion tests
module simple_jpg_tester
use, intrinsic :: iso_c_binding, only: c_associated, c_char, c_f_pointer, c_int, c_null_char, c_null_ptr, &
    &c_ptr, c_signed_char
use simple_test_utils, only: assert_int, assert_real, assert_true
use simple_string,     only: string
use simple_syslib,     only: del_file, file_exists
use simple_image,      only: image
implicit none
private
public :: run_all_mrc2jpeg_tests

integer,          parameter :: WIDTH    = 32, HEIGHT = 24, NSECTIONS = 3
real,             parameter :: SMPD     = 1.5
character(len=*), parameter :: MRC_FILE = 'tmp_mrc2jpeg_input.mrc'

interface
    function stbi_load(filename, width, height, components, desired_components) result(result) &
        &bind(C, name='stbi_load')
        import :: c_char, c_int, c_ptr
        character(c_char), intent(in)        :: filename(*)
        integer(c_int),    intent(out)       :: width, height, components
        integer(c_int),    value, intent(in) :: desired_components
        type(c_ptr) :: result
    end function stbi_load

    subroutine stbi_image_free(buffer) bind(C, name='stbi_image_free')
        import :: c_ptr
        type(c_ptr), value, intent(in) :: buffer
    end subroutine stbi_image_free
end interface

contains

    subroutine run_all_mrc2jpeg_tests()
        type(image) :: img
        integer     :: section, output_count
        call cleanup_files
        call write_source_stack
        if( .not. file_exists(string(MRC_FILE)) ) return
        call img%new([WIDTH, HEIGHT, 1], SMPD, wthreads=.false.)
        do section = 1, NSECTIONS
            call img%read(string(MRC_FILE), section)
            call img%write_jpg(string(jpeg_name(section)), quality=100)
        end do
        call img%kill
        output_count = 0
        do section = 1, NSECTIONS
            if( file_exists(string(jpeg_name(section))) ) output_count = output_count + 1
        end do
        call assert_int(NSECTIONS, output_count, 'mrc2jpeg writes one JPEG for every MRC section')
        do section = 1, NSECTIONS
            call verify_jpeg(section)
        end do
        call cleanup_files
        call test_framed_tiles
    end subroutine run_all_mrc2jpeg_tests

    !> mrc2jpeg_tiled with frames: four 32x32 images on a 2x2 grid of JPEG_DIM tiles, frames red, none, orange,
    !! none; decoded, the frames' middle rows carry their colours and an unframed tile's edge stays grey
    subroutine test_framed_tiles()
        use simple_defs,      only: JPEG_DIM
        use simple_gui_utils, only: mrc2jpeg_tiled
        character(len=*), parameter :: STK = 'tmp_mrc2jpeg_framed.mrc', JPG = 'tmp_mrc2jpeg_framed.jpg'
        type(image) :: img
        character(kind=c_char, len=:), allocatable :: c_filename
        integer(c_signed_char),        pointer     :: pixels(:)
        type(c_ptr)    :: buffer
        integer(c_int) :: w, h, comps
        integer        :: k, x, y, rgb_red(3), rgb_orange(3), rgb_none(3)
        call img%new([32, 32, 1], SMPD, wthreads=.false.)
        do k = 1,4
            do y = 1,32
                do x = 1,32
                    call img%set([x, y, 1], real(mod(x + y + k, 7)))
                end do
            end do
            call img%write(string(STK), k)
        end do
        call img%kill
        call mrc2jpeg_tiled(string(STK), string(JPG), frames=[1, 0, 2, 0])
        call assert_true(file_exists(string(JPG)), 'mrc2jpeg_tiled with frames writes its JPEG')
        if( .not. file_exists(string(JPG)) ) return
        c_filename = JPG//c_null_char
        buffer     = stbi_load(c_filename, w, h, comps, 3_c_int)
        call assert_true(c_associated(buffer), 'the framed JPEG decodes')
        if( .not. c_associated(buffer) ) return
        call assert_int(2 * JPEG_DIM, int(w), 'the framed JPEG holds two tiles across')
        call assert_int(2 * JPEG_DIM, int(h), 'and two down')
        call c_f_pointer(buffer, pixels, [3 * int(w) * int(h)])
        rgb_red    = pixel(JPEG_DIM / 2, 2)               ! tile 1, the top frame's middle row
        rgb_orange = pixel(JPEG_DIM / 2, JPEG_DIM + 2)    ! tile 3 (column 1, row 2)
        rgb_none   = pixel(JPEG_DIM + JPEG_DIM / 2, 2)    ! tile 2, no frame
        call assert_true(rgb_red(1) > 180 .and. rgb_red(2) < 90 .and. rgb_red(3) < 90, 'a red frame on tile 1')
        call assert_true(rgb_orange(1) > 180 .and. rgb_orange(2) > 100 .and. rgb_orange(2) < 210 .and. &
            &rgb_orange(3) < 90, 'an orange frame on tile 3')
        call assert_true(maxval(rgb_none) - minval(rgb_none) < 30, 'no colour at the edge of the unframed tile 2')
        call stbi_image_free(buffer)
        nullify(pixels)
        call del_file(STK)
        call del_file(JPG)

      contains

        function pixel( px, py ) result( rgb )
            integer, intent(in) :: px, py
            integer :: rgb(3), c, v
            do c = 1,3
                v = int(pixels(3 * ((py - 1) * int(w) + px - 1) + c))
                if( v < 0 ) v = v + 256
                rgb(c) = v
            end do
        end function pixel

    end subroutine test_framed_tiles

    subroutine write_source_stack()
        type(image) :: img
        integer     :: section, x, y
        call img%new([WIDTH, HEIGHT, 1], SMPD, wthreads=.false.)
        do section = 1, NSECTIONS
            do y = 1, HEIGHT
                do x = 1, WIDTH
                    call img%set([x, y, 1], source_value(section, x, y))
                end do
            end do
            call img%write(string(MRC_FILE), section)
        end do
        call img%kill
        call assert_true(file_exists(string(MRC_FILE)), 'mrc2jpeg source stack is created internally')
    end subroutine write_source_stack

    subroutine verify_jpeg(section)
        integer, intent(in) :: section
        character(kind=c_char, len=:), allocatable :: c_filename
        integer(c_signed_char),        pointer     :: pixels(:)
        type(c_ptr)                                :: buffer
        integer(c_int)                             :: width_read, height_read, components
        integer                                    :: x, y, channel, idx, actual, expected, max_error, max_channel_delta, rgb(3)
        real                                       :: mean_abs_error
        character(len=:),              allocatable :: filename
        filename = jpeg_name(section)
        call assert_true(file_exists(string(filename)), 'mrc2jpeg writes section '//section_name(section))
        if( .not. file_exists(string(filename)) ) return

        c_filename = filename//c_null_char
        buffer      = c_null_ptr
        buffer      = stbi_load(c_filename, width_read, height_read, components, 3_c_int)
        call assert_true(c_associated(buffer), 'mrc2jpeg output decodes as JPEG for section '//section_name(section))
        if( .not. c_associated(buffer) ) return

        call assert_int(WIDTH,  int(width_read),  'mrc2jpeg JPEG width for section '//section_name(section))
        call assert_int(HEIGHT, int(height_read), 'mrc2jpeg JPEG height for section '//section_name(section))
        call assert_int(3, int(components), 'mrc2jpeg JPEG has three encoded components for section '//section_name(section))
        if( width_read /= WIDTH .or. height_read /= HEIGHT )then
            call stbi_image_free(buffer)
            return
        endif
        call c_f_pointer(buffer, pixels, [3 * int(width_read) * int(height_read)])
        mean_abs_error = 0.
        max_error      = 0
        max_channel_delta = 0
        do y = 1, HEIGHT
            do x = 1, WIDTH
                expected = int(255. * source_value(section, x, y))
                do channel = 1, 3
                    idx    = 3 * ((y - 1) * WIDTH + x - 1) + channel
                    actual = int(pixels(idx))
                    if( actual < 0 ) actual = actual + 256
                    rgb(channel) = actual
                    mean_abs_error = mean_abs_error + real(abs(actual - expected))
                    max_error      = max(max_error, abs(actual - expected))
                end do
                max_channel_delta = max(max_channel_delta, maxval(rgb) - minval(rgb))
            end do
        end do
        mean_abs_error = mean_abs_error / real(3 * WIDTH * HEIGHT)
        call assert_true(max_channel_delta <= 2, &
            &'mrc2jpeg decoded RGB channels differ by no more than 2 for section '//section_name(section))
        call assert_true(max_error <= 5, 'mrc2jpeg maximum decoded pixel error <= 5 for section '//section_name(section))
        call assert_real(0., mean_abs_error, 2., &
            &'mrc2jpeg mean absolute decoded pixel error <= 2 for section '//section_name(section))
        call stbi_image_free(buffer)
        nullify(pixels)
    end subroutine verify_jpeg

    pure real function source_value(section, x, y) result(value)
        integer, intent(in) :: section, x, y
        select case(section)
            case(1)
                value = real(x - 1) / real(WIDTH - 1)
            case(2)
                value = real(y - 1) / real(HEIGHT - 1)
            case default
                value = 0.5 * (real(x - 1) / real(WIDTH - 1) + real(y - 1) / real(HEIGHT - 1))
        end select
    end function source_value

    function jpeg_name(section) result(filename)
        integer, intent(in) :: section
        character(len=28) :: filename
        write(filename, '(A,I3.3,A)') 'tmp_mrc2jpeg_output_', section, '.jpeg'
    end function jpeg_name

    function section_name(section) result(name)
        integer, intent(in) :: section
        character(len=1) :: name
        write(name, '(I1)') section
    end function section_name

    subroutine cleanup_files()
        integer :: section
        call del_file(MRC_FILE)
        do section = 1, NSECTIONS
            call del_file(jpeg_name(section))
        end do
    end subroutine cleanup_files

end module simple_jpg_tester
