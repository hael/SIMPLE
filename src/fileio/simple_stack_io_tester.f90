!@descr: unit test routines for stack_io: buffered contiguous MRC stack reading and writing, float32 and float16
module simple_stack_io_tester
use, intrinsic :: iso_c_binding, only: c_float
use simple_test_utils ! assertions etc.
use simple_string,       only: string
use simple_string_utils, only: int2str
use simple_syslib,       only: del_file, file_exists
use simple_image,        only: image
use simple_stack_io,     only: stack_io
use simple_imghead,      only: MrcImgHead, MRC_MODE_FLOAT32, MRC_MODE_FLOAT16, MRC_NVERSION_20141
implicit none
private
public :: run_all_stack_io_tests

integer, parameter :: BOX   = 16   ! box size of the small synthetic stacks
integer, parameter :: NIMGS = 5    ! odd, so a buffer of two ends in a partial flush/refill
integer, parameter :: BUFSZ = 2
! 1025**2 * 2 elements exceed two float16 write buffers (FLOAT16_WRITE_BUFFER_ELEMS = 1024**2 in
! simple_imgfile) and 1025 does not divide the buffer length, so every row boundary case of the
! chunked float16 writer is exercised: two full flushes and a partial tail
integer, parameter :: CHUNK_BOX   = 1025
integer, parameter :: CHUNK_NIMGS = 2
integer, parameter :: MRC_HEADER_NBYTES = 1024
real,    parameter :: SMPD = 1.3
real,    parameter :: TOL  = 1.0e-7
character(len=*), parameter :: SOURCE_STACK  = 'tmp_stack_io_source.mrc'
character(len=*), parameter :: IMAGE_COPY    = 'tmp_stack_io_image_copy.mrc'
character(len=*), parameter :: STACK_COPY    = 'tmp_stack_io_stack_copy.mrc'
character(len=*), parameter :: FLOAT16_STACK = 'tmp_stack_io_float16.mrc'
character(len=*), parameter :: CHUNK_STACK   = 'tmp_stack_io_float16_chunked.mrc'

contains

    subroutine run_all_stack_io_tests()
        write(*,'(A)') '**** running all stack_io tests ****'
        call test_open_state()
        call test_read_image_written_stack()
        call test_image_reader_to_stack_io_writer()
        call test_stack_io_copy()
        call test_buffer_sizes_and_forward_skips()
        call test_float32_header()
        call test_float16_header_and_roundtrip()
        call test_float16_chunked_roundtrip()
        call cleanup()
    end subroutine run_all_stack_io_tests

    !---------------- open/close bookkeeping ----------------

    subroutine test_open_state()
        type(stack_io) :: writer, reader
        type(image)    :: img
        write(*,'(A)') 'test_open_state'
        call assert_false(writer%stk_is_open(), 'fresh stack_io is not open')
        call writer%open(string(SOURCE_STACK),SMPD,'write',box=BOX,bufsz=BUFSZ)
        call assert_true(writer%stk_is_open(),  'open for write reports open')
        call assert_int(0, writer%get_nptcls(), 'no particle count on write')
        call assert_true(all(writer%get_ldim() == [BOX,BOX,1]), 'ldim set from box on write')
        call img%new([BOX,BOX,1],SMPD,wthreads=.false.)
        call set_pattern(img,1)
        call writer%write(1,img)
        call writer%close
        call assert_false(writer%stk_is_open(), 'closed writer reports closed')
        call assert_true(all(writer%get_ldim() == 0), 'close resets ldim')
        call assert_true(file_exists(string(SOURCE_STACK)), 'close flushes the buffer to disk')
        call reader%open(string(SOURCE_STACK),SMPD,'read',bufsz=BUFSZ)
        call assert_true(reader%stk_is_open(), 'open for read reports open')
        call assert_int(1, reader%get_nptcls(), 'one image written')
        call assert_true(reader%same_stk(string(SOURCE_STACK),[BOX,BOX,1]),      'same_stk: same name and ldim')
        call assert_false(reader%same_stk(string(IMAGE_COPY),[BOX,BOX,1]),      'same_stk: other name')
        call assert_false(reader%same_stk(string(SOURCE_STACK),[BOX+2,BOX,1]),  'same_stk: other ldim')
        call reader%close
        call assert_false(reader%stk_is_open(), 'closed reader reports closed')
        call img%kill
    end subroutine test_open_state

    !---------------- round trips ----------------

    subroutine test_read_image_written_stack()
        write(*,'(A)') 'test_read_image_written_stack'
        call create_synthetic_stack(string(SOURCE_STACK))
        call verify_stack(string(SOURCE_STACK), BUFSZ, 'image-written stack')
    end subroutine test_read_image_written_stack

    subroutine test_image_reader_to_stack_io_writer()
        type(stack_io) :: writer
        type(image)    :: img
        integer :: iimg
        write(*,'(A)') 'test_image_reader_to_stack_io_writer'
        call img%new([BOX,BOX,1],SMPD,wthreads=.false.)
        call writer%open(string(IMAGE_COPY),SMPD,'write',box=BOX,bufsz=BUFSZ)
        do iimg = 1,NIMGS
            call img%read(string(SOURCE_STACK),iimg)
            call writer%write(iimg,img)
        enddo
        call writer%close
        call img%kill
        call verify_stack(string(IMAGE_COPY), BUFSZ, 'image reader to stack_io writer')
    end subroutine test_image_reader_to_stack_io_writer

    subroutine test_stack_io_copy()
        type(stack_io) :: reader, writer
        type(image)    :: img
        integer :: iimg
        write(*,'(A)') 'test_stack_io_copy'
        call reader%open(string(SOURCE_STACK),SMPD,'read',bufsz=BUFSZ)
        call assert_int(NIMGS, reader%get_nptcls(), 'source image count')
        call writer%open(string(STACK_COPY),SMPD,'write',box=BOX,bufsz=BUFSZ)
        call img%new([BOX,BOX,1],SMPD,wthreads=.false.)
        do iimg = 1,NIMGS
            call reader%read(iimg,img)
            call writer%write(iimg,img)
        enddo
        call reader%close
        call writer%close
        call img%kill
        call verify_stack(string(STACK_COPY), BUFSZ, 'stack_io reader to stack_io writer')
    end subroutine test_stack_io_copy

    ! the read buffer is refilled on demand: a buffer of three leaves a partial last window,
    ! a buffer larger than the stack is clamped to it, and forward reads may skip windows
    subroutine test_buffer_sizes_and_forward_skips()
        type(stack_io) :: reader
        type(image)    :: img
        integer :: iimg
        write(*,'(A)') 'test_buffer_sizes_and_forward_skips'
        call verify_stack(string(SOURCE_STACK), 3,         'buffer of three')
        call verify_stack(string(SOURCE_STACK), NIMGS,     'buffer holding the whole stack')
        call verify_stack(string(SOURCE_STACK), NIMGS + 7, 'buffer larger than the stack')
        call img%new([BOX,BOX,1],SMPD,wthreads=.false.)
        ! every other image, buffer of two: 3 skips the rest of window [1,2]
        call reader%open(string(SOURCE_STACK),SMPD,'read',bufsz=BUFSZ)
        do iimg = 1,NIMGS,2
            call reader%read(iimg,img)
            call assert_pattern(img,iimg,'forward skip, buffer of two')
        enddo
        call reader%close
        ! last image first, buffer of one: windows [1],[2],[3],[4] are skipped in one call
        call reader%open(string(SOURCE_STACK),SMPD,'read',bufsz=1)
        call reader%read(NIMGS,img)
        call assert_pattern(img,NIMGS,'jump to the last image, buffer of one')
        call reader%close
        ! get_image re-extracts an image that is already in the buffer
        call reader%open(string(SOURCE_STACK),SMPD,'read',bufsz=BUFSZ)
        call reader%read(2,img)
        call img%zero
        call reader%get_image(1,img)
        call assert_pattern(img,1,'get_image from the current buffer')
        call reader%close
        call img%kill
    end subroutine test_buffer_sizes_and_forward_skips

    !---------------- MRC modes ----------------

    subroutine test_float32_header()
        write(*,'(A)') 'test_float32_header'
        call assert_mrc_header(string(STACK_COPY), MRC_MODE_FLOAT32, 4, 'float32')
    end subroutine test_float32_header

    subroutine test_float16_header_and_roundtrip()
        type(stack_io) :: writer
        type(image)    :: img
        integer :: iimg
        write(*,'(A)') 'test_float16_header_and_roundtrip'
        call writer%open(string(FLOAT16_STACK),SMPD,'write',box=BOX,bufsz=BUFSZ,wfloat16=.true.)
        call img%new([BOX,BOX,1],SMPD,wthreads=.false.)
        do iimg = 1,NIMGS
            call set_pattern(img,iimg)
            call writer%write(iimg,img)
        enddo
        call writer%close
        call img%kill
        call assert_mrc_header(string(FLOAT16_STACK), MRC_MODE_FLOAT16, 2, 'float16')
        ! the pattern values are exactly representable in float16, so the round trip is exact
        call verify_stack(string(FLOAT16_STACK), BUFSZ, 'float16 stack')
    end subroutine test_float16_header_and_roundtrip

    subroutine test_float16_chunked_roundtrip()
        type(stack_io) :: writer, reader
        type(image)    :: source, actual
        real(kind=c_float), pointer :: pixels(:,:,:) => null()
        real(kind=c_float) :: expected(CHUNK_BOX,CHUNK_BOX)
        integer :: iimg, ldim(3)
        write(*,'(A)') 'test_float16_chunked_roundtrip'
        call source%new([CHUNK_BOX,CHUNK_BOX,1],SMPD,wthreads=.false.)
        call actual%new([CHUNK_BOX,CHUNK_BOX,1],SMPD,wthreads=.false.)
        call writer%open(string(CHUNK_STACK),SMPD,'write',box=CHUNK_BOX,bufsz=CHUNK_NIMGS,wfloat16=.true.)
        do iimg = 1,CHUNK_NIMGS
            call source%get_rmat_ptr(pixels)
            call fill_chunk_pattern(pixels(1:CHUNK_BOX,1:CHUNK_BOX,1),iimg)
            call writer%write(iimg,source)
        enddo
        call writer%close
        call assert_mrc_header(string(CHUNK_STACK), MRC_MODE_FLOAT16, 2, 'chunked float16', ldim=[CHUNK_BOX,CHUNK_BOX,CHUNK_NIMGS])
        call reader%open(string(CHUNK_STACK),SMPD,'read',bufsz=1)
        ldim = reader%get_ldim()
        call assert_int(CHUNK_NIMGS, reader%get_nptcls(), 'chunked float16: image count')
        call assert_true(all(ldim == [CHUNK_BOX,CHUNK_BOX,1]), 'chunked float16: dimensions')
        do iimg = 1,CHUNK_NIMGS
            call reader%read(iimg,actual)
            call actual%get_rmat_ptr(pixels)
            call fill_chunk_pattern(expected,iimg)
            call assert_true(maxval(abs(pixels(1:CHUNK_BOX,1:CHUNK_BOX,1)-expected)) <= TOL,&
                &'chunked float16: image '//int2str(iimg)//' survives the buffer flushes exactly')
        enddo
        call reader%close
        call source%kill
        call actual%kill
    end subroutine test_float16_chunked_roundtrip

    !---------------- helpers ----------------

    subroutine create_synthetic_stack(fname)
        class(string), intent(in) :: fname
        type(image) :: img
        integer :: iimg
        call img%new([BOX,BOX,1],SMPD,wthreads=.false.)
        do iimg = 1,NIMGS
            call set_pattern(img,iimg)
            call img%write(fname,iimg,del_if_exists=(iimg == 1))
        enddo
        call img%kill
    end subroutine create_synthetic_stack

    ! reads the whole stack with the given buffer size and checks count, dimensions and every pixel
    subroutine verify_stack(fname, bufsz, label)
        class(string),    intent(in) :: fname
        integer,          intent(in) :: bufsz
        character(len=*), intent(in) :: label
        type(stack_io) :: reader
        type(image)    :: img
        integer :: iimg, ldim(3)
        call reader%open(fname,SMPD,'read',bufsz=bufsz)
        ldim = reader%get_ldim()
        call assert_int(NIMGS, reader%get_nptcls(), label//': image count')
        call assert_true(all(ldim == [BOX,BOX,1]), label//': dimensions')
        call img%new(ldim,SMPD,wthreads=.false.)
        do iimg = 1,NIMGS
            call reader%read(iimg,img)
            call assert_pattern(img,iimg,label)
        enddo
        call reader%close
        call img%kill
    end subroutine verify_stack

    subroutine set_pattern(img,iimg)
        type(image), intent(inout) :: img
        integer,     intent(in)    :: iimg
        real(kind=c_float), pointer :: pixels(:,:,:) => null()
        call img%get_rmat_ptr(pixels)
        call fill_pattern(pixels(1:BOX,1:BOX,1),iimg)
    end subroutine set_pattern

    subroutine assert_pattern(img,iimg,label)
        type(image),      intent(inout) :: img
        integer,          intent(in)    :: iimg
        character(len=*), intent(in)    :: label
        real(kind=c_float), pointer :: pixels(:,:,:) => null()
        real(kind=c_float) :: expected(BOX,BOX)
        call fill_pattern(expected,iimg)
        call img%get_rmat_ptr(pixels)
        call assert_true(maxval(abs(pixels(1:BOX,1:BOX,1)-expected)) <= TOL,&
            &label//': pixel values of image '//int2str(iimg))
    end subroutine assert_pattern

    ! image index as background, first row: signs, halves, a float16 subnormal, the float16 maximum, +-0
    ! (all exactly representable in float16, so the same pattern verifies float32 and float16 stacks)
    subroutine fill_pattern(pixels,iimg)
        real(kind=c_float), intent(out) :: pixels(:,:)
        integer,            intent(in)  :: iimg
        pixels        = real(iimg,kind=c_float)
        pixels(1,1)   = 2.0_c_float
        pixels(2,1)   = -2.0_c_float
        pixels(3,1)   = 1.0_c_float
        pixels(4,1)   = -1.0_c_float
        pixels(5,1)   = 0.5_c_float
        pixels(6,1)   = 0.333251953125_c_float
        pixels(7,1)   = 6.103515625e-5_c_float
        pixels(8,1)   = 65504.0_c_float
        pixels(9,1)   = 0.0_c_float
        pixels(10,1)  = -0.0_c_float
    end subroutine fill_pattern

    ! position-dependent multiples of 1/32 below 64 (exact in float16), offset per image,
    ! so a pixel displaced across a buffer flush is detected
    subroutine fill_chunk_pattern(pixels,iimg)
        real(kind=c_float), intent(out) :: pixels(:,:)
        integer,            intent(in)  :: iimg
        integer :: i, j
        do j = 1,size(pixels,2)
            do i = 1,size(pixels,1)
                pixels(i,j) = real(mod(i + 7*j + 13*iimg, 2048),kind=c_float) / 32.0_c_float
            enddo
        enddo
    end subroutine fill_chunk_pattern

    subroutine assert_mrc_header(fname, mode, bytes_per_pixel, label, ldim)
        class(string),     intent(in) :: fname
        integer,           intent(in) :: mode, bytes_per_pixel
        character(len=*),  intent(in) :: label
        integer, optional, intent(in) :: ldim(3)
        type(MrcImgHead) :: header
        integer :: funit, io_stat, dims(3), expected_dims(3)
        integer(kind=8) :: file_nbytes
        expected_dims = [BOX,BOX,NIMGS]
        if( present(ldim) ) expected_dims = ldim
        call header%new(expected_dims)
        open(newunit=funit,file=fname%to_char(),access='stream',form='unformatted',action='read',status='old',iostat=io_stat)
        call assert_int(0, io_stat, label//' header: file opens')
        if( io_stat /= 0 ) return
        call header%read(funit)
        close(funit)
        dims = header%getDims()
        call assert_int(mode, header%getMode(), label//' header: MRC mode')
        if( mode == MRC_MODE_FLOAT16 )then
            call assert_int(MRC_NVERSION_20141, header%nversion, label//' header: MRC 2014 version stamp')
        endif
        call assert_true(all(dims == expected_dims), label//' header: dimensions')
        inquire(file=fname%to_char(),size=file_nbytes,iostat=io_stat)
        call assert_int(0, io_stat, label//' header: file size readable')
        call assert_true(file_nbytes == int(MRC_HEADER_NBYTES,kind=8) + int(bytes_per_pixel,kind=8)*product(int(expected_dims,kind=8)),&
            &label//' header: file size is header plus '//int2str(bytes_per_pixel)//' bytes per pixel')
        call header%kill
    end subroutine assert_mrc_header

    subroutine cleanup()
        call del_file(string(SOURCE_STACK))
        call del_file(string(IMAGE_COPY))
        call del_file(string(STACK_COPY))
        call del_file(string(FLOAT16_STACK))
        call del_file(string(CHUNK_STACK))
    end subroutine cleanup

end module simple_stack_io_tester
