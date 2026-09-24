!@descr: unit test routines for stack_io and dstack_io: buffered and discrete MRC stack reading and writing, float32, int16 and float16
! The discrete reader (simple_discrete_stack_io) opens a window of stacks and reads them
! concurrently, one OpenMP thread per open stack (three, also in the one-thread gate: the
! concurrent reads are the point). The float16 checks pin the rounding of values that are not
! representable (round half to even), the bit patterns on disk, the image layer's mode
! inheritance and header statistics, and the subnormal and signed-zero boundaries.
module simple_stack_io_tester
use, intrinsic :: iso_c_binding, only: c_float
use, intrinsic :: iso_fortran_env, only: int16, int32
use simple_test_utils ! assertions etc.
use simple_string,       only: string
use simple_string_utils, only: int2str
use simple_syslib,       only: del_file, file_exists
use simple_image,        only: image
use simple_stack_io,     only: stack_io
use simple_discrete_stack_io, only: dstack_io
use simple_imghead,      only: MrcImgHead, MRC_MODE_FLOAT32, MRC_MODE_FLOAT16, MRC_NVERSION_20141, find_ldim_nptcls
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
character(len=*), parameter :: ROUND_STACK   = 'tmp_stack_io_float16_rounding.mrc'
character(len=*), parameter :: F16_IMAGE     = 'tmp_stack_io_float16_image.mrc'
character(len=*), parameter :: F16_BOUNDS    = 'tmp_stack_io_float16_boundaries.mrc'
! discrete reads: twelve stacks of four images, three stacks open and read concurrently
integer, parameter :: DSTK_NSTKS = 12, DSTK_NIMGS = 4, DSTK_WINDOW = 3

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
        call test_float16_rounding_and_payload()
        call test_float16_image_roundtrip()
        call test_float16_encoder_boundaries()
        call test_dstack_parallel_read_float32()
        call test_dstack_parallel_read_int16()
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


    !> values that float16 cannot hold are rounded half to even; the bits on disk are the IEEE half encodings
    subroutine test_float16_rounding_and_payload()
        integer, parameter :: NTEST_IMGS = 5
        type(stack_io) :: writer, reader
        type(image)    :: source, actual, expected_img
        real(kind=c_float), pointer :: actual_rmat(:,:,:) => null(), expected_rmat(:,:,:) => null()
        integer :: iimg
        write(*,'(A)') 'test_float16_rounding_and_payload'
        call source%new([BOX,BOX,1],SMPD,wthreads=.false.)
        call actual%new([BOX,BOX,1],SMPD,wthreads=.false.)
        call expected_img%new([BOX,BOX,1],SMPD,wthreads=.false.)
        call writer%open(string(ROUND_STACK),SMPD,'write',box=BOX,bufsz=2,wfloat16=.true.)
        do iimg = 1,NTEST_IMGS
            call set_f16_rounding_pattern(source,iimg,quantized=.false.)
            call writer%write(iimg,source)
        enddo
        call writer%close
        call assert_mrc_header(string(ROUND_STACK), MRC_MODE_FLOAT16, 2, 'float16 rounding stack', ldim=[BOX,BOX,NTEST_IMGS])
        call assert_f16_payload(string(ROUND_STACK), 'float16 rounding stack')
        call reader%open(string(ROUND_STACK),SMPD,'read',bufsz=3)
        call assert_int(NTEST_IMGS, reader%get_nptcls(), 'float16 rounding stack: image count')
        do iimg = 1,NTEST_IMGS
            call reader%read(iimg,actual)
            call set_f16_rounding_pattern(expected_img,iimg,quantized=.true.)
            call actual%get_rmat_ptr(actual_rmat)
            call expected_img%get_rmat_ptr(expected_rmat)
            call assert_true(maxval(abs(actual_rmat(1:BOX,1:BOX,1)-expected_rmat(1:BOX,1:BOX,1))) <= TOL, &
                &'float16 rounding stack: image '//int2str(iimg)//' reads back rounded half to even')
        enddo
        call reader%close
        call source%kill
        call actual%kill
        call expected_img%kill
    end subroutine test_float16_rounding_and_payload

    !> the image layer: float16 on the first write, inherited by the next, kept by a header-statistics update
    subroutine test_float16_image_roundtrip()
        type(image) :: source, actual, expected_img
        real(kind=c_float), pointer :: actual_rmat(:,:,:) => null(), expected_rmat(:,:,:) => null()
        real    :: stats(4)
        integer :: iimg
        write(*,'(A)') 'test_float16_image_roundtrip'
        call source%new([BOX,BOX,1],SMPD,wthreads=.false.)
        call actual%new([BOX,BOX,1],SMPD,wthreads=.false.)
        call expected_img%new([BOX,BOX,1],SMPD,wthreads=.false.)
        do iimg = 1,2
            call set_f16_rounding_pattern(source,iimg,quantized=.false.)
            if( iimg == 1 )then
                call source%write(string(F16_IMAGE),iimg,del_if_exists=.true.,wfloat16=.true.)
            else
                call source%write(string(F16_IMAGE),iimg)
            endif
        enddo
        do iimg = 1,2
            call actual%read(string(F16_IMAGE),iimg)
            call set_f16_rounding_pattern(expected_img,iimg,quantized=.true.)
            call actual%get_rmat_ptr(actual_rmat)
            call expected_img%get_rmat_ptr(expected_rmat)
            call assert_true(maxval(abs(actual_rmat(1:BOX,1:BOX,1)-expected_rmat(1:BOX,1:BOX,1))) <= TOL, &
                &'float16 image: image '//int2str(iimg)//' round trip, the second write inheriting the mode')
        enddo
        stats = [-1.0,65504.0,0.0,1.0]
        call source%update_header_stats(string(F16_IMAGE),stats)
        call assert_mrc_header(string(F16_IMAGE), MRC_MODE_FLOAT16, 2, 'float16 image after a statistics update', &
            &ldim=[BOX,BOX,2])
        call assert_f16_payload(string(F16_IMAGE), 'float16 image')
        call source%kill
        call actual%kill
        call expected_img%kill
    end subroutine test_float16_image_roundtrip

    !> subnormals are exact, half the smallest subnormal rounds to a signed zero, 3/2 of it to twice it
    subroutine test_float16_encoder_boundaries()
        integer(int32), parameter :: EXPECTED_BITS(8) = [1_int32,32769_int32,256_int32,33024_int32, &
            &0_int32,32768_int32,2_int32,32770_int32]
        type(image) :: source, actual
        real(kind=c_float), pointer :: pixels(:,:,:) => null(), actual_pixels(:,:,:) => null()
        integer(int32) :: bits(8)
        write(*,'(A)') 'test_float16_encoder_boundaries'
        call source%new([BOX,BOX,1],SMPD,wthreads=.false.)
        call actual%new([BOX,BOX,1],SMPD,wthreads=.false.)
        call source%get_rmat_ptr(pixels)
        pixels = 1.0_c_float
        pixels(1,1,1) = 2.0_c_float**(-24)
        pixels(2,1,1) = -2.0_c_float**(-24)
        pixels(3,1,1) = 2.0_c_float**(-16)
        pixels(4,1,1) = -2.0_c_float**(-16)
        pixels(5,1,1) = 2.0_c_float**(-25)
        pixels(6,1,1) = -2.0_c_float**(-25)
        pixels(7,1,1) = 3.0_c_float*2.0_c_float**(-25)
        pixels(8,1,1) = -3.0_c_float*2.0_c_float**(-25)
        call source%write(string(F16_BOUNDS),1,del_if_exists=.true.,wfloat16=.true.)
        call read_f16_first_row(string(F16_BOUNDS), bits)
        call assert_true(all(bits == EXPECTED_BITS), 'float16 boundaries: subnormal, underflow and sign bits on disk')
        call actual%read(string(F16_BOUNDS),1)
        call actual%get_rmat_ptr(actual_pixels)
        pixels(5,1,1) = 0.0_c_float
        pixels(6,1,1) = -0.0_c_float
        pixels(7,1,1) = 2.0_c_float**(-23)
        pixels(8,1,1) = -2.0_c_float**(-23)
        call assert_true(all(actual_pixels(1:BOX,1:BOX,1) == pixels(1:BOX,1:BOX,1)), &
            &'float16 boundaries: zero and subnormal values decode exactly')
        call assert_true(transfer(actual_pixels(6,1,1),0_int32) == not(huge(0_int32)), &
            &'float16 boundaries: negative zero keeps its sign')
        call source%kill
        call actual%kill
    end subroutine test_float16_encoder_boundaries

    !---------------- discrete stack reads ----------------

    !> float32 stacks, three open and read concurrently through dstack_io
    subroutine test_dstack_parallel_read_float32()
        type(string) :: stknames(DSTK_NSTKS)
        type(image)  :: img
        real(kind=c_float), pointer :: rmat(:,:,:) => null()
        integer :: istk, iimg
        write(*,'(A)') 'test_dstack_parallel_read_float32'
        call img%new([BOX,BOX,1], SMPD, wthreads=.false.)
        do istk = 1,DSTK_NSTKS
            stknames(istk) = 'tmp_dstack_io_f32_'//int2str(istk)//'.mrc'
            do iimg = 1,DSTK_NIMGS
                call img%get_rmat_ptr(rmat)
                rmat = real(100*istk + iimg, kind=c_float)
                call img%write(stknames(istk), iimg, del_if_exists=(iimg == 1))
            enddo
        enddo
        call img%kill
        call read_stacks_concurrently(stknames, 'float32 discrete read')
    end subroutine test_dstack_parallel_read_float32

    !> 16-bit integer stacks (MRC mode 1) through the same concurrent reader
    subroutine test_dstack_parallel_read_int16()
        type(string) :: stknames(DSTK_NSTKS)
        integer :: istk
        write(*,'(A)') 'test_dstack_parallel_read_int16'
        do istk = 1,DSTK_NSTKS
            stknames(istk) = 'tmp_dstack_io_i16_'//int2str(istk)//'.mrc'
            call write_int16_stack(stknames(istk), istk)
        enddo
        call read_stacks_concurrently(stknames, 'int16 discrete read')
    end subroutine test_dstack_parallel_read_int16

    !> opens DSTK_WINDOW stacks at a time, reads their images on one thread per stack, checks
    !! every image (constant 100*istk + iimg) and deletes the stacks
    subroutine read_stacks_concurrently( stknames, label )
        type(string),     intent(inout) :: stknames(DSTK_NSTKS)
        character(len=*), intent(in)    :: label
        type(dstack_io) :: dstkios(DSTK_WINDOW)
        type(image)     :: read_imgs(DSTK_NSTKS,DSTK_NIMGS)
        real(kind=c_float), pointer :: rmat(:,:,:) => null()
        integer :: istk, iimg, ldim(3), nptcls, stk_from, stk_to, iopen, nopen, nwrong
        do stk_from = 1,DSTK_NSTKS,DSTK_WINDOW
            stk_to = min(stk_from + DSTK_WINDOW - 1, DSTK_NSTKS)
            nopen  = stk_to - stk_from + 1
            do iopen = 1,nopen
                istk = stk_from + iopen - 1
                call find_ldim_nptcls(stknames(istk), ldim, nptcls)
                call dstkios(iopen)%new(SMPD, BOX)
                call dstkios(iopen)%cache_stack_info(stknames(istk), ldim, nptcls)
                call dstkios(iopen)%open(stknames(istk))
                do iimg = 1,DSTK_NIMGS
                    call read_imgs(istk,iimg)%new([BOX,BOX,1], SMPD, wthreads=.false.)
                enddo
            enddo
            !$omp parallel do default(shared) private(iopen,istk,iimg) schedule(static) proc_bind(close) num_threads(DSTK_WINDOW)
            do iopen = 1,nopen
                istk = stk_from + iopen - 1
                do iimg = 1,DSTK_NIMGS
                    call dstkios(iopen)%read(stknames(istk), iimg, read_imgs(istk,iimg))
                enddo
            enddo
            !$omp end parallel do
            do iopen = 1,nopen
                call dstkios(iopen)%kill
            enddo
        enddo
        nwrong = 0
        do istk = 1,DSTK_NSTKS
            do iimg = 1,DSTK_NIMGS
                call read_imgs(istk,iimg)%get_rmat_ptr(rmat)
                if( maxval(abs(rmat(1:BOX,1:BOX,1) - real(100*istk + iimg))) > 1.e-6 ) nwrong = nwrong + 1
                call read_imgs(istk,iimg)%kill
            enddo
            call del_file(stknames(istk))
        enddo
        call assert_int(0, nwrong, label//': every image of every stack reads back as written')
    end subroutine read_stacks_concurrently

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

    !> image index as background; first row: exactly representable values (as in fill_pattern), and at
    !! 9 and 10 values float16 cannot hold: 1 + 2^-11 lies halfway and rounds to even (1), 1.0006 rounds to
    !! 1 + 2^-10; quantized=.true. gives what must be read back
    subroutine set_f16_rounding_pattern(img,iimg,quantized)
        type(image), intent(inout) :: img
        integer,     intent(in)    :: iimg
        logical,     intent(in)    :: quantized
        real(kind=c_float), pointer :: pixels(:,:,:) => null()
        call img%get_rmat_ptr(pixels)
        pixels = real(iimg,kind=c_float)
        pixels(1,1,1)  = 2.0_c_float
        pixels(2,1,1)  = -2.0_c_float
        pixels(3,1,1)  = 1.0_c_float
        pixels(4,1,1)  = -1.0_c_float
        pixels(5,1,1)  = 0.5_c_float
        pixels(6,1,1)  = 0.333251953125_c_float
        pixels(7,1,1)  = 6.103515625e-5_c_float
        pixels(8,1,1)  = 65504.0_c_float
        pixels(11,1,1) = 0.0_c_float
        pixels(12,1,1) = -0.0_c_float
        if( quantized )then
            pixels(9,1,1)  = 1.0_c_float
            pixels(10,1,1) = 1.0009765625_c_float
        else
            pixels(9,1,1)  = 1.00048828125_c_float
            pixels(10,1,1) = 1.0006_c_float
        endif
    end subroutine set_f16_rounding_pattern

    !> the first row of the rounding pattern as IEEE half-precision bit patterns on disk
    subroutine assert_f16_payload(fname, label)
        class(string),    intent(in) :: fname
        character(len=*), intent(in) :: label
        integer(int32), parameter :: EXPECTED_BITS(12) = [16384_int32,49152_int32,15360_int32,48128_int32, &
            &14336_int32,13653_int32,1024_int32,31743_int32,15360_int32,15361_int32,0_int32,32768_int32]
        integer(int32) :: bits(12)
        call read_f16_first_row(fname, bits)
        call assert_true(all(bits == EXPECTED_BITS), label//': half-precision bit patterns on disk')
    end subroutine assert_f16_payload

    !> the leading pixels of the first image of a float16 MRC file as unsigned 16-bit patterns
    subroutine read_f16_first_row(fname, bits)
        class(string),  intent(in)  :: fname
        integer(int32), intent(out) :: bits(:)
        integer(int16) :: plane(BOX,BOX)
        integer :: funit, io_stat, ipixel
        bits = -1
        open(newunit=funit,file=fname%to_char(),access='stream',form='unformatted',action='read',status='old',iostat=io_stat)
        call assert_int(0, io_stat, fname%to_char()//': opens for the payload check')
        if( io_stat /= 0 ) return
        read(unit=funit,pos=MRC_HEADER_NBYTES+1,iostat=io_stat) plane
        close(funit)
        call assert_int(0, io_stat, fname%to_char()//': first image payload readable')
        if( io_stat /= 0 ) return
        do ipixel = 1,size(bits)
            bits(ipixel) = iand(int(plane(ipixel,1),int32),65535_int32)
        enddo
    end subroutine read_f16_first_row

    !> a mode-1 (int16) MRC stack of DSTK_NIMGS constant images 100*istk + iimg, written by hand
    subroutine write_int16_stack(stkname, istk)
        class(string), intent(in) :: stkname
        integer,       intent(in) :: istk
        type(MrcImgHead) :: header
        integer(int16) :: plane(BOX,BOX)
        integer :: funit, io_stat, iimg
        integer(kind=8) :: first_byte, image_nbytes
        call header%new([BOX,BOX,DSTK_NIMGS])
        call header%setMode(1)
        call header%setPixSz(SMPD)
        call header%setMinPixVal(real(100*istk + 1))
        call header%setMaxPixVal(real(100*istk + DSTK_NIMGS))
        call header%setMean(real(100*istk) + real(DSTK_NIMGS + 1) / 2.)
        open(newunit=funit, file=stkname%to_char(), access='stream', form='unformatted', &
            &action='readwrite', status='replace', iostat=io_stat)
        call assert_int(0, io_stat, 'int16 test stack opens for writing')
        if( io_stat /= 0 ) return
        call header%write(funit)
        first_byte   = int(header%firstDataByte(),kind=8)
        image_nbytes = int(BOX * BOX * 2,kind=8)
        do iimg = 1,DSTK_NIMGS
            plane = int(100*istk + iimg, int16)
            write(unit=funit, pos=first_byte + int(iimg - 1,kind=8) * image_nbytes, iostat=io_stat) plane
        enddo
        close(funit)
        call header%kill
    end subroutine write_int16_stack

    subroutine cleanup()
        call del_file(string(SOURCE_STACK))
        call del_file(string(IMAGE_COPY))
        call del_file(string(STACK_COPY))
        call del_file(string(FLOAT16_STACK))
        call del_file(string(CHUNK_STACK))
        call del_file(string(ROUND_STACK))
        call del_file(string(F16_IMAGE))
        call del_file(string(F16_BOUNDS))
    end subroutine cleanup

end module simple_stack_io_tester
