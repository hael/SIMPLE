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
use simple_syslib,       only: del_file, file_exists, simple_rename
use simple_image,        only: image
use simple_stack_io,     only: stack_io
use simple_discrete_stack_io, only: dstack_io
use simple_stack_read_session, only: stack_read_session, stack_read_counts, plan_spans, SPAN_GAP_BUDGET, &
    &PREFETCH_MAX_STACKS
use simple_syslib,       only: io_set_faults, open_fd_count, set_open_files_limit, ro_open, ro_close, fadvise_willneed
use, intrinsic :: iso_fortran_env, only: int64
use simple_forked_process, only: forked_process, FORK_STATUS_FAILED
use simple_image_slab,   only: image_slab, SLAB_OK, SLAB_OVER_BUDGET, SLAB_NO_MEMORY
!$ use omp_lib,          only: omp_get_max_threads, omp_get_num_procs, omp_set_num_threads
use simple_cmdline,      only: cmdline
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
! read session: stacks of 1, 5 and 200 images; a large-box stack for the gap budget
integer, parameter :: SESS_NIMGS(3) = [1, 5, 200]
integer, parameter :: NMANY_STKS = 20   ! more stacks than the read streams of any device (16 at most)

! a child whose read must stop it: a stack of another box, a missing or truncated stack, or descriptors exhausted
integer, parameter :: FAIL_BOX = 1, FAIL_MISSING = 2, FAIL_TRUNCATED = 3, FAIL_NOFILES = 4
type, extends(forked_process) :: failing_reader
    integer :: what = FAIL_BOX
contains
    procedure :: execute => read_and_fail
end type failing_reader
character(len=*), parameter :: TRUNC_STACK  = 'tmp_session_truncated.mrc'
character(len=*), parameter :: EXTHDR_STACK = 'tmp_session_exthdr.mrc'
character(len=*), parameter :: BIG_STACK    = 'tmp_session_sparse_4GiB.mrc'
character(len=*), parameter :: FT_STACK     = 'tmp_session_mode4.mrc'

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
        call test_session_bitexact(2)
        call test_session_bitexact(MRC_MODE_FLOAT16)
        call test_session_bitexact(0)
        call test_session_bitexact(1)
        call test_session_bitexact(6)
        call test_session_many_stacks()
        call test_session_one_open_per_stack()
        call test_session_spans()
        call test_session_gap_budget()
        call test_session_dimension_mismatch_stops()
        call test_session_prefetch()
        call test_posix_extended_header()
        call test_posix_offsets_beyond_4GiB()
        call test_posix_missing_and_truncated_stop()
        call test_posix_short_reads_eintr_and_advice()
        call test_posix_fallback_mode()
        call test_posix_parallel_opens_no_leak()
        call test_posix_descriptor_exhaustion_stops()
        call test_slab_roundtrip_reuse_and_immutable()
        call test_slab_budget_and_failed_allocation()
        call test_slab_float32_nan_payload()
        call test_slab_get_many_threads()
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
    ! a buffer larger than the stack is clamped to it, forward reads may skip windows, and
    ! backward and out-of-order reads reload the window that holds the image
    subroutine test_buffer_sizes_and_forward_skips()
        integer, parameter :: OUT_OF_ORDER(7) = [4, 1, 5, 2, 3, 3, 1]
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
        ! backward, buffer of two: 5 (partial window [5]), 4 and 3 ([3,4]), 2 and 1 ([1,2])
        call reader%open(string(SOURCE_STACK),SMPD,'read',bufsz=BUFSZ)
        do iimg = NIMGS,1,-1
            call reader%read(iimg,img)
            call assert_pattern(img,iimg,'backward read, buffer of two')
        enddo
        call reader%close
        ! out of order, buffer of two: every window is left and re-entered
        call reader%open(string(SOURCE_STACK),SMPD,'read',bufsz=BUFSZ)
        do iimg = 1,size(OUT_OF_ORDER)
            call reader%read(OUT_OF_ORDER(iimg),img)
            call assert_pattern(img,OUT_OF_ORDER(iimg),'out-of-order read, buffer of two')
        enddo
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
        integer :: istk, iimg, stk_from, stk_to, iopen, nopen, nwrong
        do stk_from = 1,DSTK_NSTKS,DSTK_WINDOW
            stk_to = min(stk_from + DSTK_WINDOW - 1, DSTK_NSTKS)
            nopen  = stk_to - stk_from + 1
            do iopen = 1,nopen
                istk = stk_from + iopen - 1
                call dstkios(iopen)%new(SMPD, BOX)
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

    !---------------- the read session ----------------

    !> Images delivered by the session are byte for byte those of single reads, in the order asked,
    !! for sorted, unsorted and repeated indices over stacks of 1, 5 and 200 images, read in batches
    !! of four so that stacks are split across batch boundaries
    subroutine test_session_bitexact( mode )
        integer, intent(in) :: mode
        type(string) :: stks(3)
        character(len=8) :: tag
        integer, allocatable :: istk(:), inds(:)
        tag = 'mode '//int2str(mode)
        write(*,'(A)') 'test_session_bitexact '//trim(tag)
        call write_session_stacks(stks, mode)
        ! sorted: whole small stacks, then runs and gaps in the large one
        istk = [1, 2,2,2,2,2, 3,3,3,3,3,3,3,3,3,3,3,3,3,3]
        inds = [1, 1,2,3,4,5, 1,2,3,4,7,8,20,21,22,100,150,151,199,200]
        call check_session_against_single_reads(stks, istk, inds, trim(tag)//' sorted')
        ! unsorted, the stacks interleaved
        istk = [3,2,3,1,2,3,3,2,3]
        inds = [150,3,2,1,1,151,1,5,77]
        call check_session_against_single_reads(stks, istk, inds, trim(tag)//' unsorted')
        ! repeated indices
        istk = [3,3,2,3,2,3,1,1]
        inds = [10,10,2,11,2,10,1,1]
        call check_session_against_single_reads(stks, istk, inds, trim(tag)//' duplicates')
        call delete_stacks(stks)
    end subroutine test_session_bitexact

    !> A pass over more stacks than there are read streams opens every stack once and reads it
    !! exactly, the stacks split across batches of four
    subroutine test_session_many_stacks()
        type(string) :: stks(NMANY_STKS)
        type(string) :: one(3)
        integer :: istk(5*NMANY_STKS), inds(5*NMANY_STKS), i, k, nopens
        write(*,'(A)') 'test_session_many_stacks'
        do i = 1, NMANY_STKS
            call write_session_stacks(one, 2)
            stks(i) = 'tmp_session_many_'//int2str(i)//'.mrc'
            call simple_rename(one(2), stks(i))
            call delete_stacks(one)
        end do
        k = 0
        do i = 1, NMANY_STKS
            istk(k+1:k+5) = i
            inds(k+1:k+5) = [1,2,3,4,5]
            k = k + 5
        end do
        call check_session_against_single_reads(stks, istk, inds, 'many stacks', nopens)
        call assert_int(NMANY_STKS, nopens, 'many stacks: every stack opened once')
        do i = 1, NMANY_STKS
            call del_file(stks(i))
        end do
    end subroutine test_session_many_stacks

    !> One open per stack per pass: a stack split across two batches is opened once; a stack
    !! reached again later in a batch (two stack rows sharing one file) is one group, one open
    subroutine test_session_one_open_per_stack()
        type(stack_read_session) :: sess
        type(stack_read_counts)  :: c
        type(string) :: stks(3)
        type(image)  :: imgs(8)
        integer :: i
        write(*,'(A)') 'test_session_one_open_per_stack'
        call write_session_stacks(stks, 2)
        do i = 1, size(imgs)
            call imgs(i)%new([BOX,BOX,1], SMPD, wthreads=.false.)
        end do
        call sess%new(BOX, SMPD, report=.false.)
        ! batch 1 ends inside the 200-image stack, batch 2 continues it
        call sess%read([stks(2),stks(2),stks(2),stks(2),stks(2),stks(3),stks(3)], [1,2,3,4,5,1,2], imgs(1:7))
        call sess%read([stks(3),stks(3),stks(3),stks(3),stks(1)], [3,4,5,6,1], imgs(1:5))
        c = sess%get_counts()
        call assert_int(3, int(c%nopens), 'three stacks, one split across batches: three opens')
        call sess%kill
        ! one file reached by two logical stack rows within a batch
        call sess%new(BOX, SMPD, report=.false.)
        call sess%read([stks(3),stks(2),stks(3)], [1,1,2], imgs(1:3))
        c = sess%get_counts()
        call assert_int(2, int(c%nopens), 'a file reached twice in a batch is opened once')
        call assert_int(2, int(c%nrequests), 'its two adjacent images are one request')
        call sess%kill
        do i = 1, size(imgs)
            call imgs(i)%kill
        end do
        call delete_stacks(stks)
    end subroutine test_session_one_open_per_stack

    !> Spans: adjacent and near images share a request; the over-read of a span never exceeds the
    !! bytes it needs, so distant images stay separate requests
    subroutine test_session_spans()
        type(string) :: stks(3)
        integer(kind=8), parameter :: IMG_BYTES = int(BOX*BOX*4,8)
        write(*,'(A)') 'test_session_spans'
        call write_session_stacks(stks, 2)
        ! one gap of one image inside a run of four: one request over five images
        call check_span_counts(stks(3), [1,2,3,5], 1, 5*IMG_BYTES, 4*IMG_BYTES, 'run with one gap')
        ! every other image: one request, the over-read (4) within the needed bytes (5)
        call check_span_counts(stks(3), [1,3,5,7,9], 1, 9*IMG_BYTES, 5*IMG_BYTES, 'every other image')
        ! distant groups: [1-3], [10], [50-51], [120] are four requests and no over-read
        call check_span_counts(stks(3), [1,2,3,10,50,51,120], 4, 7*IMG_BYTES, 7*IMG_BYTES, 'distant groups')
        ! repeated and unsorted indices: one span over 10..12, each image read once
        call check_span_counts(stks(3), [12,10,10,11,12], 1, 3*IMG_BYTES, 5*IMG_BYTES, 'repeated indices')
        call delete_stacks(stks)
    end subroutine test_session_spans

    !> A span bridges a gap only within SPAN_GAP_BUDGET bytes and within the image cap of a span;
    !! the session's cap depends on its stream count, so the planner is checked with explicit caps
    subroutine test_session_gap_budget()
        integer(kind=8) :: ibytes
        integer :: nsp, first(16), last(16), mfrom(16), mto(16)
        write(*,'(A)') 'test_session_gap_budget'
        ! 576 KiB images: a gap of three (1.7 MiB) is within the budget, a gap of four (2.3 MiB) is not
        ibytes = 384_8 * 384_8 * 4_8
        call assert_true(3_8*ibytes <= SPAN_GAP_BUDGET .and. 4_8*ibytes > SPAN_GAP_BUDGET, 'the images straddle the gap budget')
        call plan_spans([1,2,3,4,5,9], ibytes, 64, nsp, first, last, mfrom, mto)
        call assert_int(1, nsp, 'a gap of three large images is bridged')
        call assert_true(first(1) == 1 .and. last(1) == 9 .and. mfrom(1) == 1 .and. mto(1) == 6, 'one span over 1..9')
        call plan_spans([1,2,3,4,5,10], ibytes, 64, nsp, first, last, mfrom, mto)
        call assert_int(2, nsp, 'a gap of four large images is not')
        call assert_true(last(1) == 5 .and. first(2) == 10 .and. mto(1) == 5 .and. mfrom(2) == 6, 'spans 1..5 and 10')
        ! the image cap: the bridged case split when a span may hold only eight images
        call plan_spans([1,2,3,4,5,9], ibytes, 8, nsp, first, last, mfrom, mto)
        call assert_int(2, nsp, 'a span never exceeds its image cap')
        ! repeats stay in their span and the over-read never exceeds the needed images
        call plan_spans([3,3,4,6,6,9,20], 1024_8, 64, nsp, first, last, mfrom, mto)
        call assert_int(2, nsp, 'repeats and a near gap make one span, a far image another')
        call assert_true(first(1) == 3 .and. last(1) == 9 .and. mto(1) == 6 .and. first(2) == 20, 'spans 3..9 and 20')
    end subroutine test_session_gap_budget

    !> The session checks every stack's header against its box: a mismatch stops the process
    subroutine test_session_dimension_mismatch_stops()
        type(string) :: stks(3)
        write(*,'(A)') 'test_session_dimension_mismatch_stops'
        call write_session_stacks(stks, 2)
        call expect_child_stops(FAIL_BOX, 'reading a stack of another box stops the process')
        call delete_stacks(stks)
    end subroutine test_session_dimension_mismatch_stops

    !> Reading ahead changes no byte: a pass read in batches, each followed by the prefetch of the next,
    !! equals single reads; a correct guess serves every later batch from the window, a wrong one nothing
    subroutine test_session_prefetch()
        type(string) :: stks(3)
        integer :: istk(23), inds(23), k
        integer(kind=8) :: pre, used
        write(*,'(A)') 'test_session_prefetch'
        call write_session_stacks(stks, MRC_MODE_FLOAT16)
        istk = [2,2,2,2,2, (3, k=1,18)]
        inds = [1,2,3,4,5, 1,2,3,4,7,8,20,21,22,40,41,42,100,150,151,152,199,200]
        call check_session_against_single_reads(stks, istk, inds, 'prefetched pass', prefetched=pre, used=used, ahead=1)
        call assert_true(pre > 0_8 .and. used == pre, 'prefetch: the window serves the later batches and nothing else')
        call check_session_against_single_reads(stks, istk, inds, 'wrong guesses', prefetched=pre, used=used, ahead=2)
        call assert_true(pre > 0_8 .and. used == 0_8, 'prefetch: a wrong guess is read ahead and never used')
        call delete_stacks(stks)
    end subroutine test_session_prefetch

    !> An extended header (nsymbt > 0) moves the first image: the descriptor path finds it as the unit path does
    subroutine test_posix_extended_header()
        type(string) :: stks(1)
        integer :: k_all
        write(*,'(A)') 'test_posix_extended_header'
        stks(1) = EXTHDR_STACK
        call write_raw_float_stack(stks(1), 2, 5, 384, [(k_all, k_all=1,5)])
        call check_session_against_single_reads(stks, [1,1,1,1], [5,1,3,4], 'extended header')
        call check_pixel_value(stks(1), 3, 'extended header: image 3 holds its pattern')
        call del_file(stks(1))
    end subroutine test_posix_extended_header

    !> A stack whose last image lies beyond 4 GiB (sparse file): 64-bit offsets on both paths
    subroutine test_posix_offsets_beyond_4GiB()
        type(string) :: stks(1)
        integer, parameter :: NBIG = 4300000   ! 1 KiB images: image NBIG starts beyond 4 GiB
        write(*,'(A)') 'test_posix_offsets_beyond_4GiB'
        stks(1) = BIG_STACK
        call write_raw_float_stack(stks(1), 2, NBIG, 0, [NBIG-1, NBIG])
        call check_session_against_single_reads(stks, [1,1,1], [1,NBIG-1,NBIG], 'offsets beyond 4 GiB')
        call check_pixel_value(stks(1), NBIG, 'beyond 4 GiB: the last image holds its pattern')
        call del_file(stks(1))
    end subroutine test_posix_offsets_beyond_4GiB

    subroutine test_posix_missing_and_truncated_stop()
        write(*,'(A)') 'test_posix_missing_and_truncated_stop'
        call expect_child_stops(FAIL_MISSING,   'a missing stack stops the process')
        call write_raw_float_stack(string(TRUNC_STACK), 2, 5, 0, [1,2,3])
        call expect_child_stops(FAIL_TRUNCATED, 'a truncated stack (end of file before the image) stops the process')
        call del_file(string(TRUNC_STACK))
    end subroutine test_posix_missing_and_truncated_stop

    !> Reads of at most 100 bytes per call with an EINTR every third call, and advice the platform does not
    !! give, change no byte
    subroutine test_posix_short_reads_eintr_and_advice()
        type(string) :: stks(3)
        integer :: fd, rc
        write(*,'(A)') 'test_posix_short_reads_eintr_and_advice'
        call write_session_stacks(stks, 2)
        fd = ro_open(stks(3))
        call assert_true(fd >= 0, 'a stack opens on a plain descriptor')
        rc = fadvise_willneed(fd, 0_int64, 4096_int64)
        call assert_true(rc == 0 .or. rc == 1, 'readahead advice is given, or reported absent where the platform has none')
        call io_set_faults(100_int64, 3, 1)
        call assert_int(1, fadvise_willneed(fd, 0_int64, 4096_int64), 'advice reported unsupported under the hook')
        rc = ro_close(fd)
        call check_session_against_single_reads(stks, [3,3,3,2,2,3,1], [1,2,3,1,5,200,1], 'short reads and EINTR', ahead=1)
        call io_set_faults(0_int64, 0, 0)
        call delete_stacks(stks)
    end subroutine test_posix_short_reads_eintr_and_advice

    !> A mode the descriptor path does not decode (4, complex) goes to the image file reader, unchanged
    subroutine test_posix_fallback_mode()
        type(string) :: stks(1)
        write(*,'(A)') 'test_posix_fallback_mode'
        stks(1) = FT_STACK
        call write_raw_float_stack(stks(1), 4, 3, 0, [1,2,3])
        call check_session_against_single_reads(stks, [1,1,1], [3,1,2], 'fallback mode 4')
        call del_file(stks(1))
    end subroutine test_posix_fallback_mode

    !> Ten passes over more stacks than streams, the stacks opened in parallel each time: no descriptor
    !! outlives its session, and a session holds at most the window's stacks and the boundary stack
    subroutine test_posix_parallel_opens_no_leak()
        type(stack_read_session) :: sess
        type(string) :: stks(NMANY_STKS), one(3)
        type(string) :: names(5*NMANY_STKS)
        type(image)  :: imgs(5*NMANY_STKS)
        integer :: inds(5*NMANY_STKS), i, k, ipass, nfd0, nfd_held
        write(*,'(A)') 'test_posix_parallel_opens_no_leak'
        do i = 1, NMANY_STKS
            call write_session_stacks(one, 2)
            stks(i) = 'tmp_session_many_'//int2str(i)//'.mrc'
            call simple_rename(one(2), stks(i))
            call delete_stacks(one)
        end do
        k = 0
        do i = 1, NMANY_STKS
            names(k+1:k+5) = stks(i)
            inds(k+1:k+5)  = [1,2,3,4,5]
            k = k + 5
        end do
        do i = 1, size(imgs)
            call imgs(i)%new([BOX,BOX,1], SMPD, wthreads=.false.)
        end do
        nfd0 = open_fd_count()
        nfd_held = 0
        do ipass = 1, 10
            call sess%new(BOX, SMPD, report=.false.)
            call sess%read(names(1:50), inds(1:50), imgs(1:50))
            call sess%prefetch(names(51:), inds(51:))
            nfd_held = max(nfd_held, open_fd_count() - nfd0)
            call sess%read(names(51:), inds(51:), imgs(51:))
            call sess%kill
        end do
        call assert_true(nfd_held >= 1 .and. nfd_held <= PREFETCH_MAX_STACKS + 1, &
            &'a session holds the window and the boundary stack open, no more')
        call assert_int(nfd0, open_fd_count(), 'no descriptor outlives its session')
        call assert_true(abs(imgs(100)%get_rmat_at(2,2,1) - 5.) < TOL, 'the last image of the last stack is its fifth')
        do i = 1, size(imgs)
            call imgs(i)%kill
        end do
        do i = 1, NMANY_STKS
            call del_file(stks(i))
        end do
    end subroutine test_posix_parallel_opens_no_leak

    !> Out of descriptors, an open fails and the session stops the process with the reason
    subroutine test_posix_descriptor_exhaustion_stops()
        type(string) :: stks(3)
        write(*,'(A)') 'test_posix_descriptor_exhaustion_stops'
        call write_session_stacks(stks, 2)
        call expect_child_stops(FAIL_NOFILES, 'running out of descriptors stops the process')
        call delete_stacks(stks)
    end subroutine test_posix_descriptor_exhaustion_stops

    !> The particle buffer's slab: images read, put and sealed come back bit for bit (padding zeroed, real space)
    !! into junk-filled images, as often as asked (one copy per state of a weighted reconstruction), and
    !! changing a copy leaves the slab unchanged
    subroutine test_slab_roundtrip_reuse_and_immutable()
        type(stack_read_session) :: sess
        type(image_slab) :: slab
        type(dstack_io)  :: single
        type(string)     :: stks(3)
        type(image)      :: imgs(6), ref, got
        real(kind=c_float), pointer :: a(:,:,:) => null(), b(:,:,:) => null()
        integer :: inds(6), i, status, nwrong, irep
        write(*,'(A)') 'test_slab_roundtrip_reuse_and_immutable'
        call write_session_stacks(stks, MRC_MODE_FLOAT16)
        inds = [200, 1, 77, 2, 150, 3]
        do i = 1, 6
            call imgs(i)%new([BOX,BOX,1], SMPD, wthreads=.false.)
        end do
        call sess%new(BOX, SMPD, report=.false.)
        call sess%read([(stks(3), i=1,6)], inds, imgs)
        call sess%kill
        call slab%new(BOX, 6, 6_int64*BOX*BOX*4, status)
        call assert_int(SLAB_OK, status, 'a slab within its budget is allocated')
        do i = 1, 6
            call slab%put(i, imgs(i))
            call imgs(i)%kill
        end do
        call slab%seal
        call single%new(SMPD, BOX)
        call ref%new([BOX,BOX,1], SMPD, wthreads=.false.)
        call got%new([BOX,BOX,1], SMPD, wthreads=.false.)
        nwrong = 0
        do irep = 1, 3
            do i = 1, 6
                call fill_with_junk(ref)
                call single%read(stks(3), inds(i), ref)
                call fill_with_junk(got)
                call slab%get(i, got)
                call ref%get_rmat_ptr_padded(a)
                call got%get_rmat_ptr_padded(b)
                if( any(transfer(a, 0_int32, size(a)) /= transfer(b, 0_int32, size(b))) ) nwrong = nwrong + 1
                if( got%is_ft() ) nwrong = nwrong + 1
                ! a consumer changes its copy in place; the slab must not see it
                b = 2. * b + 1.
            end do
        end do
        call assert_int(0, nwrong, 'slab copies equal single reads bit for bit, every time, after their copies changed')
        call single%kill
        call ref%kill
        call got%kill
        call slab%kill
        call delete_stacks(stks)
    end subroutine test_slab_roundtrip_reuse_and_immutable

    !> the parallel copy of many slots (the batches served from the particle buffer): at one thread and 20 times at
    !! every processor, every image bit for bit its slot, padding zeroed, real space, whatever the images held
    subroutine test_slab_get_many_threads()
        integer, parameter :: NSLOTS = 97, NGET = 61, SBOX = 40
        type(image_slab) :: slab
        type(image)      :: put_img, ref(NGET), got(NGET)
        real(kind=c_float), pointer :: a(:,:,:) => null(), b(:,:,:) => null()
        integer :: slots(NGET), i, k, irep, status, nwrong, nthr_prev, nthr_full
        write(*,'(A)') 'test_slab_get_many_threads'
        nthr_prev = 1
        nthr_full = 1
        !$ nthr_prev = omp_get_max_threads()
        !$ nthr_full = omp_get_num_procs()
        call slab%new(SBOX, NSLOTS, int(NSLOTS,int64)*SBOX*SBOX*4, status)
        call assert_int(SLAB_OK, status, 'the slab of the parallel copy is allocated')
        call put_img%new([SBOX,SBOX,1], SMPD, wthreads=.false.)
        do i = 1, NSLOTS
            call put_img%get_rmat_ptr(a)
            a(1:SBOX,1:SBOX,1) = reshape([(real(mod(i*7919 + 31*k, 1000)) / 7., k=1,SBOX*SBOX)], [SBOX,SBOX])
            call slab%put(i, put_img)
        end do
        call slab%seal
        ! slots out of order and repeated
        slots = [(mod(37*i, NSLOTS) + 1, i=1,NGET)]
        slots(NGET) = slots(1)
        do i = 1, NGET
            call ref(i)%new([SBOX,SBOX,1], SMPD, wthreads=.false.)
            call got(i)%new([SBOX,SBOX,1], SMPD, wthreads=.false.)
            call fill_with_junk(ref(i))
        end do
        !$ call omp_set_num_threads(1)
        call slab%get_many(slots, ref)
        !$ call omp_set_num_threads(nthr_full)
        nwrong = 0
        do irep = 1, 20
            do i = 1, NGET
                call fill_with_junk(got(i))
            end do
            call slab%get_many(slots, got)
            do i = 1, NGET
                call ref(i)%get_rmat_ptr_padded(a)
                call got(i)%get_rmat_ptr_padded(b)
                if( any(transfer(a, 0_int32, size(a)) /= transfer(b, 0_int32, size(b))) ) nwrong = nwrong + 1
                if( got(i)%is_ft() ) nwrong = nwrong + 1
            end do
        end do
        !$ call omp_set_num_threads(nthr_prev)
        call assert_int(0, nwrong, 'the parallel slab copies equal the copies at one thread, bit for bit')
        nwrong = 0
        do i = 1, NGET
            call slab%get(slots(i), put_img)
            call put_img%get_rmat_ptr_padded(a)
            call ref(i)%get_rmat_ptr_padded(b)
            if( any(transfer(a, 0_int32, size(a)) /= transfer(b, 0_int32, size(b))) ) nwrong = nwrong + 1
        end do
        call assert_int(0, nwrong, 'a copy of many slots equals the copy of each slot')
        do i = 1, NGET
            call ref(i)%kill
            call got(i)%kill
        end do
        call put_img%kill
        call slab%kill
    end subroutine test_slab_get_many_threads

    !> float32 pixels, a signalling NaN with a payload among them, pass the read and the slab bit for bit
    subroutine test_slab_float32_nan_payload()
        integer(int32), parameter :: NAN_BITS = int(z'7FA00001', int32)
        type(stack_read_session) :: sess
        type(image_slab) :: slab
        type(dstack_io)  :: single
        type(string)     :: stks(1)
        type(image)      :: imgs(3), ref, got
        real(kind=c_float), pointer :: a(:,:,:) => null(), b(:,:,:) => null()
        integer :: i, status, funit, nwrong
        write(*,'(A)') 'test_slab_float32_nan_payload'
        stks(1) = 'tmp_slab_nan.mrc'
        call write_raw_float_stack(stks(1), 2, 3, 0, [1,2,3])
        open(newunit=funit, file=stks(1)%to_char(), access='stream', form='unformatted', action='readwrite', status='old')
        write(funit, pos=MRC_HEADER_NBYTES + BOX*BOX*4 + 2*BOX*4 + 2*4 + 1) NAN_BITS   ! image 2, pixel (3,3)
        close(funit)
        do i = 1, 3
            call imgs(i)%new([BOX,BOX,1], SMPD, wthreads=.false.)
        end do
        call sess%new(BOX, SMPD, report=.false.)
        call sess%read([(stks(1), i=1,3)], [1,2,3], imgs)
        call sess%kill
        call slab%new(BOX, 3, 3_int64*BOX*BOX*4, status)
        do i = 1, 3
            call slab%put(i, imgs(i))
            call imgs(i)%kill
        end do
        call slab%seal
        call single%new(SMPD, BOX)
        call ref%new([BOX,BOX,1], SMPD, wthreads=.false.)
        call got%new([BOX,BOX,1], SMPD, wthreads=.false.)
        nwrong = 0
        do i = 1, 3
            call fill_with_junk(ref)
            call single%read(stks(1), i, ref)
            call fill_with_junk(got)
            call slab%get(i, got)
            call ref%get_rmat_ptr_padded(a)
            call got%get_rmat_ptr_padded(b)
            if( any(transfer(a, 0_int32, size(a)) /= transfer(b, 0_int32, size(b))) ) nwrong = nwrong + 1
            if( i == 2 .and. transfer(b(3,3,1), 0_int32) /= NAN_BITS ) nwrong = nwrong + 1
        end do
        call assert_int(0, nwrong, 'float32 and a NaN payload pass the read and the slab bit for bit')
        call single%kill
        call ref%kill
        call got%kill
        call slab%kill
        call del_file(stks(1))
    end subroutine test_slab_float32_nan_payload

    !> A slab above its budget, or one whose allocation fails, stays empty and reports why
    subroutine test_slab_budget_and_failed_allocation()
        type(image_slab) :: slab
        integer :: status
        write(*,'(A)') 'test_slab_budget_and_failed_allocation'
        call slab%new(BOX, 10, 10_int64*BOX*BOX*4 - 1_int64, status)
        call assert_int(SLAB_OVER_BUDGET, status, 'a slab one byte above its budget is refused')
        call assert_false(slab%exists(), 'a refused slab stays empty')
        call slab%new(4096, 1000000000, huge(1_int64), status)
        call assert_int(SLAB_NO_MEMORY, status, 'a slab that cannot be allocated reports it')
        call assert_false(slab%exists(), 'a failed slab stays empty')
        call slab%kill
    end subroutine test_slab_budget_and_failed_allocation

    subroutine expect_child_stops( what, label )
        integer,          intent(in) :: what
        character(len=*), intent(in) :: label
        type(failing_reader) :: child
        child%what = what
        call child%start(name=string('SESSION_FAILURE'))
        call child%await_final_status()
        call assert_int(FORK_STATUS_FAILED, child%status(), label)
        call child%destroy
    end subroutine expect_child_stops

    subroutine read_and_fail( self, cline )
        class(failing_reader), intent(inout) :: self
        class(cmdline),        intent(inout) :: cline
        type(stack_read_session) :: sess
        type(image) :: imgs(1)
        integer :: rc, fd
        select case(self%what)
            case(FAIL_BOX)
                call imgs(1)%new([2*BOX,2*BOX,1], SMPD, wthreads=.false.)
                call sess%new(2*BOX, SMPD, report=.false.)
                call sess%read([string(session_stack_name(2, 2))], [1], imgs)
            case(FAIL_MISSING)
                call imgs(1)%new([BOX,BOX,1], SMPD, wthreads=.false.)
                call sess%new(BOX, SMPD, report=.false.)
                call sess%read([string('tmp_session_no_such_stack.mrc')], [1], imgs)
            case(FAIL_TRUNCATED)
                call imgs(1)%new([BOX,BOX,1], SMPD, wthreads=.false.)
                call sess%new(BOX, SMPD, report=.false.)
                call sess%read([string(TRUNC_STACK)], [5], imgs)
            case(FAIL_NOFILES)
                call imgs(1)%new([BOX,BOX,1], SMPD, wthreads=.false.)
                call sess%new(BOX, SMPD, report=.false.)
                ! an open takes the lowest free descriptor, so a soft limit at that number refuses the next open
                fd = ro_open(string(session_stack_name(3, 2)))
                rc = ro_close(fd)
                rc = set_open_files_limit(fd)
                call sess%read([string(session_stack_name(3, 2))], [1], imgs)
        end select
        ! reached only when the read did not stop the process: a clean exit fails the test
        call sess%kill
        call imgs(1)%kill
    end subroutine read_and_fail

    !> an MRC stack of nimgs float32 (mode 2) or complex-flagged (mode 4) images of box BOX, nsymbt bytes of
    !! extended header, with only the images listed written (the rest left a hole: a sparse file)
    subroutine write_raw_float_stack( stkname, mode, nimgs, nsymbt, written )
        class(string), intent(in) :: stkname
        integer,       intent(in) :: mode, nimgs, nsymbt, written(:)
        type(MrcImgHead)   :: header
        real(kind=c_float) :: plane(BOX,BOX)
        integer(kind=1)    :: ext(max(1,nsymbt))
        integer :: funit, io_stat, k, i, j
        integer(kind=8) :: first_byte, image_nbytes
        call header%new([BOX,BOX,nimgs])
        call header%setMode(mode)
        call header%setPixSz(SMPD)
        open(newunit=funit, file=stkname%to_char(), access='stream', form='unformatted', &
            &action='readwrite', status='replace', iostat=io_stat)
        call assert_int(0, io_stat, 'raw test stack opens for writing')
        if( io_stat /= 0 ) return
        call header%write(funit)
        if( nsymbt > 0 )then
            write(unit=funit, pos=93) int(nsymbt, int32)
            ext = 7_1
            write(unit=funit, pos=MRC_HEADER_NBYTES+1) ext(1:nsymbt)
        endif
        first_byte   = int(MRC_HEADER_NBYTES + nsymbt + 1, kind=8)
        image_nbytes = int(BOX * BOX * 4, kind=8)
        do k = 1, size(written)
            do j = 1, BOX
                do i = 1, BOX
                    plane(i,j) = session_pixel(i, j, written(k), 1)
                end do
            end do
            write(unit=funit, pos=first_byte + int(written(k) - 1,kind=8) * image_nbytes) plane
        end do
        close(funit)
        call header%kill
    end subroutine write_raw_float_stack

    !> pixel (2,2) of image iimg read through a session holds the image index (session_pixel)
    subroutine check_pixel_value( stkname, iimg, label )
        type(string),     intent(in) :: stkname
        integer,          intent(in) :: iimg
        character(len=*), intent(in) :: label
        type(stack_read_session) :: sess
        type(image) :: imgs(1)
        call imgs(1)%new([BOX,BOX,1], SMPD, wthreads=.false.)
        call sess%new(BOX, SMPD, report=.false.)
        call sess%read([stkname], [iimg], imgs)
        call sess%kill
        call assert_true(abs(imgs(1)%get_rmat_at(2,2,1) - real(iimg)) < 1.e-3, label)
        call imgs(1)%kill
    end subroutine check_pixel_value

    !> reads the list in batches of four through one session and compares every image, all of its
    !! storage padding included, bit for bit with a single read of the same image
    !! ahead=1: each batch is followed by the prefetch of the next; ahead=2: by the prefetch of a wrong guess
    subroutine check_session_against_single_reads( stks, istk, inds, label, nopens, prefetched, used, ahead )
        type(string),              intent(in)  :: stks(:)
        integer,                   intent(in)  :: istk(:), inds(:)
        character(len=*),          intent(in)  :: label
        integer,         optional, intent(out) :: nopens
        integer(kind=8), optional, intent(out) :: prefetched, used
        integer,         optional, intent(in)  :: ahead
        type(stack_read_session) :: sess
        type(stack_read_counts)  :: c
        type(dstack_io) :: single
        type(image), allocatable :: imgs(:)
        type(image)  :: ref
        type(string), allocatable :: names(:)
        real(kind=c_float), pointer :: a(:,:,:) => null(), b(:,:,:) => null()
        integer :: n, i, ib, ie, nwrong
        n = size(inds)
        allocate(imgs(n), names(n))
        do i = 1, n
            names(i) = stks(istk(i))
            call imgs(i)%new([BOX,BOX,1], SMPD, wthreads=.false.)
            call fill_with_junk(imgs(i))
        end do
        call sess%new(BOX, SMPD, report=.false.)
        do ib = 1, n, 4
            ie = min(n, ib + 3)
            call sess%read(names(ib:ie), inds(ib:ie), imgs(ib:ie))
            if( present(ahead) .and. ie < n )then
                if( ahead == 1 ) call sess%prefetch(names(ie+1:min(n,ie+4)), inds(ie+1:min(n,ie+4)))
                if( ahead == 2 ) call sess%prefetch(names(ie+1:min(n,ie+4)), wrong_guess(inds(ie+1:min(n,ie+4))))
            endif
        end do
        c = sess%get_counts()
        if( present(nopens) ) nopens = int(c%nopens)
        if( present(prefetched) ) prefetched = c%bytes_prefetched
        if( present(used) )       used       = c%bytes_prefetch_used
        call sess%kill
        call single%new(SMPD, BOX)
        call ref%new([BOX,BOX,1], SMPD, wthreads=.false.)
        nwrong = 0
        do i = 1, n
            call fill_with_junk(ref)
            call single%read(names(i), inds(i), ref)
            call ref%get_rmat_ptr_padded(a)
            call imgs(i)%get_rmat_ptr_padded(b)
            if( any(transfer(a, 0_int32, size(a)) /= transfer(b, 0_int32, size(b))) ) nwrong = nwrong + 1
            if( imgs(i)%is_ft() .neqv. ref%is_ft() ) nwrong = nwrong + 1
            call imgs(i)%kill
        end do
        call single%kill
        call ref%kill
        call assert_int(0, nwrong, label//': every image equals its single read bit for bit')
    end subroutine check_session_against_single_reads

    !> padding and pixels set to a value no read gives, the image flagged as a Fourier transform: a read
    !! must leave neither behind
    subroutine fill_with_junk( img )
        type(image), intent(inout) :: img
        real(kind=c_float), pointer :: a(:,:,:) => null()
        call img%get_rmat_ptr_padded(a)
        a = -12345.
        call img%set_ft(.true.)
    end subroutine fill_with_junk

    !> indices no batch of the prefetch tests asks for: the last image of the stack, then counting down
    function wrong_guess( inds ) result( guess )
        integer, intent(in)  :: inds(:)
        integer :: guess(size(inds)), k
        guess = [(SESS_NIMGS(3) - 60 - k, k=1,size(inds))]
    end function wrong_guess

    subroutine check_span_counts( stkname, inds, nreq, bytes_read, bytes_req, label )
        type(string),     intent(in) :: stkname
        integer,          intent(in) :: inds(:), nreq
        integer(kind=8),  intent(in) :: bytes_read, bytes_req
        character(len=*), intent(in) :: label
        type(stack_read_session) :: sess
        type(stack_read_counts)  :: c
        type(image), allocatable :: imgs(:)
        integer :: i
        allocate(imgs(size(inds)))
        do i = 1, size(inds)
            call imgs(i)%new([BOX,BOX,1], SMPD, wthreads=.false.)
        end do
        call sess%new(BOX, SMPD, report=.false.)
        call sess%read([(stkname, i=1,size(inds))], inds, imgs)
        c = sess%get_counts()
        call assert_int(nreq, int(c%nrequests), label//': requests')
        call assert_true(c%bytes_read == bytes_read, label//': bytes read')
        call assert_true(c%bytes_req  == bytes_req,  label//': bytes delivered')
        do i = 1, size(inds)
            call assert_true(abs(imgs(i)%get_rmat_at(2,2,1) - session_pixel(2, 2, inds(i), 3)) < 1.e-3, &
                &label//': image '//int2str(i)//' holds its index')
            call imgs(i)%kill
        end do
        call sess%kill
    end subroutine check_span_counts

    !> stacks of 1, 5 and 200 images of distinct pixels in MRC mode 2 or 12 (written by stack_io) or
    !! 0, 1 or 6 (written raw: unsigned bytes, signed and unsigned 16-bit integers over their whole range)
    subroutine write_session_stacks( stks, mode )
        type(string), intent(out) :: stks(3)
        integer,      intent(in)  :: mode
        type(stack_io) :: writer
        type(image)    :: img
        real(kind=c_float), pointer :: pixels(:,:,:) => null()
        integer :: istk, iimg, i, j
        do istk = 1, 3
            stks(istk) = session_stack_name(istk, mode)
            select case(mode)
                case(2, MRC_MODE_FLOAT16)
                    call img%new([BOX,BOX,1], SMPD, wthreads=.false.)
                    call writer%open(stks(istk), SMPD, 'write', box=BOX, bufsz=BUFSZ, wfloat16=(mode == MRC_MODE_FLOAT16))
                    do iimg = 1, SESS_NIMGS(istk)
                        call img%get_rmat_ptr(pixels)
                        do j = 1, BOX
                            do i = 1, BOX
                                pixels(i,j,1) = session_pixel(i, j, iimg, istk)
                            end do
                        end do
                        call writer%write(iimg, img)
                    end do
                    call writer%close
                    call img%kill
                case default
                    call write_raw_integer_stack(stks(istk), mode, istk)
            end select
        end do
    end subroutine write_session_stacks

    !> an MRC stack of SESS_NIMGS(istk) integer images in mode 0, 1 or 6 with the image index at pixel (2,2)
    subroutine write_raw_integer_stack( stkname, mode, istk )
        class(string), intent(in) :: stkname
        integer,       intent(in) :: mode, istk
        type(MrcImgHead) :: header
        integer(kind=1)  :: bytes(BOX,BOX)
        integer(int16)   :: words(BOX,BOX)
        integer :: funit, io_stat, iimg, i, j, v, b
        integer(kind=8) :: first_byte, image_nbytes
        call header%new([BOX,BOX,SESS_NIMGS(istk)])
        call header%setMode(mode)
        call header%setPixSz(SMPD)
        open(newunit=funit, file=stkname%to_char(), access='stream', form='unformatted', &
            &action='readwrite', status='replace', iostat=io_stat)
        call assert_int(0, io_stat, 'integer test stack opens for writing')
        if( io_stat /= 0 ) return
        call header%write(funit)
        first_byte   = int(header%firstDataByte(),kind=8)
        image_nbytes = int(BOX * BOX * merge(1, 2, mode == 0),kind=8)
        do iimg = 1, SESS_NIMGS(istk)
            do j = 1, BOX
                do i = 1, BOX
                    v = mod(i*37 + j*101 + iimg*7919 + istk*131, 65536)
                    if( i == 2 .and. j == 2 ) v = iimg
                    b = iand(v, 255)
                    bytes(i,j) = int(merge(b - 256, b, b > 127), kind=1)
                    words(i,j) = int(merge(v - 65536, v, v > 32767), int16)
                end do
            end do
            if( mode == 0 )then
                write(unit=funit, pos=first_byte + int(iimg - 1,kind=8) * image_nbytes, iostat=io_stat) bytes
            else
                write(unit=funit, pos=first_byte + int(iimg - 1,kind=8) * image_nbytes, iostat=io_stat) words
            endif
        end do
        close(funit)
        call header%kill
    end subroutine write_raw_integer_stack

    function session_stack_name( istk, mode ) result( fname )
        integer, intent(in) :: istk, mode
        character(len=:), allocatable :: fname
        fname = 'tmp_session_m'//int2str(mode)//'_'//int2str(SESS_NIMGS(istk))//'.mrc'
    end function session_stack_name

    !> image index at pixel (2,2), elsewhere a smooth signed pattern that float16 does not hold exactly
    real function session_pixel( i, j, iimg, istk )
        integer, intent(in) :: i, j, iimg, istk
        if( i == 2 .and. j == 2 )then
            session_pixel = real(iimg)
        else
            session_pixel = 37.3 * sin(0.37*real(i) + 0.71*real(j) + 1.3*real(iimg) + 2.9*real(istk)) + 0.01*real(iimg)
        endif
    end function session_pixel

    subroutine delete_stacks( stks )
        type(string), intent(inout) :: stks(:)
        integer :: istk
        do istk = 1, size(stks)
            call del_file(stks(istk))
            call stks(istk)%kill
        end do
    end subroutine delete_stacks

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
