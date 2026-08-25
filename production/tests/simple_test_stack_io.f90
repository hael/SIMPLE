program simple_test_stack_io
use simple_core_module_api
use simple_stack_io, only: stack_io
use simple_image,    only: image
use simple_imghead,  only: MrcImgHead, MRC_MODE_FLOAT16, MRC_NVERSION_20141
implicit none
#include "simple_local_flags.inc"

integer, parameter :: BOX = 16, NIMGS = 5, BUFSZ = 2
integer, parameter :: BENCH_BOX = 256, BENCH_NIMGS = 128, BENCH_BUFSZ = 16, BENCH_REPEATS = 20
integer, parameter :: CHUNK_BOX = 1025, CHUNK_NIMGS = 2
real,    parameter :: SMPD = 1.3
character(len=*), parameter :: SOURCE_STACK  = 'simple_test_stack_io_source.mrc'
character(len=*), parameter :: IMAGE_COPY    = 'simple_test_stack_io_image_copy.mrc'
character(len=*), parameter :: STACK_COPY    = 'simple_test_stack_io_stack_copy.mrc'
character(len=*), parameter :: FLOAT16_STACK = 'simple_test_stack_io_float16.mrc'
character(len=*), parameter :: FLOAT16_CHUNK_STACK = 'simple_test_stack_io_float16_chunked.mrc'
integer :: tests_passed

tests_passed = 0

call create_synthetic_stack(string(SOURCE_STACK))
call verify_stack(string(SOURCE_STACK))
call report_pass('synthetic image stack read through stack_io')

call copy_with_image_reader(string(SOURCE_STACK),string(IMAGE_COPY))
call verify_stack(string(IMAGE_COPY))
call report_pass('image reader to stack_io writer round trip')

call copy_with_stack_io(string(SOURCE_STACK),string(STACK_COPY))
call verify_stack(string(STACK_COPY))
call report_pass('stack_io reader/writer copy round trip')

call create_float16_stack(string(FLOAT16_STACK))
call assert_float16_header(string(FLOAT16_STACK))
call report_pass('float16 MRC header mode/version/dimensions/size')
call verify_stack(string(FLOAT16_STACK))
call report_pass('float16 stack_io write/read round trip')
call test_float16_chunked_roundtrip
call benchmark_float_writes

call del_file(string(SOURCE_STACK))
call del_file(string(IMAGE_COPY))
call del_file(string(STACK_COPY))
call del_file(string(FLOAT16_STACK))

write(logfhandle,'(a,i0,a)') '=== ALL ', tests_passed, ' STACK_IO TESTS PASSED ==='
call simple_end('**** SIMPLE_TEST_STACK_IO NORMAL STOP ****')

contains

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

    subroutine copy_with_image_reader(source,dest)
        class(string), intent(in) :: source, dest
        type(stack_io) :: writer
        type(image)    :: img
        integer :: iimg
        call img%new([BOX,BOX,1],SMPD,wthreads=.false.)
        call writer%open(dest,SMPD,'write',box=BOX,bufsz=BUFSZ)
        do iimg = 1,NIMGS
            call img%read(source,iimg)
            call writer%write(iimg,img)
        enddo
        call writer%close
        call img%kill
    end subroutine copy_with_image_reader

    subroutine copy_with_stack_io(source,dest)
        class(string), intent(in) :: source, dest
        type(stack_io) :: reader, writer
        type(image)    :: img
        integer :: iimg
        call reader%open(source,SMPD,'read',bufsz=BUFSZ)
        if( reader%get_nptcls() /= NIMGS ) THROW_HARD('synthetic source image count mismatch')
        call writer%open(dest,SMPD,'write',box=BOX,bufsz=BUFSZ)
        call img%new([BOX,BOX,1],SMPD,wthreads=.false.)
        do iimg = 1,NIMGS
            call reader%read(iimg,img)
            call writer%write(iimg,img)
        enddo
        call reader%close
        call writer%close
        call img%kill
    end subroutine copy_with_stack_io

    subroutine create_float16_stack(fname)
        class(string), intent(in) :: fname
        type(stack_io) :: writer
        type(image)    :: img
        integer :: iimg
        call writer%open(fname,SMPD,'write',box=BOX,bufsz=BUFSZ,wfloat16=.true.)
        call img%new([BOX,BOX,1],SMPD,wthreads=.false.)
        do iimg = 1,NIMGS
            call set_pattern(img,iimg)
            call writer%write(iimg,img)
        enddo
        call writer%close
        call img%kill
    end subroutine create_float16_stack

    subroutine test_float16_chunked_roundtrip
        type(stack_io) :: writer, reader
        type(image)    :: source, actual
        real(kind=c_float), pointer :: source_pixels(:,:,:) => null(), actual_pixels(:,:,:) => null()
        integer :: iimg, ldim(3)
        call source%new([CHUNK_BOX,CHUNK_BOX,1],SMPD,wthreads=.false.)
        call actual%new([CHUNK_BOX,CHUNK_BOX,1],SMPD,wthreads=.false.)
        call source%get_rmat_ptr(source_pixels)
        source_pixels(1:CHUNK_BOX,1:CHUNK_BOX,1) = 1.25_c_float
        call writer%open(string(FLOAT16_CHUNK_STACK),SMPD,'write',box=CHUNK_BOX,bufsz=CHUNK_NIMGS,wfloat16=.true.)
        do iimg = 1,CHUNK_NIMGS
            call writer%write(iimg,source)
        enddo
        call writer%close
        call reader%open(string(FLOAT16_CHUNK_STACK),SMPD,'read',bufsz=1)
        ldim = reader%get_ldim()
        if( reader%get_nptcls() /= CHUNK_NIMGS ) THROW_HARD('chunked float16 image count mismatch')
        if( any(ldim /= [CHUNK_BOX,CHUNK_BOX,1]) ) THROW_HARD('chunked float16 dimensions mismatch')
        do iimg = 1,CHUNK_NIMGS
            call reader%read(iimg,actual)
            call actual%get_rmat_ptr(actual_pixels)
            if( maxval(abs(actual_pixels(1:CHUNK_BOX,1:CHUNK_BOX,1)-1.25_c_float)) > 1.e-7 )then
                THROW_HARD('chunked float16 round trip failed')
            endif
        enddo
        call reader%close
        call source%kill
        call actual%kill
        call del_file(string(FLOAT16_CHUNK_STACK))
        call report_pass('float16 chunked write/read round trip')
    end subroutine test_float16_chunked_roundtrip

    subroutine benchmark_float_writes
        character(len=*), parameter :: FLOAT32_FILE = 'simple_test_stack_io_bench_float32.mrc'
        character(len=*), parameter :: FLOAT16_FILE = 'simple_test_stack_io_bench_float16.mrc'
        type(image) :: bench_img
        real(kind=c_float), pointer :: pixels(:,:,:) => null()
        real(timer_int_kind) :: elapsed_float32, elapsed_float16
        real(dp) :: output_mib_float32, output_mib_float16, total_mvox, time_float32, time_float16
        integer :: irepeat
        call bench_img%new([BENCH_BOX,BENCH_BOX,1],SMPD,wthreads=.false.)
        call bench_img%get_rmat_ptr(pixels)
        pixels = 1.25_c_float
        elapsed_float32 = 0.0_timer_int_kind
        elapsed_float16 = 0.0_timer_int_kind
        do irepeat = 1,BENCH_REPEATS
            if( is_odd(irepeat) )then
                elapsed_float32 = elapsed_float32 + time_stack_write(string(FLOAT32_FILE),bench_img,.false.)
                elapsed_float16 = elapsed_float16 + time_stack_write(string(FLOAT16_FILE),bench_img,.true.)
            else
                elapsed_float16 = elapsed_float16 + time_stack_write(string(FLOAT16_FILE),bench_img,.true.)
                elapsed_float32 = elapsed_float32 + time_stack_write(string(FLOAT32_FILE),bench_img,.false.)
            endif
        enddo
        elapsed_float32 = elapsed_float32 / real(BENCH_REPEATS,kind=timer_int_kind)
        elapsed_float16 = elapsed_float16 / real(BENCH_REPEATS,kind=timer_int_kind)
        time_float32 = max(real(elapsed_float32,kind=dp),DTINY)
        time_float16 = max(real(elapsed_float16,kind=dp),DTINY)
        total_mvox = real(BENCH_BOX*BENCH_BOX,kind=dp) * real(BENCH_NIMGS,kind=dp) / 1.0e6_dp
        output_mib_float32 = total_mvox * 1.0e6_dp * 4.0_dp / real(1024**2,kind=dp)
        output_mib_float16 = total_mvox * 1.0e6_dp * 2.0_dp / real(1024**2,kind=dp)
        write(logfhandle,'(a)') '=== WRITE-ONLY FLOAT32 VS FLOAT16 MRC BENCHMARK ==='
        write(logfhandle,'(a,i0,a,i0,a,i0)') 'box=', BENCH_BOX, ', images=', BENCH_NIMGS, ', repeats=', BENCH_REPEATS
        write(logfhandle,'(a,f10.4,a,f10.2,a,f10.2,a)') 'float32: ', time_float32, ' s, ', &
            &total_mvox/time_float32, ' Mvox/s, ', output_mib_float32/time_float32, ' MiB/s'
        write(logfhandle,'(a,f10.4,a,f10.2,a,f10.2,a)') 'float16: ', time_float16, ' s, ', &
            &total_mvox/time_float16, ' Mvox/s, ', output_mib_float16/time_float16, ' MiB/s'
        write(logfhandle,'(a,f10.3)') 'float16 write speedup (float32 time / float16 time): ', time_float32/time_float16
        call benchmark_float_reads(string(FLOAT32_FILE),string(FLOAT16_FILE))
        call bench_img%kill
        call del_file(string(FLOAT32_FILE))
        call del_file(string(FLOAT16_FILE))
    end subroutine benchmark_float_writes

    subroutine benchmark_float_reads(float32_file,float16_file)
        class(string), intent(in) :: float32_file, float16_file
        real(timer_int_kind) :: elapsed_float32, elapsed_float16
        real(dp) :: input_mib_float32, input_mib_float16, total_mvox, time_float32, time_float16
        integer :: irepeat
        ! Warm the filesystem cache and initialize the float16 converter dispatch before timing.
        elapsed_float32 = time_stack_read(float32_file)
        elapsed_float16 = time_stack_read(float16_file)
        elapsed_float32 = 0.0_timer_int_kind
        elapsed_float16 = 0.0_timer_int_kind
        do irepeat = 1,BENCH_REPEATS
            if( is_odd(irepeat) )then
                elapsed_float32 = elapsed_float32 + time_stack_read(float32_file)
                elapsed_float16 = elapsed_float16 + time_stack_read(float16_file)
            else
                elapsed_float16 = elapsed_float16 + time_stack_read(float16_file)
                elapsed_float32 = elapsed_float32 + time_stack_read(float32_file)
            endif
        enddo
        elapsed_float32 = elapsed_float32 / real(BENCH_REPEATS,kind=timer_int_kind)
        elapsed_float16 = elapsed_float16 / real(BENCH_REPEATS,kind=timer_int_kind)
        time_float32 = max(real(elapsed_float32,kind=dp),DTINY)
        time_float16 = max(real(elapsed_float16,kind=dp),DTINY)
        total_mvox = real(BENCH_BOX*BENCH_BOX,kind=dp) * real(BENCH_NIMGS,kind=dp) / 1.0e6_dp
        input_mib_float32 = total_mvox * 1.0e6_dp * 4.0_dp / real(1024**2,kind=dp)
        input_mib_float16 = total_mvox * 1.0e6_dp * 2.0_dp / real(1024**2,kind=dp)
        write(logfhandle,'(a)') '=== CACHE-WARMED READ-ONLY FLOAT32 VS FLOAT16 MRC BENCHMARK ==='
        write(logfhandle,'(a,i0,a,i0,a,i0)') 'box=', BENCH_BOX, ', images=', BENCH_NIMGS, ', repeats=', BENCH_REPEATS
        write(logfhandle,'(a,f10.4,a,f10.2,a,f10.2,a)') 'float32: ', time_float32, ' s, ', &
            &total_mvox/time_float32, ' Mvox/s, ', input_mib_float32/time_float32, ' MiB/s'
        write(logfhandle,'(a,f10.4,a,f10.2,a,f10.2,a)') 'float16: ', time_float16, ' s, ', &
            &total_mvox/time_float16, ' Mvox/s, ', input_mib_float16/time_float16, ' MiB/s'
        write(logfhandle,'(a,f10.3)') 'float16 read speedup (float32 time / float16 time): ', time_float32/time_float16
    end subroutine benchmark_float_reads

    function time_stack_write(fname,bench_img,wfloat16) result( elapsed )
        class(string), intent(in)    :: fname
        type(image),   intent(inout) :: bench_img
        logical,       intent(in)    :: wfloat16
        real(timer_int_kind) :: elapsed
        integer(timer_int_kind) :: start_time
        type(stack_io) :: writer
        integer :: iimg
        call writer%open(fname,SMPD,'write',box=BENCH_BOX,bufsz=BENCH_BUFSZ,wfloat16=wfloat16)
        start_time = tic()
        do iimg = 1,BENCH_NIMGS
            call writer%write(iimg,bench_img)
        enddo
        call writer%close
        elapsed = toc(start_time)
    end function time_stack_write

    function time_stack_read(fname) result( elapsed )
        class(string), intent(in) :: fname
        real(timer_int_kind) :: elapsed
        integer(timer_int_kind) :: start_time
        type(stack_io) :: reader
        type(image)    :: bench_img
        real(kind=c_float), pointer :: pixels(:,:,:) => null()
        integer :: iimg
        call reader%open(fname,SMPD,'read',bufsz=BENCH_BUFSZ)
        call bench_img%new([BENCH_BOX,BENCH_BOX,1],SMPD,wthreads=.false.)
        start_time = tic()
        do iimg = 1,BENCH_NIMGS
            call reader%read(iimg,bench_img)
        enddo
        elapsed = toc(start_time)
        call reader%close
        call bench_img%get_rmat_ptr(pixels)
        if( any(pixels(1:BENCH_BOX,1:BENCH_BOX,1) /= 1.25_c_float) ) THROW_HARD('float read benchmark data mismatch')
        call bench_img%kill
    end function time_stack_read

    subroutine verify_stack(fname)
        class(string), intent(in) :: fname
        type(stack_io) :: reader
        type(image)    :: img
        integer :: iimg, ldim(3)
        call reader%open(fname,SMPD,'read',bufsz=BUFSZ)
        ldim = reader%get_ldim()
        if( reader%get_nptcls() /= NIMGS ) THROW_HARD('stack image count mismatch')
        if( any(ldim /= [BOX,BOX,1]) ) THROW_HARD('stack dimensions mismatch')
        call img%new(ldim,SMPD,wthreads=.false.)
        do iimg = 1,NIMGS
            call reader%read(iimg,img)
            call assert_pattern(img,iimg)
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

    subroutine assert_pattern(img,iimg)
        type(image), intent(inout) :: img
        integer,     intent(in)    :: iimg
        real(kind=c_float), pointer :: pixels(:,:,:) => null()
        real(kind=c_float) :: expected(BOX,BOX)
        call fill_pattern(expected,iimg)
        call img%get_rmat_ptr(pixels)
        if( maxval(abs(pixels(1:BOX,1:BOX,1)-expected)) > 1.e-7 ) THROW_HARD('stack pixel values mismatch')
    end subroutine assert_pattern

    subroutine fill_pattern(pixels,iimg)
        real(kind=c_float), intent(out) :: pixels(:,:)
        integer,            intent(in)  :: iimg
        pixels = real(iimg,kind=c_float)
        pixels(1,1) = 2.0_c_float
        pixels(2,1) = -2.0_c_float
        pixels(3,1) = 1.0_c_float
        pixels(4,1) = -1.0_c_float
        pixels(5,1) = 0.5_c_float
        pixels(6,1) = 0.333251953125_c_float
        pixels(7,1) = 6.103515625e-5_c_float
        pixels(8,1) = 65504.0_c_float
        pixels(9,1) = 0.0_c_float
        pixels(10,1) = -0.0_c_float
    end subroutine fill_pattern

    subroutine assert_float16_header(fname)
        class(string), intent(in) :: fname
        type(MrcImgHead) :: header
        integer :: funit, io_stat, ldim(3)
        integer(kind=8) :: file_nbytes, expected_nbytes
        call header%new([BOX,BOX,NIMGS])
        open(newunit=funit,file=trim(fname%to_char()),access='stream',form='unformatted', &
            &action='read',status='old',iostat=io_stat)
        if( io_stat /= 0 ) THROW_HARD('failed to open float16 MRC header')
        call header%read(funit)
        close(funit)
        ldim = header%getDims()
        if( header%getMode() /= MRC_MODE_FLOAT16 ) THROW_HARD('float16 MRC mode mismatch')
        if( header%nversion /= MRC_NVERSION_20141 ) THROW_HARD('float16 MRC version mismatch')
        if( any(ldim /= [BOX,BOX,NIMGS]) ) THROW_HARD('float16 MRC dimensions mismatch')
        inquire(file=trim(fname%to_char()),size=file_nbytes,iostat=io_stat)
        if( io_stat /= 0 ) THROW_HARD('failed to inspect float16 MRC size')
        expected_nbytes = 1024_8 + int(2*BOX*BOX*NIMGS,kind=8)
        if( file_nbytes /= expected_nbytes ) THROW_HARD('float16 MRC file size mismatch')
        call header%kill
    end subroutine assert_float16_header

    subroutine report_pass(label)
        character(len=*), intent(in) :: label
        tests_passed = tests_passed + 1
        write(logfhandle,'(a)') '>>> TEST ['//trim(label)//'] PASS'
    end subroutine report_pass

end program simple_test_stack_io
