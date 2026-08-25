program simple_test_discrete_stack_io
use, intrinsic :: iso_fortran_env, only: int16, int32
use simple_core_module_api
use simple_discrete_stack_io, only: dstack_io
use simple_image,             only: image
use simple_imghead,           only: find_ldim_nptcls, MrcImgHead, MRC_MODE_FLOAT16, MRC_NVERSION_20141
use simple_stack_io,          only: stack_io
implicit none
#include "simple_local_flags.inc"

integer, parameter :: NSTKS = 12, NIMGS = 4, BOX = 32, OPEN_WINDOW = 3
real,    parameter :: SMPD = 1.0
type(dstack_io) :: dstkios(OPEN_WINDOW)
type(image)     :: img, read_imgs(NSTKS,NIMGS)
type(string)    :: stknames(NSTKS)
real(kind=c_float), pointer :: rmat(:,:,:) => null()
integer :: istk, iimg, ldim(3), nptcls, stk_from, stk_to, iopen, nopen, tests_passed
real    :: expected

tests_passed = 0
call img%new([BOX,BOX,1], SMPD)
do istk = 1,NSTKS
    stknames(istk) = 'simple_test_discrete_stack_io_'//int2str(istk)//'.mrc'
    do iimg = 1,NIMGS
        call img%get_rmat_ptr(rmat)
        rmat = real(100*istk + iimg, kind=c_float)
        call img%write(stknames(istk), iimg, del_if_exists=(iimg == 1))
    enddo
enddo
call img%kill

do stk_from = 1,NSTKS,OPEN_WINDOW
    stk_to = min(stk_from + OPEN_WINDOW - 1, NSTKS)
    nopen  = stk_to - stk_from + 1
    do iopen = 1,nopen
        istk = stk_from + iopen - 1
        call find_ldim_nptcls(stknames(istk), ldim, nptcls)
        call dstkios(iopen)%new(SMPD, BOX)
        call dstkios(iopen)%cache_stack_info(stknames(istk), ldim, nptcls)
        call dstkios(iopen)%open(stknames(istk))
        do iimg = 1,NIMGS
            call read_imgs(istk,iimg)%new([BOX,BOX,1], SMPD, wthreads=.false.)
        enddo
    enddo

    !$omp parallel do default(shared) private(iopen,istk,iimg) schedule(static) proc_bind(close) num_threads(OPEN_WINDOW)
    do iopen = 1,nopen
        istk = stk_from + iopen - 1
        do iimg = 1,NIMGS
            call dstkios(iopen)%read(stknames(istk), iimg, read_imgs(istk,iimg))
        enddo
    enddo
    !$omp end parallel do

    do iopen = 1,nopen
        call dstkios(iopen)%kill
    enddo
enddo

do istk = 1,NSTKS
    do iimg = 1,NIMGS
        expected = real(100*istk + iimg)
        call read_imgs(istk,iimg)%get_rmat_ptr(rmat)
        if( maxval(abs(rmat(1:BOX,1:BOX,1) - expected)) > 1.e-6 )then
            write(logfhandle,*) 'istk/iimg/expected/min/max: ', istk, iimg, expected, &
                minval(rmat(1:BOX,1:BOX,1)), maxval(rmat(1:BOX,1:BOX,1))
            THROW_HARD('discrete stack io parallel read returned wrong image')
        endif
        call read_imgs(istk,iimg)%kill
    enddo
    call del_file(stknames(istk))
    call stknames(istk)%kill
enddo
call report_pass('float32 discrete-stack parallel read')

do istk = 1,NSTKS
    stknames(istk) = 'simple_test_discrete_stack_io_i16_'//int2str(istk)//'.mrc'
    call write_int16_stack(stknames(istk), istk)
enddo

do stk_from = 1,NSTKS,OPEN_WINDOW
    stk_to = min(stk_from + OPEN_WINDOW - 1, NSTKS)
    nopen  = stk_to - stk_from + 1
    do iopen = 1,nopen
        istk = stk_from + iopen - 1
        call find_ldim_nptcls(stknames(istk), ldim, nptcls)
        call dstkios(iopen)%new(SMPD, BOX)
        call dstkios(iopen)%cache_stack_info(stknames(istk), ldim, nptcls)
        call dstkios(iopen)%open(stknames(istk))
        do iimg = 1,NIMGS
            call read_imgs(istk,iimg)%new([BOX,BOX,1], SMPD, wthreads=.false.)
        enddo
    enddo

    !$omp parallel do default(shared) private(iopen,istk,iimg) schedule(static) proc_bind(close) num_threads(OPEN_WINDOW)
    do iopen = 1,nopen
        istk = stk_from + iopen - 1
        do iimg = 1,NIMGS
            call dstkios(iopen)%read(stknames(istk), iimg, read_imgs(istk,iimg))
        enddo
    enddo
    !$omp end parallel do

    do iopen = 1,nopen
        call dstkios(iopen)%kill
    enddo
enddo

do istk = 1,NSTKS
    do iimg = 1,NIMGS
        expected = real(100*istk + iimg)
        call read_imgs(istk,iimg)%get_rmat_ptr(rmat)
        if( maxval(abs(rmat(1:BOX,1:BOX,1) - expected)) > 1.e-6 )then
            write(logfhandle,*) 'i16 istk/iimg/expected/min/max: ', istk, iimg, expected, &
                minval(rmat(1:BOX,1:BOX,1)), maxval(rmat(1:BOX,1:BOX,1))
            THROW_HARD('discrete stack io parallel 16-bit read returned wrong image')
        endif
        call read_imgs(istk,iimg)%kill
    enddo
    call del_file(stknames(istk))
    call stknames(istk)%kill
enddo
call report_pass('int16 discrete-stack parallel read')

call test_float16_stack_roundtrip
call test_float16_image_roundtrip
call test_float16_encoder_boundaries

write(logfhandle,'(a,i0,a)') '=== ALL ', tests_passed, ' DISCRETE STACK I/O TESTS PASSED ==='
call simple_end('**** SIMPLE_TEST_DISCRETE_STACK_IO NORMAL STOP ****')

contains

    subroutine report_pass(label)
        character(len=*), intent(in) :: label
        tests_passed = tests_passed + 1
        write(logfhandle,'(a)') '>>> TEST ['//trim(label)//'] PASS'
    end subroutine report_pass

    subroutine test_float16_stack_roundtrip
        integer, parameter :: NTEST_IMGS = 5
        type(stack_io) :: stkio_w, stkio_r
        type(image)    :: source, actual, expected_img
        type(string)   :: stkname
        real(kind=c_float), pointer :: actual_rmat(:,:,:) => null(), expected_rmat(:,:,:) => null()
        integer :: iimg
        stkname = 'simple_test_discrete_stack_io_float16.mrc'
        call source%new([BOX,BOX,1],SMPD,wthreads=.false.)
        call actual%new([BOX,BOX,1],SMPD,wthreads=.false.)
        call expected_img%new([BOX,BOX,1],SMPD,wthreads=.false.)
        call stkio_w%open(stkname,SMPD,'write',box=BOX,bufsz=2,wfloat16=.true.)
        do iimg = 1,NTEST_IMGS
            call set_float16_pattern(source,iimg,quantized=.false.)
            call stkio_w%write(iimg,source)
        enddo
        call stkio_w%close
        call assert_float16_header(stkname,NTEST_IMGS)
        call report_pass('float16 stack header mode/version/dimensions/size')
        call assert_float16_payload(stkname)
        call report_pass('float16 stack payload encoding')
        call stkio_r%open(stkname,SMPD,'read',bufsz=3)
        if( stkio_r%get_nptcls() /= NTEST_IMGS ) THROW_HARD('float16 stack image count mismatch')
        do iimg = 1,NTEST_IMGS
            call stkio_r%read(iimg,actual)
            call set_float16_pattern(expected_img,iimg,quantized=.true.)
            call actual%get_rmat_ptr(actual_rmat)
            call expected_img%get_rmat_ptr(expected_rmat)
            if( maxval(abs(actual_rmat(1:BOX,1:BOX,1)-expected_rmat(1:BOX,1:BOX,1))) > 1.e-7 )then
                THROW_HARD('float16 stack round trip failed')
            endif
        enddo
        call report_pass('float16 stack_io buffered round trip')
        call stkio_r%close
        call source%kill
        call actual%kill
        call expected_img%kill
        call del_file(stkname)
        call stkname%kill
    end subroutine test_float16_stack_roundtrip

    subroutine test_float16_image_roundtrip
        type(image)  :: source, actual, expected_img
        type(string) :: fname
        real(kind=c_float), pointer :: actual_rmat(:,:,:) => null(), expected_rmat(:,:,:) => null()
        real :: stats(4)
        integer :: iimg
        fname = 'simple_test_image_io_float16.mrc'
        call source%new([BOX,BOX,1],SMPD,wthreads=.false.)
        call actual%new([BOX,BOX,1],SMPD,wthreads=.false.)
        call expected_img%new([BOX,BOX,1],SMPD,wthreads=.false.)
        do iimg = 1,2
            call set_float16_pattern(source,iimg,quantized=.false.)
            if( iimg == 1 )then
                call source%write(fname,iimg,del_if_exists=.true.,wfloat16=.true.)
            else
                call source%write(fname,iimg)
            endif
        enddo
        do iimg = 1,2
            call actual%read(fname,iimg)
            call set_float16_pattern(expected_img,iimg,quantized=.true.)
            call actual%get_rmat_ptr(actual_rmat)
            call expected_img%get_rmat_ptr(expected_rmat)
            if( maxval(abs(actual_rmat(1:BOX,1:BOX,1)-expected_rmat(1:BOX,1:BOX,1))) > 1.e-7 )then
                THROW_HARD('float16 image round trip failed')
            endif
        enddo
        call report_pass('float16 image write/read round trip and mode inheritance')
        stats = [-1.0,65504.0,0.0,1.0]
        call source%update_header_stats(fname,stats)
        call assert_float16_header(fname,2)
        call report_pass('float16 image header statistics preserve mode')
        call assert_float16_payload(fname)
        call report_pass('float16 image payload encoding')
        call source%kill
        call actual%kill
        call expected_img%kill
        call del_file(fname)
        call fname%kill
    end subroutine test_float16_image_roundtrip

    subroutine test_float16_encoder_boundaries
        type(image)  :: source, actual
        type(string) :: fname
        real(kind=c_float), pointer :: pixels(:,:,:) => null(), actual_pixels(:,:,:) => null()
        fname = 'simple_test_image_io_float16_boundaries.mrc'
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
        call source%write(fname,1,del_if_exists=.true.,wfloat16=.true.)
        call assert_float16_boundary_payload(fname)
        call actual%read(fname,1)
        call actual%get_rmat_ptr(actual_pixels)
        pixels(5,1,1) = 0.0_c_float
        pixels(6,1,1) = -0.0_c_float
        pixels(7,1,1) = 2.0_c_float**(-23)
        pixels(8,1,1) = -2.0_c_float**(-23)
        if( any(actual_pixels(1:BOX,1:BOX,1) /= pixels(1:BOX,1:BOX,1)) )then
            THROW_HARD('float16 zero/subnormal read back failed')
        endif
        if( transfer(actual_pixels(6,1,1),0_int32) /= not(huge(0_int32)) )then
            THROW_HARD('float16 negative zero sign was not preserved')
        endif
        call report_pass('float16 zero/subnormal exact encoding and decoding')
        call source%kill
        call actual%kill
        call del_file(fname)
        call fname%kill
    end subroutine test_float16_encoder_boundaries

    subroutine set_float16_pattern(img,iimg,quantized)
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
    end subroutine set_float16_pattern

    subroutine assert_float16_header(fname,nimgs)
        type(string), intent(in) :: fname
        integer,      intent(in) :: nimgs
        type(MrcImgHead) :: header
        integer :: funit, io_stat, dims(3)
        integer(kind=8) :: file_nbytes, expected_nbytes
        call header%new([BOX,BOX,nimgs])
        open(newunit=funit,file=trim(fname%to_char()),access='stream',form='unformatted', &
            &action='read',status='old',iostat=io_stat)
        if( io_stat /= 0 ) THROW_HARD('failed to open float16 MRC header')
        call header%read(funit)
        close(funit)
        dims = header%getDims()
        if( header%getMode() /= MRC_MODE_FLOAT16 ) THROW_HARD('float16 MRC mode mismatch')
        if( header%nversion /= MRC_NVERSION_20141 ) THROW_HARD('float16 MRC version mismatch')
        if( any(dims /= [BOX,BOX,nimgs]) ) THROW_HARD('float16 MRC dimensions mismatch')
        inquire(file=trim(fname%to_char()),size=file_nbytes,iostat=io_stat)
        if( io_stat /= 0 ) THROW_HARD('failed to inspect float16 MRC size')
        expected_nbytes = 1024_8 + int(2*BOX*BOX*nimgs,kind=8)
        if( file_nbytes /= expected_nbytes ) THROW_HARD('float16 MRC file size mismatch')
        call header%kill
    end subroutine assert_float16_header

    subroutine assert_float16_payload(fname)
        type(string), intent(in) :: fname
        integer(int32), parameter :: EXPECTED_BITS(12) = [16384_int32,49152_int32,15360_int32,48128_int32, &
            &14336_int32,13653_int32,1024_int32,31743_int32,15360_int32,15361_int32,0_int32,32768_int32]
        integer(kind=2) :: plane(BOX,BOX)
        integer(int32) :: bits
        integer :: funit, io_stat, ipixel
        open(newunit=funit,file=trim(fname%to_char()),access='stream',form='unformatted', &
            &action='read',status='old',iostat=io_stat)
        if( io_stat /= 0 ) THROW_HARD('failed to open float16 MRC payload')
        read(unit=funit,pos=1025,iostat=io_stat) plane
        close(funit)
        if( io_stat /= 0 ) THROW_HARD('failed to read float16 MRC payload')
        do ipixel = 1,size(EXPECTED_BITS)
            bits = iand(int(plane(ipixel,1),int32),65535_int32)
            if( bits /= EXPECTED_BITS(ipixel) ) THROW_HARD('float16 MRC payload bits mismatch')
        enddo
    end subroutine assert_float16_payload

    subroutine assert_float16_boundary_payload(fname)
        type(string), intent(in) :: fname
        integer(int32), parameter :: EXPECTED_BITS(8) = [1_int32,32769_int32,256_int32,33024_int32, &
            &0_int32,32768_int32,2_int32,32770_int32]
        integer(kind=2) :: plane(BOX,BOX)
        integer(int32) :: bits
        integer :: funit, io_stat, ipixel
        open(newunit=funit,file=trim(fname%to_char()),access='stream',form='unformatted', &
            &action='read',status='old',iostat=io_stat)
        if( io_stat /= 0 ) THROW_HARD('failed to open float16 boundary payload')
        read(unit=funit,pos=1025,iostat=io_stat) plane
        close(funit)
        if( io_stat /= 0 ) THROW_HARD('failed to read float16 boundary payload')
        do ipixel = 1,size(EXPECTED_BITS)
            bits = iand(int(plane(ipixel,1),int32),65535_int32)
            if( bits /= EXPECTED_BITS(ipixel) ) THROW_HARD('float16 boundary payload bits mismatch')
        enddo
    end subroutine assert_float16_boundary_payload

    subroutine write_int16_stack(stkname, istk)
        type(string), intent(in) :: stkname
        integer,      intent(in) :: istk
        type(MrcImgHead) :: header
        integer(int16), allocatable :: plane(:,:)
        integer :: funit, io_stat, iimg
        integer(kind=8) :: first_byte, image_nbytes
        allocate(plane(BOX,BOX))
        call header%new([BOX,BOX,NIMGS])
        call header%setMode(1)
        call header%setPixSz(SMPD)
        call header%setMinPixVal(real(100*istk + 1))
        call header%setMaxPixVal(real(100*istk + NIMGS))
        call header%setMean(real(100*istk) + real(NIMGS + 1) / 2.)
        open(newunit=funit, file=trim(stkname%to_char()), access='stream', form='unformatted', &
            &action='readwrite', status='replace', iostat=io_stat)
        if( io_stat /= 0 ) THROW_HARD('failed to open 16-bit test stack')
        call header%write(funit)
        first_byte   = int(header%firstDataByte(),kind=8)
        image_nbytes = int(BOX * BOX * 2,kind=8)
        do iimg = 1,NIMGS
            plane = int(100*istk + iimg, int16)
            write(unit=funit, pos=first_byte + int(iimg - 1,kind=8) * image_nbytes, iostat=io_stat) plane
            if( io_stat /= 0 ) THROW_HARD('failed to write 16-bit test stack')
        enddo
        close(funit)
        call header%kill
    end subroutine write_int16_stack

end program simple_test_discrete_stack_io
