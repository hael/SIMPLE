program simple_test_discrete_stack_io
use, intrinsic :: iso_fortran_env, only: int16
use simple_core_module_api
use simple_discrete_stack_io, only: dstack_io
use simple_image,             only: image
use simple_imghead,           only: find_ldim_nptcls, MrcImgHead
implicit none
#include "simple_local_flags.inc"

integer, parameter :: NSTKS = 12, NIMGS = 4, BOX = 32, OPEN_WINDOW = 3
real,    parameter :: SMPD = 1.0
type(dstack_io) :: dstkios(OPEN_WINDOW)
type(image)     :: img, read_imgs(NSTKS,NIMGS)
type(string)    :: stknames(NSTKS)
real(kind=c_float), pointer :: rmat(:,:,:) => null()
integer :: istk, iimg, ldim(3), nptcls, stk_from, stk_to, iopen, nopen
real    :: expected

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

call simple_end('**** SIMPLE_TEST_DISCRETE_STACK_IO NORMAL STOP ****')

contains

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
