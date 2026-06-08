program simple_test_discrete_stack_io
use simple_core_module_api
use simple_discrete_stack_io, only: dstack_io
use simple_image,             only: image
use simple_imghead,           only: find_ldim_nptcls
implicit none
#include "simple_local_flags.inc"

integer, parameter :: NSTKS = 4, NIMGS = 6, BOX = 32
real,    parameter :: SMPD = 1.0
type(dstack_io) :: dstkios(NSTKS)
type(image)     :: img, read_imgs(NSTKS,NIMGS)
type(string)    :: stknames(NSTKS)
real(kind=c_float), pointer :: rmat(:,:,:) => null()
integer :: istk, iimg, ldim(3), nptcls
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

do istk = 1,NSTKS
    call find_ldim_nptcls(stknames(istk), ldim, nptcls)
    call dstkios(istk)%new(SMPD, BOX)
    call dstkios(istk)%cache_stack_info(stknames(istk), ldim, nptcls)
    call dstkios(istk)%open(stknames(istk))
    do iimg = 1,NIMGS
        call read_imgs(istk,iimg)%new([BOX,BOX,1], SMPD, wthreads=.false.)
    enddo
enddo

!$omp parallel do default(shared) private(istk,iimg) schedule(static) proc_bind(close) num_threads(NSTKS)
do istk = 1,NSTKS
    do iimg = 1,NIMGS
        call dstkios(istk)%read(stknames(istk), iimg, read_imgs(istk,iimg))
    enddo
enddo
!$omp end parallel do

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
    call dstkios(istk)%kill
    call del_file(stknames(istk))
    call stknames(istk)%kill
enddo

call simple_end('**** SIMPLE_TEST_DISCRETE_STACK_IO NORMAL STOP ****')

end program simple_test_discrete_stack_io
