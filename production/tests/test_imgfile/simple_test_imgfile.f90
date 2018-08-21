program simple_test_imgfile
include 'simple_lib.f08'
use simple_image,   only: image
use simple_imgfile, only: imgfile
use simple_imghead
implicit none
#include "simple_local_flags.inc"
integer       :: ldim(3), i, j, cnt
real          :: smpd=2., corr, corrs(20)
type(image)   :: img, img_2
type(image)   :: imgs(20)
logical       :: ft=.false.

! SELF-CONSISTENCY TESTS

! create a square
ldim = [120,120,1]
call img%new(ldim, smpd)
call img%square(20)
call img%write('squares_mrc.mrc',1)
! write stacks of 5 squares
do i=1,5
    if( ft ) call img%fft()
    call img%write('squares_spider.spi',i)
    call img%write('squares_mrc.mrc',i)
end do
! create a cube
ldim = [120,120,120]
call img%new(ldim, smpd)
call img%square(20)
! write volume files
do i=1,5
    if( ft ) call img%fft()
    call img%write('cube_spider.spi')
    call img%write('cube_mrc.mrc')
end do
! convert the cubes from SPIDER to MRC & vice versa
do i=1,5
    call img%read('cube_spider.spi')
    if( ft ) call img%ifft()
    call img%write('cube_spider_converted.mrc')
    call img%read('cube_mrc.mrc')
    if( ft ) call img%ifft()
    call img%write('cube_mrc_converted.spi')
end do
! test SPIDER vs. MRC & converted vs. nonconverted
do i=1,4
    call imgs(i)%new(ldim, smpd)
    call imgs(i)%read('cube_spider.spi')
    call imgs(i)%read('cube_spider_converted.mrc')
    call imgs(i)%read('cube_mrc.mrc')
    call imgs(i)%read('cube_mrc_converted.spi')
    if( ft ) call imgs(i)%ifft()
end do
do i=1,3
    do j=i+1,4
        corr = imgs(i)%corr(imgs(j))
        if( corr < 0.99999 )then
            THROW_HARD('SPIDER vs. MRC & converted vs. nonconverted test failed')
        endif
    end do
end do
end program
