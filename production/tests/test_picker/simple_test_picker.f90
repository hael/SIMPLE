program simple_test_picker
include 'simple_lib.f08'
use simple_parameters
use simple_pickseg
use simple_segmentation
use simple_image,    only: image
use simple_binimage, only: binimage
use simple_tvfilter, only: tvfilter
implicit none
#include "simple_local_flags.inc"
type(image)              :: microg, img_sdevs
type(binimage)           :: micbin 
type(tvfilter)           :: tvf
type(parameters), target :: params
character(len=STDLEN),     parameter   :: filetable='filetab.txt'
character(len=LONGSTRLEN), allocatable :: micname(:)
character(len=STDLEN) :: outputfile, fbody
real, parameter       :: smpd=0.732, LAMBDA=3
integer               :: nfiles, ldim(3), ifoo, i
integer, parameter    :: winsz=24 ! this is the windows size that works for the cryoEM polymeric-NPs
params_glob => params
params_glob%pcontrast = 'black'
params_glob%lp        = 10.
params_glob%nsig      = 1.5
call read_filetable(filetable, micname)
nfiles=size(micname)

do i = 1, nfiles
    print *, trim(micname(i))
    call find_ldim_nptcls(micname(i), ldim, ifoo)
    call microg%new(ldim, smpd)
    call micbin%new_bimg(ldim, smpd)
    call microg%read(micname(i))
    fbody      = get_fbody(basename(micname(i)),'mrc')
    outputfile = trim(fbody)//'.mrc'
    call microg%write(outputfile)
    ! low pass filter micrograph
    call microg%fft()
    call microg%bp(0., params_glob%lp)
    call microg%ifft()
    ! TV denoising
    call tvf%new()
    call tvf%apply_filter(microg, LAMBDA)
    call tvf%kill()
    outputfile = trim(fbody)//'_HP_FIL.mrc'
    print *, ">>> FILTERING HIGH FREQUENCIES   ",trim(outputfile)
    call microg%write(outputfile)
    call sauvola(microg, winsz, img_sdevs)
    outputfile=trim(fbody)//'_SAU_SDEVS.mrc'
    call img_sdevs%write(outputfile)
    ! binarize sdevs image with Otsu
    call otsu_img(img_sdevs)
    call micbin%transfer2bimg(img_sdevs)
    print *, ">>> SAUVOL BINARISATION          ",trim(outputfile)
    outputfile = trim(fbody)//'_BIN_SAU_SDEVS.mrc'
    call micbin%write_bimg(outputfile)
    outputfile = trim(fbody)//'_BIN_SAU.mrc'
    call microg%write(outputfile)
    call micbin%erode()
    call micbin%erode()
    call micbin%fill_holes()
    outputfile=trim(fbody)//'_BIN_SAU_SDEVS_FH.mrc'
    call micbin%write_bimg(outputfile)
    call microg%kill()
    call img_sdevs%kill()
    call micbin%kill_bimg()
enddo 

do i = 1, nfiles
    print *, trim(micname(i))
    call find_ldim_nptcls(micname(i), ldim, ifoo)
    call microg%new(ldim, smpd)
    call microg%read(micname(i))
    fbody      = get_fbody(basename(micname(i)),'mrc')
    outputfile = trim(fbody)//'_SBACK.mrc'
    print *, ">>> SUBSTRACTING BACKGROUND  ",trim(outputfile)
    call microg%subtract_background(params_glob%lp)
    call microg%write(outputfile)
    call microg%kill()
enddo

do i = 1, nfiles
    print *, trim(micname(i))
    call find_ldim_nptcls(micname(i), ldim, ifoo)
    call microg%new(ldim, smpd)
    call micbin%new_bimg(ldim, smpd)
    call microg%read(micname(i))
    ! low pass filter micrograph
    call microg%fft()
    call microg%bp(0., params_glob%lp)
    call microg%ifft()
    ! TV denoising
    call tvf%new()
    call tvf%apply_filter(microg, LAMBDA)
    call tvf%kill()
    fbody=get_fbody(basename(micname(i)),'mrc')
    outputfile = trim(fbody)//'_filt.mrc'
    print *, ">>> FILTERING HIGH FREQUENCIES   ",trim(outputfile)
    call microg%write(outputfile)
    outputfile = trim(fbody)//'_BIN_OTS.mrc'
    print *, ">>> OTSU'S BINARISATION   ",trim(outputfile)
    call otsu_img(microg,tight=.true.)
    !call otsu_img(microg)
    call micbin%transfer2bimg(microg)
    call micbin%write_bimg(outputfile)
    call microg%kill()
    call micbin%kill_bimg()
enddo 

end program simple_test_picker