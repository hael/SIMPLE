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
type(binimage)           :: microg, img_sdevs, img_cc, micraw
type(image)              :: img_win
type(tvfilter)           :: tvf
type(parameters), target :: params
type(stats_struct)       :: sz_stats, diam_stats
character(len=STDLEN),     parameter   :: filetable='filetab.txt'
character(len=LONGSTRLEN), allocatable :: micname(:)
character(len=STDLEN) :: outputfile, fbody
real, parameter       :: smpd=0.732, LAMBDA=3
integer               :: nfiles, ldim(3), ifoo, i, j, nboxes, sz_min, sz_max, box_raw, boxcoord(2)
integer, parameter    :: winsz=24 ! this is the windows size that works for the cryoEM polymeric-NPs
integer, allocatable  :: sz(:)
real,    allocatable  :: diams(:), masscens(:,:)
real :: px(3)
logical :: outside
params_glob => params
params_glob%pcontrast = 'black'
params_glob%hp        = 10.
call read_filetable(filetable, micname)
nfiles=size(micname)

do i = 1, nfiles
    print *, trim(micname(i))
    call find_ldim_nptcls(micname(i), ldim, ifoo)
    call microg%new_bimg(ldim, smpd)
    call microg%read(micname(i))
    call micraw%new_bimg(ldim, smpd)
    call micraw%read(micname(i))
    fbody      = get_fbody(basename(micname(i)),'mrc')
    outputfile = trim(fbody)//'.mrc'
    call microg%write(outputfile)
    ! low pass filter micrograph
    call microg%fft()
    call microg%bp(0., params_glob%hp)
    call microg%ifft()
    ! TV denoising
    call tvf%new()
    call tvf%apply_filter(microg, LAMBDA)
    call tvf%kill()
    outputfile = trim(fbody)//'_HP_FIL.mrc'
    print *, ">>> FILTERING HIGH FREQUENCIES   ",trim(outputfile)
    call microg%write(outputfile)
    call sauvola(microg, winsz, img_sdevs)
    outputfile = trim(fbody)//'_SAU_SDEVS.mrc'
    call img_sdevs%write(outputfile)
    outputfile = trim(fbody)//'_BIN_SAU.mrc'
    call microg%write(outputfile)
    ! binarize sdevs image with Otsu
    call otsu_img(img_sdevs)
    call microg%transfer2bimg(img_sdevs)
    print *, ">>> SAUVOL BINARISATION          ",trim(outputfile)
    outputfile = trim(fbody)//'_BIN_SAU_SDEVS.mrc'
    call microg%write_bimg(outputfile)
    call microg%erode()
    call microg%erode()
    call microg%set_largestcc2background()
    call microg%inv_bimg()
    outputfile=trim(fbody)//'_BIN_SAU_SDEVS_FH.mrc'
    call microg%write_bimg(outputfile)

    ! identify connected components
    call microg%find_ccs(img_cc)
    call img_cc%get_nccs(nboxes)
    ! eliminate connected components that are too large or too small
    sz = img_cc%size_ccs()
    call calc_stats(real(sz), sz_stats)
    print *, 'nboxes before elimination: ', nboxes
    print *, 'avg size: ', sz_stats%avg
    print *, 'med size: ', sz_stats%med
    print *, 'sde size: ', sz_stats%sdev
    print *, 'min size: ', sz_stats%minv
    print *, 'max size: ', sz_stats%maxv
    sz_min = nint(sz_stats%avg - params_glob%ndev * sz_stats%sdev)
    sz_max = nint(sz_stats%avg + params_glob%ndev * sz_stats%sdev)
    call img_cc%elim_ccs([sz_min,sz_max])
    call img_cc%get_nccs(nboxes)
    sz = img_cc%size_ccs()
    call calc_stats(real(sz), sz_stats)
    print *, 'nboxes after elimination: ', nboxes
    print *, 'avg size: ', sz_stats%avg
    print *, 'med size: ', sz_stats%med
    print *, 'sde size: ', sz_stats%sdev
    print *, 'min size: ', sz_stats%minv
    print *, 'max size: ', sz_stats%maxv
    allocate(diams(nboxes), source=0.)
    do j = 1, nboxes
        call img_cc%diameter_cc(j, diams(j))
    enddo
    call calc_stats(diams, diam_stats)
    print *, 'CC diameter (in Angs) stadistics'
    print *, 'avg diam: ', diam_stats%avg
    print *, 'med diam: ', diam_stats%med
    print *, 'sde diam: ', diam_stats%sdev
    print *, 'min diam: ', diam_stats%minv
    print *, 'max diam: ', diam_stats%maxv
    box_raw = find_magic_box(2 * nint(diam_stats%med/smpd))
    call img_win%new([box_raw,box_raw,1], smpd)
    if( allocated(masscens) ) deallocate(masscens)
    allocate(masscens(nboxes,2), source=0.)
    do j = 1, nboxes
        px = center_mass_cc(j)
        masscens(j,:2) = px(:2)
        boxcoord = nint((real(masscens(j,:2))-real(box_raw)/2.))
        call micraw%window_slim(boxcoord, box_raw, img_win, outside)
        call img_win%write(trim(fbody)//'_extracted.mrc', j)
        outputfile=trim(fbody)//'_extracted.mrc'
    enddo

    if( allocated(sz) ) deallocate(sz)
    if( allocated(diams) ) deallocate(diams)
    call img_sdevs%kill()
    call img_win%kill()
    call micraw%kill()
    call img_cc%kill_bimg()
enddo 

do i = 1, nfiles
    print *, trim(micname(i))
    call find_ldim_nptcls(micname(i), ldim, ifoo)
    call microg%new_bimg(ldim, smpd)
    call microg%read(micname(i))
    fbody      = get_fbody(basename(micname(i)),'mrc')
    outputfile = trim(fbody)//'_SBACK.mrc'
    print *, ">>> SUBSTRACTING BACKGROUND  ",trim(outputfile)
    call microg%subtract_background(params_glob%hp)
    call microg%write(outputfile)
    call microg%kill()
enddo

do i = 1, nfiles
    print *, trim(micname(i))
    call find_ldim_nptcls(micname(i), ldim, ifoo)
    call microg%new_bimg(ldim, smpd)
    call microg%read(micname(i))
    ! low pass filter micrograph
    call microg%fft()
    call microg%bp(0., params_glob%hp)
    call microg%ifft()
    ! TV denoising
    call tvf%new()
    call tvf%apply_filter(microg, LAMBDA)
    call tvf%kill()
    fbody      = get_fbody(basename(micname(i)),'mrc')
    outputfile = trim(fbody)//'_filt.mrc'
    print *, ">>> FILTERING HIGH FREQUENCIES   ",trim(outputfile)
    call microg%write(outputfile)
    outputfile = trim(fbody)//'_BIN_OTS.mrc'
    print *, ">>> OTSU'S BINARISATION   ",trim(outputfile)
    call otsu_img(microg,tight=.true.)
    !call otsu_img(microg)
    call microg%write(outputfile)
    call microg%kill_bimg()
enddo

contains

function center_mass_cc( i_cc ) result( px )
    integer, intent(in) :: i_cc
    real :: px(3)
    integer, allocatable :: pos(:,:)
    integer, allocatable :: imat_cc(:,:,:)
    imat_cc = int(img_cc%get_rmat())
    where( imat_cc .ne. i_cc ) imat_cc = 0
    call get_pixel_pos(imat_cc,pos)
    px(1) = sum(pos(1,:))/real(size(pos,dim = 2))
    px(2) = sum(pos(2,:))/real(size(pos,dim = 2))
    px(3) = 1.
    if( allocated(imat_cc) ) deallocate(imat_cc)
end function center_mass_cc

end program simple_test_picker