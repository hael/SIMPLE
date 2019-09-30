program simple_test_nano_mask
include 'simple_lib.f08'
use simple_image
use simple_segmentation
implicit none

character(len=*), parameter :: STK='selected.spi'
real,             parameter :: SMPD=0.358
integer,          parameter :: NGROW=3, WINSZ=1, EDGE=12

real, allocatable :: ccsizes(:), diams(:)
type(image)       :: img, img_bin, img_inv_bin, cc_img
integer           :: ldim(3), nptcls, iptcl, loc(1)
real              :: thresh(3)

call find_ldim_nptcls(STK, ldim, nptcls)
allocate(diams(nptcls), source=0.)
call img%new(ldim, SMPD)
do iptcl = 1, nptcls
    call progress(iptcl, nptcls)
    call img%read(STK, iptcl)
    call img_bin%copy(img)
    call img_bin%NLmean
    call otsu_robust_fast(img_bin, is2D=.true., noneg=.true., thresh=thresh)
    call img_bin%write('binarised_otsu.mrc', iptcl)
    call img_bin%grow_bins(NGROW)
    ! find the largest connected component
    call img_bin%find_connected_comps(cc_img)
    ccsizes = cc_img%size_connected_comps()
    loc     = maxloc(ccsizes)
    ! estimate its diameter
    call cc_img%diameter_cc(loc(1), diams(iptcl))
    ! turn it into a binary image for mask creation
    call cc_img%cc2bin(loc(1))
    call cc_img%write('binarised_otsu_grown.mrc', iptcl)
    call cc_img%real_space_filter(WINSZ, 'median')
    call cc_img%write('binarised_otsu_grown_median.mrc', iptcl)
    call cc_img%cos_edge(EDGE)
    call cc_img%write('masks_otsu.mrc', iptcl)
    call img%remove_neg
    call img%mul(cc_img)
    call img%write('automasked_otsu.mrc', iptcl)
end do
call img%kill
end program simple_test_nano_mask
