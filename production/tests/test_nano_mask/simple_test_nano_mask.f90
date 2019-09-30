program simple_test_nano_mask
include 'simple_lib.f08'
use simple_image
use simple_segmentation
implicit none

character(len=*), parameter :: STK='selected.spi'
real,             parameter :: SMPD=0.358
integer,          parameter :: NGROW=4, WINSZ=3, EDGE=20

type(image) :: img, img_bin, img_inv_bin, img_copy
integer     :: ldim(3), nptcls, iptcl
real        :: thresh(3)

call find_ldim_nptcls(STK, ldim, nptcls)
call img%new(ldim, SMPD)
do iptcl = 1, nptcls
    call img%read(STK, iptcl)
    call img_copy%copy(img)
    call img_bin%copy(img)
    call img_bin%NLmean
    call otsu_robust_fast(img_bin, is2D=.true., noneg=.false., thresh=thresh)
    call img_bin%grow_bins(NGROW)
    call img_bin%real_space_filter(WINSZ, 'median')
    call img_bin%write('binarised.mrc', iptcl)
    call img_bin%cos_edge(EDGE)
    call img_bin%write('masks.mrc', iptcl)
    call img%remove_neg
    call img%mul(img_bin)
    call img%write('automasked.mrc', iptcl)
end do
call img%kill
end program simple_test_nano_mask
