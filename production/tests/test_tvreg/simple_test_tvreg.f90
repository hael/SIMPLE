program simple_test_tvreg
include 'simple_lib.f08'
use simple_tvlam_opt, only: tvlam_opt
use simple_image,     only: image
implicit none

character(len=*), parameter :: stk_e  = 'cavgs_iter015_even.mrc'
character(len=*), parameter :: stk_o  = 'cavgs_iter015_odd.mrc'
character(len=*), parameter :: stk_f  = 'tvfiltered.mrc'
integer,          parameter :: box    = 104, ncls = 90
real,             parameter :: smpd   = 3.138
real,             parameter :: msk    = 200. / smpd
type(image)     :: img_e, img_o, img
integer         :: icls
real            :: lam
type(tvlam_opt) :: tvlamfind

call img_e%new([box,box,1], smpd) 
call img_o%new([box,box,1], smpd)
call img%new([box,box,1], smpd)
call tvlamfind%new([box,box,1], smpd, msk)
do icls = 1, ncls
    call img_e%read(stk_e, icls)
    call img_o%read(stk_o, icls)
    call tvlamfind%set_img_ptrs(img_e, img_o)
    call tvlamfind%minimize(lam)
    call tvlamfind%get_tvfiltered(img)
    call img%write(stk_f, icls)
end do

end program simple_test_tvreg
