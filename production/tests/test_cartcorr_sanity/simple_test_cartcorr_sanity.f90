program simple_test_cartcorr_sanity
include 'simple_lib.f08'
use simple_image,        only: image
use simple_ft_expanded,  only: ft_expanded
use simple_ftexp_shsrch, only: ftexp_shsrch
implicit none
#include "simple_local_flags.inc"

integer, parameter   :: NITS=100
real, parameter      :: TRS=5.0
type(image)          :: img1, img2, tmp
type(ft_expanded)    :: ftexp1, ftexp2
type(ftexp_shsrch)   :: ftexp_shsrch1
integer              :: i
real, parameter      :: hp=100.0, lp=8.0
real                 :: cxy(3), x, y, diff_old, diff_recast, diff_corr, old_corr
! make a square for testing the shift alignment
call img1%new([200,200,1], 1.77)
call img1%square(40)
call img1%fft()
img2 = img1
tmp  = img1
! prepare objects for corrcalc
call ftexp1%new(img1,hp,lp,.true.)
! test
diff_old    = 0.
diff_recast = 0.
diff_corr   = 0.
do i=1,nits
    DebugPrint  'iteration: ', i
    x           = real(nint(ran3()*2.0*TRS-TRS))
    y           = real(nint(ran3()*2.0*TRS-TRS))
    DebugPrint  'shifted:     ', x, y
    img2 = img1
    call img2%shift([x,y,0.])
    cxy         = find_shift_old()
    old_corr    = cxy(1)
    DebugPrint  'neg(old):    ', -cxy(2:3)
    diff_old    = diff_old+sum(abs(-cxy(2:3)-[x,y]))
    call ftexp2%new(img2,hp,lp,.true.)
    call ftexp_shsrch1%new(ftexp1,ftexp2,TRS)
    cxy         = find_shift_recast()
    DebugPrint  'neg(recast): ', -cxy(2:3)
    diff_recast = diff_recast+sum(abs(-cxy(2:3)-[x,y]))
    diff_corr   = diff_corr+abs(old_corr-cxy(1))
end do
write(logfhandle,*) 'diff_old:    ', diff_old
write(logfhandle,*) 'diff_recast: ', diff_recast
write(logfhandle,*) 'diff_corr:   ', diff_corr
contains

    function find_shift_old() result( cxy )
        integer :: sh,xsh,ysh
        real :: cxy(3), corr, shvec(3)
        sh = nint(TRS)
        cxy(1) = -1.
        do xsh=-sh,sh
            do ysh=-sh,sh
                shvec = [real(xsh),real(ysh),0.]
                corr = img1%corr_shifted(img2, shvec, lp_dyn=lp)
                if( corr > cxy(1) )then
                    cxy(1) = corr
                    cxy(2) = real(xsh)
                    cxy(3) = real(ysh)
                endif
            end do
        end do
    end function find_shift_old

    function find_shift_recast() result( cxy )
        integer :: sh,xsh,ysh
        real :: cxy(3), corr, shvec(3)
        sh = nint(TRS)
        cxy(1) = -1.
        do xsh=-sh,sh
            do ysh=-sh,sh
                shvec = [real(xsh),real(ysh),0.]
                corr = real(ftexp_shsrch1%corr_shifted_8(dble(shvec)))! ftexp1%corr_shifted_8(ftexp2,dble(shvec)))
                if( corr > cxy(1) )then
                    cxy(1) = corr
                    cxy(2) = real(xsh)
                    cxy(3) = real(ysh)
                endif
            end do
        end do
    end function find_shift_recast

end program simple_test_cartcorr_sanity
