program simple_test_cartcorr_sanity
use simple_image,       only: image
use simple_ft_expanded, only: ft_expanded
use simple_rnd,         only: ran3
use simple_defs
implicit none
integer, parameter   :: NITS=100
real, parameter      :: TRS=5.0
logical, parameter   :: debug=.false.
type(image)          :: img1, img2, tmp
type(ft_expanded)    :: ftexp1, ftexp2
integer              :: i
real, parameter      :: hp=100.0, lp=8.0
real                 :: cxy(3), x, y, diff_old, diff_recast, diff_corr, old_corr
! make a square for testing the shift alignment
call img1%new([200,200,1], 1.77)
call img1%square(40)
img2 = img1
tmp  = img1
! prepare objects for corrcalc
call ftexp1%new(img1,hp,lp)
! test
diff_old    = 0.
diff_recast = 0.
diff_corr   = 0. 
do i=1,nits
    if( debug ) print *, 'iteration: ', i
    x           = real(nint(ran3()*2.0*TRS-TRS))
    y           = real(nint(ran3()*2.0*TRS-TRS))
    if( debug ) print *, 'shifted:     ', x, y
    call img1%shift(x,y,imgout=img2)  
    cxy         = find_shift_old()
    old_corr    = cxy(1)
    if( debug ) print *, 'neg(old):    ', -cxy(2:3)
    diff_old    = diff_old+sum(abs(-cxy(2:3)-[x,y]))
    call ftexp2%new(img2,hp,lp)
    cxy         = find_shift_recast()
    if( debug ) print *, 'neg(recast): ', -cxy(2:3)
    diff_recast = diff_recast+sum(abs(-cxy(2:3)-[x,y]))
    diff_corr   = diff_corr+abs(old_corr-cxy(1))
    
end do
write(*,*) 'diff_old:    ', diff_old
write(*,*) 'diff_recast: ', diff_recast
write(*,*) 'diff_corr:   ', diff_corr
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
    end function
    
    function find_shift_recast() result( cxy )
        integer :: sh,xsh,ysh
        real :: cxy(3), corr, shvec(3)
        sh = nint(TRS)
        cxy(1) = -1.
        do xsh=-sh,sh
            do ysh=-sh,sh
                shvec = [real(xsh),real(ysh),0.]
                corr = ftexp1%corr_shifted(ftexp2,shvec)
                if( corr > cxy(1) )then
                    cxy(1) = corr
                    cxy(2) = real(xsh)
                    cxy(3) = real(ysh)
                endif
            end do
        end do
    end function

end program simple_test_cartcorr_sanity