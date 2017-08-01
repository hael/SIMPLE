module simple_polarft_tester
use simple_image,     only: image
use simple_polarft,   only: polarft
use simple_projector, only: projector 
use simple_defs       ! singleton
implicit none

public :: init_polarft_tester, run_polarft_tests, test_polarft_shift, test_polarft_corrcalc_shift
private

type(polarft)      :: polft
type(polarft)      :: polft2rot
type(projector)    :: proj
type(image)        :: img, img2rot, img_noisy, img_copy, img_pad
real               :: angerr, corr, avg_angerr, avg_corr, msk
real               :: avgerr_of_SNR, avgcorr_of_SNR
integer, parameter :: ntests=100
integer            :: kfromto(2)

contains

    subroutine init_polarft_tester( img2test, mskin ) 
        type(image), intent(in) :: img2test
        real, intent(in)        :: mskin
        integer                 :: ldim(3)
        ! set constants
        ldim       = img2test%get_ldim()
        kfromto(1) = 2
        msk        = mskin
        ! store input image
        img        = img2test
        ldim = img%get_ldim()
        call img_pad%new([2*ldim(1),2*ldim(2),1], img%get_smpd())
        img_copy   = img      
        proj       = projector()
    end subroutine

    subroutine run_polarft_tests( lp20e, lp15e, lp10e, lp7e, snr02e, snr01e, snr005e, snr001e )
        real, intent(out) :: lp20e,lp15e,lp10e,lp7e,snr02e,snr01e,snr005e,snr001e
        avgerr_of_SNR  = 0.
        avgcorr_of_SNR = 0.
        lp20e          = 0.
        lp15e          = 0.
        lp10e          = 0.
        lp7e           = 0.
        snr02e         = 0.
        snr01e         = 0.
        snr005e        = 0.
        snr001e        = 0.
        ! snr=0.2
        call polarft_tester(20.,0.2,msk)
        lp20e = lp20e+avg_angerr
        call polarft_tester(15.,0.2,msk)
        lp15e = lp15e+avg_angerr
        call polarft_tester(10.,0.2,msk)
        lp10e = lp10e+avg_angerr
        call polarft_tester( 7.,0.2,msk)
        lp7e = lp7e+avg_angerr
        snr02e = avgerr_of_SNR/4.
        ! snr=0.1
        avgerr_of_SNR = 0.
        avgcorr_of_SNR = 0.
        call polarft_tester(20.,0.1,msk)
        lp20e = lp20e+avg_angerr
        call polarft_tester(15.,0.1,msk)
        lp15e = lp15e+avg_angerr
        call polarft_tester(10.,0.1,msk)
        lp10e = lp10e+avg_angerr
        call polarft_tester( 7.,0.1,msk)
        lp7e = lp7e+avg_angerr
        snr01e = avgerr_of_SNR/4.
        ! snr=0.05
        avgerr_of_SNR = 0.
        avgcorr_of_SNR = 0.
        call polarft_tester(20.,0.05,msk)
        lp20e = lp20e+avg_angerr
        call polarft_tester(15.,0.05,msk)
        lp15e = lp15e+avg_angerr
        call polarft_tester(10.,0.05,msk)
        lp10e = lp10e+avg_angerr
        call polarft_tester( 7.,0.05,msk)
        lp7e = lp7e+avg_angerr
        snr005e = avgerr_of_SNR/4.
        ! snr=0.01
        avgerr_of_SNR = 0.
        avgcorr_of_SNR = 0.
        call polarft_tester(20.,0.01,msk)
        lp20e = lp20e+avg_angerr
        call polarft_tester(15.,0.01,msk)
        lp15e = lp15e+avg_angerr
        call polarft_tester(10.,0.01,msk)
        lp10e = lp10e+avg_angerr
        call polarft_tester( 7.,0.01,msk)
        lp7e = lp7e+avg_angerr
        snr001e = avgerr_of_SNR/4.
        lp20e = lp20e/4.
        lp15e = lp15e/4.
        lp10e = lp10e/4.
        lp7e  = lp7e/4.
    end subroutine
    
    !>  \brief  testing w respect to snr
    subroutine polarft_tester( lp, snr, msk )
        real, intent(in) :: lp, snr, msk
        integer :: i
        kfromto(2) = img%get_find(lp)
        call polft%new(    kfromto, kfromto(2))
        call polft2rot%new(kfromto, kfromto(2))
        avg_angerr = 0.
        avg_corr   = 0.
        do i=1,ntests
            img_noisy = img
            call img_noisy%add_gauran(snr)
            call img%mask(msk,'soft')
            call img_noisy%mask(msk,'soft')
            call test_polarft(img, img_noisy, msk, lp, angerr, corr)
            avg_angerr = avg_angerr+angerr
            avg_corr   = avg_corr+corr
            ! put original image back
            img = img_copy
        end do
        avg_angerr    = avg_angerr/real(ntests)
        avg_corr      = avg_corr/real(ntests)
        avgerr_of_SNR = avgerr_of_SNR+avg_angerr
    end subroutine

    !>  \brief  polarft unit test
    subroutine test_polarft( img, img2rot, msk, lp, angerr, corr )
        use simple_rnd,   only: ran3
        class(image), intent(inout) :: img, img2rot
        real, intent(in)            :: msk, lp
        real, intent(out)           :: angerr, corr
        type(image)                 :: img_rot, img2
        real                        :: ranang, truang, ang
        integer                     :: ldim(3)
        ldim = img%get_ldim()
        call img_rot%new(ldim, img%get_smpd())
        call img2%new(   ldim, img%get_smpd())
        ranang = ran3()*360.
        truang = 360.-ranang
        call img2rot%rtsq(ranang,0.,0.,img_rot)       ! rotating the image into img_rot
        call proj%img2polarft(img,polft)              ! making the polarft
        call proj%img2polarft(img_rot,polft2rot)      ! making the rotated polarft
        call polft%rotsrch_slow(polft2rot, ang, corr) ! search for an angle using img as ref
        call img_rot%rtsq(ang,0.,0.,img2)             ! rotate back into img2 using the angle found
        call img_rot%rtsq(truang,0.,0.,img)           ! rotating back into img using the correct angle
        angerr = abs(ang-truang)                      ! this is the absolute deviation (angular error)
        corr = img2%corr(img, lp)                     ! this is the correlation
    end subroutine
    
    subroutine test_polarft_shift
        type(image)   :: img, img_shifted
        type(polarft) :: pimg, pimg_shifted
        integer, parameter :: trs=5
        integer :: x, y
        write(*,'(a)') '**info(simple_polarft_tester): testing polarft shift'
        proj = projector()
        call img%new([100,100,1], 2.)
        call img_shifted%new([100,100,1], 2.)
        call img%square(20)
        call pimg%new([2,11], 45, ptcl=.false.)
        call img%fwd_ft
        call proj%img2polarft(img, pimg)
        call pimg_shifted%new([2,11], 45, ptcl=.false.)
        call pimg_shifted%set_ldim([100,100,1])
        do x=-trs,trs
            do y=-trs,trs
                call img%shift(real(x), real(y), imgout=img_shifted)
                call proj%img2polarft(img_shifted, pimg_shifted)
                corr = pimg%corr_shifted(pimg_shifted, 1, [real(x),real(y)])
                if( corr < 0.99 )then
                    print *, 'corr = ', corr
                    stop 'failure! first part of test_polarft_shift; simple_polarft_tester'
                endif
                call pimg%shift([real(x),real(y)])
                corr = pimg%corr_shifted(pimg_shifted, 1)
                if( corr < 0.99 )then
                    print *, 'corr = ', corr
                    stop 'failure! second part of test_polarft_shift; simple_polarft_tester'
                endif
            end do
        end do
        write(*,'(a)') 'TEST_POLARFT_SHIFT COMPLETED SUCCESSFULLY ;-)'         
    end subroutine
    
    subroutine test_polarft_corrcalc_shift
        use simple_polarft_corrcalc, only: polarft_corrcalc
        type(image)            :: img, img_shifted
        type(polarft)          :: pimg, pimg_shifted
        type(polarft_corrcalc) :: pftcc
        integer, parameter     :: trs=2
        integer                :: x, y, r
        real                   :: shvec(2), cc
        write(*,'(a)') '**info(test_polarft_corrcalc_shift, part 1): testing the shifting consistency by correlation'
        proj = projector()
        call img%new([100,100,1], 2.)
        call img_shifted%new([100,100,1], 2.)
        call img%square(20)
        call img%fwd_ft
        ! make a conforming polarft_corrcalc object of 3 refs and 5 ptcls
        call pftcc%new(3, [1,5], [100,100,1], [2,11], 45, 'no')
        call proj%img2polarft(1, img, pftcc, isptcl=.false.) ! this is the reference
        do x=-trs,trs
            do y=-trs,trs
                call img%shift(real(x), real(y), imgout=img_shifted)
                call proj%img2polarft(1, img_shifted, pftcc, isptcl=.true.) ! this is the particle (shifted)
                corr = pftcc%corr(1, 1, 1, [real(x),real(y)])
                if( corr < 0.99 )then
                    print *, 'corr = ', corr
                    stop 'failure! first part of test_polarft_corrcalc_shift; simple_polarft_tester'
                endif
            end do
        end do
        write(*,'(a)') 'TEST_POLARFT_CORRCALC_SHIFT COMPLETED SUCCESSFULLY ;-)'
    end subroutine

end module simple_polarft_tester
