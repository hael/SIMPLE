module simple_ssim
!    include 'simple_lib.f08'
use simple_defs
use simple_image, only: image
implicit none
private
#include "simple_local_flags.inc"
public :: ssim, dssim
contains
    function ssim (img1, img2, ssim_map, weighting)
        real :: ssim
        type(image), intent(in) :: img1,img2
        type(image), intent(inout), optional:: ssim_map
        integer, intent(in), optional :: weighting(3)
        real ::  C1 = 6.5025
        real ::  C2 = 58.5225
        real, allocatable:: temp1(:,:,:), temp2(:,:,:), temp3(:,:,:)
        type(image) :: img1_img2, img1_temp, img2_temp,&
            img1_sq,img2_sq, mu1, mu2, mu1_sq, mu2_sq, mu1_mu2, &
            sigma1_sq, sigma2_sq, sigma12, ssim_temp
        type(image) :: luminance, contrast, structure, sigma1, sigma2
        integer :: w(3)
        integer :: ldim(3) , nChan
        nChan = 1
        w = (/ 1, 1, 1 /)
        if(present(weighting)) w = weighting
        ldim = img1%get_ldim()

        img1_sq =  img1 *  img1
        img2_sq= img2 * img2
        img1_img2=img1 * img2

        call mu1%new([ldim(1), ldim(2), nChan ], smpd=1.)
        call mu2%new([ldim(1), ldim(2), nChan ], 1.)

        call mu1_sq%new([ldim(1), ldim(2), nChan ],1.)
        call mu2_sq%new([ldim(1), ldim(2), nChan ],1.)
        call mu1_mu2%new([ldim(1), ldim(2), nChan ],1.)

        call sigma1_sq%new([ldim(1), ldim(2), nChan ],1.)
        call sigma2_sq%new([ldim(1), ldim(2), nChan ],1.)
        call sigma12%new([ldim(1), ldim(2), nChan ],1.)

        allocate(temp1(ldim(1), ldim(2), nChan ))
        allocate( temp2(ldim(1), ldim(2), nChan ))
        allocate( temp3(ldim(1), ldim(2), nChan ))

        call  ssim_temp%new([ldim(1), ldim(2), nChan ],1.)


        mu1 =img1
        call mu1%real_space_filter(11, 'average')
        mu2 =img2
        call mu2%real_space_filter(11, 'average')

        mu1_sq=mu1*mu1
        mu2_sq = mu2*mu2
        mu1_mu2=mu1*mu2

        sigma1_sq=  img1_sq
        call sigma1_sq%real_space_filter(11, 'average')
        sigma1_sq = sigma1_sq - mu1_sq
        !cvSmooth( img1_sq, sigma1_sq, CV_GAUSSIAN, 11, 11, 1.5 );
        !cvAddWeighted( sigma1_sq, 1, mu1_sq, -1, 0, sigma1_sq );
        sigma2_sq=  img2_sq
        call sigma2_sq%real_space_filter(11, 'average')
        sigma2_sq = sigma2_sq - mu2_sq
        !cvSmooth( img2_sq, sigma2_sq, CV_GAUSSIAN, 11, 11, 1.5 );
        !!cvAddWeighted( sigma2_sq, 1, mu2_sq, -1, 0, sigma2_sq );
        sigma12=img1_img2
        call sigma12%real_space_filter(11, 'average')
        !cvSmooth( img1_img2, sigma12, CV_GAUSSIAN, 11, 11, 1.5 );
        !cvAddWeighted( sigma12, 1, mu1_mu2, -1, 0, sigma12 );
        sigma12 = sigma12 - mu1_mu2
#if 0

        call sigma1%new([ldim(1), ldim(2), nChan ])
        call sigma2%new([ldim(1), ldim(2), nChan ])
        call luminance%new([ldim(1), ldim(2), nChan ])
        call contrast%new([ldim(1), ldim(2), nChan ])
        call structure%new([ldim(1), ldim(2), nChan ])

        sigma1 = sigma1_sq%sqrt()
        sigma2 = sigma2_sq%sqrt()
        luminance = (2*mu1_mu2 + C1) / (mu1_sq + mu2_sq + C1)
        contrast = (2*sigma1*sigma2 + C2 ) /( sigma1_sq + sigma2_sq + C2)
        structure = (sigma12 + 0.5*C2) / (sigma1*sigma2 + 0.5*C2)
        ssim_temp=(luminance^w(1)) * (contrast^w(2)) * (structure^w(3))
#endif

        temp1 = 2.0 * mu1_mu2%get_rmat() + C1
        temp2 = 2.0 * sigma12%get_rmat() + C2
        temp3 = temp1 * temp2
        temp1 = (mu1_sq%get_rmat() + mu2_sq%get_rmat() +(C1))
        temp2 =(sigma1_sq%get_rmat() + sigma2_sq%get_rmat() + C2)
        ! ((mu1_sq + mu2_sq + C1).*(sigma1_sq + sigma2_sq + C2))
        temp1 = temp2 * temp1
        ! ((2*mu1_mu2 + C1).*(2*sigma12 + C2))./((mu1_sq + mu2_sq + C1).*(sigma1_sq + sigma2_sq + C2))
        call ssim_map%set_rmat( temp3 / temp1)
        ssim =  ssim_map%mean()

        ! through observation, there is approximately
        ! 1% error max with the original matlab program

        VerbosePrint "SSIM index"
        VerbosePrint ssim * 100 , "%"

        call img1_img2%kill
        call img1_temp%kill
        call img2_temp%kill
        call img1_sq%kill
        call img2_sq%kill
        call mu1%kill
        call mu2%kill
        call mu1_sq%kill
        call mu2_sq%kill
        call  mu1_mu2%kill
        call sigma1_sq%kill
        call sigma2_sq%kill
        call sigma12%kill
        deallocate(temp1,temp2, temp3)

        if(present(ssim_map)) ssim_map = ssim_temp
    end function ssim

    function dssim (img1, img2)
        real :: dssim
        type(image), intent(in) :: img1,img2

        dssim = (1-ssim(img1,img2))/2

    end function dssim


end module simple_ssim
