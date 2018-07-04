module simple_ssim
include 'simple_lib.f08'
use simple_image, only: image
implicit none
private
#include "simple_local_flags.inc"
public :: ssim, dssim
contains
    function ssim (img1, img2, ssim_map, weighting, save_temps)
        !$ use omp_lib
        real :: ssim
        type(image), intent(in) :: img1,img2
        type(image), intent(inout), optional:: ssim_map
        integer, intent(in), optional :: weighting(3)
        logical, intent(in), optional :: save_temps
        real ::  C1 = 6.5025
        real ::  C2 = 58.5225
        real, allocatable:: temp1(:,:,:), temp2(:,:,:), temp3(:,:,:)
        type(image) :: img1_img2, img1_temp, img2_temp,&
            img1_sq,img2_sq, mu1, mu2, mu1_sq, mu2_sq, mu1_mu2, &
            sigma1_sq, sigma2_sq, sigma12, ssim_temp
        type(image) ::  contrast, luminance,  structure, sigma1, sigma2
        integer :: w(3), winsz
        integer :: ldim(3) , nChan
        verbose=.true.
        nChan = 1
        w = (/ 1, 1, 1 /)
        print *," In ssim "
        if(present(weighting)) w = weighting
        ldim = img1%get_ldim()
        print *," size of images ", ldim
        nChan=ldim(3)


        call mu1%new([ldim(1), ldim(2), nChan ], smpd=1.)

        call mu2%new([ldim(1), ldim(2), nChan ], 1.)

        call mu1_sq%new([ldim(1), ldim(2), nChan ],1.)

        call mu2_sq%new([ldim(1), ldim(2), nChan ],1.)

        call mu1_mu2%new([ldim(1), ldim(2), nChan ],1.)


        call sigma1_sq%new([ldim(1), ldim(2), nChan ],1.)

        call sigma2_sq%new([ldim(1), ldim(2), nChan ],1.)

        call sigma12%new([ldim(1), ldim(2), nChan ],1.)

        allocate( temp1(ldim(1), ldim(2), nChan ))

        allocate( temp2(ldim(1), ldim(2), nChan ))

        allocate( temp3(ldim(1), ldim(2), nChan ))

        call  ssim_temp%new([ldim(1), ldim(2), nChan ],1.)

        VerbosePrint "SSIM allocated "
        winsz = 3

        ! omp parallel workshare
        ! omp section
        img1_sq =  img1%pow(2)
 VerbosePrint "SSIM img1^2 "
        ! omp section
 img2_sq = img2%pow(2)
  VerbosePrint "SSIM img2^2 "
        ! omp section
        mu1 =img1
        call mu1%real_space_filter(winsz, 'average')
         VerbosePrint "SSIM mu1 filter "
        mu1_sq= mu1 * mu1
         VerbosePrint "SSIM mu1 "
        ! omp section
        mu2 = img2
        call mu2%real_space_filter(winsz, 'average')
        VerbosePrint "SSIM mu2 filter "
        mu2_sq = mu2%pow(2)
        ! omp end parallel workshare

        VerbosePrint "SSIM mu2 "
        !omp parallel sections
        !omp section
        img1_img2 = img1 * img2
        !omp section
        mu1_mu2= mu1 * mu2
        !omp end parallel sections
        VerbosePrint "SSIM mu1mu2"

        ! omp parallel sections
        ! omp section
        sigma1_sq=  img1_sq
        call sigma1_sq%real_space_filter(winsz, 'average')
        VerbosePrint "SSIM sigma1_sq filter "
        sigma1_sq = sigma1_sq - mu1_sq
        ! omp section
        sigma2_sq=  img2_sq
        call sigma2_sq%real_space_filter(winsz, 'average')
        VerbosePrint "SSIM sigma2_sq filter "
        sigma2_sq = sigma2_sq - mu2_sq
        ! omp section
        sigma12=img1_img2
        call sigma12%real_space_filter(winsz, 'average')
        VerbosePrint "SSIM sigma12 filter "
        sigma12 = sigma12 - mu1_mu2
        ! omp end parallel sections
        VerbosePrint "SSIM sigma "

        if(present(save_temps))then

            call sigma1%new([ldim(1), ldim(2), nChan ],1.)
            call sigma2%new([ldim(1), ldim(2), nChan ],1.)
            call luminance%new([ldim(1), ldim(2), nChan ],1.)
            call contrast%new([ldim(1), ldim(2), nChan ],1.)
            call structure%new([ldim(1), ldim(2), nChan ],1.)
            call mu1_mu2%mul(2.)
            sigma1 = sigma1_sq%sq_root()
            sigma2 = sigma2_sq%sq_root()
            luminance = (mu1_mu2 + C1) / ( (mu1_sq + mu2_sq) + C1)
            contrast = ((sigma1*sigma2)*2. + C2 ) / ( sigma1_sq + sigma2_sq + C2)
            structure = (sigma12 + 0.5*C2) / (sigma1*sigma2 + 0.5*C2)
            ssim_temp= luminance%pow(w(1)) * contrast%pow(w(2)) * structure%pow(w(3))

            call luminance%write('luminance.mrc')
            call contrast%write('contrast.mrc')
            call structure%write('structure.mrc')
            call ssim_temp%write('ssim_classical.mrc')
            call sigma1%kill
            call sigma2%kill
            call luminance%kill
            call contrast%kill
            call structure%kill
            call mu1_mu2%kill

        end if


        !omp parallel sections
        !omp section
        temp1 = 2.0 * mu1_mu2%get_rmat() + C1
        !omp section
        temp2 = 2.0 * sigma12%get_rmat() + C2
        !omp end parallel sections

        !omp parallel workshare
        temp3 = temp1 * temp2
        !omp end parallel workshare



        VerbosePrint "SSIM temps"
        !omp parallel sections
        !omp section
        temp1 = (mu1_sq%get_rmat() + mu2_sq%get_rmat() +(C1))
        !omp section
        temp2 =(sigma1_sq%get_rmat() + sigma2_sq%get_rmat() + C2)
        ! ((mu1_sq + mu2_sq + C1).*(sigma1_sq + sigma2_sq + C2))
        !omp end parallel sections



        !omp parallel workshare
        temp1 = temp2 * temp1
        ! ((2*mu1_mu2 + C1).*(2*sigma12 + C2))./((mu1_sq + mu2_sq + C1).*(sigma1_sq + sigma2_sq + C2))
        !omp end parallel workshare


        call ssim_temp%set_rmat( temp3 / temp1 )
        ssim =  ssim_temp%mean()

        if(present(ssim_map)) call ssim_map%copy(ssim_temp)
        call ssim_temp%write('ssim_map.mrc')

        ! through observation, there is approximately
        ! 1% error max with the original matlab program

        VerbosePrint "SSIM index", ssim * 100 , "%"
        VerbosePrint "DSSIM index", 1-(ssim * 100)/2 , "%"

        call ssim_temp%kill
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
    end function ssim

    function ssim2 (img1, img2)
        !$ use omp_lib
        type(image), intent(in) :: img1,img2
        real :: ssim2, c1, c2, uxt,uyt,sxt,syt,sxyt
        integer :: i,j,k,s,t,u,npix, winsz, ldims(3),iwinsz
        real, allocatable :: rmat1(:,:,:), rmat2(:,:,:), ssim(:,:,:),&
            ux(:,:,:),uy(:,:,:),uxy(:,:,:), sx(:,:,:), sy(:,:,:), sxy(:,:,:)
        winsz = 5
        iwinsz = 2
        if( img1%is_3d() )then
            npix = (2*winsz+1)**3
        else
            npix = (2*winsz+1)**2
        endif
        ldims =img1%get_ldim()
        allocate(rmat1(1-winsz:ldims(1)+winsz-1,1-winsz:ldims(1)+winsz-1,1-winsz:ldims(1)+winsz-1 ),&
            rmat2(1-winsz:ldims(1)+winsz-1,1-winsz:ldims(1)+winsz-1,1-winsz:ldims(1)+winsz-1 ))
        ssim =  img1%get_rmat()
        rmat1(1:ldims(1), 1:ldims(2),1:ldims(3)) = img1%get_rmat()
        rmat2(1:ldims(1), 1:ldims(2),1:ldims(3)) = img2%get_rmat()
        ux  = img1%get_rmat()
        uy  = img2%get_rmat()
        uxy = img1%get_rmat()
        sx  = img1%get_rmat()
        sy  = img2%get_rmat()
        sxy = img2%get_rmat()
        ux=0.;uy=0;uxy=0.; sx=0; sy=0; sxy=0;
        !$omp parallel do collapse(3) default(shared) private(i,j,k) &
        !$omp proc_bind(close)
        do i=1,ldims(1)
            do j=1,ldims(2)
                do k=1,ldims(3)

                    !$omp parallel do collapse(3) default(shared) private(s,t,u) proc_bind(close) &
                    !$omp reduction(+:uxt,uyt,sxt,syt,sxyt)
                    do s=i-winsz,i+winsz-1
                        do t=j-winsz,j+winsz-1
                            do u=k-winsz,k+winsz-1
                               uxt = uxt + rmat1(s,t,u)
                               uyt = uyt + rmat2(s,t,u)
                               sxt = sxt + rmat1(s,t,u)**2
                               syt = syt + rmat2(s,t,u)**2
                               sxyt = sxyt + rmat1(s,t,u)* rmat2(s,t,u)
                            end do
                        end do
                    end do
                    !$omp end parallel do

                    ux(i,j,k) = uxt/npix
                    uy(i,j,k) = uyt/npix
                    sx(i,j,k) = sqrt( ((sxt/npix) - ((uxt/npix)**2.)) )
                    sy(i,j,k) = sqrt( ((syt/npix) - ((uyt/npix)**2.)) )

                    sxy(i,j,k) = sqrt( ((sxyt/npix) - ((uxt*uyt)/(npix**2)) ) )
                end do
            end do
        end do
        !$omp end parallel do
        !$omp parallel workshare
        ssim = ( 2. * ux*uy + c1  ) * (2*sxy + c2  ) /&
            ( ux**2. + uy**2. + c1  ) * (sx**2. + sy**2. + c2 )
        !$omp end parallel workshare
        ssim2 = sum(ssim)/product(ldims)

    end function ssim2


    function dssim (img1, img2)
        real :: dssim
        type(image), intent(in) :: img1,img2

        dssim = (1-ssim(img1,img2))/2

    end function dssim


end module simple_ssim
