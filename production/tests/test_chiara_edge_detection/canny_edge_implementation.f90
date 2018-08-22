program canny_no_thresh
use simple_image, only : image
use simple_jpg
use simple_edge_detector
implicit none
type(image)  :: img_in, img_out, img_in2, img_in1, img_in05
real :: thresh(2)
thresh = reshape([0.5,0.99],[2])
call img_in%new([128,128,1],1.)
call img_in2%new([128,128,1],1.)
call img_in1%new([128,128,1],1.)
call img_in05%new([128,128,1],1.)
call img_in%read('/home/lenovoc30/Desktop/Edge Detection/one_projection.mrc')
call img_in2%read('/home/lenovoc30/Desktop/Edge Detection/one_projection2.mrc')
call img_in1%read('/home/lenovoc30/Desktop/Edge Detection/one_projection1.mrc')
call img_in05%read('/home/lenovoc30/Desktop/Edge Detection/one_projection05.mrc')
!Let's use this image and apply three different metods for edge detection on it:
!1) Sobel, threshold;
!2) Canny, threshold;
!3) Canny, threshold: NO.
!NOISELESS IMAGE
    call sobel(img_in,img_out, thresh(1))           !1)
    call img_out%write('Sobel.mrc')
    call canny(img_in,img_out, [0.1,0.5])           !2)
    call img_out%write('Canny_thresh.mrc')
    call canny(img_in,img_out)                      !3)
    call img_out%write('Canny_NO_thresh.mrc')
!IMAGE WITH SNR = 2
    call sobel(img_in2,img_out, thresh(2))          !1)
    call img_out%write('Sobel2.mrc')
    call canny(img_in2,img_out, [0.5,0.95])         !2)
    call img_out%write('Canny_thresh2.mrc')
    call canny(img_in2,img_out)                     !3)
    call img_out%write('Canny_NO_thresh2.mrc')
!IMAGE WITH SNR = 1
    call sobel(img_in1,img_out, thresh(2))          !1)
    call img_out%write('Sobel1.mrc')
    call canny(img_in1,img_out, [0.5,0.95])         !2)
    call img_out%write('Canny_thresh1.mrc')
    call canny(img_in1,img_out)                     !3)
    call img_out%write('Canny_NO_thresh1.mrc')
!IMAGE WITH SNR = 0.5
    call sobel(img_in05,img_out, thresh(2))         !1)
    call img_out%write('Sobel05.mrc')
    call canny(img_in05,img_out, [0.5,0.95])        !2)
    call img_out%write('Canny_thresh05.mrc')
    call canny(img_in05,img_out)                    !3)
    call img_out%write('Canny_NO_thresh05.mrc')
end program canny_no_thresh
