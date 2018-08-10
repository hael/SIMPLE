program simple_test_chiara_NLmean
  include 'simple_lib.f08'
  use simple_image, only : image
  implicit none
  type(image)        :: img_in
  !STILL TO WORK ON
  ! call img_in%new([1024,1024,1],1.)
  ! call img_in%read('/home/lenovoc30/Desktop/ANTERGOS/forctf/0001_forctf.mrc')
  ! call img_in%write('OriginalCryo.mrc')
  ! call img_in%NLmean()
  ! call img_in%write('DenoisedCryo.mrc')
end program simple_test_chiara_NLmean
!call img_in%new([960,928,1],1.)
!call img_in%read('/home/lenovoc30/Desktop/MassCenter/real_case/shrunken_hpassfiltered.mrc') !not working
!call img_in%read('/home/lenovoc30/Desktop/MassCenter/NegativeStaining/shrunken_hpassfiltered.mrc')
