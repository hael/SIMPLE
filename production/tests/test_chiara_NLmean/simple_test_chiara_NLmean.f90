program simple_test_chiara_NLmean
  include 'simple_lib.f08'
  use simple_image, only : image
  implicit none
  type(image)        :: img, img_small
  real ::  rmat(100,100,1)
  real, allocatable :: rmat_padded(:,:,:)
  !STILL TO WORK ON!
  !On too big images, why don't you split it into 4 images and analyse
  !seperatly each of them?
  call img%new([1024,1024,1],1.)
  call img%read('/home/lenovoc30/Desktop/PickingResults/SomeExamples/NegStainingWorking/try23Nov/bp_filtered.mrc')
  call img%NLmean()
  call img%write('NLmeanFiltered.mrc')
end program simple_test_chiara_NLmean
!call img_in%new([960,928,1],1.)
!call img_in%read('/home/lenovoc30/Desktop/MassCenter/real_case/shrunken_hpassfiltered.mrc') !TO TRY
!call img_in%read('/home/lenovoc30/Desktop/MassCenter/NegativeStaining/shrunken_hpassfiltered.mrc')
