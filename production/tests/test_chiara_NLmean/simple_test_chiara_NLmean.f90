program simple_test_chiara_NLmean
  include 'simple_lib.f08'
  use simple_image, only : image
  implicit none
  type(image)        :: img
  !STILL TO WORK ON!
  !On too big images, why don't you split it into 4 images and analyse
  !seperatly each of them?
  call img%new([960,928,1],1.)
  call img%read('shrunken_hpassfiltered.mrc')
  call img%NLmean()
  call img%write('NLmeanFilteredCryo.mrc')
end program simple_test_chiara_NLmean
