program simple_test_serialize
  include 'simple_lib.f08'
  use simple_micops
  use simple_image
  use simple_stackops
  use simple_ppca
  use simple_math
implicit none
type(image)          :: img, img_msk, img_rev
real, allocatable    :: pcavec(:)
logical, allocatable :: l_mask(:,:,:)
integer              :: ldim_shrunken(3), n_images, D, recsz, ifeat, box_shrunken, ldim(3)
real,    allocatable :: avg(:), feat(:), dat(:), matrix(:,:,:)
integer, allocatable :: coords(:,:)
integer, parameter   :: box = 256
type(ppca)           :: my_ppca
type(image)          :: mic_denoised, mic_denoised_norm
real                 :: mskrad, smpd_shrunken
call img%new([box,box,1], 1.0)
call img_rev%new([box,box,1], 1.0)
call img%square(60)
call img_msk%disc([box,box,1], 1.0, 80.)
call img_msk%vis
l_mask = img_msk%bin2logical()
call img%vis
pcavec = img%serialize(l_mask)
call img_rev%unserialize(l_mask, pcavec)
call img_rev%vis
end program simple_test_serialize
