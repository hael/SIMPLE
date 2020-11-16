program ppca_denoising
include 'simple_lib.f08'
use simple_micops
use simple_image
use simple_stackops
use simple_ppca
use simple_math
implicit none
integer              :: ldim_shrunken(3), n_images, D, recsz, ifeat, box_shrunken, ldim(3)
real,    allocatable :: avg(:), feat(:), dat(:), matrix(:,:,:)
logical, allocatable :: l_mask(:,:,:)
integer, allocatable :: coords(:,:)
! integer, parameter   :: BOX = 32, OFFSET = BOX/2-1 , BOFFSET = 3 !BOX/2-1
! real,    parameter   :: SHRINK = 2
integer, parameter   :: BOX = 256, OFFSET = BOX/4-1, BOFFSET = 20 !BOX/4-1
real,    parameter   :: SHRINK = 4.
type(ppca)           :: my_ppca
type(image)          :: img, img_rev, img_msk, mic_denoised, mic_denoised_norm
real                 :: mskrad, smpd_shrunken
logical              :: do_overlap, do_discextract

do_overlap = OFFSET /= BOX/SHRINK
!do_discextract = .true.
do_discextract = .false.

call read_micrograph( micfname = '/home/lenovoc30/Desktop/ANTERGOS/forctf/0001_forctf.mrc', smpd = 1.0)
!call read_micrograph( micfname = '/home/lenovoc30/Desktop/ANTERGOS/forctf/0002_forctf.mrc', smpd = 1.0)  !TO SEE THE EFFECTS OF CARBON PRESENCE
!call read_micrograph( micfname = '/home/lenovoc30/Desktop/PPCA/outstk.mrc', smpd = 1.0)

! SHRINK-fold downsampling (micops)
call shrink_micrograph(SHRINK, ldim_shrunken, smpd_shrunken)
! setting box + high-pass filtering (micops)
call set_box(BOX, box_shrunken)
! extracting with an offset of pixels
call extract_boxes2file(OFFSET, 'extracted_windows.mrc', n_images, BOFFSET)

! generate vectors for compression
mskrad = (real(BOX)/2.)/SHRINK
if(.not. do_overlap) mskrad = (real(BOX)/2.)/SHRINK  !Biggest radius possible
call img_msk%disc([box_shrunken,box_shrunken,1], 1., mskrad, l_mask)
if(.not. do_discextract ) then
  l_mask = .true.                 !extract every pixel in the window
 matrix = logical2bin(l_mask)
 call img_msk%set_rmat(matrix)    !update img_msk
 deallocate(matrix)
endif

call make_pattern_stack('extracted_windows.mrc', 'vecs4ppca.bin', l_mask , D, recsz, avg) !(stackops)
!call del_file('extracted_windows.mrc')

! n_images = N(# images in stack), D(dimension of the features vectors), Q(# eigenvectors)
call my_ppca%new(n_images, D, 5)
! ppca with maxits=10
call my_ppca%master('vecs4ppca.bin', recsz ,'feat_stk.bin', 10)

! sample the generative model and generate back the micrograph
call img_rev%new ([box_shrunken,box_shrunken,1],1.)       !it will be the generate sample
call mic_denoised%new(ldim_shrunken, smpd_shrunken)       !denoised version of the shrunken high pass filtered
call mic_denoised_norm%new(ldim_shrunken, smpd_shrunken)  !just to keep track of the overlapping, so that I can divide
coords = generate_sampling_coordinates(OFFSET)            !(micops)
do ifeat=1,n_images
    dat = my_ppca%generate(ifeat,avg)        !generate back
    call img_rev%read('extracted_windows.mrc', ifeat)
    call img_rev%unserialize(l_mask, dat)
    !call img_rev%write('generative_samples.mrc', ifeat)
    call mic_denoised%add_window(img_rev, coords(ifeat,:))       !(image) sum everything, then I will divide
    call mic_denoised_norm%add_window(img_msk, coords(ifeat,:))  !to build denominator
end do
call mic_denoised%div(mic_denoised_norm)
call mic_denoised%write('denoised_micrograph.mrc', 1)
call mic_denoised_norm%write('denominator.mrc', 1)
end program ppca_denoising
