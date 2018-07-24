program simple_test_chiara_understanding_ppca
  include 'simple_lib.f08'
  use simple_micops
  use simple_image
  use simple_stackops
  use simple_ppca
  use simple_math
  implicit none
  integer              :: ldim_shrunken(3), n_images, D, recsz, ifeat, box_shrunken, ldim(3), n_proj, cnt, xind, yind
  real,    allocatable :: avg(:), dat(:), rmat(:,:,:)
  logical, allocatable :: l_mask(:,:,:)
  integer, allocatable :: coord_build(:,:)
  real,    parameter   :: SHRINK = 2.
  integer, parameter   :: BOX = 256, OFFSET = BOX/SHRINK, BOFFSET = 1, eigen = 8, iter = 10 !BOX/SHRINK
  type(ppca)           :: my_ppca
  type(image)          :: img, img_rev, img_msk, mic_denoised, mic_denoised_norm, mic_original
  real                 :: smpd_shrunken

  !process to build up my micrograph
  n_proj = 64
  ! generate sampling coordinates
  cnt = 0
  do xind=1,BOX*8,BOX
      do yind=1,BOX*8,BOX
          cnt = cnt + 1
      end do
  end do
  allocate(coord_build(cnt,2))
  cnt = 0
  do xind=1,BOX*8,BOX
      do yind=1,BOX*8,BOX
          cnt = cnt + 1
          coord_build(cnt,:) = [xind,yind]
      end do
  end do
  write(*,*)'cnt:', cnt
  call mic_original%new([BOX*8, BOX*8,1], 1.) !building my micrograph
  call img%new([BOX,BOX,1],1.)
  do ifeat=1,n_proj
      call img%read('outstk.mrc', ifeat)
      call mic_original%add_window(img, coord_build(ifeat,:))
  end do
  deallocate(coord_build)
!To add "carbon" (a source of high variance)
      call img%read('outstk.mrc', 17)
      rmat = mic_original%get_rmat()
      rmat(50:60,50:60,1) = 4.
      rmat(500:600,1800:2000,1) = 3.5
      rmat(964:1300,500:599,1) = 4.
      rmat(100:150,700:860,1) = 3.5
      call mic_original%set_rmat(rmat)
      deallocate(rmat)
  ! first iteration of ppca
  call mic_original%write('original_micrograph.mrc')
  call read_micrograph( micfname = 'original_micrograph.mrc', smpd = 1.0)
  call shrink_micrograph(SHRINK, ldim_shrunken, smpd_shrunken)
  call set_box(BOX, box_shrunken, 2.)  !2 is the SNR
  call extract_boxes2file(OFFSET, 'extracted_windows2.mrc', n_images, BOFFSET)
  call img_msk%new([box_shrunken,box_shrunken,1], 1.)
  img_msk = 1.
  allocate(l_mask(box_shrunken,box_shrunken,1), source = .true.)
  call make_pattern_stack('extracted_windows2.mrc', 'vecs4ppca.bin', l_mask , D, recsz, avg) !(stackops)
  call my_ppca%new(n_images, D, eigen)
  call my_ppca%master('vecs4ppca.bin', recsz ,'feat_stk.bin', iter)
  ! sample the generative model and generate back the micrograph
  call img_rev%new ([box_shrunken,box_shrunken,1],1.)       !it will be the generate sample
  call mic_denoised%new(ldim_shrunken, smpd_shrunken)       !denoised version of the shrunken high pass filtered
  call mic_denoised_norm%new(ldim_shrunken, smpd_shrunken)  !just to keep track of the overlapping, so that I can divide
  do ifeat=1,n_images
      dat = my_ppca%generate(ifeat,avg)
      call img_rev%unserialize(l_mask, dat)
      call img_rev%write('generative_samples2.mrc', ifeat)
  end do
  ! second iteration pf ppca with different level of noise
  call read_micrograph( micfname = 'original_micrograph.mrc', smpd = 1.0)
  call shrink_micrograph(SHRINK, ldim_shrunken, smpd_shrunken)
  call set_box(BOX, box_shrunken, 1.)
  call extract_boxes2file(OFFSET, 'extracted_windows1.mrc', n_images, BOFFSET)
  call make_pattern_stack('extracted_windows1.mrc', 'vecs4ppca.bin', l_mask , D, recsz, avg) !(stackops)
  call my_ppca%new(n_images, D, eigen)
  call my_ppca%master('vecs4ppca.bin', recsz ,'feat_stk.bin', iter)
  img_rev = 0.
  mic_denoised= 0.
  mic_denoised_norm= 0.
  dat = 0.
  do ifeat=1,n_images
      dat = my_ppca%generate(ifeat,avg)
      call img_rev%unserialize(l_mask, dat)
      call img_rev%write('generative_samples1.mrc', ifeat)
  end do
 ! third iteration pf ppca with higher level of noise
  call read_micrograph( micfname = 'original_micrograph.mrc', smpd = 1.0)
  call shrink_micrograph(SHRINK, ldim_shrunken, smpd_shrunken)
  call set_box(BOX, box_shrunken, 0.1)
  call extract_boxes2file(OFFSET, 'extracted_windows01.mrc', n_images, BOFFSET)
  call make_pattern_stack('extracted_windows01.mrc', 'vecs4ppca.bin', l_mask , D, recsz, avg) !(stackops)
  call my_ppca%new(n_images, D, eigen)
  call my_ppca%master('vecs4ppca.bin', recsz ,'feat_stk.bin', iter)
  img_rev = 0.
  mic_denoised= 0.
  mic_denoised_norm= 0.
  dat = 0.
  do ifeat=1,n_images
      dat = my_ppca%generate(ifeat,avg)
      call img_rev%unserialize(l_mask, dat)
      call img_rev%write('generative_samples01.mrc', ifeat)
  end do
  deallocate(avg, dat,l_mask)
end program simple_test_chiara_understanding_ppca
