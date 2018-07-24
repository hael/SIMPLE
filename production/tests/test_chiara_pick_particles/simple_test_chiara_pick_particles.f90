module chiara_pick_particles_mod
  include 'simple_lib.f08'
  use simple_image
  implicit none
contains
  function center_edge( img_in, lp, thresh, bin_img, discard ) result( xyz )
      class(image),   intent(inout)      :: img_in
      real,           intent(in)         :: lp, thresh
      logical, intent(out)               :: discard
      type(image)       :: tmp, bin_img, imgcc
      real              :: xyz(3)
      integer           :: ldim(3), min_sz
      real, allocatable :: rmat(:,:,:)
      min_sz = 400
      ldim = img_in%get_ldim()
      call tmp%copy(img_in)
      call tmp%bp(0., lp)
      call tmp%ifft()
      rmat = tmp%get_rmat()
      where( rmat < 0. )
           rmat = 0.
      end where
      call tmp%sobel(bin_img, thresh)
      call bin_img%real_space_filter(3, 'median') !median filtering allows me to calculate cc in an easy way
      call bin_img%calc_cc(imgcc)
      call imgcc%prepare_cc(discard, min_sz)
      xyz = imgcc%masscen()
      call tmp%kill
  end function center_edge

  function is_picked( new_coord,saved_coord, part_radius, saved ) result(yes_no)
    integer, intent(in) :: new_coord(2)       !Coordinates of a new window to extract
    integer, intent(in) :: saved_coord(:,:)   !Coordinates of extracted windows
    integer, intent(in) :: part_radius        !Approximate radius of the particle
    integer, intent(in), optional :: saved    !How many particles have already been saved
    logical :: yes_no
    integer :: iwind, ssaved, s(2)

    s = shape(saved_coord)
    if(s(2) /= 2) stop 'Dimension error'
    yes_no = .false.
    ssaved = s(1)
    if(present(saved)) ssaved = saved
    do iwind = 1, ssaved
        if(  sqrt(real((new_coord(1)-saved_coord(iwind,1))**2 + (new_coord(2)-saved_coord(iwind,2)) **2)) <= real(part_radius)) then
          yes_no = .true.
          return
        endif
    enddo
  end function is_picked
end module chiara_pick_particles_mod

program simple_test_chiara_pick_particles
  include 'simple_lib.f08'
  use simple_micops
  use simple_image
  use simple_stackops
  use simple_math
  use chiara_pick_particles_mod
  integer              :: ldim_shrunken(3), n_images, ifeat, box_shrunken, ldim(3), cnt, xind, yind, part_radius
  integer, allocatable :: coord_build(:,:), coord(:,:), coord_center_edge(:,:)
  real,    parameter   :: SHRINK = 4.
  integer, parameter   :: BOX = 370, OFFSET = BOX/SHRINK-20, BOFFSET = 1
  type(image)          :: img, img_extrac, imgwin, mic_original, mic_shrunken, bin_img, imgcc
  real                 :: smpd_shrunken, lp, xyz(3), sobel_thresh
  logical              :: outside, discard, picked
  call read_micrograph( micfname = &
  & '/home/lenovoc30/Desktop/MassCenter/MicrographsMArion/micrographs/FoilHole_1574549_Data_1584496_1584497_20180703_1915-40767_intg.mrc'&
  &, smpd = 1.)
  call shrink_micrograph(SHRINK, ldim_shrunken, smpd_shrunken)
  call set_box(BOX, box_shrunken)  !Here I can decide to add noise
  call extract_boxes2file(OFFSET, 'extracted_windows.mrc', n_images, BOFFSET, coord)
  call mic_shrunken%new(ldim_shrunken, smpd_shrunken)   !Shrunken-high pass filtered version
  call mic_shrunken%read('shrunken_hpassfiltered.mrc')
  call img_extrac%new ([box_shrunken,box_shrunken,1],1.)   !it will be the extracted window
  call imgwin%new ([box_shrunken,box_shrunken,1],1.)       !it will be the centered window
  lp = 4. !I chose the threshold for low-pass filtering
  sobel_thresh = 0.1
  part_radius = 36

  ! lp = 5.5
  ! sobel_thresh = 0.1

  picked = .false.
  cnt  = 0
  allocate(coord_center_edge(n_images,2), source = 0)

  do ifeat = 1, n_images
      call img_extrac%read('extracted_windows.mrc', ifeat)

      ! call img_extrac%sobel(bin_img, 0.15)
      ! call bin_img%write('sobel_detection.mrc',ifeat)
      ! call bin_img%real_space_filter(1, 'median')
      ! call bin_img%write('median_filtered.mrc',ifeat)

      xyz = center_edge(img_extrac, lp, sobel_thresh, bin_img, discard )
      call bin_img%write('binary_images_sobel.mrc', ifeat)
      call mic_shrunken%window_slim( [int(xyz(1))+coord(ifeat,1), int(xyz(2))+coord(ifeat,2)], &
                                                &  box_shrunken, imgwin, outside )
      if(ifeat > 1) picked = is_picked([int(xyz(1))+coord(ifeat,1), int(xyz(2))+coord(ifeat,2)],coord_center_edge, part_radius, ifeat-1 )
      coord_center_edge(ifeat,:) =     [int(xyz(1))+coord(ifeat,1), int(xyz(2))+coord(ifeat,2)]
     if(.not. picked) then  !window not already selected
          if(.not. discard) then
              if( .not. outside )then
                  cnt = cnt + 1
                  call imgwin%write('centered_particles_edge.mrc', cnt)
              endif
          endif
      endif
  end do
end program simple_test_chiara_pick_particles
