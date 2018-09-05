program simple_test_chiara_pick_edge
  include 'simple_lib.f08'
  use simple_picker_chiara
  use simple_micops
  use simple_image
  use simple_stackops
  use simple_math
  use simple_edge_detector, only: automatic_thresh_sobel
  use simple_parameters,    only: parameters
  use simple_cmdline,       only: cmdline
  type(cmdline)      :: cline
  type(parameters)   :: params
  type(image)        :: mic_shrunken, mic_bin, mic_lp, mic_copy, mic_closed
  type(image)        :: imgcc, imgwin, img_thresh
  integer            :: ldim_shrunken(3), box_shrunken,  xind, yind, min_sz, max_sz
  real               :: part_radius
  real,    parameter :: SHRINK = 4.
  real               :: smpd_shrunken, lp, sobel_thresh(1), ave, sdev, maxv, minv
  logical            :: outside, discard, picked
  ! In this new approach I don't extract any window, I just work directly on the entire
  ! micrograph using edge detection.

  ! Image processing steps: 1) Shrink and high pass filtering
  !                         2) Low pass filtering
  !                         3) Edge Detection
  !                         4) Median Filtering
  !                         5) Connected components (cc) identification
  !                         6) cc filtering

  ! As a test use
  !      fname = '/home/lenovoc30/Desktop/MassCenter/NegativeStaining/16.06.34 CCD Acquire_0000_1.mrc'
  !      smpd  = 1.
  !      part_radius = 15.
  !      part_concentration = 0.2

  if( command_argument_count() < 3 )then
      write(*,'(a)',advance='no') 'simple_test_chiara_try smpd=<sampling distance(in A)> [fname = file name] [part_radius = <radius of the particle (# pixels)]'
      stop
  endif
  call cline%parse_oldschool
  call cline%checkvar('fname', 1)
  call cline%checkvar('smpd', 2)
  call cline%checkvar('part_radius', 3)
  !Set defaults
  call params%new(cline)  !<read cline parameters
  ! 0) Reading and saving original micrograph
  call read_micrograph(micfname = params%fname, smpd = params%smpd)
  ! 1) Shrink and high pass filtering
  part_radius = params%part_radius
  call shrink_micrograph(SHRINK, ldim_shrunken, smpd_shrunken)
  call set_box(int(SHRINK*(4*part_radius+10)), box_shrunken)
  call mic_shrunken%new(ldim_shrunken, smpd_shrunken)
  call mic_shrunken%read('shrunken_hpassfiltered.mrc')
  ! 2) Low pass filtering
  lp = 35.
  call mic_copy%copy(mic_shrunken) !work on a copy not to modify the original mic
  call mic_copy%bp(0.,lp)
  call mic_copy%ifft()
  call mic_lp%copy(mic_shrunken) !work on a copy not to modify the original mic
  call mic_lp%bp(0.,20.)
  call mic_lp%ifft()
  call mic_lp%write('LowPassFiltered.mrc')

  !
  ! ! 3) Edge Detection
  ! call mic_copy%stats( ave, sdev, maxv, minv )
  ! call mic_copy%bin(ave+.7*sdev)
  ! call mic_copy%write('Bin1.mrc')
  ! call mic_copy%real_space_filter(4, 'median') !median filtering allows easy calculation of cc
  ! call mic_copy%write('Bin1Median.mrc')
  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! ! call mic_bin%morpho_opening()
  ! ! call mic_bin%morpho_closing()
  ! ! call mic_bin%write('MorphoOpenedClosed.mrc')
  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! ! 5) Connected components (cc) identification
  ! call imgcc%new(ldim_shrunken, smpd_shrunken)
  ! call mic_copy%find_connected_comps(imgcc)
  ! call imgcc%write('ConnectedComponents.mrc')
  ! ! 6) cc filtering


  call imgcc%new(ldim_shrunken, smpd_shrunken)
  call imgcc%read('/home/lenovoc30/Desktop/MassCenter/try1/ConnectedComponents.mrc')

  !min_sz =  int((15/100)*3*part_radius**2)   !15% of the size of the particle (suppose it's circular)
  min_sz = 10*int(part_radius+5)
  max_sz = 70*int(part_radius)
  call imgcc%elim_cc([min_sz,max_sz])
  call imgcc%write('ConnectedComponentsElimin.mrc')
  call extract_particles_NOmasscen(mic_lp, imgcc, int(part_radius))
end program simple_test_chiara_pick_edge
