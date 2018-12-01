program simple_test_chiara_try1
  include 'simple_lib.f08'
  use simple_picker_chiara, only : extract_particles_NOmasscen
  use simple_micops
  use simple_image,         only : image
  use simple_stackops
  use simple_math
  use simple_segmentation, only: sobel, automatic_thresh_sobel, canny
  use simple_parameters,    only: parameters
  use simple_cmdline,       only: cmdline
  type(cmdline)      :: cline
  type(parameters)   :: params
  type(image)        :: mic_shrunken, mic_lp, mic_copy
  type(image)        :: imgcc, imgwin
  integer            :: ldim_shrunken(3), box_shrunken, min_sz, max_sz
  real               :: smpd_shrunken, part_radius
  real,    parameter :: SHRINK = 4.
  real               :: lp, thresh(1), ave, sdev, maxv, minv


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
  !      smpd        = 1.
  !      part_radius = 15.
  !      detector    = bin

  if( command_argument_count() < 4 )then
      write(logfhandle,'(a)',advance='no') 'simple_test_chiara_try smpd=<sampling distance(in A)> [fname = file name] [part_radius = <radius of the particle (# pixels)] [detector= <binarisation method> (sobel|canny|bin)]'
      stop
  endif
  call cline%parse_oldschool
  call cline%checkvar('fname', 1)
  call cline%checkvar('smpd', 2)
  call cline%checkvar('part_radius', 3)
  call cline%checkvar('detector', 4)
  !Set defaults
  call params%new(cline)  !<read cline parameters
  ! 0) Reading and saving original micrograph
  call read_micrograph(micfname = params%fname, smpd = params%smpd)

  ! 1) Shrink and high pass filtering
  part_radius = params%part_radius
  call shrink_micrograph(SHRINK, ldim_shrunken, smpd_shrunken)
  call set_box(int(SHRINK*( 4*part_radius+10 )), box_shrunken) !here it is high pass filt
  call mic_shrunken%new(ldim_shrunken, smpd_shrunken)
  call mic_shrunken%read('shrunken_hpassfiltered.mrc')

  ! 1.1) Standardization by Histogram stretching
  call mic_shrunken%hist_stretching(mic_copy)

  ! 1.2) Non local mean filtering
  call mic_copy%NLmean()
  call mic_copy%write('NLmean_filtered.mrc')

  ! 2) Low pass filtering
  lp = 35. !40
  call mic_copy%copy(mic_shrunken) !work on a copy not to modify the original mic
  call mic_copy%bp(0.,lp)
  call mic_copy%ifft()
  call mic_copy%write('bp_filtered.mrc')


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!just to draw on that version, which is the most visible for human eye!!!!!!!!!!!!!
  call mic_lp%copy(mic_shrunken) !work on a copy not to modify the original mic
  call mic_lp%bp(0.,20.)
  call mic_lp%ifft()
  call mic_lp%write('LowPassFiltered.mrc')
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! 3) Edge Detection
  call mic_copy%stats( ave, sdev, maxv, minv )
  if(params%detector .eq. 'sobel') then
    thresh(1) = ave+.7*sdev
    call sobel(mic_copy,thresh)
  else if (params%detector .eq. 'bin') then
    call mic_copy%bin(ave+.7*sdev) !(ave+.8*sdev)
    print *, 'threshold = ', ave+.7*sdev
  else if (params%detector .eq. 'canny') then
    call canny(mic_copy)
  endif
  call mic_copy%write('Bin1.mrc')
  call mic_copy%real_space_filter(3,'median') !median filtering allows easy calculation of cc. winsz=5??
  call mic_copy%write('Bin1Median.mrc')
  ! call automatic_thresh_sobel(mic_copy, 0.2, thresh)
  ! print *, 'Selected threshold = ', thresh
  ! 5) Connected components (cc) identification
  call imgcc%new(ldim_shrunken, smpd_shrunken)

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !TO MAKE IT FASTER
  !call mic_copy%read('Bin1Filtered.mrc')
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  call mic_copy%find_connected_comps(imgcc)
  call imgcc%write('ConnectedComponents.mrc')

  ! 6) cc filtering
  min_sz = 20*int(part_radius)  !CHOSE A PERCENTAGE OF THE of the size of the particle
  max_sz = 150*int(part_radius)
  call imgcc%elim_cc([min_sz,max_sz])
  call imgcc%write('ConnectedComponentsElimin.mrc')
  call extract_particles_NOmasscen(mic_lp, imgcc, int(part_radius))
  open(unit = 17, file = "PickerInfo.txt")
  write(unit = 17, fmt = '(a)') '>>>>>>>>>>>>>>>>>>>>PARTICLE PICKING>>>>>>>>>>>>>>>>>>'
  write(unit = 17, fmt = '(a)') ''
  write(unit = 17, fmt = "(a,f0.0)")  'Mic Shrunken by factor ', SHRINK
  write(unit = 17, fmt = "(a,i0,tr1,i0,tr1,i0)") 'Dim after shrink ', ldim_shrunken
  write(unit = 17, fmt = "(a,f0.0)")  'Smpd after shrink ', smpd_shrunken
  write(unit = 17, fmt = "(a,i0)")  'Hp box ', int(SHRINK*( 4*part_radius+10 ))
  write(unit = 17, fmt = "(a,f0.0)")  'Lp  ', lp
  write(unit = 17, fmt = "(a,i0,tr1, i0)")  'Connected Components size filtering  ', min_sz, max_sz
  write(unit = 17, fmt = '(a)') ''
  write(unit = 17, fmt = "(a)")  'SELECTED PARAMETERS '
  write(unit = 17, fmt = '(a)') ''
  write(unit = 17, fmt = "(a,tr1,f0.0)")  'smpd ', params%smpd
  write(unit = 17, fmt = "(a,tr1,f0.0)")  'part_radius  ', params%part_radius
  write(unit = 17, fmt = "(a,a)")  'detector  ', params%detector
  close(17, status = "keep")
end program simple_test_chiara_try1
