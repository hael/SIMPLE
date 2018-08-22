program simple_test_chiara_pick_edge
include 'simple_lib.f08'
use simple_picker_chiara
use simple_micops
use simple_image
use simple_stackops
use simple_math

use simple_edge_detector, only : automatic_thresh_sobel
use simple_parameters, only: parameters
use simple_cmdline,    only: cmdline
type(cmdline)      :: cline
type(parameters)   :: params
type(image)        :: mic_shrunken,  mic_bin, mic_back_no_aggreg, mic_copy
type(image)        :: imgcc, imgwin
integer            :: ldim_shrunken(3), box_shrunken,  xind, yind, min_sz, max_sz
real               :: part_radius
real,    parameter :: SHRINK = 4.
integer, parameter :: BOX = 230, OFFSET = BOX/SHRINK-20, BOFFSET = 1
real               :: smpd_shrunken, lp, sobel_thresh(1)
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
if( command_argument_count() < 3 )then
    write(*,'(a)',advance='no') 'simple_test_chiara_try smpd=<sampling distance(in A)> [fname = file name] [part_radius = <radius of the particle (# pixels)]'
    stop
endif
call cline%parse_oldschool
call cline%checkvar('fname', 1)
call cline%checkvar('smpd', 2)
call cline%checkvar('part_radius', 3)
!Set defaults
if( .not. cline%defined('part_concentration') ) call cline%set('part_concentration', 0.3)
call params%new(cline)  !<read cline parameters
! 0) Reading and saving original micrograph
call read_micrograph(micfname = params%fname, smpd = params%smpd)!, .true.)
! 1) Shrink and high pass filtering
call shrink_micrograph(SHRINK, ldim_shrunken, smpd_shrunken)
call set_box(BOX, box_shrunken)
call mic_shrunken%new(ldim_shrunken, smpd_shrunken)
call mic_shrunken%read('shrunken_hpassfiltered.mrc')
! call mic_original%new(ldim_shrunken*SHRINK, smpd_shrunken/SHRINK)
! call mic_original%read('Original_micrograph.mrc')
! 2) Low pass filtering
lp = 35.
call mic_copy%copy(mic_shrunken) !work on a copy not to modify the original mic
call mic_copy%bp(0.,lp)
call mic_copy%ifft()
call mic_copy%write('LowPassFiltered.mrc')
! 3) Edge Detection
call automatic_thresh_sobel(mic_copy, mic_bin, params%part_concentration, sobel_thresh(1))  !Automatic threshold selection
call mic_bin%write('AutomaticThresh.mrc')
! 4) Median Filtering
call mic_bin%real_space_filter(6, 'median') !median filtering allows me to calculate cc in an easy way
call mic_bin%write('MedianFiltered.mrc')
! 5) Connected components (cc) identification
call imgcc%new(ldim_shrunken, smpd_shrunken)
call mic_bin%find_connected_comps(imgcc)
! 6) cc filtering
part_radius = params%part_radius
min_sz =  7*int(part_radius)
max_sz = 56*int(part_radius)
call imgcc%elim_cc([min_sz,max_sz])
call imgcc%write('ElimSmallBigCC.mrc')
! Particle extraction
call extract_particles(mic_shrunken, imgcc, int(part_radius))
end program simple_test_chiara_pick_edge
