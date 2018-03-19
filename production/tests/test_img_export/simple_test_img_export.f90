 !! Modified by Michael Eager, Feb 2018
program simple_test_img_export
include 'simple_lib.f08'
use simple_jpg, only : test_jpg_export
use simple_test_export_jpg, only : test_jpg_image
use simple_test_export_libgd, only : test_png_io, test_jpeg_io, create_raw_png_tmp
use simple_test_export_pnm, only : test_pnm_io
use simple_cmdline, only: cmdline
implicit none
#include "simple_local_flags.inc"
type(cmdline)         :: cline
logical               :: be_verbose=.false.
character(len=8)      :: datestr
character(len=STDLEN) :: folder, command, olddir
if( command_argument_count() > 0 )then
   write(*,'(a)',advance='no') 'simple_test_img_export vol1=<volume.mrc> msk=<mask radius(in pixels)>'
   write(*,'(a)') ' smpd=<sampling distance(in A)> [nthr=<number of threads{1}>] [verbose=<yes|no{no}>]'
   call cline%parse_oldschool
   call cline%check
   if( cline%defined('vol1') )then
      if( trim(cline%get_carg('vol1')) .eq. '' )then
         print *,' No input volume for testing'
      endif
   endif
   be_verbose = .false.
   if( cline%defined('verbose') )then
      if( trim(cline%get_carg('verbose')) .eq. 'yes' )then
         be_verbose = .true.
      endif
   endif
endif
call seed_rnd
call date_and_time(date=datestr)
folder = trim('SIMPLE_TEST_IMGEXPORT_'//datestr)
call simple_mkdir(folder)
call simple_chdir(folder, olddir)

call test_pnm_io
call test_jpg_image(.true.)
call test_jpg_export

#ifdef _LIBGD
 !   call test_jpeg_io
 !   call test_png_io
#endif
call simple_chdir(olddir)

end program simple_test_img_export
