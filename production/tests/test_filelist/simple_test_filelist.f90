program simple_test_filelist
!$ use omp_lib
!$ use omp_lib_kinds
include 'simple_lib.f08'
implicit none
#include "simple_local_flags.inc"
type(string),     allocatable :: thumb_files(:)
character(len=*), parameter   :: DIR_THUMBS = 'thumbnails/'

call simple_list_files('*jpg', thumb_files)
call simple_mkdir(string(DIR_THUMBS))
call move_files2dir(string(DIR_THUMBS), thumb_files)
call thumb_files%kill

end program simple_test_filelist
