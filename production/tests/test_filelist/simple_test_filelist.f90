program simple_test_filelist
!$ use omp_lib
!$ use omp_lib_kinds
include 'simple_lib.f08'
implicit none
#include "simple_local_flags.inc"
character(len=LONGSTRLEN), allocatable :: thumb_files(:)
character(len=*), parameter :: DIR_THUMBS = 'thumbnails/'
integer :: i

call simple_list_files('*jpg', thumb_files)
call simple_mkdir(DIR_THUMBS)
call move_files2dir(DIR_THUMBS, thumb_files)


! do i = 1, size(file_list)
!     print *, trim(file_list(i))
! end do

deallocate(thumb_files)

end program simple_test_filelist
