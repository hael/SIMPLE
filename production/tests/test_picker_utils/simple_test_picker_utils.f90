program simple_test_picker_utils
!$ use omp_lib
!$ use omp_lib_kinds
include 'simple_lib.f08'
use simple_picker_utils, only: picker_utils
use simple_image,        only: image
implicit none

character(len=*), parameter :: micname = '/Users/elmlundho/Processing/relion_tut/20170629_00021_frameImage_intg.mrc'
! character(len=*), parameter :: micname = '/home/elmlundho/cache/pick_bench/optimal_movie_average.mrc'
real,             parameter :: SMPD    = 0.885, MOLDIAM = 180., MAXDIAM = MOLDIAM + 20.
! real,             parameter :: SMPD    = 2.0, MOLDIAM = 200., MAXDIAM = MOLDIAM + 20.
type(image)        :: micimg
type(picker_utils) :: putils
integer            :: ldim(3), ifoo, nthr

!$ nthr = omp_get_max_threads()
!$ call omp_set_num_threads(nthr)
nthr_glob = nthr

call find_ldim_nptcls(micname, ldim, ifoo)
call micimg%new(ldim, SMPD)
call micimg%read(micname)
call putils%set_mics(micimg, SMPD)
call putils%gauconv_mic_shrink1(MAXDIAM)
call putils%set_positions(MAXDIAM)
call putils%analyze_boximgs1(MAXDIAM)

end program simple_test_picker_utils
