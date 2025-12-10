program simple_test_star_export
include 'simple_lib.f08'
!$ use omp_lib
use simple_sp_project
use simple_starfile
use simple_parameters
use simple_cmdline
implicit none
#include "simple_local_flags.inc"
character(len=*), parameter :: projfile   = 'test.simple'
character(len=*), parameter :: opticsfile = 'optics.simple'
type(sp_project)            :: spproj
integer(timer_int_kind)     :: ms0
real(timer_int_kind)        :: ms_complete

if(.not. file_exists(string(projfile))) THROW_HARD(projfile // " does not exist")

ms0 = tic()
call spproj%read(string(projfile))
ms_complete = toc(ms0)
print *,'read project file in : ', ms_complete; call flush(6)

ms0 = tic()
call spproj%write_mics_star()
ms_complete = toc(ms0)
print *,'write_mics_star file in : ', ms_complete; call flush(6)

ms0 = tic()
call spproj%write_ptcl2D_star()
ms_complete = toc(ms0)
print *,'write_ptcl2D_table file in : ', ms_complete; call flush(6)

end program simple_test_star_export
