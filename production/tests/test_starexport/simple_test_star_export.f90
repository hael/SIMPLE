program simple_test_star_export
include 'simple_lib.f08'
!$ use omp_lib
use simple_sp_project,         only: sp_project
use simple_starproject_stream, only: starproject_stream

implicit none
#include "simple_local_flags.inc"

character(len=*), parameter :: projfile   = 'test.simple'
character(len=*), parameter :: opticsfile = 'optics.simple'

type(sp_project)            :: spproj, opticsproj
type(starproject_stream)    :: starproj_stream
integer(timer_int_kind)     :: ms0
real(timer_int_kind)        :: ms_complete

if(.not. file_exists(projfile)) THROW_HARD(projfile // " does not exist")

ms0 = tic()
call spproj%read(projfile)
ms_complete = toc(ms0)
print *,'read project file in : ', ms_complete; call flush(6)

if(file_exists(opticsfile)) then
    ms0 = tic()
    call opticsproj%read(opticsfile)
    ms_complete = toc(ms0)
    print *,'read optics file in : ', ms_complete; call flush(6)
endif

print *
print *, "######### MICROGRAPHS #########"

print *,'number micrographs   : ', spproj%os_mic%get_noris()   ; call flush(6)

ms0 = tic()
call starproj_stream%stream_export_micrographs( spproj, '.', optics_set=.false., filename="test_mics_nooptics.star")
ms_complete = toc(ms0)
print *,'micrographs exported without optics in : ', ms_complete; call flush(6)

if(file_exists(opticsfile)) then

    ms0 = tic()
    call starproj_stream%copy_optics(spproj, opticsproj)
    ms_complete = toc(ms0)
    print *,'micrographs optics copied in : ', ms_complete; call flush(6)

    ms0 = tic()
    call starproj_stream%stream_export_micrographs( spproj, '.', optics_set=.true., filename="test_mics_optics.star")
    ms_complete = toc(ms0)
    print *,'micrographs exported with optics in : ', ms_complete; call flush(6)

endif

print *
print *, "######### PARTICLES #########"

print *,'number particles     : ', spproj%os_ptcl2D%get_noris(); call flush(6)

ms0 = tic()
call starproj_stream%stream_export_particles_2D( spproj, '.', optics_set=.false., filename="test_ptcls_nooptics.star", verbose=.true.)
ms_complete = toc(ms0)
print *,'particles exported without optics in : ', ms_complete; call flush(6)

if(file_exists(opticsfile)) then

    ms0 = tic()
    call starproj_stream%copy_optics(spproj, opticsproj)
    ms_complete = toc(ms0)
    print *,'particles optics copied in : ', ms_complete; call flush(6)

    ms0 = tic()
    call starproj_stream%stream_export_particles_2D( spproj, '.', optics_set=.true., filename="test_ptcls_optics.star", verbose=.true.)
    ms_complete = toc(ms0)
    print *,'particles exported with optics in : ', ms_complete; call flush(6)

endif

end program simple_test_star_export
