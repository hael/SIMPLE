program simple_test_cavgs_selection
include 'simple_lib.f08'
use simple_image,    only: image
implicit none
#include "simple_local_flags.inc"
type(image) :: cavg
character(len=STDLEN), parameter :: cavgs_file='cavgs_iter015.mrc'
character(len=STDLEN), parameter :: filetable='filetab.txt'
!integer, allocatable :: cavgs_number(:)
real,    parameter   :: smpd=0.732
integer :: ldim(3), ncavgs, icavgs, isel
!call read_filetable(filetable, cavgs_number)
call find_ldim_nptcls(cavgs_file, ldim, ncavgs)
print *, 'Number of cavgs', ncavgs, ldim
read(5,*) isel
do icavgs = 1, ncavgs
    print *, icavgs, "cavgs"
    call cavg%new([ldim(1),ldim(2),1], smpd)
    call cavg%read(cavgs_file, icavgs)
    if( icavgs .eq. isel )then
        call cavg%write('pickrefs.mrc')
        exit
    endif
    call cavg%kill()
enddo
end program simple_test_cavgs_selection
