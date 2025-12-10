program simple_test_mrc2jpeg
include 'simple_lib.f08' 
use simple_image
use simple_jpg
use simple_cmdline,    only: cmdline
use simple_parameters, only: parameters
implicit none
type(string), allocatable :: micname(:)
type(image)               :: microg 
type(string)              :: outputfile, fbody
integer                   :: i, j, nfiles, ldim(3), ifoo, ldim_refs(3)
type(cmdline)             :: cline
type(parameters)          :: p
if( command_argument_count() < 2 )then
    write(logfhandle,'(a)') 'Usage: simple_test_mrc2jpg filetab=filetab.txt smpd=0.3'
    stop
endif
call cline%parse_oldschool
call cline%checkvar('filetab', 1)
call cline%checkvar('smpd',    2)
call cline%check
call p%new(cline)
call read_filetable(p%filetab, micname)
nfiles=size(micname)
do i=1,nfiles
    fbody = get_fbody(basename(micname(i)),'mrc')
    call find_ldim_nptcls(micname(i), ldim, ifoo)
    ldim_refs = [ldim(1), ldim(2), 1]
    if(ldim(3)==1)then
        call microg%new([ldim(1), ldim(2), 1], p%smpd)
        call microg%read(micname(i))
        outputfile = fbody//'.jpeg'
        write(logfhandle,'(a)') '>>> WRITING '//outputfile%to_char()
        call microg%write_jpg(outputfile)
        call microg%kill()
    else
        call microg%new([ldim(1), ldim(2), 1], p%smpd)
        do j=1,ldim(3)
            call microg%read(micname(i),j)
            outputfile = fbody//int2str_pad(j,3)//'.jpeg'
            write(logfhandle,'(a)') '>>> WRITING '//outputfile%to_char()
            call microg%write_jpg(outputfile)
        enddo
        call microg%kill()
    endif
enddo
end program simple_test_mrc2jpeg
