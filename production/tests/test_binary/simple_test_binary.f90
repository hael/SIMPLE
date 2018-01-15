program simple_test_binary
use simple_fileio
implicit none
character(len=32) :: str
integer :: funit
str = "hello binary world!"
call fopen(funit,'hworld.bin', access='STREAM', action='READWRITE', status='UNKNOWN', form='UNFORMATTED')
write(funit,pos=1) str
call fclose(funit)
end program simple_test_binary
