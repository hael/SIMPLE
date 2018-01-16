program simple_test_binary
use simple_fileio
implicit none
character(len=32) :: str
integer(kind=1)   :: byte_array(32)
integer :: funit
str = "hello binary world!"
call fopen(funit,'hworld.bin', access='STREAM', action='READWRITE', status='UNKNOWN', form='UNFORMATTED')
write(funit,pos=1) str
call fclose(funit)

print *, 'sizeof(byte_array): ', sizeof(byte_array)
print *, 'sizeof(str)       : ', sizeof(str)

call fopen(funit,'hworld.bin', access='STREAM', action='READWRITE', status='UNKNOWN', form='UNFORMATTED')
read(funit,pos=1) byte_array
call fclose(funit)

call fopen(funit,'hworld_as_bytearray.bin', access='STREAM', action='READWRITE', status='UNKNOWN', form='UNFORMATTED')
write(funit,pos=1) byte_array
call fclose(funit)

end program simple_test_binary
