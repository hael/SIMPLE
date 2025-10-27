program simple_test_file_rename
include 'simple_lib.f08'
implicit none

character(len=*), parameter   :: FNAME = '/Users/elmlundho/src/SIMPLE/build/CTestTestfile.cmake'
character(len=*), parameter   :: HELLO_SUFFIX = '_hello'
character(len=:), allocatable :: bname

print *, append2basename(FNAME, HELLO_SUFFIX)

end program simple_test_file_rename