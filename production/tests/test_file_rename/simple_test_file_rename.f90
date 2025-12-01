program simple_test_file_rename
include 'simple_lib.f08'
implicit none
character(len=*), parameter :: FNAME = '/Users/elmlundho/src/SIMPLE/build/CTestTestfile.cmake'
character(len=*), parameter :: HELLO_SUFFIX = '_hello'
type(string) :: new_name

new_name = append2basename(string(FNAME), string(HELLO_SUFFIX))
print *, new_name%to_char()

end program simple_test_file_rename
