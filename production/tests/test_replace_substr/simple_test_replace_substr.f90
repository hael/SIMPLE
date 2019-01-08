program simple_test_replace_substr
use simple_strings
implicit none
character(len=:), allocatable :: str1, str2
integer, allocatable :: inds(:)
allocate(str1, source='{hans{elmlund{hans_elmlund')

call replace_substr(str1, '{', '\{', str2)

print *, str2
end program simple_test_replace_substr
