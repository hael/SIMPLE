program simple_test_string
use simple_string
implicit none

type(string) :: a, b, c, d, e

a = "HELLO"
b = " WORLD"
d = ""

c = a // b

print *, "c:|", c%to_char(), "|"
print *, "Length:", c%strlen(), " vs. ",len_trim(c%to_char())
print *, "d:|", d%to_char(), "|"
print *, "Empty length:", d%strlen()
print *, "e:|", e%to_char(), "|"
print *, "Unallocated length:", d%strlen()
end program simple_test_string