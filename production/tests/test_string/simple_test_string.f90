program simple_test_string
use simple_string
implicit none

type(string) :: a, b, c

a = "HELLO"
b = " WORLD"

c = a // b

print *, "c:", c%to_char()
print *, "Length:", c%strlen()
    
end program simple_test_string