program test_foofactory
use simple_foofactory, only: foofactory
use simple_foo, only: foo
implicit none

class(foofactory), allocatable :: factory
class(foo), pointer :: mydata

allocate( foofactory :: factory )

mydata => factory%construct('foor', 5)
call mydata%print

mydata => factory%construct('fooi', 5)
call mydata%print

end program