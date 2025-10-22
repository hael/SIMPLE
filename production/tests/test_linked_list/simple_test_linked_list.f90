program simple_test_linked_list
use simple_linked_list
implicit none
type(linked_list) :: lst
type(list_iterator) :: it
class(*), allocatable :: any
integer :: i
character(:), allocatable :: s

call lst%push_back(42)              ! integer
call lst%push_back(3.14159d0)       ! real(8)
s = 'hello'; call lst%push_front(s) ! character(:)

print *, 'size =', lst%size()       ! -> 3

! random access
call lst%at(2, any)
select type(any)
type is (integer)
   print *, 'at(2) as integer =', any
class default
   print *, 'at(2) is not integer'
end select

! iterate
it = lst%begin()
do while (it%has_next())
   call it%next(any)
   select type(any)
   type is (integer)
      print *, 'int:', any
   type is (real(kind(1.0d0)))
      print *, 'real:', any
   type is (character(*))
      print *, 'char:', any
   class default
      print *, 'other type'
   end select
end do

! pop
call lst%pop_front(any)
print *, 'after pop_front, size =', lst%size()

call lst%kill
end program simple_test_linked_list
