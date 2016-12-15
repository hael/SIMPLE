module num

  type, abstract :: abstr_numtype
  contains
      procedure(generic_assign), private, deferred :: assign
      generic :: assignment(=) => assign
      procedure(generic_add), private, deferred :: add
      generic :: operator(+) => add
  end type

  abstract interface
  
      subroutine generic_assign(a,b)
          import :: abstr_numtype
          class(abstr_numtype), intent(out) :: a
          class(abstr_numtype), intent(in) :: b
      end subroutine
      
      function generic_add(a,b) result (r)
          import :: abstr_numtype
          class(abstr_numtype), intent(in) :: a,b
          class(abstr_numtype), allocatable :: r
      end function
      
  end interface

  type, extends(abstr_numtype) :: my_integer
      integer, public :: value
  contains
     procedure :: assign => assign_my_integer
     procedure :: add => add_my_integer 
  end type

  contains
  
    subroutine assign_my_integer(a,b)
        class(my_integer),    intent(out) :: a
        class(abstr_numtype), intent(in)  :: b
        select type (b)
            type is (my_integer)
                a%value = b%value
        end select
    end subroutine

    function add_my_integer(a,b) result(r)
        class(my_integer),    intent(in)  :: a
        class(abstr_numtype), intent(in)  :: b
        class(abstr_numtype), allocatable :: r
        select type (b)
            type is (my_integer)
                allocate(my_integer :: r)
                select type (r)
                  type is (my_integer)
                    r%value = a%value+b%value
                end select
        end select
    end function

end module

program main
  use num
  class(my_integer), allocatable :: a, b, c
  allocate(my_integer :: a)
  allocate(my_integer :: b)
  allocate(my_integer :: c)
  a=my_integer(1)
  b=my_integer(2)
  c = a+b
  write (*,*) c%value
end program