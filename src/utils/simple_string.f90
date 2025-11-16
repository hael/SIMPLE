module simple_string
use simple_defs
use simple_error, only: simple_exception
use, intrinsic :: iso_c_binding
implicit none
private

public :: string

type :: string
    private
    character(len=:), allocatable :: buffer
  contains
    procedure, private :: assign
    generic            :: assignment(=) => assign
    procedure, private :: append
    generic            :: operator(//)  => append
    procedure          :: strlen
    procedure          :: to_char
end type string

contains

    ! Assignment: str = "literal" or str = other_str
    subroutine assign(self, rhs)
        class(string), intent(inout) :: self
        class(*),      intent(in)    :: rhs
        select type (rhs)
        type is (string)
            call move_alloc(rhs%buffer, self%buffer)
        type is (character(*))
            self%buffer = rhs
        class default
            call simple_exception('Unsupported assignment', __FILENAME__ , __LINE__)
        end select
    end subroutine assign

    ! Append: str = str // "foo"  or str = str // other_str
    function append(lhs, rhs) result(res)
        class(string), intent(in) :: lhs
        class(*),      intent(in) :: rhs
        type(string)              :: res
        character(len=:), allocatable :: tmp
        select type (rhs)
            type is (string)
                tmp = lhs%buffer // rhs%buffer
            type is (character(*))
                tmp = lhs%buffer // rhs
            class default
                call simple_exception('Unsupported append', __FILENAME__ , __LINE__)
        end select
        if( allocated(tmp) ) call move_alloc(tmp, res%buffer)
    end function append

    ! Length
    pure integer function strlen(self)
        class(string), intent(in) :: self
        strlen = len(self%buffer)
    end function strlen

    ! Convert to plain CHARACTER for legacy code
    pure function to_char(self) result(c)
        class(string),       intent(in) :: self
        character(len=len(self%buffer)) :: c
        if( len(c) > 0) c = self%buffer
    end function to_char

end module simple_string
