! integer_set.f90
!
!     Example belonging to "Modern Fortran in Practice" by Arjen Markus
!
!     This work is licensed under the Creative Commons Attribution 3.0 Unported License.
!     To view a copy of this license, visit http://creativecommons.org/licenses/by/3.0/
!     or send a letter to:
!     Creative Commons, 444 Castro Street, Suite 900, Mountain View, California, 94041, USA.
!
!     Simple program to determine the performance differences between
!     using array-valued functions and ordinary do-loops
!
!     Deal with sets of integers, using a "declarative" or "functional" programming style
!
module integer_sets
    implicit none

    type INTEGER_OPERAND
        integer                        :: operation
        integer                        :: value
        type(INTEGER_OPERAND), pointer :: first
        type(INTEGER_OPERAND), pointer :: second
    end type INTEGER_OPERAND

    type INTEGER_RELATION
        integer                        :: relation
        logical                        :: value
        type(INTEGER_OPERAND), pointer :: first
        type(INTEGER_OPERAND), pointer :: second
    end type INTEGER_RELATION

    type INTEGER_PAIR_SET
        type(INTEGER_RELATION)         :: relation
        integer, dimension(:), pointer :: cached
        type(INTEGER_OPERAND), pointer :: first
        type(INTEGER_OPERAND), pointer :: second
    end type INTEGER_PAIR_SET

    type INTEGER_RANGE
        type(INTEGER_OPERAND), pointer :: value
        integer                        :: min_value
        integer                        :: max_value
    end type INTEGER_RANGE

    interface assignment(=)
        module procedure integer_set_value
    end interface

    interface operator(+)
        module procedure integer_add
        module procedure integer_add_v
        module procedure integer_add_vv
    end interface

    interface operator(*)
        module procedure integer_multiply
        module procedure integer_multiply_v
        module procedure integer_multiply_vv
    end interface

    interface operator(**)
        module procedure integer_exponentiate_v
    end interface

    interface operator(.eq.)
        module procedure integer_equal
        module procedure integer_equal_v
    end interface

contains

function create_pair_set( relation, xin, yin ) result(set)
    type(INTEGER_RELATION), intent(in)         :: relation
    type(INTEGER_RANGE),    intent(in)         :: xin
    type(INTEGER_RANGE),    intent(in)         :: yin
    type(INTEGER_PAIR_SET)                     :: set

    integer, dimension(:), pointer             :: values
    integer, dimension(:), pointer             :: new_values
    integer                                    :: i
    integer                                    :: j
    integer                                    :: k
    type(INTEGER_RANGE)                        :: xrange
    type(INTEGER_RANGE)                        :: yrange

    xrange = xin
    yrange = yin

    set%first  => xrange%value
    set%second => yrange%value

    allocate( values(2*(1+xrange%max_value-xrange%min_value)) )

    k = 0
    do j = yrange%min_value, yrange%max_value
        yrange%value = j

        do i = xrange%min_value, xrange%max_value
            xrange%value = i

            if ( integer_relation_eval(relation) ) then
                k = k + 1
                if ( 2*k > size(values) ) then
                    allocate( new_values(2*size(values)) )
                    new_values(1:size(values)) = values
                    deallocate( values )
                    values => new_values
                endif

                values(2*k-1) = i
                values(2*k)   = j
            endif
        enddo
    enddo

    !
    ! Store it all in a suitable array
    !
    set%relation = relation
    allocate( new_values(2*k) )
    new_values = values(1:2*k)
    deallocate( values )
    set%cached => new_values
end function create_pair_set

function pair( set, idx )
    type(INTEGER_PAIR_SET), intent(in) :: set
    integer, intent(in)                :: idx
    integer, dimension(2)              :: pair

    pair = set%cached(2*idx-1:2*idx)
end function pair

integer function number_elements( set )
    type(INTEGER_PAIR_SET), intent(in) :: set

    number_elements = size(set%cached)/2
end function number_elements

logical function has_element( set, pair )
    type(INTEGER_PAIR_SET), intent(inout) :: set
    integer, dimension(:),  intent(in)    :: pair

    set%first  = pair(1)
    set%second = pair(2)

    has_element = integer_relation_eval(set%relation)

end function has_element

function range( operand, min_value, max_value )
    type(INTEGER_OPERAND), target      :: operand
    integer, intent(in)                :: min_value
    integer, intent(in)                :: max_value
    type(INTEGER_RANGE)                :: range

    range%value     => operand
    range%min_value =  min_value
    range%max_value =  max_value
end function range

subroutine integer_set_value( x, value )
    type(INTEGER_OPERAND), intent(out) :: x
    integer, intent(in)                :: value

    x%operation = 0
    x%value     = value
    x%first     => null()
    x%second    => null()
end subroutine integer_set_value

function integer_add( x, y ) result(add)
    type(INTEGER_OPERAND), intent(in), target  :: x
    type(INTEGER_OPERAND), intent(in), target  :: y
    type(INTEGER_OPERAND), pointer             :: add

    allocate( add )

    add%operation =  1
    add%first     => x
    add%second    => y
end function integer_add

function integer_add_v( x, y ) result(add)
    type(INTEGER_OPERAND), intent(in), target  :: x
    integer, intent(in)                        :: y
    type(INTEGER_OPERAND), pointer             :: add

    type(INTEGER_OPERAND), pointer             :: yy

    allocate( yy )
    yy = y

    allocate( add )

    add%operation =  1
    add%first     => x
    add%second    => yy
end function integer_add_v

function integer_add_vv( x, y ) result(add)
    integer, intent(in)                        :: x
    type(INTEGER_OPERAND), intent(in), target  :: y
    type(INTEGER_OPERAND), pointer             :: add

    type(INTEGER_OPERAND), pointer             :: xx

    allocate( xx )
    xx = x

    allocate( add )

    add%operation =  1
    add%first     => xx
    add%second    => y
end function integer_add_vv

function integer_multiply( x, y ) result(multiply)
    type(INTEGER_OPERAND), intent(in), target  :: x
    type(INTEGER_OPERAND), intent(in), target  :: y
    type(INTEGER_OPERAND), pointer             :: multiply

    allocate( multiply )

    multiply%operation =  2
    multiply%first     => x
    multiply%second    => y
end function integer_multiply

function integer_multiply_v( x, y ) result(multiply)
    type(INTEGER_OPERAND), intent(in), target  :: x
    integer, intent(in)                        :: y
    type(INTEGER_OPERAND), pointer             :: multiply

    type(INTEGER_OPERAND), pointer             :: yy

    allocate( yy )
    yy = y

    allocate( multiply )

    multiply%operation =  2
    multiply%first     => x
    multiply%second    => yy
end function integer_multiply_v

function integer_multiply_vv( x, y ) result(multiply)
    integer, intent(in)                        :: x
    type(INTEGER_OPERAND), intent(in), target  :: y
    type(INTEGER_OPERAND), pointer             :: multiply

    type(INTEGER_OPERAND), pointer             :: xx

    allocate( xx )
    xx = x

    allocate( multiply )

    multiply%operation =  2
    multiply%first     => xx
    multiply%second    => y
end function integer_multiply_vv

function integer_exponentiate_v( x, y ) result(exponentiate)
    type(INTEGER_OPERAND), intent(in), target  :: x
    integer, intent(in)                        :: y
    type(INTEGER_OPERAND), pointer             :: exponentiate

    type(INTEGER_OPERAND), pointer             :: yy

    allocate( yy )
    yy = y

    allocate( exponentiate )

    exponentiate%operation =  3
    exponentiate%first     => x
    exponentiate%second    => yy
end function integer_exponentiate_v

function integer_equal( x, y ) result(equal)
    type(INTEGER_OPERAND), intent(in), target   :: x
    type(INTEGER_OPERAND), intent(in), target   :: y
    type(INTEGER_RELATION), pointer             :: equal

    allocate( equal )

    equal%relation =  1
    equal%first    => x
    equal%second   => y
end function integer_equal

function integer_equal_v( x, y ) result(equal)
    type(INTEGER_OPERAND), intent(in), target   :: x
    integer, intent(in)                         :: y
    type(INTEGER_RELATION), pointer             :: equal

    type(INTEGER_OPERAND), pointer  :: yy

    allocate( equal )
    allocate( yy    )

    yy = y

    equal%relation =  1
    equal%first    => x
    equal%second   => yy
end function integer_equal_v

recursive subroutine integer_eval( x )
    type(INTEGER_OPERAND)   :: x

    if ( associated( x%first  ) ) call integer_eval( x%first  )
    if ( associated( x%second ) ) call integer_eval( x%second )

    select case( x%operation )
    case ( 0 )
        ! Nothing to be done

    case ( 1 )
        x%value = x%first%value + x%second%value

    case ( 2 )
        x%value = x%first%value * x%second%value

    case ( 3 )
        x%value = x%first%value ** x%second%value

    case default
        ! Nothing to be done
    end select

end subroutine integer_eval

function integer_relation_eval( relation ) result(value)
    type(INTEGER_RELATION)  :: relation
    logical                 :: value

    call integer_eval( relation%first  )
    call integer_eval( relation%second )

    select case( relation%relation )
        case ( 1 )
            value = relation%first%value == relation%second%value
        case default
            !
            ! An unknown or unimplemented operation:
            ! take the easy way out
            !
            stop 'Sorry, operation is not supported'
    end select

end function integer_relation_eval

end module integer_sets

program test_sets
    use integer_sets

    type(INTEGER_OPERAND), target   :: x
    type(INTEGER_OPERAND), target   :: y

    type(INTEGER_RELATION), pointer :: relation
    type(INTEGER_PAIR_SET)          :: set
    integer, dimension(2)           :: values

    relation => x + y == 0

    x = 1
    y = -1
    write(*,*) 'x, y: 1, -1 ', integer_relation_eval( relation )

    x = 2
    y = -1
    write(*,*) 'x, y: 2, -1 ', integer_relation_eval( relation )

    write(*,*) 'First few solutions to Pell''s equation y**2 = 2*x**2 + 1'

    set = create_pair_set( y**2 == 2*x**2 + 1, range(x, 0, 1000), &
                                               range(y, 0, 1000)  )

    do i = 1,number_elements(set)
        values = pair(set, i)
        write(*, '(i4,'': ('',i0,'','',i0,'')'')' ) i, values
    enddo

    write(*,*) 'First few solutions to equation y**2 = 2*x**2 - 1'

    set = create_pair_set( y**2 == 2*x**2 + (-1), range(x, 0, 1000), &
                                                  range(y, 0, 1000)  )

    do i = 1,number_elements(set)
        values = pair(set, i)
        write(*, '(i4,'': ('',i0,'','',i0,'')'')' ) i, values
    enddo

!
! Yet another one: to show that the relationship persists
!
    write(*,*) 'Does the set {(x,y) | y**2 = 3*x**2 + 1} contain (1,2) or (3,3)?'

    set = create_pair_set( y**2 == 3*x**2 + 1, range(x, 0, 0), &
                                               range(y, 0, 0)  )

    write(*,*) '(1,2) in this set? ', has_element(set, (/1,2/))
    write(*,*) '(3,3) in this set? ', has_element(set, (/3,3/))
end program
