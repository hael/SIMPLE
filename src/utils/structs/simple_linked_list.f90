module simple_linked_list
implicit none

public :: linked_list, list_iterator
private
#include "simple_local_flags.inc"

type :: node
    class(*), allocatable :: value
    type(node), pointer   :: next => null()
contains
    final :: kill_node
end type node

type :: linked_list
    private
    type(node), pointer :: head => null()
    type(node), pointer :: tail => null()
    integer             :: n    = 0
contains
    ! construction/lifecycle
    !--push/pop
    procedure          :: push_front
    procedure          :: push_back
    procedure          :: pop_front
    !--destructor
    procedure          :: kill
    ! accessors
    procedure          :: front
    procedure          :: back
    procedure          :: at
    ! checkers
    procedure          :: size
    procedure          :: is_empty
    ! list operations
    procedure, private :: assign
    generic            :: assignment(=) => assign
    procedure, private :: append
    generic            :: operator(//)  => append
    procedure          :: slice
    procedure          :: replace_with
    ! iteraton
    procedure          :: begin
    procedure          :: end_iter
    ! type-safe API (wrappers)
    procedure          :: push_back_int
    procedure          :: push_back_real
    procedure          :: push_back_complex
    procedure          :: push_back_logical
    procedure          :: push_back_char
    procedure          :: front_int
    procedure          :: at_int
    procedure          :: front_char
    procedure          :: at_char
end type linked_list

type :: list_iterator
    private
    type(node), pointer :: current => null()
contains
    procedure :: has_value
    procedure :: next
    procedure :: get
    procedure :: set
    procedure :: equals
    procedure :: advance
    procedure :: index
end type list_iterator

contains

    !======================
    ! construction/lifecycle
    !======================

    subroutine push_front(self, x)
        class(linked_list), intent(inout) :: self
        class(*),           intent(in)    :: x
        type(node), pointer :: p
        allocate(p)
        allocate(p%value, source=x)
        p%next => self%head
        self%head => p
        if (.not. associated(self%tail)) self%tail => p
        self%n = self%n + 1
    end subroutine push_front

    subroutine push_back(self, x)
        class(linked_list), intent(inout) :: self
        class(*),           intent(in)    :: x
        type(node), pointer :: p
        allocate(p)
        allocate(p%value, source=x)
        nullify(p%next)
        if (associated(self%tail)) then
            self%tail%next => p
        else
            self%head => p
        end if
        self%tail => p
        self%n = self%n + 1
    end subroutine push_back

    subroutine pop_front(self, x, had_value)
        class(linked_list),    intent(inout) :: self
        class(*), allocatable, intent(inout) :: x
        logical, optional,     intent(out)   :: had_value
        type(node), pointer :: p
        if (.not. associated(self%head)) then
            if (present(had_value)) had_value = .false.
            return
        end if
        p => self%head
        self%head => p%next
        if (.not. associated(self%head)) nullify(self%tail)
        if (allocated(p%value)) then
            if( allocated(x) ) deallocate(x)
            allocate(x, source=p%value)
            if (present(had_value)) had_value = .true.
        else
            if (present(had_value)) had_value = .false.
        end if
        call kill_node(p)
        deallocate(p)
        self%n = self%n - 1
    end subroutine pop_front

    ! helper
    pure subroutine kill_node(self)
        type(node), intent(inout) :: self
        if (allocated(self%value)) deallocate(self%value)
        nullify(self%next)
    end subroutine kill_node

    elemental subroutine kill(self)
        class(linked_list), intent(inout) :: self
        type(node), pointer :: p, q
        p => self%head
        do while (associated(p))
            q => p%next
            call kill_node(p)  ! dealloc value & nullify next
            deallocate(p)
            p => q
        end do
        nullify(self%head, self%tail)
        self%n = 0
    end subroutine kill

    !======================
    ! Accessors
    !======================

    subroutine front(self, x)
        class(linked_list),    intent(in)    :: self
        class(*), allocatable, intent(inout) :: x
        if (.not. associated(self%head)) stop "front(): list is empty"
        if( allocated(x) ) deallocate(x)
        allocate(x, source=self%head%value)
    end subroutine front

    subroutine back(self, x)
        class(linked_list),     intent(in)    :: self
        class(*), allocatable,  intent(inout) :: x
        if (.not. associated(self%tail)) stop "back(): list is empty"
        if( allocated(x) ) deallocate(x)
        allocate(x, source=self%tail%value)
    end subroutine back

    subroutine at(self, idx, x)
        class(linked_list),    intent(in)    :: self
        integer,               intent(in)    :: idx   ! 1-based index
        class(*), allocatable, intent(inout) :: x
        type(node), pointer :: p
        integer :: i
        if (idx < 1 .or. idx > self%n) stop "at(): index out of range"
        p => self%head
        do i = 2, idx
            p => p%next
        end do
        if( allocated(x) ) deallocate(x)
        allocate(x, source=p%value)
    end subroutine at

    !======================
    ! Checkers
    !======================

    pure elemental integer function size(self) result(k)
        class(linked_list), intent(in) :: self
        k = self%n
    end function size

    pure logical function is_empty(self) result(tf)
        class(linked_list), intent(in) :: self
        tf = (self%n == 0)
    end function is_empty

    !======================
    ! List operations
    !======================

    ! deep polymorphic copy
    subroutine assign(self, other)
        class(linked_list), intent(inout) :: self
        class(linked_list), intent(in)    :: other
        type(node), pointer :: p
        if (associated(self%head)) call self%kill()
        p => other%head
        do while (associated(p))
            call self%push_back(p%value)
            p => p%next
        end do
    end subroutine assign

    function append(self, other) result(res)
        class(linked_list), intent(in) :: self, other
        type(linked_list)              :: res
        class(*),   allocatable        :: val
        type(node), pointer            :: p
        ! Copy first list
        p => self%head
        do while (associated(p))
            call res%push_back(p%value)
            p => p%next
        end do
        ! Append second
        p => other%head
        do while (associated(p))
            call res%push_back(p%value)
            p => p%next
        end do
    end function append

    
    !======================
    ! Slice/append
    !======================

    subroutine slice(self, i1, i2, this)
        class(linked_list), intent(in)    :: self
        integer,             intent(in)   :: i1, i2
        type(linked_list),   intent(out)  :: this
        type(node), pointer :: p
        integer :: i
        if (i1 < 1 .or. i2 > self%n .or. i1 > i2) stop "slice(): invalid range"
        call this%kill()
        p => self%head
        do i = 1, i1 - 1
            p => p%next
        end do
        do i = i1, i2
            call this%push_back(p%value)
            p => p%next
        end do
    end subroutine slice

    !======================
    ! Move semantics (steal nodes from other)
    !======================

    subroutine replace_with(self, other)
        class(linked_list), intent(inout) :: self
        class(linked_list), intent(inout) :: other
        ! Move nodes from other into self, leaving other empty.
        ! This is O(1) (pointer steal).
        if (associated(self%head)) call self%kill()
        if (.not. associated(other%head)) then
            ! other empty -> leave self empty
            nullify(self%head, self%tail)
            self%n = 0
            return
        end if
        ! steal descriptor
        self%head => other%head
        self%tail => other%tail
        self%n = other%n
        ! empty other
        nullify(other%head, other%tail)
        other%n = 0
    end subroutine replace_with

    !======================
    ! Iteration
    !======================

    function begin(self) result(it)
        class(linked_list), intent(in) :: self
        type(list_iterator) :: it
        it%current => self%head
    end function begin

    function end_iter(self) result(it)
        class(linked_list), intent(in) :: self
        type(list_iterator) :: it
        nullify(it%current)  ! matches C++ end(): a special sentinel
    end function end_iter

    pure logical function has_value(self) result(tf)
        class(list_iterator), intent(in) :: self
        tf = associated(self%current)
    end function has_value

    subroutine next(self)
        class(list_iterator), intent(inout) :: self
        if (.not. associated(self%current)) return
        self%current => self%current%next
    end subroutine next

    subroutine get(self, x)
        class(list_iterator),  intent(in)    :: self
        class(*), allocatable, intent(inout) :: x
        if (.not. associated(self%current)) stop "iterator: dereferencing end()"
        if( allocated(x) ) deallocate(x)
        allocate(x, source=self%current%value)
    end subroutine get

     subroutine set(self, x)
        class(list_iterator), intent(inout) :: self
        class(*),             intent(in)    :: x
        if (.not. associated(self%current)) stop "iterator: dereferencing end()"
        if( allocated(self%current%value) ) deallocate(self%current%value)
        allocate(self%current%value, source=x)
    end subroutine set

    pure logical function equals(self, other) result(same)
        class(list_iterator), intent(in) :: self, other
        same = associated(self%current, other%current)
    end function equals

    subroutine advance(self, n)
        class(list_iterator), intent(inout) :: self
        integer, intent(in) :: n
        integer :: i
        do i = 1, n
            if (.not. associated(self%current)) exit
            self%current => self%current%next
        end do
    end subroutine advance

    integer function index(self, lst) result(i)
        class(list_iterator), intent(in) :: self
        class(linked_list),   intent(in) :: lst
        type(node), pointer :: p
        i = 1
        p => lst%head
        do while (associated(p))
            if (associated(p, self%current)) return
            i = i + 1
            p => p%next
        end do
        i = -1  ! not found
    end function index

    !======================
    ! Type-safe wrapper APIs
    ! (thin wrappers calling polymorphic methods)
    !======================

    subroutine push_back_int(self, v)
        class(linked_list), intent(inout) :: self
        integer,            intent(in)    :: v
        call self%push_back(v)
    end subroutine push_back_int

    subroutine push_back_real(self, v)
        class(linked_list), intent(inout) :: self
        real,               intent(in)    :: v
        call self%push_back(v)
    end subroutine push_back_real

    subroutine push_back_complex(self, v)
        class(linked_list), intent(inout) :: self
        complex,            intent(in)    :: v
        call self%push_back(v)
    end subroutine push_back_complex

    subroutine push_back_logical(self, v)
        class(linked_list), intent(inout) :: self
        logical,            intent(in)    :: v
        call self%push_back(v)
    end subroutine push_back_logical

    subroutine push_back_char(self, v)
        class(linked_list), intent(inout) :: self
        character(len=*),   intent(in)    :: v
        character(:), allocatable :: tmp
        tmp = v
        call self%push_back(tmp)
    end subroutine push_back_char

    subroutine front_int(self, val)
        class(linked_list), intent(in)  :: self
        integer,            intent(out) :: val
        class(*), allocatable :: tmp
        if (.not. associated(self%head)) stop "front_int: empty"
        call self%front(tmp)
        select type(t => tmp)
        type is (integer)
            val = t
        class default
            stop "front_int: stored element not integer"
        end select
        if (allocated(tmp)) deallocate(tmp)
    end subroutine front_int

    subroutine at_int(self, idx, val)
        class(linked_list), intent(in)  :: self
        integer,            intent(in)  :: idx
        integer,            intent(out) :: val
        class(*), allocatable :: tmp
        call self%at(idx, tmp)
        select type(t => tmp)
        type is (integer)
            val = t
        class default
            stop "at_int: stored element not integer"
        end select
        if (allocated(tmp)) deallocate(tmp)
    end subroutine at_int

    subroutine front_char(self, val)
        class(linked_list),        intent(in)  :: self
        character(:), allocatable, intent(out) :: val
        class(*), allocatable :: tmp
        call self%front(tmp)
        select type(t => tmp)
        type is (character(*))
            val = t
        class default
            stop "front_char: stored element not character"
        end select
        if (allocated(tmp)) deallocate(tmp)
    end subroutine front_char

    subroutine at_char(self, idx, val)
        class(linked_list),        intent(in)  :: self
        integer,                   intent(in)  :: idx
        character(:), allocatable, intent(out) :: val
        class(*), allocatable :: tmp
        call self%at(idx, tmp)
        select type(t => tmp)
        type is (character(*))
            val = t
        class default
            stop "at_char: stored element not character"
        end select
        if (allocated(tmp)) deallocate(tmp)
    end subroutine at_char

end module simple_linked_list
