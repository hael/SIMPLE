module simple_linked_list
include 'simple_lib.f08'
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
    procedure          :: push_back_int
    procedure          :: push_back_real
    procedure          :: push_back_complex
    procedure          :: push_back_logical
    procedure          :: push_back_char
    procedure          :: pop_front
    !--destructor
    procedure          :: kill
    ! accessors
    procedure          :: front
    procedure          :: back
    procedure          :: at
    procedure          :: at_int
    procedure          :: at_real
    procedure          :: at_complex
    procedure          :: at_logical
    procedure          :: at_char
    ! checkers
    procedure          :: size
    procedure          :: is_empty
    ! list operations
    procedure, private :: assign
    generic            :: assignment(=) => assign
    procedure, private :: append
    generic            :: operator(//)  => append
    procedure          :: slice
    procedure          :: replace_at
    procedure          :: replace_iterator
    procedure          :: replace_with
    ! iteraton
    procedure          :: begin
    procedure          :: end_iter
end type linked_list

type :: list_iterator
    private
    type(node), pointer :: current => null()
contains
    procedure :: has_value
    procedure :: next
    procedure :: getter
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
        if (.not. associated(self%head)) THROW_HARD('list is empty')
        if( allocated(x) ) deallocate(x)
        allocate(x, source=self%head%value)
    end subroutine front

    subroutine back(self, x)
        class(linked_list),     intent(in)    :: self
        class(*), allocatable,  intent(inout) :: x
        if (.not. associated(self%tail)) THROW_HARD('list is empty')
        if( allocated(x) ) deallocate(x)
        allocate(x, source=self%tail%value)
    end subroutine back

    subroutine at(self, idx, x)
        class(linked_list),    intent(in)    :: self
        integer,               intent(in)    :: idx   ! 1-based index
        class(*), allocatable, intent(inout) :: x
        type(node), pointer :: p
        integer :: i
        if (idx < 1 .or. idx > self%n) THROW_HARD('index out of range')
        p => self%head
        do i = 2, idx
            p => p%next
        end do
        if( allocated(x) ) deallocate(x)
        allocate(x, source=p%value)
    end subroutine at

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
            THROW_HARD('not an integer type')
        end select
        if (allocated(tmp)) deallocate(tmp)
    end subroutine at_int

    subroutine at_real(self, idx, val)
        class(linked_list), intent(in)  :: self
        integer,            intent(in)  :: idx
        real,               intent(out) :: val
        class(*), allocatable :: tmp
        call self%at(idx, tmp)
        select type(t => tmp)
        type is (real)
            val = t
        class default
            THROW_HARD('not a real type')
        end select
        if (allocated(tmp)) deallocate(tmp)
    end subroutine at_real

    subroutine at_complex(self, idx, val)
        class(linked_list), intent(in)  :: self
        integer,            intent(in)  :: idx
        complex,            intent(out) :: val
        class(*), allocatable :: tmp
        call self%at(idx, tmp)
        select type(t => tmp)
        type is (complex)
            val = t
        class default
            THROW_HARD('not a complex type')
        end select
        if (allocated(tmp)) deallocate(tmp)
    end subroutine at_complex

    subroutine at_logical(self, idx, val)
        class(linked_list), intent(in)  :: self
        integer,            intent(in)  :: idx
        logical,            intent(out) :: val
        class(*), allocatable :: tmp
        call self%at(idx, tmp)
        select type(t => tmp)
        type is (logical)
            val = t
        class default
            THROW_HARD('not a logical type')
        end select
        if (allocated(tmp)) deallocate(tmp)
    end subroutine at_logical

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
            THROW_HARD('not a character type')
        end select
        if (allocated(tmp)) deallocate(tmp)
    end subroutine at_char

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
        if (i1 < 1 .or. i2 > self%n .or. i1 > i2) THROW_HARD('invalid range')
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

    !----------------------------------------------------------------
    ! Replace node at 1-based index `idx` with a new node containing `x`.
    ! Preserves list length and updates head/tail as necessary.
    ! Safe for changing dynamic type of the polymorphic allocatable value.
    !----------------------------------------------------------------
    subroutine replace_at(self, idx, x)
        class(linked_list), intent(inout) :: self
        integer,            intent(in)    :: idx
        class(*),           intent(in)    :: x
        type(node), pointer :: old, prev, succ, pnew, p
        integer :: i
        if (idx < 1 .or. idx > self%n) then
            THROW_HARD('replace_at: index out of range')
        end if
        ! locate node to replace and previous node
        prev => null()
        old  => self%head
        do i = 2, idx
            prev => old
            old  => old%next
        end do
        ! now old => node at idx, prev => previous (null if idx==1)
        ! allocate new node and set its next pointer (will point to old%next)
        allocate(pnew)
        pnew%next => null()
        ! capture successor before we destroy old
        succ => old%next
        pnew%next => succ
        ! allocate polymorphic storage with correct type
        allocate(pnew%value, source=x)
        ! splice new node into the list in place of old
        if (associated(prev)) then
            prev%next => pnew
        else
            ! replacing head
            self%head => pnew
        end if
        ! update tail if we replaced the previous tail node
        if (associated(self%tail, old)) then
            self%tail => pnew
        end if
        ! destroy old node (value and pointer) then deallocate it
        call kill_node(old)   ! this will deallocate old%value and nullify old%next
        deallocate(old)
        ! list length unchanged
    end subroutine replace_at

    !----------------------------------------------------------------
    ! Replace the element referenced by an iterator with a new value.
    ! The iterator will remain valid and now point to the replacement.
    !----------------------------------------------------------------
    subroutine replace_iterator(self, it, x)
        class(linked_list),  intent(inout) :: self
        type(list_iterator), intent(inout) :: it
        class(*),            intent(in)    :: x
        type(node), pointer :: old, prev, succ, pnew, p
        if (.not. associated(it%current)) then
            THROW_HARD("replace_iterator: iterator does not reference a node")
        end if
        old => it%current
        ! locate previous node for relinking
        prev => null()
        p    => self%head
        if (associated(p, old)) then
            prev => null()      ! old == head
        else
            do while (associated(p) .and. .not.associated(p%next, old))
                p => p%next
            end do
            if (.not.associated(p)) then
                THROW_HARD("replace_iterator: node not found in list (corrupted iterator)")
            end if
            prev => p
        end if
        ! successor before old is destroyed
        succ => old%next
        ! allocate new node
        allocate(pnew)
        pnew%next => succ
        ! allocate polymorphic storage with correct type
        allocate(pnew%value, source=x)
        ! link new node into list
        if (associated(prev)) then
            prev%next => pnew
        else
            self%head => pnew
        end if
        ! update tail if needed
        if (associated(self%tail, old)) then
            self%tail => pnew
        end if
        ! destroy old node and free memory
        call kill_node(old)
        deallocate(old)
        ! iterator now moves to the replacement node
        it%current => pnew
    end subroutine replace_iterator

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

    subroutine getter(self, x)
        class(list_iterator),  intent(in)    :: self
        class(*), allocatable, intent(inout) :: x
        if (.not. associated(self%current)) THROW_HARD('iterator: dereferencing end()')
        if( allocated(x) ) deallocate(x)
        allocate(x, source=self%current%value)
    end subroutine getter

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

end module simple_linked_list
