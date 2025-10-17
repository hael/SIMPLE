module simple_linked_list
implicit none

public :: linked_list, list_iterator
private

type :: node
    class(*), allocatable :: value
    type(node), pointer   :: next => null()
contains
    final :: finalize_node
end type node

type :: linked_list
    private
    type(node), pointer :: head => null()
    type(node), pointer :: tail => null()
    integer             :: n    = 0
contains
    procedure :: push_front
    procedure :: push_back
    procedure :: pop_front
    procedure :: front
    procedure :: back
    procedure :: at
    procedure :: clear
    procedure :: size
    procedure :: is_empty
    procedure :: begin ! iterator to first
    final     :: finalize_list
end type linked_list

type :: list_iterator
    private
    type(node), pointer :: current => null()
contains
    procedure :: has_next
    procedure :: next  => iter_next
end type list_iterator

contains

    pure subroutine finalize_node(self)
        type(node), intent(inout) :: self
        if (allocated(self%value)) deallocate(self%value)
        nullify(self%next)
    end subroutine finalize_node

    pure elemental integer function size(self) result(k)
        class(linked_list), intent(in) :: self
        k = self%n
    end function size

    pure logical function is_empty(self) result(tf)
        class(linked_list), intent(in) :: self
        tf = (self%n == 0)
    end function is_empty

    elemental subroutine clear(self)
        class(linked_list), intent(inout) :: self
        type(node), pointer :: p, q
        p => self%head
        do while (associated(p))
            q => p%next
            call finalize_node(p)  ! dealloc value & nullify next
            deallocate(p)
            p => q
        end do
        nullify(self%head, self%tail)
        self%n = 0
    end subroutine clear

    subroutine finalize_list(self)
        type(linked_list), intent(inout) :: self
        call self%clear()
    end subroutine finalize_list

    !======================
    ! Push operations
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
        class(*), intent(in)             :: x
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

    !======================
    ! Pop / access
    !======================

    subroutine pop_front(self, x, had_value)
        class(linked_list),    intent(inout) :: self
        class(*), allocatable, intent(out)   :: x
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
            allocate(x, source=p%value)
            if (present(had_value)) had_value = .true.
        else
            if (present(had_value)) had_value = .false.
        end if
        call finalize_node(p)
        deallocate(p)
        self%n = self%n - 1
    end subroutine pop_front

    subroutine front(self, x)
        class(linked_list),    intent(in)  :: self
        class(*), allocatable, intent(out) :: x
        if (.not. associated(self%head)) stop "front(): list is empty"
        allocate(x, source=self%head%value)
    end subroutine front

    subroutine back(self, x)
        class(linked_list),    intent(in)  :: self
        class(*), allocatable, intent(out) :: x
        if (.not. associated(self%tail)) stop "back(): list is empty"
        allocate(x, source=self%tail%value)
    end subroutine back

    subroutine at(self, idx, x)
        class(linked_list),    intent(in)  :: self
        integer,               intent(in)  :: idx   ! 1-based index
        class(*), allocatable, intent(out) :: x
        type(node), pointer :: p
        integer :: i
        if (idx < 1 .or. idx > self%n) stop "at(): index out of range"
        p => self%head
        do i = 2, idx
            p => p%next
        end do
        allocate(x, source=p%value)
    end subroutine at

    !======================
    ! Iteration
    !======================

    function begin(self) result(it)
        class(linked_list), intent(in) :: self
        type(list_iterator) :: it
        it%current => self%head
    end function begin

    pure logical function has_next(self) result(tf)
        class(list_iterator), intent(in) :: self
        tf = associated(self%current)
    end function has_next

    subroutine iter_next(self, x)
        class(list_iterator),  intent(inout) :: self
        class(*), allocatable, intent(out)   :: x
        if (.not. associated(self%current)) stop "iterator: no more elements"
        allocate(x, source=self%current%value)
        self%current => self%current%next
    end subroutine iter_next

end module simple_linked_list
        