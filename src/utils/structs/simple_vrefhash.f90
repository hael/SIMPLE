!@descr: string-key, polymorphic *reference* hash (stores pointers; updates are visible)
!
! example usage of vrefhash:
!
! use simple_vrefhash
! use simple_ui_program, only: ui_program
!
! type(vrefhash) :: h
! type(ui_program), target :: abinitio2D
!
! class(*), pointer :: p
! logical :: ok
!
! call h%init()
!
! call h%set_ref("abinitio2D", abinitio2D)

! call h%get_ref("abinitio2D", p, ok)
! if (ok) then
!     select type(u => p)
!     type is (ui_program)
!         ! u is the SAME instance as abinitio2D (updates visible)
!         ! (call u%whatever ...)
!     class default
!         call simple_exception("wrong dynamic type in hash", __FILENAME__, __LINE__)
!     end select
! end if
!
module simple_vrefhash
use simple_error,  only: simple_exception
use simple_string, only: string
implicit none

public :: vrefhash

!-------------------------
! Internal entry node
!-------------------------
type :: vrefhash_node
    type(string)                 :: key
    class(*),            pointer :: p => null()     ! <--- reference to external TARGET
    type(vrefhash_node), pointer :: next => null()
contains
    final :: finalize_node
end type vrefhash_node

type :: vrefhash_bucket
    type(vrefhash_node), pointer :: head => null()
end type vrefhash_bucket

!-------------------------
! Hash table
!-------------------------
type :: vrefhash
    private
    type(vrefhash_bucket), allocatable :: buckets(:)
    integer :: nbuckets = 0
    integer :: n        = 0
    logical :: exists   = .false.
contains
    ! lifecycle
    procedure          :: init
    procedure          :: destroy
    procedure          :: clear
    ! set/get/del (reference semantics)
    procedure, private :: set_ref_char
    procedure, private :: set_ref_str
    generic            :: set_ref => set_ref_char, set_ref_str
    procedure, private :: get_ref_char
    procedure, private :: get_ref_str
    generic            :: get_ref => get_ref_char, get_ref_str
    procedure, private :: del_char
    procedure, private :: del_str
    generic            :: del => del_char, del_str
    ! queries
    procedure, private :: has_key_char
    procedure, private :: has_key_str
    generic            :: has_key => has_key_char, has_key_str
    procedure          :: count
end type vrefhash

contains

    !=====================
    ! Finalizer / helpers
    !=====================

    subroutine finalize_node(self)
        type(vrefhash_node), intent(inout) :: self
        nullify(self%p)
        nullify(self%next)
        call self%key%kill()
    end subroutine finalize_node

    integer function fnv1a_hash_bytes(s) result(h)
        character(len=*), intent(in) :: s
        integer :: i, c
        integer(kind=8) :: x
        ! 64-bit FNV-1a reduced to default integer via mod
        x = int(Z'CBF29CE484222325', kind=8)
        do i = 1, len_trim(s)
            c = iachar(s(i:i))
            x = ieor(x, int(c,kind=8))
            x = x * int(Z'00000100000001B3', kind=8)
        end do
        ! map to positive default integer range
        h = int( iand(x, int(Z'7FFFFFFF', kind=8)) )
    end function fnv1a_hash_bytes

    integer function bucket_index(self, k) result(b)
        class(vrefhash),  intent(in) :: self
        character(len=*), intent(in) :: k
        integer :: hv
        if (.not. self%exists) then
            call simple_exception('vrefhash not initialized', __FILENAME__, __LINE__)
        end if
        hv = fnv1a_hash_bytes(k)
        b  = 1 + mod(hv, self%nbuckets)
    end function bucket_index

    subroutine find_node(self, b, k, prev, cur)
        class(vrefhash),               intent(in)  :: self
        integer,                       intent(in)  :: b
        character(len=*),              intent(in)  :: k
        type(vrefhash_node), pointer,  intent(out) :: prev
        type(vrefhash_node), pointer,  intent(out) :: cur
        prev => null()
        cur  => self%buckets(b)%head
        do while (associated(cur))
            if (cur%key == k) return
            prev => cur
            cur  => cur%next
        end do
    end subroutine find_node

    !=====================
    ! lifecycle
    !=====================

    subroutine init(self, nbuckets)
        class(vrefhash), intent(inout) :: self
        integer, optional, intent(in)  :: nbuckets
        integer :: nloc, i
        call self%destroy()
        if (present(nbuckets)) then
            nloc = max(8, nbuckets)
        else
            nloc = 256
        end if
        allocate(self%buckets(nloc))
        do i = 1, nloc
            nullify(self%buckets(i)%head)
        end do
        self%nbuckets = nloc
        self%n        = 0
        self%exists   = .true.
    end subroutine init

    subroutine clear(self)
        class(vrefhash), intent(inout) :: self
        integer :: i
        type(vrefhash_node), pointer :: p, q
        if (.not. self%exists) return
        do i = 1, self%nbuckets
            p => self%buckets(i)%head
            do while (associated(p))
                q => p%next
                call finalize_node(p)
                deallocate(p)
                p => q
            end do
            nullify(self%buckets(i)%head)
        end do
        self%n = 0
    end subroutine clear

    subroutine destroy(self)
        class(vrefhash), intent(inout) :: self
        if (.not. self%exists) return
        call self%clear()
        if (allocated(self%buckets)) deallocate(self%buckets)
        self%nbuckets = 0
        self%n        = 0
        self%exists   = .false.
    end subroutine destroy

    !=====================
    ! queries
    !=====================

    integer function count(self) result(k)
        class(vrefhash), intent(in) :: self
        k = self%n
    end function count

    logical function has_key_char(self, key) result(tf)
        class(vrefhash),     intent(in) :: self
        character(len=*),    intent(in) :: key
        type(vrefhash_node), pointer    :: prev, cur
        character(len=:), allocatable   :: k
        integer :: b
        if (.not. self%exists) then
            tf = .false.; return
        end if
        k = trim(adjustl(key))
        b = bucket_index(self, k)
        call find_node(self, b, k, prev, cur)
        tf = associated(cur)
    end function has_key_char

    logical function has_key_str(self, key) result(tf)
        class(vrefhash), intent(in) :: self
        class(string),   intent(in) :: key
        character(len=:), allocatable :: k
        type(vrefhash_node), pointer :: prev, cur
        integer :: b
        if (.not. self%exists) then
            tf = .false.; return
        end if
        k = adjustl(key%to_char())
        b = bucket_index(self, k)
        call find_node(self, b, k, prev, cur)
        tf = associated(cur)
    end function has_key_str

    !=====================
    ! set_ref (store pointer)
    !=====================

    subroutine set_ref_char(self, key, val)
        class(vrefhash),               intent(inout) :: self
        character(len=*),              intent(in)    :: key
        class(*), target,              intent(inout) :: val
        character(len=:), allocatable :: k
        integer :: b
        type(vrefhash_node), pointer :: prev, cur, pnew
        if (.not. self%exists) call self%init()
        k = trim(adjustl(key))
        b = bucket_index(self, k)
        call find_node(self, b, k, prev, cur)
        if (associated(cur)) then
            ! replace pointer (reference semantics)
            cur%p => val
            return
        end if
        ! insert new node at head (cheap)
        allocate(pnew)
        pnew%key = k
        pnew%p   => val
        pnew%next => self%buckets(b)%head
        self%buckets(b)%head => pnew
        self%n = self%n + 1
    end subroutine set_ref_char

    subroutine set_ref_str(self, key, val)
        class(vrefhash),               intent(inout) :: self
        class(string),                 intent(in)    :: key
        class(*), target,              intent(inout) :: val
        character(len=:), allocatable :: k
        integer :: b
        type(vrefhash_node), pointer :: prev, cur, pnew
        if (.not. self%exists) call self%init()
        k = adjustl(key%to_char())
        b = bucket_index(self, k)
        call find_node(self, b, k, prev, cur)
        if (associated(cur)) then
            cur%p => val
            return
        end if
        allocate(pnew)
        pnew%key = k
        pnew%p   => val
        pnew%next => self%buckets(b)%head
        self%buckets(b)%head => pnew
        self%n = self%n + 1
    end subroutine set_ref_str

    !=====================
    ! get_ref (return pointer)
    !=====================

    subroutine get_ref_char(self, key, p, found)
        class(vrefhash),              intent(in)  :: self
        character(len=*),             intent(in)  :: key
        class(*), pointer,            intent(out) :: p
        logical, optional,            intent(out) :: found
        character(len=:), allocatable :: k
        integer :: b
        type(vrefhash_node), pointer :: prev, cur
        nullify(p)
        if (present(found)) found = .false.
        if (.not. self%exists) return
        k = trim(adjustl(key))
        b = bucket_index(self, k)
        call find_node(self, b, k, prev, cur)
        if (associated(cur)) then
            p => cur%p
            if (present(found)) found = associated(p)
        end if
    end subroutine get_ref_char

    subroutine get_ref_str(self, key, p, found)
        class(vrefhash),              intent(in)  :: self
        class(string),                intent(in)  :: key
        class(*), pointer,            intent(out) :: p
        logical, optional,            intent(out) :: found
        character(len=:), allocatable :: k
        integer :: b
        type(vrefhash_node), pointer :: prev, cur
        nullify(p)
        if (present(found)) found = .false.
        if (.not. self%exists) return
        k = adjustl(key%to_char())
        b = bucket_index(self, k)
        call find_node(self, b, k, prev, cur)
        if (associated(cur)) then
            p => cur%p
            if (present(found)) found = associated(p)
        end if
    end subroutine get_ref_str

    !=====================
    ! del (remove key)
    !=====================

    subroutine del_char(self, key)
        class(vrefhash),  intent(inout) :: self
        character(len=*), intent(in)    :: key
        character(len=:), allocatable :: k
        integer :: b
        type(vrefhash_node), pointer :: prev, cur
        if (.not. self%exists) return
        k = trim(adjustl(key))
        b = bucket_index(self, k)
        call find_node(self, b, k, prev, cur)
        if (.not. associated(cur)) return
        if (associated(prev)) then
            prev%next => cur%next
        else
            self%buckets(b)%head => cur%next
        end if
        call finalize_node(cur)
        deallocate(cur)
        self%n = self%n - 1
    end subroutine del_char

    subroutine del_str(self, key)
        class(vrefhash), intent(inout) :: self
        class(string),   intent(in)    :: key
        character(len=:), allocatable :: k
        integer :: b
        type(vrefhash_node), pointer :: prev, cur
        if (.not. self%exists) return
        k = adjustl(key%to_char())
        b = bucket_index(self, k)
        call find_node(self, b, k, prev, cur)
        if (.not. associated(cur)) return
        if (associated(prev)) then
            prev%next => cur%next
        else
            self%buckets(b)%head => cur%next
        end if
        call finalize_node(cur)
        deallocate(cur)
        self%n = self%n - 1
    end subroutine del_str

end module simple_vrefhash
