module simple_rec_list
include 'simple_lib.f08'
use simple_linked_list
implicit none

public :: rec, project_rec, process_rec, chunk_rec
public :: rec_list, rec_iterator
private
#include "simple_local_flags.inc"

! ========= Base and derived record types =========

type, abstract :: rec
    integer :: id = 0 ! unique ID
contains
    procedure(rec_print), deferred :: print
end type rec

abstract interface
    subroutine rec_print(self)
        import rec
        class(rec), intent(in) :: self
    end subroutine rec_print
end interface

! Convenience type to hold information about individual project files
type, extends(rec) :: project_rec
    type(string) :: projname            ! project file name
    integer      :: micind     = 0      ! index of micrograph in project
    integer      :: nptcls     = 0      ! # of particles
    integer      :: nptcls_sel = 0      ! # of particles (state=1)
    logical      :: included   = .false.! whether the record has been imported
contains
    procedure :: print => print_proj
end type project_rec

! Convenience type to hold information about processes
type, extends(rec) :: process_rec
    type(string) :: str_id              ! unique string ID
    type(string) :: folder              ! location
    type(string) :: projfile            ! project filename
    type(string) :: volume              ! volume filename
    type(string) :: alnvolume           ! aligned volume filename
    logical      :: submitted = .false. ! process has been submitted (running)
    logical      :: completed = .false. ! has completed
    logical      :: included  = .false. ! whether the record has been post-processed/analyzed
contains
    procedure :: print => print_proc
end type process_rec

! Convenience type to keep track of converged chunks
type, extends(rec) :: chunk_rec
    type(string) :: projfile            ! project filename
    logical      :: busy      = .false. ! true after submission and until completion is detected
    logical      :: processed = .false. ! chunk: has converged; set: has been clustered/selected/matched
    logical      :: included  = .false. ! whether the set has been imported into the pool
contains
    procedure :: print => print_file
end type chunk_rec

! ========= Type-safe container wrapper =========

type :: rec_iterator
    private
    type(list_iterator) :: it
contains
    procedure          :: valid        => iter_valid
    procedure          :: next         => iter_next
    procedure          :: advance      => iter_advance
    procedure, private :: iter_get_rec, iter_get_project_rec, iter_get_process_rec, iter_get_chunk_rec
    generic            :: get => iter_get_project_rec, iter_get_process_rec, iter_get_chunk_rec
    procedure          :: replace      => iter_replace
    procedure          :: equals       => iter_equals
    generic            :: operator(==) => equals
end type rec_iterator

type :: rec_list
    private
    type(linked_list) :: list
contains
    ! Modifiers
    procedure          :: push_back
    procedure          :: push2project_list
    procedure          :: push2chunk_list
    procedure          :: replace_at
    procedure          :: replace_iterator
    procedure          :: replace_with
    procedure, private :: append
    generic            :: operator(//) => append
    procedure          :: subset_to
    procedure          :: slice
    ! Accessors
    procedure          :: size
    procedure          :: at_rec
    procedure, private :: at_project_rec, at_process_rec, at_chunk_rec
    generic            :: at => at_project_rec, at_process_rec, at_chunk_rec
    procedure          :: begin
    procedure          :: end_iter
    ! Lifecycle
    procedure, private :: assign
    generic            :: assignment(=) => assign
    procedure          :: kill
    ! Specialized methods
    procedure          :: get_ids
    procedure          :: get_included_flags
    procedure          :: get_processed_flags
    procedure          :: get_busy_flags
    procedure          :: set_included_flags
    procedure          :: get_nptcls_tot
    procedure          :: get_nptcls_sel_tot
    procedure          :: get_projfiles
end type rec_list

contains

    ! -------------------- Printing --------------------

    subroutine print_proj(self)
        class(project_rec), intent(in) :: self
        print '("project_rec(id=",i0,", name=",a,")")', self%id, self%projname%to_char()
    end subroutine print_proj

    subroutine print_proc(self)
        class(process_rec), intent(in) :: self
        print '("process_rec(id=",i0,", id=",a,")")', self%id, self%str_id%to_char()
    end subroutine print_proc

    subroutine print_file(self)
        class(chunk_rec), intent(in) :: self
        print '("chunk_rec(id=",i0,", file=",a,")")', self%id, self%projfile%to_char()
    end subroutine print_file

    ! -------------------- Rec Iterator API --------------------

    logical function iter_valid(self)
        class(rec_iterator), intent(in) :: self
        iter_valid = self%it%has_value()
    end function iter_valid

    subroutine iter_next(self)
        class(rec_iterator), intent(inout) :: self
        call self%it%next()
    end subroutine iter_next

    subroutine iter_advance(self, n)
        class(rec_iterator), intent(inout) :: self
        integer,             intent(in)    :: n
        call self%it%advance(n)
    end subroutine iter_advance

    subroutine iter_replace(self, lst, r)
        class(rec_iterator), intent(inout) :: self
        class(rec),          intent(in)    :: r
        class(rec_list),     intent(inout) :: lst
        call lst%list%replace_iterator(self%it, r)
    end subroutine iter_replace

    subroutine iter_get_rec(self, x)
        class(rec_iterator),     intent(in)  :: self
        class(rec), allocatable, intent(out) :: x
        class(*),   allocatable :: tmp
        call self%it%getter(tmp)
        select type(tmp)
        class is(rec)
            allocate(x, source=tmp)
        class default
            THROW_HARD("ERROR: non-rec stored in rec_list")
        end select
        if (allocated(tmp)) deallocate(tmp)
    end subroutine iter_get_rec

    subroutine iter_get_project_rec(self, x)
        class(rec_iterator), intent(in)  :: self
        type(project_rec),   intent(out) :: x
        class(rec), allocatable :: tmp
        call self%iter_get_rec(tmp)
        select type(t => tmp)
        type is (project_rec)
            x = t
        class default
            THROW_HARD('not project_rec type')
        end select
        if (allocated(tmp)) deallocate(tmp)
    end subroutine iter_get_project_rec

    subroutine iter_get_process_rec(self, x)
        class(rec_iterator), intent(in)  :: self
        type(process_rec),   intent(out) :: x
        class(rec), allocatable :: tmp
        call self%iter_get_rec(tmp)
        select type(t => tmp)
        type is (process_rec)
            x = t
        class default
            THROW_HARD('not process_rec type')
        end select
        if (allocated(tmp)) deallocate(tmp)
    end subroutine iter_get_process_rec

    subroutine iter_get_chunk_rec(self, x)
        class(rec_iterator), intent(in)  :: self
        type(chunk_rec),     intent(out) :: x
        class(rec), allocatable :: tmp
        call self%iter_get_rec(tmp)
        select type(t => tmp)
        type is (chunk_rec)
            x = t
        class default
            THROW_HARD('not chunk_rec type')
        end select
        if (allocated(tmp)) deallocate(tmp)
    end subroutine iter_get_chunk_rec

    logical function iter_equals(self, other)
        class(rec_iterator), intent(in) :: self, other
        iter_equals = self%it%equals(other%it)
    end function iter_equals

    ! -------------------- rec_list API --------------------

    subroutine push_back(self, x)
        class(rec_list), intent(inout) :: self
        class(rec),      intent(in)    :: x
        call self%list%push_back(x)
    end subroutine push_back

    subroutine push2chunk_list( chunk_list, fname, id, processed )
        class(rec_list), intent(inout) :: chunk_list
        class(string),   intent(in)    :: fname
        integer,         intent(in)    :: id
        logical,         intent(in)    :: processed
        type(chunk_rec) :: crec
        crec%id         = id
        crec%projfile   = fname
        crec%busy       = .false.
        crec%processed  = processed
        crec%included   = .false.
        call chunk_list%push_back(crec)
    end subroutine push2chunk_list

    subroutine push2project_list( project_list, id, folder, projfile )
        class(rec_list), intent(inout) :: project_list
        class(string),   intent(in)    :: id, folder, projfile
        type(process_rec) :: record
        record%str_id    = id
        record%folder    = folder
        record%projfile  = projfile
        record%volume    = ''
        record%alnvolume = ''
        record%submitted = .false.
        record%completed = .false.
        record%included  = .false.
        call project_list%push_back(record)
    end subroutine push2project_list

    subroutine replace_at(self, idx, x)
        class(rec_list), intent(inout) :: self
        integer,         intent(in)    :: idx
        class(rec),      intent(in)    :: x
        call self%list%replace_at(idx, x)
    end subroutine replace_at

    subroutine replace_iterator(self, it, x)
        class(rec_list),    intent(inout) :: self
        type(rec_iterator), intent(inout) :: it
        class(rec),         intent(in)    :: x
        call self%list%replace_iterator(it%it, x)
    end subroutine replace_iterator

    pure integer function size(self)
        class(rec_list), intent(in) :: self
        size = self%list%size()
    end function size

    subroutine at_rec(self, idx, x)
        class(rec_list), intent(in) :: self
        integer,         intent(in) :: idx
        class(rec), allocatable, intent(out) :: x
        class(*), allocatable :: tmp
        call self%list%at(idx, tmp)
        select type(tmp)
        class is(rec)
            allocate(x, source=tmp)
        class default
            THROW_HARD("rec_list: non-rec type stored")
        end select
        if (allocated(tmp)) deallocate(tmp)
    end subroutine at_rec

    subroutine at_project_rec(self, idx, x)
        class(rec_list),   intent(in) :: self
        integer,           intent(in) :: idx
        type(project_rec), intent(out) :: x
        class(rec), allocatable :: tmp
        call self%at_rec(idx, tmp)
        select type(t => tmp)
        type is (project_rec)
            x = t
        class default
            THROW_HARD('not project_rec type')
        end select
        if (allocated(tmp)) deallocate(tmp)
    end subroutine at_project_rec

    subroutine at_process_rec(self, idx, x)
        class(rec_list),   intent(in) :: self
        integer,           intent(in) :: idx
        type(process_rec), intent(out) :: x
        class(rec), allocatable :: tmp
        call self%at_rec(idx, tmp)
        select type(t => tmp)
        type is (process_rec)
            x = t
        class default
            THROW_HARD('not process_rec type')
        end select
        if (allocated(tmp)) deallocate(tmp)
    end subroutine at_process_rec

    subroutine at_chunk_rec(self, idx, x)
        class(rec_list),   intent(in) :: self
        integer,           intent(in) :: idx
        type(chunk_rec), intent(out) :: x
        class(rec), allocatable :: tmp
        call self%at_rec(idx, tmp)
        select type(t => tmp)
        type is (chunk_rec)
            x = t
        class default
            THROW_HARD('not chunk_rec type')
        end select
        if (allocated(tmp)) deallocate(tmp)
    end subroutine at_chunk_rec

    ! ---------------- Iterators ----------------

    function begin(self) result(it)
        class(rec_list), intent(in) :: self
        type(rec_iterator) :: it
        it%it = self%list%begin()
    end function begin

    function end_iter(self) result(it)
        class(rec_list), intent(in) :: self
        type(rec_iterator) :: it
        it%it = self%list%end_iter()
    end function end_iter

    ! ---------------- Higher-level operations ----------------

    function append(self, other) result( res )
        class(rec_list), intent(in) :: self
        class(rec_list), intent(in) :: other
        type(rec_list) :: res
        res%list = self%list // other%list
    end function append

    subroutine subset_to(self, from, i1, i2)
        class(rec_list), intent(inout) :: self
        class(rec_list), intent(in)    :: from
        integer,         intent(in)    :: i1, i2
        type(linked_list) :: tmp
        call from%list%slice(i1, i2, tmp)
        call self%list%replace_with(tmp)   ! move semantics
    end subroutine subset_to

    subroutine slice(self, i1, i2, this)
        class(rec_list), intent(in)  :: self
        integer,         intent(in)  :: i1, i2
        type(rec_list),  intent(out) :: this
        call self%list%slice(i1, i2, this%list)
    end subroutine slice

    subroutine replace_with(self, other)
        class(rec_list), intent(inout) :: self
        class(rec_list), intent(inout) :: other
        call self%list%replace_with(other%list)
    end subroutine replace_with

    ! ---------------- Lifecycle ----------------

    subroutine assign(self, other)
        class(rec_list), intent(inout) :: self
        class(rec_list), intent(in)    :: other
        self%list = other%list
    end subroutine assign

    subroutine kill(self)
        class(rec_list), intent(inout) :: self
        call self%list%kill()
    end subroutine kill

    ! ---------------- Specialized methods ----------------
    
    function get_ids( self ) result( ids )
        class(rec_list), intent(in) :: self
        integer,    allocatable     :: ids(:)
        class(rec), allocatable     :: r
        type(rec_iterator)          :: it
        integer :: n, i
        n = self%size()
        if (n == 0) return
        allocate(ids(n), source=0)
        it = self%begin()
        i = 0
        do while (it%valid())
            i = i + 1
            call it%iter_get_rec(r)
            ids(i) = r%id
            call it%next()
        end do
    end function get_ids 

    !=====================================================================
    ! Return logical mask of "included" flags for all stored records
    !=====================================================================
    function get_included_flags(self) result(mask)
        class(rec_list), intent(in) :: self
        logical, allocatable        :: mask(:)
        type(rec_iterator)          :: it
        class(rec), allocatable     :: r
        integer :: n, i
        n = self%size()
        if (n == 0) return
        allocate(mask(n))
        it = self%begin()
        i = 0
        do while (it%valid())
            i = i + 1
            call it%iter_get_rec(r)
            select type(t=>r)
            type is(project_rec)
                mask(i) = t%included
            type is(process_rec)
                mask(i) = t%included
            type is(chunk_rec)
                mask(i) = t%included
            class default
                mask(i) = .false.
            end select
            call it%next()
        end do
    end function get_included_flags

    function get_processed_flags(self) result(mask)
        class(rec_list), intent(in) :: self
        logical, allocatable        :: mask(:)
        type(rec_iterator)          :: it
        class(rec), allocatable     :: r
        integer :: n, i
        n = self%size()
        if (n == 0) return
        allocate(mask(n))
        it = self%begin()
        i = 0
        do while (it%valid())
            i = i + 1
            call it%iter_get_rec(r)
            select type(t=>r)
            type is(chunk_rec)
                mask(i) = t%processed
            class default
                mask(i) = .false.
            end select
            call it%next()
        end do
    end function get_processed_flags

    function get_busy_flags(self) result(mask)
        class(rec_list), intent(in) :: self
        logical, allocatable        :: mask(:)
        type(rec_iterator)          :: it
        class(rec), allocatable     :: r
        integer :: n, i
        n = self%size()
        if (n == 0) return
        allocate(mask(n))
        it = self%begin()
        i = 0
        do while (it%valid())
            i = i + 1
            call it%iter_get_rec(r)
            select type(t=>r)
            type is(chunk_rec)
                mask(i) = t%busy
            class default
                mask(i) = .false.
            end select
            call it%next()
        end do
    end function get_busy_flags

    !=====================================================================
    ! Mark included=.true. for entire list or subset range
    ! Uses iterator replacement mechanics.
    !=====================================================================
    subroutine set_included_flags(self, fromto)
        class(rec_list),   intent(inout) :: self
        integer, optional, intent(in)    :: fromto(2)
        type(rec_iterator) :: it
        class(rec), allocatable :: r
        integer :: i, stop_idx
        if (self%size() == 0) return
        it = self%begin()
        i = 1
        stop_idx = self%size()
        if (present(fromto)) then
            call it%advance(fromto(1)-1)
            i        = fromto(1)
            stop_idx = fromto(2)
        end if
        do while (it%valid())
            call it%iter_get_rec(r)
            select type(t => r)
            type is(project_rec)
                t%included = .true.
                call it%replace(self, t)
            type is(process_rec)
                t%included = .true.
                call it%replace(self, t)
            type is(chunk_rec)
                t%included = .true.
                call it%replace(self, t)
            end select
            if (i == stop_idx) exit
            i = i + 1
            call it%next()
        end do
    end subroutine set_included_flags

    !=====================================================================
    ! Sum total particles, optionally only count NOT included ones
    !=====================================================================
    integer function get_nptcls_tot(self, l_not_included) result(total)
        class(rec_list),   intent(in) :: self
        logical, optional, intent(in) :: l_not_included
        type(rec_iterator) :: it
        class(rec), allocatable :: r
        logical :: not_include_flag
        total = 0
        if (self%size() == 0) return
        not_include_flag = .false.
        if (present(l_not_included)) not_include_flag = l_not_included
        it = self%begin()
        do while (it%valid())
            call it%iter_get_rec(r)
            select type(t => r)
            type is(project_rec)
                if (not_include_flag) then
                    if (.not. t%included) total = total + t%nptcls
                else
                    total = total + t%nptcls
                end if
            class default
                THROW_HARD('requires project_rec only')
            end select
            call it%next()
        end do
    end function get_nptcls_tot

    !=====================================================================
    ! Count selected particles
    !=====================================================================
    integer function get_nptcls_sel_tot(self) result(total)
        class(rec_list), intent(in) :: self
        type(rec_iterator) :: it
        class(rec), allocatable :: r
        total = 0
        if (self%size() == 0) return
        it = self%begin()
        do while (it%valid())
            call it%iter_get_rec(r)
            select type(t => r)
            type is(project_rec)
                total = total + t%nptcls_sel
            class default
                THROW_HARD('requires project_rec only')
            end select
            call it%next()
        end do
    end function get_nptcls_sel_tot

    function get_projfiles(self, fromto) result( projfiles )
        class(rec_list), intent(in) :: self
        integer,         intent(in) :: fromto(2)
        type(string),   allocatable :: projfiles(:)
        class(rec),     allocatable :: r
        type(rec_list)     :: tmplst
        type(rec_iterator) :: it
        integer            :: n, i
        call self%slice(fromto(1), fromto(2), tmplst)
        n  = tmplst%size()
        allocate(projfiles(n))
        it = tmplst%begin()
        i  = 0
        do while (it%valid())
            i = i + 1
            call it%iter_get_rec(r)
            select type(t => r)
            type is(chunk_rec)
                projfiles(i) = t%projfile
            class default
                THROW_HARD('requires chunk_rec only')
            end select
            call it%next()
        enddo
        call tmplst%kill
    end function get_projfiles

end module simple_rec_list
