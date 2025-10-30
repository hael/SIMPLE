module simple_chatgpt_string_class
use iso_fortran_env, only: dp => real64
implicit none

public :: string
private

interface string
    module procedure :: string_from_char_, string_from_int_, string_from_real_
end interface string

type :: string
    private
    character(:), allocatable :: data
contains
    ! Core
    procedure, private :: string_assign
    procedure, private :: string_concat
    procedure, private :: string_equal
    generic   :: assignment(=)=> string_assign
    generic   :: operator(//) => string_concat
    generic   :: operator(==) => string_equal
    procedure :: append       => string_append
    ! Transformations
    procedure :: trimmed      => string_trimmed
    procedure :: strip        => string_strip
    procedure :: lstrip       => string_lstrip
    procedure :: rstrip       => string_rstrip
    procedure :: to_lower     => string_to_lower
    procedure :: to_upper     => string_to_upper
    procedure :: replace      => string_replace
    procedure :: split        => string_split
    procedure :: join         => string_join
    procedure :: contains     => string_contains
    procedure :: is_empty     => string_is_empty

    ! Numeric conversions
    procedure :: to_int       => string_to_int
    procedure :: to_real      => string_to_real
    procedure :: from_int     => string_from_int
    procedure :: from_real    => string_from_real

    ! Utility
    procedure :: pad_left     => string_pad_left
    procedure :: pad_right    => string_pad_right
end type string

contains

!===============================================================
! Constructors
!===============================================================
    pure function string_from_char_(s) result(str)
        character(len=*), intent(in) :: s
        type(string) :: str
        str%data = s
    end function string_from_char_

    pure function string_from_int_(i) result(str)
        integer, intent(in) :: i
        type(string) :: str
        character(len=64) :: tmp
        write(tmp, '(I0)') i
        str%data = trim(tmp)
    end function string_from_int_

    pure function string_from_real_(x) result(str)
        real(dp), intent(in) :: x
        type(string) :: str
        character(len=64) :: tmp
        write(tmp, '(G0)') x
        str%data = trim(tmp)
    end function string_from_real_

!===============================================================
! Basic operations
!===============================================================
    elemental subroutine string_assign(lhs, rhs)
        class(string), intent(inout) :: lhs
        character(len=*), intent(in) :: rhs
        lhs%data = rhs
    end subroutine string_assign

    function string_concat(a, b) result(c)
        class(string), intent(in) :: a, b
        type(string) :: c
        c%data = a%data // b%data
    end function string_concat

    pure logical function string_equal(a, b)
        class(string), intent(in) :: a, b
        string_equal = trim(a%data) == trim(b%data)
    end function string_equal

    subroutine string_append(this, txt)
        class(string), intent(inout) :: this
        character(len=*), intent(in) :: txt
        if (.not. allocated(this%data)) then
            this%data = txt
        else
            this%data = this%data // txt
        end if
    end subroutine string_append

!===============================================================
! Whitespace and trimming
!===============================================================
    pure function string_trimmed(this) result(out)
        class(string), intent(in) :: this
        type(string) :: out
        out%data = trim(adjustl(this%data))
    end function string_trimmed

    pure function string_strip(this) result(out)
        class(string), intent(in) :: this
        type(string) :: out
        character(len=:), allocatable :: tmp
        integer :: i, j, n
        tmp = this%data
        n = len(tmp)
        i = 1
        do while (i <= n .and. tmp(i:i) == ' ')
            i = i + 1
        end do
        j = n
        do while (j >= i .and. tmp(j:j) == ' ')
            j = j - 1
        end do
        if (i <= j) then
            out%data = tmp(i:j)
        else
            out%data = ''
        end if
    end function string_strip

    pure function string_lstrip(this) result(out)
        class(string), intent(in) :: this
        type(string) :: out
        integer :: i
        i = 1
        do while (i <= len(this%data) .and. this%data(i:i) == ' ')
            i = i + 1
        end do
        out%data = this%data(i:)
    end function string_lstrip

    pure function string_rstrip(this) result(out)
        class(string), intent(in) :: this
        type(string) :: out
        integer :: j
        j = len(this%data)
        do while (j >= 1 .and. this%data(j:j) == ' ')
            j = j - 1
        end do
        if (j > 0) then
            out%data = this%data(:j)
        else
            out%data = ''
        end if
    end function string_rstrip

!===============================================================
! Transformations
!===============================================================
    pure function string_to_lower(this) result(out)
        class(string), intent(in) :: this
        type(string) :: out
        integer :: i, n
        character(len=:), allocatable :: tmp
        tmp = this%data
        do i = 1, len(tmp)
            n = index('ABCDEFGHIJKLMNOPQRSTUVWXYZ', tmp(i:i))
            if (n > 0) tmp(i:i) = 'abcdefghijklmnopqrstuvwxyz'(n:n)
        end do
        out%data = tmp
    end function string_to_lower

    pure function string_to_upper(this) result(out)
        class(string), intent(in) :: this
        type(string) :: out
        integer :: i, n
        character(len=:), allocatable :: tmp
        tmp = this%data
        do i = 1, len(tmp)
            n = index('abcdefghijklmnopqrstuvwxyz', tmp(i:i))
            if (n > 0) tmp(i:i) = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'(n:n)
        end do
        out%data = tmp
    end function string_to_upper

    function string_replace(this, target, repl) result(out)
        class(string), intent(in) :: this
        character(len=*), intent(in) :: target, repl
        type(string) :: out
        integer :: pos
        character(len=:), allocatable :: s
        s = this%data
        pos = index(s, target)
        do while (pos > 0)
            s = s(1:pos-1) // repl // s(pos+len(target):)
            pos = index(s, target)
        end do
        out%data = s
    end function string_replace

!===============================================================
! Search / splitting
!===============================================================
    pure logical function string_contains(this, substr)
        class(string), intent(in) :: this
        character(len=*), intent(in) :: substr
        string_contains = index(this%data, substr) > 0
    end function string_contains

    function string_split(this, delim) result(parts)
        class(string), intent(in) :: this
        character(len=*), intent(in) :: delim
        type(string), allocatable :: parts(:)
        integer :: i, start, pos, n
        character(len=:), allocatable :: s
        s = this%data
        n = 0
        start = 1
        do
            pos = index(s(start:), delim)
            if (pos == 0) exit
            n = n + 1
            call add_part(parts, s(start:start+pos-2))
            start = start + pos + len(delim) - 1
        end do
        if (start <= len(s)) call add_part(parts, s(start:))
    contains
        subroutine add_part(arr, chunk)
            type(string), allocatable, intent(inout) :: arr(:)
            character(len=*), intent(in) :: chunk
            integer :: m
            type(string), allocatable :: tmp(:)
            if (.not. allocated(arr)) then
                allocate(arr(1))
                arr(1)%data = chunk
            else
                m = size(arr)
                allocate(tmp(m))
                tmp = arr
                deallocate(arr)
                allocate(arr(m+1))
                arr(1:m) = tmp
                arr(m+1)%data = chunk
            end if
        end subroutine add_part
    end function string_split

    function string_join(this, array, delim) result(out)
        class(string), intent(in) :: this
        type(string), intent(in) :: array(:)
        character(len=*), intent(in), optional :: delim
        type(string) :: out
        integer :: i
        character(len=*), parameter :: default_delim = ' '
        character(len=:), allocatable :: d
        d = merge(delim, default_delim, present(delim))
        out%data = ''
        do i = 1, size(array)
            out%data = out%data // array(i)%data
            if (i < size(array)) out%data = out%data // d
        end do
    end function string_join

!===============================================================
! Numeric conversion
!===============================================================
    function string_to_int(this, istat) result(val)
        class(string), intent(in) :: this
        integer, intent(out), optional :: istat
        integer :: val, ios
        read(this%data, *, iostat=ios) val
        if (present(istat)) istat = ios
        if (ios /= 0) val = 0
    end function string_to_int

    function string_to_real(this, istat) result(val)
        class(string), intent(in) :: this
        integer, intent(out), optional :: istat
        real(dp) :: val
        integer :: ios
        read(this%data, *, iostat=ios) val
        if (present(istat)) istat = ios
        if (ios /= 0) val = 0.0_dp
    end function string_to_real

    pure function string_from_int(this, i) result(out)
        class(string), intent(in) :: this
        integer, intent(in) :: i
        type(string) :: out
        character(len=64) :: tmp
        write(tmp, '(I0)') i
        out%data = trim(tmp)
    end function string_from_int

    pure function string_from_real(this, x) result(out)
        class(string), intent(in) :: this
        real(dp), intent(in) :: x
        type(string) :: out
        character(len=64) :: tmp
        write(tmp, '(G0)') x
        out%data = trim(tmp)
    end function string_from_real

!===============================================================
! Utility
!===============================================================
    pure logical function string_is_empty(this)
        class(string), intent(in) :: this
        string_is_empty = (.not. allocated(this%data)) .or. (len_trim(this%data) == 0)
    end function string_is_empty

    function string_pad_left(this, width) result(out)
        class(string), intent(in) :: this
        integer, intent(in) :: width
        type(string) :: out
        integer :: n
        n = len_trim(this%data)
        if (n >= width) then
            out%data = this%data
        else
            out%data = repeat(' ', width - n) // trim(this%data)
        end if
    end function string_pad_left

    function string_pad_right(this, width) result(out)
        class(string), intent(in) :: this
        integer, intent(in) :: width
        type(string) :: out
        integer :: n
        n = len_trim(this%data)
        if (n >= width) then
            out%data = this%data
        else
            out%data = trim(this%data) // repeat(' ', width - n)
        end if
    end function string_pad_right

    !===============================================================
! File I/O helpers
!===============================================================
    function read_file(filename, iostat) result(out)
        character(len=*), intent(in) :: filename
        integer, intent(out), optional :: iostat
        type(string) :: out
        character(len=512) :: line
        integer :: unit, ios
        character(len=:), allocatable :: buffer

        open(newunit=unit, file=filename, status='old', action='read', iostat=ios)
        if (ios /= 0) then
            if (present(iostat)) iostat = ios
            out%data = ''
            return
        end if

        buffer = ''
        do
            read(unit, '(A)', iostat=ios) line
            if (ios /= 0) exit
            buffer = buffer // trim(line) // new_line('a')
        end do
        close(unit)

        out%data = buffer
        if (present(iostat)) iostat = 0
    end function read_file

    function read_lines(filename, iostat) result(lines)
        character(len=*), intent(in) :: filename
        integer, intent(out), optional :: iostat
        type(string), allocatable :: lines(:)
        character(len=512) :: line
        integer :: unit, ios

        open(newunit=unit, file=filename, status='old', action='read', iostat=ios)
        if (ios /= 0) then
            if (present(iostat)) iostat = ios
            allocate(lines(0))
            return
        end if

        allocate(lines(0))
        do
            read(unit, '(A)', iostat=ios) line
            if (ios /= 0) exit
            call append_line(lines, line)
        end do
        close(unit)

        if (present(iostat)) iostat = 0
    contains
        subroutine append_line(arr, text)
            type(string), allocatable, intent(inout) :: arr(:)
            character(len=*), intent(in) :: text
            integer :: n
            type(string), allocatable :: tmp(:)
            if (.not. allocated(arr)) then
                allocate(arr(1))
                arr(1)%data = text
            else
                n = size(arr)
                allocate(tmp(n))
                tmp = arr
                deallocate(arr)
                allocate(arr(n+1))
                arr(1:n) = tmp
                arr(n+1)%data = text
            end if
        end subroutine append_line
    end function read_lines

    subroutine write_file(filename, this, append_mode, iostat)
        character(len=*), intent(in) :: filename
        class(string), intent(in) :: this
        logical, intent(in), optional :: append_mode
        integer, intent(out), optional :: iostat
        integer :: unit, ios
        character(len=*), parameter :: form = 'formatted'
        character(len=:), allocatable :: status
        logical :: do_append

        do_append = .false.
        if (present(append_mode)) do_append = append_mode
        status = merge('unknown', 'replace', do_append)

        open(newunit=unit, file=filename, status=status, action='write', position= &
             merge('append', 'rewind', do_append), form=form, iostat=ios)
        if (ios /= 0) then
            if (present(iostat)) iostat = ios
            return
        end if

        write(unit, '(A)') trim(this%data)
        close(unit)
        if (present(iostat)) iostat = 0
    end subroutine write_file

    subroutine write_lines(filename, arr, append_mode, iostat)
        character(len=*), intent(in) :: filename
        type(string), intent(in) :: arr(:)
        logical, intent(in), optional :: append_mode
        integer, intent(out), optional :: iostat
        integer :: unit, ios, i
        character(len=*), parameter :: form = 'formatted'
        character(len=:), allocatable :: status
        logical :: do_append

        do_append = .false.
        if (present(append_mode)) do_append = append_mode
        status = merge('unknown', 'replace', do_append)

        open(newunit=unit, file=filename, status=status, action='write', position= &
             merge('append', 'rewind', do_append), form=form, iostat=ios)
        if (ios /= 0) then
            if (present(iostat)) iostat = ios
            return
        end if

        do i = 1, size(arr)
            write(unit, '(A)') trim(arr(i)%data)
        end do
        close(unit)
        if (present(iostat)) iostat = 0
    end subroutine write_lines

end module simple_chatgpt_string_class
