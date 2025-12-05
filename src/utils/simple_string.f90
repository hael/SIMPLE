module simple_string
use simple_error, only: simple_exception
use simple_defs_string
use simple_defs
implicit none

public :: string
private

interface string
    module procedure char2str
    module procedure int2str
    module procedure real2str
    module procedure dble2str
end interface string

type :: string
    private
    character(len=:), allocatable :: buffer
  contains
    procedure, private :: assign
    generic            :: assignment(=) => assign
    procedure, private :: append
    generic            :: operator(//)  => append
    procedure          :: strlen
    procedure          :: strlen_trim
    procedure          :: raw
    procedure, private :: to_char_1, to_char_2
    generic            :: to_char => to_char_1, to_char_2
    procedure, private :: to_static_1, to_static_2
    generic            :: to_static => to_static_1, to_static_2
    procedure          :: to_int
    procedure          :: to_real
    procedure          :: to_dble
    procedure, private :: has_substr_1, has_substr_2
    generic            :: has_substr => has_substr_1, has_substr_2
    procedure          :: substr_remove
    procedure, private :: ends_with_substr_1, ends_with_substr_2
    generic            :: ends_with_substr => ends_with_substr_1, ends_with_substr_2
    procedure, private :: substr_ind_1, substr_ind_2
    generic            :: substr_ind => substr_ind_1, substr_ind_2
    procedure          :: substr_replace
    procedure          :: is_blank
    procedure          :: is_allocated
    procedure, private :: is_equal_1, is_equal_2
    generic            :: operator(.eq.) => is_equal_1, is_equal_2
    generic            :: operator(==) => is_equal_1, is_equal_2
    procedure, private :: is_not_equal_1, is_not_equal_2
    generic            :: operator(.ne.) => is_not_equal_1, is_not_equal_2
    generic            :: operator(/=) => is_not_equal_1, is_not_equal_2
    procedure          :: readline
    procedure          :: writeline
    procedure, private :: set_real, set_dble, set_int
    generic            :: set => set_real, set_dble, set_int
    procedure          :: kill
end type string

contains

    function char2str( rhs ) result( self )
        character(len=*), intent(in) :: rhs
        type(string) :: self
        call self%assign(rhs)
    end function char2str

    function int2str(rhs) result(self)
        integer, intent(in) :: rhs
        type(string)        :: self
        call self%set(rhs)
    end function int2str

    function real2str(rhs) result(self)
        real, intent(in) :: rhs
        type(string)     :: self
        call self%set(rhs)
    end function real2str

    function dble2str(rhs) result(self)
        real(dp), intent(in) :: rhs
        type(string)         :: self
        call self%set(rhs)
    end function dble2str

    subroutine assign(self, rhs)
        class(string), intent(inout) :: self
        class(*),      intent(in)    :: rhs
        character(len=:), allocatable :: tmp
        call self%kill
        select type (rhs)
            type is (string)
                tmp = rhs%buffer
            type is (character(*))
                tmp = rhs
            class default
                call simple_exception('Unsupported assignment', __FILENAME__ , __LINE__)
        end select
        if( allocated(tmp) ) call move_alloc(tmp, self%buffer)
    end subroutine assign

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

    elemental integer function strlen(self)
        class(string), intent(in) :: self
        strlen = 0
        if( allocated(self%buffer) ) strlen = len(self%buffer)
    end function strlen

    elemental integer function strlen_trim(self)
        class(string), intent(in) :: self
        strlen_trim = 0
        if( allocated(self%buffer) ) strlen_trim = len_trim(self%buffer)
    end function strlen_trim

    pure function raw(self) result(c)
        class(string), intent(in) :: self
        character(:), allocatable :: c
        if (allocated(self%buffer)) then
            allocate(character(len(self%buffer)) :: c)
            if (len(c) > 0) c = self%buffer
        else
            allocate(character(0) :: c)
        end if
    end function raw

    pure function to_char_1(self) result(c)
        class(string), intent(in) :: self
        character(:), allocatable :: c
        if (.not. allocated(self%buffer)) then
            c = ''
        else
            allocate(character(len_trim(self%buffer)) :: c)
            if (len(c) > 0) c = trim(self%buffer)
        end if
    end function to_char_1

    pure function to_char_2(self, fromto) result(c)
        class(string), intent(in) :: self
        integer,       intent(in) :: fromto(2)
        character(:), allocatable :: c
        integer :: i1, i2
        i1 = fromto(1)
        i2 = fromto(2)
        if (.not. allocated(self%buffer) .or. i1 < 1 .or. i2 < i1 .or. i2 > len(self%buffer)) then
            c = ''  ! return empty
        else
            allocate(character(len(self%buffer(i1:i2))) :: c)
            c = self%buffer(i1:i2)
        end if
    end function to_char_2

    ! this needs to be used in OpenMP sections
    pure subroutine to_static_1( self, str_static )
        class(string),    intent(in)  :: self
        character(len=*), intent(out) :: str_static
        if( allocated(self%buffer) )then
            str_static = self%buffer
        else
            str_static = ''
        endif
    end subroutine to_static_1

    ! this needs to be used in OpenMP sections
    pure subroutine to_static_2( self, fromto, str_static )
        class(string),    intent(in)  :: self
        integer,          intent(in)  :: fromto(2)
        character(len=*), intent(out) :: str_static
        integer :: i1, i2
        i1 = fromto(1)
        i2 = fromto(2)
        if (.not. allocated(self%buffer) .or. i1 < 1 .or. i2 < i1 .or. i2 > len(self%buffer)) then
            str_static = '' ! return empty
        else
            str_static = self%buffer(i1:i2)
        end if
    end subroutine to_static_2

    integer function to_int( self, io_stat )
        class(string),     intent(in)  :: self
        integer, optional, intent(out) :: io_stat 
        integer :: ios
        to_int = 0
        if( allocated(self%buffer) )then
            read(self%buffer,*,iostat=ios) to_int
        else
            ios = -1
        endif
        if( present(io_stat) ) io_stat = ios
    end function to_int

    real function to_real( self, io_stat )
        class(string),     intent(in)  :: self
        integer, optional, intent(out) :: io_stat 
        integer :: ios
        to_real = 0.
        if( allocated(self%buffer) )then
            read(self%buffer,*,iostat=ios) to_real
        else
            ios = -1
        endif
        if( present(io_stat) ) io_stat = ios
    end function to_real

    real(dp) function to_dble( self, io_stat )
        class(string),     intent(in)  :: self
        integer, optional, intent(out) :: io_stat 
        integer :: ios
        to_dble = 0.d0
        if( allocated(self%buffer) )then
            read(self%buffer,*,iostat=ios) to_dble
        else
            ios = -1
        endif
        if( present(io_stat) ) io_stat = ios
    end function to_dble

    elemental logical function has_substr_1( self, substr )
        class(string), intent(in) :: self, substr
        has_substr_1 = .false.
        if( allocated(self%buffer) .and. allocated(substr%buffer) ) has_substr_1 = .not. (index(self%buffer, substr%buffer) == 0)
    end function has_substr_1

    elemental logical function has_substr_2( self, substr )
        class(string),    intent(in) :: self
        character(len=*), intent(in) :: substr
        has_substr_2 = .false.
        if( allocated(self%buffer) ) has_substr_2 = .not. (index(self%buffer, substr) == 0)
    end function has_substr_2

    function substr_remove( self, substr ) result( self_sub )
        class(string), intent(in) :: self, substr
        type(string) :: self_sub
        if( allocated(self%buffer) .and. allocated(substr%buffer) )then
            call rm_substr(self%buffer, substr%buffer, self_sub%buffer)
        endif

        contains

            subroutine rm_substr( str, substr, str_out )
                character(len=*),              intent(in)    :: str, substr
                character(len=:), allocatable, intent(inout) :: str_out
                integer :: l, lout, pos
                str_out = str
                l = len(substr)
                if( l == 0 )return
                if( index(str_out, substr) == 0 )return
                pos = index(str_out, substr)
                do while( pos > 0 )
                    lout = len(str_out)
                    if( (pos==1) .and. (lout==l) )then
                        str_out = ''
                        return
                    else if( pos == 1)then
                        str_out = str_out(pos+l:lout)
                    else
                        if( pos+l-1 == lout )then
                            str_out = str_out(1:pos-1)
                        else
                            str_out = str_out(1:pos-1)//str_out(pos+l:lout)
                        endif
                    endif
                    pos = index(str_out, substr)
                end do
            end subroutine rm_substr

    end function substr_remove

    elemental logical function ends_with_substr_1( str, substr )
        class(string), intent(in) :: str, substr
        integer :: start_substr, lenstr
        ends_with_substr_1  = .false.
        if( allocated(str%buffer) .and. allocated(substr%buffer) )then
            lenstr              = len(str%buffer)
            start_substr        = lenstr - len(substr%buffer) + 1
            if( start_substr >= 1 )then
                ends_with_substr_1 = str%buffer(start_substr:lenstr) == substr%buffer
            endif
        endif
    end function ends_with_substr_1

    elemental logical function ends_with_substr_2( str, substr )
        class(string),    intent(in) :: str
        character(len=*), intent(in) :: substr
        integer :: start_substr, lenstr
        ends_with_substr_2  = .false.
        if( allocated(str%buffer) )then
            lenstr              = len(str%buffer)
            start_substr        = lenstr - len(substr) + 1
            if( start_substr >= 1 )then
                ends_with_substr_2 = str%buffer(start_substr:lenstr) == substr
            endif
        endif
    end function ends_with_substr_2

    elemental integer function substr_ind_1( self, substr, back )
        class(string),     intent(in) :: self
        character(len=*),  intent(in) :: substr
        logical, optional, intent(in) :: back 
        substr_ind_1 = 0
        if( allocated(self%buffer) )then
            if( present(back) )then
                substr_ind_1 = index(self%buffer, substr, back=back)
            else
                substr_ind_1 = index(self%buffer, substr)
            endif
        endif
    end function substr_ind_1

    integer function substr_ind_2( self, substr, back )
        class(string),     intent(in) :: self, substr
        logical, optional, intent(in) :: back 
        substr_ind_2 = 0
        if( allocated(self%buffer) .and. allocated(substr%buffer) )then
            if( present(back) )then
                substr_ind_2 = index(self%buffer, substr%buffer, back=back)
            else
                substr_ind_2 = index(self%buffer, substr%buffer)
            endif
        endif
    end function substr_ind_2

    pure subroutine substr_replace( self, s1, s2, one, back )
        class(string),     intent(inout) :: self
        character(len=*),  intent(in)    :: s1, s2
        logical, optional, intent(in)    :: one, back
        if( allocated(self%buffer) )then
            call rpl_substr(self%buffer, s1, s2, one, back)
        endif

      contains

        !>  Replace all occurrences of s1 in str with s2
        !!  Optionally one occurence (one=.true.) and
        !!  optionally starting from end of string (back=.true.)
        pure subroutine rpl_substr( str, s1, s2, one, back )
            character(len=:), allocatable, intent(inout) :: str
            character(len=*),              intent(in)    :: s1, s2
            logical,             optional, intent(in)    :: one, back
            character(len=:), allocatable :: tmp
            integer :: i, n, ilen, ilen1
            logical :: l_one, l_back
            l_one  = .false.
            l_back = .false.
            if( present(one) ) l_one = one
            if( l_one )then
                if( present(back) ) l_back = back
            endif
            if ( len(str) > 0 ) then
                tmp   = ''
                ilen1 = len(s1)
                do
                    ilen = len(str)
                    i    = index(str,s1,back=l_back)
                    if (i>0) then
                        if (i>1) tmp = tmp//str(1:i-1)
                        tmp = tmp//s2
                        n   = i+ilen1
                        if (n<=ilen) then
                            str = str(n:ilen)
                        else
                            exit
                        end if
                        if( l_one )then
                            tmp = tmp//str
                            exit
                        endif
                    else
                        tmp = tmp//str
                        exit
                    end if
                end do
                str = tmp
            end if
        end subroutine rpl_substr
    
    end subroutine substr_replace

    logical function is_blank( self )
        class(string), intent(in) :: self
        is_blank = .true.
        if( allocated(self%buffer) ) is_blank = len_trim(self%buffer) == 0        
    end function is_blank

    elemental logical function is_allocated( self )
        class(string), intent(in) :: self
        is_allocated = allocated(self%buffer)
    end function is_allocated

    ! invariant to leading and trailing spaces
    elemental logical function is_equal_1( self1, self2 )
        class(string), intent(in) :: self1, self2
        is_equal_1 = .false.
        if( allocated(self1%buffer) .and. allocated(self2%buffer) )then
            is_equal_1 = trim(adjustl(self1%buffer)) .eq. trim(adjustl(self2%buffer))
        endif
    end function is_equal_1

    ! invariant to leading and trailing spaces
    elemental logical function is_equal_2( self, str_char )
        class(string),    intent(in) :: self
        character(len=*), intent(in) :: str_char
        is_equal_2 = .false.
        if( allocated(self%buffer) )then
            is_equal_2 = trim(adjustl(self%buffer)) .eq. trim(adjustl(str_char))
        endif
    end function is_equal_2

    ! invariant to leading and trailing spaces
    elemental logical function is_not_equal_1( self1, self2 )
        class(string), intent(in) :: self1, self2
        is_not_equal_1 = .true.
        if( allocated(self1%buffer) .and. allocated(self2%buffer) )then
            is_not_equal_1 = trim(adjustl(self1%buffer)) .ne. trim(adjustl(self2%buffer))
        endif
    end function is_not_equal_1

    ! invariant to leading and trailing spaces
    elemental logical function is_not_equal_2( self, str_char )
        class(string),    intent(in) :: self
        character(len=*), intent(in) :: str_char
        is_not_equal_2 = .true.
        if( allocated(self%buffer) )then
            is_not_equal_2 = trim(adjustl(self%buffer)) .ne. trim(adjustl(str_char))
        endif
    end function is_not_equal_2

    elemental subroutine set_real(self, val)
        class(string), intent(inout) :: self
        real,          intent(in)    :: val
        character(len=STDLEN) :: tmp ! temporary buffer for WRITE
        call self%kill
        write(tmp, '(G0)') val       ! G0 = automatic-width formatting
        self%buffer = trim(tmp)      ! assign trimmed text to allocatable string
    end subroutine set_real

    elemental subroutine set_dble(self, val)
        class(string), intent(inout) :: self
        real(dp),      intent(in)    :: val
        character(len=STDLEN) :: tmp ! temporary buffer for WRITE
        call self%kill
        write(tmp, '(G0)') val       ! G0 = automatic-width formatting
        self%buffer = trim(tmp)      ! assign trimmed text to allocatable string
    end subroutine set_dble

    elemental subroutine set_int(self, val)
        class(string), intent(inout) :: self
        integer,      intent(in)     :: val
        character(len=STDLEN) :: tmp ! temporary buffer for WRITE
        call self%kill
        write(tmp, '(G0)') val       ! G0 = automatic-width formatting
        self%buffer = trim(tmp)      ! assign trimmed text to allocatable string
    end subroutine set_int
    
    ! I/O

    subroutine readline( self, fhandle, ios )
        class(string), intent(inout) :: self
        integer,       intent(in)    :: fhandle
        integer,       intent(out)   :: ios
        character(len=XLONGSTRLEN)   :: buffer_static 
        call self%kill
        read(fhandle, '(A)', iostat=ios) buffer_static
        if( ios == 0 ) self%buffer = trim(adjustl(buffer_static))
    end subroutine readline

    subroutine writeline( self, fhandle, ios )
        class(string), intent(inout) :: self
        integer,       intent(in)    :: fhandle
        integer,       intent(out)   :: ios
        if( allocated(self%buffer) )then
            write(fhandle, '(A)', iostat=ios) trim(adjustl(self%buffer))
        else
            ios = -1
        endif
    end subroutine writeline
        
    elemental subroutine kill( self )
        class(string), intent(inout) :: self
        if( allocated(self%buffer) ) deallocate(self%buffer)
    end subroutine kill

end module simple_string
