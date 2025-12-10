! contains various string manipulation subroutines and functions
! adapted from George Benthien's string module http://gbenthien.net/strings/str-index.html
! modification and additions by Cyril Reboul, Michael Eager & Hans Elmlund
module simple_string_utils
use simple_error,  only: simple_exception
use simple_string, only: string
use simple_defs
use simple_defs_string
use, intrinsic :: iso_c_binding
implicit none

interface str2format
    module procedure str2format_1
    module procedure str2format_2
end interface str2format

interface str2int
    module procedure str2int_1
    module procedure str2int_2
    module procedure str2int_3
    module procedure str2int_4
end interface str2int

interface str2real
    module procedure str2real_1
    module procedure str2real_2
end interface str2real

interface real2str
    module procedure realsp2str
    module procedure realdp2str
end interface real2str

interface split
    module procedure split_1
    module procedure split_2
end interface split

interface findloc_str
    module procedure findloc_str_1
    module procedure findloc_str_2
end interface findloc_str

contains

    pure function spaces( n ) result( str )
        integer, intent(in) :: n
        character(len=:), allocatable :: str
        if( n <= 0 ) return
        allocate(character(len=n) :: str)
        str(:) = ' '
    end function spaces

    function cast_str_types( str ) result( str_tmp )
        class(*), intent(in) :: str
        character(len=:), allocatable :: str_tmp
        select type(str)
            type is (string)
                str_tmp = str%to_char()
            type is (character(*))
                str_tmp = trim(str)
            class default
                call simple_exception('Unsupported type', __FILENAME__ , __LINE__)
        end select
    end function cast_str_types

    !> \brief  assesses whether a string represents a filename
    subroutine str2format_1( str, format, rvar, ivar )
        character(len=*),              intent(in)  :: str
        character(len=:), allocatable, intent(out) :: format
        real,                          intent(out) :: rvar
        integer,                       intent(out) :: ivar
        real(dp) :: rvar_dp
        call str2format_2( str, format, rvar_dp, ivar )
        rvar = real(rvar_dp)
    end subroutine str2format_1

    !> \brief  assesses whether a string represents a filename
    subroutine str2format_2( str, format, rvar, ivar )
        character(len=*),              intent(in)  :: str
        character(len=:), allocatable, intent(out) :: format
        real(kind=dp),                 intent(out) :: rvar
        integer,                       intent(out) :: ivar
        integer :: iostat, i
        logical :: str_is_file
        i = index(str, '.')
        ! check if first symbol after . is a character
        str_is_file = .false.
        if( i /= 0 .and. i < len_trim(str)  )then
            if( char_is_a_letter(str(i+1:i+1)) )str_is_file = .true.
        endif
        ! file name
        if( str_is_file )then
            allocate( format, source='file' )
            return
        endif
        ! check if string contains / -> directory
        i = index(str, PATH_SEPARATOR)
        if( i /= 0 )then
            allocate(format, source='dir')
            return
        endif
        i = index(str,',')
        if( i /= 0 )then
            allocate(format,source='char')
            return
        endif
        ! real
        read(str,*,iostat=iostat) rvar
        if( iostat == 0 )then
            allocate( format, source='real' )
            return
        endif
        ! integer
        read(str,*, iostat=iostat) ivar
        if( iostat == 0 )then
            allocate( format, source='int' )
            return
        endif
        ! char is last resort
        allocate( format, source='char' )
    end subroutine str2format_2

    logical function char_is_a_letter( c )
        character(len=1), intent(in) :: c
        char_is_a_letter = (scan(LOWER_CASE_LETTERS,c)>0) .or. (scan(UPPER_CASE_LETTERS,c)>0)
    end function char_is_a_letter

    logical function char_is_a_number( c )
        character(len=1), intent(in) :: c
        char_is_a_number = scan(INTEGERS,c) > 0
    end function char_is_a_number

    !> \brief  parses the string 'str' into arguments args(1), ..., args(nargs) based on
    !!         the delimiters contained in the string 'delim'. The integer output variable
    !!         nargs contains the number of arguments found.
    subroutine parsestr( str, delim, args, nargs )
        character(len=*), intent(inout) :: str
        character(len=*), intent(in)    :: delim
        character(len=*)                :: args(:)
        integer,          intent(out)   :: nargs
        character(len=len_trim(str))    :: strsav
        integer                         :: na, i, lenstr
        strsav=str
        call compact(str)
        na=size(args)
        do i=1,na
            args(i)=''
        end do
        nargs=0
        lenstr=len_trim(str)
        if(lenstr==0) return
        do
            if(len_trim(str) == 0) exit
            nargs=nargs+1
            call split(str, delim, args(nargs))
        end do
        str=strsav
    end subroutine parsestr

    !> \brief  converts multiple spaces and tabs to single spaces;
    !!         deletes control characters; removes initial spaces.
    subroutine compact( str )
        character(len=*), intent(inout) :: str
        character(len=1):: ch
        character(len=len_trim(str)):: outstr
        integer :: lenstr, isp, k, i, ich
        str=adjustl(str)
        lenstr=len_trim(str)
        outstr=''
        isp=0
        k=0
        do i=1,lenstr
            ch=str(i:i)
            ich=iachar(ch)
            select case(ich)
                case(9,32) ! space or tab character
                    if(isp==0) then
                        k=k+1
                        outstr(k:k)=''
                    end if
                    isp=1
                case(33:) ! not a space, quote, or control character
                    k=k+1
                    outstr(k:k)=ch
                    isp=0
            end select
        end do
        str=adjustl(outstr)
    end subroutine compact

    !> \brief  finds the first instance of a character 'delim' in the
    !!         the string 'str'. The characters before the found delimiter are
    !!         output in 'before'. The characters after the found delimiter are
    !!         output in 'str'.
    subroutine split_str( str, delim, before )
        character(len=*), intent(inout) :: str
        character(len=1), intent(in)    :: delim
        character(len=*), intent(out)   :: before
        integer :: index
        str    = trim(str)
        index  = scan(str, delim)
        before = adjustl(trim(str(1:index-1)))
        str    = adjustl(str(index+1:))
    end subroutine split_str

    !> \brief  finds the first instance of a character 'delim' in the
    !!         the string 'str'. The characters before the found delimiter are
    !!         output in 'before'. The characters after the found delimiter are
    !!         output in 'str'.
    subroutine split_1( str, delim, before )
        character(len=*), intent(inout) :: str
        character(len=*), intent(in)    :: delim
        character(len=*), intent(out)   :: before
        character :: ch,cha
        integer   :: lenstr, k, iposa, i, ipos
        str=adjustl(str)
        call compact(str)
        lenstr=len_trim(str)
        if(lenstr == 0) return ! string str is empty
        k=0
        before=''
        do i=1,lenstr
            ch=str(i:i)
            ipos=index(delim,ch)
            if(ipos == 0)then ! character is not a delimiter
                k=k+1
                before(k:k)=ch
                cycle
            end if
            if(ch /= ' ')then ! character is a delimiter that is not a space
                str=str(i+1:)
                exit
            end if
            cha=str(i+1:i+1)  ! character is a space delimiter
            iposa=index(delim,cha)
            if(iposa > 0)then ! next character is a delimiter
                str=str(i+2:)
                exit
            else
                str=str(i+1:)
                exit
            end if
        end do
        if(i >= lenstr) str=''
        str=adjustl(str)      ! remove initial spaces
        return
    end subroutine split_1

    subroutine split_2( str_, delim, before_ )
        class(string),    intent(inout) :: str_
        character(len=*), intent(in)    :: delim
        class(string),    intent(out)   :: before_
        character             :: ch,cha
        integer               :: lenstr, k, iposa, i, ipos
        character(len=STDLEN) :: str, before
        str = str_%to_char()
        str=adjustl(str)
        call compact(str)
        lenstr=len_trim(str)
        if(lenstr == 0) return ! string str is empty
        k=0
        before=''
        do i=1,lenstr
            ch=str(i:i)
            ipos=index(delim,ch)
            if(ipos == 0)then ! character is not a delimiter
                k=k+1
                before(k:k)=ch
                cycle
            end if
            if(ch /= ' ')then ! character is a delimiter that is not a space
                str=str(i+1:)
                exit
            end if
            cha=str(i+1:i+1)  ! character is a space delimiter
            iposa=index(delim,cha)
            if(iposa > 0)then ! next character is a delimiter
                str=str(i+2:)
                exit
            else
                str=str(i+1:)
                exit
            end if
        end do
        if(i >= lenstr) str=''
        str     = adjustl(str)      ! remove initial spaces
        str_    = trim(adjustl(str))
        before_ = trim(adjustl(before))
        return
    end subroutine split_2

    function list_of_ints2arr( listofints ) result( iarr )
        character(len=*), intent(in) :: listofints
        integer, allocatable :: iarr(:)
        character(len=:), allocatable :: str, before
        integer :: index, cnt
        str   = adjustl(trim(listofints))
        index = scan(str, ',')
        if( index == 0 )then
            allocate(iarr(1), source=0)
            iarr(1) = str2int(str)
            return
        endif
        ! first, count commas
        cnt = 0
        do
            index  = scan(str, ',')
            if( index == 0 ) exit
            cnt    = cnt + 1
            before = adjustl(trim(str(1:index-1)))
            str    = adjustl(str(index+1:))
        end do
        ! allocate array
        allocate(iarr(cnt + 1), source=0)
        str = adjustl(trim(listofints))
        cnt = 0
        ! parse integers
        do
            cnt = cnt + 1
            index = scan(str, ',')
            if( index == 0 )then
                str       = adjustl(str(index+1:))
                iarr(cnt) = str2int(str)
                exit
            else
                before    = adjustl(trim(str(1:index-1)))
                iarr(cnt) = str2int(before)
                str       = adjustl(str(index+1:))
            endif    
        end do
    end function list_of_ints2arr

    !> \brief  removes punctuation (except comma) characters in string str
    subroutine removepunct( str )
        character(len=*), intent(inout) :: str
        character(len=1):: ch
        character(len=len_trim(str))::outstr
        integer :: lenstr, k, i, ich
        str = adjustl(str)
        lenstr = len_trim(str)
        outstr = ''
        k = 0
        do i=1,lenstr
            ch = str(i:i)
            ich = iachar(ch)
            select case(ich)
                case(32:43)   ! exclamation double-quote hash dollar percent ampersand quote ( ) * +  characters
                    cycle
                case(45:47)   ! . - fwd-slash characters
                    cycle
                case(58:64)   ! : ; < = > question-mark at characters
                    cycle
                case(91:94)   !  _ ^ [ ] characters
                    cycle
                case(123:126) ! { | }
                    cycle
                case default
                    k = k + 1
                    outstr(k:k) = ch
            end select
        end do
        str = adjustl(outstr)
    end subroutine removepunct

    !> \brief  removes spaces, tabs, and control characters in string str
    pure subroutine removesp( str )
        character(len=*), intent(inout) :: str
        character(len=1):: ch
        character(len=len_trim(str))::outstr
        integer :: lenstr, k, i, ich
        str = adjustl(str)
        lenstr = len_trim(str)
        outstr = ''
        k = 0
        do i=1,lenstr
            ch = str(i:i)
            ich = iachar(ch)
            select case(ich)
                case(0:32)  ! space, tab, or control character
                    cycle
                case(33:)
                    k = k + 1
                    outstr(k:k) = ch
            end select
        end do
        str=adjustl(outstr)
    end subroutine removesp

    function str2real_1( str ) result( rval )
        character(len=*), intent(in) :: str
        type(string) :: str_cast
        real         :: rval
        str_cast = trim(adjustl(str))
        rval = str_cast%to_real()
        call str_cast%kill
    end function str2real_1

    function str2real_2( str ) result( rval )
        class(string), intent(in) :: str
        real :: rval
        rval = str%to_real()
    end function str2real_2

    !> \brief  converts a real number to a string
    function realsp2str(rval) result(str)
        real, intent(in)  :: rval
        character(len=32) :: str
        type(string)      :: str_cast
        str_cast = string(rval)
        str = str_cast%to_char()
        call str_cast%kill
    end function realsp2str

    !> \brief  converts a real number to a string
    function realdp2str( rval ) result( str )
        real(kind=dp), intent(in) :: rval
        character(len=32)         :: str
        type(string)              :: str_cast
        str_cast = string(rval)
        str = str_cast%to_char()
        call str_cast%kill
    end function realdp2str

    !>  \brief  turn integer variable into character variable
    function int2str( intg ) result( str )
        integer, intent(in) :: intg
        character(len=:), allocatable :: str
        type(string) :: str_cast
        str_cast = string(intg)
        str = str_cast%to_char()
        call str_cast%kill
    end function int2str

    !>  \brief  turn integer variable into zero padded character variable
    function int2str_pad( intg, numlen ) result( str )
        integer, intent(in)           :: intg, numlen
        character(len=:), allocatable :: str, str_tmp, str_tmp2
        integer :: slen
        if( numlen == 0 )then
            str = int2str(intg)
            return
        endif
        str_tmp = int2str(intg)
        slen    = len(str_tmp)
        if( slen >= numlen )then
            allocate(str, source=str_tmp)
        else
            do while( len(str_tmp) < numlen )
                allocate(str_tmp2, source='0'//str_tmp)
                if( allocated(str) ) deallocate(str)
                allocate(str, source=str_tmp2)
                deallocate(str_tmp)
                allocate(str_tmp, source=str_tmp2)
                deallocate(str_tmp2)
            end do
        endif
        deallocate(str_tmp)
    end function int2str_pad

    integer function str2int_1( str )
        character(len=*), intent(in) :: str
        type(string) :: str_cast
        str_cast  = trim(adjustl(str))
        str2int_1 = str_cast%to_int()
        call str_cast%kill
    end function str2int_1

    integer function str2int_2( str )
        class(string), intent(in) :: str
        str2int_2 = str%to_int()
    end function str2int_2

    integer function str2int_3( str, io_stat )
        class(string), intent(in)  :: str
        integer,       intent(out) :: io_stat
        str2int_3 = str%to_int(io_stat)
    end function str2int_3

    integer function str2int_4( str, io_stat )
        character(len=*), intent(in)  :: str
        integer,          intent(out) :: io_stat
        type(string) :: str_cast
        str_cast  = trim(adjustl(str))
        str2int_4 = str_cast%to_int(io_stat)
        call str_cast%kill
    end function str2int_4

    !>  \brief  pad string with spaces
    function str_pad( str_in, len_padded ) result( str_out )
        character(len=*), intent(in)  :: str_in
        integer,          intent(in)  :: len_padded
        character(len=:), allocatable :: str_out
        integer :: slen
        slen = len(str_in)
        if( slen >= len_padded )then
            str_out = str_in
        else
            str_out = str_in // spaces(len_padded - slen)
        endif
    end function str_pad

    !>  \brief  maps out the number characters in a string
    !!          .true. means is number, .false. means is not
    function map_str_nrs( str ) result( map )
        character(len=*), intent(in) :: str
        logical, allocatable :: map(:)
        integer :: lstr, i, j, ind
        lstr = len(str)
        allocate(map(lstr))
        map = .false.
        do j=1,lstr
            do i=0,9
                ind = index(str(j:j), int2str(i))
                if( ind == 0 )then
                    cycle
                else
                    map(j) = .true.
                    exit
                endif
            end do
         end do
    end function map_str_nrs

    !>  \brief  check for presence of substring
    logical function str_has_substr( str, substr )
        character(len=*), intent(in) :: str, substr
        str_has_substr = .not. (index(str, substr) == 0)
    end function str_has_substr

    !>  \brief  Convert string to lower case
    pure function lowercase( str )
        character(len=*), intent(in) :: str
        character(len=(len(str)))    :: lowercase
        integer :: i, n
        lowercase = str
        do i = 1, len(str)
            n = index(UPPER_CASE_LETTERS,lowercase(i:i))
            if( n /= 0 ) lowercase(i:i) = LOWER_CASE_LETTERS(n:n)
        end do
    end function lowercase

    !>  \brief  Convert string to upper case
    pure function uppercase( str )
        character(len=*), intent(in)  ::  str
        character(len=len(str))       :: uppercase
        integer :: i, n
        uppercase = str
        do i=1,len(str)
            n = index(LOWER_CASE_LETTERS,uppercase(i:i))
            if( n /= 0 ) uppercase(i:i) = UPPER_CASE_LETTERS(n:n)
        enddo
    end function uppercase

    !>  \brief works out whether a character string is a comment line
    logical function str_is_comment( line )
        character(len=*), intent(in)  ::  line
        integer ::  pos1,linelen
        str_is_comment = .false.
        if(str_is_blank(line)) return
        pos1 = first_non_blank(line)
        ! if we didn't find a non-blank character, then the line must be blank!
        if( pos1 .eq. 0 ) return
        linelen=len_trim(line)
        if( scan(line(pos1:),'#!;', back=.false.) .eq. 1 )then
            ! the first non-blank character is a comment character.
            str_is_comment = .true.
        endif
    end function str_is_comment

    !>  \brief works out whether a character string is blank
    pure logical function str_is_blank( line )
        character(len=*),   intent(in)  ::  line
        if( len_trim(line) == 0 )then
             str_is_blank = .true.
         else
             str_is_blank = .false.
         end if
    end function str_is_blank

    !>  \brief  Find the first non-blank character in a string and return its position.
    pure integer function first_non_blank( str, back )
        character(len=*), intent(in)  :: str !< input string
        logical, optional, intent(in) :: back !<  if .true., we'll look for the last non-blank character in the string
        logical ::  bback
        ! reverse?
        bback = .false.
        if (present(back)) bback = back
        first_non_blank = verify(str,BLANKS_AND_NEW_LINES,back=bback)
    end function first_non_blank

    !>  \brief  find the first blank character in a string and return its position.
    pure integer function first_blank( str, back )
        character(len=*), intent(in)  :: str !< input string
        logical, optional, intent(in) :: back !<  if .true., we'll look for the last blank character in the string
        logical ::  bback
        ! reverse?
        if (present(back)) then
            bback = back
        else
            bback = .false.
        endif
        first_blank = scan(str,BLANKS_AND_NEW_LINES,back=bback)
    end function first_blank

    !>  \brief  test whether two strings are equal, ignoring blank characters
    logical function strings_are_equal( instr1, instr2, case_sensitive )
        character(len=*), intent(in)  :: instr1 !< input string
        character(len=*), intent(in)  :: instr2 !< input string
        logical, optional, intent(in) :: case_sensitive  !< is the comparison case-sensitive
        integer :: first_non_blank_1, first_non_blank_2
        integer :: last_non_blank_1, last_non_blank_2
        logical :: ccase_sensitive
        character(len=:), allocatable ::  local_str1, local_str2
        ! Will we be case sensitive?
        ccase_sensitive = .true.
        if( present(case_sensitive) ) ccase_sensitive = case_sensitive
        ! Copy to local string, with or without capitalization
        if( ccase_sensitive )then
            local_str1 = trim(adjustl(instr1))
            local_str2 = trim(adjustl(instr2))
        else
            local_str1 = upperCase(instr1)
            local_str2 = upperCase(instr2)
        endif
        ! Find positions of first and last non-blank characters
        first_non_blank_1 = verify(local_str1,blank_characters)
        last_non_blank_1  = verify(local_str1,blank_characters,back=.true.)
        first_non_blank_2 = verify(local_str2,blank_characters)
        last_non_blank_2  = verify(local_str2,blank_characters,back=.true.)
        if( first_non_blank_1 .eq. 0 .and. first_non_blank_2 .eq. 0 )then
            ! both strings are blank
            strings_are_equal = .true.
        else if( first_non_blank_1 .eq. 0 .or. first_non_blank_2 .eq. 0 )then
            ! one string is blank, the other isn't
            strings_are_equal = .false.
        else if( index(local_str1(first_non_blank_1:last_non_blank_1), &
                       local_str2(first_non_blank_1:last_non_blank_2)) .eq. 1  &
                 .and. last_non_blank_1-first_non_blank_1 .eq. last_non_blank_2-first_non_blank_2 )then
            ! neither string is blank, and the strings match
            strings_are_equal = .true.
        else
            ! all other cases
            strings_are_equal = .false.
        endif
    end function strings_are_equal

    !>  \brief  count the number of records in a line
    function cnt_recs_per_line( line, separators ) result( nrecs )
        character(len=*)            :: line        !<  line to be split
        character(len=*), optional  :: separators  !<  characters which separate words, if not present, default is blank characters (space, tabs...)
        character(len=LINE_MAX_LEN) :: buffer
        integer :: nrecs, pos1, pos2
        if (str_is_blank(line) .or. str_is_comment(line)) then
            nrecs = 0
        else
            buffer = trim(line)
            nrecs = 0; pos2 = 0; pos1 = 0
            do
                if (present(separators)) then
                    ! find the first non-separator character
                    pos1 = verify(buffer,separators,back=.false.)
                    if( pos1 .ne. 0 )then
                        ! find the last non-separator character after that
                        pos2 = pos1+scan(buffer(pos1:),separators,back=.false.)-2
                    endif
                else
                    ! find the first non-blank character
                    pos1 = first_non_blank(buffer)
                    if (pos1 .ne. 0) then
                        ! find the last non-blank character after that
                        pos2 = pos1+first_blank(buffer(pos1:))-2
                    endif
                endif
                if( pos1 .ne. 0 .and. pos2 .ne. 0 )then ! if we found a word
                    nrecs = nrecs + 1
                endif
                if( pos2+1 .ge. len_trim(buffer) )then ! if we reached end of the line
                    exit
                else
                    buffer = buffer(pos2+1:)
                endif
            enddo
        endif
    end function cnt_recs_per_line

    pure function findloc_str_1( str_arr, str ) result( pos )
        class(string), intent(in) :: str_arr(:)
        class(string), intent(in) :: str 
        integer :: i, pos
        pos = 0
        do i = 1, size(str_arr)
            if( str_arr(i) .eq. str )then
                pos = i
                return
            endif
        end do
    end function findloc_str_1 

    pure function findloc_str_2( str_arr, str ) result( pos )
        class(string),    intent(in) :: str_arr(:)
        character(len=*), intent(in) :: str 
        integer :: i, pos
        pos = 0
        do i = 1, size(str_arr)
            if( str_arr(i) .eq. trim(adjustl(str)) )then
                pos = i
                return
            endif
        end do
    end function findloc_str_2

    !> \brief  Lexographical sort.
    !> \param strArr is a one-dimensional array of character strings to be  sorted in ascending lexical order.
    !>   the sorted array. The characters of
    !>         the elements of the string array are not modified. If blanks or punctuation characters are
    !>         to be ignored, this needs to be taken care of before calling.
    subroutine lex_sort( strArr, CaseSens, inds )
        character(len=*),               intent(inout) :: strArr(:)
        logical, optional,              intent(in)    :: CaseSens  !< case-sensitive sorting
        integer, optional, allocatable, intent(out)   :: inds(:)
        integer, allocatable :: indexarray(:)
        integer              :: low, high, k
        logical              :: LCaseSens
        LCaseSens    = .false.
        if( present(CaseSens) ) LCaseSens = CaseSens
        low          = 1
        high         = size(strArr)
        allocate(indexarray(high))
        indexarray = (/(k,k=low,high)/)
        call quicksort(low, high)
        strArr = strArr(indexarray)
        if( present(inds) )then
            if( allocated(inds) ) deallocate(inds)
            allocate(inds(high), source=indexarray)
        endif
        deallocate(indexarray)

      contains

          recursive subroutine quicksort( low, high )
              integer, intent (in) :: low, high
              integer :: pivotlocation
              if( low < high )then
                  call partition(low, high, pivotlocation)
                  call quicksort(low, pivotlocation-1)
                  call quicksort(pivotlocation+1, high)
              end if
          end subroutine quicksort

          subroutine partition( low, high, pivotlocation )
              integer, intent (in) :: low, high
              integer, intent(out) :: pivotlocation
              integer :: k, lastsmall
              call swap(indexarray(low), indexarray((low+high)/2))
              lastsmall = low
              do k=low+1,high
                  if(strComp(strArr(indexarray(k)), strArr(indexarray(low))))then
                      lastsmall = lastsmall + 1
                      call swap(indexarray(lastsmall), indexarray(k))
                  end if
              end do
              call swap(indexarray(low), indexarray(lastsmall))
              pivotlocation = lastsmall
          end subroutine partition

          subroutine swap( m, n )
              integer, intent (inout) :: m, n
              integer :: temp
              temp = m
              m = n
              n = temp
          end subroutine swap

          function strComp( p, q ) result( lexLess )
              character(len=*), intent(in) :: p, q
              logical :: lexLess
              integer :: kq, k
              if( LCaseSens )then
                  lexLess = p < q
              else
                  kq = 1
                  do k=1,max(len_trim(p),len_trim(q))
                      if( UpperCaseFirst(p(k:k)) == UpperCaseFirst(q(k:k)) )then
                          cycle
                      else
                          kq = k
                          exit
                      end if
                  end do
                  lexLess = UpperCaseFirst(p(kq:kq)) < UpperCaseFirst(q(kq:kq))
              end if
          end function strComp

          function UpperCaseFirst( letter ) result( L )
              character(len=*), intent (in) :: letter
              character(len=1)              :: L
              character(len=26), parameter  :: Lower = "abcdefghijklmnopqrstuvwxyz"
              character(len=26), parameter  :: Upper = "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
              integer :: k
              k = index(Lower, letter)
              if (k > 0) then
                  L = Upper(k:k)
              else
                  L = letter
              end if
          end function UpperCaseFirst

    end subroutine lex_sort

    !>  \brief  converst to C string
    pure function to_cstring( f_string )result( c_string )
        character(len=*),intent(in)  :: f_string
        character(len=1,kind=c_char) :: c_string(len_trim(f_string)+1)
        integer :: i, n
        n = len_trim(f_string)
        do i=1,n
            c_string(i) = f_string(i:i)
        enddo
        c_string(n+1) = C_NULL_CHAR
    end function to_cstring

end module simple_string_utils
