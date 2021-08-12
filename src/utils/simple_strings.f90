! contains various string manipulation subroutines and functions
! adapted from George Benthien's string module http://gbenthien.net/strings/str-index.html
! modification and additions by Cyril Reboul, Michael Eager & Hans Elmlund
module simple_strings
use simple_defs   ! singleton
use simple_error, only: simple_exception
use, intrinsic :: iso_c_binding
implicit none

character(len=*), parameter :: LOWER_CASE_LETTERS = 'abcdefghijklmnopqrstuvwxyz'
character(len=*), parameter :: UPPER_CASE_LETTERS = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
character(len=*), parameter :: INTEGERS = '1234567890'
character(kind=c_char,len=*), parameter :: NEW_LINES_C = C_FORM_FEED // &
C_NEW_LINE // C_CARRIAGE_RETURN // C_VERTICAL_TAB
character(kind=c_char,len=*), parameter :: BLANK_C_CHARACTERS = C_NULL_CHAR // C_HORIZONTAL_TAB
character(len=*), parameter :: BLANK_CHARACTERS =  ' '//BLANK_C_CHARACTERS
character(len=*), parameter :: BLANKS_AND_NEW_LINES = BLANK_CHARACTERS // NEW_LINES_C
character(len=*), parameter :: INTEGER_DIGITS = '10' !< Maximum number of digits expected when reading an integer

private :: LOWER_CASE_LETTERS, UPPER_CASE_LETTERS, INTEGERS, NEW_LINES_C, BLANK_C_CHARACTERS,&
BLANK_CHARACTERS, BLANKS_AND_NEW_LINES, INTEGER_DIGITS
public

contains

    pure function spaces( n ) result( str )
        integer, intent(in) :: n
        character(len=:), allocatable :: str
        if( n <= 0 ) return
        allocate(character(len=n) :: str)
        str(:) = ' '
    end function spaces

    !> \brief  assesses whether a string represents a filename
    subroutine str2format( str, format, rvar, ivar )
        character(len=*),              intent(in)  :: str
        character(len=:), allocatable, intent(out) :: format
        real,                          intent(out) :: rvar
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
    end subroutine str2format

    function str2latex( str2convert ) result( str_latex )
        character(len=*), intent(in)  :: str2convert
        character(len=:), allocatable :: str_latex, str
        allocate(str_latex, source=adjustl(trim(str2convert)))
        call replace_substr(str_latex, '{', '\{', str)
        deallocate(str_latex)
        call replace_substr(str, '}', '\}', str_latex)
        deallocate(str)
        call replace_substr(str_latex, '_', '\_', str)
        deallocate(str_latex)
        call replace_substr(str, '%', '\%', str_latex)
        deallocate(str)
        call replace_substr(str_latex, '#', '\#', str)
        ! deallocate(str_latex)
        ! call replace_substr(str, 'in A', 'in \AA{}', str_latex)
        ! THE ISSUE HERE IS THAT THE SUBSTRING CANNOT HAVE SPACES
        str_latex = str
    end function str2latex

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
    subroutine parsestr(str,delim,args,nargs)
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
            call split(str,delim,args(nargs))
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
    subroutine split_str(str, delim, before)
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
    subroutine split(str,delim,before)
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
    end subroutine split

    !> \brief  removes punctuation (except comma) characters in string str
    subroutine removepunct(str)
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
    pure subroutine removesp(str)
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

    !> \brief  converts number string to a single precision real number
    function str2real(str) result(rval)
        character(len=*), intent(in) :: str
        real    :: rval
        integer :: io_stat
        character(len=100) :: io_msg
        read(str, *, iostat=io_stat, iomsg=io_msg) rval
        if( io_stat .ne. 0 )then
            call simple_exception('I/O: '//trim(io_msg)//'; str2real', __FILENAME__ , __LINE__)
        endif
    end function str2real

    !> \brief  converts a real number to a string
    pure function real2str(rval) result(str)
        real, intent(in)  :: rval
        character(len=32) :: str
        write(str,*) rval
        call removesp(str)
    end function real2str

    !> \brief  converts a real number to a string
    pure function dbl2str(rval) result(str)
        real(dp), intent(in)  :: rval
        character(len=32) :: str
        write(str,*) rval
        call removesp(str)
    end function dbl2str

    !>  \brief  turn integer variable into character variable
    !! - now supports negative values
    pure function int2str(intg) result(string)
        integer, intent(in)           :: intg
        character(len=:), allocatable :: string
        integer :: ndigs_int, intg_this
        logical :: isnegative
        isnegative=.false.
        intg_this=intg
        if( intg < 0 )then
            isnegative=.true.
            intg_this =  -intg
        end if
        if (intg_this .eq. 0) then
            ndigs_int = 1
        else
            ndigs_int = int(log10(real(intg_this))) + 1
        endif
        if (isnegative) ndigs_int = ndigs_int + 1
        allocate(character(len=ndigs_int) :: string)
        if (isnegative)then
            write(string,'(a1,i0)') '-',intg_this
        else
            write(string,'(i0)') intg_this
        end if
    end function int2str

    !>  \brief  turn integer variable into zero padded character variable
    function int2str_pad( intg, numlen ) result(string)
        integer, intent(in)           :: intg, numlen
        character(len=:), allocatable :: string, str_tmp, str_tmp2
        integer :: slen
        if( numlen == 0 )then
            string = int2str(intg)
            return
        endif
        str_tmp = int2str(intg)
        slen    = len(str_tmp)
        if( slen >= numlen )then
            allocate(string, source=str_tmp)
        else
            do while( len(str_tmp) < numlen )
                allocate(str_tmp2, source='0'//str_tmp)
                if( allocated(string) ) deallocate(string)
                allocate(string, source=str_tmp2)
                deallocate(str_tmp)
                allocate(str_tmp, source=str_tmp2)
                deallocate(str_tmp2)
            end do
        endif
        deallocate(str_tmp)
    end function int2str_pad

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

    !>  \brief  turns a string into an integer variable
    pure subroutine str2int( string, io_stat, ivar )
        character(len=*), intent(in)  :: string
        integer,          intent(out) :: io_stat, ivar
        read(string,*, iostat=io_stat) ivar
    end subroutine str2int

    !>  \brief  maps out the number characters in a string
    !!          .true. means is number, .false. means is not
    pure function map_str_nrs( str ) result( map )
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

    !> \brief copies the trimmed and justified version of the instr into outstr
    subroutine cpStr(instr, outstr)
        character(len=*),              intent(in)    :: instr
        character(len=:), allocatable, intent(inout) :: outstr
        outstr = trim(adjustl(instr))
    end subroutine cpStr

    !>  \brief  check for presence of substring
    logical function str_has_substr( str, substr )
        character(len=*), intent(in) :: str, substr
        str_has_substr = .not. (index(str, substr) == 0)
    end function str_has_substr

    !>  \brief  removes occurrences of substr in str into str_out
    subroutine void_substr( str, substr, str_out)
        character(len=*), intent(in)                 :: str, substr
        character(len=:), allocatable, intent(inout) :: str_out
        integer :: l, lout, pos
        str_out = trim(str)
        l = len_trim(substr)
        if( l == 0 )return
        if( .not.str_has_substr(str_out, substr) )return
        pos = index(str_out, substr)
        do while( pos > 0 )
            lout = len_trim(str_out)
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
    end subroutine

    subroutine replace_substr( str, substr, repl, str_out )
        character(len=*), intent(inout) :: str
        character(len=*), intent(in)    :: substr, repl
        character(len=:), allocatable   :: str_out
        integer :: nargs, iarg, iback, ls
        character(len=LONGSTRLEN) :: args(32)
        logical :: back_occur
        iback      = index(str, substr, back=.true.)
        ls         = len(substr)
        back_occur = len(str) - ls + 1 == iback
        call parsestr(str,substr,args,nargs)
        if( nargs == 0 )then
            allocate(str_out, source=str)
            return
        else
            str_out = ''
            do iarg=1,nargs
                if( iarg == nargs .and. .not. back_occur )then
                    str_out = str_out // trim(args(iarg))
                else
                    str_out = str_out // trim(args(iarg)) // repl
                endif
            end do
        endif
    end subroutine replace_substr

    !>  \brief  Convert string to all upper case
    !!  Adapted from
    !!  Figure 3.5B, pg 80, "Upgrading to Fortran 90", by Cooper Redwine,
    !!  1995 Springer-Verlag, New York.
    function upperCase( str )
        character(len=*),   intent(in)  ::  str
        character(len=len(str))  ::  upperCase
        integer :: i, n
        upperCase = str
        do i=1,len(str)
            n = index(LOWER_CASE_LETTERS,upperCase(i:i))
            if( n .ne. 0 ) upperCase(i:i) = UPPER_CASE_LETTERS(n:n)
        enddo
    end function upperCase

    !>  \brief works out whether a character string is a comment line
    logical function strIsComment( line )
        character(len=*), intent(in)  ::  line
        integer ::  pos1,linelen
        strIsComment = .false.
        if(strIsEmpty(line)) return
        pos1 = firstNonBlank(line)
        ! if we didn't find a non-blank character, then the line must be blank!
        if( pos1 .eq. 0 ) return
        linelen=len_trim(line)
        if( scan(line(pos1:),'#!;', back=.false.) .eq. 1 )then
            ! the first non-blank character is a comment character.
            strIsComment = .true.
        endif
    end function strIsComment

    !>  \brief works out whether a character string is blank
    pure logical function strIsEmpty(line)
        character(len=*),   intent(in)  ::  line
        if( len_trim(line) == 0 )then
             strIsEmpty = .true.
         else
             strIsEmpty = .false.
         end if
    end function strIsEmpty

    !>  \brief works out whether a character string is blank
    pure logical function strIsBlank(line)
        character(len=*),   intent(in)  ::  line
        if( trim(line) .eq. '' )then
             strIsBlank = .true.
        else if( firstNonBlank(line) .eq. 0 )then
             strIsBlank = .true.
        else
             strIsBlank = .false.
        endif
    end function strIsBlank

    !>  \brief  Find the first non-blank character in a string and return its position.
    pure integer function firstNonBlank( string, back )
        character(len=*), intent(in)  :: string !< input string
        logical, optional, intent(in) :: back !<  if .true., we'll look for the last non-blank character in the string
        logical ::  bback
        ! reverse?
        bback = .false.
        if (present(back)) bback = back
        firstNonBlank = verify(string,blanks_and_new_lines,back=bback)
    end function firstNonBlank

    !>  \brief  find the first blank character in a string and return its position.
    pure integer function firstBlank(string, back)
        character(len=*), intent(in)  :: string !< input string
        logical, optional, intent(in) :: back !<  if .true., we'll look for the last blank character in the string
        logical ::  bback
        ! reverse?
        if (present(back)) then
            bback = back
        else
            bback = .false.
        endif
        firstBlank = scan(string,blanks_and_new_lines,back=bback)
    end function firstBlank

    !>  \brief  test whether two strings are equal, ignoring blank characters
    logical function stringsAreEqual(instr1, instr2, case_sensitive)
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
            call cpStr(instr1,local_str1)
            call cpStr(instr2,local_str2)
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
            stringsAreEqual = .true.
        else if( first_non_blank_1 .eq. 0 .or. first_non_blank_2 .eq. 0 )then
            ! one string is blank, the other isn't
            stringsAreEqual = .false.
        else if( index(local_str1(first_non_blank_1:last_non_blank_1), &
                       local_str2(first_non_blank_1:last_non_blank_2)) .eq. 1  &
                 .and. last_non_blank_1-first_non_blank_1 .eq. last_non_blank_2-first_non_blank_2 )then
            ! neither string is blank, and the strings match
            stringsAreEqual = .true.
        else
            ! all other cases
            stringsAreEqual = .false.
        endif
    end function stringsAreEqual

    !>  \brief  count the number of records in a line
    function cntRecsPerLine( line, separators ) result( nrecs )
        character(len=*)            :: line        !<  line to be split
        character(len=*), optional  :: separators  !<  characters which separate words, if not present, default is blank characters (space, tabs...)
        character(len=LINE_MAX_LEN) :: buffer
        integer :: nrecs, pos1, pos2
        if (strIsBlank(line) .or. strIsComment(line)) then
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
                    pos1 = firstNonBlank(buffer)
                    if (pos1 .ne. 0) then
                        ! find the last non-blank character after that
                        pos2 = pos1+firstBlank(buffer(pos1:))-2
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
    end function cntRecsPerLine

    !! reverse string
    recursive function strflip (string) result (res)
        implicit none
        character (*), intent (in) :: string
        character (len (string)) :: res
        if (len (string) == 0) then
            res = ''
        else
            res = string(len(string):)//strflip(string(:len(string)-1))
        end if
    end function strflip

    !> \brief  Lexographical sort.
    !> \param strArr is a one-dimensional array of character strings to be  sorted in ascending lexical order.
    !>   the sorted array. The characters of
    !>         the elements of the string array are not modified. If blanks or punctuation characters are
    !>         to be ignored, this needs to be taken care of before calling.
    subroutine lexSort( strArr, CaseSens, inds )
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
        allocate(indexarray(high), stat=alloc_stat)
        if( alloc_stat /= 0 ) then
            write(logfhandle,'(a)') 'ERROR: Allocation failure!'
            write(logfhandle,'(a)') 'In: lexSort; simple_strings'
            stop
        endif
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

    end subroutine lexSort

    !>  \brief  converst to C string
    pure function toCstring( f_string )result( c_string )
        character(len=*),intent(in)  :: f_string
        character(len=1,kind=c_char) :: c_string(len_trim(f_string)+1)
        integer :: i, n
        n = len_trim(f_string)
        do i=1,n
            c_string(i) = f_string(i:i)
        enddo
        c_string(n+1) = C_NULL_CHAR
    end function toCstring

end module simple_strings
