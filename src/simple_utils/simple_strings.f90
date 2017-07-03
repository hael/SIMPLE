module simple_strings
use simple_defs ! singleton
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
character(len=*), parameter :: INTEGER_DIGITS = '10' ! Maximum number of digits expected when reading an integer

private :: LOWER_CASE_LETTERS, UPPER_CASE_LETTERS, INTEGERS, NEW_LINES_C, BLANK_C_CHARACTERS,&
BLANK_CHARACTERS, BLANKS_AND_NEW_LINES, INTEGER_DIGITS
public

contains

    !> \brief  assesses whether a string represents a filename
    function str2format( str) result( format )
        character(len=*), intent(in)  :: str
        character(len=:), allocatable :: format
        integer :: io_stat, ifoo
        logical :: str_is_file
        str_is_file = .true.
        if( index(str, '.spi') /= 0 )then
        else if( index(str, '.mrc') /= 0 )then
        else if( index(str, '.map') /= 0 )then
        else if( index(str, '.ctf') /= 0 )then
        else if( index(str, '.bin') /= 0 )then
        else if( index(str, '.txt') /= 0 )then
        else if( index(str, '.asc') /= 0 )then
        else if( index(str, '.log') /= 0 )then
        else if( index(str, '.dat') /= 0 )then
        else if( index(str, '.box') /= 0 )then
        else if( index(str, '.raw') /= 0 )then
        else if( index(str, '.hed') /= 0 )then
        else if( index(str, '.img') /= 0 )then
        else if( index(str, '.ccp') /= 0 )then
        else
            str_is_file = .false.
        endif
        if( str_is_file )then
            allocate( format, source='file' )
            return
        endif
        if( index(str, '.') /= 0 )then
            allocate( format, source='real' )
            return
        endif
        call str2int(str, io_stat, ifoo )
        if( io_stat==0 )then
            allocate( format, source='int' )
            return
        endif
        allocate( format, source='char' )
    end function str2format

    !> \brief  parses the string 'str' into arguments args(1), ..., args(nargs) based on
    !!         the delimiters contained in the string 'delim'. The integer output variable
    !!         nargs contains the number of arguments found.
    subroutine parse(str,delim,args,nargs)
        character(len=*), intent(inout) :: str
        character(len=*), intent(in)    :: delim
        character(len=*)                :: args(:)
        integer,          intent(out)   :: nargs
        character(len=len_trim(str))    :: strsav
        integer                         :: na, i, lenstr, k
        strsav=str
        call compact(str)
        na=size(args)
        do i=1,na
            args(i)=''
        end do
        nargs=0
        lenstr=len_trim(str)
        if(lenstr==0) return
        k=0
        do
            if(len_trim(str) == 0) exit
            nargs=nargs+1
            call split(str,delim,args(nargs))
        end do
        str=strsav
    end subroutine parse

    !> \brief  converts multiple spaces and tabs to single spaces;
    !!         deletes control characters; removes initial spaces.
    subroutine compact(str)
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

    !> \brief  replace any tokens in string (str) with char rch
    subroutine replace(str, tokens, rch)
    character(len=*), intent(inout) :: str
    character(len=*), intent(inout) :: tokens     ! array of searchable chars
    character(len=1), intent(in) :: rch     ! replace any element in tokens with char
    character(len=1):: ch,tok
    character(len=len_trim(str))::outstr
    integer :: lenstr,lensstr, k, i,j, ich, itoken
    str = adjustl(str)
    tokens = adjustl(tokens)
    lensstr = len_trim(tokens)
    outstr = adjustl(str) ! initialise outstr
    do j=1,lensstr
        str = adjustl(outstr)  ! update str for each element in sstr in case rch == ''
        lenstr = len_trim(str) ! re-adjust length
        outstr = ''            ! reset outstr
        k = 0
        tok = tokens(j:j)
        itoken = iachar(tok)
        do i=1,lenstr
            ch = str(i:i)
            ich = iachar(ch)
            if(ich == itoken) then  ! character in schs
                if ( len_trim(rch) == 0 ) then
                    cycle
                else
                    k = k + 1
                    outstr(k:k) = rch
                end if
            else
                k = k+1
                outstr(k:k) = ch
            end if
        end do
    end do
    str = adjustl(outstr)
    end subroutine replace

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
        case(32:43)  ! ! " # $ % & ' ( ) * +  characters
            cycle
        case(45:47)  ! . - / characters
            cycle
        case(58:64)  ! : ; < = > ? @ characters
            cycle
        case(91:94)  !  _ ^ [ ] characters
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
    subroutine removesp(str)
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
        character(len=*), intent(in) ::str
        real    :: rval
        integer :: io_stat
        character(len=100) :: io_msg
        read(str, *, iostat=io_stat, iomsg=io_msg) rval
        if( io_stat .ne. 0 )then
            write(*,*) 'iomsg: ', trim(io_msg)
            stop 'ERROR(I/O); simple_strings :: str2real'
        endif
    end function str2real

    !> \brief  converts a real number to a string
    function real2str(rval) result(str)
        real, intent(in)  :: rval
        character(len=32) :: str
        write(str,*) rval
        call removesp(str)
    end function real2str

    !>  \brief  turn integer variable into character variable
    function int2str(intg) result(string)
        integer, intent(in)           :: intg
        character(len=:), allocatable :: string
        integer :: ndigs_int
        if( intg < 0 )then
            stop 'cannot convert negative intg 2 str; simple_strings :: int2str'
        else if (intg .eq. 0) then
            ndigs_int = 1
        else
            ndigs_int = int(log10(real(intg))) + 1
        endif
        allocate(character(len=ndigs_int) :: string)
        write(string,'(i0)') intg
    end function int2str

    !>  \brief  turn integer variable into zero padded character variable
    function int2str_pad(intg, numlen) result(string)
        integer, intent(in)           :: intg, numlen
        character(len=:), allocatable :: string, str_tmp, str_tmp2
        integer :: str_tmp_len, slen
        str_tmp = int2str(intg)
        slen    = len(str_tmp)
        if( slen > numlen )then
            print *, 'len(int2str(intg)): ', slen
            print *, 'numlen            : ', numlen
            stop 'length of number > than desired length; simple_strings :: in2str_pad'
        else if( len(str_tmp) == numlen )then
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

    !>  \brief  turns a string into an integer variable
    subroutine str2int( string, io_stat, ivar )
        character(len=*), intent(in)  :: string
        integer,          intent(out) :: io_stat, ivar
        read(string,*, iostat=io_stat) ivar
    end subroutine str2int

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
        integer ::  pos1
        pos1 = firstNonBlank(line)
        ! if we didn't find a non-blank character, then the line must be blank!
        if( pos1 .eq. 0 )then
            strIsComment = .false.
        else if( scan(line(pos1:),'#!c;', back=.false.) .eq. 1 )then
            ! the first non-blank character is a comment character. if that's a c, we need to check that the following character is a space.
            if( stringsAreEqual(line(pos1:pos1),'c') )then
                if( stringsAreEqual(line(pos1+1:pos1+1),' ') )then
                    strIsComment = .true.
                else
                    strIsComment = .false.
                endif
            else
                strIsComment = .true.
            endif
        else
            strIsComment = .false.
        endif
    end function strIsComment

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
        character(len=*), intent(in)  :: string
        logical, optional, intent(in) :: back !<  if .true., we'll look for the last non-blank character in the string
        logical ::  bback
        ! reverse?
        bback = .false.
        if (present(back)) bback = back
        firstNonBlank = verify(string,blanks_and_new_lines,back=bback)
    end function firstNonBlank

    !>  \brief  find the first blank character in a string and return its position.
    pure integer function firstBlank(string, back)
        character(len=*), intent(in)  :: string
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
        character(len=*), intent(in)  :: instr1
        character(len=*), intent(in)  :: instr2
        logical, optional, intent(in) :: case_sensitive
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
        character(len=line_max_len) :: buffer
        integer :: nrecs, pos1, pos2
        if (strIsBlank(line) .or. strIsComment(line)) then
            nrecs = 0
        else
            buffer = trim(line)
            nrecs = 0
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

    !> \brief  Lexographical sort. On input, strArr is a one-dimensional array of character strings to be
    !!         sorted in ascending lexical order. On output, strArr is the sorted array. The characters of
    !!         the elements of the string array are not modified. If blanks or punctuation characters are
    !!         to be ignored, this needs to be taken care of before calling.
    subroutine lexSort( strArr, CaseSens )
        character(len=*),  intent(inout) :: strArr(:)
        logical, optional, intent(in)    :: CaseSens
        integer, allocatable :: indexarray(:)
        integer              :: low, high, alloc_stat, k
        logical              :: LCaseSens
        LCaseSens = .false.
        if( present(CaseSens) ) LCaseSens = CaseSens
        low  = 1
        high = size(strArr)
        allocate(indexarray(high), stat=alloc_stat)
        if( alloc_stat /= 0 ) then
            write(*,'(a)') 'ERROR: Allocation failure!'
            write(*,'(a)') 'In: lexSort; simple_strings'
            stop
        endif
        indexarray = (/(k,k=low,high)/)
        call quicksort(low, high)
        strArr = strArr(indexarray)
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
                      if( UpperCase(p(k:k)) == UpperCase(q(k:k)) )then
                          cycle
                      else
                          kq = k
                          exit
                      end if
                  end do
                  lexLess = UpperCase(p(kq:kq)) < UpperCase(q(kq:kq))
              end if
          end function strComp

          function UpperCase( letter ) result( L )
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
          end function UpperCase

    end subroutine lexSort

end module simple_strings
