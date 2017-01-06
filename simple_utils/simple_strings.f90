module simple_strings
implicit none

contains

    !> \brief  parses the string 'str' into arguments args(1), ..., args(nargs) based on
    !!         the delimiters contained in the string 'delim'. The integer output variable 
    !!         nargs contains the number of arguments found.
    subroutine parse(str,delim,args,nargs)
        character(len=*), intent(inout) :: str
        character(len=*), intent(in)    :: delim
        character(len=*)                :: args(:)
        integer, intent(out)            :: nargs
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
    
    !> \brief  removes spaces, tabs, and control characters in string str
    subroutine removesp(str)
        character(len=*), intent(inout) :: str
        character(len=1):: ch
        character(len=len_trim(str))::outstr
        integer :: lenstr, k, i, ich
        str=adjustl(str)
        lenstr=len_trim(str)
        outstr=''
        k=0
        do i=1,lenstr
            ch=str(i:i)
            ich=iachar(ch)
            select case(ich)    
                case(0:32)  ! space, tab, or control character
                    cycle       
                case(33:)  
                    k=k+1
                    outstr(k:k)=ch
            end select
        end do
        str=adjustl(outstr)
    end subroutine removesp
    
    !> \brief  converts number string to a single precision real number
    function str2real(str) result(rval)
        character(len=*), intent(in) ::str
        real    :: rval
        integer :: io_stat
        read(str, *, iostat=io_stat) rval
        if( io_stat .ne. 0 )then
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
    
end module simple_strings
