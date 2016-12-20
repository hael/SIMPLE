! ============================================================================
! Name        : simple_yaml_strings
! Author      : Frederic Bonnet
! Version     : 1.0
! Date        : 10th of February 2016
! Description : yaml module for printing output in a yaml format.
!             : This modules manages the string in question.
! ============================================================================
!
module simple_yaml_strings
  use simple_defs
  use simple_cuda_defs
  use simple_file_defs

  implicit none

  private
  !Not a parameter in order to be used by C bindings but constant
  integer, parameter :: max_value_length=95
  !default format used
  character(len=*), parameter :: yaml_int_fmt ='(i0)'
  character(len=*), parameter :: yaml_real_fmt='(1pe18.9)'
  character(len=*), parameter :: yaml_dble_fmt='(1pg26.16e3)'!'(1pe25.17)'
  character(len=*), parameter :: yaml_char_fmt='(a)'
  !operator oevrload via interface
  interface operator(//)
     module procedure string_and_integer,string_and_double,string_and_long
  end interface operator(//)
  !Convert into a character string yaml_toa(xxx,fmt)
  interface yaml_toa       
     module procedure yaml_itoa,yaml_litoa,yaml_dtoa,yaml_rtoa, yaml_ctoa
     !module procedure yaml_ttoa
     module procedure yaml_ltoa
  end interface yaml_toa
  !Give the default format corresponding to the nature of the data
  interface cnv_fmt 
     module procedure fmt_i,fmt_r,fmt_d,fmt_a,fmt_li
  end interface cnv_fmt

  !pulbic methods
  public :: file_strcpy,buffer_string
  public :: is_atoi
  public :: read_fraction_string
  !String writters
  public :: yaml_date_and_time_toa, yaml_date_toa
  public :: yaml_toa
  public :: align_message
  public :: shiftstr
  public :: operator(//)
  
contains
  
  pure function fmt_i(data)
    implicit none
    integer, intent(in) :: data
    character(len=len(yaml_int_fmt)) :: fmt_i
    fmt_i=yaml_int_fmt
  end function fmt_i
  pure function fmt_r(data)
    implicit none
    real(f_simple), intent(in) :: data
    character(len=len(yaml_real_fmt)) :: fmt_r
    fmt_r=yaml_real_fmt
  end function fmt_r
  pure function fmt_d(data)
    implicit none
    real(f_double), intent(in) :: data
    character(len=len(yaml_dble_fmt)) :: fmt_d
    fmt_d=yaml_dble_fmt
  end function fmt_d
  pure function fmt_a(data)
    implicit none
    character(len=*), intent(in) :: data
    character(len=len(yaml_char_fmt)) :: fmt_a
    fmt_a=yaml_char_fmt
  end function fmt_a
  pure function fmt_li(data)
    implicit none
    integer(kind=8), intent(in) :: data
    character(len=len(yaml_int_fmt)) :: fmt_li
    fmt_li=yaml_int_fmt
  end function fmt_li

  !routine to buffer a string
  pure subroutine buffer_string(string,string_length,buffer, &
                                string_pos,back,istat)
    implicit none
    integer, intent(in) :: string_length
    character(len=string_length), intent(inout) :: string
    integer, intent(inout) :: string_pos
    character(len=*), intent(in) :: buffer
    logical, optional, intent(in) :: back
    integer, optional, intent(out) :: istat
    !local variables
    integer :: length_add
    
    if(present(istat)) istat = 0
    length_add = len(buffer)
    if (length_add+string_pos > string_length) then
       length_add = len_trim(buffer)
    end if

    if (length_add+string_pos > string_length) then
       if (present(istat)) then
          istat = -1
          return
       else if (length_add+string_pos > string_length) then
          length_add = string_length-string_pos-1
       end if
    end if

    if (length_add==0) return
    if (present(back)) then
       if (back) then
          call shiftstr(string,length_add)
          string(1:length_add)=buffer(1:length_add)
       else
          string(string_pos+1:string_pos+length_add)=buffer(1:length_add)
       end if
    else
       string(string_pos+1:string_pos+length_add)=buffer(1:length_add)
    end if

    string_pos=string_pos+length_add

    return
  end subroutine buffer_string

  !routine used to align the message
  subroutine align_message(rigid,maxlen,tabval,anchor,message)
    implicit none
    logical, intent(in) :: rigid
    integer, intent(in) :: maxlen
    integer, intent(in) :: tabval
    character(len=*), intent(in) :: anchor
    character(len=maxlen), intent(inout) :: message
    !local variables
    integer :: iscpos,ishift
    !cannot align, tabular too far
    if (tabval>maxlen) return

    iscpos=index(message,anchor)      
    ishift=tabval-iscpos
    if (rigid) then
       call shiftstr(message,ishift)
    else
       message=message(1:iscpos-1)//repeat(' ',ishift)//anchor//&
            message(iscpos+1:maxlen-ishift)  ! shift right 
    end if
    return
  end subroutine align_message
  
  !> Write the strings as they were written by write
  pure subroutine file_strcpy(dest,src)
    implicit none
    character(len=*), intent(out) :: dest
    character(len=*), intent(in) :: src
    !local variables
    integer :: i

    dest=repeat(' ',len(dest))
    
    do i=1,min(len(src),len(dest))
       dest(i:i)=src(i:i)
    end do

    return
  end subroutine file_strcpy
  !Shifts characters in in the string 'str' n positions (positive values
  !denote a right shift and negative values denote a left shift). Characters
  !that are shifted off the end are lost. Positions opened up by the shift 
  !are replaced by spaces.
  !This routine has been downloaded from the website
  !http://gbenthien.net/strings/index.html
  pure subroutine shiftstr(str,n)
    implicit none
    integer, intent(in) :: n
    character(len=*), intent(inout) :: str
    !local variables
    integer :: lenstr,nabs
    
    lenstr=len(str)
    nabs=iabs(n)
    if(nabs>=lenstr) then
       str=repeat(' ',lenstr)
       return
    end if
    if(n<0) str=str(nabs+1:)//repeat(' ',nabs)  ! shift left
    if(n>0) str=repeat(' ',nabs)//str(:lenstr-nabs)  ! shift right 
    return
  end subroutine shiftstr

  pure function is_atoi(str) result(yes)
    implicit none
    character(len=*), intent(in) :: str
    logical :: yes
    !local variables
    !TODO: need to implement the function 
  end function is_atoi

  !Read a real or real/real, real:real 
  !Here the fraction is indicated by the ':' or '/'
  !The problem is that / is a separator for Fortran
  pure subroutine read_fraction_string(string,var,ierror)
    implicit none
    !Arguments
    character(len=*), intent(in) :: string
    double precision, intent(out) :: var
    integer, intent(out) :: ierror
    !Local variables
    character(len=256) :: tmp
    integer :: num,den,pfr,psp

    !First look at the first blank after trim
    tmp(1:len(tmp))=trim(adjustl(string))
    psp = scan(tmp,' ')
    !see whether there is a fraction in the string
    if(psp==0) psp=len(tmp)
    pfr = scan(tmp(1:psp),':')
    if (pfr == 0) pfr = scan(tmp(1:psp),'/')
    !It is not a fraction
    if (pfr == 0) then
       read(tmp(1:psp),*,iostat=ierror) var
    else 
       read(tmp(1:pfr-1),*,iostat=ierror) num
       read(tmp(pfr+1:psp),*,iostat=ierror) den
       if (ierror == 0) var=dble(num)/dble(den)
    end if
    !Value by defaut
    if (ierror /= 0) var = huge(1.d0) 
  end subroutine read_fraction_string

  pure function string_and_integer(a,num) result(c)
    implicit none
    integer(f_integer), intent(in) :: num
    character(len=*), intent(in) :: a
    character(len=len_trim(adjustl(a))+len_trim(yaml_itoa(num))) :: c
    c=a//trim(yaml_toa(num))
  end function string_and_integer
  pure function string_and_long(a,num) result(c)
    implicit none
    integer(f_long), intent(in) :: num
    character(len=*), intent(in) :: a
    character(len=len_trim(adjustl(a))+len_trim(yaml_litoa(num))) :: c
    c=a//trim(yaml_toa(num))
  end function string_and_long
  pure function string_and_double(a,num) result(c)
    implicit none
    real(f_double), intent(in) :: num
    character(len=*), intent(in) :: a
    character(len=len_trim(adjustl(a))+len_trim(yaml_dtoa(num))) :: c
    c=a//trim(yaml_toa(num))
  end function string_and_double

  !implementation of the interface yaml_toa
  !Convert integer to character
  pure function yaml_itoa(data,fmt) result(str)
    implicit none
    integer, intent(in) :: data
    include 'simple_yaml_toa-incFl.f90'
  end function yaml_itoa

  !Convert long integer to character
  pure function yaml_litoa(data,fmt) result(str)
    implicit none
    integer(kind=8), intent(in) :: data
    include 'simple_yaml_toa-incFl.f90'
  end function yaml_litoa

  !Convert double to character
  pure function yaml_dtoa(data,fmt) result(str)
    implicit none
    real(f_double), intent(in) :: data
    include 'simple_yaml_toa-incFl.f90'
  end function yaml_dtoa
  
  !Convert real to character
  pure function yaml_rtoa(data,fmt) result(str)
    implicit none
    real(f_simple), intent(in) :: data
    include 'simple_yaml_toa-incFl.f90'
  end function yaml_rtoa

  !character to character
  pure function yaml_ctoa(d,fmt)
    implicit none
    character(len=*), intent(in) :: d
    character(len=len(d)) :: yaml_ctoa
    character(len=*),optional, intent(in) :: fmt
    if (present(fmt)) then
       write(yaml_ctoa,fmt) trim(d)
    else
       call file_strcpy(src=d,dest=yaml_ctoa)
    end if
  end function yaml_ctoa
  ! Convert logical to character
  pure function yaml_ltoa(l,fmt)
    implicit none
    logical, intent(in) :: l
    character(len=max_value_length) :: yaml_ltoa
    character(len=*), optional, intent(in) :: fmt

    yaml_ltoa=repeat(' ',max_value_length)

    if (present(fmt)) then
       write(yaml_ltoa,fmt) l
    else
       if (l) then
          write(yaml_ltoa,'(a3)')'Yes'
       else
          write(yaml_ltoa,'(a3)')'No'
       end if
    end if

    yaml_ltoa=yaml_adjust(yaml_ltoa)
  end function yaml_ltoa

  !Convert logical to character
  !TODO: need to take care of the 
  !  pure function yaml_ttoa(data,fmt) result(str)
!    implicit none
!    logical, intent(in) :: data
!    include 'simple_yaml_toa-incFl.f90'
!  end function yaml_ttoa
 
  !adjusting method
  pure function yaml_adjust(str,clean)
    implicit none
    character(len=*), intent(in) :: str
    logical, intent(in), optional :: clean
    character(len=max_value_length) :: yaml_adjust
    !local variables
    logical :: clean0

    clean0=.true.
    if (present(clean)) clean0=clean

    yaml_adjust=adjustl(str)

    if (clean0) yaml_adjust=clean_zero(yaml_adjust)

    !put a space if there is no sign
    if (yaml_adjust(1:1)/='-') then
       call shiftstr(yaml_adjust,1)
    else
       call shiftstr(yaml_adjust,0)
    end if

  end function yaml_adjust

  !clean to zero the string
  pure function clean_zero(str)
    implicit none
    character(len=*), intent(in) :: str
    character(len=max_value_length) :: clean_zero
    !local variables
    integer :: idot,iexpo,i

    !first fill with all the values up to the dot if it exist
    idot=scan(str,'.')
    if (idot==0) then
       !no dot, nothing to clean
       clean_zero(1:max_value_length)=str
    else
       !first find the position of the end of the string
!       iend=len_trim(str)
       !then search for the position of the exponent or of the space if present
       iexpo=scan(str(idot+2:),'eE ')+idot+1
       !print *,'there',iexpo,'str',str(idot+2:)
       if (iexpo==idot+1) iexpo=len(str)+1
       i=iexpo
       find_last_zero: do while(i > idot+1) !first digit after comma always stays
          i=i-1
          if (str(i:i) /= '0') exit find_last_zero
       end do find_last_zero
       clean_zero(1:i)=str(1:i)
       !print *,'here',i,clean_zero(1:i),'iexpo',iexpo,str(iexpo:)
       !then copy the exponent
       if (str(iexpo:) /= 'E+00' .and. str(iexpo:) /= 'e+00' .and. str(iexpo:) /= 'E+000' .and. &
            str(iexpo:) /= 'e+000') then
          clean_zero(i+1:max_value_length)=str(iexpo:)
       else
          clean_zero(i+1:max_value_length)=' '
       end if
       !try to put at the old position a termination character
       !clean_zero(iend:iend)=char(0)
    end if
  end function clean_zero
  
  !Yaml Spaced format for Date and Time
  function yaml_date_and_time_toa(values,zone)
    implicit none
    logical, optional, intent(in) :: zone
    integer, dimension(8), optional, intent(in) :: values
    character(len=max_value_length) :: yaml_date_and_time_toa
    !local variables
    character(len=*), parameter :: &
         deffmt='i4.4,"-",i2.2,"-",i2.2," ",i2.2,":",i2.2,":",i2.2,".",i3.3'
    logical :: zon
    integer :: zonhrs,zonmin
    integer, dimension(8) :: vals
    character(len=4) :: sgn

    zon=.false.
    if (present(zone)) zon=zone

    if (present(values)) then
       vals=values
    else
       call date_and_time(values=vals)
    end if

    if (zon) then
       zonmin=abs(mod(vals(4),60))
       zonhrs=abs(vals(4)/60)
       if (vals(4) < 0) then
          sgn='" -"'
       else
          sgn='" +"'
       end if
       write(yaml_date_and_time_toa,'('//deffmt//','//sgn//',i2.2,":",i2.2)')&
            vals(1:3),vals(5:8),zonhrs,zonmin

    else
       write(yaml_date_and_time_toa,'('//deffmt//')')vals(1:3),vals(5:8)
    end if

    !There is no - sign so we skip this step (TD)
    !yaml_date_and_time_toa=yaml_adjust(yaml_date_and_time_toa)

  end function yaml_date_and_time_toa

  !Yaml Spaced format for Date
  function yaml_date_toa(values,zone)
    implicit none
    logical, optional, intent(in) :: zone
    integer, dimension(8), optional, intent(in) :: values
    character(len=max_value_length) :: yaml_date_toa
    !local variables
    character(len=*), parameter :: deffmt='i4.4,"-",i2.2,"-",i2.2'
    logical :: zon
    integer :: zonhrs,zonmin
    integer, dimension(8) :: vals
    character(len=4) :: sgn

    zon=.false.
    if (present(zone)) zon=zone

    if (present(values)) then
       vals=values
    else
       call date_and_time(values=vals)
    end if

    if (zon) then
       zonmin=abs(mod(vals(4),60))
       zonhrs=abs(vals(4)/60)
       if (vals(4) < 0) then
          sgn='" -"'
       else
          sgn='" +"'
       end if
       write(yaml_date_toa,'('//deffmt//','//sgn)vals(1:3)

    else
       write(yaml_date_toa,'('//deffmt//')')vals(1:3)
    end if

    !There is no - sign so we skip this step (TD)
    !yaml_date_and_time_toa=yaml_adjust(yaml_date_and_time_toa)

  end function yaml_date_toa

  

end module simple_yaml_strings
