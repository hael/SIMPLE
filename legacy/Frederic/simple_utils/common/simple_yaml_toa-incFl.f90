!> @file
!! Include file used in yaml_strings.f90
!! Body of the yaml_toa template.
!! yaml: Yet Another Markup Language (ML for Human)
!! @author
!!    Copyright (C) 2013-2013 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS


  integer :: ipos
  character(len=max_value_length) :: str
  character(len=*), optional, intent(in) :: fmt

  !print *,'here',data,fmt
  str=repeat(' ',max_value_length)
  if (present(fmt)) then
     write(str,fmt) data
     !error recovery (is iostat /=0 in these cases?)
     !if the format has failed the result is full of stars
     if (trim(str) == repeat('*',len_trim(str))) then
        write(str,cnv_fmt(data)) data
     else
        !another source of failure is the exponent too short like 10+300 instead of 10E+300
        ipos=max(index(str,'+',back=.true.),index(str,'-',back=.true.)) !only one of these should be present
        if (ipos > 2) then
           if (str(ipos-1:ipos-1) /= 'E' .and. str(ipos-1:ipos-1) /= 'e')&
                write(str,cnv_fmt(data)) data
        end if
     end if
  else
     write(str,cnv_fmt(data)) data
  end if
  !otherwise write it in free format
  if (trim(str) == repeat('*',len_trim(str))) write(str,*) data
  !print *,'hereagain',str,data,fmt
  str=yaml_adjust(str,clean=.not. present(fmt))
