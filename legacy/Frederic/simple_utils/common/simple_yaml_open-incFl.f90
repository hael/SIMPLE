  !general types used overall
! character(len=*), optional, intent(in) :: advance
  !some local variables
  
  logical :: doflow
  integer :: msg_lgt
  integer :: unt,strm
  integer :: tb,ipos
  character(len=3) :: adv
  character(len=tot_max_record_length) :: towrite
  
  unt=0
  if (present(unit)) unt=unit
  call get_stream(unt,strm)

  doflow=streams(strm)%flowrite
  !override if already active
  if (present(flow)) doflow=flow .or. doflow

  !Position of the cursor
  ipos=max(streams(strm)%icursor,streams(strm)%indent)

  msg_lgt=0
  !put the message
  if (present(mapname) .and. len_trim(mapname)>0) then
     call buffer_string(towrite,len(towrite),trim(mapname),msg_lgt)
     !add some spaces if required
     if (present(tabbing)) then
        ipos=ipos+msg_lgt
        tb=max(tabbing-ipos-1,1)
        call buffer_string(towrite,len(towrite),repeat(' ',tb),msg_lgt)
        ipos=ipos+tb
     end if
     !put the semicolon
     call buffer_string(towrite,len(towrite),':',msg_lgt)
  end if
  !put the optional tag
  if (present(tag).and. len_trim(tag)>0) then
     call buffer_string(towrite,len(towrite),' !',msg_lgt)
     call buffer_string(towrite,len(towrite),trim(tag),msg_lgt)
  end if
  !put the optional name
  if (present(label).and. len_trim(label)>0) then
     call buffer_string(towrite,len(towrite),' &',msg_lgt)
     call buffer_string(towrite,len(towrite),trim(label),msg_lgt)
  end if

  call open_level(streams(strm),doflow)

  if (doflow .or. msg_lgt==0) then
     adv='no '
  else
     adv='yes'
     if (present(advance)) adv = advance
  end if
