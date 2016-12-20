! ============================================================================
! Name        : simple_yaml_output
! Author      : Frederic Bonnet
! Version     : 1.0
! Date        : 10th of February 2016
! Description : yaml module for printing output in a yaml format
!             :
! ============================================================================
!
module simple_yaml_output
  use simple_defs
  use simple_cuda_defs
  use simple_yaml_strings
  use simple_eglossary
  use simple_error_handling
  use simple_eglossary_lowlev
  implicit none

  private
  !> Yaml events for dump routine
  integer, parameter :: NONE                   = -1000
  integer, parameter :: DOCUMENT_START         = -1001
  integer, parameter :: DOCUMENT_END           = -1002
  integer, parameter :: MAPPING_START          = -1003
  integer, parameter :: MAPPING_END            = -1004
  integer, parameter :: SEQUENCE_START         = -1005
  integer, parameter :: SEQUENCE_END           = -1006
  integer, parameter :: SCALAR                 = -1007
  integer, parameter :: COMMENT                = -1008
  integer, parameter :: MAPPING                = -1009
  integer, parameter :: SEQUENCE_ELEM          = -1010
  integer, parameter :: NEWLINE                = -1011
  integer, parameter :: COMMA_TO_BE_PUT        =  10
  integer, parameter :: DEFAULT_STREAM_ID      =  0

  integer, parameter :: tot_max_record_length=95 !Max record length by default
  integer, parameter :: tot_max_flow_events=500  !Max flow events
  integer, parameter :: tot_streams=10           !Max total number of streams
  integer, parameter :: tab=5                    !Default number for tabbing

  integer, parameter :: n_streams = 10           !Total number of streams

  integer :: active_streams=0!N streams (stdout always active after init)
  integer :: default_stream=1!Id of the default stream

  type yaml_stream !stream for the document 1,..,n 
     logical :: document_closed=.true. !start closed doc
     logical :: pp_allowed=.true.      !Pretty printing allowed
     integer :: indent = 1             !spaces for yaml output
     integer :: indent_previous=0      !Indent level prior to flow writing
     integer :: indent_step=2          !Indentation level
     integer :: unit = 6               !strandard output
     integer :: icursor=1              !Running cursor position on the line
     integer :: tabref=40              !Position of tabular in scalar assignm
     integer :: ilevel=0               !Number of opened levels
     integer :: iflowlevel=0           !flowrite levels simoultaneously enabled
     integer :: ilast=0                !Last level with flow==.false.
     integer :: Wall=-1                !Warning messages of level Wall(-1:none)
     logical :: flowrite=.false.       !Write in flow (.false.=no .true.=yes)
     integer :: ievt_flow=0            !Track events in flowrite
     integer :: icommentline=0         !Active if commented line being written
     integer :: itab=0                 !Tabbing to have a look on
     integer :: itab_active=0          !# active tabbings for line in flowrite
     integer :: max_record_length=tot_max_record_length
     integer, dimension(tot_max_flow_events) :: flow_events=0 !set flow event
     integer, dimension(tot_max_record_length/tab) :: linetab=0!Value of the tabbing in the line
     type(eglossary), pointer :: egloss_warning=>null() !eglossary warnings
  end type yaml_stream

  ! Global variables used through out
  type(yaml_stream), dimension(n_streams), save :: streams
  integer, dimension(n_streams) :: stream_units = 6 !default unit
  type(eglossary), pointer :: stream_files

  logical :: module_initialized=.false.!module already referenced or not
  ! Error ids
  integer :: YAML_STREAM_ALREADY_PRESENT !Trying to create a stream already present
  integer :: YAML_STREAM_NOT_FOUND       !Trying to seach for a absent unit
  integer :: YAML_UNIT_INCONSISTENCY     !Internal error, unit inconsistency
  integer :: YAML_INVALID                !Invalid action, unit inconsistency

  !interfaces
  interface yaml_map
     module procedure yaml_map, yaml_map_egloss
     !general types
     module procedure yaml_map_li,yaml_map_i
     module procedure yaml_map_f, yaml_map_d
     module procedure yaml_map_l
  end interface yaml_map
  
  !Here all the public routines used
  public :: yaml_new_document,yaml_release_document
  public :: yaml_map, yaml_map_open, yaml_map_close
  public :: yaml_sequence, yaml_sequence_open, yaml_sequence_close
  public :: yaml_scalar
  !yaml setters and getters
  public :: yaml_set_stream
  !warning and error checkers
  public :: yaml_warning, yaml_comment
  !dumpers
  public :: yaml_write
  public :: yaml_egloss_write
  public :: yaml_output_errors

contains
  
  !Initialize the yaml output errors, error messages
  subroutine yaml_output_errors()
    implicit none
    !initialize error messages
    call file_err_define('YAML_INVALID','Generic error of yaml module, invalid operation',&
         YAML_INVALID)
    call file_err_define('YAML_STREAM_ALREADY_PRESENT','The stream is already present',&
         YAML_STREAM_ALREADY_PRESENT)
    call file_err_define('YAML_STREAM_NOT_FOUND','The stream has not been found',&
         YAML_STREAM_NOT_FOUND)
    call file_err_define('YAML_UNIT_INCONSISTENCY',&
         'The array of the units is not in agreement with the array of the streams',&
         YAML_UNIT_INCONSISTENCY,&
         err_action='This is an internal error of yaml_output module, contact developers')
    !the module is ready for usage
    call egloss_init(stream_files)
    module_initialized=.true.
  end subroutine yaml_output_errors

  !TODO: move back to section of the setters and getters
  !yaml_set_stream
  subroutine yaml_set_stream(unit,filename,istat,tabbing,record_length,&
                             position,setdefault)
    use simple_file_utils
    implicit none
    integer, optional, intent(in) :: unit !unit user specified(default=6)
    integer, optional, intent(in) :: tabbing!Indicate a tabbing for the stream
    integer, optional, intent(in) :: record_length!Max # columns stream
    character(len=*), optional, intent(in) :: filename
    character(len=*), optional, intent(in) :: position !unit position
    integer, optional, intent(out) :: istat      !Status
    logical, optional, intent(in) :: setdefault  !new stream default or not
    !local variables
    character(len=2) :: char_out
    integer :: unt,recl_file
    logical :: unit_is_open,set_default
    integer :: err,rc
    character(len=15) :: pos
    !indexers
    integer :: istream
    !function calls
    integer :: convert_int2char_pos_c
    integer :: strlen

    !TODO: to be implenmented, sets the stream
    if (present(unit)) then
       unt=unit
    else
       if (present(filename)) then
          if (has_key(stream_files,trim(filename))) then
             unt = stream_files // trim(filename)
          else
             unt = file_get_free_unit()
          end if
       else
          unt=6
       end if
    end if

    !write(*,*) "yaml_set_stream unt: ",unt
    !open a fortran unit if needed
    recl_file = 0
    if (present(filename) .and. unt /= 6 ) then
       !TODO: need to implement if this condition occurs
       !write(*,*) "in present filename condition: ",unt

       !inquire whether file already exists
       inquire(unit=unt,opened=unit_is_open,iostat=err)
       !TODO: throw a yaml error using file_err in error_handling module
       rc = convert_int2char_pos_c(char_out,abs(err))
       if (file_err(err /=0,'error in unit inquiring, err='//trim(char_out),&
            YAML_INVALID)) return
       if (unit_is_open) then
          if (present(istat)) then
             istat=YAML_STREAM_ALREADY_PRESENT
          else
             rc = convert_int2char_pos_c(char_out,abs(unt))
             call file_err_throw('The unit '//trim(char_out)//' is already present',&
                  YAML_STREAM_ALREADY_PRESENT)
          end if          
       end if
       if (present(position))then
          pos(1:strlen(pos))=position
       else
          pos(1:strlen(pos))='append'
       end if

       if (.not. unit_is_open) then

          inquire(file=trim(filename),opened=unit_is_open,iostat=err)
          rc = convert_int2char_pos_c(char_out,abs(err))
          if (file_err(err /=0,'error in file inquiring, ierr='//trim(char_out),&
               YAML_INVALID)) return
          if (unit_is_open) then
             if(present(istat)) then
                istat=YAML_STREAM_ALREADY_PRESENT
             else
                call file_err_throw('The file '//trim(filename)//' is already connected',&
                     YAML_STREAM_ALREADY_PRESENT)
             end if
          end if
          call file_err_open_try()
          call file_open(trim(filename),unt,position=trim(pos))
          err=file_get_last_error()
          call file_err_close_try()
          if (present(istat)) istat=err

       end if

       if (err == 0 .and. .not. unit_is_open) then
          !inquire the record length for the unit
          !inquire(unit=unt,recl=recl_file)
          if (present(record_length)) &
               call file_recl(unt,record_length,recl_file)
          !throw to yaml format
          rc = convert_int2char_pos_c(char_out,abs(unt))
          !write(*,*)"char_out: ", char_out," unt: ",unt
          !TODO: need to fix the segmentaion fault due to uninitialized pointer
          !write(*,*)"in set_stream stream_files%data%key:",stream_files%data%key
          !write(*,*)"in set_stream stream_files%data%value:",stream_files%data%value
          call set(stream_files//trim(filename), trim(char_out))
       end if
    end if

    !write(*,*)"line 225 in set_stream_files unt: ",unt," char_out: ",char_out
    
    if (present(setdefault)) then
       set_default=setdefault
    else
       set_default = .true.
    end if
    !check if the unit has already been assigned
    do istream=1,active_streams
       if (unt==stream_units(istream)) then
          !throw an error
          rc = convert_int2char_pos_c(char_out,abs(unt))
          if (file_err(.not. present(istat),'Unit: '//trim(char_out)&
               //' already present',err_id=YAML_STREAM_ALREADY_PRESENT)) then
             return
          else
             istat=YAML_STREAM_ALREADY_PRESENT
             return
          end if
       end if
    end do
    !If there are no active streams set default cannot be false 
    if (.not. set_default .and. active_streams==0) then
       active_streams=active_streams+1
       !initialise the stream
       streams(active_streams)=stream_null()
       streams(active_streams)%unit=6
       stream_units(active_streams)=6
       streams(active_streams)%max_record_length=92
    end if

    active_streams=active_streams+1
    !initialise the stream
    streams(active_streams)=stream_null()
    streams(active_streams)%unit=unt
    stream_units(active_streams)=unt
    !set last opened stream as default stream
    if (set_default) default_stream=active_streams

    if (present(tabbing)) then
       streams(active_streams)%tabref=tabbing
       if(tabbing==0) streams(active_streams)%pp_allowed=.false.
    end if

    !keep record length to be lower than the maximum allowed by the processor
    if (present(record_length)) then
       if (recl_file<=0) recl_file=record_length
       streams(active_streams)%max_record_length=recl_file
    end if

    return
  end subroutine yaml_set_stream

  !Nullifyer of the stream
  function stream_null() result(strm)
    implicit none
    type(yaml_stream) :: strm
    !TODO: need to add here entries as yaml_stream grows 
    strm%document_closed=.true. !start closed doc
    strm%pp_allowed=.true.
    strm%indent = 1             !spaces for yaml output
    strm%unit = 6               !strandard output
    strm%icursor=1              !Running cursor position on the line
    strm%max_record_length=tot_max_record_length
    strm%flow_events=0 !set flow event
    strm%tabref=40
  end function stream_null
  
  !comments
  subroutine yaml_comment(message,advance,unit_in,hfill,tabbing)
    implicit none
    character(len=*), intent(in) :: message !comment (with #)
    character(len=*), optional, intent(in) :: advance
    integer, optional, intent(in) :: unit_in
    character(len=*), optional, intent(in) :: hfill !fill line wth character
    integer, optional, intent(in) :: tabbing   !Number of space for tabbing
    !local variables
    integer :: unit,strm,ipos,msg_lgt,tb
    integer :: lstart,lmsg,lend,lspace,hmax
    character(len=3) :: adv
    character(len=tot_max_record_length) :: towrite
    integer :: strlen
    !start of the execution commands
    unit = DEFAULT_STREAM_ID
    if (present(unit_in)) unit = unit_in
    call get_stream(unit,strm)

    !comment to be written
    if(present(advance)) then
       adv=advance
    else
       adv='yes'
    end if

    !Beginning the message
    lstart = 1
    !length of message
    lmsg = len_trim(message)
    !write(*,*) "message:",message," , length of message: ",lmsg

    !split the message if too long
    do
       !position the cursor
       ipos = max(streams(strm)%icursor,streams(strm)%indent)
       msg_lgt = 0
       if (present(tabbing)) then
          tb=max(tabbing-ipos-1,1)
          call buffer_string(towrite,len(towrite),repeat(' ',tb),msg_lgt)
          !call buffer_string(towrite,len(towrite), &
          !                   message(lstart:lstart+lend-1),msg_lgt)
          ipos=ipos+tb
       end if
       !look for the last character
       lend=len_trim(message(lstart:))
       if  (lend+msg_lgt+2 > streams(strm)%max_record_length ) then
          lend = streams(strm)%max_record_length-msg_lgt-2
          lspace = index(message(lstart:lstart+lend-1),' ',back=.true.)
          if (lspace /= 0) then
             lend = lspace
          end if
       end if
       call buffer_string(towrite,len(towrite),message(lstart:lstart+lend-1),msg_lgt)

       !write(*,*)'there: ',trim(towrite),' ,',lstart,lend

       !check if possible to hfill
       hmax = max(streams(strm)%max_record_length-ipos-len_trim(message)-3,0)
       if (present(hfill)) hmax=hmax/len(hfill)
       if (present(hfill) .and. hmax > 0) then
          call yaml_write(streams(strm),repeat(hfill,hmax)//' '//towrite(1:msg_lgt),advance=adv,event=COMMENT)
       else
          call yaml_write(streams(strm),towrite(1:msg_lgt),advance=adv,event=COMMENT)
       end if

       lstart=lstart+lend
       if (lstart>lmsg) then
          exit
       end if
    end do
    return
  end subroutine yaml_comment
  
  !warners and error checkers
  subroutine yaml_warning(message,level,unit_in)
    implicit none
    character(len=*),intent(in) :: message
    integer, optional, intent(in) :: level
    integer, optional, intent(in) :: unit_in
    !local variables
    integer :: istream
    integer :: unt
    integer :: idx
    type(eglossary),pointer :: egloss_tmp
    
    unt = DEFAULT_STREAM_ID
    if (present(unit_in)) unt = unit_in
    call get_stream(unt,istream)

    call yaml_comment('WARNING: '//trim(message),'no',unit_in=unt)
    !TODO: implement warning message for generic message
    if (.not. streams(istream)%document_closed) then

       if (.not. associated(streams(istream)%egloss_warning)) &
            call egloss_init(streams(istream)%egloss_warning)
       !add warning to list
       !TODO: Can't convert TYPE(list_container) to TYPE(eglossary) at (1)
       !egloss_tmp = streams(istream)%egloss_warning .get. 'WARNING'
       idx=egloss_tmp .index. trim(message)
       !TODO: need to implement the add interface in lowlev with type method
    end if
    if (present(level)) then
       if (level <= streams(istream)%Wall) then
          call yaml_write(streams(istream),' Critical warning level reached, aborting...')
          !TODO: need to implement the method below
          call yaml_release_document(unit=unit_in)
          stop
       end if
    end if

    return
  end subroutine yaml_warning
  
  !Method to release document from the warning method
  subroutine yaml_release_document(unit)
    implicit none
    integer, optional, intent(in) :: unit
    !TODO: need to implement the method
    !yaml_new_line needs to be implemented
    !yaml_egloss_dump and yaml_flush_document

    return
  end subroutine yaml_release_document
  
  !writters and dumpers
  subroutine dump()
    implicit none
    !TODO: implement the dump method, traditional yaml_dump
    return
  end subroutine dump
  !writing out the eglossary
  subroutine yaml_egloss_write(egloss,unit,flow,verbatim)
    implicit none
    type(eglossary),intent(in),pointer :: egloss
    integer,intent(in) , optional :: unit
    logical,intent(in),optional :: flow
    logical,intent(in),optional :: verbatim
    !local variables
    logical :: flowrite,verb,default_flow
    integer :: unt

    !TODO: implement the dumper method for the eglossary
    !write(*,*)
    !write(*,*) "******Need to implement the yaml_egloss_dump*****"
    !write(*,*)

    flowrite=.false.
    default_flow=.true.
    if (present(flow))then
       flowrite=flow
       default_flow=.false.
    end if
    unt=DEFAULT_STREAM_ID
    if(present(unit))unt=unit
    verb=.false.
    if(present(verbatim))verb=verbatim

    if (.not.associated(egloss))then
       call scalar('null')
    elseif(associated(egloss%child))then
       call yaml_egloss_dump_(egloss%child)
    else
       call scalar(egloss_value(egloss))
    end if
    
  contains

    recursive subroutine yaml_egloss_dump_(egloss)
      implicit none
      type(eglossary),pointer,intent(in) :: egloss
      !local variables
      
      !TODO: implement the dumper

      if(.not.associated(egloss)) return
      if(associated(egloss%child)) then
         if(egloss_item(egloss)>=0) call sequence(adv='no')

         if(egloss_len(egloss)>0) then
            if(switch_flow(egloss)) flowrite=.true.
            call open_seq(egloss_key(egloss))
            call yaml_egloss_dump_(egloss%child)
            call close_seq()
            !restoe normal flowriting if it has been changed before
            if (switch_flow(egloss)) flowrite=.false.
         else
            if (egloss_item(egloss) >= 0) then
               call yaml_egloss_dump_(egloss%child)
            else
               call open_map(egloss_key(egloss))
               call yaml_egloss_dump_(egloss%child)
               call close_map()
            end if
         end if
      else
         if (egloss_item(egloss) >= 0) then
            call sequence(val=egloss_value(egloss))
         else
            call map(egloss_key(egloss),egloss_value(egloss))
         end if
      end if
      !applying the recursion
      call yaml_egloss_dump_(egloss_next(egloss))
      
    end subroutine yaml_egloss_dump_
    !wrappers code to the yaml_map_open
    subroutine open_map(key)
      implicit none
      character(len=*), intent(in) :: key
      if (verb) then
         call yaml_comment('call yaml_mapping_open("'//trim(key)//&
              '",flow='//trim(flw(flowrite))//&
              ',unit='//trim(adjustl(yaml_toa(unt)))//')')
      else
         call yaml_map_open(trim(key),flow=flowrite,unit=unt)
      end if
    end subroutine open_map
    !wrappers code to the yaml_map_open
    subroutine close_map()
      implicit none
      if (verb) then
         call yaml_comment('call yaml_mapping_close('//&
              'unit='//trim(adjustl(yaml_toa(unt)))//')')
      else
         call yaml_map_close(unit=unt)
      end if
    end subroutine close_map
    !wrappers code to the yaml_map
    subroutine map(key,val)
      implicit none
      character(len=*), intent(in) :: key,val
      if (verb) then
         call yaml_comment('call yaml_map("'//trim(key)//'","'//trim(val)//&
              '",unit='//trim(adjustl(yaml_toa(unt)))//'")')
      else
         call yaml_map(trim(key),trim(val),unit=unt)
      end if
    end subroutine map
    !wrappers code to the yaml_sequence_open
    subroutine open_seq(key)
      implicit none
      character(len=*), intent(in) :: key
      if (verb) then
         call yaml_comment('call yaml_sequence_open("'//trim(key)//&
              '",flow='//trim(flw(flowrite))//&
              ',unit='//trim(adjustl(yaml_toa(unt)))//')')
      else
         call yaml_sequence_open(trim(key),flow=flowrite,unit=unt)
      end if
    end subroutine open_seq
    !wrappers code to the yaml_sequence_close
    subroutine close_seq()
      implicit none
      if (verb) then
         call yaml_comment('call yaml_sequence_close('//&
              'unit='//trim(adjustl(yaml_toa(unt)))//')')
      else
         call yaml_sequence_close(unit=unt)
      end if
    end subroutine close_seq

    !wrappers code to the yaml_sequence
    subroutine sequence(adv,val)
      implicit none
      character(len=*), intent(in), optional :: val,adv

      if (present(val) .and. present(adv)) then
         if (verb) then
            call yaml_comment('call yaml_sequence("'//trim(val)//&
                 '",advance="'//trim(adv)//&
                 ',unit='//trim(adjustl(yaml_toa(unt)))//'")')
         else
            call yaml_sequence(trim(val),advance=adv,unit=unt)
         end if
      else if (present(adv)) then
         if (verb) then
            call yaml_comment('call yaml_sequence(advance="'//trim(adv)//&
                 ',unit='//trim(adjustl(yaml_toa(unt)))//'")')
         else
            call yaml_sequence(advance=adv,unit=unt)
         end if
      else if (present(val)) then
         if (verb) then
            call yaml_comment('call yaml_sequence("'//trim(val)//&
                 ',unit='//trim(adjustl(yaml_toa(unt)))//'")')
         else
            call yaml_sequence(trim(val),unit=unt)
         end if
      end if
    end subroutine sequence

    !wrapper code over yaml_scalar method
    subroutine scalar(val)
      implicit none
      character(len=*), intent(in) :: val
      if (verb) then
         call yaml_comment('call yaml_scalar("'//trim(val)//'",advance="'//trim(advc(flowrite))//&
              '",unit='//trim(adjustl(yaml_toa(unt)))//')')
      else
         call yaml_scalar(trim(val),advance=advc(flowrite),unit=unt)
      end if
    end subroutine scalar

    !handlers
    function switch_flow(egloss)
      implicit none
      type(eglossary), pointer, intent(in) :: egloss
      logical :: switch_flow
      
      switch_flow=default_flow .and. last_level(egloss) &
           .and. egloss_len(egloss) <=5 .and. egloss_len(egloss) > 1
      
    end function switch_flow

    !determine if we are at the last level of the eglossary
    function last_level(egloss)
      implicit none
      type(eglossary), pointer, intent(in) :: egloss
      logical :: last_level
      !local variables
      type(eglossary), pointer :: egloss_tmp
      
      !we should have a sequence of only scalar values
      last_level = egloss_len(egloss) > 0
      if (last_level) then
         egloss_tmp=>egloss_next(egloss%child)
         do while(associated(egloss_tmp))
            if (egloss_len(egloss_tmp) > 0 .or. egloss_size(egloss_tmp) > 0) then
               last_level=.false.
               nullify(egloss_tmp)
            else
               egloss_tmp=>egloss_next(egloss_tmp)
            end if
         end do
      end if
    end function last_level

    function advc(flow_tmp)
      implicit none
      logical, intent(in) :: flow_tmp
      character(len=3) :: advc

      if (flow_tmp) then
         advc='no '
      else
         advc='yes'
      end if

    end function advc

    function flw(flow_tmp)
      implicit none
      logical, intent(in) :: flow_tmp
      character(len=7) :: flw

      if (flow_tmp) then
         flw='.true.'
      else
         flw='.false.'
      end if
    end function flw

  end subroutine yaml_egloss_write

  !yaml writter
  subroutine yaml_write(stream,message,advance,event,istat)
    implicit none
    type(yaml_stream), intent(inout) :: stream          !< Stream to handle
    character(len=*), intent(in) :: message             !< Message to dump
    character(len=*), intent(in), optional :: advance   !< Advance option
    integer, intent(in), optional :: event              !< Event to handle
    integer, intent(out), optional :: istat             !< Status error
    !llocal variables
    integer :: evt,indent_lgt,msg_lgt,prefix_lgt,shift_lgt
    integer :: towrite_lgt
    logical :: ladv,reset_line,change_line,pretty_print,reset_tabbing
    logical :: comma_postponed,extra_line
    character(len=3) :: adv
    character(len=5) :: prefix
    character(len=1) :: anchor
    character(len=stream%max_record_length) :: towrite
    !TODO: need to implement the writer method
    if (present(istat)) istat=0
    if (present(event)) then
       evt=event
    else
       evt=SCALAR
    end if

    ladv=.not.stream%flowrite
    if (present(advance)) then
       if (trim(advance)=='no' .or. trim(advance)=='NO') then
          ladv=.false.
       else if (trim(advance)=='yes' .or. trim(advance)=='YES') then
          ladv=.true.          
       end if
    end if
    if (ladv) then
       adv='yes'
    else
       adv='no'
    end if

    !reset line {yes|no(no)}
    reset_line=.false.
    !line changer {yes|no(no)}
    change_line=.false.
    !indentation
    indent_lgt=indent_value(stream,evt)
    !get number of active objects to be written
    towrite=repeat(' ',len(towrite))
    msg_lgt=0
    if(.not. present(istat)) then
       if(len_trim(message)>0) call buffer_string(towrite,len(towrite),message,msg_lgt)
    else
       call buffer_string(towrite,len(towrite),message,msg_lgt,istat=istat)
    end if
    prefix_lgt=0
    prefix=repeat(' ',len(prefix))
    if (put_comma(stream,evt)) then
       call buffer_string(prefix,len(prefix),', ',prefix_lgt)
    end if
    !postpoing the next comma
    comma_postponed=comma_not_needed(evt) .or. (flow_is_ending(evt) .and. stream%iflowlevel==1)
    !default no pretty printing 
    pretty_print=.false.
    shift_lgt=0
    !reset the tabbing
    reset_tabbing=.false.

    !use an event select data structure 
    select case(evt)
    case(SEQUENCE_START)
       if (.not.stream%flowrite) then
          call open_indent_level(stream)
       else
          call buffer_string(towrite,len(towrite),' [',msg_lgt)
          stream%flowrite=.true.
          !added for prettty printing
          reset_tabbing=.true.
       end if
    case(SEQUENCE_END)
       if (.not.stream%flowrite) then
          call close_indent_level(stream)
       else
          if (stream%iflowlevel>1 .and. ladv) then
             call buffer_string(prefix,len(prefix),' ]',prefix_lgt,back=.true.)
             stream%flowrite=.true.
          else
             call buffer_string(prefix,len(prefix),' ]',prefix_lgt)
          end if
          reset_line=ladv
       end if
    case(MAPPING_START)
       if (.not.stream%flowrite) then
          call open_indent_level(stream)
       else
          call buffer_string(towrite,len(towrite),' {',msg_lgt)
          stream%flowrite=.true.
          reset_tabbing=.true.
       end if
    case(MAPPING_END)
       if (.not.stream%flowrite) then
          call close_indent_level(stream)
       else
          if (stream%iflowlevel>1 .and. ladv) then
             call buffer_string(prefix,len(prefix),' }',prefix_lgt,back=.true.)
             reset_line=.true.
          else
             call buffer_string(prefix,len(prefix),' }',prefix_lgt)
          end if
          reset_line=ladv
       end if
    case(COMMENT)
       if (stream%icommentline==0) then
          call buffer_string(prefix,len(prefix),' #',prefix_lgt)
       end if
       if (.not. ladv) then
          stream%icommentline=1
       else
          reset_line=.true.
       end if
    case(MAPPING)
       pretty_print=.true. .and. stream%pp_allowed
       anchor=':'
    case(SEQUENCE_ELEM)
       if (.not.stream%flowrite) then
          indent_lgt=indent_lgt-2
          call buffer_string(prefix,len(prefix),'- ',prefix_lgt)
       else
          if (msg_lgt>0)comma_postponed=.false.
          pretty_print=.true. .and. stream%pp_allowed
          anchor='.'
       end if
    case(SCALAR)
       !TODO: implement if required. not needed yet
    case(NEWLINE)

       if (stream%flowrite) then
          change_line=.true.
          stream%flowrite=.true.
          reset_line=ladv
          msg_lgt=0
       else
          change_line=.true.
          reset_line=.true.
          msg_lgt=0
       end if
    end select
    !adjust the towrite string to match with the closest tabular
    if (pretty_print) then
       call pretty_printing(stream%flowrite,anchor,towrite,&
            stream%icursor,indent_lgt,prefix_lgt,&
            msg_lgt,stream%max_record_length,shift_lgt,change_line)
    end if

    extra_line=.false.
    if (change_line) then
       !first write prefix, if needed
       if (prefix_lgt>0) then
          write(stream%unit,'(a)')prefix(1:prefix_lgt)
       else if (msg_lgt >0 .or. evt == NEWLINE) then
          !change line, only if istat is not present
          if (.not. present(istat)) then
             write(stream%unit,*)
          else
             extra_line=.true.
          end if
       end if
       stream%icursor=1
       towrite_lgt=msg_lgt+shift_lgt
    else
       call shiftstr(towrite,prefix_lgt)
       if (prefix_lgt > 0)towrite(1:prefix_lgt)=prefix(1:prefix_lgt)
       towrite_lgt=prefix_lgt+msg_lgt+shift_lgt
    end if
    !print *,'adv',trim(adv),towrite_lgt,stream%icursor,extra_line,msg_lgt,towrite_lgt
    !here we should check whether the size of the string exceeds the maximum length
    if (towrite_lgt > 0) then
       if (towrite_lgt > stream%max_record_length) then
          if (present(istat)) then
             istat=-1
             return
          else
             !crop the writing
             towrite_lgt=stream%max_record_length
             !print *,'towrite', repeat(' ',max(indent_lgt,0))//towrite(1:towrite_lgt),' end'
             !stop 'ERROR (dump): writing exceeds record size'
          end if
       else
          if (extra_line) write(stream%unit,*)
!!!!!!!!!!!!!!!!!!!!!!!!!!!! TODO: need to fix copmments !!!!!!!!!!!!!!!!!!!!
!          write(*,fmt='(a,i0,a)',advance="no") '(indent_lgt ',indent_lgt,')'
!TODO: here need to fix the stream unit number
!          write(*,*)"------------In yaml_write--------------"
!          write(*,'(a)',advance=trim(adv))&
!               repeat(' ',max(indent_lgt,0))//towrite(1:towrite_lgt)
!          write(*,*)"---------------------------------------"
!          write(*,*) "in yaml_write, stream%unit: ", stream%unit

          write(1,'(a)',advance=trim(adv))&
               repeat(' ',max(indent_lgt,0))//towrite(1:towrite_lgt)

          write(stream%unit,'(a)',advance=trim(adv))&
               repeat(' ',max(indent_lgt,0))//towrite(1:towrite_lgt)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       end if
    end if

    !if advancing i/o cursor is again one
    if (ladv) then
       stream%icursor=1
    else
       !cursor after writing
       stream%icursor=stream%icursor+indent_lgt+towrite_lgt
    end if

    if (reset_tabbing) then
       stream%itab_active=0
       stream%itab=0
    end if

    if (reset_line) call carriage_return(stream)

    !keep history of the event for a flowrite
    !needed for the comma
    if (stream%flowrite) then
       stream%ievt_flow=modulo(stream%ievt_flow,tot_max_flow_events)+1 !to avoid boundary problems
       if (comma_postponed) then
          stream%flow_events(stream%ievt_flow)=evt
       else
          stream%flow_events(stream%ievt_flow)=COMMA_TO_BE_PUT
       end if
    else
       stream%ievt_flow=0
    end if
    
  contains
    !TODO: include the helpers

    subroutine pretty_printing(rigid,anchor,message,icursor,  &
                               indent_lgt,prefix_lgt,msg_lgt, &
                               max_lgt,shift_lgt,change_line)
      implicit none
      logical, intent(in) :: rigid
      integer, intent(in) :: icursor,prefix_lgt,msg_lgt,max_lgt
      integer, intent(inout) :: indent_lgt
      character(len=*), intent(in) :: anchor
      character(len=*), intent(inout) :: message
      logical, intent(out) :: change_line
      integer, intent(out) :: shift_lgt
      !local variables
      integer :: iscpos,ianchor_pos,tabeff

      change_line=.false.
      iscpos=index(message,anchor)
      shift_lgt=0

      ianchor_pos=icursor+prefix_lgt+indent_lgt+iscpos-1
      call closest_tab(ianchor_pos,tabeff)
      !try to see if lines enters
      shift_lgt=tabeff-ianchor_pos

      if (icursor+msg_lgt+prefix_lgt+indent_lgt+shift_lgt >= max_lgt) then
         change_line=.true.
         if (stream%itab==stream%itab_active .and. stream%itab > 1)&
              stream%itab_active=max(stream%itab_active-1,0)
         stream%itab=1

         if (indent_lgt==0) indent_lgt=1
         ianchor_pos=indent_lgt+iscpos
         call closest_tab(ianchor_pos,tabeff)
         shift_lgt=tabeff-ianchor_pos
      end if

      !here the size of the message is known, message adjusted to anchor
      call align_message(rigid,len(message),shift_lgt+iscpos,anchor,message)

      return
    end subroutine pretty_printing

    !Calculate the reference tabular value
    subroutine closest_tab(ianchor_pos,tabeff)
      implicit none
      integer, intent(in) :: ianchor_pos
      integer, intent(out) :: tabeff

      if (stream%flowrite) then
         !check that the tabbing is already done, otherwise add another tab
         if (stream%itab < stream%itab_active) then
            !realign the value to the tabbing
            do
               if (ianchor_pos <= stream%linetab(stream%itab) .or. &
                    stream%itab==stream%itab_active) exit
               stream%itab=modulo(stream%itab,tot_max_record_length/tab)+1
            end do
         end if

         if (stream%itab < stream%itab_active .and. stream%itab>0) then
            tabeff=stream%linetab(stream%itab)
         else
            tabeff=ianchor_pos
            stream%itab=modulo(stream%itab,tot_max_record_length/tab)+1
            stream%itab_active=modulo(stream%itab_active,tot_max_record_length/tab)+1
            stream%linetab(stream%itab_active)=tabeff
         end if
      else
         !for the moment do not check compatibility of the line
         tabeff=max(stream%tabref,ianchor_pos)
      end if
      return
    end subroutine closest_tab

    !comma insertion method
    function put_comma(stream,evt)
      implicit none
      type(yaml_stream),intent(in) :: stream
      integer, intent(inout) :: evt
      logical :: put_comma
      put_comma=stream%flowrite .and. stream%ievt_flow>0
      if (stream%ievt_flow>0) then
         put_comma=stream%flow_events(stream%ievt_flow)==COMMA_TO_BE_PUT
         !TODO: maybe a case when the comma is not needed
      end if
      if (flow_is_ending(evt))put_comma=.false.
    end function put_comma
    
  end subroutine yaml_write

  !When the flow is starting
  function flow_is_starting(evt)
    implicit none
    integer,intent(in) :: evt
    logical :: flow_is_starting
    flow_is_starting = (evt==MAPPING_START .or. evt==SEQUENCE_START)
  end function flow_is_starting
  !When the flow is ending
  function flow_is_ending(evt)
    implicit none
    integer,intent(in) :: evt
    logical :: flow_is_ending
    flow_is_ending = (evt==MAPPING_END .or. evt==SEQUENCE_END)
  end function flow_is_ending
  !when comma is not needed
  pure function comma_not_needed(evt)
    implicit none
    integer, intent(in) :: evt
    logical :: comma_not_needed

    comma_not_needed=evt==NONE           .or. &
                     evt==MAPPING_START  .or. &
                     evt==SEQUENCE_START .or. &
                     evt==SCALAR         .or. &
                     evt==COMMENT        .or. &
                     evt==SEQUENCE_ELEM  .or. &
                     evt==NEWLINE
  end function comma_not_needed
  
  !Increase indentation of stream by changing flow level
  subroutine open_indent_level(stream)
    implicit none
    type(yaml_stream), intent(inout) :: stream !Stream to handle
    stream%indent=stream%indent+stream%indent_step
    return
  end subroutine open_indent_level
  !Decrease indentation of strean without changing flow level
  subroutine close_indent_level(stream)
    implicit none
    type(yaml_stream), intent(inout) :: stream
    stream%indent=max(stream%indent-stream%indent_step,0) !to prevent bugs
  end subroutine close_indent_level
  !Reset line control quantities, and reset the indentation
  subroutine carriage_return(stream)
    implicit none
    type(yaml_stream), intent(inout) :: stream
    !if a yaml_comment is called put the has in front
    stream%icommentline=0
    !beginining of the line
    stream%icursor=1
    !no tabbing decided yet
    stream%itab_active=0
    stream%itab=0
    !all needed commas are placed in the previous line
  end subroutine carriage_return
    
  !Indenting function
  function indent_value(stream,evt)
    implicit none
    type(yaml_stream), intent(in) :: stream    !Stream to handle
    integer, intent(in), optional :: evt       !Event to handle
    !local variables
    integer :: indent_value
    if (.not.stream%flowrite .and. stream%icursor==1) then
       indent_value=stream%indent!max(stream%indent,0) !to prevent bugs
       !if first time in the flow recuperate the saved indent
    else if (stream%icursor==1 .and. stream%iflowlevel==1 &
         .and. stream%ievt_flow==0) then
       indent_value=stream%indent_previous
    else
       indent_value=0!1 !TODO: check indenting value here
       if (stream%icursor==1) indent_value=1
    end if
    if (evt==DOCUMENT_START) indent_value=0
  end function indent_value
  
  !new document methods
  subroutine yaml_new_document(unit_in)
    implicit none
    !global variables
    integer, optional, intent(in) :: unit_in
    !local variables
    integer :: unit,istream

    unit = DEFAULT_STREAM_ID
    if (present(unit_in)) unit = unit_in
    call get_stream(unit,istream)

    if (streams(istream)%document_closed) then
       if (streams(istream)%indent /= 1 ) then
          call yaml_warning("Indentation error",unit_in=stream_units(istream))
          streams(istream)%indent=1
       end if
       call dump()!TODO: remove temporary call dev
       call yaml_write(streams(istream),'---',event=DOCUMENT_START)
       streams(istream)%flow_events = NONE
       streams(istream)%document_closed = .false.
    end if
    
    return
  end subroutine yaml_new_document

  !SETTERS AND GETTERS

  !get the stream
  subroutine get_stream(unit_in,stream_out,istat)
    implicit none
    !global variables
    integer, intent(in) :: unit_in
    integer, intent(out) :: stream_out
    integer, optional,intent(out) :: istat !return code
    !local variables
    logical :: stream_found
    integer :: istream,prev_def,err

    if (present(istat)) istat = 0

    if (unit_in == DEFAULT_STREAM_ID) then
       !if there are no stream activated then do
       if ( active_streams == 0 ) call yaml_set_stream(record_length=92,istat=err)
       stream_out = default_stream
    else
       stream_found = .false.
       do istream = 1, active_streams
          if (stream_units(istream) == unit_in) then
             stream_out = istream
             stream_found = .true.
             exit
          end if
       end do

       if (.not. stream_found) then
          if (present(istat)) then
             istat = YAML_STREAM_NOT_FOUND
          else
             !otherwise activate it
             if (active_streams == 0 .and. unit_in /=6 ) &
                  call yaml_set_stream(record_length=92,istat=err)
             prev_def = default_stream
             call yaml_set_stream(unit=unit_in,tabbing=0)
             stream_out = default_stream
             default_stream = prev_def
          end if
       end if
    end if
    
    return
  end subroutine get_stream
  
  ! MAPPERS METHODS
  
  !mapping a scalar
  subroutine yaml_scalar(message,advance,unit,hfill)
    implicit none
    character(len=*), optional, intent(in) :: hfill
    character(len=*), intent(in) :: message
    integer, optional, intent(in) :: unit
    character(len=*), intent(in), optional :: advance
    !local variables
    integer :: unt,istream,hmax
    character(len=3) :: adv

    unt=DEFAULT_STREAM_ID
    if (present(unit)) unt=unit
    call get_stream(unt,istream)

    !comment to be written
    if (present(advance)) then
       adv=advance
    else
       adv='yes'
    end if
    if (present(hfill)) then
       hmax = max(streams(istream)%max_record_length-&
            max(streams(istream)%icursor,streams(istream)%indent)-&
            len_trim(message)-3,0)
       hmax=hmax/len(hfill)

       call yaml_write(streams(istream),&
            repeat(hfill,hmax)//' '//trim(message),&
            advance=adv,event=COMMENT)
    else
       call yaml_write(streams(istream),trim(message),advance=adv,event=SCALAR)
    end if

    return
  end subroutine yaml_scalar
  !subroutine that opens the yaml map
  subroutine yaml_map_open(mapname,label,tag,flow,tabbing,advance,unit)
    implicit none
    character(len=*), optional, intent(in) :: mapname
    character(len=*), optional, intent(in) :: label
    character(len=*), optional, intent(in) :: tag
    logical, optional, intent(in) :: flow
    character(len=*), optional, intent(in) :: advance
    integer, optional, intent(in) :: unit
    integer, optional, intent(in) :: tabbing
    include 'simple_yaml_open-incFl.f90'
    call yaml_write(streams(strm),towrite(1:msg_lgt),advance=trim(adv),event=MAPPING_START)
    return
  end subroutine yaml_map_open
  !method that closes the yaml map
  subroutine yaml_map_close(advance,unit)
    implicit none
    !TODO: need to implement the new document method
    integer, optional, intent(in) :: unit !< @copydoc doc::unit
    character(len=*), optional, intent(in) :: advance !<@copydoc doc::advance
    !local variables
    integer :: unt,strm
    character(len=3) :: adv
    logical :: doflow

    unt=DEFAULT_STREAM_ID
    if (present(unit)) unt=unit
    call get_stream(unt,strm)

    if (streams(strm)%iflowlevel > 1) then
       adv='no'
    else
       adv='yes'
    end if
    if (present(advance)) adv=advance

    call yaml_write(streams(strm),' ',advance=trim(adv),event=MAPPING_END)

    doflow = (streams(strm)%flowrite)
    call close_level(streams(strm),doflow)
    return
  end subroutine yaml_map_close

  subroutine yaml_map(mapname, mapvalue,label,tag,advance,unit)
    implicit none
    character(len=*), intent(in) :: mapname
    character(len=*), intent(in) :: mapvalue

    character(len=*), optional, intent(in) :: label
    character(len=*), optional, intent(in) :: tag
    character(len=*), optional, intent(in) :: advance
    integer, optional, intent(in) :: unit
    !local varaiables
    integer :: unt, istream
    logical :: cut,redo_line
    integer :: err
    integer :: msg_lgt,icut,istr,msg_lgt_ck,idbg
    character(len=max_field_length) :: lbl
    character(len=3) :: adv
    character(len=tot_max_record_length) :: towrite

    !TODO: need to implement the new document method
    !write(*,*) "yaml_map starts..."
    !write(*,*)"mapname: ", mapname, ", mapvalue: ",mapvalue

    unt=DEFAULT_STREAM_ID
    if (present(unit)) unt=unit
    call get_stream(unt,istream)

    lbl(1:len(lbl))=' '
    if (present(label))lbl(1:len(lbl))=label

    msg_lgt = 0
    !put the message
    call buffer_string(towrite,len(towrite),trim(mapname),msg_lgt)
    !put the semicolon
    call buffer_string(towrite,len(towrite),':',msg_lgt)
    !put the optional tag
    if (present(tag) .and. len_trim(tag) > 0) then
       call buffer_string(towrite,len(towrite),' !',msg_lgt)
       call buffer_string(towrite,len(towrite),trim(tag),msg_lgt)
    end if
    !put the optional name
    if (present(label) .and. len_trim(label) > 0) then
       call buffer_string(towrite,len(towrite),' &',msg_lgt)
       call buffer_string(towrite,len(towrite),trim(label),msg_lgt)
    end if
    !put a space
    call buffer_string(towrite,len(towrite),' ',msg_lgt)

    !while putting the message verify that the string is not too long
    msg_lgt_ck=msg_lgt
    !write(*,*) 'here'
    if (len_trim(mapvalue) == 0) then
       call buffer_string(towrite,len(towrite),"null",msg_lgt,istat=err)
    else
       call buffer_string(towrite,len(towrite),trim(mapvalue),msg_lgt,istat=err)
    end if
    !write(*,*) 'here2',err
    if (err ==0) then
       call yaml_write(streams(istream),towrite(1:msg_lgt),&
                       advance=trim(adv),event=MAPPING,istat=err)
    end if

    !print *, 'here2b',err
    redo_line=err/=0
    if (redo_line) then
       !print *, 'here3',err,msg_lgt_ck,msg_lgt
       if (streams(istream)%flowrite) then
          call yaml_write(streams(istream),towrite(1:msg_lgt_ck),advance=trim(adv),event=SCALAR)
       else
          if (present(label)) then
             call yaml_map_open(mapname,label=label,unit=unt)
          else
             call yaml_map_open(mapname,unit=unt)
          end if
       end if
!       if (streams(istream)%flowrite) call yaml_newline(unit=unt)
       !first, if the cursor is already gone, carriage return
       if (streams(istream)%icursor >= streams(istream)%max_record_length) &
            call yaml_write(streams(istream),' ',advance='yes',event=SCALAR)
       icut=len_trim(mapvalue)
       istr=1
       cut=.true.
       msg_lgt=0
       idbg=0
       cut_line: do while(cut)
          idbg=idbg+1
          !print *,'hereOUTPU',cut,icut,idbg,streams(istream)%icursor,streams(istream)%max_record_length
       !verify where the message can be cut
          !print *,'test2',index(trim((mapvalue(istr:istr+icut-1))),' ',back=.true.)
          cut=.false.
          cut_message :do while(icut > streams(istream)%max_record_length - &
               max(streams(istream)%icursor,streams(istream)%indent))
             icut=index(trim((mapvalue(istr:istr+icut-1))),' ',back=.true.)
             !print *,'test',icut,streams(istream)%max_record_length,&
             !     max(streams(istream)%icursor,streams(istream)%indent),&
             !     streams(istream)%max_record_length - &
             !     max(streams(istream)%icursor,streams(istream)%indent),istr
             cut=.true.
          end do cut_message
          !if the first line is too long cut it abruptly
          if (icut == 0) icut = streams(istream)%max_record_length - &
               max(streams(istream)%icursor,streams(istream)%indent)+1
          call buffer_string(towrite,len(towrite),mapvalue(istr:istr+icut-1),msg_lgt)
          if (streams(istream)%flowrite .and. .not. cut) &
               call buffer_string(towrite,len(towrite),',',msg_lgt)
          call yaml_write(streams(istream),towrite(1:msg_lgt),advance='yes',event=SCALAR)
          istr=istr+icut
          icut=len_trim(mapvalue)-istr+1
          !print *,'icut',istr,icut,mapvalue(istr:istr+icut-1),cut,istr+icut-1,len_trim(mapvalue)
          msg_lgt=0
         if (idbg==1000) exit cut_line !to avoid infinite loops
       end do cut_line
       if (.not.streams(istream)%flowrite) call yaml_map_close(unit=unt)
    end if
   
    return
  end subroutine yaml_map
  !mapps a eglossary
  subroutine yaml_map_egloss(mapname,mapvalue,label,unit,flow)
    implicit none
    character(len=*), intent(in) :: mapname
    type(eglossary), pointer, intent(in) :: mapvalue
    character(len=*), optional, intent(in) :: label
    integer, optional, intent(in) :: unit
    logical, optional, intent(in) :: flow
    !local variables
    integer :: unt, istream
    character(len=max_field_length) :: lbl

    unt=DEFAULT_STREAM_ID
    if (present(unit)) unt=unit
    call get_stream(unt,istream)

    lbl(1:len(lbl))=' '
    if (present(label)) lbl(1:len(lbl))=label

    if (associated(mapvalue)) then
       if (present(flow)) then
          call yaml_map_open(mapname,label=lbl,flow=flow,unit=unt)
          !TODO: the egloss dumper need to be tested
          call yaml_egloss_write(mapvalue,unit=unt,flow=flow)
          !write(*,*) "TODO: the egloss dumper need to be tested"
!          stop
       else
          call yaml_map_open(mapname,label=lbl,unit=unt)
          !TODO: the egloss dumper need to be tested
          call yaml_egloss_write(mapvalue,unit=unt)
          !write(*,*) "TODO: the egloss dumper need to be tested"
!          stop
       end if
       call yaml_map_close(unit=unt)
    else
       call yaml_map(mapname,'<nullified eglossary>',label=lbl,unit=unt)
    end if

    return
  end subroutine yaml_map_egloss
  
  subroutine yaml_map_li(mapname,mapvalue,label,advance,unit,fmt)
    implicit none
!    character(len=*), intent(in) :: mapname
    integer(kind=8), intent(in) :: mapvalue
!    character(len=*), optional, intent(in) :: label
!    character(len=*), optional, intent(in) :: advance
!    character(len=*), optional, intent(in) :: fmt
!    integer, optional, intent(in) :: unit
    include 'simple_yaml_map-incFl.f90'
  end subroutine yaml_map_li

  subroutine yaml_map_i(mapname,mapvalue,label,advance,unit,fmt)
    implicit none
    integer(kind=4), intent(in) :: mapvalue

!    character(len=*), intent(in) :: mapname
!    character(len=*), optional, intent(in) :: label
!    character(len=*), optional, intent(in) :: advance
!    character(len=*), optional, intent(in) :: fmt
!    integer, optional, intent(in) :: unit
    !TODO: need to throw to file simple_yaml_map-incFl.f90
    include 'simple_yaml_map-incFl.f90'
    !    include 'simple_yaml_open-incFl.f90'
!    integer :: msg_lgt
!    integer :: unt,strm
!    integer :: tb,ipos
!    character(len=3) :: adv
!    character(len=tot_max_record_length) :: towrite

!    unt = 0
!    if (present(unit))unt=unit
!    call get_stream(unt,strm)
    
!    adv='def' !default value
!    if (present(advance)) adv=advance

!    msg_lgt = 0
    !put the message
!    call buffer_string(towrite,len(towrite),trim(mapname),msg_lgt)
    !put the semicolon
!    call buffer_string(towrite,len(towrite),': ',msg_lgt)
    !put the optional name
!    if (present(label)) then
!       call buffer_string(towrite,len(towrite),' &',msg_lgt)
!       call buffer_string(towrite,len(towrite),trim(label)//' ',msg_lgt)
!    end if
    !put the value
!    if (present(fmt)) then
!       call buffer_string(towrite,len(towrite),trim(yaml_toa(mapvalue,fmt=fmt)),msg_lgt)
!    else
!       call buffer_string(towrite,len(towrite),trim(yaml_toa(mapvalue)),msg_lgt)
!    end if
    !TODO: need to implement the dumping method called the yaml_write
    !    call dump(streams(strm),towrite(1:msg_lgt),advance=trim(adv),event=MAPPING)

    return
  end subroutine yaml_map_i

  subroutine yaml_map_f(mapname,mapvalue,label,advance,unit,fmt)
    implicit none
    real(sp), intent(in) :: mapvalue
!    character(len=*), intent(in) :: mapname
!    character(len=*), optional, intent(in) :: label
!    character(len=*), optional, intent(in) :: advance
!    character(len=*), optional, intent(in) :: fmt
!    integer, optional, intent(in) :: unit
    !TODO: need to implement here
    include 'simple_yaml_map-incFl.f90'
!    include 'simple_yaml_open-incFl.f90'
  end subroutine yaml_map_f
  subroutine yaml_map_d(mapname,mapvalue,label,advance,unit,fmt)
    implicit none
    real(dp), intent(in) :: mapvalue
!    character(len=*), intent(in) :: mapname
!    character(len=*), optional, intent(in) :: label
!    character(len=*), optional, intent(in) :: advance
!    character(len=*), optional, intent(in) :: fmt
!    integer, optional, intent(in) :: unit
    !TODO: need to implement here
    include 'simple_yaml_map-incFl.f90'
!    include 'simple_yaml_open-incFl.f90'
  end subroutine yaml_map_d
  subroutine yaml_map_l(mapname,mapvalue,label,advance,unit,fmt)
    implicit none
    logical, intent(in) :: mapvalue
!    character(len=*), intent(in) :: mapname
!    character(len=*), optional, intent(in) :: label
!    character(len=*), optional, intent(in) :: advance
!    character(len=*), optional, intent(in) :: fmt
!    integer, optional, intent(in) :: unit
    !TODO: need to implement here
    include 'simple_yaml_map-incFl.f90'
    !    include 'simple_yaml_open-incFl.f90'
  end subroutine yaml_map_l

  ! SEQUENCERS METHODS

  subroutine yaml_close_stream(unit,istat)
    implicit none
    integer, optional, intent(in) :: unit
    integer, optional, intent(out) :: istat
    !local variables
    integer :: unt,istatus,strm,funt
    type(eglossary), pointer :: iter
    unt=DEFAULT_STREAM_ID
    if (present(unit)) unt=unit
    call get_stream(unt,strm,istat=istatus)

    !TODO: need to implement the new document method
    return
  end subroutine yaml_close_stream
  
  subroutine yaml_sequence(seqvalue,label,advance,unit,padding)
    implicit none
    character(len=*), optional, intent(in) :: seqvalue
    character(len=*), optional, intent(in) :: label
    character(len=*), optional, intent(in) :: advance
    integer, optional, intent(in) :: unit
    integer, intent(in), optional :: padding
    !local variables
    integer :: unt, istream
    integer :: msg_lgt,tb,istat,ipos,jpos,kpos
    character(len=3) :: adv
    character(len=tot_max_record_length) :: towrite

    unt=DEFAULT_STREAM_ID
    if (present(unit)) unt=unit
    call get_stream(unt,istream)

    adv='def' !default value
    if (present(advance)) adv=advance

    msg_lgt=0
    !put the optional name
    if (present(label)) then
       call buffer_string(towrite,len(towrite),' &',msg_lgt)
       call buffer_string(towrite,len(towrite),trim(label)//' ',msg_lgt)
    end if
    !put the value
    if (present(seqvalue)) &
         call buffer_string(towrite,len(towrite),trim(seqvalue),msg_lgt)

    if (present(padding)) then
       tb=padding-len_trim(seqvalue)
       if (tb > 0) call buffer_string(towrite,len(towrite),repeat(' ',tb),msg_lgt)
    end if
    !try to see if the line is too long
    call yaml_write(streams(istream),towrite(1:msg_lgt),advance=trim(adv),event=SEQUENCE_ELEM,istat=istat)
    if (istat /=0 .and. .not. streams(istream)%flowrite) then
       ipos=1
       jpos=msg_lgt
       kpos=jpos
       loop_seq: do
          call yaml_write(streams(istream),towrite(ipos:kpos),advance=trim(adv),event=SCALAR,istat=istat)
          if (istat /=0) then
             !continue searching
             kpos=index(towrite(ipos:jpos),' ',back=.true.)
             if (kpos == 0) kpos=jpos-1
          else
             ipos=kpos+1
             jpos=msg_lgt
             kpos=jpos
          end if
          if (ipos > msg_lgt) exit loop_seq
       end do loop_seq
    end if

    return
  end subroutine yaml_sequence

  !Open a level
  subroutine open_level(stream,doflow)
    implicit none
    type(yaml_stream), intent(inout) :: stream
    logical, intent(in) :: doflow
    stream%ilevel = stream%ilevel + 1
    if(doflow) then
       call open_flow_level(stream)
    else
       stream%ilast = stream%ilevel
    end if
  end subroutine open_level

  !Open a flow level (Indent more)
  subroutine open_flow_level(stream)
    implicit none
    type(yaml_stream), intent(inout) :: stream
    if (.not.stream%flowrite) then
       if (stream%iflowlevel==0) stream%indent_previous=stream%indent
       stream%indent=1
    end if
    stream%iflowlevel=stream%iflowlevel+1
    if (.not.stream%flowrite) stream%flowrite=.true. !start to write
  end subroutine open_flow_level

  !close a level
  subroutine close_level(stream,doflow)
    implicit none
    type(yaml_stream), intent(inout) :: stream
    logical, intent(in) :: doflow
    !local variables
    stream%ilevel = stream%ilevel - 1
    if(doflow) then
       call close_flow_level(stream)
    else
       stream%ilast = min(stream%ilevel,stream%ilast)
    end if
    return
  end subroutine close_level

  !close a flow level
  subroutine close_flow_level(stream)
    implicit none
    type(yaml_stream), intent(inout) :: stream
    !local variables
    stream%iflowlevel=stream%iflowlevel-1
    if (stream%iflowlevel==0) then
       stream%indent=stream%indent_previous
       stream%flowrite=.false.
       !reset the events in the flow
       stream%flow_events=NONE
       stream%ievt_flow=0
    else
       stream%indent=1
       stream%flowrite=.true.
    end if
    return
  end subroutine close_flow_level

  !sequence opener
  subroutine yaml_sequence_open(mapname,label,tag,flow,tabbing,advance,unit)
    implicit none
    character(len=*), optional, intent(in) :: mapname !< Key of the sequence. @copydoc doc::mapname
    character(len=*), optional, intent(in) :: label   !< @copydoc doc::label
    character(len=*), optional, intent(in) :: tag     !< @copydoc doc::tag
    logical, optional, intent(in) :: flow             !< @copydoc doc::flow
    character(len=*), optional, intent(in) :: advance !< @copydoc doc::advance
    integer, optional, intent(in) :: unit             !< @copydoc doc::unit
    integer, optional, intent(in) :: tabbing          !< @copydoc doc::tabbing
    include 'simple_yaml_open-incFl.f90'
    !TODO: need to implement the new document method
    !write(*,*) "in yaml_sequence_open: "
    !write(*,*) "mapname: ",mapname
    !write(*,*) "label: ",label
    !write(*,*) "tag: ",tag
    !write(*,*) "flow: ",flow
    !write(*,*) "tabbing: ",tabbing
    !write(*,*) "advance: ",advance
    !write(*,*) "unit: ",unit

    call yaml_write(streams(strm),towrite(1:msg_lgt),advance=trim(adv),&
                    event=SEQUENCE_START)
    return
  end subroutine yaml_sequence_open

  subroutine yaml_sequence_close(advance,unit)
    implicit none
    character(len=*), optional, intent(in) :: advance
    integer, optional, intent(in) :: unit
    !local variables
    integer :: unt,strm
    character(len=3) :: adv
    logical :: doflow

    unt=DEFAULT_STREAM_ID
    if (present(unit)) unt=unit
    call get_stream(unt,strm)

    if (streams(strm)%iflowlevel>1) then
       adv='no'
    else
       adv='yes'
       if (present(advance)) adv=advance
    end if

    call yaml_write(streams(strm),' ',advance=trim(adv),event=SEQUENCE_END)
    doflow = (streams(strm)%flowrite)
    call close_level(streams(strm),doflow)

    return
  end subroutine yaml_sequence_close


end module simple_yaml_output
