! ============================================================================
! Name        : simple_file_utils
! Author      : Frederic Bonnet
! Version     : 1.0
! Date        : 12th of February 2016
! Description : Module to handle the files and the like for simple library.
!             : 
! ============================================================================
!
module simple_file_utils
  use simple_err_defs
  use simple_eglossary
  use simple_error_handling
  use simple_file_defs
  use simple_yaml_strings

  implicit none

  private
  !serial handlers
  public :: file_utils_errors
  public :: file_time
  public :: file_unit,file_exists
  public :: file_open, file_close, file_write
  public :: file_recl
  !getters
  public :: file_get_free_unit
  !parralel IO handlers
  public :: file_Parrallel_open
  public :: file_Parrallel_write
  public :: file_Parrallel_close
  public :: file_Parrallel_read
  
contains
!*******************************************************************************
!    Serial Data IO
!
!*******************************************************************************
!
  subroutine file_utils_errors()
    implicit none
    !TODO: throw to yaml error handler
    !Opening
    call file_err_define('OPEN_FILE_SUCESS', 'File open has succeeded',&
         ERR_OPEN_FILE_SUCCESS, err_action='No action')
    call file_err_define('OPEN_FILE_FAIL', 'File open has fail',&
         ERR_OPEN_FILE_FAIL, err_action='stop')
    !close
    call file_err_define('CLOSE_FILE_SUCCESS', 'File close has succeed',&
         ERR_CLOSE_FILE_SUCCESS, err_action='No action')
    call file_err_define('CLOSE_FILE_SUCCESS', 'File close has fail',&
         ERR_CLOSE_FILE_FAIL, err_action='stop')
    !Inquire
    call file_err_define('INQUIRE_FILE_SUCESS', 'File inquiry has succeeded',&
         ERR_INQUIRE_SUCCESS, err_action='No action')
    call file_err_define('INQUIRE_FILE_FAIL', 'File inquiry has fail',&
         ERR_INQUIRE_FAIL, err_action='stop')
    
    return
  end subroutine file_utils_errors
  
  subroutine file_time()
    implicit none
    !TODO: times the time taken for IO using the simple_timing.c for accuracy
    return
  end subroutine file_time

!  subroutine file_err_throw()
!    implicit none
!
!    return
!  end subroutine file_err_throw
  !subroutine to get the record length
  subroutine file_recl(unit, recl_max,recl)
    implicit none
    integer, intent(in)  :: unit !unit to check the record length upon
    integer, intent(in)  :: recl_max!max val for record length
    integer, intent(out) :: recl !return value
    !local variables
    integer :: err_recl
    integer :: err,unt
    integer(kind=recl_kind) :: recl_file
    logical :: unit_is_open
    !start of the execution commands
    recl=recl_max
    err_recl=-1
    recl_file=int(-1234567891,kind=recl_kind)
    unt = unit
    inquire(unit=unt,opened=unit_is_open,iostat=err)
    if (err==0 .and. .not. unit_is_open) then
       inquire(unit=unt,recl=recl_file,iostat=err_recl)
    end if
    if (err_recl==0) then
       recl=int(min(int(recl_max,kind=recl_kind),recl_file))
    end if
    if (recl<=0) recl=recl_max
    return
  end subroutine file_recl

  !subroutine to 
  
  subroutine file_exists(file,exists)
    implicit none
    character(len=*), intent(in) :: file
    logical, intent(out) :: exists
    !local variable
    integer :: err,rc
    character(len=1) :: char_out
    !function calls
    integer :: convert_int2char_pos_c
    exists = .false.
    inquire(file=trim(file), exist=exists, iostat=err)
    if ( err /= 0 ) then
       !TODO: insert the yaml to a something to yaml file
       rc = convert_int2char_pos_c(char_out,abs(err))
       call file_err_throw('Error in inquiring file='//&
            trim(file)//' for number, iostat='//"-"//trim(char_out),&
            err_id=INPUT_OUTPUT_ERROR)
    end if
    exists = exists .and. err ==0
    return
  end subroutine file_exists

  !gets a free unit
  function file_get_free_unit(unit) result(unit_out)
    implicit none
    integer, intent(in),optional :: unit
    integer :: unit_out
    !integer :: file_get_free_unit
    !local variables
    logical :: unit_is_open
    integer :: unt, err, rc
    character(len=1) :: char_out_err !TODO: need to check that the value len ok
    character(len=1) :: char_out_unt !TODO: need to check that the value len ok
    !function calls
    integer :: convert_int2char_pos_c
    
    unit_is_open=.true.
    unt=7
    if (present(unit)) unt=unit
    inquire(unit=unt,opened=unit_is_open,iostat=err)
    do while(unit_is_open .and. err==0)
       unt=unt+1
       inquire(unit=unt,opened=unit_is_open,iostat=err)
    end do
    if (err /=0) then
       !TODO: need to implement the yaml output for the throw
       rc = convert_int2char_pos_c(char_out_err,abs(err))
       rc = convert_int2char_pos_c(char_out_unt,abs(unt))
       call file_err_throw('Error in inquiring unit='//&
            trim(char_out_unt)//', iostat='//trim(char_out_err),&
            err_id=INPUT_OUTPUT_ERROR)
    end if
    unit_out = unt    
  end function file_get_free_unit

  !gets a file unit
  subroutine file_unit(file,unit)
    implicit none
    character(len=*), intent(in) :: file
    integer, intent(out) :: unit
    !local varaiables
    logical :: exists
    integer :: err, unt , rc
    character(len=1) :: char_out
    !function calls
    integer :: convert_int2char_pos_c
    
    unit = -1
    call file_exists(file,exists)
    if (exists) then
       inquire(file=trim(file), number=unt, iostat=err)
       if (err/=0) then
          !TODO: insert the yaml to a something to yaml file
          rc = convert_int2char_pos_c(char_out,abs(err))
          call file_err_throw('Error in inquiring file='//&
               trim(file)//' for number, iostat='//"-"//trim(char_out),&
               err_id=INPUT_OUTPUT_ERROR)
       else
          unit=unt
       end if
    end if

    return
  end subroutine file_unit
  
  subroutine file_open(file,unit,status,position,action,binary)
    use simple_yaml_strings, only: file_strcpy 
    implicit none
    
    character(len=*), intent(in) :: file
    integer, intent(inout) :: unit
    character(len=*), intent(in), optional :: status
    character(len=*), intent(in), optional :: position
    character(len=*), intent(in), optional :: action
    logical, intent(in), optional :: binary
    !local varaibles
    character(len=7) :: f_status
    character(len=11):: f_form
    character(len=6) :: f_position
    character(len=9) :: f_action
    integer :: unt,err,rc
    character(len=1) :: char_out
    !function calls
    integer :: convert_int2char_pos_c

    !check if the file is already open
    call file_unit(file,unt)
    
    !TODO: need to check the logic here fo file_nit
    if (unt /= -1 ) then
       unit = unt
    else
       unt = unit

       !useful open specifiers
       call file_strcpy(src='unknown',dest=f_status)
       if (present(status)) call file_strcpy(src=status,dest=f_status)

       call file_strcpy(src='formatted',dest=f_form)
       if (present(binary)) then
          if (binary) call file_strcpy(src='unformatted',dest=f_form)
       end if

       call file_strcpy(src='asis',dest=f_position)
       if (present(position)) call file_strcpy(src=position,dest=f_position)

       call file_strcpy(src='readwrite',dest=f_action)
       if (present(action)) call file_strcpy(src=action,dest=f_action)

       write(*,*) unt," ",trim(file)," ", f_status," ", f_form," ",f_position," ",f_action
       
       !then open the file with the given unit
       open(unit=unt,file=trim(file),status=f_status,form=f_form,&
            position=f_position,action=f_action,iostat=err)
       if (err /= 0) then
          !write(*,*) "files cannnot be open!!"
          !TODO: throw an yaml output error message to yaml file
          rc = convert_int2char_pos_c(char_out,abs(err))
          call file_err_throw('Error in opening file='//&
               trim(file)//' for number, iostat='//"-"//trim(char_out),&
               err_id=FILE_OPENING_ERROR)
       else
          unit=unt
       end if

    end if

    return
  end subroutine file_open
  
  ! subroutine to close the files for serial IO

  subroutine file_close(unit)
    implicit none
    integer, intent(in) :: unit
    !local variables
    integer :: err,rc
    character(len=1) :: char_out
    character(len=4) :: unit_out
    !function calls
    integer :: convert_int2char_pos_c
    
    if (unit>0) then
       close(unit,iostat=err)
       if (err/=0)  then
          !TODO: throw an yaml output error message to yaml file
          rc = convert_int2char_pos_c(char_out,abs(err))
          rc = convert_int2char_pos_c(unit_out,abs(unit))
          call file_err_throw('Error in closing file unit='//&
               trim(unit_out)//' for number, iostat='//"-"//trim(char_out),&
               err_id=ERR_CLOSE_FILE_FAIL)
       end if
    end if
    
    return
  end subroutine file_close
  
  subroutine file_write()
    implicit none
    !TODO: writes to file
    return
  end subroutine file_write
!*******************************************************************************
!    Parrallel Data IO
!
!*******************************************************************************
!
!*******************************************************************************
! DESCRIPTION
! subroutine to open the files for parrallel IO
!
!*******************************************************************************
! SOURCE
  subroutine file_Parrallel_open(nnodes,file,unit,status,position,action,bin,&
                                 filenumbr,hstD,fpioD)
    use simple_yaml_strings, only: file_strcpy
    use simple_file_defs
    use simple_cuda_defs
    implicit none
    integer                                           :: nnodes
    character(len=*), intent(in)                      :: file
    integer, intent(inout)                            :: unit
    character(len=*), intent(in), optional            :: status
    character(len=*), intent(in), optional            :: position
    character(len=*), intent(in), optional            :: action
    logical, intent(in), optional                     :: bin
    integer,optional,intent(in)                       :: filenumbr
    type(systemDetails),intent(in),optional           :: hstD
    type(fileDetails_Parallel_IO),intent(in),optional :: fpioD
    !local varaibles
    integer            :: unt,err
    character(len=7)   :: f_status
    character(len=11)  :: f_form
    character(len=6)   :: f_position
    character(len=9)   :: f_action
    !local parrallel variables
    integer            :: myfile, procnum
    character(len=120) :: myfilename, hostfilename
    character(len=150) :: systemcall
    character(len=50)  :: nodename
    character(len=5)   :: filesuffix
    integer            :: rc = RC_SUCCESS
    !check if the file is already open
    call file_unit(file,unt)

    write(*,*) file,unit,unt

    !TODO: need to check the logic here fo file_nit
    if (unt == -1 ) then
       unit = unt
    else

       write(*,*) file,unit,unt
       write(*,*) fpioD%file,fpioD%unit,unt

       !useful open specifiers
       call file_strcpy(src='unknown',dest=f_status)
       if (present(status)) call file_strcpy(src=status,dest=f_status)

       call file_strcpy(src='formatted',dest=f_form)
       if (present(bin)) then
          if (bin) call file_strcpy(src='unformatted',dest=f_form)
       end if

       call file_strcpy(src='asis',dest=f_position)
       if (present(position)) call file_strcpy(src=position,dest=f_position)

       call file_strcpy(src='readwrite',dest=f_action)
       if (present(action)) call file_strcpy(src=action,dest=f_action)

       !then open the file with the given unit
       unit=unt

    end if
    !if filenumber is present then we proceeding parrallel opens
    if (present(filenumbr)) then
       write(*,*) "filenumbr is pesent, filenumbr: ",filenumbr
       !TODO: outside the pseudo code need to replace with proper node finder
       if ( present(hstD)) then

          procnum = hstD%nCPUcores ! TODO: find the local processor number
          myfile = 1000 + procnum + filenumbr*MAXPROCS
          write(*,*) "myfile: ",myfile

          !Making a unique file name for scratch space for each processors to
          !be opened writting the suffix. Each files are written to local disk
          !on each nodes.
          write(filesuffix,fmt='(a2,i3.3)')'.p',procnum
          myfilename = '/scratch/' // trim(adjustl(file)) // filesuffix 
          open(myfile,file=myfilename,status='replace',form='unformatted',action='write')

          hostfilename = trim(adjustl(file)) // '.node' // filesuffix
          write(*,*) hostfilename
          systemcall = 'hostname > ' // hostfilename
          write(*,*) systemcall
          rc = system(systemcall)
          write(*,*) rc
          if (rc .ne. RC_SUCCESS) then
             write(*,*) 'ERROR: failed to output node name for processor ', procnum
             rc = RC_FAIL
          end if
          
       end if

       
    end if

    return
  end subroutine file_Parrallel_open
!*******************************************************************************
! DESCRIPTION
! subroutine to write data files in parrallel IO
!
!*******************************************************************************
! SOURCE
  subroutine file_Parrallel_write(m,n,data,filenumbr,hstD)
    use simple_defs
    use simple_file_defs
    implicit none
    !data to be written to IO
    integer                                 :: m,n
    complex(sp)                             :: data(m,*)
    integer,optional,intent(in)             :: filenumbr
    type(systemDetails),intent(in),optional :: hstD
    !local variables
    !local parrallel variables
    integer            :: myfile, procnum

    !start of the execution commands
    if (present(filenumbr)) then
       write(*,*) "file_Parrallel_write filenumbr: ",filenumbr
       procnum = hstD%nCPUcores ! find the local processor number
       myfile = 1000 + procnum + filenumbr*MAXPROCS
       write(*,*) "file_Parrallel_write myfile: ",myfile

       write(myfile) data(1:m,1:n)

    end if

    return
  end subroutine file_Parrallel_write
!*******************************************************************************
! DESCRIPTION
! subroutine to write data files in parrallel IO
!
!*******************************************************************************
! SOURCE
  subroutine file_Parrallel_read(m,n,data,filenumbr,hstD)
    use simple_defs
    use simple_file_defs
    implicit none
    !data to be written to IO
    integer                                 :: m,n
    complex(sp)                             :: data(m,*)
    integer,optional,intent(in)             :: filenumbr
    type(systemDetails),intent(in),optional :: hstD
    !local variables
    !local parrallel variables
    integer            :: myfile, procnum

    !start of the execution commands
    if (present(filenumbr)) then
       write(*,*) "file_Parrallel_write filenumbr: ",filenumbr
       procnum = hstD%nCPUcores ! find the local processor number
       myfile = 1000 + procnum + filenumbr*MAXPROCS
       write(*,*) "file_Parrallel_write myfile: ",myfile

       read(myfile) data(1:m,1:n)

    end if

    return
  end subroutine file_Parrallel_read
!*******************************************************************************
! DESCRIPTION
! subroutine to close data files in parrallel IO
!
!*******************************************************************************
! SOURCE
  subroutine file_Parrallel_close(file,filenumbr,hstD)
    use simple_defs
    use simple_cuda_defs
    use simple_file_defs
    implicit none
    !data to be written to IO
    character(len=*), intent(in)                      :: file
    integer,optional,intent(in)             :: filenumbr
    type(systemDetails),intent(in),optional :: hstD
    !local variables
    !local parrallel variables
    integer            :: myfile, procnum
    character(len=120) :: myfilename, hostfilename
    character(len=150) :: systemcall
    character(len=5)   :: filesuffix
    integer            :: rc = RC_SUCCESS

    !start of the execution commands
    if (present(filenumbr)) then
       write(*,*) "file_Parrallel_write filenumbr: ",filenumbr
       procnum = hstD%nCPUcores ! find the local processor number
       myfile = 1000 + procnum + filenumbr*MAXPROCS
       write(*,*) "file_Parrallel_write myfile: ",myfile

       close(myfile)

       !check the size of the files
       write(filesuffix,fmt='(a2,i3.3)')'.p',procnum
       myfilename = '/scratch/' // trim(adjustl(file)) // filesuffix 
       hostfilename = trim(adjustl(file)) // '.size' // filesuffix
       systemcall = '/bin/ls -l '// trim(adjustl(myfilename)) // ' > ' // hostfilename
       rc = system(systemcall)
       if (rc .ne. RC_SUCCESS) then
          write(*,*) 'ERROR: failed to output node name for processor ', procnum
          rc = RC_FAIL
      end if
      
    end if

    return
  end subroutine file_Parrallel_close
    
  
end module simple_file_utils
