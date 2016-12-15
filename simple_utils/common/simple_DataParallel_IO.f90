! ============================================================================
! Name        : simple_DataParallel_IO
! Author      : Frederic Bonnet
! Version     : 1.0
! Date        : 19th of March 2016
! Description : Module to handle parallel IO
!             : files and the like for simple library.
!             : 
! ============================================================================
!
module simple_DataParallel_IO
  use simple_defs
  use simple_file_utils
  use simple_file_defs

  implicit none

  private

  public :: hello_DataParallel_IO
  public :: writeData
  
contains
  subroutine hello_DataParallel_IO()
    implicit none
    !TODO: says hello first routine
    write(*,*)"Hello from parrallel IO"
    return
  end subroutine hello_DataParallel_IO

  subroutine writeData(m,n,data,nnodes,hstD,fpioD)!,file,unit,status,position,action,binary)
    implicit none
    !data to be written to IO
    integer                                 :: m,n
    complex(sp)                             :: data(m,*)
    !system details
    integer,intent(in),optional             :: nnodes
    type(systemDetails),intent(in),optional :: hstD
    type(fileDetails_Parallel_IO),intent(in),optional :: fpioD
    !file details
    character(len=80)                       :: file
    integer                                 :: unit
    character(len=80)                       :: status
    character(len=80)                       :: position
    character(len=80)                       :: action
    logical                                 :: binary
    !local variables
    integer                                 :: linkfile
    !start of the execution commands

    !first open the file
    linkfile = 1
    if (present(fpioD)) then
       file = fpioD%file
       unit = fpioD%unit
       status = fpioD%status
       position = fpioD%position
       action = fpioD%action
       binary = fpioD%binary
       call file_Parrallel_open(nnodes, file, unit, status, position, action, &
                                binary, linkfile, hstD, fpioD)
    end if

    call file_Parrallel_write(m,n,data,linkfile,hstD)

    call file_Parrallel_close(file,linkfile,hstD)
    
    !closing the file
    call file_close(unit)
    
    return
  end subroutine writeData

  
end module simple_DataParallel_IO
