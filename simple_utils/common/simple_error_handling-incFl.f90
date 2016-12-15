! ============================================================================
! Name        : simple_error_handling-inc
! Author      : Frederic Bonnet
! Version     : 1.0
! Date        : 08th of July 2016
! Description : File that inlcudes convoluted subroutines between
!             : simple_eglossary and simple_error_handling modules.
! ============================================================================
!
  !function to get the last error
  function file_get_last_error(add_msg)
    implicit none
    character(len=*), intent(out), optional :: add_msg
    integer :: file_get_last_error
    !local variables
    integer :: err
    !TODO: implement the last error
    err=egloss_len(egloss_present_error)-1
    if (err>=0) then
       file_get_last_error=egloss_present_error//err//ERRID
       if (present(add_msg)) add_msg=egloss_present_error//err//ERR_ADD_INFO
    else
       file_get_last_error=0
       if (present(add_msg)) add_msg=repeat(' ',len(add_msg))
    end if

  end function file_get_last_error
