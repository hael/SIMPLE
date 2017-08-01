! ============================================================================
! Name        : simple_error_handling
! Author      : Frederic Bonnet
! Version     : 1.0
! Date        : 27th of June 2016
! Description : Module to define a glossary of fucntionalities, thisd includes
!             : the simple_error_handling module, the simple_yaml module and 
!             : the like.
! ============================================================================
!
module simple_error_handling
  !  use simple_eglossary
  use simple_eglossary_lowlev
  use simple_yaml_strings
  use simple_err_defs

  implicit none
  
  private 
  !control the error environment in simple_error_handling
  logical,public :: try_environment=.false.
  type(eglossary), pointer :: egloss_errors=>null()  !the global eglossary
  type(eglossary), pointer :: egloss_present_error=>null()

  !Stack of egloss_present_error for nested try (opne and close)
  type, private :: error_stack
     !egloss_present_error point to here.
     type(eglossary), pointer :: current => null()
     type(error_stack), pointer :: previous => null() !previous error
  end type error_stack

  !Stack of errors for try clause
  type(error_stack), pointer :: error_pipelines=>null() 
  !public functions
  public :: file_err_initialise, file_err_finalise
  public :: file_err_throw, file_err_define
  public :: file_err,file_err_open_try,file_err_close_try
  public :: file_err_clean
  !pipelines methods
  public :: error_pipelines_clean
  
contains

  !subrouitne to handle the err initialisation
  subroutine file_err_initialise()
    implicit none
    !TODO: need to implemetn the initialiser
    if (associated(error_pipelines)) then
       call error_pipelines_clean()
    end if
    allocate(error_pipelines)
    call egloss_init(error_pipelines%current)
    egloss_present_error=>error_pipelines%current
    call egloss_init(egloss_errors)
    return
  end subroutine file_err_initialise
  !cleaner for the error pipeline
  subroutine error_pipelines_clean()
    implicit none
    type(error_stack), pointer :: stack
    nullify(egloss_present_error)
    do while(associated(error_pipelines))
      call egloss_free(error_pipelines%current)
      stack=>error_pipelines%previous
      deallocate(error_pipelines)
      error_pipelines=>stack
    end do
  end subroutine error_pipelines_clean
  
  !finalising the hnadlers and free resources
  subroutine file_err_finalise
    implicit none
    !TODO: need to implement all the releasers
    call egloss_free(egloss_errors)
    return
  end subroutine file_err_finalise
  
  !cleaner
  subroutine file_err_clean()
    implicit none
    !TODO: implement this method
    nullify(egloss_present_error)
    call egloss_free(error_pipelines%current)
    call egloss_init(error_pipelines%current)
    egloss_present_error=>error_pipelines%current
    return
  end subroutine file_err_clean
  !Openners and closers
  !subroutine to try openning the file stream
  subroutine file_err_open_try()
    !use simple_eglossary    
    implicit none
    type(error_stack), pointer :: stack
    try_environment=.true.
    allocate(stack)
    stack%previous=>error_pipelines
    error_pipelines=>stack
    call egloss_init(error_pipelines%current)
    egloss_present_error=>error_pipelines%current
    return
  end subroutine file_err_open_try
  !subroutine to do a try close 
  subroutine file_err_close_try()
    implicit none
    type(error_stack), pointer :: stack
    if (associated(error_pipelines%previous)) then
       nullify(egloss_present_error)
       call egloss_free(error_pipelines%current)
       stack=>error_pipelines%previous
       deallocate(error_pipelines)
       error_pipelines=>stack
       egloss_present_error=>error_pipelines%current
       try_environment=associated(error_pipelines%previous)
    else
       call file_err_clean()
       try_environment=.false.
    end if
    return
  end subroutine file_err_close_try

  !error handlers

  !subroutine to bring up the arror
  recursive function file_err(condition,err_msg,err_id,err_name)
    use simple_yaml_strings
    implicit none
    logical, intent(in), optional :: condition !the condition which raise the error
    integer, intent(in), optional :: err_id !code of the error to be thrown
    !should already have been defined by file_err_define
    character(len=*), intent(in), optional :: err_name !error name
    character(len=*), intent(in), optional :: err_msg  !error message
    logical :: file_err
    !local variables
    character(len=max_field_length) :: message
    integer :: strlen
    
    if (present(condition)) then
       file_err=condition
    else
       file_err=.true.
    end if
    !TODO: need to finish the method here and include the callback functions
    if (file_err) then
       if (present(err_msg)) then
          message(1:strlen(message))=err_msg
       else
          message(1:strlen(message))='UNKOWN'
       end if

       if (present(err_id)) then
          call file_err_throw(message,err_id=err_id)
       else if (present(err_name)) then
          call file_err_throw(message,err_name=err_name)
       else
          call file_err_throw(message)
       end if

    end if
    
  end function file_err

  !handler to define the errors
  subroutine file_err_define(err_name,err_msg, err_id,err_action)
    implicit none
    character(len=*), intent(in) :: err_name               !name of the error
    character(len=*), intent(in) :: err_msg                !error message
    integer, intent(out) :: err_id                         !code of the error

    character(len=*), intent(in), optional :: err_action   !not really sure
    !local variable
    type(eglossary), pointer :: egloss_error

    !assure initialization of the library in case of misuse
    if (.not. associated(egloss_errors)) then !call f_err_initialize()
       write(*,*) "simple_error_handling library not initialized"
       write(*,*) "simple_file_lib_initialized should be called"
       stop
    end if
    
    err_id=ERR_GENERIC
    err_id=egloss_len(egloss_errors)

    !TODO: need to implement the rest of the method to define the errors
    call egloss_init(egloss_error)
    
    return
  end subroutine file_err_define

  !method to throw the gerror message
  subroutine file_err_throw(err_msg, err_id, err_name)
    use simple_yaml_strings
    implicit none
    integer, intent(in), optional :: err_id      !error code to be raised.
                                                 ! already defined f_err_define
    character(len=*), intent(in), optional :: err_name  ! error name
    character(len=*), intent(in), optional :: err_msg   ! error message
    !local variables
    integer :: new_errcode
!    if (.not.associated(egloss_present_error)) then
    !TODO: need to fix the pointer from the simple_eglosssary module
!    end if
    new_errcode=ERR_GENERIC
    !start of the excution cammands
    !TODO: insert the handler for the message output
    write(*,*) err_msg, err_id, err_name

    if (present(err_name)) then
       !TODO: need to implement the .index. operator
    end if

    
    return
  end subroutine file_err_throw

end module simple_error_handling
