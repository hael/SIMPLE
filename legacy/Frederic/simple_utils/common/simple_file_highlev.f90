! 
!*******************************************************************************
!     Author: Frederic D.R. Bonnet date: 28th of June 2016.
!
! Name:
! simple_file_highlev - High level routine to drive the utilities methods.
!
!*******************************************************************************
!
module simple_file_highlev
  implicit none

contains
  !subroutine to initialise the file lib errors
  subroutine simple_file_lib_errors()
    use simple_eglossary
    use simple_file_utils
    use simple_yaml_output
    use simple_dynamic_memory
    use simple_timing
    implicit none
    !initialising the errors for the lib
    call glossary_errors                     !DONE, more need to be added
    call file_utils_errors                   !DONE
    call yaml_output_errors                  !DONE
    call simple_timing_errors                !DONE
    !dynamic memory errors handlers
    call dynamic_memory_errors               !DONE
    return
  end subroutine simple_file_lib_errors

  !Subroutine to initialise the 
  subroutine simple_file_lib_initialise()
    use simple_eglossary
    use simple_error_handling
    use simple_timing
    use simple_dynamic_memory
    implicit none

    !file error initialisors
    call file_err_initialise                  !DONE
    call simple_file_lib_errors
    !dynamic memory initialisation
    call file_malloc_initialise
    !file timing initiliasors
    call simple_file_timing_initialise
    !TODO: need to initialise the categories

    return
  end subroutine simple_file_lib_initialise

end module simple_file_highlev
