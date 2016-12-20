!*******************************************************************************
!     Author: Frederic D.R. Bonnet date: 4th of March 2013.
!
! Name:
! greeting_version - Various utilities for version and greeting.
!
! Description:
! greeting_version provides version numbering rouintes for the codes and initial
! greeting routines to make the user know which technology that has been used.
!*******************************************************************************
!
module greeting_version

  use simple_textHandler

  implicit none

  interface
  end interface

contains
!*******************************************************************************
! DESCRIPTION
! subroutine to say hello GPU using MATH as a check initial code to see if
! modules is responding
!
!*******************************************************************************
! SOURCE
  subroutine hello_cpu_discrepancy_checker()
    implicit none

    !start of the execution commands

!    write(*,*)
    write(*,*) "hello CPU world with CPU vs GPU discrepancy checker"
    write(*,*)

    return    
  end subroutine hello_cpu_discrepancy_checker
  !byeer
  subroutine bye_cpu_discrepancy_checker()
    implicit none

    !start of the execution commands

!    write(*,*)
    write(*,*) "bye CPU world with CPU vs GPU discrepancy checker"
    write(*,*)

    return    
  end subroutine bye_cpu_discrepancy_checker
!*******************************************************************************
! DESCRIPTION
! subroutine to say hello GPU using MATH as a check initial code to see if
! modules is responding
!
!*******************************************************************************
! SOURCE
  subroutine hello_gpu_math()
    implicit none

    !start of the execution commands

!    write(*,*)
    write(*,*) "hello GPU world with simple_math_gpu"
    write(*,*)

    return    
  end subroutine hello_gpu_math
  !byeer
  subroutine bye_gpu_math()
    implicit none

    !start of the execution commands

!    write(*,*)
    write(*,*) "bye GPU world with simple_math_gpu"
    write(*,*)

    return    
  end subroutine bye_gpu_math
!*******************************************************************************
! DESCRIPTION
! subroutine to say hello GPU using MAGMA as a check initial code to see if
! modules is responding
!
!*******************************************************************************
! SOURCE
  subroutine hello_gpu_magma()
    implicit none

    !start of the execution commands

!    write(*,*)
    write(*,*) "hello GPU world with MAGMA"
    write(*,*)

    return    
  end subroutine hello_gpu_magma
  !byeer
  subroutine bye_gpu_magma()
    implicit none

    !start of the execution commands

!    write(*,*)
    write(*,*) "bye GPU world with MAGMA"
    write(*,*)

    return    
  end subroutine bye_gpu_magma
!*******************************************************************************
! DESCRIPTION
! subroutine to stamp Simple environment.
!
!*******************************************************************************
! SOURCE
  subroutine simple_version()
    implicit none
    write(*,*) "* Simple version:                                           *"
    write(*,*) "* May 16th 2015 - version 2.1                               *"
    return
  end subroutine simple_version
  subroutine hello_simple()
    implicit none

    !start of the execution commands

    write(*,*)"#Logo for the application of Simple                  "
    write(*,*)"                                                     "
    write(*,*)"   ####      #    #    #  #####   #       ######     "
    write(*,*)"  #          #    ##  ##  #    #  #       #          "
    write(*,*)"   ####      #    # ## #  #    #  #       #####      "
    write(*,*)"       #     #    #    #  #####   #       #          "
    write(*,*)"  #    #     #    #    #  #       #       #          "
    write(*,*)"   ####      #    #    #  #       ######  ######     "
    write(*,*)"                                                     "
    write(*,*)"By Frederic Bonnet 4th of March 2015                 "
    write(*,*)

    return
  end subroutine hello_simple
  !byeer
  subroutine bye_simple()
    implicit none

    !start of the execution commands

!    write(*,*)
    write(*,*)"#Logo for the application of Simple                  "
    write(*,*)"                                                     "
    write(*,*)"   ####      #    #    #  #####   #       ######     "
    write(*,*)"  #          #    ##  ##  #    #  #       #          "
    write(*,*)"   ####      #    # ## #  #    #  #       #####      "
    write(*,*)"       #     #    #    #  #####   #       #          "
    write(*,*)"  #    #     #    #    #  #       #       #          "
    write(*,*)"   ####      #    #    #  #       ######  ######     "
    write(*,*)"                                                     "
    write(*,*)"By Frederic Bonnet 4th of March 2015                 "
    write(*,*)

    return
  end subroutine bye_simple
!*******************************************************************************
! DESCRIPTION
! subroutine to say hello GPU as a check initial code to see if modules is
! responding
!
!*******************************************************************************
! SOURCE
  subroutine hello_gpu()
    implicit none

    !start of the execution commands

!    write(*,*)
    write(*,*) "hello GPU world"
    write(*,*)

    return
  end subroutine hello_gpu

  subroutine bye_gpu()
    implicit none

    !start of the execution commands

!    write(*,*)
    write(*,*) "bye GPU world"
    write(*,*)

    return
  end subroutine bye_gpu
!*******************************************************************************
! DESCRIPTION
! subroutine to say hello CPU as a check initial code to see if modules is
! responding
!
!*******************************************************************************
! SOURCE
  subroutine hello_cpu()
    implicit none

    !start of the execution commands

!    write(*,*)
    write(*,*) "hello CPU world"
    write(*,*)

    return
  end subroutine hello_cpu

  subroutine bye_cpu()
    implicit none

    !start of the execution commands

!    write(*,*)
    write(*,*) "Bye CPU world"
    write(*,*)

    return
  end subroutine bye_cpu
!*******************************************************************************
! DESCRIPTION
! subroutine to say hello MPI as a check initial code to see if modules is
! responding
!
!*******************************************************************************
! SOURCE

  subroutine hello_mpi(my_rank,n_proc)
    implicit none
    !global variables
    integer                :: my_rank      !rank of process
    integer                :: n_proc       !number of processes

    !start of the execution commands

!    write(*,*)
    write(*,'(x,a,i5,x,a,i5,x,a)') "Hello MPI world, I am Rank: ",my_rank,"/",n_proc," of procs"
!    write(*,*)

    return
  end subroutine hello_mpi

  subroutine bye_mpi()
    implicit none

    !start of the execution commands

!    write(*,*)
    write(*,*) "bye MPI world"
    write(*,*)

    return
  end subroutine bye_mpi

!*******************************************************************************
! DESCRIPTION
! subroutine to create space and empty lines in the display!
!*******************************************************************************
! SOURCE
  subroutine spacer(nlines,nlength,in_line,type_lines)
    implicit none
    !global variables
    integer                     :: nlines,nlength
    character(len=*),intent(in) :: type_lines
    logical,intent(in)          :: in_line
    !local variables
    integer                     :: ilines,ilength
    integer,allocatable         :: integer_array(:)
    !start of the execution commands
    allocate(integer_array(1:nlength))

    select case (trim(lower_case(type_lines)))

    case ("lines")
       do ilines = 1,nlines
          write(*,*)
       end do
    case ("*")
       do ilines = 1,nlines
          do ilength = 1,nlength
             integer_array(ilength) = ichar("*")
          end do
          write(*,*)achar(integer_array(1:nlength))
       end do
    case ("=")
       do ilines = 1,nlines
          do ilength = 1,nlength
             integer_array(ilength) = ichar("=")
          end do
          write(*,*)achar(integer_array(1:nlength))
       end do

    case default

    end select !end select case (trim(lower_case(type_lines)))

    !freeing the ressources
    deallocate(integer_array)

    return
  end subroutine spacer

end module greeting_version
