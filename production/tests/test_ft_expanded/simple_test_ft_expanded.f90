!------------------------------------------------------------------------------!
! SIMPLE v2.5          Elmlund & Elmlund Lab         simplecryoem.com          !
!------------------------------------------------------------------------------!
!> Test program for simple_ft_expanded module
!
! @author
! 
!
! DESCRIPTION: 
!> Brief description of moduletest.
!
! REVISION HISTORY:
! DD Mmm YYYY - Initial Version
! TODO_dd_mmm_yyyy - TODO_describe_appropriate_changes - TODO_name
!------------------------------------------------------------------------------
program simple_test_ft_expanded
use simple_ft_expanded_tester
use simple_defs
implicit none
integer(8) :: version(3)
print *, 'SIMPLE VERSION:', trim(adjustl(SIMPLE_LIB_VERSION)), " (", trim(adjustl(build_descr)), ")"
print *, 'SIMPLE_PATH: ', trim(adjustl(SIMPLE_PATH))
print *, 'COMPILER: ', trim(adjustl(COMPILER))
version(1) =  INT(COMPILER_VERSION,8)/INT(1000000,8)
version(2) =  (COMPILER_VERSION - version(1))/1000
write(*,'(A,I0,A,I0)') 'COMPILER VERSION: ', version(1), '.',version(2)
call exec_ft_expanded_test
end program simple_test_ft_expanded
