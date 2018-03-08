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
include 'simple_lib.f08'
use simple_ft_expanded_tester
implicit none
call exec_ft_expanded_test
end program simple_test_ft_expanded
