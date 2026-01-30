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
use simple_core_module_api
use simple_ftexp_shsrch
implicit none
call seed_rnd
call test_ftexp_shsrch
call test_ftexp_shsrch2
end program simple_test_ft_expanded
