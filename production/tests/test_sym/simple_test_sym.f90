program simple_test_sym
include 'simple_lib.f08'
use simple_sym,  only: sym, sym_tester
use simple_oris, only: oris
use simple_ori,  only: ori
use simple_timer
implicit none
type(sym)                     :: se
integer                       :: nsubgrps, i
character(len=:), allocatable :: pgrp_str
call se%new('d6')
nsubgrps = se%get_nsubgrp()
! print *, '# subgroups: ', nsubgrps
! do i=1,nsubgrps
!     print *, se%get_subgrp_descr(i)
! end do
call sym_tester('c1')
call sym_tester('c4')
call sym_tester('c5')
call sym_tester('d2')
call sym_tester('d7')
call sym_tester('t')
call sym_tester('o')
call sym_tester('i')
end program simple_test_sym
