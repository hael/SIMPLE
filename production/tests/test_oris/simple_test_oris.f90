program simple_test_oris
include 'simple_lib.f08'
use simple_oris
use simple_ori
implicit none
type(oris) :: o, o1, o2, o_subset, o1_subset
logical    :: test_passed
integer    :: inds(10)=[1,2,3,4,5,6,7,8,9,10]
#include "simple_local_flags.inc"
test_passed=.true.
print *,'>>> CREATE ORIS '
call o%new(100, is_ptcl=.false.)
o1 = oris(100, is_ptcl=.false.)
print *,'>>> REALLOCATE '
call o1%reallocate(200)
if( o1%get_noris() .ne. 200 ) test_passed = .false.
if( .not. test_passed ) THROW_HARD('reallocate oris failed!')
print *,'>>> EXTRACT SUBSET'
o_subset=o%extract_subset(1, 10)
o1_subset=o1%extract_subset(inds)
print *,'>>> ORIS ELEMENT EXISTS ', o%exists(1)
if(.not. o%exists(1)) test_passed=.false.
call o%rnd_oris(5.)
print *,'>>> ORIS WRITE'
call o%write(string('test_oris_rndoris.txt'))
print *,'>>> ORIS READ'
call o2%read(string('test_oris_rndoris.txt'))
print *,'>>> ORIS WRITE 2'
call o2%write(string('test_oris_rndoris_copy.txt'))
call o%kill()
call o1%kill()
call o_subset%kill()
call o1_subset%kill()
if( test_passed )then
   print *, '>>> TEST PASSED'
else
   THROW_HARD('>>> TEST FAILED')
endif
end program simple_test_oris
