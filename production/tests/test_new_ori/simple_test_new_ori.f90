program simple_test_new_ori
use simple_oris, only: oris
use simple_ori,  only: ori
use simple_syscalls
use simple_image, only: image
implicit none
type(ori)   :: o
type(oris)  :: os
integer     :: i
!real        :: rval
!type(image) :: img
integer, parameter :: NORIS=20000  !, NTST=1000
o = ori()
call o%set('movie', 'mymovie.mrc')
call o%set('movie_ctf', 'mymovie_forctf.mrc')
call os%new(NORIS)
do i=1,NORIS
    call os%set_ori(i,o)
end do
call os%write('oris_test.txt')
call os%read('oris_test.txt')
call os%write('oris_test.txt')
end program simple_test_new_ori
