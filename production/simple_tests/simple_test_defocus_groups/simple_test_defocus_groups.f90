program simple_test_defocus_groups
use simple_oris, only: oris
use simple_ori,  only: ori
implicit none
type(oris)           :: os_rnd_ctf, os
integer, parameter   :: group_pops(5) = [5,7,11,13,10]
integer, allocatable :: defgroups(:)
real,    allocatable :: ctfparams(:,:)
integer              :: cnt, i, j
type(ori)            :: o
call os_rnd_ctf%new(5)
call os_rnd_ctf%rnd_ctf(300., 2.7, 0.7, 3.0, 2.0, 1.0)
call os%new(sum(group_pops))
cnt = 0
do i=1,5
    o = os_rnd_ctf%get_ori(i)
    do j=1,group_pops(i)
        cnt = cnt + 1
        call os%set_ori(cnt, o)
    end do
end do
call os%write('before.txt')
call os%get_defocus_groups([1,cnt], defgroups, ctfparams)
call os%new(cnt)
do i=1,5
    print *, ctfparams(i,:)
end do
do i=1,cnt
    call os%set(i, 'kv',     ctfparams(defgroups(i),1))
    call os%set(i, 'cs',     ctfparams(defgroups(i),2))
    call os%set(i, 'fraca',  ctfparams(defgroups(i),3))
    call os%set(i, 'dfx',    ctfparams(defgroups(i),4))
    call os%set(i, 'dfy',    ctfparams(defgroups(i),5))
    call os%set(i, 'angast', ctfparams(defgroups(i),6))
end do
call os%write('after.txt')
end program simple_test_defocus_groups
