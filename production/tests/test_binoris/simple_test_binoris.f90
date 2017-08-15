program simple_test_binoris
use simple_oris,       only: oris
use simple_map_reduce, only: split_nobjs_even
use simple_strings,    only: int2str
use simple_binoris,    only: binoris
implicit none
type(oris)           :: a
integer, allocatable :: parts(:,:)
type(binoris)        :: bos
integer              :: ipart, iptcl
call a%new(10)
call a%rnd_oris
call a%write('rndoris.txt')
parts = split_nobjs_even(10,5)
do ipart=1,5
    call bos%new(a, parts(ipart,:))
    call bos%open('bos'//int2str(ipart)//'.bin', del_if_exists=.true.)
    call bos%write_header
    do iptcl=parts(ipart,1),parts(ipart,2)
        call bos%write_record(a, iptcl)
    end do
    call bos%kill
end do
call a%kill
call a%new(10)
do ipart=1,5
    call bos%open('bos'//int2str(ipart)//'.bin')
    do iptcl=parts(ipart,1),parts(ipart,2)
        call bos%read_record(a, iptcl)
    end do
    call bos%kill
end do
call a%write('rndoris_read.txt')
end program simple_test_binoris
