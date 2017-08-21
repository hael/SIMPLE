program simple_test_binoris
use simple_oris,       only: oris
use simple_map_reduce, only: split_nobjs_even
use simple_strings,    only: int2str
use simple_binoris,    only: binoris
implicit none
type(oris)           :: a, os_peak
integer, allocatable :: parts(:,:)
type(binoris)        :: bos
integer              :: ipart, iptcl, n_recs, n_peaks
call a%new(10)
call a%rnd_oris
call a%write('oris1_facit.txt')
parts = split_nobjs_even(10,5)
do ipart=1,5
    call bos%new(a, parts(ipart,:))
    call bos%open('bos'//int2str(ipart)//'.bin', del_if_exists=.true.)
    call bos%write_header
    do iptcl=parts(ipart,1),parts(ipart,2)
        call bos%write_record(iptcl, a)
    end do
    call bos%kill
end do
call a%kill
call a%new(10)
do ipart=1,5
    call bos%open('bos'//int2str(ipart)//'.bin')
    do iptcl=parts(ipart,1),parts(ipart,2)
        call bos%read_record(iptcl, a)
    end do
    call bos%kill
end do
call a%write('oris1.txt')
call a%new(10)
call a%rnd_oris
call a%write('oris2_facit.txt')
call os_peak%new(3)
call os_peak%rnd_oris()
call os_peak%rnd_states(2)
call os_peak%rnd_corrs()
call os_peak%set_all2single('ow', 0.5)
call os_peak%set_all2single('proj', 5.0)
call os_peak%write('os_peak_facit.txt')
do ipart=1,5
    call bos%new(a, parts(ipart,:), os_peak)
    call bos%open('bos'//int2str(ipart)//'.bin', del_if_exists=.true.)
    call bos%write_header
    do iptcl=parts(ipart,1),parts(ipart,2)
        call bos%write_record(iptcl, a, os_peak)
    end do
    call bos%kill
end do
call a%kill
call a%new(10)
do ipart=1,5
    call bos%open('bos'//int2str(ipart)//'.bin')
    do iptcl=parts(ipart,1),parts(ipart,2)
        call bos%read_record(iptcl, a, os_peak)
        call os_peak%write('os_peak'//int2str(ipart)//'.txt')
    end do
    call bos%kill
end do
call a%write('oris2.txt')
end program simple_test_binoris
