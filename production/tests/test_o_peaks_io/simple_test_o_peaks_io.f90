program simple_test_o_peaks_io
use simple_o_peaks_io
use simple_oris, only: oris
implicit none
type(oris) :: o_peaks(3), even_projs, o_peaks_read(3)
integer    :: i, n_nozero
! generate o_peaks
call o_peaks(1)%new(3)
call o_peaks(2)%new(3)
call o_peaks(3)%new(3)
call even_projs%new(100)
call even_projs%spiral
call o_peaks(1)%rnd_oris_discrete_from(even_projs)
call o_peaks(2)%rnd_oris_discrete_from(even_projs)
call o_peaks(3)%rnd_oris_discrete_from(even_projs)
call o_peaks(1)%set_all2single('ow', 1.)
call o_peaks(2)%set_all2single('ow', 1.)
call o_peaks(3)%set_all2single('ow', 1.)
! write o_peaks
call open_o_peaks_io('mypeaks.bin')
do i=1,3
    call write_o_peak(o_peaks(i), [1,3], i)
end do
call close_o_peaks_io
! read o_peaks
call open_o_peaks_io('mypeaks.bin')
do i=1,3
    call read_o_peak(o_peaks_read(i), [1,3], i, n_nozero)
end do
call close_o_peaks_io
do i=1,3
    call o_peaks(i)%print_(i)
    call o_peaks_read(i)%print_(i)
    print *, '******************'
end do
end program simple_test_o_peaks_io
