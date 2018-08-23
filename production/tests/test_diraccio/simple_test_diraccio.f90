program simple_test_diraccio
use simple_fileio
implicit none

type o_peak
    real    :: eul(3)   = 0.
    real    :: shift(2) = 0.
    real    :: corr     = 0.
    real    :: ow       = 0.
    integer :: iptcl    = 0
    integer :: proj     = 0
    integer :: ipeak    = 0
    integer :: state    = 0
end type o_peak

type(o_peak) :: o_peak2disc(3), o_peak2read(3)
integer :: recsz, funit, io_stat, i
real    :: arr(5)

o_peak2disc(1) = o_peak(eul=[1.,2.,3.], shift=[1.,2.], corr=1., ow=1., iptcl=1, proj=1, ipeak=1, state=1)
o_peak2disc(2) = o_peak(eul=[4.,5.,6.], shift=[2.,3.], corr=0.8, ow=0.8, iptcl=2, proj=2, ipeak=2, state=2)
o_peak2disc(3) = o_peak(eul=[5.,6.,7.], shift=[3.,4.], corr=0.6, ow=0.6, iptcl=3, proj=3, ipeak=3, state=3)
inquire(iolength=recsz) o_peak2disc
call fopen(funit,'test_draccio.bin','replace','unknown',io_stat,'direct','unformatted',recl=recsz)
call fileiochk('fopen failed for test_draccio.bin',io_stat)
write(funit,rec=1) o_peak2disc
write(funit,rec=3) o_peak2disc
write(funit,rec=5) o_peak2disc
call fclose(funit)

call fopen(funit,'test_draccio.bin','read','old',io_stat,'direct','unformatted',recl=recsz)
call fileiochk('fopen failed for test_draccio.bin',io_stat)
do i=1,5
    read(funit,rec=i) o_peak2read
    print *, o_peak2read
end do
call fclose(funit)

end program simple_test_diraccio
