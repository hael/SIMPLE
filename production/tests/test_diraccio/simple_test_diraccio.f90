program simple_test_diraccio
use simple_fileio
implicit none
integer :: recsz, funit, io_stat, i
real    :: arr(5)
inquire(iolength=recsz) arr
arr(1:3) = 1.
arr(4:5) = 2.
call fopen(funit,'test_draccio.bin','replace','unknown',io_stat,'direct','unformatted',recl=recsz)
call fileiochk('fopen failed for test_draccio.bin',io_stat)
write(funit,rec=1) arr
write(funit,rec=3) arr
write(funit,rec=5) arr
call fclose(funit)

call fopen(funit,'test_draccio.bin','read','old',io_stat,'direct','unformatted',recl=recsz)
call fileiochk('fopen failed for test_draccio.bin',io_stat)
do i=1,5
    read(funit,rec=i) arr
    print *, i, arr
end do
call fclose(funit)
end program simple_test_diraccio
