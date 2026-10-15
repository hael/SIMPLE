!@descr: benchmark report files: one per family (and partition) in BENCH/, one record appended per report
module simple_bench_io
use simple_string,       only: string
use simple_string_utils, only: int2str_pad
use simple_fileio,       only: fopen
use simple_syslib,       only: dir_exists, file_exists, simple_mkdir
implicit none

public :: BENCH_DIR, bench_report_fname, open_bench_report
private

character(len=*), parameter :: BENCH_DIR = 'BENCH'

contains

    !> BENCH/<FAMILY>_BENCH.txt, or BENCH/<FAMILY>_BENCH_PARTppp.txt for the reports of one partition
    function bench_report_fname( family, part, numlen ) result( fname )
        character(len=*),  intent(in) :: family
        integer, optional, intent(in) :: part, numlen
        type(string) :: fname
        integer :: nlen
        if( present(part) )then
            nlen = 1
            if( present(numlen) ) nlen = max(1, numlen)
            fname = string(BENCH_DIR//'/'//trim(family)//'_BENCH_PART'//int2str_pad(part, nlen)//'.txt')
        else
            fname = string(BENCH_DIR//'/'//trim(family)//'_BENCH.txt')
        endif
    end function bench_report_fname

    !> Opens the family's report file for appending one record, creating BENCH/ when it is missing, and writes
    !! the record header: the family, and the iteration and partition when given. Each partition has a file of
    !! its own, so concurrent workers never append to the same file.
    subroutine open_bench_report( fnr, family, iter, part, numlen )
        integer,           intent(out) :: fnr
        character(len=*),  intent(in)  :: family
        integer, optional, intent(in)  :: iter, part, numlen
        type(string) :: fname
        logical      :: l_exists
        if( .not. dir_exists(BENCH_DIR) ) call simple_mkdir(BENCH_DIR, verbose=.false.)
        fname = bench_report_fname(family, part, numlen)
        l_exists = file_exists(fname)
        call fopen(fnr, FILE=fname, STATUS='UNKNOWN', ACTION='WRITE', POSITION='APPEND')
        if( l_exists ) write(fnr,'(a)') ''
        write(fnr,'(a)') '*** BENCHMARK RECORD ***'
        write(fnr,'(a,a)') 'record family    : ', trim(family)
        if( present(iter) ) write(fnr,'(a,i0)') 'record iteration : ', iter
        if( present(part) ) write(fnr,'(a,i0)') 'record partition : ', part
        call fname%kill
    end subroutine open_bench_report

end module simple_bench_io
