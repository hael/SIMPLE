!@descr: unit tests for image file headers (simple_imghead): SPIDER header geometry, SPIDER and MRC header round trips, header probes
! Replaces the in-module self-test test_imghead (one box of 120, dimensions only, a file left behind).
! A SPIDER header is labrec records of lenbyt = 4*nx bytes, labrec = ceiling(1024/lenbyt), so
! labbyt = labrec*lenbyt is at least 1024 bytes and the data start after it. getLabbyt returned
! lenbyt (fixed 2026-09-25): headers were written short (nx words, so a box below 43 lost the pixel
! size) and read short, and a box below 43 read past the 43 header fields (the first build of the
! image tester stopped there, on a SPIDER volume of box 32).
module simple_imghead_tester
use simple_test_utils
use simple_defs
use simple_string,       only: string
use simple_string_utils, only: int2str
use simple_syslib,       only: del_file
use simple_imghead,      only: SpiImgHead, MrcImgHead, MRC_MODE_FLOAT32, find_ldim_nptcls, find_img_smpd
implicit none
private
public :: run_all_imghead_tests

character(len=*), parameter :: TMP_SPI = 'tmp_imghead_tester.spi'
character(len=*), parameter :: TMP_MRC = 'tmp_imghead_tester.mrc'
real,             parameter :: SMPD    = 1.25

contains

    subroutine run_all_imghead_tests()
        write(*,'(A)') '**** running all imghead tests ****'
        call test_spider_geometry()
        call test_spider_roundtrip()
        call test_spider_probe()
        call test_mrc_roundtrip()
        call del_file(TMP_SPI)
        call del_file(TMP_MRC)
    end subroutine run_all_imghead_tests

    ! boxes 32 and 120 need several records, 300 fits in one
    subroutine test_spider_geometry()
        integer, parameter :: NXS(3) = [32, 120, 300], LABRECS(3) = [8, 3, 1], LABBYTS(3) = [1024, 1440, 1200]
        type(SpiImgHead) :: hed
        integer :: i
        write(*,'(A)') 'test_spider_geometry'
        do i = 1,3
            call hed%new([NXS(i),NXS(i),1])
            call assert_int(4*NXS(i),     hed%getLenbyt(),     'SPIDER box '//int2str(NXS(i))//': a record is 4*nx bytes')
            call assert_int(LABRECS(i),   hed%getLabrec(),     'SPIDER box '//int2str(NXS(i))//': ceiling(1024/lenbyt) records')
            call assert_int(LABBYTS(i),   hed%getLabbyt(),     'SPIDER box '//int2str(NXS(i))//': getLabbyt is labrec*lenbyt')
            call assert_int(LABBYTS(i)+1, hed%firstDataByte(), 'SPIDER box '//int2str(NXS(i))//': the data follow the header')
            call hed%kill
        end do
    end subroutine test_spider_geometry

    ! write a header, check the file holds exactly the header, read it back into a fresh object
    subroutine test_spider_roundtrip()
        integer, parameter :: LDIMS(3,3) = reshape([32,32,1, 120,120,1, 32,32,32], [3,3])
        type(SpiImgHead) :: hed, back
        integer(kind=8)  :: fsize
        integer          :: i, funit, ios, maxim
        character(len=:), allocatable :: tag
        write(*,'(A)') 'test_spider_roundtrip'
        do i = 1,3
            tag   = 'SPIDER '//int2str(LDIMS(1,i))//'x'//int2str(LDIMS(2,i))//'x'//int2str(LDIMS(3,i))//': '
            maxim = merge(7, 1, LDIMS(3,i) == 1)
            call del_file(TMP_SPI)
            call hed%new(LDIMS(:,i))
            call hed%setPixSz(SMPD)
            call hed%setMaxim(maxim)
            open(newunit=funit, file=TMP_SPI, access='stream', form='unformatted', status='replace', action='write', iostat=ios)
            call assert_int(0, ios, tag//'the file opens for writing')
            if( ios /= 0 ) cycle
            call hed%write(funit)
            close(funit)
            inquire(file=TMP_SPI, size=fsize)
            call assert_true(fsize == int(hed%getLabbyt(), kind=8), tag//'write puts the whole header, labbyt bytes')
            call back%new(LDIMS(:,i))
            open(newunit=funit, file=TMP_SPI, access='stream', form='unformatted', status='old', action='read', iostat=ios)
            call assert_int(0, ios, tag//'the file opens for reading')
            if( ios /= 0 ) cycle
            call back%read(funit)
            close(funit)
            call assert_true(all(back%getDims() == LDIMS(:,i)),   tag//'the dimensions read back')
            call assert_real(SMPD, back%getPixSz(), 0.,            tag//'the pixel size reads back')
            call assert_int(merge(1, 3, LDIMS(3,i) == 1), back%getIform(), tag//'iform reads back (1 image, 3 volume)')
            call assert_int(maxim, back%getMaxim(),                tag//'maxim reads back')
            call assert_int(hed%getLabbyt(), back%getLabbyt(),     tag//'labbyt reads back')
            call hed%kill
            call back%kill
        end do
        call del_file(TMP_SPI)
    end subroutine test_spider_roundtrip

    ! a header is read by an object made for another box (as when a file is probed), and the
    ! probes of the image layer read the same fields
    subroutine test_spider_probe()
        type(SpiImgHead) :: hed, probe
        integer :: funit, ios, ldim(3), nptcls
        write(*,'(A)') 'test_spider_probe'
        call del_file(TMP_SPI)
        call hed%new([120,120,1])
        call hed%setPixSz(SMPD)
        call hed%setMaxim(7)
        open(newunit=funit, file=TMP_SPI, access='stream', form='unformatted', status='replace', action='readwrite', iostat=ios)
        call assert_int(0, ios, 'SPIDER probe: the file opens')
        if( ios /= 0 ) return
        call hed%write(funit)
        call probe%new([32,32,1])
        call probe%read(funit)
        close(funit)
        call assert_true(all(probe%getDims() == [120,120,1]), 'SPIDER probe: a box-32 header reads a box-120 file')
        call assert_int(1440, probe%getLabbyt(), 'SPIDER probe: and its header length')
        call find_ldim_nptcls(string(TMP_SPI), ldim, nptcls)
        call assert_true(all(ldim == [120,120,1]), 'find_ldim_nptcls: SPIDER dimensions')
        call assert_int(7, nptcls, 'find_ldim_nptcls: SPIDER stack size is maxim')
        call assert_real(SMPD, find_img_smpd(string(TMP_SPI)), 0., 'find_img_smpd: SPIDER pixel size')
        call hed%kill
        call probe%kill
        call del_file(TMP_SPI)
    end subroutine test_spider_probe

    ! an MRC header is 1024 bytes (no extended header) and reads back into an object made without dimensions
    subroutine test_mrc_roundtrip()
        type(MrcImgHead) :: hed, back
        integer(kind=8)  :: fsize
        integer :: funit, ios
        write(*,'(A)') 'test_mrc_roundtrip'
        call del_file(TMP_MRC)
        call hed%new([40,30,5])
        call hed%setPixSz(SMPD)
        open(newunit=funit, file=TMP_MRC, access='stream', form='unformatted', status='replace', action='write', iostat=ios)
        call assert_int(0, ios, 'MRC: the file opens for writing')
        if( ios /= 0 ) return
        call hed%write(funit)
        close(funit)
        inquire(file=TMP_MRC, size=fsize)
        call assert_true(fsize == 1024_8, 'MRC: the header is 1024 bytes')
        call back%new()
        open(newunit=funit, file=TMP_MRC, access='stream', form='unformatted', status='old', action='read', iostat=ios)
        call assert_int(0, ios, 'MRC: the file opens for reading')
        if( ios /= 0 ) return
        call back%read(funit)
        close(funit)
        call assert_true(all(back%getDims() == [40,30,5]), 'MRC: the dimensions read back')
        call assert_real(SMPD, back%getPixSz(), 1.e-6,      'MRC: the pixel size reads back')
        call assert_int(MRC_MODE_FLOAT32, back%getMode(),   'MRC: the default mode is float32')
        call hed%kill
        call back%kill
        call del_file(TMP_MRC)
    end subroutine test_mrc_roundtrip

end module simple_imghead_tester
