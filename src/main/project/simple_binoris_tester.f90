!@descr: unit test routines for the binary orientation file (binoris)
! Header bookkeeping, segment round trips, in-place segment rewrites, legacy particle records, and the
! binoris_io / sp_project front doors.
module simple_binoris_tester
use simple_test_utils    ! assertions etc.
use simple_defs_ori,     only: N_PTCL_ORIPARAMS
use simple_type_defs,    only: MIC_SEG, STK_SEG, PTCL2D_SEG
use simple_string,       only: string
use simple_string_utils, only: int2str
use simple_syslib,       only: del_file, file_exists
use simple_fileio,       only: nlines
use simple_ori,          only: ori
use simple_oris,         only: oris
use simple_binoris,      only: binoris, binoris_seginfo
use simple_binoris_io,   only: binread_oritab, binread_ctfparams_state_eo, binread_nlines, binwrite_oritab
use simple_sp_project,   only: sp_project
implicit none
private
public :: run_all_binoris_tests

! the file header: MAX_N_SEGMENTS (20) segment records of five 8-byte integers (simple_binoris)
integer, parameter :: HEADER_NBYTES     = 20 * 5 * 8
integer, parameter :: PTCL_REC_NBYTES   = N_PTCL_ORIPARAMS * 4
integer, parameter :: NLEGACY_VALUES    = 40 ! a particle record width from before the last slot additions
integer, parameter :: NMICS = 3, NSTKS = 2, NPTCLS = 6
real,    parameter :: TOL = 1.0e-5
character(len=*), parameter :: BIN_FILE = 'tmp_binoris_test.simple'
character(len=*), parameter :: TXT_FILE = 'tmp_binoris_test.txt'
character(len=*), parameter :: PROJ_FILE  = 'tmp_binoris_project.simple'
character(len=*), parameter :: PROJ_FILE2 = 'tmp_binoris_project_missing.simple'

contains

    subroutine run_all_binoris_tests()
        write(*,'(A)') '**** running all binoris tests ****'
        call test_open_close_and_empty_header()
        call test_string_segment_roundtrip()
        call test_particle_segment_roundtrip()
        call test_legacy_narrow_particle_records()
        call test_write_segment_inside()
        call test_write_segment_inside_strings()
        call test_project_write_segment_inside()
        call test_binoris_io_text_dispatch()
        call test_binoris_io_project_dispatch()
        call cleanup()
    end subroutine run_all_binoris_tests

    !---------------- header bookkeeping ----------------

    subroutine test_open_close_and_empty_header()
        type(binoris) :: bos
        integer(kind=8) :: nbytes
        write(*,'(A)') 'test_open_close_and_empty_header'
        call del_file(string(BIN_FILE))
        call assert_false(bos%is_opened(),      'fresh binoris is not opened')
        call assert_int(0, bos%get_n_segments(), 'fresh binoris has no segments')
        call assert_true(bos%get_n_bytes_tot() == int(HEADER_NBYTES,kind=8), 'empty file is header only')
        call bos%open(string(BIN_FILE), del_if_exists=.true.)
        call assert_true(bos%is_opened(),       'open on a new name reports opened')
        call bos%write_header
        call bos%close
        call assert_false(bos%is_opened(),      'close reports closed')
        call assert_true(file_exists(string(BIN_FILE)), 'open + write_header creates the file')
        inquire(file=BIN_FILE, size=nbytes)
        call assert_true(nbytes == int(HEADER_NBYTES,kind=8), 'header-only file is '//int2str(HEADER_NBYTES)//' bytes')
        call bos%open(string(BIN_FILE))
        call assert_true(bos%is_opened(),        'reopen of a header-only file')
        call assert_int(0, bos%get_n_segments(), 'header-only file has no populated segments')
        call bos%close
    end subroutine test_open_close_and_empty_header

    !---------------- string-backed segment (mic) ----------------

    subroutine test_string_segment_roundtrip()
        type(binoris) :: bos
        type(oris)    :: os, os2
        type(ori)     :: o, o2
        type(string)  :: s, s2
        type(string),          allocatable :: sarr(:)
        type(binoris_seginfo), allocatable :: infos(:)
        integer,               allocatable :: seginds(:)
        integer(kind=8) :: nbytes
        integer :: i, nbpr
        write(*,'(A)') 'test_string_segment_roundtrip'
        call make_mic_oris(os, NMICS)
        nbpr = os%max_ori_strlen_trim()
        call bos%open(string(BIN_FILE), del_if_exists=.true.)
        call bos%write_segment(MIC_SEG, os)
        call bos%write_header
        call assert_int(int(MIC_SEG), bos%get_n_segments(),                 'n_segments after writing segment 1')
        call assert_int(NMICS, bos%get_n_records(MIC_SEG),                  'n_records of the mic segment')
        call assert_true(all(bos%get_fromto(MIC_SEG) == [1,NMICS]),        'fromto of the mic segment')
        call assert_int(nbpr, bos%get_n_bytes_per_record(MIC_SEG),          'record width is the longest trimmed record')
        call assert_true(bos%get_first_data_byte(MIC_SEG) == int(HEADER_NBYTES+1,kind=8), 'data starts right after the header')
        call assert_true(bos%get_n_bytes_tot() == int(HEADER_NBYTES + NMICS*nbpr,kind=8), 'total bytes: header + records')
        call bos%close
        inquire(file=BIN_FILE, size=nbytes)
        call assert_true(nbytes == int(HEADER_NBYTES + NMICS*nbpr,kind=8), 'file size equals the header total')
        ! reopen: header read back
        call bos%open(string(BIN_FILE))
        call assert_int(int(MIC_SEG), bos%get_n_segments(),         'reopen: n_segments')
        call assert_int(NMICS, bos%get_n_records(MIC_SEG),          'reopen: n_records')
        call assert_int(nbpr,  bos%get_n_bytes_per_record(MIC_SEG), 'reopen: record width')
        call bos%get_segments_info(seginds, infos)
        call assert_int(1, size(seginds),                           'segments_info: one populated segment')
        call assert_int(int(MIC_SEG), seginds(1),                   'segments_info: it is the mic segment')
        call assert_true(infos(1)%n_records == NMICS,               'segments_info: n_records')
        call assert_true(infos(1)%n_bytes_per_record == nbpr,       'segments_info: record width')
        call assert_true(infos(1)%first_data_byte == HEADER_NBYTES+1, 'segments_info: first data byte')
        ! records
        call os2%new(NMICS, is_ptcl=.false.)
        call bos%read_segment(MIC_SEG, os2)
        call assert_records_equal(os, os2, 1, NMICS, 'mic segment round trip')
        call bos%read_first_segment_record(MIC_SEG, o)
        call os%get_ori(1, o2)
        s  = o2%ori2str()
        s2 = o%ori2str()
        call assert_string_eq(s%to_char(), s2, 'read_first_segment_record returns record 1')
        ! raw string read (what merge_algndocs uses)
        allocate(sarr(NMICS))
        call bos%read_segment(MIC_SEG, sarr)
        do i = 1,NMICS
            s = os%ori2str(i)
            call assert_string_eq(s%to_char(), string(trim(sarr(i)%to_char())), 'string read of mic record '//int2str(i))
        enddo
        call bos%close
        call os%kill
        call os2%kill
        call o%kill
        call o2%kill
    end subroutine test_string_segment_roundtrip

    !---------------- fixed-width particle segment (ptcl2D) ----------------

    subroutine test_particle_segment_roundtrip()
        type(binoris) :: bos
        type(oris)    :: os, os2, os3
        integer, parameter :: FROMTO(2) = [2,5]
        integer :: i
        write(*,'(A)') 'test_particle_segment_roundtrip'
        call make_ptcl_oris(os, NPTCLS)
        call os%set(2, 'kv', 300.) ! not a fixed slot: particle records do not carry it
        call bos%open(string(BIN_FILE), del_if_exists=.true.)
        call bos%write_segment(PTCL2D_SEG, os, fromto=FROMTO)
        call bos%write_header
        call bos%close
        call bos%open(string(BIN_FILE))
        call assert_int(int(PTCL2D_SEG), bos%get_n_segments(),                  'n_segments after writing segment 3 only')
        call assert_int(0, bos%get_n_records(MIC_SEG),                          'segment 1 stays empty')
        call assert_int(FROMTO(2)-FROMTO(1)+1, bos%get_n_records(PTCL2D_SEG),   'n_records of the partial particle segment')
        call assert_true(all(bos%get_fromto(PTCL2D_SEG) == FROMTO),             'fromto of the partial particle segment')
        call assert_int(PTCL_REC_NBYTES, bos%get_n_bytes_per_record(PTCL2D_SEG), 'particle records are N_PTCL_ORIPARAMS reals')
        call assert_true(bos%get_first_data_byte(PTCL2D_SEG) == int(HEADER_NBYTES+1,kind=8), 'empty segments take no bytes')
        ! whole segment: records land at their absolute indices
        call os2%new(NPTCLS, is_ptcl=.true.)
        call bos%read_segment(PTCL2D_SEG, os2)
        do i = FROMTO(1),FROMTO(2)
            call assert_ptcl_equal(os, os2, i, 'particle round trip')
        enddo
        call assert_int(0, os2%get_state(1),      'record outside the written range is untouched (state)')
        call assert_real(0., os2%get(1,'corr'), TOL, 'record outside the written range is untouched (corr)')
        call assert_false(os2%isthere(2, 'kv'), 'non-slot key is not carried by a particle record')
        ! sub-range read
        call os3%new(NPTCLS, is_ptcl=.true.)
        call bos%read_segment(PTCL2D_SEG, os3, fromto=[3,4])
        call assert_ptcl_equal(os, os3, 3, 'sub-range read')
        call assert_ptcl_equal(os, os3, 4, 'sub-range read')
        call assert_int(0, os3%get_state(2), 'sub-range read leaves record 2 untouched')
        call assert_int(0, os3%get_state(5), 'sub-range read leaves record 5 untouched')
        call bos%close
        call os%kill
        call os2%kill
        call os3%kill
    end subroutine test_particle_segment_roundtrip

    ! records written by an older build with fewer slots read back with zeros in the new slots
    subroutine test_legacy_narrow_particle_records()
        type(binoris) :: bos
        type(oris)    :: os, os2
        type(ori)     :: o
        real    :: prec(N_PTCL_ORIPARAMS)
        integer :: i, funit, io_stat
        integer(kind=8) :: pos
        write(*,'(A)') 'test_legacy_narrow_particle_records'
        call make_ptcl_oris(os, NPTCLS)
        ! header claims NLEGACY_VALUES reals per record; the payload is written by hand
        call bos%open(string(BIN_FILE), del_if_exists=.true.)
        call bos%add_segment(PTCL2D_SEG, [1,NPTCLS], NLEGACY_VALUES*4)
        call bos%update_byte_ranges
        call bos%write_header
        pos = bos%get_first_data_byte(PTCL2D_SEG)
        call bos%close
        open(newunit=funit, file=BIN_FILE, access='stream', form='unformatted', action='readwrite', status='old', iostat=io_stat)
        call assert_int(0, io_stat, 'legacy file reopened for the raw payload')
        if( io_stat /= 0 ) return
        do i = 1,NPTCLS
            call os%get_ori(i, o)
            call o%ori2prec(prec)
            write(unit=funit, pos=pos) prec(1:NLEGACY_VALUES)
            pos = pos + NLEGACY_VALUES*4
        enddo
        close(funit)
        call bos%open(string(BIN_FILE))
        call assert_int(NLEGACY_VALUES*4, bos%get_n_bytes_per_record(PTCL2D_SEG), 'legacy record width read from the header')
        call os2%new(NPTCLS, is_ptcl=.true.)
        call bos%read_segment(PTCL2D_SEG, os2)
        call bos%close
        do i = 1,NPTCLS
            call assert_real(os%get(i,'corr'),   os2%get(i,'corr'),   TOL, 'legacy record: slot 3 (corr)')
            call assert_int(os%get_state(i),     os2%get_state(i),         'legacy record: slot 23 (state)')
            call assert_real(os%get(i,'lp_est'), os2%get(i,'lp_est'), TOL, 'legacy record: slot 40 (lp_est), the last one present')
            call assert_real(0., os2%get(i,'sampled'), TOL,               'legacy record: slot 45 (sampled) reads as zero')
            call assert_real(0., os2%get(i,'cluster'), TOL,               'legacy record: slot 46 (cluster) reads as zero')
        enddo
        call os%kill
        call os2%kill
        call o%kill
    end subroutine test_legacy_narrow_particle_records

    !---------------- in-place rewrite of one segment ----------------

    ! three segments; the middle one is rewritten longer, then shorter: the segments after it
    ! move but keep their bytes, the header total tracks the file size
    subroutine test_write_segment_inside()
        type(binoris) :: bos
        type(oris)    :: os_mic, os_stk, os_ptcl, os_stk_long, os_stk_short, rd_mic, rd_stk, rd_ptcl
        integer(kind=8) :: nbytes
        integer :: nbpr_mic, nbpr_stk, i
        write(*,'(A)') 'test_write_segment_inside'
        call make_mic_oris(os_mic, NMICS)
        call make_stk_oris(os_stk, NSTKS, 'stack')
        call make_ptcl_oris(os_ptcl, NPTCLS)
        call make_stk_oris(os_stk_long,  NSTKS+2, 'a_much_longer_stack_file_name')
        call make_stk_oris(os_stk_short, 1,       'stk')
        nbpr_mic = os_mic%max_ori_strlen_trim()
        call bos%open(string(BIN_FILE), del_if_exists=.true.)
        call bos%write_segment(MIC_SEG,    os_mic)
        call bos%write_segment(STK_SEG,    os_stk)
        call bos%write_segment(PTCL2D_SEG, os_ptcl)
        call bos%write_header
        call bos%close
        ! grow the middle segment
        call bos%open(string(BIN_FILE))
        call bos%write_segment_inside(STK_SEG, os_stk_long)
        call bos%close
        nbpr_stk = os_stk_long%max_ori_strlen_trim()
        call bos%open(string(BIN_FILE))
        call assert_int(NSTKS+2,  bos%get_n_records(STK_SEG),         'grown segment: n_records')
        call assert_int(nbpr_stk, bos%get_n_bytes_per_record(STK_SEG), 'grown segment: record width')
        call assert_true(bos%get_first_data_byte(PTCL2D_SEG) == int(HEADER_NBYTES + NMICS*nbpr_mic + (NSTKS+2)*nbpr_stk + 1,kind=8),&
            &'grown segment: the particle segment moved down')
        inquire(file=BIN_FILE, size=nbytes)
        call assert_true(nbytes == bos%get_n_bytes_tot(), 'grown segment: file size equals the header total')
        call read_three(bos, rd_mic, rd_stk, rd_ptcl, NSTKS+2)
        call assert_records_equal(os_mic,      rd_mic, 1, NMICS,   'grown segment: mic records intact')
        call assert_records_equal(os_stk_long, rd_stk, 1, NSTKS+2, 'grown segment: new stk records')
        do i = 1,NPTCLS
            call assert_ptcl_equal(os_ptcl, rd_ptcl, i, 'grown segment: particle records intact')
        enddo
        call bos%close
        ! shrink it
        call bos%open(string(BIN_FILE))
        call bos%write_segment_inside(STK_SEG, os_stk_short)
        call bos%close
        nbpr_stk = os_stk_short%max_ori_strlen_trim()
        call bos%open(string(BIN_FILE))
        call assert_int(1, bos%get_n_records(STK_SEG), 'shrunk segment: n_records')
        call assert_true(bos%get_first_data_byte(PTCL2D_SEG) == int(HEADER_NBYTES + NMICS*nbpr_mic + nbpr_stk + 1,kind=8),&
            &'shrunk segment: the particle segment moved up')
        inquire(file=BIN_FILE, size=nbytes)
        call assert_true(nbytes >= bos%get_n_bytes_tot(), 'shrunk segment: file holds at least the header total')
        call read_three(bos, rd_mic, rd_stk, rd_ptcl, 1)
        call assert_records_equal(os_mic,       rd_mic, 1, NMICS, 'shrunk segment: mic records intact')
        call assert_records_equal(os_stk_short, rd_stk, 1, 1,     'shrunk segment: new stk record')
        do i = 1,NPTCLS
            call assert_ptcl_equal(os_ptcl, rd_ptcl, i, 'shrunk segment: particle records intact')
        enddo
        call bos%close
        call os_mic%kill
        call os_stk%kill
        call os_ptcl%kill
        call os_stk_long%kill
        call os_stk_short%kill
        call rd_mic%kill
        call rd_stk%kill
        call rd_ptcl%kill
    end subroutine test_write_segment_inside

    ! the string-array variant (merge_algndocs): records arrive as strings with a caller-supplied width
    subroutine test_write_segment_inside_strings()
        type(binoris) :: bos
        type(oris)    :: os_mic, os_stk, os_ptcl, rd_mic, rd_stk, rd_ptcl
        type(string), allocatable :: recs(:)
        integer :: i, strlen_max
        write(*,'(A)') 'test_write_segment_inside_strings'
        call make_mic_oris(os_mic, NMICS)
        call make_stk_oris(os_stk, NSTKS, 'stack')
        call make_ptcl_oris(os_ptcl, NPTCLS)
        call bos%open(string(BIN_FILE), del_if_exists=.true.)
        call bos%write_segment(MIC_SEG,    os_mic)
        call bos%write_segment(STK_SEG,    os_stk)
        call bos%write_segment(PTCL2D_SEG, os_ptcl)
        call bos%write_header
        call bos%close
        ! new mic records as strings, one more than before and wider
        call os_mic%kill
        call make_mic_oris(os_mic, NMICS+1, prefix='renamed_micrograph')
        allocate(recs(NMICS+1))
        strlen_max = 0
        do i = 1,NMICS+1
            recs(i)    = os_mic%ori2str(i)
            strlen_max = max(strlen_max, recs(i)%strlen_trim())
        enddo
        call bos%open(string(BIN_FILE))
        call bos%write_segment_inside(MIC_SEG, recs, [1,NMICS+1], strlen_max)
        call bos%close
        call bos%open(string(BIN_FILE))
        call assert_int(NMICS+1,    bos%get_n_records(MIC_SEG),         'string rewrite: n_records')
        call assert_int(strlen_max, bos%get_n_bytes_per_record(MIC_SEG), 'string rewrite: record width as supplied')
        call read_three(bos, rd_mic, rd_stk, rd_ptcl, NSTKS, nmic_recs=NMICS+1)
        call assert_records_equal(os_mic, rd_mic, 1, NMICS+1, 'string rewrite: new mic records')
        call assert_records_equal(os_stk, rd_stk, 1, NSTKS,   'string rewrite: stk records intact')
        do i = 1,NPTCLS
            call assert_ptcl_equal(os_ptcl, rd_ptcl, i, 'string rewrite: particle records intact')
        enddo
        call bos%close
        call os_mic%kill
        call os_stk%kill
        call os_ptcl%kill
        call rd_mic%kill
        call rd_stk%kill
        call rd_ptcl%kill
    end subroutine test_write_segment_inside_strings

    !---------------- sp_project front door (the former inside_write test) ----------------

    subroutine test_project_write_segment_inside()
        type(sp_project) :: proj, before, after, fallback
        type(oris)       :: stk_new
        integer :: i
        write(*,'(A)') 'test_project_write_segment_inside'
        call del_file(string(PROJ_FILE))
        call del_file(string(PROJ_FILE2))
        call make_project(proj, string(PROJ_FILE))
        call proj%write(string(PROJ_FILE))
        call before%read(string(PROJ_FILE))
        call assert_int(NMICS,  before%os_mic%get_noris(),    'baseline project: mics')
        call assert_int(NSTKS,  before%os_stk%get_noris(),    'baseline project: stacks')
        call assert_int(NPTCLS, before%os_ptcl2D%get_noris(), 'baseline project: ptcl2D')
        call assert_int(NPTCLS, before%os_ptcl3D%get_noris(), 'baseline project: ptcl3D')
        call assert_int(1,      before%projinfo%get_noris(),  'baseline project: projinfo')
        ! rewrite only the stk segment, in place (stk_new is a twin of the replacement, for the comparison)
        call make_stk_oris(proj%os_stk, NSTKS+1, 'replacement_stack_with_a_longer_name')
        call make_stk_oris(stk_new,     NSTKS+1, 'replacement_stack_with_a_longer_name')
        call proj%write_segment_inside('stk', string(PROJ_FILE))
        call after%read(string(PROJ_FILE))
        call assert_int(NSTKS+1, after%os_stk%get_noris(), 'inside write: stk segment replaced')
        call assert_records_equal(stk_new, after%os_stk, 1, NSTKS+1, 'inside write: new stk records')
        call assert_records_equal(before%os_mic,    after%os_mic,    1, NMICS,  'inside write: mic segment untouched')
        call assert_records_equal(before%projinfo,  after%projinfo,  1, 1,      'inside write: projinfo untouched')
        do i = 1,NPTCLS
            call assert_ptcl_equal(before%os_ptcl2D, after%os_ptcl2D, i, 'inside write: ptcl2D untouched')
            call assert_ptcl_equal(before%os_ptcl3D, after%os_ptcl3D, i, 'inside write: ptcl3D untouched')
        enddo
        ! a missing project file falls back to a full write
        call proj%write_segment_inside('stk', string(PROJ_FILE2))
        call assert_true(file_exists(string(PROJ_FILE2)), 'inside write on a missing file creates it')
        call fallback%read(string(PROJ_FILE2))
        call assert_int(NMICS,   fallback%os_mic%get_noris(),    'fallback full write: mics')
        call assert_int(NSTKS+1, fallback%os_stk%get_noris(),    'fallback full write: stacks')
        call assert_int(NPTCLS,  fallback%os_ptcl2D%get_noris(), 'fallback full write: ptcl2D')
        call proj%kill
        call before%kill
        call after%kill
        call fallback%kill
        call stk_new%kill
    end subroutine test_project_write_segment_inside

    !---------------- binoris_io: format dispatch ----------------

    subroutine test_binoris_io_text_dispatch()
        type(sp_project) :: spproj
        type(oris)       :: os, os2, os3
        type(string)     :: s
        integer :: i
        write(*,'(A)') 'test_binoris_io_text_dispatch'
        call del_file(string(TXT_FILE))
        call make_mic_oris(os, NMICS)
        call binwrite_oritab(string(TXT_FILE), spproj, os, [1,NMICS])
        call assert_true(file_exists(string(TXT_FILE)), '.txt: binwrite_oritab writes the text table')
        call assert_int(NMICS, nlines(string(TXT_FILE)), '.txt: one line per record')
        call assert_int(NMICS, binread_nlines(string(TXT_FILE), int(MIC_SEG)), '.txt: binread_nlines counts lines')
        call os2%new(NMICS, is_ptcl=.false.)
        call binread_oritab(string(TXT_FILE), spproj, os2, [1,NMICS])
        call assert_records_equal(os, os2, 1, NMICS, '.txt: binread_oritab round trip')
        ! ctf parameters merge into an existing table, the other keys stay
        call os3%new(NMICS, is_ptcl=.false.)
        do i = 1,NMICS
            call os3%set(i, 'movie', 'other_'//int2str(i)//'.mrc')
            call os3%set(i, 'boxfile', 'boxes_'//int2str(i)//'.box')
            call os3%set(i, 'dfx', 9.)
            call os3%set(i, 'state', 0)
        enddo
        call binread_ctfparams_state_eo(string(TXT_FILE), spproj, os3)
        do i = 1,NMICS
            call assert_real(os%get(i,'dfx'),     os3%get(i,'dfx'),     TOL, '.txt: ctfparams merge updates dfx')
            call assert_real(os%get(i,'dfy'),     os3%get(i,'dfy'),     TOL, '.txt: ctfparams merge updates dfy')
            call assert_real(os%get(i,'angast'),  os3%get(i,'angast'),  TOL, '.txt: ctfparams merge updates angast')
            call assert_real(os%get(i,'kv'),      os3%get(i,'kv'),      TOL, '.txt: ctfparams merge updates kv')
            call assert_real(os%get(i,'phshift'), os3%get(i,'phshift'), TOL, '.txt: ctfparams merge updates phshift')
            call assert_int(os%get_state(i),      os3%get_state(i),          '.txt: ctfparams merge updates state')
            s = os3%get_str(i, 'boxfile')
            call assert_char('boxes_'//int2str(i)//'.box', s%to_char(), '.txt: ctfparams merge keeps keys not in the file')
            s = os3%get_str(i, 'movie')
            call assert_char('other_'//int2str(i)//'.mrc', s%to_char(), '.txt: ctfparams merge does not touch other char keys')
        enddo
        call os%kill
        call os2%kill
        call os3%kill
    end subroutine test_binoris_io_text_dispatch

    subroutine test_binoris_io_project_dispatch()
        type(sp_project) :: proj, rd, rd_ctf
        type(oris)       :: dummy
        integer :: i
        write(*,'(A)') 'test_binoris_io_project_dispatch'
        call del_file(string(PROJ_FILE))
        call make_project(proj, string(PROJ_FILE))
        ! one segment written through the front door
        call binwrite_oritab(string(PROJ_FILE), proj, dummy, [1,NMICS], isegment=MIC_SEG)
        call assert_true(file_exists(string(PROJ_FILE)), '.simple: binwrite_oritab writes the project file')
        call assert_int(NMICS, binread_nlines(string(PROJ_FILE), int(MIC_SEG)),  '.simple: binread_nlines is the segment record count')
        call assert_int(0,     binread_nlines(string(PROJ_FILE), int(STK_SEG)),  '.simple: a segment not written has no records')
        ! the whole project
        call proj%write(string(PROJ_FILE))
        call assert_int(NPTCLS, binread_nlines(string(PROJ_FILE), int(PTCL2D_SEG)), '.simple: particle segment record count')
        call binread_oritab(string(PROJ_FILE), rd, dummy, [1,NPTCLS])
        call assert_int(NMICS,  rd%os_mic%get_noris(),    '.simple: binread_oritab reads every segment (mics)')
        call assert_int(NSTKS,  rd%os_stk%get_noris(),    '.simple: binread_oritab reads every segment (stacks)')
        call assert_int(NPTCLS, rd%os_ptcl2D%get_noris(), '.simple: binread_oritab reads every segment (ptcl2D)')
        do i = 1,NSTKS
            call assert_real(proj%os_stk%get(i,'dfx'), rd%os_stk%get(i,'dfx'), TOL, '.simple: stk dfx')
            call assert_real(proj%os_stk%get(i,'box'), rd%os_stk%get(i,'box'), TOL, '.simple: stk box')
        enddo
        do i = 1,NPTCLS
            call assert_ptcl_equal(proj%os_ptcl2D, rd%os_ptcl2D, i, '.simple: ptcl2D through binread_oritab')
        enddo
        ! ctf/state/eo only
        call binread_ctfparams_state_eo(string(PROJ_FILE), rd_ctf, dummy)
        call assert_int(NSTKS, rd_ctf%os_stk%get_noris(), '.simple: ctfparams read sizes the stk segment')
        do i = 1,NSTKS
            call assert_real(proj%os_stk%get(i,'dfx'),    rd_ctf%os_stk%get(i,'dfx'),    TOL, '.simple: ctfparams read restores dfx')
            call assert_real(proj%os_stk%get(i,'dfy'),    rd_ctf%os_stk%get(i,'dfy'),    TOL, '.simple: ctfparams read restores dfy')
            call assert_real(proj%os_stk%get(i,'angast'), rd_ctf%os_stk%get(i,'angast'), TOL, '.simple: ctfparams read restores angast')
            call assert_real(proj%os_stk%get(i,'kv'),     rd_ctf%os_stk%get(i,'kv'),     TOL, '.simple: ctfparams read restores kv')
            call assert_real(proj%os_stk%get(i,'cs'),     rd_ctf%os_stk%get(i,'cs'),     TOL, '.simple: ctfparams read restores cs')
            call assert_real(proj%os_stk%get(i,'fraca'),  rd_ctf%os_stk%get(i,'fraca'),  TOL, '.simple: ctfparams read restores fraca')
            call assert_real(proj%os_stk%get(i,'smpd'),   rd_ctf%os_stk%get(i,'smpd'),   TOL, '.simple: ctfparams read restores smpd')
            call assert_int(proj%os_stk%get_state(i),     rd_ctf%os_stk%get_state(i),         '.simple: ctfparams read restores state')
        enddo
        call proj%kill
        call rd%kill
        call rd_ctf%kill
    end subroutine test_binoris_io_project_dispatch

    !---------------- fixtures ----------------

    ! hash-backed micrograph records with record lengths that differ
    subroutine make_mic_oris(os, n, prefix)
        type(oris),                 intent(inout) :: os
        integer,                    intent(in)    :: n
        character(len=*), optional, intent(in)    :: prefix
        character(len=:), allocatable :: pfx
        integer :: i
        pfx = 'movie'
        if( present(prefix) ) pfx = prefix
        call os%new(n, is_ptcl=.false.)
        do i = 1,n
            call os%set(i, 'movie',   pfx//'_'//int2str(i)//'.mrc')
            call os%set(i, 'intg',    pfx//'_'//int2str(i)//'_intg.mrc')
            call os%set(i, 'smpd',    1.1)
            call os%set(i, 'kv',      300.)
            call os%set(i, 'cs',      2.7)
            call os%set(i, 'fraca',   0.1)
            call os%set(i, 'dfx',     1.5 + 0.01*real(i))
            call os%set(i, 'dfy',     1.6 + 0.01*real(i))
            call os%set(i, 'angast',  10.*real(i))
            call os%set(i, 'phshift', 0.25*real(i))
            call os%set(i, 'state',   1)
            call os%set(i, 'nptcls',  100 + i)
            if( i == n ) call os%set(i, 'ctfres', 3.75) ! last record longer than the others
        enddo
    end subroutine make_mic_oris

    subroutine make_stk_oris(os, n, stkname)
        type(oris),       intent(inout) :: os
        integer,          intent(in)    :: n
        character(len=*), intent(in)    :: stkname
        integer :: i
        call os%new(n, is_ptcl=.false.)
        do i = 1,n
            call os%set(i, 'stk',     stkname//'_'//int2str(i)//'.mrcs')
            call os%set(i, 'ctf',     'yes')
            call os%set(i, 'smpd',    1.25)
            call os%set(i, 'kv',      200.)
            call os%set(i, 'cs',      1.4)
            call os%set(i, 'fraca',   0.07)
            call os%set(i, 'dfx',     2.0 + 0.1*real(i))
            call os%set(i, 'dfy',     2.1 + 0.1*real(i))
            call os%set(i, 'angast',  5.*real(i))
            call os%set(i, 'phshift', 0.)
            call os%set(i, 'box',     64 + 32*i)
            call os%set(i, 'nptcls',  NPTCLS)
            call os%set(i, 'fromp',   1)
            call os%set(i, 'top',     NPTCLS)
            call os%set(i, 'state',   1)
        enddo
    end subroutine make_stk_oris

    ! particle records: fixed slots only, deterministic values, phase shift inside [0,2pi)
    subroutine make_ptcl_oris(os, n)
        type(oris), intent(inout) :: os
        integer,    intent(in)    :: n
        integer :: i
        call os%new(n, is_ptcl=.true.)
        do i = 1,n
            call os%set_euler(i, [30.*real(i), 10.*real(i), 5.*real(i)])
            call os%set_shift(i, [0.5*real(i), -0.25*real(i)])
            call os%set_stkind(i, 1)
            call os%set(i, 'indstk',    i)
            call os%set_state(i, mod(i,2))
            call os%set_class(i, 1 + mod(i,3))
            call os%set_dfx(i, 1.50 + 0.01*real(i))
            call os%set_dfy(i, 1.60 + 0.01*real(i))
            call os%set(i, 'angast',    10.*real(i))
            call os%set(i, 'phshift',   0.1*real(i))
            call os%set(i, 'corr',      0.5 + 0.01*real(i))
            call os%set(i, 'eo',        mod(i,2))
            call os%set(i, 'lp_est',    8. + real(i))
            call os%set(i, 'sampled',   i)
            call os%set(i, 'cluster',   10 + i)
            call os%set(i, 'updatecnt', i + 10)
        enddo
    end subroutine make_ptcl_oris

    subroutine make_project(proj, projfile)
        type(sp_project), intent(inout) :: proj
        type(string),     intent(in)    :: projfile
        call proj%kill
        call make_mic_oris(proj%os_mic, NMICS)
        call make_stk_oris(proj%os_stk, NSTKS, 'stack')
        call make_ptcl_oris(proj%os_ptcl2D, NPTCLS)
        call make_ptcl_oris(proj%os_ptcl3D, NPTCLS)
        call proj%update_projinfo(projfile)
    end subroutine make_project

    !---------------- comparisons ----------------

    ! dummy names must not coincide with the module constants (NMICS, NSTKS): Fortran is case-insensitive
    ! and a dummy hides the host constant, so an absent optional would be read instead of it
    subroutine read_three(bos, rd_mic, rd_stk, rd_ptcl, nstk_recs, nmic_recs)
        type(binoris),     intent(inout) :: bos
        type(oris),        intent(inout) :: rd_mic, rd_stk, rd_ptcl
        integer,           intent(in)    :: nstk_recs
        integer, optional, intent(in)    :: nmic_recs
        integer :: nm
        nm = NMICS
        if( present(nmic_recs) ) nm = nmic_recs
        call rd_mic%new(nm,        is_ptcl=.false.)
        call rd_stk%new(nstk_recs, is_ptcl=.false.)
        call rd_ptcl%new(NPTCLS, is_ptcl=.true.)
        call bos%read_segment(MIC_SEG,    rd_mic)
        call bos%read_segment(STK_SEG,    rd_stk)
        call bos%read_segment(PTCL2D_SEG, rd_ptcl)
    end subroutine read_three

    ! hash-backed records: the text form is the record, so it must come back verbatim
    subroutine assert_records_equal(expected, actual, ifrom, ito, label)
        type(oris),       intent(in) :: expected, actual
        integer,          intent(in) :: ifrom, ito
        character(len=*), intent(in) :: label
        type(string) :: s
        integer :: i
        do i = ifrom,ito
            s = expected%ori2str(i)
            call assert_string_eq(s%to_char(), actual%ori2str(i), label//': record '//int2str(i))
        enddo
    end subroutine assert_records_equal

    ! particle records: slot by slot for the slots the fixtures set
    subroutine assert_ptcl_equal(expected, actual, i, label)
        type(oris),       intent(in) :: expected, actual
        integer,          intent(in) :: i
        character(len=*), intent(in) :: label
        real :: e1(3), e2(3), sh1(2), sh2(2)
        integer :: j
        e1  = expected%get_euler(i)
        e2  = actual%get_euler(i)
        sh1 = expected%get_2Dshift(i)
        sh2 = actual%get_2Dshift(i)
        do j = 1,3
            call assert_real(e1(j), e2(j), TOL, label//': Euler angle of record '//int2str(i))
        enddo
        call assert_real(sh1(1), sh2(1), TOL, label//': shift x of record '//int2str(i))
        call assert_real(sh1(2), sh2(2), TOL, label//': shift y of record '//int2str(i))
        call assert_int(expected%get_state(i), actual%get_state(i), label//': state of record '//int2str(i))
        call assert_int(expected%get_class(i), actual%get_class(i), label//': class of record '//int2str(i))
        call assert_real(expected%get_dfx(i),         actual%get_dfx(i),         TOL, label//': dfx of record '//int2str(i))
        call assert_real(expected%get_dfy(i),         actual%get_dfy(i),         TOL, label//': dfy of record '//int2str(i))
        call assert_real(expected%get(i,'angast'),    actual%get(i,'angast'),    TOL, label//': angast of record '//int2str(i))
        call assert_real(expected%get(i,'phshift'),   actual%get(i,'phshift'),   TOL, label//': phshift of record '//int2str(i))
        call assert_real(expected%get(i,'corr'),      actual%get(i,'corr'),      TOL, label//': corr of record '//int2str(i))
        call assert_real(expected%get(i,'eo'),        actual%get(i,'eo'),        TOL, label//': eo of record '//int2str(i))
        call assert_real(expected%get(i,'lp_est'),    actual%get(i,'lp_est'),    TOL, label//': lp_est of record '//int2str(i))
        call assert_real(expected%get(i,'sampled'),   actual%get(i,'sampled'),   TOL, label//': sampled of record '//int2str(i))
        call assert_real(expected%get(i,'cluster'),   actual%get(i,'cluster'),   TOL, label//': cluster of record '//int2str(i))
        call assert_real(expected%get(i,'updatecnt'), actual%get(i,'updatecnt'), TOL, label//': updatecnt of record '//int2str(i))
        call assert_real(expected%get(i,'indstk'),    actual%get(i,'indstk'),    TOL, label//': indstk of record '//int2str(i))
    end subroutine assert_ptcl_equal

    subroutine cleanup()
        call del_file(string(BIN_FILE))
        call del_file(string(TXT_FILE))
        call del_file(string(PROJ_FILE))
        call del_file(string(PROJ_FILE2))
    end subroutine cleanup

end module simple_binoris_tester
