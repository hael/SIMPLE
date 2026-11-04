!@descr: dose manifest round trips, interval selection, recoverable metadata errors and actual MRC header validation
module simple_dose_series_manifest_tester
use simple_defs,                 only: dp
use simple_string,               only: string
use simple_image,                only: image
use simple_dose_series_manifest, only: dose_series_manifest
use simple_syslib,               only: simple_mkdir, simple_rename, file_exists, del_file
use simple_test_utils,           only: assert_true, assert_false, assert_int, enter_fixture, leave_fixture, tests_failed
use, intrinsic :: iso_fortran_env, only: int8
use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan
implicit none
private
public :: run_all_dose_series_manifest_tests
integer, parameter :: BOX = 8
real(dp), parameter :: SMPD = 1.5_dp, DPF = 1.25_dp

contains

    !> Run hermetic manifest checks in the reconstruction area.
    subroutine run_all_dose_series_manifest_tests()
        type(string) :: saved, root
        integer :: failures
        failures = tests_failed
        call enter_fixture('dose_series_manifest', saved, root)
        call make_maps()
        call test_roundtrip()
        call test_bad_metadata()
        call test_bad_headers()
        call test_publication()
        call test_final_line()
        call leave_fixture(saved, root, failures)
    end subroutine run_all_dose_series_manifest_tests

    !> Write six distinct half maps on one known grid; the manifest carries their dose coordinates.
    subroutine make_maps()
        class(image), allocatable :: vol
        integer :: i
        type(string) :: name
        allocate(vol)
        call simple_mkdir('data')
        call simple_mkdir('data/volume_dir')
        call vol%new([BOX,BOX,BOX], real(SMPD))
        do i = 1,6
            vol = real(i)
            name = string('data/half')//string(i)//'.mrc'
            call vol%write(name, wfloat16=i == 6)
        enddo
        call vol%kill
    end subroutine make_maps

    !> Construct a shifted-origin series whose last interval overlaps only its immediate predecessor.
    subroutine make_series(series)
        class(dose_series_manifest), intent(inout) :: series
        call series%new(BOX,SMPD,200._dp,DPF,10,5,5,string('0123456789ABCDEF'))
        call series%add_window(4,6,3.75_dp,7.5_dp,string('half1.mrc'),string('half2.mrc'))
        call series%add_window(7,9,7.5_dp,11.25_dp,string('half3.mrc'),string('half4.mrc'))
        call series%add_window(8,10,8.75_dp,12.5_dp,string('half5.mrc'),string('half6.mrc'))
    end subroutine make_series

    !> Round-trip physical doses, select the final overlap and relocate the complete input directory.
    subroutine test_roundtrip()
        type(dose_series_manifest) :: original, loaded
        type(string) :: message, even, odd
        integer :: box, a, b
        real(dp) :: smpd, kv, lo, hi
        logical :: ok
        call make_series(original)
        ok = original%write(string('data/dose_series.txt'),message)
        call assert_true(ok,'manifest write: '//message%to_char())
        ok = loaded%read(string('data/dose_series.txt'),message)
        call assert_true(ok,'manifest read: '//message%to_char())
        if( .not. ok ) return
        call assert_int(3,loaded%get_nwindows(),'publication keeps every actual interval')
        call loaded%get_geometry(box,smpd,kv)
        call assert_int(BOX,box,'manifest box')
        call assert_true(abs(smpd-SMPD) < 1.e-12_dp,'manifest sampling round trip')
        call assert_true(abs(kv-200._dp) < 1.e-12_dp,'manifest voltage remains physical')
        call loaded%get_window(1,a,b,lo,hi,even,odd)
        call assert_int(4,a,'first original fraction is retained')
        call assert_true(abs(lo-3.75_dp) < 1.e-12_dp,'dose before fromf is retained')
        call assert_true(abs(hi-7.5_dp) < 1.e-12_dp,'group width multiplies saved decoded-fraction dose once')
        ok = loaded%select_windows(message)
        call assert_true(ok,'consumer selects disjoint windows')
        call assert_int(2,loaded%get_nwindows(),'final predecessor is excluded')
        call loaded%get_window(2,a,b,lo,hi,even,odd)
        call assert_int(8,a,'full-width last interval survives')
        call assert_int(10,b,'last fraction survives')
        call simple_rename('data','moved data')
        ok = loaded%read(string('moved data/dose_series.txt'),message)
        call assert_true(ok,'relative paths survive relocation into a directory with spaces: '//message%to_char())
        if( ok )then
            ok = loaded%select_windows(message)
            call assert_true(ok,'selection accepts resolved paths containing a parent-directory space')
            call loaded%get_window(1,a,b,lo,hi,even,odd)
            call assert_true(file_exists(even) .and. file_exists(odd),'relocated paths resolve to existing halves')
        endif
        call simple_rename('moved data','data')
        call original%kill
        call loaded%kill
    end subroutine test_roundtrip

    !> Exercise parse failures and scientific invariants through recoverable reader errors.
    subroutine test_bad_metadata()
        type(dose_series_manifest) :: series
        type(string) :: message
        integer :: which
        logical :: ok
        do which = 1,16
            call write_case(which)
            ok = series%read(string('data/case.txt'),message)
            call assert_false(ok,'invalid manifest case')
            call assert_true(message%strlen() > 0,'invalid manifest explains rejection')
            call assert_int(0,series%get_nwindows(),'failed read leaves no partial series')
        enddo
        call series%kill
    end subroutine test_bad_metadata

    !> Supply malformed versions, metadata, doses, overlap, paths and missing files independently of the writer.
    subroutine write_case(which)
        integer, intent(in) :: which
        integer :: unit
        open(newunit=unit,file='data/case.txt',status='replace')
        if( which == 1 )then
            write(unit,'(A)') '# SIMPLE dose series 10'
        else
            write(unit,'(A)') '# SIMPLE dose series 1'
        endif
        write(unit,'(A)') '# program refine_motion_model'
        if( which == 2 )then
            write(unit,'(A)') '# state 2'
        else
            write(unit,'(A)') '# state 1'
        endif
        write(unit,'(A)') '# box 8'
        if( which == 3 )then
            write(unit,'(A)') '# smpd NaN'
        else
            write(unit,'(A)') '# smpd 1.5'
        endif
        if( which /= 4 ) write(unit,'(A)') '# kv 200'
        if( which == 5 )then
            write(unit,'(A)') '# nptcls 11'
        else
            write(unit,'(A)') '# nptcls 10'
        endif
        write(unit,'(A)') '# neven 5'
        write(unit,'(A)') '# nodd 5'
        if( which == 6 )then
            write(unit,'(A)') '# eo_digest not-a-digest'
        else
            write(unit,'(A)') '# eo_digest 0123456789abcdef'
        endif
        if( which == 15 )then
            write(unit,'(A)') '# dose_per_fraction 1.0E308'
        else
            write(unit,'(A)') '# dose_per_fraction 1.25'
        endif
        if( which == 7 ) write(unit,'(A)') '# box 8'
        select case(which)
            case(8)
                write(unit,'(A)') '4 6 0.0 7.5 half1.mrc half2.mrc'
            case(9)
                write(unit,'(A)') '4 6 3.75 7.5 half1.mrc half1.mrc'
            case(10)
                write(unit,'(A)') '4 6 3.75 7.5 missing.mrc half2.mrc'
            case(11)
                write(unit,'(A)') '4 6 3.75 7.5 half1.mrc half2.mrc unexpected'
            case(12)
                write(unit,'(A)') '4 6 3.75 NaN half1.mrc half2.mrc'
            case(13)
                write(unit,'(A)') '4 6 3.75 7.5 '//repeat('a',4200)//' half2.mrc'
            case(16)
                write(unit,'(A)') '4 6 3.75 7.5 volume_dir half2.mrc'
            case default
                write(unit,'(A)') '4 6 3.75 7.5 half1.mrc half2.mrc'
        end select
        if( which == 14 )then
            write(unit,'(A)') '6 8 6.25 10.0 half3.mrc half4.mrc'
        else
            write(unit,'(A)') '7 9 7.5 11.25 half3.mrc half4.mrc'
        endif
        close(unit)
    end subroutine write_case

    !> Modify actual header words and payload size, ensuring preallocated image geometry cannot hide errors.
    subroutine test_bad_headers()
        type(dose_series_manifest) :: series
        type(string) :: message
        integer(int8) :: original(1024)
        integer :: unit, which
        logical :: ok
        real :: nan
        open(newunit=unit,file='data/half1.mrc',status='old',access='stream',form='unformatted')
        read(unit,pos=1) original
        close(unit)
        nan = ieee_value(0.,ieee_quiet_nan)
        call make_series(series)
        do which = 1,8
            open(newunit=unit,file='data/half1.mrc',status='old',access='stream',form='unformatted')
            write(unit,pos=1) original
            select case(which)
                case(1); write(unit,pos=5) BOX+2
                case(2); write(unit,pos=13) 4
                case(3); write(unit,pos=45) 2.*real(BOX)*real(SMPD)
                case(4); write(unit,pos=65) 2
                case(5); write(unit,pos=197) 1.
                case(6); write(unit,pos=197) nan
                case(7); write(unit,pos=93) -1
                case(8); write(unit,pos=41) nan
            end select
            close(unit)
            ok = series%write(string('data/bad_header.txt'),message)
            call assert_false(ok,'actual MRC header case')
            call assert_false(file_exists('data/bad_header.txt'),'invalid header never publishes a manifest')
        enddo
        open(newunit=unit,file='data/half1.mrc',status='old',access='stream',form='unformatted')
        write(unit,pos=1) original
        close(unit)
        call simple_rename('data/half1.mrc','data/saved_half1.mrc')
        open(newunit=unit,file='data/half1.mrc',status='new',access='stream',form='unformatted')
        write(unit) original(:100)
        close(unit)
        ok = series%write(string('data/bad_header.txt'),message)
        call assert_false(ok,'truncated MRC header returns an error')
        open(newunit=unit,file='data/half1.mrc',status='replace',access='stream',form='unformatted')
        write(unit) original
        close(unit)
        ok = series%write(string('data/bad_header.txt'),message)
        call assert_false(ok,'truncated MRC payload returns an error')
        call del_file('data/half1.mrc')
        call simple_rename('data/saved_half1.mrc','data/half1.mrc')
        call series%kill
    end subroutine test_bad_headers

    !> Failed publication does not expose a partial series or overwrite somebody else's temporary file.
    subroutine test_publication()
        type(dose_series_manifest) :: series
        type(string) :: message
        integer :: unit
        logical :: ok
        character(len=32) :: content
        call make_series(series)
        open(newunit=unit,file='data/new.txt.tmp',status='new')
        write(unit,'(A)') 'owned by another attempt'
        close(unit)
        ok = series%write(string('data/new.txt'),message)
        call assert_false(ok,'an occupied temporary path prevents publication')
        call assert_false(file_exists('data/new.txt'),'no incomplete final manifest is created')
        open(newunit=unit,file='data/new.txt.tmp',status='old')
        read(unit,'(A)') content
        close(unit)
        call assert_true(trim(content) == 'owned by another attempt','foreign temporary contents are preserved')
        call series%kill
    end subroutine test_publication


    !> A complete last row without a terminating newline is still a complete manifest.
    subroutine test_final_line()
        type(dose_series_manifest) :: series
        type(string) :: message
        character(len=:), allocatable :: content
        integer :: unit, length
        logical :: ok
        call write_case(0)
        inquire(file='data/case.txt',size=length)
        allocate(character(len=length) :: content)
        open(newunit=unit,file='data/case.txt',access='stream',form='unformatted',status='old')
        read(unit) content
        close(unit)
        open(newunit=unit,file='data/no_newline.txt',access='stream',form='unformatted',status='new')
        write(unit) content(:length-1)
        close(unit)
        ok = series%read(string('data/no_newline.txt'),message)
        call assert_true(ok,'final row without newline: '//message%to_char())
        call assert_int(2,series%get_nwindows(),'unterminated final row is retained')
        call series%kill
    end subroutine test_final_line

end module simple_dose_series_manifest_tester
