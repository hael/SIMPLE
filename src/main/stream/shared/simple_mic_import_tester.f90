!@descr: unit tests for appending upstream micrograph segments (simple_mic_import)
module simple_mic_import_tester
use simple_test_utils
use simple_string,       only: string
use simple_string_utils, only: int2str
use simple_fileio,       only: del_file
use simple_syslib,       only: get_process_id
use simple_oris,         only: oris
use simple_sp_project,   only: sp_project
use simple_mic_import,   only: append_mics_from_projects
implicit none
private
public :: run_all_mic_import_tests

contains

    subroutine run_all_mic_import_tests()
        write(*,'(A)') '**** running all micrograph import tests ****'
        call test_append_all_in_file_order()
        call test_append_accepted_only()
        call test_append_to_existing_and_empty_list()
    end subroutine run_all_mic_import_tests

    subroutine test_append_all_in_file_order()
        type(oris)   :: os_mic
        type(string) :: fnames(2)
        integer      :: nappended
        write(*,'(A)') 'test_append_all_in_file_order'
        call write_fixtures(fnames)
        call append_mics_from_projects(os_mic, fnames, .false., nappended)
        call assert_int(6, nappended,           'every micrograph of both files is appended')
        call assert_int(6, os_mic%get_noris(),  'the segment holds six micrographs')
        if( os_mic%get_noris() == 6 )then
            call assert_true(all(nint(os_mic%get_all('importind')) == [1,2,3,4,5,6]), 'file order is kept')
            call assert_int(0, os_mic%get_state(2), 'a rejected micrograph keeps its state')
        endif
        call os_mic%kill
        call delete_fixtures(fnames)
    end subroutine test_append_all_in_file_order

    subroutine test_append_accepted_only()
        type(oris)   :: os_mic
        type(string) :: fnames(2)
        integer      :: nappended
        write(*,'(A)') 'test_append_accepted_only'
        call write_fixtures(fnames)
        call append_mics_from_projects(os_mic, fnames, .true., nappended)
        call assert_int(4, nappended, 'only the accepted micrographs are appended')
        if( os_mic%get_noris() == 4 )then
            call assert_true(all(nint(os_mic%get_all('importind')) == [1,3,4,5]), 'the accepted ones, in file order')
        endif
        call os_mic%kill
        call delete_fixtures(fnames)
    end subroutine test_append_accepted_only

    subroutine test_append_to_existing_and_empty_list()
        type(oris)                :: os_mic
        type(string)              :: fnames(2)
        type(string), allocatable :: none(:)
        integer                   :: nappended
        write(*,'(A)') 'test_append_to_existing_and_empty_list'
        call write_fixtures(fnames)
        call append_mics_from_projects(os_mic, fnames(1:1), .false., nappended)
        call append_mics_from_projects(os_mic, fnames(2:2), .false., nappended)
        call assert_int(3, nappended,          'the second call appends the second file')
        call assert_int(6, os_mic%get_noris(), 'the segment grows across calls')
        if( os_mic%get_noris() == 6 )then
            call assert_true(all(nint(os_mic%get_all('importind')) == [1,2,3,4,5,6]), 'earlier micrographs stay first')
        endif
        allocate(none(0))
        call append_mics_from_projects(os_mic, none, .false., nappended)
        call assert_int(0, nappended,          'an empty file list appends nothing')
        call assert_int(6, os_mic%get_noris(), 'an empty file list leaves the segment alone')
        call os_mic%kill
        call delete_fixtures(fnames)
    end subroutine test_append_to_existing_and_empty_list

    ! two projects of three micrographs, import indices 1-3 and 4-6; micrographs 2 and 6 rejected
    subroutine write_fixtures( fnames )
        type(string), intent(inout) :: fnames(2)
        integer, parameter :: STATES(6) = [1, 0, 1, 1, 1, 0]
        type(sp_project) :: proj
        integer          :: ifile, imic, iglob
        do ifile = 1,2
            fnames(ifile) = 'mic_import_test_'//int2str(get_process_id())//'_'//int2str(ifile)//'.simple'
            call proj%os_mic%new(3, is_ptcl=.false.)
            do imic = 1,3
                iglob = (ifile - 1) * 3 + imic
                call proj%os_mic%set(imic, 'importind', real(iglob))
                call proj%os_mic%set_state(imic, STATES(iglob))
            enddo
            call proj%update_projinfo(fnames(ifile))
            call proj%write(fnames(ifile))
            call proj%kill
        enddo
    end subroutine write_fixtures

    subroutine delete_fixtures( fnames )
        type(string), intent(in) :: fnames(2)
        call del_file(fnames(1))
        call del_file(fnames(2))
    end subroutine delete_fixtures

end module simple_mic_import_tester
