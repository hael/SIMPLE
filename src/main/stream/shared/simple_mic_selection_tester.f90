!@descr: unit tests for the micrograph threshold rejection rule (simple_mic_selection)
module simple_mic_selection_tester
use simple_test_utils
use simple_string,        only: string
use simple_string_utils,  only: int2str
use simple_fileio,        only: del_file, simple_touch
use simple_syslib,        only: get_process_id
use simple_oris,          only: oris
use simple_mic_selection, only: reject_mics_by_thresholds, reject_mics_without_particles
implicit none
private
public :: run_all_mic_selection_tests

real, parameter :: CTFRES_THRES  = 10.0
real, parameter :: ICEFRAC_THRES = 1.0
real, parameter :: ASTIG_THRES   = 10.0

contains

    subroutine run_all_mic_selection_tests()
        write(*,'(A)') '**** running all micrograph threshold selection tests ****'
        call test_each_threshold_rejects()
        call test_absent_thresholds_and_missing_keys()
        call test_rejected_mics_stay_rejected_and_uncounted()
        call test_reject_mics_without_particles()
    end subroutine run_all_mic_selection_tests

    !> 1 has particles and a box file; 2 has no particles; 3 no particle count; 4 a missing box
    !! file; 5 no box file entry; 6 is already rejected
    subroutine test_reject_mics_without_particles()
        type(oris)   :: os
        type(string) :: boxfile
        integer      :: nrejected, imic
        write(*,'(A)') 'test_reject_mics_without_particles'
        boxfile = 'mic_selection_test_'//int2str(get_process_id())//'.box'
        call simple_touch(boxfile%to_char())
        call os%new(6, is_ptcl=.false.)
        do imic = 1,6
            call os%set_state(imic, 1)
            call os%set(imic, 'nptcls',  3.)
            call os%set(imic, 'boxfile', boxfile)
        enddo
        call os%set(2, 'nptcls', 0.)
        call os%delete_entry(3, 'nptcls')
        call os%set(4, 'boxfile', string('mic_selection_test_missing.box'))
        call os%delete_entry(5, 'boxfile')
        call os%set_state(6, 0)
        call reject_mics_without_particles(os, nrejected)
        call assert_int(4, nrejected,       'four micrographs lack particles or a box file')
        call assert_int(1, os%get_state(1), 'a micrograph with particles and a box file is kept')
        do imic = 2,5
            call assert_int(0, os%get_state(imic), 'micrograph '//int2str(imic)//' is rejected')
        enddo
        call assert_int(0, os%get_state(6), 'an already rejected micrograph stays rejected')
        call os%kill
        call del_file(boxfile)
    end subroutine test_reject_mics_without_particles

    subroutine test_each_threshold_rejects()
        type(oris) :: os
        integer    :: nrejected
        write(*,'(A)') 'test_each_threshold_rejects'
        call make_mics(os)
        call reject_mics_by_thresholds(os, nrejected, ctfres=CTFRES_THRES, icefrac=ICEFRAC_THRES, astig=ASTIG_THRES)
        call assert_int(3, nrejected,        'three micrographs exceed one threshold each')
        call assert_int(1, os%get_state(1),  'a micrograph under every threshold is kept')
        call assert_int(0, os%get_state(2),  'ctfres above its threshold rejects')
        call assert_int(0, os%get_state(3),  'icefrac above its threshold rejects')
        call assert_int(0, os%get_state(4),  'astig above its threshold rejects')
        call assert_int(1, os%get_state(5),  'a value equal to its threshold is kept')
        call os%kill
    end subroutine test_each_threshold_rejects

    subroutine test_absent_thresholds_and_missing_keys()
        type(oris) :: os
        integer    :: nrejected
        write(*,'(A)') 'test_absent_thresholds_and_missing_keys'
        call make_mics(os)
        call reject_mics_by_thresholds(os, nrejected, ctfres=CTFRES_THRES)
        call assert_int(1, nrejected,       'only the supplied threshold is applied')
        call assert_int(1, os%get_state(3), 'icefrac is not applied when not supplied')
        call assert_int(1, os%get_state(4), 'astig is not applied when not supplied')
        call os%kill
        ! a micrograph without the key is not rejected on it
        call os%new(1, is_ptcl=.false.)
        call os%set_state(1, 1)
        call os%set(1, 'icefrac', 0.2)
        call reject_mics_by_thresholds(os, nrejected, ctfres=1.0)
        call assert_int(0, nrejected,       'a missing key never rejects')
        call assert_int(1, os%get_state(1), 'a micrograph without ctfres is kept')
        call os%kill
    end subroutine test_absent_thresholds_and_missing_keys

    subroutine test_rejected_mics_stay_rejected_and_uncounted()
        type(oris) :: os
        integer    :: nrejected
        write(*,'(A)') 'test_rejected_mics_stay_rejected_and_uncounted'
        call make_mics(os)
        call os%set_state(2, 0)
        call reject_mics_by_thresholds(os, nrejected, ctfres=CTFRES_THRES, icefrac=ICEFRAC_THRES, astig=ASTIG_THRES)
        call assert_int(2, nrejected,       'an already rejected micrograph is not counted again')
        call assert_int(0, os%get_state(2), 'an already rejected micrograph stays rejected')
        call reject_mics_by_thresholds(os, nrejected, ctfres=CTFRES_THRES, icefrac=ICEFRAC_THRES, astig=ASTIG_THRES)
        call assert_int(0, nrejected,       'a second pass rejects nothing new')
        call os%kill
    end subroutine test_rejected_mics_stay_rejected_and_uncounted

    ! five accepted micrographs: 1 passes every threshold, 2 fails ctfres, 3 fails icefrac,
    ! 4 fails astig, and 5 sits exactly on every threshold
    subroutine make_mics( os )
        type(oris), intent(inout) :: os
        integer :: imic
        call os%new(5, is_ptcl=.false.)
        do imic = 1,5
            call os%set_state(imic, 1)
            call os%set(imic, 'ctfres',  5.0)
            call os%set(imic, 'icefrac', 0.5)
            call os%set(imic, 'astig',   1.0)
        enddo
        call os%set(2, 'ctfres',  12.0)
        call os%set(3, 'icefrac', 1.5)
        call os%set(4, 'astig',   20.0)
        call os%set(5, 'ctfres',  CTFRES_THRES)
        call os%set(5, 'icefrac', ICEFRAC_THRES)
        call os%set(5, 'astig',   ASTIG_THRES)
    end subroutine make_mics

end module simple_mic_selection_tester
