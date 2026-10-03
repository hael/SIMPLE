!@descr: unit tests for optics-group assignment from beam-image shifts (simple_optics_groups)
module simple_optics_groups_tester
use simple_test_utils
use simple_string,        only: string
use simple_sp_project,    only: sp_project
use simple_optics_groups, only: assign_optics_groups
implicit none
private
public :: run_all_optics_groups_tests

integer, parameter :: NMICS      = 5
real,    parameter :: TILT_THRES = 0.5
! two shift clusters, near (0,0) twice and near (5,5) three times, one per tilt group
real,    parameter :: SHIFT_X(NMICS)    = [0.0, 0.1, 5.0, 5.1, 4.9]
real,    parameter :: SHIFT_Y(NMICS)    = [0.0,-0.1, 5.0, 4.9, 5.1]
integer, parameter :: TILT_GROUP(NMICS) = [1, 1, 2, 2, 2]
real,    parameter :: SMPD = 1.3, CS = 2.7, KV = 300.0, FRACA = 0.1

contains

    subroutine run_all_optics_groups_tests()
        write(*,'(A)') '**** running all optics group assignment tests ****'
        call test_two_shift_clusters_with_beamtilt()
        call test_two_shift_clusters_without_beamtilt()
        call test_group_offset()
        call test_ctf_constants_from_first_accepted()
    end subroutine run_all_optics_groups_tests

    subroutine test_two_shift_clusters_with_beamtilt()
        type(sp_project) :: spproj
        write(*,'(A)') 'test_two_shift_clusters_with_beamtilt'
        call make_mics(spproj)
        call assign_optics_groups(spproj, TILT_THRES, .true., 0)
        call check_two_groups(spproj, 0, 'with beam tilt')
        call spproj%kill
    end subroutine test_two_shift_clusters_with_beamtilt

    !> without beam tilt all micrographs form one tilt group; the shifts still separate the clusters
    subroutine test_two_shift_clusters_without_beamtilt()
        type(sp_project) :: spproj
        write(*,'(A)') 'test_two_shift_clusters_without_beamtilt'
        call make_mics(spproj)
        call assign_optics_groups(spproj, TILT_THRES, .false., 0)
        call check_two_groups(spproj, 0, 'without beam tilt')
        call spproj%kill
    end subroutine test_two_shift_clusters_without_beamtilt

    subroutine test_group_offset()
        type(sp_project) :: spproj
        type(string)     :: ogname
        write(*,'(A)') 'test_group_offset'
        call make_mics(spproj)
        call assign_optics_groups(spproj, TILT_THRES, .true., 10)
        call check_two_groups(spproj, 10, 'with offset 10')
        if( spproj%os_optics%get_noris() == 2 )then
            call assert_int(11, spproj%os_optics%get_int(1, 'ogid'), 'offset: the first group id')
            call assert_int(12, spproj%os_optics%get_int(2, 'ogid'), 'offset: the second group id')
            ogname = spproj%os_optics%get_str(2, 'ogname')
            call assert_char('opticsgroup12', ogname%to_char(), 'offset: the group name carries the id')
        endif
        call spproj%kill
    end subroutine test_group_offset

    !> the optics rows take smpd from the first accepted micrograph, or from the first when none is
    subroutine test_ctf_constants_from_first_accepted()
        type(sp_project) :: spproj
        integer          :: imic
        write(*,'(A)') 'test_ctf_constants_from_first_accepted'
        call make_mics(spproj)
        call spproj%os_mic%set_state(1, 0)
        call spproj%os_mic%set(1, 'smpd', 9.9)
        call assign_optics_groups(spproj, TILT_THRES, .true., 0)
        call assert_real(SMPD, spproj%os_optics%get(1, 'smpd'), 1.e-6, 'smpd of the first accepted micrograph')
        call assert_int(2, spproj%os_optics%get_noris(), 'a rejected micrograph is still grouped')
        do imic = 1,NMICS
            call spproj%os_mic%set_state(imic, 0)
        enddo
        call assign_optics_groups(spproj, TILT_THRES, .true., 0)
        call assert_real(9.9, spproj%os_optics%get(1, 'smpd'), 1.e-6, 'smpd of the first micrograph when none is accepted')
        call spproj%kill
    end subroutine test_ctf_constants_from_first_accepted

    ! the expectations of the stream optics test (simple_stream_tester): two groups of 2 and 3
    ! micrographs with centroids (0.05,-0.05) and (5,5), ids offset+1 and offset+2
    subroutine check_two_groups( spproj, offset, label )
        type(sp_project), intent(inout) :: spproj
        integer,          intent(in)    :: offset
        character(len=*), intent(in)    :: label
        integer :: group_a, group_b
        logical :: l_groups
        call assert_int(2, spproj%os_optics%get_noris(), label//': two optics groups')
        if( spproj%os_optics%get_noris() /= 2 ) return
        group_a  = spproj%os_mic%get_int(1, 'ogid') - offset
        group_b  = spproj%os_mic%get_int(3, 'ogid') - offset
        l_groups = group_a /= group_b .and. all([group_a, group_b] >= 1) .and. all([group_a, group_b] <= 2)
        call assert_true(l_groups, label//': the two shift clusters are in different groups')
        if( .not. l_groups ) return
        call assert_int(group_a + offset, spproj%os_mic%get_int(2, 'ogid'), label//': micrograph 2 in the (0,0) group')
        call assert_int(group_b + offset, spproj%os_mic%get_int(4, 'ogid'), label//': micrograph 4 in the (5,5) group')
        call assert_int(group_b + offset, spproj%os_mic%get_int(5, 'ogid'), label//': micrograph 5 in the (5,5) group')
        call assert_int(2, nint(spproj%os_optics%get(group_a, 'pop')), label//': the (0,0) group holds two')
        call assert_int(3, nint(spproj%os_optics%get(group_b, 'pop')), label//': the (5,5) group holds three')
        call assert_real( 0.05, spproj%os_optics%get(group_a, 'opcx'), 0.01, label//': (0,0) centroid x')
        call assert_real(-0.05, spproj%os_optics%get(group_a, 'opcy'), 0.01, label//': (0,0) centroid y')
        call assert_real( 5.00, spproj%os_optics%get(group_b, 'opcx'), 0.01, label//': (5,5) centroid x')
        call assert_real( 5.00, spproj%os_optics%get(group_b, 'opcy'), 0.01, label//': (5,5) centroid y')
    end subroutine check_two_groups

    subroutine make_mics( spproj )
        type(sp_project), intent(inout) :: spproj
        integer :: imic
        call spproj%os_mic%new(NMICS, is_ptcl=.false.)
        do imic = 1,NMICS
            call spproj%os_mic%set_state(imic, 1)
            call spproj%os_mic%set(imic, 'smpd',    SMPD)
            call spproj%os_mic%set(imic, 'cs',      CS)
            call spproj%os_mic%set(imic, 'kv',      KV)
            call spproj%os_mic%set(imic, 'fraca',   FRACA)
            call spproj%os_mic%set(imic, 'shiftx',  SHIFT_X(imic))
            call spproj%os_mic%set(imic, 'shifty',  SHIFT_Y(imic))
            call spproj%os_mic%set(imic, 'tiltgrp', real(TILT_GROUP(imic)))
        enddo
    end subroutine make_mics

end module simple_optics_groups_tester
