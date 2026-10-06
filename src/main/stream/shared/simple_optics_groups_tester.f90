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
        call test_equal_shifts_one_group()
        call test_chain_is_one_group()
        call test_ids_kept_across_regroupings()
    end subroutine run_all_optics_groups_tests

    !> equal shifts (no dir_meta: all 0) form one group, quickly and without a distance matrix
    subroutine test_equal_shifts_one_group()
        integer, parameter :: NMANY = 50000
        type(sp_project) :: spproj
        real, allocatable :: xs(:)
        write(*,'(A)') 'test_equal_shifts_one_group'
        allocate(xs(NMANY), source=0.)
        call make_line(spproj, xs)
        call assign_optics_groups(spproj, TILT_THRES, .false., 0)
        call assert_int(1, spproj%os_optics%get_noris(), '50,000 equal shifts: one group')
        call assert_int(NMANY, nint(spproj%os_optics%get(1, 'pop')), 'holding every micrograph')
        call spproj%kill
    end subroutine test_equal_shifts_one_group

    !> single linkage: shifts each within the threshold of the next form one group, a gap wider
    !! than it another
    subroutine test_chain_is_one_group()
        type(sp_project) :: spproj
        integer          :: i
        write(*,'(A)') 'test_chain_is_one_group'
        call make_line(spproj, [0.0, 0.4, 0.8, 1.2, 3.0])
        call assign_optics_groups(spproj, TILT_THRES, .false., 0)
        call assert_int(2, spproj%os_optics%get_noris(), 'a chain and an outlier: two groups')
        call assert_true(all([(spproj%os_mic%get_int(i, 'ogid'), i=2,4)] == spproj%os_mic%get_int(1, 'ogid')),&
            &'the chain is one group')
        call assert_true(spproj%os_mic%get_int(5, 'ogid') /= spproj%os_mic%get_int(1, 'ogid'), 'the outlier is another')
        call spproj%kill
    end subroutine test_chain_is_one_group

    !> regrouped with last_ogid, a group keeps the id most of its micrographs had: adding a
    !! micrograph keeps both ids; when two groups merge the larger keeps its id; a new group takes
    !! the next id never given; os_optics rows follow the ids
    subroutine test_ids_kept_across_regroupings()
        type(sp_project) :: spproj
        integer :: last_ogid, id_a, id_b
        write(*,'(A)') 'test_ids_kept_across_regroupings'
        ! group A at 0 (two), group B at 0.9 (three): apart, as 0.9 exceeds the threshold
        call make_line(spproj, [0.0, 0.0, 0.9, 0.9, 0.9])
        last_ogid = 0
        call assign_optics_groups(spproj, TILT_THRES, .false., 0, last_ogid=last_ogid)
        call assert_int(2, spproj%os_optics%get_noris(), 'two groups')
        call assert_int(2, last_ogid, 'ids 1 and 2 given')
        id_a = spproj%os_mic%get_int(1, 'ogid')
        id_b = spproj%os_mic%get_int(3, 'ogid')
        ! one more micrograph in B
        call add_mic(spproj, 0.9)
        call assign_optics_groups(spproj, TILT_THRES, .false., 0, last_ogid=last_ogid)
        call assert_int(id_a, spproj%os_mic%get_int(1, 'ogid'), 'A keeps its id')
        call assert_int(id_b, spproj%os_mic%get_int(3, 'ogid'), 'B keeps its id')
        call assert_int(id_b, spproj%os_mic%get_int(6, 'ogid'), 'the new micrograph joins B')
        ! a bridge at 0.45 joins A and B: the merged group keeps B's id, B being larger
        call add_mic(spproj, 0.45)
        call assign_optics_groups(spproj, TILT_THRES, .false., 0, last_ogid=last_ogid)
        call assert_int(1, spproj%os_optics%get_noris(), 'merged: one group')
        call assert_int(id_b, spproj%os_mic%get_int(1, 'ogid'), 'the merged group keeps the larger group''s id')
        call assert_int(id_b, spproj%os_optics%get_int(1, 'ogid'), 'and its row has it')
        ! far away: a new group with the next id never given
        call add_mic(spproj, 20.0)
        call assign_optics_groups(spproj, TILT_THRES, .false., 0, last_ogid=last_ogid)
        call assert_int(2, spproj%os_optics%get_noris(), 'a new group')
        call assert_int(3, spproj%os_mic%get_int(8, 'ogid'), 'with id 3, not the merged-away id')
        call assert_int(3, last_ogid, 'the highest id given')
        call assert_true(spproj%os_optics%get_int(1, 'ogid') < spproj%os_optics%get_int(2, 'ogid'), 'rows in id order')
        call spproj%kill
    end subroutine test_ids_kept_across_regroupings

    ! micrographs with shifts (@p xs, 0), all accepted, one tilt group
    subroutine make_line( spproj, xs )
        type(sp_project), intent(inout) :: spproj
        real,             intent(in)    :: xs(:)
        integer :: imic
        call spproj%os_mic%new(size(xs), is_ptcl=.false.)
        do imic = 1,size(xs)
            call set_mic(spproj, imic, xs(imic))
        enddo
    end subroutine make_line

    ! one more micrograph at shift (@p x, 0), without a group
    subroutine add_mic( spproj, x )
        type(sp_project), intent(inout) :: spproj
        real,             intent(in)    :: x
        integer :: n
        n = spproj%os_mic%get_noris()
        call spproj%os_mic%reallocate(n + 1)
        call set_mic(spproj, n + 1, x)
    end subroutine add_mic

    subroutine set_mic( spproj, imic, x )
        type(sp_project), intent(inout) :: spproj
        integer,          intent(in)    :: imic
        real,             intent(in)    :: x
        call spproj%os_mic%set_state(imic, 1)
        call spproj%os_mic%set(imic, 'smpd',    SMPD)
        call spproj%os_mic%set(imic, 'cs',      CS)
        call spproj%os_mic%set(imic, 'kv',      KV)
        call spproj%os_mic%set(imic, 'fraca',   FRACA)
        call spproj%os_mic%set(imic, 'shiftx',  x)
        call spproj%os_mic%set(imic, 'shifty',  0.0)
        call spproj%os_mic%set(imic, 'tiltgrp', 1.0)
    end subroutine set_mic


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
