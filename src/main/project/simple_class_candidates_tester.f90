!@descr: unit tests of class_candidates: the round trip of a published artifact, its merge from a distributed run's part files, and its refusals (none registered, another class count, another particle layout, a missing or truncated file)
submodule (simple_class_candidates) simple_class_candidates_tester
use simple_test_utils
use simple_syslib,       only: simple_getcwd, del_file
use simple_string_utils, only: int2str
use simple_fileio,       only: files_identical
implicit none

integer,          parameter :: NP    = 12
integer,          parameter :: NCLS  = 7
integer,          parameter :: KMAX  = 4
character(len=*), parameter :: PROJ  = 'class_candidates_test.simple'
character(len=*), parameter :: FNAME = 'class_candidates_test_cands.bin'

contains

    module subroutine run_all_class_candidates_tests()
        type(string) :: cwd_saved, root
        integer      :: nfail0
        write(*,'(A)') '**** running all class candidates tests ****'
        nfail0 = tests_failed
        call enter_fixture('class_candidates', cwd_saved, root)
        call test_round_trip()
        call test_parts()
        call test_refusals()
        call leave_fixture(cwd_saved, root, nfail0)
    end subroutine run_all_class_candidates_tests

    !> rows 3..10 predicted: the even ones seeded with candidates (row mod NCLS + 1, then the next classes), the odd
    !! ones dense with one candidate
    subroutine prediction_table( pinds, routes, ncands, cands )
        integer, allocatable, intent(out) :: pinds(:), routes(:), ncands(:), cands(:,:)
        integer :: i, j, n
        n = 8
        allocate(pinds(n), routes(n), ncands(n), cands(KMAX,n), source=0)
        do i = 1,n
            pinds(i) = i + 2
            if( mod(pinds(i),2) == 0 )then
                routes(i) = 1
                ncands(i) = 1 + mod(i, KMAX)
            else
                routes(i) = 2
                ncands(i) = 1
            endif
            do j = 1,ncands(i)
                cands(j,i) = 1 + mod(pinds(i) + j - 1, NCLS)
            enddo
        enddo
    end subroutine prediction_table

    subroutine test_round_trip()
        class(sp_project), allocatable :: project
        type(class_candidates) :: pub, back
        integer, allocatable   :: pinds(:), routes(:), ncands(:), cands(:,:)
        character(len=:), allocatable :: msg
        integer :: got(KMAX), n, i, status
        logical :: l_same, l_unseeded
        write(*,'(A)') 'test_round_trip'
        allocate(project)
        call make_project(project, 'cands_particles')
        call prediction_table(pinds, routes, ncands, cands)
        call pub%publish(project, project%os_ptcl2D, string(PROJ), string(FNAME), pinds, routes, ncands, cands, NCLS)
        call project%kill
        call project%read(string(PROJ))
        call back%new(project, project%os_ptcl2D, NCLS, status, msg)
        call assert_int(CANDS_OK, status, 'the published artifact opens through the project''s out segment')
        l_same = .true.
        do i = 1,size(pinds)
            if( back%is_seeded(pinds(i)) .neqv. (routes(i) == 1) ) l_same = .false.
            if( routes(i) /= 1 ) cycle
            call back%get_cands(pinds(i), n, got)
            if( n /= ncands(i) ) l_same = .false.
            if( n == ncands(i) )then
                if( any(got(1:n) /= cands(1:n,i)) ) l_same = .false.
            endif
        enddo
        call assert_true(l_same, 'the seeded rows and their candidates read back in order')
        l_unseeded = .not. (back%is_seeded(1) .or. back%is_seeded(2) .or. back%is_seeded(NP) .or. back%is_seeded(0) &
            &.or. back%is_seeded(NP + 1))
        call assert_true(l_unseeded, 'unpredicted rows and rows outside the particles are not seeded')
        call assert_int(count(routes == 1), back%get_nseeded([(i, i = 1,NP)]), 'the seeded count')
        call pub%kill
        call back%kill
        call project%kill
    end subroutine test_round_trip

    !> three parts (rows 1-2 with no predicted row, 3-7, 8-NP), each writing its range, merge into the bytes the
    !! whole prediction publishes
    subroutine test_parts()
        integer, parameter :: NPARTS = 3
        integer, parameter :: RANGES(2,NPARTS) = reshape([1,2, 3,7, 8,NP], [2,NPARTS])
        character(len=*), parameter :: FNAME_PARTS = 'class_candidates_test_cands_parts.bin'
        class(sp_project), allocatable :: project
        type(class_candidates) :: pub
        type(string)           :: part_fnames(NPARTS)
        integer, allocatable   :: pinds(:), routes(:), ncands(:), cands(:,:)
        logical, allocatable   :: in_part(:)
        integer, allocatable   :: cols(:)
        integer :: ipart, i
        write(*,'(A)') 'test_parts'
        allocate(project)
        call make_project(project, 'cands_particles')
        call prediction_table(pinds, routes, ncands, cands)
        call pub%publish(project, project%os_ptcl2D, string(PROJ), string(FNAME), pinds, routes, ncands, cands, NCLS)
        do ipart = 1,NPARTS
            part_fnames(ipart) = 'class_candidates_test_part'//int2str(ipart)//'.bin'
            in_part = pinds >= RANGES(1,ipart) .and. pinds <= RANGES(2,ipart)
            cols    = pack([(i, i = 1,size(pinds))], in_part)
            call pub%write_part(NP, RANGES(1,ipart), RANGES(2,ipart), part_fnames(ipart), pinds(cols), routes(cols), &
                &ncands(cols), cands(:,cols), NCLS)
        enddo
        call pub%publish_parts(project, project%os_ptcl2D, string(PROJ), string(FNAME_PARTS), part_fnames, NCLS)
        call assert_true(files_identical(string(FNAME), string(FNAME_PARTS)), 'the parts merge into the whole prediction''s artifact, byte for byte')
        call pub%kill
        call project%kill
        do ipart = 1,NPARTS
            call del_file(part_fnames(ipart))
        enddo
        call del_file(FNAME_PARTS)
    end subroutine test_parts

    subroutine test_refusals()
        class(sp_project), allocatable :: project, other
        type(class_candidates) :: pub, back
        integer, allocatable   :: pinds(:), routes(:), ncands(:), cands(:,:)
        character(len=:), allocatable :: msg
        integer :: status, u
        write(*,'(A)') 'test_refusals'
        allocate(project, other)
        ! no artifact registered
        call make_project(project, 'cands_particles')
        call back%new(project, project%os_ptcl2D, NCLS, status, msg)
        call assert_int(CANDS_NONE, status, 'a project without the artifact: none')
        call prediction_table(pinds, routes, ncands, cands)
        call pub%publish(project, project%os_ptcl2D, string(PROJ), string(FNAME), pinds, routes, ncands, cands, NCLS)
        ! another class count
        call back%new(project, project%os_ptcl2D, NCLS + 1, status, msg)
        call assert_int(CANDS_STALE, status, 'another class count: stale')
        ! another particle layout: the same rows from another stack, the artifact's out entry carried over
        call make_project(other, 'cands_other_particles')
        other%os_out = project%os_out
        call back%new(other, other%os_ptcl2D, NCLS, status, msg)
        call assert_int(CANDS_STALE, status, 'another particle layout (the digest): stale')
        call assert_true(.not. back%is_seeded(4), 'a stale artifact seeds nobody')
        ! a truncated file
        open(newunit=u, file=FNAME, access='stream', form='unformatted', status='replace', action='write')
        write(u) 'CLSCAND1'
        close(u)
        call back%new(project, project%os_ptcl2D, NCLS, status, msg)
        call assert_int(CANDS_UNREADABLE, status, 'a truncated file: unreadable')
        ! a missing file
        call del_file(FNAME)
        call back%new(project, project%os_ptcl2D, NCLS, status, msg)
        call assert_int(CANDS_UNREADABLE, status, 'a missing file: unreadable')
        call pub%kill
        call back%kill
        call project%kill
        call other%kill
    end subroutine test_refusals

    !> NP particle rows of one stack named stkname, all selected, in both particle fields
    subroutine make_project( project, stkname )
        class(sp_project), intent(inout) :: project
        character(len=*),  intent(in)    :: stkname
        type(string) :: cwd
        integer :: i
        call project%kill
        call project%projinfo%new(1, is_ptcl=.false.)
        call project%projinfo%set(1, 'projname', PROJ(1:len(PROJ)-7))
        call project%projinfo%set(1, 'projfile', PROJ)
        call simple_getcwd(cwd)
        call project%projinfo%set(1, 'cwd', cwd%to_char())
        call project%os_stk%new(1, is_ptcl=.false.)
        call project%os_stk%set(1, 'stk',        '/data/'//stkname//'.mrcs')
        call project%os_stk%set(1, 'fromp',      1)
        call project%os_stk%set(1, 'top',        NP)
        call project%os_stk%set(1, 'nptcls_stk', NP)
        call project%os_stk%set(1, 'box',        16)
        call project%os_stk%set(1, 'smpd',       1.5)
        call project%os_ptcl2D%new(NP, is_ptcl=.true.)
        do i = 1,NP
            call project%os_ptcl2D%set(i, 'stkind', 1)
            call project%os_ptcl2D%set(i, 'indstk', i)
            call project%os_ptcl2D%set_state(i, 1)
        enddo
        project%os_ptcl3D = project%os_ptcl2D
        call project%write(string(PROJ))
    end subroutine make_project

end submodule simple_class_candidates_tester
